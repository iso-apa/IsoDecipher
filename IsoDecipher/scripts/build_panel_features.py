#!/usr/bin/env python3
"""
IsoDecipher: Build feature panel for isoform quantification
-----------------------------------------------------------
From an Ensemble GTF + list of genes, collapse isoforms into
polyA groups (shared 3'ends) and emit features for transcript 
assignment from Cell Ranger BAMs.

Key design decisions:
 - retained_intron / NMD transcripts contribute APA site COORDINATES
   but are EXCLUDED from UTR length calculations (utr_source field)
 - utr_source field tags each group: protein_coding / mixed / non_coding
 - Downstream analysis should filter on utr_source for UTR trajectory

Features:
 - polyA_group: merge window around transcript ends
 - avg_spliced_utr: mean 3'UTR length (protein_coding transcripts only)
 - avg_genomic_utr: mean genomic UTR distance (protein_coding only)
 - utr_source: protein_coding | mixed | non_coding

Usage:
python IsoDecipher/scripts/build_panel_features.py \
    --gtf data/Homo_sapiens.GRCh38.115.gtf \
    --genes data/gene_list.txt \
    --out results/panel_features.csv \
    [--tolerance 10] \
    [--include-nmd] \
    [--include-retained-intron] \
    [--no-filter]
"""

import argparse
import gffutils
import pandas as pd
from collections import defaultdict
import os


# Biotypes that don't have CDS by definition — don't penalize them
NO_CDS_BIOTYPES = {
    "lncRNA",
    "processed_transcript",
    "non_stop_decay",
    "sense_intronic",
    "sense_overlapping",
    "antisense",
}

# Biotypes that contribute coordinates but NOT UTR length
NON_CODING_BIOTYPES = {
    "retained_intron",
    "nonsense_mediated_decay",
    "processed_transcript",
    "lncRNA",
    "non_stop_decay",
    "sense_intronic",
    "sense_overlapping",
    "antisense",
}

IG_WHITELIST = {"IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build cleavage-centered polyA panel (distance-based model)"
    )
    parser.add_argument("--gtf",    required=True)
    parser.add_argument("--genes",  required=True)
    parser.add_argument("--db",     help="Optional: specific path for gffutils DB")
    parser.add_argument("--out",    required=True)
    parser.add_argument("--custom_params", help="Custom tolerance TSV/CSV")
    parser.add_argument("--tolerance", type=int, default=10)
    parser.add_argument("--no-filter", action="store_true")
    parser.add_argument("--include-nmd", action="store_true",
                        help="Include NMD transcripts (coordinates only, no UTR)")
    parser.add_argument("--include-retained-intron", action="store_true",
                        help="Include retained intron transcripts (coordinates only, no UTR)")
    return parser.parse_args()


def load_custom_parameters(file):
    if not file:
        return {}
    df = pd.read_csv(file, sep=None, engine="python")
    if not {"gene", "end_tolerance"} <= set(df.columns):
        raise ValueError("Custom parameter must contain columns: gene, end_tolerance")
    return df.set_index("gene")["end_tolerance"].to_dict()


def load_or_build_db(gtf):
    db_path = gtf + ".db"
    build = False
    if os.path.exists(db_path):
        if os.path.getmtime(gtf) > os.path.getmtime(db_path):
            print("[IsoDecipher] GTF newer than DB. Rebuilding...")
            build = True
        else:
            print("[IsoDecipher] Using existing DB")
    else:
        print("[IsoDecipher] No DB found. Building...")
    if build:
        gffutils.create_db(
            gtf, dbfn=db_path, force=True, keep_order=True,
            merge_strategy="merge",
            disable_infer_genes=True, disable_infer_transcripts=True
        )
    return gffutils.FeatureDB(db_path, keep_order=True)


def collect_transcript_end(db, gene_list, skip_biotypes):
    """
    Collect strand-aware transcript 3' boundaries.

    Key change from v1:
    - retained_intron / NMD transcripts are NO LONGER skipped entirely
      (if --include-retained-intron / --include-nmd flags are set)
    - Instead, their coordinates are collected but
      spliced_utr_length and genomic_utr_length are set to None
    - biotype field is recorded for downstream utr_source tagging
    """
    gene_data   = defaultdict(list)
    all_genes   = {
        g.attributes.get("gene_name", [g.id])[0].upper(): g
        for g in db.features_of_type("gene")
    }
    missing_genes    = []
    skipped_biotype  = 0
    skipped_no_cds   = 0

    for gene_name in gene_list:
        gene = all_genes.get(gene_name)
        if gene is None:
            missing_genes.append(gene_name)
            continue

        chrom       = gene.seqid
        transcripts = db.children(gene, featuretype="transcript")

        for tx in transcripts:
            tx_id   = tx.attributes.get("transcript_id",   [""])[0]
            tx_name = tx.attributes.get("transcript_name", [""])[0]
            strand  = tx.strand
            biotype = tx.attributes.get("transcript_biotype", [""])[0]

            # Hard skip — truly useless biotypes
            if biotype in skip_biotypes:
                skipped_biotype += 1
                continue

            # CDS check — only for protein_coding
            cds = list(db.children(tx, featuretype="CDS"))
            if not cds and biotype not in NO_CDS_BIOTYPES:
                skipped_no_cds += 1
                continue

            exons = list(db.children(tx, featuretype="exon"))
            if exons:
                if strand == "+":
                    last_exon = max(exons, key=lambda e: e.end)
                    coord     = last_exon.end
                else:
                    last_exon = min(exons, key=lambda e: e.start)
                    coord     = last_exon.start
            else:
                coord = tx.end if strand == "+" else tx.start

            # UTR length — only computed for protein_coding
            # Non-coding biotypes get None → excluded from avg_spliced_utr
            is_protein_coding = (biotype == "protein_coding")

            genomic_utr_dist = None
            if is_protein_coding and cds:
                if strand == "+":
                    cds_end          = max(c.end for c in cds)
                    genomic_utr_dist = coord - cds_end
                else:
                    cds_start        = min(c.start for c in cds)
                    genomic_utr_dist = cds_start - coord

            spliced_utr_len = None
            if is_protein_coding:
                utr_features = list(db.children(tx, featuretype="three_prime_UTR"))
                if utr_features:
                    spliced_utr_len = sum(u.end - u.start + 1 for u in utr_features)
                elif genomic_utr_dist is not None:
                    spliced_utr_len = genomic_utr_dist

            gene_data[gene_name].append({
                "coord":              coord,
                "biotype":            biotype,
                "transcript_id":      tx_id,
                "transcript_name":    tx_name,
                "genomic_utr_length": genomic_utr_dist,  # None for non-coding
                "spliced_utr_length": spliced_utr_len,   # None for non-coding
                "chrom":              chrom,
                "strand":             strand,
            })

        strands = set(t['strand'] for t in gene_data[gene_name])
        if len(strands) > 1:
            print(f"[ERROR] Gene {gene_name} has transcripts on multiple strands: {strands}")

    if missing_genes:
        print(f"\n[WARN] The following {len(missing_genes)} genes were not found in the GTF:")
        print(f"       {', '.join(missing_genes)}")
        print(f"       Please check your gene list spelling or GTF version.\n")

    print(f"[FILTER] Skipped {skipped_biotype} transcripts by biotype")
    print(f"[FILTER] Skipped {skipped_no_cds} transcripts with no CDS annotation")

    return gene_data


def cluster_transcript_ends(transcript_data, gene_name, param_dict, default_tolerance=10):
    tolerance  = param_dict.get(gene_name.upper(), default_tolerance)
    sorted_tx  = sorted(transcript_data, key=lambda x: x['coord'])
    if not sorted_tx:
        return []

    clusters       = []
    current_cluster = [sorted_tx[0]]
    for i in range(1, len(sorted_tx)):
        dist = sorted_tx[i]["coord"] - sorted_tx[i-1]["coord"]
        if dist <= tolerance:
            current_cluster.append(sorted_tx[i])
        else:
            clusters.append(current_cluster)
            current_cluster = [sorted_tx[i]]
    clusters.append(current_cluster)

    if gene_name.upper() in param_dict:
        print(f"  [CLUSTER] {gene_name}: custom tolerance {tolerance}bp")

    return clusters


def get_utr_source(cluster):
    """
    Classify the UTR data source for a cluster of transcripts.

    Returns:
        'protein_coding' — all transcripts are protein_coding
        'mixed'          — mix of protein_coding + non-coding
        'non_coding'     — no protein_coding transcripts (retained_intron, NMD, etc.)
    """
    biotypes = set(tx['biotype'] for tx in cluster)
    if biotypes == {"protein_coding"}:
        return "protein_coding"
    elif "protein_coding" in biotypes:
        return "mixed"
    else:
        return "non_coding"


def filter_panel_features(df_panel):
    """
    Remove poorly annotated polyA groups.

    Rules:
    1. Drop groups where avg_spliced_utr == 0 AND num_transcripts == 1
       AND utr_source == 'protein_coding'
       (true zero-UTR singleton from protein_coding = likely annotation artifact)
    2. Keep non_coding singletons only if they're the dominant group
    3. Always keep the dominant group per gene
    """
    filtered_rows = []
    removed       = 0

    for gene, gdf in df_panel.groupby('gene'):
        dominant_idx = gdf['num_transcripts'].idxmax()

        for idx, row in gdf.iterrows():
            is_zero_utr      = row['avg_spliced_utr'] == 0
            is_singleton     = row['num_transcripts'] == 1
            is_dominant      = (idx == dominant_idx)
            is_protein_coding = row['utr_source'] == 'protein_coding'

            if is_dominant:
                filtered_rows.append(row)
            elif is_zero_utr and is_singleton and is_protein_coding:
                removed += 1  # artifact — drop
            else:
                filtered_rows.append(row)  # keep non_coding, mixed, multi-tx

    df_filtered = pd.DataFrame(filtered_rows).reset_index(drop=True)
    df_filtered['polyA_group'] = df_filtered.groupby('gene').cumcount()

    # Re-assign IG labels
    ig_mask = df_filtered['gene'].str.upper().isin(IG_WHITELIST)
    df_filtered.loc[ig_mask, 'user_label'] = df_filtered.loc[ig_mask, 'polyA_group'].apply(
        lambda i: "Secreted" if i == 0 else "Membrane"
    )

    print(f"\n[FILTER] Removed {removed} zero-UTR singleton groups (protein_coding only)")
    print(f"[FILTER] Kept {len(df_filtered)} groups ({len(df_panel)} before filter)")

    return df_filtered


def print_panel_summary(df):
    groups_per_gene = df.groupby('gene')['polyA_group'].count()
    print(f"\n[SUMMARY] Group count distribution:")
    print(f"  1 group  (no analysis): {(groups_per_gene == 1).sum()} genes")
    print(f"  2 groups (PUI):         {(groups_per_gene == 2).sum()} genes")
    print(f"  3+ groups (Entropy):    {(groups_per_gene >= 3).sum()} genes")

    print(f"\n[SUMMARY] UTR source distribution:")
    for src, count in df['utr_source'].value_counts().items():
        print(f"  {src:<20} {count:>6} groups")


def main():
    args = parse_args()

    # Hard skip biotypes — truly useless regardless of flags
    skip_biotypes = {"misc_RNA", "pseudogene"}

    # Retained intron — skip unless --include-retained-intron
    if not args.include_retained_intron:
        skip_biotypes.add("retained_intron")
        print("[INFO] Retained intron transcripts excluded (use --include-retained-intron to include)")
    else:
        print("[INFO] Retained intron transcripts INCLUDED (coordinates only, UTR=None)")

    # NMD — skip unless --include-nmd
    if not args.include_nmd:
        skip_biotypes.add("nonsense_mediated_decay")
        print("[INFO] NMD transcripts excluded (use --include-nmd to include)")
    else:
        print("[INFO] NMD transcripts INCLUDED (coordinates only, UTR=None)")

    db = load_or_build_db(args.gtf)

    if not os.path.exists(args.genes):
        raise FileNotFoundError(f"Gene list not found: {args.genes}")
    with open(args.genes, "r") as f:
        gene_list = [line.strip().upper()
                     for line in f
                     if line.strip() and not line.strip().startswith("#")]

    custom_params = load_custom_parameters(args.custom_params)

    print(f"[IsoDecipher] Collecting transcript ends for {len(gene_list)} genes...")
    print(f"[IsoDecipher] Global clustering tolerance: {args.tolerance}bp")

    raw_gene_data = collect_transcript_end(db, gene_list, skip_biotypes)

    panel_rows = []

    for gene_name in gene_list:
        transcripts = raw_gene_data.get(gene_name, [])
        if not transcripts:
            continue

        clusters = cluster_transcript_ends(
            transcripts, gene_name, custom_params,
            default_tolerance=args.tolerance
        )

        strand = transcripts[0]['strand']
        if strand == "+":
            clusters.sort(key=lambda c: min(tx['coord'] for tx in c))
        else:
            clusters.sort(key=lambda c: max(tx['coord'] for tx in c), reverse=True)

        for i, cluster in enumerate(clusters):
            coords    = [tx['coord'] for tx in cluster]
            rep_coord = int(sum(coords) / len(coords))

            # UTR — only from protein_coding transcripts
            spliced_utr_lengths = [tx['spliced_utr_length'] for tx in cluster
                                   if tx['spliced_utr_length'] is not None]
            avg_spliced_utr = (sum(spliced_utr_lengths) / len(spliced_utr_lengths)
                               if spliced_utr_lengths else 0)

            genomic_utr_lengths = [tx['genomic_utr_length'] for tx in cluster
                                   if tx['genomic_utr_length'] is not None]
            avg_genomic_utr = (sum(genomic_utr_lengths) / len(genomic_utr_lengths)
                               if genomic_utr_lengths else 0)

            # UTR source tag
            utr_source = get_utr_source(cluster)

            row = {
                "gene":             gene_name,
                "polyA_group":      i,
                "rep_coord":        rep_coord,
                "coord_min":        min(coords),
                "coord_max":        max(coords),
                "coord_spread":     max(coords) - min(coords),
                "strand":           strand,
                "chrom":            cluster[0]['chrom'],
                "avg_spliced_utr":  round(avg_spliced_utr, 2),
                "avg_genomic_utr":  round(avg_genomic_utr, 2),
                "utr_source":       utr_source,       # ← new
                "num_transcripts":  len(cluster),     # fixed typo: was num_transcirpts
                "transcript_ids":   ";".join([tx['transcript_id']   for tx in cluster]),
                "transcript_names": ";".join([tx['transcript_name'] for tx in cluster]),
            }

            if gene_name.upper() in IG_WHITELIST:
                row["user_label"] = "Secreted" if i == 0 else "Membrane"
            else:
                row["user_label"] = "N/A"

            panel_rows.append(row)

    df_panel = pd.DataFrame(panel_rows)

    if not args.no_filter:
        df_panel = filter_panel_features(df_panel)
    else:
        print("\n[FILTER] Skipping filter (--no-filter flag set)")

    print_panel_summary(df_panel)

    df_panel.to_csv(args.out, index=False)
    print(f"\n[SUCCESS] IsoDecipher Panel complete!")
    print(f" - Output file:    {args.out}")
    print(f" - Total features: {len(df_panel)}")
    print(f" - Genes processed: {df_panel['gene'].nunique()}")


if __name__ == "__main__":
    main()
