#!/usr/bin/env python3
# ==============================================================================
# Copyright (c) 2026 IsoMatrix Suite. All rights reserved.
# ==============================================================================
"""
IsoDecipher: build_panel_features.py
------------------------------------
Phase 2/2 of Panel Building: Dynamic Clustering & Analysis
Loads pre-built reference indices (.pkl) to perform:
  - PolyASite v2.0 gene assignment
  - Two-phase clustering (v2.0 first, then GTF absorption)
  - Auto-calibrated site_type via Exon Tree
  - Cross-group UTR imputation with Intron Shield

Usage:
  python build_panel_features.py \
      --indices-pkl reference/hg38_iso_indices.pkl \
      --polyasite   /path/to/polyasite2_hg38.bed \
      --out         results/panel/panel_features_global_v1.1.csv \
      --tolerance   150
"""

import argparse
import os
import pickle
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from intervaltree import IntervalTree
from _constants import (CODING_BIOTYPES, IG_BIOTYPES, TR_BIOTYPES,
                        IG_PSEUDO_BIOTYPES, TR_PSEUDO_BIOTYPES,
                        get_pas_motif_strength)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

INCLUDE_SITE_CLASSES = {"TE", "EX", "IN", "DS"}
DS_MAX_DIST     = 1000
MAX_BORROW_DIST = 3000


# ---------------------------------------------------------------------------
# Gene assignment
# ---------------------------------------------------------------------------

def _pick_best(hits, coord, strand):
    pc  = [g for g in hits if g["biotype"] in CODING_BIOTYPES]
    lnc = [g for g in hits if "lncRNA" in g["biotype"]]
    pool = pc or lnc or hits

    if len(pool) == 1:
        return pool[0]

    # Tiebreak: closest gene 3' end to coord, then longest gene
    def sort_key(g):
        gene_3end = g["end"] if strand == "+" else g["start"]
        return (abs(gene_3end - coord), -(g["end"] - g["start"]))

    return min(pool, key=sort_key)


def assign_gene(chrom, coord, strand, gene_index):
    key  = (chrom, strand)
    tree = gene_index.get(key)
    if not tree:
        return None, None, None

    # Step 1: gene body overlap (IntervalTree O(log N))
    overlaps = tree[coord]
    if overlaps:
        hits = [iv.data for iv in overlaps]
        best = _pick_best(hits, coord, strand)
        return best["name"], best["biotype"], "overlap"

    # Step 2: DS fallback — nearest gene 3' end within DS_MAX_DIST
    ds_cands = []
    for iv in tree:
        g         = iv.data
        gene_3end = g["end"] if strand == "+" else g["start"]
        dist      = (coord - gene_3end) if strand == "+" else (gene_3end - coord)
        if 0 < dist <= DS_MAX_DIST:
            ds_cands.append((dist, g))

    if ds_cands:
        pc   = [(d, g) for d, g in ds_cands if g["biotype"] in CODING_BIOTYPES]
        lnc  = [(d, g) for d, g in ds_cands if "lncRNA" in g["biotype"]]
        pool = pc or lnc or ds_cands
        best_dist, best_g = min(pool, key=lambda x: x[0])
        bio = ("protein_coding" if pool is pc else
               "non_coding"     if pool is lnc else best_g["biotype"])
        return best_g["name"], bio, "ds_fallback"

    return None, None, None


# ---------------------------------------------------------------------------
# Per-gene clustering (two-phase)
# ---------------------------------------------------------------------------

def build_gene_clusters(gene_name, strand, v2_sites, gtf_3ends, tolerance):
    groups = []

    # Phase A: cluster v2.0 sites (150bp consecutive gap)
    if v2_sites:
        v2_sorted = sorted(v2_sites, key=lambda x: x["coord"])
        v2_clusters = []
        current     = [v2_sorted[0]]
        for site in v2_sorted[1:]:
            if abs(site["coord"] - current[-1]["coord"]) <= tolerance:
                current.append(site)
            else:
                v2_clusters.append(current)
                current = [site]
        v2_clusters.append(current)

        for cluster in v2_clusters:
            # rep_coord = highest fraction, tie-break by usage
            rep = max(cluster, key=lambda x: (x["fraction"], x["usage"]))
            groups.append({
                "type":          "PA2",
                "rep_coord":     rep["coord"],
                "coord_min":     min(s["coord"] for s in cluster),
                "coord_max":     max(s["coord"] for s in cluster),
                "v2_sites":      cluster,
                "v2_fraction":   rep["fraction"],
                "v2_usage":      rep["usage"],
                "v2_site_class": rep["site_class"],
                "v2_signal":     rep.get("signal"),
                "num_polyAsites":len(cluster),
                "gtf_txs":       [],
            })

    # Phase B: absorb GTF 3' ends into NEAREST v2.0 cluster (nearest-first, not first-match)
    gtf_standalone = []
    for tx in gtf_3ends:
        coord     = tx["coord"]
        best_grp  = None
        best_dist = float('inf')
        for grp in groups:
            if (grp["coord_min"] - tolerance) <= coord <= (grp["coord_max"] + tolerance):
                dist = min(
                    abs(coord - grp["coord_min"]),
                    abs(coord - grp["coord_max"]),
                    abs(coord - grp["rep_coord"]),
                )
                if dist < best_dist:
                    best_dist = dist
                    best_grp  = grp
        if best_grp is not None:
            best_grp["gtf_txs"].append(tx)
            best_grp["type"]      = "GTF+PA2"
            best_grp["coord_min"] = min(best_grp["coord_min"], coord)
            best_grp["coord_max"] = max(best_grp["coord_max"], coord)
        else:
            gtf_standalone.append(tx)

    # Cluster standalone GTF 3' ends
    if gtf_standalone:
        gtf_sorted = sorted(gtf_standalone, key=lambda x: x["coord"])
        gtf_clusters = []
        current      = [gtf_sorted[0]]
        for tx in gtf_sorted[1:]:
            if abs(tx["coord"] - current[-1]["coord"]) <= tolerance:
                current.append(tx)
            else:
                gtf_clusters.append(current)
                current = [tx]
        gtf_clusters.append(current)

        for cluster in gtf_clusters:
            rep_coord = int(np.round(np.mean([t["coord"] for t in cluster])))
            groups.append({
                "type":          "GTF",
                "rep_coord":     rep_coord,
                "coord_min":     min(t["coord"] for t in cluster),
                "coord_max":     max(t["coord"] for t in cluster),
                "v2_sites":      [],
                "v2_fraction":   np.nan,
                "v2_usage":      np.nan,
                "v2_site_class": None,
                "v2_signal":     None,
                "num_polyAsites":0,
                "gtf_txs":       cluster,
            })

    # Sort proximal → distal (strand-aware)
    reverse = (strand == "-")
    groups.sort(key=lambda g: g["rep_coord"], reverse=reverse)
    return groups


# ---------------------------------------------------------------------------
# Row builder (Phase A — base metadata, no imputation yet)
# ---------------------------------------------------------------------------

def build_row(gene_name, group_idx, grp, strand, chrom,
              gene_meta, gene_merged_exons, gene_merged_term_exons,
              nmd_transcripts=None):
    panel_source    = grp["type"]
    coord           = grp["rep_coord"]
    gtf_txs         = grp["gtf_txs"]
    num_transcripts = len(gtf_txs)
    transcript_ids  = ";".join(t["transcript_id"]   for t in gtf_txs)
    transcript_names= ";".join(t["transcript_name"] for t in gtf_txs)

    gene_biotype   = gene_meta.get(gene_name, {}).get("biotype", "unknown")
    gene_is_coding = gene_biotype in CODING_BIOTYPES

    # pc_fraction and _coding_support (internal use only, not in CSV output)
    if gtf_txs:
        pc_count    = sum(1 for t in gtf_txs if t["is_coding"])
        pc_fraction = round(pc_count / num_transcripts, 3)
        has_coding    = any(t["is_coding"]     for t in gtf_txs)
        has_noncoding = any(not t["is_coding"] for t in gtf_txs)
        _coding_support = ("mixed"          if has_coding and has_noncoding else
                      "protein_coding" if has_coding else
                      "non_coding")
    else:
        pc_fraction = 1.0 if gene_is_coding else 0.0
        _coding_support  = "protein_coding" if gene_is_coding else "non_coding"

    # nmd / ig / tr / pseudogene fractions
    if gtf_txs:
        total_tx = len(gtf_txs)
        nmd_set  = nmd_transcripts or set()
        nmd_fraction       = sum(1 for t in gtf_txs if t["transcript_id"] in nmd_set)              / total_tx
        ig_fraction        = sum(1 for t in gtf_txs if t.get("biotype") in IG_BIOTYPES)            / total_tx
        tr_fraction        = sum(1 for t in gtf_txs if t.get("biotype") in TR_BIOTYPES)            / total_tx
        ig_pseudo_fraction = sum(1 for t in gtf_txs if t.get("biotype") in IG_PSEUDO_BIOTYPES)    / total_tx
        tr_pseudo_fraction = sum(1 for t in gtf_txs if t.get("biotype") in TR_PSEUDO_BIOTYPES)    / total_tx
    else:
        nmd_fraction = ig_fraction = tr_fraction = 0.0
        ig_pseudo_fraction = tr_pseudo_fraction = 0.0

    # cds_end_coord: stop codon position (mode across coding transcripts)
    # Sites sharing the same cds_end_coord truly compete in the same CPA pool
    coding_txs = [t for t in gtf_txs
                  if t.get("is_coding") and t.get("cds_end") is not None]
    if coding_txs:
        cds_end_coord = Counter(t["cds_end"] for t in coding_txs).most_common(1)[0][0]
        cds_group_id  = f"{gene_name}_CDS{cds_end_coord}"
    else:
        cds_end_coord = None
        cds_group_id  = f"{gene_name}_noCDS"

    # avg_utr_bp — Phase A: GTF/GTF+PA2 only
    avg_utr_bp = np.nan
    utr_confidence  = None

    if _coding_support in ("protein_coding", "mixed") and gtf_txs:
        coding_utrs = []
        tx_confidences = []
        for t in gtf_txs:
            if not t["is_coding"]:
                continue
            if t.get("spliced_utr") is not None:
                coding_utrs.append(t["spliced_utr"])
                tx_confidences.append(t.get("utr_confidence"))
        if coding_utrs:
            avg_utr_bp = round(float(np.mean(coding_utrs)), 2)
            # "spliced" only when every contributing transcript has a true
            # three_prime_utr annotation; any genomic fallback downgrades to "genomic".
            utr_confidence = ("spliced"
                              if all(c == "spliced" for c in tx_confidences)
                              else "genomic")

    # site_type: Auto-Calibration via General + Terminal Exon Trees
    # Applied to ALL groups regardless of GTF transcript support — a site with
    # GTF txs can still be ALE (non-terminal exon), IPA (intronic), or 3'Ext.
    tree      = gene_merged_exons.get(gene_name)
    term_tree = gene_merged_term_exons.get(gene_name)

    if tree:
        if not tree.at(coord):
            # Not in any exon: downstream of all exons → 3'Ext, otherwise → IPA
            if strand == "+":
                max_exon_end = max(iv.end for iv in tree)
                site_type = "3'Ext" if coord >= max_exon_end else "IPA"
            else:
                min_exon_start = min(iv.begin for iv in tree)
                site_type = "3'Ext" if coord <= min_exon_start else "IPA"
        else:
            # In an exon: terminal exon → Tandem_UTR, non-terminal → ALE
            if term_tree and term_tree.at(coord):
                site_type = "Tandem_UTR"
            else:
                site_type = "ALE"
    elif gtf_txs:
        # No exon tree for this gene but has transcript support — assume terminal
        site_type = "Tandem_UTR"
    else:
        site_type = "Unknown"

    # avg_utr_bp is only meaningful for Tandem_UTR and 3'Ext (UTR read-through).
    # IPA/ALE/Unknown inherit transcript UTR lengths from the gene but those
    # describe the canonical isoform, not this premature/alternative cleavage site.
    if site_type not in ("Tandem_UTR", "3'Ext"):
        avg_utr_bp     = np.nan
        utr_confidence = None

    # pa2_coords: PA2/GTF+PA2 → all PA2 cleavage site coords
    # GTF-only → rep_coord + most distal tx end (two anchors to avoid penalizing
    # reads far from the mean when group span is large)
    if grp["v2_sites"]:
        _pa2_coords = ";".join(str(s["coord"]) for s in grp["v2_sites"])
    elif grp["gtf_txs"]:
        rep    = int(coord)
        if strand == "-":
            distal = max(grp["coord_min"], rep - 420)
        else:
            distal = min(grp["coord_max"], rep + 420)
        _pa2_coords = str(rep) if distal == rep else f"{rep};{distal}"
    else:
        _pa2_coords = str(int(coord))

    return {
        "gene":             gene_name,
        "polyA_group":      group_idx,
        "rep_coord":        coord,
        "strand":           strand,
        "chrom":            chrom,
        "coord_min":        grp["coord_min"],
        "coord_max":        grp["coord_max"],
        "coord_spread":     grp["coord_max"] - grp["coord_min"],
        "avg_utr_bp":  avg_utr_bp,
        "utr_confidence":   utr_confidence,
        "_coding_support":  _coding_support,  # internal: coding/mixed/non_coding
        "pc_fraction":      pc_fraction,
        "panel_source":     panel_source,
        "num_transcripts":  num_transcripts,
        "num_polyAsites":   grp["num_polyAsites"],
        "pa2_coords":       _pa2_coords,
        "transcript_ids":   transcript_ids,
        "transcript_names": transcript_names,
        "v2_fraction":      grp["v2_fraction"],
        "v2_usage":         grp["v2_usage"],
        "v2_site_class":    grp["v2_site_class"],
        "PAS_motif":        grp.get("v2_signal"),
        "site_type":        site_type,
        # UTR intron flag: True if any coding transcript in this group
        # has a non-contiguous three_prime_utr (intron within 3'UTR)
        # Biological significance: EJC downstream of stop codon → NMD susceptibility
        "utr_has_intron":   any(t.get("utr_has_intron", False)
                                for t in gtf_txs if t["is_coding"]),
        "nmd_fraction":        round(nmd_fraction,       3),
        "ig_fraction":         round(ig_fraction,        3),
        "tr_fraction":         round(tr_fraction,        3),
        "ig_pseudo_fraction":  round(ig_pseudo_fraction, 3),
        "tr_pseudo_fraction":  round(tr_pseudo_fraction, 3),
        "cds_end_coord":       cds_end_coord,
        "cds_group_id":        cds_group_id,
        "is_singleton":     False,  # updated after all genes processed in main()
        # reserved for future integration
        "polyaDB_matched":  None,
        "isocape_validated":None,
        "cancer_relevant":  None,
        "functional_class": None,
    }


# ---------------------------------------------------------------------------
# Phase B: Cross-group UTR imputation with Intron Shield
# ---------------------------------------------------------------------------

def impute_pa2_utrs(gene_rows, gene_name, strand, gene_merged_exons):
    """
    For PA2-only protein_coding groups with NaN UTR:
    Find nearest GTF anchor (spliced UTR known), compute offset,
    validate via Exon Tree, write estimated UTR.

    Key design decisions:
    - Iterate ALL anchors (not just upstream), pick nearest VALID one
    - Skip anchors that yield negative utr_est (wrong direction)
      PA2-only sites are always >150bp from any GTF anchor after clustering
    - Exon Tree rules: no intron crossing, no upstream jump into intron
    """
    anchors = [r for r in gene_rows
               if r["utr_confidence"] == "spliced"
               and not np.isnan(r["avg_utr_bp"])]

    if not anchors:
        return

    tree = gene_merged_exons.get(gene_name)

    for r in gene_rows:
        if not (r["panel_source"] == "PA2" and
                r["_coding_support"] in ("protein_coding", "mixed") and
                r["site_type"] in ("Tandem_UTR", "3'Ext") and
                np.isnan(r["avg_utr_bp"])):
            continue

        pa2_coord      = r["rep_coord"]
        best_utr_est   = None
        best_dist      = float("inf")

        for a in anchors:
            anchor_coord = a["rep_coord"]
            dist         = abs(pa2_coord - anchor_coord)

            if dist > MAX_BORROW_DIST:
                continue

            # Compute offset (positive = PA2 is downstream of anchor)
            offset  = ((pa2_coord - anchor_coord) if strand == "+"
                       else (anchor_coord - pa2_coord))
            utr_est = a["avg_utr_bp"] + offset

            # Skip anchors that would yield deeply negative UTR
            # negative utr_est → anchor is on wrong side; skip
            # (PA2-only sites are always >150bp from GTF after clustering)
            if utr_est < 0:
                continue

            # Exon Tree validation
            is_valid = True
            if tree:
                start_c = min(anchor_coord, pa2_coord)
                end_c   = max(anchor_coord, pa2_coord)

                # Rule 1: no intron crossing (overlap > 1 exon block)
                if len(tree.overlap(start_c, end_c + 1)) > 1:
                    is_valid = False
                else:
                    # Rule 2: PA2 must not jump upstream into previous intron
                    anchor_blocks = tree.at(anchor_coord)
                    if anchor_blocks:
                        ablock = list(anchor_blocks)[0]
                        if strand == "+" and pa2_coord < ablock.begin:
                            is_valid = False
                        elif strand == "-" and pa2_coord >= ablock.end:
                            is_valid = False

            if is_valid and dist < best_dist:
                best_dist    = dist
                best_utr_est = utr_est

        # Write result
        if best_utr_est is not None:
            r["avg_utr_bp"] = round(float(best_utr_est), 2)
            r["utr_confidence"]  = "estimated"


# ---------------------------------------------------------------------------
# Phase C: CDS end coord imputation for PA2-only groups
# ---------------------------------------------------------------------------

def impute_cds_end_coords(gene_rows):
    """
    For PA2-only groups with no cds_end_coord (no GTF transcripts → no CDS info),
    inherit cds_end_coord from the nearest group within the same gene that has one.
    Preserves the Dual-Track framework: PA2-only sites join the correct CDS anchor.
    """
    have_cds    = [r for r in gene_rows if r.get('cds_end_coord') is not None]
    missing_cds = [r for r in gene_rows if r.get('cds_end_coord') is None]

    if not have_cds or not missing_cds:
        return

    for r in missing_cds:
        best = min(have_cds, key=lambda x: abs(x['rep_coord'] - r['rep_coord']))
        r['cds_end_coord'] = best['cds_end_coord']
        r['cds_group_id']  = best['cds_group_id']


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="IsoDecipher Dynamic Panel Builder (v2.0-first architecture)"
    )
    p.add_argument("--indices-pkl",  required=True,
                   help="Pickled reference indices from build_reference_indices.py")
    p.add_argument("--polyasite",    required=True,
                   help="PolyASite v2.0 BED file (polyasite2_hg38.bed)")
    p.add_argument("--out",          required=True,
                   help="Output panel CSV")
    p.add_argument("--tolerance",    type=int, default=150,
                   help="Clustering tolerance in bp (default: 150)")
    p.add_argument("--min-fraction", type=float, default=0.0,
                   help="Min v2.0 fraction to include a site (default: 0, keep all)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    # Load pre-built reference indices
    print(f"[panel] Loading reference indices: {args.indices_pkl}")
    with open(args.indices_pkl, "rb") as f:
        indices = pickle.load(f)
    gene_index        = indices["gene_index"]
    gene_meta         = indices["gene_meta"]
    transcript_3ends  = indices["transcript_3ends"]
    gene_merged_exons      = indices["gene_merged_exons"]
    gene_merged_term_exons = indices["gene_merged_term_exons"]
    nmd_transcripts        = indices.get("nmd_transcripts", set())
    print(f"[panel]   gene_index:        {sum(len(t) for t in gene_index.values()):,} genes")
    print(f"[panel]   nmd_transcripts:   {len(nmd_transcripts):,}")
    print(f"[panel]   transcript_3ends:  {len(transcript_3ends):,} genes")
    print(f"[panel]   gene_merged_exons: {len(gene_merged_exons):,} genes")

    # Load PolyASite v2.0
    print(f"[panel] Loading PolyASite v2.0: {args.polyasite}")
    v2 = pd.read_csv(
        args.polyasite, sep="\t", header=None, low_memory=False,
        names=["chrom","start","end","name","usage","strand",
               "fraction","count","score2","site_class","signal"],
    )
    v2["chrom"]     = v2["chrom"].astype(str)
    v2["rep_coord"] = v2["name"].str.split(":").str[1].astype(int)
    v2 = v2[v2["site_class"].isin(INCLUDE_SITE_CLASSES)].copy()
    print(f"[panel]   {len(v2):,} sites after class filter ({', '.join(sorted(INCLUDE_SITE_CLASSES))})")

    if args.min_fraction > 0:
        v2 = v2[v2["fraction"] >= args.min_fraction]
        print(f"[panel]   {len(v2):,} sites after fraction >= {args.min_fraction}")

    # Assign v2.0 sites to genes
    print(f"[panel] Assigning {len(v2):,} v2.0 sites to genes...")
    gene_v2_sites   = defaultdict(list)
    skipped_no_gene = 0

    for i, (_, row) in enumerate(v2.iterrows()):
        if i % 50000 == 0 and i > 0:
            print(f"[panel]   {i:,}/{len(v2):,} sites assigned...")
        chrom  = str(row["chrom"])
        coord  = int(row["rep_coord"])
        strand = row["strand"]
        gname, _, _ = assign_gene(chrom, coord, strand, gene_index)
        if gname is None:
            skipped_no_gene += 1
            continue
        gene_v2_sites[gname].append({
            "coord":      coord,
            "fraction":   float(row["fraction"]),
            "usage":      float(row["usage"]),
            "site_class": row["site_class"],
            "signal":     str(row["signal"]) if pd.notna(row["signal"]) else None,
            "chrom":      chrom,
            "strand":     strand,
        })

    print(f"[panel]   Assigned: {len(gene_v2_sites):,} genes "
          f"(skipped {skipped_no_gene:,} unassigned sites)")

    # Per-gene clustering and row building
    all_genes = set(gene_v2_sites.keys()) | set(transcript_3ends.keys())
    print(f"[panel] Building clusters for {len(all_genes):,} genes...")

    panel_rows = []
    total      = len(all_genes)

    for gi, gene_name in enumerate(sorted(all_genes)):
        if gi % 10000 == 0 and gi > 0:
            print(f"[panel]   Progress: {gi:,}/{total:,} genes...")

        meta   = gene_meta.get(gene_name, {})
        strand = meta.get("strand", "+")
        chrom  = meta.get("chrom",  "")
        if not chrom:
            continue

        v2_sites     = gene_v2_sites.get(gene_name, [])
        all_gene_txs = transcript_3ends.get(gene_name, [])

        groups = build_gene_clusters(
            gene_name, strand, v2_sites, all_gene_txs, args.tolerance
        )

        # Phase A: build base rows
        gene_rows = []
        for group_idx, grp in enumerate(groups):
            row = build_row(
                gene_name, group_idx, grp, strand, chrom,
                gene_meta, gene_merged_exons, gene_merged_term_exons,
                nmd_transcripts=nmd_transcripts,
            )
            gene_rows.append(row)

        # Phase B: cross-group UTR imputation with Intron Shield
        impute_pa2_utrs(gene_rows, gene_name, strand, gene_merged_exons)
        # Phase C: inherit cds_end_coord for PA2-only groups without GTF transcripts
        impute_cds_end_coords(gene_rows)

        panel_rows.extend(gene_rows)

    # Build DataFrame
    print(f"[panel] Building output DataFrame...")
    df = pd.DataFrame(panel_rows)

    # Singleton genes (only 1 polyA group) are KEPT in panel.
    # Used for insert size KDE calibration in validate_insert_size.py.
    # assign_reads.py filters them out before read assignment.
    groups_per_gene  = df.groupby("gene")["polyA_group"].max() + 1
    singleton_genes  = set(groups_per_gene[groups_per_gene == 1].index)
    n_singleton_feat = df["gene"].isin(singleton_genes).sum()
    df["is_singleton"] = df["gene"].isin(singleton_genes)

    df["pas_motif_strength"] = df.apply(
        lambda r: get_pas_motif_strength(r["PAS_motif"], r["panel_source"]), axis=1
    )

    # Summary
    print(f"\n{'='*60}")
    print(f"[SUMMARY] IsoDecipher Panel v1.1")
    print(f"{'='*60}")
    print(f"  Total features:  {len(df):,}")
    print(f"  Genes:           {df['gene'].nunique():,}")
    print(f"  Tolerance:       {args.tolerance} bp")

    print(f"\n  panel_source:")
    for src, cnt in df["panel_source"].value_counts().items():
        print(f"    {src:<12} {cnt:>8,}  ({cnt/len(df)*100:.1f}%)")

    print(f"\n  site_type (Auto-Calibrated via Exon Tree):")
    for st, cnt in df["site_type"].value_counts().items():
        print(f"    {st:<15} {cnt:>8,}")

    print(f"\n  utr_confidence:")
    for conf, cnt in df["utr_confidence"].value_counts(dropna=False).items():
        print(f"    {str(conf):<18} {cnt:>8,}")

    print(f"\n  Group count per gene:")
    print(f"    1 group  (singleton, KDE only): {len(singleton_genes):>8,} genes "
          f"({n_singleton_feat:,} features — kept for insert size calibration)")
    print(f"    2 groups:           {(groups_per_gene == 2).sum():>8,} genes")
    print(f"    3+ groups:          {(groups_per_gene >= 3).sum():>8,} genes")

    print(f"\n  pc_fraction distribution:")
    print(f"    == 1.0 (all coding):  {(df['pc_fraction'] == 1.0).sum():>8,}")
    print(f"    >= 0.5 (majority):    {(df['pc_fraction'] >= 0.5).sum():>8,}")
    print(f"    == 0.0 (no coding):   {(df['pc_fraction'] == 0.0).sum():>8,}")

    # Drop internal-use columns before saving
    df = df.drop(columns=["_coding_support"], errors="ignore")
    df.to_csv(args.out, index=False)
    print(f"\n[SUCCESS] Panel saved → {args.out}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
