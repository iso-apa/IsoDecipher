#!/usr/bin/env python3
# ==============================================================================
# Copyright (c) 2026 IsoMatrix Suite. All rights reserved.
# ==============================================================================
"""
IsoDecipher: build_reference_indices.py
---------------------------------------
Phase 1/2 of Panel Building: Reference Indexer
Parses the GTF database to build static IntervalTrees and exact Spliced UTR records.
Outputs a pickled dictionary containing:
  - gene_index
  - gene_meta
  - transcript_3ends
  - gene_merged_exons (Intron Shield)

Usage:
  python build_reference_indices.py \
      --gtf-db  /path/to/Homo_sapiens.GRCh38.115.gtf.db \
      --out-pkl reference/hg38_iso_indices.pkl
"""

import argparse
import os
import pickle
from collections import defaultdict
import gffutils
from intervaltree import IntervalTree
from _constants import CODING_BIOTYPES


def parse_args():
    p = argparse.ArgumentParser(description="IsoDecipher Reference Indexer")
    p.add_argument("--gtf-db",  required=True, help="gffutils SQLite DB (.db)")
    p.add_argument("--out-pkl", required=True, help="Output pickle file path (.pkl)")
    return p.parse_args()


def build_gene_index(db):
    print("[ref_builder] Building gene index (IntervalTree)...")
    gene_index = defaultdict(IntervalTree)
    gene_meta  = {}
    n = 0
    for g in db.features_of_type("gene"):
        chrom   = str(g.seqid)
        gname   = g.attributes.get("gene_name", [g.id])[0]
        biotype = g.attributes.get("gene_biotype", ["unknown"])[0]
        gene_index[(chrom, g.strand)][g.start:g.end + 1] = {
            "start":   g.start,
            "end":     g.end,
            "name":    gname,
            "biotype": biotype,
        }
        gene_meta[gname] = {
            "biotype": biotype,
            "chrom":   chrom,
            "strand":  g.strand,
            "start":   g.start,
            "end":     g.end,
        }
        n += 1
    print(f"[ref_builder]   {n:,} genes loaded into IntervalTree")
    return gene_index, gene_meta


def build_transcript_index(db, gene_meta):
    print("[ref_builder] Building transcript + UTR index...")

    # Pass 1: tx_id → gene_name mapping + retained_intron blacklist + canonical transcript selection
    print("[ref_builder]   Pass 1/5: mapping transcripts to genes...")
    tx_to_gene          = {}
    # Biotypes excluded from Exon Tree:
    # - retained_intron: their "exons" span true introns, destroying gap structure
    # - protein_coding_CDS_not_defined: large unverified exons that swallow intron regions
    EXON_TREE_SKIP_BIOTYPES = {"retained_intron", "protein_coding_CDS_not_defined"}

    # Canonical transcript priority (lower = better).
    # MANE_Select: single NCBI+Ensembl-agreed canonical per gene — matches IGV default.
    # Ensembl_canonical: Ensembl's own pick when no MANE (e.g. non-human genes).
    # appris_principal_*: APPRIS structural scores, used as tertiary fallback.
    # Anything without a priority tag gets 99 (last resort).
    _CANONICAL_PRIORITY = {
        "MANE_Select":        0,
        "Ensembl_canonical":  1,
        "appris_principal_1": 2,
        "appris_principal_2": 3,
        "appris_principal_3": 4,
        "appris_principal_4": 5,
        "appris_principal_5": 6,
    }
    canonical_tx_per_gene = {}  # gene_name → (priority, tx_id)

    retained_intron_txs = set()
    nmd_transcripts     = set()
    for tx in db.features_of_type("transcript"):
        tx_id    = tx.attributes.get("transcript_id", [tx.id])[0]
        gname    = tx.attributes.get("gene_name", [None])[0]
        tx_bio   = tx.attributes.get("transcript_biotype", [""])[0]
        tags     = tx.attributes.get("tag", [])
        if gname:
            tx_to_gene[tx_id] = gname
        if tx_bio in EXON_TREE_SKIP_BIOTYPES:
            retained_intron_txs.add(tx_id)
        if tx_bio == "nonsense_mediated_decay":
            nmd_transcripts.add(tx_id)

        # Canonical selection: skip excluded biotypes; pick highest-priority tag.
        if gname and tx_bio not in EXON_TREE_SKIP_BIOTYPES:
            best = 99
            for tag in tags:
                p = _CANONICAL_PRIORITY.get(tag)
                if p is not None and p < best:
                    best = p
            existing = canonical_tx_per_gene.get(gname)
            if existing is None or best < existing[0]:
                canonical_tx_per_gene[gname] = (best, tx_id)

    print(f"[ref_builder]   Transcripts excluded from Exon Tree: {len(retained_intron_txs):,} "
          f"(retained_intron + protein_coding_CDS_not_defined)")
    print(f"[ref_builder]   NMD transcripts: {len(nmd_transcripts):,}")
    print(f"[ref_builder]   Canonical transcripts selected: {len(canonical_tx_per_gene):,} genes")

    # Pass 2: CDS per transcript
    print("[ref_builder]   Pass 2/5: loading CDS features...")
    cds_by_tx = defaultdict(list)
    for c in db.features_of_type("CDS"):
        tx_id = c.attributes.get("transcript_id", [None])[0]
        if tx_id:
            cds_by_tx[tx_id].append(c)

    # Pass 3: three_prime_utr lengths + intron-in-UTR detection
    print("[ref_builder]   Pass 3/5: loading three_prime_utr features...")
    utr_lens        = defaultdict(int)
    utr_feats_by_tx = defaultdict(list)

    for u in db.features_of_type("three_prime_utr"):  # Ensembl GTF uses lowercase
        tx_id = u.attributes.get("transcript_id", [None])[0]
        if not tx_id:
            continue
        utr_lens[tx_id] += (u.end - u.start + 1)
        utr_feats_by_tx[tx_id].append((u.start, u.end))

    # detect transcripts whose 3'UTR contains an intron
    # (multiple non-contiguous three_prime_utr features)
    utr_has_intron = set()
    for tx_id, feats in utr_feats_by_tx.items():
        if len(feats) > 1:
            sorted_feats = sorted(feats)
            for i in range(1, len(sorted_feats)):
                if sorted_feats[i][0] > sorted_feats[i-1][1] + 1:
                    utr_has_intron.add(tx_id)
                    break

    print(f"[ref_builder]   Transcripts with intron in 3\'UTR: {len(utr_has_intron):,}")

    # Pass 4: build transcript_3ends
    print("[ref_builder]   Pass 4/5: building transcript 3' end index...")
    transcript_3ends = defaultdict(list)
    n_tx = 0
    n_spliced   = 0
    n_genomic   = 0
    n_no_utr    = 0

    for tx in db.features_of_type("transcript"):
        n_tx += 1
        gname = tx.attributes.get("gene_name", [None])[0]
        if not gname:
            continue
        meta = gene_meta.get(gname)
        if not meta:
            continue

        gstrand      = meta["strand"]
        gene_biotype = meta["biotype"]
        tx_biotype   = tx.attributes.get("transcript_biotype", [gene_biotype])[0]
        tx_id        = tx.attributes.get("transcript_id",   [""])[0]
        tx_name      = tx.attributes.get("transcript_name", [tx_id])[0]

        tx_3end = tx.end if gstrand == "+" else tx.start

        # CDS end (stop codon position) — computed for ALL transcripts upfront
        # Used for CDS-anchored APA grouping regardless of UTR confidence path
        cds_feats = cds_by_tx.get(tx_id, [])
        if cds_feats:
            cds_end = (max(c.end   for c in cds_feats) if gstrand == "+"
                       else min(c.start for c in cds_feats))
        else:
            cds_end = None

        # Spliced UTR: three_prime_UTR features (most accurate)
        if tx_id in utr_lens:
            spliced_utr     = utr_lens[tx_id]
            utr_confidence  = "spliced"
            n_spliced += 1
        else:
            # Fallback: genomic distance from tx 3' end to CDS end
            # Valid only when no intron lies between CDS stop and transcript end;
            # labeled "genomic" to distinguish from truly spliced three_prime_utr features.
            if cds_end is not None:
                est = (tx_3end - cds_end) if gstrand == "+" else (cds_end - tx_3end)
                if est >= 0:
                    spliced_utr    = est
                    utr_confidence = "genomic"
                    n_genomic += 1
                else:
                    spliced_utr    = None
                    utr_confidence = None
                    n_no_utr += 1
            else:
                spliced_utr    = None
                utr_confidence = None
                n_no_utr += 1

        is_coding = (
            tx_biotype in CODING_BIOTYPES or
            (tx_biotype == "retained_intron" and gene_biotype in CODING_BIOTYPES)
        )

        transcript_3ends[gname].append({
            "coord":           tx_3end,
            "transcript_id":   tx_id,
            "transcript_name": tx_name,
            "biotype":         tx_biotype,
            "spliced_utr":     spliced_utr,
            "utr_confidence":  utr_confidence,
            "is_coding":       is_coding,
            "utr_has_intron":  tx_id in utr_has_intron,
            "cds_end":         cds_end,
        })

    print(f"[ref_builder]   {n_tx:,} transcripts processed")
    print(f"[ref_builder]   spliced UTR from three_prime_UTR: {n_spliced:,}")
    print(f"[ref_builder]   spliced UTR from genomic fallback: {n_genomic:,}")
    print(f"[ref_builder]   no UTR (non-coding or no CDS):    {n_no_utr:,}")

    # Pass 5: General Exon Tree (Intron Shield) + Terminal Exon Tree (ALE Calibration)
    print("[ref_builder]   Pass 5/5: building merged exon trees (Intron Shield + ALE)...")
    gene_exons      = defaultdict(IntervalTree)
    gene_term_exons = defaultdict(IntervalTree)  # terminal exon only

    for tx in db.features_of_type("transcript"):
        tx_id = tx.attributes.get("transcript_id", [""])[0]
        if tx_id not in tx_to_gene:
            continue
        if tx_id in retained_intron_txs:
            continue

        gname  = tx_to_gene[tx_id]
        meta   = gene_meta.get(gname)
        if not meta:
            continue
        gstrand = meta["strand"]

        exons = sorted(db.children(tx, featuretype="exon"), key=lambda e: e.start)
        if not exons:
            continue

        # General Exon Tree
        for ex in exons:
            gene_exons[gname].addi(ex.start, ex.end + 1)

        # Terminal Exon Tree
        # - strand: terminal exon = smallest coord (first in ascending sort)
        # + strand: terminal exon = largest coord (last in ascending sort)
        terminal = exons[0] if gstrand == "-" else exons[-1]
        gene_term_exons[gname].addi(terminal.start, terminal.end + 1)

    gene_merged_exons      = {}
    gene_merged_term_exons = {}

    for gname in gene_exons:
        tree_copy = gene_exons[gname].copy()
        tree_copy.merge_overlaps()
        gene_merged_exons[gname] = tree_copy

        if gname in gene_term_exons:
            term_copy = gene_term_exons[gname].copy()
            term_copy.merge_overlaps()
            gene_merged_term_exons[gname] = term_copy

    print(f"[ref_builder]   General exon trees: {len(gene_merged_exons):,} genes")
    print(f"[ref_builder]   Terminal exon trees: {len(gene_merged_term_exons):,} genes")

    # Canonical exon structure: single isoform per gene (MANE_Select > Ensembl_canonical > APPRIS P1).
    # Used by IsoDecipher Gaussian fit to walk upstream exons without phantom exons from union merging.
    # Stored as a sorted list of (start, end) tuples in ascending genomic order.
    # Caller uses gene_meta[gene]['strand'] to determine which end is terminal vs upstream.
    print("[ref_builder]   Building canonical exon structure (MANE_Select > Ensembl_canonical > APPRIS)...")
    gene_canonical_exons = {}
    _prio_counts = [0, 0, 0, 0]  # MANE, Ensembl_canonical, APPRIS, fallback

    for gname, (prio, tx_id) in canonical_tx_per_gene.items():
        try:
            tx_feat = db[tx_id]
        except Exception:
            continue
        exons = sorted(db.children(tx_feat, featuretype="exon"), key=lambda e: e.start)
        if not exons:
            continue
        gene_canonical_exons[gname] = [(e.start, e.end) for e in exons]
        if prio == 0:   _prio_counts[0] += 1
        elif prio == 1: _prio_counts[1] += 1
        elif prio <= 6: _prio_counts[2] += 1
        else:           _prio_counts[3] += 1

    print(f"[ref_builder]   Canonical exon genes: {len(gene_canonical_exons):,} "
          f"(MANE={_prio_counts[0]:,} | Ensembl_can={_prio_counts[1]:,} | "
          f"APPRIS={_prio_counts[2]:,} | fallback={_prio_counts[3]:,})")
    print(f"[ref_builder]   Transcript 3' end index: {len(transcript_3ends):,} genes")
    return transcript_3ends, gene_merged_exons, gene_merged_term_exons, nmd_transcripts, gene_canonical_exons


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out_pkl)), exist_ok=True)

    print(f"[ref_builder] Loading GTF DB: {args.gtf_db}")
    db = gffutils.FeatureDB(args.gtf_db, keep_order=True)

    gene_index, gene_meta             = build_gene_index(db)
    transcript_3ends, gene_merged_exons, gene_merged_term_exons, nmd_transcripts, gene_canonical_exons = build_transcript_index(db, gene_meta)

    indices = {
        "gene_index":             gene_index,
        "gene_meta":              gene_meta,
        "transcript_3ends":       transcript_3ends,
        "gene_merged_exons":      gene_merged_exons,
        "gene_merged_term_exons": gene_merged_term_exons,
        "nmd_transcripts":        nmd_transcripts,
        "gene_canonical_exons":   gene_canonical_exons,
    }

    print(f"[ref_builder] Serializing indices to {args.out_pkl}...")
    with open(args.out_pkl, "wb") as f:
        pickle.dump(indices, f)
    print(f"[SUCCESS] Reference indices saved → {args.out_pkl}")


if __name__ == "__main__":
    main()
