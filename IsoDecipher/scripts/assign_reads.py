# ==============================================================================
# Copyright (c) 2026 IsoMatrix Suite. All rights reserved.
# 
# IsoDecipher is dual-licensed:
# 1. For academic and non-commercial use, it is licensed under the AGPLv3.
# 2. For commercial and enterprise use, a Commercial License is required.
# 
# See the LICENSE file in the project root for more details.
# ==============================================================================

import pysam
import pandas as pd
import gzip
from collections import defaultdict
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Assign reads to isoforms per sample")
    parser.add_argument("--bam",      required=True, help="Path to input BAM file")
    parser.add_argument("--panel",    default="results/panel_features.csv")
    parser.add_argument("--out",      required=True, help="Path to save output CSV")
    parser.add_argument("--barcodes", default=None,
                        help="Path to filtered barcodes (TSV, CSV, or .gz)")
    parser.add_argument("--barcode-col", type=int, default=None,
                        help="0-based column index for barcodes. Auto-detected if not set.")
    parser.add_argument("--barcode-sep", default=None,
                        help="Delimiter. Auto-detected from extension if not set.")
    parser.add_argument("--chroms", default=None,
                        help="Comma-separated chromosomes to include. "
                             "Default: standard chr1-22,chrX,chrY,chrM only. "
                             "Use 'all' to include all contigs.")
    parser.add_argument("--window", type=int, default=350,
                        help="Upstream read scatter window around PA site. "
                             "Default: 350 (based on 10x max insert 600bp - backbone 200bp - read2 90bp + buffer)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Flexible barcode loader (same logic as IsoCAPE)
# ---------------------------------------------------------------------------

def detect_sep(path):
    p = path.lower().replace('.gz', '')
    return '\t' if p.endswith('.tsv') else ','


def load_barcodes(barcodes_path, barcode_col=None, sep=None):
    """
    Load valid cell barcodes from Cell Ranger output.
    Handles:
    - Single-column TSV/CSV (standard barcodes.tsv.gz)
    - Two-column CSV with genome prefix (GRCh38,BARCODE-1)
    - Gzipped or plain text
    - Auto-detects barcode column
    """
    if sep is None:
        sep = detect_sep(barcodes_path)

    opener = gzip.open if barcodes_path.endswith('.gz') else open
    with opener(barcodes_path, 'rt') as fh:
        first_line = fh.readline().strip()

    fields = first_line.split(sep)
    n_cols = len(fields)

    if barcode_col is None:
        barcode_col = 0
        for i, f in enumerate(fields):
            cleaned = f.replace('-1', '').split('-')[0]
            if all(c in 'ACGTacgt' for c in cleaned) and len(cleaned) >= 12:
                barcode_col = i
                break

    print(f"[filter] Barcode file: {n_cols} column(s), sep='{sep}', using column {barcode_col}")

    df = pd.read_csv(barcodes_path, header=None, sep=sep)
    barcodes = set(df[barcode_col].astype(str))

    print(f"[filter] Loaded {len(barcodes):,} filtered barcodes from {barcodes_path}")
    return barcodes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    panel = pd.read_csv(args.panel)

    # Standard chromosomes filter
    STANDARD_CHROMS = set(
        [f"chr{i}" for i in range(1, 23)] +
        ["chrX", "chrY", "chrM"] +
        [str(i) for i in range(1, 23)] +
        ["X", "Y", "MT"]
    )

    if args.chroms == 'all':
        allowed_chroms = None  # no filter
        print("[filter] Chromosome filter: OFF (all contigs included)")
    elif args.chroms:
        allowed_chroms = set(c.strip() for c in args.chroms.split(','))
        print(f"[filter] Chromosome filter: {len(allowed_chroms)} user-specified chroms")
    else:
        allowed_chroms = STANDARD_CHROMS
        print("[filter] Chromosome filter: standard chr1-22,chrX,chrY,chrM "
              "(use --chroms all to include contigs)")

    targets = defaultdict(list)
    skipped_chroms = set()

    for row in panel.itertuples(index=False):
        chrom = row.chrom
        if allowed_chroms is not None and chrom not in allowed_chroms:
            # try with/without chr prefix
            alt = chrom.lstrip('chr') if chrom.startswith('chr') else f"chr{chrom}"
            if alt in allowed_chroms:
                chrom = alt
            else:
                skipped_chroms.add(row.chrom)
                continue
        targets[chrom].append({
            'gene':        row.gene,
            'pos':         row.rep_coord,
            'strand':      row.strand,
            'group':       row.polyA_group,
            'label':       row.user_label,
            'spliced_utr': row.avg_spliced_utr,
            'genomic_utr': row.avg_genomic_utr,
        })

    if skipped_chroms:
        print(f"[filter] Skipped {len(skipped_chroms)} non-standard contigs: "
              f"{', '.join(sorted(skipped_chroms)[:5])}{'...' if len(skipped_chroms) > 5 else ''}")

    # Load barcodes with flexible loader
    valid_barcodes = None
    if args.barcodes:
        valid_barcodes = load_barcodes(
            args.barcodes,
            barcode_col=args.barcode_col,
            sep=args.barcode_sep,
        )

    bam    = pysam.AlignmentFile(args.bam, "rb")
    
    # Use dynamic window from arguments
    window = args.window

    counts        = defaultdict(set)
    umi_best_match = {}

    for chrom, sites in targets.items():
        current_chrom = chrom
        if current_chrom not in bam.references:
            if f"chr{current_chrom}" in bam.references:
                current_chrom = f"chr{current_chrom}"
            else:
                print(f"⚠️  Contig {chrom} not found in BAM. Skipping...")
                continue

        print(f"Current Chromosome: {current_chrom} | "
              f"Total unique UMIs captured so far: {len(umi_best_match)}")

        for site in sites:
            # Fetch upstream only — downstream reads belong to neighbor genes
            if site['strand'] == '+':
                fetch_start = max(0, site['pos'] - window)
                fetch_end   = site['pos']
            else:
                fetch_start = site['pos']
                fetch_end   = site['pos'] + window
            for read in bam.fetch(current_chrom, fetch_start, fetch_end):
                if not (read.has_tag("CB") and read.has_tag("UB")):
                    continue
                if (site["strand"] == "+" and read.is_reverse) or \
                   (site["strand"] == "-" and not read.is_reverse):
                    continue

                cb = read.get_tag("CB")

                if valid_barcodes is not None and cb not in valid_barcodes:
                    continue

                read_3_prime = (read.reference_end
                                if site['strand'] == '+'
                                else read.reference_start)
                dist = abs(read_3_prime - site['pos'])

                if dist <= window:
                    ub      = read.get_tag("UB")
                    umi_key = (cb, ub)
                    label   = site['label']
                    if not label or label == 'N/A' or str(label) == 'nan':
                        feature = f"{site['gene']}_G{site['group']}"
                    else:
                        feature = f"{site['gene']}_G{site['group']}_{label}"

                    if umi_key not in umi_best_match or dist < umi_best_match[umi_key][1]:
                        umi_best_match[umi_key] = (feature, dist)

    for (cb, ub), (feature, dist) in umi_best_match.items():
        counts[(cb, feature)].add(ub)

    results = []
    for (cb, feature), umis in counts.items():
        results.append({
            'cell_barcode': cb,
            'feature':      feature,
            'count':        len(umis),
        })

    if results:
        df     = pd.DataFrame(results)
        matrix = df.pivot(index='cell_barcode', columns='feature',
                          values='count').fillna(0)
        matrix.to_csv(args.out)
        print(f"✅ Saved: {args.out}")
        print(f"   Cells:    {matrix.shape[0]:,}")
        print(f"   Features: {matrix.shape[1]:,}")
        print(f"   UMIs:     {int(matrix.values.sum()):,}")
    else:
        print(f"⚠️  No reads assigned for {args.bam}")


if __name__ == "__main__":
    main()
