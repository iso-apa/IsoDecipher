# ==============================================================================
# Copyright (c) 2026 IsoMatrix Suite. All rights reserved.
# 
# IsoDecipher is dual-licensed:
# 1. For academic and non-commercial use, it is licensed under the AGPLv3.
# 2. For commercial and enterprise use, a Commercial License is required.
# 
# See the LICENSE file in the project root for more details.
# ==============================================================================

import os

# Default lookup bundled with IsoDecipher
_SCRIPT_DIR    = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_LOOKUP = os.path.join(
    _SCRIPT_DIR, '..', 'reference',
    'insert_size_lookup_10x_3p_v3.json'
)

import pysam
import pandas as pd
import gzip
from collections import defaultdict
import argparse

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
    parser.add_argument("--window", type=int, default=420,
                        help="Upstream read scatter window around PA site. "
                             "Default: 420bp (P95 of empirical insert size distribution + buffer). "
                             "Previously 350bp.")
    parser.add_argument("--insert-size-lookup", default=_DEFAULT_LOOKUP,
                        help="KDE lookup table for probabilistic read assignment. "
                             "Default: bundled 10x 3' v3 reference "
                             "(IsoDecipher/reference/insert_size_lookup_10x_3p_v3.json). "
                             "Generate custom lookup with validate_insert_size.py "
                             "if using a different library prep. "
                             "Pass 'none' to use shortest-distance assignment instead.")
    parser.add_argument("--debug-gene", default=None,
                        help="Gene name to output per-read assignment details for IGV validation. "
                             "Outputs a TSV with read coordinates, site assigned, offset, and score.")
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

    # Load KDE insert size lookup for probabilistic assignment
    prob_array = None
    lookup_path = args.insert_size_lookup
    if lookup_path and lookup_path.lower() != 'none':
        if os.path.exists(lookup_path):
            import json
            import numpy as np
            lookup    = json.load(open(lookup_path))
            prob_array = np.array(lookup['probs'])
            print(f"[assign] Probabilistic assignment: KDE lookup loaded")
            print(f"         {lookup_path}")
            print(f"         median_offset={lookup.get('offset_mu', 'N/A')}bp | "
                  f"window={args.window}bp")
        else:
            print(f"[assign] WARNING: lookup file not found: {lookup_path}")
            print(f"         Falling back to shortest-distance assignment.")
            print(f"         Run validate_insert_size.py to generate a lookup table.")
    else:
        print(f"[assign] Shortest-distance assignment (--insert-size-lookup none).")

    def score_read(dist):
        """
        Return assignment score for a read at given offset distance.
        Higher score = better assignment.
        KDE mode: probability from empirical distribution
        Fallback:  negative distance (closer = higher score)
        """
        if prob_array is not None:
            if dist < 0 or dist > len(prob_array) - 1:
                return 0.0
            return float(prob_array[int(dist)])
        else:
            return -dist  # negative distance: closer = less negative = higher

    # Use dynamic window from arguments
    window = args.window

    # Debug mode: if debug-gene specified, only process that gene's chromosome
    debug_gene = args.debug_gene
    if debug_gene:
        gene_rows = panel[panel['gene'] == debug_gene]
        if gene_rows.empty:
            print(f"[debug] ERROR: gene '{debug_gene}' not found in panel")
            return
        debug_chrom = gene_rows['chrom'].iloc[0]
        # normalize chrom prefix
        if not str(debug_chrom).startswith('chr'):
            debug_chrom = f'chr{debug_chrom}'
        print(f"[debug] Gene {debug_gene} → only processing {debug_chrom}")
        # filter targets to only this chrom
        targets = {k: v for k, v in targets.items()
                   if (k.startswith('chr') and k == debug_chrom) or
                      (not k.startswith('chr') and f'chr{k}' == debug_chrom)}

    debug_rows   = []
    umi_debug_info = {}

    # Global tracker for the final output matrix
    global_soft_counts = defaultdict(float)  # (cb, feature) -> fractional count
    global_umi_count = 0                     # Total UMIs processed

    for chrom, sites in targets.items():
        current_chrom = chrom
        if current_chrom not in bam.references:
            if f"chr{current_chrom}" in bam.references:
                current_chrom = f"chr{current_chrom}"
            else:
                print(f"⚠️  Contig {chrom} not found in BAM. Skipping...")
                continue

        print(f"Current Chromosome: {current_chrom} | "
              f"Features processed so far: {len(global_soft_counts)}")

        # Initialize tracking dict locally per chromosome to prevent memory bloat
        umi_best_match = {}  # (cb, ub) -> {feature: score, ...}

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
                    # user_label (e.g. Secreted/Membrane for IGH genes) is stored
                    # in panel CSV as metadata only — NOT embedded in feature name.
                    # Keeping feature names as Gene_G0 / Gene_G1 ensures that
                    # integrate_samples.py's panel join (keyed on gene + '_G' + group)
                    # correctly matches IGH features and populates utr_source /
                    # avg_spliced_utr from the panel CSV.
                    feature = f"{site['gene']}_G{site['group']}"

                    s = score_read(dist)

                    # Soft assignment: accumulate scores for ALL sites
                    # each UMI -> dict of {feature: score}
                    if umi_key not in umi_best_match:
                        umi_best_match[umi_key] = {}
                    umi_best_match[umi_key][feature] = max(
                        umi_best_match[umi_key].get(feature, 0), s
                    )

                    # Debug tracking (keep highest-score site for display)
                    if debug_gene and site['gene'] == debug_gene:
                        prev = umi_debug_info.get(umi_key, {})
                        if s > prev.get('score', -1):
                            umi_debug_info[umi_key] = {
                                'read_id':    read.query_name,
                                'cb':         cb,
                                'ub':         ub,
                                'chrom':      current_chrom,
                                'read_start': read.reference_start,
                                'read_end':   read.reference_end,
                                'read_3p':    read_3_prime,
                                'site':       feature,
                                'site_pos':   site['pos'],
                                'offset_bp':  dist,
                                'score':      round(s, 6),
                                'method':     'KDE' if prob_array is not None else 'shortest_dist',
                                'strand':     site['strand'],
                            }

        # Soft assignment for the current chromosome: normalize scores per UMI
        for umi_key, site_scores in umi_best_match.items():
            cb = umi_key[0]
            total_score = sum(site_scores.values())
            if total_score <= 0:
                continue
            for feature, score in site_scores.items():
                global_soft_counts[(cb, feature)] += score / total_score

        # Track total UMIs and explicitly flush memory for this chromosome
        global_umi_count += len(umi_best_match)
        umi_best_match.clear()
        import gc
        gc.collect()

    # Debug output
    if debug_gene and umi_debug_info:
        debug_path = args.out.replace('.csv', f'_debug_{debug_gene}.tsv')
        debug_df   = pd.DataFrame(list(umi_debug_info.values()))
        debug_df   = debug_df.sort_values('offset_bp')
        debug_df.to_csv(debug_path, sep='\t', index=False)
        print(f"\n[debug] {debug_gene} assignment details → {debug_path}")
        print(f"[debug] {len(debug_df)} UMIs assigned to {debug_gene}")
        print(f"\n[debug] Site distribution:")
        print(debug_df['site'].value_counts().to_string())
        print(f"\n[debug] Offset distribution per site:")
        print(debug_df.groupby('site')['offset_bp'].describe().round(1).to_string())
        print(f"\n[debug] Sample reads (sorted by offset):")
        cols = ['read_id','chrom','read_start','read_end','read_3p',
                'site','site_pos','offset_bp','score','method']
        print(debug_df[cols].head(20).to_string(index=False))
        print(f"\n[debug] IGV coordinates for top reads:")
        for _, row in debug_df.head(10).iterrows():
            igv = f"{row['chrom']}:{row['read_start']}-{row['read_end']}"
            print(f"  {igv:<35} offset={row['offset_bp']:>4}bp  → {row['site']}")

    results = []
    for (cb, feature), frac_count in global_soft_counts.items():
        results.append({
            'cell_barcode': cb,
            'feature':      feature,
            'count':        round(frac_count, 4),  # 4 decimal places
        })

    if results:
        df     = pd.DataFrame(results)
        matrix = df.pivot(index='cell_barcode', columns='feature',
                          values='count').fillna(0)
        matrix.to_csv(args.out)
        print(f"✅ Saved: {args.out}")
        print(f"   Cells:    {matrix.shape[0]:,}")
        print(f"   Features: {matrix.shape[1]:,}")
        print(f"   UMIs:     {global_umi_count:,} (soft-assigned → fractional counts)")
    else:
        print(f"⚠️  No reads assigned for {args.bam}")


if __name__ == "__main__":
    main()