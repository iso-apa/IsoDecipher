# ==============================================================================
# Copyright (c) 2026 IsoMatrix Suite. All rights reserved.
# 
# IsoDecipher is dual-licensed:
# 1. For academic and non-commercial use, it is licensed under the AGPLv3.
# 2. For commercial and enterprise use, a Commercial License is required.
# 
# See the LICENSE file in the project root for more details.
# ==============================================================================

import argparse
import gc
import gzip
import json
import os
import pickle
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

# Default lookup bundled with IsoDecipher
_SCRIPT_DIR    = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_LOOKUP = os.path.join(
    _SCRIPT_DIR, '..', 'reference',
    'insert_size_lookup_10x_3p_v3.json'
)

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
    parser.add_argument("--window", type=int, default=None,
                        help="Fetch window upstream of PA site (bp). Default: 420bp.")
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
    parser.add_argument("--indices-pkl", default=None,
                        help="Path to hg38_iso_indices.pkl from build_reference_indices.py. "
                             "Enables exon-space read fetching and scoring for Tandem_UTR and ALE sites, "
                             "correcting for short terminal exons where reads pile up in upstream exons.")
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
# Exon-space helpers (Tandem_UTR / ALE short-exon correction)
# ---------------------------------------------------------------------------

def _get_exon_list_for_site(gene, pa_coord, site_type, gene_canonical_exons, gene_merged_exons_sorted):
    """
    Returns exon list [(start, end), ...] with 0-based coords sorted ascending.
    Tandem_UTR → MANE_Select canonical exons (matches IGV, no phantom exons).
    ALE → canonical if pa_coord falls inside a canonical exon; otherwise merged union
          (ALE exon M may not be present in the canonical isoform).
    Returns None if no exon data available for this gene.
    """
    canonical = gene_canonical_exons.get(gene)

    if site_type == 'Tandem_UTR':
        if canonical:
            return canonical
        # Fallback: merged union exons when canonical not yet in pkl (pre-rebuild).
        # Phantom exons rare for Tandem_UTR terminal regions; exon walk still correct.
        return gene_merged_exons_sorted.get(gene)

    if site_type == 'ALE':
        if canonical and any(s <= pa_coord <= e for s, e in canonical):
            return canonical
        # ALE exon not in canonical transcript → use union merged exons
        return gene_merged_exons_sorted.get(gene)

    return None


def _pa_exon_upstream_bp(pa_coord, strand, exon_list):
    """
    Return how many bp are available upstream of pa_coord within its own exon.
    If pa_coord falls in a long exon, this will be >= window → no exon walk needed.
    """
    for ex_s, ex_e in exon_list:
        if ex_s <= pa_coord <= ex_e:
            return (pa_coord - ex_s) if strand == '+' else (ex_e - pa_coord)
    return 0  # PA not in any annotated exon


def _get_exon_space_intervals(pa_coord, strand, exon_list, window=420):
    """
    Walk upstream from pa_coord through the exon chain, collecting up to `window` bp of
    spliced exon space.  Continues past exon N-1, N-2, ... until the budget is spent —
    handles genes (e.g. RPL/RPS) where multiple consecutive exons are shorter than `window`.

    exon_list : list of (start, end) 0-based inclusive, sorted ascending.
    Returns   : list of (start, end) 0-based inclusive genomic intervals, sorted ascending.
                These are the pysam fetch targets (caller adds +1 to end for half-open fetch).
    """
    if not exon_list:
        return []

    intervals = []
    remaining = window

    if strand == '+':
        # Upstream = lower genomic coords.  Walk the exon list in reverse (high→low).
        # Clip each exon's right boundary to pa_coord so we stay upstream.
        valid = [(s, min(e, pa_coord)) for s, e in exon_list if s <= pa_coord]
        for ex_s, ex_e in reversed(valid):
            if remaining <= 0:
                break
            span = ex_e - ex_s + 1
            if span <= 0:
                continue
            if span >= remaining:
                intervals.append((ex_e - remaining + 1, ex_e))
                remaining = 0
            else:
                intervals.append((ex_s, ex_e))
                remaining -= span

    else:  # '-'
        # Upstream on - strand = higher genomic coords.  Walk exon list forward (low→high).
        # Clip each exon's left boundary to pa_coord so we stay upstream.
        valid = [(max(s, pa_coord), e) for s, e in exon_list if e >= pa_coord]
        for ex_s, ex_e in valid:
            if remaining <= 0:
                break
            span = ex_e - ex_s + 1
            if span <= 0:
                continue
            if span >= remaining:
                intervals.append((ex_s, ex_s + remaining - 1))
                remaining = 0
            else:
                intervals.append((ex_s, ex_e))
                remaining -= span

    return sorted(intervals)


def _exon_space_dist(read_3p, pa_coord, exon_list, strand):
    """
    Spliced (exon-space) distance from read_3p to pa_coord.
    Counts only exonic bases between the two positions — intronic gaps are not counted.
    Returns the same unit as the current KDE lookup (insert size in spliced bp).

    exon_list : list of (start, end) 0-based inclusive, sorted ascending.
    Returns   : distance in bp (non-negative).  Returns genomic |dist| if exon_list empty.
    """
    if not exon_list:
        return abs(pa_coord - read_3p)

    lo = read_3p if strand == '+' else pa_coord
    hi = pa_coord if strand == '+' else read_3p   # half-open upper bound

    if lo >= hi:
        return 0

    dist = 0
    for ex_s, ex_e in exon_list:
        # Overlap of [lo, hi) with [ex_s, ex_e] (inclusive end → ex_e+1 exclusive)
        ov_s = max(ex_s, lo)
        ov_e = min(ex_e + 1, hi)
        if ov_e > ov_s:
            dist += ov_e - ov_s
    return dist


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    panel = pd.read_csv(args.panel)

    # Panel is expected to be pre-filtered by prep_sample_panel.py
    # (singleton removal + expressed gene filter + chrom filter)
    print(f"[assign] Panel loaded: {len(panel):,} features, "
          f"{panel['gene'].nunique():,} genes")

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
            'cmin':        getattr(row, 'coord_min', row.rep_coord),
            'cmax':        getattr(row, 'coord_max', row.rep_coord),
            'pa2_coords':  [int(float(c)) for c in str(getattr(row, 'pa2_coords', '') or '').split(';')
                            if c.strip() and c.strip().lower() != 'nan'],
            'strand':      row.strand,
            'group':       row.polyA_group,
            'spliced_utr': getattr(row, 'avg_utr_bp', None),
            'site_type':   getattr(row, 'site_type', 'Unknown'),
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
    lookup     = None
    lookup_path = args.insert_size_lookup
    if lookup_path and lookup_path.lower() != 'none':
        if os.path.exists(lookup_path):
            lookup     = json.load(open(lookup_path))
            prob_array = np.array(lookup['probs'])
            print(f"[assign] Probabilistic assignment: KDE lookup loaded")
            print(f"         {lookup_path}")
            print(f"         median_offset={lookup.get('offset_median', 'N/A')}bp")
        else:
            print(f"[assign] WARNING: lookup file not found: {lookup_path}")
            print(f"         Falling back to shortest-distance assignment.")
            print(f"         Run validate_insert_size.py to generate a lookup table.")
    else:
        print(f"[assign] Shortest-distance assignment (--insert-size-lookup none).")

    # Load exon-space indices (optional — enables Tandem_UTR/ALE correction)
    gene_canonical_exons   = {}
    gene_merged_exons_sorted = {}
    if args.indices_pkl:
        print(f"[assign] Loading exon-space indices: {args.indices_pkl}")
        with open(args.indices_pkl, 'rb') as _f:
            _idx = pickle.load(_f)
        gene_canonical_exons = _idx.get('gene_canonical_exons', {})
        # Convert IntervalTree → sorted list of (start, end_inclusive) for each gene.
        # IntervalTree stores [iv.begin, iv.end) half-open where iv.end = ex.end + 1.
        for _g, _tree in _idx.get('gene_merged_exons', {}).items():
            gene_merged_exons_sorted[_g] = sorted(
                (iv.begin, iv.end - 1) for iv in _tree
            )
        print(f"[assign]   Canonical exon sets: {len(gene_canonical_exons):,} genes")
        print(f"[assign]   Merged exon sets:    {len(gene_merged_exons_sorted):,} genes")
        print(f"[assign]   Exon-space fetch active for Tandem_UTR / ALE sites")
    else:
        print(f"[assign] Exon-space fetch: disabled (no --indices-pkl provided)")

    # Resolve fetch window: explicit --window overrides default 420bp
    if args.window is not None:
        window = args.window
        print(f"[assign] Window: {window}bp (user-specified)")
    else:
        window = 420
        print(f"[assign] Window: {window}bp (default)")

    def score_read(dist):
        """
        Return assignment score for a read at given offset distance.
        Higher score = better assignment.
        KDE mode: probability from empirical distribution
        Fallback:  negative distance (closer = higher score)
        """
        if prob_array is not None:
            if dist > len(prob_array) - 1:
                return 0.0
            return float(prob_array[int(dist)])
        else:
            return -dist  # negative distance: closer = less negative = higher

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

        print(f"[assign] Progress: {current_chrom} | "
              f"Features processed so far: {len(global_soft_counts)}")

        # Initialize tracking dict locally per chromosome to prevent memory bloat
        umi_best_match = {}  # (cb, ub) -> {feature: score, ...}

        for site in sites:
            # ------------------------------------------------------------------
            # Fetch strategy: exon-space walk for Tandem_UTR / ALE;
            # original genomic window for all other site types.
            #
            # Exon-space walk: collects exactly `window` bp of spliced exon
            # space upstream of the PA coord, spanning as many exons as needed.
            # This correctly captures reads that pile up in exon N-1, N-2, etc.
            # when the terminal exon is shorter than the KDE mode (~137 bp).
            # ------------------------------------------------------------------
            site_type = site.get('site_type', 'Unknown')
            exon_list      = None
            use_exon_space = False
            pa2_list       = site.get('pa2_coords') or [site['pos']]

            if site_type in ('Tandem_UTR', 'ALE') and (gene_canonical_exons or gene_merged_exons_sorted):
                exon_list = _get_exon_list_for_site(
                    site['gene'], site['pos'], site_type,
                    gene_canonical_exons, gene_merged_exons_sorted,
                )
                if exon_list:
                    # Fast-path: if the PA exon already has >= window bp upstream of
                    # coord_max (the rightmost PA), no walk needed for any site in the group.
                    pa_cmax = max(pa2_list) if site['strand'] == '+' else min(pa2_list)
                    upstream_in_pa_exon = _pa_exon_upstream_bp(
                        pa_cmax, site['strand'], exon_list
                    )
                    use_exon_space = upstream_in_pa_exon < window

            if use_exon_space:
                # Walk 420bp upstream from the leftmost PA (min on + strand) so all
                # PA sites in the group receive equal upstream coverage.
                # Then add the full PA cluster interval [min_pa, max_pa] so reads
                # between coord_min and coord_max are also fetched.
                walk_anchor = min(pa2_list) if site['strand'] == '+' else max(pa2_list)
                exon_ivs    = _get_exon_space_intervals(walk_anchor, site['strand'], exon_list, window)
                cluster_iv  = (min(pa2_list), max(pa2_list))
                combined    = sorted(exon_ivs + [cluster_iv])
                fetch_intervals = []
                for _s, _e in combined:
                    if fetch_intervals and _s <= fetch_intervals[-1][1] + 1:
                        fetch_intervals[-1] = (fetch_intervals[-1][0], max(fetch_intervals[-1][1], _e))
                    else:
                        fetch_intervals.append((_s, _e))
                if not fetch_intervals:
                    use_exon_space = False  # no exons found → fall through to genomic

            if not use_exon_space:
                cmin = site.get('cmin', site['pos'])
                cmax = site.get('cmax', site['pos'])
                if site['strand'] == '+':
                    fetch_intervals = [(max(0, cmin - window), cmax)]
                else:
                    fetch_intervals = [(cmin, cmax + window)]

            for _fetch_start, _fetch_end in fetch_intervals:
                # exon_list intervals use 0-based inclusive end → pysam needs exclusive end
                _pysam_end = _fetch_end + 1 if use_exon_space else _fetch_end
                for read in bam.fetch(current_chrom, _fetch_start, _pysam_end):
                    if read.is_secondary or read.is_supplementary or read.is_unmapped:
                        continue
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
                    dists = []
                    for p in pa2_list:
                        if site['strand'] == '-':
                            if read_3_prime >= p:    # read upstream on neg strand (higher coord)
                                dists.append(
                                    _exon_space_dist(read_3_prime, p, exon_list, '-')
                                    if use_exon_space else read_3_prime - p
                                )
                        else:
                            if read_3_prime <= p:    # read upstream on pos strand (lower coord)
                                dists.append(
                                    _exon_space_dist(read_3_prime, p, exon_list, '+')
                                    if use_exon_space else p - read_3_prime
                                )

                    if not dists:
                        continue                     # read is downstream of all PA sites

                    dist = min(dists)

                    if dist <= window:
                        ub      = read.get_tag("UB")
                        umi_key = (cb, ub)
                        # user_label (e.g. Secreted/Membrane for IGH genes) is stored
                        # in panel CSV as metadata only — NOT embedded in feature name.
                        # Keeping feature names as Gene_G0 / Gene_G1 ensures that
                        # integrate_samples.py's panel join (keyed on gene + '_G' + group)
                        # correctly matches IGH features and populates utr_source /
                        # avg_utr_bp from the panel CSV.
                        feature = f"{site['gene']}_G{site['group']}"

                        # max KDE over all PA sites in group: a read 150bp upstream of
                        # any PA site scores the same regardless of which site it came from
                        s = max(score_read(d) for d in dists)

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
                                    'pa2_coords': site.get('pa2_coords', [site['pos']]),
                                    'offset_bp':  dist,
                                    'score':      round(s, 6),
                                    'method':     ('KDE_exon' if use_exon_space else 'KDE')
                                                  if prob_array is not None
                                                  else ('dist_exon' if use_exon_space else 'shortest_dist'),
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
        gc.collect()

    # Debug output
    if debug_gene and umi_debug_info:
        debug_path = f"{args.out}_debug_{debug_gene}.tsv"
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

    if global_soft_counts:
        import scipy.sparse

        obs_names = sorted(set(cb   for cb,   _ in global_soft_counts))
        var_names = sorted(set(feat for _,  feat in global_soft_counts))

        obs_idx = {cb:   i for i, cb   in enumerate(obs_names)}
        var_idx = {feat: j for j, feat in enumerate(var_names)}

        rows, cols, data = [], [], []
        for (cb, feat), val in global_soft_counts.items():
            rows.append(obs_idx[cb])
            cols.append(var_idx[feat])
            data.append(val)

        mat = scipy.sparse.csr_matrix(
            (data, (rows, cols)),
            shape=(len(obs_names), len(var_names)),
            dtype='float32',
        )

        out_base = args.out
        scipy.sparse.save_npz(f"{out_base}_matrix.npz", mat)
        with open(f"{out_base}_obs.txt", 'w') as f:
            f.write('\n'.join(obs_names))
        with open(f"{out_base}_var.txt", 'w') as f:
            f.write('\n'.join(var_names))

        print(f"✅ Saved: {out_base}_matrix.npz / _obs.txt / _var.txt")
        print(f"   Cells:    {len(obs_names):,}")
        print(f"   Features: {len(var_names):,}")
        print(f"   UMIs:     {global_umi_count:,} (soft-assigned → fractional counts)")

        # Panel composition summary (one JSON per sample for pie charts / QC)
        # panel_exp: expressed features only, ordered to match mat columns
        panel['feature_id'] = panel['gene'] + '_G' + panel['polyA_group'].astype(str)
        panel_exp = panel.set_index('feature_id').loc[var_names].copy()
        pcols = set(panel_exp.columns)

        has_all = all(c in pcols for c in (
            'pc_fraction', 'ig_fraction', 'tr_fraction',
            'ig_pseudo_fraction', 'tr_pseudo_fraction',
        ))

        def _feat(mask):
            return int(mask.sum())

        def _umi(mask):
            return float(mat[:, mask.values].sum())

        if has_all:
            pc_mask  = ((panel_exp['pc_fraction'] > 0) &
                        (panel_exp['ig_fraction'] == 0) &
                        (panel_exp['tr_fraction'] == 0) &
                        (panel_exp['ig_pseudo_fraction'] == 0) &
                        (panel_exp['tr_pseudo_fraction'] == 0))
            ig_mask  = panel_exp['ig_fraction']        > 0
            tr_mask  = panel_exp['tr_fraction']        > 0
            igp_mask = panel_exp['ig_pseudo_fraction'] > 0
            trp_mask = panel_exp['tr_pseudo_fraction'] > 0
            nc_mask  = panel_exp['pc_fraction'] == 0

            biotype_features = {
                'protein_coding': _feat(pc_mask),
                'ig':             _feat(ig_mask),
                'tr':             _feat(tr_mask),
                'ig_pseudogene':  _feat(igp_mask),
                'tr_pseudogene':  _feat(trp_mask),
                'non_coding':     _feat(nc_mask),
            }
            biotype_umis = {
                'protein_coding': _umi(pc_mask),
                'ig':             _umi(ig_mask),
                'tr':             _umi(tr_mask),
                'ig_pseudogene':  _umi(igp_mask),
                'tr_pseudogene':  _umi(trp_mask),
                'non_coding':     _umi(nc_mask),
            }
        else:
            biotype_features = {}
            biotype_umis     = {}
            if 'pc_fraction' in pcols:
                biotype_features['protein_coding'] = _feat(panel_exp['pc_fraction'] > 0)
                biotype_features['non_coding']     = _feat(panel_exp['pc_fraction'] == 0)
                biotype_umis['protein_coding']     = _umi(panel_exp['pc_fraction'] > 0)
                biotype_umis['non_coding']         = _umi(panel_exp['pc_fraction'] == 0)

        if 'nmd_fraction' in pcols:
            nmd_mask = panel_exp['nmd_fraction'] > 0
            biotype_features['nmd'] = _feat(nmd_mask)
            biotype_umis['nmd']     = _umi(nmd_mask)

        summary = {
            'sample':     os.path.basename(out_base),
            'cells':      len(obs_names),
            'features':   len(var_names),
            'total_umis': float(mat.sum()),

            'site_type_features': (panel_exp['site_type'].value_counts().to_dict()
                                   if 'site_type' in pcols else {}),
            'site_type_umis': ({
                k: float(mat[:, (panel_exp['site_type'] == k).values].sum())
                for k in panel_exp['site_type'].unique()
            } if 'site_type' in pcols else {}),

            'biotype_features': biotype_features,
            'biotype_umis':     biotype_umis,

            'pas_motif_strength_features': (
                panel_exp['pas_motif_strength'].value_counts().to_dict()
                if 'pas_motif_strength' in pcols else {}),
            'pas_motif_strength_umis': ({
                k: float(mat[:, (panel_exp['pas_motif_strength'] == k).values].sum())
                for k in panel_exp['pas_motif_strength'].unique()
            } if 'pas_motif_strength' in pcols else {}),

            'utr_confidence_features': (
                panel_exp['utr_confidence'].value_counts(dropna=False)
                .rename(lambda x: str(x)).to_dict()
                if 'utr_confidence' in pcols else {}),
            'utr_confidence_umis': ({
                str(k): float(mat[:, (panel_exp['utr_confidence'].fillna('nan') == str(k)).values].sum())
                for k in panel_exp['utr_confidence'].fillna('nan').unique()
            } if 'utr_confidence' in pcols else {}),
        }

        summary_path = f"{out_base}_summary.json"
        with open(summary_path, 'w') as fh:
            json.dump(summary, fh, indent=2)
        print(f"[assign] Summary saved → {summary_path}")
    else:
        print(f"⚠️  No reads assigned for {args.bam}")


if __name__ == "__main__":
    main()