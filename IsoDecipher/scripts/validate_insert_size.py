#!/usr/bin/env python3
# ==============================================================================
# Copyright (c) 2026 IsoMatrix Suite. All rights reserved.
# ==============================================================================
"""
IsoDecipher: validate_insert_size.py
-------------------------------------
Estimate empirical insert size KDE from singleton polyA sites in the panel.

Singleton sites (genes with exactly 1 polyA group) are ideal calibration
anchors because their reads cannot be confused with a competing nearby site.
This matches the original PolyASite v2.0 single-PAS gene design, but now
uses the panel's own rep_coord (which is the v2.0 coordinate for PA2/GTF+PA2
sites), ensuring coord-system consistency with assign_reads.py.

Usage:
  python validate_insert_size.py \\
      --bam   /path/to/sample.bam \\
      --panel results/panel/panel_features_global_v1.1.csv \\
      --out   results/exp97_insert_size.png \\
      --n-reads 1000000
"""

import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


def parse_args():
    parser = argparse.ArgumentParser(
        description="Estimate insert size KDE from singleton panel sites"
    )
    parser.add_argument("--bam",        required=True)
    parser.add_argument("--panel",      required=True,
                        help="panel_features_global_v1.1.csv (with is_singleton column)")
    parser.add_argument("--out",        required=True,
                        help="Output PNG path (lookup JSON saved alongside)")
    parser.add_argument("--n-reads",    type=int, default=1_000_000,
                        help="Max reads to sample (default: 1,000,000)")
    parser.add_argument("--window",     type=int, default=500,
                        help="Fetch window around each site (default: 500bp)")
    parser.add_argument("--n-sites",    type=int, default=200,
                        help="Top N sites by avg_utr_bp (default: 200)")
    parser.add_argument("--no-singleton", action="store_true", default=False,
                        help="Disable singleton filter (default: singleton ON)")
    parser.add_argument("--max-coord-spread", type=int, default=10,
                        help="Exclude anchor sites whose panel group internally spans "
                             "more than this many bp (coord_max - coord_min). Prevents "
                             "within-group micro-heterogeneity (multiple pa2_coords merged "
                             "into one group) from contaminating the KDE kernel with extra "
                             "spread that isn't true fragment-length variance. "
                             "Default: 10bp. Set to -1 to disable.")
    parser.add_argument("--min-reads", type=int, default=2000,
                        help="Minimum reads needed for a stable KDE. If the initial "
                             "--n-sites yields fewer, automatically retries with a larger "
                             "n_sites (relaxing top-avg_utr_bp selectivity) until this is "
                             "met or --max-n-sites is reached. Low-depth samples (e.g. "
                             "exp105) can collect too few reads at n_sites=50 for a stable "
                             "KDE; but going past ~n_sites=150 risks pulling in "
                             "tissue-restricted/off-target genes (empirically: keratin/LCE/ "
                             "skin genes start appearing around rank 150-200) that distort "
                             "the kernel shape for non-matching cell types. Default: 2000.")
    parser.add_argument("--smooth-factor", type=float, default=1.5,
                        help="Multiplier applied on top of Scott's-rule bandwidth. "
                             "1.0 = pure Scott's rule, which under-smooths for "
                             "low-read-count samples (n_reads ~3000-5000 gives a visibly "
                             "wiggly curve vs the clean single-mode shape seen in IGV pileups). "
                             "Default: 1.5.")
    parser.add_argument("--max-n-sites", type=int, default=150,
                        help="Ceiling for automatic --n-sites escalation when read count "
                             "is below --min-reads. Mostly a formality now: select_sites() "
                             "requires avg_utr_bp >= --window, which structurally caps the "
                             "eligible pool (~106 sites at window=500bp) regardless of how "
                             "high n_sites is requested — the earlier failure mode (keratin/ "
                             "LCE/skin genes at rank 150-200 distorting the KDE) was these "
                             "sites' UTRs being shorter than the fetch window, not their "
                             "rank per se. Default: 150.")
    return parser.parse_args()


def select_sites(panel_path, n_sites, singleton_only=False, max_coord_spread=10,
                  window=500):
    """
    Select the top N Tandem_UTR sites by avg_utr_bp for KDE calibration.

    avg_utr_bp is NOT a proxy for housekeeping/broad expression — it is purely the
    annotated 3'UTR length, which has no inherent relationship to how widely a gene
    is expressed (rank 1-50 by this sort includes olfactory receptors, testis-only
    genes, developmental homeobox genes — the opposite of "housekeeping"). What it
    actually buys us: requiring avg_utr_bp >= window (see below) guarantees the
    fetch window used to collect offsets stays entirely inside the annotated,
    spliced 3'UTR — never spilling into an intron, a neighboring gene, or (for
    compact paralog-dense families like keratins) a homologous nearby gene that
    would cross-contaminate the offset measurement. Requiring Tandem_UTR +
    PA2/GTF+PA2 ensures reads are real 3'-end signal at verified cleavage sites.

    avg_utr_bp >= window filter: empirically, ranks beyond ~107 (at window=500bp)
    drop below the window size — exactly where tissue-restricted/paralog-dense
    genes (keratin/LCE/CST family members) start appearing and distorting the KDE
    on low-depth samples (verified: n_sites=200 on exp105 caused a +57bp median
    jump vs n_sites<=150, driven by these short-UTR sites). Filtering directly on
    avg_utr_bp >= window replaces the earlier empirical n_sites<=150 heuristic
    with a structural guarantee, and naturally caps how large n_sites can safely
    go for a given window.

    coord_spread filter (default max 10bp): is_singleton only guarantees no OTHER
    polyA_group exists for this gene — it says nothing about whether THIS group
    itself is an internal merge of multiple nearby pa2_coords (up to 150bp apart).
    An anchor site built from a merged multi-coordinate group would bake that
    within-group micro-heterogeneity into the KDE as if it were fragment-length
    noise, inflating/distorting the kernel used for matched-filter peak calling.
    Sorted by avg_utr_bp descending, fully deterministic across samples.
    """
    panel = pd.read_csv(panel_path, low_memory=False)

    mask = (
        (panel['panel_source'].isin(['PA2', 'GTF+PA2'])) &
        (panel['site_type'] == 'Tandem_UTR') &
        (panel['avg_utr_bp'].notna()) &
        (panel['avg_utr_bp'] >= window)
    )
    print(f"[validate] avg_utr_bp >= window ({window}bp) filter: "
          f"{mask.sum():,} eligible sites (guarantees fetch window stays inside UTR)")

    if singleton_only:
        if 'is_singleton' not in panel.columns:
            print("[validate] is_singleton column not found — computing from panel")
            groups_per_gene = panel.groupby('gene')['polyA_group'].max() + 1
            singleton_genes = set(groups_per_gene[groups_per_gene == 1].index)
            panel['is_singleton'] = panel['gene'].isin(singleton_genes)
        mask = mask & (panel['is_singleton'] == True)
        print(f"[validate] Singleton filter: ON")
    else:
        print(f"[validate] Singleton filter: OFF (all Tandem_UTR PA2 sites)")

    if max_coord_spread >= 0:
        if 'coord_spread' in panel.columns:
            before = mask.sum()
            mask = mask & (panel['coord_spread'] <= max_coord_spread)
            print(f"[validate] coord_spread filter: <= {max_coord_spread}bp "
                  f"(excluded {before - mask.sum():,} internally-merged groups)")
        else:
            print(f"[validate] WARNING: coord_spread column not found — skipping filter")

    sites = panel[mask].sort_values('avg_utr_bp', ascending=False).head(n_sites)

    print(f"[validate] Using {len(sites):,} Tandem_UTR sites (top by avg_utr_bp)")
    if len(sites) < 50:
        print(f"[validate] WARNING: only {len(sites)} sites — results may be noisy")
    return sites


def collect_offsets(sites, bam_path, n_reads, window):
    """Collect read-to-site offsets from BAM."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    offsets = []
    n_total = 0
    max_per_site = max(1, n_reads // len(sites))  # equal budget — prevents any single site dominating KDE

    for row in sites.itertuples(index=False):
        if n_total >= n_reads:
            break

        chrom  = str(row.chrom)
        if chrom not in bam.references:
            chrom = f'chr{chrom}'
        if chrom not in bam.references:
            continue

        pos    = int(row.rep_coord)
        strand = row.strand

        if strand == '+':
            fetch_start = max(0, pos - window)
            fetch_end   = pos + 50
        else:
            fetch_start = max(0, pos - 50)
            fetch_end   = pos + window

        n_this_site = 0
        for read in bam.fetch(chrom, fetch_start, fetch_end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not (read.has_tag("CB") and read.has_tag("UB")):
                continue
            read_strand = '-' if read.is_reverse else '+'
            if read_strand != strand:
                continue

            if strand == '+':
                cleavage = read.reference_end
                offset   = pos - cleavage
            else:
                cleavage = read.reference_start
                offset   = cleavage - pos

            if 0 <= offset <= window:
                offsets.append(offset)
                n_total    += 1
                n_this_site += 1

            if n_this_site >= max_per_site or n_total >= n_reads:
                break

    bam.close()
    return np.array(offsets)


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    # Select calibration sites + collect offsets, escalating n_sites if too few
    # reads come back (low-depth samples). The avg_utr_bp >= window filter in
    # select_sites() structurally caps how far this can safely go — once the
    # eligible pool (sites with UTR long enough to contain the fetch window) is
    # exhausted, further escalation just returns the same site set.
    n_sites_try = args.n_sites
    while True:
        sites = select_sites(args.panel, n_sites_try, singleton_only=not args.no_singleton,
                              max_coord_spread=args.max_coord_spread, window=args.window)
        print(f"[validate] Sampling up to {args.n_reads:,} reads from BAM...")
        arr = collect_offsets(sites, args.bam, args.n_reads, args.window)

        if len(arr) >= args.min_reads or n_sites_try >= args.max_n_sites:
            if len(arr) < args.min_reads:
                print(f"[validate] WARNING: only {len(arr):,} reads even at "
                      f"n_sites={n_sites_try} (target {args.min_reads:,}). Sample is "
                      f"likely too shallow for a stable KDE — inspect the offset "
                      f"distribution below before trusting this lookup.")
            elif n_sites_try > args.n_sites:
                print(f"[validate] NOTE: n_sites escalated from {args.n_sites} to "
                      f"{n_sites_try} to reach {args.min_reads:,} reads. Sites beyond "
                      f"~rank 150 can be tissue-restricted/off-target for this cell "
                      f"type — check the site list if the KDE shape looks unusual.")
            break

        n_sites_next = min(args.max_n_sites, int(n_sites_try * 1.5) + 25)
        print(f"[validate] Only {len(arr):,} reads (< {args.min_reads:,} threshold) at "
              f"n_sites={n_sites_try} — retrying with n_sites={n_sites_next}")
        n_sites_try = n_sites_next

    print(f"\n[validate] Collected {len(arr):,} reads")
    print(f"[validate] Offset distribution:")
    print(f"  min:    {arr.min():.0f}bp")
    print(f"  p5:     {np.percentile(arr, 5):.0f}bp")
    print(f"  median: {np.median(arr):.0f}bp")
    print(f"  mean:   {arr.mean():.0f}bp")
    print(f"  p95:    {np.percentile(arr, 95):.0f}bp")
    print(f"  max:    {arr.max():.0f}bp")
    print(f"  Reads within 420bp: {(arr <= 420).sum():,} ({(arr<=420).mean():.1%})")

    # KDE fit — Scott's rule, optionally widened for extra smoothing (low-read-count
    # samples give a visibly wiggly curve at bw_method='scott' alone).
    kde = gaussian_kde(arr, bw_method='scott')
    if args.smooth_factor != 1.0:
        kde.set_bandwidth(bw_method=kde.factor * args.smooth_factor)
    x_eval     = np.arange(0, 501, 1)
    kde_probs  = np.maximum(kde.evaluate(x_eval), 0)

    # Save lookup JSON (0–420bp, normalized)
    prob_420 = kde_probs[:421].copy()
    prob_420 = prob_420 / prob_420.sum()

    # Summit = argmax of the (smoothed) KDE, not the median. For a right-skewed
    # distribution (hard cutoff at 0bp, long tail) median > mode — downstream
    # consumers that want "the most likely offset" (e.g. master_peak_merge.py's
    # pa_coord = rep_coord - KDE_mode) need the summit, not the median.
    offset_summit = float(np.argmax(prob_420))
    print(f"  summit: {offset_summit:.0f}bp  (mode of KDE — use this, not median, "
          f"as the single-point offset correction)")

    lookup_path = args.out.replace('.png', '_lookup.json')
    lookup = {
        'offsets':             list(range(421)),
        'probs':               prob_420.tolist(),
        'offset_summit':       offset_summit,
        'offset_median':       float(np.median(arr)),
        'offset_p5':           float(np.percentile(arr, 5)),
        'offset_p95':          float(np.percentile(arr, 95)),
        'n_reads':             int(len(arr)),
        'n_sites':             int(len(sites)),
        'bw_method':           'scott',
        'bw_factor':           float(kde.factor),
        'smooth_factor':       float(args.smooth_factor),
    }
    with open(lookup_path, 'w') as f:
        json.dump(lookup, f, indent=2)
    print(f"\n[validate] Lookup saved → {lookup_path}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    site_label = "singleton sites" if not args.no_singleton else "all Tandem_UTR PA2 sites"
    fig.suptitle(
        f'Insert Size Distribution ({site_label}, n={len(arr):,} reads)',
        fontsize=12
    )

    # Left: histogram + KDE
    ax = axes[0]
    ax.hist(arr, bins=100, range=(0, 500), color='steelblue',
            alpha=0.5, density=True, label='Empirical')
    ax.plot(x_eval, kde_probs, color='purple', linewidth=2.5,
            label=f'KDE (bw={kde.factor:.4f})')
    ax.axvline(offset_summit, color='black', linestyle='--',
               label=f'Summit: {offset_summit:.0f}bp')
    ax.axvline(np.median(arr), color='dimgray', linestyle='-.', alpha=0.7,
               label=f'Median: {np.median(arr):.0f}bp')
    ax.axvline(np.percentile(arr, 95), color='gray', linestyle=':',
               label=f'P95: {np.percentile(arr, 95):.0f}bp')
    ax.set_xlabel('Offset from polyA site (bp)')
    ax.set_ylabel('Density')
    ax.set_title('Histogram + KDE')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 500)

    # Right: CDF
    ax = axes[1]
    sorted_arr = np.sort(arr)
    cdf_emp = np.arange(1, len(sorted_arr)+1) / len(sorted_arr)
    ax.plot(sorted_arr, cdf_emp, 'steelblue', linewidth=2, label='Empirical CDF')
    kde_cdf = np.cumsum(kde_probs) / kde_probs.sum()
    ax.plot(x_eval, kde_cdf, color='purple', linewidth=1.5,
            linestyle='--', label='KDE CDF')
    ax.axvline(420, color='gray', linestyle=':', label='420bp window')
    ax.set_xlabel('Offset (bp)')
    ax.set_ylabel('Cumulative fraction')
    ax.set_title('CDF')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    plt.close()
    print(f"[validate] Plot saved → {args.out}")

    print(f"  --insert-size-lookup {lookup_path}")


if __name__ == '__main__':
    main()
