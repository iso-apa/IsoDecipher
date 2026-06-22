#!/usr/bin/env python3
"""
Validate 10x library insert size distribution from known canonical polyA sites.

Uses known 3'end sites from IsoDecipher panel to measure actual
cleavage_coord offset distribution. This validates:
1. The 350bp window assumption
2. The actual peak of the offset distribution
3. Whether probabilistic assignment would help

Usage:
python validate_insert_size.py \
    --bam   /Volumes/Lexar/bam/Breast_Cancer_3p_possorted_genome_bam.bam \
    --panel results/panel_features.csv \
    --out   ~/Desktop/insert_size_dist.png \
    --n-reads 200000
"""

import argparse
import numpy as np
import pysam
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam",     required=True)
    parser.add_argument("--panel",   required=True)
    parser.add_argument("--out",     default=os.path.expanduser("~/Desktop/insert_size_dist.png"))
    parser.add_argument("--n-reads", type=int, default=200_000,
                        help="Max reads to sample (default: 200,000)")
    parser.add_argument("--window",  type=int, default=500,
                        help="Fetch window around each site (default: 500bp)")
    return parser.parse_args()


def main():
    import os
    args = parse_args()

    panel = pd.read_csv(args.panel)

    if 'avg_spliced_utr' in panel.columns:
        # If using PolyASite2 single-PAS panel (all avg_spliced_utr == 999)
        # use random sample to avoid bias
        unique_vals = panel['avg_spliced_utr'].nunique()
        if unique_vals == 1:
            top_sites = panel.sample(n=min(200, len(panel)), random_state=42)
            print(f"[validate] PolyASite2 mode: random sample {len(top_sites)} sites")
        else:
            top_sites = panel.sort_values('avg_spliced_utr', ascending=False).head(200)
            print(f"[validate] Panel mode: top {len(top_sites)} sites by UTR length")
    else:
        op_sites = panel.head(200)

    print(f"[validate] Using {len(top_sites)} canonical sites")
    print(f"[validate] Sampling up to {args.n_reads:,} reads")

    bam    = pysam.AlignmentFile(args.bam, "rb")
    offsets_plus  = []  # +strand offsets
    offsets_minus = []  # -strand offsets
    n_total = 0

    for row in top_sites.itertuples(index=False):
        if n_total >= args.n_reads:
            break

        chrom  = row.chrom
        # Handle chr prefix mismatch (panel: '19', BAM: 'chr19')
        if chrom not in bam.references:
            chrom = f'chr{chrom}'
        if chrom not in bam.references:
            continue

        pos    = int(row.rep_coord)
        strand = row.strand

        # Fetch wider window to see full distribution
        if strand == '+':
            fetch_start = max(0, pos - args.window)
            fetch_end   = pos + 50  # small buffer on right
        else:
            fetch_start = pos - 50
            fetch_end   = pos + args.window

        for read in bam.fetch(chrom, fetch_start, fetch_end):
            if read.is_unmapped or read.is_secondary:
                continue
            if not (read.has_tag("CB") and read.has_tag("UB")):
                continue

            read_strand = '-' if read.is_reverse else '+'
            if read_strand != strand:
                continue

            # cleavage_coord (same as bam_to_parquet)
            if strand == '+':
                cleavage = read.reference_end
                offset   = pos - cleavage  # how far UPSTREAM of polyA site
            else:
                cleavage = read.reference_start
                offset   = cleavage - pos  # how far DOWNSTREAM of polyA site

            # Only keep reads within reasonable range
            if 0 <= offset <= args.window:
                if strand == '+':
                    offsets_plus.append(offset)
                else:
                    offsets_minus.append(offset)
                n_total += 1

            if n_total >= args.n_reads:
                break

    bam.close()

    print(f"\n[validate] Collected {n_total:,} reads")
    print(f"  +strand: {len(offsets_plus):,}")
    print(f"  -strand: {len(offsets_minus):,}")

    all_offsets = offsets_plus + offsets_minus
    arr = np.array(all_offsets)

    print(f"\n[validate] Offset distribution (cleavage_coord distance from polyA site):")
    print(f"  min:    {arr.min():.0f}bp")
    print(f"  p5:     {np.percentile(arr, 5):.0f}bp")
    print(f"  p25:    {np.percentile(arr, 25):.0f}bp")
    print(f"  median: {np.median(arr):.0f}bp")
    print(f"  mean:   {arr.mean():.0f}bp")
    print(f"  p75:    {np.percentile(arr, 75):.0f}bp")
    print(f"  p95:    {np.percentile(arr, 95):.0f}bp")
    print(f"  max:    {arr.max():.0f}bp")
    print(f"\n  Reads within 350bp: {(arr <= 350).sum():,} ({(arr<=350).mean():.1%})")
    print(f"  Reads within 200bp: {(arr <= 200).sum():,} ({(arr<=200).mean():.1%})")
    print(f"  Reads within 100bp: {(arr <= 100).sum():,} ({(arr<=100).mean():.1%})")

    # Distribution fitting comparison
    print(f"\n[validate] Fitting distributions...")

    from scipy.stats import lognorm, gamma, norm, kstest
    import warnings

    arr_pos = arr[arr > 0]  # log-normal needs positive values

    fits = {}

    # 1. Normal
    mu_n, sigma_n = norm.fit(arr)
    ks_n, p_n = kstest(arr, 'norm', args=(mu_n, sigma_n))
    fits['Normal'] = dict(params=(mu_n, sigma_n), ks=ks_n, p=p_n,
                          pdf=lambda x: norm.pdf(x, mu_n, sigma_n))

    # 2. Log-normal
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sigma_ln, loc_ln, scale_ln = lognorm.fit(arr_pos, floc=0)
    ks_ln, p_ln = kstest(arr_pos, 'lognorm',
                         args=(sigma_ln, loc_ln, scale_ln))
    fits['Log-normal'] = dict(
        params=(sigma_ln, loc_ln, scale_ln), ks=ks_ln, p=p_ln,
        pdf=lambda x: lognorm.pdf(x, sigma_ln, loc_ln, scale_ln))

    # 3. Gamma
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a_g, loc_g, scale_g = gamma.fit(arr_pos, floc=0)
    ks_g, p_g = kstest(arr_pos, 'gamma', args=(a_g, loc_g, scale_g))
    fits['Gamma'] = dict(
        params=(a_g, loc_g, scale_g), ks=ks_g, p=p_g,
        pdf=lambda x: gamma.pdf(x, a_g, loc_g, scale_g))

    # KDE fit (empirical, no distribution assumption)
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(arr, bw_method='scott')
    x_eval = np.arange(0, 501, 1)
    kde_probs = kde.evaluate(x_eval)
    kde_probs = np.maximum(kde_probs, 0)  # no negative

    # KS stat for KDE (manual)
    from scipy.stats import kstest
    ks_kde = kstest(arr, lambda x: np.interp(
        x, x_eval, np.cumsum(kde_probs) / kde_probs.sum())).statistic
    print(f"  {'KDE':<15} {ks_kde:<12.4f} {'N/A':<12} bw={kde.factor:.4f}")
    print(f"\n  Best fit: KDE (empirical, no assumption)")

    # Save lookup table
    lookup_path = args.out.replace('.png', '_lookup.json')
    import json
    # Normalize to probability (sum to 1 over 0-420bp)
    prob_420 = kde_probs[:421].copy()
    prob_420 = prob_420 / prob_420.sum()
    lookup = {
        'offsets':    list(range(421)),
        'probs':      prob_420.tolist(),
        'offset_mu':  float(np.median(arr)),
        'offset_p5':  float(np.percentile(arr, 5)),
        'offset_p95': float(np.percentile(arr, 95)),
        'recommended_window': int(np.percentile(arr, 95)) + 20,
        'n_reads':    len(arr),
    }
    with open(lookup_path, 'w') as f:
        json.dump(lookup, f, indent=2)
    print(f"  Lookup table saved → {lookup_path}")
    print(f"  Usage in assign_reads.py:")
    print(f"    import json, numpy as np")
    print(f"    lookup = json.load(open('{lookup_path}'))")
    print(f"    prob_array = np.array(lookup['probs'])")
    print(f"    def insert_prob(offset):")
    print(f"        if offset < 0 or offset > 420: return 0.0")
    print(f"        return float(prob_array[int(offset)])")

    # Plot with fits
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'10x 3\' Library Insert Size Distribution\n'
                 f'(n={len(all_offsets):,} reads from {len(top_sites)} canonical sites)',
                 fontsize=13)

    x = np.linspace(0, 500, 500)
    colors = {'Normal': 'red', 'Log-normal': 'green', 'Gamma': 'orange'}

    # Left: histogram + fits
    ax = axes[0]
    ax.hist(arr, bins=100, range=(0, 500), color='steelblue',
            alpha=0.5, density=True, label='Empirical')
    for name, fit in fits.items():
        ax.plot(x, fit['pdf'](x), color=colors[name],
                linewidth=1.5, linestyle='--', label=f"{name} (KS={fit['ks']:.3f})")
    # KDE
    ax.plot(x_eval, kde_probs, color='purple', linewidth=2.5,
            label=f'KDE (KS={ks_kde:.3f})')
    ax.axvline(np.median(arr), color='black', linestyle='--',
               label=f'Median: {np.median(arr):.0f}bp')
    ax.set_xlabel('Offset from polyA site (bp)')
    ax.set_ylabel('Density')
    ax.set_title('Histogram + Distribution Fits')
    ax.legend(fontsize=7)
    ax.set_xlim(0, 500)

    # Middle: Q-Q plots
    ax = axes[1]
    from scipy.stats import probplot
    for name, fit in fits.items():
        if name == 'Normal':
            res = probplot(arr, dist='norm', sparams=fit['params'])
        elif name == 'Log-normal':
            res = probplot(arr_pos, dist='lognorm', sparams=fit['params'])
        else:
            res = probplot(arr_pos, dist='gamma', sparams=fit['params'])
        ax.plot(res[0][0], res[0][1], '.', alpha=0.3,
                color=colors[name], markersize=2, label=name)
        ax.plot(res[0][0], res[1][1] + res[1][0]*res[0][0],
                '-', color=colors[name], linewidth=1.5)
    ax.set_xlabel('Theoretical quantiles')
    ax.set_ylabel('Sample quantiles')
    ax.set_title('Q-Q Plot (closer to diagonal = better fit)')
    ax.legend(fontsize=8)

    # Right: CDF comparison
    ax = axes[2]
    sorted_arr = np.sort(arr)
    cdf_emp = np.arange(1, len(sorted_arr)+1) / len(sorted_arr)
    ax.plot(sorted_arr, cdf_emp, 'steelblue', linewidth=2, label='Empirical')
    for name, fit in fits.items():
        if name == 'Normal':
            from scipy.stats import norm as norm_dist
            cdf_fit = norm_dist.cdf(x, *fit['params'])
        elif name == 'Log-normal':
            cdf_fit = lognorm.cdf(x, *fit['params'])
        else:
            cdf_fit = gamma.cdf(x, *fit['params'])
        ax.plot(x, cdf_fit, color=colors[name], linewidth=1.5,
                linestyle='--', label=name)
    ax.axvline(350, color='gray', linestyle=':', label='350bp window')
    ax.set_xlabel('Offset (bp)')
    ax.set_ylabel('Cumulative fraction')
    ax.set_title('CDF Comparison')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    print(f"\n[validate] Plot saved: {args.out}")
    plt.close()

    # Final recommendation
    p95 = np.percentile(arr, 95)
    print(f"\n[validate] Recommendation:")
    print(f"  Use KDE empirical distribution (best fit, no assumption)")
    print(f"  P95 offset = {p95:.0f}bp → recommended window = {int(p95)+20}bp")
    print(f"  Lookup table: {lookup_path}")
    print(f"  Load in assign_reads.py with --insert-size-lookup {lookup_path}")


if __name__ == '__main__':
    import os
    main()
