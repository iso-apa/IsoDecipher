#!/usr/bin/env python3
# ==============================================================================
# Copyright (c) 2026 IsoMatrix Suite. All rights reserved.
# ==============================================================================
"""
IsoDecipher: prep_sample_panel.py
----------------------------------
Pre-filter the global panel for single or multiple samples before assign_reads.py.

Three filters applied:
  1. Remove singleton genes (only 1 polyA group — no APA analysis)
  2. Keep only genes expressed in ANY of the provided samples (Union)
  3. Keep only standard chromosomes (chr1-22, chrX, chrY, chrM)

Output: a batch-specific active panel CSV ready for assign_reads.py.
"""

import argparse
import os
import pandas as pd
import scanpy as sc

STANDARD_CHROMS = {
    str(i) for i in range(1, 23)
} | {'X', 'Y', 'MT', 'M'} | {
    f'chr{i}' for i in range(1, 23)
} | {'chrX', 'chrY', 'chrM', 'chrMT'}


def parse_args():
    p = argparse.ArgumentParser(
        description="Prepare batch-specific active panel for assign_reads.py"
    )
    p.add_argument("--global-panel", required=True,
                   help="Global panel CSV (panel_features_global_v1.1.csv)")
    p.add_argument("--matrix-dirs",  required=True, nargs='+',
                   help="List of Cell Ranger filtered_feature_bc_matrix/ directories")
    p.add_argument("--out",          required=True,
                   help="Output batch-specific panel CSV")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    # Load global panel
    print(f"[prep] Loading global panel: {args.global_panel}")
    panel = pd.read_csv(args.global_panel, low_memory=False)
    print(f"[prep]   {len(panel):,} features, {panel['gene'].nunique():,} genes")

    # Filter 1: remove singleton genes
    print("[prep] Filter 1: removing singleton genes...")
    if 'is_singleton' in panel.columns:
        before = len(panel)
        panel  = panel[panel['is_singleton'] == False].copy()
        print(f"[prep]   Removed {before - len(panel):,} singleton features "
              f"→ {len(panel):,} remaining")
    else:
        groups_per_gene   = panel.groupby('gene')['polyA_group'].max() + 1
        multi_group_genes = set(groups_per_gene[groups_per_gene >= 2].index)
        before = len(panel)
        panel  = panel[panel['gene'].isin(multi_group_genes)].copy()
        print(f"[prep]   Removed {before - len(panel):,} singleton features "
              f"→ {len(panel):,} remaining")

    # Filter 2: union of expressed genes across all samples
    print(f"[prep] Filter 2: finding expressed genes across {len(args.matrix_dirs)} samples...")
    all_expressed_genes = set()

    for i, matrix_dir in enumerate(args.matrix_dirs, 1):
        print(f"[prep]   Reading sample {i}/{len(args.matrix_dirs)}: {matrix_dir}")
        adata = sc.read_10x_mtx(
            matrix_dir,
            var_names='gene_symbols',
            cache=False,
            gex_only=True,
        )
        expressed = set(adata.var_names[adata.X.sum(axis=0).A1 > 0])
        print(f"[prep]     → {len(expressed):,} expressed genes")
        all_expressed_genes.update(expressed)

    print(f"[prep]   Total unique expressed genes across all samples: {len(all_expressed_genes):,}")

    before = len(panel)
    panel  = panel[panel['gene'].str.upper().isin(
        {g.upper() for g in all_expressed_genes}
    )].copy()
    print(f"[prep]   Removed {before - len(panel):,} unexpressed features "
          f"→ {len(panel):,} remaining")

    # Filter 3: standard chromosomes
    print("[prep] Filter 3: standard chromosomes only...")
    before = len(panel)
    panel  = panel[panel['chrom'].astype(str).isin(STANDARD_CHROMS)].copy()
    print(f"[prep]   Removed {before - len(panel):,} non-standard chrom features "
          f"→ {len(panel):,} remaining")

    # Summary
    print(f"\n[prep] Final active batch panel:")
    print(f"  Features: {len(panel):,}")
    print(f"  Genes:    {panel['gene'].nunique():,}")
    if 'site_type' in panel.columns:
        print(f"  site_type:")
        for st, cnt in panel['site_type'].value_counts().items():
            print(f"    {st:<15} {cnt:>8,}")

    panel.to_csv(args.out, index=False)
    print(f"\n[SUCCESS] Active batch panel saved → {args.out}")


if __name__ == "__main__":
    main()
