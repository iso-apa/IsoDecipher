#!/usr/bin/env python3
# ==============================================================================
# Copyright (c) 2026 IsoMatrix Suite. All rights reserved.
# 
# IsoDecipher is dual-licensed:
# 1. For academic and non-commercial use, it is licensed under the AGPLv3.
# 2. For commercial and enterprise use, a Commercial License is required.
# 
# See the LICENSE file in the project root for more details.
# ==============================================================================
"""
IsoDecipher: Integrate GEX, Isoform, and CITE-seq data across samples
"""
import scanpy as sc
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import argparse
import gc
import os

def main():
    parser = argparse.ArgumentParser(description="Integrate GEX, Isoforms, and CITE-seq across samples")
    parser.add_argument("--samples", nargs="+", required=True, help="List of sample IDs")
    parser.add_argument("--data_dir", required=True, help="Local directory containing experiment folders")
    parser.add_argument("--iso_dir", required=True, help="Directory containing isoform CSVs")
    parser.add_argument('--suffix', default='_isoform_count.csv', help='Isoform count file suffix')
    parser.add_argument("--out", required=True, help="Path for output H5AD")
    args = parser.parse_args()

    exps = args.samples
    print(f"[integrate] Processing {len(exps)} samples: {exps}")
    adatas = []

    for exp in exps:
        print(f"\n--- Processing {exp} ---")

        # 1. 
        h5_path = os.path.join(args.data_dir, exp, "filtered_feature_bc_matrix.h5")
        if not os.path.exists(h5_path):
            print(f"  [SKIP] H5 not found: {h5_path}")
            continue

        adata_h5 = sc.read_10x_h5(h5_path, gex_only=False)
        adata_h5.var_names_make_unique()

        # Split GEX and ADT
        genes = adata_h5[:, adata_h5.var['feature_types'] == 'Gene Expression'].copy()
        proteins = adata_h5[:, adata_h5.var['feature_types'] == 'Antibody Capture'].copy()
        del adata_h5
        gc.collect()
        print(f"  GEX: {genes.n_obs} cells × {genes.n_vars} genes")

        # 2. 
        iso_path = os.path.join(args.iso_dir, f"{exp}{args.suffix}")

        if os.path.exists(iso_path):
            df_iso = pd.read_csv(iso_path, index_col=0, engine='c')
            df_iso = df_iso.astype('float32')
            print(f"  Isoform: {df_iso.shape[0]} cells × {df_iso.shape[1]} features")

            adata_iso = ad.AnnData(X=csr_matrix(df_iso.values))
            adata_iso.obs_names = df_iso.index
            adata_iso.var_names = df_iso.columns
            adata_iso.var['feature_types'] = 'Isoform'

            # Normalize barcodes
            genes.obs_names = genes.obs_names.str.replace("-1", "", regex=False)
            adata_iso.obs_names = adata_iso.obs_names.str.replace("-1", "", regex=False)

            common_cells = genes.obs_names.intersection(adata_iso.obs_names)
            print(f"  Common cells (GEX ∩ Isoform): {len(common_cells)}")

            genes = genes[common_cells].copy()
            adata_iso = adata_iso[common_cells].copy()

            # Combine GEX + Isoform
            combined_sample = ad.concat([genes, adata_iso], axis=1, merge="first")
        else:
            print(f"  [WARN] Isoform CSV not found: {iso_path} — using GEX only")
            combined_sample = genes

        # 3. 加入 ADT
        if proteins.n_vars > 0:
            print(f"  ADT: {proteins.n_vars} proteins")
            proteins.obs_names = proteins.obs_names.str.replace("-1", "", regex=False)
            common_pro = combined_sample.obs_names.intersection(proteins.obs_names)
            proteins = proteins[common_pro].copy()
            proteins.var_names = [f"prot_{n}" for n in proteins.var_names]
            proteins.var['feature_types'] = 'ADT'

            combined_sample = combined_sample[common_pro].copy()
            combined_sample = ad.concat([combined_sample, proteins], axis=1, merge="first")

        combined_sample.obs['batch']  = exp
        combined_sample.obs['sample'] = exp
        combined_sample.obs_names = [f"{exp}_{bc}" for bc in combined_sample.obs_names]

        if not isinstance(combined_sample.X, csr_matrix):
            combined_sample.X = csr_matrix(combined_sample.X)

        print(f"  Final: {combined_sample.n_obs} cells × {combined_sample.n_vars} features")
        adatas.append(combined_sample)
        gc.collect()

    if not adatas:
        print("[ERROR] No samples loaded successfully.")
        return

    # --- Merge all samples ---
    print(f"\n[integrate] Merging {len(adatas)} samples...")
    adata_final = ad.concat(adatas, join='outer', fill_value=0)

    adata_final.var['feature_types'] = 'Gene Expression'
    adata_final.var.loc[adata_final.var_names.str.contains(r'_G\d+', regex=True), 'feature_types'] = 'Isoform'
    adata_final.var.loc[adata_final.var_names.str.startswith('prot_'), 'feature_types'] = 'ADT'

    clean_names = [name[:-4] if name.endswith('_nan') else name for name in adata_final.var_names]
    adata_final.var_names = clean_names

    n_gex = (adata_final.var['feature_types'] == 'Gene Expression').sum()
    n_iso = (adata_final.var['feature_types'] == 'Isoform').sum()
    n_adt = (adata_final.var['feature_types'] == 'ADT').sum()

    print(f"\n[SUCCESS] Combined AnnData:")
    print(f"  Cells:          {adata_final.n_obs:,}")
    print(f"  GEX features:   {n_gex:,}")
    print(f"  Isoform features: {n_iso:,}")
    print(f"  ADT features:   {n_adt:,}")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    adata_final.write(args.out)
    print(f"\n✅ Saved to: {args.out}")

if __name__ == "__main__":
    main()