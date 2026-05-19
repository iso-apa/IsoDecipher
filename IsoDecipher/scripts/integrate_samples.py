#!/usr/bin/env python3
"""
IsoDecipher: Integrate GEX, Isoform, and CITE-seq data across samples
----------------------------------------------------------------------
Loads Cell Ranger H5 files from GCS, merges with IsoDecipher isoform
counts, and produces a combined AnnData object.

Usage:
python IsoDecipher/scripts/integrate_samples.py \
    --exp_list data/samples.txt \
    --gcs_dir gs://isodecipher-bam/samples/samples \
    --iso_dir results \
    -suffix _isoform_count_expanded.csv \
    --out results/master_mosaic_combined.h5ad
"""
import scanpy as sc
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import argparse
import subprocess
import gc
import os


def download_from_gcs(gcs_path, local_path):
    """Download a file from GCS to local path."""
    result = subprocess.run(
        ["gsutil", "-o", "GSUtil:parallel_process_count=1", "cp", gcs_path, local_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  [WARN] Failed to download {gcs_path}: {result.stderr.strip()}")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Integrate GEX, Isoforms, and CITE-seq across samples"
    )
    parser.add_argument("--exp_list", default=None,
                        help="Text file with one experiment ID per line")
    parser.add_argument("--samples", nargs="+", default=None,
                        help="List of sample IDs (alternative to --exp_list)")
    parser.add_argument("--gcs_dir",
                        default="gs://isodecipher-bam/samples/samples",
                        help="GCS directory containing experiment folders")
    parser.add_argument("--data_dir", default=None,
                        help="Local directory containing experiment folders (overrides GCS)")
    parser.add_argument("--iso_dir", default="results",
                        help="Directory containing isoform CSVs (results/counts/)")
    parser.add_argument('--suffix', default='_isoform_count.csv',
                    help='Isoform count file suffix')
    parser.add_argument("--out", default="results/master_mosaic_combined.h5ad",
                        help="Path for output H5AD")
    args = parser.parse_args()

    # Load sample list from either --samples or --exp_list
    if args.samples:
        exps = args.samples
    elif args.exp_list:
        with open(args.exp_list, 'r') as f:
            exps = [line.strip() for line in f if line.strip()
                    and not line.strip().startswith('#')]
    else:
        raise ValueError("Must provide either --samples or --exp_list")

    print(f"[integrate] Processing {len(exps)} samples: {exps}")
    adatas = []

    for exp in exps:
        print(f"\n--- Processing {exp} ---")

        # --- Load H5 (from GCS or local) ---
        if args.data_dir:
            h5_path = os.path.join(args.data_dir, exp, "filtered_feature_bc_matrix.h5")
            if not os.path.exists(h5_path):
                print(f"  [SKIP] H5 not found: {h5_path}")
                continue
        else:
            # Download from GCS to /tmp
            gcs_h5 = f"{args.gcs_dir}/{exp}/filtered_feature_bc_matrix.h5"
            h5_path = f"/tmp/{exp}_features.h5"
            print(f"  Downloading H5 from GCS...")
            if not download_from_gcs(gcs_h5, h5_path):
                continue

        adata_h5 = sc.read_10x_h5(h5_path, gex_only=False)
        adata_h5.var_names_make_unique()

        # Clean up /tmp after loading
        if not args.data_dir and os.path.exists(h5_path):
            os.remove(h5_path)

        # Split GEX and ADT
        genes = adata_h5[:, adata_h5.var['feature_types'] == 'Gene Expression'].copy()
        proteins = adata_h5[:, adata_h5.var['feature_types'] == 'Antibody Capture'].copy()
        del adata_h5
        gc.collect()

        print(f"  GEX: {genes.n_obs} cells × {genes.n_vars} genes")

        # --- Load Isoform Counts ---
        # Fixed path: results/counts/{exp}_isoform_count.csv
        iso_path = os.path.join(args.iso_dir, "counts", f"{exp}{args.suffix}")

        if os.path.exists(iso_path):
            df_iso = pd.read_csv(iso_path, index_col=0, engine='c')
            df_iso = df_iso.astype('float32')
            print(f"  Isoform: {df_iso.shape[0]} cells × {df_iso.shape[1]} features")

            adata_iso = ad.AnnData(X=csr_matrix(df_iso.values))
            adata_iso.obs_names = df_iso.index
            adata_iso.var_names = df_iso.columns
            adata_iso.var['feature_types'] = 'Isoform'

            # Normalize barcodes (strip -1 suffix if present)
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

        # --- Add ADT if exists ---
        if proteins.n_vars > 0:
            print(f"  ADT: {proteins.n_vars} proteins")
            proteins.obs_names = proteins.obs_names.str.replace("-1", "", regex=False)
            common_pro = combined_sample.obs_names.intersection(proteins.obs_names)
            proteins = proteins[common_pro].copy()
            proteins.var_names = [f"prot_{n}" for n in proteins.var_names]
            proteins.var['feature_types'] = 'ADT'

            # Align cells
            combined_sample = combined_sample[common_pro].copy()
            combined_sample = ad.concat(
                [combined_sample, proteins], axis=1, merge="first"
            )

        # --- Finalize sample metadata ---
        combined_sample.obs['batch'] = exp
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

    # Re-apply feature_types after concat
    adata_final.var['feature_types'] = 'Gene Expression'
    adata_final.var.loc[
        adata_final.var_names.str.contains(r'_G\d+', regex=True),
        'feature_types'
    ] = 'Isoform'
    adata_final.var.loc[
        adata_final.var_names.str.startswith('prot_'),
        'feature_types'
    ] = 'ADT'

    # Summary
    n_gex = (adata_final.var['feature_types'] == 'Gene Expression').sum()
    n_iso = (adata_final.var['feature_types'] == 'Isoform').sum()
    n_adt = (adata_final.var['feature_types'] == 'ADT').sum()

    print(f"\n[SUCCESS] Combined AnnData:")
    print(f"  Cells:          {adata_final.n_obs:,}")
    print(f"  GEX features:   {n_gex:,}")
    print(f"  Isoform features: {n_iso:,}")
    print(f"  ADT features:   {n_adt:,}")

    # Save
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    adata_final.write(args.out)
    print(f"\n✅ Saved to: {args.out}")


if __name__ == "__main__":
    main()
