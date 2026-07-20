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
import warnings
import scanpy as sc
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import scipy.sparse
import argparse
import gc
import os

def main():
    parser = argparse.ArgumentParser(description="Integrate GEX, Isoforms, and CITE-seq across samples")
    parser.add_argument("--samples", nargs="+", required=True, help="List of sample IDs")
    parser.add_argument("--data_dir", required=True, help="Local directory containing experiment folders")
    parser.add_argument("--iso_dir", required=True, help="Directory containing isoform CSVs")
    parser.add_argument('--suffix', default='{sample}_isoform_count',
                        help='Isoform base path pattern (no extension). {sample} is replaced '
                             'with the sample ID, e.g. {sample}/counts/{sample}_isoform_count')
    parser.add_argument("--out", required=True, help="Path for output H5AD")
    parser.add_argument("--panel", default=None,
                        help="Path to panel_features.csv for iso.var metadata enrichment")
    parser.add_argument("--h5-paths", nargs="+", default=None,
                        help="Per-sample h5 overrides in sample=path form. "
                             "Supersedes data_dir-based discovery.")
    args = parser.parse_args()

    exps = args.samples
    print(f"[integrate] Processing {len(exps)} samples: {exps}")

    h5_map = {}
    if args.h5_paths:
        for item in args.h5_paths:
            k, v = item.split("=", 1)
            h5_map[k] = v

    adatas = []

    for exp in exps:
        print(f"\n--- Processing {exp} ---")

        # 1.
        h5_path = h5_map.get(exp) or os.path.join(args.data_dir, exp, "filtered_feature_bc_matrix.h5")
        if not os.path.exists(h5_path):
            print(f"  [SKIP] H5 not found: {h5_path}")
            continue

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Variable names are not unique")
            adata_h5 = sc.read_10x_h5(h5_path, gex_only=False)
        adata_h5.var_names_make_unique()

        # Split GEX and ADT
        genes = adata_h5[:, adata_h5.var['feature_types'] == 'Gene Expression'].copy()
        proteins = adata_h5[:, adata_h5.var['feature_types'] == 'Antibody Capture'].copy()
        del adata_h5
        gc.collect()
        print(f"  GEX: {genes.n_obs} cells × {genes.n_vars} genes")

        # 2.
        iso_base = os.path.join(args.iso_dir, args.suffix.replace('{sample}', exp))
        npz_path = f"{iso_base}_matrix.npz"
        obs_path = f"{iso_base}_obs.txt"
        var_path = f"{iso_base}_var.txt"

        if os.path.exists(npz_path):
            mat           = scipy.sparse.load_npz(npz_path)
            obs_names_iso = open(obs_path).read().splitlines()
            var_names_iso = open(var_path).read().splitlines()
            print(f"  Isoform: {mat.shape[0]} cells × {mat.shape[1]} features")

            adata_iso = ad.AnnData(X=mat)
            adata_iso.obs_names = obs_names_iso
            adata_iso.var_names = var_names_iso
            adata_iso.var['feature_types'] = 'Isoform'
            del mat
            gc.collect()

            # Normalize barcodes
            genes.obs_names     = genes.obs_names.str.replace("-1", "", regex=False)
            adata_iso.obs_names = pd.Index(obs_names_iso).str.replace("-1", "", regex=False)

            common_cells = genes.obs_names.intersection(adata_iso.obs_names)
            print(f"  Common cells (GEX ∩ Isoform): {len(common_cells)}")

            genes     = genes[common_cells].copy()
            adata_iso = adata_iso[common_cells].copy()

            combined_sample = ad.concat([genes, adata_iso], axis=1, merge="first")
            del genes, adata_iso
            gc.collect()
        else:
            print(f"  [WARN] Isoform NPZ not found: {npz_path} — using GEX only")
            combined_sample = genes
            del genes
            gc.collect()

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

    # --- Join panel metadata into iso.var ---
    if args.panel and os.path.exists(args.panel):
        print(f"\n[integrate] Joining panel metadata from {args.panel}...")
        panel = pd.read_csv(args.panel)
        panel['feature_id'] = panel['gene'] + '_G' + panel['polyA_group'].astype(str)
        # adata.var columns for isoform features (panel v1.1)
        VAR_COLS = [c for c in [
            'gene',               # gene symbol
            'rep_coord',          # representative cleavage coordinate
            'avg_utr_bp',         # spliced UTR length (bp); NaN for IPA/ALE
            'utr_confidence',     # spliced | estimated | None
            'pc_fraction',        # fraction of coding transcripts (0-1)
            'panel_source',       # GTF | PA2 | GTF+PA2
            'v2_fraction',        # PolyASite v2.0 sample support (0-1)
            'v2_site_class',      # TE | IN | EX | DS
            'site_type',          # Tandem_UTR | IPA | ALE | 3'Ext | Unknown
            'utr_has_intron',     # True if 3'UTR contains intron (NMD flag)
            'nmd_fraction',       # fraction of NMD transcripts at this site
            'ig_fraction',        # fraction of IG gene transcripts
            'tr_fraction',        # fraction of TR gene transcripts
            'ig_pseudo_fraction', # fraction of IG pseudogene transcripts
            'tr_pseudo_fraction', # fraction of TR pseudogene transcripts
            'pas_motif_strength', # strong | moderate | weak | unknown_gtf | unknown_pa2
            'cds_end_coord',      # stop codon genomic coordinate (CDS anchor)
            'cds_group_id',       # {gene}_CDS{coord} or {gene}_noCDS
        ] if c in panel.columns]

        panel_indexed = panel.set_index('feature_id')[VAR_COLS].copy()
        panel_indexed = panel_indexed.rename(columns={'gene': 'gene_name'})
        panel_indexed = panel_indexed[~panel_indexed.index.duplicated(keep='first')]
        adata_final.var = adata_final.var.join(panel_indexed, how='left')
        n_joined = adata_final.var['panel_source'].notna().sum() if 'panel_source' in adata_final.var.columns else adata_final.var['avg_utr_bp'].notna().sum()
        print(f"  Joined metadata for {n_joined:,} isoform features")
    else:
        if args.panel:
            print(f"[WARN] Panel file not found: {args.panel} — skipping metadata join")

    n_gex = (adata_final.var['feature_types'] == 'Gene Expression').sum()
    n_iso = (adata_final.var['feature_types'] == 'Isoform').sum()
    n_adt = (adata_final.var['feature_types'] == 'ADT').sum()

    print(f"\n[SUCCESS] Combined AnnData:")
    print(f"  Cells:          {adata_final.n_obs:,}")
    print(f"  GEX features:   {n_gex:,}")
    print(f"  Isoform features: {n_iso:,}")
    print(f"  ADT features:   {n_adt:,}")

    # Fix var dtypes before writing — panel join leaves NaN in bool/string columns
    # for GEX/ADT features, causing h5py TypeError on write
    if 'utr_has_intron' in adata_final.var.columns:
        adata_final.var['utr_has_intron'] = (
            adata_final.var['utr_has_intron'] == True
        )
    for col in ['utr_confidence', 'panel_source', 'v2_site_class', 'site_type', 'gene_name']:
        if col in adata_final.var.columns:
            adata_final.var[col] = adata_final.var[col].fillna('').astype(str)

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    adata_final.write(args.out)
    print(f"\n✅ Saved to: {args.out}")

if __name__ == "__main__":
    main()