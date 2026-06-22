#!/usr/bin/env python3
"""
make_gene_list.py
-----------------
Extract expressed gene names from a 10x Genomics filtered feature
barcode matrix (Cell Ranger or Space Ranger output).

Compatible with:
  - Cell Ranger (scRNA-seq): filtered_feature_bc_matrix/
  - Space Ranger (Visium):   filtered_feature_bc_matrix/
  - Any 10x pipeline with matrix.mtx.gz + features.tsv.gz + barcodes.tsv.gz

Usage:
  # scRNA-seq (Cell Ranger)
  python make_gene_list.py \
      --matrix /path/to/filtered_feature_bc_matrix/ \
      --out    data/gene_list.txt \
      --min-cells 3

  # Visium (Space Ranger)
  python make_gene_list.py \
      --matrix /path/to/filtered_feature_bc_matrix/ \
      --out    data/gene_list_visium.txt \
      --min-cells 1

Copyright (c) 2026 IsoMatrix Suite. All rights reserved.
"""

import argparse
import gzip
import os
import scipy.io
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract expressed gene names from 10x filtered matrix. "
                    "Compatible with Cell Ranger (scRNA-seq) and "
                    "Space Ranger (Visium)."
    )
    parser.add_argument("--matrix", required=True,
                        help="Path to filtered_feature_bc_matrix/ directory")
    parser.add_argument("--out",    required=True,
                        help="Output gene list text file (one gene per line)")
    parser.add_argument("--min-cells", type=int, default=3,
                        help="Min cells/spots expressing a gene (default: 3). "
                             "For scRNA-seq: number of cells. "
                             "For Visium: number of spots.")
    return parser.parse_args()


def load_features(matrix_dir):
    """
    Load features.tsv.gz from 10x matrix directory.
    Returns list of (gene_id, gene_name, feature_type).
    """
    features_path = os.path.join(matrix_dir, "features.tsv.gz")
    if not os.path.exists(features_path):
        # Try without .gz
        features_path = os.path.join(matrix_dir, "genes.tsv")
        if not os.path.exists(features_path):
            raise FileNotFoundError(
                f"features.tsv.gz or genes.tsv not found in {matrix_dir}")
        # Old format: gene_id, gene_name only
        features = []
        with open(features_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                gene_id   = parts[0] if len(parts) > 0 else ''
                gene_name = parts[1] if len(parts) > 1 else parts[0]
                features.append((gene_id, gene_name, 'Gene Expression'))
        return features

    features = []
    with gzip.open(features_path, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            gene_id      = parts[0] if len(parts) > 0 else ''
            gene_name    = parts[1] if len(parts) > 1 else parts[0]
            feature_type = parts[2] if len(parts) > 2 else 'Gene Expression'
            features.append((gene_id, gene_name, feature_type))
    return features


def load_matrix(matrix_dir):
    """Load matrix.mtx.gz as sparse CSR matrix."""
    matrix_path = os.path.join(matrix_dir, "matrix.mtx.gz")
    if not os.path.exists(matrix_path):
        matrix_path = os.path.join(matrix_dir, "matrix.mtx")
    mat = scipy.io.mmread(matrix_path)
    return mat.tocsr()


def main():
    args = parse_args()

    print(f"[make_gene_list] Loading matrix from: {args.matrix}")
    features = load_features(args.matrix)
    mat      = load_matrix(args.matrix)

    n_barcodes = mat.shape[1]
    data_type  = "spots" if n_barcodes < 10_000 else "cells"

    print(f"[make_gene_list] Matrix: {mat.shape[0]:,} features × "
          f"{n_barcodes:,} {data_type}")
    print(f"[make_gene_list] Total features: {len(features):,}")

    # Filter to Gene Expression features only (not Antibody Capture etc.)
    gene_indices = [
        i for i, (_, _, ftype) in enumerate(features)
        if ftype == 'Gene Expression'
    ]
    print(f"[make_gene_list] Gene Expression features: {len(gene_indices):,}")

    # Count cells/spots per gene
    gene_mat          = mat[gene_indices, :]
    cells_per_gene    = np.array((gene_mat > 0).sum(axis=1)).flatten()

    # Filter by min_cells
    expressed_mask  = cells_per_gene >= args.min_cells
    expressed_genes = [
        features[gene_indices[i]][1]
        for i in range(len(gene_indices))
        if expressed_mask[i]
    ]

    # Remove duplicates, sort
    expressed_genes = sorted(set(expressed_genes))

    print(f"[make_gene_list] Expressed genes "
          f"(≥{args.min_cells} {data_type}): {len(expressed_genes):,}")

    # Save
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, 'w') as f:
        f.write('\n'.join(expressed_genes) + '\n')

    print(f"[make_gene_list] Saved → {args.out}")
    print(f"\n  Next step: build_panel_features.py")
    print(f"  python IsoDecipher/scripts/build_panel_features.py \\")
    print(f"      --gtf   /Volumes/Lexar/reference/hg38/Homo_sapiens.GRCh38.115.gtf \\")
    print(f"      --genes {args.out} \\")
    print(f"      --out   results/panel/panel_features.csv \\")
    print(f"      --tolerance 10")


if __name__ == '__main__':
    main()
