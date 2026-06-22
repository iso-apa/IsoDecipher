"""
Entropy Velocity Calculation Script
Calculates the structural remodeling rate of isoform entropy 
along a pseudotemporal trajectory.
"""

import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline
import scanpy as ad

def calculate_entropy_velocity(adata, df_entropy, pseudotime_key='iso_pseudotime', s_factor=0.01):
    """
    Computes d(Entropy)/d(Pseudotime) using smoothing splines.
    """
    print(f"⚡ Starting alignment and velocity calculation using {pseudotime_key}...")

    # 1. Align cells between entropy matrix and AnnData object
    common_cells = df_entropy.index.intersection(adata.obs_names)
    dpt_full = adata.obs.loc[common_cells, pseudotime_key].values
    df_entropy_aligned = df_entropy.loc[common_cells]

    df_velocity = pd.DataFrame(index=common_cells, columns=df_entropy_aligned.columns)

    # 2. Iterative calculation per gene
    for gene in df_entropy_aligned.columns:
        entropy_vals = df_entropy_aligned[gene].values
        
        # Filter NaN values for stable spline fitting
        mask = ~np.isnan(entropy_vals) & ~np.isnan(dpt_full)
        dpt_clean = dpt_full[mask]
        ent_clean = entropy_vals[mask]
        
        # Minimum requirement for spline interpolation
        if len(dpt_clean) < 10:
            df_velocity.loc[common_cells[mask], gene] = 0
            continue
            
        # Sort by pseudotime
        idx = np.argsort(dpt_clean)
        dpt_sorted = dpt_clean[idx]
        ent_sorted = ent_clean[idx]
        
        # Fit smoothing spline
        # s_factor controls smoothness; higher = more rigid, lower = more sensitive to noise
        spline = UnivariateSpline(dpt_sorted, ent_sorted, s=len(ent_sorted) * s_factor)
        
        # Calculate first derivative (velocity)
        velocity = spline.derivative()(dpt_clean)
        df_velocity.loc[common_cells[mask], gene] = velocity

    # 3. Store results back into AnnData
    # Filling NaN with 0 for cells outside the high-confidence window
    adata.obsm['isoform_velocity'] = df_velocity.reindex(adata.obs_names).fillna(0).values
    
    print(f"✅ Calculation complete! Velocity matrix shape: {adata.obsm['isoform_velocity'].shape}")
    return adata

if __name__ == "__main__":
    # Integration point: Replace these with your actual objects
    # adata_gex = ...
    # df_gene_entropy = ...
    
    # Example call:
    # adata_gex = calculate_entropy_velocity(adata_gex, df_gene_entropy)
    pass