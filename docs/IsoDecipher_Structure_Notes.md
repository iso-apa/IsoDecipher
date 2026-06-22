# IsoDecipher — Architecture & Figure Planning Notes

> Last updated: June 2026
> Author: Rene Cheng

---

## Paper Figure Structure

### Figure 1 (drawn LAST): IsoDecipher Platform Overview
> To be designed after Figure 2-5 are complete.
> Functions as graphical abstract / system overview.

- **Panel A**: The Biological Blindspot
  - GEX: sees only abundance change
  - IsoDecipher: sees structural rewiring (UTR length, isoform switching)

- **Panel B**: Core Mathematical Engine
  - Isoform Entropy: per-gene Shannon entropy of isoform usage
  - Global Decision Entropy: cell-level mean across all valid genes
  - C_state formula: [UTR matrix + Isoform Usage + Entropy] → true cell state

- **Panel C**: Computational Workflow
  - BAM → Panel Build → Insert Size → Soft Assignment → Features → C_state → Velocity

- **Panel D**: Performance Benchmark (numbers from Figure 2-5)
  - vs GEX-based trajectory
  - vs scAPA / Sierra / MAAPER
  - Metrics TBD after biological figures complete

---

### Figure 2: IsoDecipher Benchmark
> First biological figure. Prove isoform-based trajectory works independently.

- Dataset: GSE212138 (exp97) + GSE229042 (exp105, exp106)
  - n=~18,500 cells after QC + IFN filter (ifn_score>0 & bcell_score>1.3)
- GEX HVG UMAP vs Isoform Entropy HVG UMAP (side by side)
- DPT pseudotime on both
- Canonical gene trajectories (GEX vs Isoform side by side):
  - Non-IG: MS4A1, PAX5, PRDM1, XBP1, CD59, ITGB2, TNFRSF13C
  - IG (isotype-masked): IGHM, IGHG1, IGHG3, IGHA1 (G0 Secreted vs G1 Membrane)
- GEX vs Isoform pseudotime correlation: Spearman R=0.656

---

### Figure 3: APA Biology — UTR Shortening & Stressed PC
> The biological discovery figure.

- **Panel A**: Global UTR shortening
  - Weighted mean UTR along B cell score
  - Waddington landscape (Before/Switching/Terminal)

- **Panel B**: Global APA Entropy
  - Inverted-U pattern along pseudotime
  - Peak at switching window
  - B cell score vs entropy scatter (hexbin)

- **Panel C**: Stressed PC finding
  - High APA entropy PC = UPR/apoptosis enrichment
  - Scatter: B cell score vs entropy, colored by UPR stress score
  - Module scores: UPR_Stress, Secretory_Stress, RNA_Splicing, Apoptosis
  - SHAP analysis: MALAT1, SRSF5, LINC-PINT, XAF1, SSR3, EIF1

- **Panel D**: Gene-level exceptions
  - PRDM1 paradox: long UTR in terminal PC (counter to global shortening)
  - SRSF5 isoform switching
  - Other interesting genes (TMBIM6, CD59, DERL3, DNAJB6)

---

### Figure 4: C_state Manifold
> Methodological leap — multi-modal cell state definition.

- **Panel A**: C_state formula / schematic
  - [UTR Matrix + Isoform Usage Matrix + Global Decision Entropy]
  - Normalization strategy (Z-score per modality)
  - Integration: weighted PCA or VAE latent space

- **Panel B**: C_state reconstructed UMAP
  - Normal PC vs Stressed PC cleanly separated
  - Compare with Figure 2 UMAP

- **Panel C**: Benchmark
  - GEX UMAP vs Isoform UMAP vs C_state UMAP
  - Metric: separation of Normal PC vs Stressed PC
  - TBD: Silhouette score / cluster purity

> Status: Design phase. Implementation TBD.

---

### Figure 5: Structural Dynamics & Velocity
> The culmination — IsoDecipher predicts cell fate direction.

- **Panel A**: KNN entropy velocity vector field on C_state UMAP
- **Panel B**: Manifold displacement (structural mutation hotspots)
- **Panel C**: 1D spline velocity for key genes (PRDM1 isoforms etc.)

> Status: Design phase. Implementation TBD.

---

## IsoDecipher Pipeline Architecture

### 1. Panel Build (`build_panel_features.py`)
- Input: GTF + gene list
- Consecutive gap-based clustering, tolerance=125bp (IGV validated)
- Group label: `{gene}_G{polyA_group}_{label}`
- Output: `panel_features.csv`

### 2. Insert Size Validation (`validate_insert_size.py`)
- Per-sample KDE lookup (no parametric assumption)
- Results:
  - exp97: median=158bp, P5=28bp, P95=401bp
  - exp105: median=157bp, P5=28bp, P95=400bp
  - exp106: median=169bp, P5=31bp, P95=409bp
- Output: `{sample}_insert_size_dist_lookup.json`

### 3. Soft Assignment (`assign_reads.py`)
- Fractional counting: `soft_count[(cb, feature)] += score / total_score`
- Sum per UMI = 1.0
- Output: UMI-deduplicated isoform matrix (float, 4 decimal places)

---

## Features

### Per-Cell
| Feature | Description |
|---|---|
| `weighted_utr` | Expression-weighted mean 3' UTR (bp) |
| `Global_APA_Entropy` | Mean per-gene Shannon entropy (bits), non-IGH genes only |
| `bcell_score` | naive_score − plasma_score |
| `iso_pseudotime` | DPT on isoform entropy HVG manifold |
| `gex_pseudotime` | DPT on GEX manifold |

### Per-Gene
| Feature | Description |
|---|---|
| Per-gene APA entropy | Shannon entropy of isoform usage per gene per cell |
| Per-gene weighted UTR | Expression-weighted UTR for single gene |
| Proximal site score | Σ(expr_i × exp(−i)) / Σ(expr_i) |

### Planned (C_state components)
| Feature | Description | Status |
|---|---|---|
| UTR matrix | Per-gene weighted UTR per cell | Implemented |
| Isoform usage matrix | Per-gene proximal site score per cell | Implemented |
| Global Decision Entropy | Mean per-gene entropy | Implemented |
| C_state latent space | Multi-modal integration | Design |

---

## Analytical Tools

### Implemented
| Tool | Key Result |
|---|---|
| Entropy trajectory | Inverted-U along pseudotime |
| UTR trajectory | 450bp → 150bp shortening |
| Waddington landscape | 3 distinct valleys |
| B cell score | Naive/Plasma cell identity axis |
| RF + SHAP | Library size rank#1; MALAT1 rank#2 after correction |
| Stress module scores | UPR/Apoptosis enriched in high entropy PC |

### Planned
| Tool | Description | Status |
|---|---|---|
| KNN entropy velocity | APA entropy direction in KNN graph | Design |
| DPT-based entropy velocity | 1D spline derivative | Design |
| Manifold displacement | APA-based cell displacement | Design |
| C_state manifold | Multi-modal UMAP | Design |
| IsoDecipher-GPT | Foundation model for APA profiles | Active dev |

---

## Key Decisions & Rationale

| Decision | Rationale |
|---|---|
| Tolerance=125bp consecutive gap | IGV validated: CD59 G0/G1/G2 = same peak |
| Soft assignment (not hard) | Bhattacharyya BC=0.94 for 35bp apart sites |
| IGH excluded from entropy HVG | IGH dominates entropy variance; independent benchmark |
| IGH excluded from isoform UMAP | Compositional bias; IG is output not driver of differentiation |
| IFN filter (ifn_score>0 & bcell_score>1.3) | exp97 activated B cells contaminate trajectory |
| Library size as RF covariate | R=0.429 between entropy and library size; biological coupling confirmed |
| Entropy computed on raw counts | Proportion-based; scale-invariant; decoupled from GEX normalization |
| .raw = post-normalize pre-scale | Scanpy convention; true raw counts in layers['raw_counts'] and adata_true_raw |

---

## Sensitivity Analysis Results

| Test | R | Conclusion |
|---|---|---|
| Global entropy vs library size | 0.429 | Biological coupling, not pure artifact |
| Stable gene entropy (≥50% cells) vs library size | 0.453 | Consistent with original |
| Original vs stable entropy | 0.792 | Robust across gene sets |
| MALAT1 RF rank after library size covariate | #2 | Independent biological signal |

---

## Code To Add Later
- [ ] `build_panel_features.py` core loop
- [ ] `assign_reads.py` soft assignment loop
- [ ] Entropy computation
- [ ] Weighted UTR computation
- [ ] RF + SHAP analysis
- [ ] C_state integration code (TBD)

