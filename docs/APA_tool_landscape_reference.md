# APA Tool Landscape — Verified Reference Summary
*Compiled from web search verification, 2026-06-19*
*For use in IsoDecipher/IsoCAPE related-work writing*

---

## 1. Classification Framework (from Bi et al. 2026, NAR — independent 3rd-party benchmark)

Bi X, Chen Z, Ye M, Zhang T, He D, Wu X. "Benchmarking computational methods for
identifying and quantifying polyadenylation sites from 3' tag-based single-cell
RNA-seq data." Nucleic Acids Research, Volume 54, Issue 9, 2026, gkag490.
https://doi.org/10.1093/nar/gkag490

This is the most comprehensive, independent (non-tool-author) benchmark found.
10 methods, 9 simulated + 25 real datasets, 4 protocols, 3 species.

### Category A: De novo pA identification (peak-calling based)
Tools: scAPA, polyApipe, Sierra, scAPAtrap, SCAPE

### Category B: Annotation-based pA identification
Subcategory B1 "read-support filtering of annotated pAs": MAAPER, scUTRquant
Subcategory B2 "peak calling followed by annotation filtering": SCAPTURE, Infernape, scraps

**IsoDecipher belongs conceptually in Category B** (annotation/GTF-anchored),
closest in spirit to scUTRquant and MAAPER.

---

## 2. Individual Tool Profiles

### SCAPE / SCAPE-APA
- Zhou et al. 2022, Nucleic Acids Research, e66 (original SCAPE)
- Updated version: SCAPE-APA, bioRxiv 2024.03.12.584547 (Taichi-accelerated, binned reads)
- GitHub: LuChenLab/SCAPE
- **Method**: Bayesian mixture model. Each pA site = one Gaussian component
  (mean mu, std delta, weight pi) + one uniform noise component. Insert size (fragment
  length f) and poly(A) tail length (s) are BOTH modeled with uncertainty;
  pA site position theta = x + (f - s) - 1, inferred via EM + BIC for component count.
- **De novo**: pA site position is an ESTIMATED PARAMETER, not a fixed coordinate.
- Self-reported (in original paper): outperforms Sierra/scAPA/scAPAtrap/SCAPTURE/MAAPER
  on simulated data (precision/recall/F-score).
- **Independent benchmark (Bi et al. 2026) — less favorable**:
  - "regardless of window size, SCAPE and scAPA exhibited relatively poor
    prediction precision"
  - >62% of SCAPE/scAPA unique pAs lacked any poly(A) signal motif; only ~20%
    predicted as true positive by DeepPASTA
  - Performed poorly on CEL-seq/Drop-seq (insert-size params tuned for 10x/Microwell-seq)
  - BUT: retained high sensitivity at low depth (15M reads) — robust to downsampling,
    similar to scUTRquant/MAAPER, unlike Sierra/scraps
- **Key difference from IsoDecipher**: SCAPE uses insert size to INFER unknown
  site location (de novo). IsoDecipher uses insert size (KDE) to assign reads
  PROBABILISTICALLY among already-known candidate sites (GTF-anchored). Same raw
  signal (insert size), fundamentally different use (estimation vs. classification).

### Sierra
- GitHub: VCCRI/Sierra
- De novo, peak-calling based, coverage-fitting approach
- Recommended for high-depth data (per Bi et al. 2026 conclusions)
- Highly depth-sensitive: "Sierra's coverage-fitting approach... led to sharp
  pA yield drops at low coverage" (Bi et al. 2026)

### scAPAtrap
- Identification + quantification from scRNA-seq
- Can detect intergenic pAs (broader genomic region coverage than most tools)
- Recommended for "moderate read lengths including intergenic sites" (benchmark consensus)
- Highly sensitive (detects more sites than Sierra in some comparisons) but
  lower precision/localization accuracy in independent benchmark

### scDaPars (+ DaPars / DaPars2, bulk precursor)
- DaPars: Xia et al. 2014, Nature Communications (original bulk tool)
- scDaPars: Gao YP et al., Genome Research 2021. GitHub: YiPeng-Gao/scDaPars
- **Method**: Computes PDUI (Percentage of Distal pA site Usage Index) =
  distal isoform abundance / (distal + proximal isoform abundance). ONLY
  handles two-endpoint (proximal vs. distal) APA, not multi-site.
- Step 1 (DaPars2): linear regression to find proximal pA location from coverage.
  Step 2 (scDaPars): KNN graph on sparse PDUI matrix + NNLS regression to
  IMPUTE dropout genes using neighboring cells.
- **No evidence found of biotype (protein_coding vs non_coding) filtering**
  in any reviewed scDaPars/DaPars documentation or methods description.
- **Key structural difference from IsoDecipher**: PDUI is a strictly two-point
  ratio (proximal/distal only); cannot represent genes with >2 isoform groups
  (e.g., the ST6GAL1 9-group, DDX5 7-group cases IsoDecipher panel handles).
  IsoDecipher's weighted UTR / entropy framework is N-point by design.
- Uses imputation (NNLS borrowing from neighbors) — different philosophy from
  IsoDecipher's no-imputation, soft-probabilistic-assignment approach.

### scUTRquant (Mayrlab)
- Bioinformatics pipeline, Bioconductor SingleCellExperiment output, AnnData support (v0.5.0+)
- GitHub: Mayrlab/scUTRquant
- Manuscript: "Quantification of alternative 3'UTR isoforms from single cell
  RNA-seq data with scUTRquant" (bioRxiv 2021.11.22.469635); follow-up
  "Comprehensive annotation of 3'UTRs from primary cells..."
- **Method**: Built on kallisto bus (pseudo-alignment). Augments GENCODE
  protein-coding transcripts with high-confidence cleavage sites called from
  Mouse Cell Atlas / Human Cell Landscape aggregate data. Truncates transcripts
  to last 500nt, deduplicates. Uses cleanUpdTSeq (Bioconductor) to filter
  internal-priming-derived false cleavage sites.
- **Explicitly states**: "augmenting the PROTEIN CODING transcripts that have
  verified 3' ends" — i.e., scUTRquant by design ONLY uses protein-coding
  transcripts. This is the closest published precedent to IsoDecipher's
  utr_source filtering, BUT:
  - No evidence found that they specifically address IG_C_gene/TR_C_gene
    misclassification (the immunoglobulin biotype edge case IsoDecipher fixed)
  - No evidence found of a quantified demonstration of "what happens if you
    DON'T filter" (i.e., no equivalent of IsoDecipher's R=-0.624 noncoding-
    fraction-vs-differentiation finding)
- Uses txcutr (Bioconductor pkg) for generating compatible indexes from
  arbitrary GTF — flexible annotation-anchored design, philosophically closest
  to IsoDecipher's panel construction approach.
- Robust to depth: retained high sensitivity at 15M-read downsampling (Bi et al. 2026)

### MAAPER
- Li WV, Zheng D, Wang R et al. 2021, Genome Biology, 22:222
- Annotation-based, "read-support filtering of annotated pAs"
- Independent benchmark: "MAAPER's strict adherence to annotations maximized
  precision but limited novel pA discovery" — highest overlap rate with
  reference pAs among ALL methods tested; near-100% precision at 30nt window
- Outputs only POPULATION-level (not single-cell resolution) pA expression —
  important limitation vs. IsoDecipher's true single-cell soft-assignment output
- Robust to depth (retains high sensitivity at low coverage, like scUTRquant/SCAPE)

### SCAPTURE
- Genome Biology 2021. "a deep learning-embedded pipeline"
- De novo peak calling + DeepPASS (embedded deep learning model) for PAS validation
- **Directly relevant observation**: explicitly notes "a bias of mapped reads
  to annotated... PASs in 10x Chromium data for BOTH protein-coding (e.g. GAPDH)
  AND noncoding (e.g. NORAD, GAS5) genes" — i.e., SCAPTURE authors OBSERVED
  that non-coding genes also show 3'-end read pileup in 10x data.
  **However**: this was reported as a descriptive technical observation about
  PAS-calling behavior, NOT followed up with any quantification of how this
  affects downstream summary statistics (e.g., weighted UTR, PDUI) or any
  correlation with differentiation/cell state. This is the key gap IsoDecipher's
  finding (R=-0.624, noncoding fraction vs. bcell_score) fills.
- Robust to depth ("relatively smaller declines" in pA yield at low coverage)

### Infernape, scraps, polyApipe, scAPA
- Lower priority / less directly comparable; scraps and MAAPER had near-100%
  precision at 30nt window but identified far fewer unique/novel pAs (most
  conservative methods in the benchmark)
- polyApipe: de novo but achieved comparable annotation-overlap rate to
  annotation-based methods; among the most reliable de novo methods tested

---

## 3. Bulk RNA-seq precursors (for IsoCAPE / general context)
- DaPars (2014, Nat Commun) — first de novo bulk APA tool, proximal/distal comparison
- DaPars2 — updated version, used as upstream step for scDaPars
- TAPAS, Aptardi — other bulk de novo tools
- APAtrap — bulk precursor to scAPAtrap

---

## 4. Key Cross-Cutting Findings Relevant to IsoDecipher's Methodology

### 4.1 Depth-dependence is a known, benchmarked, universal problem
Bi et al. 2026 directly tested depth-dependence via downsampling a 146M-read
mouse sperm BAM to 10-90%:
- **Depth-robust** (retain sensitivity at low depth): scUTRquant, SCAPE, MAAPER
  — all are annotation-anchored / pre-defined-site methods
- **Depth-sensitive** (sharp decline at low depth): Sierra, scraps
- This independently validates IsoDecipher's design choice (GTF-anchored,
  fixed coordinates) as conferring depth robustness — but does NOT mean
  IsoDecipher's downstream summary statistics (entropy, weighted UTR) are
  automatically immune to depth confounds; that requires separate validation
  (which IsoDecipher's own analysis already demonstrated: entropy IS depth-
  confounded at the per-cell scalar level even with GTF-anchored sites,
  while weighted UTR is depth-robust).

### 4.2 Aggregation reduces quantification error (precedent for "consensus" strategies)
From simulated-data benchmarking: "focusing on multi-method consensus sites
reduced group-level MAPE to near barcode levels... elevated group-level error
stems from low-abundance, low-stability true pAs identified by few methods."
This is independent validation of the general principle behind pseudobulk/
microbulk aggregation and consensus-feature-selection strategies for stabilizing
low-count APA quantification — directly relevant to the entropy correction
strategies explored for IsoDecipher (Chao-Shen, consensus filtering, KNN-based
microbulk for Figure 5 velocity).

### 4.3 Annotation-based vs. de novo: precision/discovery trade-off
"annotation-dependent methods (scraps, MAAPER, scUTRquant) exhibit superior pA
localization accuracy than de novo methods" but at the cost of limited novel
site discovery. This is the central trade-off IsoDecipher should explicitly
acknowledge: choosing GTF-anchored panel construction sacrifices de novo
discovery capability in exchange for precision, depth-robustness, and
interpretability tied to known transcript annotation.

---

## 5. Gaps / Original Contribution Candidates for IsoDecipher (verified, conservative framing)

1. **IG_C_gene/TR_C_gene biotype misclassification fix**: No evidence found in
   any reviewed tool's documentation of explicit handling for immunoglobulin/
   TCR gene biotypes being mistakenly excluded from "protein_coding" classification.
   This appears to be a genuinely under-addressed annotation pitfall.

2. **Quantified demonstration that un-filtered non-coding RNA fraction
   confounds weighted UTR / APA summary statistics, with correlation to
   differentiation state**: SCAPTURE observed the underlying phenomenon
   (non-coding genes show 3' read pileup) but did not quantify its downstream
   effect on summary statistics or its relationship to cell state. This appears
   to be IsoDecipher's most defensible original contribution from today's review.

3. **N-point (not just proximal/distal) isoform group framework**: scDaPars/
   DaPars's PDUI is fundamentally a 2-point design choice. IsoDecipher's
   panel/group structure natively supports genes with arbitrary numbers of
   APA sites (e.g., ST6GAL1 with 9 groups) — should be framed as a structural
   capability difference, not a "superiority" claim, since 2-point PDUI is a
   deliberate simplification with its own justification in those papers.

## 6. Cautions for Writing
- Do NOT claim "first" or "only" without exhaustive literature search beyond
  today's queries — these are promising leads, not exhaustive proof of novelty.
- SCAPE's own paper's self-reported superiority claims should be cited as
  "SCAPE reported..." and cross-referenced against the independent Bi et al.
  2026 benchmark, which found notably weaker precision for SCAPE.
- Avoid marketing language ("complete victory", "sole immunity", etc.) — state
  findings as specific, falsifiable technical claims with citations.

---

## 8. Entropy / Depth-Confound Problem — Literature Context and Open Strategies
*Connects the entropy instability findings from earlier in this session
(library-size confound, Chao-Shen/Miller-Madow limitations, DonorD depth
outlier) to today's literature review. Two candidate strategies remain open
and unimplemented.*

### 8.1 This is a known, benchmarked, field-wide problem — not unique to IsoDecipher
Section 4.1 above already establishes that pA/isoform identification sensitivity
is depth-dependent across essentially all tools (Bi et al. 2026 downsampling
experiment). The IsoDecipher-specific finding goes one level deeper: even with
a fixed, GTF-anchored set of candidate sites (i.e., even when site IDENTIFICATION
is not the bottleneck), the downstream entropy SUMMARY STATISTIC computed from
soft-assigned read counts remains depth-confounded at the per-cell level
(library size vs. entropy: R=0.38-0.44 in shallow day7 samples; DonorD's median
library size of 13584 vs. 7000-9000 in other Day13 batches drives a
disproportionate share of its elevated entropy). This is a distinct failure
mode from "missing the site entirely" — it is "detecting the right sites but
mis-estimating their relative usage due to small per-cell, per-gene sample size."
No tool reviewed today (SCAPE, scUTRquant, MAAPER, scDaPars, SCAPTURE) was found
to explicitly address this estimator-level bias; the closest related concept
found was the Bi et al. 2026 observation that low-abundance/low-stability sites
detected by few methods inflate group-level MAPE — i.e., independent confirmation
that low-count APA quantification is inherently noisy across the field, motivating
aggregation-based correction strategies (see 4.2).

### 8.2 Strategy 1 — Stratified (per-experiment) consensus feature selection
*Discussed but not yet validated end-to-end.*

Core idea: rather than selecting entropy-valid genes globally (which lets the
deepest-sequenced experiment, exp106/DonorD, dominate the ranking), compute
entropy independently within each of the 3 sequencing runs (exp97, exp105,
exp106), take the top-N genes by entropy variance from each run separately,
and use the union (or intersection) as a "consensus" gene set before any
depth correction is applied.

Rationale: this is conceptually the single-cell analogue of Seurat's
`SelectIntegrationFeatures` (selecting features that are informative within
each batch rather than letting one batch's scale dominate) — a legitimate,
established strategy category (stratified, not supervised, since no cell-type
label is used, only experiment/batch identity).

**Open concerns not yet resolved**:
- Per-experiment top-N selection may still itself be library-size-driven
  *within* each experiment (needs to be checked: does entropy ranking within
  a single exp still correlate with that exp's internal library size
  distribution?).
- Union strategy risks including genes that are only reliably estimable in
  one of the three experiments — being in the "consensus" list doesn't
  automatically make a gene's entropy estimate reliable in the experiments
  where it wasn't independently selected.
- Still requires a downstream correction step (Chao-Shen or downsampling) on
  the resulting consensus set — itself only partially effective per the
  Chao-Shen validation already performed (real reduction in DonorC/DonorD gap,
  but not full elimination; works well for very high-count genes like MALAT1,
  poorly for typical mid/low-count genes like TMBIM6).

### 8.3 Strategy 2 — Micro-pseudobulk pooling
*Discussed but not yet validated end-to-end. See Section 7.2 for the
Figure-5-velocity-specific version of this idea (KNN-neighbor pooling);
this section covers the more general, non-velocity application — i.e.,
using pseudobulk purely to get a usable Global APA Entropy value per
cluster/region, independent of any vector-field/velocity goal.*

Core idea: aggregate counts from a pool of k cells (e.g., k neighboring
cells in a validated, batch-corrected embedding, or k cells within the
same cluster/leiden group) into a single "micro-bulk" pseudo-sample before
computing entropy. This directly targets the root statistical problem
(per-cell, per-gene total counts of ~4-7 for most entropy-valid genes are
too low for any single-cell-resolution entropy estimator, corrected or not,
to be reliable) by trading single-cell resolution for estimator reliability.

**Relationship to literature**: Section 4.2 (Bi et al. 2026) provides indirect
support — their simulated-data result showing that restricting to
multi-method consensus (i.e., more reliably-detected, implicitly higher-count)
sites brings group-level error down to barcode-level error is consistent with
the general principle that aggregation reduces estimator variance for
low-abundance features. No tool reviewed today implements pseudobulk pooling
specifically for entropy estimation in scRNA-seq APA analysis; scDaPars's NNLS
neighbor-borrowing imputation (Section 2, scDaPars profile) is the closest
related idea but is an imputation strategy for filling dropout, not an
aggregation strategy for stabilizing variance, and operates on PDUI
(2-point) rather than N-point entropy.

**Open concerns not yet resolved** (shared with Section 7.2):
- Choice of k, grouping strategy (KNN vs. cluster-based vs. other), and
  avoiding re-introduction of batch confound via the pooling step itself
  — all unresolved design questions, see 7.2 for full discussion.
- Unlike the Figure-5 velocity use case, a non-velocity micro-pseudobulk
  entropy value would most naturally be reported at the cluster/region level
  (e.g., "Global APA Entropy per Leiden cluster, pseudobulked") rather than
  per-cell — this needs to be made explicit when used, since it changes the
  resolution/interpretation of any entropy comparison (cluster-level claims,
  not single-cell claims).

### 8.4 Status: both strategies remain candidate solutions, not yet decided
As of this session, Figure 3's entropy treatment is planned to be the
conservative, already-validated framing: report entropy's depth-confound
as a methodological finding/caveat, anchor the main differentiation-axis
narrative on protein-coding-only weighted UTR (robust, validated, IGH-fix
applied), and defer any of Sections 8.2/8.3's correction strategies to
future work or a clearly-labeled supplementary/exploratory analysis —
unless one of the two strategies above is implemented and validated
(showing convincing removal of the library-size correlation AND preservation
of genuine biological signal) before the manuscript is finalized.
*Discussion notes from 2026-06-19, for later implementation after Figure 4
(C_state / robust cell ordering) is complete. Not yet implemented — record only.*

### 7.1 Candidate velocity metrics

**(a) UTR velocity — ROBUST, highest priority**
- d(weighted_utr) / d(differentiation axis)
- Uses the already-validated protein-coding-only weighted UTR metric
  (robust to depth, unaffected by non-coding RNA confound after the
  utr_source/IG_C_gene fix)
- Should form the core/anchor panel of Figure 5

**(b) Non-coding RNA fraction velocity — ROBUST, second priority**
- d(noncoding_fraction) / d(differentiation axis)
- Motivated by the R=-0.624 correlation (noncoding fraction vs. bcell_score)
  discovered while debugging the weighted-UTR confound
- Independent axis from (a): captures global transcriptome shift toward
  non-coding-RNA-heavy programs (MALAT1/NEAT1 upregulation in terminally
  differentiated/secretory cells is consistent with known biology), not
  structural 3'UTR shortening per se
- Should NOT be highly correlated with (a) if it's adding independent signal;
  worth checking this empirically before/after building the vector field

**(c) Entropy velocity — RISKY, requires fix before use**
- Naive per-cell entropy is NOT stable enough to differentiate directly
  (noise gets amplified, not smoothed, by taking a derivative)
- Two candidate fixes discussed (both deferred to Figure 5 implementation):
  - **Pseudobulk/microbulk pooling**: aggregate k neighboring cells'
    isoform counts before computing entropy, raising per-gene total count
    into a statistically reliable range (c.f. the Chao-Shen test showing
    MALAT1 at count~141-220 corrects reliably vs. TMBIM6 at count~4-7 not)
  - **KNN-connectivity-based local gradient (preferred over pseudotime)**:
    instead of relying on a single global ordering (DPT was abandoned due
    to batch contamination in the PCA space), define "velocity direction"
    as the local gradient of (microbulk-corrected) entropy across a cell's
    k-nearest-neighbor graph. This avoids needing one canonical pseudotime
    axis and is closer in spirit to how scVelo embeds velocity vectors
    onto a UMAP via neighbor relationships.

### 7.2 Microbulk pooling — design considerations (not yet implemented)

Open design questions to resolve before implementation:
1. **Grouping cells for pooling**:
   - Random grouping: simple but destroys real cell-to-cell heterogeneity,
     loses resolution along the trajectory — not suitable for velocity.
   - Sliding window along pseudotime/differentiation score: preserves
     trajectory resolution, but reintroduces dependency on having a
     reliable ordering axis (circular with the DPT problem this is meant
     to avoid, unless the axis used is bcell_score or protein-coding
     weighted_utr, which are already validated as robust).
   - **KNN-neighbor-based grouping (preferred)**: pool each cell with its
     k nearest neighbors in a chosen embedding. Avoids needing a single
     global axis.
2. **Avoiding re-introducing batch confound via pooling**:
   - Must stratify by batch BEFORE pooling, or use neighbor relationships
     from a BBKNN-corrected graph (already batch-aware by construction),
     to avoid mixing high-depth (DonorD) and low-depth (DonorA/B) cells
     into the same microbulk, which would create new, harder-to-interpret
     confounds.
   - Do NOT build the KNN graph from the ISO Entropy UMAP itself (the
     embedding used to select entropy-valid genes) — this would be
     circular (using entropy-derived structure to "fix" entropy).
     Use the GEX UMAP's BBKNN graph instead (validated to have more even
     batch distribution).
3. **Choosing k (pool size)**:
   - Too small: doesn't solve the underlying sparsity problem.
   - Too large: loses the time/state resolution needed for "velocity"
     to mean anything (a microbulk spanning early-to-late cells can't
     represent an instantaneous rate).
   - Rough estimate: per-gene single-cell median count for entropy-valid
     genes is ~4-7; to reach a Chao-Shen-reliable range (~140-220, per the
     MALAT1 test) would suggest k on the order of 20-30 cells, though this
     should be empirically tuned per-gene or per-gene-class rather than
     fixed globally.

### 7.3 "Decision entropy" — undefined, needs clarification
Mentioned in discussion as a candidate Figure 5 metric but not yet defined.
If it refers to "uncertainty in isoform assignment for a given read/cell,"
it reduces to a variant of entropy and will face the same low-count
reliability issues. If it refers to something structurally different (e.g.
a classification confidence margin), it needs a precise definition before
being treated as implementable.

### 7.4 Priority ordering for Figure 5 implementation
1. UTR velocity (robust, ready once Figure 4 ordering/embedding is settled)
2. Non-coding fraction velocity (robust, same dependency)
3. KNN-based microbulk entropy velocity (exploratory; depends on resolving
   open design questions in 7.2, and ideally validated against the UTR
   velocity field for consistency before being presented as a primary result)
4. "Decision entropy" — deferred until defined


## 9. Major Update — Comprehensive Review Paper + PASTA + scPAISO + scUTRquant clustering details
*Added 2026-06-19, second research session. Source: Fahmi, Saha, Song, Lou,
Yong, Zhang. "Computational methods for alternative polyadenylation and
splicing in post-transcriptional gene regulation." Experimental & Molecular
Medicine 57(8), 2025. https://doi.org/10.1038/s12276-025-01496-z (Open Access,
CC BY 4.0) — plus PASTA, scPAISO, and scUTRquant/mca-utrome follow-ups.*

### 9.1 Review paper Table 1 — Bulk RNA-seq APA methods (complete list, 27 tools)
For IsoCAPE (de novo) reference. Columns: Name / APA location / Approach /
Language / Novel APA / Differential APA / Year.

| Name | Location | Approach | Lang | Novel | Diff | Year |
|---|---|---|---|---|---|---|
| 3'T-fill | UTR-APA | New experimental technique | — | Yes | No | 2013 |
| Naive Bayes classifier | UTR-APA | ML site prediction | Perl | Yes | No | 2013 |
| change-point | UTR-APA | Read density change detection | Java | No | Yes | 2014 |
| DaPars | UTR-APA | Modeling read density changes | Python | Yes | Yes | 2014 |
| GETUTR | UTR-APA | Modeling read density changes | Python | Yes | No | 2015 |
| IsoSCM | UTR-APA | Modeling read density changes | Java | Yes | Yes | 2015 |
| Roar | UTR-APA | Read density change (2-site only) | R | No | Yes | 2016 |
| Omni-PolyA | site ID | ML prediction | Webtool | Yes | No | 2017 |
| QAPA | UTR-APA | Annotated-site read density | Python, R | No | Yes | 2018 |
| PAQR_KAPAC | UTR-APA | Annotated-site read density | Python, R | No | Yes | 2018 |
| APAtrap | UTR-APA | Sliding window + MSE model | R, Perl | Yes | Yes | 2018 |
| IntMAP | UTR-APA | Modeling read density changes | MATLAB | Yes | Yes | 2018 |
| TAPAS | UTR-APA, IPA | Pruned Exact Linear Time + read coverage | Shell, R | Yes | Yes | 2018 |
| APARENT | site ID | DL (CNN) site prediction | Python | Yes | No | 2019 |
| DeeReCT-PolyA | site ID | DL (CNN) site prediction | Python | Yes | No | 2019 |
| mountainClimber | UTR-APA, IPA | Read density change | Python | Yes | Yes | 2019 |
| SANPolyA | site ID | DL (self-attention) | Python | Yes | No | 2020 |
| APAlyzer | UTR-APA, IPA | Annotated-site, Fisher/t-test | R | No | Yes | 2020 |
| PolyA-Miner | UTR-APA | 3'-end-seq peak modeling | Python | Yes | Yes | 2020 |
| Aptardi | UTR-APA | DL (biLSTM) site prediction | Shell, Python | Yes | Yes | 2021 |
| MAAPER | UTR-APA, IPA | Modeling read density changes | R | Yes | Yes | 2021 |
| IPAFinder | IPA | Read density MSE ratio | Python | Yes | Yes | 2021 |
| DeeRect-APA | site ID | DL site usage estimation | Python | Yes | No | 2022 |
| APA-Scan | UTR-APA | Read density change | Python | Yes | Yes | 2022 |
| APAIQ | UTR-APA, IPA | DL site prediction | Python | Yes | Yes | 2023 |
| InPACT | IPA | CNN + Salmon + Dirichlet-multinomial (DRIMSeq) | Python | Yes | Yes | 2024 |
| PolyAMiner-Bulk | UTR-APA | BERT (C/PAS-BERT) + beta-binomial | Python | Yes | Yes | 2024 |

Key mechanistic notes from this table:
- **Roar and change-point both assume only TWO poly(A) sites per 3'UTR** —
  this is the same 2-point limitation discussed re: DaPars/PDUI (Section 2,
  scDaPars profile). Confirms 2-point assumption is a widespread, not isolated,
  design pattern in this field — useful context for framing IsoDecipher's
  N-point design as a structural choice rather than an isolated quirk others missed.
- **InPACT** uses Dirichlet-multinomial (via DRIMSeq) specifically for IPA
  differential testing — third confirmed use of Dirichlet-multinomial in this
  space (alongside PASTA, Section 9.4, and scraps, Section 9.5), suggesting
  this is becoming a standard statistical choice for polyA/isoform count
  overdispersion modeling.

### 9.2 Review paper Table 2 — Single-cell APA methods (complete list, 8 tools)
This is the authoritative, independently-compiled list to check IsoDecipher's
positioning against. Columns: Name / APA location / Approach / Language /
Novel APA / Differential APA / Year.

| Name | Location | Approach | Lang | Novel | Diff | Year |
|---|---|---|---|---|---|---|
| scAPA | UTR-APA, IPA | Peak calling | R | Yes | Yes | 2019 |
| Sierra | UTR-APA, IPA | Peak calling | R | Yes | Yes | 2020 |
| scAPAtrap | UTR-APA, IPA | Peak calling | R | Yes | Yes | 2021 |
| scDAPA | UTR-APA | Read density change (histogram-based, see 9.3) | Shell, R | Yes | Yes | 2020 |
| SCAPTURE | UTR-APA, IPA | Peak calling + DL (DeepPASS) | Shell, R, Python | Yes | Yes | 2021 |
| SAPAS | UTR-APA, IPA | Peak calling | Python, R | Yes | Yes | 2021 |
| scDaPars | UTR-APA | Read density modeling (PDUI + NNLS imputation) | R | Yes | Yes | 2021 |
| SCAPE | UTR-APA | Peak calling (probabilistic mixture model) | Python, R | Yes | Yes | 2022 |

**Important gap check**: this review's Table 2 does NOT include scUTRquant,
MAAPER (listed separately in Table 1 as bulk-capable but applicable to
single-cell too per Bi et al. 2026), scPAISO (too new, post-dates this review),
or scAPAmod. The review's single-cell table appears less complete than the
Bi et al. 2026 benchmark's method list — cross-reference both, don't rely on
either alone as the canonical list.

**Also notable**: every single tool in Table 2 is marked "Yes" for both Novel
APA and Differential APA — i.e., the entire published single-cell APA
landscape (per this review) is de novo-capable de facto. This sharpens
IsoDecipher's positioning: an annotation-anchored, non-de-novo design is
the less common choice in the single-cell space specifically (more common
in bulk, per Table 1's "No" entries for Novel APA: Roar, change-point,
QAPA, PAQR_KAPAC, APAlyzer). This is a meaningful framing point — IsoDecipher
imports a bulk-RNA-seq-style annotation-anchored design philosophy into the
single-cell domain, where it is comparatively rare.

### 9.3 scDAPA — newly identified, distinct from scDaPars (note capitalization)
"scDAPA" (Table 2) is a DIFFERENT tool from "scDaPars" (Section 2 profile) —
easy to confuse due to near-identical names. scDAPA (2020): "uses a
histogram-based method to bin 3'-ends and calculates the site distribution
difference index to identify APA changes between cell groups... results
visualized through smooth density plots and isoform profiles." Distinct
algorithmic approach (histogram binning + distribution difference index)
vs. scDaPars (regression + PDUI + NNLS imputation). Not yet deep-dived;
flag for follow-up if a detailed comparison is needed.

### 9.4 PASTA (Satija Lab / Seurat ecosystem)
Kowalski*, Wessels*, Linder* et al. Cell 2024 (DOI confirmed via cell.com).
GitHub: satijalab/PASTA. "PolyA Site analysis using relative Transcript Abundance."

- **Core method**: Dirichlet-multinomial distribution to model a BACKGROUND
  distribution of polyA site usage (in control/NT cells), explicitly
  **"controlling for gene expression."** Computes "polyA-residuals" — per-cell,
  per-site over/under-utilization relative to this background model.
- **Key statistical innovation directly relevant to IsoDecipher's entropy
  problem**: "parameterize overdispersion estimates individually for each
  polyA site but then REGULARIZE these estimates ACROSS SIMILAR SITES" —
  i.e., a formal Bayesian/empirical-Bayes-style shrinkage approach to
  stabilize low-count site-level estimates by borrowing statistical strength
  from similar sites. This is a more principled (if more complex to
  implement) alternative to the ad hoc Chao-Shen/consensus-filtering/
  microbulk strategies discussed in Section 8 — explicitly designed to solve
  the same underlying problem (unreliable low-count single-cell APA estimates)
  that Sections 8.2/8.3 are trying to address.
- Analogized in the source text to sctransform's per-gene-then-regularized
  overdispersion modeling for standard gene expression — i.e., this is
  philosophically the "do for APA what sctransform does for gene counts"
  approach.
- Fully integrated into Seurat (R ecosystem) — not directly portable to
  IsoDecipher's Python/Scanpy stack without reimplementation, but the
  *statistical strategy* (per-feature overdispersion + cross-feature
  regularization) is implementable independently of the R/Seurat machinery.
- **Action item for Section 8 (entropy strategies)**: consider PASTA's
  regularization approach as a third candidate strategy (alongside
  stratified-consensus and microbulk), potentially superior in principle
  since it doesn't require discarding single-cell resolution (unlike
  pseudobulk) or making ad hoc per-experiment splits (unlike stratified
  consensus) — but would require new implementation work, not adaptation of
  an existing Python tool.

### 9.5 scraps — additional detail
(Brief tool already listed in Section 2; expanding here.)
- Supports **three input modes**, in order of preference: (1) paired
  alignments (preferred — analogous in spirit to scPAISO's Read1+Read2
  approach, though scraps does not appear to mandate it), (2) existing Cell
  Ranger BAMs, (3) Read2-only alignment. This flexibility (graceful
  degradation to Read2-only) is more accommodating of varied input data than
  scPAISO's apparent hard requirement for sufficiently long Read1 cDNA content.
- Demonstrated workflow includes UMAP panels for a single gene (NR1D1) showing
  proximal/distal site counts AND their Dirichlet-multinomial residuals
  side-by-side — again confirming Dirichlet-multinomial as a recurring
  statistical choice (4th confirmed instance: PASTA, InPACT/DRIMSeq, scraps,
  and likely related to BANDITS in the AS literature, Table 3).

### 9.6 scPAISO (2025, bioRxiv 2025.08.20.669565) — Read1-based approach
GitHub: yongjieliu/scPAISO. Not yet in any benchmark (too new — postdates
Bi et al. 2026's data collection and the Fahmi et al. 2025 review).

**Core innovation — uses Read1, not Read2, for PAS identification**:
- Standard 10x Chromium 3' library: Read1 = 16bp barcode + 10/12bp UMI +
  30bp polyT + ~90bp cDNA (under PE150 sequencing); Read2 = 150bp cDNA.
  Read1's cDNA portion is conventionally DISCARDED after barcode/UMI
  extraction by virtually all other tools.
- Read1's 5' end sits directly at the true cleavage site, yielding "sharp
  peaks" (half of detected PASs span <20bp; 95% <70bp) — i.e., near-
  nucleotide-resolution site calls, in contrast to Read2-based methods'
  diffuse, statistically-inferred site locations (which is exactly the
  problem IsoDecipher's KDE/insert-size approach and SCAPE's Gaussian
  mixture model are both built to handle).
- Method: STAR-align Read1 separately (with specific clip/filter params),
  extract 5' end positions, MACS3 peak calling, internal-priming filter
  (6-mer AAAAAA check, same general strategy as IsoDecipher's
  `--filter-polya-tract`), merge sites within 50bp, classify as tPA
  (terminal/3'UTR) / iPA (intronic) / ePA (exonic). Read2 is then assigned
  to the nearest PAS using a fitted distance model (best-fit among Weibull/
  Gamma/Log-normal, chosen empirically per dataset) — conceptually similar
  in spirit to IsoDecipher's per-sample KDE calibration, but using a fitted
  parametric distribution rather than a non-parametric KDE, and using
  Read1-derived PASs as the reference rather than GTF-derived ones.

**Critical limitation (important for honest comparison)**: this approach
fundamentally REQUIRES a sequencing configuration where Read1 contains
sufficient cDNA beyond the barcode/UMI/polyT (e.g., PE150 with ~90bp of
usable cDNA in Read1). Many standard 10x configurations (e.g., shorter
Read1, such as 26bp covering only barcode+UMI) cannot provide this signal at
all. This is a hard data-availability constraint, not a tunable parameter —
should be framed as "scPAISO exploits a sequencing-configuration-dependent
information source unavailable in many standard datasets," not as a strict
superiority claim, when discussed in IsoDecipher's related work.

**Classification scheme (tPA/iPA/ePA)** — more granular than IsoDecipher's
current G0/G1/G2... index, which encodes relative distance-from-gene-body
ordering but does not explicitly label genomic context (3'UTR vs intron vs
exon). Worth considering as a potential `iso.var` metadata extension
(e.g., a `site_type` column) for a future IsoDecipher version, separate
from the existing `utr_source` (coding/non-coding) column.

**"Rank score" metric**: 0 (exclusive proximal) to 1 (exclusive distal),
computed as a weighted sum using the PAS's ORDINAL POSITION (not physical
bp distance) as the weight: "multiplied the estimate by the ordinal number
of the PAS." This is a meaningful methodological contrast with IsoDecipher's
`weighted_utr`, which weights by actual physical UTR length in bp. Ordinal
weighting treats unevenly-spaced sites as equally spaced; physical-length
weighting (IsoDecipher's choice) preserves true spacing information. Worth
stating explicitly as a design trade-off in any related-work comparison
(physical-distance weighting is more information-preserving when site
spacing is uneven, which is common — e.g., recall the IsoDecipher panel
data showing wildly uneven inter-site spacing within single genes like
ST6GAL1 or PRDM1).

**Notable validated biological findings from scPAISO** (useful as cross-
validation precedent for IsoDecipher's own B-cell/PC differentiation
findings): in hematopoietic differentiation, PA-isoform-level differential
expression captured changes invisible at the gene level (e.g., HAT1, TAL1,
TEN1, WDR33 — genes with no gene-level DE but significant PA-isoform-level
DE), directly paralleling IsoDecipher's repeated finding that ISO-space
structure reveals signal absent from GEX-space (e.g., the IFN+ plasma cell
subpopulation, MALAT1 isoform-level signal). This is useful supporting
precedent (independent confirmation of the general phenomenon "APA reveals
state invisible to gene expression") even though scPAISO is methodologically
very different from IsoDecipher.

**Uses pseudobulk + k-means** for global UTR-length-pattern discovery across
cell types — another independent confirmation that pseudobulk aggregation is
a standard, accepted strategy for stabilizing low-count APA signal at the
cost of single-cell resolution (further supporting Section 8.3's micro-
pseudobulk strategy as methodologically conventional, not ad hoc).

### 9.7 scUTRquant / mca-utrome — clustering parameter details (resolves
open question from Section 2)
GitHub: Mayrlab/mca-utrome (companion repo to scUTRquant, builds the
reference UTRome index). Manuscript: Fansler, Mitschka & Mayr, Nat Commun
15:4050 (2024).

Confirmed: scUTRquant/mca-utrome DOES perform a clustering/merging step
analogous to IsoDecipher's consecutive-gap clustering, via TWO separate
tunable parameters:
- **epsilon (e)**: "distance for initial merging of cleavage sites" — used
  when building the reference UTRome database (one-time, applies to all
  downstream samples using that reference). Observed default in a released
  index: e=30 (i.e., 30bp), from the filename `utrome.e30.t5.gc25.pas3.f0.9999.w500`.
- **mergeDistance (m)**: "distance within which to merge isoforms when
  quantifying in scUTRquant (kallisto)" — a second, separate merge step
  applied at quantification time.
- **truncation width (w)**: "maximum transcript distance to cleavage site"
  — observed default w=500 (500bp), conceptually adjacent to (but distinct
  in purpose from) IsoDecipher's 420bp read-fetch window (IsoDecipher's
  420bp is derived from the empirical P95 insert-size offset; scUTRquant's
  500bp appears to be a transcript-truncation parameter for index-building,
  not an empirically-derived window).

**Key comparison point (now resolvable)**: no documentation found describing
HOW the e=30bp default was chosen — i.e., no equivalent of IsoDecipher's
Bhattacharyya Coefficient validation across 6 independent datasets that
empirically grounds the 125bp tolerance choice. scUTRquant's epsilon appears
to be a fixed, reference-database-level choice (set once when building the
UTRome, not re-derived per new dataset), whereas IsoDecipher's per-sample
KDE calibration and BC-validated tolerance are explicitly designed to be
dataset-specific and re-derivable. This remains a legitimate, specific
differentiator: **IsoDecipher's tolerance threshold is empirically validated
and dataset-adaptive; scUTRquant's merge parameters are fixed values
embedded in a precompiled reference index.**

### 9.8 Isoform naming convention — field survey
Surveyed naming conventions across reviewed tools:
- **scDaPars/DaPars, 3UTR-Quant, Roar, change-point**: binary `proximal`/
  `distal` labels (column names like `proximal_id`, `distal_id`). Inherently
  assumes (or is restricted to) exactly 2 sites per gene.
- **scPAISO**: `tPA`/`iPA`/`ePA` — categorical by genomic context (terminal
  3'UTR / intronic / exonic), not by ordinal position; a gene can have
  multiple PASs within the same category.
- **IsoDecipher**: `G0, G1, G2...` — numeric, ordinal, by increasing distance
  from gene body (G0 = most proximal). Natively supports arbitrary N sites
  per gene without naming ambiguity (e.g., ST6GAL1's 9 groups, DDX5's 7
  groups render unambiguously; a binary proximal/distal scheme cannot
  represent these without an awkward ad hoc extension).
- No other tool surveyed today was found using a comparable N-ary ordinal
  naming scheme — this appears to be a genuine, if modest, design difference
  worth noting in Methods (not necessarily framed as a major contribution,
  but as a deliberate structural choice enabling the N-point analyses
  IsoDecipher performs that 2-point-designed tools like PDUI/RUD cannot).

### 9.9 Updated synthesis — where IsoDecipher sits in the landscape
Combining all findings (Sections 1-9):
- IsoDecipher is annotation-anchored (GTF-based), placing it among the
  minority of single-cell APA tools using this design (per 9.2, all 8
  reviewed single-cell tools are de novo-capable; closest annotation-anchored
  precedents — scUTRquant, MAAPER — are in the broader literature but not
  in this particular review's single-cell table).
- IsoDecipher's N-point group structure (vs. 2-point PDUI/RUD/rank-score
  designs used by DaPars/scDaPars/Roar/change-point/scPAISO) is a genuine,
  consistently-observed structural differentiator across every 2-point tool
  surveyed.
- IsoDecipher's empirically-validated (BC-based, cross-dataset), per-sample-
  adaptive tolerance threshold contrasts with scUTRquant's fixed, reference-
  index-embedded merge distance — a real, specific, citable difference.
- IsoDecipher's IG_C_gene/TR_C_gene biotype fix and the quantified
  non-coding-RNA-contamination-of-weighted-UTR finding (R=-0.624) remain
  the most original, well-evidenced contributions identified in two full
  sessions of literature review — no tool surveyed (including scUTRquant,
  which does filter to protein-coding only) was found to address the
  immunoglobulin-biotype edge case or to quantify the downstream
  differentiation-correlated bias from not filtering.
- Three independent precedents for entropy/usage-estimate stabilization
  strategies now identified for Section 8's open problem: Bi et al. 2026's
  consensus-site aggregation (Section 4.2), PASTA's per-site overdispersion
  + cross-site regularization (Section 9.4), and scPAISO's pseudobulk +
  k-means for pattern discovery (Section 9.6) — none directly reusable
  off-the-shelf for IsoDecipher's Python stack, but all confirming that some
  form of statistical pooling/regularization is the field's accepted answer
  to single-cell APA sparsity, supporting further investment in Section 8.2/
  8.3/PASTA-style strategies as methodologically mainstream rather than ad hoc.

## 10. Per-Sample / Per-BAM Fragment-Size Adaptation — Cross-Tool Comparison
*Added 2026-06-19, third pass. Motivated by IsoCAPE's breast cancer FACS-sorted
dataset, which showed a dramatic short-fragment-enriched insert size
distribution distinct from other datasets (see Figure 2B six-dataset
comparison). Question: does any existing tool re-calibrate its distance/
insert-size model independently for each input BAM/sample, the way
IsoDecipher's per-sample KDE does?*

### 10.1 scDAPA — confirmed NO per-sample adaptation
Explicitly uses a **fixed, global bin width of 100bp** ("divide the dispersed
3' ends into distinct bins with the same width (default 100 bp)") to compute
the Site Distribution Difference (SDD) index. No mechanism found for
adjusting bin width per BAM/sample based on that sample's actual fragment
size distribution or sequencing/library quality. Bi et al. 2026 (k-mer/file)
text confirms.

### 10.2 scUTRquant / mca-utrome — confirmed NO per-sample adaptation
All three relevant parameters (epsilon=30bp merge distance, mergeDistance(m),
truncation width w=500bp — Section 9.7) are **fixed values baked into the
reference UTRome index at build time**. The same precompiled index/kallisto
target is reused across all downstream samples regardless of each sample's
own fragment size distribution or library quality. No per-BAM recalibration step found.

### 10.3 DaPars/scDaPars, QAPA, PAQR_KAPAC — no evidence of per-sample
insert-size modeling found
These rely on read-density/coverage-based regression to infer site locations
and proportions; no description found of any insert-size-distribution
estimation step (per-sample or otherwise) in any documentation reviewed.
Coverage-based regression methods are, by construction, more directly
exposed to fragment-size-driven artifacts: a sample with an enriched short-
fragment population (e.g., from FACS sorting / degraded RNA, as observed in
IsoCAPE's breast cancer dataset) would directly distort the read-density
profile these methods regress against, with no described correction mechanism.

### 10.4 SCAPE — has per-dataset insert-size modeling, but likely protocol-
level rather than fully per-sample
SCAPE's Bayesian mixture model explicitly estimates the insert size (f)
distribution "from the read pairs... mapped to large constitutive regions
such as 3' UTR" — i.e., it IS estimated from the data at hand, not hardcoded.
**However**, the independent Bi et al. 2026 benchmark found SCAPE "performed
poorly on CEL-seq and Drop-seq data, likely due to mismatched empirical
insert-size parameters (optimized for 10x Chromium/Microwell-seq)" — this
strongly suggests SCAPE's insert-size calibration is tuned/validated at the
**sequencing-protocol level**, not independently re-derived for every input
BAM regardless of sample-specific quality variation (e.g., it is unclear
whether SCAPE would correctly adapt to a degraded/short-fragment-enriched
10x sample the way IsoDecipher's BC-validated, per-BAM KDE approach is
explicitly designed to). This is an important nuance: "estimated from data"
is not the same guarantee as "re-validated per individual sample," and the
cross-protocol failure mode observed empirically suggests the latter,
stronger property may not hold for SCAPE.

### 10.5 scPAISO — closest precedent found, but conditional on a different
data requirement
Explicitly fits multiple candidate distributions (Weibull/Gamma/Log-normal)
to the Read2-to-PAS distance using a positive-control dataset constructed
**from that same sample's own paired Read1/Read2 data**, and selects the
best-fitting model per analysis ("the model with the highest global accuracy
was subsequently used"). This is methodologically the closest precedent
found to IsoDecipher's per-sample KDE recalibration — both re-derive a
distance model from each dataset's own empirical read behavior rather than
applying a fixed global/protocol-level default.
**However**, scPAISO's approach is entirely contingent on having sufficient
Read1 cDNA content (Section 9.6's critical limitation) — a data-availability
constraint that does not apply universally across 10x configurations.
IsoDecipher's per-sample KDE works directly from standard Read2-based insert
size information without requiring any special Read1 length, making it
applicable to a broader range of standard 10x library preparations,
including (as directly demonstrated in Figure 2B) a degraded, short-fragment-
enriched FACS-sorted breast cancer sample.

### 10.6 Summary table

| Tool | Per-sample/per-BAM distance recalibration? | Granularity | Notes |
|---|---|---|---|
| scDAPA | No | Fixed global (100bp bins) | — |
| scUTRquant/mca-utrome | No | Fixed, baked into reference index | epsilon/mergeDistance/w all fixed at index-build time |
| DaPars/scDaPars | No evidence found | — | Coverage-regression-based; directly exposed to fragment-size artifacts |
| QAPA/PAQR_KAPAC | No evidence found | — | Annotated-site read-density based |
| Sierra | No (implied) | Coverage-fitting, depth-sensitive | Per Bi et al. 2026, sharp pA yield drop at low coverage — consistent with no robust per-sample adaptation |
| SCAPE | Partial — estimated from data | Likely protocol-level, not confirmed per-BAM | Independent benchmark shows cross-protocol failure (CEL-seq/Drop-seq), suggesting tuning is not fully sample-adaptive |
| scPAISO (2025) | Yes — closest precedent | Per-sample (uses that sample's own Read1/Read2 pairs) | BUT requires long-Read1 sequencing configuration; not applicable to standard short-Read1 10x runs |
| **IsoDecipher** | **Yes** | **Per-sample (per-BAM independent KDE)** | Validated specifically on a degraded/short-fragment FACS-sorted sample (breast cancer, Figure 2B); works from standard Read2 insert-size info, no special sequencing config required |

### 10.7 Conclusion — this is a genuine, citable differentiator
Across all tools surveyed in this and the prior research session, **no tool
other than IsoDecipher and the newly-identified scPAISO (2025) was found to
re-derive its core distance/insert-size model independently for each input
BAM/sample**, and scPAISO's version is gated behind a non-universal
sequencing-configuration requirement (long Read1) that IsoDecipher's
approach does not need. The remaining tools (the clear majority: scDAPA,
scUTRquant, DaPars-family, QAPA-family, and likely Sierra) all rely on
either a fixed global parameter or a fixed reference-index-embedded
parameter, making them structurally vulnerable to sample-to-sample
variation in fragment size distribution — precisely the kind of variation
demonstrated empirically in IsoCAPE's breast cancer dataset (large excess
of sub-50bp offsets from FACS-sorting-induced degradation) and in
IsoDecipher's own six-dataset BC validation (Figure 2B), where peak insert
sizes ranged from ~25bp (breast cancer) to ~200-250bp (B cell differentiation/
PBMC) across otherwise-comparable 10x Chromium 3' datasets.

This should be framed precisely: not "every other tool fails on degraded
samples" (untested, overly broad), but rather "every other tool's published
methodology relies on a distance model that is fixed in advance (globally
or per sequencing protocol) rather than independently calibrated per input
sample, which is a structural difference from IsoDecipher's approach and
would be expected, on theoretical grounds, to behave less robustly on
samples whose fragment size distribution deviates substantially from
whatever assumption is baked into that fixed parameter." The empirical
breast cancer validation (Figure 2B) demonstrates IsoDecipher's approach
handles this specific, real deviation correctly; a rigorous head-to-head
benchmark against other tools on the same degraded sample (not yet
performed) would be needed to make a direct comparative robustness claim.

## 11. Output Format Ecosystem — AnnData/Scanpy vs R/Bioconductor/Seurat
*Added 2026-06-19, fourth pass.*

### 11.1 Confirmed output formats across reviewed tools

| Tool | Native output format | AnnData/Python support |
|---|---|---|
| scUTRquant | Bioconductor `SingleCellExperiment` (R) | Added later (v0.5.0+) as a secondary option — explicitly described as making it "easier for users who prefer Python and working with the scverse ecosystem," implying it was not the original design target |
| scDaPars | R-style PDUI matrix (genes x cells) | No AnnData support found in reviewed documentation |
| PASTA | Fully Seurat-integrated (R) | None — "interfaces directly with Seurat" |
| scAPA, Sierra, scAPAtrap, SCAPTURE | R (GitHub language profile predominantly R) | Not confirmed; likely Seurat/SCE objects given R-based implementation |
| scPAISO | Uses Seurat v4.0.1 for downstream analysis per Methods | Implies R/Seurat-centric workflow despite some Python components (STAR alignment scripting) |
| SCAPE/SCAPE-APA | Implemented in Python | Specific matrix output format (AnnData vs. plain CSV/TSV) not confirmed in documentation reviewed — would need direct inspection of output files/code to verify |
| **IsoDecipher** | **Native AnnData (.h5ad)**, integrated GEX+ADT+Isoform in one object | **Yes — designed-in from the start** |

### 11.2 Interpretation
The single-cell APA tool landscape is historically rooted in the R/
Bioconductor/Seurat ecosystem (consistent with APA analysis methods'
origins in bulk RNA-seq tools like DaPars, also R-based). Python/AnnData
support, where present at all (scUTRquant), was added as a secondary
accommodation after the fact, not built in from the start. SCAPE is the
only Python-native tool identified, but its precise output data structure
was not confirmed to be AnnData-formatted.

IsoDecipher's AnnData-native design — producing a single object that
directly merges GEX (Cell Ranger), ADT (CITE-seq), and Isoform (soft-
assigned fractional counts) as three feature types in one `.X` matrix
with shared `.obs`/`.var` — is a genuine ecosystem-positioning difference
from the majority of existing tools. This matters practically: a
Scanpy-based analyst working with GEX/ADT data already in AnnData format
can incorporate IsoDecipher's isoform output directly into their existing
pipeline without an R-to-Python conversion step, which most other reviewed
tools would require (e.g., via `sceasy`, `zellkonverter`, or manual export/
import) before integration with a Scanpy-based downstream workflow.

This should be stated as a practical/ecosystem-compatibility advantage,
not a methodological superiority claim — the underlying statistical
methods of R-based tools are not inferior because of their language/
format choice.

---

## 12. Open Question — Do Existing Tools Demonstrate End-to-End Biological
Validation on Their Own Data?
*Added 2026-06-19. Raised as a hypothesis by the user: tool papers may focus
on benchmarking against ground truth (simulated data, matched bulk/3'-seq)
without necessarily presenting a complete, novel biological discovery
story using their own tool end-to-end (upstream BAM processing through
downstream trajectory/differentiation biology). This is a testable claim
requiring direct verification per tool, not an assumption — status of each
tool reviewed below is graded by what was actually found in today's searches,
not inferred.*

### 12.1 What was actually observed per tool (evidence-graded)

**SCAPE** — Has a genuine biological application: "illustrating the global
APA landscape in mouse cell atlas and human glioblastoma," and the original
paper discusses APA dynamics "during erythropoiesis and induced pluripotent
stem cell (iPSC) differentiation." This IS a differentiation/trajectory
biology application using the tool's own pipeline. Not merely a benchmark-
only paper.

**scDaPars** — Demonstrated on "primary breast cancer," distinguishing
"tumor-specific and immune-cell-type-specific APA landscape," with a
described comparison to scAPA/Sierra (Figure 4, A-C: UMAP clustering
comparison) and tumor vs. normal cell PDUI scatter plots. This is a real
biological application (cancer-immune-cell APA heterogeneity), not just
a synthetic benchmark.

**PASTA** — Demonstrated on "a dataset of circulating human peripheral
blood mononuclear cells" and the companion Cell 2024 paper title is
"Multiplexed single-cell characterization of alternative polyadenylation
regulators" — i.e., the paper's PRIMARY contribution appears to be a
perturbation/regulator-screening biological study (likely CRISPR-based,
given "regulators" framing and "NT cells" control terminology seen in
Section 9.4), with PASTA as the supporting statistical tool. This is the
strongest example found of a tool built specifically to enable and report
a novel biological finding, not a benchmarking exercise.

**scPAISO (2025)** — The most extensively validated on real biological
applications among all tools reviewed: hematopoiesis (bone marrow,
GSE196676), systemic sclerosis/SSc (skin, GSE249279), and a multi-tissue
mouse atlas (brain/liver/muscle/skin/hematopoietic). Reports specific
novel findings (CBLB isoform switching in ITP, RBP-3'UTR-length coupling
across tissues, etc.). This is extensive, multi-dataset biological
validation — arguably more extensive than IsoDecipher's current single
B-cell-differentiation dataset scope.

**scUTRquant** — Companion paper title "Quantifying 3'UTR length from
scRNA-seq data reveals changes independent of gene expression" — applied
across "474 cell types and 2,134 perturbations" per the comprehensive
annotation paper (Section 2 profile) — this is also a large-scale, real
biological discovery paper (3'UTR length as an independent regulatory axis
from gene expression), not benchmark-only.

**scDAPA** — Applied to "neuroretinas from mouse," identifying genes with
dynamic APA across cell groups, including a specific example gene
(AI607873) with cell-type-specific proximal/distal site preference. Real
biological application, smaller in scope than the above.

### 12.2 Revised assessment — the original hypothesis is only partially
supported
Based on actual evidence gathered, **most of the major tools reviewed DO
present genuine biological discovery applications using their own
pipelines end-to-end**, not purely benchmark-only papers. PASTA and
scPAISO in particular appear to have built substantial, multi-finding
biological narratives as their primary contribution (not an afterthought).
**The hypothesis that "tool builders don't validate biologically" is not
well supported by what was found today** — it should not be used as a
talking point for IsoDecipher's positioning without significant caveats,
or it risks being an easily-falsifiable claim if a reviewer checks these
same papers.

### 12.3 Where a genuine, evidence-supported gap DOES appear to exist
What today's review does NOT show any other tool doing, specifically:
1. **End-to-end ownership from raw BAM through a fully-integrated
   multi-modal (GEX+ADT+Isoform) AnnData object, in a single coherent
   pipeline with Scanpy-native downstream analysis** — most tools either
   (a) require external GEX quantification (Cell Ranger) and only handle
   the isoform layer themselves, output to R, with multi-modal integration
   left to the user, or (b) are R/Seurat-native and don't natively handle
   CITE-seq/ADT integration alongside isoform-level APA in one object.
2. **Demonstrated robustness to a degraded/short-fragment clinical
   sample via per-sample distance-model recalibration** (Section 10) —
   this remains a genuinely distinctive, evidenced claim.
3. **The specific combination of (a) GTF-anchored panel construction with
   empirically-validated tolerance (BC method), (b) per-sample KDE
   calibration, (c) AnnData-native multi-omic integration, and (d)
   explicit immunoglobulin/biotype-aware UTR filtering** — no single
   other tool was found combining all four of these properties; each
   individual property has at least partial precedent elsewhere (Section 9.9),
   but the combination, plus the IG_C_gene fix specifically, remains
   IsoDecipher's most defensible distinguishing claim.

### 12.4 Implication for positioning strategy
Given 12.2, IsoDecipher's competitive narrative should NOT lean on "we are
the only tool that validates on real biology" (false/risky claim). It
should instead lean on the verified, narrower claims from Sections 9-11:
annotation-anchored design in a field dominated by de novo single-cell
methods; N-point group structure vs. prevailing 2-point PDUI/RUD designs;
empirically-derived (not fixed) tolerance/distance parameters re-validated
per sample; AnnData-native multi-omic (GEX+ADT+Isoform) integration in an
ecosystem otherwise dominated by R/Bioconductor; and the specific,
quantified IG_C_gene/non-coding-RNA biotype contamination finding. A
useful framing for a forthcoming B-cell-differentiation results section
is not "no one else does biology" but rather "this specific combination of
design choices, applied here to plasma cell differentiation, surfaces a
mechanistic finding (e.g., stress/UPR-associated entropy elevation, the
non-coding-RNA UTR confound) that the field's existing toolkit — built on
fixed-parameter or 2-point-restricted designs — would not straightforwardly
have caught."

## 13. SCAPE Deep-Dive — Insert Size Modeling Mechanics (Full Methods Detail)
*Added 2026-06-19, fifth pass. Source: Zhou et al. 2022, NAR 50(11):e66,
full Methods section (direct fetch, not abstract-level).*

### 13.1 How SCAPE actually estimates its insert size distribution
Critical detail not visible from the abstract/intro: SCAPE's insert size
distribution is NOT estimated from reads at the genes of interest. It is
estimated from a separate, restricted subset of the genome:

> "We select genes with only **one exon** and pick **cleavage reads**
> mapped to these one-exon genes, which excludes the biases brought by
> splicing in insert size calculation... Distances between R1 and R2 of
> these reads, plus the length of poly(A) part in R1, provide the insert
> size distribution."

- "Cleavage reads" = reads where **both R1 and R2 uniquely map and directly
  span the cleavage site** — described in the paper as "a small proportion
  of reads."
- Restricting to single-exon genes is necessary to avoid intron-driven
  distance artifacts contaminating the insert-size estimate.
- **Practical implication**: SCAPE's insert-size distribution is estimated
  from a genome-wide pool of single-exon genes' cleavage reads — NOT from
  the specific genes under analysis, and NOT (necessarily) from a sample-
  representative cross-section of all gene/transcript types. This is a
  meaningfully different sampling strategy from IsoDecipher's, which derives
  its tolerance-validation KDE from PolyASite v2.0's curated single-PAS
  genes (n=5,728, chosen for unambiguous ground-truth status) and its
  operational KDE directly from the panel's own annotated sites — both
  IsoDecipher KDEs are computed from BAM-wide read behavior at the actual
  panel genes (or a principled proxy set), not restricted to a structurally
  convenient subset (single-exon genes) that may not be representative of
  the genes of biological interest.

### 13.2 Parametric (single Gaussian) vs. non-parametric (KDE) — the
critical mathematical difference
SCAPE's insert size model is a **single Gaussian distribution** with two
scalar parameters:
> "p(x_n | s_n, θ_nk) = N(x_n | θ_nk + s_n + 1 − μ_f, σ_f²)"
> "where μ_f and σ_f are the mean and standard deviation of fragment
> length distribution (Gaussian)."

Empirically estimated value reported in the paper: **"cDNA fragments are
around 300bp with 50bp standard deviation"** (a single point estimate,
not a full empirical distribution) — and in the simulation studies, fragment
length mean/std were drawn from ranges [250,350] and [20,40] respectively,
i.e., the entire insert-size behavior is collapsed to two numbers (μ_f, σ_f)
per analysis run.

**This is the single most important mechanistic contrast with IsoDecipher**:
IsoDecipher's insert-size model is a **non-parametric KDE** — the full
empirical shape of the insert size distribution is preserved, including
multi-modality or skew. SCAPE's single-Gaussian assumption cannot represent:
- **Bimodal distributions** — exactly what was empirically observed in
  IsoCAPE's breast cancer (FACS-sorted) dataset (Figure 2B), where a large
  population of sub-50bp fragments (degradation/short-fragment artifact)
  coexists with a normal ~180-250bp population. A single Gaussian fit to
  this data would either average across both modes (badly mis-locating
  both populations) or be dominated by whichever mode has more reads,
  systematically mis-assigning reads from the minority mode.
- **Skewed/heavy-tailed distributions** — any departure from symmetric
  Gaussian shape, which the empirical insert-size curves in Figure 2B show
  to varying degrees across all six validation datasets even in the
  "normal" (non-degraded) cases.

This directly and mechanistically explains the independent Bi et al. 2026
benchmark finding (Section 2, SCAPE profile) that **SCAPE "performed poorly
on CEL-seq and Drop-seq data, likely due to mismatched empirical insert-size
parameters (optimized for 10x Chromium/Microwell-seq)"** — a single Gaussian
(μ_f, σ_f) tuned to behave well for one protocol's typical insert-size shape
has no mechanism to adapt to a differently-shaped distribution in another
protocol, because the model family itself (Gaussian) cannot represent
arbitrary shapes regardless of how the two parameters are re-fit.

**This is now a fully evidenced, mechanistically-explained, citable
differentiator** — not merely an assumption. IsoDecipher's non-parametric
KDE is structurally capable of representing the kind of degraded/bimodal
insert-size distribution observed in the breast cancer dataset; SCAPE's
single-Gaussian model family is not, by mathematical construction,
regardless of how its two parameters are estimated.

### 13.3 BIC-driven component selection — confirms depth-dependent
complexity detection is a known, expected statistical property
SCAPE uses **Bayesian Information Criterion (BIC)** to automatically select
the number of isoform components (K) per gene:
> "The number of components is automatically selected using Bayesian
> Information Criterion (BIC)."

BIC penalizes model complexity against log-likelihood improvement, and
log-likelihood scales with the number of reads (i.e., sequencing depth) at
a given gene. This means: **for the same true underlying biology, a
higher-depth sample will, by the mathematical structure of BIC model
selection, more readily justify additional components (detect more pA
sites) than a lower-depth sample** — this is a direct, formal confirmation
(from a different tool's explicit methodology) of the same depth-confound
phenomenon independently observed in IsoDecipher's entropy analysis
(Section 8.1): DonorD's higher library size was associated with apparently
higher "complexity" (entropy) not solely attributable to underlying biology.

This is useful supporting evidence: the depth-confound is not a flaw
specific to IsoDecipher's entropy calculation — it is a structural property
of *any* method (including a fundamentally different, de novo, BIC-based
one like SCAPE) that lets sequencing depth influence how much "complexity"
(number of isoforms, or entropy magnitude) gets detected/inferred per gene
per cell. This supports framing the Section 8 entropy-depth problem as a
field-wide statistical phenomenon worth explicitly addressing, not an
IsoDecipher-specific weakness.

### 13.4 SCAPE's "Expected pA Length" formula — confirms the ordinal/
relative-position weighting pattern
$$\bar\theta = \sum_{k=1}^{K} \pi_k \cdot \frac{\theta_k - \theta_1}{\theta_K - \theta_1}$$

This normalizes every gene's pA site positions to a **relative [0,1] scale**
(proximal-most = 0, distal-most = 1) and computes a usage-weighted average
on that relative scale — **not on absolute physical bp distance**. This is
mechanistically the same category of design choice as scPAISO's ordinal
rank score (Section 9.6): both normalize away the actual physical spacing
between sites, treating a gene with sites 50bp apart and a gene with sites
5000bp apart identically once rescaled to [0,1].

This reinforces the Section 9.6 observation as a genuine, now twice-
independently-confirmed pattern: **the field's dominant convention is to
summarize APA usage via relative/ordinal positioning** (SCAPE's normalized
θ̄, scPAISO's ordinal rank score, and by extension DaPars/scDaPars's binary
proximal/distal PDUI, which is the degenerate K=2 case of the same idea).
IsoDecipher's `weighted_utr` — using actual physical UTR length in bp as
the weight — is a structurally different and more information-preserving
choice when site spacing is biologically uneven (which IsoDecipher's own
panel data directly demonstrates: e.g., ST6GAL1's 9 sites span from
~0bp to ~3132bp average UTR length in highly uneven increments; PRDM1's
4 sites span 25bp to ~2456bp). This trade-off — relative/ordinal
(simpler, scale-invariant, used by SCAPE/scPAISO/PDUI-family) vs.
physical-distance-weighted (more information-preserving but requires
reliable UTR-length annotation per site, as IsoDecipher's `avg_spliced_utr`
column provides) — is worth stating explicitly as a Methods-level design
decision with cited precedent on both sides, not framed as one approach
being simply "better."

### 13.5 Updated entry for Section 10 (per-sample fragment-size adaptation table)
Given 13.1-13.2, Section 10.4's SCAPE entry should be sharpened: SCAPE's
insert-size estimate IS recomputed from each input BAM's own cleavage
reads (so in a narrow sense it is "per-sample"), but the model family it
fits to that data (a single Gaussian, 2 parameters) cannot represent the
multi-modal, degraded-sample fragment-size behavior IsoDecipher's KDE
approach is built to handle, and is estimated only from a structurally
restricted subset of the genome (single-exon genes) rather than from the
panel/genes of direct interest. The distinction is not "per-sample vs.
not per-sample" but **"per-sample point-parameter (mean/std) re-estimation
within a fixed Gaussian model family" (SCAPE) vs. "per-sample full
non-parametric distributional re-estimation" (IsoDecipher)** — the latter
is strictly more flexible and is the property that matters for handling
the breast cancer dataset's bimodal short-fragment contamination correctly.

## 14. Corrections, Confirmations, and IsoCAPE/IsoFormer-Relevant Findings
*Added 2026-06-19, sixth pass. Includes a correction of an earlier error in
this document, confirmed publication details for PASTA, peak-caller
alternatives relevant to IsoCAPE, and the SCAP foundational paper.*

### 14.1 CORRECTION — PASTA does not use Read1-based detection
An earlier note in this document (Section 12, implicitly) conflated PASTA
with scPAISO's Read1-based cleavage-site detection. **This was incorrect.**
To be precise:
- **scPAISO** (Section 9.6) is the tool that uses Read1's cDNA content for
  near-nucleotide-resolution cleavage site detection.
- **PASTA** (Section 9.4) does NOT process raw reads/BAMs at all. It is a
  purely downstream statistical analysis package (R) that consumes a
  pre-computed polyA site count matrix — by default from **polyApipe**
  (Section 14.3), a separate, Python+R peak-calling pipeline. PASTA itself
  contains no read-level or insert-size logic.

### 14.2 PASTA — confirmed publication details
- **Journal**: Cell, Volume 187, Issue 16, pp. 4408-4425.e23 (2024)
- **DOI**: 10.1016/j.cell.2024.06.005
- **Title**: "Multiplexed single-cell characterization of alternative
  polyadenylation regulators"
- **Authors**: Kowalski MH*, Wessels HH*, Linder J*, Dalgarno C, Mascio I,
  Choudhary S, Hartman A, Hao Y, Kundaje A, Satija R (* = co-first authors)
- **Institution**: New York Genome Center / NYU / Stanford
- **Received Feb 9 2023; published online June 25 2024** (long review cycle,
  bioRxiv preprint from Feb 2023 — relevant if citing the earlier preprint
  version vs. final Cell version)
- **Core dataset**: introduces **CPA-Perturb-seq**, a multiplexed CRISPR
  perturbation screen of **42 cleavage-and-polyadenylation (CPA) regulator
  genes**, with a 3' scRNA-seq readout, performed in HEK293FT and K562 cell
  lines (data publicly released via Zenodo as Seurat .Rds objects).
- **Important clarification on the "background model"**: per the abstract
  and Figure 2 legend ("Average usage of 5,335 proximal polyA sites in
  **NT [non-targeting] cells**... vs. CSTF3-perturbed cells"), the
  background distribution in the PASTA manuscript's primary use case is
  computed from **non-targeting (NT) control cells within the same
  CRISPR screen** — i.e., this is a designed experimental control group,
  not merely "all cells in the dataset" by default. The PBMC vignette
  (Section 12) is a secondary demonstration where no NT-control structure
  exists, and there the background is computed from all cells in that
  dataset (as confirmed by the "Calculating background distribution / Using
  all cells" log message in the GitHub issue, Section 9.4/12). **This
  resolves the open question from Section 12**: PASTA's background model is
  flexible — it uses whatever reference group is specified/available
  (ideally a designed control, e.g., NT cells in a perturbation screen),
  defaulting to "all cells" when no such designed control exists. This
  reinforces the Section 12 conclusion: for a cancer-only dataset with no
  built-in normal/NT reference, PASTA's "all cells" fallback would compute
  a background from the tumor population itself, which would not detect
  tumor-vs-normal-level deviation — the same external-reference design
  (analogous to IsoCAPE's existing `delta_ce_external` healthy-BAM
  reference framework) remains necessary for a cancer diagnostic use case.

### 14.3 polyApipe — the actual peak-calling layer beneath PASTA
GitHub: MonashBioinformaticsPlatform/polyApipe.
- **Architecture**: Python 3 component (`polyApipe.py`) performs **peak
  calling and UMI counting** from 10x Genomics scRNA-seq BAMs; a companion
  R component (`polyApiper`) handles downstream analysis of the resulting
  UMI count matrix. This Python(upstream)+R(downstream) split mirrors
  exactly the kind of cross-language pipeline architecture under
  consideration for linking IsoDecipher/IsoCAPE (Python) outputs with
  R-based statistical tools (PASTA, etc.) — confirms this hybrid pattern
  is an established, working precedent in the field, not a novel risk.
- De novo peak caller, not annotation-anchored — same category as Sierra,
  scAPA, scAPAtrap (Section 1).

### 14.4 SCAP foundational paper — confirmed citation and direct relevance
**Cheng LC**, Zheng D, Baljinnyam E, Sun F, Ogami K, Yeung PL, Hoque M, Lu
C-W, Manley JL, Tian B. "Widespread transcript shortening through
alternative polyadenylation in secretory cell differentiation." *Nature
Communications* 11, 3182 (2020). DOI: 10.1038/s41467-020-16959-2.
PMID: 32576858.

- **Note on authorship**: first author is **Larry C. Cheng** (Tian lab,
  Rutgers), unrelated to Rene Yu-Hong Cheng (IsoDecipher) — same surname,
  different researcher.
- **Directly establishes the SCAP (secretion-coupled APA) concept**
  motivating IsoDecipher's biological framing: primary finding is in
  syncytiotrophoblast (placental) differentiation, but the paper explicitly
  states **"this mechanism, named secretion-coupled APA (SCAP), is also
  executed in B cell differentiation to plasma cells"** — i.e., the
  B-cell-to-plasma-cell SCAP phenomenon IsoDecipher's case study examines
  was already identified (if less exhaustively characterized at single-cell
  resolution) in this 2020 paper. This should be cited as the primary
  motivating/foundational reference for IsoDecipher's biological narrative,
  with IsoDecipher framed as providing single-cell-resolution, GTF-anchored,
  depth-robust quantification of a phenomenon whose existence in this
  specific cell system (B cell → plasma cell) was already established here.
- Key mechanistic claim from this paper directly relevant to IsoDecipher's
  Methods/Discussion: SCAP-driven shortening is explicitly **"unrelated to
  cell proliferation"** (a previously known APA-shortening driver) and
  instead **"accompanies increased secretory functions"**, with shortened
  3'UTRs conferring **higher mRNA stability** — useful mechanistic framing
  for IsoDecipher's own UTR-shortening-vs-differentiation findings.
- This paper is independently cited by the PASTA Cell 2024 paper (Section
  9.4/14.2) in the same biological context (plasma cell 3'UTR usage) —
  confirming both groups are working from the same foundational reference,
  positioning IsoDecipher's B-cell case study within an actively-cited,
  recognized area of the field rather than an isolated application.

### 14.5 Peak-caller alternatives for IsoCAPE (de novo detection layer)
Direct response to the question of whether a better-established peak
caller exists that IsoCAPE could adopt instead of (or alongside) a custom
implementation, given IsoCAPE's stated focus on short-read data requiring
sufficient per-sample depth (as opposed to a long-read approach):

| Tool | Peak-calling method | Notes for IsoCAPE consideration |
|---|---|---|
| **SCAPTURE** | HOMER (Poisson-based) peak calling + DeepPASS (CNN) validation; transcriptomic-annotation-aware; **performs peak identification across all cells collectively** (per Section 14.6's newly-found 2024 review) | Most directly comparable target for IsoCAPE's AluCE-type novel-site discovery use case, given DeepPASS's explicit purpose of distinguishing true PAS from false positives — relevant to IsoCAPE's internal-priming/AluCE disambiguation problem |
| **Sierra** | Gaussian curve fitting on coverage, splice-aware | R-based; per Bi et al. 2026, most depth-sensitive of all tools tested — likely unsuitable if IsoCAPE needs to remain usable at moderate depth |
| **scAPA** | HOMER (Poisson), **per-cell-type peak identification** (separate calling per cell type, unlike SCAPTURE's pooled approach) | The per-cell-type design could be directly useful if IsoCAPE wants cell-type-resolved novel site discovery rather than a single pooled call set |
| **polyApipe** | Custom Python peak caller + UMI counting | Already proven to interoperate with an R downstream package (PASTA) — useful precedent for cross-language design, less proven on accuracy specifically |
| **scAPAtrap** | "Genome-wide sensitive peak calling... can accurately locate pAs without using prior genome annotation, even for very low read coverage" (self-described) | Explicitly markets itself as more low-coverage-robust than coverage-fitting approaches — worth direct benchmarking against IsoCAPE's own approach if low-depth robustness is a priority |
| **LAPA** (Section 14.7, new) | Cluster-extension with "patience parameter," followed by Gaussian-kernel-smoothed peak picking | Designed for **long-read RNA-seq**, not directly applicable to IsoCAPE's stated short-read focus, but its clustering algorithm (extend while read-end counts exceed cutoff, terminate via patience parameter) is conceptually adjacent to IsoDecipher's own consecutive-gap clustering and could be worth examining purely as an algorithmic reference even if not adopted wholesale |

**Recommendation for IsoCAPE evaluation**: SCAPTURE (pooled, DeepPASS-validated)
and scAPAtrap (claimed low-coverage robustness) appear to be the two most
relevant existing tools to directly benchmark IsoCAPE's custom peak caller
against, given IsoCAPE's specific concerns (sufficient per-sample depth
requirement, AluCE/internal-priming disambiguation). Neither has been
deep-dived in this research session yet — flagged for follow-up when
returning to IsoCAPE development.

### 14.6 New finding — 2024 review confirms pooled-vs-per-cell-type peak
calling as a real design axis
"Guidelines for alternative polyadenylation identification tools using
single-cell and spatial transcriptomics data" (bioRxiv, Dec 2024) explicitly
contrasts: **"SCAPTURE performs identification across all cells
collectively, while scAPA identifies peaks separately for each cell
type."** This is a previously unrecorded design axis (pooled vs.
per-cell-type peak calling) relevant to IsoCAPE's design choices — worth
a full read of this review when IsoCAPE work resumes, as it appears to be
the most current (Dec 2024) and most granular methodological comparison
found across both research sessions.

### 14.7 New finding — LAPA (long-read-focused, noted for future reference
only, not directly applicable to IsoCAPE's short-read focus)
LAPA (bioRxiv 2022.11.08.515683): "a computational toolkit to study APA
from diverse data sources such as LR-RNA-seq [long-read] and 3'-seq...
generic enough to analyze any 3'-seq or long-read RNA-seq protocol."
Available via PyPI (Python). Notable capability not seen elsewhere:
**"correction of annotated transcript ends"** — i.e., LAPA can use
empirical peak-calling results to retroactively correct/refine the
reference GTF's annotated 3' end coordinates, rather than only trusting
the GTF as ground truth. This is explicitly noted as relevant only if
IsoCAPE's scope expands to long-read data in the future; **per the user's
clarification, IsoCAPE is currently focused on short-read sequencing
specifically because it requires sufficient per-sample depth for reliable
quantification** (a constraint long-read methods, with their typically
lower per-sample throughput, would not as easily satisfy at single-cell
scale) — so LAPA is recorded here as a reference algorithm/concept
(particularly its GTF-correction capability) rather than a tool under
active consideration for adoption.

### 14.8 IsoFormer — design concept clarification (recorded from user
discussion, not yet implemented)
IsoFormer's intended scope, as clarified: a (non-LLM) model relating APA
"grammar"/patterns to cell-type/cell-state identity, intended to work
downstream of and in combination with IsoCAPE (novel/de novo site
discovery, short-read-focused) and IsoDecipher (annotation-anchored
quantification). Explicitly NOT conceived as an LLM-based architecture.
No further design details discussed yet in this session — flagged for a
dedicated future design session once IsoCAPE work resumes, at which point
this document's tool survey (especially Section 9.6's scPAISO biological
findings — e.g., PA-isoform-level DE invisible at the gene level, directly
relevant to an APA-to-cell-identity mapping model) should be revisited as
relevant prior art / validation precedent.

## 15. Funding/Fellowship Positioning Notes (Activate Fellowship, NSF SBIR)
*Added 2026-06-19, seventh pass. Strategic notes for grant/fellowship
applications — not yet drafted into application text, recorded here as a
reference checklist for when that writing happens.*

### 15.1 Scope check — what today's competitive landscape actually covers
Nearly all tools surveyed across this document (scAPA, Sierra, scAPAtrap,
SCAPTURE, SAPAS, scDaPars, SCAPE, scDAPA, scUTRquant, PASTA/polyApipe,
scPAISO) are **single-cell 3' scRNA-seq tools** (10x Chromium-style protocols),
i.e., directly the same niche as IsoDecipher. Two edge cases are NOT direct
competitors: LAPA (long-read focused, Section 14.7) and MAAPER (originally
bulk RNA-seq, included in single-cell benchmarks but not natively
single-cell-designed). **Conclusion: the single-cell 3' APA quantification
space is genuinely crowded** (10+ academic tools, including groups at
Satija lab/NYGC, Mayr lab/MSKCC, and others) — an application centered
purely on "our quantification algorithm is more accurate" is a weak,
likely-contested position given this density.

### 15.2 Recommended reframing — from "better algorithm" to "missing
infrastructure for translation"
Given 15.1, the strategic recommendation is to shift the core funding
narrative away from algorithmic novelty/accuracy claims (a crowded,
hard-to-defend space) and toward **translational/infrastructure gaps**
that this session's research repeatedly confirmed are real and largely
unaddressed by the existing academic toolkit:

1. **Ecosystem accessibility (AnnData-native, Scanpy-integrated)**
   - Verified (Section 11): nearly the entire existing toolkit is
     R/Bioconductor/Seurat-native (scUTRquant's SCE, scDaPars' R matrices,
     PASTA's Seurat integration, scPAISO's Seurat-based downstream
     analysis). AnnData/Python support, where present at all, was added
     as an afterthought (scUTRquant v0.5.0).
   - This is a genuine, demonstrable engineering/usability gap for any
     user building modern ML/DL pipelines (which predominantly use
     Python/PyTorch/scverse), directly relevant to a translational/SBIR
     "can this be deployed and integrated quickly" evaluation criterion,
     as opposed to a purely academic novelty criterion.

2. **Per-sample, depth/quality-adaptive calibration — validated on a
   degraded clinical-type sample**
   - Verified (Section 10): of all tools surveyed, only IsoDecipher and
     the very recently published scPAISO (2025) re-derive their core
     distance model independently per input BAM/sample. scPAISO's version
     requires a non-universal sequencing configuration (long Read1);
     IsoDecipher's does not.
   - IsoDecipher's per-sample KDE approach has been empirically validated
     on a real degraded/short-fragment-enriched sample (FACS-sorted breast
     cancer, Figure 2B) — directly relevant to the kind of variable-quality
     clinical sample heterogeneity a diagnostic/translational reviewer
     would expect to be a major practical obstacle.
   - This is a concrete "robustness to real-world sample variability"
     claim, mechanistically explained (Section 13.2: SCAPE's single-Gaussian
     model family cannot represent the bimodal fragment-size distortion
     this sample exhibits; IsoDecipher's non-parametric KDE can) — a much
     stronger, falsifiable technical claim than a vague accuracy comparison.

3. **Demonstrated rigor on a specific, non-obvious methodological pitfall**
   (IG_C_gene/TR_C_gene biotype misclassification + quantified non-coding
   RNA contamination of weighted UTR, R=-0.624 with differentiation state)
   - This level of detail-oriented validation (catching and quantifying a
     subtle annotation/biotype bug that silently corrupts a commonly-used
     summary statistic) is a credibility signal especially relevant to
     reviewers evaluating rigor for eventual clinical/regulatory use —
     this kind of issue is exactly the type of subtle confound a
     diagnostic-grade tool needs to have caught and documented.

4. **Long-term architecture story: IsoCAPE + IsoDecipher + IsoFormer as a
   three-layer pipeline**
   - No tool surveyed in this document combines (a) de novo/novel site
     discovery (IsoCAPE's role), (b) annotation-anchored, depth-robust,
     biotype-aware quantification (IsoDecipher's current role), and (c) a
     model relating APA "grammar"/patterns to cell identity/state
     (IsoFormer's proposed, not-yet-LLM-based role) in one coherent,
     Python-native pipeline. This three-layer story is well-suited to a
     "Phase II / scaling potential" narrative (a key SBIR evaluation axis)
     even though IsoFormer itself is not yet implemented — it can be
     presented as the natural extension of validated current capability
     (IsoDecipher) plus an already-partially-built discovery layer
     (IsoCAPE), not pure speculation.

### 15.3 What NOT to claim (avoid overreach, per Section 12's findings)
- Do NOT claim other tools lack real biological validation — Section 12
  found this is largely false (SCAPE, scDaPars, PASTA, scPAISO, and
  scUTRquant all have substantial, real biological-discovery applications
  using their own pipelines). This claim is easily falsifiable by a
  knowledgeable reviewer and would damage credibility.
- Do NOT claim outright superiority on quantification accuracy without a
  head-to-head benchmark (not yet performed) — the per-sample-calibration
  advantage (15.2.2) is currently a *mechanistically well-explained,
  theoretically grounded* claim, not yet a directly-measured comparative
  accuracy claim against another specific tool on the same dataset.
- Frame all comparative claims at the level of "structural design choice"
  (annotation-anchored vs. de novo; fixed vs. per-sample-calibrated
  parameters; N-point vs. 2-point group structure; R/Seurat-native vs.
  AnnData-native) rather than "IsoDecipher is better" — this is both more
  defensible and maps more directly onto concrete review criteria
  (engineering reproducibility, translational readiness, robustness to
  real sample heterogeneity) than a generic superiority claim would.

### 15.4 Action items for application drafting (not yet done)
- Draft a "translational gap analysis" paragraph using 15.2's four points
  as the core structure, citing the specific verified findings (and tool
  names/papers) from Sections 9-14 of this document for each claim.
- Consider explicitly citing the Bi et al. 2026 NAR benchmark and the
  Fahmi et al. 2025 EMM review (Section 9.1-9.2) as evidence that this is
  a recognized, actively-surveyed field — useful for establishing that the
  application is grounded in current literature awareness, not operating
  in isolation.
- The SCAP foundational citation (Cheng LC et al. 2020 Nat Commun, Section
  14.4) should anchor the biological motivation section of any application,
  independent of the tool-comparison sections.

## 16. Panel Expansion Strategy — PolyASite Versions, Cryptic Site Databases,
## and the IsoCAPE-to-IsoDecipher Feedback Loop
*Added 2026-06-19, eighth pass. IMPORTANT CONTEXT CORRECTION: the panel is
no longer a curated/focused gene list (the earlier ~367-gene B-cell panel
has been abandoned) — IsoDecipher's panel construction is now intended to
be genome-wide (any expressed gene is included) to support global/
unsupervised analysis. All panel-expansion discussion below assumes a
genome-wide panel, not a focused gene list.*

### 16.1 No general-purpose "cryptic exon/intronic polyA database" exists
Searched specifically for a curated, downloadable, version-controlled
database of cryptic exons / intronic polyadenylation sites analogous to
PolyASite for canonical sites. **Conclusion: no such general-purpose
resource was found.** The cryptic-APA knowledge that exists in the
literature is fragmented across disease-specific compendia (e.g., the
TDP-43/ALS-FTD cryptic APA Supplementary Tables from Bryce-Smith et al.
2025, Section 16.3) — these are one-off supplementary materials tied to a
specific perturbation/disease context, not a maintained, programmatically
queryable, cross-disease atlas. **Implication**: there is no shortcut —
expanding IsoDecipher's panel to include non-GTF-annotated/cryptic sites
requires either (a) systematic discovery via IsoCAPE's own de novo pipeline,
or (b) manual curation of individual well-characterized disease-specific
sites from the literature (Section 16.4). There is no single database to
simply download and merge in.

### 16.2 PolyASite version comparison — which to use for what
Two distinct, both-still-relevant PolyASite versions now confirmed:

**PolyASite v2.0** (Herrmann et al. 2020, NAR; already used in IsoDecipher's
Figure 2B tolerance validation):
- Built from **bulk 3' end sequencing** data (>12 distinct protocols,
  primarily cell lines), human/mouse/worm, GRCh38/Ensembl 96.
- Fully automated, containerized pipeline (zavolanlab/polyAsite_workflow,
  GitHub) — reproducible, versioned.
- Appropriate for: high-confidence "ground truth" validation use cases
  (as IsoDecipher already uses it) — clean protocol, well-characterized
  internal-priming flagging.

**PolyASite v3.0** (Moon, Herrmann, Mironov, Zavolan 2025, NAR Database
Issue 53(D1):D197-D204, DOI 10.1093/nar/gkae1043 — received Sept 14 2024,
published Jan 6 2025):
- Built directly from **single-cell RNA-seq data** (10x Genomics V2/V3 3'
  kits) using the SCINPAS workflow, aggregated to **tissue-level** (not
  cell-level — authors explicitly state "the resolution of scRNA-seq was
  too low for robust detection of cell-level differences in PAS usage").
- Data sources: Human Cell Atlas, NEMO brain repository (human); Tabula
  Muris Senis (mouse); C. elegans embryonic 10x data.
- **Directly methodologically aligned with IsoDecipher's own input data
  type** (10x Chromium 3' scRNA-seq) — same library prep/capture biology
  as IsoDecipher's own BAMs, unlike v2.0's heterogeneous bulk protocols.
- Provides **three stringency levels** (20% / 62% / 87% canonical
  polyadenylation-signal-motif presence) — allows the user to trade off
  sensitivity vs. specificity explicitly, rather than using a single fixed
  confidence threshold.
- **Validation against v2.0**: at 62% stringency, ~75% of v3.0 **terminal
  exon** sites perfectly match v2.0 sites within a matching-distance
  threshold; **intronic and intergenic sites show substantially more
  discrepancy** between v2.0 and v3.0. This is an important caution: any
  IPA-type (intronic) site sourced from v3.0 should be treated with lower
  confidence than a terminal-exon-type site, consistent with IPA's general
  reputation across this whole document as the hardest-to-reliably-call
  APA category (Section 9.6, scPAISO findings; Section 16.3, TDP-43 paper's
  explicit statement that IPA events were the hardest category to validate).
- **PAS-to-gene assignment rule** (useful as a reference design choice for
  IsoDecipher's own panel logic): when a PAS overlaps multiple genes,
  priority is (1) protein-coding gene, (2) lncRNA gene, then tie-break by
  (i) closest 3' end border, (ii) longest overlapping gene, (iii) random.

**Recommendation**: keep v2.0 as the tolerance-validation reference
(Figure 2B, unchanged); use v3.0 (filtered to a conservative stringency,
e.g. 62% or 87%) as the primary source for genome-wide panel expansion,
given its direct technical alignment with IsoDecipher's own 10x 3'
scRNA-seq data type. Treat any v3.0-sourced intronic/IPA-type sites with
extra caution given the documented v2.0/v3.0 discrepancy in that category.

### 16.3 TDP-43/ALS-FTD cryptic APA compendium — not directly usable, but
provides the operational classification framework and a concrete novel-site
validation rule
Full review of Bryce-Smith et al. 2025, Nat Neurosci 28:2190-2200 (DOI
10.1038/s41593-025-02050-w) confirms: this is a **disease/perturbation-
specific** compendium (227 cryptic APA events induced specifically by
TDP-43 depletion in iPSC-derived neurons), not a general reference —
directly merging its site list into a B-cell/cancer-focused panel would
not be appropriate without independent justification for each site's
relevance to the target biology.

**What IS directly reusable from this paper**:
1. **Operational three-category classification** (ALE / 3'Ext / IPA),
   precisely defined:
   - **ALE** (Alternative Last Exon): defined by an **upstream novel splice
     junction** that creates a new terminal exon.
   - **3'Ext** (3'UTR Extension): **independent of splicing** — occurs
     downstream of an annotated distal 3'UTR, extending it. This is the
     category structurally closest to IsoDecipher's existing G0/G1/G2
     framework (a new, more-distal site beyond the previously-known
     furthest site).
   - **IPA** (Intronic Polyadenylation): occurs within an intron, **in the
     absence of upstream alternative splicing** — distinguishes IPA from
     ALE specifically by the absence of a novel splice junction. (Note:
     this refines the simpler tPA/iPA/ePA scheme from scPAISO, Section 9.6 —
     the TDP-43 paper's ALE/3'Ext/IPA scheme is more mechanistically
     precise and should be preferred as the basis for any `site_type`
     metadata column added to IsoDecipher's panel.)
2. **Concrete, reusable novel-site validation rule**: a putative novel
   last exon is accepted as credible if its PAS is either **(a) within
   100nt of a PolyASite-annotated PAS**, or **(b) contains a conserved
   poly(A) signal hexamer** in the terminal 100nt. This is a directly
   implementable filter for any IsoCAPE-discovered candidate site before
   it is considered for inclusion in IsoDecipher's panel (Section 16.5).
3. Cryptic APA events were defined quantitatively as: **<10% mean usage in
   controls AND >10% usage change after perturbation** — a reusable
   pattern for defining "this site is biologically interesting/disease-
   associated" in a differential context, separate from the question of
   whether the site exists/should be in the panel at all.
4. Confirms (consistent with this document's repeated finding) that **IPA
   events are the hardest category to reliably detect and validate** — in
   postmortem bulk tissue specifically, the paper notes "IPA detection is
   further complicated by the fact that normal pre-mRNA reads also map to
   IPA regions, creating significant noise" — directly relevant caution
   for any IsoDecipher/IsoCAPE work expanding into IPA-type site detection.

### 16.4 AR-V7 (androgen receptor splice variant 7) — mechanism clarification
and as a concrete worked example
**Important mechanistic correction**: AR-V7 is **not** a simple 3'UTR-
length/proximal-distal APA event (i.e., not directly analogous to
IsoDecipher's existing G0/G1/G2 weighted-UTR framework). The actual
mechanism:
- AR-V7 transcript = normal exons 1-3 (encoding the N-terminal domain and
  DNA-binding domain), **spliced via a novel splice junction to a "cryptic
  exon 3" (CE3) located within intron 3** of the AR gene.
- This is, in the TDP-43 paper's terminology (Section 16.3), structurally
  an **ALE event**: CE3 is an alternative last exon defined by an upstream
  novel splice junction, not a simple distal/proximal polyA choice within
  an existing terminal exon.
- Splicing into CE3 causes premature translation termination, producing
  a truncated, constitutively-active AR protein lacking the ligand-binding
  domain (16 novel C-terminal amino acids replace the normal LBD) — this
  is the well-established mechanism driving castration-resistant prostate
  cancer (CRPC) and resistance to AR-targeted therapies (abiraterone,
  enzalutamide). AR-V7 is an FDA-relevant, clinically-used liquid-biopsy
  biomarker (multiple cited clinical studies, e.g. PMC6307949, Clin Cancer
  Res 2022).
- The CE3 polyA site itself is precisely characterized in the literature
  (multiple independent papers using 3'RACE, IGV read-pileup inspection
  distinguishing true Ex3-CE3 splice junction reads from background
  intronic reads, etc.) — i.e., this is a **well-validated, high-confidence,
  disease-specific cryptic ALE site with a known genomic coordinate**, not
  a speculative or low-confidence candidate.

**Why this matters for panel design**: AR-V7/CE3 is a strong worked
example of "Line B" sites (Section 16.5) — a site that will never appear
in a healthy-tissue atlas like PolyASite (because it's essentially absent
in untreated primary prostate cancer, <1% expression, only emerging under
treatment-induced selective pressure — Section search result PMC6307949),
but is extremely well-characterized and clinically important. This is
exactly the kind of site that justifies a manually-curated, disease-
specific supplementary site list, separate from genome-wide healthy-tissue
atlas-derived expansion.

### 16.5 Recommended panel expansion architecture — two separate tracks
Given 16.1-16.4, the panel expansion strategy should explicitly separate
two tracks rather than merging everything into one undifferentiated
process:

**Track A — Genome-wide, healthy-tissue-baseline site expansion**
- Source: PolyASite v3.0 (Section 16.2), filtered to a conservative
  stringency (e.g., 62% or 87% canonical-motif presence).
- Goal: improve baseline panel completeness/coverage genome-wide, catching
  real sites that GTF transcript-end annotation may have missed, using a
  resource that is technically well-matched to IsoDecipher's own 10x 3'
  scRNA-seq data type.
- Treat IPA/intronic-type v3.0 sites with extra caution (lower confidence
  category per the v2.0/v3.0 discrepancy noted in 16.2).
- This track does NOT depend on IsoCAPE; it can be implemented directly
  against the existing GTF-anchored panel construction pipeline.

**Track B — Disease/biology-specific curated cryptic site supplement**
- Source: (i) individual well-characterized literature sites (e.g., AR-V7/
  CE3, and any future analogous biomarkers as they are identified), each
  manually verified and documented with a PMID/coordinate/evidence source;
  (ii) IsoCAPE's own de novo discoveries, once IsoCAPE work resumes,
  **filtered through the TDP-43 paper's novel-site validation rule**
  (within 100nt of a PolyASite PAS, or containing a conserved polyA signal
  hexamer — Section 16.3, point 2) before being promoted into this list.
- Maintained as a **separate CSV/table** from the genome-wide GTF-anchored
  panel (e.g., `known_disease_sites.csv`), using a schema compatible with
  but distinct from the main panel: gene, coordinate, `site_type` (using
  the ALE/3'Ext/IPA classification from 16.3), `evidence_source` (PMID or
  "IsoCAPE-discovered, validated [date]"), and any relevant disease/
  biomarker context fields.
- **This is the concrete mechanism for the IsoCAPE→IsoDecipher feedback
  loop the user proposed**: IsoCAPE's future de novo discoveries get
  validated against the 100nt/hexamer rule, then appended to this Track B
  table, which is unioned with the Track A/GTF-anchored panel at
  `integrate_samples.py` time (or at panel-build time) to produce the
  final, expanding panel used for quantification. Track B sites accumulate
  over time as IsoCAPE is applied to more datasets — i.e., IsoDecipher's
  panel becomes incrementally more complete with each IsoCAPE discovery
  cycle, without requiring IsoCAPE to be re-run on every new IsoDecipher
  analysis.
