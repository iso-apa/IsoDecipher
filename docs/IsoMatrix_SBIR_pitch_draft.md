# IsoMatrix Suite — SBIR / Fellowship Application Draft
*Working document — compiled 2026-06-22*
*Sources: APA_tool_landscape_reference.md (today's research session)*
*Adaptable for: NIH SBIR PA-24-185 / NSF SBIR BT5 / Activate Fellowship Cohort 2027*

---

## SECTION 1: THE PROBLEM (Significance)

### 1.1 Alternative Polyadenylation is an underexplored layer of gene regulation
with direct clinical relevance

The vast majority of human genes produce multiple mRNA isoforms through
alternative polyadenylation (APA) — the selection of different cleavage and
polyadenylation sites within the 3' end of a pre-mRNA transcript. APA
determines the length of the 3' untranslated region (UTR), which in turn
controls mRNA stability, translation efficiency, and subcellular localization.
APA dysregulation is now recognized as a hallmark of cancer and immune
dysfunction:

- Global 3'UTR shortening in secretory cells (including B cell-to-plasma cell
  differentiation) is mechanistically linked to increased secretory output
  (Cheng LC et al. 2020, Nat Commun — the SCAP mechanism).
- Disease-specific APA events produce protein isoforms with fundamentally
  altered function: AR-V7 (androgen receptor splice variant 7), an APA/ALE
  event generating a constitutively active truncated AR protein, drives
  castration-resistant prostate cancer (CRPC) and is an FDA-relevant
  liquid-biopsy biomarker.
- Intronic polyadenylation (IPA) events regulated by TDP-43 and other
  RNA-binding proteins are directly implicated in ALS/FTD neurodegeneration.

Despite this, APA remains systematically understudied in clinical genomics.
The reason is not scientific — it is infrastructural.

### 1.2 The existing single-cell APA toolkit has critical gaps that prevent
clinical translation

Eight single-cell APA tools have been published (2019-2022, per Fahmi et al.
2025 review); none has reached clinical deployment. The key barriers are:

**Gap 1 — Ecosystem incompatibility**
All eight published tools (scAPA, Sierra, scAPAtrap, SCAPTURE, SAPAS,
scDaPars, scDAPA, SCAPE) are R/Bioconductor/Seurat-native. Modern
computational biology and machine learning pipelines are predominantly
Python/AnnData/Scanpy-based. No tool in the reviewed literature provides
native AnnData output integrating GEX, ADT (protein), and APA isoform data
in a single object ready for downstream ML/DL workflows.

**Gap 2 — Fixed parameters fail on real clinical samples**
Of all tools reviewed, only IsoDecipher and the very recently published
scPAISO (2025) re-derive their distance model per input BAM. All other tools
use parameters fixed at reference-index build time (e.g. scUTRquant's
epsilon=30bp, scDAPA's 100bp bins). IsoDecipher has empirically demonstrated
per-sample KDE recalibration on a degraded FACS-sorted breast cancer sample
(bimodal insert-size distribution incompatible with any single-Gaussian model,
as used by SCAPE). This robustness to real-world sample quality variation is
a prerequisite for clinical deployment.

**Gap 3 — Biotype annotation errors silently corrupt summary statistics**
(a) Immunoglobulin constant region genes (IGHM, IGHG1-4, IGHA1-2, IGHE) are
tagged in Ensembl/GENCODE GTF with biotype "IG_C_gene" rather than
"protein_coding," causing all reviewed tools to treat these genes as
non-coding — silently zeroing their UTR length annotation and excluding them
from APA analysis. This is particularly critical for B cell and plasma cell
biology where IGH expression is the dominant transcriptional feature.
(b) Non-coding RNA (e.g. MALAT1, NEAT1) produces 3'-end-enriched reads in
10x Chromium data (observed but not quantified in SCAPTURE 2021). When
included in weighted UTR calculations without biotype filtering, their
increasing expression fraction during terminal differentiation (Spearman
R=-0.624 with differentiation score across 16,116 cells) artificially inflates
apparent UTR shortening — producing a spurious biological signal. No published
tool was found to explicitly quantify or correct for this confound.

**Gap 4 — No clinical-grade on-premise deployment pathway**
Hospital and pharmaceutical data governance requirements (HIPAA, GxP, EU GDPR)
prevent uploading patient genomic data to public cloud infrastructure. No
existing APA tool provides an enterprise-grade, containerized, on-premise
deployment architecture suitable for clinical or pharmaceutical use.

**Gap 5 — No end-to-end learning architecture connecting APA patterns to
cell identity and disease state**
Existing tools output isoform count matrices. None provides a trained model
that maps these APA patterns to cell state, disease risk, or drug target
relevance — the step required to convert a quantification tool into a clinical
decision-support product.

**Gap 6 — Static panels cannot capture disease-specific APA events**
All annotation-based tools (scUTRquant, MAAPER) rely on fixed reference panels
derived from GTF transcript annotations and/or healthy-tissue bulk 3' end
sequencing databases (e.g. PolyASite v2.0, 221 libraries from normal tissues).
These panels systematically miss disease-specific cleavage sites that are:
(a) absent from healthy-tissue atlases — e.g. AR-V7/CE3 in castration-resistant
    prostate cancer, present at <1% frequency in normal tissue and emerging only
    under treatment-induced selective pressure;
(b) intronic or alternative-last-exon (ALE) events not captured by standard
    GTF transcript 3' end annotation;
(c) novel sites discovered only through de novo analysis of disease samples.
No existing tool provides a mechanism for users to incorporate curated
disease-specific sites into their quantification panel, or to systematically
expand the panel as new discoveries are made. This means that any tool with
a fixed panel is structurally incapable of detecting the most clinically
actionable APA events — those unique to disease states.

---

## SECTION 2: THE SOLUTION (Innovation)

### 2.1 IsoMatrix — a three-layer, end-to-end APA intelligence platform

IsoMatrix is a continuously evolving suite of three complementary tools,
each addressing a distinct layer of the APA analysis problem:

```
LAYER 1: IsoCAPE (Discovery)
  Input:  BAM files (10x Chromium 3' scRNA-seq, standard short-read)
  Task:   De novo identification of cleavage/polyadenylation sites,
          including cryptic/disease-specific sites absent from GTF
          annotation (AluCE, AR-V7/CE3-type ALE events, IPA sites)
  Output: Site coordinates + ALE/3'Ext/IPA classification +
          confidence score (100nt/PolyASite proximity or hexamer rule)
          → feeds Track B panel expansion in IsoDecipher

LAYER 2: IsoDecipher (Quantification)
  Input:  BAM files + expanding panel
          (GTF-anchored [Track A] + IsoCAPE-validated cryptic sites [Track B])
  Task:   Per-cell, per-BAM-calibrated probabilistic read assignment →
          isoform count matrix with full metadata
  Output: AnnData (.h5ad) — GEX + ADT + Isoform in one object
          iso.var: gene, G0/G1/G2... group index, avg_spliced_utr,
          utr_source (protein_coding/non_coding/mixed),
          site_type (ALE/3'Ext/IPA), coord_spread, IG_C_gene-aware

LAYER 3: IsoFormer (Semantic / Pattern Recognition)
  Input:  IsoDecipher AnnData (per-cell isoform usage vectors)
  Task:   Learn "APA grammar" — which isoform usage patterns
          co-occur, characterize specific cell states, or deviate
          from a healthy reference distribution
  Output: APA-based cell state embedding (C_state) +
          Anomaly score (deviation from healthy reference) +
          IsoScore (single-number clinical summary)
          → attention weights feed back to IsoCAPE as
            "priority regions for next-sample discovery"
```

### 2.2 The IsoCAPE → IsoDecipher → IsoFormer feedback loop

The three layers form a self-improving system:

1. IsoCAPE discovers novel sites in new datasets
2. Sites passing the validation threshold (within 100nt of PolyASite v3.0
   PAS OR containing a conserved polyA signal hexamer) are promoted to
   IsoDecipher's Track B panel
3. IsoDecipher quantifies these sites across all cells, including the newly
   discovered ones
4. IsoFormer learns which of these sites — old and new — carry the most
   discriminative signal for cell state and disease classification
5. IsoFormer's learned feature importance is used as a prior to guide
   IsoCAPE's site-discovery sensitivity in the next dataset: high-weight
   genes receive more intensive investigation; low-weight genes can be
   processed more efficiently

This loop means each new dataset makes the entire platform more accurate,
more sensitive to clinically-relevant sites, and more computationally
efficient — a compounding return on data investment not achievable with
static, one-pass pipelines.

### 2.3 Key technical differentiators (verified against published literature)

| Property | IsoMatrix | Closest competing tool | Difference |
|---|---|---|---|
| Output ecosystem | AnnData-native, GEX+ADT+ISO in one object | scUTRquant (SCE-native, AnnData added v0.5.0) | IsoMatrix designed AnnData-first; others added Python support as afterthought |
| Insert-size calibration | Per-BAM non-parametric KDE, validated on degraded clinical sample | SCAPE (single Gaussian μ/σ per dataset) | KDE can represent bimodal/skewed distributions; single Gaussian cannot — directly explains SCAPE's failure on CEL-seq/Drop-seq (Bi et al. 2026) |
| Isoform group structure | N-point (G0/G1/G2...), arbitrary # sites per gene | scDaPars/DaPars/PDUI family: 2-point only (proximal/distal) | IsoMatrix can represent ST6GAL1 (9 groups), DDX5 (7 groups); PDUI cannot |
| Biotype handling | CODING_BIOTYPES set: includes IG_C_gene/TR_C_gene/IG_V_gene etc. | No tool reviewed handles this | IGH genes correctly classified as protein-coding; UTR length correctly computed |
| Non-coding RNA filter | Explicit utr_source filter; non-coding RNA fraction tracked as independent metric | SCAPTURE noted phenomenon but did not quantify downstream effect | R=-0.624 confound quantified; protein-coding-only weighted UTR used for cross-sample comparisons |
| Site type annotation | ALE / 3'Ext / IPA (per Bryce-Smith et al. 2025 framework) | scPAISO: tPA/iPA/ePA (similar but less mechanistically precise) | ALE/3'Ext/IPA maps directly to disease literature (TDP-43, AR-V7) |
| Deployment | On-premise Docker container (HIPAA-compliant, no data egress) | None found in literature | Clinical/pharma deployable without data governance concerns |
| Panel architecture | Versioned, multi-source, expandable: GTF + PolyASite v2.0 + disease-specific tracks + IsoCAPE feedback loop | All tools use fixed panels (scUTRquant: precompiled UTRome index; MAAPER: fixed PolyASite reference) | Users can incorporate known disease sites; panel grows automatically with each IsoCAPE discovery cycle |

### 2.4 IsoFormer — architecture roadmap

**Phase I (current proposal scope): APA-VAE**
A variational autoencoder trained on IsoDecipher's per-cell isoform usage
vectors. Analogous to scVI for gene expression, but APA-specific. Produces
a low-dimensional C_state embedding that captures APA-based cell identity
independent of GEX. Computable on existing hardware (GPU workstation);
immediate application to B cell differentiation and cancer datasets.

**Phase II (18-24 months): Isoform-level Transformer**
Input tokenization: each gene's multi-isoform usage vector (G0/G1/G2...)
as a single "token" — APA grammar = co-usage patterns across genes.
Enables: (a) cross-dataset transfer learning, (b) perturbation prediction
("if you knock out RBP X, which isoform patterns change?"),
(c) systematic comparison against a healthy reference distribution
(the prerequisite for an anomaly-based IsoScore).
Requires: multi-dataset training corpus, ~100-200 donor-level samples;
GPU cluster (target: Taiwan NCHC/TWCC collaboration for training compute).

**Phase III (3+ years): APA Foundation Model**
Pretrained across tissues, diseases, and species on IsoDecipher-quantified
APA profiles. Analogous to scGPT/Geneformer but APA-specific. Enables
zero-shot inference for novel cell types or diseases. Requires: large-scale
data accumulation from IsoCAPE+IsoDecipher pipeline runs; AI inference
chip optimization for on-premise deployment (Taiwan AI chip partnership
target: inference acceleration for on-premise IsoMatrix appliance).

---

## SECTION 3: COMMERCIAL OPPORTUNITY (Significance / Commercialization)

### 3.1 Target markets

**Primary (Year 1-2): Pharmaceutical / Biotech Translational Bioinformatics**
- Customers: Computational biology teams at mid-to-large pharma
  (Seattle area: Pfizer/Seagen, BMS; remote: Genentech, AstraZeneca,
  Regeneron, Novo Nordisk)
- Pain point: Single-cell datasets exist; APA layer is unanalyzed due to
  lack of validated, Python-native, data-governance-compliant tooling
- Entry point: B2B consulting contract → software license
- Pricing model: Annual on-premise Docker license (per-site or per-dataset)
  + consulting/analysis services
- Revenue: $50K-$200K per engagement (comparable to existing single-cell
  analysis consulting rates)

**Secondary (Year 2-3): CRO (Contract Research Organizations)**
- Customers: Charles River, Covance, WuXi, Pacific Biosciences service arms
- Pain point: Need to offer APA analysis as a premium service to pharma clients
- Entry point: Become the standard APA analysis vendor for their single-cell
  service offering
- Advantage: One CRO contract → access to multiple pharma clients

**Tertiary (Year 3+): Clinical / Diagnostic Applications**
- Customers: Hospital genomics cores (Fred Hutch, UW Medicine, SCCA in Seattle)
- Pain point: APA biomarkers (AR-V7, IGHM/IGHA switching) have clinical
  relevance but no validated quantification tool approved for on-premise use
- Entry point: Research use only (RUO) → IVD pathway (LDT first)
- Note: Regulatory pathway (FDA LDT/IVD) is a long-term consideration;
  on-premise Docker architecture was explicitly chosen to minimize data
  governance complexity in early commercial deployments

### 3.2 Taiwan strategic partnership (hardware acceleration layer)

IsoFormer Phase III requires AI inference chip optimization for on-premise
deployment at clinical scale. This creates a concrete, technically-grounded
collaboration opportunity with Taiwan's semiconductor ecosystem:

- Algorithm layer (IsoMatrix): developed and owned by IsoMatrix LLC (Seattle)
- Silicon layer (inference acceleration): target collaboration with Taiwan
  AI chip design companies (MediaTek, Faraday Technology, SiGe Semiconductor)
  and ITRI (Industrial Technology Research Institute) International Cooperation
- Compute layer (foundation model training): Taiwan NCHC/TWCC GPU cluster
  (National Center for High-Performance Computing / Taiwan Computing Cloud)
  — publicly accessible to researchers, with priority allocation for
  AI-in-health projects

This hardware partnership is not speculative: it is the natural convergence
of IsoMatrix's algorithmic core (developed in isolation of any hardware
dependency) with Taiwan's established AI chip design capability, motivated
by the specific requirement for on-premise clinical deployment at a cost
point acceptable to hospitals and smaller pharmaceutical sites. The
partnership model follows the precedent of the NSF-Taiwan Deep Tech
Innovation and Partnerships Workshop (NSF Award #2414276, 2024), which
explicitly targeted US-Taiwan collaboration in biotechnology and AI.

### 3.2b Disease-Specific Panel Tracks — the recurring revenue asset

IsoMatrix's panel architecture (Section 2, differentiator table) separates the
baseline quantification engine (IsoDecipher, open-source AGPLv3) from
disease-specific panel data (proprietary, subscription-based). This creates a
compounding data asset:

**Track A — Healthy baseline panel (open, included with IsoDecipher):**
- GTF-derived sites (Ensembl GRCh38, ~52,000 features across ~15,000 genes)
- PolyASite v2.0 validated extensions (TE/IN/DS sites, fraction>10% of 221
  independent bulk 3' end sequencing libraries)
- Sufficient for standard single-cell APA analysis in non-disease contexts

**Track B — Disease-specific panel tracks (proprietary, subscription):**
Pre-built, continuously-updated cleavage site panels for specific cancer and
disease contexts, distributed as versioned CSV files included in the commercial
Docker image. Each track is derived from IsoCAPE's de novo discovery pipeline
applied to publicly available or partnered disease-specific scRNA-seq datasets,
validated against PolyASite v2.0 coordinates and polyA signal hexamers:

| Track | Key sites | Data source | Clinical relevance |
|---|---|---|---|
| Prostate cancer / CRPC | AR-V7/CE3 (ALE event, intron 3 of AR) | Literature-curated | FDA liquid biopsy biomarker; resistance to enzalutamide/abiraterone |
| Multiple myeloma | IsoCAPE-discovered from GSE193531, GSE223060/61 | scRNA-seq (3' 10x) | Plasma cell identity escape; anti-BCMA therapy resistance |
| Breast cancer | IsoCAPE-discovered from FACS-sorted samples | Proprietary + public | Triple-negative subtype; UTR shortening in aggressive disease |
| ALS/FTD | TDP-43-regulated IPA/ALE events | Bryce-Smith et al. 2025 | Neurodegeneration biomarker panel |
| B cell lymphoma | IsoCAPE-discovered (planned) | Public scRNA-seq | DLBCL subtyping; treatment response prediction |

**Subscription model:**
- Annual license includes: disease track updates + new tracks added each year
- Each new IsoCAPE discovery cycle (processing a new cancer dataset) expands
  the relevant disease track, which is distributed to subscribers as a panel
  update within their existing Docker deployment
- No data egress: the panel CSV is included in the Docker image; subscriber's
  patient data never leaves their environment

**Why this creates a durable moat:**
Unlike algorithm improvements (which competitors can replicate), validated
disease-specific APA site panels accumulate over time and require real patient
data, experimental validation, and expert curation to produce. Each new
panel track added increases the value of the subscription without requiring
subscribers to change their analysis workflow — they simply pull the latest
Docker image.

### 3.3 Why on-premise Docker (not cloud SaaS)

This is a deliberate strategic and technical choice, not a limitation:

1. HIPAA / GxP compliance: Patient genomic data cannot be uploaded to
   public cloud infrastructure under most hospital and pharmaceutical
   data governance policies. On-premise deployment eliminates this barrier.
2. Founder preference and focus: IsoMatrix's founding team has deep
   algorithm expertise, not clinical data operations experience. On-premise
   architecture means IsoMatrix never handles patient data — the algorithm
   runs inside the customer's secure environment. This keeps IsoMatrix
   focused on what it does uniquely well (APA algorithm development) and
   avoids building a data operations infrastructure that is neither a
   core competency nor a defensible moat.
3. Compounding moat: Each Docker release includes the latest IsoFormer
   weights and the latest IsoCAPE-discovered site panel. Customers on
   annual licenses automatically benefit from the IsoCAPE→IsoDecipher→
   IsoFormer feedback loop — the platform improves with every new dataset
   processed anywhere in the network, and customers receive those
   improvements as Docker image updates, not as data sharing obligations.

---

## SECTION 4: PRELIMINARY DATA (Track Record / Feasibility)

### 4.1 IsoDecipher — validated, running pipeline

- Successfully processed 4-donor, 6-batch B cell differentiation dataset
  (GSE212138 + GSE229042, 21,670 cells, 10x Chromium 3')
- Panel: genome-wide GTF-anchored, ~52,105 isoform features post-filtering
- Per-sample KDE calibration validated across 6 independent datasets
  (including FACS-sorted breast cancer with pathological short-fragment
  insert-size distribution), using Bhattacharyya Coefficient (BC>0.6)
  as the empirical tolerance-selection criterion
- Key biological finding: APA subclustering of GEX-defined plasma cell
  cluster resolves 5 distinct subpopulations (including interferon-stimulated
  plasma cells, IFIT1/IFIT2/IFITM2) invisible to gene expression alone
- Key methodological finding: non-coding RNA fraction increases with
  differentiation (Spearman R=-0.624), confounding naive weighted-UTR
  calculations; protein-coding-only filter resolves this
- IG_C_gene/TR_C_gene biotype fix implemented and validated
- TMBIM6 identified as top APA gene in plasma cell differentiation —
  independently confirmed by PASTA (Cell 2024) in a separate PBMC dataset
  using a completely different statistical method (Dirichlet-multinomial
  residuals vs. Shannon entropy), providing strong cross-tool biological
  validation

### 4.2 IsoCAPE — architecture established

- MGA AluCE signal characterized as tumor-specific antisense Alu
  transcription (cancer hypomethylation) via strand flag analysis
- Internal priming filter designed; PAS_WINDOW=350bp derived from
  first principles using 10x library anatomy
- KDE clustering (bandwidth 100bp CE/PA, 50bp AluCE) adopted to replace
  sliding window for site aggregation
- CE_candidate / PA_candidate site type classification implemented

### 4.3 IsoDecipher-GPT (IsoFormer precursor) — architecture designed

- 100bp window, 6-mer sliding window tokenizer (stride=1, 95 tokens, vocab=4,096)
- Tier 1 Transformer with span masking MLM
- Tier 2 per-gene histogram pooling (over mean pooling to preserve
  proximal/distal site distribution)
- Files completed: bam_to_parquet.py, kmer_tokenizer.py, vocab.py
- GitHub: renegibson/isodecipher-gpt

---

## SECTION 5: SPECIFIC AIMS (NIH SBIR format, 1-page target)

*[Draft — to be tightened to exactly 1 page for NIH submission]*

**Overall goal**: Develop IsoMatrix into a validated, commercially deployable,
on-premise platform for single-cell alternative polyadenylation (APA) analysis,
enabling pharmaceutical and clinical users to access APA-based drug target
discovery and biomarker quantification within their existing data governance
frameworks.

**Aim 1: Complete and validate IsoDecipher v1.0 for pharmaceutical deployment**
(Months 1-6)
Deliverables:
- Finalize IsoDecipher manuscript (Figure 3: protein-coding-only UTR analysis;
  Figure 4: C_state latent space; Figure 5: UTR velocity + non-coding fraction
  velocity; Figure 6: IPA/ALE site type analysis)
- Implement site_type (ALE/3'Ext/IPA) annotation in build_panel_features.py
- Implement PolyASite v3.0 Track A panel expansion (genome-wide, 62%+ stringency)
- Package IsoDecipher as on-premise Docker container (HIPAA-compliant,
  no data egress, versioned panel updates)
- Success criterion: Docker deployment validated at 2 external sites
  (target: UW Medicine Genomics Core + 1 pharma partner)

**Aim 2: Train and validate IsoFormer Phase I (APA-VAE)**
(Months 3-12)
Deliverables:
- Train APA-VAE on B cell differentiation + cancer datasets
  (B cell: 21,670 cells [current]; cancer: target 2 additional datasets
  via pharma consulting partnerships)
- Demonstrate C_state embedding outperforms GEX-only embedding for
  APA-relevant cell state discrimination (plasma cell subpopulation
  resolution; interferon response identification)
- Implement IsoScore v0.1: deviation of a cell's APA profile from a
  healthy reference distribution (healthy reference: PolyASite v3.0
  tissue-level profiles as prior)
- Success criterion: IsoScore v0.1 discriminates cancer from normal cells
  in ≥1 independent validation dataset

**Aim 3: Demonstrate IsoCAPE→IsoDecipher feedback loop on a disease-relevant
cryptic site**
(Months 6-18)
Deliverables:
- Process ≥1 prostate cancer or CRPC dataset through IsoCAPE
- Validate detection of AR-V7/CE3-type ALE events (known ground truth,
  well-characterized coordinates in literature) as a positive control
  for IsoCAPE's cryptic site discovery capability
- Promote validated novel sites to IsoDecipher Track B panel
- Demonstrate IsoDecipher quantification of newly-discovered site
  across single cells in the same dataset
- Success criterion: AR-V7/CE3 detection sensitivity ≥80% at ≤5% FDR
  in a dataset with known AR-V7+ and AR-V7- cells

---

## SECTION 6: FORMAT ADAPTATION NOTES

### For NIH SBIR (PA-24-185, September 5 2026 deadline):
- Section 5 Specific Aims → tighten to exactly 1 page (750-800 words)
- Section 1 → Significance section (3-4 pages of Research Strategy)
- Section 2 → Innovation section (1-2 pages)
- Aims 1-3 → Approach section (6-7 pages) with timeline/milestones
- Section 3 → Commercialization Plan (separate section, 12 pages max)
- Section 4 → Preliminary Data (weave into Approach; also add to Biosketch)
- Key citation to include: Cheng LC et al. 2020 Nat Commun (SCAP mechanism);
  Bi et al. 2026 NAR (benchmark); Fahmi et al. 2025 EMM (landscape review);
  Kowalski et al. 2024 Cell (PASTA, confirms TMBIM6 finding cross-validates)

### For NSF SBIR (BT5 Life Science Research Tools, November 4 2026):
- Project Pitch (to submit first, ~4 sections, 1500-3500 chars each):
  1. Problem: Gap 1-4 from Section 1 above (ecosystem, calibration,
     biotype, deployment)
  2. Solution: Three-layer IsoMatrix architecture (Section 2)
  3. Team: Rene Cheng PhD — unique intersection of APA biology,
     single-cell engineering, and clinical domain knowledge
  4. Milestones: Aims 1-3 condensed to 3 bullet points each
- Full proposal adapts Section 1-4 above into NSF format
  (15 pages total vs NIH's longer format)

### For Activate Fellowship Cohort 2027 (apply ~October 2026):
- Hardware angle (if confirmed eligible): IsoMatrix on-premise Docker +
  Taiwan AI chip inference acceleration = hardware appliance pathway
  (Section 2.4 Phase III + Section 3.2 Taiwan partnership)
- If hardware eligibility confirmed: Section 3.2 becomes the core
  "hardware vision" narrative, with software-first as the de-risking phase
- If not eligible: redirect effort to NIH/NSF, skip Activate
- ACTION REQUIRED BEFORE WRITING: email apply@activate.org to confirm
  whether "on-premise hardware appliance (software-defined, with planned
  AI chip co-development with Taiwan partner)" qualifies as hardware-based

---

## APPENDIX: KEY REFERENCES (verified, with DOIs)

1. Cheng LC, Zheng D, Baljinnyam E, et al. "Widespread transcript shortening
   through alternative polyadenylation in secretory cell differentiation."
   Nat Commun 11:3182 (2020). DOI: 10.1038/s41467-020-16959-2
   [SCAP mechanism — primary motivating reference]

2. Kowalski MH*, Wessels HH*, Linder J*, et al. "Multiplexed single-cell
   characterization of alternative polyadenylation regulators."
   Cell 187(16):4408-4425.e23 (2024). DOI: 10.1016/j.cell.2024.06.005
   [PASTA — confirms TMBIM6 cross-validation; Dirichlet-multinomial approach]

3. Bi X, Chen Z, Ye M, et al. "Benchmarking computational methods for
   identifying and quantifying polyadenylation sites from 3' tag-based
   single-cell RNA-seq data."
   Nucleic Acids Research 54(9):gkag490 (2026). DOI: 10.1093/nar/gkag490
   [Independent benchmark — validates depth-robustness of annotation-anchored tools]

4. Fahmi M, Saha A, Song S, et al. "Computational methods for alternative
   polyadenylation and splicing in post-transcriptional gene regulation."
   Exp Mol Med 57(8) (2025). DOI: 10.1038/s12276-025-01496-z
   [Landscape review — Table 1 (bulk, 27 tools) + Table 2 (single-cell, 8 tools)]

5. Moon J, Herrmann CJ, Mironov V, Zavolan M. "PolyASite v3.0."
   Nucleic Acids Research 53(D1):D197-D204 (2025). DOI: 10.1093/nar/gkae1043
   [PolyASite v3.0 — scRNA-seq-derived, recommended for Track A panel expansion]

6. Bryce-Smith S, et al. "TDP-43 regulates cryptic splicing and
   polyadenylation in ALS-FTD."
   Nat Neurosci 28:2190-2200 (2025). DOI: 10.1038/s41593-025-02050-w
   [ALE/3'Ext/IPA classification framework + novel site validation rule]

7. Zhou B, et al. "SCAPE: a mixture model revealing single-cell
   polyadenylation diversity and cellular dynamics during cell differentiation
   and reprogramming." Nucleic Acids Res 50(11):e66 (2022).
   DOI: 10.1093/nar/gkac167
   [SCAPE — single Gaussian insert-size model; Section 13 technical deep-dive]

8. Fansler MM, Mitschka S, Mayr C. "scUTRquant."
   Nat Commun 15:4050 (2024). [scUTRquant — closest annotation-anchored precedent]
