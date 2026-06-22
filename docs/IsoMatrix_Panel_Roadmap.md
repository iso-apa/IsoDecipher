# IsoMatrix Panel Expansion Roadmap
*Compiled 2026-06-22 — based on APA tool landscape research + panel validation session*
*For use in: build_panel_features.py development, manuscript methods, SBIR application*

---

## 1. Design Philosophy

The IsoMatrix panel is NOT a static reference file.
It is a versioned, source-annotated, continuously-expanding asset that:
- Has a clear provenance for every feature (which source, which evidence)
- Supports multiple confidence tiers that users can filter by
- Grows systematically as IsoCAPE discovers novel sites
- Is fully reproducible (panel_config stored in adata.uns)

Core principle:
"A feature in the panel must have at least one of:
 (a) a GTF transcript 3' end annotation, or
 (b) independent bulk 3' end sequencing support (PolyASite v2.0), or
 (c) de novo discovery + validation in real data (IsoCAPE), or
 (d) curated literature evidence with precise coordinates (e.g. AR-V7/CE3)"

---

## 2. Panel Source Architecture

Every feature in the panel carries a `panel_source` tag.
Current and planned values:

| panel_source | Description | Confidence | Status |
|---|---|---|---|
| GTF | Ensembl GRCh38.115 transcript 3' end clustering (125bp tolerance) | Medium | Current (v1.0) |
| PA2 | PolyASite v2.0 site (bulk 3' end seq, not in GTF panel) | High | Planned (v1.1) |
| GTF+PA2 | GTF site confirmed by v2.0 within 125bp | Highest | Added via annotation (v1.1) |
| IsoCAPE | De novo discovery via IsoCAPE pipeline, internal validation only | Low-Medium | Future (v2.0) |
| IsoCAPE+PA2 | IsoCAPE discovery confirmed by v2.0 within 100nt | Medium-High | Future (v2.0) |
| IsoCAPE+hexamer | IsoCAPE discovery with conserved polyA signal hexamer | Medium | Future (v2.0) |
| Curated | Literature-curated known disease site (e.g. AR-V7/CE3) | Verified | Ongoing |

---

## 3. Panel Version History and Roadmap

### v1.0 (Current — completed 2026-06-19)
Source: GTF only (Ensembl GRCh38.115)
Features: 52,105 isoform groups across 15,102 genes
Key fixes in this version:
- CODING_BIOTYPES expansion: IG_C_gene/TR_C_gene/IG_V/D/J_gene now
  correctly classified as protein_coding (fixes IGH zero-UTR bug)
- assign_reads.py: user_label (Secreted/Membrane) removed from feature name,
  now stored only in panel metadata column
- utr_source column: protein_coding / non_coding / mixed
- avg_spliced_utr column: physical UTR length in bp (now correct for IGH)
- v2_matched, v2_usage, v2_site_class columns: PolyASite v2.0 annotation
  (63.7% of panel features matched within 125bp)

Validated findings:
- protein_coding v2_matched rate: 70.1%
- non_coding v2_matched rate: 37.8% (lower, consistent with annotation quality)
- IGH G0 (Secreted) UTR << G1 (Membrane) UTR — SCAP mechanism confirmed
- weighted_utr_pc_with_igh vs bcell_score: R=+0.574 (SCAP signal, correct
  direction after IGH fix)
- Non-coding RNA fraction vs bcell_score: R=0.299 (much weaker than the
  previously-reported R=-0.624, which was an artifact of IGH misclassification
  as non_coding)

### v1.1 (Planned — next sprint)
Source: GTF + PolyASite v2.0 Track A
New features: ~15,000-25,000 additional sites (TE/IN/DS, fraction>0.1)
Key changes:
- Add Track A sites from v2.0 (see Section 4)
- Add site_type column: ALE / 3'Ext / IPA (see Section 5)
- Add functional_class column (see Section 6)
- Add panel_source column to all features
- Update panel_config stored in adata.uns (see Section 8)
Target: panel size ~65,000-75,000 features

### v2.0 (Future — after IsoCAPE development complete)
Source: GTF + PA2 + IsoCAPE discoveries + disease tracks
New features: IsoCAPE-validated novel sites per cancer/disease type
Key changes:
- Disease tracks as separate CSV files, merged at build time
- IsoCAPE feedback loop operational (Section 7)
- panel_source fully populated for all features
- scAPAdb cross-validation layer (if B cell/plasma cell data confirmed)
Target: panel size variable by user configuration

---

## 4. Track A — PolyASite v2.0 Baseline Expansion

### 4.1 Data source
File: /Volumes/Lexar/reference/hg38/polyasite2_hg38.bed
Format: chrom, start, end, name (chrom:coord:strand), usage (TPM),
        strand, fraction (% of 221 libraries supporting this site),
        count, score2, site_class (IN/TE/AE/EX/IG/AI/DS/AU), signal
Total sites: 569,005
Data source: 221 bulk 3' end sequencing libraries, 9 protocols, >1.1B reads
Key column: fraction = percentage of 221 libraries detecting this site
  (NOT a percentage of cells; this is a cross-library support metric)
  fraction > 0.1 means at least 22 independent libraries confirmed this site

### 4.2 Filtering criteria (validated 2026-06-22)
Include:
  site_class IN ['TE', 'IN', 'DS']  — terminal-exon-related sites only
  fraction > 0.1                     — confirmed in >22 independent libraries
  NOT already in GTF panel (125bp tolerance check)
Exclude:
  IG (intergenic), AI (antisense intronic), AU (antisense UTR) — not gene-body APA
  AE (alternative exon) — lower confidence; move to Track B/Curated
  fraction <= 0.1 — insufficient cross-library support

Note on usage (TPM) vs fraction:
  fraction is the correct confidence metric here, NOT usage (TPM average).
  A site can have high usage in one cancer sample (inflating TPM average)
  but low fraction (few libraries detect it) — indicating it is disease-specific,
  NOT a healthy baseline site. The converse also holds.
  Cancer-specific sites with low fraction should go to Track B, not Track A.

### 4.3 Quantification (to be run)
After applying filters:
  Tier 1 (TE/IN/DS, usage>0.1 TPM): 39,711 candidate sites
  After 125bp clustering: 39,147 representative sites
  Sites with fraction>0.1 (to be re-quantified): TBD

  Note: Re-run filtering with fraction>0.1 as primary criterion,
  usage>0.1 TPM as secondary criterion. Expected final count: ~20,000-30,000.

### 4.4 Implementation steps
1. Filter v2.0 bed file to TE/IN/DS sites, fraction>0.1, not in panel
2. Run 125bp consecutive-gap clustering (same logic as build_panel_features.py)
3. Assign each cluster to a gene using GTF gene body overlap (gffutils)
   Priority: protein_coding > lncRNA > other
   Discard if no gene found (truly intergenic)
4. Compute avg_spliced_utr:
   Use GTF CDS 3' end of the assigned gene as reference
   UTR length = abs(site_coord - CDS_3prime_end)
   Flag as 0 if no CDS found (non-coding gene)
5. Merge with existing panel:
   Match to existing GTF features within 125bp — if match found,
   update v2_matched=True, v2_usage, v2_site_class, panel_source=GTF+PA2
   If no match — add as new row with panel_source=PA2
6. Assign site_type and functional_class (see Sections 5-6)
7. Update build_panel_features.py with --polyasite-v2 flag
   (Track A expansion is optional, off by default for minimal panel)

### 4.5 Version notes: v2.0 vs v2.0 of PolyASite
PolyASite v3.0 (Moon et al. 2025, NAR) exists and is derived from scRNA-seq
(10x Chromium 3', tissue-level). It is methodologically aligned with
IsoDecipher's input data type, but:
  - Site coordinates are LESS precise than v2.0 (inferred from read pile-up)
  - ~25% of v3.0 terminal exon sites have no v2.0 match
  - Intronic/intergenic sites show much higher v2.0/v3.0 discrepancy
Recommendation: Keep v2.0 as Track A source (highest precision);
use v3.0 only for cross-validation of IsoCAPE-discovered sites (Section 7).

---

## 5. site_type Column — ALE / 3'Ext / IPA Classification

### 5.1 Definition (Bryce-Smith et al. 2025, Nat Neurosci framework)
ALE (Alternative Last Exon):
  Defined by an UPSTREAM novel splice junction creating a new terminal exon.
  The polyA site is within this novel exon.
  Example: AR-V7/CE3 (intron 3 of AR gene)
  v2.0 site_class mapping: AE (primary), sometimes EX

3'Ext (3'UTR Extension):
  INDEPENDENT of splicing. Occurs downstream of an annotated distal 3'UTR.
  The polyA site extends the known terminal exon.
  Example: A more distal site than the GTF-annotated end of a terminal exon.
  v2.0 site_class mapping: TE (if within annotated terminal exon boundary),
                            DS (if downstream of annotated end)

IPA (Intronic Polyadenylation):
  Occurs within an intron, IN THE ABSENCE of upstream alternative splicing.
  Distinguished from ALE specifically by the lack of a novel splice junction.
  Example: IGHM secreted vs membrane isoform selection
           (uses an intronic PAS to terminate before the transmembrane exon)
  v2.0 site_class mapping: IN (when intron-located)

### 5.2 Implementation in build_panel_features.py
Function classify_site_type(rep_coord, transcript_exons, strand) already
written (session 2026-06-19). Key logic:
  - rep_coord in terminal exon → check if in last exon: tPA (→ 3'Ext/Canonical)
  - rep_coord in non-terminal exon → ePA (→ ALE candidate, check for splice junction)
  - rep_coord in intron → iPA (→ IPA)
  - rep_coord downstream of all exons → DS (→ 3'Ext, novel extension)

Note: v2.0 site_class used as prior/guide, but GTF-based determination
is the authoritative classification.

### 5.3 v2.0 site_class → ALE/3'Ext/IPA approximate mapping
TE → 3'Ext or Canonical_3UTR (within terminal exon)
IN → IPA (within non-terminal region, intronic location)
DS → 3'Ext (downstream extension beyond annotated end)
AE → ALE (alternative exon, check for splice junction evidence)
EX → ALE or ePA (exonic non-terminal)
IG → Novel/unknown (intergenic, exclude from Track A)
AI, AU → Antisense (exclude from panel entirely)

---

## 6. functional_class Column

Purpose: human-readable biological interpretation for downstream analysis,
         cell state modeling (IsoFormer), and clinical reporting.
         Maps directly to what the APA switch DOES, not just where it is.

Values and assignment logic:

| functional_class | Assignment rule | Biological meaning |
|---|---|---|
| Canonical_3UTR | protein_coding, site_type=3'Ext, G_mid (neither first nor last group) | Normal proximal/distal switching within terminal exon |
| UTR_shortening | protein_coding, site_type=3'Ext, G0 (most proximal site) | Shift to shorter UTR → increased mRNA stability (SCAP-associated) |
| UTR_lengthening | protein_coding, site_type=3'Ext, G_last (most distal site) | Shift to longer UTR → more regulatory elements, reduced stability |
| IPA_escape | site_type=IPA | Premature termination → truncated protein (cancer-associated) |
| ALE_switch | site_type=ALE | Novel last exon → functionally distinct protein isoform |
| NonCoding_APA | utr_source=non_coding | APA in lncRNA/non-coding context (tracked but excluded from UTR metrics) |

Note: Secretory_switch (IGH G0 Secreted vs G1 Membrane) is intentionally
NOT a separate functional_class. Reason: the G0/G1 naming within the
standard IPA/3'Ext framework already captures this, and the Secreted/Membrane
distinction is captured by user_label in the panel CSV (not the primary
analysis label). The UTR length difference (G0=~130bp, G1=~1400bp for IGHG)
is captured in avg_spliced_utr and is the key analytical feature.

---

## 7. Track B — Disease/Cancer-Specific Site Expansion

### 7.1 Design principles
Track B sites are stored as SEPARATE CSV files from the main panel,
merged at build_panel time via the --disease-tracks argument.
Each disease track has its own version and provenance.

Rationale for separation:
- Healthy baseline (Track A/GTF) ≠ disease-specific sites (Track B)
- Users doing normal biology studies should NOT have cancer-specific sites
  in their panel by default (adds noise, not signal)
- Cancer researchers can opt-in specific disease tracks
- Different diseases have different site lists — one panel file per disease

### 7.2 Disease track CSV schema
Columns (same as main panel, plus disease-specific metadata):
  gene, chrom, coord, strand, site_type, functional_class,
  avg_spliced_utr (estimated), utr_source,
  panel_source (Curated / IsoCAPE / IsoCAPE+PA2 / IsoCAPE+hexamer),
  evidence_pmid, cancer_type, clinical_significance, date_added,
  isocape_dataset_gse (if discovered by IsoCAPE)

Validation rule for IsoCAPE-discovered sites to be promoted to Track B:
  EITHER: within 100nt of a PolyASite v2.0 annotated PAS
  OR: contains a conserved polyA signal hexamer in terminal 100nt
  (per Bryce-Smith et al. 2025, Nat Neurosci — same rule they used for TDP-43
  cryptic APA validation)

### 7.3 Current Track B files (to be created)

known_sites_prostate_cancer.csv:
  AR-V7/CE3 — ALE event, intron 3 of AR gene
  Coordinates: chr X, intron 3 of AR (need to confirm precise CE3 coordinates)
  Mechanism: novel splice junction to cryptic exon 3 → constitutively active AR
  Clinical significance: FDA-relevant liquid biopsy biomarker for CRPC
  Evidence: multiple PMIDs (AR-V7 resistance literature)
  panel_source: Curated

known_sites_myeloma.csv:
  To be populated from:
  (a) GSE193531 (Nature Commun 2022) processed through IsoCAPE
  (b) GSE223060/GSE223061 (Cancer Research 2023, MM single-cell)
  Currently empty — placeholder for IsoCAPE output

known_sites_als_ftd.csv:
  227 TDP-43-regulated cryptic APA events from Bryce-Smith et al. 2025
  Note: these are CONDITION-SPECIFIC (TDP-43 depletion in iPSC neurons),
  NOT general cancer sites — use only for neurodegeneration research track

### 7.4 Future disease tracks (planned, not started)
  known_sites_breast_cancer.csv (IsoCAPE on breast cancer BAMs)
  known_sites_lung_cancer.csv
  known_sites_lymphoma.csv (B cell lymphoma scRNA-seq data)

### 7.5 IsoCAPE → IsoDecipher feedback loop
This is the core self-improving mechanism of IsoMatrix:

1. IsoCAPE processes a new cancer BAM dataset (e.g. GSE193531 MM data)
2. IsoCAPE outputs: novel site candidates with coordinates + site_type
3. Validation filter applied:
   - within 100nt of v2.0 PAS: panel_source = IsoCAPE+PA2
   - hexamer present: panel_source = IsoCAPE+hexamer
   - neither: stays as IsoCAPE (low confidence, Track B only after manual review)
4. Validated sites appended to appropriate disease track CSV
5. User runs build_panel with --disease-tracks flag
6. New sites integrated into panel with correct panel_source and metadata
7. IsoDecipher re-run with expanded panel → quantifies new sites
8. IsoFormer learns which of the new sites carry diagnostic signal
9. IsoFormer attention weights fed back to IsoCAPE as priors for next dataset

Result: each new cancer dataset makes the whole platform more sensitive,
without requiring re-discovery of already-known sites.

---

## 8. panel_config in adata.uns

Every h5ad produced by IsoDecipher must store the full panel provenance:

```python
adata.uns['panel_config'] = {
    'isodecipher_version': '1.1.0',
    'panel_version': '1.1',
    'build_date': '2026-06-22',
    'gtf': 'Homo_sapiens.GRCh38.115.gtf',
    'gtf_source': 'Ensembl release 115, GRCh38',
    'tolerance_bp': 125,
    'tolerance_validation': 'BC>0.6 across 6 independent datasets (Figure 2B)',
    'polyasite_v2': {
        'used': True,
        'file': 'polyasite2_hg38.bed',
        'filters': {
            'site_class': ['TE', 'IN', 'DS'],
            'fraction_threshold': 0.1,
            'note': 'fraction = % of 221 libraries supporting site'
        }
    },
    'disease_tracks': [],          # list of disease track CSV filenames used
    'isocape_sites': None,         # GSE accession(s) if IsoCAPE sites included
    'coding_biotypes': [           # CODING_BIOTYPES set used in build_panel
        'protein_coding', 'IG_C_gene', 'IG_V_gene', 'IG_D_gene', 'IG_J_gene',
        'IG_LV_gene', 'TR_C_gene', 'TR_V_gene', 'TR_D_gene', 'TR_J_gene'
    ],
    'n_features_total': 52105,
    'n_features_by_source': {
        'GTF': 52105,
        'PA2': 0,                  # will be updated in v1.1
        'GTF+PA2': 0,              # matched features get updated label
        'IsoCAPE': 0,
        'Curated': 0
    },
    'n_genes': 15102,
    'utr_source_distribution': {
        'protein_coding': 41820,
        'non_coding': 10259,
        'mixed': 26
    }
}
```

This makes any IsoDecipher output FULLY REPRODUCIBLE and auditable.
Clinicians, regulators, or collaborators receiving an h5ad can inspect
exactly what went into the panel that produced every feature in the matrix.

---

## 9. Paper Methods Statement (draft)

"The IsoMatrix panel is a versioned, multi-source reference of
cleavage/polyadenylation sites, designed to expand systematically as
new data becomes available. The current version (v1.1) integrates
three sources with distinct confidence levels: (i) GTF-derived sites from
Ensembl GRCh38 release 115, clustered using a 125bp tolerance window
validated by Bhattacharyya Coefficient analysis across six independent
datasets; (ii) PolyASite v2.0 sites (Herrmann et al. 2020, NAR) supported
by >10% of 221 independent bulk 3' end sequencing libraries (fraction>0.1),
restricted to terminal-exon-proximal classes (TE/IN/DS); and (iii) curated
disease-specific sites from the literature, including clinically-validated
APA/ALE events (e.g., AR-V7/CE3 in castration-resistant prostate cancer).
Each feature carries a panel_source annotation indicating its evidence basis.
Disease-specific sites are distributed as separate track files, allowing
users to include only cancer types relevant to their analysis.
The panel is designed to expand continuously as IsoCAPE (IsoMatrix's de novo
cleavage site discovery module) identifies and validates novel sites in
new datasets, creating a closed feedback loop between discovery (IsoCAPE)
and quantification (IsoDecipher) (Figure 1)."

---

## 10. Immediate Action Items (ordered by priority)

### This week:
1. Re-filter v2.0 candidates using fraction>0.1 (not just usage>0.1 TPM)
   → Run quantification to see how many sites remain
2. Implement gene assignment for Track A candidates (gffutils, GTF gene body)
3. Implement avg_spliced_utr calculation for new sites
4. Implement site_type classification (classify_site_type() function, already written)
5. Implement functional_class assignment
6. Add panel_source column to all features
7. Update build_panel_features.py with Track A expansion logic
8. Re-run build_panel → assign_reads → integrate (final v1.1 pipeline run)
9. Verify panel_config stored in adata.uns after integrate

### Next sprint:
10. Create known_sites_prostate_cancer.csv (AR-V7/CE3 coordinates)
11. Download GSE193531 (MM precursor stages scRNA-seq)
12. Run IsoCAPE on GSE193531 to populate known_sites_myeloma.csv (first entries)
13. Add --disease-tracks CLI argument to build_panel_features.py
14. Add --polyasite-v2 CLI argument to build_panel_features.py

### Future:
15. scAPAdb cross-validation (check if B cell/plasma cell datasets exist)
16. TC3A download → identify B cell/lymphoma relevant APA genes
    → find in v2.0 → tag as cancer_relevant=True in panel
17. IsoFormer Phase I (APA-VAE) training once panel v1.1 is stable

---

## 11. Key References

1. PolyASite v2.0: Herrmann et al. 2020, NAR. DOI: 10.1093/nar/gkz918
2. PolyASite v3.0: Moon et al. 2025, NAR. DOI: 10.1093/nar/gkae1043
3. ALE/3'Ext/IPA framework: Bryce-Smith et al. 2025, Nat Neurosci.
   DOI: 10.1038/s41593-025-02050-w
4. SCAP mechanism (foundational): Cheng LC et al. 2020, Nat Commun.
   DOI: 10.1038/s41467-020-16959-2
5. AR-V7 mechanism: multiple PMIDs (see known_sites_prostate_cancer.csv)
6. scAPAdb: Zhu et al. 2022, NAR. DOI: 10.1093/nar/gkab795
7. TC3A: Feng et al. 2018, Nat Commun
8. TREND-DB: Marini et al. 2021
9. GSE193531 (MM precursor stages): Nature Commun 2022
10. GSE223060/61 (MM single-cell): Cancer Research 2023
