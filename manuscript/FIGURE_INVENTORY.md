# Figure & Table Inventory — submission reorganisation (round 15 / Brain-Stroke build)

Mapping of the round-14 internal figures (Fig0–Fig41 + aux) into **8 main figures**
(each main figure = a set of *independent*, Nature-styled subpanels, no composite) and
**supplementary figures**, per §10 of the round-14 manuscript. All figures re-rendered:
Arial + Okabe-Ito, height 5 cm, width 5 cm (square) or 6.5 cm (landscape), **no notes on
the image** (in-figure suptitles/titles stripped; captions live in the manuscript).

Source file → panel name convention: each internal figure becomes one or more panels
named `FigN_<a/b/c>.png` (main) or `FigSN_<a/b/c>.png` (supplementary).

## MAIN FIGURES

- **Figure 1 — Overview / HRR framework (schematic).**
  - 1a HRR→RBC/heme-stress→ferroptosis hypothesis schematic  (rebuild from `Fig10_axis_schematic`)
  - 1b Multi-layer study design + cohorts schematic           (rebuild from `Fig0` design elements / new)

- **Figure 2 — Acute whole-blood heme–immune remodelling is largely composition-driven.**
  - 2a Volcano, GSE36791 DEGs                       (from `Fig1_volcano_GSE36791`)
  - 2b Candidate heatmap (DEG ∩ red-cell signature) (from `Fig3_heatmap_candidates`)
  - 2c Enrichment of 32 candidates                  (from `Fig4_enrichment` panel A)
  - 2d Hallmark/Reactome GSEA NES                   (from `Fig4_enrichment` panel B / `Fig5`)
  - 2e Marker deconvolution, immune shift           (from `Fig14_bulk_deconvolution`)
  - 2f Composition-adjusted DE (heme genes)         (from `Fig19_composition_adjusted`)

- **Figure 3 — Single-cell brain atlas localises heme handling to myeloid, ferroptosis to neurons.**
  - 3a Heme/iron module score by cell type          (from `Fig8_scrna_summary`)
  - 3b Cell-type-specific DE (neuronal ferroptosis) (from `Fig12_scrna_deep`)
  - 3c Temporal dynamics Sham→24h→72h               (from `Fig13_scrna_temporal`)
  - 3d Cell-state transitions                       (from `Fig20_cellstate`)

- **Figure 4 — HP→CD163 refines into a heme-handling program (MRC1 + DAMP).**
  - 4a CellChat-style ligand–receptor test          (from `Fig16_cellchat_style`)
  - 4b NicheNet-style ligand→target ranking         (from `Fig25_nichenet_style`)
  - 4c CITE-seq CD163 protein vs RNA                (from `Fig34_cite_seq_cd163`)

- **Figure 5 — HDFI: a conditional heme-detoxification-failure state (leakage disclosed).**
  - 5a NRF2/BACH1/HMOX1 regulon by cell type        (from `Fig15_regulon_bach1_nrf2`)
  - 5b HDFI vs neuronal ferroptosis-execution       (from `Fig18_hdfi_centrepiece`)
  - 5c Pathway-specificity + PCA factor audit       (from `Fig36_hdfi_latent_axis`)
  - 5d Robustness + GPX4 leakage audit              (from `Fig41_hdfi_robustness`)

- **Figure 6 — Regulatory (NRF2/ARE) support for HMOX1 activation.**
  - 6a ARE/NRF2 promoter density                    (from `Fig21_regulatory_motif`)
  - 6b ENCODE NFE2L2 ChIP enrichment at HMOX1       (from `Fig26_nrf2_chip_binding`)
  - 6c Macrophage ATAC: HMOX1 openness + ARE peaks  (from `Fig32_scatac_hmox1`)

- **Figure 7 — Ferroptotic execution lipids measured; external-validation limits stated.**
  - 7a ox-PE execution lipidomics                   (from `Fig33_lipidomics_ferroptosis`)
  - 7b CSF proteome axis-protein detection          (from `Fig37_proteomics_hdfi`)
  - 7c CSF ferritin ↔ neuronal-injury association   (from `Fig40_csf_hdfi_neuronal_injury`)
  - 7d HDFI external ICH validation (null)          (from `Fig23_hdfi_external`)

- **Figure 8 — Integrated model (associated-with, not causal).**
  - 8a Redrawn graphical abstract                   (from `Fig0_graphical_abstract`, state framing)

## SUPPLEMENTARY FIGURES

- **Figure S1** — Signed WGCNA network, SAH module GS/MM                (`Fig7_wgcna`)
- **Figure S2** — ML consensus hub rankings (LASSO/RF/SVM-RFE)          (`Fig5_ml_ranking`)
- **Figure S3** — Molecular score internal CV + external transfer AUC   (`Fig6_score_and_validation`)
- **Figure S4** — XGBoost feature-set comparison                        (`Fig9_xgb`)
- **Figure S5** — Microglia/macrophage DPT pseudotime (negative)        (`Fig11_scrna_trajectory`)
- **Figure S6** — FerrDb-categorised ferroptosis GSVA execution         (`Fig17_ferroptosis_gsva`)
- **Figure S7** — Two-sample MR (HP→SAH) forest + sensitivity           (`Fig22_mr_forest`)
- **Figure S8** — pQTL–MR panel, heme-axis proteins × SAH/ICH           (`Fig24_mr_panel`)
- **Figure S9** — Spatial co-localisation null (GSE254161 Visium)       (`Fig27_spatial_coloc`)
- **Figure S10** — Drug-target availability shortlist (OpenTargets)     (`Fig28_drugrepurposing_targets`)
- **Figure S11** — Rat stroke plasma lipidomics (arachidonate)          (`Fig29_lipidomics_stroke`)
- **Figure S12** — LINCS L1000 signature-reversal ranking               (`Fig35_drug_reversal_lincs`)
- **Figure S13** — Human SAH-CSF lipidome (verbatim PUFA/ether-PE)      (`Fig38_sah_csf_lipidomics`)
- **Figure S14** — Mouse SAH immune-atlas HDFI temporal program         (`Fig39_mouse_sah_atlas`)
- **Figure S15** — scRNA UMAP cell type / heme score / marker dotplot   (a `umap_celltype`, b `umap_heme`, c `dotplot__markers`)

## TABLES

MAIN (three-line):
- **Table 1** — Data and design (cohorts / accessions / tissue / comparison / platform / n)

SUPPLEMENTARY (three-line):
- **Table S1** — 32 candidate genes (log2FC, P, adj-P, signature membership)
- **Table S2** — Hub genes and ML ranking (LASSO/RF/SVM-RFE/SVM-Perm)
- **Table S3** — GSEA pathways (NES, FDR q, leading-edge)
- **Table S4** — WGCNA module–trait correlation + hub GS/MM
- **Table S5** — MR instruments + estimates (HP→SAH; IVW/Egger/WM/Q)
- **Table S6** — pQTL–MR panel results (10 pairs)
- **Table S7** — Cell-type-specific DE, selected myeloid/neuronal genes
- **Table S8** — Regulon / detox-failure index by cell type
- **Table S9** — FerrDb ferroptosis gene categories used
- **Table S10** — CSF proteome axis proteins (PXD030593) + neuronal-injury correlations
- **Table S11** — Ferroptosis execution lipids (ST003337/ST004928 ox-PE)
- **Table S12** — HDFI robustness & leakage audit numbers
- **Table S13** — Public data-accession inventory (GEO / PRIDE / GWAS-Catalog / Metabolomics-Workbench / ENCODE / JASPAR)

## DROPPED FROM SUBMISSION (planning / reviewer-round artefacts — "notes")
- Round-14 header ("Working manuscript… round 14"), the §10 submission-architecture
  table, the §10-B reviewer-attack-map table, the §7 "honest evaluation of the plan"
  table, and all "the reviewer asked / author-must-produce / this round / Round N"
  meta-commentary. Their substantive content is folded into Results, Limitations and
  Future directions as normal prose.
