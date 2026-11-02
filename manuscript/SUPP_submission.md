# Supplementary Information

# A heme-detoxification-failure state links red-cell dysfunction to ferroptotic vulnerability after subarachnoid hemorrhage

## Supplementary Figures

**Figure S1. Co-expression network.** (A) module-eigengene correlation with SAH status; (B) gene-significance / module-membership for heme–iron genes in the SAH-correlated module.

**Figure S2. Machine-learning hub selection.** (A) consensus machine-learning hub selection across LASSO, random-forest and SVM-RFE.

**Figure S3. Molecular-score validation.** (A–C) molecular-score internal cross-validation and external transfer to the remote-aSAH cohort.

**Figure S4. Gradient-boosting comparison.** (A) gradient-boosting comparison of four feature sets (internal versus external AUC).

**Figure S5. Single-cell reference panels.** (A) UMAP of the SAH brain atlas coloured by cell type; (B) UMAPs coloured by the heme/iron module score and by Hp, Cd163 and Hmox1 expression; (C) canonical-marker dotplot across cell types.

**Figure S6. Mouse SAH immune-atlas HDFI program.** (A) immune composition across sham, 24 h and 72 h; (B) HDFI by cell type and timepoint; (C) cell-type × timepoint score heatmap.

**Figure S7. Pseudotime trajectory (negative).** (A–C) microglia/macrophage diffusion-pseudotime trajectory (homeostatic-rooted).

**Figure S8. FerrDb ferroptosis scoring.** (A–C) FerrDb-categorised ferroptosis Driver/Suppressor/Iron and execution scores.

**Figure S9. Stroke plasma lipidomics.** (A) rat post-stroke plasma lipidomics (arachidonate and related metabolites).

**Figure S10. Human SAH-CSF lipidome (verbatim).** (A) human SAH-CSF lipidome PUFA/ether-PE elevations (reported from the source study) with the mouse execution-lipid fold-changes as quantitative anchor.

**Figure S11. Mendelian-randomization of haptoglobin on SAH.** (A) two-sample Mendelian-randomization of haptoglobin on SAH (instrument and summary forest).

**Figure S12. pQTL–Mendelian-randomization panel.** (A) Mendelian-randomization panel of heme-axis pQTLs on SAH and intracerebral hemorrhage.

**Figure S13. Spatial co-localisation (null).** (A–C) spatial co-localisation test on an independent brain-injury Visium dataset (null).

**Figure S14. Drug-target availability.** (A) drug-target availability shortlist for the axis proteins (OpenTargets).

**Figure S15. LINCS signature reversal.** (A, B) LINCS L1000 direction-aware signature-reversal ranking.

## Supplementary Tables

**Table S1.** Public data-accession inventory.

| Type | Accession | Content |
|---|---|---|
| GEO | GSE36791 | acute SAH whole blood (GPL10558); discovery |
| GEO | GSE73378 | remote/past aSAH whole blood (GPL10558); test |
| GEO | GSE311852 | mouse SAH single-cell brain atlas (sham/24 h/72 h) |
| GEO | GSE163256 | ICH blood monocytes (external validation) |
| GEO | GSE296792 | ICH whole blood (external validation) |
| GEO | GSE269478 | human monocyte RNA+ADT CITE-seq |
| GEO | GSE172116 | human macrophage bulk ATAC-seq (LPS contrast) |
| GEO | GSE254161 | mouse cortical-injury 10x Visium (spatial) |
| PRIDE | PXD030593 | human SAH-CSF label-free DDA proteome |
| GWAS | GCST90019428 | haptoglobin plasma pQTL (Eur, N = 10,708) |
| GWAS | GCST90473583 | SAH GWAS (Eur, N = 458,440) |
| GWAS | GCST90473584 | ICH GWAS (Non-Finnish Eur, N = 458,440) |
| GWAS | GCST90270865 / GCST90451683 / GCST90469465 / GCST90470840 | ferritin / hepcidin / HMOX1 / TFRC pQTLs |
| ENCODE | ENCSR584GHV / ENCSR707IUN / ENCSR488EES / ENCSR197WGI | NFE2L2 ChIP-seq (A549/HeLa/HepG2/IMR-90) |
| JASPAR | MA0150.3 | NFE2L2 PWM |
| Metabolomics-Workbench | ST002968 / ST003337 / ST004928 | stroke plasma / MEF / MDSC lipid panels |

**Table S2.** The 32 candidate genes at the intersection of SAH differentially expressed genes and the red-cell/HRR signature.

| Gene | log2FC | P | adj. P | Signature source |
|---|---|---|---|---|
| HSPA9 | -0.51 | 9.6e-06 | 0.00057 | - |
| MAPK14 | 0.58 | 1e-05 | 0.0006 | - |
| CDK6 | -0.39 | 1.1e-05 | 0.0006 | - |
| HP | 1.3 | 1.3e-05 | 0.00067 | curated_iron_handling_oxidative |
| RPS19 | -0.58 | 3.2e-05 | 0.001 | - |
| BCL11A | -0.37 | 4.3e-05 | 0.001 | curated_erythropoiesis_maturation |
| PRMT1 | -0.37 | 4.4e-05 | 0.001 | - |
| STAT5B | 0.4 | 7.5e-05 | 0.002 | - |
| HMGB2 | 0.55 | 0.0001 | 0.002 | - |
| SLC40A1 | 0.67 | 0.00013 | 0.003 | curated_iron_handling_oxidative |
| HSPA1A | 0.61 | 0.00021 | 0.004 | - |
| RPS14 | -0.52 | 0.00035 | 0.005 | - |
| TGFBR3 | -0.53 | 0.00036 | 0.006 | - |
| SPI1 | 0.41 | 0.00049 | 0.007 | - |
| GTPBP4 | -0.31 | 0.00065 | 0.008 | - |
| CFLAR | 0.5 | 0.00069 | 0.009 | - |
| NFE2 | 0.43 | 0.00081 | 0.010 | curated_erythropoiesis_maturation |
| HSPA1B | 0.4 | 0.001 | 0.012 | - |
| STAT3 | 0.36 | 0.001 | 0.012 | - |
| FTH1 | 0.29 | 0.001 | 0.013 | curated_iron_handling_oxidative |
| SLC11A1 | 0.7 | 0.001 | 0.015 | - |
| ACVR1B | 0.27 | 0.001 | 0.015 | - |
| SLC25A37 | 0.48 | 0.002 | 0.016 | curated_heme_biosynthesis |
| SLC25A5 | -0.33 | 0.003 | 0.024 | - |
| UBA5 | -0.12 | 0.003 | 0.027 | - |
| JAK2 | 0.43 | 0.005 | 0.034 | - |
| CDK5RAP3 | -0.23 | 0.005 | 0.034 | - |
| LYN | 0.31 | 0.005 | 0.034 | - |
| LYAR | -0.24 | 0.006 | 0.039 | - |
| GYPB | 0.85 | 0.007 | 0.043 | curated_RBC_membrane_structure |
| SH2B3 | -0.25 | 0.007 | 0.047 | - |
| ETS1 | -0.97 | 2.7e-07 | 9.2e-05 | - |

**Table S3.** Gene-set enrichment analysis (preranked) of the discovery cohort, selected pathways.

| Pathway | NES | Nominal P | FDR q |
|---|---|---|---|
| Neutrophil Degranulation | 2.22 | <0.001 | <0.001 |
| Complement and coagulation cascades | 2.11 | <0.001 | 0.001 |
| TNF-alpha Signaling via NF-kB | 1.92 | <0.001 | 0.013 |
| Reactive Oxygen Species Pathway | 1.89 | <0.001 | 0.017 |
| heme Metabolism | 1.89 | <0.001 | 0.018 |
| Hypoxia | 1.45 | 0.004 | 0.230 |
| Complement | 1.43 | 0.007 | 0.229 |

**Table S4.** Whole-blood heme outcomes before and after adjustment for neutrophil, monocyte, and T pseudo-fractions.

| Outcome | P (unadjusted) | P (adjusted) | q (adjusted) |
|---|---|---|---|
| HP | 0.001 | 0.868 | 0.954 |
| HMOX1 | 0.237 | 0.411 | 0.831 |
| FTH1 | 0.005 | 0.694 | 0.954 |
| FTL | 0.079 | 0.439 | 0.831 |
| SLC40A1 | 0.002 | 0.865 | 0.954 |
| MARC1 | 2.6e-05 | 0.453 | 0.831 |
| GPX4 | 0.089 | 0.297 | 0.831 |
| SLC25A37 | 0.002 | 0.279 | 0.831 |
| NFE2 | 0.00062 | 0.569 | 0.893 |
| HEME_MODULE | 6.2e-08 | 0.171 | 0.831 |
| HDFI | 0.576 | 0.954 | 0.954 |

**Table S5.** WGCNA module-eigengene correlation with SAH status (top modules).

| Module | Genes | r with SAH | P |
|---|---|---|---|
| ME133 | 1068 | 0.55 | 3.8e-06 |
| ME60 | 61 | 0.55 | 5.2e-06 |
| ME65 | 760 | 0.52 | 1.4e-05 |
| ME124 | 121 | 0.38 | 0.002 |
| ME132 | 108 | 0.36 | 0.004 |
| ME131 | 36 | 0.33 | 0.010 |
| ME52 | 387 | 0.29 | 0.023 |
| ME58 | 105 | 0.24 | 0.059 |
| ME94 | 53 | 0.17 | 0.185 |
| ME36 | 79 | 0.15 | 0.248 |

**Table S6.** Consensus machine-learning hub genes (top-ranked by at least three of LASSO / random-forest / SVM-RFE criteria).

| Gene | log2FC | adj. P | Ranking votes |
|---|---|---|---|
| CST7 | 1.26 | 3.7e-05 | 3 |
| MARC1 | 1.18 | 8.1e-05 | 3 |
| NFE2 | 0.43 | 0.010 | 3 |
| IL2RB | -1.16 | 2.2e-05 | 3 |

**Table S7.** Cell-type-specific differential expression (SAH vs sham), selected heme-handling and ferroptosis genes.

| Cell type | Gene | log2FC | adj. P |
|---|---|---|---|
| Neutrophil | Hp | 2.01 | 2.7e-13 |
| Neutrophil | Hmox1 | 6.57 | 1.2e-24 |
| Neutrophil | Fth1 | 0.43 | 4.4e-18 |
| Neutrophil | Ftl1 | 0.71 | 1.7e-21 |
| Neutrophil | Lcn2 | 2.23 | 7.4e-16 |
| Neutrophil | Mmp9 | 1.75 | 3.3e-10 |
| Neuron | Slc7a11 | 3.63 | 8.4e-21 |
| Neuron | Gpx4 | -0.24 | 0.030 |
| Neuron | Acsl4 | 1.53 | 1.7e-23 |
| Neuron | Ptgs2 | 5.72 | 5.5e-24 |
| Neuron | Hmox1 | 3.79 | 1.9e-100 |

**Table S8.** Regulon and heme-detoxification-failure activity, SAH vs sham, by cell type.

| Cell type | Metric | Delta | Cohen d | P |
|---|---|---|---|---|
| Microglia/Mac | NRF2_activity | 0.31 | 1.14 | <0.001 |
| Microglia/Mac | BACH1_activity | -0.1 | -0.09 | 0.226 |
| Microglia/Mac | BACH1_expr | 0.33 | 0.33 | 5.4e-42 |
| Microglia/Mac | HMOX1 | 0.87 | 1.15 | <0.001 |
| Microglia/Mac | DetoxFailure | 0.48 | 0.54 | 4.3e-110 |
| Neutrophil | NRF2_activity | 0.76 | 1.81 | 2.9e-28 |
| Neutrophil | BACH1_activity | 0.07 | 0.06 | 0.535 |
| Neutrophil | BACH1_expr | 0.88 | 0.99 | 1.1e-15 |
| Neutrophil | HMOX1 | 1.36 | 1.8 | 1.2e-27 |
| Neutrophil | DetoxFailure | -0.67 | -0.49 | 0.00026 |
| Neuron | NRF2_activity | 0.54 | 1.33 | 4.2e-113 |
| Neuron | BACH1_activity | -0.06 | -0.06 | 0.535 |
| Neuron | BACH1_expr | 0.69 | 0.75 | 6.7e-47 |
| Neuron | HMOX1 | 1.13 | 1.32 | 2.5e-115 |
| Neuron | DetoxFailure | -0.02 | -0.02 | 0.504 |
| Astrocyte | NRF2_activity | 0.11 | 0.35 | 2.4e-39 |
| Astrocyte | BACH1_activity | 0.72 | 0.67 | 5.1e-92 |
| Astrocyte | BACH1_expr | 1.01 | 1.06 | 3.6e-229 |
| Astrocyte | HMOX1 | 0.61 | 1.01 | 1.3e-216 |
| Astrocyte | DetoxFailure | 0.74 | 0.67 | 6e-92 |
| Oligodendrocyte | NRF2_activity | 0.32 | 0.73 | 9.4e-07 |
| Oligodendrocyte | BACH1_activity | 0.28 | 0.22 | 0.450 |
| Oligodendrocyte | BACH1_expr | 0.54 | 0.56 | 2.1e-08 |
| Oligodendrocyte | HMOX1 | 0.8 | 0.89 | 2.5e-51 |
| Oligodendrocyte | DetoxFailure | 0.64 | 0.39 | 0.020 |

**Table S9.** HDFI robustness and GPX4-leakage audit on the neuronal association (raw-count-z pipeline).

| Test | Result |
|---|---|
| Canonical neuronal association | rho 0.337, P 9.2e-44 |
| Leakage-free (no shared gene) | rho 0.071, P 0.005 |
| GPX4 removed from outcome only | rho 0.067, P 0.007 |
| Sign-reversed index | rho -0.337 |
| Random-weighting null (n=5000) | one-sided empirical P 0.071; best random rho 0.457 |
| PCA-independence | partial rho 0.295, P 1.5e-33 |
| Leave-one-out (drop Hmox1) | rho 0.399, P 1.8e-62 |
| Leave-one-out (drop Fth1) | rho 0.274, P 4.5e-29 |
| Leave-one-out (drop Ftl1) | rho 0.278, P 8.3e-30 |
| Leave-one-out (drop Gpx4) | rho 0.095, P 0.00015 |

**Table S10.** Oxidised-PE ferroptosis execution markers across two public targeted metabolomics datasets.

| Analyte | Dataset | Class | Fold change | P | Cliff's delta |
|---|---|---|---|---|---|
| PE(16:0/20:4+2[O])T4_1I | ST003337 | ox-PE (oxidised phosphatidylethanolamine) | 17.1 | 0.057 | 1 |
| PE(18:0/20:4+1[O])T3_2I | ST003337 | ox-PE (oxidised phosphatidylethanolamine) | 1.69 | 0.229 | 0.67 |
| PE(18:0/20:4+2[O])T4_1I | ST003337 | ox-PE (oxidised phosphatidylethanolamine) | 24.1 | 0.057 | 1 |
| PE(18:0/20:4+2[O])T4_2I | ST003337 | ox-PE (oxidised phosphatidylethanolamine) | 18.7 | 0.057 | 1 |
| PE(18:0/20:4+2[O])T4_3I | ST003337 | ox-PE (oxidised phosphatidylethanolamine) | 27.8 | 0.057 | 1 |
| FA(20:4) | ST004928 | free arachidonate (AA) | 2.44 | 0.100 | 1 |
| PE(16:0_20:4) | ST004928 | arachidonoyl-PE substrate (16:0/20:4) | 0.791 | 0.100 | -1 |
| PE(18:0_20:4) | ST004928 | arachidonoyl-PE substrate (18:0/20:4) | 0.848 | 0.100 | -1 |
| PE(18:0_18:2) | ST004928 | linoleoyl-PE substrate (18:0/18:2 = LA esterified) | 0.802 | 0.100 | -1 |

**Table S11.** Human SAH cerebrospinal-fluid proteome: axis proteins detected, SAH vs control.

| Source | Protein | Role | log2FC | Mann-Whitney P |
|---|---|---|---|---|
| PXD030593 | HP | heme scavenger (haptoglobin) | 0.64 | 0.057 |
| PXD030593 | CD163 | heme-Hb endocytosis receptor | -0.53 | 0.114 |
| PXD030593 | FTH1 | iron storage (ferritin heavy) | -0.04 | 0.886 |
| PXD030593 | FTL | iron storage (ferritin light) | 0.5 | 0.486 |
| PXD030593 | TFRC | iron uptake (CD71) | 0.45 | 0.486 |
| PMC11267790 (verbatim) | HP | independent SAH-CSF cohort (authors' Table 1) | 1.62 | NA |
| PMC11267790 (verbatim) | TFRC | independent SAH-CSF cohort (authors' Table 1) | 0.83 | NA |

**Table S12.** Within-CSF association of axis proteins with the neuronal-injury protein score (PXD030593, n = 8 samples).

| Protein | rho vs injury score | P |
|---|---|---|
| HP | 0.02 | 0.955 |
| CD163 | -0.48 | 0.233 |
| FTH1 | 0.86 | 0.007 |
| FTL | 0.76 | 0.028 |
| TFRC | -0.24 | 0.570 |

**Table S13.** Two-sample Mendelian randomization: haptoglobin plasma level (GCST90019428) on SAH (GCST90473583).

| Estimate | Value |
|---|---|
| Instruments retained | 9 |
| IVW beta (log OR per SD) | -0.023 |
| IVW SE | 0.046 |
| IVW P | 0.625 |
| IVW OR (95% CI) | 0.978 [0.893, 1.07] |
| Cochran Q (df), P | 4.84 (8), 0.774 |
| MR-Egger beta, P | -0.021, 0.633 |
| MR-Egger intercept P | 0.943 |
| Weighted-median beta, P | 0.001, 1.000 |
| Mean / min instrument F | 335 / 30 |

**Table S14.** Mendelian-randomization panel: heme-axis pQTLs on SAH and intracerebral hemorrhage (IVW odds ratios).

| Exposure | Outcome | IVs | OR | 95% CI | P |
|---|---|---|---|---|---|
| Haptoglobin | Subarachnoid haemorrhage (ICD10 I60) | 9 | 0.978 | [0.893, 1.07] | 0.625 |
| Haptoglobin | Intracerebral haemorrhage (ICD10 I61) | 9 | 1 | [0.923, 1.09] | 0.961 |
| Ferritin (total) | Subarachnoid haemorrhage (ICD10 I60) | 79 | 1.1 | [0.854, 1.41] | 0.464 |
| Ferritin (total) | Intracerebral haemorrhage (ICD10 I61) | 79 | 0.943 | [0.751, 1.18] | 0.609 |
| Heme oxygenase 1 (HMOX1) | Subarachnoid haemorrhage (ICD10 I60) | 45 | 1.03 | [0.861, 1.23] | 0.743 |
| Heme oxygenase 1 (HMOX1) | Intracerebral haemorrhage (ICD10 I61) | 45 | 1.07 | [0.913, 1.26] | 0.390 |
| Transferrin receptor (TFRC) | Subarachnoid haemorrhage (ICD10 I60) | 37 | 1.13 | [0.943, 1.35] | 0.187 |
| Transferrin receptor (TFRC) | Intracerebral haemorrhage (ICD10 I61) | 37 | 1.17 | [0.991, 1.37] | 0.064 |
| Hepcidin (HAMP) | Subarachnoid haemorrhage (ICD10 I60) | 17 | 1.09 | [0.723, 1.64] | 0.686 |
| Hepcidin (HAMP) | Intracerebral haemorrhage (ICD10 I61) | 17 | 1.01 | [0.698, 1.46] | 0.966 |
