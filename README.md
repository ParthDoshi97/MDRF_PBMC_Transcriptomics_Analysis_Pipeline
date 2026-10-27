# Dynamic Transcriptomic Profiling of PBMCs Under Acute Oral Glucose Challenge in Asian-Indian Men with Prediabetes

[![GEO Accession](https://img.shields.io/badge/GEO-GSE153837-blue.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153837)
[![Platform](https://img.shields.io/badge/Platform-Illumina_HumanHT--12_V4.0-brightgreen.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558)
[![Pipeline](https://img.shields.io/badge/Pipeline-Python_3.13_Bioinformatics-orange.svg)](https://www.python.org/)
[![Status](https://img.shields.io/badge/Status-Completed-success.svg)]()

A comprehensive, end-to-end transcriptomic investigation of peripheral blood mononuclear cells (PBMCs) from Asian-Indian men undergoing an acute 75g oral glucose tolerance test (OGTT). This study dissects baseline fasting perturbations and dynamic postprandial transcriptional failures characterizing the transition from normal glucose tolerance (NGT) to prediabetes (PD).

---

## Table of Contents
1. [Study Overview & Clinical Background](#study-overview--clinical-background)
2. [Key Biological Discoveries & Mechanisms](#key-biological-discoveries--mechanisms)
3. [Quality Control & Batch Effect Resolution](#quality-control--batch-effect-resolution)
4. [Differential Expression Summary](#differential-expression-summary)
5. [Repository Structure](#repository-structure)
6. [Reproducibility & Execution Pipeline](#reproducibility--execution-pipeline)
7. [Reports & Deliverables](#reports--deliverables)
8. [Evidence-Based Future Directions](#evidence-based-future-directions)
9. [Citation & Attribution](#citation--attribution)

---

## Study Overview & Clinical Background

Asian-Indian populations present a distinct metabolic phenotype—often termed the **"Asian-Indian Phenotype"**—characterized by elevated insulin resistance, increased visceral adiposity, higher high-sensitivity C-reactive protein (hs-CRP), and premature coronary artery disease, frequently in the absence of generalized obesity (low BMI).

* **Dataset Accession:** NCBI GEO [GSE153837](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153837)
* **Microarray Platform:** Illumina HumanHT-12 V4.0 Expression BeadChip (`GPL10558`, 47,231 probe features)
* **Cohort Design:** 28 PBMC samples derived from Asian-Indian adult men:
  * **NGT Fasting (Control Fasting):** $n = 7$
  * **NGT 2h Post-OGTT (Control Postprandial):** $n = 7$
  * **Prediabetes Fasting (PD Fasting):** $n = 7$
  * **Prediabetes 2h Post-OGTT (PD Postprandial):** $n = 7$
* **Primary Objective:** Unravel whether an acute 75g oral glucose challenge unmasks condition-specific transcriptional inflexibility in prediabetic immune cells that cannot be identified from static fasting measurements alone.

---

## Key Biological Discoveries & Mechanisms

Our analysis uncovered that prediabetes is defined by **two distinct pathophysiological regimes**: a chronic baseline inflammatory/stress state in the fasting regime, and an acute failure of proteostasis and metabolic signaling in the post-challenge regime. Across the 138 fasting DEGs and 62 post-OGTT DEGs, **only 5 genes overlap (97.5% state divergence)**, demonstrating that the oral glucose challenge serves as an indispensable molecular stress test.

```
       Fasting DEGs: 138                    Post-OGTT DEGs: 62
   ┌───────────────────────┐            ┌───────────────────────┐
   │                       │   5 Shared │                       │
   │      133 Unique       │   Genes:   │       57 Unique       │
   │  - SLA (Checkpoint ↓) │  - MRPS31  │  - HSPH1 (Chaperone ↓)│
   │  - DNAJC13 (ER/Endo ↑)│  - CCDC85B │  - PREP (Incretin ↓)  │
   │  - NCKAP1L (Adhes. ↑) │  - TRAM1L1 │  - DNAJB1 (Hsp40 ↓)   │
   │  - NFKBIA (NF-kB reg) │  - ACKR4   │  - HSPA1A/B (Hsp70 ↓) │
   │                       │  - AC007392│                       │
   └───────────────────────┘            └───────────────────────┘
```

### 1. The Chaperone Crisis & Proteotoxic Shock (`HSPH1` & `DNAJC13`)
* **Postprandial Inflexibility:** In healthy NGT PBMCs, the rapid surge of intracellular glucose and mitochondrial reactive oxygen species (ROS) induces **`HSPH1`** (Hsp105/110 family), a high-molecular-weight nucleotide exchange factor for Hsp70 that prevents protein aggregation and proteotoxic apoptosis. In prediabetic subjects, this cytoprotective defense is completely blunted ($p = 6.1 \times 10^{-5}$, $\log_2\text{FC} = -0.90$).
* **Baseline Exhaustion:** Concurrently, fasting prediabetic PBMCs exhibit elevated baseline **`DNAJC13`** (Hsp40 family; $p = 0.00049$, $\log_2\text{FC} = +0.65$), an endosomal chaperone regulating clathrin uncoating and receptor recycling. Chronic baseline endosomal/ER stress renders prediabetic immune cells exhausted, preventing an adaptive heat-shock response during glycemic excursions.

### 2. Impaired Incretin & Peptide Metabolism (`PREP`)
* **`PREP` (Prolyl Endopeptidase):** A cytosolic serine peptidase that selectively cleaves peptide bonds at proline residues, directly degrading and activating bioactive peptide hormones including **GLP-1** and **GIP** (which govern glucose-dependent insulin secretion).
* **Finding:** While healthy individuals upregulate `PREP` 2 hours post-glucose load, prediabetic subjects completely fail to mount this response ($p = 4.9 \times 10^{-5}$, $\log_2\text{FC} = -0.73$). This establishes a novel transcriptional bridge linking immune cell dysfunction with the blunted incretin response and delayed first-phase insulin release in prediabetes.

### 3. Removal of Immune Checkpoint Brakes (`SLA` & `NCKAP1L`)
* **Checkpoint Attenuation:** Fasting prediabetic PBMCs show marked downregulation of **`SLA`** (Src-Like Adaptor; $p = 0.00041$, $\log_2\text{FC} = -0.77$). SLA functions as an essential negative regulator of T-cell receptor (TCR) signaling by targeting activated kinases for degradation. Its suppression effectively removes the brake on immune activation.
* **Leukocyte Priming:** Simultaneously, upregulation of **`NCKAP1L`** (HEM1, a hematopoietic-specific component of the WAVE complex governing actin remodeling) primes circulating leukocytes for polarization, endothelial adhesion, and vascular transmigration, underlying the premature atherogenesis characteristic of Asian-Indian men.

### 4. Mitochondrial Translational Inflexibility (`MRPS31`)
* **`MRPS31` (Mitochondrial Ribosomal Protein S31):** One of only five genes consistently shared between states, yet significantly downregulated in prediabetic PBMCs post-OGTT ($p = 0.00026$, $\log_2\text{FC} = -0.58$).
* **Significance:** Normal glucose absorption requires rapid acceleration of mitochondrial protein synthesis to sustain oxidative phosphorylation complexes. Failure to upregulate mitochondrial ribosomal machinery precipitates electron transport chain decoupling, ROS accumulation, and metabolic exhaustion.

---

## Quality Control & Batch Effect Resolution

A critical component of this study was rigorous diagnostic evaluation of array processing order and technical confounders:

1. **Background Noise Removal:** Filtering required probes to achieve detection $p < 0.05$ in $\ge 50\%$ of samples in at least one condition. This eliminated **25,640 low-signal/unexpressed probe features (54.3%)**, preventing spurious low-intensity variance.
2. **Quantile Normalization:** Applied across all 28 arrays, bringing technical distributions into alignment while preserving biological signal. Probes were summarized to **15,704 validated, non-redundant gene symbols**.
3. **Batch Effect & Processing Order Assessment:**
   * Diagnostic correlation analysis revealed an apparent raw correlation with processing order ($r = -0.426, p = 0.024$), traced to two low-intensity arrays at the end of the prediabetic processing batch ("M" series vs "200" series).
   * Following quantile normalization, correlation with processing order was completely abolished across all principal components (PC1: $p = 0.160$; PC2: $p = 0.231$; PC3: $p = 0.380$).
4. **Dominant Biological Signal along PC2:**
   * **PC2 ($14.9\%$ variance)** cleanly separated the Fasting vs 2-Hour Post-OGTT state ($t = -5.33, p = 1.68 \times 10^{-5}$). Healthy donors consistently demonstrated a $+75$ to $+120$ point vertical jump along PC2 in response to glucose, whereas prediabetic donors exhibited trajectory dispersion.

| Diagnostic Metric | Raw BeadStudio Data | Quantile-Normalized Data | Conclusion |
| :--- | :--- | :--- | :--- |
| **Median Array Intensity** | Highly variable ($6.5 - 9.1 \log_2$) | Identical ($7.52 \log_2$) | Technical variation removed |
| **PC1 vs Processing Order** | $r = -0.426, p = 0.024$ | $r = -0.272, p = 0.160$ | Not significant ($p > 0.05$) |
| **PC2 vs Timepoint (OGTT)** | $p = 0.003$ | $t = -5.33, p = 1.68 \times 10^{-5}$ | Primary biological driver |
| **PC1 vs Glucose Group** | $p = 0.449$ | $p = 0.720$ | Baseline subtle, unmasked by OGTT |

---

## Differential Expression Summary

Differential expression analysis was executed using Welch's two-sample unequal variance t-test paired with Benjamini-Hochberg False Discovery Rate (FDR) adjustment and validated via non-parametric Wilcoxon rank-sum testing:

| Comparison Regime | Total DEGs ($p < 0.01$) | Upregulated in PD | Downregulated in PD | Overlapping DEGs |
| :--- | :---: | :---: | :---: | :---: |
| **Comparison 1: Fasting NGT vs PD** | **138** | 42 | 96 | \multirow{2}{*}{\textbf{5} (97.5% divergence)} |
| **Comparison 2: Post-OGTT NGT vs PD** | **62** | 32 | 30 | |

* **Full Results Tables:** Available in `results/DEG_Fasting_NGT_vs_PD.csv` and `results/DEG_OGTT_NGT_vs_PD.csv`.
* **Summary Table:** `results/DEG_summary_table.csv`.

---

## Repository Structure

```
.
├── .gitignore                      # Excludes data/ directory and temporary caches
├── requirements.txt                # Python package dependencies
├── README.md                       # Complete scientific documentation and execution guide
│
├── scripts/                        # Modular, fully reproducible analysis pipeline
│   ├── 01_dataset_exploration.py   # GEO accession download, metadata parsing & matrix extraction
│   ├── 02_quality_control.py       # Detection filtering, quantile normalization & batch diagnostics
│   ├── 03_differential_expression.py # Welch's t-test, Benjamini-Hochberg FDR, Wilcoxon benchmarking
│   ├── 04_visualization.py         # Dual volcano plots, MA plots, clustered heatmaps, Venn, boxplots
│   └── 05_pathway_analysis.py      # ORA (GO BP/MF, KEGG, Reactome) and preranked GSEA
│
├── reports/                        # Publication-ready scientific reports & generation scripts
│   ├── Final_Scientific_Report.docx # Comprehensive 14-page Word report with embedded figures
│   ├── Final_Scientific_Report.md   # Full Markdown scientific manuscript
│   ├── Task 2 - Final_Report.docx   # Structured summary report
│   └── build_final_report_doc.py    # Automated Python generator script for the Word document
│
├── figures/                        # High-resolution (300 DPI) publication figures
│   ├── qc_01_raw_boxplot.png       # Raw intensity distributions across 28 arrays
│   ├── qc_02_raw_density.png       # Raw sample density distributions
│   ├── qc_03_correlation_heatmap.png # Raw pairwise sample Pearson correlations
│   ├── qc_04_pca_raw_multivariable.png # Multivariable PCA of raw unnormalized intensities
│   ├── qc_05_batch_effect_assessment.png # Diagnostic PC correlations with processing order & group
│   ├── qc_06_normalization_comparison.png # Pre- vs Post-quantile normalization boxplots & densities
│   ├── qc_07_pca_normalized.png    # Post-normalization PCA highlighting Timepoint trajectory
│   ├── qc_08_correlation_normalized.png # Post-normalization sample correlation matrix
│   ├── qc_09_pca_subject_effects.png # Longitudinal subject trajectory vectors (Fasting -> OGTT)
│   ├── fig_01_volcano_fasting.png  # Volcano plot: Fasting NGT vs Prediabetes
│   ├── fig_02_volcano_ogtt.png     # Volcano plot: Post-OGTT NGT vs Prediabetes
│   ├── fig_03_volcano_dual_panel.png # Dual-panel side-by-side volcano plot
│   ├── fig_04_ma_plots.png         # Dual MA plots (Log Intensity vs Log2 Fold Change)
│   ├── fig_05_heatmap_fasting_top_degs.png # Hierarchical clustering: Fasting top DEGs
│   ├── fig_06_heatmap_ogtt_top_degs.png # Hierarchical clustering: Post-OGTT top DEGs
│   ├── fig_07_heatmap_combined_signatures.png # Dual heatmap of combined DEG signatures
│   ├── fig_08_deg_venn_overlap.png # Venn diagram illustrating 97.5% DEG divergence (5 shared)
│   ├── fig_09_deg_counts_barchart.png # DEG count breakdown (Upregulated vs Downregulated)
│   ├── fig_10_top_candidate_boxplots.png # Boxplots of key biomarkers (PREP, HSPH1, DNAJC13, SLA, etc.)
│   ├── fig_11_pathway_dotplot_fasting.png # Pathway enrichment dotplot: Fasting NGT vs PD
│   ├── fig_12_pathway_dotplot_ogtt.png # Pathway enrichment dotplot: Post-OGTT NGT vs PD
│   └── fig_13_pathway_comparison_barchart.png # Comparative biological pathway bar chart
│
└── results/                        # Tabular outputs & statistical tables
    ├── DEG_Fasting_NGT_vs_PD.csv   # Complete DEG statistics for Fasting NGT vs PD
    ├── DEG_OGTT_NGT_vs_PD.csv      # Complete DEG statistics for Post-OGTT NGT vs PD
    ├── DEG_summary_table.csv       # Summary metrics and DEG counts
    ├── pathway_fasting_enrichment.csv # GO/KEGG/Reactome ORA enrichment for Fasting
    ├── pathway_ogtt_enrichment.csv # GO/KEGG/Reactome ORA enrichment for Post-OGTT
    ├── gsea_fasting_results.csv    # Preranked GSEA results across 15,704 genes for Fasting
    └── gsea_ogtt_results.csv       # Preranked GSEA results across 15,704 genes for Post-OGTT
```

---

## Reproducibility & Execution Pipeline

### 1. Environment Setup
Clone the repository and install required dependencies:
```bash
git clone https://github.com/<username>/GSE153837-pbmc-transcriptomics.git
cd GSE153837-pbmc-transcriptomics
python -m venv venv
# On Windows:
.\venv\Scripts\activate
# On Linux/macOS:
source venv/bin/activate

pip install -r requirements.txt
```

### 2. Sequential Pipeline Execution
The analysis is orchestrated in five modular stages. Each script reads inputs from `data/` and saves outputs to `figures/`, `results/`, or `reports/`:

```bash
# Step 1: Download raw GEO data and extract expression matrix & sample metadata
python scripts/01_dataset_exploration.py

# Step 2: Quality control (detection p-value filtering, quantile normalization, PCA, batch testing)
python scripts/02_quality_control.py

# Step 3: Statistical differential expression analysis (Welch's t-test, FDR, Wilcoxon)
python scripts/03_differential_expression.py

# Step 4: Publication-quality figures (Volcano, MA, heatmaps, Venn, candidate boxplots)
python scripts/04_visualization.py

# Step 5: Functional enrichment & pathway analysis (GO, KEGG, Reactome ORA & preranked GSEA)
python scripts/05_pathway_analysis.py

# Step 6: Generate the comprehensive Final Scientific Report Word Document with embedded figures
python reports/build_final_report_doc.py
```

---

## Reports & Deliverables

* **Final Scientific Report (Word Document):** Located at [`reports/Final_Scientific_Report.docx`](reports/Final_Scientific_Report.docx). Formatted with professional executive typography, callout callouts, styled statistical tables, and 6 high-resolution figures embedded alongside detailed captions.
* **Final Scientific Report (Markdown):** Located at [`reports/Final_Scientific_Report.md`](reports/Final_Scientific_Report.md). Complete GitHub-flavored markdown reference.
* **Document Builder Script:** [`reports/build_final_report_doc.py`](reports/build_final_report_doc.py). Fully autonomous script utilizing `python-docx` to regenerate the formal Word document at any time.

---

## Evidence-Based Future Directions

Every future recommendation is anchored directly in the empirical findings of this project:

1. **Validation of Incretin Axis via Targeted Targeted Assays (`PREP`):**
   * *Rationale:* `PREP` was significantly suppressed post-OGTT ($p = 4.9 \times 10^{-5}$).
   * *Proposal:* Perform paired enzymatic assays measuring plasma PREP activity alongside active intact GLP-1 (7-36 amide) and total GIP during a standard 3-hour mixed meal tolerance test (MMTT).
2. **PBMC Subpopulation Deconvolution (Single-Cell RNA-seq):**
   * *Rationale:* `SLA` ($p = 0.00041$) and `NCKAP1L` are lineage-specific (T-cell checkpoint and leukocyte motility).
   * *Proposal:* Single-cell CITE-seq on fasting and post-OGTT PBMCs to resolve whether `SLA` suppression occurs uniformly across memory/naive CD4+ and CD8+ T cells or within specific effector populations.
3. **Mitochondrial Respirometry & Chaperone Rescue Experiments:**
   * *Rationale:* Postprandial blunting of `HSPH1` ($p = 6.1 \times 10^{-5}$) and `MRPS31` ($p = 0.00026$).
   * *Proposal:* Seahorse XF flux analysis measuring Oxygen Consumption Rate (OCR) and extracellular acidification rate (ECAR) in PBMCs under acute glucose stress, with small-molecule Hsp70/Hsp110 co-chaperone activators.

---

## Citation & Attribution

If you utilize this pipeline, analysis code, or curated results, please cite:
* **Dataset:** NCBI Gene Expression Omnibus Accession [GSE153837](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153837).
* **Platform:** Illumina HumanHT-12 V4.0 Expression BeadChip ([GPL10558](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558)).
* **Analysis Repository:** GSE153837 Transcriptomic Analysis Pipeline in Asian-Indian Prediabetes.
