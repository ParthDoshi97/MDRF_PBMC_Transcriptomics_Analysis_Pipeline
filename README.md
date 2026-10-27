# MDRF Use Case 2: PBMC Transcriptomics Analysis Pipeline

This repository contains the end-to-end Python bioinformatics pipeline for analyzing **PBMC Gene Expression Profiling Under Acute Oral Glucose Challenge (75g OGTT) in Asian-Indian Men with Prediabetes** from NCBI GEO:
* **`GSE153837`**: Asian-Indian Cohort ($n = 28$ samples across $4$ experimental groups)
* **Platform:** Illumina HumanHT-12 V4.0 Expression BeadChip (`GPL10558`; $47,311$ probe features)
* **Original Study:** Mallu et al., *Clinical Nutrition* (2021) | **Cohorts:** Normal Glucose Tolerance (NGT) vs. Prediabetes (PD) in Fasting ($0\text{ hr}$) and Post-OGTT ($2\text{ hr}$) states

---

## Folder Structure

```text
MDRF_PBMC_Transcriptomics_Analysis_Pipeline/
├── README.md                      # Quick-start guide and pipeline overview
├── requirements.txt               # Required Python libraries
├── .gitignore                     # Excludes raw data/ (~50 MB) and temporary caches
├── scripts/                       # Step-by-step analysis scripts (Tasks 1–5)
│   ├── 01_dataset_exploration.py     # Task 1: Downloads GSE153837, parses metadata, extracts expression
│   ├── 02_quality_control.py         # Task 2: Detection p-value filtering, quantile normalization, PCA & batch diagnostics
│   ├── 03_differential_expression.py # Task 3: Welch's t-test, Benjamini-Hochberg FDR, Wilcoxon benchmarking
│   ├── 04_visualization.py           # Task 4: Volcano, MA plots, clustered heatmaps, Venn, & candidate boxplots
│   └── 05_pathway_analysis.py        # Task 5: Functional enrichment (GO, KEGG, Reactome) & preranked GSEA
├── reports/                       # Scientific research reports & generation scripts
│   ├── Comprehensive_Scientific_Report_GSE153837.docx # Publication-style research article report (.docx)
│   ├── Final_Scientific_Report.docx  # Master Word scientific report
│   ├── Final_Scientific_Report.md    # Markdown version of the scientific report
│   └── build_academic_report_doc.py  # Script to generate clean academic Word report
├── figures/                       # Generated 300 DPI publication figures (QC & Tasks 1–5)
│   ├── qc_05_batch_effect_assessment.png
│   ├── qc_07_pca_normalized.png
│   ├── fig_03_volcano_dual_panel.png
│   ├── fig_07_heatmap_combined_signatures.png
│   ├── fig_08_deg_venn_overlap.png
│   ├── fig_10_top_candidate_boxplots.png
│   └── fig_13_pathway_comparison_barchart.png
└── results/                       # Output summary tables (CSV) from Tasks 1–5
    ├── DEG_Fasting_NGT_vs_PD.csv
    ├── DEG_OGTT_NGT_vs_PD.csv
    ├── DEG_summary_table.csv
    ├── pathway_fasting_enrichment.csv
    ├── pathway_ogtt_enrichment.csv
    ├── gsea_fasting_results.csv
    └── gsea_ogtt_results.csv
```

---

## How to Run the Pipeline

### 1. Install Dependencies
```bash
pip install -r requirements.txt
```

### 2. Run Each Step in Order
```bash
# Task 1: Dataset Exploration & Raw Data Parsing (outputs data/metadata.csv & raw_expression.csv)
python scripts/01_dataset_exploration.py

# Task 2: Quality Control, Quantile Normalization & Batch Assessment (outputs data/gene_expression.csv & figures/qc_*.png)
python scripts/02_quality_control.py

# Task 3: Differential Expression Analysis (outputs results/DEG_*.csv)
python scripts/03_differential_expression.py

# Task 4: Generate Publication Figures (outputs figures/fig_01 through fig_10)
python scripts/04_visualization.py

# Task 5: Functional Enrichment & GSEA (outputs results/pathway_*.csv, results/gsea_*.csv, & figures/fig_11–13)
python scripts/05_pathway_analysis.py

# Optional: Generate Clean Academic Word Report (outputs reports/Comprehensive_Scientific_Report_GSE153837.docx)
python reports/build_academic_report_doc.py
```

---

## Summary of Key Outputs

| Task | Script | Key Output File(s) | Main Result |
| :--- | :--- | :--- | :--- |
| **Task 1** | `01_dataset_exploration.py` | `data/metadata.csv`, `data/raw_expression.csv` | 28 samples, Asian-Indian males, pre-matched for age & BMI across 4 experimental groups (NGT/PD × Fasting/OGTT). |
| **Task 2** | `02_quality_control.py` | `data/gene_expression.csv`, `figures/qc_05_batch_effect_assessment.png`, `figures/qc_07_pca_normalized.png` | Filtered 25,713 unexpressed probes (54.3%); quantile normalized to 15,704 genes; resolved batch order artifacts ($p > 0.05$); Timepoint identified as dominant biological driver on PC2 ($p = 1.68 \times 10^{-5}$). |
| **Task 3** | `03_differential_expression.py` | `results/DEG_Fasting_NGT_vs_PD.csv`, `results/DEG_OGTT_NGT_vs_PD.csv`, `results/DEG_summary_table.csv` | Identified 138 Fasting DEGs (42 up, 96 down) and 62 Post-OGTT DEGs (32 up, 30 down); validated via Welch's t-test and Wilcoxon benchmarking. |
| **Task 4** | `04_visualization.py` | `figures/fig_01`–`fig_10` (`.png`) | Volcano plots, MA plots, clustered heatmaps, Venn diagram (only 5 shared genes / 97.5% state divergence: the acute glucose stress test), and candidate boxplots (*PREP*, *HSPH1*, *DNAJC13*, *SLA*, *NCKAP1L*, *MRPS31*). |
| **Task 5** | `05_pathway_analysis.py` | `results/pathway_*.csv`, `results/gsea_*.csv`, `figures/fig_11`–`fig_13` | Disrupted biological pathways: Chaperone Crisis (*HSPH1* / *DNAJC13*), Incretin Metabolism (*PREP*), Immune Checkpoint Attenuation (*SLA* / *NCKAP1L*), and Mitochondrial Translation (*MRPS31*). |
