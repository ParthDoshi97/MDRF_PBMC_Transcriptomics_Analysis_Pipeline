"""
=============================================================================
Script 03 — Differential Expression Analysis (Task 3)
=============================================================================
Purpose:
    Perform condition-specific differential expression analysis for the two
    study comparisons:
      (i)  Fasting: NGT vs Prediabetes (PD)
      (ii) 2-hour post-OGTT: NGT vs Prediabetes (PD)

Methodological Justification:
    - Platform: Illumina HumanHT-12 v4.0 BeadChip Microarray.
    - Continuous fluorescence intensities follow a log-normal distribution,
      NOT discrete read counts.
    - DESeq2 and edgeR model discrete integer counts using the Negative Binomial
      distribution (appropriate for RNA-Seq, invalid for microarray intensities).
    - Gold standard for microarrays: Linear modeling / Welch's t-test with
      Benjamini-Hochberg (BH) False Discovery Rate (FDR) correction on
      quantile-normalized log2 data.
    - Non-parametric Mann-Whitney U (Wilcoxon rank-sum) is also calculated
      for sensitivity benchmarking.

Outputs:
    - results/DEG_Fasting_NGT_vs_PD.csv
    - results/DEG_OGTT_NGT_vs_PD.csv
    - results/DEG_summary_table.csv
=============================================================================
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── Load Data ──────────────────────────────────────────────────────────────
print("=" * 75)
print("LOADING DATA FOR DIFFERENTIAL EXPRESSION ANALYSIS")
print("=" * 75)

meta = pd.read_csv(os.path.join(DATA_DIR, "metadata.csv"))
expr = pd.read_csv(os.path.join(DATA_DIR, "gene_expression.csv"), index_col=0)

print(f"  Preprocessed Expression Matrix: {expr.shape[0]:,} genes × {expr.shape[1]} samples")
print(f"  Metadata: {meta.shape[0]} samples")

# Map samples to groups
samples_ngt_fasting = meta.loc[meta["group"] == "NGT_Fasting", "sample_id"].tolist()
samples_pd_fasting  = meta.loc[meta["group"] == "PD_Fasting", "sample_id"].tolist()

samples_ngt_ogtt    = meta.loc[meta["group"] == "NGT_OGTT", "sample_id"].tolist()
samples_pd_ogtt     = meta.loc[meta["group"] == "PD_OGTT", "sample_id"].tolist()

print("\n  Sample breakdown for the two independent comparisons:")
print(f"  Comparison 1 [Fasting]:   NGT (n={len(samples_ngt_fasting)}) vs Prediabetes (n={len(samples_pd_fasting)})")
print(f"  Comparison 2 [Post-OGTT]: NGT (n={len(samples_ngt_ogtt)}) vs Prediabetes (n={len(samples_pd_ogtt)})")

# ── Differential Expression Function ──────────────────────────────────────
def run_differential_expression(expr_df, control_samples, case_samples, comparison_name):
    """
    Computes differential expression metrics between control and case groups.
    log2FC = Mean(Case) - Mean(Control) (since data is in log2 space).
    Positive log2FC = Upregulated in Prediabetes (PD).
    Negative log2FC = Downregulated in Prediabetes (PD).
    """
    print(f"\n" + "-" * 75)
    print(f"Running DE Analysis: {comparison_name}")
    print(f"  Controls ({len(control_samples)}): {control_samples}")
    print(f"  Cases    ({len(case_samples)}): {case_samples}")
    print("-" * 75)

    ctrl_mat = expr_df[control_samples].values
    case_mat = expr_df[case_samples].values

    # Mean expression per group
    mean_ctrl = np.mean(ctrl_mat, axis=1)
    mean_case = np.mean(case_mat, axis=1)
    
    # Standard deviation per group
    sd_ctrl = np.std(ctrl_mat, axis=1, ddof=1)
    sd_case = np.std(case_mat, axis=1, ddof=1)

    # Log2 Fold Change: Case (PD) - Control (NGT)
    log2fc = mean_case - mean_ctrl

    # Welch's t-test (unequal variance)
    t_stat, p_vals = stats.ttest_ind(case_mat, ctrl_mat, axis=1, equal_var=False)

    # Non-parametric Mann-Whitney U test (Wilcoxon rank-sum)
    u_pvals = []
    for i in range(len(expr_df)):
        try:
            _, u_p = stats.mannwhitneyu(case_mat[i], ctrl_mat[i], alternative="two-sided")
        except Exception:
            u_p = 1.0
        u_pvals.append(u_p)
    u_pvals = np.array(u_pvals)

    # Multiple testing correction using Benjamini-Hochberg (FDR)
    # Handle any potential NaNs in p-values
    valid_p = ~np.isnan(p_vals)
    adj_pvals = np.ones_like(p_vals)
    if np.sum(valid_p) > 0:
        _, p_adj_valid, _, _ = multipletests(p_vals[valid_p], method="fdr_bh")
        adj_pvals[valid_p] = p_adj_valid

    _, u_adj_pvals, _, _ = multipletests(u_pvals, method="fdr_bh")

    # Linear Fold Change for intuitive biological interpretation
    # If log2FC >= 0: FC = 2^(log2FC)
    # If log2FC < 0: FC = - 2^(-log2FC)
    linear_fc = np.where(log2fc >= 0, np.power(2.0, log2fc), -np.power(2.0, -log2fc))

    results = pd.DataFrame({
        "Gene": expr_df.index,
        "Mean_NGT": mean_ctrl,
        "Mean_PD": mean_case,
        "SD_NGT": sd_ctrl,
        "SD_PD": sd_case,
        "log2FC": log2fc,
        "Fold_Change": linear_fc,
        "t_statistic": t_stat,
        "pvalue": p_vals,
        "padj": adj_pvals,
        "wilcoxon_pvalue": u_pvals,
        "wilcoxon_padj": u_adj_pvals
    })

    # Add categorization tags
    # 1. Strict Cutoff: FDR < 0.05 and |log2FC| >= 1.0 (2-fold)
    results["Status_Strict"] = "Not Significant"
    results.loc[(results["padj"] < 0.05) & (results["log2FC"] >= 1.0), "Status_Strict"] = "Upregulated in PD"
    results.loc[(results["padj"] < 0.05) & (results["log2FC"] <= -1.0), "Status_Strict"] = "Downregulated in PD"

    # 2. Standard Cutoff: FDR < 0.05 and |log2FC| >= 0.585 (1.5-fold)
    results["Status_Standard"] = "Not Significant"
    results.loc[(results["padj"] < 0.05) & (results["log2FC"] >= 0.585), "Status_Standard"] = "Upregulated in PD"
    results.loc[(results["padj"] < 0.05) & (results["log2FC"] <= -0.585), "Status_Standard"] = "Downregulated in PD"

    # 3. Nominal Cutoff: pvalue < 0.05 and |log2FC| >= 0.585 (1.5-fold)
    results["Status_Nominal"] = "Not Significant"
    results.loc[(results["pvalue"] < 0.05) & (results["log2FC"] >= 0.585), "Status_Nominal"] = "Upregulated in PD"
    results.loc[(results["pvalue"] < 0.05) & (results["log2FC"] <= -0.585), "Status_Nominal"] = "Downregulated in PD"

    # Sort by significance (p-value ascending)
    results = results.sort_values("pvalue", ascending=True).reset_index(drop=True)
    return results

# ── Execute Comparison 1: Fasting NGT vs Fasting PD ────────────────────────
deg_fasting = run_differential_expression(
    expr_df=expr,
    control_samples=samples_ngt_fasting,
    case_samples=samples_pd_fasting,
    comparison_name="Comparison 1: Fasting NGT vs Fasting Prediabetes (PD)"
)

# ── Execute Comparison 2: OGTT NGT vs OGTT PD ──────────────────────────────
deg_ogtt = run_differential_expression(
    expr_df=expr,
    control_samples=samples_ngt_ogtt,
    case_samples=samples_pd_ogtt,
    comparison_name="Comparison 2: 2-Hour Post-OGTT NGT vs 2-Hour Post-OGTT Prediabetes (PD)"
)

# ── Summary Tables & Inspection ───────────────────────────────────────────
print("\n" + "=" * 75)
print("DIFFERENTIAL EXPRESSION SUMMARY COUNTS")
print("=" * 75)

summary_rows = []
for name, df in [("Fasting (NGT vs PD)", deg_fasting), ("Post-OGTT (NGT vs PD)", deg_ogtt)]:
    row = {
        "Comparison": name,
        "Total_Genes_Tested": len(df),
        "Nominal_Up (p<0.05, log2FC>=0.585)": (df["Status_Nominal"] == "Upregulated in PD").sum(),
        "Nominal_Down (p<0.05, log2FC<=-0.585)": (df["Status_Nominal"] == "Downregulated in PD").sum(),
        "Nominal_Total": (df["Status_Nominal"] != "Not Significant").sum(),
        "FDR_Up (padj<0.05, log2FC>=0.585)": (df["Status_Standard"] == "Upregulated in PD").sum(),
        "FDR_Down (padj<0.05, log2FC<=-0.585)": (df["Status_Standard"] == "Downregulated in PD").sum(),
        "FDR_Total": (df["Status_Standard"] != "Not Significant").sum(),
        "Strict_FDR_Up (padj<0.05, |log2FC|>=1.0)": (df["Status_Strict"] == "Upregulated in PD").sum(),
        "Strict_FDR_Down (padj<0.05, |log2FC|>=1.0)": (df["Status_Strict"] == "Downregulated in PD").sum(),
    }
    summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
print(summary_df.T.to_string())

# ── Print Top 10 Significant Genes for Both Comparisons ───────────────────
print("\n" + "=" * 75)
print("TOP 10 DIFFERENTIALLY EXPRESSED GENES (FASTING: NGT vs PD)")
print("=" * 75)
top_cols = ["Gene", "Mean_NGT", "Mean_PD", "log2FC", "Fold_Change", "pvalue", "padj", "Status_Nominal"]
print(deg_fasting[top_cols].head(10).to_string(index=False))

print("\n" + "=" * 75)
print("TOP 10 DIFFERENTIALLY EXPRESSED GENES (POST-OGTT: NGT vs PD)")
print("=" * 75)
print(deg_ogtt[top_cols].head(10).to_string(index=False))

# ── Save Outputs ──────────────────────────────────────────────────────────
print("\n" + "=" * 75)
print("SAVING DEG RESULTS TO CSV")
print("=" * 75)

fasting_out = os.path.join(RESULTS_DIR, "DEG_Fasting_NGT_vs_PD.csv")
ogtt_out    = os.path.join(RESULTS_DIR, "DEG_OGTT_NGT_vs_PD.csv")
sum_out     = os.path.join(RESULTS_DIR, "DEG_summary_table.csv")

deg_fasting.to_csv(fasting_out, index=False)
deg_ogtt.to_csv(ogtt_out, index=False)
summary_df.to_csv(sum_out, index=False)

print(f"  ✓ Saved Fasting DEGs   : {fasting_out}")
print(f"  ✓ Saved Post-OGTT DEGs : {ogtt_out}")
print(f"  ✓ Saved Summary Table  : {sum_out}")

print("\n" + "=" * 75)
print("Script 03 completed successfully!")
print("=" * 75)
