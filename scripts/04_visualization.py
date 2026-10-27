"""
=============================================================================
Script 04 — Data Visualization (Task 4)
=============================================================================
Purpose:
    Generate publication-quality visualizations to effectively communicate
    the transcriptomic signatures identified in:
      (i)  Fasting: NGT vs Prediabetes (PD)
      (ii) 2-hour post-OGTT: NGT vs Prediabetes (PD)

Visualizations Generated:
    1. Volcano Plots:
       - fig_01_volcano_fasting.png
       - fig_02_volcano_ogtt.png
       - fig_03_volcano_dual_panel.png
    2. MA Plots:
       - fig_04_ma_plots.png
    3. Hierarchical Clustered Heatmaps:
       - fig_05_heatmap_fasting_top_degs.png
       - fig_06_heatmap_ogtt_top_degs.png
       - fig_07_heatmap_combined_signatures.png
    4. Venn Diagram & DEG Overlap:
       - fig_08_deg_venn_overlap.png
    5. DEG Counts Summary Bar Chart:
       - fig_09_deg_counts_barchart.png
    6. Candidate Biomarker Profile Plots:
       - fig_10_top_candidate_boxplots.png

Outputs:
    - All figures saved into figures/ directory (high-resolution PNGs).
=============================================================================
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2

warnings.filterwarnings("ignore")

# ── Configure Plot Styling ──────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.titlesize": 14,
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "savefig.bbox": "tight"
})

# Color palette
PALETTE = {
    "NGT_Fasting": "#1976D2",  # Deep Blue
    "PD_Fasting":  "#E53935",  # Strong Red
    "NGT_OGTT":     "#388E3C",  # Green
    "PD_OGTT":      "#F57C00",  # Orange
    "Up":           "#D32F2F",  # Red for upregulated in PD
    "Down":         "#1976D2",  # Blue for downregulated in PD
    "NS":           "#BDBDBD"   # Gray for not significant
}

# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIG_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ── Load Data ──────────────────────────────────────────────────────────────
print("=" * 75)
print("LOADING DATA FOR VISUALIZATION")
print("=" * 75)

meta = pd.read_csv(os.path.join(DATA_DIR, "metadata.csv"))
meta["short_label"] = meta["title"].apply(lambda x: x.split("_", 1)[1] if "_" in str(x) else str(x))
expr = pd.read_csv(os.path.join(DATA_DIR, "gene_expression.csv"), index_col=0)
deg_fasting = pd.read_csv(os.path.join(RESULTS_DIR, "DEG_Fasting_NGT_vs_PD.csv"))
deg_ogtt    = pd.read_csv(os.path.join(RESULTS_DIR, "DEG_OGTT_NGT_vs_PD.csv"))

print(f"  Fasting DEG Table : {deg_fasting.shape[0]:,} genes")
print(f"  Post-OGTT DEG Table: {deg_ogtt.shape[0]:,} genes")
print(f"  Expression Matrix : {expr.shape[0]:,} genes × {expr.shape[1]} samples")

# Extract sample sets
samples_ngt_fasting = meta.loc[meta["group"] == "NGT_Fasting", "sample_id"].tolist()
samples_pd_fasting  = meta.loc[meta["group"] == "PD_Fasting", "sample_id"].tolist()
samples_ngt_ogtt    = meta.loc[meta["group"] == "NGT_OGTT", "sample_id"].tolist()
samples_pd_ogtt     = meta.loc[meta["group"] == "PD_OGTT", "sample_id"].tolist()

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 1 & 2: INDIVIDUAL VOLCANO PLOTS
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("GENERATING VOLCANO PLOTS")
print("=" * 75)

def plot_volcano(deg_df, title, filename, top_n_labels=10):
    df = deg_df.copy()
    df["-log10(pvalue)"] = -np.log10(df["pvalue"])
    
    # Define status
    conditions = [
        (df["pvalue"] < 0.05) & (df["log2FC"] >= 0.585),
        (df["pvalue"] < 0.05) & (df["log2FC"] <= -0.585)
    ]
    choices = ["Upregulated in PD", "Downregulated in PD"]
    df["Category"] = np.select(conditions, choices, default="Not Significant")

    fig, ax = plt.subplots(figsize=(8, 6.5))
    
    # Plot points
    sns.scatterplot(
        data=df[df["Category"] == "Not Significant"],
        x="log2FC", y="-log10(pvalue)",
        color=PALETTE["NS"], alpha=0.35, s=20, ax=ax, label="Not Significant"
    )
    sns.scatterplot(
        data=df[df["Category"] == "Downregulated in PD"],
        x="log2FC", y="-log10(pvalue)",
        color=PALETTE["Down"], alpha=0.85, s=40, ax=ax, label=f"Down in PD ({sum(df['Category'] == 'Downregulated in PD')})"
    )
    sns.scatterplot(
        data=df[df["Category"] == "Upregulated in PD"],
        x="log2FC", y="-log10(pvalue)",
        color=PALETTE["Up"], alpha=0.85, s=40, ax=ax, label=f"Up in PD ({sum(df['Category'] == 'Upregulated in PD')})"
    )

    # Threshold guidelines
    ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axvline(0.585, color="black", linestyle=":", linewidth=0.8, alpha=0.7)
    ax.axvline(-0.585, color="black", linestyle=":", linewidth=0.8, alpha=0.7)

    # Label top genes
    sig_genes = df[df["Category"] != "Not Significant"].sort_values("pvalue").head(top_n_labels)
    for _, row in sig_genes.iterrows():
        ax.annotate(
            row["Gene"],
            xy=(row["log2FC"], row["-log10(pvalue)"]),
            xytext=(5, 5), textcoords="offset points",
            fontsize=8, fontweight="bold",
            arrowprops=dict(arrowstyle="-", color="black", lw=0.5, alpha=0.6)
        )

    ax.set_title(title, pad=12, fontweight="bold")
    ax.set_xlabel("log2 Fold Change (PD vs NGT)")
    ax.set_ylabel("-log10 (p-value)")
    ax.legend(frameon=True, facecolor="white", edgecolor="#E0E0E0", loc="upper right")
    ax.grid(True, linestyle="--", alpha=0.3)

    plt.savefig(filename)
    plt.close()
    print(f"  ✓ Saved: {filename}")

plot_volcano(
    deg_fasting,
    title="Differential Expression: Fasting (NGT vs Prediabetes)",
    filename=os.path.join(FIG_DIR, "fig_01_volcano_fasting.png")
)

plot_volcano(
    deg_ogtt,
    title="Differential Expression: 2-Hour Post-OGTT (NGT vs Prediabetes)",
    filename=os.path.join(FIG_DIR, "fig_02_volcano_ogtt.png")
)

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 3: DUAL-PANEL PUBLICATION VOLCANO PLOT
# ══════════════════════════════════════════════════════════════════════════
print("  Generating dual-panel volcano plot...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5), sharey=True)

for ax, df, title, comp_tag in zip([ax1, ax2], [deg_fasting, deg_ogtt],
                                   ["(A) Fasting: NGT vs Prediabetes", "(B) 2-Hour Post-OGTT: NGT vs Prediabetes"],
                                   ["Fasting", "Post-OGTT"]):
    plot_df = df.copy()
    plot_df["-log10(p)"] = -np.log10(plot_df["pvalue"])
    
    # Categories
    conds = [
        (plot_df["pvalue"] < 0.05) & (plot_df["log2FC"] >= 0.585),
        (plot_df["pvalue"] < 0.05) & (plot_df["log2FC"] <= -0.585)
    ]
    plot_df["Cat"] = np.select(conds, ["Up", "Down"], default="NS")
    
    ax.scatter(plot_df.loc[plot_df["Cat"] == "NS", "log2FC"], plot_df.loc[plot_df["Cat"] == "NS", "-log10(p)"],
               c=PALETTE["NS"], alpha=0.3, s=16, label="Not Significant")
    ax.scatter(plot_df.loc[plot_df["Cat"] == "Down", "log2FC"], plot_df.loc[plot_df["Cat"] == "Down", "-log10(p)"],
               c=PALETTE["Down"], alpha=0.85, s=35, label=f"Down in PD ({sum(plot_df['Cat'] == 'Down')})")
    ax.scatter(plot_df.loc[plot_df["Cat"] == "Up", "log2FC"], plot_df.loc[plot_df["Cat"] == "Up", "-log10(p)"],
               c=PALETTE["Up"], alpha=0.85, s=35, label=f"Up in PD ({sum(plot_df['Cat'] == 'Up')})")

    ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axvline(0.585, color="black", linestyle=":", linewidth=0.8, alpha=0.7)
    ax.axvline(-0.585, color="black", linestyle=":", linewidth=0.8, alpha=0.7)

    # Label top 8
    top_candidates = plot_df[plot_df["Cat"] != "NS"].sort_values("pvalue").head(8)
    for _, r in top_candidates.iterrows():
        ax.annotate(r["Gene"], xy=(r["log2FC"], r["-log10(p)"]),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=8, fontweight="bold")

    ax.set_title(title, fontweight="bold", pad=10)
    ax.set_xlabel("log2 Fold Change (Prediabetes vs NGT)")
    if ax == ax1:
        ax.set_ylabel("-log10 (p-value)")
    ax.legend(frameon=True, facecolor="white", edgecolor="#E0E0E0", loc="upper right")
    ax.grid(True, linestyle="--", alpha=0.3)

plt.suptitle("Condition-Specific Differential Expression Signatures in PBMCs", y=1.02, fontweight="bold", fontsize=15)
dual_volcano_path = os.path.join(FIG_DIR, "fig_03_volcano_dual_panel.png")
plt.savefig(dual_volcano_path)
plt.close()
print(f"  ✓ Saved: {dual_volcano_path}")

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 4: MA PLOTS (Mean Expression vs log2FC)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("GENERATING MA PLOTS")
print("=" * 75)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), sharey=True)

for ax, df, title in zip([ax1, ax2], [deg_fasting, deg_ogtt],
                         ["(A) MA Plot: Fasting (NGT vs PD)", "(B) MA Plot: 2-Hour Post-OGTT (NGT vs PD)"]):
    avg_expr = (df["Mean_NGT"] + df["Mean_PD"]) / 2.0
    sig = (df["pvalue"] < 0.05) & (np.abs(df["log2FC"]) >= 0.585)
    
    ax.scatter(avg_expr[~sig], df.loc[~sig, "log2FC"], c=PALETTE["NS"], alpha=0.3, s=15, label="Not Significant")
    ax.scatter(avg_expr[sig], df.loc[sig, "log2FC"], c=PALETTE["Up"], alpha=0.85, s=35, label=f"Significant DEG ({sum(sig)})")
    
    ax.axhline(0, color="black", linestyle="-", linewidth=0.9)
    ax.axhline(0.585, color="black", linestyle=":", linewidth=0.8, alpha=0.7)
    ax.axhline(-0.585, color="black", linestyle=":", linewidth=0.8, alpha=0.7)
    
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Average log2 Expression")
    if ax == ax1:
        ax.set_ylabel("log2 Fold Change (PD vs NGT)")
    ax.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#E0E0E0")
    ax.grid(True, linestyle="--", alpha=0.3)

plt.suptitle("MA Plots: Expression Level vs Fold Change", y=1.02, fontweight="bold", fontsize=14)
ma_path = os.path.join(FIG_DIR, "fig_04_ma_plots.png")
plt.savefig(ma_path)
plt.close()
print(f"  ✓ Saved: {ma_path}")

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 5, 6, 7: HIERARCHICAL CLUSTERED HEATMAPS
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("GENERATING CLUSTERED HEATMAPS")
print("=" * 75)

def plot_heatmap(expr_sub, col_meta, title, filename):
    # Z-score normalize across rows (genes)
    z_expr = expr_sub.apply(lambda x: (x - x.mean()) / (x.std() + 1e-8), axis=1)

    # Color bar for columns
    group_colors = {
        "NGT_Fasting": PALETTE["NGT_Fasting"],
        "PD_Fasting":  PALETTE["PD_Fasting"],
        "NGT_OGTT":     PALETTE["NGT_OGTT"],
        "PD_OGTT":      PALETTE["PD_OGTT"]
    }
    col_colors = col_meta["group"].map(group_colors)
    col_colors.index = col_meta["sample_id"]

    g = sns.clustermap(
        z_expr,
        cmap="vlag",
        center=0,
        vmin=-2.5, vmax=2.5,
        col_colors=col_colors,
        col_cluster=True,
        row_cluster=True,
        linewidths=0.3,
        figsize=(11, 9),
        yticklabels=True,
        xticklabels=col_meta["short_label"].values
    )
    g.fig.suptitle(title, y=1.02, fontweight="bold", fontsize=13)
    g.savefig(filename)
    plt.close("all")
    print(f"  ✓ Saved: {filename}")

# Fasting Heatmap (Top 35 DEGs)
fasting_sig = deg_fasting.loc[(deg_fasting["pvalue"] < 0.05) & (np.abs(deg_fasting["log2FC"]) >= 0.585)].sort_values("pvalue")
top_fasting_genes = fasting_sig["Gene"].head(35).tolist()
meta_fasting = meta[meta["timepoint"] == "Fasting"].sort_values("group")
expr_fasting = expr.loc[top_fasting_genes, meta_fasting["sample_id"]]

plot_heatmap(
    expr_fasting,
    meta_fasting,
    title="Top Differentially Expressed Genes in Fasting State (NGT vs PD)",
    filename=os.path.join(FIG_DIR, "fig_05_heatmap_fasting_top_degs.png")
)

# Post-OGTT Heatmap (Top 35 DEGs)
ogtt_sig = deg_ogtt.loc[(deg_ogtt["pvalue"] < 0.05) & (np.abs(deg_ogtt["log2FC"]) >= 0.585)].sort_values("pvalue")
top_ogtt_genes = ogtt_sig["Gene"].head(35).tolist()
meta_ogtt = meta[meta["timepoint"] == "OGTT"].sort_values("group")
expr_ogtt = expr.loc[top_ogtt_genes, meta_ogtt["sample_id"]]

plot_heatmap(
    expr_ogtt,
    meta_ogtt,
    title="Top Differentially Expressed Genes in 2-Hour Post-OGTT State (NGT vs PD)",
    filename=os.path.join(FIG_DIR, "fig_06_heatmap_ogtt_top_degs.png")
)

# Combined Signature Heatmap (All 4 conditions for prioritized genes)
prioritized_genes = list(dict.fromkeys(top_fasting_genes[:20] + top_ogtt_genes[:20]))
meta_all = meta.sort_values(["timepoint", "condition"])
expr_combined = expr.loc[prioritized_genes, meta_all["sample_id"]]

plot_heatmap(
    expr_combined,
    meta_all,
    title="Dynamic PBMC Transcriptomic Signatures Across Fasting & Post-OGTT",
    filename=os.path.join(FIG_DIR, "fig_07_heatmap_combined_signatures.png")
)

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 8: VENN DIAGRAM & DEG OVERLAP
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("GENERATING VENN DIAGRAM")
print("=" * 75)

fasting_genes_set = set(fasting_sig["Gene"])
ogtt_genes_set    = set(ogtt_sig["Gene"])

shared_genes = fasting_genes_set.intersection(ogtt_genes_set)
fasting_only = fasting_genes_set - shared_genes
ogtt_only    = ogtt_genes_set - shared_genes

print(f"  Fasting DEGs Total       : {len(fasting_genes_set)}")
print(f"  Post-OGTT DEGs Total     : {len(ogtt_genes_set)}")
print(f"  Shared across conditions : {len(shared_genes)}")
print(f"  Fasting-unique DEGs      : {len(fasting_only)}")
print(f"  Post-OGTT-unique DEGs    : {len(ogtt_only)}")
if len(shared_genes) > 0:
    print(f"  Shared genes list        : {list(shared_genes)}")

fig, ax = plt.subplots(figsize=(7, 6))
v = venn2(
    subsets=(len(fasting_only), len(ogtt_only), len(shared_genes)),
    set_labels=("Fasting DEGs\n(NGT vs PD)", "Post-OGTT DEGs\n(NGT vs PD)"),
    set_colors=("#1976D2", "#F57C00"),
    alpha=0.6,
    ax=ax
)

# Styling labels
for text in v.set_labels:
    text.set_fontsize(11)
    text.set_fontweight("bold")
for text in v.subset_labels:
    if text:
        text.set_fontsize(13)
        text.set_fontweight("bold")

ax.set_title("Overlap of Differentially Expressed Genes\n(Fasting vs 2-Hour Post-OGTT)", pad=15, fontweight="bold")
venn_path = os.path.join(FIG_DIR, "fig_08_deg_venn_overlap.png")
plt.savefig(venn_path)
plt.close()
print(f"  ✓ Saved: {venn_path}")

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 9: DEG COUNTS SUMMARY BAR CHART
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("GENERATING DEG COUNTS BAR CHART")
print("=" * 75)

counts_data = pd.DataFrame({
    "Comparison": ["Fasting (NGT vs PD)", "Fasting (NGT vs PD)", "Post-OGTT (NGT vs PD)", "Post-OGTT (NGT vs PD)"],
    "Direction": ["Upregulated in PD", "Downregulated in PD", "Upregulated in PD", "Downregulated in PD"],
    "Count": [
        sum((deg_fasting["pvalue"] < 0.05) & (deg_fasting["log2FC"] >= 0.585)),
        sum((deg_fasting["pvalue"] < 0.05) & (deg_fasting["log2FC"] <= -0.585)),
        sum((deg_ogtt["pvalue"] < 0.05) & (deg_ogtt["log2FC"] >= 0.585)),
        sum((deg_ogtt["pvalue"] < 0.05) & (deg_ogtt["log2FC"] <= -0.585)),
    ]
})

fig, ax = plt.subplots(figsize=(8, 5.5))
bars = sns.barplot(
    data=counts_data,
    x="Comparison", y="Count", hue="Direction",
    palette={"Upregulated in PD": PALETTE["Up"], "Downregulated in PD": PALETTE["Down"]},
    ax=ax, edgecolor="black", linewidth=0.7
)

for p in ax.patches:
    h = p.get_height()
    if h > 0:
        ax.annotate(f"{int(h)}", (p.get_x() + p.get_width() / 2., h),
                    ha="center", va="bottom", fontsize=10, fontweight="bold",
                    xytext=(0, 3), textcoords="offset points")

ax.set_title("Differentially Expressed Genes by Direction (p < 0.05, |log2FC| >= 0.585)", pad=12, fontweight="bold")
ax.set_ylabel("Number of Genes")
ax.set_xlabel("")
ax.set_ylim(0, max(counts_data["Count"]) * 1.15)
ax.legend(title="Direction of Regulation", frameon=True, facecolor="white", edgecolor="#E0E0E0")
ax.grid(axis="y", linestyle="--", alpha=0.3)

barchart_path = os.path.join(FIG_DIR, "fig_09_deg_counts_barchart.png")
plt.savefig(barchart_path)
plt.close()
print(f"  ✓ Saved: {barchart_path}")

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 10: CANDIDATE BIOMARKER PROFILE BOXPLOTS
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("GENERATING CANDIDATE BIOMARKER BOXPLOTS")
print("=" * 75)

candidate_genes = ["PREP", "HSPH1", "DNAJC13", "SLA", "NCKAP1L", "MRPS31"]
fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True)
axes = axes.flatten()

group_order = ["NGT_Fasting", "PD_Fasting", "NGT_OGTT", "PD_OGTT"]
group_palette = {
    "NGT_Fasting": PALETTE["NGT_Fasting"],
    "PD_Fasting":  PALETTE["PD_Fasting"],
    "NGT_OGTT":     PALETTE["NGT_OGTT"],
    "PD_OGTT":      PALETTE["PD_OGTT"]
}

for i, gene in enumerate(candidate_genes):
    ax = axes[i]
    if gene in expr.index:
        plot_df = pd.DataFrame({
            "Expression": expr.loc[gene].values,
            "Group": meta["group"].values
        })
        
        # Boxplot with overlay strip
        sns.boxplot(
            data=plot_df, x="Group", y="Expression", order=group_order,
            palette=group_palette, ax=ax, width=0.5, boxprops=dict(alpha=0.75),
            fliersize=0
        )
        sns.stripplot(
            data=plot_df, x="Group", y="Expression", order=group_order,
            color="black", size=5, alpha=0.6, jitter=0.2, ax=ax
        )
        
        ax.set_title(f"{gene}", fontweight="bold", fontsize=12)
        ax.set_ylabel("Normalized log2 Expression")
        ax.set_xlabel("")
        ax.grid(axis="y", linestyle="--", alpha=0.3)
    else:
        ax.text(0.5, 0.5, f"{gene} not found", ha="center", va="center")

for ax in axes[-3:]:
    ax.set_xticklabels(["NGT\nFasting", "PD\nFasting", "NGT\nOGTT", "PD\nOGTT"], rotation=0)

plt.suptitle("Expression Dynamics of Key Prioritized Biomarker Genes Across Conditions", y=1.02, fontweight="bold", fontsize=14)
plt.tight_layout()
candidate_path = os.path.join(FIG_DIR, "fig_10_top_candidate_boxplots.png")
plt.savefig(candidate_path)
plt.close()
print(f"  ✓ Saved: {candidate_path}")

print("\n" + "=" * 75)
print("Script 04 completed successfully! All publication figures generated.")
print("=" * 75)
