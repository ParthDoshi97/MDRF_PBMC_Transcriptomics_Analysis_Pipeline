"""
=============================================================================
Script 02 — Data Quality Assessment & Preprocessing (Task 2)
=============================================================================
Purpose:
    Assess data quality, check for batch effects, perform PCA to understand
    variability, preprocess the data (normalization, probe filtering,
    probe-to-gene mapping).

Inputs:
    - data/metadata.csv
    - data/raw_expression.csv
    - GEO soft file (for detection p-values and platform annotation)

Outputs:
    - data/normalized_expression.csv   : Quantile-normalized probe-level data
    - data/gene_expression.csv         : Gene-level collapsed expression
    - data/detection_pvals.csv         : Detection p-value matrix
    - figures/qc_*.png                 : All QC figures
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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy import stats
import GEOparse

warnings.filterwarnings("ignore")
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "font.size": 10,
    "figure.figsize": (12, 8),
    "savefig.bbox": "tight"
})

# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR  = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────────
print("=" * 70)
print("LOADING DATA")
print("=" * 70)

meta = pd.read_csv(os.path.join(DATA_DIR, "metadata.csv"))
expr = pd.read_csv(os.path.join(DATA_DIR, "raw_expression.csv"), index_col=0)
expr = expr.apply(pd.to_numeric, errors="coerce")

print(f"  Expression matrix : {expr.shape[0]:,} probes × {expr.shape[1]} samples")
print(f"  Metadata          : {meta.shape[0]} samples × {meta.shape[1]} columns")

# Ensure sample order matches
samples = [s for s in expr.columns if s in meta["sample_id"].values]
expr = expr[samples]
meta = meta[meta["sample_id"].isin(samples)].set_index("sample_id").loc[samples].reset_index()

# Short labels for plots
meta["short_label"] = meta["title"].apply(lambda x: x.split("_", 1)[1] if "_" in str(x) else str(x))

print(f"  Matched samples   : {len(samples)}")

# ── Load detection p-values from GEO soft file ─────────────────────────────
print("\n" + "=" * 70)
print("EXTRACTING DETECTION P-VALUES")
print("=" * 70)

gse = GEOparse.get_GEO(geo="GSE153837", destdir=DATA_DIR, silent=True)

det_pval_dict = {}
for gsm_name, gsm in gse.gsms.items():
    tbl = gsm.table
    if "Detection Pval" in tbl.columns:
        det_pval_dict[gsm_name] = tbl.set_index("ID_REF")["Detection Pval"].astype(float)

det_pvals = pd.DataFrame(det_pval_dict)
det_pvals = det_pvals.loc[expr.index, expr.columns]
det_pvals.to_csv(os.path.join(DATA_DIR, "detection_pvals.csv"))
print(f"  ✓ Detection p-values extracted: {det_pvals.shape}")

# ══════════════════════════════════════════════════════════════════════════
# PART A: RAW DATA QUALITY ASSESSMENT
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART A: RAW DATA QUALITY ASSESSMENT")
print("=" * 70)

# ── A1. Box plot — raw expression distributions ──────────────────────────
print("\n  [A1] Box plot of raw expression distributions...")

fig, ax = plt.subplots(figsize=(16, 7))
colors = {"NGT_Fasting": "#2196F3", "NGT_OGTT": "#4CAF50",
          "PD_Fasting": "#F44336", "PD_OGTT": "#FF9800"}
sample_colors = [colors.get(g, "gray") for g in meta["group"]]

bp = ax.boxplot([expr[s].dropna().values for s in expr.columns],
                patch_artist=True, showfliers=False, widths=0.6)
for patch, color in zip(bp["boxes"], sample_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax.set_xticklabels(meta["short_label"].values, rotation=90, fontsize=7)
ax.set_ylabel("Expression (log2)")
ax.set_title("Raw Expression Distributions per Sample")

# Legend
from matplotlib.patches import Patch
legend_items = [Patch(facecolor=c, label=g, alpha=0.7) for g, c in colors.items()]
ax.legend(handles=legend_items, loc="upper right", fontsize=8)

plt.savefig(os.path.join(FIG_DIR, "qc_01_raw_boxplot.png"))
plt.close()
print("    ✓ Saved: figures/qc_01_raw_boxplot.png")

# ── A2. Density plots — raw expression ───────────────────────────────────
print("  [A2] Density plot of raw expression distributions...")

fig, ax = plt.subplots(figsize=(14, 7))
for i, s in enumerate(expr.columns):
    group = meta.loc[meta["sample_id"] == s, "group"].values[0]
    ax.hist(expr[s].dropna().values, bins=100, alpha=0.3, density=True,
            color=colors.get(group, "gray"), label=group if i < 4 else "")

ax.set_xlabel("Expression (log2)")
ax.set_ylabel("Density")
ax.set_title("Raw Expression Density Distribution")
ax.legend(handles=legend_items, fontsize=8)
plt.savefig(os.path.join(FIG_DIR, "qc_02_raw_density.png"))
plt.close()
print("    ✓ Saved: figures/qc_02_raw_density.png")

# ── A3. Sample-sample correlation heatmap ────────────────────────────────
print("  [A3] Sample correlation heatmap...")

corr = expr.corr(method="pearson")
fig, ax = plt.subplots(figsize=(14, 12))

# Row/col colors for annotation
group_lut = {g: c for g, c in colors.items()}
row_colors = meta.set_index("sample_id")["group"].map(group_lut)

g = sns.clustermap(corr, cmap="RdYlBu_r", vmin=corr.min().min(), vmax=1,
                   figsize=(14, 12),
                   row_colors=row_colors, col_colors=row_colors,
                   xticklabels=meta["short_label"].values,
                   yticklabels=meta["short_label"].values,
                   linewidths=0.5)
g.fig.suptitle("Sample-Sample Pearson Correlation (Raw Data)", y=1.02, fontsize=14)
g.savefig(os.path.join(FIG_DIR, "qc_03_correlation_heatmap.png"))
plt.close("all")
print("    ✓ Saved: figures/qc_03_correlation_heatmap.png")

# ── A4. PCA — raw data, multiple variable colorings ─────────────────────
print("  [A4] PCA on raw data (checking variability sources)...")

# Transpose: samples as rows, probes as features
X_raw = expr.T.values
X_scaled = StandardScaler().fit_transform(X_raw)

pca_raw = PCA(n_components=min(10, len(samples)))
pcs_raw = pca_raw.fit_transform(X_scaled)
var_explained = pca_raw.explained_variance_ratio_ * 100

pca_df = pd.DataFrame({
    "PC1": pcs_raw[:, 0], "PC2": pcs_raw[:, 1], "PC3": pcs_raw[:, 2],
    "sample_id": expr.columns
})
pca_df = pca_df.merge(meta[["sample_id", "condition", "timepoint", "group"]], on="sample_id")

# Plot PCA colored by different variables
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# PC1 vs PC2 — by Group
ax = axes[0, 0]
for grp, color in colors.items():
    mask = pca_df["group"] == grp
    ax.scatter(pca_df.loc[mask, "PC1"], pca_df.loc[mask, "PC2"],
               c=color, label=grp, s=80, edgecolors="black", linewidths=0.5, alpha=0.8)
ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
ax.set_title("PCA: Colored by Group (Condition × Timepoint)")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# PC1 vs PC2 — by Condition only
ax = axes[0, 1]
cond_colors = {"NGT": "#2196F3", "PD": "#F44336"}
for cond, color in cond_colors.items():
    mask = pca_df["condition"] == cond
    ax.scatter(pca_df.loc[mask, "PC1"], pca_df.loc[mask, "PC2"],
               c=color, label=cond, s=80, edgecolors="black", linewidths=0.5, alpha=0.8)
ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
ax.set_title("PCA: Colored by Condition (NGT vs PD)")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# PC1 vs PC2 — by Timepoint only
ax = axes[1, 0]
tp_colors = {"Fasting": "#9C27B0", "OGTT": "#FF9800"}
for tp, color in tp_colors.items():
    mask = pca_df["timepoint"] == tp
    ax.scatter(pca_df.loc[mask, "PC1"], pca_df.loc[mask, "PC2"],
               c=color, label=tp, s=80, edgecolors="black", linewidths=0.5, alpha=0.8)
ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
ax.set_title("PCA: Colored by Timepoint (Fasting vs OGTT)")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Scree plot — variance explained
ax = axes[1, 1]
ax.bar(range(1, len(var_explained) + 1), var_explained, color="#607D8B", alpha=0.7)
ax.plot(range(1, len(var_explained) + 1), np.cumsum(var_explained), "ro-", markersize=6)
ax.set_xlabel("Principal Component")
ax.set_ylabel("Variance Explained (%)")
ax.set_title("Scree Plot — Variance Explained per PC")
ax.set_xticks(range(1, len(var_explained) + 1))
ax.grid(True, alpha=0.3)

# Add cumulative % labels
for i, cum in enumerate(np.cumsum(var_explained)):
    ax.annotate(f"{cum:.0f}%", (i + 1, cum), textcoords="offset points",
                xytext=(0, 8), ha="center", fontsize=7, color="red")

plt.suptitle("PCA Analysis — Raw Data (Before Preprocessing)", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "qc_04_pca_raw_multivariable.png"))
plt.close()
print("    ✓ Saved: figures/qc_04_pca_raw_multivariable.png")

# ── A5. Batch effect assessment ──────────────────────────────────────────
print("  [A5] Batch effect assessment...")

# Extract potential batch variables from sample titles
# Sample IDs often encode processing batch info
meta["sample_number"] = meta["title"].apply(lambda x: str(x).split("_")[0] if "_" in str(x) else "")

# Check if sample ordering correlates with PC1/PC2
# (if samples were processed in order, this could indicate a batch effect)
pca_df_merged = pca_df.merge(meta[["sample_id", "sample_number"]], on="sample_id")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# PC1 vs sample order
ax = axes[0]
sample_order = range(len(pca_df_merged))
scatter = ax.scatter(sample_order, pca_df_merged["PC1"],
                     c=[colors.get(g, "gray") for g in pca_df_merged["group"]],
                     s=80, edgecolors="black", linewidths=0.5)
r_val, p_val = stats.pearsonr(list(sample_order), pca_df_merged["PC1"])
ax.set_xlabel("Sample Processing Order")
ax.set_ylabel("PC1")
ax.set_title(f"PC1 vs Processing Order\n(r={r_val:.3f}, p={p_val:.3f})")
ax.grid(True, alpha=0.3)

# PC2 vs sample order
ax = axes[1]
scatter = ax.scatter(sample_order, pca_df_merged["PC2"],
                     c=[colors.get(g, "gray") for g in pca_df_merged["group"]],
                     s=80, edgecolors="black", linewidths=0.5)
r_val, p_val = stats.pearsonr(list(sample_order), pca_df_merged["PC2"])
ax.set_xlabel("Sample Processing Order")
ax.set_ylabel("PC2")
ax.set_title(f"PC2 vs Processing Order\n(r={r_val:.3f}, p={p_val:.3f})")
ax.grid(True, alpha=0.3)

# Inter-group vs intra-group correlation comparison (batch indicator)
ax = axes[2]
intra_corrs = []
inter_corrs = []
for i, s1 in enumerate(expr.columns):
    g1 = meta.loc[meta["sample_id"] == s1, "group"].values[0]
    for j, s2 in enumerate(expr.columns):
        if i >= j:
            continue
        g2 = meta.loc[meta["sample_id"] == s2, "group"].values[0]
        r = corr.loc[s1, s2]
        if g1 == g2:
            intra_corrs.append(r)
        else:
            inter_corrs.append(r)

ax.boxplot([intra_corrs, inter_corrs], labels=["Within-Group", "Between-Group"],
           patch_artist=True,
           boxprops=dict(facecolor="#E3F2FD", alpha=0.7))
t_stat, t_pval = stats.ttest_ind(intra_corrs, inter_corrs)
ax.set_ylabel("Pearson Correlation")
ax.set_title(f"Intra vs Inter-Group Correlation\n(t-test p={t_pval:.2e})")
ax.grid(True, alpha=0.3)

plt.suptitle("Batch Effect Assessment", fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "qc_05_batch_effect_assessment.png"))
plt.close()
print("    ✓ Saved: figures/qc_05_batch_effect_assessment.png")

# Statistical tests: does condition/timepoint associate with top PCs?
print("\n  ── Statistical association of variables with PCs ──")
print(f"  {'Variable':<15s} {'PC':<5s} {'Test':<12s} {'Statistic':>10s} {'P-value':>12s} {'Significant?':>14s}")
print("  " + "-" * 70)

for var in ["condition", "timepoint", "group"]:
    for pc_name in ["PC1", "PC2", "PC3"]:
        groups_data = []
        group_labels = pca_df[var].unique()
        for gl in group_labels:
            groups_data.append(pca_df.loc[pca_df[var] == gl, pc_name].values)
        
        if len(groups_data) == 2:
            stat, pval = stats.ttest_ind(groups_data[0], groups_data[1])
            test = "t-test"
        else:
            stat, pval = stats.f_oneway(*groups_data)
            test = "ANOVA"
        
        sig = "YES ***" if pval < 0.001 else "YES **" if pval < 0.01 else "YES *" if pval < 0.05 else "No"
        print(f"  {var:<15s} {pc_name:<5s} {test:<12s} {stat:>10.3f} {pval:>12.2e} {sig:>14s}")

# ══════════════════════════════════════════════════════════════════════════
# PART B: PREPROCESSING
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART B: PREPROCESSING")
print("=" * 70)

# ── B1. Filter low-expression probes using detection p-values ────────────
print("\n  [B1] Filtering low-expression probes...")

# A probe is considered "detected" in a sample if detection p-value < 0.05
# Keep probes detected in at least 50% of samples in any ONE group
det_mask = det_pvals < 0.05  # True = detected

probe_keep = pd.Series(False, index=expr.index)
for grp in meta["group"].unique():
    grp_samples = meta.loc[meta["group"] == grp, "sample_id"].values
    grp_det = det_mask[grp_samples]
    # Probe detected in >= 50% of samples in this group
    detected_in_grp = grp_det.mean(axis=1) >= 0.5
    probe_keep = probe_keep | detected_in_grp

n_before = len(expr)
expr_filt = expr.loc[probe_keep]
n_after = len(expr_filt)
print(f"    Probes before filtering : {n_before:,}")
print(f"    Probes after filtering  : {n_after:,}")
print(f"    Probes removed          : {n_before - n_after:,} ({(n_before - n_after)/n_before*100:.1f}%)")

# ── B2. Quantile normalization ───────────────────────────────────────────
print("\n  [B2] Quantile normalization...")

def quantile_normalize(df):
    """Quantile normalize a DataFrame (probes × samples)."""
    rank_mean = df.stack().groupby(
        df.rank(method="first").stack().astype(int)
    ).mean()
    normalized = df.rank(method="min").stack().astype(int).map(rank_mean).unstack()
    return normalized

expr_norm = quantile_normalize(expr_filt)
print(f"    ✓ Quantile normalization complete.")
print(f"    Shape: {expr_norm.shape}")

# ── B3. Remove zero-variance probes ──────────────────────────────────────
print("\n  [B3] Removing zero/near-zero variance probes...")

probe_var = expr_norm.var(axis=1)
var_threshold = probe_var.quantile(0.01)  # bottom 1%
expr_norm = expr_norm.loc[probe_var > var_threshold]
print(f"    Probes after variance filter: {expr_norm.shape[0]:,}")

# ── B4. Probe-to-gene mapping ───────────────────────────────────────────
print("\n  [B4] Probe-to-gene mapping...")

# Get platform annotation
gpl = list(gse.gpls.values())[0]
annot = gpl.table[["ID", "Symbol", "Entrez_Gene_ID", "ILMN_Gene", "Definition"]].copy()
annot = annot.rename(columns={"ID": "Probe_ID"})
annot = annot.dropna(subset=["Symbol"])
annot = annot[annot["Symbol"].str.strip() != ""]
annot = annot.set_index("Probe_ID")

print(f"    Annotated probes in platform: {len(annot):,}")

# Map probes to genes
common_probes = expr_norm.index.intersection(annot.index)
expr_annot = expr_norm.loc[common_probes].copy()
expr_annot["Symbol"] = annot.loc[common_probes, "Symbol"]
print(f"    Probes with gene annotation : {len(expr_annot):,}")

# Collapse probes to genes: keep probe with highest mean expression
expr_annot["mean_expr"] = expr_annot[samples].mean(axis=1)
expr_annot = expr_annot.sort_values("mean_expr", ascending=False)
expr_gene = expr_annot.drop_duplicates(subset="Symbol", keep="first")
gene_symbols = expr_gene["Symbol"].values
expr_gene = expr_gene[samples].copy()
expr_gene.index = gene_symbols
expr_gene.index.name = "Gene"

print(f"    Unique genes after collapse  : {expr_gene.shape[0]:,}")

# ══════════════════════════════════════════════════════════════════════════
# PART C: POST-NORMALIZATION QC
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART C: POST-NORMALIZATION QC")
print("=" * 70)

# ── C1. Box plot — normalized expression ─────────────────────────────────
print("\n  [C1] Box plot — normalized expression...")

fig, axes = plt.subplots(2, 1, figsize=(16, 12))

# Before (raw, filtered)
ax = axes[0]
bp = ax.boxplot([expr_filt[s].dropna().values for s in expr_filt.columns],
                patch_artist=True, showfliers=False, widths=0.6)
for patch, color in zip(bp["boxes"], sample_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_xticklabels(meta["short_label"].values, rotation=90, fontsize=7)
ax.set_ylabel("Expression (log2)")
ax.set_title("Before Normalization (After Probe Filtering)")

# After normalization
ax = axes[1]
bp = ax.boxplot([expr_norm[s].dropna().values for s in expr_norm.columns],
                patch_artist=True, showfliers=False, widths=0.6)
for patch, color in zip(bp["boxes"], sample_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_xticklabels(meta["short_label"].values, rotation=90, fontsize=7)
ax.set_ylabel("Expression (log2)")
ax.set_title("After Quantile Normalization")

legend_items = [Patch(facecolor=c, label=g, alpha=0.7) for g, c in colors.items()]
axes[0].legend(handles=legend_items, loc="upper right", fontsize=7)
axes[1].legend(handles=legend_items, loc="upper right", fontsize=7)

plt.suptitle("Expression Distribution — Before vs After Normalization", fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "qc_06_normalization_comparison.png"))
plt.close()
print("    ✓ Saved: figures/qc_06_normalization_comparison.png")

# ── C2. PCA — post-normalization (gene-level) ───────────────────────────
print("  [C2] PCA on normalized gene-level data...")

X_gene = expr_gene.T.values
X_gene_scaled = StandardScaler().fit_transform(X_gene)

pca_norm = PCA(n_components=min(10, len(samples)))
pcs_norm = pca_norm.fit_transform(X_gene_scaled)
var_norm = pca_norm.explained_variance_ratio_ * 100

pca_norm_df = pd.DataFrame({
    "PC1": pcs_norm[:, 0], "PC2": pcs_norm[:, 1], "PC3": pcs_norm[:, 2],
    "sample_id": expr_gene.columns
})
pca_norm_df = pca_norm_df.merge(meta[["sample_id", "condition", "timepoint", "group", "short_label"]], on="sample_id")

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# PC1 vs PC2 — by Group
ax = axes[0, 0]
for grp, color in colors.items():
    mask = pca_norm_df["group"] == grp
    ax.scatter(pca_norm_df.loc[mask, "PC1"], pca_norm_df.loc[mask, "PC2"],
               c=color, label=grp, s=100, edgecolors="black", linewidths=0.5, alpha=0.85)
    # Add sample labels
    for _, row in pca_norm_df.loc[mask].iterrows():
        ax.annotate(row["short_label"], (row["PC1"], row["PC2"]),
                    fontsize=5, alpha=0.7, xytext=(4, 4), textcoords="offset points")
ax.set_xlabel(f"PC1 ({var_norm[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_norm[1]:.1f}%)")
ax.set_title("Normalized PCA: by Group")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# PC1 vs PC2 — by Condition
ax = axes[0, 1]
for cond, color in cond_colors.items():
    mask = pca_norm_df["condition"] == cond
    ax.scatter(pca_norm_df.loc[mask, "PC1"], pca_norm_df.loc[mask, "PC2"],
               c=color, label=cond, s=100, edgecolors="black", linewidths=0.5, alpha=0.85)
ax.set_xlabel(f"PC1 ({var_norm[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_norm[1]:.1f}%)")
ax.set_title("Normalized PCA: by Condition (NGT vs PD)")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# PC1 vs PC2 — by Timepoint
ax = axes[1, 0]
for tp, color in tp_colors.items():
    mask = pca_norm_df["timepoint"] == tp
    ax.scatter(pca_norm_df.loc[mask, "PC1"], pca_norm_df.loc[mask, "PC2"],
               c=color, label=tp, s=100, edgecolors="black", linewidths=0.5, alpha=0.85)
ax.set_xlabel(f"PC1 ({var_norm[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_norm[1]:.1f}%)")
ax.set_title("Normalized PCA: by Timepoint (Fasting vs OGTT)")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Scree plot — normalized
ax = axes[1, 1]
ax.bar(range(1, len(var_norm) + 1), var_norm, color="#607D8B", alpha=0.7)
ax.plot(range(1, len(var_norm) + 1), np.cumsum(var_norm), "ro-", markersize=6)
ax.set_xlabel("Principal Component")
ax.set_ylabel("Variance Explained (%)")
ax.set_title("Scree Plot — After Normalization")
ax.set_xticks(range(1, len(var_norm) + 1))
ax.grid(True, alpha=0.3)
for i, cum in enumerate(np.cumsum(var_norm)):
    ax.annotate(f"{cum:.0f}%", (i + 1, cum), textcoords="offset points",
                xytext=(0, 8), ha="center", fontsize=7, color="red")

plt.suptitle("PCA Analysis — After Normalization (Gene-Level)", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "qc_07_pca_normalized.png"))
plt.close()
print("    ✓ Saved: figures/qc_07_pca_normalized.png")

# ── C3. PC association tests — post-normalization ────────────────────────
print("\n  ── Post-normalization: Statistical association of variables with PCs ──")
print(f"  {'Variable':<15s} {'PC':<5s} {'Test':<12s} {'Statistic':>10s} {'P-value':>12s} {'Significant?':>14s}")
print("  " + "-" * 70)

for var in ["condition", "timepoint", "group"]:
    for pc_name in ["PC1", "PC2", "PC3"]:
        groups_data = []
        group_labels = pca_norm_df[var].unique()
        for gl in group_labels:
            groups_data.append(pca_norm_df.loc[pca_norm_df[var] == gl, pc_name].values)
        
        if len(groups_data) == 2:
            stat, pval = stats.ttest_ind(groups_data[0], groups_data[1])
            test = "t-test"
        else:
            stat, pval = stats.f_oneway(*groups_data)
            test = "ANOVA"
        
        sig = "YES ***" if pval < 0.001 else "YES **" if pval < 0.01 else "YES *" if pval < 0.05 else "No"
        print(f"  {var:<15s} {pc_name:<5s} {test:<12s} {stat:>10.3f} {pval:>12.2e} {sig:>14s}")

# ── C4. Sample correlation heatmap — post normalization ─────────────────
print("\n  [C4] Post-normalization correlation heatmap...")

corr_norm = expr_gene.corr(method="pearson")
row_colors_norm = meta.set_index("sample_id")["group"].map(group_lut)

g = sns.clustermap(corr_norm, cmap="RdYlBu_r", vmin=corr_norm.min().min(), vmax=1,
                   figsize=(14, 12),
                   row_colors=row_colors_norm, col_colors=row_colors_norm,
                   xticklabels=meta["short_label"].values,
                   yticklabels=meta["short_label"].values,
                   linewidths=0.5)
g.fig.suptitle("Sample-Sample Correlation (After Normalization, Gene-Level)", y=1.02, fontsize=14)
g.savefig(os.path.join(FIG_DIR, "qc_08_correlation_normalized.png"))
plt.close("all")
print("    ✓ Saved: figures/qc_08_correlation_normalized.png")

# ── C5. Check for outlier samples ────────────────────────────────────────
print("\n  [C5] Outlier detection...")

# Method: samples with mean correlation < mean - 2*SD are potential outliers
mean_corr_per_sample = corr_norm.mean()
overall_mean = mean_corr_per_sample.mean()
overall_std = mean_corr_per_sample.std()
outlier_threshold = overall_mean - 2 * overall_std

print(f"    Mean inter-sample correlation : {overall_mean:.4f}")
print(f"    Std                            : {overall_std:.4f}")
print(f"    Outlier threshold (mean-2SD)   : {outlier_threshold:.4f}")

outliers = mean_corr_per_sample[mean_corr_per_sample < outlier_threshold]
if len(outliers) > 0:
    print(f"\n    ⚠ Potential outlier samples:")
    for s, r in outliers.items():
        grp = meta.loc[meta["sample_id"] == s, "group"].values[0]
        print(f"      {s} ({grp}): mean corr = {r:.4f}")
else:
    print(f"    ✓ No outlier samples detected.")

# ══════════════════════════════════════════════════════════════════════════
# SAVE OUTPUTS
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SAVING PROCESSED DATA")
print("=" * 70)

# Save normalized probe-level expression
norm_path = os.path.join(DATA_DIR, "normalized_expression.csv")
expr_norm.to_csv(norm_path)
print(f"  ✓ Normalized expression (probe-level): {norm_path}")

# Save gene-level expression
gene_path = os.path.join(DATA_DIR, "gene_expression.csv")
expr_gene.to_csv(gene_path)
print(f"  ✓ Gene-level expression: {gene_path}")

# ══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SUMMARY — QC & Preprocessing Results")
print("=" * 70)
print(f"""
  ┌──────────────────────────────────────────────────────────┐
  │ PREPROCESSING PIPELINE SUMMARY                           │
  ├──────────────────────────────────────────────────────────┤
  │ Raw probes                : {expr.shape[0]:>8,}                   │
  │ After detection filtering : {expr_filt.shape[0]:>8,}                   │
  │ After variance filtering  : {expr_norm.shape[0]:>8,}                   │
  │ Annotated probes          : {len(common_probes):>8,}                   │
  │ Final unique genes        : {expr_gene.shape[0]:>8,}                   │
  │ Samples                   : {expr_gene.shape[1]:>8}                   │
  ├──────────────────────────────────────────────────────────┤
  │ Data scale     : log2 (already transformed)              │
  │ Normalization  : Quantile normalization                  │
  │ Gene collapse  : Max mean expression probe per gene      │
  │ Outliers       : {'None detected' if len(outliers)==0 else f'{len(outliers)} potential outlier(s)'}                            │
  └──────────────────────────────────────────────────────────┘

  Figures generated:
    - qc_01_raw_boxplot.png
    - qc_02_raw_density.png
    - qc_03_correlation_heatmap.png
    - qc_04_pca_raw_multivariable.png
    - qc_05_batch_effect_assessment.png
    - qc_06_normalization_comparison.png
    - qc_07_pca_normalized.png
    - qc_08_correlation_normalized.png

  Next Step: Run 03_differential_expression.py
""")

print("=" * 70)
print("Script 02 completed successfully!")
print("=" * 70)
