"""
=============================================================================
Script 01 — Dataset Exploration (Task 1)
=============================================================================
Purpose:
    Download GEO dataset GSE153837 and explore its structure, experimental
    design, sample information, and expression matrix characteristics.

Dataset:
    GSE153837 — Effects of oral glucose load on gene expression of PBMCs
    in Asian-Indian men.
    Platform: Illumina HumanHT-12 v4.0 Expression BeadChip (GPL10558)

Outputs:
    - data/metadata.csv          : Cleaned sample metadata
    - data/raw_expression.csv    : Raw expression matrix (probes × samples)
    - Console: Comprehensive dataset summary
=============================================================================
"""

import os
import sys
import warnings
import pandas as pd
import numpy as np
import GEOparse

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)

# ── 1. Download the dataset ───────────────────────────────────────────────
print("=" * 70)
print("STEP 1: Downloading GSE153837 from GEO")
print("=" * 70)

gse = GEOparse.get_GEO(geo="GSE153837", destdir=DATA_DIR, silent=True)
print(f"✓ Dataset downloaded successfully.\n")

# ── 2. General dataset information ────────────────────────────────────────
print("=" * 70)
print("STEP 2: General Dataset Information")
print("=" * 70)

metadata_keys = ["title", "summary", "overall_design", "type", "platform_id"]
for key in metadata_keys:
    value = gse.metadata.get(key, ["N/A"])
    if isinstance(value, list):
        value = " ".join(value)
    print(f"  {key:>20s}: {value}")

print(f"\n  {'Number of samples':>20s}: {len(gse.gsms)}")
print(f"  {'Number of platforms':>20s}: {len(gse.gpls)}")
for gpl_name, gpl in gse.gpls.items():
    print(f"  {'Platform':>20s}: {gpl_name} — {gpl.metadata.get('title', ['N/A'])[0]}")
    print(f"  {'Technology':>20s}: {gpl.metadata.get('technology', ['N/A'])[0]}")
    print(f"  {'Organism':>20s}: {gpl.metadata.get('organism', ['N/A'])[0]}")

# ── 3. Extract sample phenotype data ─────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 3: Extracting Sample Phenotype Data")
print("=" * 70)

pheno_records = []
for gsm_name, gsm in gse.gsms.items():
    record = {"sample_id": gsm_name, "title": gsm.metadata.get("title", [""])[0]}
    
    # Extract all characteristics
    chars = gsm.metadata.get("characteristics_ch1", [])
    for ch in chars:
        if ":" in ch:
            key, val = ch.split(":", 1)
            record[key.strip().lower().replace(" ", "_")] = val.strip()
        else:
            record[ch.strip()] = ch.strip()
    
    # Also capture source and description
    record["source"] = gsm.metadata.get("source_name_ch1", [""])[0]
    record["description"] = gsm.metadata.get("description", [""])[0] if gsm.metadata.get("description") else ""
    
    pheno_records.append(record)

pheno_df = pd.DataFrame(pheno_records)
print(f"\n  Phenotype table shape: {pheno_df.shape}")
print(f"\n  Columns found: {list(pheno_df.columns)}")
print(f"\n  Full phenotype table:\n")
print(pheno_df.to_string(index=False, max_colwidth=60))

# ── 4. Classify samples into analysis groups ─────────────────────────────
print("\n\n" + "=" * 70)
print("STEP 4: Classifying Samples into Analysis Groups")
print("=" * 70)

# Attempt to identify condition (NGT/PD) and timepoint (Fasting/OGTT)
# from available metadata columns
print("\n  Inspecting sample titles and metadata for group classification...\n")

# Print unique values per column to identify grouping variables
for col in pheno_df.columns:
    unique_vals = pheno_df[col].unique()
    if len(unique_vals) <= 10:
        print(f"  {col}: {list(unique_vals)}")

# Try to parse groups from title or characteristics
# Common patterns: "NGT_Fasting_1", "PD_OGTT_2", etc.
def classify_condition(row):
    """Classify sample condition (NGT or PD) from available metadata."""
    text = " ".join(str(v).lower() for v in row.values)
    if "ngt" in text or "normal glucose tolerance" in text or "normal" in text:
        return "NGT"
    elif "pd" in text or "prediabet" in text or "pre-diabet" in text:
        return "PD"
    return "Unknown"

def classify_timepoint(row):
    """Classify timepoint (Fasting or OGTT) from available metadata."""
    text = " ".join(str(v).lower() for v in row.values)
    if "fasting" in text or "fast" in text or "baseline" in text or "0h" in text or "0 h" in text:
        return "Fasting"
    elif "ogtt" in text or "2h" in text or "2 h" in text or "post" in text or "glucose" in text:
        return "OGTT"
    return "Unknown"

pheno_df["condition"] = pheno_df.apply(classify_condition, axis=1)
pheno_df["timepoint"] = pheno_df.apply(classify_timepoint, axis=1)
pheno_df["group"] = pheno_df["condition"] + "_" + pheno_df["timepoint"]

print(f"\n  Group assignments:")
print(pheno_df[["sample_id", "title", "condition", "timepoint", "group"]].to_string(index=False))

print(f"\n  Group counts:")
print(pheno_df["group"].value_counts().to_string())

# Validate: we expect 4 groups
n_groups = pheno_df["group"].nunique()
print(f"\n  Total unique groups: {n_groups}")
if n_groups == 4:
    print("  ✓ All 4 expected groups identified (NGT_Fasting, NGT_OGTT, PD_Fasting, PD_OGTT)")
else:
    print(f"  ⚠ Expected 4 groups, found {n_groups}. Manual inspection may be needed.")

# ── 5. Extract and summarize expression matrix ──────────────────────────
print("\n" + "=" * 70)
print("STEP 5: Expression Matrix Summary")
print("=" * 70)

# Get the expression data from the first GPL
first_gsm = list(gse.gsms.values())[0]
sample_table = first_gsm.table
print(f"\n  Sample table columns (first sample): {list(sample_table.columns)}")
print(f"  Sample table shape (first sample): {sample_table.shape}")

# Build full expression matrix
expr_dict = {}
for gsm_name, gsm in gse.gsms.items():
    tbl = gsm.table
    if "ID_REF" in tbl.columns and "VALUE" in tbl.columns:
        expr_dict[gsm_name] = tbl.set_index("ID_REF")["VALUE"]

expr_df = pd.DataFrame(expr_dict)
expr_df.index.name = "Probe_ID"

print(f"\n  Expression matrix shape: {expr_df.shape}")
print(f"  Probes: {expr_df.shape[0]:,}")
print(f"  Samples: {expr_df.shape[1]}")
print(f"\n  Data types:\n{expr_df.dtypes.value_counts().to_string()}")

# Convert to numeric if needed
expr_df = expr_df.apply(pd.to_numeric, errors="coerce")

print(f"\n  Missing values: {expr_df.isna().sum().sum()}")
print(f"  Zero values: {(expr_df == 0).sum().sum()}")

print(f"\n  Summary statistics (across all probes and samples):")
print(f"    Min    : {expr_df.min().min():.4f}")
print(f"    Max    : {expr_df.max().max():.4f}")
print(f"    Mean   : {expr_df.mean().mean():.4f}")
print(f"    Median : {expr_df.median().median():.4f}")
print(f"    Std    : {expr_df.std().mean():.4f}")

# Check if data is likely log-transformed
max_val = expr_df.max().max()
if max_val > 100:
    print(f"\n  ⓘ  Data appears to be on RAW scale (max = {max_val:.2f})")
    print("     → Log2 transformation will be needed in QC step")
else:
    print(f"\n  ⓘ  Data appears to be LOG-transformed (max = {max_val:.2f})")
    print("     → No log transformation needed")

# ── 6. Check for detection p-values ──────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 6: Detection P-values")
print("=" * 70)

first_gsm = list(gse.gsms.values())[0]
if "Detection Pval" in first_gsm.table.columns:
    print("  ✓ Detection p-values available — will use for probe filtering in QC")
    pval_col = "Detection Pval"
elif "DETECTION_PVAL" in first_gsm.table.columns:
    print("  ✓ Detection p-values available — will use for probe filtering in QC")
    pval_col = "DETECTION_PVAL"
else:
    extra_cols = [c for c in first_gsm.table.columns if c not in ["ID_REF", "VALUE"]]
    print(f"  Additional columns in sample table: {extra_cols}")
    if extra_cols:
        print("  Will inspect these for detection p-values in QC step")
    else:
        print("  ⚠ No detection p-values found — will use expression-based filtering")

# ── 7. Platform annotation (probe-to-gene mapping) ──────────────────────
print("\n" + "=" * 70)
print("STEP 7: Platform Annotation Preview")
print("=" * 70)

for gpl_name, gpl in gse.gpls.items():
    gpl_table = gpl.table
    print(f"\n  Platform: {gpl_name}")
    print(f"  Annotation table shape: {gpl_table.shape}")
    print(f"  Columns: {list(gpl_table.columns)}")
    
    # Look for gene symbol column
    gene_cols = [c for c in gpl_table.columns if "gene" in c.lower() or "symbol" in c.lower()]
    print(f"  Gene-related columns: {gene_cols}")
    
    if gene_cols:
        print(f"\n  Preview of annotation (first 5 rows):")
        preview_cols = ["ID"] + gene_cols[:3]
        preview_cols = [c for c in preview_cols if c in gpl_table.columns]
        print(gpl_table[preview_cols].head().to_string(index=False))

# ── 8. Save outputs ─────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 8: Saving Outputs")
print("=" * 70)

# Save metadata
meta_path = os.path.join(DATA_DIR, "metadata.csv")
pheno_df.to_csv(meta_path, index=False)
print(f"  ✓ Metadata saved to: {meta_path}")

# Save raw expression
expr_path = os.path.join(DATA_DIR, "raw_expression.csv")
expr_df.to_csv(expr_path)
print(f"  ✓ Raw expression matrix saved to: {expr_path}")

# ── 9. Summary of observations ──────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY — Key Observations for Downstream Analysis")
print("=" * 70)

print(f"""
  1. Dataset: GSE153837 (Illumina HumanHT-12 v4.0)
  2. Samples: {len(gse.gsms)} total
  3. Probes: {expr_df.shape[0]:,}
  4. Groups: {pheno_df['group'].value_counts().to_dict()}
  5. Expression range: [{expr_df.min().min():.2f}, {expr_df.max().max():.2f}]
  6. Missing values: {expr_df.isna().sum().sum()}
  
  Next Step: Run 02_quality_control.py for data quality assessment
""")

print("=" * 70)
print("Script 01 completed successfully!")
print("=" * 70)
