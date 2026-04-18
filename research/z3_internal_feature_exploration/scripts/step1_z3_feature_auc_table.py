#!/usr/bin/env python3
"""Step 1: Z3 Internal Feature AUC Scan.

Goal: test whether any candidate feature separates TP from FP within Z3
      (TO mode, LOH + extreme AF + HPFineNGroups <= 1).

Output: per-sample x feature AUC table + heatmap.
"""
import os
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore')

# ── Config ─────────────────────────────────────────────────────────
MASTER = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
          "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz")
OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/z3_internal_feature_exploration"
FIG_DIR = os.path.join(OUT_DIR, "figures")
DATA_DIR = os.path.join(OUT_DIR, "data")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

SAMPLE_ORDER = ["H2009", "H1437", "HCC1395", "HCC1395_DORADO",
                "HCC1937", "HCC1954", "COLO829"]

# Candidate features (abs-value AUC: higher = more TP-distinguishing)
FEATURES = [
    ("caller_af",             "caller"),
    ("AlleleDelta",           "methylation-HP"),
    ("HPFineF",               "methylation-HP"),
    ("PairwiseMedianDist",    "methylation-HP-free"),
    ("HeuristicScore",        "methylation-HP-free"),
    ("HPMergedDelta",         "methylation-HP-free"),
    ("NumCpGs",               "context"),
    ("NumReads",              "coverage"),
    ("Coverage_Multiple",     "CN"),
    ("LabelAllelePermanovaF", "methylation-HP"),
    ("HPFineNGroups",         "control"),
    ("nhp_ratio",             "self-phasing"),
]

COLS_NEEDED = [
    "sample", "mode", "truth_label",
    "core_loh_like", "to_loh_bed_hit",
    "caller_af", "HPFineNGroups", "NumReads", "NumCpGs",
    "Coverage_Multiple", "AlleleDelta", "HPFineF",
    "PairwiseMedianDist", "HeuristicScore", "HPMergedDelta",
    "LabelAllelePermanovaF", "NHP3", "NHP0",
]

# ── Load Data ──────────────────────────────────────────────────────
print("Loading master dataset...")
df = pd.read_csv(MASTER, sep='\t', usecols=COLS_NEEDED, low_memory=False)
print(f"Loaded {len(df):,} rows")

# Numeric conversions
num_cols = [c for c in COLS_NEEDED if c not in ("sample", "mode", "truth_label",
                                                "core_loh_like", "to_loh_bed_hit")]
for c in num_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

for col in ["core_loh_like", "to_loh_bed_hit"]:
    df[col] = df[col].map({"True": True, "False": False, True: True, False: False}).fillna(False)

df["is_tp"] = df["truth_label"] == "TP"
df["loh_flag"] = np.where(df["mode"] == "paired", df["core_loh_like"], df["to_loh_bed_hit"])


# Composite feature: self-phasing artifact proxy — (NHP3 + NHP0) / NumReads
df["nhp_ratio"] = (df["NHP3"].fillna(0) + df["NHP0"].fillna(0)) / df["NumReads"].replace(0, np.nan)


# ── Zone Z3 Assignment (TO mode only) ──────────────────────────────
def is_extreme_af(af):
    if pd.isna(af):
        return False
    return af < 0.1 or af > 0.9


df_to = df[df["mode"] == "to"].copy()
df_to["af_extreme"] = df_to["caller_af"].apply(is_extreme_af)
df_to["in_z3"] = df_to["loh_flag"] & df_to["af_extreme"] & (df_to["HPFineNGroups"] <= 1)

print(f"\nTO Z3 total regions: {df_to['in_z3'].sum():,}")

# ── Per-sample AUC Scan ────────────────────────────────────────────
print("\n=== Per-Sample Z3 Feature AUC ===")

results = []
for sample in SAMPLE_ORDER:
    ds = df_to[(df_to["sample"] == sample) & df_to["in_z3"]]
    n = len(ds)
    n_tp = ds["is_tp"].sum()
    n_fp = n - n_tp

    if n < 20 or n_tp < 5 or n_fp < 5:
        print(f"  {sample:16s}: insufficient N=({n}, TP={n_tp}, FP={n_fp}) — skipping")
        for feat, family in FEATURES:
            results.append({
                "sample": sample, "feature": feat, "family": family,
                "n": n, "n_tp": n_tp, "n_fp": n_fp,
                "auc": np.nan, "auc_abs": np.nan, "note": "insufficient_n"
            })
        continue

    for feat, family in FEATURES:
        values = ds[feat].values
        labels = ds["is_tp"].values
        mask = ~np.isnan(values)

        if mask.sum() < 10 or labels[mask].sum() < 3 or (~labels[mask]).sum() < 3:
            auc = np.nan
        else:
            try:
                auc = roc_auc_score(labels[mask], values[mask])
            except ValueError:
                auc = np.nan

        auc_abs = max(auc, 1 - auc) if not np.isnan(auc) else np.nan
        results.append({
            "sample": sample, "feature": feat, "family": family,
            "n": n, "n_tp": int(n_tp), "n_fp": int(n_fp),
            "auc": auc, "auc_abs": auc_abs,
            "note": "ok" if not np.isnan(auc) else "nan_feature"
        })

    # Print top-3 features for this sample
    sdf = pd.DataFrame([r for r in results if r["sample"] == sample and not np.isnan(r.get("auc_abs", np.nan))])
    if len(sdf):
        top3 = sdf.nlargest(3, "auc_abs")
        print(f"  {sample:16s} N={n:5d} (TP={n_tp}, FP={n_fp}): "
              + ", ".join([f"{r['feature']}={r['auc_abs']:.3f}" for _, r in top3.iterrows()]))

# ── Save ───────────────────────────────────────────────────────────
res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(DATA_DIR, "z3_feature_auc.tsv"), sep='\t', index=False)

# ── Heatmap ────────────────────────────────────────────────────────
pivot = res_df.pivot_table(index="feature", columns="sample", values="auc_abs")
pivot = pivot.reindex([f[0] for f in FEATURES])
pivot = pivot[SAMPLE_ORDER]

fig, ax = plt.subplots(figsize=(12, 8))
data = pivot.values
im = ax.imshow(data, cmap='RdYlGn', aspect='auto', vmin=0.5, vmax=0.75)

ax.set_xticks(range(len(SAMPLE_ORDER)))
ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha='right', fontsize=9)
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=9)

for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        val = data[i, j]
        if not np.isnan(val):
            color = 'white' if val < 0.55 or val > 0.68 else 'black'
            weight = 'bold' if val >= 0.60 else 'normal'
            ax.text(j, i, f'{val:.3f}', ha='center', va='center',
                    fontsize=8, color=color, fontweight=weight)

ax.set_title("Z3 Internal Feature AUC (|AUC| for TP vs FP in TO Z3 per sample)\n"
             "Green ≥0.60 = potentially useful; grey <0.55 = no signal",
             fontweight='bold')
plt.colorbar(im, ax=ax, label='|AUC|', shrink=0.8)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "z3_feature_auc_heatmap.png"), dpi=150, bbox_inches='tight')
plt.close()

# ── Summary ────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY — Step 1 Decision Gate")
print("=" * 70)

valid = res_df[res_df["auc_abs"].notna()]

# Count features where AUC >= 0.60 across >=3 samples
print("\nFeatures with |AUC| ≥ 0.60 in ≥3 samples:")
gate_hits = []
for feat, _ in FEATURES:
    fdf = valid[valid["feature"] == feat]
    n_hits = (fdf["auc_abs"] >= 0.60).sum()
    n_samples = len(fdf)
    mean_auc = fdf["auc_abs"].mean()
    if n_hits >= 3:
        samples_hit = fdf[fdf["auc_abs"] >= 0.60]["sample"].tolist()
        print(f"  {feat:26s} — {n_hits}/{n_samples} samples, mean={mean_auc:.3f}, "
              f"hits: {samples_hit}")
        gate_hits.append(feat)

if not gate_hits:
    print("  (none — Step 1 NEGATIVE gate)")

# Global top-5 features by max across-samples AUC
print("\nTop-5 features by max |AUC| across samples:")
top5 = valid.nlargest(5, "auc_abs")[["sample", "feature", "auc_abs", "n_tp", "n_fp"]]
print(top5.to_string(index=False))

print(f"\nOutputs: {DATA_DIR}/z3_feature_auc.tsv")
print(f"         {FIG_DIR}/z3_feature_auc_heatmap.png")
print("\nDecision:")
if gate_hits:
    print(f"  → POSITIVE gate: proceed to Step 2 on features {gate_hits}")
else:
    print("  → NEGATIVE gate: Step 2 skipped; run Step 2.5 independently (H-Z3d)")
