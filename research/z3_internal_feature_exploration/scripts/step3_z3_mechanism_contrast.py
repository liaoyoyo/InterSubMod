#!/usr/bin/env python3
"""Step 3: HCC1954 vs HCC1395 Z3 Mechanism Contrast.

Context: Both are breast cancer cell lines but Z3 TP rate differs 11x
    (HCC1954 ~0.050, HCC1395 ~0.551). Step 1 showed HCC1954 Z3
    NumReads/Coverage_Multiple AUC=0.673 (only single-sample signal).
    This step documents the HCC1954-specific mechanism.

Analyses:
    1. chr-level Z3 FP/TP distribution (germline cluster hypothesis)
    2. PairwiseMedianDist distribution (germline-like reads proxy)
    3. NumCpGs distribution (CpG island bias)
    4. NumReads distribution (NR confound / coverage effect)
    5. Coverage_Multiple distribution (CN artifact)
"""
import os
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings('ignore')

MASTER = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
          "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz")
OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/z3_internal_feature_exploration"
FIG_DIR = os.path.join(OUT_DIR, "figures")
DATA_DIR = os.path.join(OUT_DIR, "data")

COLS = [
    "sample", "mode", "truth_label", "Chr",
    "core_loh_like", "to_loh_bed_hit",
    "caller_af", "HPFineNGroups", "NumReads", "NumCpGs",
    "Coverage_Multiple", "AlleleDelta", "PairwiseMedianDist",
    "HeuristicScore",
]

print("Loading master dataset...")
df = pd.read_csv(MASTER, sep='\t', usecols=COLS, low_memory=False)

num_cols = ["caller_af", "HPFineNGroups", "NumReads", "NumCpGs",
            "Coverage_Multiple", "AlleleDelta", "PairwiseMedianDist",
            "HeuristicScore"]
for c in num_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

for col in ["core_loh_like", "to_loh_bed_hit"]:
    df[col] = df[col].map({"True": True, "False": False,
                           True: True, False: False}).fillna(False)

df["is_tp"] = df["truth_label"] == "TP"
df["loh_flag"] = np.where(df["mode"] == "paired", df["core_loh_like"],
                          df["to_loh_bed_hit"])

df_to = df[df["mode"] == "to"].copy()
df_to["af_extreme"] = df_to["caller_af"].apply(
    lambda x: (not pd.isna(x)) and (x < 0.1 or x > 0.9))
df_to["in_z3"] = (df_to["loh_flag"] & df_to["af_extreme"]
                  & (df_to["HPFineNGroups"] <= 1))

z3_1954 = df_to[(df_to["sample"] == "HCC1954") & df_to["in_z3"]].copy()
z3_1395 = df_to[(df_to["sample"] == "HCC1395") & df_to["in_z3"]].copy()

print(f"\nHCC1954 Z3: n={len(z3_1954)}, TP={z3_1954['is_tp'].sum()}, "
      f"FP={(~z3_1954['is_tp']).sum()}, TP rate={z3_1954['is_tp'].mean():.3f}")
print(f"HCC1395 Z3: n={len(z3_1395)}, TP={z3_1395['is_tp'].sum()}, "
      f"FP={(~z3_1395['is_tp']).sum()}, TP rate={z3_1395['is_tp'].mean():.3f}")

# ── 1. chr-level Z3 FP/TP distribution ─────────────────────────────
print("\n=== 1. chr-level Z3 distribution ===")

chr_rows = []
for sample, ds in [("HCC1954", z3_1954), ("HCC1395", z3_1395)]:
    tot = len(ds)
    for chrom, g in ds.groupby("Chr"):
        n = len(g)
        n_tp = g["is_tp"].sum()
        n_fp = n - n_tp
        chr_rows.append({
            "sample": sample, "chr": chrom, "n": n,
            "n_tp": int(n_tp), "n_fp": int(n_fp),
            "pct_of_sample_z3": n / tot if tot else 0,
            "tp_rate": n_tp / n if n else np.nan,
        })

chr_df = pd.DataFrame(chr_rows)
chr_df = chr_df.sort_values(["sample", "n"], ascending=[True, False])
chr_df.to_csv(os.path.join(DATA_DIR, "z3_mechanism_chr_dist.tsv"),
              sep='\t', index=False)

# Top-5 chromosomes by volume per sample
for sample in ["HCC1954", "HCC1395"]:
    sdf = chr_df[chr_df["sample"] == sample].head(5)
    print(f"\n  {sample} top-5 chr (by Z3 count):")
    for _, r in sdf.iterrows():
        print(f"    {r['chr']:8s} n={r['n']:5d} "
              f"({r['pct_of_sample_z3']*100:5.1f}%) "
              f"TP={r['n_tp']:5d} FP={r['n_fp']:5d} "
              f"TP_rate={r['tp_rate']:.3f}")

# ── 2-5. Feature distribution contrast (TP vs FP within each sample) ───
print("\n=== 2-5. Feature distribution contrast (TP vs FP) ===")

features = ["PairwiseMedianDist", "NumCpGs", "NumReads", "Coverage_Multiple"]
dist_rows = []
for sample, ds in [("HCC1954", z3_1954), ("HCC1395", z3_1395)]:
    tp = ds[ds["is_tp"]]
    fp = ds[~ds["is_tp"]]
    for feat in features:
        tp_v = tp[feat].dropna()
        fp_v = fp[feat].dropna()
        if len(tp_v) >= 5 and len(fp_v) >= 5:
            mw_stat, mw_p = stats.mannwhitneyu(tp_v, fp_v, alternative='two-sided')
        else:
            mw_stat, mw_p = np.nan, np.nan
        dist_rows.append({
            "sample": sample, "feature": feat,
            "tp_n": len(tp_v), "tp_median": tp_v.median(),
            "tp_q25": tp_v.quantile(0.25), "tp_q75": tp_v.quantile(0.75),
            "fp_n": len(fp_v), "fp_median": fp_v.median(),
            "fp_q25": fp_v.quantile(0.25), "fp_q75": fp_v.quantile(0.75),
            "mw_p": mw_p,
        })

dist_df = pd.DataFrame(dist_rows)
dist_df.to_csv(os.path.join(DATA_DIR, "z3_mechanism_feature_dist.tsv"),
               sep='\t', index=False)

print(dist_df.to_string(index=False))

# ── Figure: grid of violin plots HCC1954 vs HCC1395, 4 features × TP/FP ──
fig, axes = plt.subplots(2, 4, figsize=(16, 8))

for i, (sample, ds) in enumerate([("HCC1954", z3_1954), ("HCC1395", z3_1395)]):
    tp = ds[ds["is_tp"]]
    fp = ds[~ds["is_tp"]]
    for j, feat in enumerate(features):
        ax = axes[i, j]
        tp_v = tp[feat].dropna().values
        fp_v = fp[feat].dropna().values

        if len(tp_v) and len(fp_v):
            parts = ax.violinplot([tp_v, fp_v], positions=[1, 2],
                                  showmedians=True, widths=0.7)
            for pc, color in zip(parts['bodies'], ['#2ecc71', '#e74c3c']):
                pc.set_facecolor(color)
                pc.set_alpha(0.6)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"TP\n(n={len(tp_v)})",
                            f"FP\n(n={len(fp_v)})"], fontsize=8)
        ax.set_title(f"{sample} — {feat}", fontsize=9)
        if feat == "PairwiseMedianDist":
            ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)

fig.suptitle("Z3 Mechanism Contrast — HCC1954 (top, TP rate=0.050) vs "
             "HCC1395 (bottom, TP rate=0.551)",
             fontsize=11, fontweight='bold', y=1.00)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "z3_hcc1954_vs_hcc1395_mechanism.png"),
            dpi=150, bbox_inches='tight')
plt.close()

# ── Figure: chromosomal Z3 density (HCC1954 vs HCC1395) ────────────
# Order chromosomes: 1-22, X, Y
chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

fig, axes = plt.subplots(2, 1, figsize=(14, 6), sharex=True)
for ax, sample in zip(axes, ["HCC1954", "HCC1395"]):
    sdf = chr_df[chr_df["sample"] == sample].set_index("chr")
    sdf = sdf.reindex([c for c in chr_order if c in sdf.index])
    x = np.arange(len(sdf))
    ax.bar(x, sdf["n_tp"], color="#2ecc71", label=f"TP (n={sdf['n_tp'].sum()})")
    ax.bar(x, sdf["n_fp"], bottom=sdf["n_tp"],
           color="#e74c3c", label=f"FP (n={sdf['n_fp'].sum()})")
    ax.set_xticks(x)
    ax.set_xticklabels(sdf.index, rotation=60, fontsize=7, ha='right')
    ax.set_ylabel("Z3 region count")
    ax.set_title(f"{sample} — Z3 regions per chromosome "
                 f"(overall TP rate={sdf['n_tp'].sum()/(sdf['n_tp'].sum()+sdf['n_fp'].sum()):.3f})",
                 fontsize=10)
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "z3_chr_distribution_contrast.png"),
            dpi=150, bbox_inches='tight')
plt.close()

# ── Summary ────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY — HCC1954 vs HCC1395 Z3 Mechanism")
print("=" * 70)

for sample, ds in [("HCC1954", z3_1954), ("HCC1395", z3_1395)]:
    chr_top = chr_df[chr_df["sample"] == sample].head(3)
    chr_str = ", ".join([f"{r['chr']}({r['pct_of_sample_z3']*100:.0f}%)"
                         for _, r in chr_top.iterrows()])
    print(f"\n{sample}:")
    print(f"  Z3 n={len(ds)}, TP rate={ds['is_tp'].mean():.3f}")
    print(f"  Top-3 chr: {chr_str}")
    print(f"  NumReads median (TP/FP): "
          f"{ds[ds['is_tp']]['NumReads'].median():.0f} / "
          f"{ds[~ds['is_tp']]['NumReads'].median():.0f}")
    print(f"  Coverage_Multiple median (TP/FP): "
          f"{ds[ds['is_tp']]['Coverage_Multiple'].median():.2f} / "
          f"{ds[~ds['is_tp']]['Coverage_Multiple'].median():.2f}")

print(f"\nOutputs:")
print(f"  {DATA_DIR}/z3_mechanism_chr_dist.tsv")
print(f"  {DATA_DIR}/z3_mechanism_feature_dist.tsv")
print(f"  {FIG_DIR}/z3_hcc1954_vs_hcc1395_mechanism.png")
print(f"  {FIG_DIR}/z3_chr_distribution_contrast.png")
