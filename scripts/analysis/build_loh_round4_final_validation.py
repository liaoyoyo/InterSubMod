#!/usr/bin/env python3
"""
LOH Round 4 Final Validation:
  Step 1 — Baseline F1 確認（從 benchmark_comparison.tsv 讀取完整 TP/FP/FN，含 FN 的真實 F1）
  Step 2 — LOH 覆蓋率確認（LOH 數據 vs 完整 benchmark 的 TP 覆蓋比例）
  Step 3 — LOH + HPMergedSig 深度驗證
            （80 個 paired FP 的 per-sample 分解、chr 分佈、HPMergedDelta 分佈、TO 端驗證）
  Step 4 — 全量 Benchmark Filter Simulation（含 FN 的正確 F1 計算）
  Step 5 — A≥30 vs 30≤A<50 vs A≥50 雙層 Tier 定義

生成時間: 2026-03-28
"""

import os

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

# ============================================================
# 路徑設定
# ============================================================
ROUND1_WS   = Path(
    os.environ.get(
        "LOH_ROUND1_WS",
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit",
    )
)
CANONICAL   = Path("/big7_disk/liaoyoyo2001/big7_disk_output/canonical")
OUT_DIR     = Path(
    os.environ.get(
        "LOH_ROUND4_OUT_DIR",
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260328_loh_round4_final_validation",
    )
)
FIG_DIR     = OUT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 樣本 benchmark 路徑（以 InterSubMod row 為主，fallback to LongPhase-S）
BENCHMARK_FILES = {
    "HCC1395":         CANONICAL/"HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/benchmark_comparison.tsv",
    "COLO829":         CANONICAL/"COLO829/paired_full/20260315_COLO829_paired_full_full_complete_matrix/benchmark_comparison.tsv",
    "H1437":           CANONICAL/"H1437/paired_full/20260315_H1437_paired_full_full_complete_matrix/benchmark_comparison.tsv",
    "H2009":           CANONICAL/"H2009/paired_full/20260315_H2009_paired_full_full_complete_matrix/benchmark_comparison.tsv",
    "HCC1937":         CANONICAL/"HCC1937/paired_full/20260315_HCC1937_paired_full_full_complete_matrix/benchmark_comparison.tsv",
    "HCC1954":         CANONICAL/"HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix/benchmark_comparison.tsv",
    "HCC1395_DORADO":  CANONICAL/"HCC1395_DORADO/paired_full/20260315_HCC1395_DORADO_paired_full_full_complete_matrix/benchmark_comparison.tsv",
}
PAIRED_SAMPLES = list(BENCHMARK_FILES.keys())
PREFERRED_METHOD = "LongPhase-S"  # 使用 LongPhase-S baseline（未加 InterSubMod 規則）

# Tier 定義
def assign_tier(n):
    if pd.isna(n): return "C0"
    n = int(n)
    if n == 0:    return "C0"
    elif n < 10:  return "C"
    elif n < 30:  return "B"
    else:         return "A"

def assign_tier_dual(n):
    """雙層 A 定義：A+(≥50) vs A(30-49) vs B vs C vs C0"""
    if pd.isna(n): return "C0"
    n = int(n)
    if n == 0:     return "C0"
    elif n < 10:   return "C"
    elif n < 30:   return "B"
    elif n < 50:   return "A"
    else:          return "A+"

def compute_f1_full(tp, fp, fn):
    """F1 with FN: F1 = 2*TP / (2*TP + FP + FN)"""
    denom = 2*tp + fp + fn
    if denom == 0: return np.nan
    return 2*tp / denom

# ============================================================
# 載入 LOH 資料
# ============================================================
print("[1/8] 載入 LOH 資料 ...")
df_raw = pd.read_csv(ROUND1_WS/"all_region_rows.tsv.gz", sep="\t", low_memory=False)
num_cols = ["effective_hp_reads","hp_ratio_core","hp0_ratio","hp3_ratio",
            "NHP0","NHP3","NumReads","UnassignedAffinity",
            "HPMergedDelta","HPMergedP","HPFineF","HPFineP",
            "AlleleDelta","AlleleP","CramersV","PairwiseMedianDist",
            "HeuristicScore","GlobalP","caller_af"]
for c in num_cols:
    if c in df_raw.columns:
        df_raw[c] = pd.to_numeric(df_raw[c], errors="coerce")
bool_cols = ["HPMergedSig","HPFineSig","AlleleSig","core_loh_like","PassedGating","Significant"]
for c in bool_cols:
    if c in df_raw.columns:
        df_raw[c] = df_raw[c].astype(str).str.strip().str.lower().isin(["true","1","yes"])

df_raw["hp_support_tier"]  = df_raw["effective_hp_reads"].apply(assign_tier)
df_raw["tier_dual"]        = df_raw["effective_hp_reads"].apply(assign_tier_dual)
df = df_raw[df_raw["truth_label"].isin(["TP","FP"])].copy()

print(f"  總計 {len(df):,} TP/FP rows (TP={df['truth_label'].eq('TP').sum():,}, FP={df['truth_label'].eq('FP').sum():,})")

# ============================================================
# Step 1: Baseline F1 確認（從 benchmark_comparison.tsv）
# ============================================================
print("\n[2/8] Step 1: Baseline F1 確認 ...")
baseline_rows = []
for sample, fpath in BENCHMARK_FILES.items():
    if not fpath.exists():
        print(f"  WARNING: {sample} benchmark not found at {fpath}")
        continue
    bdf = pd.read_csv(fpath, sep="\t", comment="#")
    # 優先取 InterSubMod 行，沒有就取 LongPhase-S
    row_ism = bdf[bdf["method"]=="InterSubMod"]
    row_lps = bdf[bdf["method"]=="LongPhase-S"]
    row = row_ism.iloc[0] if len(row_ism) > 0 else (row_lps.iloc[0] if len(row_lps) > 0 else None)
    if row is None:
        continue
    truth_total = int(row["truth_total"])
    tp  = int(row["tp"])
    fp  = int(row["fp"])
    fn  = int(row["fn"])
    pre = float(row["precision"])
    rec = float(row["recall"])
    f1  = float(row["f1"])
    method = str(row["method"])
    baseline_rows.append({
        "sample": sample, "method": method,
        "truth_total": truth_total, "tp": tp, "fp": fp, "fn": fn,
        "precision": pre, "recall": rec, "f1_true": f1
    })
    print(f"  {sample} [{method}]: truth={truth_total:,}, TP={tp:,}, FP={fp:,}, FN={fn:,}, F1={f1:.4f}")

baseline_df = pd.DataFrame(baseline_rows)
baseline_df.to_csv(OUT_DIR/"step1_baseline_f1.tsv", sep="\t", index=False)

# ============================================================
# Step 2: LOH 覆蓋率確認
# ============================================================
print("\n[3/8] Step 2: LOH 覆蓋率確認 ...")
coverage_rows = []
for sample, grp in df[df["mode"]=="paired"].groupby("sample"):
    tp_loh = grp[grp["truth_label"]=="TP"].shape[0]
    fp_loh = grp[grp["truth_label"]=="FP"].shape[0]
    brow = baseline_df[baseline_df["sample"]==sample]
    if len(brow) == 0:
        print(f"  {sample}: no benchmark row, skip")
        continue
    brow = brow.iloc[0]
    cov_tp = tp_loh / brow["tp"] if brow["tp"] > 0 else np.nan
    cov_fp = fp_loh / brow["fp"] if brow["fp"] > 0 else np.nan
    print(f"  {sample}: LOH_TP={tp_loh:,} / bench_TP={brow['tp']:,} = {cov_tp:.3f} | LOH_FP={fp_loh:,} / bench_FP={brow['fp']:,} = {cov_fp:.3f}")
    coverage_rows.append({"sample": sample, "loh_tp": tp_loh, "bench_tp": brow["tp"],
                          "loh_fp": fp_loh, "bench_fp": brow["fp"],
                          "tp_coverage": cov_tp, "fp_coverage": cov_fp})

coverage_df = pd.DataFrame(coverage_rows)
coverage_df.to_csv(OUT_DIR/"step2_loh_coverage.tsv", sep="\t", index=False)

# ============================================================
# Step 3: LOH + HPMergedSig 深度驗證
# ============================================================
print("\n[4/8] Step 3: LOH + HPMergedSig 深度驗證 ...")

# Tier A LOH-like region（paired）
df_paired_A = df[(df["mode"]=="paired") &
                 (df["hp_support_tier"]=="A") &
                 df["core_loh_like"]].copy()

print(f"  Paired Tier A LOH-like: {len(df_paired_A):,} rows "
      f"(TP={df_paired_A['truth_label'].eq('TP').sum():,}, FP={df_paired_A['truth_label'].eq('FP').sum():,})")

# 3a: 2×2 LOH×HPSig 分解（全局 + per-sample）
joint_rows = []
for mode in ["paired","to"]:
    df_A = df[(df["mode"]==mode) & (df["hp_support_tier"]=="A")].copy()
    for loh in [True, False]:
        for sig in [True, False]:
            sub = df_A[(df_A["core_loh_like"]==loh) & (df_A["HPMergedSig"]==sig)]
            tp_n = (sub["truth_label"]=="TP").sum()
            fp_n = (sub["truth_label"]=="FP").sum()
            joint_rows.append({
                "mode": mode, "loh_like": loh, "hp_merged_sig": sig,
                "n_tp": tp_n, "n_fp": fp_n, "n_total": len(sub),
                "fp_frac": fp_n/len(sub) if len(sub) > 0 else np.nan
            })

joint_df = pd.DataFrame(joint_rows)
joint_df.to_csv(OUT_DIR/"step3_joint_loh_hpsig_global.tsv", sep="\t", index=False)

print("\n  Joint 2×2 Table (all samples):")
for _, row in joint_df.iterrows():
    print(f"    {row['mode']}/LOH={row['loh_like']}/HPSig={row['hp_merged_sig']}: "
          f"TP={row['n_tp']:,}, FP={row['n_fp']:,}, FP%={100*row['fp_frac']:.2f}%")

# 3b: 80 FPs per-sample 分解
fp_joint_paired = df_paired_A[(df_paired_A["HPMergedSig"]==True) &
                               (df_paired_A["truth_label"]=="FP")]
print(f"\n  Paired Tier A LOH=Y+HPSig=Y FPs: {len(fp_joint_paired):,} total")
per_sample_fp = fp_joint_paired.groupby("sample").size().reset_index(name="fp_count")
per_sample_tp = df_paired_A[(df_paired_A["HPMergedSig"]==True) & (df_paired_A["truth_label"]=="TP")].groupby("sample").size().reset_index(name="tp_count")
sample_decomp = per_sample_tp.merge(per_sample_fp, on="sample", how="outer").fillna(0)
sample_decomp["fp_count"] = sample_decomp["fp_count"].astype(int)
sample_decomp["tp_count"] = sample_decomp["tp_count"].astype(int)
sample_decomp["fp_frac"] = sample_decomp["fp_count"] / (sample_decomp["tp_count"] + sample_decomp["fp_count"])
sample_decomp = sample_decomp.sort_values("fp_count", ascending=False)
print("  Per-sample (LOH=Y+HPSig=Y):")
for _, row in sample_decomp.iterrows():
    print(f"    {row['sample']}: TP={row['tp_count']}, FP={row['fp_count']}, FP%={100*row['fp_frac']:.1f}%")
sample_decomp.to_csv(OUT_DIR/"step3_fp_per_sample_decomp.tsv", sep="\t", index=False)

# 3c: Consistency test（排除 HCC1395 + HCC1954 後的 enrichment）
excl = ["HCC1395", "HCC1954"]
for label, mask_fn in [
    ("All samples", lambda d: d),
    ("Excl HCC1395+HCC1954", lambda d: d[~d["sample"].isin(excl)])
]:
    sub = mask_fn(df[(df["mode"]=="paired") & (df["hp_support_tier"]=="A")])
    tp_total = (sub["truth_label"]=="TP").sum()
    fp_total = (sub["truth_label"]=="FP").sum()
    # LOH=Y+HPSig=Y vs LOH=Y+HPSig=N
    for sig in [True, False]:
        g = sub[sub["core_loh_like"]==True]
        g = g[g["HPMergedSig"]==sig]
        tp_n = (g["truth_label"]=="TP").sum()
        fp_n = (g["truth_label"]=="FP").sum()
        total = tp_n + fp_n
        fp_pct = 100*fp_n/total if total > 0 else 0
        print(f"  [{label}] LOH=Y HPSig={sig}: TP={tp_n:,}, FP={fp_n:,}, FP%={fp_pct:.2f}%")

# 3d: Chr 分佈
print("\n  Chromosomal distribution of paired Tier A LOH=Y+HPSig=Y FPs:")
if "Chr" in fp_joint_paired.columns:
    chr_dist = fp_joint_paired.groupby("Chr").size().reset_index(name="n_fp")
    chr_dist_tp = df_paired_A[df_paired_A["HPMergedSig"]==True].groupby("Chr").size().reset_index(name="n_total")
    chr_dist = chr_dist.merge(chr_dist_tp, on="Chr", how="outer").fillna(0)
    chr_dist["fp_count"] = chr_dist["n_fp"].astype(int)
    chr_dist["fp_pct"] = 100 * chr_dist["n_fp"] / chr_dist["n_fp"].sum()
    chr_dist = chr_dist.sort_values("n_fp", ascending=False)
    for _, row in chr_dist.head(10).iterrows():
        print(f"    {row['Chr']}: {row['fp_count']} FPs ({row['fp_pct']:.1f}%)")
    chr_dist.to_csv(OUT_DIR/"step3_chr_distribution_fp.tsv", sep="\t", index=False)

# 3e: HPMergedDelta 分佈（TP vs FP in Tier A, split by LOH×HPSig）
print("\n  HPMergedDelta abs stats (Tier A paired):")
for loh, sig, label in [
    (True, True, "LOH=Y+HPSig=Y"),
    (True, False, "LOH=Y+HPSig=N"),
    (False, True, "LOH=N+HPSig=Y"),
    (False, False, "LOH=N+HPSig=N"),
]:
    sub = df[(df["mode"]=="paired") & (df["hp_support_tier"]=="A") &
             (df["core_loh_like"]==loh) & (df["HPMergedSig"]==sig)]
    for truth in ["TP","FP"]:
        g = sub[sub["truth_label"]==truth]
        vals = g["HPMergedDelta"].abs().dropna()
        if len(vals) > 0:
            print(f"    {label}/{truth}: n={len(vals):,}, mean={vals.mean():.4f}, "
                  f"med={vals.median():.4f}, p75={vals.quantile(0.75):.4f}")

# 3f: TO validation（同樣分析 TO 端）
print("\n  TO Tier A LOH×HPSig joint:")
for _, row in joint_df[joint_df["mode"]=="to"].iterrows():
    print(f"    LOH={row['loh_like']}/HPSig={row['hp_merged_sig']}: "
          f"TP={row['n_tp']:,}, FP={row['n_fp']:,}, FP%={100*row['fp_frac']:.2f}%")

# ============================================================
# Step 4: 全量 Benchmark Filter Simulation（含 FN 的正確 F1）
# ============================================================
print("\n[5/8] Step 4: 全量 Benchmark Filter Simulation ...")

# Filter 定義
def apply_filter(df_paired_sample, filter_name):
    """回傳 removed_tp, removed_fp"""
    if filter_name == "TierA_LOH":
        mask = (df_paired_sample["hp_support_tier"]=="A") & df_paired_sample["core_loh_like"]
    elif filter_name == "TierA_LOH_HP0low":
        mask = ((df_paired_sample["hp_support_tier"]=="A") &
                df_paired_sample["core_loh_like"] &
                (df_paired_sample["hp0_ratio"] < 0.05))
    elif filter_name == "TierA_LOH_HPSig":
        mask = ((df_paired_sample["hp_support_tier"]=="A") &
                df_paired_sample["core_loh_like"] &
                df_paired_sample["HPMergedSig"])
    elif filter_name == "TierAplus_LOH":
        mask = (df_paired_sample["tier_dual"]=="A+") & df_paired_sample["core_loh_like"]
    elif filter_name == "TierAplus_LOH_HPSig":
        mask = ((df_paired_sample["tier_dual"]=="A+") &
                df_paired_sample["core_loh_like"] &
                df_paired_sample["HPMergedSig"])
    else:
        raise ValueError(f"Unknown filter: {filter_name}")
    removed = df_paired_sample[mask]
    rm_tp = (removed["truth_label"]=="TP").sum()
    rm_fp = (removed["truth_label"]=="FP").sum()
    return rm_tp, rm_fp

FILTERS = ["TierA_LOH", "TierA_LOH_HP0low", "TierA_LOH_HPSig",
           "TierAplus_LOH", "TierAplus_LOH_HPSig"]

sim_rows = []
for sample, grp in df[df["mode"]=="paired"].groupby("sample"):
    brow = baseline_df[baseline_df["sample"]==sample]
    if len(brow) == 0:
        continue
    brow = brow.iloc[0]
    tp_bench   = int(brow["tp"])
    fp_bench   = int(brow["fp"])
    fn_bench   = int(brow["fn"])
    f1_base    = compute_f1_full(tp_bench, fp_bench, fn_bench)

    for filt in FILTERS:
        rm_tp, rm_fp = apply_filter(grp, filt)
        tp_new = tp_bench - rm_tp
        fp_new = fp_bench - rm_fp
        fn_new = fn_bench + rm_tp  # removed TPs become FN
        f1_new = compute_f1_full(tp_new, fp_new, fn_new)
        f1_delta = f1_new - f1_base
        sim_rows.append({
            "sample": sample, "filter": filt,
            "bench_tp": tp_bench, "bench_fp": fp_bench, "bench_fn": fn_bench,
            "f1_baseline": f1_base,
            "rm_tp": rm_tp, "rm_fp": rm_fp,
            "rm_tp_pct": rm_tp/tp_bench if tp_bench > 0 else 0,
            "rm_fp_pct": rm_fp/fp_bench if fp_bench > 0 else 0,
            "tp_new": tp_new, "fp_new": fp_new, "fn_new": fn_new,
            "f1_new": f1_new, "f1_delta": f1_delta
        })
    print(f"  {sample}: baseline F1={f1_base:.4f} (TP={tp_bench}, FP={fp_bench}, FN={fn_bench})")

sim_df = pd.DataFrame(sim_rows)
sim_df.to_csv(OUT_DIR/"step4_full_benchmark_simulation.tsv", sep="\t", index=False)

print("\n  Filter simulation results (correct F1 with FN):")
for sample in PAIRED_SAMPLES:
    sdf = sim_df[sim_df["sample"]==sample]
    if len(sdf) == 0: continue
    print(f"\n  {sample}:")
    for _, row in sdf.iterrows():
        print(f"    {row['filter']}: rm_TP={row['rm_tp']}({100*row['rm_tp_pct']:.1f}%), "
              f"rm_FP={row['rm_fp']}({100*row['rm_fp_pct']:.1f}%), "
              f"F1={row['f1_new']:.4f} (delta={row['f1_delta']:+.4f})")

# ============================================================
# Step 5: A≥30 vs 30≤A<50 vs A≥50 雙層 Tier 分析
# ============================================================
print("\n[6/8] Step 5: 雙層 Tier 分析（A vs A+）...")

def compute_enrichment_tier(subdf, tier_mask, mode):
    sub = subdf[subdf["mode"]==mode]
    tp_total = (sub["truth_label"]=="TP").sum()
    fp_total = (sub["truth_label"]=="FP").sum()
    sub_tier = sub[tier_mask(sub)]
    tp_loh = (sub_tier[sub_tier["truth_label"]=="TP"]["core_loh_like"]).sum()
    fp_loh = (sub_tier[sub_tier["truth_label"]=="FP"]["core_loh_like"]).sum()
    tp_f = tp_loh / tp_total if tp_total > 0 else np.nan
    fp_f = fp_loh / fp_total if fp_total > 0 else np.nan
    enr = fp_f / tp_f if (tp_f and tp_f > 0) else np.nan
    a, b, c, d = int(fp_loh), int(fp_total - fp_loh), int(tp_loh), int(tp_total - tp_loh)
    try:
        _, pval = stats.fisher_exact([[a,b],[c,d]])
    except:
        pval = np.nan
    return {"mode": mode, "tp_total": tp_total, "fp_total": fp_total,
            "tp_loh": tp_loh, "fp_loh": fp_loh, "tp_frac": tp_f, "fp_frac": fp_f,
            "enrichment": enr, "fisher_p": pval}

tier_rows = []
for tier_label, tier_fn in [
    ("A(30-49)", lambda d: (d["tier_dual"]=="A")),
    ("A+(>=50)", lambda d: (d["tier_dual"]=="A+")),
    ("A_combined(>=30)", lambda d: (d["tier_dual"].isin(["A","A+"]))),
]:
    for mode in ["paired","to"]:
        row = compute_enrichment_tier(df, tier_fn, mode)
        row["tier_label"] = tier_label
        tier_rows.append(row)

tier_df = pd.DataFrame(tier_rows)
tier_df.to_csv(OUT_DIR/"step5_dual_tier_enrichment.tsv", sep="\t", index=False)

print("\n  Dual-tier enrichment:")
for _, row in tier_df.iterrows():
    print(f"  [{row['mode']}] {row['tier_label']}: enrichment={row['enrichment']:.4f}, "
          f"p={row['fisher_p']:.2e}, "
          f"TP_LOH_n={row['tp_loh']}, FP_LOH_n={row['fp_loh']}")

# ============================================================
# 視覺化
# ============================================================
print("\n[7/8] 視覺化 ...")
sns.set_style("whitegrid")
PALETTE_PAIRED = "#2166ac"
PALETTE_TO     = "#d73027"

# --- Fig 01: Baseline F1 per sample（含 TP/FP/FN bar + F1 點）---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax_left, ax_right = axes

samples_sorted = baseline_df.sort_values("f1_true", ascending=False)["sample"].tolist()
x = np.arange(len(samples_sorted))
w = 0.25
for i, (col, label, color) in enumerate([
    ("tp", "TP", "#4575b4"),
    ("fp", "FP", "#d73027"),
    ("fn", "FN", "#fdae61"),
]):
    vals = [baseline_df[baseline_df["sample"]==s][col].values[0]
            for s in samples_sorted]
    ax_left.bar(x + i*w, vals, width=w, label=label, color=color, alpha=0.8)

ax_left.set_xticks(x + w)
ax_left.set_xticklabels(samples_sorted, rotation=30, ha="right", fontsize=9)
ax_left.set_ylabel("Count")
ax_left.set_title("Complete Benchmark: TP / FP / FN per Sample\n(LongPhase-S / InterSubMod baseline)")
ax_left.legend()
ax_left.set_yscale("log")

f1_vals = [baseline_df[baseline_df["sample"]==s]["f1_true"].values[0]
           for s in samples_sorted]
ax_right.barh(samples_sorted, f1_vals, color=PALETTE_PAIRED, alpha=0.8)
ax_right.axvline(0.9, color="red", linestyle="--", alpha=0.5, label="F1=0.9")
for i, (s, v) in enumerate(zip(samples_sorted, f1_vals)):
    ax_right.text(v + 0.002, i, f"{v:.4f}", va="center", fontsize=8)
ax_right.set_xlabel("F1 (with FN)")
ax_right.set_title("True Baseline F1 (including FN)")
ax_right.set_xlim(0.8, 1.01)
ax_right.legend()
plt.tight_layout()
fig.savefig(FIG_DIR/"fig01_baseline_f1.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 01 saved")

# --- Fig 02: 80 FPs per-sample decomposition (stacked bar) ---
fig, ax = plt.subplots(figsize=(10, 5))
samples_ord = sample_decomp.sort_values("fp_count", ascending=False)["sample"].tolist()
tp_vals = [sample_decomp[sample_decomp["sample"]==s]["tp_count"].values[0]
           if s in sample_decomp["sample"].values else 0 for s in samples_ord]
fp_vals = [sample_decomp[sample_decomp["sample"]==s]["fp_count"].values[0]
           if s in sample_decomp["sample"].values else 0 for s in samples_ord]

x = np.arange(len(samples_ord))
ax.bar(x, tp_vals, label="TP (LOH=Y+HPSig=Y)", color="#4575b4", alpha=0.85)
ax.bar(x, fp_vals, bottom=tp_vals, label="FP (LOH=Y+HPSig=Y)", color="#d73027", alpha=0.85)

for i, (tp, fp) in enumerate(zip(tp_vals, fp_vals)):
    total = tp + fp
    if total > 0:
        fp_pct = 100*fp/total
        ax.text(i, total + 1, f"FP%={fp_pct:.0f}%", ha="center", fontsize=7.5, color="black")

ax.set_xticks(x)
ax.set_xticklabels(samples_ord, rotation=30, ha="right", fontsize=9)
ax.set_ylabel("Count (Tier A LOH=Y + HPMergedSig=Y)")
ax.set_title("Paired Tier A LOH=Y+HPMergedSig=Y: Per-Sample TP/FP Breakdown\n(7.4× FP enrichment vs LOH alone)")
ax.legend()
plt.tight_layout()
fig.savefig(FIG_DIR/"fig02_fp80_per_sample.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 02 saved")

# --- Fig 03: HPMergedDelta 分佈（TP vs FP, Tier A LOH-like） ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (loh, sig, title) in zip(axes, [
    (True, True, "LOH=Y + HPMergedSig=Y\n(7.4× FP enrichment group)"),
    (True, False, "LOH=Y + HPMergedSig=N\n(baseline reference)"),
]):
    sub = df[(df["mode"]=="paired") & (df["hp_support_tier"]=="A") &
             (df["core_loh_like"]==loh) & (df["HPMergedSig"]==sig)]
    for truth, color in [("TP","#4575b4"),("FP","#d73027")]:
        g = sub[sub["truth_label"]==truth]["HPMergedDelta"].abs().dropna()
        if len(g) > 0:
            ax.hist(g, bins=40, alpha=0.55, color=color,
                    label=f"{truth} (n={len(g):,})", density=True)
    ax.set_xlabel("|HPMergedDelta|")
    ax.set_ylabel("Density")
    ax.set_title(title)
    ax.legend(fontsize=8)
plt.suptitle("HPMergedDelta Distribution: TP vs FP (Paired Tier A LOH-like)", y=1.02)
plt.tight_layout()
fig.savefig(FIG_DIR/"fig03_hp_merged_delta_dist.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 03 saved")

# --- Fig 04: Chr distribution of 80 FPs ---
if "Chr" in fp_joint_paired.columns and len(fp_joint_paired) > 0:
    fig, ax = plt.subplots(figsize=(12, 4))
    chr_order = [f"chr{i}" for i in list(range(1,23))+["X","Y"]]
    chr_data = fp_joint_paired["Chr"].value_counts().reindex(chr_order).fillna(0)
    ax.bar(range(len(chr_order)), chr_data.values, color="#fdae61", alpha=0.85)
    ax.set_xticks(range(len(chr_order)))
    ax.set_xticklabels(chr_order, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Number of FPs (LOH=Y+HPSig=Y)")
    ax.set_title(f"Chromosomal Distribution of {len(fp_joint_paired)} Paired Tier A LOH+HPSig FPs")
    plt.tight_layout()
    fig.savefig(FIG_DIR/"fig04_chr_distribution_fp.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Fig 04 saved")

# --- Fig 05: 2×2 LOH×HPSig FP% heatmap (paired vs TO) ---
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, mode in zip(axes, ["paired","to"]):
    sub = joint_df[joint_df["mode"]==mode].copy()
    heat = sub.pivot_table(index="loh_like", columns="hp_merged_sig",
                           values="fp_frac", aggfunc="first")
    heat.index = [f"LOH={'Y' if v else 'N'}" for v in heat.index]
    heat.columns = [f"HPSig={'Y' if v else 'N'}" for v in heat.columns]
    sns.heatmap(heat, ax=ax, annot=True, fmt=".3f", cmap="YlOrRd",
                vmin=0, vmax=heat.values.max()*1.1, linewidths=0.5)
    ax.set_title(f"{mode.upper()}: FP Rate (Tier A)")
plt.suptitle("Joint LOH × HPMergedSig FP Rate Matrix", y=1.02)
plt.tight_layout()
fig.savefig(FIG_DIR/"fig05_joint_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 05 saved")

# --- Fig 06: Full benchmark filter simulation F1 delta ---
fig, ax = plt.subplots(figsize=(12, 5))
sim_df_ordered = sim_df[sim_df["filter"].isin(["TierA_LOH","TierA_LOH_HPSig","TierAplus_LOH","TierAplus_LOH_HPSig"])]
pivot = sim_df_ordered.pivot_table(index="sample", columns="filter", values="f1_delta", aggfunc="first")
pivot = pivot.reindex(PAIRED_SAMPLES)

x = np.arange(len(PAIRED_SAMPLES))
w = 0.22
colors = {"TierA_LOH":"#d73027","TierA_LOH_HPSig":"#fc8d59",
          "TierAplus_LOH":"#4575b4","TierAplus_LOH_HPSig":"#74add1"}
for i, filt in enumerate(["TierA_LOH","TierA_LOH_HPSig","TierAplus_LOH","TierAplus_LOH_HPSig"]):
    if filt in pivot.columns:
        vals = [pivot.loc[s, filt] if s in pivot.index else 0 for s in PAIRED_SAMPLES]
        ax.bar(x + i*w, vals, width=w, label=filt, color=colors[filt], alpha=0.85)

ax.axhline(0, color="black", linewidth=0.8)
ax.set_xticks(x + 1.5*w)
ax.set_xticklabels(PAIRED_SAMPLES, rotation=30, ha="right", fontsize=9)
ax.set_ylabel("F1 Delta (correct F1 with FN)")
ax.set_title("Full Benchmark Filter Simulation: F1 Delta\n(Correct calculation including FN)")
ax.legend(fontsize=8, ncol=2)
plt.tight_layout()
fig.savefig(FIG_DIR/"fig06_full_benchmark_f1_delta.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 06 saved")

# --- Fig 07: Dual tier enrichment ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, mode in zip(axes, ["paired","to"]):
    sub = tier_df[tier_df["mode"]==mode]
    tiers = sub["tier_label"].tolist()
    enr_vals = sub["enrichment"].tolist()
    colors_bar = ["#74add1" if e < 1 else "#d73027" for e in enr_vals]
    bars = ax.barh(tiers, enr_vals, color=colors_bar, alpha=0.8)
    ax.axvline(1.0, color="black", linewidth=1.2, linestyle="--", label="enrichment=1 (no effect)")
    for bar, row in zip(bars, sub.itertuples()):
        ax.text(row.enrichment + 0.01, bar.get_y() + bar.get_height()/2,
                f"{row.enrichment:.3f} (p={row.fisher_p:.1e})", va="center", fontsize=8)
    ax.set_xlabel("FP/TP LOH-like Enrichment")
    ax.set_title(f"{mode.upper()} — Tier A vs A+ Enrichment")
    ax.legend()
plt.suptitle("Dual Tier Definition: A(30-49) vs A+(≥50) LOH-like Enrichment", y=1.02)
plt.tight_layout()
fig.savefig(FIG_DIR/"fig07_dual_tier_enrichment.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 07 saved")

# --- Fig 08: F1 baseline comparison (Round 3 wrong vs Round 4 correct) ---
fig, ax = plt.subplots(figsize=(10, 5))
r3_f1 = {
    "H2009":0.9997,"HCC1937":0.9922,"HCC1954":0.9992,
    "COLO829":0.9689,"HCC1395":0.9896,"HCC1395_DORADO":0.9960,"H1437":0.9999
}
r4_f1 = {r["sample"]: r["f1_true"] for _, r in baseline_df.iterrows()}
samples_fig8 = [s for s in PAIRED_SAMPLES if s in r4_f1]
x = np.arange(len(samples_fig8))
w = 0.38
bars1 = ax.bar(x - w/2, [r3_f1.get(s, np.nan) for s in samples_fig8],
               width=w, label="Round 3 F1 (no FN, subset only)", color="#fdae61", alpha=0.85)
bars2 = ax.bar(x + w/2, [r4_f1.get(s, np.nan) for s in samples_fig8],
               width=w, label="Round 4 F1 (true F1, with FN)", color="#4575b4", alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(samples_fig8, rotation=30, ha="right", fontsize=9)
ax.set_ylabel("F1")
ax.set_ylim(0.8, 1.01)
ax.axhline(1.0, color="gray", linewidth=0.5, linestyle="--")
ax.set_title("Round 3 vs Round 4 F1 Baseline Comparison\n(Round 3 was overestimated due to missing FN)")
ax.legend()
plt.tight_layout()
fig.savefig(FIG_DIR/"fig08_f1_baseline_comparison.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 08 saved")

# ============================================================
# 彙整 Key Numbers
# ============================================================
print("\n[8/8] 彙整 key numbers ...")
print("\n=== Round 4 Key Numbers ===")

print("\nStep 1 — Complete Benchmark Baseline:")
for _, row in baseline_df.iterrows():
    print(f"  {row['sample']} [{row['method']}]: TP={row['tp']:,}, FP={row['fp']:,}, FN={row['fn']:,}, F1={row['f1_true']:.4f}")

print("\nStep 2 — LOH Data Coverage Ratio:")
for _, row in coverage_df.iterrows():
    print(f"  {row['sample']}: TP coverage={row['tp_coverage']:.4f}, FP coverage={row['fp_coverage']:.4f}")

print("\nStep 3 — LOH+HPMergedSig Deep Validation:")
print(f"  Paired Tier A LOH=Y+HPSig=Y total: {len(fp_joint_paired)} FPs")
print("  Per-sample FP decomposition:")
for _, row in sample_decomp.iterrows():
    print(f"    {row['sample']}: FP={row['fp_count']}, TP={row['tp_count']}, FP%={100*row['fp_frac']:.1f}%")

print("\nStep 4 — Full Benchmark Filter Simulation (with FN):")
for sample in PAIRED_SAMPLES:
    sdf = sim_df[sim_df["sample"]==sample]
    if len(sdf) == 0: continue
    row_hpsig = sdf[sdf["filter"]=="TierA_LOH_HPSig"]
    row_aplus = sdf[sdf["filter"]=="TierAplus_LOH_HPSig"]
    if len(row_hpsig) > 0:
        row = row_hpsig.iloc[0]
        print(f"  {sample}/TierA_LOH_HPSig: F1 {row['f1_baseline']:.4f}→{row['f1_new']:.4f} (delta={row['f1_delta']:+.4f})")
    if len(row_aplus) > 0:
        row = row_aplus.iloc[0]
        print(f"  {sample}/TierAplus_LOH_HPSig: F1 {row['f1_baseline']:.4f}→{row['f1_new']:.4f} (delta={row['f1_delta']:+.4f})")

print("\nStep 5 — Dual Tier Enrichment:")
for _, row in tier_df.iterrows():
    print(f"  [{row['mode']}] {row['tier_label']}: enrichment={row['enrichment']:.4f} (p={row['fisher_p']:.2e})")

print(f"\n完成！輸出位置：{OUT_DIR}")
