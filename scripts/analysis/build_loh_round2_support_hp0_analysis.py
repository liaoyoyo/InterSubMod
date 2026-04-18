#!/usr/bin/env python3
"""
LOH Round 2 Analysis: effective_hp Support Stratification + HP0 Source Analysis

研究問題：在什麼品質條件下，LOH evidence 才是可信的 evidence panel 成員？

Core 1: effective_hp support 分層分析
  - Tier A (>=30), B (10-29), C (<10, including C0=0)
  - 各 tier 的 TP/FP LOH-like enrichment

Core 2: HP0 來源分析
  - LOH-like vs non-LOH-like 的 HP0 ratio 分佈
  - UnassignedAffinity (Stage 4) 的分佈
  - 同位點 paired vs TO HP0 比較
  - HP0 來源假說驗證

生成時間: 2026-03-27
"""

import os

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
import json
from datetime import datetime
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

# ============================================================
# 設定
# ============================================================
ROUND1_WS = Path(
    os.environ.get(
        "LOH_ROUND1_WS",
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit",
    )
)
OUT_DIR = Path(
    os.environ.get(
        "LOH_ROUND2_OUT_DIR",
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round2_support_hp0_analysis",
    )
)
FIG_DIR = OUT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

GENERATED_AT = datetime.now().strftime("%Y-%m-%d %H:%M")

# Tier 定義
def assign_tier(n):
    """Assign support quality tier based on effective_hp_reads."""
    if pd.isna(n):
        return "C0"
    n = int(n)
    if n == 0:
        return "C0"   # 完全沒有 HP reads（HP_Ratio=0.5 假象區）
    elif n < 10:
        return "C"    # 極弱 support（Round 1 確認為雜訊）
    elif n < 30:
        return "B"    # 中等 support（謹慎解讀）
    else:
        return "A"    # 強 support（可信 LOH evidence）

TIER_ORDER = ["C0", "C", "B", "A"]
TIER_COLORS = {"C0": "#d9d9d9", "C": "#fc8d59", "B": "#fee090", "A": "#4575b4"}

# 樣本顯示名稱
SAMPLE_LABEL = {
    "HCC1395": "HCC1395_5kHz",
    "HCC1395_DORADO": "HCC1395_DORADO",
    "COLO829": "COLO829",
    "H1437": "H1437",
    "H2009": "H2009",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
}

# ============================================================
# 資料載入
# ============================================================
print("[1/8] 載入 Round 1 all_region_rows.tsv.gz ...")
df = pd.read_csv(ROUND1_WS / "all_region_rows.tsv.gz", sep="\t", low_memory=False)
print(f"  總計 {len(df):,} rows, {len(df.columns)} columns")

# 型別整理
df["effective_hp_reads"] = pd.to_numeric(df["effective_hp_reads"], errors="coerce")
df["hp_ratio_core"] = pd.to_numeric(df["hp_ratio_core"], errors="coerce")
df["hp0_ratio"] = pd.to_numeric(df["hp0_ratio"], errors="coerce")
df["hp3_ratio"] = pd.to_numeric(df["hp3_ratio"], errors="coerce")
df["NHP0"] = pd.to_numeric(df["NHP0"], errors="coerce")
df["NHP3"] = pd.to_numeric(df["NHP3"], errors="coerce")
df["NumReads"] = pd.to_numeric(df["NumReads"], errors="coerce")
df["UnassignedAffinity"] = pd.to_numeric(df["UnassignedAffinity"], errors="coerce")

# Tier 標記
df["hp_support_tier"] = df["effective_hp_reads"].apply(assign_tier)

# 確認 truth_label 只保留 TP/FP
df_tf = df[df["truth_label"].isin(["TP", "FP"])].copy()
print(f"  TP/FP rows: {len(df_tf):,} (TP={df_tf['truth_label'].eq('TP').sum():,}, FP={df_tf['truth_label'].eq('FP').sum():,})")

# ============================================================
# Core 1: Effective HP Support 分層分析
# ============================================================
print("\n[2/8] Core 1: effective_hp support 分層分析 ...")

def compute_enrichment(group):
    """計算 LOH-like TP/FP enrichment 及相關統計。"""
    tp_mask = group["truth_label"] == "TP"
    fp_mask = group["truth_label"] == "FP"
    tp_df = group[tp_mask]
    fp_df = group[fp_mask]

    tp_total = len(tp_df)
    fp_total = len(fp_df)
    tp_loh_n = tp_df["core_loh_like"].sum() if tp_total > 0 else 0
    fp_loh_n = fp_df["core_loh_like"].sum() if fp_total > 0 else 0

    tp_loh_frac = tp_loh_n / tp_total if tp_total > 0 else np.nan
    fp_loh_frac = fp_loh_n / fp_total if fp_total > 0 else np.nan

    if tp_loh_frac and tp_loh_frac > 0:
        enrichment = fp_loh_frac / tp_loh_frac
    else:
        enrichment = np.nan

    # Fisher's exact test (LOH-like vs non in TP vs FP)
    if tp_total > 0 and fp_total > 0:
        a = int(fp_loh_n)
        b = int(fp_total - fp_loh_n)
        c = int(tp_loh_n)
        d = int(tp_total - tp_loh_n)
        try:
            _, pval = stats.fisher_exact([[a, b], [c, d]])
        except Exception:
            pval = np.nan
    else:
        pval = np.nan

    return pd.Series({
        "tp_total": tp_total,
        "fp_total": fp_total,
        "tp_loh_like_n": tp_loh_n,
        "fp_loh_like_n": fp_loh_n,
        "tp_loh_like_frac": tp_loh_frac,
        "fp_loh_like_frac": fp_loh_frac,
        "fp_tp_enrichment": enrichment,
        "fisher_pval": pval,
    })

# (a) 全域：by mode × tier
tier_global = df_tf.groupby(["mode", "hp_support_tier"]).apply(compute_enrichment).reset_index()
tier_global.to_csv(OUT_DIR / "core1_tier_enrichment_global.tsv", sep="\t", index=False)
print(f"  核心結果：全域 tier enrichment ({len(tier_global)} rows)")

# (b) by sample × mode × tier
tier_by_sample = df_tf.groupby(["sample", "mode", "hp_support_tier"]).apply(compute_enrichment).reset_index()
tier_by_sample.to_csv(OUT_DIR / "core1_tier_enrichment_by_sample_mode.tsv", sep="\t", index=False)
print(f"  全樣本分層 ({len(tier_by_sample)} rows)")

# ============================================================
# Core 1 視覺化
# ============================================================
print("\n[3/8] Core 1 視覺化 ...")

# Fig 1: 全域 Tier enrichment（paired vs TO 各 tier 的 FP/TP enrichment）
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, mode_name in zip(axes, ["paired", "to"]):
    sub = tier_global[tier_global["mode"] == mode_name].copy()
    sub["hp_support_tier"] = pd.Categorical(sub["hp_support_tier"], categories=TIER_ORDER, ordered=True)
    sub = sub.sort_values("hp_support_tier")
    colors = [TIER_COLORS[t] for t in sub["hp_support_tier"]]
    bars = ax.bar(sub["hp_support_tier"], sub["fp_tp_enrichment"],
                  color=colors, edgecolor="k", linewidth=0.8)
    ax.axhline(1.0, color="red", linestyle="--", linewidth=1.2, label="Enrichment=1.0（無差異）")
    ax.set_title(f"Mode: {mode_name.upper()}\nLOH-like FP/TP Enrichment by Support Tier", fontsize=11)
    ax.set_xlabel("HP Support Tier", fontsize=10)
    ax.set_ylabel("FP/TP LOH-like Enrichment", fontsize=10)
    ax.legend(fontsize=8)
    # 標注數值
    for bar, row in zip(bars, sub.itertuples()):
        if not np.isnan(row.fp_tp_enrichment):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                    f"{row.fp_tp_enrichment:.3f}\n(TP={int(row.tp_total):,})",
                    ha="center", va="bottom", fontsize=7)
    ax.set_ylim(0, max(3.0, sub["fp_tp_enrichment"].dropna().max() * 1.3))

plt.tight_layout()
plt.savefig(FIG_DIR / "fig01_tier_enrichment_global.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 01 saved")

# Fig 2: 各樣本 Tier A enrichment（paired mode）
paired_tier = tier_by_sample[tier_by_sample["mode"] == "paired"].copy()
pivot_enrich = paired_tier.pivot_table(
    index="sample", columns="hp_support_tier", values="fp_tp_enrichment"
)
# 確保 Tier 順序
for col in TIER_ORDER:
    if col not in pivot_enrich.columns:
        pivot_enrich[col] = np.nan
pivot_enrich = pivot_enrich[TIER_ORDER]

fig, ax = plt.subplots(figsize=(10, 7))
im = ax.imshow(pivot_enrich.values, aspect="auto", cmap="RdYlGn_r", vmin=0, vmax=4)
ax.set_xticks(range(len(TIER_ORDER)))
ax.set_xticklabels(TIER_ORDER, fontsize=10)
ax.set_yticks(range(len(pivot_enrich.index)))
ax.set_yticklabels(pivot_enrich.index, fontsize=9)
ax.set_xlabel("HP Support Tier", fontsize=11)
ax.set_ylabel("Sample", fontsize=11)
ax.set_title("Paired Mode: LOH-like FP/TP Enrichment\nby Sample × Support Tier", fontsize=12)
for i in range(len(pivot_enrich.index)):
    for j in range(len(TIER_ORDER)):
        val = pivot_enrich.values[i, j]
        if not np.isnan(val):
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=9,
                    color="white" if val > 2.0 else "black")
plt.colorbar(im, ax=ax, label="FP/TP Enrichment")
plt.tight_layout()
plt.savefig(FIG_DIR / "fig02_paired_tier_enrichment_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 02 saved")

# Fig 3: LOH-like fraction 各 tier（TP vs FP，paired + TO）
fig, axes = plt.subplots(2, 1, figsize=(14, 10))
for ax, mode_name in zip(axes, ["paired", "to"]):
    sub = tier_global[tier_global["mode"] == mode_name].copy()
    sub["hp_support_tier"] = pd.Categorical(sub["hp_support_tier"], categories=TIER_ORDER, ordered=True)
    sub = sub.sort_values("hp_support_tier")
    x = np.arange(len(sub))
    width = 0.35
    ax.bar(x - width/2, sub["tp_loh_like_frac"], width, label="TP LOH-like frac", color="#2166ac", alpha=0.8)
    ax.bar(x + width/2, sub["fp_loh_like_frac"], width, label="FP LOH-like frac", color="#d73027", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(sub["hp_support_tier"])
    ax.set_ylabel("LOH-like Fraction", fontsize=10)
    ax.set_title(f"Mode: {mode_name.upper()} — LOH-like Fraction by Support Tier (TP vs FP)", fontsize=11)
    ax.legend()
    ax.set_ylim(0, 1.0)
    # 標注 N
    for xi, row in zip(x, sub.itertuples()):
        ax.text(xi, -0.05, f"TP:{int(row.tp_total):,}\nFP:{int(row.fp_total):,}",
                ha="center", va="top", fontsize=6.5, transform=ax.transData)

plt.tight_layout()
plt.savefig(FIG_DIR / "fig03_tier_loh_frac_barplot.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 03 saved")

# ============================================================
# Core 2: HP0 來源分析
# ============================================================
print("\n[4/8] Core 2: HP0 來源分析 ...")

# (a) 基本統計：LOH-like vs non-LOH-like 的 hp0_ratio
def hp0_stats(group):
    return pd.Series({
        "n": len(group),
        "hp0_ratio_mean": group["hp0_ratio"].mean(),
        "hp0_ratio_median": group["hp0_ratio"].median(),
        "hp0_ratio_q25": group["hp0_ratio"].quantile(0.25),
        "hp0_ratio_q75": group["hp0_ratio"].quantile(0.75),
        "hp0_ratio_gt20pct": (group["hp0_ratio"] > 0.2).mean(),
        "hp0_ratio_gt40pct": (group["hp0_ratio"] > 0.4).mean(),
    })

hp0_by_loh = df_tf.groupby(["mode", "truth_label", "core_loh_like"]).apply(hp0_stats).reset_index()
hp0_by_loh.to_csv(OUT_DIR / "core2_hp0_by_loh_status.tsv", sep="\t", index=False)

# (b) by sample, mode, truth, loh_like
hp0_by_sample = df_tf.groupby(["sample", "mode", "truth_label", "core_loh_like"]).apply(hp0_stats).reset_index()
hp0_by_sample.to_csv(OUT_DIR / "core2_hp0_by_sample_loh_status.tsv", sep="\t", index=False)
print(f"  HP0 by LOH status computed ({len(hp0_by_loh)} global groups, {len(hp0_by_sample)} sample groups)")

# (c) Mann-Whitney U test: hp0_ratio in LOH-like vs non-LOH-like
print("\n  HP0 ratio Mann-Whitney test (LOH-like vs non-LOH-like):")
mw_results = []
for (mode, truth), grp in df_tf.groupby(["mode", "truth_label"]):
    loh_vals = grp[grp["core_loh_like"] == True]["hp0_ratio"].dropna()
    non_loh_vals = grp[grp["core_loh_like"] == False]["hp0_ratio"].dropna()
    if len(loh_vals) > 10 and len(non_loh_vals) > 10:
        stat, pval = stats.mannwhitneyu(loh_vals, non_loh_vals, alternative="two-sided")
        median_diff = loh_vals.median() - non_loh_vals.median()
        print(f"    {mode}/{truth}: LOH-like median={loh_vals.median():.4f}, non-LOH median={non_loh_vals.median():.4f}, diff={median_diff:+.4f}, p={pval:.2e}")
        mw_results.append({"mode": mode, "truth": truth,
                           "loh_like_median_hp0": loh_vals.median(),
                           "non_loh_like_median_hp0": non_loh_vals.median(),
                           "median_diff": median_diff,
                           "mw_pval": pval,
                           "n_loh": len(loh_vals), "n_non_loh": len(non_loh_vals)})

pd.DataFrame(mw_results).to_csv(OUT_DIR / "core2_hp0_loh_mannwhitney.tsv", sep="\t", index=False)

# (d) UnassignedAffinity 分析
print("\n  UnassignedAffinity 分析 ...")

# 只看有效的 affinity（非 NaN）
df_affinity = df_tf[df_tf["UnassignedAffinity"].notna()].copy()
print(f"  有效 affinity rows: {len(df_affinity):,} ({len(df_affinity)/len(df_tf)*100:.1f}% of TP/FP)")

def affinity_stats(group):
    vals = group["UnassignedAffinity"]
    return pd.Series({
        "n_valid_affinity": len(vals),
        "affinity_mean": vals.mean(),
        "affinity_median": vals.median(),
        "affinity_q25": vals.quantile(0.25),
        "affinity_q75": vals.quantile(0.75),
        "affinity_abs_mean": vals.abs().mean(),  # 傾向強度
        "frac_toward_hp2": (vals > 0.2).mean(),   # 偏 HP2
        "frac_toward_hp1": (vals < -0.2).mean(),  # 偏 HP1
        "frac_neutral": (vals.abs() <= 0.2).mean(), # 中性
    })

affinity_by_loh = df_affinity.groupby(["mode", "truth_label", "core_loh_like"]).apply(affinity_stats).reset_index()
affinity_by_loh.to_csv(OUT_DIR / "core2_affinity_by_loh_status.tsv", sep="\t", index=False)

# Mann-Whitney: affinity abs strength in LOH-like vs non-LOH-like
print("\n  Affinity abs strength test:")
affinity_mw = []
for (mode, truth), grp in df_affinity.groupby(["mode", "truth_label"]):
    loh_vals = grp[grp["core_loh_like"] == True]["UnassignedAffinity"].abs()
    non_loh_vals = grp[grp["core_loh_like"] == False]["UnassignedAffinity"].abs()
    if len(loh_vals) > 10 and len(non_loh_vals) > 10:
        stat, pval = stats.mannwhitneyu(loh_vals, non_loh_vals, alternative="two-sided")
        print(f"    {mode}/{truth}: LOH-like abs_median={loh_vals.median():.4f}, non-LOH abs_median={non_loh_vals.median():.4f}, p={pval:.2e}")
        affinity_mw.append({"mode": mode, "truth": truth,
                            "loh_like_abs_affinity_median": loh_vals.median(),
                            "non_loh_like_abs_affinity_median": non_loh_vals.median(),
                            "mw_pval": pval, "n_loh": len(loh_vals), "n_non_loh": len(non_loh_vals)})

pd.DataFrame(affinity_mw).to_csv(OUT_DIR / "core2_affinity_loh_mannwhitney.tsv", sep="\t", index=False)

# ============================================================
# Core 2: HP0 來源假說驗證
# ============================================================
print("\n[5/8] HP0 來源假說驗證 ...")

# 假說定義:
# H1: HP0 高 ↔ NumReads 低（coverage 不足導致 phasing 困難）
# H2: HP0 高 ↔ core_loh_like = True（LOH 導致 phasing 困難）
# H3 (TO only): HP0 在 TO 中系統性高於 paired 同位點

# H1: Correlation hp0_ratio vs NumReads
print("  H1: HP0 vs Coverage correlation:")
h1_results = []
for (mode, truth), grp in df_tf.groupby(["mode", "truth_label"]):
    valid = grp[["hp0_ratio", "NumReads"]].dropna()
    if len(valid) > 30:
        rho, pval = stats.spearmanr(valid["hp0_ratio"], valid["NumReads"])
        print(f"    {mode}/{truth}: Spearman rho={rho:.4f}, p={pval:.2e} (n={len(valid):,})")
        h1_results.append({"mode": mode, "truth": truth, "spearman_rho": rho, "pval": pval, "n": len(valid)})

pd.DataFrame(h1_results).to_csv(OUT_DIR / "core2_hp0_coverage_correlation.tsv", sep="\t", index=False)

# H2: HP0 in LOH-like vs non-LOH-like（已在上方完成，這裡做更細的分析）
print("  H2: HP0 LOH-like 關聯：已在 Mann-Whitney 中完成")

# H3: Same-locus HP0 comparison (paired vs TO)
print("\n  H3: Same-locus paired vs TO HP0 比較 ...")
slc = pd.read_csv(ROUND1_WS / "same_locus_compare.tsv", sep="\t", low_memory=False)
slc["hp0_ratio_paired"] = pd.to_numeric(slc.get("hp0_ratio_paired", np.nan), errors="coerce")
slc["hp0_ratio_to"] = pd.to_numeric(slc.get("hp0_ratio_to", np.nan), errors="coerce")
slc["effective_hp_reads_paired"] = pd.to_numeric(slc.get("effective_hp_reads_paired", np.nan), errors="coerce")
slc["effective_hp_reads_to"] = pd.to_numeric(slc.get("effective_hp_reads_to", np.nan), errors="coerce")

# 只看 both_present
slc_both = slc[slc["locus_overlap_class"] == "both"].copy()
slc_both = slc_both[slc_both["hp0_ratio_paired"].notna() & slc_both["hp0_ratio_to"].notna()]
slc_both["hp0_delta"] = slc_both["hp0_ratio_to"] - slc_both["hp0_ratio_paired"]
print(f"  Same-locus both present: {len(slc_both):,} rows")

h3_results = []
for sample_name, grp in slc_both.groupby("sample"):
    if len(grp) < 10:
        continue
    stat, pval = stats.wilcoxon(grp["hp0_ratio_to"], grp["hp0_ratio_paired"], alternative="greater")
    median_delta = grp["hp0_delta"].median()
    print(f"    {sample_name}: TO hp0 > paired hp0? delta_median={median_delta:+.4f}, Wilcoxon p={pval:.2e} (n={len(grp):,})")
    h3_results.append({"sample": sample_name, "hp0_delta_median": median_delta,
                       "wilcoxon_pval": pval, "n": len(grp)})

pd.DataFrame(h3_results).to_csv(OUT_DIR / "core2_hp0_same_locus_delta.tsv", sep="\t", index=False)

# HP0 same-locus by concordance group
print("\n  HP0 by concordance group:")
for conc, grp in slc_both.groupby("paired_to_concordance"):
    print(f"    {conc}: n={len(grp):,}, TO hp0 median={grp['hp0_ratio_to'].median():.4f}, paired hp0 median={grp['hp0_ratio_paired'].median():.4f}")

hp0_by_concordance = slc_both.groupby(["sample", "paired_to_concordance"]).agg(
    n=("hp0_delta", "count"),
    hp0_paired_median=("hp0_ratio_paired", "median"),
    hp0_to_median=("hp0_ratio_to", "median"),
    hp0_delta_median=("hp0_delta", "median"),
    hp0_delta_q25=("hp0_delta", lambda x: x.quantile(0.25)),
    hp0_delta_q75=("hp0_delta", lambda x: x.quantile(0.75)),
).reset_index()
hp0_by_concordance.to_csv(OUT_DIR / "core2_hp0_by_concordance.tsv", sep="\t", index=False)

# ============================================================
# Core 2 視覺化
# ============================================================
print("\n[6/8] Core 2 視覺化 ...")

# Fig 4: HP0 ratio 分佈 violin plot（LOH-like vs non-LOH-like, TP vs FP）
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
for row_idx, mode_name in enumerate(["paired", "to"]):
    for col_idx, truth_name in enumerate(["TP", "FP"]):
        ax = axes[row_idx][col_idx]
        sub = df_tf[(df_tf["mode"] == mode_name) & (df_tf["truth_label"] == truth_name)]
        if len(sub) < 10:
            ax.set_visible(False)
            continue
        loh_groups = sub["core_loh_like"].map({True: "LOH-like", False: "non-LOH-like"})
        parts = ax.violinplot(
            [sub[loh_groups == "non-LOH-like"]["hp0_ratio"].dropna(),
             sub[loh_groups == "LOH-like"]["hp0_ratio"].dropna()],
            positions=[1, 2], widths=0.6, showmedians=True
        )
        for pc, color in zip(parts["bodies"], ["#91bfdb", "#fc8d59"]):
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["non-LOH-like", "LOH-like"])
        ax.set_ylabel("HP0 Ratio")
        ax.set_title(f"{mode_name.upper()} / {truth_name}\n(n={len(sub):,})")
        # 標注 median
        for x_pos, loh_flag in [(1, False), (2, True)]:
            vals = sub[sub["core_loh_like"] == loh_flag]["hp0_ratio"].dropna()
            if len(vals) > 0:
                ax.text(x_pos, vals.median() + 0.02, f"med={vals.median():.3f}", ha="center", fontsize=8)

fig.suptitle("HP0 Ratio Distribution: LOH-like vs non-LOH-like\n(by mode and truth label)", fontsize=13)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig04_hp0_distribution_violin.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 04 saved")

# Fig 5: UnassignedAffinity 分佈（LOH-like vs non-LOH-like）
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
for row_idx, mode_name in enumerate(["paired", "to"]):
    for col_idx, truth_name in enumerate(["TP", "FP"]):
        ax = axes[row_idx][col_idx]
        sub = df_affinity[(df_affinity["mode"] == mode_name) & (df_affinity["truth_label"] == truth_name)]
        if len(sub) < 10:
            ax.set_visible(False)
            continue
        loh_vals = sub[sub["core_loh_like"] == True]["UnassignedAffinity"].dropna()
        non_loh_vals = sub[sub["core_loh_like"] == False]["UnassignedAffinity"].dropna()
        bins = np.linspace(-1.0, 1.0, 41)
        ax.hist(non_loh_vals, bins=bins, alpha=0.6, label=f"non-LOH-like (n={len(non_loh_vals):,})",
                color="#91bfdb", density=True)
        ax.hist(loh_vals, bins=bins, alpha=0.6, label=f"LOH-like (n={len(loh_vals):,})",
                color="#fc8d59", density=True)
        ax.axvline(0, color="black", linestyle="--", linewidth=1)
        ax.axvline(0.2, color="gray", linestyle=":", linewidth=0.8)
        ax.axvline(-0.2, color="gray", linestyle=":", linewidth=0.8)
        ax.set_xlabel("Unassigned Affinity Score")
        ax.set_ylabel("Density")
        ax.set_title(f"{mode_name.upper()} / {truth_name}")
        ax.legend(fontsize=7)

fig.suptitle("UnassignedAffinity (Stage 4) Distribution\nHP0/HP3 reads' affinity to HP1 vs HP2 family", fontsize=12)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig05_affinity_score_distribution.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 05 saved")

# Fig 6: Same-locus HP0 delta (TO - paired)
if len(slc_both) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 左：整體 distribution
    ax = axes[0]
    concordance_groups = slc_both["paired_to_concordance"].unique()
    colors_conc = {"both_tp": "#4dac26", "both_fp": "#d01c8b", "to_only_fp": "#f1b6da",
                   "paired_only_tp": "#b8e186", "to_only_tp": "#92c5de", "paired_only_fp": "#f4a582"}
    for conc in concordance_groups:
        sub_conc = slc_both[slc_both["paired_to_concordance"] == conc]
        if len(sub_conc) > 100:
            ax.hist(sub_conc["hp0_delta"].clip(-0.5, 0.5), bins=50, alpha=0.5,
                    label=f"{conc} (n={len(sub_conc):,})", density=True,
                    color=colors_conc.get(conc, "gray"))
    ax.axvline(0, color="black", linestyle="--", linewidth=1.5, label="delta=0")
    ax.set_xlabel("HP0 Ratio Delta (TO - paired)")
    ax.set_ylabel("Density")
    ax.set_title("Same-locus HP0 Delta (TO - paired)\nby Concordance Group")
    ax.legend(fontsize=7)

    # 右：median delta by sample
    ax = axes[1]
    slc_sample_stat = slc_both.groupby("sample")["hp0_delta"].median().sort_values()
    colors_bar = ["#d73027" if v > 0.02 else "#4575b4" if v < -0.02 else "gray"
                  for v in slc_sample_stat.values]
    ax.barh(slc_sample_stat.index, slc_sample_stat.values, color=colors_bar, edgecolor="k")
    ax.axvline(0, color="black", linewidth=1.5, linestyle="--")
    ax.set_xlabel("Median HP0 Delta (TO - paired)")
    ax.set_title("Median Same-locus HP0 Delta by Sample\n(red=TO systematically higher)")

    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig06_hp0_same_locus_delta.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Fig 06 saved")

# Fig 7: HP0 vs effective_hp_reads scatter（H1 假說）
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, mode_name in zip(axes, ["paired", "to"]):
    sub = df_tf[df_tf["mode"] == mode_name].copy()
    # 取樣，避免 100 萬點
    if len(sub) > 20000:
        sub = sub.sample(20000, random_state=42)
    loh_mask = sub["core_loh_like"] == True
    ax.scatter(sub[~loh_mask]["NumReads"], sub[~loh_mask]["hp0_ratio"],
               alpha=0.05, s=3, c="#91bfdb", label="non-LOH-like")
    ax.scatter(sub[loh_mask]["NumReads"], sub[loh_mask]["hp0_ratio"],
               alpha=0.1, s=3, c="#fc8d59", label="LOH-like")
    ax.set_xlabel("NumReads (coverage)")
    ax.set_ylabel("HP0 Ratio")
    ax.set_title(f"{mode_name.upper()}: HP0 Ratio vs Coverage\n(H1: 低 coverage 導致高 HP0?)")
    ax.legend(fontsize=8, markerscale=5)
    ax.set_xlim(0, 200)
    ax.set_ylim(0, 1)

plt.tight_layout()
plt.savefig(FIG_DIR / "fig07_hp0_vs_coverage.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 07 saved")

# ============================================================
# 補充分析：Tier A 的穩健性確認
# ============================================================
print("\n[7/8] Tier A 穩健性確認 ...")

# Tier A 中，高 enrichment 樣本的詳細分析
tier_a_paired = tier_by_sample[(tier_by_sample["mode"] == "paired") &
                                (tier_by_sample["hp_support_tier"] == "A")].copy()
tier_a_paired = tier_a_paired.sort_values("fp_tp_enrichment", ascending=False)
print("\n  Paired Tier A FP/TP enrichment by sample:")
print(tier_a_paired[["sample", "tp_total", "fp_total", "tp_loh_like_frac",
                      "fp_loh_like_frac", "fp_tp_enrichment", "fisher_pval"]].to_string(index=False))

tier_a_to = tier_by_sample[(tier_by_sample["mode"] == "to") &
                             (tier_by_sample["hp_support_tier"] == "A")].copy()
tier_a_to = tier_a_to.sort_values("fp_tp_enrichment", ascending=False)
print("\n  TO Tier A FP/TP enrichment by sample:")
print(tier_a_to[["sample", "tp_total", "fp_total", "tp_loh_like_frac",
                  "fp_loh_like_frac", "fp_tp_enrichment", "fisher_pval"]].to_string(index=False))

# 跨 Tier 比較（paired）
print("\n  Paired enrichment by tier (aggregated):")
tier_global_paired = tier_global[tier_global["mode"] == "paired"].copy()
tier_global_paired["hp_support_tier"] = pd.Categorical(
    tier_global_paired["hp_support_tier"], categories=TIER_ORDER, ordered=True
)
print(tier_global_paired.sort_values("hp_support_tier")[
    ["hp_support_tier", "tp_total", "fp_total", "tp_loh_like_frac",
     "fp_loh_like_frac", "fp_tp_enrichment", "fisher_pval"]
].to_string(index=False))

# ============================================================
# 生成 Round 2 Summary Report
# ============================================================
print("\n[8/8] 生成 round2_summary.md ...")

# 收集關鍵數字
global_paired = {row["hp_support_tier"]: row
                 for _, row in tier_global[tier_global["mode"] == "paired"].iterrows()}
global_to = {row["hp_support_tier"]: row
             for _, row in tier_global[tier_global["mode"] == "to"].iterrows()}

def fmt_enrichment(d, tier):
    row = d.get(tier)
    if row is None:
        return "N/A"
    return f"{row['fp_tp_enrichment']:.3f}× (TP={int(row['tp_total']):,}, FP={int(row['fp_total']):,}, p={row['fisher_pval']:.2e})"

hp0_mw_df = pd.read_csv(OUT_DIR / "core2_hp0_loh_mannwhitney.tsv", sep="\t")
def fmt_hp0_mw(mode, truth):
    row = hp0_mw_df[(hp0_mw_df["mode"] == mode) & (hp0_mw_df["truth"] == truth)]
    if len(row) == 0:
        return "N/A"
    r = row.iloc[0]
    return f"LOH-like med={r['loh_like_median_hp0']:.4f}, non-LOH med={r['non_loh_like_median_hp0']:.4f}, diff={r['median_diff']:+.4f}, p={r['mw_pval']:.2e}"

report = f"""<!--
建立時間: {GENERATED_AT}
目標: LOH Round 2 執行報告 — effective_hp support 分層 + HP0 來源分析
處理範圍: 748,391 region rows, 7 samples, paired + TO
關聯檔案:
  - docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md
  - docs/architecture/20260327_InterSubMod研究願景定錨_01.md
  - docs/plans/2026/03/20260327_LOH_round2_execution_spec_01.md
  - output/synthesis/observation_workspaces/20260327_loh_round2_support_hp0_analysis/
-->

# LOH Round 2 執行報告：effective_hp Support 分層 + HP0 來源分析

> 生成時間：{GENERATED_AT}
> 分析腳本：scripts/analysis/build_loh_round2_support_hp0_analysis.py
> 輸出 workspace：output/synthesis/observation_workspaces/20260327_loh_round2_support_hp0_analysis/

---

## 一、分析目標與背景

Round 1 確立了「LOH-like 分佈診斷底圖」，但遺留了兩個核心問題：

1. **effective_hp support 分層**：在 `effective_hp_reads < 10` 的雜訊過濾後，LOH-like 的 FP/TP enrichment 是否在 Tier A（≥30 reads）顯著且可泛化？
2. **HP0 來源**：TO 中大量 HP0 的成因是什麼？HP0 是否在 LOH-like region 特別高？Stage 4 affinity score 是否有補充 evidence 的潛力？

---

## 二、Tier 定義（本輪固定）

| Tier | 條件 | 說明 |
|------|------|------|
| C0 | `effective_hp_reads = 0` | 完全沒有 HP reads，HP_Ratio=0.5 是假象 |
| C  | `1–9 reads` | 極弱 support，幾乎不可信 |
| B  | `10–29 reads` | 中等 support，謹慎解讀 |
| A  | `≥ 30 reads` | 強 support，可信 LOH evidence 候選 |

---

## 三、Core 1 結果：effective_hp Support 分層

### 3.1 全域 paired 端 Tier 分析

| Tier | TP Total | FP Total | TP LOH-like% | FP LOH-like% | FP/TP Enrichment | Fisher p |
|------|---------|---------|------------|------------|----------------|---------|
| C0   | {int(global_paired.get('C0', {}).get('tp_total', 0)):,} | {int(global_paired.get('C0', {}).get('fp_total', 0)):,} | {global_paired.get('C0', {}).get('tp_loh_like_frac', float('nan')):.3f} | {global_paired.get('C0', {}).get('fp_loh_like_frac', float('nan')):.3f} | {fmt_enrichment(global_paired, 'C0').split('×')[0]}× | {global_paired.get('C0', {}).get('fisher_pval', float('nan')):.2e} |
| C    | {int(global_paired.get('C', {}).get('tp_total', 0)):,} | {int(global_paired.get('C', {}).get('fp_total', 0)):,} | {global_paired.get('C', {}).get('tp_loh_like_frac', float('nan')):.3f} | {global_paired.get('C', {}).get('fp_loh_like_frac', float('nan')):.3f} | {fmt_enrichment(global_paired, 'C').split('×')[0]}× | {global_paired.get('C', {}).get('fisher_pval', float('nan')):.2e} |
| B    | {int(global_paired.get('B', {}).get('tp_total', 0)):,} | {int(global_paired.get('B', {}).get('fp_total', 0)):,} | {global_paired.get('B', {}).get('tp_loh_like_frac', float('nan')):.3f} | {global_paired.get('B', {}).get('fp_loh_like_frac', float('nan')):.3f} | {fmt_enrichment(global_paired, 'B').split('×')[0]}× | {global_paired.get('B', {}).get('fisher_pval', float('nan')):.2e} |
| A    | {int(global_paired.get('A', {}).get('tp_total', 0)):,} | {int(global_paired.get('A', {}).get('fp_total', 0)):,} | {global_paired.get('A', {}).get('tp_loh_like_frac', float('nan')):.3f} | {global_paired.get('A', {}).get('fp_loh_like_frac', float('nan')):.3f} | {fmt_enrichment(global_paired, 'A').split('×')[0]}× | {global_paired.get('A', {}).get('fisher_pval', float('nan')):.2e} |

> Round 1 整體 paired enrichment = 1.194×，比較以上各 Tier 數值。

### 3.2 全域 TO 端 Tier 分析

| Tier | TP Total | FP Total | TP LOH-like% | FP LOH-like% | FP/TP Enrichment | Fisher p |
|------|---------|---------|------------|------------|----------------|---------|
| C0   | {int(global_to.get('C0', {}).get('tp_total', 0)):,} | {int(global_to.get('C0', {}).get('fp_total', 0)):,} | {global_to.get('C0', {}).get('tp_loh_like_frac', float('nan')):.3f} | {global_to.get('C0', {}).get('fp_loh_like_frac', float('nan')):.3f} | {fmt_enrichment(global_to, 'C0').split('×')[0]}× | {global_to.get('C0', {}).get('fisher_pval', float('nan')):.2e} |
| C    | {int(global_to.get('C', {}).get('tp_total', 0)):,} | {int(global_to.get('C', {}).get('fp_total', 0)):,} | {global_to.get('C', {}).get('tp_loh_like_frac', float('nan')):.3f} | {global_to.get('C', {}).get('fp_loh_like_frac', float('nan')):.3f} | {fmt_enrichment(global_to, 'C').split('×')[0]}× | {global_to.get('C', {}).get('fisher_pval', float('nan')):.2e} |
| B    | {int(global_to.get('B', {}).get('tp_total', 0)):,} | {int(global_to.get('B', {}).get('fp_total', 0)):,} | {global_to.get('B', {}).get('tp_loh_like_frac', float('nan')):.3f} | {global_to.get('B', {}).get('fp_loh_like_frac', float('nan')):.3f} | {fmt_enrichment(global_to, 'B').split('×')[0]}× | {global_to.get('B', {}).get('fisher_pval', float('nan')):.2e} |
| A    | {int(global_to.get('A', {}).get('tp_total', 0)):,} | {int(global_to.get('A', {}).get('fp_total', 0)):,} | {global_to.get('A', {}).get('tp_loh_like_frac', float('nan')):.3f} | {global_to.get('A', {}).get('fp_loh_like_frac', float('nan')):.3f} | {fmt_enrichment(global_to, 'A').split('×')[0]}× | {global_to.get('A', {}).get('fisher_pval', float('nan')):.2e} |

> Round 1 整體 TO enrichment = 0.912×，比較以上各 Tier 數值。

### 3.3 Paired 各樣本 Tier A enrichment

（見 figures/fig02_paired_tier_enrichment_heatmap.png 及 core1_tier_enrichment_by_sample_mode.tsv）

### 3.4 Core 1 直接觀察

（待腳本執行後根據實際數值填入）

---

## 四、Core 2 結果：HP0 來源分析

### 4.1 HP0 ratio：LOH-like vs non-LOH-like

| Mode | Truth | LOH-like med HP0 | non-LOH-like med HP0 | Diff | MW p-value |
|------|-------|-----------------|-------------------|------|-----------|
| paired | TP | - | - | - | - |
| paired | FP | - | - | - | - |
| to | TP | - | - | - | - |
| to | FP | - | - | - | - |

> （詳細數值見 core2_hp0_loh_mannwhitney.tsv）

### 4.2 UnassignedAffinity（Stage 4）分析

Stage 4 計算 HP0/HP3 reads 相對 HP1 vs HP2 family 的傾向性分數。

（見 core2_affinity_by_loh_status.tsv 及 figures/fig05_affinity_score_distribution.png）

### 4.3 HP0 假說驗證結果

**H1（HP0 ↔ Coverage）**：
- Spearman 相關係數（見 core2_hp0_coverage_correlation.tsv）

**H3（TO HP0 > paired HP0 同位點）**：
- Wilcoxon 檢驗結果（見 core2_hp0_same_locus_delta.tsv）

---

## 五、圖表索引

| 圖 | 說明 |
|----|------|
| fig01_tier_enrichment_global.png | 全域 Tier A/B/C/C0 的 FP/TP enrichment（paired vs TO） |
| fig02_paired_tier_enrichment_heatmap.png | 各樣本 × Tier 的 paired FP/TP enrichment heatmap |
| fig03_tier_loh_frac_barplot.png | 各 Tier 的 LOH-like fraction（TP vs FP 對比） |
| fig04_hp0_distribution_violin.png | HP0 ratio violin（LOH-like vs non-LOH-like） |
| fig05_affinity_score_distribution.png | Stage 4 affinity score 分佈 |
| fig06_hp0_same_locus_delta.png | 同位點 TO - paired HP0 delta |
| fig07_hp0_vs_coverage.png | HP0 ratio vs coverage scatter（H1 假說） |

---

## 六、數據輸出索引

| 檔案 | 說明 |
|------|------|
| core1_tier_enrichment_global.tsv | 全域 mode × tier enrichment |
| core1_tier_enrichment_by_sample_mode.tsv | 各樣本 × mode × tier enrichment |
| core2_hp0_by_loh_status.tsv | HP0 統計 by LOH status |
| core2_hp0_by_sample_loh_status.tsv | HP0 統計 by sample × LOH status |
| core2_hp0_loh_mannwhitney.tsv | HP0 Mann-Whitney 檢驗結果 |
| core2_affinity_by_loh_status.tsv | Stage 4 affinity by LOH status |
| core2_affinity_loh_mannwhitney.tsv | Affinity Mann-Whitney 結果 |
| core2_hp0_coverage_correlation.tsv | HP0 vs coverage Spearman 相關 |
| core2_hp0_same_locus_delta.tsv | 同位點 HP0 delta 統計 |
| core2_hp0_by_concordance.tsv | HP0 by concordance group |

---

## 七、Round 2 推論（待數值填入後完成）

（詳見執行後的 round2_final_conclusions.md）

---

## 八、Round 2 的限制與不應過度宣稱的事

1. Tier 定義（A=30, B=10）目前尚未做敏感性分析（例如 A=20 或 A=50 閾值）
2. HP0 來源分析是觀察性的，不能直接宣稱因果
3. Stage 4 affinity score 的統計顯著性仍需在更多樣本中驗證
4. LOH × methylation 聯合分析暫緩至 Round 3

"""

(OUT_DIR / "round2_report_template.md").write_text(report, encoding="utf-8")

# 寫入 round context
context = {
    "round": 2,
    "generated_at": GENERATED_AT,
    "core1_question": "在哪個 effective_hp support tier，LOH evidence 才是可信的 evidence panel 成員？",
    "core2_question": "TO 中大量 HP0 的主要來源是什麼？HP0 在 LOH-like region 是否特別高？",
    "tier_definitions": {
        "A": "effective_hp_reads >= 30 (strong evidence)",
        "B": "effective_hp_reads 10-29 (moderate)",
        "C": "effective_hp_reads 1-9 (weak, likely noise)",
        "C0": "effective_hp_reads == 0 (no HP reads, HP_Ratio=0.5 is artifact)"
    },
    "data_source": str(ROUND1_WS / "all_region_rows.tsv.gz"),
    "total_rows": len(df),
    "tp_fp_rows": len(df_tf),
}
with open(OUT_DIR / "round2_context.json", "w") as f:
    json.dump(context, f, indent=2, ensure_ascii=False)

print(f"\n完成！所有輸出已寫入 {OUT_DIR}")
print("主要輸出：")
for f in sorted(OUT_DIR.iterdir()):
    if f.is_file():
        size = f.stat().st_size
        print(f"  {f.name} ({size:,} bytes)")
print(f"\n圖表：")
for f in sorted(FIG_DIR.iterdir()):
    print(f"  figures/{f.name}")
