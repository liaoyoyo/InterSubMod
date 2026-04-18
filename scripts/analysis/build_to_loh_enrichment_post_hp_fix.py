#!/usr/bin/env python3
"""
TO LOH Enrichment Analysis — Post HP Integer Tag Fix (2026-03-30)

Uses corrected TO significance_summary.csv files (updated 2026-03-30 after
HP:i:11/21/33 mapping bug fix) to compute LOH enrichment across all TO samples.

Key questions:
  1. After fix, how many TO regions are LOH-eligible (eff_hp >= 30)?
  2. What is the TO LOH enrichment ratio (FP% / TP%) by Tier?
  3. How does COLO829/H1437/H2009 specifically look (previously most impacted)?
  4. Compare corrected TO vs paired LOH enrichment patterns.

Output:
  /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
  20260330_to_loh_enrichment_post_hp_fix/
    loh_eligibility_before_after.tsv
    loh_enrichment_by_tier.tsv
    loh_enrichment_by_sample_tier.tsv
    cross_mode_comparison.tsv
    summary.json
    figures/
"""

from __future__ import annotations

import csv
import json
import warnings
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")

# ============================================================
# 路徑設定
# ============================================================
RESEARCH_ROUNDS = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds")
ROUND1_WS = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit")
OUT_DIR = RESEARCH_ROUNDS / "20260330_to_loh_enrichment_post_hp_fix"
FIG_DIR = OUT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

GENERATED_AT = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# TO 樣本設定（run_id -> significance_summary.csv 位置）
TO_SAMPLES = [
    dict(sample="HCC1395",       label="HCC1395 5kHz",    run_id="20260315_hcc1395_to_pilot"),
    dict(sample="HCC1395_DORADO",label="HCC1395 DORADO",  run_id="20260315_hcc1395_dorado_to_pilot"),
    dict(sample="COLO829",       label="COLO829",          run_id="20260317_colo829_to_pilot_1"),
    dict(sample="H1437",         label="H1437",            run_id="20260318_h1437_to_pilot_fastresume"),
    dict(sample="H2009",         label="H2009",            run_id="20260318_h2009_to_pilot_fastresume"),
    dict(sample="HCC1937",       label="HCC1937",          run_id="20260318_hcc1937_to_pilot_fastresume"),
    dict(sample="HCC1954",       label="HCC1954",          run_id="20260318_hcc1954_to_pilot_fastresume"),
]

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]

# Round 1 舊版 LOH eligible 數字（HP fix 前，從 memory 記錄）
ROUND1_OLD_LOH_ELIG = {
    # (sample, subset): (count, pct)
    ("HCC1395",  "TP"): (16665, 58.5), ("HCC1395",  "FP"): (6926, 59.7),
    ("COLO829",  "TP"): (232,   0.7),  ("COLO829",  "FP"): (227,  1.3),
    ("H1437",    "TP"): (11144, 24.5), ("H1437",    "FP"): (3319, 24.7),
    ("H2009",    "TP"): (81023, 64.5), ("H2009",    "FP"): (8579, 71.6),
    ("HCC1937",  "TP"): (10948, 86.7), ("HCC1937",  "FP"): (10027, 83.3),
    ("HCC1954",  "TP"): (11944, 70.0), ("HCC1954",  "FP"): (35074, 69.8),
    ("HCC1395_DORADO", "TP"): (13521, 46.9), ("HCC1395_DORADO", "FP"): (3797, 32.8),
}

# Paired LOH enrichment (from Round 1, paired samples, for cross-mode comparison)
ROUND1_PAIRED_ENRICHMENT = {}

# ============================================================
# 輔助函數
# ============================================================

def get_tp_fp_paths(run_id: str) -> Tuple[Path, Path]:
    base = RESEARCH_ROUNDS / run_id / "step05_intersubmod"
    tp = base / "intersubmod_tp" / "significance_summary.csv"
    fp = base / "intersubmod_fp" / "significance_summary.csv"
    return tp, fp


def load_summary(path: Path) -> Optional[pd.DataFrame]:
    if not path.exists():
        print(f"  [MISSING] {path}")
        return None
    df = pd.read_csv(path, low_memory=False)
    df["effective_hp"] = pd.to_numeric(df["HP1FamilyN"], errors="coerce").fillna(0) + \
                         pd.to_numeric(df["HP2FamilyN"], errors="coerce").fillna(0)
    df["effective_hp"] = df["effective_hp"].astype(int)
    df["loh_like"] = df["Potential_LOH"].astype(str).str.upper().isin(["TRUE", "1", "YES"])
    df["NHP0"] = pd.to_numeric(df["NHP0"], errors="coerce").fillna(0)
    df["NHP3"] = pd.to_numeric(df["NHP3"], errors="coerce").fillna(0)
    df["NumReads"] = pd.to_numeric(df["NumReads"], errors="coerce").fillna(0)
    df["hp0_ratio"] = df["NHP0"] / df["NumReads"].clip(lower=1)
    return df


def assign_tier(n: int) -> str:
    """Tier 定義：C0(0), C(<10), B(10-29), A(30-49), A+(≥50)"""
    if n == 0:    return "C0"
    elif n < 10:  return "C"
    elif n < 30:  return "B"
    elif n < 50:  return "A"
    else:         return "A+"


def compute_enrichment(tp_df: pd.DataFrame, fp_df: pd.DataFrame,
                       mask_tp=None, mask_fp=None,
                       label: str = "all") -> dict:
    """計算 LOH enrichment 基本統計"""
    if mask_tp is None:
        mask_tp = pd.Series(True, index=tp_df.index)
    if mask_fp is None:
        mask_fp = pd.Series(True, index=fp_df.index)
    sub_tp = tp_df[mask_tp]
    sub_fp = fp_df[mask_fp]
    tp_n = len(sub_tp)
    fp_n = len(sub_fp)
    tp_loh = sub_tp["loh_like"].sum()
    fp_loh = sub_fp["loh_like"].sum()
    tp_frac = tp_loh / tp_n if tp_n > 0 else 0.0
    fp_frac = fp_loh / fp_n if fp_n > 0 else 0.0
    enrichment = fp_frac / tp_frac if tp_frac > 0 else float("nan")

    # Fisher's exact test
    table = [[tp_loh, tp_n - tp_loh], [fp_loh, fp_n - fp_loh]]
    try:
        odds, pval = stats.fisher_exact(table)
    except Exception:
        odds, pval = float("nan"), float("nan")

    return {
        "subset": label,
        "tp_n": tp_n, "fp_n": fp_n,
        "tp_loh_n": int(tp_loh), "fp_loh_n": int(fp_loh),
        "tp_loh_frac": round(tp_frac, 4),
        "fp_loh_frac": round(fp_frac, 4),
        "enrichment": round(enrichment, 4) if not np.isnan(enrichment) else None,
        "fisher_p": round(pval, 6) if not np.isnan(pval) else None,
        "median_eff_hp_tp": int(sub_tp["effective_hp"].median()) if tp_n > 0 else None,
        "median_eff_hp_fp": int(sub_fp["effective_hp"].median()) if fp_n > 0 else None,
    }


# ============================================================
# Step 1: 載入所有 TO 數據
# ============================================================
print("[Step 1] 載入 TO significance_summary.csv（HP fix 後，2026-03-30）")

sample_data: Dict[str, dict] = {}  # sample -> {tp_df, fp_df, label, run_id, status}

for cfg in TO_SAMPLES:
    sample = cfg["sample"]
    run_id = cfg["run_id"]
    tp_path, fp_path = get_tp_fp_paths(run_id)
    tp_df = load_summary(tp_path)
    fp_df = load_summary(fp_path)

    ok = tp_df is not None and fp_df is not None
    modtime = ""
    if ok:
        import os
        modtime = datetime.fromtimestamp(os.path.getmtime(tp_path)).strftime("%Y-%m-%d %H:%M")

    print(f"  {sample}: TP={len(tp_df) if tp_df is not None else 'MISSING'} "
          f"FP={len(fp_df) if fp_df is not None else 'MISSING'} "
          f"(updated {modtime})")
    sample_data[sample] = dict(tp_df=tp_df, fp_df=fp_df,
                               label=cfg["label"], run_id=run_id, ok=ok)

# ============================================================
# Step 2: LOH Eligibility 比較（修正前 vs 修正後）
# ============================================================
print("\n[Step 2] LOH Eligibility 比較（eff_hp >= 30）")

eligibility_rows = []
LOH_ELIG_THRESH = 30

for sample in SAMPLE_ORDER:
    if sample not in sample_data or not sample_data[sample]["ok"]:
        continue
    tp_df = sample_data[sample]["tp_df"]
    fp_df = sample_data[sample]["fp_df"]
    label = sample_data[sample]["label"]

    for subset, df in [("TP", tp_df), ("FP", fp_df)]:
        new_n = len(df)
        new_elig = (df["effective_hp"] >= LOH_ELIG_THRESH).sum()
        new_pct = 100 * new_elig / new_n if new_n > 0 else 0

        old = ROUND1_OLD_LOH_ELIG.get((sample, subset), (None, None))
        old_elig, old_pct = old

        row = {
            "sample": sample, "label": label, "subset": subset,
            "new_total": new_n, "new_loh_elig": int(new_elig), "new_loh_elig_pct": round(new_pct, 1),
            "old_loh_elig": old_elig, "old_loh_elig_pct": old_pct,
            "delta_pct": round(new_pct - old_pct, 1) if old_pct is not None else None,
            "fold_change": round(new_elig / old_elig, 1) if old_elig and old_elig > 0 else None,
        }
        eligibility_rows.append(row)
        print(f"  {sample} {subset}: {old_pct}% → {new_pct:.1f}% "
              f"(+{new_pct - old_pct:.1f}pp, {new_elig:,} regions eligible)")

elig_df = pd.DataFrame(eligibility_rows)
elig_df.to_csv(OUT_DIR / "loh_eligibility_before_after.tsv", sep="\t", index=False)

# ============================================================
# Step 3: TO LOH Enrichment by Tier（全樣本 aggregated）
# ============================================================
print("\n[Step 3] TO LOH Enrichment by Tier（全樣本聚合）")

# Tier 定義清單
TIER_DEFS = [
    ("All",      lambda df: pd.Series(True, index=df.index)),
    ("C0",       lambda df: df["effective_hp"] == 0),
    ("C",        lambda df: (df["effective_hp"] > 0) & (df["effective_hp"] < 10)),
    ("B",        lambda df: (df["effective_hp"] >= 10) & (df["effective_hp"] < 30)),
    ("A",        lambda df: (df["effective_hp"] >= 30) & (df["effective_hp"] < 50)),
    ("A+",       lambda df: df["effective_hp"] >= 50),
    ("A≥30",     lambda df: df["effective_hp"] >= 30),
    ("A≥10",     lambda df: df["effective_hp"] >= 10),
]

# 聚合所有 TO 樣本
all_tp_list = []
all_fp_list = []
for sample in SAMPLE_ORDER:
    if sample not in sample_data or not sample_data[sample]["ok"]:
        continue
    tp_df = sample_data[sample]["tp_df"].copy()
    fp_df = sample_data[sample]["fp_df"].copy()
    tp_df["sample"] = sample
    fp_df["sample"] = sample
    all_tp_list.append(tp_df)
    all_fp_list.append(fp_df)

all_tp = pd.concat(all_tp_list, ignore_index=True)
all_fp = pd.concat(all_fp_list, ignore_index=True)

tier_rows = []
for tier_name, tier_fn in TIER_DEFS:
    row = compute_enrichment(all_tp, all_fp,
                             mask_tp=tier_fn(all_tp), mask_fp=tier_fn(all_fp),
                             label=tier_name)
    row["mode"] = "TO"
    tier_rows.append(row)
    enrich_str = f"{row['enrichment']:.3f}×" if row["enrichment"] else "N/A"
    p_str = f"p={row['fisher_p']:.2e}" if row["fisher_p"] else ""
    print(f"  Tier {tier_name:6s}: TP={row['tp_n']:6d} ({row['tp_loh_frac']:.1%} LOH)  "
          f"FP={row['fp_n']:6d} ({row['fp_loh_frac']:.1%} LOH)  "
          f"enrich={enrich_str}  {p_str}")

tier_df = pd.DataFrame(tier_rows)
tier_df.to_csv(OUT_DIR / "loh_enrichment_by_tier.tsv", sep="\t", index=False)

# ============================================================
# Step 4: Per-Sample LOH Enrichment by Tier
# ============================================================
print("\n[Step 4] Per-Sample LOH Enrichment by Tier")

per_sample_rows = []
for sample in SAMPLE_ORDER:
    if sample not in sample_data or not sample_data[sample]["ok"]:
        continue
    tp_df = sample_data[sample]["tp_df"]
    fp_df = sample_data[sample]["fp_df"]
    label = sample_data[sample]["label"]

    for tier_name, tier_fn in TIER_DEFS:
        row = compute_enrichment(tp_df, fp_df,
                                 mask_tp=tier_fn(tp_df), mask_fp=tier_fn(fp_df),
                                 label=tier_name)
        row["sample"] = sample
        row["sample_label"] = label
        row["mode"] = "TO"
        per_sample_rows.append(row)

per_sample_df = pd.DataFrame(per_sample_rows)
per_sample_df.to_csv(OUT_DIR / "loh_enrichment_by_sample_tier.tsv", sep="\t", index=False)

# 印出重點樣本 A≥30 tier 結果
print("\n  [重點樣本 eff_hp ≥ 30 tier]")
focus = ["COLO829", "H1437", "H2009"]
for sample in focus:
    row = per_sample_df[(per_sample_df["sample"]==sample) & (per_sample_df["subset"]=="A≥30")]
    if not row.empty:
        r = row.iloc[0]
        print(f"  {sample:20s}: TP={r['tp_n']:5d} ({r['tp_loh_frac']:.1%})  "
              f"FP={r['fp_n']:5d} ({r['fp_loh_frac']:.1%})  "
              f"enrich={r['enrichment']:.3f}×  p={r['fisher_p']:.2e}")
    else:
        print(f"  {sample}: no A≥30 data")

# ============================================================
# Step 5: Cross-mode comparison（TO vs Paired，Round 1 paired enrichment）
# ============================================================
print("\n[Step 5] Cross-mode 比較（TO new vs TO old vs paired）")

# 載入 Round 1 paired enrichment（from 已有 loh_enrichment_by_sample_mode.tsv）
paired_enrich_path = ROUND1_WS / "loh_enrichment_by_sample_mode.tsv"
paired_old_to_enrich = {}
if paired_enrich_path.exists():
    round1_enrich = pd.read_csv(paired_enrich_path, sep="\t")
    for _, row in round1_enrich.iterrows():
        key = (row["sample"], row["mode"])
        paired_old_to_enrich[key] = {
            "tp_loh_frac": row["tp_core_loh_like_frac"],
            "fp_loh_frac": row["fp_core_loh_like_frac"],
            "enrichment": row["core_loh_fp_vs_tp_enrichment"],
            "median_tp_eff_hp": row["median_tp_effective_hp_reads"],
        }

cross_rows = []
for sample in SAMPLE_ORDER:
    if sample not in sample_data or not sample_data[sample]["ok"]:
        continue
    tp_df = sample_data[sample]["tp_df"]
    fp_df = sample_data[sample]["fp_df"]

    # TO new: All 和 A≥30 tier
    for tier_name, tier_fn in [("All", lambda df: pd.Series(True, index=df.index)),
                                ("A≥30", lambda df: df["effective_hp"] >= 30)]:
        row_new = compute_enrichment(tp_df, fp_df,
                                     mask_tp=tier_fn(tp_df), mask_fp=tier_fn(fp_df),
                                     label=tier_name)
        row_new["sample"] = sample
        row_new["data_version"] = "TO_new_post_fix"
        row_new["tier"] = tier_name
        cross_rows.append(row_new)

    # TO old enrichment（Round 1 全域，無 tier filter）
    old = paired_old_to_enrich.get((sample, "to"))
    if old:
        cross_rows.append({
            "sample": sample, "subset": "All", "data_version": "TO_old_pre_fix",
            "tier": "All",
            "tp_n": None, "fp_n": None,
            "tp_loh_frac": round(old["tp_loh_frac"], 4),
            "fp_loh_frac": round(old["fp_loh_frac"], 4),
            "enrichment": round(old["enrichment"], 4),
            "fisher_p": None,
            "median_eff_hp_tp": old["median_tp_eff_hp"],
            "median_eff_hp_fp": None,
        })

    # Paired enrichment（from Round 1）
    old_paired = paired_old_to_enrich.get((sample, "paired"))
    if old_paired:
        cross_rows.append({
            "sample": sample, "subset": "All", "data_version": "paired_round1",
            "tier": "All",
            "tp_n": None, "fp_n": None,
            "tp_loh_frac": round(old_paired["tp_loh_frac"], 4),
            "fp_loh_frac": round(old_paired["fp_loh_frac"], 4),
            "enrichment": round(old_paired["enrichment"], 4),
            "fisher_p": None,
            "median_eff_hp_tp": old_paired["median_tp_eff_hp"],
            "median_eff_hp_fp": None,
        })

cross_df = pd.DataFrame(cross_rows)
cross_df.to_csv(OUT_DIR / "cross_mode_comparison.tsv", sep="\t", index=False)

# ============================================================
# Step 6: 可視化
# ============================================================
print("\n[Step 6] 生成圖表")

# Figure 1: LOH Eligibility Before vs After（TP only）
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("TO LOH Eligibility Change After HP Integer Tag Fix\n(eff_hp ≥ 30 threshold)",
             fontsize=14, fontweight="bold")

elig_tp = elig_df[elig_df["subset"] == "TP"].copy()
x = np.arange(len(elig_tp))
width = 0.35
ax = axes[0]
bars_old = ax.bar(x - width/2, elig_tp["old_loh_elig_pct"].fillna(0), width,
                  label="Before Fix (pre 2026-03-30)", color="#ff7f7f", edgecolor="black", linewidth=0.5)
bars_new = ax.bar(x + width/2, elig_tp["new_loh_elig_pct"], width,
                  label="After Fix (2026-03-30)", color="#7fbf7f", edgecolor="black", linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels(elig_tp["sample"].tolist(), rotation=30, ha="right")
ax.set_ylabel("LOH Eligible TP Regions (%)")
ax.set_title("TP Regions — LOH Eligible %")
ax.legend()
ax.set_ylim(0, 110)
ax.axhline(90, color="gray", linestyle="--", alpha=0.5, label="90%")
for bar in bars_new:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, h + 1, f"{h:.0f}%",
            ha="center", va="bottom", fontsize=7)

elig_fp = elig_df[elig_df["subset"] == "FP"].copy()
ax2 = axes[1]
bars_old2 = ax2.bar(x - width/2, elig_fp["old_loh_elig_pct"].fillna(0), width,
                    label="Before Fix", color="#ff7f7f", edgecolor="black", linewidth=0.5)
bars_new2 = ax2.bar(x + width/2, elig_fp["new_loh_elig_pct"], width,
                    label="After Fix", color="#7fbf7f", edgecolor="black", linewidth=0.5)
ax2.set_xticks(x)
ax2.set_xticklabels(elig_fp["sample"].tolist(), rotation=30, ha="right")
ax2.set_ylabel("LOH Eligible FP Regions (%)")
ax2.set_title("FP Regions — LOH Eligible %")
ax2.legend()
ax2.set_ylim(0, 110)
for bar in bars_new2:
    h = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2, h + 1, f"{h:.0f}%",
             ha="center", va="bottom", fontsize=7)

plt.tight_layout()
plt.savefig(FIG_DIR / "fig1_loh_eligibility_before_after.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved fig1_loh_eligibility_before_after.png")

# Figure 2: LOH Enrichment by Tier（全樣本聚合）
tier_order = ["C0", "C", "B", "A", "A+", "A≥30"]
tier_plot = tier_df[tier_df["subset"].isin(tier_order)].set_index("subset")
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(tier_order))
tp_fracs = [tier_plot.loc[t, "tp_loh_frac"] * 100 if t in tier_plot.index else 0 for t in tier_order]
fp_fracs = [tier_plot.loc[t, "fp_loh_frac"] * 100 if t in tier_plot.index else 0 for t in tier_order]
enrichments = [tier_plot.loc[t, "enrichment"] if t in tier_plot.index else 1.0 for t in tier_order]

width = 0.35
bars_tp = ax.bar(x - width/2, tp_fracs, width, label="TP LOH-like%", color="#4488cc", alpha=0.8, edgecolor="black", linewidth=0.5)
bars_fp = ax.bar(x + width/2, fp_fracs, width, label="FP LOH-like%", color="#cc4444", alpha=0.8, edgecolor="black", linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels(tier_order)
ax.set_ylabel("LOH-like Region %")
ax.set_xlabel("eff_hp Tier")
ax.set_title("TO LOH Enrichment by eff_hp Tier\n(All samples aggregated, post HP-fix)")
ax.legend()

# Enrichment text annotation
ax2b = ax.twinx()
ax2b.set_ylabel("FP/TP Enrichment Ratio", color="darkgreen")
ax2b.plot(x, enrichments, "D-", color="darkgreen", linewidth=2, markersize=8, label="FP/TP enrichment")
ax2b.axhline(1.0, color="gray", linestyle="--", alpha=0.5)
for xi, e in zip(x, enrichments):
    ax2b.annotate(f"{e:.2f}×", (xi, e), textcoords="offset points",
                  xytext=(0, 8), ha="center", color="darkgreen", fontsize=9, fontweight="bold")
ax2b.tick_params(axis="y", colors="darkgreen")

plt.tight_layout()
plt.savefig(FIG_DIR / "fig2_loh_enrichment_by_tier.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved fig2_loh_enrichment_by_tier.png")

# Figure 3: Per-sample enrichment heatmap（A≥30 tier）
fig, ax = plt.subplots(figsize=(12, 5))
focus_tiers = ["All", "B", "A", "A+", "A≥30"]
heatmap_data = []
heatmap_pvals = []
for sample in SAMPLE_ORDER:
    row_e = []
    row_p = []
    for t in focus_tiers:
        sub = per_sample_df[(per_sample_df["sample"] == sample) & (per_sample_df["subset"] == t)]
        if not sub.empty:
            row_e.append(sub.iloc[0]["enrichment"] or 1.0)
            row_p.append(sub.iloc[0]["fisher_p"] or 1.0)
        else:
            row_e.append(1.0)
            row_p.append(1.0)
    heatmap_data.append(row_e)
    heatmap_pvals.append(row_p)

hm = pd.DataFrame(heatmap_data, index=SAMPLE_ORDER, columns=focus_tiers)
sns.heatmap(hm, annot=True, fmt=".2f", cmap="RdYlGn_r", center=1.0,
            vmin=0.5, vmax=2.5, ax=ax, linewidths=0.5,
            annot_kws={"size": 10})
ax.set_title("TO LOH FP/TP Enrichment Ratio by Sample × eff_hp Tier\n(>1 = FP enriched, <1 = TP enriched)")
ax.set_xlabel("eff_hp Tier")
ax.set_ylabel("Sample")
plt.tight_layout()
plt.savefig(FIG_DIR / "fig3_per_sample_enrichment_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved fig3_per_sample_enrichment_heatmap.png")

# Figure 4: eff_hp distribution shift（COLO829, H1437, H2009）
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle("TO eff_hp Distribution After HP Fix\n(COLO829 / H1437 / H2009)", fontsize=13)
focus_samples_fig = ["COLO829", "H1437", "H2009"]
for i, sample in enumerate(focus_samples_fig):
    ax = axes[i]
    if sample not in sample_data or not sample_data[sample]["ok"]:
        ax.set_title(f"{sample} (no data)")
        continue
    tp_df = sample_data[sample]["tp_df"]
    fp_df = sample_data[sample]["fp_df"]
    bins = [0, 1, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300]
    ax.hist(tp_df["effective_hp"].clip(upper=300), bins=bins, alpha=0.5,
            color="#4488cc", label="TP", density=True)
    ax.hist(fp_df["effective_hp"].clip(upper=300), bins=bins, alpha=0.5,
            color="#cc4444", label="FP", density=True)
    ax.axvline(30, color="orange", linestyle="--", label="eff_hp=30")
    tp_elig_pct = 100 * (tp_df["effective_hp"] >= 30).mean()
    fp_elig_pct = 100 * (fp_df["effective_hp"] >= 30).mean()
    ax.set_title(f"{sample}\nTP ≥30: {tp_elig_pct:.1f}%  FP ≥30: {fp_elig_pct:.1f}%", fontsize=10)
    ax.set_xlabel("effective_hp reads")
    ax.set_ylabel("Density")
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(FIG_DIR / "fig4_eff_hp_distribution_focus.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved fig4_eff_hp_distribution_focus.png")

# ============================================================
# Step 7: 整合摘要
# ============================================================
print("\n[Step 7] 整合摘要統計")

# 重點樣本 All tier 整合
focus_all = per_sample_df[per_sample_df["subset"] == "All"][
    ["sample", "tp_n", "fp_n", "tp_loh_frac", "fp_loh_frac",
     "enrichment", "fisher_p", "median_eff_hp_tp", "median_eff_hp_fp"]]
print("\n  === TO LOH Enrichment Summary (All regions, post HP-fix) ===")
for _, r in focus_all.iterrows():
    enrich_str = f"{r['enrichment']:.3f}×" if r["enrichment"] else "N/A"
    print(f"  {r['sample']:20s}: TP={r['tp_n']:6d} LOH={r['tp_loh_frac']:.1%}  "
          f"FP={r['fp_n']:6d} LOH={r['fp_loh_frac']:.1%}  "
          f"enrich={enrich_str}")

# 重點樣本 A≥30 tier
focus_a30 = per_sample_df[per_sample_df["subset"] == "A≥30"][
    ["sample", "tp_n", "fp_n", "tp_loh_frac", "fp_loh_frac",
     "enrichment", "fisher_p"]]
print("\n  === TO LOH Enrichment (eff_hp ≥ 30, post HP-fix) ===")
for _, r in focus_a30.iterrows():
    enrich_str = f"{r['enrichment']:.3f}×" if r["enrichment"] else "N/A"
    p_str = f"p={r['fisher_p']:.2e}" if r["fisher_p"] else ""
    print(f"  {r['sample']:20s}: TP={r['tp_n']:6d} LOH={r['tp_loh_frac']:.1%}  "
          f"FP={r['fp_n']:6d} LOH={r['fp_loh_frac']:.1%}  "
          f"enrich={enrich_str} {p_str}")

# 全域聚合 A≥30 vs paired 對比
agg_a30 = tier_df[tier_df["subset"] == "A≥30"]
paired_round1 = cross_df[cross_df["data_version"] == "paired_round1"]
if not agg_a30.empty:
    r = agg_a30.iloc[0]
    print(f"\n  TO All-samples A≥30: TP={r['tp_n']} ({r['tp_loh_frac']:.1%} LOH)  "
          f"FP={r['fp_n']} ({r['fp_loh_frac']:.1%} LOH)  "
          f"enrichment={r['enrichment']:.3f}× p={r['fisher_p']:.2e}")

# summary.json
summary = {
    "generated_at": GENERATED_AT,
    "data_version": "post_hp_fix_2026-03-30",
    "samples_analyzed": SAMPLE_ORDER,
    "total_to_tp_regions": int(all_tp.shape[0]),
    "total_to_fp_regions": int(all_fp.shape[0]),
    "loh_eligibility": {
        row["sample"] + "_" + row["subset"]: {
            "old_pct": row["old_loh_elig_pct"], "new_pct": row["new_loh_elig_pct"],
            "delta_pct": row["delta_pct"], "fold_change": row["fold_change"]
        }
        for _, row in elig_df.iterrows()
    },
    "aggregated_enrichment_by_tier": {
        row["subset"]: {"tp_loh_frac": row["tp_loh_frac"], "fp_loh_frac": row["fp_loh_frac"],
                        "enrichment": row["enrichment"], "fisher_p": row["fisher_p"]}
        for _, row in tier_df.iterrows()
    },
    "output_dir": str(OUT_DIR),
}
with open(OUT_DIR / "summary.json", "w") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

print(f"\n[Done] 輸出目錄：{OUT_DIR}")
print(f"  loh_eligibility_before_after.tsv")
print(f"  loh_enrichment_by_tier.tsv")
print(f"  loh_enrichment_by_sample_tier.tsv")
print(f"  cross_mode_comparison.tsv")
print(f"  summary.json")
print(f"  figures/ (4 figures)")
