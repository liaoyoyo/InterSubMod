#!/usr/bin/env python3
"""
LOH Round 3 Analysis:
  Core 1 — TO LOH-like × HP0 Stratification
            （低 HP0 的 TO LOH-like 是否是更可信的 TP/FP 分離訊號？）
  Core 2 — LOH × Methylation Joint Analysis
            （Tier A LOH-like + HPMergedSig 是否比單獨 LOH-like 更有判別力？）
  Core 3 — Tier Threshold Sensitivity
            （Tier A=20/30/40/50 的 enrichment 是否穩健？）
  Core 4 — Paired Tier A LOH-like FP Filter Simulation
            （對 H2009、HCC1937 的 filter 效果與 F1 影響估算）

生成時間: 2026-03-27
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
import json
from datetime import datetime
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
OUT_DIR   = Path(
    os.environ.get(
        "LOH_ROUND3_OUT_DIR",
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round3_methyl_hp0_filter",
    )
)
FIG_DIR   = OUT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

GENERATED_AT = datetime.now().strftime("%Y-%m-%d %H:%M")

# Round 2 確認的 Tier 定義（固定）
def assign_tier(n):
    if pd.isna(n): return "C0"
    n = int(n)
    if n == 0: return "C0"
    elif n < 10: return "C"
    elif n < 30: return "B"
    else: return "A"

TIER_ORDER = ["C0","C","B","A"]
TIER_COLORS = {"C0":"#d9d9d9","C":"#fc8d59","B":"#fee090","A":"#4575b4"}

# HP0 strata（Round 3 新增）
def assign_hp0_stratum(r):
    if pd.isna(r): return "NA"
    if r < 0.05:  return "Low (<5%)"
    if r < 0.15:  return "Mid (5-15%)"
    return "High (≥15%)"

HP0_STRATA_ORDER = ["Low (<5%)","Mid (5-15%)","High (≥15%)"]
HP0_COLORS = {"Low (<5%)":"#2166ac","Mid (5-15%)":"#abdda4","High (≥15%)":"#d73027"}

def compute_enrichment_row(tp_total, fp_total, tp_loh_n, fp_loh_n):
    tp_frac = tp_loh_n / tp_total if tp_total > 0 else np.nan
    fp_frac = fp_loh_n / fp_total if fp_total > 0 else np.nan
    enrich  = fp_frac / tp_frac if (tp_frac and tp_frac > 0) else np.nan
    if tp_total > 0 and fp_total > 0:
        a,b,c,d = int(fp_loh_n), int(fp_total-fp_loh_n), int(tp_loh_n), int(tp_total-tp_loh_n)
        try: _, pval = stats.fisher_exact([[a,b],[c,d]])
        except: pval = np.nan
    else: pval = np.nan
    return tp_frac, fp_frac, enrich, pval

def enrichment_df(grp_df):
    tp = grp_df[grp_df["truth_label"]=="TP"]
    fp = grp_df[grp_df["truth_label"]=="FP"]
    tp_n, fp_n = len(tp), len(fp)
    tp_loh = tp["core_loh_like"].sum()
    fp_loh = fp["core_loh_like"].sum()
    tp_f, fp_f, enr, p = compute_enrichment_row(tp_n, fp_n, tp_loh, fp_loh)
    return pd.Series({"tp_total":tp_n,"fp_total":fp_n,
                      "tp_loh_n":tp_loh,"fp_loh_n":fp_loh,
                      "tp_loh_frac":tp_f,"fp_loh_frac":fp_f,
                      "fp_tp_enrichment":enr,"fisher_pval":p})

# ============================================================
# 載入資料
# ============================================================
print("[1/9] 載入資料 ...")
df = pd.read_csv(ROUND1_WS/"all_region_rows.tsv.gz", sep="\t", low_memory=False)
num_cols = ["effective_hp_reads","hp_ratio_core","hp0_ratio","hp3_ratio",
            "NHP0","NHP3","NumReads","UnassignedAffinity",
            "HPMergedDelta","HPMergedP","HPFineF","HPFineP",
            "AlleleDelta","AlleleP","CramersV","PairwiseMedianDist",
            "HeuristicScore","GlobalP","caller_af"]
for c in num_cols:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

bool_cols = ["HPMergedSig","HPFineSig","AlleleSig","core_loh_like",
             "PassedGating","Significant"]
for c in bool_cols:
    if c in df.columns:
        df[c] = df[c].astype(str).str.strip().str.lower().isin(["true","1","yes"])

df["hp_support_tier"]  = df["effective_hp_reads"].apply(assign_tier)
df["hp0_stratum"]      = df["hp0_ratio"].apply(assign_hp0_stratum)
df_tf = df[df["truth_label"].isin(["TP","FP"])].copy()
print(f"  總計 {len(df_tf):,} TP/FP rows (TP={df_tf['truth_label'].eq('TP').sum():,}, FP={df_tf['truth_label'].eq('FP').sum():,})")

# Tier A only subset
df_tA   = df_tf[df_tf["hp_support_tier"] == "A"].copy()
df_tA_L = df_tA[df_tA["core_loh_like"] == True].copy()
print(f"  Tier A LOH-like: {len(df_tA_L):,} rows")

# ============================================================
# Core 1: TO LOH-like × HP0 Stratification
# ============================================================
print("\n[2/9] Core 1: TO LOH-like × HP0 Stratification ...")

to_L = df_tA_L[df_tA_L["mode"] == "to"].copy()
pa_L = df_tA_L[df_tA_L["mode"] == "paired"].copy()

# 1A: 全域 TO Tier A LOH-like by HP0 stratum
def hp0_strat_enrichment(sub, mode_label):
    rows = []
    for stratum in HP0_STRATA_ORDER:
        grp = sub[sub["hp0_stratum"] == stratum]
        if len(grp) < 5: continue
        tp = grp[grp["truth_label"]=="TP"]
        fp = grp[grp["truth_label"]=="FP"]
        _, _, enr, p = compute_enrichment_row(len(tp),len(fp),len(tp),len(fp))  # all are LOH-like
        tp_f = len(tp)/len(grp) if len(grp) > 0 else np.nan
        fp_f = len(fp)/len(grp) if len(grp) > 0 else np.nan
        # enrichment = FP fraction / TP fraction in stratum
        total_tp_in_mode = (sub["truth_label"]=="TP").sum()
        total_fp_in_mode = (sub["truth_label"]=="FP").sum()
        if total_tp_in_mode > 0 and total_fp_in_mode > 0:
            enr_val = (len(fp)/total_fp_in_mode) / (len(tp)/total_tp_in_mode) if len(tp)>0 else np.nan
        else:
            enr_val = np.nan
        rows.append({"mode":mode_label,"hp0_stratum":stratum,
                     "n_tp":len(tp),"n_fp":len(fp),"n_total":len(grp),
                     "tp_frac_in_stratum":tp_f,"fp_frac_in_stratum":fp_f,
                     "hp0_mean":grp["hp0_ratio"].mean(),
                     "hp0_median":grp["hp0_ratio"].median()})
    return rows

c1_rows = hp0_strat_enrichment(to_L, "to") + hp0_strat_enrichment(pa_L, "paired")
c1_global = pd.DataFrame(c1_rows)
c1_global.to_csv(OUT_DIR/"core1_to_hp0_strat_global.tsv", sep="\t", index=False)

# 1B: TO 各樣本 HP0 stratification 的 TP/FP ratio
print("  TO Tier A LOH-like HP0 stratification per sample:")
c1_sample_rows = []
for sample, sgrp in to_L.groupby("sample"):
    for stratum in HP0_STRATA_ORDER:
        sg = sgrp[sgrp["hp0_stratum"]==stratum]
        if len(sg) < 5: continue
        n_tp = (sg["truth_label"]=="TP").sum()
        n_fp = (sg["truth_label"]=="FP").sum()
        frac_tp = n_tp / len(sg) if len(sg) > 0 else np.nan
        c1_sample_rows.append({"sample":sample,"hp0_stratum":stratum,
                                "n_tp":n_tp,"n_fp":n_fp,"n_total":len(sg),
                                "frac_tp":frac_tp,"frac_fp":1-frac_tp if not np.isnan(frac_tp) else np.nan,
                                "hp0_median":sg["hp0_ratio"].median()})
        print(f"    {sample}/{stratum}: TP={n_tp}, FP={n_fp}, TP%={frac_tp:.3f}")

c1_sample = pd.DataFrame(c1_sample_rows)
c1_sample.to_csv(OUT_DIR/"core1_to_hp0_strat_by_sample.tsv", sep="\t", index=False)

# 1C: HP0 threshold sweep (TO Tier A LOH-like)
print("\n  HP0 threshold sweep:")
hp0_thresholds = [0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30]
sweep_rows = []
for thresh in hp0_thresholds:
    sub_low  = to_L[to_L["hp0_ratio"] < thresh]
    sub_high = to_L[to_L["hp0_ratio"] >= thresh]
    for label, sub in [("low_hp0", sub_low), ("high_hp0", sub_high)]:
        n_tp = (sub["truth_label"]=="TP").sum()
        n_fp = (sub["truth_label"]=="FP").sum()
        frac_tp = n_tp/(n_tp+n_fp) if (n_tp+n_fp) > 0 else np.nan
        sweep_rows.append({"hp0_thresh":thresh,"group":label,
                           "n_tp":n_tp,"n_fp":n_fp,"frac_tp":frac_tp,
                           "n_total":len(sub)})
    print(f"  thresh={thresh:.2f}: low(n={len(sub_low):,} TP%={(to_L[to_L['hp0_ratio']<thresh]['truth_label']=='TP').mean():.3f}) | high(n={len(sub_high):,} TP%={(to_L[to_L['hp0_ratio']>=thresh]['truth_label']=='TP').mean():.3f})")

pd.DataFrame(sweep_rows).to_csv(OUT_DIR/"core1_hp0_threshold_sweep.tsv", sep="\t", index=False)

# ============================================================
# Core 2: LOH × Methylation Joint Analysis
# ============================================================
print("\n[3/9] Core 2: LOH × Methylation Joint Analysis ...")

# 2A: Tier A LOH-like TP vs FP 的 HPMergedSig / AlleleSig 分佈
print("  Tier A LOH-like — HPMergedSig / AlleleSig TP vs FP:")
methyl_sig_rows = []
for mode, mgrp in df_tA_L.groupby("mode"):
    for truth, tgrp in mgrp.groupby("truth_label"):
        n = len(tgrp)
        hp_sig_frac = tgrp["HPMergedSig"].mean()
        al_sig_frac = tgrp["AlleleSig"].mean() if "AlleleSig" in tgrp.columns else np.nan
        hp_delta_med = tgrp["HPMergedDelta"].median()
        hp_delta_abs = tgrp["HPMergedDelta"].abs().median()
        allele_delta_med = tgrp["AlleleDelta"].median() if "AlleleDelta" in tgrp.columns else np.nan
        methyl_sig_rows.append({"mode":mode,"truth":truth,"n":n,
                                 "HPMergedSig_frac":hp_sig_frac,
                                 "AlleleSig_frac":al_sig_frac,
                                 "HPMergedDelta_median":hp_delta_med,
                                 "HPMergedDelta_abs_median":hp_delta_abs,
                                 "AlleleDelta_median":allele_delta_med})
        print(f"    {mode}/{truth}: HPMergedSig={hp_sig_frac:.3f}, AlleleSig={al_sig_frac:.3f}, |HP_delta|_med={hp_delta_abs:.4f}")

pd.DataFrame(methyl_sig_rows).to_csv(OUT_DIR/"core2_methyl_sig_by_mode_truth.tsv", sep="\t", index=False)

# 2B: Mann-Whitney — HPMergedDelta abs 差異 (TP vs FP in Tier A LOH-like)
print("\n  Mann-Whitney: HPMergedDelta abs (TP vs FP, Tier A LOH-like):")
mw_methyl = []
for mode, mgrp in df_tA_L.groupby("mode"):
    tp_vals = mgrp[mgrp["truth_label"]=="TP"]["HPMergedDelta"].abs().dropna()
    fp_vals = mgrp[mgrp["truth_label"]=="FP"]["HPMergedDelta"].abs().dropna()
    if len(tp_vals)>10 and len(fp_vals)>10:
        stat, pval = stats.mannwhitneyu(tp_vals, fp_vals, alternative="two-sided")
        print(f"    {mode}: TP abs_median={tp_vals.median():.4f}, FP abs_median={fp_vals.median():.4f}, p={pval:.2e}")
        mw_methyl.append({"mode":mode,"tp_abs_median":tp_vals.median(),"fp_abs_median":fp_vals.median(),
                          "mw_pval":pval,"n_tp":len(tp_vals),"n_fp":len(fp_vals)})
pd.DataFrame(mw_methyl).to_csv(OUT_DIR/"core2_methyl_mw_test.tsv", sep="\t", index=False)

# 2C: Joint feature — Tier A LOH-like × HPMergedSig 的 TP/FP enrichment
print("\n  Joint feature: Tier A LOH-like × HPMergedSig:")
joint_rows = []
for mode, mgrp in df_tA.groupby("mode"):
    for loh_l in [True, False]:
        for hp_sig in [True, False]:
            sg = mgrp[(mgrp["core_loh_like"]==loh_l) & (mgrp["HPMergedSig"]==hp_sig)]
            n_tp = (sg["truth_label"]=="TP").sum()
            n_fp = (sg["truth_label"]=="FP").sum()
            frac_fp = n_fp/(n_tp+n_fp) if (n_tp+n_fp)>0 else np.nan
            label = f"LOH={loh_l}/HPSig={hp_sig}"
            print(f"    {mode}/{label}: TP={n_tp:,}, FP={n_fp:,}, FP%={frac_fp:.4f}")
            joint_rows.append({"mode":mode,"loh_like":loh_l,"hp_merged_sig":hp_sig,
                                "n_tp":n_tp,"n_fp":n_fp,"fp_frac":frac_fp,
                                "n_total":len(sg)})
pd.DataFrame(joint_rows).to_csv(OUT_DIR/"core2_joint_loh_methyl_sig.tsv", sep="\t", index=False)

# 2D: Per-sample (paired Tier A) 中 HPMergedSig TP vs FP
print("\n  Paired Tier A LOH-like HPMergedSig per sample:")
pa_tA_L = df_tA_L[df_tA_L["mode"]=="paired"]
per_sample_methyl = []
for sample, sgrp in pa_tA_L.groupby("sample"):
    tp = sgrp[sgrp["truth_label"]=="TP"]
    fp = sgrp[sgrp["truth_label"]=="FP"]
    tp_hpsig = tp["HPMergedSig"].mean() if len(tp)>0 else np.nan
    fp_hpsig = fp["HPMergedSig"].mean() if len(fp)>0 else np.nan
    tp_aldelta = tp["AlleleDelta"].abs().median() if len(tp)>0 else np.nan
    fp_aldelta = fp["AlleleDelta"].abs().median() if len(fp)>0 else np.nan
    print(f"    {sample}: TP HPSig={tp_hpsig:.3f}, FP HPSig={fp_hpsig:.3f} | TP |AlleleDelta|={tp_aldelta:.4f}, FP |AlleleDelta|={fp_aldelta:.4f}")
    per_sample_methyl.append({"sample":sample,
                               "tp_n":len(tp),"fp_n":len(fp),
                               "tp_HPMergedSig_frac":tp_hpsig,"fp_HPMergedSig_frac":fp_hpsig,
                               "tp_AlleleDelta_abs_med":tp_aldelta,"fp_AlleleDelta_abs_med":fp_aldelta})
pd.DataFrame(per_sample_methyl).to_csv(OUT_DIR/"core2_methyl_sig_paired_per_sample.tsv", sep="\t", index=False)

# 2E: TO Tier A LOH-like HP0 stratum × HPMergedSig
print("\n  TO Tier A LOH-like: HP0 stratum × HPMergedSig:")
hp0_methyl_rows = []
for stratum in HP0_STRATA_ORDER:
    sg = to_L[to_L["hp0_stratum"]==stratum]
    if len(sg) < 10: continue
    tp = sg[sg["truth_label"]=="TP"]
    fp = sg[sg["truth_label"]=="FP"]
    tp_hpsig = tp["HPMergedSig"].mean() if len(tp)>0 else np.nan
    fp_hpsig = fp["HPMergedSig"].mean() if len(fp)>0 else np.nan
    tp_abs = tp["HPMergedDelta"].abs().median() if len(tp)>0 else np.nan
    fp_abs = fp["HPMergedDelta"].abs().median() if len(fp)>0 else np.nan
    print(f"    {stratum}: TP HPSig={tp_hpsig:.3f}(n={len(tp)}), FP HPSig={fp_hpsig:.3f}(n={len(fp)})")
    print(f"             TP |delta|_med={tp_abs:.4f}, FP |delta|_med={fp_abs:.4f}")
    hp0_methyl_rows.append({"hp0_stratum":stratum,"n_tp":len(tp),"n_fp":len(fp),
                             "tp_HPMergedSig_frac":tp_hpsig,"fp_HPMergedSig_frac":fp_hpsig,
                             "tp_HPMergedDelta_abs_med":tp_abs,"fp_HPMergedDelta_abs_med":fp_abs})
pd.DataFrame(hp0_methyl_rows).to_csv(OUT_DIR/"core2_hp0_stratum_methyl.tsv", sep="\t", index=False)

# ============================================================
# Core 3: Tier Threshold Sensitivity
# ============================================================
print("\n[4/9] Core 3: Tier Threshold Sensitivity ...")

# Sweep Tier A threshold: 10, 15, 20, 25, 30, 40, 50
a_thresholds = [10, 15, 20, 25, 30, 40, 50]
sens_rows = []
for thresh in a_thresholds:
    for mode in ["paired", "to"]:
        sub = df_tf[(df_tf["mode"]==mode) & (df_tf["effective_hp_reads"] >= thresh) & (df_tf["effective_hp_reads"] > 0)]
        tp_n = (sub["truth_label"]=="TP").sum()
        fp_n = (sub["truth_label"]=="FP").sum()
        tp_loh = sub[sub["truth_label"]=="TP"]["core_loh_like"].sum()
        fp_loh = sub[sub["truth_label"]=="FP"]["core_loh_like"].sum()
        tp_f, fp_f, enr, p = compute_enrichment_row(tp_n, fp_n, tp_loh, fp_loh)
        sens_rows.append({"tier_a_thresh":thresh,"mode":mode,
                          "tp_total":tp_n,"fp_total":fp_n,
                          "tp_loh_frac":tp_f,"fp_loh_frac":fp_f,
                          "fp_tp_enrichment":enr,"fisher_pval":p})
        print(f"  A≥{thresh:2d}/{mode}: enrichment={enr:.4f}, p={p:.2e} (TP={tp_n:,}, FP={fp_n:,})")

sens_df = pd.DataFrame(sens_rows)
sens_df.to_csv(OUT_DIR/"core3_tier_threshold_sensitivity.tsv", sep="\t", index=False)

# ============================================================
# Core 4: Paired Tier A LOH-like FP Filter Simulation
# ============================================================
print("\n[5/9] Core 4: Paired Tier A LOH-like Filter Simulation ...")

# Focus on H2009 and HCC1937（Round 2 確認有強 enrichment 且顯著）
# Additional: simulate across all samples
filter_sim_rows = []

for sample in df_tf["sample"].unique():
    sdf = df_tf[(df_tf["mode"]=="paired") & (df_tf["sample"]==sample)]
    if len(sdf) < 10: continue

    tp_all  = (sdf["truth_label"]=="TP").sum()
    fp_all  = (sdf["truth_label"]=="FP").sum()
    if fp_all == 0: continue

    # 基線 precision / recall
    baseline_prec = tp_all / (tp_all + fp_all)
    baseline_recall = 1.0  # 全部 TP 都保留（baseline，相對自身）

    for filter_label, mask_fn in [
        ("TierA_LOH",     lambda d: (d["hp_support_tier"]=="A") & (d["core_loh_like"]==True)),
        ("TierA_LOH_HP0low", lambda d: (d["hp_support_tier"]=="A") & (d["core_loh_like"]==True) & (d["hp0_ratio"] < 0.05)),
        ("TierA_LOH_HPSig", lambda d: (d["hp_support_tier"]=="A") & (d["core_loh_like"]==True) & (d["HPMergedSig"]==True)),
    ]:
        removed = sdf[mask_fn(sdf)]
        removed_tp = (removed["truth_label"]=="TP").sum()
        removed_fp = (removed["truth_label"]=="FP").sum()
        remain_tp = tp_all - removed_tp
        remain_fp = fp_all - removed_fp
        new_prec = remain_tp / (remain_tp + remain_fp) if (remain_tp + remain_fp) > 0 else np.nan
        new_recall = remain_tp / tp_all if tp_all > 0 else np.nan
        f1 = 2*new_prec*new_recall / (new_prec+new_recall) if (new_prec and new_recall and (new_prec+new_recall)>0) else np.nan
        base_f1 = 2*baseline_prec*1.0 / (baseline_prec+1.0) if baseline_prec else np.nan
        filter_sim_rows.append({
            "sample": sample, "filter": filter_label,
            "tp_all": tp_all, "fp_all": fp_all,
            "removed_tp": removed_tp, "removed_fp": removed_fp,
            "removed_tp_pct": removed_tp/tp_all if tp_all>0 else np.nan,
            "removed_fp_pct": removed_fp/fp_all if fp_all>0 else np.nan,
            "new_precision": new_prec,
            "new_recall": new_recall,
            "new_f1": f1,
            "baseline_f1": base_f1,
            "f1_delta": f1-base_f1 if (f1 and base_f1) else np.nan,
        })

filter_sim = pd.DataFrame(filter_sim_rows)
filter_sim.to_csv(OUT_DIR/"core4_filter_simulation.tsv", sep="\t", index=False)

print("\n  Filter Simulation (paired, H2009 & HCC1937):")
key_samples = ["H2009","HCC1937","HCC1954","COLO829","HCC1395","HCC1395_DORADO"]
for s in key_samples:
    sub = filter_sim[filter_sim["sample"]==s]
    if len(sub)==0: continue
    print(f"\n  {s}:")
    for _, row in sub.iterrows():
        print(f"    {row['filter']}: rm_TP={int(row['removed_tp'])}({row['removed_tp_pct']:.2%}), rm_FP={int(row['removed_fp'])}({row['removed_fp_pct']:.2%}), F1_delta={row['f1_delta']:+.4f} ({row['baseline_f1']:.4f}→{row['new_f1']:.4f})")

# ============================================================
# 視覺化
# ============================================================
print("\n[6/9] 視覺化 ...")

# Fig 1: TO Tier A LOH-like — HP0 threshold sweep (TP% at different thresholds)
sweep_df = pd.read_csv(OUT_DIR/"core1_hp0_threshold_sweep.tsv", sep="\t")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax = axes[0]
low_tp  = sweep_df[sweep_df["group"]=="low_hp0"][["hp0_thresh","frac_tp","n_total"]].set_index("hp0_thresh")
high_tp = sweep_df[sweep_df["group"]=="high_hp0"][["hp0_thresh","frac_tp","n_total"]].set_index("hp0_thresh")
ax.plot(low_tp.index, low_tp["frac_tp"], "o-", color="#2166ac", linewidth=2, label="low HP0 (< thresh)")
ax.plot(high_tp.index, high_tp["frac_tp"], "s--", color="#d73027", linewidth=2, label="high HP0 (≥ thresh)")
for x, (_, row) in zip(low_tp.index, low_tp.iterrows()):
    ax.annotate(f"n={int(row['n_total'])}", (x, row['frac_tp']), xytext=(0,6),
                textcoords="offset points", ha="center", fontsize=6.5, color="#2166ac")
for x, (_, row) in zip(high_tp.index, high_tp.iterrows()):
    ax.annotate(f"n={int(row['n_total'])}", (x, row['frac_tp']), xytext=(0,-12),
                textcoords="offset points", ha="center", fontsize=6.5, color="#d73027")
ax.set_xlabel("HP0 Threshold")
ax.set_ylabel("TP Fraction")
ax.set_title("TO Tier A LOH-like: TP Fraction\nby HP0 Threshold (low vs high split)")
ax.legend()
ax.set_ylim(0, 1)
ax.axhline(to_L["truth_label"].eq("TP").mean(), color="gray", linestyle=":", label="Overall TP%")

ax = axes[1]
for sample in c1_sample["sample"].unique():
    sgrp = c1_sample[c1_sample["sample"]==sample]
    if len(sgrp) < 2: continue
    sgrp_ord = sgrp.set_index("hp0_stratum").reindex(HP0_STRATA_ORDER)
    ax.plot(range(len(HP0_STRATA_ORDER)), sgrp_ord["frac_tp"], "o-", label=sample, alpha=0.8)
ax.set_xticks(range(len(HP0_STRATA_ORDER)))
ax.set_xticklabels(HP0_STRATA_ORDER, fontsize=9)
ax.set_ylabel("TP Fraction in Stratum")
ax.set_title("TO Tier A LOH-like: TP Fraction\nby HP0 Stratum (per sample)")
ax.legend(fontsize=7)
ax.set_ylim(0, 1)

plt.tight_layout()
plt.savefig(FIG_DIR/"fig01_to_loh_hp0_strat.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 01 saved")

# Fig 2: LOH × HPMergedSig joint heatmap (FP fraction)
joint_df = pd.read_csv(OUT_DIR/"core2_joint_loh_methyl_sig.tsv", sep="\t")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, mode in zip(axes, ["paired","to"]):
    sub = joint_df[joint_df["mode"]==mode].copy()
    sub["label"] = sub.apply(lambda r: f"LOH={'Y' if r['loh_like'] else 'N'}\nHPSig={'Y' if r['hp_merged_sig'] else 'N'}", axis=1)
    sub["fp_pct"] = sub["fp_frac"] * 100
    colors = sub["fp_pct"].apply(lambda v: plt.cm.RdYlGn_r(v/20 if not np.isnan(v) else 0))
    bars = ax.bar(sub["label"], sub["fp_pct"],
                  color=[plt.cm.RdYlGn_r(v/20) for v in sub["fp_pct"].fillna(0)],
                  edgecolor="k", linewidth=0.8)
    for bar, row in zip(bars, sub.itertuples()):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.1,
                f"TP={int(row.n_tp):,}\nFP={int(row.n_fp):,}",
                ha="center", va="bottom", fontsize=7)
    ax.set_ylabel("FP% in group")
    ax.set_title(f"Mode: {mode.upper()}\nFP Fraction by LOH × HPMergedSig (Tier A)")
    ax.set_ylim(0, max(20, sub["fp_pct"].max()*1.3))

plt.tight_layout()
plt.savefig(FIG_DIR/"fig02_loh_methyl_joint.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 02 saved")

# Fig 3: HPMergedDelta distribution (TP vs FP, Tier A LOH-like, paired)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, mode in zip(axes, ["paired","to"]):
    sg = df_tA_L[df_tA_L["mode"]==mode]
    tp_delta = sg[sg["truth_label"]=="TP"]["HPMergedDelta"].dropna()
    fp_delta = sg[sg["truth_label"]=="FP"]["HPMergedDelta"].dropna()
    bins = np.linspace(-0.6, 0.6, 41)
    ax.hist(tp_delta.clip(-0.6,0.6), bins=bins, alpha=0.6, density=True,
            label=f"TP (n={len(tp_delta):,})", color="#2166ac")
    ax.hist(fp_delta.clip(-0.6,0.6), bins=bins, alpha=0.6, density=True,
            label=f"FP (n={len(fp_delta):,})", color="#d73027")
    ax.axvline(0, color="k", linewidth=1, linestyle="--")
    ax.axvline(tp_delta.median(), color="#2166ac", linewidth=1.5, linestyle=":")
    ax.axvline(fp_delta.median(), color="#d73027", linewidth=1.5, linestyle=":")
    ax.set_xlabel("HPMergedDelta (HP1-family vs HP2-family methylation diff)")
    ax.set_ylabel("Density")
    ax.set_title(f"{mode.upper()} Tier A LOH-like\nHPMergedDelta Distribution (TP vs FP)")
    ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(FIG_DIR/"fig03_hp_merged_delta_loh.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 03 saved")

# Fig 4: Tier threshold sensitivity
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, mode in zip(axes, ["paired","to"]):
    sub = sens_df[sens_df["mode"]==mode].sort_values("tier_a_thresh")
    ax.plot(sub["tier_a_thresh"], sub["fp_tp_enrichment"], "o-", color="#4575b4", linewidth=2)
    for _, row in sub.iterrows():
        ax.annotate(f"{row['fp_tp_enrichment']:.3f}", (row["tier_a_thresh"], row["fp_tp_enrichment"]),
                    xytext=(0,8), textcoords="offset points", ha="center", fontsize=8)
    ax.axhline(1.0, color="red", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xlabel("Tier A Minimum effective_hp_reads Threshold")
    ax.set_ylabel("LOH-like FP/TP Enrichment")
    ax.set_title(f"Mode: {mode.upper()}\nLOH Enrichment vs Tier A Threshold")
    ax.set_ylim(0.5, max(2.5, sub["fp_tp_enrichment"].max()*1.2))
    ax.set_xticks(sub["tier_a_thresh"])

plt.tight_layout()
plt.savefig(FIG_DIR/"fig04_tier_threshold_sensitivity.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 04 saved")

# Fig 5: Filter simulation — F1 delta by sample and filter type
fig, ax = plt.subplots(figsize=(14, 6))
filt_pivot = filter_sim.pivot_table(index="sample", columns="filter", values="f1_delta")
filt_pivot_sorted = filt_pivot.reindex(
    filt_pivot.fillna(0).abs().sum(axis=1).sort_values(ascending=False).index
)
x = np.arange(len(filt_pivot_sorted))
width = 0.25
cols = ["TierA_LOH","TierA_LOH_HP0low","TierA_LOH_HPSig"]
filt_colors = {"TierA_LOH":"#4575b4","TierA_LOH_HP0low":"#2ca25f","TierA_LOH_HPSig":"#d73027"}
for i, col in enumerate(cols):
    vals = [filt_pivot_sorted.loc[s, col] if col in filt_pivot_sorted.columns and s in filt_pivot_sorted.index else 0
            for s in filt_pivot_sorted.index]
    ax.bar(x + (i-1)*width, vals, width,
           label=col.replace("TierA_",""), color=filt_colors[col], alpha=0.8, edgecolor="k")
ax.axhline(0, color="k", linewidth=1.5)
ax.set_xticks(x)
ax.set_xticklabels(filt_pivot_sorted.index, rotation=15, fontsize=9)
ax.set_ylabel("F1 Score Delta (new - baseline)")
ax.set_title("Paired Mode Filter Simulation: F1 Delta by Sample\n(positive = F1 improves)")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(FIG_DIR/"fig05_filter_simulation_f1_delta.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 05 saved")

# Fig 6: Paired Tier A LOH-like HPMergedSig fraction (TP vs FP) per sample
per_smp_df = pd.read_csv(OUT_DIR/"core2_methyl_sig_paired_per_sample.tsv", sep="\t")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
x = np.arange(len(per_smp_df))
ax = axes[0]
w = 0.35
ax.bar(x-w/2, per_smp_df["tp_HPMergedSig_frac"], w, label="TP", color="#2166ac", alpha=0.8)
ax.bar(x+w/2, per_smp_df["fp_HPMergedSig_frac"], w, label="FP", color="#d73027", alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(per_smp_df["sample"], rotation=20, fontsize=9)
ax.set_ylabel("HPMergedSig Fraction")
ax.set_title("Paired Tier A LOH-like\nHPMergedSig Fraction (TP vs FP)")
ax.legend()
ax.set_ylim(0, 1)
for xi, row in zip(x, per_smp_df.itertuples()):
    ax.text(xi-w/2, row.tp_HPMergedSig_frac+0.02, f"n={row.tp_n}", ha="center", fontsize=6, color="#2166ac")
    ax.text(xi+w/2, row.fp_HPMergedSig_frac+0.02, f"n={row.fp_n}", ha="center", fontsize=6, color="#d73027")

ax = axes[1]
ax.bar(x-w/2, per_smp_df["tp_AlleleDelta_abs_med"], w, label="TP", color="#2166ac", alpha=0.8)
ax.bar(x+w/2, per_smp_df["fp_AlleleDelta_abs_med"], w, label="FP", color="#d73027", alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(per_smp_df["sample"], rotation=20, fontsize=9)
ax.set_ylabel("|AlleleDelta| Median")
ax.set_title("Paired Tier A LOH-like\n|AlleleDelta| Median (TP vs FP)")
ax.legend()

plt.tight_layout()
plt.savefig(FIG_DIR/"fig06_methyl_sig_per_sample.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Fig 06 saved")

# Fig 7: HP0 stratum × HPMergedSig (TO Tier A LOH-like)
hp0_meth_df = pd.read_csv(OUT_DIR/"core2_hp0_stratum_methyl.tsv", sep="\t")
if len(hp0_meth_df) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    x = np.arange(len(hp0_meth_df))
    ax = axes[0]
    w = 0.35
    ax.bar(x-w/2, hp0_meth_df["tp_HPMergedSig_frac"], w, label="TP", color="#2166ac", alpha=0.8)
    ax.bar(x+w/2, hp0_meth_df["fp_HPMergedSig_frac"], w, label="FP", color="#d73027", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(hp0_meth_df["hp0_stratum"], fontsize=9)
    ax.set_ylabel("HPMergedSig Fraction")
    ax.set_title("TO Tier A LOH-like\nHPMergedSig by HP0 Stratum (TP vs FP)")
    ax.legend()
    ax.set_ylim(0, 1)

    ax = axes[1]
    ax.bar(x-w/2, hp0_meth_df["tp_HPMergedDelta_abs_med"], w, label="TP", color="#2166ac", alpha=0.8)
    ax.bar(x+w/2, hp0_meth_df["fp_HPMergedDelta_abs_med"], w, label="FP", color="#d73027", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(hp0_meth_df["hp0_stratum"], fontsize=9)
    ax.set_ylabel("|HPMergedDelta| Median")
    ax.set_title("TO Tier A LOH-like\n|HPMergedDelta| by HP0 Stratum (TP vs FP)")
    ax.legend()

    plt.tight_layout()
    plt.savefig(FIG_DIR/"fig07_hp0_stratum_methyl.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Fig 07 saved")

# ============================================================
# 彙整 Round 3 Summary 數字
# ============================================================
print("\n[7/9] 彙整 key numbers ...")

print("\n=== Round 3 Key Numbers ===")

# Core 1
print("\nCore 1 — TO Tier A LOH-like HP0 threshold sweep:")
sw = pd.read_csv(OUT_DIR/"core1_hp0_threshold_sweep.tsv", sep="\t")
for t in [0.05, 0.10, 0.15]:
    lo = sw[(sw["hp0_thresh"]==t)&(sw["group"]=="low_hp0")]
    hi = sw[(sw["hp0_thresh"]==t)&(sw["group"]=="high_hp0")]
    if len(lo) > 0 and len(hi) > 0:
        lo_tp = lo.iloc[0]["frac_tp"]; lo_n = lo.iloc[0]["n_total"]
        hi_tp = hi.iloc[0]["frac_tp"]; hi_n = hi.iloc[0]["n_total"]
        print(f"  thresh={t:.2f}: low(n={lo_n:.0f}, TP%={lo_tp:.3f}) vs high(n={hi_n:.0f}, TP%={hi_tp:.3f}), diff={lo_tp-hi_tp:+.3f}")

# Core 2
print("\nCore 2 — Joint LOH × HPMergedSig (paired, Tier A):")
jd = pd.read_csv(OUT_DIR/"core2_joint_loh_methyl_sig.tsv", sep="\t")
for _, row in jd[jd["mode"]=="paired"].iterrows():
    tag = f"LOH={'Y' if row['loh_like'] else 'N'}/HPSig={'Y' if row['hp_merged_sig'] else 'N'}"
    print(f"  {tag}: TP={int(row['n_tp']):,}, FP={int(row['n_fp']):,}, FP%={row['fp_frac']:.4f}")

# Core 3
print("\nCore 3 — Tier A threshold sensitivity (paired):")
for _, row in sens_df[sens_df["mode"]=="paired"].iterrows():
    print(f"  A≥{int(row['tier_a_thresh'])}: enrichment={row['fp_tp_enrichment']:.4f}, p={row['fisher_pval']:.2e}")

# Core 4
print("\nCore 4 — Filter simulation H2009 & HCC1937:")
for s in ["H2009","HCC1937"]:
    sub = filter_sim[filter_sim["sample"]==s]
    print(f"\n  {s}:")
    for _, row in sub.iterrows():
        print(f"    {row['filter']}: rm_TP={int(row['removed_tp'])}({row['removed_tp_pct']:.1%}) rm_FP={int(row['removed_fp'])}({row['removed_fp_pct']:.1%}) F1={row['new_f1']:.4f}({row['f1_delta']:+.4f})")

print(f"\n完成！輸出位置：{OUT_DIR}")
