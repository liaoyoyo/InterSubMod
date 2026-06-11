#!/usr/bin/env python3
"""LR ablation: drop methylation features — does Cycle 1 ΔF1=+0.02236 survive without ISM methylation?

Pre-registered hypothesis:
  H_ABL_1: 5-feature LR (no methylation) ΔF1 ≈ +0.020 → methylation 沒貢獻 → ISM 對 F1 沒新加價值
  H_ABL_2: 5-feature LR ΔF1 ≈ +0.015 → methylation 貢獻約 +0.005 (與 5th-rank coef 一致)
  H_ABL_3: 5-feature LR ΔF1 ≪ +0.01 → methylation dominant (與 today's coef analysis 矛盾)

Test both V6_off_NG and baseline_off_NG.

Mirror Cycle 1 final_filter_and_verdict.py pipeline:
  - L2 C=1.0, lbfgs, 5-fold StratifiedKFold (seed=42 primary + 5 seeds multi-seed)
  - OOF predict_proba + τ sweep [0.10, 0.95] step 0.01
  - F1 formula: TP/(TP+FP) precision, TP/(TP+FN) recall, FN = 19288 + (30490 - TP_kept)
  - ΔF1 = F1 - CALLER_F1_BASELINE (0.7166)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
STEP1 = REPO / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way"
STEP5 = REPO / "research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot"
MASTER_4WAY = STEP1 / "step1_master_four_way.tsv"
MASTER_AUG = STEP5 / "step5_master_augmented.tsv"
FIG_DIR = REPO / "research/v6_bam_tpfp_hp_loh_cn/figures/baseline_vs_v6"
FIG_DIR.mkdir(parents=True, exist_ok=True)

CALLER_F1_BASELINE = 0.7166
FN_CALLER = 19288
PRIMARY_SEED = 42
MULTI_SEEDS = [42, 7, 13, 2026, 1395]
TAU_GRID = np.round(np.arange(0.10, 0.96, 0.01), 2)
CHOSEN_C = 1.0

# Cycle 1 full feature set (10 features)
FULL_FEATURES = [
    "NG_feature",
    "caller_af",
    "loh_inner_flag",
    "V6_off_meth_HPMergedDelta",
    "V6_off_meth_HPFineF",
    "V6_off_meth_NME_imbalance",
    "V6_off_meth_Epipoly_Delta",
    "V6_off_meth_ClusterPermanovaF",
    "chr8_flag",
    "Coverage_Multiple_imp",
]

# Ablation 5-feature subset (drop 5 methylation features)
NOMETH_FEATURES = [
    "NG_feature",
    "caller_af",
    "loh_inner_flag",
    "chr8_flag",
    "Coverage_Multiple_imp",
]


def compute_metrics(tp_kept, fp_kept):
    tp_removed = 30490 - tp_kept
    new_fn = FN_CALLER + tp_removed
    if (tp_kept + fp_kept) == 0:
        return -CALLER_F1_BASELINE, 0.0, 0.0
    precision = tp_kept / (tp_kept + fp_kept)
    recall = tp_kept / (tp_kept + new_fn)
    if precision + recall == 0:
        return -CALLER_F1_BASELINE, precision, recall
    f1 = 2 * precision * recall / (precision + recall)
    return f1 - CALLER_F1_BASELINE, precision, recall


def prep(df, ng_col):
    df = df.copy()
    df["loh_inner_flag"] = (df["loh_side"] == "Inner").astype(int)
    cov_median = df["Coverage_Multiple"].median()
    df["Coverage_Multiple_imp"] = df["Coverage_Multiple"].fillna(cov_median)
    df["chr8_flag"] = (df["chr"] == "chr8").astype(int)
    df["y"] = (df["label"] == "TP").astype(int)
    df["NG_feature"] = df[ng_col]
    return df


def impute(df, features):
    sub = df.copy()
    for c in features:
        if c not in sub.columns:
            sub[c] = np.nan
        med = float(sub[c].median()) if sub[c].notna().any() else 0.0
        if sub[c].isna().any():
            sub[c] = sub[c].fillna(med)
    return sub


def fit_oof(X, y, C, seed):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof = np.zeros(len(X))
    coefs = []
    for tr, va in skf.split(X, y):
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr])
        Xva = sc.transform(X[va])
        clf = LogisticRegression(max_iter=5000, C=C, solver="lbfgs",
                                 penalty="l2", random_state=seed)
        clf.fit(Xtr, y[tr])
        oof[va] = clf.predict_proba(Xva)[:, 1]
        coefs.append(clf.coef_[0].tolist())
    return oof, np.array(coefs)


def sweep_tau(oof, y):
    best = {"tau": 0.0, "dF1": -1e9, "tp_kept": 0, "fp_kept": 0, "P": 0, "R": 0}
    for tau in TAU_GRID:
        keep = oof >= tau
        tp_kept = int(((y == 1) & keep).sum())
        fp_kept = int(((y == 0) & keep).sum())
        dF1, P, R = compute_metrics(tp_kept, fp_kept)
        if dF1 > best["dF1"]:
            best.update(tau=float(tau), dF1=float(dF1),
                        tp_kept=tp_kept, fp_kept=fp_kept,
                        P=float(P), R=float(R))
    return best


def run_pipeline(df, ng_col, features, label):
    print(f"\n=== {label} (NG={ng_col}, {len(features)} features) ===", flush=True)
    df_prep = prep(df, ng_col)
    sub = impute(df_prep, features)
    X = sub[features].values.astype(float)
    y = sub["y"].values.astype(int)

    # Primary
    oof, coefs = fit_oof(X, y, CHOSEN_C, PRIMARY_SEED)
    primary = sweep_tau(oof, y)
    primary["features"] = features
    primary["mean_coefs"] = list(zip(features, coefs.mean(0).tolist()))
    print(f"  Primary (seed=42): tau*={primary['tau']:.2f}, ΔF1={primary['dF1']:+.5f}", flush=True)

    # Multi-seed
    multi = []
    for s in MULTI_SEEDS:
        oof_s, _ = fit_oof(X, y, CHOSEN_C, s)
        sw = sweep_tau(oof_s, y)
        multi.append({"seed": s, "tau": sw["tau"], "dF1": sw["dF1"]})
        print(f"  seed={s}: tau*={sw['tau']:.2f}, ΔF1={sw['dF1']:+.5f}", flush=True)
    multi_df = pd.DataFrame(multi)
    primary["multi_mean"] = float(multi_df["dF1"].mean())
    primary["multi_std"] = float(multi_df["dF1"].std(ddof=1))
    print(f"  Multi-seed: {primary['multi_mean']:+.5f} ± {primary['multi_std']:.5f}", flush=True)

    return primary, multi_df


def main():
    print("Loading data ...", flush=True)
    four = pd.read_csv(MASTER_4WAY, sep="\t", low_memory=False)
    aug = pd.read_csv(MASTER_AUG, sep="\t", low_memory=False)
    merged = aug.merge(four[["region_id", "baseline_off_NG"]], on="region_id", how="inner")
    print(f"  merged: {merged.shape}", flush=True)

    # 4 runs: 2 NG × 2 feature sets
    results = {}
    multi_records = {}
    for ng_col in ["V6_off_NG", "baseline_off_NG"]:
        for fset_name, fset in [("FULL_10feat", FULL_FEATURES), ("NOMETH_5feat", NOMETH_FEATURES)]:
            tag = f"{ng_col}__{fset_name}"
            sw, md = run_pipeline(merged, ng_col, fset, tag)
            results[tag] = sw
            multi_records[tag] = md

    # Ablation deltas
    print(f"\n=== ABLATION ANALYSIS ===")
    print(f"  V6 FULL 10feat        : {results['V6_off_NG__FULL_10feat']['multi_mean']:+.5f}")
    print(f"  V6 NOMETH 5feat       : {results['V6_off_NG__NOMETH_5feat']['multi_mean']:+.5f}")
    v6_meth_contrib = results['V6_off_NG__FULL_10feat']['multi_mean'] - results['V6_off_NG__NOMETH_5feat']['multi_mean']
    print(f"  V6 methylation contrib: {v6_meth_contrib:+.5f}")
    print()
    print(f"  baseline FULL 10feat  : {results['baseline_off_NG__FULL_10feat']['multi_mean']:+.5f}")
    print(f"  baseline NOMETH 5feat : {results['baseline_off_NG__NOMETH_5feat']['multi_mean']:+.5f}")
    base_meth_contrib = results['baseline_off_NG__FULL_10feat']['multi_mean'] - results['baseline_off_NG__NOMETH_5feat']['multi_mean']
    print(f"  baseline meth contrib : {base_meth_contrib:+.5f}")

    # Verdict
    print(f"\n=== PRE-REGISTERED VERDICT ===")
    nometh_v6 = results['V6_off_NG__NOMETH_5feat']['multi_mean']
    if nometh_v6 >= 0.020:
        verdict = "H_ABL_1: methylation USELESS — 5-feat 仍達 +0.020+ (ISM 對 F1 沒新加價值)"
    elif nometh_v6 >= 0.015:
        verdict = "H_ABL_2: methylation 貢獻 ~+0.005 (與 5th-rank coef +0.75 一致)"
    elif nometh_v6 >= 0.01:
        verdict = "H_ABL_2-partial: methylation 貢獻 +0.005-+0.012"
    else:
        verdict = "H_ABL_3: methylation DOMINANT — 5-feat 跌至 +0.01- (與 coef analysis 矛盾，可能有 interaction effects)"
    print(f"  {verdict}")

    # Save
    out = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "results": {k: {kk: vv for kk, vv in v.items() if kk not in ['mean_coefs']}
                    for k, v in results.items()},
        "feature_coefs": {k: results[k]["mean_coefs"] for k in results},
        "ablation_deltas": {
            "V6_methylation_contribution": v6_meth_contrib,
            "baseline_methylation_contribution": base_meth_contrib,
        },
        "verdict": verdict,
    }
    out_path = STEP1 / "step1_lr_ablation_methylation_findings.json"
    with out_path.open("w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"\n→ {out_path}")

    # TSV summary
    rows = []
    for tag, sw in results.items():
        ng, fset = tag.split("__")
        rows.append({
            "NG_column": ng, "feature_set": fset,
            "primary_tau": sw["tau"], "primary_dF1": sw["dF1"],
            "primary_P": sw["P"], "primary_R": sw["R"],
            "multi_mean_dF1": sw["multi_mean"], "multi_std_dF1": sw["multi_std"],
            "TP_kept": sw["tp_kept"], "FP_kept": sw["fp_kept"],
        })
    pd.DataFrame(rows).to_csv(STEP1 / "step1_lr_ablation_methylation.tsv", sep="\t", index=False)
    print(f"→ {STEP1 / 'step1_lr_ablation_methylation.tsv'}")

    # F11 figure — 4 bars with multi-seed std
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['DejaVu Sans', 'Droid Sans Fallback']
    fig, ax = plt.subplots(figsize=(11, 6.5))
    labels = ["V6\nFULL (10 feat)", "V6\nNoMeth (5 feat)",
              "baseline\nFULL (10 feat)", "baseline\nNoMeth (5 feat)"]
    means = [results['V6_off_NG__FULL_10feat']['multi_mean'],
             results['V6_off_NG__NOMETH_5feat']['multi_mean'],
             results['baseline_off_NG__FULL_10feat']['multi_mean'],
             results['baseline_off_NG__NOMETH_5feat']['multi_mean']]
    stds = [results['V6_off_NG__FULL_10feat']['multi_std'],
            results['V6_off_NG__NOMETH_5feat']['multi_std'],
            results['baseline_off_NG__FULL_10feat']['multi_std'],
            results['baseline_off_NG__NOMETH_5feat']['multi_std']]
    colors = ['#1f77b4', '#88c4f0', '#d62728', '#f0a0a0']
    bars = ax.bar(labels, means, yerr=stds, color=colors, capsize=6, alpha=0.85)
    ax.axhline(0.02236, color='#15803D', linestyle='--', linewidth=1.5,
               label='Cycle 1 V6 LR target +0.02236')
    ax.axhline(0.01, color='orange', linestyle=':', linewidth=1.2,
               label='H_C1_3 threshold +0.01')
    ax.axhline(0.00242, color='gray', linestyle=':', linewidth=1.2,
               label='v1.0 baseline +0.00242')
    ax.axhline(0, color='#222', linewidth=1)
    ax.set_ylabel("ΔF1 vs caller-alone (multi-seed mean)")
    ax.set_title(f"F11 — LR Ablation: drop methylation features\n"
                 f"V6 meth contrib = {v6_meth_contrib:+.5f} | baseline meth contrib = {base_meth_contrib:+.5f}\n"
                 f"VERDICT: {verdict[:80]}", fontsize=11)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, axis='y', alpha=0.3)
    for bar, mean in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width()/2, mean + 0.0008,
                f"{mean:+.5f}", ha='center', fontsize=9, fontweight='bold')
    plt.tight_layout()
    out_fig = FIG_DIR / "F11_lr_ablation_methylation.png"
    plt.savefig(out_fig, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"→ {out_fig}")


if __name__ == "__main__":
    main()
