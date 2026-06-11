#!/usr/bin/env python3
"""Apply Phase 2 Cycle 1 multi-axis LR to baseline BAM — does the +0.02236 mechanism work for baseline?

Approach:
  Cycle 1 uses 10 features including `V6_off_NG`. To test on baseline, we substitute
  `baseline_off_NG` (from step1_master_four_way) for `V6_off_NG`, keep the 9 other features
  (caller_af, loh_inner_flag, methylation×5, chr8_flag, Coverage_Multiple_imp).

  Methylation features (HPMergedDelta / HPFineF / NME_imbalance / Epipoly_Delta / ClusterPermanovaF)
  are inherited from V6 ISM output. Methylation analysis is largely HP-tag-INDEPENDENT
  (per-read base modification calls), so V6's methylation features are a reasonable
  approximation for baseline. The only BAM-specific substitution is NG.

  This isolates the question: "If we run the Cycle 1 LR with baseline_off_NG, does the
  ΔF1 = +0.02236 mechanism survive?"

  Procedure mirrors `final_filter_and_verdict.py`:
    - L2 C=1.0, lbfgs, 5-fold StratifiedKFold (seed=42 primary)
    - Out-of-fold probability + τ sweep [0.10, 0.95] step 0.01
    - Best ΔF1 = best F1 across τ - 0.7166 (caller-alone canonical)
    - Multi-seed (5 seeds) for stability check

Outputs:
  - apply_lr_baseline_results.json (full results)
  - apply_lr_baseline_tau_sweep.tsv (per τ TP/FP/dF1)
  - figures/baseline_vs_v6/F10_lr_apply_baseline.png (visualization)
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

# Constants from Cycle 1 (final_filter_and_verdict.py lines 43-44)
CALLER_F1_BASELINE = 0.7166
FN_CALLER = 19288
PRIMARY_SEED = 42
MULTI_SEEDS = [42, 7, 13, 2026, 1395]
TAU_GRID = np.round(np.arange(0.10, 0.96, 0.01), 2)
CHOSEN_C = 1.0


def compute_metrics(tp_kept, fp_kept):
    """Same formula as Cycle 1 final_filter_and_verdict.py compute_metrics()."""
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
    # Alias the NG column being tested to a generic "NG_feature" name
    df["NG_feature"] = df[ng_col]
    return df


def impute(df, features):
    sub = df.copy()
    medians = {}
    for c in features:
        if c not in sub.columns:
            sub[c] = np.nan
        med = float(sub[c].median()) if sub[c].notna().any() else 0.0
        if sub[c].isna().any():
            sub[c] = sub[c].fillna(med)
        medians[c] = med
    return sub, medians


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


def feat_names(ng_col):
    """LR features with NG_feature aliased."""
    return [
        "NG_feature",   # baseline_off_NG OR V6_off_NG depending on run
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


def run_lr_pipeline(df_merged, ng_col, label):
    """Single-seed primary + multi-seed ensemble for a given NG column."""
    print(f"\n=== Run LR with NG column = {ng_col} ({label}) ===", flush=True)
    df = prep(df_merged, ng_col)
    features = feat_names(ng_col)
    sub, medians = impute(df, features)
    X = sub[features].values.astype(float)
    y = sub["y"].values.astype(int)

    # Primary fit
    oof, coefs = fit_oof(X, y, CHOSEN_C, PRIMARY_SEED)
    sweep = sweep_tau(oof, y)
    sweep["features"] = features
    sweep["feature_coef_mean"] = list(zip(features, coefs.mean(0).tolist()))

    # Tau sweep table
    rows = []
    for tau in TAU_GRID:
        keep = oof >= tau
        tp_k = int(((y == 1) & keep).sum())
        fp_k = int(((y == 0) & keep).sum())
        dF1, P, R = compute_metrics(tp_k, fp_k)
        rows.append({"tau": float(tau), "tp_kept": tp_k, "fp_kept": fp_k,
                     "delta_F1": dF1, "P": P, "R": R})
    tau_df = pd.DataFrame(rows)

    # Multi-seed
    print(f"  Primary (seed=42): tau*={sweep['tau']:.2f}, ΔF1={sweep['dF1']:+.5f}", flush=True)
    multi = []
    for s in MULTI_SEEDS:
        oof_s, _ = fit_oof(X, y, CHOSEN_C, s)
        sw_s = sweep_tau(oof_s, y)
        multi.append({"seed": s, "tau": sw_s["tau"], "dF1": sw_s["dF1"]})
        print(f"  seed={s}: tau*={sw_s['tau']:.2f}, ΔF1={sw_s['dF1']:+.5f}", flush=True)
    multi_df = pd.DataFrame(multi)
    sweep["multi_mean_dF1"] = float(multi_df["dF1"].mean())
    sweep["multi_std_dF1"] = float(multi_df["dF1"].std(ddof=1))

    print(f"  Multi-seed mean ΔF1 = {sweep['multi_mean_dF1']:+.5f} ± {sweep['multi_std_dF1']:.5f}", flush=True)
    return sweep, tau_df, multi_df, oof


def main():
    print("Loading data ...", flush=True)
    four = pd.read_csv(MASTER_4WAY, sep="\t", low_memory=False)
    aug = pd.read_csv(MASTER_AUG, sep="\t", low_memory=False)
    print(f"  four_way: {four.shape} | aug: {aug.shape}", flush=True)

    # Join: master_augmented has V6 features + methylation; four_way has baseline_off_NG
    merged = aug.merge(four[["region_id", "baseline_off_NG"]], on="region_id", how="inner")
    print(f"  merged: {merged.shape} (inner join on region_id)", flush=True)

    # Run 1: V6 (regression test — should reproduce +0.02236)
    sweep_v6, tau_v6, multi_v6, _ = run_lr_pipeline(merged, "V6_off_NG", "V6 — regression test")

    # Run 2: baseline
    sweep_base, tau_base, multi_base, _ = run_lr_pipeline(merged, "baseline_off_NG", "baseline — V6→baseline NG substitution")

    # Save results
    out = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "ng_substitution_test": True,
        "v6_regression": {
            "primary_tau": sweep_v6["tau"],
            "primary_dF1": sweep_v6["dF1"],
            "primary_P": sweep_v6["P"],
            "primary_R": sweep_v6["R"],
            "primary_tp_kept": sweep_v6["tp_kept"],
            "primary_fp_kept": sweep_v6["fp_kept"],
            "multi_mean_dF1": sweep_v6["multi_mean_dF1"],
            "multi_std_dF1": sweep_v6["multi_std_dF1"],
            "expected_cycle1_dF1": 0.02236,
            "regression_pass": abs(sweep_v6["dF1"] - 0.02236) < 0.005,
        },
        "baseline_substitution": {
            "primary_tau": sweep_base["tau"],
            "primary_dF1": sweep_base["dF1"],
            "primary_P": sweep_base["P"],
            "primary_R": sweep_base["R"],
            "primary_tp_kept": sweep_base["tp_kept"],
            "primary_fp_kept": sweep_base["fp_kept"],
            "multi_mean_dF1": sweep_base["multi_mean_dF1"],
            "multi_std_dF1": sweep_base["multi_std_dF1"],
            "feature_coef_mean": sweep_base["feature_coef_mean"],
        },
        "delta_v6_minus_baseline_dF1": sweep_v6["dF1"] - sweep_base["dF1"],
        "verdict": (
            f"V6 NG: ΔF1={sweep_v6['dF1']:+.5f}; "
            f"baseline NG: ΔF1={sweep_base['dF1']:+.5f}; "
            f"V6 advantage at LR level: {sweep_v6['dF1'] - sweep_base['dF1']:+.5f}"
        ),
    }

    out_path = STEP1 / "apply_lr_baseline_results.json"
    with out_path.open("w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"\n→ {out_path}")

    # tau sweep TSV combined
    tau_v6["bam"] = "V6"
    tau_base["bam"] = "baseline"
    combined = pd.concat([tau_v6, tau_base], ignore_index=True)
    tau_path = STEP1 / "apply_lr_baseline_tau_sweep.tsv"
    combined.to_csv(tau_path, sep="\t", index=False)
    print(f"→ {tau_path}")

    # Plot F10 — tau sweep comparison
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['DejaVu Sans', 'Droid Sans Fallback']
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Left: τ sweep ΔF1
    ax = axes[0]
    ax.plot(tau_v6["tau"], tau_v6["delta_F1"], 'o-', color="#1f77b4", label=f"V6 LR (best ΔF1={sweep_v6['dF1']:+.4f} @ τ={sweep_v6['tau']:.2f})", markersize=4)
    ax.plot(tau_base["tau"], tau_base["delta_F1"], 's-', color="#d62728", label=f"baseline LR (best ΔF1={sweep_base['dF1']:+.4f} @ τ={sweep_base['tau']:.2f})", markersize=4)
    ax.axhline(0, color='#222', linewidth=1, label='caller-alone F1 = 0.7166')
    ax.axhline(0.02236, color='#15803D', linestyle='--', linewidth=1.5, label='Cycle 1 V6 LR target +0.02236')
    ax.set_xlabel("Threshold τ on LR P(TP)")
    ax.set_ylabel("ΔF1 vs caller-alone")
    ax.set_title("F10a — Multi-axis LR τ sweep (V6 NG vs baseline NG substitution)")
    ax.legend(loc='lower left', fontsize=9)
    ax.grid(True, alpha=0.3)

    # Right: multi-seed bar
    ax = axes[1]
    x_seeds = np.arange(len(MULTI_SEEDS))
    width = 0.4
    ax.bar(x_seeds - width/2, multi_v6["dF1"], width, color="#1f77b4", label="V6 NG", alpha=0.85)
    ax.bar(x_seeds + width/2, multi_base["dF1"], width, color="#d62728", label="baseline NG", alpha=0.85)
    ax.axhline(0.02236, color='#15803D', linestyle='--', linewidth=1.5, label='Cycle 1 V6 LR target')
    ax.set_xticks(x_seeds)
    ax.set_xticklabels([str(s) for s in MULTI_SEEDS])
    ax.set_xlabel("Random seed")
    ax.set_ylabel("ΔF1 (best τ)")
    ax.set_title(f"F10b — Multi-seed robustness\n"
                 f"V6 mean={sweep_v6['multi_mean_dF1']:+.5f}±{sweep_v6['multi_std_dF1']:.5f} | "
                 f"baseline mean={sweep_base['multi_mean_dF1']:+.5f}±{sweep_base['multi_std_dF1']:.5f}")
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, axis='y', alpha=0.3)
    for i, (v6, b) in enumerate(zip(multi_v6["dF1"], multi_base["dF1"])):
        ax.text(i - width/2, v6 + 0.0005, f"{v6:+.4f}", ha='center', fontsize=8)
        ax.text(i + width/2, b + 0.0005, f"{b:+.4f}", ha='center', fontsize=8)

    fig.suptitle("F10 — Phase 2 Cycle 1 LR applied to baseline NG (HCC1395, single-sample)\n"
                 "Test: does the +0.02236 mechanism survive V6_NG → baseline_NG substitution?",
                 fontsize=12)
    plt.tight_layout()
    out_fig = FIG_DIR / "F10_lr_apply_baseline.png"
    plt.savefig(out_fig, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"→ {out_fig}")

    # Verdict summary
    print(f"\n=== VERDICT ===")
    print(f"V6 LR (regression):     primary ΔF1 = {sweep_v6['dF1']:+.5f}, multi-seed mean = {sweep_v6['multi_mean_dF1']:+.5f}")
    print(f"baseline LR (NG subst): primary ΔF1 = {sweep_base['dF1']:+.5f}, multi-seed mean = {sweep_base['multi_mean_dF1']:+.5f}")
    print(f"V6 vs Cycle 1 published +0.02236: {'PASS ✓' if abs(sweep_v6['dF1'] - 0.02236) < 0.005 else 'DRIFT ⚠'} (drift = {sweep_v6['dF1'] - 0.02236:+.5f})")
    print(f"V6 - baseline LR ΔF1 gap: {sweep_v6['dF1'] - sweep_base['dF1']:+.5f}")
    if sweep_base["dF1"] > 0.01:
        print(f"  → baseline LR ALSO works! mechanism generalize beyond V6 NG")
    elif sweep_base["dF1"] > 0.00242:
        print(f"  → baseline LR weakly works (above v1.0 baseline but below H_C1_3)")
    else:
        print(f"  → baseline LR essentially fails — V6 NG is critical to the +0.02236 result")


if __name__ == "__main__":
    main()
