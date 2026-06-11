#!/usr/bin/env python3
"""H_NEW_4 LOSO: 9-feature LR (cycle 1 features minus caller_af) - sanity check.

Pre-reg: HCC1395 LOSO drop |ΔF1| < 0.005 confirms caller_af was the LOSO root cause.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
CYCLE2_SCRIPTS = REPO / "research/methyl_augmented_filter_phase2/cycle2/scripts"
sys.path.insert(0, str(CYCLE2_SCRIPTS))

from cross_sample_apply import (  # noqa: E402
    CYCLE1_FEATURES, SAMPLE_CALLER_F1, SAMPLE_VARIANT_TPFP,
    HCC1395_TP_TOTAL, HCC1395_FP_TOTAL, HCC1395_FN_CALLER,
    PRIMARY_SEED, compute_metrics, impute_with_medians, load_sample, reverse_solve_fn, sweep_tau,
)

import logging
LOG = logging.getLogger("h_new_4")
LOG.setLevel(logging.INFO)
_h = logging.StreamHandler(sys.stdout)
_h.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S"))
LOG.handlers.clear()
LOG.addHandler(_h)

LOSO_DIR = REPO / "research/methyl_augmented_filter_phase2/cycle4/loso_validation"
OUT_TSV = LOSO_DIR / "data" / "loso_hnew4_results.tsv"
OUT_FIG = LOSO_DIR / "figures" / "loso_hnew4_5sample.png"
OUT_SUMMARY = LOSO_DIR / "intermediate" / "loso_hnew4_summary.json"

# H_NEW_4: drop caller_af
HNEW4_FEATURES = [f for f in CYCLE1_FEATURES if f != "caller_af"]
ALL_SAMPLES = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009"]


def get_tpfp_baseline(sample, df):
    caller_f1 = SAMPLE_CALLER_F1[sample]
    if sample == "HCC1395":
        return HCC1395_TP_TOTAL, HCC1395_FP_TOTAL, HCC1395_FN_CALLER, caller_f1
    n_tp = int((df["y"] == 1).sum())
    n_fp = int((df["y"] == 0).sum())
    var_tp, var_fp, var_fn = SAMPLE_VARIANT_TPFP[sample]
    if var_fn is None:
        var_fn = reverse_solve_fn(var_tp, var_fp, caller_f1)
    return n_tp, n_fp, var_fn, caller_f1


def loso_fold(test_sample, sample_dfs, features):
    train_samples = [s for s in ALL_SAMPLES if s != test_sample]
    train_combined = pd.concat([sample_dfs[s] for s in train_samples], ignore_index=True)
    train_imp, train_medians = impute_with_medians(train_combined, features, medians=None)
    X_train = train_imp[features].values.astype(float)
    y_train = train_imp["y"].values.astype(int)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    clf = LogisticRegression(max_iter=5000, C=1.0, solver="lbfgs", penalty="l2", random_state=PRIMARY_SEED)
    clf.fit(X_train_scaled, y_train)

    test_df = sample_dfs[test_sample]
    test_imp, _ = impute_with_medians(test_df, features, medians=train_medians)
    X_test = scaler.transform(test_imp[features].values.astype(float))
    p_test = clf.predict_proba(X_test)[:, 1]
    y_test = test_imp["y"].values.astype(int)

    tp_total, fp_total, fn_caller, caller_f1 = get_tpfp_baseline(test_sample, test_df)
    best = sweep_tau(p_test, y_test, tp_total=tp_total, fp_total=fp_total,
                     fn_caller=fn_caller, caller_f1=caller_f1)

    coefs = {f: float(c) for f, c in zip(features, clf.coef_[0])}
    return {
        "test_sample": test_sample,
        "best_tau": best["tau_star"],
        "delta_F1": best["delta_F1"],
        "delta_precision": best["delta_precision"],
        "delta_recall": best["delta_recall"],
        "F1_post": best["F1_post"],
        "precision_post": best["precision"],
        "recall_post": best["recall"],
        "tp_kept": best["tp_kept"],
        "fp_kept": best["fp_kept"],
        "caller_F1": caller_f1,
        "tp_total_used": tp_total,
        "fp_total_used": fp_total,
        "coefs_json": json.dumps(coefs),
        "intercept": float(clf.intercept_[0]),
    }


def main():
    LOG.info(f"H_NEW_4 LOSO: 9-feature LR (drop caller_af) features={HNEW4_FEATURES}")
    sample_dfs = {s: load_sample(s, LOG) for s in ALL_SAMPLES}
    rows = [loso_fold(s, sample_dfs, HNEW4_FEATURES) for s in ALL_SAMPLES]
    df = pd.DataFrame(rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False)

    delta_f1 = df["delta_F1"].values
    hcc1395_loso = float(df[df["test_sample"]=="HCC1395"]["delta_F1"].iloc[0])
    n_pos = int((delta_f1 > 0.002).sum())
    mean_dF1 = float(np.mean(delta_f1))

    # Sanity: HCC1395 ΔF1 should still be ≈0
    if abs(hcc1395_loso) < 0.005:
        sanity = "CONFIRMED (caller_af was root cause of HCC1395 in-distribution overfit)"
    else:
        sanity = "VIOLATED (caller_af removal restored HCC1395 cross-sample signal — unexpected)"

    summary = {
        "method": "H_NEW_4 LOSO 9-feature (drop caller_af)",
        "features": HNEW4_FEATURES,
        "samples": ALL_SAMPLES,
        "per_sample_delta_F1": {r["test_sample"]: float(r["delta_F1"]) for _, r in df.iterrows()},
        "per_sample_best_tau": {r["test_sample"]: float(r["best_tau"]) for _, r in df.iterrows()},
        "summary": {
            "mean_delta_F1": mean_dF1,
            "n_pos_strict_above_0_002": n_pos,
            "hcc1395_LOSO_delta_F1": hcc1395_loso,
        },
        "sanity_verdict": sanity,
        "comparison": {
            "baseline_10_feature_LOSO_HCC1395": -0.00012,
            "hnew4_9_feature_LOSO_HCC1395": hcc1395_loso,
            "hcc1395_in_distribution_cycle1": +0.02236,
        },
        "pre_reg_prior_prob": 0.80,  # high confidence HCC1395 仍 ≈ 0
    }
    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    # Figure
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["#15803D" if v > 0.001 else ("#C2410C" if v < -0.001 else "#9CA3AF") for v in delta_f1]
    ax.bar(df["test_sample"], delta_f1, color=colors, edgecolor="black", linewidth=0.5)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.axhline(0.005, color="#A16207", linestyle=":", linewidth=1, label="Sanity boundary ±0.005")
    ax.axhline(-0.005, color="#A16207", linestyle=":", linewidth=1)
    for s, v in zip(df["test_sample"], delta_f1):
        ax.text(s, v + (0.0005 if v >= 0 else -0.0005), f"{v:+.5f}",
                ha="center", va="bottom" if v >= 0 else "top", fontsize=10, fontweight="bold")
    ax.set_ylabel("LOSO ΔF1 (held-out sample)")
    ax.set_title(f"H_NEW_4 LOSO — 9-feature LR (drop caller_af) - sanity check\n"
                 f"HCC1395 ΔF1={hcc1395_loso:+.5f} | mean={mean_dF1:+.5f} | {sanity[:50]}",
                 fontsize=11, fontweight="bold")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=130, bbox_inches="tight")
    plt.close(fig)

    print()
    print(f"=== H_NEW_4 LOSO Sanity ===")
    for _, r in df.iterrows():
        print(f"  {r['test_sample']:10s}: ΔF1 = {r['delta_F1']:+.5f}  (τ = {r['best_tau']:.2f})")
    print(f"  Mean ΔF1 = {mean_dF1:+.5f}")
    print(f"  HCC1395 in-distribution: +0.02236 | HCC1395 LOSO (10f baseline): -0.00012 | HCC1395 LOSO (9f drop caller_af): {hcc1395_loso:+.5f}")
    print(f"  Sanity: {sanity}")


if __name__ == "__main__":
    main()
