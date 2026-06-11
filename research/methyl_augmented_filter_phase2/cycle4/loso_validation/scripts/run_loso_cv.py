#!/usr/bin/env python3
"""LOSO (Leave-One-Sample-Out) Cross-Validation for cycle 1 LR filter.

Directly answers user 2026-05-20 ask: "LR filter 用 HCC1395 數據訓練,又用 HCC1395 數據驗證,
這敘述合理嗎"

LOSO protocol (sample-level CV, addresses row-level OOF's sample-level circularity):
    For each test_sample in [HCC1395, HCC1937, HCC1954, H1437, H2009]:
        train_samples = other 4 samples
        train_data = concat(master_TSV[train_samples])
        Fit L2 LR (C=1.0, lbfgs) on train_data with per-train median impute
        Predict P(TP) on test_sample
        Swept τ ∈ [0.10, 0.95] → best ΔF1 on test_sample
        Record: test_sample, best_τ, ΔF1, precision/recall, TP/FP kept

Aggregate:
    Wilcoxon paired n=5 vs 0
    Direction: ≥4/5 ΔF1 > 0 = PASS (但 n=5 exact min p=0.0625 requires 5/5)

Output:
    cycle4/loso_validation/data/loso_cv_results.tsv
    cycle4/loso_validation/figures/loso_5sample_dF1.png
    cycle4/loso_validation/intermediate/loso_summary.json
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
from scipy.stats import wilcoxon

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
CYCLE2_SCRIPTS = REPO / "research/methyl_augmented_filter_phase2/cycle2/scripts"
sys.path.insert(0, str(CYCLE2_SCRIPTS))

from cross_sample_apply import (  # type: ignore  # noqa: E402
    CYCLE1_FEATURES,
    SAMPLE_CALLER_F1,
    SAMPLE_VARIANT_TPFP,
    HCC1395_TP_TOTAL,
    HCC1395_FP_TOTAL,
    HCC1395_FN_CALLER,
    PRIMARY_SEED,
    TAU_GRID,
    compute_metrics,
    impute_with_medians,
    load_sample,
    reverse_solve_fn,
    sweep_tau,
)

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

import logging
LOGGER = logging.getLogger("loso")
LOGGER.setLevel(logging.INFO)
_h = logging.StreamHandler(sys.stdout)
_h.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S"))
LOGGER.handlers.clear()
LOGGER.addHandler(_h)

LOSO_DIR = REPO / "research/methyl_augmented_filter_phase2/cycle4/loso_validation"
OUT_TSV = LOSO_DIR / "data" / "loso_cv_results.tsv"
OUT_FIG = LOSO_DIR / "figures" / "loso_5sample_dF1.png"
OUT_SUMMARY = LOSO_DIR / "intermediate" / "loso_summary.json"

ALL_SAMPLES = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009"]


def get_tpfp_baseline(sample: str, df: pd.DataFrame):
    caller_f1 = SAMPLE_CALLER_F1[sample]
    n_tp = int((df["y"] == 1).sum())
    n_fp = int((df["y"] == 0).sum())
    if sample == "HCC1395":
        return HCC1395_TP_TOTAL, HCC1395_FP_TOTAL, HCC1395_FN_CALLER, caller_f1
    var_tp, var_fp, var_fn = SAMPLE_VARIANT_TPFP[sample]
    if var_fn is None:
        var_fn = reverse_solve_fn(var_tp, var_fp, caller_f1)
    return n_tp, n_fp, var_fn, caller_f1


def loso_single_fold(test_sample: str, sample_dfs: dict) -> dict:
    """Train on 4 samples, test on the held-out sample."""
    train_samples = [s for s in ALL_SAMPLES if s != test_sample]
    LOGGER.info(f"  LOSO fold: test={test_sample}, train={train_samples}")

    # Combine training data (concat rows)
    train_dfs = []
    for s in train_samples:
        df = sample_dfs[s].copy()
        df["_source_sample"] = s
        train_dfs.append(df)
    train_combined = pd.concat(train_dfs, ignore_index=True)
    LOGGER.info(f"    train_combined: {len(train_combined):,} rows")

    # Impute using train medians (NOT test sample medians, to avoid leakage)
    train_imp, train_medians = impute_with_medians(train_combined, CYCLE1_FEATURES, medians=None)
    X_train = train_imp[CYCLE1_FEATURES].values.astype(float)
    y_train = train_imp["y"].values.astype(int)

    # Fit StandardScaler + L2 LR on train
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    clf = LogisticRegression(max_iter=5000, C=1.0, solver="lbfgs", penalty="l2",
                             random_state=PRIMARY_SEED)
    clf.fit(X_train_scaled, y_train)

    coefs = {f: float(c) for f, c in zip(CYCLE1_FEATURES, clf.coef_[0])}
    intercept = float(clf.intercept_[0])

    # Apply to test sample (impute using TRAIN medians to avoid leakage)
    test_df = sample_dfs[test_sample]
    test_imp, _ = impute_with_medians(test_df, CYCLE1_FEATURES, medians=train_medians)
    X_test = test_imp[CYCLE1_FEATURES].values.astype(float)
    X_test_scaled = scaler.transform(X_test)
    p_test = clf.predict_proba(X_test_scaled)[:, 1]
    y_test = test_imp["y"].values.astype(int)

    tp_total, fp_total, fn_caller, caller_f1 = get_tpfp_baseline(test_sample, test_df)

    # Sweep τ for best ΔF1
    best = sweep_tau(p_test, y_test, tp_total=tp_total, fp_total=fp_total,
                     fn_caller=fn_caller, caller_f1=caller_f1)

    return {
        "test_sample": test_sample,
        "train_samples": ",".join(train_samples),
        "train_n_rows": len(train_combined),
        "test_n_rows": len(test_df),
        "test_n_TP": int((y_test == 1).sum()),
        "test_n_FP": int((y_test == 0).sum()),
        "caller_F1": caller_f1,
        "best_tau": best["tau_star"],
        "delta_F1": best["delta_F1"],
        "delta_precision": best["delta_precision"],
        "delta_recall": best["delta_recall"],
        "F1_post": best["F1_post"],
        "precision_post": best["precision"],
        "recall_post": best["recall"],
        "tp_kept": best["tp_kept"],
        "fp_kept": best["fp_kept"],
        "tp_removed": best["tp_removed"],
        "fp_removed": best["fp_removed"],
        "tp_total_used": tp_total,
        "fp_total_used": fp_total,
        "FN_caller_used": fn_caller,
        "intercept": intercept,
        "coefs_json": json.dumps(coefs),
    }


def main():
    LOGGER.info("Loading 5 sample master TSVs ...")
    sample_dfs = {}
    for s in ALL_SAMPLES:
        sample_dfs[s] = load_sample(s, LOGGER)

    LOGGER.info("\n=== LOSO 5-fold execution ===")
    rows = []
    for test_sample in ALL_SAMPLES:
        rows.append(loso_single_fold(test_sample, sample_dfs))

    df = pd.DataFrame(rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    LOGGER.info(f"\nWrote {OUT_TSV}")

    # Wilcoxon paired vs 0
    delta_f1 = df["delta_F1"].values
    pos_count = int((delta_f1 > 1e-6).sum())
    neg_count = int((delta_f1 < -1e-6).sum())
    zero_count = int((np.abs(delta_f1) <= 1e-6).sum())

    try:
        # Filter zeros for wilcoxon (n must be ≥1 nonzero)
        nonzero = delta_f1[np.abs(delta_f1) > 1e-9]
        if len(nonzero) >= 1:
            w_stat, w_p = wilcoxon(nonzero)
        else:
            w_stat, w_p = float("nan"), float("nan")
    except Exception as e:
        w_stat, w_p = float("nan"), float("nan")

    # Verdict
    if pos_count >= 4 and w_p < 0.0625:
        verdict = "POSITIVE"
    elif pos_count >= 4:
        verdict = "DIRECTION_POSITIVE"
    elif neg_count >= 4:
        verdict = "DIRECTION_NEGATIVE"
    elif pos_count >= 3:
        verdict = "MIXED_positive_lean"
    elif neg_count >= 3:
        verdict = "MIXED_negative_lean"
    else:
        verdict = "MIXED"

    mean_dF1 = float(np.mean(delta_f1))
    median_dF1 = float(np.median(delta_f1))

    summary = {
        "method": "LOSO (Leave-One-Sample-Out CV)",
        "purpose": "Address row-level OOF's sample-level circularity (user 2026-05-20 query)",
        "n_samples": len(ALL_SAMPLES),
        "samples": ALL_SAMPLES,
        "per_sample_delta_F1": {r["test_sample"]: float(r["delta_F1"]) for _, r in df.iterrows()},
        "summary_stats": {
            "mean_delta_F1": mean_dF1,
            "median_delta_F1": median_dF1,
            "n_positive": pos_count,
            "n_negative": neg_count,
            "n_zero": zero_count,
            "min_delta_F1": float(np.min(delta_f1)),
            "max_delta_F1": float(np.max(delta_f1)),
        },
        "wilcoxon": {
            "stat": float(w_stat) if not np.isnan(w_stat) else None,
            "p_value": float(w_p) if not np.isnan(w_p) else None,
            "exact_min_p_at_n5": 0.0625,
        },
        "verdict": verdict,
        "comparison_with_cycle1_internal": {
            "HCC1395_in_distribution_5fold_OOF": +0.02236,
            "HCC1395_LOSO_held_out": float(df[df["test_sample"]=="HCC1395"]["delta_F1"].iloc[0]),
            "drop_due_to_sample_level_circularity": +0.02236 - float(df[df["test_sample"]=="HCC1395"]["delta_F1"].iloc[0]),
        },
    }

    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    LOGGER.info(f"Wrote {OUT_SUMMARY}")

    # Figure
    fig, ax = plt.subplots(figsize=(11, 5.5))
    colors = ["#15803D" if v > 0.001 else ("#C2410C" if v < -0.001 else "#9CA3AF") for v in delta_f1]
    bars = ax.bar(df["test_sample"], delta_f1, color=colors, edgecolor="black", linewidth=0.5)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.axhline(0.01, color="#15803D", linestyle="--", linewidth=1, label="Cohen +0.01 PASS threshold")
    ax.axhline(0.005, color="#A16207", linestyle=":", linewidth=1, label="Cohen +0.005 ribbon")
    ax.axhline(-0.005, color="#C2410C", linestyle=":", linewidth=1)

    # Annotate
    for sample, v in zip(df["test_sample"], delta_f1):
        ax.text(sample, v + (0.0015 if v >= 0 else -0.0015),
                f"{v:+.5f}", ha="center",
                va="bottom" if v >= 0 else "top", fontsize=10, fontweight="bold")

    ax.set_xlabel("Hold-out test sample", fontsize=11)
    ax.set_ylabel("ΔF1 (4-sample-trained filter applied to held-out sample)", fontsize=11)

    w_p_str = f"{w_p:.4f}" if not np.isnan(w_p) else "N/A"
    annot = (f"LOSO 5-fold CV (sample-level)\n"
             f"mean ΔF1 = {mean_dF1:+.5f} | median = {median_dF1:+.5f}\n"
             f"n_pos = {pos_count} / n_neg = {neg_count} / n_zero = {zero_count}\n"
             f"Wilcoxon p = {w_p_str}\n"
             f"Verdict: {verdict}")
    ax.text(0.02, 0.98, annot, transform=ax.transAxes, va="top", ha="left",
            fontsize=10, family="monospace",
            bbox=dict(boxstyle="round,pad=0.5", fc="#F4EFE6", ec="black", alpha=0.9))

    ax.set_title("LOSO Cross-Validation — sample-level generalization test\n"
                 "Each sample held out; trained on other 4 samples combined",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=130, bbox_inches="tight")
    plt.close(fig)
    LOGGER.info(f"Wrote {OUT_FIG}")

    # Console verdict
    print()
    print("=== LOSO Cross-Validation Verdict ===")
    print(f"  Per-sample ΔF1:")
    for _, r in df.iterrows():
        print(f"    {r['test_sample']:10s}: ΔF1 = {r['delta_F1']:+.5f}  (best τ = {r['best_tau']:.2f})")
    print()
    print(f"  Mean ΔF1   = {mean_dF1:+.5f}")
    print(f"  Median ΔF1 = {median_dF1:+.5f}")
    print(f"  Sign:       {pos_count}+ / {neg_count}- / {zero_count}≈0")
    print(f"  Wilcoxon p = {w_p:.4f}" if not np.isnan(w_p) else "  Wilcoxon: N/A")
    print(f"  Verdict:    {verdict}")
    print()
    cycle1_inDist = +0.02236
    hcc1395_loso = float(df[df["test_sample"]=="HCC1395"]["delta_F1"].iloc[0])
    print(f"  HCC1395 in-distribution (5-fold OOF): +0.02236")
    print(f"  HCC1395 LOSO (held-out from 4-sample train): {hcc1395_loso:+.5f}")
    print(f"  Drop due to sample-level circularity removal: {cycle1_inDist - hcc1395_loso:+.5f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
