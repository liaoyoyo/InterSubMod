#!/usr/bin/env python3
"""Phase A4 — Multi-algorithm LOSO benchmark (LR vs DT vs RF vs XGBoost).

Question
========
PI 必問: "為何不用 DT / RF / XGBoost?" — 用戶要求驗證 cycle 4 LOSO sample-level
circularity 失敗（LR mean ΔF1 = -0.00004, 4/5 negative）**非 LR-specific** —
非線性 ML 在 sample-level held-out 也同樣失敗才能下結論 "sample-level circularity
exhausts all 10 features regardless of model class"。

Design
======
- 5-sample LOSO (held-out × 5): HCC1395 / HCC1937 / HCC1954 / H1437 / H2009
- 10 features (cycle 1 canonical, drop NumReads VIF=217)
- 4 algorithms:
    LR  — LogisticRegression L2 C=1.0 lbfgs (cycle 1 baseline)
    DT  — DecisionTreeClassifier max_depth=5, min_samples_leaf=100
    RF  — RandomForestClassifier n_estimators=200, max_depth=8, min_samples_leaf=50
    XGB — XGBClassifier n_estimators=200, max_depth=6, learning_rate=0.1
- 5 random seeds: 42, 7, 100, 999, 2026 → mean ± std (single-seed artifact 排除)
- Total: 4 algo × 5 seed × 5 sample = 100 folds
- Sanity: train F1 vs test F1 gap (overfit symptom)

Inputs (reused from cycle 2 / cycle 4)
======================================
- cross_sample_apply.py — CYCLE1_FEATURES, SAMPLE_CALLER_F1, SAMPLE_VARIANT_TPFP,
  HCC1395_TP/FP/FN constants, sweep_tau, impute_with_medians, load_sample,
  reverse_solve_fn

Outputs
=======
- phase2_completeness_audit/A4_LR_DT_RF_XGBoost_LOSO_5sample.tsv
- phase2_completeness_audit/figures/A4_4algorithm_LOSO_dF1.png
- phase2_completeness_audit/figures/A4_train_test_gap.png

Wall-time
=========
~30-90 min (XGBoost / RF dominate; XGB n_jobs=-1)
"""
from __future__ import annotations

import json
import logging
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
    TAU_GRID,
    compute_metrics,
    impute_with_medians,
    load_sample,
    reverse_solve_fn,
    sweep_tau,
)

from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler

# XGBoost (optional)
try:
    import xgboost as xgb
    HAVE_XGB = True
except ImportError:
    HAVE_XGB = False

# ─── Paths ───
AUDIT_DIR = REPO / "research/methyl_augmented_filter_phase2/phase2_completeness_audit"
OUT_TSV = AUDIT_DIR / "A4_LR_DT_RF_XGBoost_LOSO_5sample.tsv"
OUT_FIG_DF1 = AUDIT_DIR / "figures" / "A4_4algorithm_LOSO_dF1.png"
OUT_FIG_GAP = AUDIT_DIR / "figures" / "A4_train_test_gap.png"
OUT_SUMMARY = AUDIT_DIR / "data" / "A4_summary.json"

# ─── Settings ───
ALL_SAMPLES = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009"]
SEEDS = [42, 7, 100, 999, 2026]
matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

# ─── Logger ───
LOGGER = logging.getLogger("A4_multi_algo")
LOGGER.setLevel(logging.INFO)
_h = logging.StreamHandler(sys.stdout)
_h.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S"))
LOGGER.handlers.clear()
LOGGER.addHandler(_h)


# ─── Algorithm factory ───
def make_model(algo: str, seed: int):
    """Construct model from algo name + seed."""
    if algo == "LR":
        return LogisticRegression(
            max_iter=5000, C=1.0, solver="lbfgs", penalty="l2",
            random_state=seed,
        )
    if algo == "DT":
        return DecisionTreeClassifier(
            max_depth=5, min_samples_leaf=100, random_state=seed,
        )
    if algo == "RF":
        return RandomForestClassifier(
            n_estimators=200, max_depth=8, min_samples_leaf=50,
            n_jobs=-1, random_state=seed,
        )
    if algo == "XGB":
        if not HAVE_XGB:
            raise RuntimeError("xgboost not installed")
        return xgb.XGBClassifier(
            n_estimators=200, max_depth=6, learning_rate=0.1,
            random_state=seed, n_jobs=-1,
            use_label_encoder=False, eval_metric="logloss",
            verbosity=0,
        )
    raise ValueError(f"unknown algo: {algo}")


def needs_scaling(algo: str) -> bool:
    """LR needs StandardScaler; tree-based don't (DT/RF/XGB are scale-invariant)."""
    return algo == "LR"


# ─── TP/FP/FN baseline for held-out sample ───
def get_tpfp_baseline(sample: str, df: pd.DataFrame):
    caller_f1 = SAMPLE_CALLER_F1[sample]
    if sample == "HCC1395":
        return HCC1395_TP_TOTAL, HCC1395_FP_TOTAL, HCC1395_FN_CALLER, caller_f1
    var_tp, var_fp, var_fn = SAMPLE_VARIANT_TPFP[sample]
    if var_fn is None:
        var_fn = reverse_solve_fn(var_tp, var_fp, caller_f1)
    return var_tp, var_fp, var_fn, caller_f1


# ─── Single LOSO fold for one (algo, seed, test_sample) ───
def loso_single_fold(algo: str, seed: int, test_sample: str,
                     sample_dfs: dict) -> dict:
    train_samples = [s for s in ALL_SAMPLES if s != test_sample]

    # Combine training data
    train_dfs = []
    for s in train_samples:
        train_dfs.append(sample_dfs[s])
    train_combined = pd.concat(train_dfs, ignore_index=True)

    # Impute using TRAIN medians (no test leakage)
    train_imp, train_medians = impute_with_medians(
        train_combined, CYCLE1_FEATURES, medians=None,
    )
    X_train = train_imp[CYCLE1_FEATURES].values.astype(float)
    y_train = train_imp["y"].values.astype(int)

    # Apply same medians to test
    test_df = sample_dfs[test_sample]
    test_imp, _ = impute_with_medians(
        test_df, CYCLE1_FEATURES, medians=train_medians,
    )
    X_test = test_imp[CYCLE1_FEATURES].values.astype(float)
    y_test = test_imp["y"].values.astype(int)

    # Scale only for LR
    if needs_scaling(algo):
        sc = StandardScaler()
        X_train_use = sc.fit_transform(X_train)
        X_test_use = sc.transform(X_test)
    else:
        X_train_use = X_train
        X_test_use = X_test

    # Fit + predict
    clf = make_model(algo, seed)
    t0 = time.time()
    clf.fit(X_train_use, y_train)
    fit_sec = time.time() - t0

    p_train = clf.predict_proba(X_train_use)[:, 1]
    p_test = clf.predict_proba(X_test_use)[:, 1]

    # Sample-specific baselines
    tp_total, fp_total, fn_caller, caller_f1 = get_tpfp_baseline(test_sample, test_df)

    # Sweep τ on test → best ΔF1, best τ, F1_post
    best_test = sweep_tau(
        p_test, y_test, tp_total=tp_total, fp_total=fp_total,
        fn_caller=fn_caller, caller_f1=caller_f1,
    )

    # Train F1 at SAME τ as test (for gap calculation, using training-population
    # baselines: aggregate train TP/FP, and reverse-solve FN by weighted caller F1).
    # Approximation: aggregate train TP/FP/FN by sum across train samples.
    tp_train_total = 0
    fp_train_total = 0
    fn_train_total = 0
    f1_weighted_num = 0.0
    f1_weighted_den = 0
    for s in train_samples:
        if s == "HCC1395":
            tp_s, fp_s, fn_s = HCC1395_TP_TOTAL, HCC1395_FP_TOTAL, HCC1395_FN_CALLER
            cf1_s = SAMPLE_CALLER_F1[s]
        else:
            v_tp, v_fp, v_fn = SAMPLE_VARIANT_TPFP[s]
            cf1_s = SAMPLE_CALLER_F1[s]
            if v_fn is None:
                v_fn = reverse_solve_fn(v_tp, v_fp, cf1_s)
            tp_s, fp_s, fn_s = v_tp, v_fp, v_fn
        tp_train_total += tp_s
        fp_train_total += fp_s
        fn_train_total += fn_s
        f1_weighted_num += cf1_s * (tp_s + fp_s + fn_s)
        f1_weighted_den += (tp_s + fp_s + fn_s)
    caller_f1_train = f1_weighted_num / f1_weighted_den if f1_weighted_den else 0.0

    tau_best = best_test["tau_star"]
    keep_train = p_train >= tau_best
    tp_kept_tr = int(((y_train == 1) & keep_train).sum())
    fp_kept_tr = int(((y_train == 0) & keep_train).sum())
    # Use train aggregates (sample TP/FP totals) as baseline for train F1
    # Note: y_train counts == sample row counts; use sample-level TP/FP totals
    # since train_combined rows are region-level (consistent with cycle 1).
    n_tp_train = int((y_train == 1).sum())
    n_fp_train = int((y_train == 0).sum())
    # train_combined is region-level; use these as totals (FN approximation by aggregate)
    train_metrics = compute_metrics(
        tp_kept_tr, fp_kept_tr, n_tp_train, n_fp_train, fn_train_total, caller_f1_train,
    )
    train_F1 = train_metrics["F1_post"]
    test_F1 = best_test["F1_post"]
    train_test_gap = train_F1 - test_F1

    return {
        "algorithm": algo,
        "seed": seed,
        "test_sample": test_sample,
        "best_tau": best_test["tau_star"],
        "delta_F1": best_test["delta_F1"],
        "train_F1": train_F1,
        "test_F1": test_F1,
        "train_test_gap": train_test_gap,
        "n_train": len(train_combined),
        "n_test": len(test_df),
        "caller_F1": caller_f1,
        "fit_sec": fit_sec,
        "test_n_TP": int((y_test == 1).sum()),
        "test_n_FP": int((y_test == 0).sum()),
    }


# ─── Main ───
def main():
    LOGGER.info("=== Phase A4 — multi-algorithm LOSO ===")
    LOGGER.info(f"  XGBoost available: {HAVE_XGB}")
    algorithms = ["LR", "DT", "RF"]
    if HAVE_XGB:
        algorithms.append("XGB")
    else:
        LOGGER.warning("  XGBoost NOT installed — skipping XGB")
    LOGGER.info(f"  algorithms: {algorithms}")
    LOGGER.info(f"  seeds: {SEEDS}")
    LOGGER.info(f"  samples: {ALL_SAMPLES}")
    LOGGER.info(f"  total folds: {len(algorithms) * len(SEEDS) * len(ALL_SAMPLES)}")

    # Load samples once (reused across all folds)
    LOGGER.info("\nLoading 5 sample master TSVs ...")
    sample_dfs = {}
    for s in ALL_SAMPLES:
        sample_dfs[s] = load_sample(s, LOGGER)

    # Run grid: algo × seed × sample
    rows = []
    total_start = time.time()
    fold_idx = 0
    total_folds = len(algorithms) * len(SEEDS) * len(ALL_SAMPLES)
    for algo in algorithms:
        for seed in SEEDS:
            for test_sample in ALL_SAMPLES:
                fold_idx += 1
                t0 = time.time()
                LOGGER.info(f"[{fold_idx}/{total_folds}] algo={algo} seed={seed} "
                            f"test={test_sample}")
                try:
                    r = loso_single_fold(algo, seed, test_sample, sample_dfs)
                    rows.append(r)
                    LOGGER.info(
                        f"    ΔF1 = {r['delta_F1']:+.5f}  "
                        f"train_F1 = {r['train_F1']:.4f}  test_F1 = {r['test_F1']:.4f}  "
                        f"gap = {r['train_test_gap']:+.4f}  "
                        f"({time.time()-t0:.1f}s)"
                    )
                except Exception as e:
                    LOGGER.error(f"    FAILED: {e}")
                    rows.append({
                        "algorithm": algo, "seed": seed, "test_sample": test_sample,
                        "best_tau": np.nan, "delta_F1": np.nan,
                        "train_F1": np.nan, "test_F1": np.nan, "train_test_gap": np.nan,
                        "n_train": np.nan, "n_test": np.nan, "caller_F1": np.nan,
                        "fit_sec": np.nan, "test_n_TP": np.nan, "test_n_FP": np.nan,
                        "error": str(e),
                    })

    LOGGER.info(f"\nTotal wall-time: {(time.time()-total_start)/60:.1f} min")

    # ─── Save TSV ───
    df = pd.DataFrame(rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    LOGGER.info(f"Wrote {OUT_TSV}")

    # ─── Summary stats: mean ± std across 5 seeds per (algo, sample) ───
    df_clean = df.dropna(subset=["delta_F1"])
    grp = df_clean.groupby(["algorithm", "test_sample"]).agg(
        delta_F1_mean=("delta_F1", "mean"),
        delta_F1_std=("delta_F1", "std"),
        gap_mean=("train_test_gap", "mean"),
        gap_std=("train_test_gap", "std"),
        n_seed=("seed", "nunique"),
    ).reset_index()

    # Cross-sample mean per algorithm (overall ΔF1)
    overall = df_clean.groupby("algorithm").agg(
        overall_mean_dF1=("delta_F1", "mean"),
        overall_std_dF1=("delta_F1", "std"),
        overall_gap_mean=("train_test_gap", "mean"),
        n_positive=("delta_F1", lambda v: int((v > 1e-6).sum())),
        n_negative=("delta_F1", lambda v: int((v < -1e-6).sum())),
        n_total=("delta_F1", "count"),
    ).reset_index()

    summary = {
        "design": {
            "algorithms": algorithms,
            "seeds": SEEDS,
            "samples": ALL_SAMPLES,
            "features": CYCLE1_FEATURES,
            "n_features": len(CYCLE1_FEATURES),
            "total_folds": int(df_clean.shape[0]),
        },
        "per_algo_per_sample": grp.to_dict(orient="records"),
        "per_algo_overall": overall.to_dict(orient="records"),
        "verdict_test": "sample-level circularity model-specific? "
                        "Yes if non-LR ΔF1 > +0.01 and LR < 0; No if all algo ≤ 0",
    }
    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    LOGGER.info(f"Wrote {OUT_SUMMARY}")

    # ─── Figure 1: ΔF1 grouped bar (algo × sample, error bar = std across seeds) ───
    fig, ax = plt.subplots(figsize=(13, 6.5))
    n_algo = len(algorithms)
    width = 0.8 / n_algo
    x_base = np.arange(len(ALL_SAMPLES))
    colors = {"LR": "#2563EB", "DT": "#16A34A", "RF": "#D97706", "XGB": "#DC2626"}
    for i, algo in enumerate(algorithms):
        sub = grp[grp["algorithm"] == algo].set_index("test_sample").reindex(ALL_SAMPLES)
        means = sub["delta_F1_mean"].values
        stds = sub["delta_F1_std"].fillna(0).values
        x = x_base + (i - (n_algo - 1) / 2) * width
        ax.bar(x, means, width, yerr=stds, capsize=4,
               color=colors.get(algo, "#888"), label=algo,
               edgecolor="black", linewidth=0.5)
        # Annotate mean
        for xi, m in zip(x, means):
            if not np.isnan(m):
                va = "bottom" if m >= 0 else "top"
                off = 0.0008 if m >= 0 else -0.0008
                ax.text(xi, m + off, f"{m:+.4f}", ha="center", va=va,
                        fontsize=7, rotation=90)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.axhline(0.01, color="#15803D", linestyle="--", linewidth=1,
               label="Cohen +0.01 PASS")
    ax.axhline(0.005, color="#A16207", linestyle=":", linewidth=1,
               label="Cohen +0.005 ribbon")
    ax.axhline(-0.005, color="#C2410C", linestyle=":", linewidth=1)
    ax.set_xticks(x_base)
    ax.set_xticklabels(ALL_SAMPLES)
    ax.set_xlabel("Hold-out test sample", fontsize=11)
    ax.set_ylabel("ΔF1 (mean ± std across 5 seeds)", fontsize=11)
    ax.set_title("Phase A4 — Multi-algorithm LOSO benchmark\n"
                 "(5 algorithms × 5 seeds × 5 samples = held-out generalization)",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper left", fontsize=9, ncol=2)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    OUT_FIG_DF1.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG_DF1, dpi=130, bbox_inches="tight")
    plt.close(fig)
    LOGGER.info(f"Wrote {OUT_FIG_DF1}")

    # ─── Figure 2: Train-test gap (overfit symptom) ───
    fig, ax = plt.subplots(figsize=(11, 5.5))
    gap_by_algo = []
    labels = []
    for algo in algorithms:
        sub = df_clean[df_clean["algorithm"] == algo]
        gap_by_algo.append(sub["train_test_gap"].values)
        labels.append(algo)
    bp = ax.boxplot(gap_by_algo, labels=labels, patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.2))
    for patch, algo in zip(bp["boxes"], algorithms):
        patch.set_facecolor(colors.get(algo, "#888"))
        patch.set_alpha(0.6)
    # Overlay individual seeds as dots
    for i, algo in enumerate(algorithms):
        sub = df_clean[df_clean["algorithm"] == algo]
        x_jit = np.random.RandomState(42).normal(i + 1, 0.04, len(sub))
        ax.scatter(x_jit, sub["train_test_gap"].values, s=18, alpha=0.5,
                   color="black", zorder=3)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.axhline(0.05, color="#C2410C", linestyle="--", linewidth=1,
               label="Overfit warning (gap > 0.05)")
    ax.set_ylabel("train F1 − test F1 (positive = overfitting on training samples)",
                  fontsize=10)
    ax.set_title("Phase A4 — Train-test F1 gap (overfit symptom)\n"
                 "Box = 25/50/75% across (5 seeds × 5 samples) = 25 folds per algo",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_FIG_GAP, dpi=130, bbox_inches="tight")
    plt.close(fig)
    LOGGER.info(f"Wrote {OUT_FIG_GAP}")

    # ─── Console verdict ───
    print()
    print("=== Phase A4 — Multi-algorithm LOSO Verdict ===")
    print(f"{'algorithm':<8s} {'overall_mean':>12s} {'overall_std':>12s} "
          f"{'gap_mean':>10s} {'n_pos':>5s}/{'n_neg':>5s}/{'n_tot':>5s}")
    for _, row in overall.iterrows():
        print(f"{row['algorithm']:<8s} {row['overall_mean_dF1']:+12.5f} "
              f"{row['overall_std_dF1']:>12.5f} {row['overall_gap_mean']:>+10.4f} "
              f"{int(row['n_positive']):>5d}/{int(row['n_negative']):>5d}/"
              f"{int(row['n_total']):>5d}")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
