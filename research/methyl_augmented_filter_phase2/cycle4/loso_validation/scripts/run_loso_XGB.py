#!/usr/bin/env python3
"""LOSO 5-fold CV using XGBoost (n_estimators=200, max_depth=6, lr=0.1)."""
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

try:
    import xgboost as xgb
except ImportError as e:
    print(f"ERROR: xgboost not installed. pip install xgboost. ({e})", file=sys.stderr)
    sys.exit(2)

from cross_sample_apply import (  # type: ignore  # noqa: E402
    CYCLE1_FEATURES,
    SAMPLE_CALLER_F1,
    SAMPLE_VARIANT_TPFP,
    HCC1395_TP_TOTAL,
    HCC1395_FP_TOTAL,
    HCC1395_FN_CALLER,
    PRIMARY_SEED,
    impute_with_medians,
    load_sample,
    reverse_solve_fn,
    sweep_tau,
)

import logging
LOGGER = logging.getLogger("loso_XGB")
LOGGER.setLevel(logging.INFO)
_h = logging.StreamHandler(sys.stdout)
_h.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S"))
LOGGER.handlers.clear()
LOGGER.addHandler(_h)

LOSO_DIR = REPO / "research/methyl_augmented_filter_phase2/cycle4/loso_validation"
OUT_TSV = LOSO_DIR / "data" / "loso_XGB_results.tsv"
OUT_FIG = LOSO_DIR / "figures" / "loso_XGB_5sample_dF1.png"
OUT_SUMMARY = LOSO_DIR / "intermediate" / "loso_XGB_summary.json"

ALL_SAMPLES = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009"]


def get_tpfp_baseline(sample: str, df: pd.DataFrame):
    caller_f1 = SAMPLE_CALLER_F1[sample]
    if sample == "HCC1395":
        return HCC1395_TP_TOTAL, HCC1395_FP_TOTAL, HCC1395_FN_CALLER, caller_f1
    var_tp, var_fp, var_fn = SAMPLE_VARIANT_TPFP[sample]
    if var_fn is None:
        var_fn = reverse_solve_fn(var_tp, var_fp, caller_f1)
    return var_tp, var_fp, var_fn, caller_f1


def loso_fold(test_sample, sample_dfs):
    train_samples = [s for s in ALL_SAMPLES if s != test_sample]
    LOGGER.info(f"  LOSO fold: test={test_sample}, train={train_samples}")
    train_combined = pd.concat([sample_dfs[s] for s in train_samples], ignore_index=True)
    train_imp, train_medians = impute_with_medians(train_combined, CYCLE1_FEATURES, medians=None)
    X_train = train_imp[CYCLE1_FEATURES].values.astype(float)
    y_train = train_imp["y"].values.astype(int)

    clf = xgb.XGBClassifier(
        n_estimators=200, max_depth=6, learning_rate=0.1,
        random_state=PRIMARY_SEED, n_jobs=-1,
        use_label_encoder=False, eval_metric="logloss", verbosity=0,
    )
    clf.fit(X_train, y_train)

    test_df = sample_dfs[test_sample]
    test_imp, _ = impute_with_medians(test_df, CYCLE1_FEATURES, medians=train_medians)
    X_test = test_imp[CYCLE1_FEATURES].values.astype(float)
    p_test = clf.predict_proba(X_test)[:, 1]
    y_test = test_imp["y"].values.astype(int)

    tp_total, fp_total, fn_caller, caller_f1 = get_tpfp_baseline(test_sample, test_df)
    best = sweep_tau(p_test, y_test, tp_total=tp_total, fp_total=fp_total,
                     fn_caller=fn_caller, caller_f1=caller_f1)
    return {
        "test_sample": test_sample,
        "best_tau": best["tau_star"],
        "delta_F1": best["delta_F1"],
        "F1_post": best["F1_post"],
        "caller_F1": caller_f1,
        "n_train": len(train_combined),
        "n_test": len(test_df),
    }


def main():
    LOGGER.info("Loading 5 sample master TSVs ...")
    sample_dfs = {s: load_sample(s, LOGGER) for s in ALL_SAMPLES}
    rows = [loso_fold(s, sample_dfs) for s in ALL_SAMPLES]
    df = pd.DataFrame(rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False)

    delta_f1 = df["delta_F1"].values
    nonzero = delta_f1[np.abs(delta_f1) > 1e-9]
    try:
        w_stat, w_p = wilcoxon(nonzero) if len(nonzero) >= 1 else (float("nan"),) * 2
    except ValueError:
        w_stat, w_p = float("nan"), float("nan")
    summary = {
        "algorithm": "XGBClassifier",
        "hyperparameters": {"n_estimators": 200, "max_depth": 6, "learning_rate": 0.1,
                             "random_state": PRIMARY_SEED, "n_jobs": -1},
        "per_sample_dF1": dict(zip(df["test_sample"], df["delta_F1"].astype(float))),
        "mean_dF1": float(np.mean(delta_f1)),
        "median_dF1": float(np.median(delta_f1)),
        "n_positive": int((delta_f1 > 1e-6).sum()),
        "n_negative": int((delta_f1 < -1e-6).sum()),
        "wilcoxon_p": float(w_p) if not np.isnan(w_p) else None,
    }
    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["#15803D" if v > 0.001 else ("#C2410C" if v < -0.001 else "#9CA3AF") for v in delta_f1]
    ax.bar(df["test_sample"], delta_f1, color=colors, edgecolor="black")
    ax.axhline(0, color="black", linewidth=0.7)
    ax.axhline(0.01, color="#15803D", linestyle="--", linewidth=1)
    ax.set_title("LOSO — XGBoost (n=200, depth=6, lr=0.1, seed=42)")
    ax.set_ylabel("ΔF1")
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
