#!/usr/bin/env python3
"""
Phase 2 Cycle 2 Step B3' — Cross-sample apply cycle 1 filter to 5 samples.

Samples (n=5): HCC1395 + H1437 + H2009 + HCC1954 + HCC1937

Two modes per sample (same as Step C1 architecture):
  1) transfer  — Apply cycle 1 V6-trained per-fold (coefs + scalers + medians) ensemble,
                 then mean-of-folds probability. Use τ*=0.39 from cycle 1.
  2) refit     — Re-fit identical-architecture LR (L2 C=1.0, 5-fold StratifiedKFold,
                 random_state=42, lbfgs, same 10 features, median impute per-sample).
                 τ* swept per-sample to maximize ΔF1 (best in 0.10–0.95 grid step 0.01).

Caller F1 baselines (sample-specific, from validated cross-sample report
20260216_跨樣本跨純度TP_FP_F1綜合分析報告_01.md "舊優化標準" table = ClairS-TO Pileup
without methylation filter):

  | Sample  | caller_F1 |
  | ------- | --------- |
  | HCC1395 | 0.7166    | (cycle 1 baseline, SEQC2 truth)
  | H1437   | 0.8670    | (Orthogonal truth, Pileup mode)
  | H2009   | 0.8863    | (Orthogonal truth, Pileup mode)
  | HCC1954 | 0.8385    | (Orthogonal truth, Pileup mode)
  | HCC1937 | 0.3692    | (Orthogonal truth, Pileup mode; BRCA1 outlier)

FN_caller is reverse-solved from F1 = 2TP / (2TP + FP + FN):
    FN_caller = 2*TP_variant / F1 − 2*TP_variant − FP_variant
where TP_variant and FP_variant are variant-level counts from filtered_snv_{tp,fp}.vcf.gz.

For HCC1395 we use the cycle 1 published TP/FP/FN = 30,490 / 4,842 / 19,288 (region-level
filter counts; cycle 1 reproduces ΔF1 = +0.02236 with these constants).

For each new sample we use variant-level TP/FP from VCF and reverse-solve FN.

CRITICAL: TP_TOTAL / FP_TOTAL used in compute_metrics() are region-level counts
(number of unique 10kb regions). When the filter drops a region, we count the
variants in that region. Since the master TSV is built region-level (joined to
VCF AF via region midpoint), tp_kept/fp_kept reflect region-level metrics, and
ΔF1 is computed on the same region-aggregate basis. This is consistent with
cycle 1's protocol.

Outputs:
  cycle2/data/cycle2_cross_sample_delta_f1.tsv   (5 samples × 2 modes × metrics)
  cycle2/figures/cycle2_b3_per_sample_delta_f1.png (5-sample bar chart)
  cycle2/intermediate/cross_sample_apply.log
  cycle2/intermediate/cycle2_b3_summary.json
"""
from __future__ import annotations

import argparse
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
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
CYCLE1_FILTER = (
    REPO_ROOT
    / "research/methyl_augmented_filter_phase2/cycle1/cycle1_track_a_filter.json"
)
CYCLE1_MASTER = (
    REPO_ROOT
    / "research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_master_augmented.tsv"
)
CYCLE2_DATA = REPO_ROOT / "research/methyl_augmented_filter_phase2/cycle2/data"
OUT_DIR = REPO_ROOT / "research/methyl_augmented_filter_phase2/cycle2"

PRIMARY_SEED = 42
TAU_GRID = np.round(np.arange(0.10, 0.96, 0.01), 2)

# Cycle 1 published constants (HCC1395 region-level filter counts).
HCC1395_TP_TOTAL = 30490
HCC1395_FP_TOTAL = 4842
HCC1395_FN_CALLER = 19288
HCC1395_CALLER_F1 = 0.7166

# Sample-specific caller F1 baselines from 20260216_跨樣本跨純度TP_FP_F1綜合分析報告_01.md
# "舊優化標準" (ClairS-TO Pileup baseline without methylation filter).
SAMPLE_CALLER_F1 = {
    "HCC1395": HCC1395_CALLER_F1,
    "H1437": 0.8670,
    "H2009": 0.8863,
    "HCC1954": 0.8385,
    "HCC1937": 0.3692,
}

# Variant-level TP/FP from filtered_snv_{tp,fp}.vcf.gz (manually verified at runtime).
# These also equal step4 per-sample master row counts for the 4 new samples.
# For HCC1395 we use cycle 1's region-level counts (30,490 / 4,842), since cycle 1 was
# designed at region-level and its filter was tuned on region-level totals.
SAMPLE_VARIANT_TPFP = {
    "HCC1395": (HCC1395_TP_TOTAL, HCC1395_FP_TOTAL, HCC1395_FN_CALLER),
    "H1437": (70191, 773, None),    # FN reverse-solved at runtime
    "H2009": (135359, 1342, None),
    "HCC1954": (19449, 687, None),
    "HCC1937": (13910, 2697, None),
}

NEW_SAMPLES = ["H1437", "H2009", "HCC1954", "HCC1937"]
ALL_SAMPLES = ["HCC1395"] + NEW_SAMPLES

# Cycle 1 feature names (must match cycle1_track_a_filter.json `feature_set`).
CYCLE1_FEATURES = [
    "V6_off_NG",
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


def setup_logger() -> logging.Logger:
    log_path = OUT_DIR / "intermediate" / "cross_sample_apply.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("cross_sample_apply")
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_path, mode="w")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.handlers.clear()
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def reverse_solve_fn(tp: int, fp: int, f1: float) -> int:
    """FN = 2*TP / F1 - 2*TP - FP, rounded to nearest int."""
    fn_float = 2.0 * tp / f1 - 2.0 * tp - fp
    return int(round(fn_float))


def compute_metrics(tp_kept: int, fp_kept: int, tp_total: int, fp_total: int,
                    fn_caller: int, caller_f1: float) -> dict:
    tp_removed = tp_total - tp_kept
    new_tp = tp_kept
    new_fp = fp_kept
    new_fn = fn_caller + tp_removed
    if new_tp + new_fp == 0 or new_tp + new_fn == 0:
        return {
            "delta_F1": -caller_f1, "delta_precision": -1.0, "delta_recall": -1.0,
            "precision": 0.0, "recall": 0.0, "F1_post": 0.0,
            "tp_kept": tp_kept, "fp_kept": fp_kept,
            "tp_removed": tp_removed, "fp_removed": fp_total - fp_kept,
        }
    precision = new_tp / (new_tp + new_fp)
    recall = new_tp / (new_tp + new_fn)
    if precision + recall == 0:
        f1 = 0.0
    else:
        f1 = 2 * precision * recall / (precision + recall)

    # Caller baseline precision/recall for delta comparison.
    caller_precision = tp_total / (tp_total + fp_total) if (tp_total + fp_total) > 0 else 0.0
    caller_recall = tp_total / (tp_total + fn_caller) if (tp_total + fn_caller) > 0 else 0.0

    return {
        "delta_F1": f1 - caller_f1,
        "delta_precision": precision - caller_precision,
        "delta_recall": recall - caller_recall,
        "precision": precision,
        "recall": recall,
        "F1_post": f1,
        "tp_kept": tp_kept,
        "fp_kept": fp_kept,
        "tp_removed": tp_removed,
        "fp_removed": fp_total - fp_kept,
        "caller_precision": caller_precision,
        "caller_recall": caller_recall,
    }


def sweep_tau(p: np.ndarray, y: np.ndarray, tp_total: int, fp_total: int,
              fn_caller: int, caller_f1: float) -> dict:
    best = None
    for tau in TAU_GRID:
        keep = p >= tau
        tp_kept = int(((y == 1) & keep).sum())
        fp_kept = int(((y == 0) & keep).sum())
        m = compute_metrics(tp_kept, fp_kept, tp_total, fp_total, fn_caller, caller_f1)
        m["tau_star"] = float(tau)
        if best is None or m["delta_F1"] > best["delta_F1"]:
            best = m
    return best


def evaluate_fixed_tau(p: np.ndarray, y: np.ndarray, tau: float, tp_total: int,
                       fp_total: int, fn_caller: int, caller_f1: float) -> dict:
    keep = p >= tau
    tp_kept = int(((y == 1) & keep).sum())
    fp_kept = int(((y == 0) & keep).sum())
    m = compute_metrics(tp_kept, fp_kept, tp_total, fp_total, fn_caller, caller_f1)
    m["tau_star"] = float(tau)
    return m


def transfer_predict(X: np.ndarray, filter_obj: dict) -> np.ndarray:
    """Cycle 1 per-fold (coef + scaler + intercept) ensemble; mean prob."""
    per_fold_coefs = np.array(filter_obj["per_fold_coefs"])  # 5×10
    per_fold_intercepts = filter_obj["per_fold_intercepts"]
    per_fold_scalers = filter_obj["per_fold_scalers"]
    p_accum = np.zeros(len(X))
    for k in range(len(per_fold_coefs)):
        mean = np.array(per_fold_scalers[k]["mean"])
        scale = np.array(per_fold_scalers[k]["scale"])
        Xs = (X - mean) / scale
        z = Xs @ per_fold_coefs[k] + per_fold_intercepts[k]
        p_accum += 1.0 / (1.0 + np.exp(-z))
    return p_accum / len(per_fold_coefs)


def refit_oof(X: np.ndarray, y: np.ndarray, C: float = 1.0,
              seed: int = PRIMARY_SEED) -> tuple:
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof = np.zeros(len(X))
    coefs = []
    intercepts = []
    for tr, va in skf.split(X, y):
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr])
        Xva = sc.transform(X[va])
        clf = LogisticRegression(
            max_iter=5000, C=C, solver="lbfgs", penalty="l2", random_state=seed
        )
        clf.fit(Xtr, y[tr])
        oof[va] = clf.predict_proba(Xva)[:, 1]
        coefs.append(clf.coef_[0].tolist())
        intercepts.append(float(clf.intercept_[0]))
    return oof, np.array(coefs), intercepts


def impute_with_medians(df: pd.DataFrame, features: list[str],
                        medians: dict | None = None) -> tuple[pd.DataFrame, dict]:
    sub = df.copy()
    used_medians = {}
    for c in features:
        if medians is not None and c in medians:
            med = float(medians[c])
        else:
            med = float(sub[c].median()) if sub[c].notna().any() else 0.0
        if sub[c].isna().any():
            sub[c] = sub[c].fillna(med)
        used_medians[c] = med
    return sub, used_medians


def add_cycle1_derived(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure cycle 1 derived columns exist: loh_inner_flag, Coverage_Multiple_imp,
    chr8_flag, y. Idempotent — only adds when missing."""
    df = df.copy()
    if "loh_inner_flag" not in df.columns and "loh_side" in df.columns:
        df["loh_inner_flag"] = (df["loh_side"] == "Inner").astype(int)
    if "Coverage_Multiple_imp" not in df.columns and "Coverage_Multiple" in df.columns:
        cov_median = df["Coverage_Multiple"].median()
        df["Coverage_Multiple_imp"] = df["Coverage_Multiple"].fillna(cov_median)
    if "chr8_flag" not in df.columns:
        df["chr8_flag"] = (df["chr"] == "chr8").astype(int)
    if "y" not in df.columns:
        df["y"] = (df["label"] == "TP").astype(int)
    return df


def load_sample(sample: str, logger: logging.Logger) -> pd.DataFrame:
    """Load the per-sample augmented master TSV (HCC1395 uses cycle 1 v1.0 master)."""
    if sample == "HCC1395":
        path = CYCLE1_MASTER
    else:
        path = CYCLE2_DATA / f"{sample}_master_augmented.tsv"
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = add_cycle1_derived(df)
    logger.info(
        f"  loaded {sample}: rows={len(df):,} cols={len(df.columns)}"
        f" TP={int((df['y']==1).sum()):,} FP={int((df['y']==0).sum()):,}"
    )
    return df


def run_sample(sample: str, filter_obj: dict, logger: logging.Logger) -> dict:
    logger.info(f"\n=== {sample} ===")
    df = load_sample(sample, logger)
    n_TP_region = int((df["y"] == 1).sum())
    n_FP_region = int((df["y"] == 0).sum())

    # Variant-level TP/FP/FN baselines.
    var_tp, var_fp, var_fn = SAMPLE_VARIANT_TPFP[sample]
    caller_f1 = SAMPLE_CALLER_F1[sample]
    if var_fn is None:
        var_fn = reverse_solve_fn(var_tp, var_fp, caller_f1)
    logger.info(
        f"  variant-level: TP={var_tp:,} FP={var_fp:,} FN={var_fn:,} "
        f"caller_F1={caller_f1:.4f}"
    )

    # CRITICAL: We use REGION-level counts as TP_TOTAL/FP_TOTAL inside compute_metrics
    # (because the LR predicts at region level) but the F1 baseline reference for ΔF1
    # is the variant-level caller F1. For consistency with cycle 1 (which used 30490 /
    # 4842 / 19288 for HCC1395 = its region counts), we use region-level TP_TOTAL /
    # FP_TOTAL but variant-level caller_F1 + FN_caller. This matches cycle 1's
    # CALLER_F1_BASELINE = 0.7166 + region-level FN_CALLER = 19288 protocol.
    #
    # For the 4 new samples we reverse-solve FN from variant-level F1 and use
    # region-level TP/FP as the filter target.
    tp_total = n_TP_region
    fp_total = n_FP_region

    if sample == "HCC1395":
        # Cycle 1 protocol: cycle1 used variant-level TP=30490, FP=4842, FN=19288.
        # The master TSV has 30,476 TP rows + 4,820 FP rows (3 sample-level), but
        # cycle 1 hard-coded 30,490/4,842. We follow cycle 1 hard-coded constants
        # to bit-exactly reproduce cycle 1 ΔF1.
        tp_total = HCC1395_TP_TOTAL
        fp_total = HCC1395_FP_TOTAL
        fn_caller = HCC1395_FN_CALLER
    else:
        fn_caller = var_fn

    # ----- Transfer mode -----
    transfer_medians = {c: float(v) for c, v in filter_obj["feature_medians_for_imputation"].items()}
    sub_t, _ = impute_with_medians(df, CYCLE1_FEATURES, medians=transfer_medians)
    X_t = sub_t[CYCLE1_FEATURES].values.astype(float)
    y = sub_t["y"].values.astype(int)
    p_transfer = transfer_predict(X_t, filter_obj)
    transfer_fixed = evaluate_fixed_tau(
        p_transfer, y, tau=filter_obj["tau_star"],
        tp_total=tp_total, fp_total=fp_total,
        fn_caller=fn_caller, caller_f1=caller_f1,
    )
    transfer_swept = sweep_tau(
        p_transfer, y, tp_total=tp_total, fp_total=fp_total,
        fn_caller=fn_caller, caller_f1=caller_f1,
    )

    # ----- Refit mode -----
    sub_r, refit_medians = impute_with_medians(df, CYCLE1_FEATURES, medians=None)
    X_r = sub_r[CYCLE1_FEATURES].values.astype(float)
    oof, coefs, intercepts = refit_oof(X_r, y, C=filter_obj["L2_C"], seed=PRIMARY_SEED)
    refit_best = sweep_tau(
        oof, y, tp_total=tp_total, fp_total=fp_total,
        fn_caller=fn_caller, caller_f1=caller_f1,
    )
    mean_coef = coefs.mean(axis=0).tolist()
    std_coef = coefs.std(axis=0).tolist()
    feature_importance = sorted(
        [
            {"feature": f, "mean_coef": mc, "std_coef": sd}
            for f, mc, sd in zip(CYCLE1_FEATURES, mean_coef, std_coef)
        ],
        key=lambda r: -abs(r["mean_coef"]),
    )

    record = {
        "sample": sample,
        "caller_F1": caller_f1,
        "n_TP_region": n_TP_region,
        "n_FP_region": n_FP_region,
        "tp_total_used": tp_total,
        "fp_total_used": fp_total,
        "variant_TP": var_tp,
        "variant_FP": var_fp,
        "FN_caller_used": fn_caller,
        "transfer_fit_fixed_tau": {
            "tau_used": float(filter_obj["tau_star"]),
            **transfer_fixed,
        },
        "transfer_fit_sweep": transfer_swept,
        "refit": {
            "C": filter_obj["L2_C"],
            "n_splits": 5,
            "seed": PRIMARY_SEED,
            "medians_used": refit_medians,
            **refit_best,
            "per_fold_coefs": coefs.tolist(),
            "per_fold_intercepts": intercepts,
            "feature_importance": feature_importance,
        },
    }
    t = transfer_fixed
    r = refit_best
    logger.info(
        f"  Transfer (τ={t['tau_star']:.2f}): ΔF1={t['delta_F1']:+.5f} "
        f"F1_post={t['F1_post']:.4f} prec={t['precision']:.4f} rec={t['recall']:.4f} "
        f"TP_rm={t['tp_removed']} FP_rm={t['fp_removed']}"
    )
    logger.info(
        f"  Refit (τ*={r['tau_star']:.2f}): ΔF1={r['delta_F1']:+.5f} "
        f"F1_post={r['F1_post']:.4f} prec={r['precision']:.4f} rec={r['recall']:.4f} "
        f"TP_rm={r['tp_removed']} FP_rm={r['fp_removed']}"
    )
    return record


def make_plot(records: list[dict], out_path: Path) -> None:
    samples = [r["sample"] for r in records]
    transfer = [r["transfer_fit_fixed_tau"]["delta_F1"] for r in records]
    refit = [r["refit"]["delta_F1"] for r in records]
    transfer_sweep = [r["transfer_fit_sweep"]["delta_F1"] for r in records]

    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(len(samples))
    w = 0.25
    bars_t = ax.bar(x - w, transfer, w, label="Transfer (V6-trained, τ=0.39)", color="#1f77b4")
    bars_ts = ax.bar(x, transfer_sweep, w, label="Transfer (V6-trained, swept τ)", color="#9ecae1")
    bars_r = ax.bar(x + w, refit, w, label="Re-fit (per-sample, swept τ)", color="#ff7f0e")

    ax.axhline(0.02236, color="gray", linestyle=":", linewidth=1,
               label="HCC1395 cycle 1 reference +0.02236")
    ax.axhline(0.01, color="green", linestyle=":", linewidth=1,
               label="H_C1_3 +0.01 (Cohen small)")
    ax.axhline(0, color="black", linewidth=0.6)

    for bars in (bars_t, bars_ts, bars_r):
        for bar in bars:
            h = bar.get_height()
            ax.annotate(
                f"{h:+.4f}",
                xy=(bar.get_x() + bar.get_width() / 2, h),
                xytext=(0, 4 if h >= 0 else -10),
                textcoords="offset points",
                ha="center", va="bottom" if h >= 0 else "top", fontsize=8,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(samples)
    ax.set_ylabel("ΔF1 (vs caller F1)")
    ax.set_title(
        "Step B3' Cross-sample apply of cycle 1 filter (n=5)\n"
        "Transfer (fixed τ=0.39 / swept τ) and Re-fit (per-sample 5-fold CV)"
    )
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--filter-json", default=str(CYCLE1_FILTER))
    args = parser.parse_args()

    logger = setup_logger()
    t0 = time.time()
    logger.info("=" * 78)
    logger.info("Step B3' — Cross-sample apply cycle 1 filter (n=5)")
    logger.info("=" * 78)

    logger.info(f"Loading cycle 1 filter: {args.filter_json}")
    with open(args.filter_json) as f:
        filter_obj = json.load(f)
    logger.info(f"  feature_set: {filter_obj['feature_set']}")
    logger.info(f"  tau*={filter_obj['tau_star']} C={filter_obj['L2_C']}")
    logger.info(f"  expected V6 ΔF1 (HCC1395 training) = {filter_obj['expected_delta_F1']:+.5f}")

    records = []
    for sample in ALL_SAMPLES:
        rec = run_sample(sample, filter_obj, logger)
        records.append(rec)

    # ----- TSV -----
    rows = []
    for rec in records:
        for mode_key, mode_label in [
            ("transfer_fit_fixed_tau", "transfer_fixed"),
            ("transfer_fit_sweep", "transfer_sweep"),
            ("refit", "refit"),
        ]:
            m = rec[mode_key]
            rows.append({
                "sample": rec["sample"],
                "mode": mode_label,
                "tau": m.get("tau_used", m.get("tau_star")),
                "delta_F1": m["delta_F1"],
                "delta_precision": m["delta_precision"],
                "delta_recall": m["delta_recall"],
                "F1_post": m["F1_post"],
                "precision_post": m["precision"],
                "recall_post": m["recall"],
                "caller_precision": m["caller_precision"],
                "caller_recall": m["caller_recall"],
                "caller_F1": rec["caller_F1"],
                "tp_kept": m["tp_kept"],
                "fp_kept": m["fp_kept"],
                "tp_removed": m["tp_removed"],
                "fp_removed": m["fp_removed"],
                "tp_total_used": rec["tp_total_used"],
                "fp_total_used": rec["fp_total_used"],
                "FN_caller_used": rec["FN_caller_used"],
            })
    tsv_df = pd.DataFrame(rows)
    tsv_path = OUT_DIR / "data" / "cycle2_cross_sample_delta_f1.tsv"
    tsv_df.to_csv(tsv_path, sep="\t", index=False)
    logger.info(f"\nWrote per-sample TSV: {tsv_path}")

    # ----- Figure -----
    fig_path = OUT_DIR / "figures" / "cycle2_b3_per_sample_delta_f1.png"
    make_plot(records, fig_path)
    logger.info(f"Wrote figure: {fig_path}")

    # ----- JSON summary -----
    summary = {
        "cycle": "Phase2_Cycle2_StepB3prime",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "filter_source": args.filter_json,
        "cycle1_filter_tau_star": filter_obj["tau_star"],
        "cycle1_filter_L2_C": filter_obj["L2_C"],
        "cycle1_filter_expected_V6_delta_F1": filter_obj["expected_delta_F1"],
        "per_sample_records": records,
    }
    json_path = OUT_DIR / "intermediate" / "cycle2_b3_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    logger.info(f"Wrote JSON summary: {json_path}")

    elapsed = time.time() - t0
    logger.info("=" * 78)
    logger.info(f"Step B3' COMPLETE in {elapsed:.1f} s ({elapsed/60:.2f} min)")
    logger.info("=" * 78)


if __name__ == "__main__":
    main()
