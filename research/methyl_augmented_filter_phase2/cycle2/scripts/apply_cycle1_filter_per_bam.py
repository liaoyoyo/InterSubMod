#!/usr/bin/env python3
"""
Phase 2 Cycle 2 Step C1 — Apply Cycle 1 filter to V3F / V5 / V6 per-BAM masters.

Two modes per BAM:
  1) transfer  — reuse cycle1_track_a_filter.json per-fold coefs + scalers + medians
                 (V6-trained → V3F/V5/V6 apply; expected ΔF1_V6 ≈ cycle 1 +0.02236
                 because V6 is the original training BAM).
  2) refit     — re-fit identical-architecture LR (L2 C=1.0, 5-fold StratKFold,
                 same features, median impute) on this BAM's data.

H_C1_6 verdict (cross-binary sanity):
    max_var = max(ΔF1_{V3F,V5,V6}) − min(ΔF1_{V3F,V5,V6})
    PASS         if max_var < 0.005
    BORDERLINE   if 0.005 ≤ max_var ≤ 0.010
    FAIL         if max_var > 0.010

The verdict is reported separately for the transfer-fit and the re-fit modes.

Outputs:
  cycle2/data/cycle2_c1_per_bam_delta_f1.tsv  (3 BAM × 2 modes × metrics)
  cycle2/figures/cycle2_c1_per_bam_delta_f1.png  (grouped bar chart)
  cycle2/intermediate/cycle2_c1_summary.json     (machine-readable verdicts)
"""

from __future__ import annotations

import argparse
import json
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

# CJK font fallback (cycle 1 convention).
matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

CALLER_F1_BASELINE = 0.7166
FN_CALLER = 19288
TP_TOTAL = 30490
FP_TOTAL = 4842
PRIMARY_SEED = 42
TAU_GRID = np.round(np.arange(0.10, 0.96, 0.01), 2)

BAM_TAGS = ["V3F", "V5", "V6"]

# Canonical feature names produced by build_per_bam_master.py.
CANONICAL_FEATURES = [
    "bam_off_NG",
    "caller_af",
    "loh_inner_flag",
    "bam_off_meth_HPMergedDelta",
    "bam_off_meth_HPFineF",
    "bam_off_meth_NME_imbalance",
    "bam_off_meth_Epipoly_Delta",
    "bam_off_meth_ClusterPermanovaF",
    "chr8_flag",
    "Coverage_Multiple_imp",
]

# Mapping from canonical → cycle1 filter feature names (same ordering!).
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
assert len(CANONICAL_FEATURES) == len(CYCLE1_FEATURES) == 10


# ---------- utilities -------------------------------------------------------


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def compute_metrics(tp_kept: int, fp_kept: int) -> dict:
    tp_removed = TP_TOTAL - tp_kept
    new_tp = tp_kept
    new_fp = fp_kept
    new_fn = FN_CALLER + tp_removed
    if new_tp + new_fp == 0 or new_tp + new_fn == 0:
        return {
            "delta_F1": 0.0, "precision": 0.0, "recall": 0.0,
            "F1_post": 0.0, "tp_kept": tp_kept, "fp_kept": fp_kept,
            "tp_removed": tp_removed, "fp_removed": FP_TOTAL - fp_kept,
        }
    precision = new_tp / (new_tp + new_fp)
    recall = new_tp / (new_tp + new_fn)
    if precision + recall == 0:
        return {
            "delta_F1": -CALLER_F1_BASELINE, "precision": precision,
            "recall": recall, "F1_post": 0.0, "tp_kept": tp_kept,
            "fp_kept": fp_kept, "tp_removed": tp_removed,
            "fp_removed": FP_TOTAL - fp_kept,
        }
    f1 = 2 * precision * recall / (precision + recall)
    return {
        "delta_F1": f1 - CALLER_F1_BASELINE,
        "precision": precision, "recall": recall,
        "F1_post": f1, "tp_kept": tp_kept, "fp_kept": fp_kept,
        "tp_removed": tp_removed, "fp_removed": FP_TOTAL - fp_kept,
    }


def sweep_tau(p: np.ndarray, y: np.ndarray) -> dict:
    best = None
    for tau in TAU_GRID:
        keep = p >= tau
        tp_kept = int(((y == 1) & keep).sum())
        fp_kept = int(((y == 0) & keep).sum())
        m = compute_metrics(tp_kept, fp_kept)
        m["tau_star"] = float(tau)
        if best is None or m["delta_F1"] > best["delta_F1"]:
            best = m
    return best


def evaluate_fixed_tau(p: np.ndarray, y: np.ndarray, tau: float) -> dict:
    keep = p >= tau
    tp_kept = int(((y == 1) & keep).sum())
    fp_kept = int(((y == 0) & keep).sum())
    m = compute_metrics(tp_kept, fp_kept)
    m["tau_star"] = float(tau)
    return m


# ---------- transfer-fit mode ----------------------------------------------


def transfer_predict(X: np.ndarray, filter_obj: dict) -> np.ndarray:
    """Apply cycle1 per-fold (coef + scaler) ensemble, then mean over folds."""
    per_fold_coefs = np.array(filter_obj["per_fold_coefs"])  # 5×10
    per_fold_intercepts = filter_obj["per_fold_intercepts"]
    per_fold_scalers = filter_obj["per_fold_scalers"]

    n = len(X)
    p_accum = np.zeros(n)
    for k in range(len(per_fold_coefs)):
        mean = np.array(per_fold_scalers[k]["mean"])
        scale = np.array(per_fold_scalers[k]["scale"])
        Xs = (X - mean) / scale
        z = Xs @ per_fold_coefs[k] + per_fold_intercepts[k]
        p_accum += 1.0 / (1.0 + np.exp(-z))
    return p_accum / len(per_fold_coefs)


# ---------- refit mode ------------------------------------------------------


def refit_oof(X: np.ndarray, y: np.ndarray, C: float = 1.0,
              seed: int = PRIMARY_SEED) -> tuple:
    """5-fold StratifiedKFold OOF logistic regression (L2)."""
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof = np.zeros(len(X))
    coefs = []
    intercepts = []
    for tr, va in skf.split(X, y):
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr])
        Xva = sc.transform(X[va])
        clf = LogisticRegression(
            max_iter=5000, C=C, solver="lbfgs", penalty="l2", random_state=seed,
        )
        clf.fit(Xtr, y[tr])
        oof[va] = clf.predict_proba(Xva)[:, 1]
        coefs.append(clf.coef_[0].tolist())
        intercepts.append(float(clf.intercept_[0]))
    return oof, np.array(coefs), intercepts


# ---------- per-BAM driver --------------------------------------------------


def impute_with_medians(df: pd.DataFrame, features: list[str],
                        medians: dict | None = None) -> tuple[pd.DataFrame, dict]:
    """If medians is provided, impute with those (transfer mode);
    otherwise compute medians from the data (refit mode)."""
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


def run_bam(bam: str, master_path: Path, filter_obj: dict) -> dict:
    log(f"  Loading per-BAM master: {master_path}")
    df = pd.read_csv(master_path, sep="\t", low_memory=False)
    if "y" not in df.columns:
        df["y"] = (df["label"] == "TP").astype(int)

    # --- Transfer mode ---
    canonical_to_cycle1 = dict(zip(CANONICAL_FEATURES, CYCLE1_FEATURES))
    transfer_medians = {
        canonical: float(filter_obj["feature_medians_for_imputation"]
                         [canonical_to_cycle1[canonical]])
        for canonical in CANONICAL_FEATURES
    }
    sub_t, _ = impute_with_medians(df, CANONICAL_FEATURES,
                                   medians=transfer_medians)
    X_t = sub_t[CANONICAL_FEATURES].values.astype(float)
    y = sub_t["y"].values.astype(int)
    p_transfer = transfer_predict(X_t, filter_obj)

    transfer_tau_fixed = evaluate_fixed_tau(p_transfer, y,
                                            tau=filter_obj["tau_star"])
    transfer_tau_sweep = sweep_tau(p_transfer, y)

    # --- Refit mode ---
    sub_r, refit_medians = impute_with_medians(df, CANONICAL_FEATURES,
                                               medians=None)
    X_r = sub_r[CANONICAL_FEATURES].values.astype(float)
    oof, coefs, intercepts = refit_oof(X_r, y, C=filter_obj["L2_C"],
                                       seed=PRIMARY_SEED)
    refit_best = sweep_tau(oof, y)
    mean_coef = coefs.mean(axis=0).tolist()
    std_coef = coefs.std(axis=0).tolist()
    feature_importance = sorted(
        [
            {"feature": fname, "mean_coef": mc, "std_coef": sd}
            for fname, mc, sd in zip(CANONICAL_FEATURES, mean_coef, std_coef)
        ],
        key=lambda r: -abs(r["mean_coef"]),
    )

    return {
        "bam": bam,
        "n_rows": int(len(df)),
        "n_TP": int((y == 1).sum()),
        "n_FP": int((y == 0).sum()),
        "transfer_fit": {
            "tau_used": float(filter_obj["tau_star"]),
            **transfer_tau_fixed,
            "tau_sweep_best": transfer_tau_sweep,
        },
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


# ---------- verdict + plotting ---------------------------------------------


def verdict_H_C1_6(values: dict[str, float]) -> dict:
    arr = np.array(list(values.values()))
    max_var = float(arr.max() - arr.min())
    if max_var < 0.005:
        v = "PASS"
    elif max_var <= 0.010:
        v = "BORDERLINE"
    else:
        v = "FAIL"
    return {"max_var": max_var, "verdict": v, "values": values}


def make_plot(records: dict, out_path: Path) -> None:
    transfer_vals = [records[b]["transfer_fit"]["delta_F1"] for b in BAM_TAGS]
    refit_vals = [records[b]["refit"]["delta_F1"] for b in BAM_TAGS]

    fig, ax = plt.subplots(figsize=(9, 6))
    x = np.arange(len(BAM_TAGS))
    w = 0.38
    bars_t = ax.bar(x - w / 2, transfer_vals, w,
                    label="Transfer (V6-trained coefs)", color="#1f77b4")
    bars_r = ax.bar(x + w / 2, refit_vals, w,
                    label="Re-fit (per-BAM)", color="#ff7f0e")

    ax.axhline(0.02236, color="gray", linestyle=":", linewidth=1,
               label="Cycle 1 V6 reference +0.02236")
    ax.axhline(0.01, color="green", linestyle=":", linewidth=1,
               label="H_C1_3 +0.01 (Cohen small)")
    ax.axhline(0, color="black", linewidth=0.6)

    for bars in (bars_t, bars_r):
        for bar in bars:
            h = bar.get_height()
            ax.annotate(f"{h:+.5f}",
                        xy=(bar.get_x() + bar.get_width() / 2, h),
                        xytext=(0, 4 if h >= 0 else -12), textcoords="offset points",
                        ha="center", va="bottom" if h >= 0 else "top", fontsize=8)

    transfer_var = max(transfer_vals) - min(transfer_vals)
    refit_var = max(refit_vals) - min(refit_vals)
    ax.set_xticks(x)
    ax.set_xticklabels(BAM_TAGS)
    ax.set_ylabel("ΔF1 (vs caller F1=0.7166)")
    ax.set_title(
        "H_C1_6 cross-binary sanity — Cycle 1 filter applied to V3F / V5 / V6\n"
        f"max_var transfer = {transfer_var:.5f}  |  max_var refit = {refit_var:.5f}  "
        "(PASS <0.005, BORDERLINE 0.005-0.010, FAIL >0.010)"
    )
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------- main ------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter-json",
        default="/big7_disk/liaoyoyo2001/InterSubMod/research/"
                "methyl_augmented_filter_phase2/cycle1/cycle1_track_a_filter.json",
    )
    parser.add_argument(
        "--master-dir",
        default="/big7_disk/liaoyoyo2001/InterSubMod/research/"
                "methyl_augmented_filter_phase2/cycle2/data",
    )
    parser.add_argument(
        "--output-dir",
        default="/big7_disk/liaoyoyo2001/InterSubMod/research/"
                "methyl_augmented_filter_phase2/cycle2",
    )
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    (outdir / "data").mkdir(parents=True, exist_ok=True)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    (outdir / "intermediate").mkdir(parents=True, exist_ok=True)

    log(f"Loading cycle 1 filter: {args.filter_json}")
    with open(args.filter_json) as f:
        filter_obj = json.load(f)
    log(f"  feature_set: {filter_obj['feature_set']}")
    log(f"  tau* = {filter_obj['tau_star']} | C = {filter_obj['L2_C']}")
    log(f"  expected V6 ΔF1 (training BAM) = {filter_obj['expected_delta_F1']:+.5f}")

    records = {}
    for bam in BAM_TAGS:
        log(f"\n=== BAM = {bam} ===")
        master_path = Path(args.master_dir) / f"per_bam_master_{bam}.tsv"
        rec = run_bam(bam, master_path, filter_obj)
        records[bam] = rec
        t = rec["transfer_fit"]
        r = rec["refit"]
        log(
            f"  Transfer (τ={t['tau_used']:.2f}): ΔF1={t['delta_F1']:+.5f} "
            f"prec={t['precision']:.4f} rec={t['recall']:.4f} "
            f"TP_rm={t['tp_removed']} FP_rm={t['fp_removed']}"
        )
        log(
            f"  Transfer best-tau sweep: ΔF1={t['tau_sweep_best']['delta_F1']:+.5f} "
            f"@ τ*={t['tau_sweep_best']['tau_star']:.2f}"
        )
        log(
            f"  Refit (τ*={r['tau_star']:.2f}): ΔF1={r['delta_F1']:+.5f} "
            f"prec={r['precision']:.4f} rec={r['recall']:.4f} "
            f"TP_rm={r['tp_removed']} FP_rm={r['fp_removed']}"
        )

    # ---------- TSV summary ----------
    rows = []
    for bam in BAM_TAGS:
        rec = records[bam]
        for mode in ("transfer_fit", "refit"):
            m = rec[mode]
            rows.append({
                "BAM": bam,
                "mode": mode,
                "tau": m.get("tau_used", m.get("tau_star")),
                "delta_F1": m["delta_F1"],
                "F1_post": m["F1_post"],
                "precision": m["precision"],
                "recall": m["recall"],
                "tp_kept": m["tp_kept"],
                "fp_kept": m["fp_kept"],
                "tp_removed": m["tp_removed"],
                "fp_removed": m["fp_removed"],
            })
    tsv_df = pd.DataFrame(rows)
    tsv_path = outdir / "data" / "cycle2_c1_per_bam_delta_f1.tsv"
    tsv_df.to_csv(tsv_path, sep="\t", index=False)
    log(f"\nWrote per-BAM metrics: {tsv_path}")

    # ---------- Verdict ----------
    transfer_values = {b: records[b]["transfer_fit"]["delta_F1"] for b in BAM_TAGS}
    refit_values = {b: records[b]["refit"]["delta_F1"] for b in BAM_TAGS}
    transfer_verdict = verdict_H_C1_6(transfer_values)
    refit_verdict = verdict_H_C1_6(refit_values)

    log("\n=== H_C1_6 verdict ===")
    log(
        f"  Transfer  max_var = {transfer_verdict['max_var']:.5f} -> "
        f"{transfer_verdict['verdict']}  (V3F={transfer_values['V3F']:+.5f}, "
        f"V5={transfer_values['V5']:+.5f}, V6={transfer_values['V6']:+.5f})"
    )
    log(
        f"  Re-fit    max_var = {refit_verdict['max_var']:.5f} -> "
        f"{refit_verdict['verdict']}  (V3F={refit_values['V3F']:+.5f}, "
        f"V5={refit_values['V5']:+.5f}, V6={refit_values['V6']:+.5f})"
    )

    # Sanity bonus: V6 refit ΔF1 should be ≈ +0.02236.
    v6_refit_drift = refit_values["V6"] - filter_obj["expected_delta_F1"]
    log(
        f"  Reproducibility sanity: V6 re-fit ΔF1 = {refit_values['V6']:+.5f}, "
        f"cycle 1 expected = {filter_obj['expected_delta_F1']:+.5f}, "
        f"drift = {v6_refit_drift:+.5f}"
    )

    # ---------- Figure ----------
    fig_path = outdir / "figures" / "cycle2_c1_per_bam_delta_f1.png"
    make_plot(records, fig_path)
    log(f"Wrote figure: {fig_path}")

    # ---------- JSON summary ----------
    summary = {
        "cycle": "Phase2_Cycle2_StepC1_H_C1_6",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "filter_source": args.filter_json,
        "cycle1_filter_tau_star": filter_obj["tau_star"],
        "cycle1_filter_L2_C": filter_obj["L2_C"],
        "cycle1_filter_expected_V6_delta_F1": filter_obj["expected_delta_F1"],
        "per_bam_records": records,
        "verdict_transfer": transfer_verdict,
        "verdict_refit": refit_verdict,
        "v6_refit_reproducibility_drift": v6_refit_drift,
    }
    json_path = outdir / "intermediate" / "cycle2_c1_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    log(f"Wrote summary JSON: {json_path}")


if __name__ == "__main__":
    main()
