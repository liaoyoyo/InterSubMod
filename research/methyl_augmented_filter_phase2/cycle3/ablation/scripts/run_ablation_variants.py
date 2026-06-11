#!/usr/bin/env python3
"""Cycle 3 Step 1.5 — Ablation variants (minimum viable evaluation).

Pre-registered hypotheses (per plan v2.1):
    H_A1   caller_af shrinkage → HCC1954 transfer ΔF1 improves > +0.30 from -0.377
    H_M1a  Drop ISM 5 methyl features → HCC1395 refit ΔF1 drop ≥ +0.003 (Cohen)
    H_M2   ISM block partial F-test p < 0.05 (full ablation only; not done here)
    H_M0   Residualized HPFineF AUC > 0.55 (full ablation only; not done here)

This script answers H_A1 + H_M1a with 4 variants × 5 samples × 2 modes = 40 combinations.

Variants:
    full          10 features (cycle 1 baseline)
    no-caller-af  9 features (drop caller_af)         → A1 isolated
    no-methyl     5 features (drop 5 ISM methyl)      → M1a isolated
    no-both       4 features (drop caller_af + methyl)→ A1+M1a joint

Modes:
    transfer  Apply cycle 1 per-fold coefs with TARGET coefs forced to 0 (coef shrinkage).
              For "no-methyl" mock the 5 methyl coef indices to 0; "no-caller-af"
              mock index 1 to 0; "no-both" mock all 6 to 0.
              τ = cycle 1 τ*=0.39 fixed (per H_A1 spec).
    refit     Re-train L2 LR on this feature subset, 5-fold StratifiedKFold OOF,
              swept τ ∈ [0.10, 0.95] step 0.01. Per-sample medians for imputation.

Reuses functions from cycle2/scripts/cross_sample_apply.py via sys.path import.

Outputs:
    cycle3/ablation/data/cycle3_step1_5_min_ablation.tsv  (40 rows × metrics)
    cycle3/ablation/figures/cycle3_step1_5_min_ablation.png (grouped bar 5×4 refit)
    cycle3/ablation/intermediate/cycle3_step1_5_summary.json (verdict + coef snapshots)
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

REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
CYCLE2_SCRIPTS = REPO_ROOT / "research/methyl_augmented_filter_phase2/cycle2/scripts"
sys.path.insert(0, str(CYCLE2_SCRIPTS))

# Reuse cycle 2 functions (verified line ranges in plan v2.1)
from cross_sample_apply import (  # type: ignore  # noqa: E402
    CYCLE1_FEATURES,
    SAMPLE_CALLER_F1,
    SAMPLE_VARIANT_TPFP,
    HCC1395_TP_TOTAL,
    HCC1395_FP_TOTAL,
    HCC1395_FN_CALLER,
    PRIMARY_SEED,
    TAU_GRID,
    add_cycle1_derived,
    compute_metrics,
    evaluate_fixed_tau,
    impute_with_medians,
    load_sample,
    refit_oof,
    reverse_solve_fn,
    sweep_tau,
    transfer_predict,
)

import logging
LOGGER = logging.getLogger("ablation")
LOGGER.setLevel(logging.INFO)
_h = logging.StreamHandler(sys.stdout)
_h.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S"))
LOGGER.handlers.clear()
LOGGER.addHandler(_h)

ABLATION_DIR = REPO_ROOT / "research/methyl_augmented_filter_phase2/cycle3/ablation"
OUT_TSV = ABLATION_DIR / "data" / "cycle3_step1_5_min_ablation.tsv"
OUT_FIG = ABLATION_DIR / "figures" / "cycle3_step1_5_min_ablation.png"
OUT_SUMMARY = ABLATION_DIR / "intermediate" / "cycle3_step1_5_summary.json"

CYCLE1_FILTER_PATH = (
    REPO_ROOT
    / "research/methyl_augmented_filter_phase2/cycle1/cycle1_track_a_filter.json"
)

ALL_SAMPLES = ["HCC1395", "H1437", "H2009", "HCC1954", "HCC1937"]

# Feature indices in CYCLE1_FEATURES (line 116-127 of cross_sample_apply.py):
#   0 V6_off_NG
#   1 caller_af                              ← A1 drop
#   2 loh_inner_flag
#   3 V6_off_meth_HPMergedDelta              ← M1a drop
#   4 V6_off_meth_HPFineF                    ← M1a drop
#   5 V6_off_meth_NME_imbalance              ← M1a drop
#   6 V6_off_meth_Epipoly_Delta              ← M1a drop
#   7 V6_off_meth_ClusterPermanovaF          ← M1a drop
#   8 chr8_flag
#   9 Coverage_Multiple_imp
CALLER_AF_IDX = 1
METHYL_IDXS = [3, 4, 5, 6, 7]

VARIANTS = {
    "full":         {"drop_idx": [],                            "label": "Full (10 features)"},
    "no-caller-af": {"drop_idx": [CALLER_AF_IDX],               "label": "No caller_af (9 features)"},
    "no-methyl":    {"drop_idx": METHYL_IDXS,                   "label": "No methyl (5 features)"},
    "no-both":      {"drop_idx": [CALLER_AF_IDX] + METHYL_IDXS, "label": "No both (4 features)"},
}


def get_feature_subset(drop_idx: list[int]) -> tuple[list[str], list[int]]:
    keep_idx = [i for i in range(len(CYCLE1_FEATURES)) if i not in drop_idx]
    keep_features = [CYCLE1_FEATURES[i] for i in keep_idx]
    return keep_features, keep_idx


def shrink_filter_coefs(filter_obj: dict, drop_idx: list[int]) -> dict:
    """Return a copy of cycle 1 filter with target coef indices forced to 0
    (transfer-mode coef shrinkage ablation)."""
    new_filter = json.loads(json.dumps(filter_obj))  # deep copy
    coefs = np.array(new_filter["per_fold_coefs"])  # shape 5×10
    for idx in drop_idx:
        coefs[:, idx] = 0.0
    new_filter["per_fold_coefs"] = coefs.tolist()
    return new_filter


def get_tpfp_totals(sample: str, df: pd.DataFrame) -> tuple[int, int, int, float]:
    caller_f1 = SAMPLE_CALLER_F1[sample]
    if sample == "HCC1395":
        return HCC1395_TP_TOTAL, HCC1395_FP_TOTAL, HCC1395_FN_CALLER, caller_f1

    n_tp = int((df["y"] == 1).sum())
    n_fp = int((df["y"] == 0).sum())
    var_tp, var_fp, var_fn = SAMPLE_VARIANT_TPFP[sample]
    if var_fn is None:
        var_fn = reverse_solve_fn(var_tp, var_fp, caller_f1)
    return n_tp, n_fp, var_fn, caller_f1


def run_one(sample: str, variant: str, df: pd.DataFrame, filter_obj: dict,
            cycle1_tau_star: float) -> list[dict]:
    """Run one (sample, variant) combination producing both transfer + refit rows."""
    drop_idx = VARIANTS[variant]["drop_idx"]
    keep_features, keep_idx = get_feature_subset(drop_idx)
    tp_total, fp_total, fn_caller, caller_f1 = get_tpfp_totals(sample, df)
    y = df["y"].values.astype(int)

    rows = []

    # ---- TRANSFER mode (coef shrinkage on full 10-feature filter) ----
    transfer_medians = {c: float(v) for c, v in filter_obj["feature_medians_for_imputation"].items()}
    sub_t, _ = impute_with_medians(df, CYCLE1_FEATURES, medians=transfer_medians)
    X_t = sub_t[CYCLE1_FEATURES].values.astype(float)
    shrunk_filter = shrink_filter_coefs(filter_obj, drop_idx)
    p_t = transfer_predict(X_t, shrunk_filter)
    transfer_fixed = evaluate_fixed_tau(
        p_t, y, tau=cycle1_tau_star,
        tp_total=tp_total, fp_total=fp_total,
        fn_caller=fn_caller, caller_f1=caller_f1,
    )
    rows.append({
        "sample": sample,
        "variant": variant,
        "mode": "transfer_fixed",
        "tau": cycle1_tau_star,
        "n_features_used": 10 - len(drop_idx),
        "drop_idx": ",".join(str(i) for i in drop_idx) if drop_idx else "",
        **{k: v for k, v in transfer_fixed.items() if k != "tau_star"},
        "caller_F1": caller_f1,
        "tp_total_used": tp_total,
        "fp_total_used": fp_total,
        "FN_caller_used": fn_caller,
    })

    # ---- REFIT mode (truly drop features and retrain LR) ----
    if len(keep_features) >= 2:
        sub_r, _ = impute_with_medians(df, keep_features, medians=None)
        X_r = sub_r[keep_features].values.astype(float)
        oof, coefs, _ = refit_oof(X_r, y, C=filter_obj["L2_C"], seed=PRIMARY_SEED)
        refit_best = sweep_tau(
            oof, y, tp_total=tp_total, fp_total=fp_total,
            fn_caller=fn_caller, caller_f1=caller_f1,
        )
        mean_coef = coefs.mean(axis=0).tolist()
        coef_summary = {f: round(mc, 4) for f, mc in zip(keep_features, mean_coef)}
    else:
        refit_best = {
            "delta_F1": np.nan, "delta_precision": np.nan, "delta_recall": np.nan,
            "precision": np.nan, "recall": np.nan, "F1_post": np.nan,
            "tp_kept": -1, "fp_kept": -1, "tp_removed": -1, "fp_removed": -1,
            "tau_star": np.nan,
        }
        coef_summary = {}

    rows.append({
        "sample": sample,
        "variant": variant,
        "mode": "refit",
        "tau": refit_best.get("tau_star", np.nan),
        "n_features_used": len(keep_features),
        "drop_idx": ",".join(str(i) for i in drop_idx) if drop_idx else "",
        **{k: v for k, v in refit_best.items() if k != "tau_star"},
        "caller_F1": caller_f1,
        "tp_total_used": tp_total,
        "fp_total_used": fp_total,
        "FN_caller_used": fn_caller,
        "refit_coef_json": json.dumps(coef_summary),
    })

    return rows


def main():
    filter_obj = json.loads(Path(CYCLE1_FILTER_PATH).read_text())
    cycle1_tau_star = float(filter_obj["tau_star"])
    LOGGER.info(f"Loaded cycle 1 filter; τ*={cycle1_tau_star}, L2_C={filter_obj['L2_C']}")

    all_rows = []
    sample_dfs = {}
    for sample in ALL_SAMPLES:
        LOGGER.info(f"--- loading {sample} ---")
        sample_dfs[sample] = load_sample(sample, LOGGER)

    for sample in ALL_SAMPLES:
        df = sample_dfs[sample]
        for variant in VARIANTS:
            LOGGER.info(f"  running {sample} × {variant}")
            rows = run_one(sample, variant, df, filter_obj, cycle1_tau_star)
            all_rows.extend(rows)

    out = pd.DataFrame(all_rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    LOGGER.info(f"Wrote {OUT_TSV} ({len(out)} rows)")

    # ---- Pre-reg verdicts ----
    # H_A1: caller_af shrinkage → HCC1954 transfer ΔF1 improves > +0.30 from -0.377
    full_hcc1954_transfer = float(
        out[(out["sample"] == "HCC1954") & (out["variant"] == "full")
            & (out["mode"] == "transfer_fixed")]["delta_F1"].iloc[0]
    )
    nocaller_hcc1954_transfer = float(
        out[(out["sample"] == "HCC1954") & (out["variant"] == "no-caller-af")
            & (out["mode"] == "transfer_fixed")]["delta_F1"].iloc[0]
    )
    h_a1_improvement = nocaller_hcc1954_transfer - full_hcc1954_transfer
    h_a1_verdict = "PASS" if h_a1_improvement > 0.30 else (
        "FAIL" if h_a1_improvement <= 0.10 else "MARGINAL")

    # H_M1a: HCC1395 refit no-methyl ΔF1 drop ≥ +0.003 from full
    full_hcc1395_refit = float(
        out[(out["sample"] == "HCC1395") & (out["variant"] == "full")
            & (out["mode"] == "refit")]["delta_F1"].iloc[0]
    )
    nomethyl_hcc1395_refit = float(
        out[(out["sample"] == "HCC1395") & (out["variant"] == "no-methyl")
            & (out["mode"] == "refit")]["delta_F1"].iloc[0]
    )
    h_m1a_drop = full_hcc1395_refit - nomethyl_hcc1395_refit
    h_m1a_verdict = "PASS" if h_m1a_drop >= 0.003 else (
        "FAIL" if h_m1a_drop < 0.001 else "MARGINAL")

    # Cross-sample sanity: 5-sample mean refit ΔF1 per variant
    refit_means = (
        out[out["mode"] == "refit"]
        .groupby("variant")["delta_F1"]
        .agg(["mean", "min", "max"])
        .to_dict("index")
    )

    summary = {
        "cycle_id": "20260519-0031-cycle3-caller-f1-headroom",
        "ablation_step": "1.5_min",
        "tau_grid": [float(TAU_GRID[0]), float(TAU_GRID[-1])],
        "variants_run": list(VARIANTS),
        "samples_run": ALL_SAMPLES,
        "cycle1_tau_star": cycle1_tau_star,
        "verdicts": {
            "H_A1_caller_af_shrinkage": {
                "full_HCC1954_transfer_dF1": full_hcc1954_transfer,
                "no_caller_HCC1954_transfer_dF1": nocaller_hcc1954_transfer,
                "improvement": h_a1_improvement,
                "threshold_pass": 0.30,
                "threshold_fail": 0.10,
                "verdict": h_a1_verdict,
            },
            "H_M1a_drop_methyl_block": {
                "full_HCC1395_refit_dF1": full_hcc1395_refit,
                "no_methyl_HCC1395_refit_dF1": nomethyl_hcc1395_refit,
                "drop": h_m1a_drop,
                "threshold_pass": 0.003,
                "threshold_fail": 0.001,
                "verdict": h_m1a_verdict,
            },
        },
        "refit_5sample_mean_dF1_per_variant": refit_means,
    }
    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    LOGGER.info(f"Wrote {OUT_SUMMARY}")

    # ---- Figure: grouped bar (refit mode) ----
    fig, ax = plt.subplots(figsize=(11, 5.5))
    variants_order = ["full", "no-caller-af", "no-methyl", "no-both"]
    colors = {"full": "#15803D", "no-caller-af": "#C2410C",
              "no-methyl": "#A16207", "no-both": "#6B7280"}
    width = 0.20
    x = np.arange(len(ALL_SAMPLES))
    for i, variant in enumerate(variants_order):
        vals = []
        for sample in ALL_SAMPLES:
            v = out[(out["sample"] == sample) & (out["variant"] == variant)
                   & (out["mode"] == "refit")]["delta_F1"].iloc[0]
            vals.append(float(v))
        offset = (i - 1.5) * width
        ax.bar(x + offset, vals, width, label=VARIANTS[variant]["label"],
               color=colors[variant], edgecolor="black", linewidth=0.4)
        for j, v in enumerate(vals):
            ax.text(x[j] + offset, v + 0.0008 if v >= 0 else v - 0.0008,
                    f"{v:+.4f}", ha="center",
                    va="bottom" if v >= 0 else "top", fontsize=7)

    ax.axhline(0, color="black", linewidth=0.5)
    ax.axhline(0.01, color="#15803D", linestyle="--", linewidth=0.8, label="Cohen +0.01")
    ax.set_xticks(x)
    ax.set_xticklabels(ALL_SAMPLES)
    ax.set_ylabel("ΔF1 (refit per-sample OOF)")
    ax.set_title("Cycle 3 Step 1.5 Min Ablation — refit ΔF1 per variant\n"
                 "(H_A1 + H_M1a feature ablation)", fontsize=11)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=120)
    plt.close(fig)
    LOGGER.info(f"Wrote {OUT_FIG}")

    # ---- Console verdict ----
    print()
    print("=== Cycle 3 Step 1.5 Min Ablation Verdicts ===")
    print(f"  H_A1 caller_af shrinkage:")
    print(f"    full HCC1954 transfer ΔF1 = {full_hcc1954_transfer:+.5f}")
    print(f"    no-caller-af HCC1954 transfer ΔF1 = {nocaller_hcc1954_transfer:+.5f}")
    print(f"    improvement = {h_a1_improvement:+.5f}  →  H_A1: {h_a1_verdict}")
    print(f"  H_M1a drop methylation block:")
    print(f"    full HCC1395 refit ΔF1 = {full_hcc1395_refit:+.5f}")
    print(f"    no-methyl HCC1395 refit ΔF1 = {nomethyl_hcc1395_refit:+.5f}")
    print(f"    drop = {h_m1a_drop:+.5f}  →  H_M1a: {h_m1a_verdict}")
    print()
    print("  Refit 5-sample mean ΔF1 per variant:")
    for v in variants_order:
        s = refit_means[v]
        print(f"    {v:14s}  mean={s['mean']:+.5f}  range=[{s['min']:+.5f}, {s['max']:+.5f}]")

    return 0


if __name__ == "__main__":
    sys.exit(main())
