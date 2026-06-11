"""Step 5d — Robustness Audit (C.1 / C.4 / C.7)

Performs three audits on Step 5 methylation-augmented FP filter pilot:

  C.1 Multi-seed CV variance
      Re-run Step 3 ΔF1 pipeline with 5 different KFold seeds (42 / 123 / 456 / 789 / 1000).
      For each seed:
        - 5-fold CV out-of-fold P(TP) per kept cell (Model B on V6_off)
        - aggregated multi-cell pool sweep across tau in [0.10, 0.95]
        - record Strategy 1 max ΔF1, winning cell, winning τ*
      Verdict: std > 0.001 RED, 0.0005-0.001 YELLOW, ≤0.0005 GREEN.

  C.4 H1 LRT unique cell count
      From step1_lrt_per_cell.tsv (138 rows = 3 BAMs × 2 flags × 23 cells),
      filter q<0.05 rows, count UNIQUE cell_id (collapsing BAM × flag),
      and report consistency (in how many BAM × flag combos each unique cell hit).
      Verdict: ≥4 unique = GREEN, ≤2 = YELLOW.

  C.7 NaN rate per cell × covariate
      For (4 step3 aggregated winning cells) + (step1 LRT q<0.05 unique cells),
      compute Model B effective n / n_total / NaN rate across 5 methyl covariates
      per BAM × flag.
      Verdict: ≥0.7 = GREEN, 0.5-0.7 = YELLOW, <0.5 = RED.

Outputs:
  step5d_robustness_audit.md
  intermediate/step5d_c1_multiseed_variance.tsv
  intermediate/step5d_c4_unique_cells.tsv
  intermediate/step5d_c7_nan_rate.tsv
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold

from _common_step1 import (
    MASTER_AUG,
    METHYL_COVARIATES,
    STEP5_DIR,
    annotate_v6off_cell,
    log_msg,
)
from augmented_lr import assemble_design, fit_lr
from filter_sweep import (
    compute_cov_proxy,
    find_fp_rich_cells,
    select_cell_rows,
    FP_DENSITY_FILE,
    PRIMARY_BAM,
    PRIMARY_FLAG,
)

warnings.filterwarnings("ignore")

# Caller baseline constants (mirrors delta_f1.py)
TP_CALLER = 30490
FP_CALLER = 4842
FN_CALLER = 19288
F1_CALLER = 0.7166

TAU_GRID = np.round(np.arange(0.10, 0.9501, 0.01), 2)
CV_FOLDS = 5
SEEDS = [42, 123, 456, 789, 1000]

MIN_N_EFF_B = 100
MIN_N_TP = 5

OUT_C1 = STEP5_DIR / "intermediate" / "step5d_c1_multiseed_variance.tsv"
OUT_C4 = STEP5_DIR / "intermediate" / "step5d_c4_unique_cells.tsv"
OUT_C7 = STEP5_DIR / "intermediate" / "step5d_c7_nan_rate.tsv"
OUT_MD = STEP5_DIR / "step5d_robustness_audit.md"
LOG_PATH = STEP5_DIR / "intermediate" / "step5d_log.txt"


def s5d_log(msg: str) -> None:
    import datetime as _dt
    ts = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(line + "\n")


def f1_from_counts(tp: float, fp: float, fn: float) -> float:
    if tp + fp == 0 or tp + fn == 0:
        return 0.0
    p = tp / (tp + fp)
    r = tp / (tp + fn)
    if p + r == 0:
        return 0.0
    return 2 * p * r / (p + r)


# ---------------------------------------------------------------------------
# C.1 multi-seed CV variance
# ---------------------------------------------------------------------------

def cv_oof_for_seed(df_cell: pd.DataFrame, bam: str, flag: str, seed: int):
    """Replicate delta_f1.cv_out_of_fold_predictions with configurable seed.
    Returns (region_id_arr, y, p_oof, status). Skipped folds backfilled
    with in-sample predictions (same as Step 3 behaviour)."""
    X, _, mask = assemble_design(df_cell, bam, flag, "B")
    y_full = (df_cell["label"] == "TP").astype(int).values
    y = y_full[mask]
    X_use = X[mask]
    n_use = mask.sum()
    if n_use < CV_FOLDS * 2 or len(np.unique(y)) < 2:
        return np.array([], dtype=object), np.array([]), np.array([]), "insufficient_data"

    kf = KFold(n_splits=CV_FOLDS, shuffle=True, random_state=seed)
    p_oof = np.full(n_use, np.nan)
    n_fold_ok = 0
    for tr_idx, te_idx in kf.split(X_use):
        y_tr = y[tr_idx]
        if len(np.unique(y_tr)) < 2:
            continue
        r, _ = fit_lr(y_tr, X_use[tr_idx])
        if r is None:
            continue
        Xc_te = sm.add_constant(X_use[te_idx], has_constant="add")
        try:
            p_oof[te_idx] = r.predict(Xc_te)
            n_fold_ok += 1
        except Exception:
            continue
    if n_fold_ok == 0:
        return np.array([], dtype=object), np.array([]), np.array([]), "all_folds_failed"

    if np.isnan(p_oof).any():
        r_all, _ = fit_lr(y, X_use)
        if r_all is not None:
            Xc_all = sm.add_constant(X_use, has_constant="add")
            try:
                p_all = r_all.predict(Xc_all)
                nan_mask = np.isnan(p_oof)
                p_oof[nan_mask] = p_all[nan_mask]
            except Exception:
                pass
        if np.isnan(p_oof).any():
            keep = ~np.isnan(p_oof)
            df_idx = np.where(mask)[0][keep]
            return df_cell["region_id"].values[df_idx], y[keep], p_oof[keep], f"ok_partial_n_ok={n_fold_ok}"

    df_idx = np.where(mask)[0]
    return df_cell["region_id"].values[df_idx], y, p_oof, f"ok_n_ok={n_fold_ok}"


def run_c1_seed(seed: int, df: pd.DataFrame, fp_rich: list[dict], keep_cell_ids: set) -> dict:
    """Re-run Step 3 aggregated + per-cell sweep for one seed.
    Returns dict with per-seed Strategy 1 results."""
    bam, flag = PRIMARY_BAM, PRIMARY_FLAG

    per_cell_results = {}
    agg_records = []
    seen_region_ids: set = set()

    for cell in fp_rich:
        cell_id = cell["cell_id"]
        if cell_id not in keep_cell_ids:
            continue
        idx = select_cell_rows(df, cell)
        df_cell = df.loc[idx].copy()
        if len(df_cell) < 30:
            continue
        rg_ids, y, p_oof, status = cv_oof_for_seed(df_cell, bam, flag, seed)
        if len(p_oof) == 0:
            continue
        per_cell_results[cell_id] = (rg_ids, y, p_oof)
        # accumulate aggregated pool (dedup by region_id)
        for rg, y_val, p_val in zip(rg_ids, y, p_oof):
            if rg in seen_region_ids:
                continue
            seen_region_ids.add(rg)
            agg_records.append((rg, int(y_val), float(p_val)))

    # Strategy 1: max ΔF1 across all per-cell + aggregated configurations
    best_delta = -np.inf
    best_cell = ""
    best_tau = np.nan

    def evaluate_pool(label: str, y_arr: np.ndarray, p_arr: np.ndarray):
        nonlocal best_delta, best_cell, best_tau
        if len(y_arr) == 0:
            return
        tp_total = int((y_arr == 1).sum())
        fp_total = int((y_arr == 0).sum())
        for tau in TAU_GRID:
            kept = (p_arr >= tau)
            tp_kept = int(((y_arr == 1) & kept).sum())
            fp_kept = int(((y_arr == 0) & kept).sum())
            tp_removed = tp_total - tp_kept
            fp_removed = fp_total - fp_kept
            tp_post = TP_CALLER - tp_removed
            fp_post = FP_CALLER - fp_removed
            fn_post = FN_CALLER + tp_removed
            f1_post = f1_from_counts(tp_post, fp_post, fn_post)
            delta = f1_post - F1_CALLER
            if delta > best_delta:
                best_delta = delta
                best_cell = label
                best_tau = float(tau)

    # per-cell
    for cell_id, (_, y, p) in per_cell_results.items():
        evaluate_pool(cell_id, y, p)
    # aggregated
    if agg_records:
        agg_df = pd.DataFrame(agg_records, columns=["region_id", "y", "p_oof"])
        evaluate_pool("AGGREGATED", agg_df["y"].values, agg_df["p_oof"].values)

    return {
        "seed": seed,
        "max_delta_F1": float(best_delta) if best_delta != -np.inf else float("nan"),
        "winning_cell": best_cell,
        "winning_tau": best_tau,
        "n_cells_used": len(per_cell_results),
        "n_agg_rows": len(agg_records),
    }


def run_c1(df: pd.DataFrame, fp_rich: list[dict], keep_cell_ids: set) -> dict:
    s5d_log("=== C.1 START ===")
    rows = []
    for seed in SEEDS:
        s5d_log(f"  seed {seed}...")
        rec = run_c1_seed(seed, df, fp_rich, keep_cell_ids)
        s5d_log(
            f"  seed {seed}: max_ΔF1 = {rec['max_delta_F1']:.6f}, "
            f"winner = {rec['winning_cell']}, τ* = {rec['winning_tau']:.2f}"
        )
        rows.append(rec)
    df_c1 = pd.DataFrame(rows)
    df_c1.to_csv(OUT_C1, sep="\t", index=False, float_format="%.6g")
    s5d_log(f"Wrote {OUT_C1}")

    deltas = df_c1["max_delta_F1"].values
    mean_d = float(np.nanmean(deltas))
    std_d = float(np.nanstd(deltas, ddof=1))
    if std_d > 0.001:
        verdict = "RED"
    elif std_d > 0.0005:
        verdict = "YELLOW"
    else:
        verdict = "GREEN"

    unique_cells = df_c1["winning_cell"].nunique()
    unique_taus = df_c1["winning_tau"].nunique()
    consistent = unique_cells == 1 and unique_taus <= 2  # τ may jiggle 1 step

    s5d_log(
        f"C.1 verdict: mean = {mean_d:.6f} ± {std_d:.6f}  -> {verdict}  "
        f"(unique winners = {unique_cells}, unique τ = {unique_taus})"
    )

    return {
        "mean": mean_d,
        "std": std_d,
        "verdict": verdict,
        "df": df_c1,
        "unique_winners": int(unique_cells),
        "unique_taus": int(unique_taus),
        "consistent": bool(consistent),
    }


# ---------------------------------------------------------------------------
# C.4 unique cell count
# ---------------------------------------------------------------------------

def run_c4() -> dict:
    s5d_log("=== C.4 START ===")
    step1 = pd.read_csv(STEP5_DIR / "step1_lrt_per_cell.tsv", sep="\t")
    sig = step1[step1["LRT_q"] < 0.05].copy()
    total_rows = len(sig)

    # cell_id collapses BAM × flag
    grouped = (
        sig.groupby("cell_id")
        .apply(
            lambda g: pd.Series(
                {
                    "n_BAM_flag_combos_significant": len(g),
                    "BAM_flag_combos": ",".join(
                        sorted(g["BAM"].astype(str) + "_" + g["flag"].astype(str))
                    ),
                    "min_LRT_q": float(g["LRT_q"].min()),
                    "max_LRT_q": float(g["LRT_q"].max()),
                }
            )
        )
        .reset_index()
        .sort_values("min_LRT_q")
    )

    grouped.to_csv(OUT_C4, sep="\t", index=False, float_format="%.6g")
    s5d_log(f"Wrote {OUT_C4}")

    unique_cells = len(grouped)
    s5d_log(f"q<0.05 rows: {total_rows}, unique cells: {unique_cells}")
    for _, r in grouped.iterrows():
        s5d_log(
            f"  {r['cell_id']}: {int(r['n_BAM_flag_combos_significant'])}/6 combos, "
            f"min q = {r['min_LRT_q']:.3e}"
        )

    if unique_cells >= 4:
        verdict = "GREEN"
    elif unique_cells <= 2:
        verdict = "YELLOW"
    else:
        verdict = "YELLOW_LIGHT"

    s5d_log(f"C.4 verdict: unique cells = {unique_cells} -> {verdict}")
    return {
        "total_sig_rows": int(total_rows),
        "unique_cells": int(unique_cells),
        "df": grouped,
        "verdict": verdict,
    }


# ---------------------------------------------------------------------------
# C.7 NaN rate per cell × covariate × BAM × flag
# ---------------------------------------------------------------------------

def run_c7(df: pd.DataFrame, fp_rich: list[dict]) -> dict:
    s5d_log("=== C.7 START ===")

    # Cells in scope:
    #   1. 4 Step 3 winning per-cell from intermediate/step3_optimal_tau.tsv
    #   2. C.4 unique cells (step1 LRT q<0.05) — partition different (V6_off canonical),
    #      so we resolve them on df['cell_id'] directly.
    step3_opt = pd.read_csv(STEP5_DIR / "intermediate" / "step3_optimal_tau.tsv", sep="\t")
    step3_cells = step3_opt[
        (step3_opt["strategy"] == "max_DeltaF1") & (step3_opt["scope"] == "per_cell")
    ]["cell_id"].unique().tolist()

    step1 = pd.read_csv(STEP5_DIR / "step1_lrt_per_cell.tsv", sep="\t")
    step1_sig_cells = step1[step1["LRT_q"] < 0.05]["cell_id"].unique().tolist()

    rows = []

    # 1) Step 3 cells (use FP-rich cell descriptors → select_cell_rows)
    fp_rich_by_id = {c["cell_id"]: c for c in fp_rich}
    for cell_id in step3_cells:
        cell = fp_rich_by_id.get(cell_id)
        if cell is None:
            s5d_log(f"  WARN: step3 cell {cell_id} not in fp_rich list (skip)")
            continue
        idx = select_cell_rows(df, cell)
        df_cell = df.loc[idx].copy()
        n_total = len(df_cell)
        for bam in ("V3F", "V5", "V6"):
            for flag in ("off", "on"):
                _, _, mask = assemble_design(df_cell, bam, flag, "B")
                n_eff = int(mask.sum())
                # also per-covariate NaN rate
                cov_nans = {}
                for cov in METHYL_COVARIATES:
                    col = f"{bam}_{flag}_meth_{cov}"
                    series = pd.to_numeric(df_cell[col], errors="coerce")
                    cov_nans[cov] = float(series.isna().mean()) if n_total > 0 else float("nan")
                rows.append(
                    {
                        "scope": "step3_winning",
                        "cell_id": cell_id,
                        "BAM": bam,
                        "flag": flag,
                        "n_total": int(n_total),
                        "n_eff_B": n_eff,
                        "n_eff_over_total": float(n_eff / n_total) if n_total > 0 else float("nan"),
                        "nan_rate_HPMergedDelta": cov_nans["HPMergedDelta"],
                        "nan_rate_HPFineF": cov_nans["HPFineF"],
                        "nan_rate_NME_imbalance": cov_nans["NME_imbalance"],
                        "nan_rate_Epipoly_Delta": cov_nans["Epipoly_Delta"],
                        "nan_rate_ClusterPermanovaF": cov_nans["ClusterPermanovaF"],
                    }
                )

    # 2) Step 1 LRT-significant cells (V6_off canonical cell_id partition; select via df['cell_id'])
    for cell_id in step1_sig_cells:
        sub = df[df["cell_id"] == cell_id]
        n_total = len(sub)
        if n_total == 0:
            s5d_log(f"  WARN: step1 cell {cell_id} not in df (skip)")
            continue
        for bam in ("V3F", "V5", "V6"):
            for flag in ("off", "on"):
                _, _, mask = assemble_design(sub, bam, flag, "B")
                n_eff = int(mask.sum())
                cov_nans = {}
                for cov in METHYL_COVARIATES:
                    col = f"{bam}_{flag}_meth_{cov}"
                    series = pd.to_numeric(sub[col], errors="coerce")
                    cov_nans[cov] = float(series.isna().mean()) if n_total > 0 else float("nan")
                rows.append(
                    {
                        "scope": "step1_lrt_sig",
                        "cell_id": cell_id,
                        "BAM": bam,
                        "flag": flag,
                        "n_total": int(n_total),
                        "n_eff_B": n_eff,
                        "n_eff_over_total": float(n_eff / n_total) if n_total > 0 else float("nan"),
                        "nan_rate_HPMergedDelta": cov_nans["HPMergedDelta"],
                        "nan_rate_HPFineF": cov_nans["HPFineF"],
                        "nan_rate_NME_imbalance": cov_nans["NME_imbalance"],
                        "nan_rate_Epipoly_Delta": cov_nans["Epipoly_Delta"],
                        "nan_rate_ClusterPermanovaF": cov_nans["ClusterPermanovaF"],
                    }
                )

    df_c7 = pd.DataFrame(rows)
    df_c7 = df_c7.sort_values(["scope", "cell_id", "BAM", "flag"]).reset_index(drop=True)
    df_c7.to_csv(OUT_C7, sep="\t", index=False, float_format="%.6g")
    s5d_log(f"Wrote {OUT_C7}")

    # Verdict based on aggregate n_eff/n_total across all rows
    valid = df_c7.dropna(subset=["n_eff_over_total"])
    mean_eff = float(valid["n_eff_over_total"].mean()) if len(valid) else float("nan")
    median_eff = float(valid["n_eff_over_total"].median()) if len(valid) else float("nan")
    p25_eff = float(valid["n_eff_over_total"].quantile(0.25)) if len(valid) else float("nan")
    worst = valid.loc[valid["n_eff_over_total"].idxmin()] if len(valid) else None

    # Verdict uses median (most representative central tendency)
    if median_eff < 0.5:
        verdict = "RED"
    elif median_eff < 0.7:
        verdict = "YELLOW"
    else:
        verdict = "GREEN"

    s5d_log(
        f"C.7 verdict: median n_eff/n_total = {median_eff:.3f}  -> {verdict}  "
        f"(mean = {mean_eff:.3f}, p25 = {p25_eff:.3f})"
    )
    if worst is not None:
        s5d_log(
            f"  worst row: {worst['cell_id']} {worst['BAM']}_{worst['flag']} "
            f"n_eff/n_total = {worst['n_eff_over_total']:.3f} "
            f"(n_eff = {int(worst['n_eff_B'])}/{int(worst['n_total'])})"
        )

    return {
        "df": df_c7,
        "mean_eff": mean_eff,
        "median_eff": median_eff,
        "p25_eff": p25_eff,
        "worst": worst.to_dict() if worst is not None else None,
        "verdict": verdict,
    }


# ---------------------------------------------------------------------------
# Synthesis markdown
# ---------------------------------------------------------------------------

def overall_verdict(c1_v: str, c4_v: str, c7_v: str) -> tuple[str, str]:
    """Combine 3 verdicts into overall."""
    weight = {"GREEN": 0, "YELLOW_LIGHT": 1, "YELLOW": 2, "RED": 3}
    s = weight[c1_v] + weight[c4_v] + weight[c7_v]
    if s == 0:
        return "GREEN", "All three audits pass. ΔF1 result is statistically robust."
    if s <= 2:
        return (
            "GREEN_with_caveats",
            "Two of three audits pass strictly; remaining audit shows a soft margin. "
            "Effect size is small but the underlying signal is consistent.",
        )
    if s <= 4:
        return (
            "YELLOW",
            "Mixed evidence. At least one audit raises a concern that should be flagged "
            "in any external claim about ΔF1.",
        )
    return (
        "RED",
        "Multiple audits fail. Effect size is too fragile to support an external claim.",
    )


def write_markdown(c1: dict, c4: dict, c7: dict, master_dF1: float) -> None:
    overall, overall_msg = overall_verdict(c1["verdict"], c4["verdict"], c7["verdict"])

    lines: list[str] = []
    lines.append("# Step 5d — Robustness Audit (C.1 / C.4 / C.7)")
    lines.append("")
    lines.append("> Generated by `scripts/robustness_audit.py`.")
    lines.append(
        "> Audits the Step 5 Methylation-Augmented FP Filter Pilot ΔF1 = "
        f"+{master_dF1:.5f} result against three failure modes: "
        "(1) CV-fold randomness, (2) LRT signal concentration, (3) effective-n loss from NaN."
    )
    lines.append("")
    lines.append("## TL;DR")
    lines.append("")
    lines.append("| Audit | Verdict | Headline metric |")
    lines.append("|-------|---------|------------------|")
    lines.append(
        f"| C.1 Multi-seed CV variance | **{c1['verdict']}** | "
        f"ΔF1 = {c1['mean']:+.6f} ± {c1['std']:.6f} across 5 seeds |"
    )
    lines.append(
        f"| C.4 H1 LRT unique cell count | **{c4['verdict']}** | "
        f"{c4['unique_cells']} unique cells across {c4['total_sig_rows']} q<0.05 rows |"
    )
    lines.append(
        f"| C.7 NaN rate per cell × covariate | **{c7['verdict']}** | "
        f"median n_eff/n_total = {c7['median_eff']:.3f} |"
    )
    lines.append(f"| **Overall** | **{overall}** | {overall_msg} |")
    lines.append("")

    # --- §1 C.1 ---
    lines.append("## §1 C.1 — Multi-seed CV variance")
    lines.append("")
    lines.append("**Question:** Is ΔF1 = +0.00242 produced by random KFold splits rather than by genuine signal?")
    lines.append("")
    lines.append("**Method:** Re-run Step 3's CV-based ΔF1 sweep with 5 different KFold seeds "
                 f"({', '.join(str(s) for s in SEEDS)}). For each seed, record the Strategy 1 (max ΔF1) "
                 "result across per-cell + aggregated scopes.")
    lines.append("")
    lines.append("**Per-seed results:**")
    lines.append("")
    lines.append("| seed | max ΔF1 | winning cell | winning τ* | n_cells | n_agg_rows |")
    lines.append("|------|---------|--------------|-------------|---------|------------|")
    for _, r in c1["df"].iterrows():
        lines.append(
            f"| {int(r['seed'])} | {r['max_delta_F1']:+.6f} | {r['winning_cell']} | "
            f"{r['winning_tau']:.2f} | {int(r['n_cells_used'])} | {int(r['n_agg_rows'])} |"
        )
    lines.append("")
    lines.append(f"- mean ± std (n=5) = **{c1['mean']:+.6f} ± {c1['std']:.6f}**")
    lines.append(f"- unique winning cells across seeds = **{c1['unique_winners']}**")
    lines.append(f"- unique winning τ values across seeds = **{c1['unique_taus']}**")
    lines.append("")
    lines.append("**Decision thresholds:**")
    lines.append("- std ≤ 0.0005 → GREEN (robust)")
    lines.append("- 0.0005 < std ≤ 0.001 → YELLOW (marginal)")
    lines.append("- std > 0.001 → RED (effect not robust)")
    lines.append("")
    lines.append(f"**§1 Verdict: {c1['verdict']}**")
    lines.append("")

    # --- §2 C.4 ---
    lines.append("## §2 C.4 — H1 LRT unique cell count")
    lines.append("")
    lines.append("**Question:** The 16 q<0.05 rows in `step1_lrt_per_cell.tsv` — are these 16 unique cells, "
                 "or fewer cells repeated across BAM × flag combinations (3 BAMs × 2 flags = 6 combos per cell)?")
    lines.append("")
    lines.append("**Method:** Collapse `BAM × flag` dimension; count unique `cell_id` with at least one "
                 "q<0.05 hit. Report how many combos (out of 6) each unique cell hit.")
    lines.append("")
    lines.append(f"- Total q<0.05 rows: **{c4['total_sig_rows']}**")
    lines.append(f"- Unique cells: **{c4['unique_cells']}**")
    lines.append("")
    lines.append("**Per-cell consistency (combos firing q<0.05):**")
    lines.append("")
    lines.append("| cell_id | combos significant (of 6) | min LRT q | max LRT q | combos |")
    lines.append("|---------|----------------------------|-----------|-----------|--------|")
    for _, r in c4["df"].iterrows():
        lines.append(
            f"| {r['cell_id']} | {int(r['n_BAM_flag_combos_significant'])}/6 | "
            f"{r['min_LRT_q']:.3e} | {r['max_LRT_q']:.3e} | {r['BAM_flag_combos']} |"
        )
    lines.append("")
    lines.append("**Decision thresholds:**")
    lines.append("- ≥ 4 unique cells → GREEN (signal distributed across cells)")
    lines.append("- 3 unique cells → YELLOW_LIGHT")
    lines.append("- ≤ 2 unique cells → YELLOW (signal concentrated in few cells)")
    lines.append("")
    lines.append(f"**§2 Verdict: {c4['verdict']}**")
    lines.append("")

    # --- §3 C.7 ---
    lines.append("## §3 C.7 — NaN rate per cell × covariate")
    lines.append("")
    lines.append("**Question:** With Step 0 sanity showing 19-40% NaN in methylation covariates, what is "
                 "the effective n / n_total for Model B per (cell × BAM × flag)? If < 50%, the LRT signal "
                 "may be a sample-selection artefact.")
    lines.append("")
    lines.append("**Method:** For (Step 3 winning per-cell cells) ∪ (Step 1 LRT q<0.05 cells) × (3 BAMs) × "
                 "(2 flags), compute n_eff_B (Model B design matrix non-NaN rows) and per-covariate NaN rate.")
    lines.append("")
    lines.append(f"- median n_eff/n_total = **{c7['median_eff']:.3f}** (across all rows)")
    lines.append(f"- mean n_eff/n_total = {c7['mean_eff']:.3f}")
    lines.append(f"- 25th percentile n_eff/n_total = {c7['p25_eff']:.3f}")
    if c7["worst"] is not None:
        w = c7["worst"]
        lines.append(
            f"- worst row: `{w['cell_id']}` {w['BAM']}_{w['flag']} → n_eff/n_total = "
            f"{w['n_eff_over_total']:.3f} (n_eff = {int(w['n_eff_B'])} / n_total = {int(w['n_total'])})"
        )
    lines.append("")
    lines.append("**Full table:** see `intermediate/step5d_c7_nan_rate.tsv`.")
    lines.append("")
    lines.append("**Top 10 worst rows (lowest n_eff/n_total):**")
    lines.append("")
    lines.append("| scope | cell_id | BAM | flag | n_total | n_eff_B | n_eff/n_total | NaN(HPMergedDelta) | NaN(HPFineF) | NaN(NME_imb) | NaN(Epi_Delta) | NaN(PermF) |")
    lines.append("|-------|---------|-----|------|---------|---------|---------------|---------------------|---------------|---------------|------------------|-------------|")
    worst10 = c7["df"].sort_values("n_eff_over_total").head(10)
    for _, r in worst10.iterrows():
        lines.append(
            f"| {r['scope']} | {r['cell_id']} | {r['BAM']} | {r['flag']} | "
            f"{int(r['n_total'])} | {int(r['n_eff_B'])} | {r['n_eff_over_total']:.3f} | "
            f"{r['nan_rate_HPMergedDelta']:.3f} | {r['nan_rate_HPFineF']:.3f} | "
            f"{r['nan_rate_NME_imbalance']:.3f} | {r['nan_rate_Epipoly_Delta']:.3f} | "
            f"{r['nan_rate_ClusterPermanovaF']:.3f} |"
        )
    lines.append("")
    lines.append("**Decision thresholds:**")
    lines.append("- median ≥ 0.7 → GREEN (most regions retained)")
    lines.append("- 0.5 ≤ median < 0.7 → YELLOW (notable drop, possible selection bias)")
    lines.append("- median < 0.5 → RED (majority of regions dropped, LR signal sample-selection)")
    lines.append("")
    lines.append(f"**§3 Verdict: {c7['verdict']}**")
    lines.append("")

    # --- §4 Overall ---
    lines.append("## §4 Overall robustness")
    lines.append("")
    lines.append(f"- C.1 (multi-seed): **{c1['verdict']}**")
    lines.append(f"- C.4 (unique cells): **{c4['verdict']}**")
    lines.append(f"- C.7 (n_eff/n_total): **{c7['verdict']}**")
    lines.append("")
    lines.append(f"**Combined verdict: {overall}**")
    lines.append("")
    lines.append(overall_msg)
    lines.append("")

    # --- §5 impact on Step 5 verdict ---
    lines.append("## §5 Impact on Step 5 verdict / paper claim")
    lines.append("")
    if overall == "GREEN":
        lines.append(
            "- Step 5 ⭐3 candidate verdict can stand as written: ΔF1 = +0.00242 is robust "
            "to CV-fold randomness, distributed across multiple cells, and not driven by sample selection."
        )
    elif overall == "GREEN_with_caveats":
        lines.append(
            "- Step 5 ⭐3 candidate verdict **stands** but the report MUST disclose:"
        )
        if c1["verdict"] != "GREEN":
            lines.append(
                f"  - ΔF1 has cross-seed std = {c1['std']:.6f} "
                "(comparable magnitude to the point estimate itself)."
            )
        if c4["verdict"] != "GREEN":
            lines.append(
                f"  - H1 LRT signal is concentrated in {c4['unique_cells']} unique cells "
                "(not 16 independent cells)."
            )
        if c7["verdict"] != "GREEN":
            lines.append(
                f"  - Effective sample size is {c7['median_eff']:.0%} of total — "
                "Model B fit is on a NaN-filtered subset."
            )
        lines.append("- No tier upgrade until 7-sample cross-validation (P4) gates this on independent samples.")
    elif overall == "YELLOW":
        lines.append(
            "- Step 5 ⭐3 candidate is **degraded to ⭐2 (mechanism candidate only)** until concerns resolved."
        )
        lines.append(
            "- External claim of 'methylation-augmented filter improves F1' is NOT supportable on these data alone."
        )
        lines.append(
            "- Recommended next step: 7-sample P4 generalisation. If concern persists across samples → mark NEGATIVE."
        )
    else:
        lines.append("- Step 5 ⭐3 candidate is **withdrawn**. ΔF1 effect is too fragile to claim.")
        lines.append("- Recommended: pivot away from this filter design.")
    lines.append("")
    lines.append("## Files")
    lines.append(f"- `step5d_robustness_audit.md` (this file)")
    lines.append(f"- `intermediate/step5d_c1_multiseed_variance.tsv` — 5 seeds × ΔF1/winner/τ*")
    lines.append(f"- `intermediate/step5d_c4_unique_cells.tsv` — q<0.05 cell consistency")
    lines.append(f"- `intermediate/step5d_c7_nan_rate.tsv` — per (cell × BAM × flag) NaN diagnostics")
    lines.append(f"- `intermediate/step5d_log.txt` — execution log")
    lines.append("")
    lines.append("## Constraints honoured")
    lines.append("- Read-only on Step 0-3 artefacts (`step1_lrt_per_cell.tsv`, `step3_delta_f1.tsv`, "
                 "`step3_optimal_tau.tsv`, `step5_master_augmented.tsv`).")
    lines.append("- BH-FDR not re-computed (uses Step 1 LRT q as-is).")
    lines.append("- No emoji in output files.")
    lines.append("- Pandas + scipy + numpy + statsmodels + sklearn only (matches Step 1-3 stack).")

    OUT_MD.write_text("\n".join(lines))
    s5d_log(f"Wrote {OUT_MD}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    LOG_PATH.unlink(missing_ok=True)
    s5d_log("=== Step 5d robustness_audit START ===")

    # Load master TSV once (same as delta_f1.py)
    s5d_log(f"Loading {MASTER_AUG}")
    df = pd.read_csv(MASTER_AUG, sep="\t", low_memory=False)
    sig_cols = [
        c
        for c in df.columns
        if c.endswith("_meth_HPMergedSig")
        or c.endswith("_meth_HPFineSig")
        or c.endswith("_meth_AlleleP")
    ]
    for c in sig_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = annotate_v6off_cell(df)
    df["cov_proxy"] = compute_cov_proxy(
        df["V6_off_n_reads"].fillna(0).astype(float)
    )
    fp_density = pd.read_csv(FP_DENSITY_FILE, sep="\t")
    df = df.merge(
        fp_density[["region_id", "fp_density", "in_Z_AUTO"]],
        on="region_id",
        how="left",
    )
    df["in_Z_AUTO"] = df["in_Z_AUTO"].fillna(False).astype(bool)
    s5d_log(f"Master rows: {len(df):,}")

    fp_rich = find_fp_rich_cells()
    step2_summary = pd.read_csv(
        STEP5_DIR / "intermediate" / "step2_cell_summary.tsv", sep="\t"
    )
    step2_summary["passes_gate"] = (
        (step2_summary["model_used"] == "B")
        & (step2_summary["n_fit"] >= MIN_N_EFF_B)
        & (step2_summary["n_TP_fit"] >= MIN_N_TP)
    )
    keep_cell_ids = set(
        step2_summary.loc[step2_summary["passes_gate"], "cell_id"].tolist()
    )
    s5d_log(f"Cells passing Step 3 gate: {len(keep_cell_ids)}")
    for cid in sorted(keep_cell_ids):
        s5d_log(f"  KEEP: {cid}")

    # Step 3 master ΔF1 baseline
    step3_opt = pd.read_csv(
        STEP5_DIR / "intermediate" / "step3_optimal_tau.tsv", sep="\t"
    )
    s1 = step3_opt[step3_opt["strategy"] == "max_DeltaF1"]
    master_dF1 = float(s1["delta_F1"].max())
    s5d_log(f"Step 3 master max ΔF1 = {master_dF1:+.6f}")

    # Run audits
    c1 = run_c1(df, fp_rich, keep_cell_ids)
    c4 = run_c4()
    c7 = run_c7(df, fp_rich)

    write_markdown(c1, c4, c7, master_dF1)

    s5d_log("=== Step 5d robustness_audit END ===")
    print("\n=== STEP 5d ROBUSTNESS AUDIT SUMMARY ===")
    print(f"C.1 multi-seed:  mean ± std = {c1['mean']:+.6f} ± {c1['std']:.6f}  -> {c1['verdict']}")
    print(f"C.4 unique cells: {c4['unique_cells']} unique  -> {c4['verdict']}")
    print(f"C.7 n_eff/n_total: median = {c7['median_eff']:.3f}  -> {c7['verdict']}")
    overall, msg = overall_verdict(c1["verdict"], c4["verdict"], c7["verdict"])
    print(f"OVERALL: {overall}  ({msg})")


if __name__ == "__main__":
    main()
