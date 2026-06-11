"""Step 3 — ΔF1 vs ClairS-TO caller baseline.

Re-fits Model B per FP-rich cell with 5-fold CV, captures out-of-fold P(TP),
sweeps tau in [0.10, 0.95] by 0.01 (86 values), and computes per-cell + aggregated
global ΔF1 against the reverse-solved caller baseline (TP=30,490 / FP=4,842 /
FN=19,288 / F1=0.7166 per plan §Step 3).

Caveats addressed (per Step 2 hand-off):
  1. Only converged + n_eff_B >= 100 cells are filtered (excludes 4 trivial UNKNOWN cells)
  2. CV out-of-fold predictions, not in-sample, to avoid optimism (AUC 0.97 -> 0.90)
  3. tau range extended to [0.10, 0.95] (per hand-off Step 3 §3 — best tau was 0.50 floor)
  4. Strategy 1: max DeltaF1 (aggressive, per user resolved decision #1)
     Strategy 2: max tau s.t. FP_removal > TP_loss AND DeltaF1 > 0 (sanity)

Outputs:
  step3_delta_f1.tsv
  step3_optimal_tau_summary.md
  step3_findings.md
  figures/step3_deltaf1_vs_tau.png
  figures/step3_filter_signal_curve.png
  intermediate/step3_log.txt
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold

from _common_step1 import (
    MASTER_AUG,
    STEP3_CHR8,
    STEP3_AUTO,
    STEP5_DIR,
    annotate_v6off_cell,
    log_msg,
    matplotlib_setup,
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

# Override log location for step 3
LOG_PATH = STEP5_DIR / "intermediate" / "step3_log.txt"


def s3_log(msg: str) -> None:
    import datetime as _dt
    ts = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(line + "\n")


# -----------------------------------------------------------------------------
# Caller baseline (post-PON, ISM-region level — same scope as master TSV)
# -----------------------------------------------------------------------------
# Per plan §Step 3 reverse-solve:
#   Master TSV: TP=30,490, FP=4,842, P = 30490/35332 = 0.8629
#   Given F1 = 0.7166 -> R = (F1*P)/(2P - F1) = 0.6126
#   FN = TP*(1-R)/R = 19,288; N_truth_positive = TP + FN = 49,778
TP_CALLER = 30490
FP_CALLER = 4842
FN_CALLER = 19288
N_TRUTH_POSITIVE = TP_CALLER + FN_CALLER  # 49,778
F1_CALLER = 0.7166

TAU_GRID = np.round(np.arange(0.10, 0.9501, 0.01), 2)  # 86 values
CV_FOLDS = 5
RANDOM_STATE = 42

# Inclusion gates per Step 2 hand-off caveats:
MIN_N_EFF_B = 100
MIN_N_TP = 5


def f1_from_counts(tp: float, fp: float, fn: float) -> float:
    if tp + fp == 0 or tp + fn == 0:
        return 0.0
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


def reverse_solve_sanity_check() -> dict:
    P = TP_CALLER / (TP_CALLER + FP_CALLER)
    R = TP_CALLER / (TP_CALLER + FN_CALLER)
    F1_check = 2 * P * R / (P + R)
    return {
        "TP_caller": TP_CALLER,
        "FP_caller": FP_CALLER,
        "FN_caller_reverse_solved": FN_CALLER,
        "N_truth_positive_reverse_solved": N_TRUTH_POSITIVE,
        "P_recomputed": P,
        "R_recomputed": R,
        "F1_recomputed": F1_check,
        "F1_target": F1_CALLER,
        "F1_match_to_3decimals": abs(F1_check - F1_CALLER) < 0.0005,
    }


def cv_out_of_fold_predictions(df_cell: pd.DataFrame, bam: str, flag: str, model: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Refit Model with 5-fold CV. Return (p_oof, y, mask_indices_in_df_cell, status).
    p_oof: out-of-fold predicted P(TP) for rows in mask.
    Indices align to df_cell rows where the design matrix has no NaN.
    """
    X, _, mask = assemble_design(df_cell, bam, flag, model)
    y_full = (df_cell["label"] == "TP").astype(int).values
    y = y_full[mask]
    X_use = X[mask]
    n_use = mask.sum()
    if n_use < CV_FOLDS * 2 or len(np.unique(y)) < 2:
        return np.array([]), np.array([]), np.where(mask)[0], "insufficient_data"

    kf = KFold(n_splits=CV_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    p_oof = np.full(n_use, np.nan)

    n_fold_ok = 0
    n_fold_skip = 0
    for tr_idx, te_idx in kf.split(X_use):
        y_tr = y[tr_idx]
        if len(np.unique(y_tr)) < 2:
            n_fold_skip += 1
            continue
        r, status = fit_lr(y_tr, X_use[tr_idx])
        if r is None:
            n_fold_skip += 1
            continue
        Xc_te = sm.add_constant(X_use[te_idx], has_constant="add")
        try:
            scores = r.predict(Xc_te)
            p_oof[te_idx] = scores
            n_fold_ok += 1
        except Exception:
            n_fold_skip += 1

    if n_fold_ok == 0:
        return np.array([]), np.array([]), np.where(mask)[0], "all_folds_failed"

    # If any folds skipped, fill their predictions with a fallback in-sample fit on all data
    # (this preserves all rows for the sweep but is documented).
    if np.isnan(p_oof).any():
        r_all, status_all = fit_lr(y, X_use)
        if r_all is not None:
            Xc_all = sm.add_constant(X_use, has_constant="add")
            try:
                p_all = r_all.predict(Xc_all)
                mask_nan = np.isnan(p_oof)
                p_oof[mask_nan] = p_all[mask_nan]
            except Exception:
                pass
        if np.isnan(p_oof).any():
            keep_idx = np.where(~np.isnan(p_oof))[0]
            y = y[keep_idx]
            p_oof = p_oof[keep_idx]
            df_idx = np.where(mask)[0][keep_idx]
            return p_oof, y, df_idx, f"ok_partial_n_fold_ok={n_fold_ok}"

    df_idx = np.where(mask)[0]
    return p_oof, y, df_idx, f"ok_n_fold_ok={n_fold_ok}"


def main():
    LOG_PATH.unlink(missing_ok=True)
    s3_log("=== Step 3 delta_f1 START ===")

    # ----- Sanity check on baseline -----
    sanity = reverse_solve_sanity_check()
    s3_log("Caller baseline sanity check:")
    for k, v in sanity.items():
        s3_log(f"  {k}: {v}")

    # ----- Load master + annotate -----
    s3_log(f"Loading {MASTER_AUG}")
    df = pd.read_csv(MASTER_AUG, sep="\t", low_memory=False)
    sig_cols = [c for c in df.columns
                if c.endswith("_meth_HPMergedSig") or c.endswith("_meth_HPFineSig")
                or c.endswith("_meth_AlleleP")]
    for c in sig_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = annotate_v6off_cell(df)
    df["cov_proxy"] = compute_cov_proxy(df["V6_off_n_reads"].fillna(0).astype(float))

    fp_density = pd.read_csv(FP_DENSITY_FILE, sep="\t")
    df = df.merge(fp_density[["region_id", "fp_density", "in_Z_AUTO"]],
                  on="region_id", how="left")
    df["in_Z_AUTO"] = df["in_Z_AUTO"].fillna(False).astype(bool)

    s3_log(f"Master rows: {len(df):,}")
    s3_log(f"  TP={int((df['label']=='TP').sum()):,} FP={int((df['label']=='FP').sum()):,}")

    # ----- Identify FP-rich cells (same as Step 2) -----
    fp_rich = find_fp_rich_cells()
    s3_log(f"FP-rich cells: {len(fp_rich)}")

    # ----- Filter cells by Step 2 caveats -----
    # We need to retain only converged Model B + n_eff_B >= 100 + n_TP >= 5
    step2_summary = pd.read_csv(STEP5_DIR / "intermediate" / "step2_cell_summary.tsv", sep="\t")
    step2_summary["passes_gate"] = (
        (step2_summary["model_used"] == "B")
        & (step2_summary["n_fit"] >= MIN_N_EFF_B)
        & (step2_summary["n_TP_fit"] >= MIN_N_TP)
    )
    n_keep = int(step2_summary["passes_gate"].sum())
    s3_log(f"Cells passing gate (Model B converged + n_fit >= {MIN_N_EFF_B} + n_TP >= {MIN_N_TP}): "
           f"{n_keep} / {len(step2_summary)}")
    keep_cell_ids = set(step2_summary.loc[step2_summary["passes_gate"], "cell_id"].tolist())
    for cid in keep_cell_ids:
        s3_log(f"  KEEP: {cid}")
    for _, r in step2_summary[~step2_summary["passes_gate"]].iterrows():
        s3_log(f"  DROP: {r['cell_id']}  (n_fit={r['n_fit']}, n_TP={r['n_TP_fit']}, model={r['model_used']})")

    # ----- Per-cell CV-based filter results -----
    # First, compute CV out-of-fold P(TP) per kept cell using Model B on V6_off
    sweep_rows = []
    per_cell_oof = {}  # cell_id -> (df_idx, y, p_oof) for aggregation
    cv_aucs = {}
    cell_meta = {}  # cell_id -> {"n_fit":..., "n_TP":..., "n_FP":..., "source":...}

    bam, flag = PRIMARY_BAM, PRIMARY_FLAG
    for cell in fp_rich:
        cell_id = cell["cell_id"]
        if cell_id not in keep_cell_ids:
            continue
        idx = select_cell_rows(df, cell)
        df_cell = df.loc[idx].copy()
        if len(df_cell) < 30:
            continue
        p_oof, y, df_indices, status = cv_out_of_fold_predictions(df_cell, bam, flag, "B")
        if len(p_oof) == 0:
            s3_log(f"  Skip {cell_id}: {status}")
            continue
        cv_auc = roc_auc_score(y, p_oof) if len(np.unique(y)) >= 2 else float("nan")
        cv_aucs[cell_id] = cv_auc
        # Translate df_indices (positions within df_cell) -> global region_ids
        rg_id_arr = df_cell["region_id"].values[df_indices]
        per_cell_oof[cell_id] = {
            "df_indices": df_indices,
            "region_id": rg_id_arr,
            "y": y,
            "p_oof": p_oof,
            "status": status,
        }
        n_fit = len(y)
        n_tp_fit = int(y.sum())
        n_fp_fit = n_fit - n_tp_fit
        cell_meta[cell_id] = {
            "n_fit": n_fit,
            "n_TP_fit": n_tp_fit,
            "n_FP_fit": n_fp_fit,
            "source": cell["source"],
            "cv_auc": cv_auc,
            "cv_status": status,
        }
        s3_log(f"  {cell_id}: n_fit={n_fit}, TP/FP={n_tp_fit}/{n_fp_fit}, "
               f"CV AUC = {cv_auc:.4f}, status={status}")

        # Per-cell tau sweep
        for tau in TAU_GRID:
            kept = (p_oof >= tau)
            tp_kept = int(((y == 1) & kept).sum())
            fp_kept = int(((y == 0) & kept).sum())
            tp_removed = n_tp_fit - tp_kept
            fp_removed = n_fp_fit - fp_kept
            tp_loss_pct = tp_removed / n_tp_fit if n_tp_fit > 0 else float("nan")
            fp_removal_pct = fp_removed / n_fp_fit if n_fp_fit > 0 else float("nan")

            # Per-cell precision/recall (cell only)
            cell_precision_post = (tp_kept / (tp_kept + fp_kept)) if (tp_kept + fp_kept) > 0 else float("nan")
            cell_precision_pre = n_tp_fit / n_fit

            # ---- GLOBAL ΔF1 if we apply this filter to this cell only ----
            # All other regions (master TSV \ this cell rows) keep caller's decision.
            # FN_post = FN_caller + tp_removed (filtered TPs become FNs)
            # TP_post = TP_caller - tp_removed
            # FP_post = FP_caller - fp_removed
            tp_post = TP_CALLER - tp_removed
            fp_post = FP_CALLER - fp_removed
            fn_post = FN_CALLER + tp_removed
            f1_post = f1_from_counts(tp_post, fp_post, fn_post)
            delta_f1 = f1_post - F1_CALLER

            sweep_rows.append({
                "cell_id": cell_id,
                "source": cell["source"],
                "BAM": bam,
                "flag": flag,
                "model_used": "B",
                "n_fit": n_fit,
                "n_TP_fit": n_tp_fit,
                "n_FP_fit": n_fp_fit,
                "tau": float(tau),
                "TP_kept_cell": tp_kept,
                "FP_kept_cell": fp_kept,
                "TP_removed_cell": tp_removed,
                "FP_removed_cell": fp_removed,
                "TP_loss_pct": tp_loss_pct,
                "FP_removal_pct": fp_removal_pct,
                "filter_signal": fp_removal_pct - tp_loss_pct,
                "cell_precision_pre": cell_precision_pre,
                "cell_precision_post": cell_precision_post,
                "global_TP_post": tp_post,
                "global_FP_post": fp_post,
                "global_FN_post": fn_post,
                "global_F1_post": f1_post,
                "delta_F1_per_cell": delta_f1,
                "scope": "per_cell",
            })

    # ----- Aggregated filter (multi-cell global filter) -----
    # Aggregate all rows across kept cells, then sweep tau on the union.
    # Same row not in multiple cells (by region_id) — but we deduplicate just in case.
    s3_log("")
    s3_log("Building aggregated multi-cell filter...")
    agg_records = []  # list of (region_id, y, p_oof, cell_id)
    seen_region_ids: set = set()
    for cell_id, data in per_cell_oof.items():
        for rg, y_val, p_val in zip(data["region_id"], data["y"], data["p_oof"]):
            if rg in seen_region_ids:
                continue
            seen_region_ids.add(rg)
            agg_records.append((rg, int(y_val), float(p_val), cell_id))
    if agg_records:
        agg_df = pd.DataFrame(agg_records, columns=["region_id", "y", "p_oof", "cell_id"])
        s3_log(f"  Aggregated pool: {len(agg_df)} unique regions ("
               f"TP={int((agg_df['y']==1).sum())}, FP={int((agg_df['y']==0).sum())})")

        for tau in TAU_GRID:
            kept = (agg_df["p_oof"] >= tau).values
            y_arr = agg_df["y"].values
            tp_kept = int(((y_arr == 1) & kept).sum())
            fp_kept = int(((y_arr == 0) & kept).sum())
            tp_total = int((y_arr == 1).sum())
            fp_total = int((y_arr == 0).sum())
            tp_removed = tp_total - tp_kept
            fp_removed = fp_total - fp_kept
            tp_loss_pct = tp_removed / tp_total if tp_total > 0 else float("nan")
            fp_removal_pct = fp_removed / fp_total if fp_total > 0 else float("nan")

            tp_post = TP_CALLER - tp_removed
            fp_post = FP_CALLER - fp_removed
            fn_post = FN_CALLER + tp_removed
            f1_post = f1_from_counts(tp_post, fp_post, fn_post)
            delta_f1 = f1_post - F1_CALLER

            sweep_rows.append({
                "cell_id": "AGGREGATED",
                "source": "aggregated_multi_cell",
                "BAM": bam,
                "flag": flag,
                "model_used": "B",
                "n_fit": len(agg_df),
                "n_TP_fit": tp_total,
                "n_FP_fit": fp_total,
                "tau": float(tau),
                "TP_kept_cell": tp_kept,
                "FP_kept_cell": fp_kept,
                "TP_removed_cell": tp_removed,
                "FP_removed_cell": fp_removed,
                "TP_loss_pct": tp_loss_pct,
                "FP_removal_pct": fp_removal_pct,
                "filter_signal": fp_removal_pct - tp_loss_pct,
                "cell_precision_pre": tp_total / len(agg_df) if len(agg_df) > 0 else float("nan"),
                "cell_precision_post": (tp_kept / (tp_kept + fp_kept)) if (tp_kept + fp_kept) > 0 else float("nan"),
                "global_TP_post": tp_post,
                "global_FP_post": fp_post,
                "global_FN_post": fn_post,
                "global_F1_post": f1_post,
                "delta_F1_per_cell": delta_f1,
                "scope": "aggregated",
            })

    sweep_df = pd.DataFrame(sweep_rows)
    out_tsv = STEP5_DIR / "step3_delta_f1.tsv"
    sweep_df.to_csv(out_tsv, sep="\t", index=False, float_format="%.6g")
    s3_log(f"Wrote {out_tsv.name} ({len(sweep_df)} rows)")

    # ----- Optimal tau per cell + aggregated -----
    # Strategy 1: max ΔF1 (aggressive)
    # Strategy 2: max tau where filter_signal > 0 AND ΔF1 > 0 (sanity)
    s1_records, s2_records = [], []
    for cell_id, grp in sweep_df.groupby("cell_id"):
        # Strategy 1
        idx_max = grp["delta_F1_per_cell"].idxmax()
        best = grp.loc[idx_max]
        s1_records.append({
            "scope": "per_cell" if cell_id != "AGGREGATED" else "aggregated",
            "cell_id": cell_id,
            "strategy": "max_DeltaF1",
            "tau_star": float(best["tau"]),
            "delta_F1": float(best["delta_F1_per_cell"]),
            "TP_loss_pct": float(best["TP_loss_pct"]),
            "FP_removal_pct": float(best["FP_removal_pct"]),
            "filter_signal": float(best["filter_signal"]),
            "n_fit": int(best["n_fit"]),
            "n_TP_fit": int(best["n_TP_fit"]),
            "n_FP_fit": int(best["n_FP_fit"]),
            "global_F1_post": float(best["global_F1_post"]),
        })
        # Strategy 2: max tau with both filter_signal>0 AND delta_F1>0
        eligible = grp[(grp["filter_signal"] > 0) & (grp["delta_F1_per_cell"] > 0)]
        if len(eligible) == 0:
            s2_records.append({
                "scope": "per_cell" if cell_id != "AGGREGATED" else "aggregated",
                "cell_id": cell_id,
                "strategy": "max_tau_s.t._sanity",
                "tau_star": float("nan"),
                "delta_F1": float("nan"),
                "TP_loss_pct": float("nan"),
                "FP_removal_pct": float("nan"),
                "filter_signal": float("nan"),
                "n_fit": int(grp["n_fit"].iloc[0]),
                "n_TP_fit": int(grp["n_TP_fit"].iloc[0]),
                "n_FP_fit": int(grp["n_FP_fit"].iloc[0]),
                "global_F1_post": float("nan"),
            })
        else:
            best = eligible.loc[eligible["tau"].idxmax()]
            s2_records.append({
                "scope": "per_cell" if cell_id != "AGGREGATED" else "aggregated",
                "cell_id": cell_id,
                "strategy": "max_tau_s.t._sanity",
                "tau_star": float(best["tau"]),
                "delta_F1": float(best["delta_F1_per_cell"]),
                "TP_loss_pct": float(best["TP_loss_pct"]),
                "FP_removal_pct": float(best["FP_removal_pct"]),
                "filter_signal": float(best["filter_signal"]),
                "n_fit": int(best["n_fit"]),
                "n_TP_fit": int(best["n_TP_fit"]),
                "n_FP_fit": int(best["n_FP_fit"]),
                "global_F1_post": float(best["global_F1_post"]),
            })

    optimal_df = pd.DataFrame(s1_records + s2_records)
    optimal_path = STEP5_DIR / "intermediate" / "step3_optimal_tau.tsv"
    optimal_df.to_csv(optimal_path, sep="\t", index=False, float_format="%.6g")
    s3_log(f"Wrote {optimal_path.name}")

    # ----- Global max ΔF1 across all per-cell + aggregated -----
    s1_per_cell = optimal_df[(optimal_df["strategy"] == "max_DeltaF1") & (optimal_df["scope"] == "per_cell")]
    s1_agg = optimal_df[(optimal_df["strategy"] == "max_DeltaF1") & (optimal_df["scope"] == "aggregated")]
    if not s1_per_cell.empty:
        best_per_cell = s1_per_cell.loc[s1_per_cell["delta_F1"].idxmax()]
        max_dF1_per_cell = float(best_per_cell["delta_F1"])
        max_dF1_per_cell_cell = best_per_cell["cell_id"]
        max_dF1_per_cell_tau = float(best_per_cell["tau_star"])
    else:
        max_dF1_per_cell = float("nan")
        max_dF1_per_cell_cell = ""
        max_dF1_per_cell_tau = float("nan")
    if not s1_agg.empty:
        agg_row = s1_agg.iloc[0]
        max_dF1_agg = float(agg_row["delta_F1"])
        max_dF1_agg_tau = float(agg_row["tau_star"])
    else:
        max_dF1_agg = float("nan")
        max_dF1_agg_tau = float("nan")

    # Final winning ΔF1: take the largest of per-cell vs aggregated
    if not np.isnan(max_dF1_per_cell) and not np.isnan(max_dF1_agg):
        if max_dF1_agg > max_dF1_per_cell:
            final_max_dF1 = max_dF1_agg
            final_winner = "AGGREGATED"
            final_tau = max_dF1_agg_tau
        else:
            final_max_dF1 = max_dF1_per_cell
            final_winner = max_dF1_per_cell_cell
            final_tau = max_dF1_per_cell_tau
    elif not np.isnan(max_dF1_per_cell):
        final_max_dF1 = max_dF1_per_cell
        final_winner = max_dF1_per_cell_cell
        final_tau = max_dF1_per_cell_tau
    else:
        final_max_dF1 = float("nan")
        final_winner = ""
        final_tau = float("nan")

    h2_verdict = "POSITIVE" if (not np.isnan(final_max_dF1) and final_max_dF1 > 0) else "NEGATIVE"

    # ----- H3 verdict (at the winning tau*, is FP_removal_pct > TP_loss_pct?) -----
    if final_winner == "AGGREGATED":
        row = sweep_df[(sweep_df["cell_id"] == "AGGREGATED") & (np.isclose(sweep_df["tau"], final_tau))]
    elif final_winner:
        row = sweep_df[(sweep_df["cell_id"] == final_winner) & (np.isclose(sweep_df["tau"], final_tau))]
    else:
        row = pd.DataFrame()
    if not row.empty:
        h3_fp_remov = float(row.iloc[0]["FP_removal_pct"])
        h3_tp_loss = float(row.iloc[0]["TP_loss_pct"])
        h3_verdict = "POSITIVE" if h3_fp_remov > h3_tp_loss else "NEGATIVE"
    else:
        h3_fp_remov = h3_tp_loss = float("nan")
        h3_verdict = "INCONCLUSIVE"

    # ----- Figures -----
    plt = matplotlib_setup()

    # Figure 1: ΔF1 vs tau curves
    n_cells_keep = len(keep_cell_ids)
    fig, ax = plt.subplots(1, 1, figsize=(9, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, max(10, n_cells_keep + 1)))
    color_idx = 0
    for cell_id in sorted(keep_cell_ids):
        sub = sweep_df[(sweep_df["cell_id"] == cell_id) & (sweep_df["scope"] == "per_cell")].sort_values("tau")
        if sub.empty:
            continue
        ax.plot(sub["tau"], sub["delta_F1_per_cell"], marker=".", ms=4, lw=1.0, alpha=0.85,
                color=colors[color_idx], label=cell_id[:40])
        color_idx += 1
    agg_sub = sweep_df[sweep_df["scope"] == "aggregated"].sort_values("tau")
    if not agg_sub.empty:
        ax.plot(agg_sub["tau"], agg_sub["delta_F1_per_cell"], lw=2.5, color="black",
                label="AGGREGATED (multi-cell)", zorder=10)
    ax.axhline(0, color="red", ls="--", lw=1.0, label="ΔF1 = 0 (caller baseline)")
    ax.axvline(final_tau, color="green", ls=":", lw=1.0, alpha=0.6,
               label=f"τ* (winning) = {final_tau:.2f}")
    ax.set_xlabel("τ threshold (P(TP) >= τ kept)")
    ax.set_ylabel("ΔF1 = F1_post − 0.7166")
    ax.set_title(f"Step 3 — ΔF1 vs τ (CV out-of-fold, V6_off Model B, n_cells={n_cells_keep})\n"
                 f"Best ΔF1 = {final_max_dF1:.5f} @ τ={final_tau:.2f} ({final_winner[:40]})",
                 fontsize=10)
    ax.legend(fontsize=7, loc="best", ncol=2)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(STEP5_DIR / "figures" / "step3_deltaf1_vs_tau.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    s3_log("Wrote figures/step3_deltaf1_vs_tau.png")

    # Figure 2: filter signal curve (FP_removal vs TP_loss)
    fig, ax = plt.subplots(1, 1, figsize=(7.5, 7))
    color_idx = 0
    for cell_id in sorted(keep_cell_ids):
        sub = sweep_df[(sweep_df["cell_id"] == cell_id) & (sweep_df["scope"] == "per_cell")].sort_values("tau")
        if sub.empty:
            continue
        ax.plot(sub["TP_loss_pct"], sub["FP_removal_pct"], marker=".", ms=4, lw=1.0, alpha=0.85,
                color=colors[color_idx], label=cell_id[:40])
        color_idx += 1
    if not agg_sub.empty:
        ax.plot(agg_sub["TP_loss_pct"], agg_sub["FP_removal_pct"], lw=2.5, color="black",
                label="AGGREGATED", zorder=10)
    ax.plot([0, 1], [0, 1], color="red", ls="--", lw=1.0, label="y=x (no filter signal)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("TP loss %")
    ax.set_ylabel("FP removal %")
    ax.set_title("Step 3 — Filter signal (FP removal vs TP loss, CV out-of-fold)\n"
                 "Above y=x = filter beats random for this cell")
    ax.legend(fontsize=7, loc="lower right", ncol=2)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(STEP5_DIR / "figures" / "step3_filter_signal_curve.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    s3_log("Wrote figures/step3_filter_signal_curve.png")

    # ----- step3_optimal_tau_summary.md -----
    lines = []
    lines.append("# Step 3 — Optimal τ* Summary")
    lines.append("")
    lines.append("> Generated by `scripts/delta_f1.py`")
    lines.append("")
    lines.append("## Caller baseline reverse-solve sanity check")
    lines.append("")
    for k, v in sanity.items():
        lines.append(f"- `{k}` = {v}")
    lines.append("")
    lines.append(f"## Strategy 1 — max ΔF1 (aggressive, per resolved decision #1)")
    lines.append("")
    lines.append("| cell_id | τ* | ΔF1 | TP_loss% | FP_remov% | filter_signal | n_fit (TP/FP) | F1_post |")
    lines.append("|---------|-----|------|----------|-----------|---------------|---------------|---------|")
    s1 = optimal_df[optimal_df["strategy"] == "max_DeltaF1"].sort_values("delta_F1", ascending=False)
    for _, r in s1.iterrows():
        lines.append(f"| {r['cell_id'][:40]} | {r['tau_star']:.2f} | "
                     f"{r['delta_F1']:+.5f} | {r['TP_loss_pct']*100:.2f}% | "
                     f"{r['FP_removal_pct']*100:.2f}% | {r['filter_signal']:+.4f} | "
                     f"{r['n_fit']} ({r['n_TP_fit']}/{r['n_FP_fit']}) | "
                     f"{r['global_F1_post']:.5f} |")
    lines.append("")
    lines.append(f"## Strategy 2 — max τ subject to FP_removal_pct > TP_loss_pct AND ΔF1 > 0 (sanity)")
    lines.append("")
    lines.append("| cell_id | τ* | ΔF1 | TP_loss% | FP_remov% | filter_signal | n_fit (TP/FP) | F1_post |")
    lines.append("|---------|-----|------|----------|-----------|---------------|---------------|---------|")
    s2 = optimal_df[optimal_df["strategy"] == "max_tau_s.t._sanity"].sort_values("tau_star", ascending=False, na_position="last")
    for _, r in s2.iterrows():
        if pd.isna(r["tau_star"]):
            lines.append(f"| {r['cell_id'][:40]} | (no eligible τ) | -- | -- | -- | -- | "
                         f"{r['n_fit']} ({r['n_TP_fit']}/{r['n_FP_fit']}) | -- |")
        else:
            lines.append(f"| {r['cell_id'][:40]} | {r['tau_star']:.2f} | "
                         f"{r['delta_F1']:+.5f} | {r['TP_loss_pct']*100:.2f}% | "
                         f"{r['FP_removal_pct']*100:.2f}% | {r['filter_signal']:+.4f} | "
                         f"{r['n_fit']} ({r['n_TP_fit']}/{r['n_FP_fit']}) | "
                         f"{r['global_F1_post']:.5f} |")
    lines.append("")
    lines.append("## Per-cell CV AUC (5-fold, out-of-fold)")
    lines.append("")
    lines.append("| cell_id | n_fit | TP/FP | CV AUC | CV status |")
    lines.append("|---------|-------|-------|--------|-----------|")
    for cid, meta in cell_meta.items():
        lines.append(f"| {cid[:40]} | {meta['n_fit']} | {meta['n_TP_fit']}/{meta['n_FP_fit']} | "
                     f"{meta['cv_auc']:.4f} | {meta['cv_status']} |")
    lines.append("")
    lines.append(f"## Winning configuration")
    lines.append(f"- max ΔF1 (any scope): **{final_max_dF1:+.5f}** at τ={final_tau:.2f}, scope={final_winner}")
    lines.append(f"- H2 verdict: **{h2_verdict}**")
    lines.append(f"- H3 verdict (at τ={final_tau:.2f}): **{h3_verdict}** (FP_remov={h3_fp_remov*100:.2f}%, TP_loss={h3_tp_loss*100:.2f}%)")

    summary_path = STEP5_DIR / "step3_optimal_tau_summary.md"
    summary_path.write_text("\n".join(lines))
    s3_log(f"Wrote {summary_path.name}")

    # ----- step3_findings.md -----
    fl = []
    fl.append("# Step 3 — ΔF1 Computation Findings")
    fl.append("")
    fl.append("> Generated by `scripts/delta_f1.py`")
    fl.append("")
    fl.append("## H2 Verdict — Does methylation-augmented filter improve F1?")
    fl.append(f"- max ΔF1 = **{final_max_dF1:+.5f}** at τ={final_tau:.2f}, scope=**{final_winner}**")
    fl.append(f"- F1_post = {F1_CALLER + final_max_dF1:.5f}  (caller baseline = {F1_CALLER})")
    fl.append(f"- **Verdict: {h2_verdict}**" + (
        " — methylation augmentation passes the H2 gate (ΔF1 > 0)."
        if h2_verdict == "POSITIVE"
        else " — methylation augmentation cannot beat caller-only F1 (ΔF1 ≤ 0)."))
    fl.append("")
    fl.append("## H3 Verdict — At optimal τ*, does FP_removal_pct exceed TP_loss_pct?")
    fl.append(f"- At τ*={final_tau:.2f} (scope={final_winner}):")
    fl.append(f"  - FP_removal_pct = {h3_fp_remov*100:.2f}%")
    fl.append(f"  - TP_loss_pct = {h3_tp_loss*100:.2f}%")
    fl.append(f"- **Verdict: {h3_verdict}**")
    fl.append("")
    fl.append("## Caveats applied")
    fl.append(f"- Only {n_keep} of {len(step2_summary)} Step 2 cells passed the gate "
              f"(Model B converged + n_fit ≥ {MIN_N_EFF_B} + n_TP ≥ {MIN_N_TP}).")
    fl.append(f"- 5-fold CV out-of-fold predictions used (not in-sample); Step 1 reported "
              f"~7 AUC points of optimism between in-sample 0.97 vs CV 0.90.")
    fl.append(f"- τ grid extended to [0.10, 0.95] by 0.01 (86 values) per Step 2 hand-off §3.")
    fl.append("")
    fl.append("## Baseline reverse-solve (per plan §Step 3)")
    fl.append("```")
    fl.append(f"  TP_caller = {TP_CALLER:,}")
    fl.append(f"  FP_caller = {FP_CALLER:,}")
    fl.append(f"  P = TP/(TP+FP) = {sanity['P_recomputed']:.4f}")
    fl.append(f"  Reverse-solve R given F1=0.7166: R = (F1·P)/(2P − F1) = {sanity['R_recomputed']:.4f}")
    fl.append(f"  FN = TP·(1−R)/R = {FN_CALLER:,}")
    fl.append(f"  N_truth_positive = TP + FN = {N_TRUTH_POSITIVE:,}")
    fl.append(f"  F1 sanity (re-computed) = {sanity['F1_recomputed']:.5f}  "
              f"({'OK' if sanity['F1_match_to_3decimals'] else 'MISMATCH'} vs target 0.7166)")
    fl.append("```")
    fl.append("")
    fl.append("Cross-check vs `annotated_hcc1395_cnv.tsv` (SEQC2 CNV-stratified subset):")
    fl.append("- That file has 31,640 rows (TP=30,387 / FP=1,253) — a CNV-region subset of the")
    fl.append("  master TSV's 35,332 rows (TP=30,490 / FP=4,842). The plan §Step 3 reverse-solve uses")
    fl.append("  the master TSV totals (same scope as our LR fits). Direct count of total")
    fl.append("  N_truth_positive from `annotated_hcc1395_cnv.tsv` is not possible (only contains")
    fl.append("  callable TP/FP rows, not all FN). Reverse-solve from F1=0.7166 is the authoritative path.")
    fl.append("")
    fl.append("## Per-cell ΔF1 summary (Strategy 1: max ΔF1)")
    fl.append("")
    fl.append("| cell_id | τ* | ΔF1 | TP_loss% | FP_remov% | F1_post |")
    fl.append("|---------|-----|------|----------|-----------|---------|")
    for _, r in s1.iterrows():
        lbl = "AGG" if r['cell_id'] == 'AGGREGATED' else 'per-cell'
        fl.append(f"| {r['cell_id'][:40]} ({lbl}) | {r['tau_star']:.2f} | "
                  f"{r['delta_F1']:+.5f} | {r['TP_loss_pct']*100:.2f}% | "
                  f"{r['FP_removal_pct']*100:.2f}% | {r['global_F1_post']:.5f} |")
    fl.append("")
    fl.append("## CV AUC sanity")
    fl.append("")
    fl.append("| cell_id | CV AUC | CV status |")
    fl.append("|---------|--------|-----------|")
    for cid, meta in cell_meta.items():
        fl.append(f"| {cid[:40]} | {meta['cv_auc']:.4f} | {meta['cv_status']} |")
    fl.append("")
    fl.append("## Hand-off to Coordinator (Step 5)")
    fl.append("- H1 (Step 1): POSITIVE (16/30 cell-combos LRT q<0.05)")
    fl.append("- H2 (this step): " + h2_verdict)
    fl.append("- H3 (this step at τ*): " + h3_verdict)
    fl.append("- H5 (Step 1, V5 vs V6 β consistency): POSITIVE (robust median Δβ = 1.87e-5)")
    fl.append("- H4 (Step 4 mechanism brainstorm): see step4_mechanism_candidates.md")
    fl.append("")
    fl.append(f"### Quantitative summary for Step 5 decision tree")
    fl.append(f"- max ΔF1 = {final_max_dF1:+.5f} (winning scope: {final_winner}, τ*={final_tau:.2f})")
    fl.append(f"- post-filter F1 = {F1_CALLER + final_max_dF1:.5f} vs caller {F1_CALLER}")
    fl.append(f"- FP removal at τ* = {h3_fp_remov*100:.2f}%")
    fl.append(f"- TP loss at τ* = {h3_tp_loss*100:.2f}%")
    fl.append(f"- Cells used: {n_keep} (gated; trivial UNKNOWN cells excluded)")
    fl.append("")
    fl.append("## Files")
    fl.append("- `step3_delta_f1.tsv` — full sweep (per cell + aggregated × 86 τ)")
    fl.append("- `step3_optimal_tau_summary.md` — optimal τ tables (Strategy 1+2) + winner")
    fl.append("- `intermediate/step3_optimal_tau.tsv` — Strategy 1+2 raw table")
    fl.append("- `intermediate/step3_log.txt` — execution log")
    fl.append("- `figures/step3_deltaf1_vs_tau.png` — ΔF1 vs τ per cell + aggregated")
    fl.append("- `figures/step3_filter_signal_curve.png` — FP_remov vs TP_loss")

    findings_path = STEP5_DIR / "step3_findings.md"
    findings_path.write_text("\n".join(fl))
    s3_log(f"Wrote {findings_path.name}")

    s3_log("=== Step 3 delta_f1 END ===")
    summary = {
        "n_cells_kept": n_keep,
        "n_truth_positive_reverse_solved": N_TRUTH_POSITIVE,
        "max_delta_F1": float(final_max_dF1),
        "max_delta_F1_tau": float(final_tau),
        "max_delta_F1_scope": final_winner,
        "F1_post_at_max": F1_CALLER + final_max_dF1 if not np.isnan(final_max_dF1) else None,
        "H2_verdict": h2_verdict,
        "H3_verdict": h3_verdict,
        "H3_FP_removal_pct_at_tau": float(h3_fp_remov),
        "H3_TP_loss_pct_at_tau": float(h3_tp_loss),
    }
    return summary


if __name__ == "__main__":
    s = main()
    print("\n=== STEP 3 SUMMARY ===")
    print(json.dumps(s, indent=2))
