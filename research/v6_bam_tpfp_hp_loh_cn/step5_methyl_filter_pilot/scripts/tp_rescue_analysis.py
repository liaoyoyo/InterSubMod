"""Step 5c — TP Rescue Analysis.

Goal: of the 35% TPs lost by the AGGREGATED τ*=0.52 filter, can we identify a
common feature signature to design a rescue rule that keeps them while not
re-introducing the removed FPs?

Pipeline (Stages 1-6 per plan):
  1. Re-run Step 3's CV out-of-fold P(TP) on the 4 AGGREGATED cells, build the
     aggregated pool with (region_id, y, p_oof, cell_id).
  2. Mann-Whitney U feature comparison lost_TP vs kept_TP + BH-FDR.
  3. Mann-Whitney U feature comparison lost_TP vs removed_FP + BH-FDR.
  4. Source taxonomy (S1-S4) classification.
  5. Rescue rule design + projected ΔF1.
  6. Falsifiability test (per-fold rescue) + cross-sample generalizability note.

Outputs:
  step5c_tp_rescue_analysis.md
  intermediate/step5c_lost_tp_features.tsv
  intermediate/step5c_mw_test_results.tsv
  intermediate/step5c_log.txt
  figures/step5c_lost_vs_kept_TP_violin.png
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold
from statsmodels.stats.multitest import multipletests

from _common_step1 import (
    MASTER_AUG,
    STEP5_DIR,
    annotate_v6off_cell,
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

LOG_PATH = STEP5_DIR / "intermediate" / "step5c_log.txt"

CV_FOLDS = 5
RANDOM_STATE = 42
TAU_STAR = 0.52
MIN_N_EFF_B = 100
MIN_N_TP = 5

# Caller baseline (mirror delta_f1)
TP_CALLER = 30490
FP_CALLER = 4842
FN_CALLER = 19288
F1_CALLER = 0.7166


def s5c_log(msg: str) -> None:
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
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


# -----------------------------------------------------------------------------
# CV out-of-fold (mirror delta_f1.cv_out_of_fold_predictions)
# -----------------------------------------------------------------------------
def cv_out_of_fold_predictions(df_cell: pd.DataFrame, bam: str, flag: str, model: str):
    X, _, mask = assemble_design(df_cell, bam, flag, model)
    y_full = (df_cell["label"] == "TP").astype(int).values
    y = y_full[mask]
    X_use = X[mask]
    n_use = mask.sum()
    if n_use < CV_FOLDS * 2 or len(np.unique(y)) < 2:
        return np.array([]), np.array([]), np.where(mask)[0], np.array([]), "insufficient_data"

    kf = KFold(n_splits=CV_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    p_oof = np.full(n_use, np.nan)
    fold_id = np.full(n_use, -1, dtype=int)

    n_fold_ok = 0
    for fid, (tr_idx, te_idx) in enumerate(kf.split(X_use)):
        y_tr = y[tr_idx]
        if len(np.unique(y_tr)) < 2:
            continue
        r, status = fit_lr(y_tr, X_use[tr_idx])
        if r is None:
            continue
        Xc_te = sm.add_constant(X_use[te_idx], has_constant="add")
        try:
            scores = r.predict(Xc_te)
            p_oof[te_idx] = scores
            fold_id[te_idx] = fid
            n_fold_ok += 1
        except Exception:
            continue

    # If any folds skipped, fill predictions with a fallback in-sample fit on all data
    if np.isnan(p_oof).any():
        r_all, _ = fit_lr(y, X_use)
        if r_all is not None:
            Xc_all = sm.add_constant(X_use, has_constant="add")
            try:
                p_all = r_all.predict(Xc_all)
                mask_nan = np.isnan(p_oof)
                p_oof[mask_nan] = p_all[mask_nan]
                fold_id[mask_nan] = -1  # mark fallback
            except Exception:
                pass
        if np.isnan(p_oof).any():
            keep_idx = np.where(~np.isnan(p_oof))[0]
            y = y[keep_idx]
            p_oof = p_oof[keep_idx]
            fold_id = fold_id[keep_idx]
            df_idx = np.where(mask)[0][keep_idx]
            return p_oof, y, df_idx, fold_id, f"ok_partial_n_fold_ok={n_fold_ok}"

    df_idx = np.where(mask)[0]
    return p_oof, y, df_idx, fold_id, f"ok_n_fold_ok={n_fold_ok}"


# -----------------------------------------------------------------------------
# Feature set for Mann-Whitney
# -----------------------------------------------------------------------------
FEATURE_COLS = [
    # caller / region
    "caller_af",
    "Coverage_Multiple",
    "AF_master",
    # HP family counts (V6_off)
    "V6_off_0",
    "V6_off_1",
    "V6_off_2",
    "V6_off_1-1",
    "V6_off_2-1",
    "V6_off_3",
    "V6_off_11",
    "V6_off_21",
    "V6_off_33",
    "V6_off_other",
    "V6_off_NG",
    "V6_off_n_reads",
    # Methylation V6_off
    "V6_off_meth_HPMergedDelta",
    "V6_off_meth_HPMergedP",
    "V6_off_meth_HPFineF",
    "V6_off_meth_HPFineP",
    "V6_off_meth_HPFineNGroups",
    "V6_off_meth_ClusterPermanovaF",
    "V6_off_meth_ClusterPermanovaP",
    "V6_off_meth_AlleleDelta",
    "V6_off_meth_AlleleP",
    "V6_off_meth_Entropy_Imbalance",
    "V6_off_meth_NME_HP1",
    "V6_off_meth_NME_HP2",
    "V6_off_meth_Epipoly_HP1",
    "V6_off_meth_Epipoly_HP2",
    "V6_off_meth_Epipoly_Delta",
    "V6_off_meth_NME_imbalance",
    "V6_off_meth_Epipoly_imbalance",
]


def mann_whitney_compare(df_a: pd.DataFrame, df_b: pd.DataFrame, label_a: str,
                         label_b: str, features: list[str]) -> pd.DataFrame:
    rows = []
    for col in features:
        if col not in df_a.columns or col not in df_b.columns:
            continue
        a_vals = pd.to_numeric(df_a[col], errors="coerce").dropna().values
        b_vals = pd.to_numeric(df_b[col], errors="coerce").dropna().values
        # NaN rate
        n_a = len(df_a)
        n_b = len(df_b)
        nan_rate_a = (n_a - len(a_vals)) / n_a if n_a > 0 else float("nan")
        nan_rate_b = (n_b - len(b_vals)) / n_b if n_b > 0 else float("nan")
        if len(a_vals) < 3 or len(b_vals) < 3:
            rows.append({
                "feature": col, "n_a": len(a_vals), "n_b": len(b_vals),
                "median_a": float("nan"), "median_b": float("nan"),
                "mean_a": float("nan"), "mean_b": float("nan"),
                "nan_rate_a": nan_rate_a, "nan_rate_b": nan_rate_b,
                "u_stat": float("nan"), "p_value": float("nan"),
                "direction": "insufficient",
                "label_a": label_a, "label_b": label_b,
            })
            continue
        # Mann-Whitney two-sided
        try:
            u, p = stats.mannwhitneyu(a_vals, b_vals, alternative="two-sided")
        except Exception:
            u, p = float("nan"), float("nan")
        med_a = float(np.median(a_vals))
        med_b = float(np.median(b_vals))
        if med_a > med_b:
            direction = f"{label_a}>{label_b}"
        elif med_a < med_b:
            direction = f"{label_a}<{label_b}"
        else:
            direction = "equal"
        rows.append({
            "feature": col, "n_a": len(a_vals), "n_b": len(b_vals),
            "median_a": med_a, "median_b": med_b,
            "mean_a": float(np.mean(a_vals)), "mean_b": float(np.mean(b_vals)),
            "nan_rate_a": nan_rate_a, "nan_rate_b": nan_rate_b,
            "u_stat": float(u), "p_value": float(p),
            "direction": direction,
            "label_a": label_a, "label_b": label_b,
        })
    out = pd.DataFrame(rows)
    valid = out["p_value"].notna().values
    q = np.full(len(out), float("nan"))
    if valid.any():
        _, q_valid, _, _ = multipletests(out.loc[valid, "p_value"].values, method="fdr_bh")
        q[valid] = q_valid
    out["q_value"] = q
    out["sig_q05"] = out["q_value"] < 0.05
    return out


def main():
    LOG_PATH.unlink(missing_ok=True)
    s5c_log("=== Step 5c TP Rescue Analysis START ===")

    # ----- Load + annotate master -----
    s5c_log(f"Loading {MASTER_AUG}")
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
    s5c_log(f"Master rows: {len(df):,}")

    # ----- Identify the 4 AGGREGATED cells -----
    step2_summary = pd.read_csv(STEP5_DIR / "intermediate" / "step2_cell_summary.tsv", sep="\t")
    step2_summary["passes_gate"] = (
        (step2_summary["model_used"] == "B")
        & (step2_summary["n_fit"] >= MIN_N_EFF_B)
        & (step2_summary["n_TP_fit"] >= MIN_N_TP)
    )
    keep_cell_ids = set(step2_summary.loc[step2_summary["passes_gate"], "cell_id"].tolist())
    s5c_log(f"4 AGGREGATED cells (Model B gate): {sorted(keep_cell_ids)}")

    fp_rich = find_fp_rich_cells()
    bam, flag = PRIMARY_BAM, PRIMARY_FLAG

    # ----- Stage 1: CV out-of-fold P(TP) per kept cell + group assignment -----
    s5c_log("")
    s5c_log("=== Stage 1: CV out-of-fold + group assignment ===")
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
        p_oof, y, df_indices, fold_id, status = cv_out_of_fold_predictions(df_cell, bam, flag, "B")
        if len(p_oof) == 0:
            s5c_log(f"  Skip {cell_id}: {status}")
            continue
        cv_auc_val = roc_auc_score(y, p_oof) if len(np.unique(y)) >= 2 else float("nan")
        s5c_log(f"  {cell_id}: n={len(y)}, TP/FP={int(y.sum())}/{int(len(y) - y.sum())}, "
                f"CV AUC = {cv_auc_val:.4f}")
        # df_indices are positions within df_cell; resolve to df rows
        df_cell_global_idx = df_cell.index.values[df_indices]
        rg_ids = df_cell["region_id"].values[df_indices]
        for gi, rg, y_val, p_val, fid in zip(df_cell_global_idx, rg_ids, y, p_oof, fold_id):
            if rg in seen_region_ids:
                continue
            seen_region_ids.add(rg)
            agg_records.append({
                "global_idx": int(gi),
                "region_id": rg,
                "y": int(y_val),
                "p_oof": float(p_val),
                "cell_id": cell_id,
                "fold_id": int(fid),
            })

    agg_df = pd.DataFrame(agg_records)
    s5c_log(f"Aggregated pool: {len(agg_df)} unique regions "
            f"(TP={int((agg_df['y']==1).sum())}, FP={int((agg_df['y']==0).sum())})")

    # Group assignment
    agg_df["kept"] = (agg_df["p_oof"] >= TAU_STAR)
    agg_df["group"] = np.where(
        (agg_df["y"] == 1) & (~agg_df["kept"]), "lost_TP",
        np.where(
            (agg_df["y"] == 1) & (agg_df["kept"]), "kept_TP",
            np.where(
                (agg_df["y"] == 0) & (~agg_df["kept"]), "removed_FP",
                "kept_FP"
            )
        )
    )

    group_sizes = agg_df["group"].value_counts().to_dict()
    s5c_log("Group sizes:")
    for g, n in group_sizes.items():
        s5c_log(f"  {g}: {n}")

    lost_tp = agg_df[agg_df["group"] == "lost_TP"]
    kept_tp = agg_df[agg_df["group"] == "kept_TP"]
    removed_fp = agg_df[agg_df["group"] == "removed_FP"]
    kept_fp = agg_df[agg_df["group"] == "kept_FP"]

    # Attach master features to each lost_TP row
    df_view = df.set_index(df.index)
    lost_tp_master = df.loc[lost_tp["global_idx"].values].copy()
    kept_tp_master = df.loc[kept_tp["global_idx"].values].copy()
    removed_fp_master = df.loc[removed_fp["global_idx"].values].copy()
    kept_fp_master = df.loc[kept_fp["global_idx"].values].copy()

    # Add p_oof + cell_id back
    for sub_master, sub_group in [
        (lost_tp_master, lost_tp),
        (kept_tp_master, kept_tp),
        (removed_fp_master, removed_fp),
        (kept_fp_master, kept_fp),
    ]:
        sub_master.reset_index(drop=False, inplace=True)
        sub_master.rename(columns={"index": "global_idx"}, inplace=True)
        sub_master["p_oof"] = sub_group["p_oof"].values
        sub_master["fold_id"] = sub_group["fold_id"].values
        sub_master["assigned_cell_id"] = sub_group["cell_id"].values

    # Write lost_tp features file
    lost_path = STEP5_DIR / "intermediate" / "step5c_lost_tp_features.tsv"
    keep_cols = ["region_id", "chr", "pos", "label", "assigned_cell_id", "p_oof", "fold_id",
                 "caller_af", "AF_master", "Coverage_Multiple", "loh_side_norm",
                 "hp_bucket_v6off", "cov_bin",
                 "V6_off_NG", "V6_off_n_reads"]
    keep_cols += [c for c in FEATURE_COLS if c not in keep_cols]
    keep_cols = [c for c in keep_cols if c in lost_tp_master.columns]
    lost_tp_master[keep_cols].to_csv(lost_path, sep="\t", index=False, float_format="%.6g")
    s5c_log(f"Wrote {lost_path.name} ({len(lost_tp_master)} rows)")

    # ----- Stage 2: Mann-Whitney lost_TP vs kept_TP -----
    s5c_log("")
    s5c_log("=== Stage 2: Mann-Whitney lost_TP vs kept_TP ===")
    mw_lost_vs_kept = mann_whitney_compare(lost_tp_master, kept_tp_master,
                                           "lost_TP", "kept_TP", FEATURE_COLS)
    sig_s2 = mw_lost_vs_kept[mw_lost_vs_kept["sig_q05"]].sort_values("q_value")
    s5c_log(f"q<0.05 features (lost vs kept TP): {len(sig_s2)}")
    for _, r in sig_s2.iterrows():
        s5c_log(f"  {r['feature']}: med_lost={r['median_a']:.4g}, med_kept={r['median_b']:.4g}, "
                f"p={r['p_value']:.3g}, q={r['q_value']:.3g}, dir={r['direction']}")

    # ----- Stage 3: Mann-Whitney lost_TP vs removed_FP -----
    s5c_log("")
    s5c_log("=== Stage 3: Mann-Whitney lost_TP vs removed_FP ===")
    mw_lost_vs_fp = mann_whitney_compare(lost_tp_master, removed_fp_master,
                                         "lost_TP", "removed_FP", FEATURE_COLS)
    sig_s3 = mw_lost_vs_fp[mw_lost_vs_fp["sig_q05"]].sort_values("q_value")
    s5c_log(f"q<0.05 features (lost TP vs removed FP): {len(sig_s3)}")
    for _, r in sig_s3.iterrows():
        s5c_log(f"  {r['feature']}: med_lost={r['median_a']:.4g}, med_FP={r['median_b']:.4g}, "
                f"p={r['p_value']:.3g}, q={r['q_value']:.3g}, dir={r['direction']}")

    # Combine MW results into one file
    mw_lost_vs_kept["comparison"] = "lost_TP_vs_kept_TP"
    mw_lost_vs_fp["comparison"] = "lost_TP_vs_removed_FP"
    mw_all = pd.concat([mw_lost_vs_kept, mw_lost_vs_fp], ignore_index=True)
    mw_path = STEP5_DIR / "intermediate" / "step5c_mw_test_results.tsv"
    mw_all.to_csv(mw_path, sep="\t", index=False, float_format="%.6g")
    s5c_log(f"Wrote {mw_path.name} ({len(mw_all)} rows)")

    # ----- Stage 4: Source taxonomy S1-S4 -----
    s5c_log("")
    s5c_log("=== Stage 4: Source taxonomy ===")

    # S1: caller_af edge - overlap with removed_FP
    lost_af = pd.to_numeric(lost_tp_master["caller_af"], errors="coerce").dropna().values
    fp_af = pd.to_numeric(removed_fp_master["caller_af"], errors="coerce").dropna().values
    # Use empirical overlap: fraction of lost_TP caller_af in [min_FP, max_FP] using
    # a quantile-based interval (5th-95th of FP) — more meaningful than full range
    fp_p05, fp_p95 = (np.percentile(fp_af, [5, 95]) if len(fp_af) > 0 else (0.0, 1.0))
    n_lost_in_fp_range = int(((lost_af >= fp_p05) & (lost_af <= fp_p95)).sum())
    s1_frac = n_lost_in_fp_range / len(lost_af) if len(lost_af) > 0 else 0.0
    s5c_log(f"S1 caller_af edge: lost_TP in FP[p05={fp_p05:.3f}, p95={fp_p95:.3f}] = "
            f"{n_lost_in_fp_range}/{len(lost_af)} = {s1_frac*100:.1f}%")

    # S2: methylation NaN rate (NME/Epipoly)
    s2_meth_cols = ["V6_off_meth_NME_imbalance", "V6_off_meth_Epipoly_Delta",
                    "V6_off_meth_HPFineF"]
    s2_lost_nan = float(lost_tp_master[s2_meth_cols].isna().any(axis=1).mean())
    s2_kept_nan = float(kept_tp_master[s2_meth_cols].isna().any(axis=1).mean())
    s2_diff = s2_lost_nan - s2_kept_nan
    s5c_log(f"S2 methyl-NaN: lost_TP any-NaN rate={s2_lost_nan*100:.1f}%, "
            f"kept_TP={s2_kept_nan*100:.1f}%, diff={s2_diff*100:+.1f}pp")

    # S3: heterogeneous sub-clone — caller_af systematically low (< 0.3 vs kept)
    lost_low_af_frac = float(np.mean(lost_af < 0.3)) if len(lost_af) > 0 else float("nan")
    kept_af = pd.to_numeric(kept_tp_master["caller_af"], errors="coerce").dropna().values
    kept_low_af_frac = float(np.mean(kept_af < 0.3)) if len(kept_af) > 0 else float("nan")
    s5c_log(f"S3 low-AF sub-clone: lost_TP frac caller_af<0.3 = {lost_low_af_frac*100:.1f}%, "
            f"kept_TP={kept_low_af_frac*100:.1f}%")

    # S4: chr enrichment — Fisher exact chr8 vs other in lost_TP vs kept_TP
    lost_chr8 = int((lost_tp_master["chr"] == "chr8").sum())
    lost_other = len(lost_tp_master) - lost_chr8
    kept_chr8 = int((kept_tp_master["chr"] == "chr8").sum())
    kept_other = len(kept_tp_master) - kept_chr8
    contingency = np.array([[lost_chr8, lost_other],
                            [kept_chr8, kept_other]])
    try:
        odds, p_fisher = stats.fisher_exact(contingency)
    except Exception:
        odds, p_fisher = float("nan"), float("nan")
    lost_chr8_frac = lost_chr8 / len(lost_tp_master) if len(lost_tp_master) > 0 else 0.0
    kept_chr8_frac = kept_chr8 / len(kept_tp_master) if len(kept_tp_master) > 0 else 0.0
    s5c_log(f"S4 chr8 enrichment in lost_TP: lost={lost_chr8}/{len(lost_tp_master)} ({lost_chr8_frac*100:.1f}%), "
            f"kept={kept_chr8}/{len(kept_tp_master)} ({kept_chr8_frac*100:.1f}%), "
            f"odds={odds:.3g}, p_fisher={p_fisher:.3g}")

    # Decide dominant source (≥50% threshold on its explained fraction)
    sources_explained = {
        "S1_caller_af_edge": s1_frac,
        "S2_methyl_NaN": s2_lost_nan,  # absolute (not differential)
        "S3_low_AF_subclone": lost_low_af_frac,
        "S4_chr8_enrichment": lost_chr8_frac,
    }
    dom_source = max(sources_explained, key=sources_explained.get)
    dom_frac = sources_explained[dom_source]
    s5c_log(f"Dominant source: {dom_source} ({dom_frac*100:.1f}% of lost_TP)")
    if dom_frac < 0.5:
        s5c_log(f"  ⚠ Below 50% threshold — no single source dominates")

    # ----- Stage 5: Rescue rule design + projected ΔF1 -----
    s5c_log("")
    s5c_log("=== Stage 5: Rescue rule + projected ΔF1 ===")

    # Strategy: pick the top distinguishing feature(s) from Stage 3 (lost vs FP)
    # that have q<0.05 AND median direction lost > FP (or < FP).
    # Then build candidate rules using that feature's percentile threshold.
    rescue_rule_candidates = []
    # Restrict to features where lost_TP and removed_FP differ AND lost_TP has signal
    # vs kept_TP is small (so we don't accidentally undo the filter for kept FP).
    if len(sig_s3) > 0:
        # take top 3 by q-value
        top_feats = sig_s3.head(3)
        for _, feat_row in top_feats.iterrows():
            feat = feat_row["feature"]
            lost_vals = pd.to_numeric(lost_tp_master[feat], errors="coerce").dropna().values
            fp_vals = pd.to_numeric(removed_fp_master[feat], errors="coerce").dropna().values
            kept_fp_vals = pd.to_numeric(kept_fp_master[feat], errors="coerce").dropna().values \
                if len(kept_fp_master) > 0 else np.array([])
            # Heuristic threshold: midpoint of lost median and FP p75 (or p25)
            direction = feat_row["direction"]
            if "lost_TP>" in direction:
                # Rescue if value >= threshold; threshold = FP p75 (above which FPs rare)
                thr = float(np.percentile(fp_vals, 75)) if len(fp_vals) > 0 else float("nan")
                if not np.isnan(thr):
                    rescue_recall = float(np.mean(lost_vals >= thr))  # fraction of lost_TP rescued
                    fp_hits = int(np.sum(fp_vals >= thr))  # FPs that would be re-introduced
                    rule = f"{feat} >= {thr:.4g}"
                else:
                    rescue_recall = 0.0
                    fp_hits = 0
                    rule = f"{feat} >= NaN"
            else:  # lost_TP <
                thr = float(np.percentile(fp_vals, 25)) if len(fp_vals) > 0 else float("nan")
                if not np.isnan(thr):
                    rescue_recall = float(np.mean(lost_vals <= thr))
                    fp_hits = int(np.sum(fp_vals <= thr))
                    rule = f"{feat} <= {thr:.4g}"
                else:
                    rescue_recall = 0.0
                    fp_hits = 0
                    rule = f"{feat} <= NaN"
            rescue_tp_count = round(rescue_recall * len(lost_tp_master))
            rescue_rule_candidates.append({
                "feature": feat,
                "direction": direction,
                "threshold": thr,
                "rule": rule,
                "rescue_recall_on_lost_TP": rescue_recall,
                "rescue_TP_count": int(rescue_tp_count),
                "FP_reintroduced_count": fp_hits,
                "kept_FP_carry": int(np.sum(kept_fp_vals >= thr)) if "lost_TP>" in direction
                                  else int(np.sum(kept_fp_vals <= thr)),
            })

    # Compute projected ΔF1 for each candidate rule
    # Start from Step 3 post-filter state:
    #   At τ*=0.52 AGGREGATED:
    #     TP_post = TP_caller - 21 (21 lost TPs from agg)
    #     FP_post = FP_caller - removed_FP_count
    #     FN_post = FN_caller + 21
    n_lost_tp_actual = len(lost_tp_master)
    n_removed_fp_actual = len(removed_fp_master)
    tp_post_step3 = TP_CALLER - n_lost_tp_actual
    fp_post_step3 = FP_CALLER - n_removed_fp_actual
    fn_post_step3 = FN_CALLER + n_lost_tp_actual
    f1_post_step3 = f1_from_counts(tp_post_step3, fp_post_step3, fn_post_step3)
    s5c_log(f"Step 3 baseline (τ*=0.52 AGG): TP={tp_post_step3}, FP={fp_post_step3}, "
            f"FN={fn_post_step3}, F1={f1_post_step3:.5f}")
    s5c_log(f"  caller F1 = {F1_CALLER}, current ΔF1 = {f1_post_step3 - F1_CALLER:+.5f}")

    for c in rescue_rule_candidates:
        tp_rescue = c["rescue_TP_count"]
        fp_reintro = c["FP_reintroduced_count"]
        tp_new = tp_post_step3 + tp_rescue
        fp_new = fp_post_step3 + fp_reintro
        fn_new = fn_post_step3 - tp_rescue
        f1_new = f1_from_counts(tp_new, fp_new, fn_new)
        c["F1_after_rescue"] = f1_new
        c["delta_F1_vs_caller"] = f1_new - F1_CALLER
        c["delta_F1_vs_step3"] = f1_new - f1_post_step3
        s5c_log(f"Rule '{c['rule']}': TP+={tp_rescue}, FP+={fp_reintro} -> "
                f"F1={f1_new:.5f}, ΔF1_vs_caller={c['delta_F1_vs_caller']:+.5f}, "
                f"ΔF1_vs_step3={c['delta_F1_vs_step3']:+.5f}")

    # Pick best rule (max delta_F1_vs_step3)
    if rescue_rule_candidates:
        best_rule = max(rescue_rule_candidates,
                        key=lambda r: r.get("delta_F1_vs_step3", -1))
    else:
        best_rule = None

    # ----- Stage 6: Falsifiability — per-fold validation -----
    s5c_log("")
    s5c_log("=== Stage 6: Falsifiability — per-fold rescue test ===")
    fold_results = []
    if best_rule is not None and not np.isnan(best_rule["threshold"]):
        thr = best_rule["threshold"]
        feat = best_rule["feature"]
        direction = best_rule["direction"]
        # For each fold, compute rescue recall and FP-reintroduce count on
        # that fold's held-out lost_TP / removed_FP.
        for fid in sorted(lost_tp_master["fold_id"].unique()):
            if fid == -1:
                continue
            lost_f = lost_tp_master[lost_tp_master["fold_id"] == fid]
            fp_f = removed_fp_master[removed_fp_master["fold_id"] == fid]
            lv = pd.to_numeric(lost_f[feat], errors="coerce").dropna().values
            fv = pd.to_numeric(fp_f[feat], errors="coerce").dropna().values
            if "lost_TP>" in direction:
                recall = float(np.mean(lv >= thr)) if len(lv) > 0 else float("nan")
                fp_hit = int(np.sum(fv >= thr)) if len(fv) > 0 else 0
            else:
                recall = float(np.mean(lv <= thr)) if len(lv) > 0 else float("nan")
                fp_hit = int(np.sum(fv <= thr)) if len(fv) > 0 else 0
            fold_results.append({
                "fold_id": int(fid),
                "n_lost_TP_fold": len(lost_f),
                "n_removed_FP_fold": len(fp_f),
                "rescue_recall_fold": recall,
                "FP_reintroduced_fold": fp_hit,
            })
            s5c_log(f"  Fold {fid}: lost={len(lost_f)}, FP={len(fp_f)}, "
                    f"recall={recall:.3f}, FP_reintro={fp_hit}")

        # Variability across folds
        recalls = [r["rescue_recall_fold"] for r in fold_results
                   if not np.isnan(r["rescue_recall_fold"])]
        recall_cv = float(np.std(recalls) / np.mean(recalls)) if recalls and np.mean(recalls) > 0 else float("nan")
        s5c_log(f"  Across folds: mean recall={np.mean(recalls):.3f}, "
                f"std={np.std(recalls):.3f}, CV={recall_cv:.3f}")

    # ----- Generate markdown report -----
    md = []
    md.append("# Step 5c — TP Rescue Analysis")
    md.append("")
    md.append("> Generated by `scripts/tp_rescue_analysis.py`")
    md.append("")
    md.append("## Question")
    md.append("")
    md.append("Step 3 ΔF1 = +0.00242 @ τ*=0.52 (AGGREGATED 4 cells) but TP_loss=35% "
              "(FP_removal=98.26%). Can we identify a common feature signature among "
              "the lost TPs to rescue them without re-introducing the removed FPs?")
    md.append("")
    md.append("## Stage 1 — 4-group assignment")
    md.append("")
    md.append("4 AGGREGATED cells (Model B converged + n_fit≥100 + n_TP≥5):")
    md.append("")
    for cid in sorted(keep_cell_ids):
        md.append(f"- `{cid}`")
    md.append("")
    md.append(f"Total unique regions in aggregated pool: **{len(agg_df)}** "
              f"(TP={int((agg_df['y']==1).sum())}, FP={int((agg_df['y']==0).sum())})")
    md.append("")
    md.append("τ* = **0.52** (Step 3 AGGREGATED winning threshold)")
    md.append("")
    md.append("| group | count | description |")
    md.append("|-------|-------|-------------|")
    md.append(f"| lost_TP | {group_sizes.get('lost_TP', 0)} | truth=TP AND P(TP) < 0.52 |")
    md.append(f"| kept_TP | {group_sizes.get('kept_TP', 0)} | truth=TP AND P(TP) >= 0.52 |")
    md.append(f"| removed_FP | {group_sizes.get('removed_FP', 0)} | truth=FP AND P(TP) < 0.52 |")
    md.append(f"| kept_FP | {group_sizes.get('kept_FP', 0)} | truth=FP AND P(TP) >= 0.52 |")
    md.append("")
    md.append(f"Expected (from Step 3): lost_TP ~21, kept_TP ~39, removed_FP ~338, kept_FP ~6")
    md.append("")
    md.append("## Stage 2 — lost_TP vs kept_TP (Mann-Whitney + BH-FDR)")
    md.append("")
    md.append(f"Tested {len(mw_lost_vs_kept)} features. Significant at q<0.05: "
              f"**{len(sig_s2)}**")
    md.append("")
    if len(sig_s2) > 0:
        md.append("| feature | n_lost | n_kept | median_lost | median_kept | p | q | direction |")
        md.append("|---------|--------|--------|-------------|-------------|---|---|-----------|")
        for _, r in sig_s2.head(15).iterrows():
            md.append(f"| {r['feature']} | {int(r['n_a'])} | {int(r['n_b'])} | "
                      f"{r['median_a']:.4g} | {r['median_b']:.4g} | "
                      f"{r['p_value']:.3g} | {r['q_value']:.3g} | {r['direction']} |")
    else:
        md.append("No features distinguish lost_TP from kept_TP at q<0.05.")
        md.append("")
        md.append("Top 5 by raw p-value (for inspection only):")
        md.append("")
        md.append("| feature | n_lost | n_kept | median_lost | median_kept | p | q | direction |")
        md.append("|---------|--------|--------|-------------|-------------|---|---|-----------|")
        for _, r in mw_lost_vs_kept.sort_values("p_value").head(5).iterrows():
            md.append(f"| {r['feature']} | {int(r['n_a'])} | {int(r['n_b'])} | "
                      f"{r['median_a']:.4g} | {r['median_b']:.4g} | "
                      f"{r['p_value']:.3g} | {r['q_value']:.3g} | {r['direction']} |")
    md.append("")
    md.append("## Stage 3 — lost_TP vs removed_FP (Mann-Whitney + BH-FDR)")
    md.append("")
    md.append(f"Significant at q<0.05: **{len(sig_s3)}**")
    md.append("")
    if len(sig_s3) > 0:
        md.append("| feature | n_lost | n_FP | median_lost | median_FP | p | q | direction |")
        md.append("|---------|--------|------|-------------|-----------|---|---|-----------|")
        for _, r in sig_s3.head(15).iterrows():
            md.append(f"| {r['feature']} | {int(r['n_a'])} | {int(r['n_b'])} | "
                      f"{r['median_a']:.4g} | {r['median_b']:.4g} | "
                      f"{r['p_value']:.3g} | {r['q_value']:.3g} | {r['direction']} |")
    else:
        md.append("No features distinguish lost_TP from removed_FP at q<0.05.")
        md.append("")
        md.append("Top 5 by raw p-value:")
        md.append("")
        md.append("| feature | n_lost | n_FP | median_lost | median_FP | p | q | direction |")
        md.append("|---------|--------|------|-------------|-----------|---|---|-----------|")
        for _, r in mw_lost_vs_fp.sort_values("p_value").head(5).iterrows():
            md.append(f"| {r['feature']} | {int(r['n_a'])} | {int(r['n_b'])} | "
                      f"{r['median_a']:.4g} | {r['median_b']:.4g} | "
                      f"{r['p_value']:.3g} | {r['q_value']:.3g} | {r['direction']} |")
    md.append("")
    md.append("## Stage 4 — Source taxonomy (S1-S4)")
    md.append("")
    md.append("| source | metric | value | interpretation |")
    md.append("|--------|--------|-------|----------------|")
    md.append(f"| S1 caller_af edge | fraction of lost_TP in FP[p05={fp_p05:.3f}, p95={fp_p95:.3f}] | "
              f"{s1_frac*100:.1f}% | "
              f"{'⚠ majority overlap with FP AF range' if s1_frac > 0.5 else 'partial overlap'} |")
    md.append(f"| S2 methyl NaN | lost_TP any-NaN(NME/Epipoly/HPFineF) | "
              f"{s2_lost_nan*100:.1f}% (vs kept {s2_kept_nan*100:.1f}%) | "
              f"{'⚠ NaN-driven loss' if s2_diff > 0.2 else 'NaN rate similar across groups'} |")
    md.append(f"| S3 low-AF sub-clone | lost_TP frac caller_af<0.3 | "
              f"{lost_low_af_frac*100:.1f}% (vs kept {kept_low_af_frac*100:.1f}%) | "
              f"{'⚠ low-AF dominant' if lost_low_af_frac > 0.5 else 'AF spread similar'} |")
    md.append(f"| S4 chr8 enrichment | lost_TP frac chr8 | "
              f"{lost_chr8_frac*100:.1f}% (vs kept {kept_chr8_frac*100:.1f}%) | "
              f"Fisher p={p_fisher:.3g}, odds={odds:.3g} |")
    md.append("")
    md.append(f"**Dominant source**: `{dom_source}` ({dom_frac*100:.1f}% of lost_TP)")
    if dom_frac < 0.5:
        md.append("")
        md.append("⚠ No single source explains ≥50% of lost_TP — heterogeneous origin.")
    md.append("")
    md.append("## Stage 5 — Rescue rule design + projected ΔF1")
    md.append("")
    md.append(f"Step 3 baseline state (τ*=0.52 AGG applied):")
    md.append("")
    md.append(f"- TP_post = {tp_post_step3}")
    md.append(f"- FP_post = {fp_post_step3}")
    md.append(f"- FN_post = {fn_post_step3}")
    md.append(f"- F1_post = {f1_post_step3:.5f}")
    md.append(f"- ΔF1 vs caller = {f1_post_step3 - F1_CALLER:+.5f}")
    md.append("")
    if rescue_rule_candidates:
        md.append("### Candidate rescue rules")
        md.append("")
        md.append("Rule form: `if P(TP)<0.52 AND <rescue condition> then KEEP`")
        md.append("")
        md.append("| feature | rule | rescue_recall on lost_TP | rescue_TP | FP_reintroduced | "
                  "F1_after_rescue | ΔF1 vs caller | ΔF1 vs Step3 |")
        md.append("|---------|------|------------|-----------|-----------------|------------------|---------------|---------------|")
        for c in rescue_rule_candidates:
            md.append(f"| {c['feature']} | `{c['rule']}` | "
                      f"{c['rescue_recall_on_lost_TP']*100:.1f}% | "
                      f"{c['rescue_TP_count']} | {c['FP_reintroduced_count']} | "
                      f"{c['F1_after_rescue']:.5f} | "
                      f"{c['delta_F1_vs_caller']:+.5f} | "
                      f"{c['delta_F1_vs_step3']:+.5f} |")
        md.append("")
        if best_rule is not None:
            md.append(f"**Best rule (max ΔF1 vs Step 3)**: `{best_rule['rule']}` "
                      f"→ ΔF1 vs Step3 = {best_rule.get('delta_F1_vs_step3', float('nan')):+.5f}")
    else:
        md.append("⚠ No Stage 3 q<0.05 features available — no statistically supported rescue "
                  "rule can be derived.")
    md.append("")
    md.append("## Stage 6 — Falsifiability & limitations")
    md.append("")
    if fold_results:
        md.append("### Per-fold rescue rule stability")
        md.append("")
        md.append(f"Best rule: `{best_rule['rule']}` evaluated per CV fold.")
        md.append("")
        md.append("| fold | n_lost_TP | n_removed_FP | rescue_recall | FP_reintroduced |")
        md.append("|------|-----------|--------------|---------------|-----------------|")
        for r in fold_results:
            md.append(f"| {r['fold_id']} | {r['n_lost_TP_fold']} | {r['n_removed_FP_fold']} | "
                      f"{r['rescue_recall_fold']:.3f} | {r['FP_reintroduced_fold']} |")
        recalls_v = [r["rescue_recall_fold"] for r in fold_results
                     if not np.isnan(r["rescue_recall_fold"])]
        if recalls_v:
            md.append("")
            md.append(f"- Mean recall = {np.mean(recalls_v):.3f}, "
                      f"std = {np.std(recalls_v):.3f}, "
                      f"CV = {np.std(recalls_v)/np.mean(recalls_v):.3f}" if np.mean(recalls_v) > 0 else "")
    else:
        md.append("Per-fold falsifiability not computed (no valid rescue rule).")
    md.append("")
    md.append("### Cross-sample generalizability")
    md.append("")
    md.append("- This analysis is single-sample (HCC1395). InterSubMod history shows "
              "single-sample → cross-sample collapse rate ~50-70% (see MEMORY.md "
              "`feedback_evidence_driven_iteration_workflow.md` and prior NEGATIVE "
              "cycles where single-sample AUC>0.6 dropped to <0.55 across 7 samples).")
    md.append("- Recommendation: any rescue rule passing single-sample falsifiability "
              "should be re-tested on COLO829 / H2009 / HCC1954 before claiming a "
              "production filter.")
    md.append("")
    md.append("## Limitations")
    md.append("")
    md.append(f"- lost_TP n = {len(lost_tp_master)} — very small; Mann-Whitney has limited power.")
    md.append(f"- removed_FP n = {len(removed_fp_master)} — substantially larger; "
              f"asymmetric power may inflate false positives in MW tests despite BH-FDR.")
    md.append("- Rescue rule thresholds derived in-fold from the same data they are evaluated "
              "on — per-fold validation in Stage 6 attempts to mitigate but is not a true "
              "held-out test.")
    md.append("- 4 AGGREGATED cells represent ~1.1% of master rows (404/35,332). "
              "Generalization to whole-genome filter is not assured.")
    md.append("")
    md.append("## Files")
    md.append("")
    md.append("- `intermediate/step5c_lost_tp_features.tsv` — per-row features for lost_TP set")
    md.append("- `intermediate/step5c_mw_test_results.tsv` — Mann-Whitney results (both comparisons)")
    md.append("- `intermediate/step5c_log.txt` — execution log")
    md.append("- `figures/step5c_lost_vs_kept_TP_violin.png` — top 5 distinguishing features violin")
    md.append("")

    # ----- Violin plot for top-5 features (lost vs kept TP) -----
    plt = matplotlib_setup()
    # Pick top 5 by p-value (regardless of q significance, so we always have 5)
    top5 = mw_lost_vs_kept.sort_values("p_value").head(5)
    if len(top5) >= 1:
        ncol = min(5, len(top5))
        fig, axes = plt.subplots(1, ncol, figsize=(ncol * 3.0, 4.5))
        if ncol == 1:
            axes = [axes]
        for i, (_, feat_row) in enumerate(top5.iterrows()):
            feat = feat_row["feature"]
            lost_vals = pd.to_numeric(lost_tp_master[feat], errors="coerce").dropna().values
            kept_vals = pd.to_numeric(kept_tp_master[feat], errors="coerce").dropna().values
            ax = axes[i]
            data = [lost_vals, kept_vals]
            positions = [1, 2]
            if all(len(d) >= 2 for d in data):
                parts = ax.violinplot(data, positions=positions, showmeans=False,
                                      showmedians=True, widths=0.7)
                for j, pc in enumerate(parts['bodies']):
                    pc.set_facecolor(['#d62728', '#1f77b4'][j])
                    pc.set_alpha(0.6)
            ax.scatter(np.full(len(lost_vals), 1) + np.random.uniform(-0.05, 0.05, len(lost_vals)),
                       lost_vals, color='#d62728', s=8, alpha=0.5)
            ax.scatter(np.full(len(kept_vals), 2) + np.random.uniform(-0.05, 0.05, len(kept_vals)),
                       kept_vals, color='#1f77b4', s=8, alpha=0.5)
            ax.set_xticks(positions)
            ax.set_xticklabels([f"lost\n(n={len(lost_vals)})", f"kept\n(n={len(kept_vals)})"],
                               fontsize=8)
            ax.set_title(f"{feat[:30]}\np={feat_row['p_value']:.3g} q={feat_row['q_value']:.3g}",
                         fontsize=8)
            ax.grid(True, alpha=0.3)
        fig.suptitle("Step 5c — Top 5 features distinguishing lost_TP vs kept_TP\n"
                     "(at τ*=0.52, AGGREGATED 4 cells, CV out-of-fold P(TP))",
                     fontsize=10)
        plt.tight_layout()
        fig_path = STEP5_DIR / "figures" / "step5c_lost_vs_kept_TP_violin.png"
        fig.savefig(fig_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        s5c_log(f"Wrote {fig_path.name}")

    out_md = STEP5_DIR / "step5c_tp_rescue_analysis.md"
    out_md.write_text("\n".join(md))
    s5c_log(f"Wrote {out_md.name}")

    summary = {
        "lost_TP": int(group_sizes.get("lost_TP", 0)),
        "kept_TP": int(group_sizes.get("kept_TP", 0)),
        "removed_FP": int(group_sizes.get("removed_FP", 0)),
        "kept_FP": int(group_sizes.get("kept_FP", 0)),
        "n_sig_lost_vs_kept_TP": int(len(sig_s2)),
        "n_sig_lost_vs_removed_FP": int(len(sig_s3)),
        "dominant_source": dom_source,
        "dominant_source_frac": float(dom_frac),
        "best_rescue_rule": best_rule["rule"] if best_rule else None,
        "best_rescue_delta_F1_vs_step3": (best_rule.get("delta_F1_vs_step3")
                                          if best_rule else None),
        "best_rescue_delta_F1_vs_caller": (best_rule.get("delta_F1_vs_caller")
                                           if best_rule else None),
        "n_folds_validated": len(fold_results),
    }
    s5c_log("=== Step 5c END ===")
    return summary


if __name__ == "__main__":
    s = main()
    print("\n=== STEP 5c SUMMARY ===")
    print(json.dumps(s, indent=2, default=str))
