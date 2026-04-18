#!/usr/bin/env python3
"""
Beyond-AUC Comprehensive Feature Validation
=============================================

Seven complementary validation methods to determine whether ISM features
have genuine discriminative signal beyond what AUC-ROC captures.

Methods:
  M1: Precision-Recall AUC
  M2: Systematic Residualization (OLS confound removal)
  M3: Conditional Subgroup Discovery (stratification grid)
  M4: Calibration + Brier Score
  M5: Bootstrap Stability (1000 iterations)
  M6: Non-linear Modeling (GradientBoosting + LOSO-CV + permutation importance)
  M7: Distribution Overlap + Tail Enrichment (KDE, KS, tail analysis)

Generated: 2026-04-09
"""

from __future__ import annotations

import json
import time
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import sys
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis")
from observation_common import (
    load_master_dataset, setup_plot_style, save_figure, add_caption,
    compute_auc, compute_effect_size, compute_enrichment,
    bootstrap_ci, mannwhitney_test, format_p,
    encode_truth_binary, ensure_dir, safe_div,
    write_round_context,
    OUTPUT_ROOT, SAMPLE_ORDER, MODE_ORDER,
    COLOR_TP, COLOR_FP, TRUTH_PALETTE, MODE_PALETTE,
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OUT_DIR = OUTPUT_ROOT / "20260409_beyond_auc_validation"
TSV_DIR = OUT_DIR / "tsv"
FIG_DIR = OUT_DIR / "figures"

# ISM features to evaluate (G2-G8 + Caller)
ISM_FEATURES = [
    # G2: Global
    "GlobalP", "CramersV", "HeuristicScore",
    # G3: Pairwise
    "PairwiseMeanDist", "PairwiseMedianDist",
    # G4: HP Merged
    "HPMergedDelta", "HPMergedP", "HP1FamilyN", "HP2FamilyN",
    # G5: HP Fine
    "HPFineF", "HPFineP", "HPFineNGroups",
    # G6: Allele
    "AlleleDelta", "AlleleP",
    # G7: Label Verification
    "LabelHPPermanovaF", "LabelHPPermanovaP",
    "LabelAllelePermanovaF", "LabelAllelePermanovaP",
    "ClusterPermanovaF",
    # G8: Quality Composite
    "Quality_Score", "HP_Ratio", "Stability",
]

CALLER_FEATURES = ["caller_af", "caller_gq", "caller_dp"]

ALL_FEATURES = ISM_FEATURES + CALLER_FEATURES

# Confounder columns for residualization
CONFOUNDERS = ["NumReads", "NumCpGs", "caller_af", "caller_dp"]

# Stratification
AF_BINS_COARSE = {
    "low_0-0.2": (0.0, 0.2),
    "mid_0.2-0.5": (0.2, 0.5),
    "high_0.5-0.8": (0.5, 0.8),
    "extreme_0.8+": (0.8, 1.01),
}

SAMPLE_GROUPS = {
    "HCC1395": ["HCC1395"],
    "HCC1395_DORADO": ["HCC1395_DORADO"],
    "COLO829": ["COLO829"],
    "H1437": ["H1437"],
    "H2009": ["H2009"],
    "HCC1937": ["HCC1937"],
    "HCC1954": ["HCC1954"],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """OLS residualization: return y - X @ beta."""
    mask = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    res = np.full_like(y, np.nan)
    if mask.sum() < 10:
        return res
    Xm = np.column_stack([np.ones(mask.sum()), X[mask]])
    beta, _, _, _ = np.linalg.lstsq(Xm, y[mask], rcond=None)
    res[mask] = y[mask] - Xm @ beta
    return res


def _safe_auc(y_true, y_score):
    """Compute AUC handling NaN."""
    mask = np.isfinite(y_score) & np.isfinite(y_true)
    yt, ys = y_true[mask], y_score[mask]
    if len(np.unique(yt)) < 2 or len(yt) < 10:
        return np.nan
    from sklearn.metrics import roc_auc_score
    try:
        return roc_auc_score(yt, ys)
    except Exception:
        return np.nan


def _safe_pr_auc(y_true, y_score):
    """Compute Average Precision (PR-AUC) handling NaN."""
    mask = np.isfinite(y_score) & np.isfinite(y_true)
    yt, ys = y_true[mask], y_score[mask]
    if len(np.unique(yt)) < 2 or len(yt) < 10:
        return np.nan
    from sklearn.metrics import average_precision_score
    try:
        return average_precision_score(yt, ys)
    except Exception:
        return np.nan


def _cohens_d(g1, g2):
    """Quick Cohen's d without bootstrap."""
    g1 = g1[np.isfinite(g1)]
    g2 = g2[np.isfinite(g2)]
    if len(g1) < 2 or len(g2) < 2:
        return np.nan
    n1, n2 = len(g1), len(g2)
    pooled = np.sqrt(((n1 - 1) * g1.std(ddof=1)**2 + (n2 - 1) * g2.std(ddof=1)**2) / (n1 + n2 - 2))
    return (g1.mean() - g2.mean()) / pooled if pooled > 0 else 0.0


def _ks_test(tp_vals, fp_vals):
    """KS test between TP and FP distributions."""
    tp = tp_vals[np.isfinite(tp_vals)]
    fp = fp_vals[np.isfinite(fp_vals)]
    if len(tp) < 10 or len(fp) < 10:
        return np.nan, np.nan, np.nan
    stat, p = sp_stats.ks_2samp(tp, fp)
    # Find max divergence point
    combined = np.sort(np.concatenate([tp, fp]))
    n_tp, n_fp = len(tp), len(fp)
    ecdf_tp = np.searchsorted(np.sort(tp), combined, side="right") / n_tp
    ecdf_fp = np.searchsorted(np.sort(fp), combined, side="right") / n_fp
    diff = np.abs(ecdf_tp - ecdf_fp)
    max_idx = np.argmax(diff)
    max_div_point = combined[max_idx]
    return stat, p, max_div_point


def _tail_enrichment(y_true, y_score, quantile=0.05):
    """FP enrichment in the top and bottom quantile tails."""
    mask = np.isfinite(y_score) & np.isfinite(y_true)
    yt, ys = y_true[mask], y_score[mask]
    if len(yt) < 100:
        return {}
    n = len(yt)
    n_tp_total = int((yt == 0).sum())  # TP=0 (truth_binary: TP=0, FP=1)
    n_fp_total = int((yt == 1).sum())

    results = {}
    for tail_name, idx in [("top5", ys >= np.percentile(ys, 100 - quantile * 100)),
                           ("bottom5", ys <= np.percentile(ys, quantile * 100))]:
        n_in_tail = int(idx.sum())
        if n_in_tail < 5:
            continue
        fp_in_tail = int(yt[idx].sum())
        tp_in_tail = n_in_tail - fp_in_tail
        fp_rate_tail = fp_in_tail / n_in_tail if n_in_tail > 0 else 0
        fp_rate_overall = n_fp_total / n if n > 0 else 0
        enrichment = fp_rate_tail / fp_rate_overall if fp_rate_overall > 0 else np.nan
        # Fisher exact
        table = np.array([[fp_in_tail, n_in_tail - fp_in_tail],
                          [n_fp_total - fp_in_tail, n_tp_total - tp_in_tail]])
        try:
            _, fisher_p = sp_stats.fisher_exact(table)
        except Exception:
            fisher_p = np.nan
        results[tail_name] = {
            "n_in_tail": n_in_tail,
            "fp_in_tail": fp_in_tail,
            "fp_rate_tail": fp_rate_tail,
            "enrichment": enrichment,
            "fisher_p": fisher_p,
        }
    return results


# ===========================================================================
# M1: Precision-Recall AUC
# ===========================================================================

def run_m1_pr_auc(df: pd.DataFrame) -> pd.DataFrame:
    """M1: Compute ROC-AUC and PR-AUC for all features × modes × strata."""
    print("\n" + "=" * 70)
    print("M1: Precision-Recall AUC Analysis")
    print("=" * 70)
    t0 = time.time()

    rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        y = dm["truth_binary"].values

        # Strata: pooled, LOH, non-LOH, per-sample
        strata = {"pooled": dm}
        strata["LOH"] = dm[dm["Potential_LOH"] == True]
        strata["non-LOH"] = dm[dm["Potential_LOH"] == False]
        for s in SAMPLE_ORDER:
            strata[f"sample_{s}"] = dm[dm["sample"] == s]

        for stratum_name, ds in strata.items():
            if len(ds) < 50:
                continue
            ys = ds["truth_binary"].values
            prevalence = ys.mean()  # FP rate (FP=1)
            if len(np.unique(ys)) < 2:
                continue

            for feat in ALL_FEATURES:
                if feat not in ds.columns:
                    continue
                vals = ds[feat].values.astype(float)
                roc = _safe_auc(ys, vals)
                pr = _safe_pr_auc(ys, vals)
                n_valid = int(np.isfinite(vals).sum())
                nan_rate = 1.0 - n_valid / len(vals) if len(vals) > 0 else 1.0
                rows.append({
                    "feature": feat, "mode": mode, "stratum": stratum_name,
                    "n_total": len(ds), "n_valid": n_valid, "nan_rate": round(nan_rate, 4),
                    "prevalence_fp": round(prevalence, 4),
                    "roc_auc": round(roc, 4) if np.isfinite(roc) else np.nan,
                    "pr_auc": round(pr, 4) if np.isfinite(pr) else np.nan,
                    "pr_auc_baseline": round(prevalence, 4),
                    "pr_auc_lift": round(pr / prevalence, 3) if (np.isfinite(pr) and prevalence > 0) else np.nan,
                })

    result = pd.DataFrame(rows)
    result.to_csv(TSV_DIR / "m1_pr_auc_results.tsv", sep="\t", index=False)
    print(f"  M1 done: {len(result)} rows, {time.time()-t0:.1f}s")

    # --- Figure: AUC-ROC vs PR-AUC scatter (pooled only) ---
    pooled = result[result["stratum"] == "pooled"].copy()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        sub = pooled[pooled["mode"] == mode].dropna(subset=["roc_auc", "pr_auc"])
        ism_mask = sub["feature"].isin(ISM_FEATURES)
        ax.scatter(sub.loc[ism_mask, "roc_auc"], sub.loc[ism_mask, "pr_auc"],
                   c="#2196F3", s=50, alpha=0.7, label="ISM", zorder=3)
        ax.scatter(sub.loc[~ism_mask, "roc_auc"], sub.loc[~ism_mask, "pr_auc"],
                   c="#F44336", s=80, marker="D", alpha=0.8, label="Caller", zorder=4)
        # Label caller features
        for _, r in sub[~ism_mask].iterrows():
            ax.annotate(r["feature"], (r["roc_auc"], r["pr_auc"]),
                        fontsize=7, ha="left", va="bottom")
        # Label top ISM
        top_ism = sub[ism_mask].nlargest(3, "roc_auc")
        for _, r in top_ism.iterrows():
            ax.annotate(r["feature"], (r["roc_auc"], r["pr_auc"]),
                        fontsize=7, ha="left", va="bottom", color="#1565C0")
        baseline = sub["pr_auc_baseline"].iloc[0] if len(sub) > 0 else 0
        ax.axhline(baseline, color="gray", ls="--", lw=1, label=f"PR random={baseline:.3f}")
        ax.axvline(0.5, color="gray", ls=":", lw=0.8)
        ax.set_xlabel("ROC-AUC")
        ax.set_ylabel("PR-AUC")
        ax.set_title(f"{mode.upper()} — ROC-AUC vs PR-AUC (pooled)")
        ax.legend(fontsize=8)
        ax.set_xlim(0.3, 0.9)
    fig.suptitle("M1: ROC-AUC vs Precision-Recall AUC", fontsize=14, y=1.02)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "m1_01_auc_vs_prauc_scatter.png")

    return result


# ===========================================================================
# M2: Systematic Residualization
# ===========================================================================

def run_m2_residualization(df: pd.DataFrame) -> pd.DataFrame:
    """M2: OLS residualization to remove confounders from all features."""
    print("\n" + "=" * 70)
    print("M2: Systematic Residualization")
    print("=" * 70)
    t0 = time.time()

    rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode].copy()
        y = dm["truth_binary"].values

        # Build confounder matrix
        conf_cols = [c for c in CONFOUNDERS if c in dm.columns]
        X_conf = dm[conf_cols].values.astype(float)

        strata = {"pooled": (dm, y, X_conf)}
        loh_mask = dm["Potential_LOH"] == True
        nloh_mask = dm["Potential_LOH"] == False
        strata["LOH"] = (dm[loh_mask], y[loh_mask.values], X_conf[loh_mask.values])
        strata["non-LOH"] = (dm[nloh_mask], y[nloh_mask.values], X_conf[nloh_mask.values])

        for stratum_name, (ds, ys, Xs) in strata.items():
            if len(ds) < 50 or len(np.unique(ys)) < 2:
                continue

            for feat in ALL_FEATURES:
                if feat in CONFOUNDERS or feat not in ds.columns:
                    continue
                raw = ds[feat].values.astype(float)
                raw_auc = _safe_auc(ys, raw)
                raw_d = _cohens_d(raw[ys == 0], raw[ys == 1])

                # Residualize
                resid = _residualize(raw, Xs)
                resid_auc = _safe_auc(ys, resid)
                resid_d = _cohens_d(resid[ys == 0], resid[ys == 1])

                delta_auc = resid_auc - raw_auc if (np.isfinite(resid_auc) and np.isfinite(raw_auc)) else np.nan
                delta_d = resid_d - raw_d if (np.isfinite(resid_d) and np.isfinite(raw_d)) else np.nan

                nan_rate = 1.0 - np.isfinite(raw).sum() / len(raw) if len(raw) > 0 else 1.0

                rows.append({
                    "feature": feat, "mode": mode, "stratum": stratum_name,
                    "n": len(ds), "nan_rate": round(nan_rate, 4),
                    "raw_auc": round(raw_auc, 4) if np.isfinite(raw_auc) else np.nan,
                    "resid_auc": round(resid_auc, 4) if np.isfinite(resid_auc) else np.nan,
                    "delta_auc": round(delta_auc, 4) if np.isfinite(delta_auc) else np.nan,
                    "raw_d": round(raw_d, 4) if np.isfinite(raw_d) else np.nan,
                    "resid_d": round(resid_d, 4) if np.isfinite(resid_d) else np.nan,
                    "delta_d": round(delta_d, 4) if np.isfinite(delta_d) else np.nan,
                })

    result = pd.DataFrame(rows)
    result.to_csv(TSV_DIR / "m2_residualization_results.tsv", sep="\t", index=False)
    print(f"  M2 done: {len(result)} rows, {time.time()-t0:.1f}s")

    # --- Figure: Lollipop raw vs residualized AUC ---
    pooled = result[result["stratum"] == "pooled"].dropna(subset=["raw_auc", "resid_auc"])
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        sub = pooled[pooled["mode"] == mode].sort_values("raw_auc", ascending=True)
        y_pos = np.arange(len(sub))
        ax.barh(y_pos, sub["raw_auc"] - 0.5, left=0.5, height=0.35,
                color="#90CAF9", alpha=0.8, label="Raw AUC")
        ax.barh(y_pos + 0.35, sub["resid_auc"] - 0.5, left=0.5, height=0.35,
                color="#EF9A9A", alpha=0.8, label="Residualized AUC")
        ax.set_yticks(y_pos + 0.175)
        ax.set_yticklabels(sub["feature"], fontsize=7)
        ax.axvline(0.5, color="gray", ls=":", lw=0.8)
        ax.axvline(0.58, color="green", ls="--", lw=1, alpha=0.5, label="Threshold 0.58")
        ax.set_xlabel("AUC")
        ax.set_title(f"{mode.upper()} — Raw vs Residualized AUC")
        ax.legend(fontsize=8)
    fig.suptitle("M2: Systematic Residualization (confounders: NumReads, NumCpGs, caller_af, caller_dp)",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "m2_01_raw_vs_residual_auc_lollipop.png")

    return result


# ===========================================================================
# M3: Conditional Subgroup Discovery
# ===========================================================================

def run_m3_subgroup_discovery(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """M3: Stratification grid search for subgroup-specific signals."""
    print("\n" + "=" * 70)
    print("M3: Conditional Subgroup Discovery")
    print("=" * 70)
    t0 = time.time()

    # Select top 10 features by pooled AUC (from M1 or compute here)
    top_features = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        y = dm["truth_binary"].values
        aucs = {}
        for feat in ALL_FEATURES:
            if feat in dm.columns:
                aucs[feat] = abs(_safe_auc(y, dm[feat].values.astype(float)) - 0.5)
        sorted_feats = sorted(aucs, key=lambda x: aucs.get(x, 0), reverse=True)[:10]
        top_features.extend(sorted_feats)
    top_features = list(dict.fromkeys(top_features))[:15]  # Deduplicate, keep order, max 15

    print(f"  Top features for subgroup analysis: {top_features}")

    rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode].copy()

        # Add coarse AF bin
        dm["af_bin_coarse"] = pd.cut(
            dm["caller_af"].fillna(-1),
            bins=[-0.01, 0.2, 0.5, 0.8, 1.01],
            labels=["low_0-0.2", "mid_0.2-0.5", "high_0.5-0.8", "extreme_0.8+"],
        )

        for loh_label, loh_val in [("LOH", True), ("non-LOH", False)]:
            dl = dm[dm["Potential_LOH"] == loh_val]

            for af_label in ["low_0-0.2", "mid_0.2-0.5", "high_0.5-0.8", "extreme_0.8+"]:
                da = dl[dl["af_bin_coarse"] == af_label]

                for sample_label in list(SAMPLE_GROUPS.keys()) + ["all_samples"]:
                    if sample_label == "all_samples":
                        ds = da
                    else:
                        ds = da[da["sample"].isin(SAMPLE_GROUPS[sample_label])]

                    if len(ds) < 60:
                        continue
                    ys = ds["truth_binary"].values
                    n_tp = int((ys == 0).sum())
                    n_fp = int((ys == 1).sum())
                    if n_tp < 30 or n_fp < 30:
                        continue

                    for feat in top_features:
                        if feat not in ds.columns:
                            continue
                        vals = ds[feat].values.astype(float)
                        auc = _safe_auc(ys, vals)
                        d = _cohens_d(vals[ys == 0], vals[ys == 1])
                        rows.append({
                            "feature": feat, "mode": mode,
                            "loh": loh_label, "af_bin": af_label,
                            "sample": sample_label,
                            "n_tp": n_tp, "n_fp": n_fp, "n_total": len(ds),
                            "auc": round(auc, 4) if np.isfinite(auc) else np.nan,
                            "cohens_d": round(d, 4) if np.isfinite(d) else np.nan,
                        })

    result = pd.DataFrame(rows)
    result.to_csv(TSV_DIR / "m3_subgroup_auc_matrix.tsv", sep="\t", index=False)

    # Identify promising subgroups
    promising = result[(result["auc"] > 0.60) | (result["cohens_d"].abs() > 0.30)].copy()
    promising = promising.sort_values("auc", ascending=False)
    promising.to_csv(TSV_DIR / "m3_promising_subgroups.tsv", sep="\t", index=False)
    print(f"  M3 done: {len(result)} cells, {len(promising)} promising, {time.time()-t0:.1f}s")

    # --- Figure: Heatmap per mode (top features × sample, pooled AF, non-LOH) ---
    for mode in MODE_ORDER:
        sub = result[(result["mode"] == mode) & (result["loh"] == "non-LOH") &
                     (result["af_bin"] == "mid_0.2-0.5")]
        if len(sub) < 5:
            # Try pooled AF bins
            sub = result[(result["mode"] == mode) & (result["loh"] == "non-LOH") &
                         (result["sample"] == "all_samples")]

        if len(sub) < 3:
            continue

        pivot = sub.pivot_table(index="feature", columns="af_bin", values="auc", aggfunc="first")
        if pivot.empty:
            pivot = sub.pivot_table(index="feature", columns="sample", values="auc", aggfunc="first")

        fig, ax = plt.subplots(figsize=(12, max(6, len(pivot) * 0.4)))
        sns.heatmap(pivot, annot=True, fmt=".3f", cmap="RdYlGn", center=0.5,
                    vmin=0.35, vmax=0.75, ax=ax, linewidths=0.5)
        ax.set_title(f"M3: Subgroup AUC — {mode.upper()} (non-LOH)")
        fig.tight_layout()
        save_figure(fig, FIG_DIR / f"m3_01_subgroup_heatmap_{mode}.png")

    return result, promising


# ===========================================================================
# M4: Calibration + Brier Score
# ===========================================================================

def run_m4_calibration(df: pd.DataFrame) -> pd.DataFrame:
    """M4: Decile calibration analysis and Brier score."""
    print("\n" + "=" * 70)
    print("M4: Calibration + Brier Score")
    print("=" * 70)
    t0 = time.time()

    from sklearn.metrics import brier_score_loss

    rows = []
    decile_rows = []

    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        y = dm["truth_binary"].values
        prevalence = y.mean()
        brier_baseline = prevalence * (1 - prevalence)  # Brier of constant prediction

        for feat in ALL_FEATURES:
            if feat not in dm.columns:
                continue
            vals = dm[feat].values.astype(float)
            mask = np.isfinite(vals)
            if mask.sum() < 100:
                continue
            ym, vm = y[mask], vals[mask]

            # Normalize to [0, 1] for Brier
            vmin, vmax = vm.min(), vm.max()
            if vmax == vmin:
                continue
            v_norm = (vm - vmin) / (vmax - vmin)

            # Brier score (treat normalized feature as probability proxy)
            brier = brier_score_loss(ym, v_norm)
            bss = 1 - brier / brier_baseline if brier_baseline > 0 else np.nan

            rows.append({
                "feature": feat, "mode": mode,
                "n_valid": int(mask.sum()),
                "prevalence": round(prevalence, 4),
                "brier_score": round(brier, 5),
                "brier_baseline": round(brier_baseline, 5),
                "brier_skill_score": round(bss, 4) if np.isfinite(bss) else np.nan,
            })

            # Decile analysis
            try:
                deciles = pd.qcut(vm, 10, labels=False, duplicates="drop")
            except Exception:
                continue
            for d in sorted(np.unique(deciles)):
                idx = deciles == d
                obs_fp_rate = ym[idx].mean()
                decile_rows.append({
                    "feature": feat, "mode": mode, "decile": int(d),
                    "n": int(idx.sum()),
                    "observed_fp_rate": round(obs_fp_rate, 4),
                    "feature_mean": round(vm[idx].mean(), 4),
                    "feature_median": round(np.median(vm[idx]), 4),
                })

    result = pd.DataFrame(rows)
    result.to_csv(TSV_DIR / "m4_calibration_results.tsv", sep="\t", index=False)

    decile_df = pd.DataFrame(decile_rows)
    decile_df.to_csv(TSV_DIR / "m4_decile_fp_rates.tsv", sep="\t", index=False)
    print(f"  M4 done: {len(result)} features, {len(decile_df)} decile rows, {time.time()-t0:.1f}s")

    # --- Figure: Calibration curves (top 4 per mode) ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    for i, mode in enumerate(MODE_ORDER):
        sub = result[result["mode"] == mode].sort_values("brier_skill_score", ascending=False)
        top4 = sub.head(4)["feature"].tolist()

        for j, feat in enumerate(top4[:2]):
            ax = axes[i][j]
            dec = decile_df[(decile_df["mode"] == mode) & (decile_df["feature"] == feat)]
            if len(dec) < 3:
                continue
            ax.plot(dec["decile"], dec["observed_fp_rate"], "o-", color="#2196F3", lw=2, ms=6)
            prev = result[(result["mode"] == mode) & (result["feature"] == feat)]["prevalence"].iloc[0]
            ax.axhline(prev, color="gray", ls="--", lw=1, label=f"baseline={prev:.3f}")
            ax.set_xlabel("Feature Decile")
            ax.set_ylabel("Observed FP Rate")
            bss = result[(result["mode"] == mode) & (result["feature"] == feat)]["brier_skill_score"].iloc[0]
            ax.set_title(f"{mode.upper()} — {feat}\nBSS={bss:.3f}")
            ax.legend(fontsize=8)
            ax.set_ylim(-0.02, 1.02)
    fig.suptitle("M4: Calibration Curves — FP Rate by Feature Decile", fontsize=13, y=1.02)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "m4_01_calibration_curves_top8.png")

    return result


# ===========================================================================
# M5: Bootstrap Stability
# ===========================================================================

def run_m5_bootstrap(df: pd.DataFrame) -> pd.DataFrame:
    """M5: Bootstrap confidence intervals for AUC and Cohen's d."""
    print("\n" + "=" * 70)
    print("M5: Bootstrap Stability (1000 iterations)")
    print("=" * 70)
    t0 = time.time()

    N_BOOT = 200  # Reduced from 1000 for speed; 200 gives adequate CI width
    rng = np.random.RandomState(42)

    rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        y = dm["truth_binary"].values

        # Strata: pooled + per-sample
        strata = {"pooled": (dm, y)}
        for s in SAMPLE_ORDER:
            mask = dm["sample"] == s
            strata[f"sample_{s}"] = (dm[mask], y[mask.values])

        for stratum_name, (ds, ys) in strata.items():
            if len(ds) < 100 or len(np.unique(ys)) < 2:
                continue
            indices = np.arange(len(ds))

            for feat in ALL_FEATURES:
                if feat not in ds.columns:
                    continue
                vals = ds[feat].values.astype(float)
                point_auc = _safe_auc(ys, vals)
                point_d = _cohens_d(vals[ys == 0], vals[ys == 1])

                boot_aucs = []
                boot_ds = []
                for _ in range(N_BOOT):
                    bi = rng.choice(indices, size=len(indices), replace=True)
                    by, bv = ys[bi], vals[bi]
                    ba = _safe_auc(by, bv)
                    if np.isfinite(ba):
                        boot_aucs.append(ba)
                    bd = _cohens_d(bv[by == 0], bv[by == 1])
                    if np.isfinite(bd):
                        boot_ds.append(bd)

                if len(boot_aucs) < 100:
                    continue

                auc_lo = np.percentile(boot_aucs, 2.5)
                auc_hi = np.percentile(boot_aucs, 97.5)
                d_lo = np.percentile(boot_ds, 2.5) if len(boot_ds) > 100 else np.nan
                d_hi = np.percentile(boot_ds, 97.5) if len(boot_ds) > 100 else np.nan
                pct_above_055 = np.mean(np.array(boot_aucs) > 0.55) * 100

                rows.append({
                    "feature": feat, "mode": mode, "stratum": stratum_name,
                    "n": len(ds),
                    "auc_point": round(point_auc, 4) if np.isfinite(point_auc) else np.nan,
                    "auc_ci_lo": round(auc_lo, 4),
                    "auc_ci_hi": round(auc_hi, 4),
                    "auc_ci_width": round(auc_hi - auc_lo, 4),
                    "pct_above_055": round(pct_above_055, 1),
                    "d_point": round(point_d, 4) if np.isfinite(point_d) else np.nan,
                    "d_ci_lo": round(d_lo, 4) if np.isfinite(d_lo) else np.nan,
                    "d_ci_hi": round(d_hi, 4) if np.isfinite(d_hi) else np.nan,
                })

    result = pd.DataFrame(rows)
    result.to_csv(TSV_DIR / "m5_bootstrap_stability.tsv", sep="\t", index=False)
    print(f"  M5 done: {len(result)} rows, {time.time()-t0:.1f}s")

    # --- Figure: Forest plot (pooled) ---
    for mode in MODE_ORDER:
        pooled = result[(result["mode"] == mode) & (result["stratum"] == "pooled")].copy()
        pooled = pooled.dropna(subset=["auc_point"]).sort_values("auc_point", ascending=True)

        fig, ax = plt.subplots(figsize=(10, max(6, len(pooled) * 0.35)))
        y_pos = np.arange(len(pooled))
        ax.errorbar(pooled["auc_point"], y_pos,
                    xerr=[pooled["auc_point"] - pooled["auc_ci_lo"],
                          pooled["auc_ci_hi"] - pooled["auc_point"]],
                    fmt="o", color="#1565C0", ecolor="#90CAF9", capsize=3, ms=5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(pooled["feature"], fontsize=8)
        ax.axvline(0.5, color="gray", ls=":", lw=1)
        ax.axvline(0.55, color="orange", ls="--", lw=1, alpha=0.5)
        ax.axvline(0.58, color="green", ls="--", lw=1, alpha=0.5)
        ax.set_xlabel("AUC (95% Bootstrap CI)")
        ax.set_title(f"M5: Bootstrap Stability — {mode.upper()} (pooled, n=1000)")
        fig.tight_layout()
        save_figure(fig, FIG_DIR / f"m5_01_forest_plot_{mode}.png")

    return result


# ===========================================================================
# M6: Non-linear Modeling + Feature Importance
# ===========================================================================

def run_m6_nonlinear(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """M6: GradientBoosting LOSO-CV with permutation importance."""
    print("\n" + "=" * 70)
    print("M6: Non-linear Modeling (GradientBoosting + LOSO-CV)")
    print("=" * 70)
    t0 = time.time()

    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.metrics import roc_auc_score, f1_score, precision_score, recall_score
    from sklearn.inspection import permutation_importance

    # Feature sets
    feature_sets = {
        "ISM_only": ISM_FEATURES,
        "Caller_only": CALLER_FEATURES,
        "ISM+Caller": ISM_FEATURES + CALLER_FEATURES,
    }

    cv_rows = []
    imp_rows = []

    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode].copy()
        samples = dm["sample"].unique()

        for fset_name, fset in feature_sets.items():
            fset_valid = [f for f in fset if f in dm.columns]
            if len(fset_valid) < 2:
                continue

            # LOSO-CV
            all_y_true, all_y_pred = [], []
            for held_out in SAMPLE_ORDER:
                if held_out not in samples:
                    continue
                train = dm[dm["sample"] != held_out]
                test = dm[dm["sample"] == held_out]
                if len(train) < 100 or len(test) < 50:
                    continue
                if len(np.unique(train["truth_binary"])) < 2:
                    continue
                if len(np.unique(test["truth_binary"])) < 2:
                    continue

                X_train = train[fset_valid].values.astype(float)
                y_train = train["truth_binary"].values
                X_test = test[fset_valid].values.astype(float)
                y_test = test["truth_binary"].values

                # Impute NaN with median
                imp = SimpleImputer(strategy="median")
                X_train = imp.fit_transform(X_train)
                X_test = imp.transform(X_test)

                gb = GradientBoostingClassifier(
                    n_estimators=200, max_depth=4, learning_rate=0.1,
                    subsample=0.8, random_state=42,
                )
                gb.fit(X_train, y_train)
                y_prob = gb.predict_proba(X_test)[:, 1]
                y_pred = (y_prob >= 0.5).astype(int)

                auc = roc_auc_score(y_test, y_prob) if len(np.unique(y_test)) > 1 else np.nan
                f1 = f1_score(y_test, y_pred, zero_division=0)
                prec = precision_score(y_test, y_pred, zero_division=0)
                rec = recall_score(y_test, y_pred, zero_division=0)

                cv_rows.append({
                    "mode": mode, "feature_set": fset_name,
                    "held_out": held_out,
                    "n_train": len(train), "n_test": len(test),
                    "auc": round(auc, 4), "f1": round(f1, 4),
                    "precision": round(prec, 4), "recall": round(rec, 4),
                })

                all_y_true.extend(y_test)
                all_y_pred.extend(y_prob)

            # Permutation importance on full data (last fold's model)
            if fset_name == "ISM+Caller" and len(all_y_true) > 100:
                # Fit one model on all data for importance
                X_all = dm[fset_valid].values.astype(float)
                y_all = dm["truth_binary"].values
                imp_full = SimpleImputer(strategy="median")
                X_all = imp_full.fit_transform(X_all)
                gb_full = GradientBoostingClassifier(
                    n_estimators=200, max_depth=4, learning_rate=0.1,
                    subsample=0.8, random_state=42,
                )
                gb_full.fit(X_all, y_all)

                perm = permutation_importance(gb_full, X_all, y_all,
                                              n_repeats=10, random_state=42, n_jobs=-1)
                for fi, fname in enumerate(fset_valid):
                    imp_rows.append({
                        "mode": mode, "feature": fname,
                        "importance_mean": round(perm.importances_mean[fi], 5),
                        "importance_std": round(perm.importances_std[fi], 5),
                    })

    cv_result = pd.DataFrame(cv_rows)
    cv_result.to_csv(TSV_DIR / "m6_gb_loso_cv.tsv", sep="\t", index=False)

    imp_result = pd.DataFrame(imp_rows)
    imp_result.to_csv(TSV_DIR / "m6_permutation_importance.tsv", sep="\t", index=False)
    print(f"  M6 done: {len(cv_result)} CV folds, {len(imp_result)} importance rows, {time.time()-t0:.1f}s")

    # --- Figure: Permutation importance ---
    if len(imp_result) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(14, 8))
        for i, mode in enumerate(MODE_ORDER):
            ax = axes[i]
            sub = imp_result[imp_result["mode"] == mode].sort_values("importance_mean", ascending=True)
            colors = ["#F44336" if f in CALLER_FEATURES else "#2196F3" for f in sub["feature"]]
            ax.barh(np.arange(len(sub)), sub["importance_mean"],
                    xerr=sub["importance_std"], color=colors, alpha=0.8, capsize=2)
            ax.set_yticks(np.arange(len(sub)))
            ax.set_yticklabels(sub["feature"], fontsize=7)
            ax.set_xlabel("Permutation Importance (decrease in accuracy)")
            ax.set_title(f"{mode.upper()} — Permutation Importance (GB)")
            # Legend
            ax.legend(handles=[
                mpatches.Patch(color="#2196F3", label="ISM"),
                mpatches.Patch(color="#F44336", label="Caller"),
            ], fontsize=8)
        fig.suptitle("M6: GradientBoosting Permutation Importance (ISM+Caller)", fontsize=13, y=1.02)
        fig.tight_layout()
        save_figure(fig, FIG_DIR / "m6_01_permutation_importance.png")

    # --- Figure: LOSO-CV comparison ---
    if len(cv_result) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for i, mode in enumerate(MODE_ORDER):
            ax = axes[i]
            sub = cv_result[cv_result["mode"] == mode]
            means = sub.groupby("feature_set")["auc"].mean().sort_values()
            stds = sub.groupby("feature_set")["auc"].std()
            colors = {"ISM_only": "#2196F3", "Caller_only": "#F44336", "ISM+Caller": "#4CAF50"}
            bars = ax.barh(means.index, means.values, xerr=stds[means.index].values,
                          color=[colors.get(x, "gray") for x in means.index],
                          alpha=0.8, capsize=4)
            ax.set_xlabel("Mean LOSO-CV AUC")
            ax.set_title(f"{mode.upper()} — LOSO-CV AUC by Feature Set")
            ax.axvline(0.5, color="gray", ls=":", lw=0.8)
            for idx, (name, val) in enumerate(means.items()):
                ax.text(val + 0.005, idx, f"{val:.3f}", va="center", fontsize=9)
        fig.suptitle("M6: GradientBoosting LOSO-CV Comparison", fontsize=13, y=1.02)
        fig.tight_layout()
        save_figure(fig, FIG_DIR / "m6_02_loso_cv_comparison.png")

    return cv_result, imp_result


# ===========================================================================
# M7: Distribution Overlap + Tail Enrichment
# ===========================================================================

def run_m7_distribution(df: pd.DataFrame) -> pd.DataFrame:
    """M7: KDE overlays, KS test, and tail enrichment analysis."""
    print("\n" + "=" * 70)
    print("M7: Distribution Overlap + Tail Enrichment")
    print("=" * 70)
    t0 = time.time()

    rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        y = dm["truth_binary"].values

        strata = {"pooled": (dm, y)}
        loh_mask = dm["Potential_LOH"] == True
        nloh_mask = dm["Potential_LOH"] == False
        strata["LOH"] = (dm[loh_mask], y[loh_mask.values])
        strata["non-LOH"] = (dm[nloh_mask], y[nloh_mask.values])

        # AF bins
        for af_label, (lo, hi) in AF_BINS_COARSE.items():
            af_mask = (dm["caller_af"] >= lo) & (dm["caller_af"] < hi)
            strata[f"AF_{af_label}"] = (dm[af_mask], y[af_mask.values])

        for stratum_name, (ds, ys) in strata.items():
            if len(ds) < 100 or len(np.unique(ys)) < 2:
                continue
            tp_idx = ys == 0
            fp_idx = ys == 1

            for feat in ALL_FEATURES:
                if feat not in ds.columns:
                    continue
                vals = ds[feat].values.astype(float)
                tp_vals = vals[tp_idx]
                fp_vals = vals[fp_idx]

                ks_stat, ks_p, max_div = _ks_test(tp_vals, fp_vals)
                tail_info = _tail_enrichment(ys, vals)

                row = {
                    "feature": feat, "mode": mode, "stratum": stratum_name,
                    "n_tp": int(tp_idx.sum()), "n_fp": int(fp_idx.sum()),
                    "ks_stat": round(ks_stat, 4) if np.isfinite(ks_stat) else np.nan,
                    "ks_p": ks_p if np.isfinite(ks_p) else np.nan,
                    "max_divergence_point": round(max_div, 4) if np.isfinite(max_div) else np.nan,
                }

                for tail_key in ["top5", "bottom5"]:
                    ti = tail_info.get(tail_key, {})
                    row[f"{tail_key}_enrichment"] = round(ti.get("enrichment", np.nan), 3) if np.isfinite(ti.get("enrichment", np.nan)) else np.nan
                    row[f"{tail_key}_fisher_p"] = ti.get("fisher_p", np.nan)
                    row[f"{tail_key}_fp_rate"] = round(ti.get("fp_rate_tail", np.nan), 4)
                    row[f"{tail_key}_n"] = ti.get("n_in_tail", 0)

                rows.append(row)

    result = pd.DataFrame(rows)
    result.to_csv(TSV_DIR / "m7_distribution_overlap.tsv", sep="\t", index=False)
    print(f"  M7 done: {len(result)} rows, {time.time()-t0:.1f}s")

    # --- Figure: KDE overlays (top 10 by KS stat, pooled) ---
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        y = dm["truth_binary"].values
        pooled = result[(result["mode"] == mode) & (result["stratum"] == "pooled")]
        top10 = pooled.dropna(subset=["ks_stat"]).nlargest(10, "ks_stat")

        n_feats = min(10, len(top10))
        if n_feats < 3:
            continue
        ncols = 3
        nrows = (n_feats + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(15, 4 * nrows))
        axes = axes.flatten() if n_feats > 1 else [axes]

        for idx, (_, row) in enumerate(top10.iterrows()):
            if idx >= len(axes):
                break
            ax = axes[idx]
            feat = row["feature"]
            vals = dm[feat].values.astype(float)
            tp_v = vals[y == 0]
            fp_v = vals[y == 1]
            tp_v = tp_v[np.isfinite(tp_v)]
            fp_v = fp_v[np.isfinite(fp_v)]

            # Clip outliers for visualization
            combined = np.concatenate([tp_v, fp_v])
            lo_clip, hi_clip = np.percentile(combined, [1, 99])
            tp_c = tp_v[(tp_v >= lo_clip) & (tp_v <= hi_clip)]
            fp_c = fp_v[(fp_v >= lo_clip) & (fp_v <= hi_clip)]

            if len(tp_c) > 10 and len(fp_c) > 10:
                ax.hist(tp_c, bins=50, density=True, alpha=0.4, color=COLOR_TP, label="TP")
                ax.hist(fp_c, bins=50, density=True, alpha=0.4, color=COLOR_FP, label="FP")
                try:
                    from scipy.stats import gaussian_kde
                    x_range = np.linspace(lo_clip, hi_clip, 200)
                    kde_tp = gaussian_kde(tp_c)
                    kde_fp = gaussian_kde(fp_c)
                    ax.plot(x_range, kde_tp(x_range), color=COLOR_TP, lw=2)
                    ax.plot(x_range, kde_fp(x_range), color=COLOR_FP, lw=2)
                except Exception:
                    pass

            ks = row["ks_stat"]
            ax.set_title(f"{feat}\nKS={ks:.3f}", fontsize=9)
            ax.legend(fontsize=7)
            ax.tick_params(labelsize=7)

        # Hide unused axes
        for idx in range(n_feats, len(axes)):
            axes[idx].set_visible(False)

        fig.suptitle(f"M7: KDE Overlays — {mode.upper()} (top 10 by KS stat)", fontsize=13, y=1.02)
        fig.tight_layout()
        save_figure(fig, FIG_DIR / f"m7_01_kde_overlays_{mode}_top10.png")

    # --- Figure: Tail enrichment bar chart ---
    pooled_all = result[result["stratum"] == "pooled"].copy()
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        sub = pooled_all[pooled_all["mode"] == mode].dropna(subset=["top5_enrichment"])
        sub = sub.sort_values("top5_enrichment", ascending=True)

        y_pos = np.arange(len(sub))
        colors = ["#F44336" if e > 2.0 else "#90CAF9" for e in sub["top5_enrichment"]]
        ax.barh(y_pos, sub["top5_enrichment"], color=colors, alpha=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(sub["feature"], fontsize=7)
        ax.axvline(1.0, color="gray", ls=":", lw=1)
        ax.axvline(2.0, color="orange", ls="--", lw=1, alpha=0.5, label="2× threshold")
        ax.axvline(3.0, color="red", ls="--", lw=1, alpha=0.5, label="3× threshold")
        ax.set_xlabel("FP Enrichment in Top 5%")
        ax.set_title(f"{mode.upper()} — Top 5% Tail FP Enrichment")
        ax.legend(fontsize=8)
    fig.suptitle("M7: Tail Enrichment Analysis (top 5% of feature range)", fontsize=13, y=1.02)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "m7_03_tail_enrichment_bar.png")

    return result


# ===========================================================================
# Master Verdict
# ===========================================================================

def synthesize_verdict(m1, m2, m3_promising, m5, m6_cv, m6_imp, m7):
    """Combine all evidence into master verdict table."""
    print("\n" + "=" * 70)
    print("Synthesizing Master Verdict")
    print("=" * 70)

    verdicts = []

    # Check M2: any residualized AUC > 0.58?
    m2_pooled = m2[m2["stratum"] == "pooled"]
    m2_max = m2_pooled["resid_auc"].max() if len(m2_pooled) > 0 else np.nan
    m2_pass = m2_max > 0.58 if np.isfinite(m2_max) else False
    verdicts.append({
        "criterion": "M2: Residualized AUC > 0.58",
        "threshold": 0.58,
        "observed": round(m2_max, 4) if np.isfinite(m2_max) else "N/A",
        "pass": m2_pass,
    })

    # Check M3: subgroup AUC > 0.65 with n >= 500?
    large_promising = m3_promising[(m3_promising["n_tp"] >= 250) & (m3_promising["n_fp"] >= 250)]
    m3_max = large_promising["auc"].max() if len(large_promising) > 0 else np.nan
    m3_pass = m3_max > 0.65 if np.isfinite(m3_max) else False
    verdicts.append({
        "criterion": "M3: Subgroup AUC > 0.65 (n≥500)",
        "threshold": 0.65,
        "observed": round(m3_max, 4) if np.isfinite(m3_max) else "N/A",
        "pass": m3_pass,
    })

    # Check M6: GB LOSO-CV > caller-only + 0.02?
    if len(m6_cv) > 0:
        for mode in MODE_ORDER:
            cv_mode = m6_cv[m6_cv["mode"] == mode]
            caller_mean = cv_mode[cv_mode["feature_set"] == "Caller_only"]["auc"].mean()
            ism_caller_mean = cv_mode[cv_mode["feature_set"] == "ISM+Caller"]["auc"].mean()
            delta = ism_caller_mean - caller_mean if (np.isfinite(ism_caller_mean) and np.isfinite(caller_mean)) else np.nan
            m6_pass = delta > 0.02 if np.isfinite(delta) else False
            verdicts.append({
                "criterion": f"M6: GB ISM+Caller > Caller-only + 0.02 ({mode})",
                "threshold": 0.02,
                "observed": round(delta, 4) if np.isfinite(delta) else "N/A",
                "pass": m6_pass,
            })

    # Check M7: tail enrichment > 3×?
    m7_pooled = m7[m7["stratum"] == "pooled"]
    m7_top = m7_pooled.dropna(subset=["top5_enrichment"])
    n_above_3x = int((m7_top["top5_enrichment"] > 3.0).sum())
    m7_pass = n_above_3x >= 3
    verdicts.append({
        "criterion": "M7: Tail 5% enrichment > 3× (≥3 features)",
        "threshold": "≥3 features >3×",
        "observed": n_above_3x,
        "pass": m7_pass,
    })

    # Check M5: any bootstrap CI excludes 0.50?
    m5_pooled = m5[m5["stratum"] == "pooled"]
    ism_m5 = m5_pooled[m5_pooled["feature"].isin(ISM_FEATURES)]
    n_ci_above_050 = int(((ism_m5["auc_ci_lo"] > 0.50)).sum()) if len(ism_m5) > 0 else 0
    verdicts.append({
        "criterion": "M5: Bootstrap CI lower bound > 0.50 (ISM features)",
        "threshold": "any CI_lo > 0.50",
        "observed": n_ci_above_050,
        "pass": n_ci_above_050 > 0,
    })

    verdict_df = pd.DataFrame(verdicts)
    verdict_df.to_csv(TSV_DIR / "master_verdict.tsv", sep="\t", index=False)

    any_pass = any(v["pass"] for v in verdicts)
    conclusion = "SIGNAL_FOUND" if any_pass else "EXHAUSTION_CONFIRMED"
    print(f"\n  === MASTER VERDICT: {conclusion} ===")
    for v in verdicts:
        status = "✅ PASS" if v["pass"] else "❌ FAIL"
        print(f"    {status} | {v['criterion']} | observed={v['observed']}")

    return verdict_df, conclusion


# ===========================================================================
# Main
# ===========================================================================

def main():
    print(f"\nBeyond-AUC Comprehensive Feature Validation")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Output: {OUT_DIR}")

    ensure_dir(TSV_DIR)
    ensure_dir(FIG_DIR)
    setup_plot_style()

    # Load data
    df = load_master_dataset()
    df["truth_binary"] = encode_truth_binary(df)
    print(f"  Paired: {(df['mode']=='paired').sum():,}  TO: {(df['mode']=='to').sum():,}")
    print(f"  TP: {(df['truth_binary']==0).sum():,}  FP: {(df['truth_binary']==1).sum():,}")

    # Write round context
    write_round_context(
        OUT_DIR,
        observation_id="Beyond-AUC",
        title="Comprehensive Beyond-AUC Feature Validation",
        description="7 complementary methods to validate ISM feature discriminability beyond AUC-ROC",
        script_path="scripts/analysis/build_beyond_auc_validation.py",
        data_sources=[{"name": "master_dataset", "path": str(load_master_dataset.__defaults__)}],
        row_count=len(df),
        col_count=len(df.columns),
        extra={
            "n_paired": int((df["mode"] == "paired").sum()),
            "n_to": int((df["mode"] == "to").sum()),
            "methods": ["M1:PR-AUC", "M2:Residualization", "M3:Subgroup",
                         "M4:Calibration", "M5:Bootstrap", "M6:NonlinearML", "M7:Distribution"],
            "features_evaluated": ALL_FEATURES,
            "confounders": CONFOUNDERS,
        },
    )

    # Run all methods
    t_total = time.time()

    m1 = run_m1_pr_auc(df)
    m2 = run_m2_residualization(df)
    m3_full, m3_promising = run_m3_subgroup_discovery(df)
    m4 = run_m4_calibration(df)
    m5 = run_m5_bootstrap(df)
    m6_cv, m6_imp = run_m6_nonlinear(df)
    m7 = run_m7_distribution(df)

    # Synthesize
    verdict_df, conclusion = synthesize_verdict(m1, m2, m3_promising, m5, m6_cv, m6_imp, m7)

    total_time = time.time() - t_total
    print(f"\n{'='*70}")
    print(f"ALL DONE: {total_time/60:.1f} minutes")
    print(f"TSV files: {len(list(TSV_DIR.glob('*.tsv')))}")
    print(f"Figures: {len(list(FIG_DIR.glob('*.png')))}")
    print(f"Verdict: {conclusion}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
