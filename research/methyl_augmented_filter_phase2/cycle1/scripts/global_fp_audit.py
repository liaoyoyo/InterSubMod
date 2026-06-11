#!/usr/bin/env python3
"""
Phase 2 Cycle 1 Step 0 - Global FP Exploration Audit
=====================================================

4 Stages in 1 script:
  Stage 1: FP Concentration Map  (caller_af bin x LOH inner/outer x Coverage_Multiple bin x chr8 vs other)
  Stage 2: High-AF FP profile    (caller_af > 0.3)
  Stage 3: Global LR sweep       (35,332 rows, no cell gate, 5-fold OOF, tau sweep 0.10..0.95)
  Stage 4: Heterogeneous threshold (top 20 FP-rich cells, cell-specific tau, aggregate)

Output:
  cycle1_step0_fp_concentration_map.tsv
  cycle1_step0_high_af_fp_profile.tsv
  cycle1_step0_global_lr_sweep.tsv
  cycle1_step0_heterogeneous_threshold.tsv
  figures/*.png

Decision criteria pre-registered:
  D1: top 10 FP-rich cells cover >= 70% all FP    -> H_C1_1
  D2: global LR max delta_F1 > +0.00242           -> H_C1_2
  D3: heterogeneous aggregate delta_F1 >= +0.01   -> H_C1_3
  D4: high-AF zone incremental delta_F1 >= +0.003 -> H_C1_4
"""

import argparse
import json
import os
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

# CJK font fallback (per feedback_matplotlib_cjk_font_rule.md)
matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

# ---------- Constants (pre-registered, hard-coded) ----------

# v1.0 baseline (from plan v2.0)
CALLER_F1_BASELINE = 0.7166
N_TRUTH = 49778  # FN reverse-solved per v1.0: paper-spec total truth count
# Per v1.0 reverse-solve: FN = 49778 - TP = 49778 - 30490 = 19288
FN_CALLER = 19288

# Decision thresholds
D1_THRESHOLD_TOP10_COVERAGE = 0.70
D2_THRESHOLD_GLOBAL_DELTA_F1 = 0.00242  # v1.0 cell-gated baseline
D3_THRESHOLD_HETERO_DELTA_F1 = 0.01
D4_THRESHOLD_HIGH_AF_DELTA_F1 = 0.003

# Stratification bins
AF_BINS = [0.0, 0.10, 0.20, 0.30, 0.50, 0.70, 1.01]  # 6 bins
AF_LABELS = ["[0.00,0.10)", "[0.10,0.20)", "[0.20,0.30)", "[0.30,0.50)", "[0.50,0.70)", "[0.70,1.00]"]
COV_BINS = [-np.inf, 0.5, 0.9, 1.1, 1.3, np.inf]  # 5 bins
COV_LABELS = ["<0.5", "[0.5,0.9)", "[0.9,1.1)", "[1.1,1.3)", ">=1.3"]

# tau sweep
TAU_GRID = np.round(np.arange(0.10, 0.96, 0.01), 2)

# LR features (no cell gate; Model B style + 5 methyl covariates)
LR_FEATURES = [
    "V6_off_NG",
    "caller_af",
    "NumReads_master",
    "loh_inner_flag",          # binary
    "Coverage_Multiple_imp",   # imputed
    "V6_off_meth_HPMergedDelta",
    "V6_off_meth_HPFineF",
    "V6_off_meth_NME_imbalance",
    "V6_off_meth_Epipoly_Delta",
    "V6_off_meth_ClusterPermanovaF",
    "chr8_flag",
]

RANDOM_SEED = 42


# ---------- Utility functions ----------

def log(msg, logf=None):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    if logf is not None:
        logf.write(line + "\n")
        logf.flush()


def compute_delta_f1(tp_kept, fp_kept, fn_added):
    """
    Compute delta_F1 vs caller baseline given filter outcome.

    Filter removes some predictions:
      - TP filtered out -> reclassified as FN
      - FP filtered out -> correctly removed

    New counts:
      TP_new = tp_kept
      FP_new = fp_kept
      FN_new = FN_CALLER + (TP_total - tp_kept)  [tp_removed becomes FN]

    Wait: per v1.0 framework, baseline TP=30490, FN=19288 (reverse-solve N_truth=49778).
    Filter rule selects subset to keep. Removed-TP becomes new-FN (we lose recall).
    Removed-FP is gain (we lose false-positives).

    Recall = TP / (TP+FN)
    Precision = TP / (TP+FP)
    F1 = 2 * P * R / (P+R)
    """
    tp_total = 30490
    fp_total = 4842
    tp_removed = tp_total - tp_kept
    fp_removed = fp_total - fp_kept

    new_tp = tp_kept
    new_fp = fp_kept
    new_fn = FN_CALLER + tp_removed

    if new_tp + new_fp == 0 or new_tp + new_fn == 0:
        return 0.0, 0.0, 0.0
    precision = new_tp / (new_tp + new_fp)
    recall = new_tp / (new_tp + new_fn)
    if precision + recall == 0:
        return 0.0, precision, recall
    f1 = 2 * precision * recall / (precision + recall)
    return f1 - CALLER_F1_BASELINE, precision, recall


def prep_dataframe(df):
    """
    Add derived columns:
      - loh_inner_flag (1 if loh_side == Inner else 0)
      - Coverage_Multiple_imp (median imputed)
      - chr8_flag (1 if chr == chr8)
      - af_bin (categorical)
      - cov_bin (categorical)
      - cell_id (4-tuple string)
    """
    df = df.copy()
    df["loh_inner_flag"] = (df["loh_side"] == "Inner").astype(int)
    cov_median = df["Coverage_Multiple"].median()
    df["Coverage_Multiple_imp"] = df["Coverage_Multiple"].fillna(cov_median)
    df["chr8_flag"] = (df["chr"] == "chr8").astype(int)
    df["af_bin"] = pd.cut(df["caller_af"], bins=AF_BINS, labels=AF_LABELS, include_lowest=True, right=False)
    df["cov_bin"] = pd.cut(df["Coverage_Multiple_imp"], bins=COV_BINS, labels=COV_LABELS, right=False)
    df["loh_str"] = np.where(df["loh_inner_flag"] == 1, "Inner", "Outer")
    df["chr8_str"] = np.where(df["chr8_flag"] == 1, "chr8", "other")
    df["cell_id"] = (
        df["af_bin"].astype(str) + " | " +
        df["loh_str"] + " | " +
        df["cov_bin"].astype(str) + " | " +
        df["chr8_str"]
    )
    df["y"] = (df["label"] == "TP").astype(int)
    return df


# ---------- Stage 1: FP concentration map ----------

def stage1_concentration_map(df, outdir, logf):
    log("=" * 70, logf)
    log("STAGE 1 - FP CONCENTRATION MAP", logf)
    log("=" * 70, logf)

    cell_stats = df.groupby("cell_id").agg(
        n_TP=("y", lambda x: int((x == 1).sum())),
        n_FP=("y", lambda x: int((x == 0).sum())),
        n_total=("y", "size"),
    ).reset_index()
    cell_stats["FP_rate"] = cell_stats["n_FP"] / cell_stats["n_total"]
    cell_stats["TP_rate"] = cell_stats["n_TP"] / cell_stats["n_total"]
    total_fp = int((df["y"] == 0).sum())
    cell_stats["FP_share"] = cell_stats["n_FP"] / total_fp

    cell_stats = cell_stats.sort_values("n_FP", ascending=False).reset_index(drop=True)

    top10 = cell_stats.head(10)
    top10_share_sum = top10["FP_share"].sum()
    log(f"Total cells: {len(cell_stats)}", logf)
    log(f"Total FP: {total_fp}", logf)
    log(f"Top 10 FP-rich cells cumulative FP_share: {top10_share_sum*100:.2f}%", logf)
    log("Top 10 FP-rich cells:", logf)
    for _, row in top10.iterrows():
        log(f"  {row['cell_id']:60s}  n_FP={row['n_FP']:4d}  n_TP={row['n_TP']:5d}  FP_rate={row['FP_rate']:.3f}  FP_share={row['FP_share']*100:.1f}%", logf)

    d1_pass = bool(top10_share_sum >= D1_THRESHOLD_TOP10_COVERAGE)
    log(f"D1 verdict: {'PASS' if d1_pass else 'FAIL'} (top10 FP_share={top10_share_sum*100:.2f}%, threshold {D1_THRESHOLD_TOP10_COVERAGE*100:.0f}%)", logf)

    cell_stats.to_csv(outdir / "data" / "cycle1_step0_fp_concentration_map.tsv", sep="\t", index=False)

    # Figure: bar of top 20 cells by FP_share
    top20 = cell_stats.head(20)
    fig, ax = plt.subplots(figsize=(12, 8))
    bars = ax.barh(range(len(top20)), top20["FP_share"].values[::-1] * 100, color="#d62728")
    ax.set_yticks(range(len(top20)))
    ax.set_yticklabels(top20["cell_id"].values[::-1], fontsize=8)
    ax.set_xlabel("FP share (%)")
    ax.set_title(f"Top 20 FP-rich cells (n=120 cells total)\nTop10 cumulative FP_share = {top10_share_sum*100:.1f}% (D1 threshold {D1_THRESHOLD_TOP10_COVERAGE*100:.0f}%, {'PASS' if d1_pass else 'FAIL'})")
    ax.axvline(0, color="black", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "cycle1_step0_fp_concentration_heatmap.png", dpi=150)
    plt.close(fig)

    return {
        "n_cells_total": len(cell_stats),
        "n_cells_nonempty": int((cell_stats["n_total"] > 0).sum()),
        "total_FP": total_fp,
        "top10_FP_share": float(top10_share_sum),
        "top10_cells": top10[["cell_id", "n_FP", "n_TP", "FP_rate", "FP_share"]].to_dict(orient="records"),
        "D1_pass": d1_pass,
        "cell_stats_top20": cell_stats.head(20),
    }


# ---------- Stage 2: High-AF FP profile ----------

def stage2_high_af_profile(df, outdir, logf, stage3_result=None):
    log("=" * 70, logf)
    log("STAGE 2 - HIGH-AF FP PROFILE (caller_af > 0.3)", logf)
    log("=" * 70, logf)

    high_af = df[df["caller_af"] > 0.3]
    high_af_fp = high_af[high_af["y"] == 0]
    high_af_tp = high_af[high_af["y"] == 1]
    log(f"High-AF rows: {len(high_af)} (TP={len(high_af_tp)}, FP={len(high_af_fp)})", logf)
    log(f"  high-AF FP fraction of all FP: {len(high_af_fp)/4842*100:.1f}%", logf)
    log(f"  high-AF TP fraction of all TP: {len(high_af_tp)/30490*100:.1f}%", logf)
    if len(high_af_fp):
        chr8_share = (high_af_fp["chr"] == "chr8").mean()
        log(f"  high-AF FP chr8 share: {chr8_share*100:.1f}%", logf)

    # Per-cell profile within high-AF
    profile = high_af.groupby("cell_id").agg(
        n_TP=("y", lambda x: int((x == 1).sum())),
        n_FP=("y", lambda x: int((x == 0).sum())),
        n_total=("y", "size"),
    ).reset_index()
    profile["FP_rate"] = profile["n_FP"] / profile["n_total"]
    profile = profile.sort_values("n_FP", ascending=False).reset_index(drop=True)
    profile.to_csv(outdir / "data" / "cycle1_step0_high_af_fp_profile.tsv", sep="\t", index=False)

    # D4: high-AF incremental delta_F1 vs full LR using stage3 results
    # Simple proxy: filter rule = "filter all FP in high-AF zone" upper-bound assuming oracle
    # Real D4 needs LR on high-AF subset; here we compute oracle upper bound + LR-based estimate
    # Oracle: remove all high-AF FP (360) and zero high-AF TP loss
    oracle_dF1, oracle_p, oracle_r = compute_delta_f1(tp_kept=30490, fp_kept=4842 - len(high_af_fp), fn_added=0)
    log(f"  Oracle upper bound (filter ALL high-AF FP, keep all TP): delta_F1 = {oracle_dF1:.5f}", logf)

    # LR-based: refit on high-AF subset, find best tau
    d4_pass, d4_best_dF1, d4_best_tau, d4_curve = run_lr_on_subset(high_af, "high_af_zone", logf)
    log(f"D4 high-AF LR best delta_F1 (incremental): {d4_best_dF1:.5f} at tau={d4_best_tau}", logf)
    log(f"D4 verdict: {'PASS' if d4_pass else 'FAIL'} (threshold {D4_THRESHOLD_HIGH_AF_DELTA_F1:+.5f})", logf)

    return {
        "high_af_n": int(len(high_af)),
        "high_af_FP": int(len(high_af_fp)),
        "high_af_TP": int(len(high_af_tp)),
        "high_af_FP_share_of_all_FP": float(len(high_af_fp) / 4842),
        "high_af_chr8_share": float((high_af_fp["chr"] == "chr8").mean()) if len(high_af_fp) else 0.0,
        "D4_oracle_dF1": float(oracle_dF1),
        "D4_lr_best_dF1": float(d4_best_dF1),
        "D4_lr_best_tau": float(d4_best_tau),
        "D4_pass": bool(d4_pass),
    }


def run_lr_on_subset(sub_df, label, logf, n_splits=5):
    """
    Fit LR on subset and find best tau within subset's own scope.
    Returns subset-incremental delta_F1 (i.e., apply rule only to this subset).

    delta_F1 here = F1_after - F1_baseline, where:
      - baseline counts use FULL data (caller F1=0.7166)
      - filter applies only within subset (other rows unchanged)
    """
    sub = sub_df.copy()
    # Drop rows with NaN in core features
    keep_mask = sub[LR_FEATURES].notna().all(axis=1)
    sub_clean = sub[keep_mask].copy()
    if len(sub_clean) < 100:
        log(f"  Subset '{label}': too few rows after dropna ({len(sub_clean)}), skip", logf)
        return False, 0.0, 0.0, pd.DataFrame()

    X = sub_clean[LR_FEATURES].values.astype(float)
    y = sub_clean["y"].values.astype(int)

    # 5-fold OOF predictions
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_SEED)
    oof_pred = np.zeros(len(sub_clean))

    if y.sum() == 0 or y.sum() == len(y):
        log(f"  Subset '{label}': degenerate label distribution, skip", logf)
        return False, 0.0, 0.0, pd.DataFrame()

    for fold_idx, (tr, va) in enumerate(skf.split(X, y)):
        scaler = StandardScaler()
        Xtr = scaler.fit_transform(X[tr])
        Xva = scaler.transform(X[va])
        clf = LogisticRegression(max_iter=2000, C=1.0, solver="lbfgs", random_state=RANDOM_SEED)
        clf.fit(Xtr, y[tr])
        oof_pred[va] = clf.predict_proba(Xva)[:, 1]

    # tau sweep: keep rows with P(TP) >= tau, filter out the rest in subset
    # Apply only to subset, other rows unchanged
    sub_idx = sub_clean.index
    tp_total = 30490
    fp_total = 4842

    # tp/fp within subset
    sub_tp_total = int((y == 1).sum())
    sub_fp_total = int((y == 0).sum())
    # tp/fp outside subset (unchanged by filter)
    out_tp = tp_total - sub_tp_total
    out_fp = fp_total - sub_fp_total

    rows = []
    for tau in TAU_GRID:
        keep_mask_local = oof_pred >= tau
        sub_tp_kept = int(((y == 1) & keep_mask_local).sum())
        sub_fp_kept = int(((y == 0) & keep_mask_local).sum())
        full_tp_kept = out_tp + sub_tp_kept
        full_fp_kept = out_fp + sub_fp_kept
        dF1, prec, rec = compute_delta_f1(full_tp_kept, full_fp_kept, 0)
        rows.append({
            "subset": label,
            "tau": tau,
            "sub_tp_kept": sub_tp_kept,
            "sub_fp_kept": sub_fp_kept,
            "sub_tp_removed": sub_tp_total - sub_tp_kept,
            "sub_fp_removed": sub_fp_total - sub_fp_kept,
            "delta_F1": dF1,
            "precision_new": prec,
            "recall_new": rec,
        })
    curve = pd.DataFrame(rows)
    best = curve.loc[curve["delta_F1"].idxmax()]
    best_dF1 = float(best["delta_F1"])
    best_tau = float(best["tau"])

    return (best_dF1 >= D4_THRESHOLD_HIGH_AF_DELTA_F1), best_dF1, best_tau, curve


# ---------- Stage 3: Global LR sweep ----------

def stage3_global_lr(df, outdir, logf):
    log("=" * 70, logf)
    log("STAGE 3 - GLOBAL LR SWEEP (all 35,332 rows, no cell gate)", logf)
    log("=" * 70, logf)

    results_combined = []

    for strategy in ["A_dropna", "B_impute"]:
        log(f"\n--- Strategy {strategy} ---", logf)
        if strategy == "A_dropna":
            keep = df[LR_FEATURES].notna().all(axis=1)
            sub = df[keep].copy()
            log(f"  Drop NaN rows: kept {len(sub)}/{len(df)}", logf)
        else:
            sub = df.copy()
            for c in LR_FEATURES:
                if sub[c].isna().any():
                    med = sub[c].median()
                    sub[c] = sub[c].fillna(med)
            log(f"  Impute median: kept all {len(sub)} rows", logf)

        X = sub[LR_FEATURES].values.astype(float)
        y = sub["y"].values.astype(int)

        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
        oof_pred = np.zeros(len(sub))
        for fold_idx, (tr, va) in enumerate(skf.split(X, y)):
            scaler = StandardScaler()
            Xtr = scaler.fit_transform(X[tr])
            Xva = scaler.transform(X[va])
            clf = LogisticRegression(max_iter=2000, C=1.0, solver="lbfgs", random_state=RANDOM_SEED)
            clf.fit(Xtr, y[tr])
            oof_pred[va] = clf.predict_proba(Xva)[:, 1]
            log(f"    fold {fold_idx+1}/5 done", logf)

        # tau sweep
        sub_tp_total = int((y == 1).sum())
        sub_fp_total = int((y == 0).sum())
        # rows dropped by Strategy A => unchanged (no filter)
        dropped_tp = 30490 - sub_tp_total
        dropped_fp = 4842 - sub_fp_total

        rows = []
        for tau in TAU_GRID:
            keep_mask = oof_pred >= tau
            tp_kept_in_sub = int(((y == 1) & keep_mask).sum())
            fp_kept_in_sub = int(((y == 0) & keep_mask).sum())
            full_tp_kept = dropped_tp + tp_kept_in_sub
            full_fp_kept = dropped_fp + fp_kept_in_sub
            dF1, prec, rec = compute_delta_f1(full_tp_kept, full_fp_kept, 0)
            rows.append({
                "strategy": strategy,
                "tau": tau,
                "tp_kept": full_tp_kept,
                "fp_kept": full_fp_kept,
                "tp_removed": 30490 - full_tp_kept,
                "fp_removed": 4842 - full_fp_kept,
                "delta_F1": dF1,
                "precision_new": prec,
                "recall_new": rec,
            })
        curve = pd.DataFrame(rows)
        best = curve.loc[curve["delta_F1"].idxmax()]
        log(f"  Best: tau={best['tau']:.2f}  delta_F1={best['delta_F1']:.5f}  TP_rm={int(best['tp_removed'])}  FP_rm={int(best['fp_removed'])}", logf)
        results_combined.append(curve)

    full_curve = pd.concat(results_combined, ignore_index=True)
    full_curve.to_csv(outdir / "data" / "cycle1_step0_global_lr_sweep.tsv", sep="\t", index=False)

    # Best across both strategies
    best_a = full_curve[full_curve["strategy"] == "A_dropna"].loc[full_curve[full_curve["strategy"] == "A_dropna"]["delta_F1"].idxmax()]
    best_b = full_curve[full_curve["strategy"] == "B_impute"].loc[full_curve[full_curve["strategy"] == "B_impute"]["delta_F1"].idxmax()]
    max_dF1 = max(float(best_a["delta_F1"]), float(best_b["delta_F1"]))

    d2_pass = bool(max_dF1 > D2_THRESHOLD_GLOBAL_DELTA_F1)
    log(f"\nD2 verdict: {'PASS' if d2_pass else 'FAIL'}  (max delta_F1 = {max_dF1:.5f}; threshold {D2_THRESHOLD_GLOBAL_DELTA_F1:+.5f})", logf)
    log(f"  Strategy A best: tau={best_a['tau']:.2f} delta_F1={best_a['delta_F1']:.5f}", logf)
    log(f"  Strategy B best: tau={best_b['tau']:.2f} delta_F1={best_b['delta_F1']:.5f}", logf)

    # Figure: delta_F1 vs tau both strategies
    fig, ax = plt.subplots(figsize=(10, 6))
    for strat, color in [("A_dropna", "#1f77b4"), ("B_impute", "#ff7f0e")]:
        sub_c = full_curve[full_curve["strategy"] == strat]
        ax.plot(sub_c["tau"], sub_c["delta_F1"], "-o", markersize=3, label=strat, color=color)
    ax.axhline(D2_THRESHOLD_GLOBAL_DELTA_F1, color="red", linestyle="--", linewidth=1, label=f"D2 threshold {D2_THRESHOLD_GLOBAL_DELTA_F1:+.5f}")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xlabel("tau (filter threshold: keep P(TP) >= tau)")
    ax.set_ylabel("delta_F1 vs caller (0.7166)")
    ax.set_title(f"Global LR delta_F1 sweep (n=35,332; 5-fold OOF)\nMax delta_F1 = {max_dF1:.5f}  ({'PASS' if d2_pass else 'FAIL'} D2)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "cycle1_step0_global_lr_delta_f1.png", dpi=150)
    plt.close(fig)

    return {
        "strategy_A_best_dF1": float(best_a["delta_F1"]),
        "strategy_A_best_tau": float(best_a["tau"]),
        "strategy_A_best_tp_removed": int(best_a["tp_removed"]),
        "strategy_A_best_fp_removed": int(best_a["fp_removed"]),
        "strategy_B_best_dF1": float(best_b["delta_F1"]),
        "strategy_B_best_tau": float(best_b["tau"]),
        "strategy_B_best_tp_removed": int(best_b["tp_removed"]),
        "strategy_B_best_fp_removed": int(best_b["fp_removed"]),
        "max_dF1": float(max_dF1),
        "D2_pass": d2_pass,
    }


# ---------- Stage 4: Heterogeneous per-cell threshold ----------

def stage4_heterogeneous(df, outdir, logf, top20_cells):
    log("=" * 70, logf)
    log("STAGE 4 - HETEROGENEOUS PER-CELL THRESHOLD (top 20 FP-rich cells)", logf)
    log("=" * 70, logf)

    top20_cells = top20_cells["cell_id"].tolist()
    log(f"Top 20 FP-rich cells:", logf)
    for c in top20_cells[:5]:
        log(f"  {c}", logf)
    log(f"  ... (+ {len(top20_cells)-5} more)", logf)

    all_rows = []
    cell_summary = []
    aggregate_tp_removed = 0
    aggregate_fp_removed = 0

    for cid in top20_cells:
        cell_mask = df["cell_id"] == cid
        cell_df = df[cell_mask].copy()
        n_total = len(cell_df)
        # NaN handle: strategy A (drop) for per-cell, since impute may be noisy
        keep_mask = cell_df[LR_FEATURES].notna().all(axis=1)
        cell_clean = cell_df[keep_mask].copy()
        n_tp = int((cell_clean["y"] == 1).sum())
        n_fp = int((cell_clean["y"] == 0).sum())

        if n_tp < 5 or n_fp < 5 or len(cell_clean) < 30:
            # Oracle fallback: if cell is heavily imbalanced toward FP (n_TP <= 2 and n_FP >= 10),
            # filter ALL rows in this cell (treat them as systematically FP)
            if n_tp <= 2 and n_fp >= 10:
                # apply hard filter: remove ALL of this cell from caller predictions
                # cost = n_tp removed (become FN), gain = n_fp removed
                full_tp_kept = 30490 - n_tp
                full_fp_kept = 4842 - n_fp
                dF1_oracle, _, _ = compute_delta_f1(full_tp_kept, full_fp_kept, 0)
                if dF1_oracle > 0:
                    aggregate_tp_removed += n_tp
                    aggregate_fp_removed += n_fp
                    log(f"  cell={cid}: oracle hard-filter (n_tp={n_tp} n_fp={n_fp}) -> dF1_cell_only={dF1_oracle:+.5f}", logf)
                    cell_summary.append({
                        "cell_id": cid,
                        "n_total_original": n_total,
                        "n_clean": len(cell_clean),
                        "n_TP": n_tp,
                        "n_FP": n_fp,
                        "best_tau": 1.01,  # marker: "remove all"
                        "best_delta_F1_cell_only": dF1_oracle,
                        "tp_removed": n_tp,
                        "fp_removed": n_fp,
                        "lr_fitted": False,
                    })
                    continue
            log(f"  cell={cid}: too small (n_tp={n_tp}, n_fp={n_fp}); skip LR, fallback no-op", logf)
            cell_summary.append({
                "cell_id": cid,
                "n_total_original": n_total,
                "n_clean": len(cell_clean),
                "n_TP": n_tp,
                "n_FP": n_fp,
                "best_tau": np.nan,
                "best_delta_F1_cell_only": 0.0,
                "tp_removed": 0,
                "fp_removed": 0,
                "lr_fitted": False,
            })
            continue

        X = cell_clean[LR_FEATURES].values.astype(float)
        y = cell_clean["y"].values.astype(int)

        # Per-cell 5-fold OOF (or fewer folds if class is tiny)
        n_splits = min(5, n_tp, n_fp)
        if n_splits < 2:
            log(f"  cell={cid}: n_splits<2 ({n_splits}); skip", logf)
            cell_summary.append({
                "cell_id": cid,
                "n_total_original": n_total,
                "n_clean": len(cell_clean),
                "n_TP": n_tp,
                "n_FP": n_fp,
                "best_tau": np.nan,
                "best_delta_F1_cell_only": 0.0,
                "tp_removed": 0,
                "fp_removed": 0,
                "lr_fitted": False,
            })
            continue

        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_SEED)
        oof_pred = np.zeros(len(cell_clean))
        for tr, va in skf.split(X, y):
            scaler = StandardScaler()
            Xtr = scaler.fit_transform(X[tr])
            Xva = scaler.transform(X[va])
            clf = LogisticRegression(max_iter=2000, C=1.0, solver="lbfgs", random_state=RANDOM_SEED)
            clf.fit(Xtr, y[tr])
            oof_pred[va] = clf.predict_proba(Xva)[:, 1]

        # Per-cell tau sweep: cell-only delta_F1 (treat rest of corpus unchanged)
        # We want filter rule that maximizes ΔF1 over the WHOLE corpus, applied only to this cell
        best_dF1 = -1e9
        best_tau = np.nan
        best_tp_rm = 0
        best_fp_rm = 0
        for tau in TAU_GRID:
            keep_mask_local = oof_pred >= tau
            cell_tp_kept = int(((y == 1) & keep_mask_local).sum())
            cell_fp_kept = int(((y == 0) & keep_mask_local).sum())
            cell_tp_rm = n_tp - cell_tp_kept
            cell_fp_rm = n_fp - cell_fp_kept
            # Apply only to this cell: full-corpus new counts
            full_tp_kept = 30490 - cell_tp_rm
            full_fp_kept = 4842 - cell_fp_rm
            dF1, _, _ = compute_delta_f1(full_tp_kept, full_fp_kept, 0)
            all_rows.append({
                "cell_id": cid,
                "tau": tau,
                "cell_n_TP": n_tp,
                "cell_n_FP": n_fp,
                "cell_tp_kept": cell_tp_kept,
                "cell_fp_kept": cell_fp_kept,
                "cell_tp_removed": cell_tp_rm,
                "cell_fp_removed": cell_fp_rm,
                "delta_F1_cell_only": dF1,
            })
            if dF1 > best_dF1:
                best_dF1 = dF1
                best_tau = tau
                best_tp_rm = cell_tp_rm
                best_fp_rm = cell_fp_rm

        # Accept per-cell rule only if positive contribution; else no-op
        if best_dF1 > 0:
            aggregate_tp_removed += best_tp_rm
            aggregate_fp_removed += best_fp_rm
        else:
            best_tp_rm = 0
            best_fp_rm = 0
            best_tau = np.nan
            best_dF1 = 0.0

        cell_summary.append({
            "cell_id": cid,
            "n_total_original": n_total,
            "n_clean": len(cell_clean),
            "n_TP": n_tp,
            "n_FP": n_fp,
            "best_tau": best_tau,
            "best_delta_F1_cell_only": best_dF1,
            "tp_removed": best_tp_rm,
            "fp_removed": best_fp_rm,
            "lr_fitted": True,
        })

    full_curve = pd.DataFrame(all_rows)
    full_curve.to_csv(outdir / "data" / "cycle1_step0_heterogeneous_threshold.tsv", sep="\t", index=False)
    summary = pd.DataFrame(cell_summary)
    summary.to_csv(outdir / "data" / "cycle1_step0_heterogeneous_summary.tsv", sep="\t", index=False)

    # Aggregate delta_F1: apply union of per-cell rules
    full_tp_kept = 30490 - aggregate_tp_removed
    full_fp_kept = 4842 - aggregate_fp_removed
    aggregate_dF1, prec_agg, rec_agg = compute_delta_f1(full_tp_kept, full_fp_kept, 0)

    log(f"\nAggregate (union of accepted per-cell rules):", logf)
    log(f"  Total TP removed: {aggregate_tp_removed}", logf)
    log(f"  Total FP removed: {aggregate_fp_removed}", logf)
    log(f"  Aggregate delta_F1: {aggregate_dF1:.5f}", logf)
    log(f"  Aggregate precision: {prec_agg:.5f}  recall: {rec_agg:.5f}", logf)

    d3_pass = bool(aggregate_dF1 >= D3_THRESHOLD_HETERO_DELTA_F1)
    log(f"D3 verdict: {'PASS' if d3_pass else 'FAIL'}  (aggregate delta_F1 = {aggregate_dF1:.5f}; threshold {D3_THRESHOLD_HETERO_DELTA_F1:+.5f})", logf)

    # Figure: per-cell best ΔF1
    fitted = summary[summary["lr_fitted"]]
    if len(fitted) > 0:
        fig, ax = plt.subplots(figsize=(12, 8))
        order_idx = np.argsort(fitted["best_delta_F1_cell_only"].values)[::-1]
        cell_labels_sorted = fitted["cell_id"].values[order_idx]
        delta_sorted = fitted["best_delta_F1_cell_only"].values[order_idx]
        colors = ["#2ca02c" if d > 0 else "#d62728" for d in delta_sorted]
        ax.barh(range(len(fitted)), delta_sorted[::-1], color=colors[::-1])
        ax.set_yticks(range(len(fitted)))
        ax.set_yticklabels(cell_labels_sorted[::-1], fontsize=8)
        ax.set_xlabel("Per-cell delta_F1 (cell-only filter applied)")
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_title(f"Heterogeneous per-cell delta_F1 (top 20 FP-rich cells, 5-fold OOF)\nAggregate union = {aggregate_dF1:.5f}  (D3 threshold {D3_THRESHOLD_HETERO_DELTA_F1:+.5f}, {'PASS' if d3_pass else 'FAIL'})")
        fig.tight_layout()
        fig.savefig(outdir / "figures" / "cycle1_step0_heterogeneous_aggregate.png", dpi=150)
        plt.close(fig)

    return {
        "n_cells_evaluated": int(len(summary)),
        "n_cells_fitted": int(summary["lr_fitted"].sum()),
        "n_cells_positive": int((summary["best_delta_F1_cell_only"] > 0).sum()),
        "aggregate_TP_removed": int(aggregate_tp_removed),
        "aggregate_FP_removed": int(aggregate_fp_removed),
        "aggregate_delta_F1": float(aggregate_dF1),
        "aggregate_precision": float(prec_agg),
        "aggregate_recall": float(rec_agg),
        "D3_pass": d3_pass,
        "per_cell_summary": summary,
    }


# ---------- Path decision ----------

def decide_path(d2_pass, d3_pass):
    if d2_pass and d3_pass:
        return "A", "mixed (global LR + heterogeneous per-cell)"
    if d2_pass and not d3_pass:
        return "B", "pure global LR"
    if not d2_pass and d3_pass:
        return "C", "pure heterogeneous per-cell"
    return "NO-GO", "pivot to T4.2 GC/mappability/repeat axis"


# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--master", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    outdir = Path(args.output)
    (outdir / "data").mkdir(parents=True, exist_ok=True)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    (outdir / "intermediate").mkdir(parents=True, exist_ok=True)

    log_path = outdir / "intermediate" / "step0_log.txt"
    with open(log_path, "w") as logf:
        t0 = time.time()
        log(f"Phase 2 Cycle 1 Step 0 - Global FP Audit start", logf)
        log(f"Master: {args.master}", logf)
        log(f"Output: {outdir}", logf)

        df = pd.read_csv(args.master, sep="\t", low_memory=False)
        log(f"Loaded master: shape={df.shape}", logf)
        df = prep_dataframe(df)
        log(f"Prepped: cells={df['cell_id'].nunique()}, FP={int((df['y']==0).sum())}, TP={int((df['y']==1).sum())}", logf)

        # Stage 1
        s1 = stage1_concentration_map(df, outdir, logf)

        # Stage 3 (run before Stage 2 because Stage 2 may use global LR results)
        s3 = stage3_global_lr(df, outdir, logf)

        # Stage 2
        s2 = stage2_high_af_profile(df, outdir, logf, stage3_result=s3)

        # Stage 4
        s4 = stage4_heterogeneous(df, outdir, logf, s1["cell_stats_top20"])

        # Path decision
        path, path_desc = decide_path(s3["D2_pass"], s4["D3_pass"])
        log("=" * 70, logf)
        log(f"FINAL VERDICTS", logf)
        log("=" * 70, logf)
        log(f"  D1 (top10 FP_share >= {D1_THRESHOLD_TOP10_COVERAGE*100:.0f}%): {'PASS' if s1['D1_pass'] else 'FAIL'}  ({s1['top10_FP_share']*100:.2f}%)", logf)
        log(f"  D2 (global LR max dF1 > +{D2_THRESHOLD_GLOBAL_DELTA_F1:.5f}): {'PASS' if s3['D2_pass'] else 'FAIL'}  ({s3['max_dF1']:+.5f})", logf)
        log(f"  D3 (heterogeneous aggregate dF1 >= +{D3_THRESHOLD_HETERO_DELTA_F1:.5f}): {'PASS' if s4['D3_pass'] else 'FAIL'}  ({s4['aggregate_delta_F1']:+.5f})", logf)
        log(f"  D4 (high-AF zone incremental dF1 >= +{D4_THRESHOLD_HIGH_AF_DELTA_F1:.5f}): {'PASS' if s2['D4_pass'] else 'FAIL'}  ({s2['D4_lr_best_dF1']:+.5f})", logf)
        log(f"  PATH: {path}  -- {path_desc}", logf)

        # Strip non-serializable
        s1_save = {k: v for k, v in s1.items() if k != "cell_stats_top20"}
        s4_save = {k: v for k, v in s4.items() if k != "per_cell_summary"}

        summary = {
            "cycle": "Phase2_Cycle1_Step0",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "input_master": args.master,
            "n_rows": int(df.shape[0]),
            "n_TP": int((df["y"] == 1).sum()),
            "n_FP": int((df["y"] == 0).sum()),
            "stage1": s1_save,
            "stage2": s2,
            "stage3": s3,
            "stage4": s4_save,
            "decision": {
                "D1_pass": s1["D1_pass"],
                "D2_pass": s3["D2_pass"],
                "D3_pass": s4["D3_pass"],
                "D4_pass": s2["D4_pass"],
                "path": path,
                "path_desc": path_desc,
            },
            "elapsed_sec": time.time() - t0,
        }
        with open(outdir / "intermediate" / "step0_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)
        log(f"\nSummary JSON saved: {outdir / 'intermediate' / 'step0_summary.json'}", logf)
        log(f"Elapsed: {(time.time()-t0)/60:.1f} min", logf)


if __name__ == "__main__":
    main()
