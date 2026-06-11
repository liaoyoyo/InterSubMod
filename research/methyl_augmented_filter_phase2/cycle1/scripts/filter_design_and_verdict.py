#!/usr/bin/env python3
"""
Phase 2 Cycle 1 Step 1 + Step 2 — Filter Design (Path B) + HCC1395 ΔF1 verdict
================================================================================

Step 1 — Filter design (Path B: pure global LR)
  Stage 1: VIF audit + L2 regularization sweep on 11 features
           - Diagnose Coverage_Multiple_imp / NumReads_master collinearity (R-Step0-2)
           - Pick L2 C value that stabilizes coefs (no sign flip across folds)
           - Lasso (L1) comparison: identify zero-out features
  Stage 2: NaN mechanism (R-Step0-1, 8.7x Strategy A vs B gap)
           - For NME_imbalance / Epipoly_Delta / HPFineF NaN rows:
             test MNAR (NaN ~ caller_af bin) vs MCAR (NaN random)
           - Logistic regression NaN_indicator ~ caller_af bin + cov bin
  Stage 3: Final filter rule design
           - Apply selected C, refit LR with L2 (Strategy B median impute)
           - 5-fold StratifiedKFold OOF, tau sweep
           - Output filter rule JSON

Step 2 — HCC1395 ΔF1 verdict
  Stage 1: Apply final filter to HCC1395 master TSV (n=35,332)
           - 5-fold OOF ΔF1
           - Compare with v1.0 +0.00242 baseline + Step 0 (+0.02637)
  Stage 2: Step 5c lost TP cross-check (R-Step0-4)
           - For 21 Step 5c lost TP region_ids, compute Cycle 1 filter P(TP)
           - rescue rate = fraction with P(TP) >= tau*
           - Verify global LR is not disproportionately removing low-AF subclone TP
  Stage 3: Robustness audit (R-Step0-5)
           - Multi-seed (5 seeds: 42,7,13,2026,1395) 5-fold CV variance
           - mean ± std ΔF1
  Stage 4: Verdict H_C1_2 / H_C1_3 / H_C1_4

Output:
  cycle1_track_a_filter.json
  cycle1_step1_filter_design.md (table fragments)
  cycle1_track_a_findings.md (table fragments)
  figures/
    cycle1_step1_vif_audit.png
    cycle1_step1_l2_regularization_sweep.png
    cycle1_step1_nan_mechanism.png
    cycle1_step2_lost_tp_predictions.png
    cycle1_step2_multiseed_variance.png
  intermediate/step1_step2_log.txt
"""

import argparse
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

# CJK font fallback
matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

# ---------- Constants (mirror Step 0) ----------

CALLER_F1_BASELINE = 0.7166
FN_CALLER = 19288  # from Step 0

# Filter feature set (Path B 11 features)
LR_FEATURES = [
    "V6_off_NG",
    "caller_af",
    "NumReads_master",
    "loh_inner_flag",
    "Coverage_Multiple_imp",
    "V6_off_meth_HPMergedDelta",
    "V6_off_meth_HPFineF",
    "V6_off_meth_NME_imbalance",
    "V6_off_meth_Epipoly_Delta",
    "V6_off_meth_ClusterPermanovaF",
    "chr8_flag",
]

# Methylation features with high NaN rate (audit MNAR/MCAR)
NAN_FEATURES = [
    "V6_off_meth_NME_imbalance",
    "V6_off_meth_Epipoly_Delta",
    "V6_off_meth_HPFineF",
]

# AF bins for NaN-mechanism stratification
AF_BINS = [0.0, 0.10, 0.20, 0.30, 0.50, 0.70, 1.01]
AF_LABELS = ["[0.00,0.10)", "[0.10,0.20)", "[0.20,0.30)", "[0.30,0.50)", "[0.50,0.70)", "[0.70,1.00]"]
COV_BINS = [-np.inf, 0.5, 0.9, 1.1, 1.3, np.inf]
COV_LABELS = ["<0.5", "[0.5,0.9)", "[0.9,1.1)", "[1.1,1.3)", ">=1.3"]

TAU_GRID = np.round(np.arange(0.10, 0.96, 0.01), 2)

C_VALUES = [0.001, 0.01, 0.1, 1.0, 10.0]  # L2 sweep

PRIMARY_SEED = 42
MULTI_SEEDS = [42, 7, 13, 2026, 1395]

# Step 5c lost TP region_ids (21 regions; see step5c_lost_tp_features.tsv)
STEP5C_LOST_TP_REGIONS = [
    "chr8_78174671_78184671", "chr8_82021565_82031565", "chr8_82067360_82077360",
    "chr8_82858123_82868123", "chr8_85434106_85444106", "chr8_85907839_85917839",
    "chr8_88286102_88296102", "chr8_122037156_122047156", "chr8_81782436_81792436",
    "chr8_82070646_82080646", "chr8_82571703_82581703", "chr8_82600457_82610457",
    "chr8_82838318_82848318", "chr8_82946427_82956427", "chr8_83076109_83086109",
    "chr8_83123414_83133414", "chr8_88419839_88429839", "chr8_94091239_94101239",
    "chr14_106265441_106275441", "chr8_87696219_87706219", "chr9_42161476_42171476",
]


# ---------- Utilities ----------

def log(msg, logf=None):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    if logf is not None:
        logf.write(line + "\n")
        logf.flush()


def compute_delta_f1(tp_kept, fp_kept):
    tp_total = 30490
    fp_total = 4842
    tp_removed = tp_total - tp_kept
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
    df = df.copy()
    df["loh_inner_flag"] = (df["loh_side"] == "Inner").astype(int)
    cov_median = df["Coverage_Multiple"].median()
    df["Coverage_Multiple_imp"] = df["Coverage_Multiple"].fillna(cov_median)
    df["chr8_flag"] = (df["chr"] == "chr8").astype(int)
    df["af_bin"] = pd.cut(df["caller_af"], bins=AF_BINS, labels=AF_LABELS,
                          include_lowest=True, right=False)
    df["y"] = (df["label"] == "TP").astype(int)
    return df


def impute_strategy_b(df, features):
    sub = df.copy()
    medians = {}
    for c in features:
        if sub[c].isna().any():
            med = sub[c].median()
            sub[c] = sub[c].fillna(med)
            medians[c] = float(med)
        else:
            medians[c] = float(sub[c].median()) if sub[c].notna().any() else 0.0
    return sub, medians


# ---------- Step 1 Stage 1: VIF + L2 sweep ----------

def compute_vif(X, feature_names):
    """
    Variance Inflation Factor: VIF_i = 1 / (1 - R^2_i)
    where R^2_i is R-squared regressing feature i on all other features.
    """
    vifs = []
    X = np.asarray(X, dtype=float)
    for i, name in enumerate(feature_names):
        y_i = X[:, i]
        X_others = np.delete(X, i, axis=1)
        # Add constant
        X_aug = np.column_stack([np.ones(len(X_others)), X_others])
        try:
            beta, *_ = np.linalg.lstsq(X_aug, y_i, rcond=None)
            y_hat = X_aug @ beta
            ss_res = np.sum((y_i - y_hat) ** 2)
            ss_tot = np.sum((y_i - y_i.mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
            vif = 1.0 / (1.0 - r2) if r2 < 1 - 1e-12 else float("inf")
        except Exception:
            vif = float("nan")
        vifs.append((name, r2, vif))
    return vifs


def stage1_vif_l2(df, outdir, logf):
    log("=" * 70, logf)
    log("STEP 1 STAGE 1 - VIF + L2 regularization sweep", logf)
    log("=" * 70, logf)

    sub, medians = impute_strategy_b(df, LR_FEATURES)
    X = sub[LR_FEATURES].values.astype(float)
    y = sub["y"].values.astype(int)

    # VIF audit
    scaler = StandardScaler()
    X_std = scaler.fit_transform(X)
    vifs = compute_vif(X_std, LR_FEATURES)
    log("VIF audit (standardized features):", logf)
    vif_table = []
    for name, r2, vif in vifs:
        log(f"  {name:35s}  R^2={r2:.4f}  VIF={vif:.3f}", logf)
        vif_table.append({"feature": name, "R2": r2, "VIF": vif})
    vif_df = pd.DataFrame(vif_table)
    vif_df.to_csv(outdir / "data" / "cycle1_step1_vif.tsv", sep="\t", index=False)

    high_vif = vif_df[vif_df["VIF"] > 5.0]
    log(f"\nFeatures with VIF > 5 (collinearity concern): {len(high_vif)}", logf)
    for _, row in high_vif.iterrows():
        log(f"  {row['feature']}: VIF={row['VIF']:.2f}", logf)

    # L2 sweep: track coefs vs C
    log("\nL2 regularization sweep (C in [0.001, 0.01, 0.1, 1.0, 10.0]):", logf)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=PRIMARY_SEED)
    coef_records = []
    df1_records = []
    for C in C_VALUES:
        per_fold_coefs = []
        oof_pred = np.zeros(len(sub))
        for fold_idx, (tr, va) in enumerate(skf.split(X, y)):
            sc = StandardScaler()
            Xtr = sc.fit_transform(X[tr])
            Xva = sc.transform(X[va])
            clf = LogisticRegression(max_iter=5000, C=C, solver="lbfgs",
                                     penalty="l2", random_state=PRIMARY_SEED)
            clf.fit(Xtr, y[tr])
            oof_pred[va] = clf.predict_proba(Xva)[:, 1]
            per_fold_coefs.append(clf.coef_[0].tolist())
        per_fold_coefs = np.array(per_fold_coefs)  # 5 x 11
        mean_coef = per_fold_coefs.mean(axis=0)
        std_coef = per_fold_coefs.std(axis=0)
        # Track sign-flip across folds
        sign_flips = []
        for j in range(per_fold_coefs.shape[1]):
            sgn = np.sign(per_fold_coefs[:, j])
            sign_flips.append(int(len(set(sgn[sgn != 0])) > 1))
        n_flips = int(sum(sign_flips))
        for j, fname in enumerate(LR_FEATURES):
            coef_records.append({
                "C": C,
                "feature": fname,
                "mean_coef": float(mean_coef[j]),
                "std_coef": float(std_coef[j]),
                "sign_flip": int(sign_flips[j]),
            })
        # Best tau ΔF1 for this C
        best_dF1 = -1e9
        best_tau = 0.0
        for tau in TAU_GRID:
            keep = oof_pred >= tau
            tp_kept = int(((y == 1) & keep).sum())
            fp_kept = int(((y == 0) & keep).sum())
            dF1, _, _ = compute_delta_f1(tp_kept, fp_kept)
            if dF1 > best_dF1:
                best_dF1 = dF1
                best_tau = float(tau)
        df1_records.append({
            "C": C,
            "best_dF1": float(best_dF1),
            "best_tau": float(best_tau),
            "n_sign_flips": n_flips,
            "max_abs_coef": float(np.max(np.abs(mean_coef))),
        })
        log(f"  C={C}: best ΔF1={best_dF1:+.5f} @ tau={best_tau}; n_sign_flips={n_flips}; max|coef|={np.max(np.abs(mean_coef)):.3f}", logf)

    coef_df = pd.DataFrame(coef_records)
    coef_df.to_csv(outdir / "data" / "cycle1_step1_l2_coefs.tsv", sep="\t", index=False)
    df1_df = pd.DataFrame(df1_records)
    df1_df.to_csv(outdir / "data" / "cycle1_step1_l2_df1.tsv", sep="\t", index=False)

    # L1 (Lasso) sanity comparison at canonical strength
    log("\nL1 (Lasso) sanity at C=0.1:", logf)
    skf2 = StratifiedKFold(n_splits=5, shuffle=True, random_state=PRIMARY_SEED)
    per_fold_l1 = []
    for tr, va in skf2.split(X, y):
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr])
        clf = LogisticRegression(max_iter=5000, C=0.1, solver="liblinear",
                                 penalty="l1", random_state=PRIMARY_SEED)
        clf.fit(Xtr, y[tr])
        per_fold_l1.append(clf.coef_[0].tolist())
    per_fold_l1 = np.array(per_fold_l1)
    mean_l1 = per_fold_l1.mean(axis=0)
    zero_fraction = (np.abs(per_fold_l1) < 1e-6).mean(axis=0)
    l1_records = []
    for j, fname in enumerate(LR_FEATURES):
        l1_records.append({
            "feature": fname,
            "mean_coef_l1": float(mean_l1[j]),
            "zero_fraction_across_folds": float(zero_fraction[j]),
        })
        log(f"  {fname:35s}  mean_coef={mean_l1[j]:+.4f}  zero_frac={zero_fraction[j]:.2f}", logf)
    l1_df = pd.DataFrame(l1_records)
    l1_df.to_csv(outdir / "data" / "cycle1_step1_l1_coefs.tsv", sep="\t", index=False)

    # Decision: pick C that minimizes sign-flips AND keeps best_dF1 within 0.5% of max (broad plateau)
    df1_df_sorted = df1_df.sort_values("best_dF1", ascending=False)
    max_dF1 = df1_df_sorted["best_dF1"].max()
    candidates = df1_df_sorted[
        (df1_df_sorted["best_dF1"] >= max_dF1 - 0.001) &
        (df1_df_sorted["n_sign_flips"] == 0)
    ]
    if len(candidates) == 0:
        candidates = df1_df_sorted[df1_df_sorted["best_dF1"] >= max_dF1 - 0.001]
    # Prefer smaller |coef| (more regularized) within candidate plateau
    chosen = candidates.sort_values("max_abs_coef").iloc[0]
    chosen_C = float(chosen["C"])
    chosen_tau = float(chosen["best_tau"])
    log(f"\nChosen L2 C={chosen_C} (best ΔF1={chosen['best_dF1']:+.5f} @ tau={chosen_tau}; max|coef|={chosen['max_abs_coef']:.3f}; sign_flips={chosen['n_sign_flips']})", logf)

    # Figure: VIF bar chart
    fig, ax = plt.subplots(figsize=(10, 5))
    order = vif_df.sort_values("VIF", ascending=False)
    bars = ax.bar(range(len(order)), order["VIF"].clip(upper=50), color=["#d62728" if v > 5 else "#1f77b4" for v in order["VIF"]])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order["feature"], rotation=45, ha="right", fontsize=8)
    ax.axhline(5, color="red", linestyle="--", linewidth=1, label="VIF=5 threshold")
    ax.set_ylabel("VIF (clip at 50)")
    ax.set_title("VIF audit — collinearity diagnosis (Step 1 Stage 1)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "cycle1_step1_vif_audit.png", dpi=150)
    plt.close(fig)

    # Figure: L2 sweep coefs heatmap + best ΔF1
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    # left: heatmap of mean_coef vs C
    pivot = coef_df.pivot(index="feature", columns="C", values="mean_coef")
    pivot = pivot.reindex(LR_FEATURES)
    im = axes[0].imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                        vmin=-max(abs(pivot.values.min()), abs(pivot.values.max())),
                        vmax=max(abs(pivot.values.min()), abs(pivot.values.max())))
    axes[0].set_yticks(range(len(LR_FEATURES)))
    axes[0].set_yticklabels(LR_FEATURES, fontsize=8)
    axes[0].set_xticks(range(len(C_VALUES)))
    axes[0].set_xticklabels(C_VALUES)
    axes[0].set_xlabel("L2 C (smaller = more regularization)")
    axes[0].set_title("Mean LR coefficient vs C (Strategy B impute, 5-fold)")
    plt.colorbar(im, ax=axes[0])
    # annotate values
    for i in range(len(LR_FEATURES)):
        for j in range(len(C_VALUES)):
            val = pivot.values[i, j]
            color = "white" if abs(val) > 3 else "black"
            axes[0].text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, color=color)

    # right: best_dF1 vs C
    axes[1].plot(df1_df["C"], df1_df["best_dF1"], "-o", color="#1f77b4")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("L2 C")
    axes[1].set_ylabel("Best ΔF1 across tau grid")
    axes[1].axhline(0.00242, color="gray", linestyle="--", linewidth=1, label="v1.0 baseline +0.00242")
    axes[1].axhline(0.01, color="red", linestyle="--", linewidth=1, label="H_C1_3 threshold +0.01")
    axes[1].set_title("Best ΔF1 vs L2 C")
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    for _, r in df1_df.iterrows():
        axes[1].annotate(f"flips={int(r['n_sign_flips'])}",
                         xy=(r["C"], r["best_dF1"]), fontsize=7,
                         xytext=(5, -10), textcoords="offset points")
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "cycle1_step1_l2_regularization_sweep.png", dpi=150)
    plt.close(fig)

    return {
        "vif_table": vif_table,
        "n_high_vif": int(len(high_vif)),
        "l2_sweep": df1_records,
        "chosen_C": chosen_C,
        "chosen_tau": chosen_tau,
        "chosen_dF1": float(chosen["best_dF1"]),
        "max_abs_coef_chosen": float(chosen["max_abs_coef"]),
        "sign_flips_chosen": int(chosen["n_sign_flips"]),
        "l1_records": l1_records,
        "medians": medians,
    }


# ---------- Step 1 Stage 2: NaN mechanism ----------

def stage2_nan_mechanism(df, outdir, logf):
    log("=" * 70, logf)
    log("STEP 1 STAGE 2 - NaN mechanism (MCAR vs MNAR)", logf)
    log("=" * 70, logf)

    rows = []
    for feat in NAN_FEATURES:
        is_nan = df[feat].isna().astype(int)
        nan_rate = float(is_nan.mean())
        log(f"\nFeature: {feat}  NaN rate = {nan_rate*100:.1f}%", logf)
        # By AF bin
        by_af = df.groupby("af_bin", observed=False)["__nan__"].agg(["mean", "count"]) if False else None
        af_grp = df.assign(__nan__=is_nan).groupby("af_bin", observed=False).agg(
            nan_rate=("__nan__", "mean"),
            n=("__nan__", "size"),
        ).reset_index()
        log("  AF-bin NaN rates:", logf)
        for _, r in af_grp.iterrows():
            log(f"    {str(r['af_bin']):20s}  nan_rate={r['nan_rate']*100:6.2f}%  n={int(r['n'])}", logf)

        # Chi-square: NaN x af_bin
        ct = pd.crosstab(df["af_bin"], is_nan)
        if ct.shape[0] >= 2 and ct.shape[1] == 2:
            chi2, p_chi, dof, _ = stats.chi2_contingency(ct.values)
        else:
            chi2, p_chi, dof = np.nan, np.nan, np.nan

        # Mann-Whitney: caller_af for NaN==1 vs NaN==0
        af1 = df.loc[is_nan == 1, "caller_af"].dropna().values
        af0 = df.loc[is_nan == 0, "caller_af"].dropna().values
        if len(af1) >= 5 and len(af0) >= 5:
            mw_u, mw_p = stats.mannwhitneyu(af1, af0, alternative="two-sided")
            med_af_nan = float(np.median(af1))
            med_af_nonnan = float(np.median(af0))
        else:
            mw_u, mw_p = np.nan, np.nan
            med_af_nan, med_af_nonnan = np.nan, np.nan

        # Effect size: range of nan_rate across AF bins
        nan_rate_range = float(af_grp["nan_rate"].max() - af_grp["nan_rate"].min())

        # Logistic regression: NaN ~ caller_af
        try:
            sc = StandardScaler()
            xa = sc.fit_transform(df["caller_af"].fillna(df["caller_af"].median()).values.reshape(-1, 1))
            clf_n = LogisticRegression(max_iter=5000, random_state=PRIMARY_SEED)
            clf_n.fit(xa, is_nan)
            logit_coef = float(clf_n.coef_[0][0])
            # Pseudo R^2 (McFadden-ish via log-likelihood)
            ll_full = -np.log(clf_n.predict_proba(xa)[range(len(is_nan)), is_nan]).sum()
            p_null = is_nan.mean()
            ll_null = -(is_nan * np.log(max(p_null, 1e-9)) + (1 - is_nan) * np.log(max(1 - p_null, 1e-9))).sum()
            mcfadden = 1 - (ll_full / ll_null) if ll_null > 0 else np.nan
        except Exception as e:
            logit_coef, mcfadden = np.nan, np.nan

        # Verdict: MNAR if (chi-square p<0.01 AND nan_rate_range>0.10) OR (MW p<0.01 AND |median diff|>0.05)
        is_mnar = (
            (not np.isnan(p_chi) and p_chi < 0.01 and nan_rate_range > 0.10)
            or (not np.isnan(mw_p) and mw_p < 0.01 and abs(med_af_nan - med_af_nonnan) > 0.05)
        )

        log(f"  Chi-square NaN x af_bin: p={p_chi:.3e} dof={dof}", logf)
        log(f"  Mann-Whitney caller_af | NaN vs non-NaN: p={mw_p:.3e}  median NaN={med_af_nan:.3f}  non-NaN={med_af_nonnan:.3f}", logf)
        log(f"  Logit NaN ~ caller_af: coef={logit_coef:+.3f}  pseudo-R^2={mcfadden:.4f}", logf)
        log(f"  Verdict: {'MNAR' if is_mnar else 'MCAR'} (effect size: range={nan_rate_range:.3f})", logf)

        rows.append({
            "feature": feat,
            "nan_rate": nan_rate,
            "chi2_p": float(p_chi) if not np.isnan(p_chi) else None,
            "mw_p": float(mw_p) if not np.isnan(mw_p) else None,
            "median_af_nan": med_af_nan,
            "median_af_nonnan": med_af_nonnan,
            "logit_coef": logit_coef if not np.isnan(logit_coef) else None,
            "mcfadden_pseudo_r2": mcfadden if not np.isnan(mcfadden) else None,
            "nan_rate_range_across_af_bins": nan_rate_range,
            "verdict": "MNAR" if is_mnar else "MCAR",
        })

    nan_df = pd.DataFrame(rows)
    nan_df.to_csv(outdir / "data" / "cycle1_step1_nan_mechanism.tsv", sep="\t", index=False)

    # Figure: NaN rate by AF bin per feature
    fig, axes = plt.subplots(1, len(NAN_FEATURES), figsize=(5 * len(NAN_FEATURES), 5), sharey=True)
    if len(NAN_FEATURES) == 1:
        axes = [axes]
    for ax, feat, rec in zip(axes, NAN_FEATURES, rows):
        is_nan = df[feat].isna().astype(int)
        af_grp = df.assign(__nan__=is_nan).groupby("af_bin", observed=False).agg(
            nan_rate=("__nan__", "mean"),
            n=("__nan__", "size"),
        ).reset_index()
        ax.bar(range(len(af_grp)), af_grp["nan_rate"].values * 100,
               color="#d62728" if rec["verdict"] == "MNAR" else "#2ca02c")
        ax.set_xticks(range(len(af_grp)))
        ax.set_xticklabels(af_grp["af_bin"], rotation=45, ha="right", fontsize=8)
        ax.set_title(f"{feat}\n{rec['verdict']} (range={rec['nan_rate_range_across_af_bins']:.2f}, p_chi={rec['chi2_p']:.1e})", fontsize=9)
        ax.set_ylabel("NaN rate (%)")
        ax.grid(alpha=0.3)
    fig.suptitle("NaN mechanism: NaN rate vs caller_af bin", fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "cycle1_step1_nan_mechanism.png", dpi=150)
    plt.close(fig)

    return rows


# ---------- Step 1 Stage 3: Final filter rule ----------

def stage3_final_filter(df, outdir, logf, chosen_C):
    log("=" * 70, logf)
    log(f"STEP 1 STAGE 3 - Final filter design (L2 C={chosen_C})", logf)
    log("=" * 70, logf)

    sub, medians = impute_strategy_b(df, LR_FEATURES)
    X = sub[LR_FEATURES].values.astype(float)
    y = sub["y"].values.astype(int)

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=PRIMARY_SEED)
    oof_pred = np.zeros(len(sub))
    per_fold_coefs = []
    per_fold_intercepts = []
    per_fold_scalers = []
    for fold_idx, (tr, va) in enumerate(skf.split(X, y)):
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr])
        Xva = sc.transform(X[va])
        clf = LogisticRegression(max_iter=5000, C=chosen_C, solver="lbfgs",
                                 penalty="l2", random_state=PRIMARY_SEED)
        clf.fit(Xtr, y[tr])
        oof_pred[va] = clf.predict_proba(Xva)[:, 1]
        per_fold_coefs.append(clf.coef_[0].tolist())
        per_fold_intercepts.append(float(clf.intercept_[0]))
        per_fold_scalers.append({
            "mean": sc.mean_.tolist(),
            "scale": sc.scale_.tolist(),
        })
    per_fold_coefs_arr = np.array(per_fold_coefs)

    # tau sweep
    rows = []
    for tau in TAU_GRID:
        keep = oof_pred >= tau
        tp_kept = int(((y == 1) & keep).sum())
        fp_kept = int(((y == 0) & keep).sum())
        dF1, prec, rec = compute_delta_f1(tp_kept, fp_kept)
        rows.append({
            "tau": float(tau),
            "tp_kept": tp_kept,
            "fp_kept": fp_kept,
            "tp_removed": 30490 - tp_kept,
            "fp_removed": 4842 - fp_kept,
            "delta_F1": float(dF1),
            "precision": float(prec),
            "recall": float(rec),
        })
    sweep = pd.DataFrame(rows)
    sweep.to_csv(outdir / "data" / "cycle1_step1_final_tau_sweep.tsv", sep="\t", index=False)
    best = sweep.loc[sweep["delta_F1"].idxmax()]
    log(f"Best tau: {best['tau']:.2f}  ΔF1={best['delta_F1']:+.5f}  TP_rm={int(best['tp_removed'])} FP_rm={int(best['fp_removed'])}", logf)

    # Save OOF predictions for Step 2
    oof_df = pd.DataFrame({
        "region_id": sub["region_id"].values,
        "label": sub["label"].values,
        "y": y,
        "caller_af": sub["caller_af"].values,
        "p_oof": oof_pred,
    })
    oof_df.to_csv(outdir / "data" / "cycle1_step1_oof_predictions.tsv", sep="\t", index=False)

    # Feature importance ranking
    mean_coef = per_fold_coefs_arr.mean(axis=0)
    importance = sorted(
        zip(LR_FEATURES, mean_coef.tolist(),
            per_fold_coefs_arr.std(axis=0).tolist()),
        key=lambda t: -abs(t[1]),
    )
    log("Feature importance (sorted by |mean coef|):", logf)
    for fname, mc, sd in importance:
        log(f"  {fname:35s}  mean_coef={mc:+.4f}  std_across_folds={sd:.4f}", logf)

    return {
        "tau_star": float(best["tau"]),
        "best_dF1": float(best["delta_F1"]),
        "best_tp_removed": int(best["tp_removed"]),
        "best_fp_removed": int(best["fp_removed"]),
        "per_fold_coefs": per_fold_coefs,
        "per_fold_intercepts": per_fold_intercepts,
        "per_fold_scalers": per_fold_scalers,
        "feature_order": LR_FEATURES,
        "medians": medians,
        "L2_C": float(chosen_C),
        "feature_importance": [
            {"feature": f, "mean_coef": mc, "std_coef": sd}
            for f, mc, sd in importance
        ],
        "oof_pred": oof_pred,
        "y_arr": y,
    }


# ---------- Step 2 Stage 1+2: Apply filter + lost TP cross-check ----------

def step2_apply_and_lost_tp(df, outdir, logf, filter_obj):
    log("=" * 70, logf)
    log("STEP 2 STAGE 1+2 - Apply filter + Step 5c lost TP cross-check", logf)
    log("=" * 70, logf)

    tau = filter_obj["tau_star"]
    oof_pred = filter_obj["oof_pred"]
    y = filter_obj["y_arr"]

    keep = oof_pred >= tau
    tp_kept = int(((y == 1) & keep).sum())
    fp_kept = int(((y == 0) & keep).sum())
    dF1, prec, rec = compute_delta_f1(tp_kept, fp_kept)
    log(f"Filter applied: tau={tau:.2f}", logf)
    log(f"  TP kept: {tp_kept} (removed {30490-tp_kept})", logf)
    log(f"  FP kept: {fp_kept} (removed {4842-fp_kept})", logf)
    log(f"  Precision: {prec:.4f}  Recall: {rec:.4f}  ΔF1: {dF1:+.5f}", logf)
    log(f"  vs v1.0 baseline +0.00242: ratio = {dF1/0.00242:.2f}x", logf)
    log(f"  vs H_C1_3 threshold +0.01: ratio = {dF1/0.01:.2f}x", logf)

    # Stage 2: Step 5c lost TP cross-check
    sub, _ = impute_strategy_b(df, LR_FEATURES)
    region_ids = sub["region_id"].values
    region_to_idx = {rid: i for i, rid in enumerate(region_ids)}

    lost_records = []
    in_master_count = 0
    rescued_count = 0
    for rid in STEP5C_LOST_TP_REGIONS:
        if rid in region_to_idx:
            in_master_count += 1
            i = region_to_idx[rid]
            p = float(oof_pred[i])
            kept_here = bool(p >= tau)
            if kept_here:
                rescued_count += 1
            lost_records.append({
                "region_id": rid,
                "in_master": True,
                "caller_af": float(sub["caller_af"].iloc[i]),
                "p_oof": p,
                "kept_by_cycle1_filter": kept_here,
            })
        else:
            lost_records.append({
                "region_id": rid,
                "in_master": False,
                "caller_af": None,
                "p_oof": None,
                "kept_by_cycle1_filter": None,
            })
    rescue_rate = rescued_count / max(in_master_count, 1)
    log(f"\nStep 5c 21 lost TP cross-check:", logf)
    log(f"  In master (out of 21): {in_master_count}", logf)
    log(f"  Rescued by cycle 1 filter (P(TP)>={tau:.2f}): {rescued_count}/{in_master_count} = {rescue_rate*100:.1f}%", logf)
    log(f"  Rescue verdict: {'rescues majority' if rescue_rate >= 0.5 else 'still removes low-AF subclone TP'}", logf)

    lost_df = pd.DataFrame(lost_records)
    lost_df.to_csv(outdir / "data" / "cycle1_step2_lost_tp_predictions.tsv", sep="\t", index=False)

    # Figure: P(TP) histogram for 21 lost TP vs all TP and FP
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(oof_pred[y == 1], bins=40, alpha=0.4, color="#2ca02c", label=f"All TP n={int((y==1).sum())}", density=True)
    ax.hist(oof_pred[y == 0], bins=40, alpha=0.4, color="#d62728", label=f"All FP n={int((y==0).sum())}", density=True)
    p_vals_lost = [r["p_oof"] for r in lost_records if r["in_master"]]
    if p_vals_lost:
        ax.scatter(p_vals_lost, [0.5] * len(p_vals_lost),
                   color="orange", s=60, marker="v", zorder=5,
                   label=f"21 Step 5c lost TP (n={in_master_count} in master)")
    ax.axvline(tau, color="black", linestyle="--", label=f"tau*={tau:.2f}")
    ax.set_xlabel("P(TP) from cycle 1 filter (OOF)")
    ax.set_ylabel("Density")
    ax.set_title(f"Step 5c lost TP cross-check\nRescue rate: {rescued_count}/{in_master_count} = {rescue_rate*100:.1f}%")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "cycle1_step2_lost_tp_predictions.png", dpi=150)
    plt.close(fig)

    return {
        "tau_star": tau,
        "tp_kept": tp_kept,
        "fp_kept": fp_kept,
        "tp_removed": 30490 - tp_kept,
        "fp_removed": 4842 - fp_kept,
        "precision": float(prec),
        "recall": float(rec),
        "delta_F1": float(dF1),
        "vs_v1_0_baseline_ratio": float(dF1 / 0.00242),
        "vs_H_C1_3_threshold_ratio": float(dF1 / 0.01),
        "step5c_lost_tp_in_master": int(in_master_count),
        "step5c_lost_tp_rescued": int(rescued_count),
        "step5c_rescue_rate": float(rescue_rate),
    }


# ---------- Step 2 Stage 3: Multi-seed robustness ----------

def step2_multiseed(df, outdir, logf, chosen_C):
    log("=" * 70, logf)
    log("STEP 2 STAGE 3 - Multi-seed robustness", logf)
    log("=" * 70, logf)

    sub, _ = impute_strategy_b(df, LR_FEATURES)
    X = sub[LR_FEATURES].values.astype(float)
    y = sub["y"].values.astype(int)

    records = []
    for seed in MULTI_SEEDS:
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
        oof_pred = np.zeros(len(sub))
        for tr, va in skf.split(X, y):
            sc = StandardScaler()
            Xtr = sc.fit_transform(X[tr])
            Xva = sc.transform(X[va])
            clf = LogisticRegression(max_iter=5000, C=chosen_C, solver="lbfgs",
                                     penalty="l2", random_state=seed)
            clf.fit(Xtr, y[tr])
            oof_pred[va] = clf.predict_proba(Xva)[:, 1]
        # Find best tau for this seed
        best_dF1 = -1e9
        best_tau = 0.0
        for tau in TAU_GRID:
            keep = oof_pred >= tau
            tp_kept = int(((y == 1) & keep).sum())
            fp_kept = int(((y == 0) & keep).sum())
            dF1, _, _ = compute_delta_f1(tp_kept, fp_kept)
            if dF1 > best_dF1:
                best_dF1 = dF1
                best_tau = float(tau)
        records.append({"seed": seed, "best_dF1": best_dF1, "best_tau": best_tau})
        log(f"  seed={seed}: best ΔF1={best_dF1:+.5f} @ tau={best_tau}", logf)

    rec_df = pd.DataFrame(records)
    rec_df.to_csv(outdir / "data" / "cycle1_step2_multiseed.tsv", sep="\t", index=False)
    mean_dF1 = float(rec_df["best_dF1"].mean())
    std_dF1 = float(rec_df["best_dF1"].std(ddof=1))
    log(f"\nMulti-seed (n={len(MULTI_SEEDS)}): ΔF1 mean={mean_dF1:+.5f} std={std_dF1:.5f}", logf)
    log(f"  Stability verdict: {'stable' if std_dF1 < 0.001 else 'marginal'} (threshold std<0.001)", logf)

    # Figure
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(range(len(records)), rec_df["best_dF1"], color="#1f77b4")
    ax.errorbar(range(len(records)), rec_df["best_dF1"],
                yerr=[0.0] * len(records), fmt="none")
    ax.axhline(mean_dF1, color="red", linestyle="--", label=f"mean={mean_dF1:+.5f}")
    ax.fill_between(range(len(records)), mean_dF1 - std_dF1, mean_dF1 + std_dF1,
                    alpha=0.2, color="red", label=f"±1σ (std={std_dF1:.5f})")
    ax.axhline(0.00242, color="gray", linestyle=":", label="v1.0 baseline +0.00242")
    ax.axhline(0.01, color="orange", linestyle=":", label="H_C1_3 threshold +0.01")
    ax.set_xticks(range(len(records)))
    ax.set_xticklabels([str(s) for s in MULTI_SEEDS])
    ax.set_xlabel("Random seed")
    ax.set_ylabel("Best ΔF1 (5-fold OOF)")
    ax.set_title(f"Multi-seed robustness (n={len(MULTI_SEEDS)} seeds)\nmean={mean_dF1:+.5f} ± {std_dF1:.5f}")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "cycle1_step2_multiseed_variance.png", dpi=150)
    plt.close(fig)

    return {
        "seeds": MULTI_SEEDS,
        "records": records,
        "mean_dF1": mean_dF1,
        "std_dF1": std_dF1,
        "stable": bool(std_dF1 < 0.001),
    }


# ---------- Step 2 Stage 4: H_C1_2/3/4 verdicts ----------

def step2_verdicts(filter_apply, multiseed_obj, df, outdir, logf, chosen_C):
    log("=" * 70, logf)
    log("STEP 2 STAGE 4 - Pre-registered hypothesis verdicts", logf)
    log("=" * 70, logf)

    dF1_primary = filter_apply["delta_F1"]
    dF1_multi_mean = multiseed_obj["mean_dF1"]

    H_C1_2_pass = bool(dF1_primary > 0.00242)
    H_C1_3_pass = bool(dF1_primary >= 0.01)

    # H_C1_4: high-AF zone incremental ΔF1
    sub, _ = impute_strategy_b(df, LR_FEATURES)
    high_af_mask = sub["caller_af"] > 0.3
    X = sub[LR_FEATURES].values.astype(float)
    y = sub["y"].values.astype(int)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=PRIMARY_SEED)
    oof_pred_full = np.zeros(len(sub))
    for tr, va in skf.split(X, y):
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr])
        Xva = sc.transform(X[va])
        clf = LogisticRegression(max_iter=5000, C=chosen_C, solver="lbfgs",
                                 penalty="l2", random_state=PRIMARY_SEED)
        clf.fit(Xtr, y[tr])
        oof_pred_full[va] = clf.predict_proba(Xva)[:, 1]
    # high-AF only filter at tau*
    tau = filter_apply["tau_star"]
    high_af_idx = np.where(high_af_mask.values)[0]
    high_af_y = y[high_af_idx]
    high_af_p = oof_pred_full[high_af_idx]
    keep_local = high_af_p >= tau
    sub_tp_kept = int(((high_af_y == 1) & keep_local).sum())
    sub_fp_kept = int(((high_af_y == 0) & keep_local).sum())
    sub_tp_total = int((high_af_y == 1).sum())
    sub_fp_total = int((high_af_y == 0).sum())
    # Outside high-AF unchanged
    out_tp = 30490 - sub_tp_total
    out_fp = 4842 - sub_fp_total
    full_tp_kept = out_tp + sub_tp_kept
    full_fp_kept = out_fp + sub_fp_kept
    dF1_highaf_only, _, _ = compute_delta_f1(full_tp_kept, full_fp_kept)
    H_C1_4_pass = bool(dF1_highaf_only >= 0.003)

    log(f"H_C1_2 (ΔF1 > +0.00242): observed {dF1_primary:+.5f} -> {'PASS' if H_C1_2_pass else 'FAIL'}", logf)
    log(f"H_C1_3 (ΔF1 >= +0.01):   observed {dF1_primary:+.5f} -> {'PASS' if H_C1_3_pass else 'FAIL'}", logf)
    log(f"H_C1_4 (high-AF zone incremental ΔF1 >= +0.003): observed {dF1_highaf_only:+.5f} -> {'PASS' if H_C1_4_pass else 'FAIL'}", logf)

    return {
        "H_C1_2": {"threshold": 0.00242, "observed": dF1_primary, "pass": H_C1_2_pass},
        "H_C1_3": {"threshold": 0.01, "observed": dF1_primary, "pass": H_C1_3_pass},
        "H_C1_4": {"threshold": 0.003, "observed": dF1_highaf_only, "pass": H_C1_4_pass,
                   "note": "high-AF (caller_af>0.3) zone-only filter incremental ΔF1"},
    }


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

    log_path = outdir / "intermediate" / "step1_step2_log.txt"
    with open(log_path, "w") as logf:
        t0 = time.time()
        log("Phase 2 Cycle 1 Step 1+2 start", logf)
        log(f"Master: {args.master}", logf)
        log(f"Output: {outdir}", logf)

        df = pd.read_csv(args.master, sep="\t", low_memory=False)
        log(f"Loaded master: shape={df.shape}", logf)
        df = prep_dataframe(df)
        log(f"Prepped: TP={int((df['y']==1).sum())}, FP={int((df['y']==0).sum())}", logf)

        s1 = stage1_vif_l2(df, outdir, logf)
        s2 = stage2_nan_mechanism(df, outdir, logf)
        chosen_C = s1["chosen_C"]
        s3 = stage3_final_filter(df, outdir, logf, chosen_C)

        # Save filter rule JSON
        filter_obj_save = {
            "cycle": "Phase2_Cycle1",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "path": "B (pure global LR)",
            "feature_set": LR_FEATURES,
            "n_features": len(LR_FEATURES),
            "L2_C": chosen_C,
            "tau_star": s3["tau_star"],
            "feature_medians_for_imputation": s3["medians"],
            "per_fold_coefs": s3["per_fold_coefs"],
            "per_fold_intercepts": s3["per_fold_intercepts"],
            "per_fold_scalers": s3["per_fold_scalers"],
            "feature_importance": s3["feature_importance"],
            "expected_delta_F1": s3["best_dF1"],
            "expected_tp_removed": s3["best_tp_removed"],
            "expected_fp_removed": s3["best_fp_removed"],
            "vif_table": s1["vif_table"],
            "n_high_vif_features": s1["n_high_vif"],
            "l2_sweep": s1["l2_sweep"],
            "nan_mechanism": s2,
            "step1_chain_position": "P3 PILOT (filter design) -> P4 generalize (Track B)",
        }
        with open(outdir / "cycle1_track_a_filter.json", "w") as f:
            json.dump(filter_obj_save, f, indent=2, default=str)
        log(f"\nFilter rule saved: {outdir/'cycle1_track_a_filter.json'}", logf)

        # Step 2
        filter_apply = step2_apply_and_lost_tp(df, outdir, logf, s3)
        multiseed = step2_multiseed(df, outdir, logf, chosen_C)
        verdicts = step2_verdicts(filter_apply, multiseed, df, outdir, logf, chosen_C)

        summary = {
            "cycle": "Phase2_Cycle1_Step1_Step2",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "step1_stage1_vif_l2": {k: v for k, v in s1.items() if k != "vif_table"},
            "step1_stage1_vif_table": s1["vif_table"],
            "step1_stage2_nan_mechanism": s2,
            "step1_stage3_filter": {
                "tau_star": s3["tau_star"],
                "best_dF1": s3["best_dF1"],
                "best_tp_removed": s3["best_tp_removed"],
                "best_fp_removed": s3["best_fp_removed"],
                "L2_C": s3["L2_C"],
                "feature_importance": s3["feature_importance"],
            },
            "step2_stage1_apply": filter_apply,
            "step2_stage3_multiseed": multiseed,
            "step2_stage4_verdicts": verdicts,
            "elapsed_sec": time.time() - t0,
        }
        with open(outdir / "intermediate" / "step1_step2_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)
        log(f"Summary JSON saved: {outdir/'intermediate'/'step1_step2_summary.json'}", logf)
        log(f"Elapsed: {(time.time()-t0)/60:.1f} min", logf)


if __name__ == "__main__":
    main()
