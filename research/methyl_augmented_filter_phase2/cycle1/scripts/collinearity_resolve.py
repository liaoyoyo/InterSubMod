#!/usr/bin/env python3
"""
Cycle 1 Step 1 Stage 1B — Resolve Coverage_Multiple_imp / NumReads_master collinearity.

VIF audit revealed VIF=215, 217 for the two -- severe collinearity.
Approach: Test 4 feature configurations and pick the one with best ΔF1 + tame |coef|.

Configurations:
  cfg_full  : all 11 features (baseline)
  cfg_drop_cov : drop Coverage_Multiple_imp, keep NumReads_master
  cfg_drop_nr  : drop NumReads_master, keep Coverage_Multiple_imp
  cfg_ratio    : replace both with reads_per_unit_cov = NumReads_master / Coverage_Multiple_imp
                 + keep one as anchor

For each: L2 sweep C in [0.01, 0.1, 1.0, 10.0], 5-fold OOF, best tau ΔF1
Pick best (cfg, C) such that:
  - ΔF1 within 1% of max
  - max|coef| minimized (interpretability + L2 not explosive)
  - All other VIFs <= 5
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
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

CALLER_F1_BASELINE = 0.7166
FN_CALLER = 19288
PRIMARY_SEED = 42
TAU_GRID = np.round(np.arange(0.10, 0.96, 0.01), 2)

BASE_FEATURES = [
    "V6_off_NG", "caller_af", "loh_inner_flag",
    "V6_off_meth_HPMergedDelta", "V6_off_meth_HPFineF",
    "V6_off_meth_NME_imbalance", "V6_off_meth_Epipoly_Delta",
    "V6_off_meth_ClusterPermanovaF", "chr8_flag",
]

AF_BINS = [0.0, 0.10, 0.20, 0.30, 0.50, 0.70, 1.01]
AF_LABELS = ["[0.00,0.10)", "[0.10,0.20)", "[0.20,0.30)", "[0.30,0.50)", "[0.50,0.70)", "[0.70,1.00]"]
COV_BINS = [-np.inf, 0.5, 0.9, 1.1, 1.3, np.inf]
COV_LABELS = ["<0.5", "[0.5,0.9)", "[0.9,1.1)", "[1.1,1.3)", ">=1.3"]


def log(msg, logf=None):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    if logf is not None:
        logf.write(line + "\n")
        logf.flush()


def compute_delta_f1(tp_kept, fp_kept):
    tp_removed = 30490 - tp_kept
    new_tp = tp_kept
    new_fp = fp_kept
    new_fn = FN_CALLER + tp_removed
    if new_tp + new_fp == 0 or new_tp + new_fn == 0:
        return 0.0
    precision = new_tp / (new_tp + new_fp)
    recall = new_tp / (new_tp + new_fn)
    if precision + recall == 0:
        return 0.0
    f1 = 2 * precision * recall / (precision + recall)
    return f1 - CALLER_F1_BASELINE


def compute_vif(X, names):
    out = []
    for i, name in enumerate(names):
        y_i = X[:, i]
        X_others = np.delete(X, i, axis=1)
        X_aug = np.column_stack([np.ones(len(X_others)), X_others])
        beta, *_ = np.linalg.lstsq(X_aug, y_i, rcond=None)
        y_hat = X_aug @ beta
        ss_res = np.sum((y_i - y_hat) ** 2)
        ss_tot = np.sum((y_i - y_i.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        vif = 1.0 / (1.0 - r2) if r2 < 1 - 1e-12 else float("inf")
        out.append((name, r2, vif))
    return out


def prep(df):
    df = df.copy()
    df["loh_inner_flag"] = (df["loh_side"] == "Inner").astype(int)
    cov_median = df["Coverage_Multiple"].median()
    df["Coverage_Multiple_imp"] = df["Coverage_Multiple"].fillna(cov_median)
    df["chr8_flag"] = (df["chr"] == "chr8").astype(int)
    df["y"] = (df["label"] == "TP").astype(int)
    df["af_bin"] = pd.cut(df["caller_af"], bins=AF_BINS, labels=AF_LABELS,
                          include_lowest=True, right=False)
    return df


def evaluate(df, feature_set, C_list, label, logf):
    sub = df.copy()
    for c in feature_set:
        if sub[c].isna().any():
            sub[c] = sub[c].fillna(sub[c].median())
    X = sub[feature_set].values.astype(float)
    y = sub["y"].values.astype(int)

    sc0 = StandardScaler()
    X_std = sc0.fit_transform(X)
    vifs = compute_vif(X_std, feature_set)
    max_vif = max(v for _, _, v in vifs if not np.isinf(v))

    results = []
    for C in C_list:
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=PRIMARY_SEED)
        oof_pred = np.zeros(len(sub))
        per_fold_coefs = []
        for tr, va in skf.split(X, y):
            sc = StandardScaler()
            Xtr = sc.fit_transform(X[tr])
            Xva = sc.transform(X[va])
            clf = LogisticRegression(max_iter=5000, C=C, solver="lbfgs",
                                     penalty="l2", random_state=PRIMARY_SEED)
            clf.fit(Xtr, y[tr])
            oof_pred[va] = clf.predict_proba(Xva)[:, 1]
            per_fold_coefs.append(clf.coef_[0].tolist())
        per_fold_coefs = np.array(per_fold_coefs)
        mean_coef = per_fold_coefs.mean(axis=0)
        max_abs = float(np.max(np.abs(mean_coef)))
        sign_flips = sum(
            int(len(set(np.sign(per_fold_coefs[:, j])[np.sign(per_fold_coefs[:, j]) != 0])) > 1)
            for j in range(per_fold_coefs.shape[1])
        )

        best_dF1 = -1e9
        best_tau = 0.0
        for tau in TAU_GRID:
            keep = oof_pred >= tau
            tp_kept = int(((y == 1) & keep).sum())
            fp_kept = int(((y == 0) & keep).sum())
            dF1 = compute_delta_f1(tp_kept, fp_kept)
            if dF1 > best_dF1:
                best_dF1 = dF1
                best_tau = float(tau)
        results.append({
            "config": label,
            "C": C,
            "best_dF1": float(best_dF1),
            "best_tau": float(best_tau),
            "max_abs_coef": max_abs,
            "n_sign_flips": int(sign_flips),
            "max_vif": float(max_vif),
            "n_features": len(feature_set),
        })
        log(f"  {label} C={C}: ΔF1={best_dF1:+.5f} @ tau={best_tau} | max|coef|={max_abs:.2f} | sign_flips={sign_flips} | max_VIF={max_vif:.2f}", logf)
    return results, vifs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--master", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    outdir = Path(args.output)
    (outdir / "data").mkdir(parents=True, exist_ok=True)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    (outdir / "intermediate").mkdir(parents=True, exist_ok=True)
    log_path = outdir / "intermediate" / "stage1b_collinearity_log.txt"

    with open(log_path, "w") as logf:
        log("Collinearity resolution — 4 config comparison", logf)
        df = pd.read_csv(args.master, sep="\t", low_memory=False)
        df = prep(df)
        log(f"Loaded: {df.shape}", logf)

        C_LIST = [0.01, 0.1, 1.0, 10.0]

        # cfg_full
        log("\n=== cfg_full (11 features) ===", logf)
        full_feat = BASE_FEATURES + ["NumReads_master", "Coverage_Multiple_imp"]
        res_full, vif_full = evaluate(df, full_feat, C_LIST, "cfg_full", logf)
        log(f"  VIFs: {[(n, round(v,1)) for n,_,v in vif_full]}", logf)

        # cfg_drop_cov
        log("\n=== cfg_drop_cov (drop Coverage_Multiple_imp; keep NumReads) ===", logf)
        feat_dropcov = BASE_FEATURES + ["NumReads_master"]
        res_dropcov, vif_dropcov = evaluate(df, feat_dropcov, C_LIST, "cfg_drop_cov", logf)
        log(f"  VIFs: {[(n, round(v,2)) for n,_,v in vif_dropcov]}", logf)

        # cfg_drop_nr
        log("\n=== cfg_drop_nr (drop NumReads_master; keep Coverage_Multiple_imp) ===", logf)
        feat_dropnr = BASE_FEATURES + ["Coverage_Multiple_imp"]
        res_dropnr, vif_dropnr = evaluate(df, feat_dropnr, C_LIST, "cfg_drop_nr", logf)
        log(f"  VIFs: {[(n, round(v,2)) for n,_,v in vif_dropnr]}", logf)

        # cfg_ratio
        log("\n=== cfg_ratio (replace both with reads/cov + keep Coverage as anchor) ===", logf)
        df["reads_per_unit_cov"] = df["NumReads_master"] / df["Coverage_Multiple_imp"].replace(0, np.nan)
        df["reads_per_unit_cov"] = df["reads_per_unit_cov"].fillna(df["reads_per_unit_cov"].median())
        feat_ratio = BASE_FEATURES + ["Coverage_Multiple_imp", "reads_per_unit_cov"]
        res_ratio, vif_ratio = evaluate(df, feat_ratio, C_LIST, "cfg_ratio", logf)
        log(f"  VIFs: {[(n, round(v,2)) for n,_,v in vif_ratio]}", logf)

        # Combine
        all_res = res_full + res_dropcov + res_dropnr + res_ratio
        df_res = pd.DataFrame(all_res)
        df_res.to_csv(outdir / "data" / "cycle1_step1_collinearity_comparison.tsv",
                      sep="\t", index=False)

        # Decision logic:
        # 1. Find config-C with best ΔF1
        # 2. Find config-C within 1% of best ΔF1 with smallest max|coef|, no inf VIF
        max_dF1 = df_res["best_dF1"].max()
        candidates = df_res[df_res["best_dF1"] >= max_dF1 - 0.002]
        log(f"\nCandidates within 0.002 of max ΔF1 ({max_dF1:+.5f}):", logf)
        for _, r in candidates.iterrows():
            log(f"  {r['config']}/C={r['C']}: ΔF1={r['best_dF1']:+.5f} max|coef|={r['max_abs_coef']:.2f} max_VIF={r['max_vif']:.1f}", logf)

        # Prefer: max_VIF < 5 (no severe collinearity)
        good_vif = candidates[candidates["max_vif"] < 5]
        if len(good_vif) > 0:
            chosen = good_vif.sort_values("max_abs_coef").iloc[0]
        else:
            # fallback: smallest max|coef| even if VIF high
            chosen = candidates.sort_values("max_abs_coef").iloc[0]

        log(f"\nChosen: config={chosen['config']} C={chosen['C']} ΔF1={chosen['best_dF1']:+.5f} max|coef|={chosen['max_abs_coef']:.2f} max_VIF={chosen['max_vif']:.1f}", logf)

        # Save decision summary
        summary = {
            "all_results": all_res,
            "chosen": chosen.to_dict(),
            "max_dF1_across_configs": float(max_dF1),
        }
        with open(outdir / "intermediate" / "stage1b_decision.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)
        log(f"\nDecision saved to {outdir/'intermediate'/'stage1b_decision.json'}", logf)

        # Figure: bar of best_dF1 across configs at C=1.0
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for cfg, color in zip(["cfg_full", "cfg_drop_cov", "cfg_drop_nr", "cfg_ratio"],
                              ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]):
            sub = df_res[df_res["config"] == cfg]
            axes[0].plot(sub["C"], sub["best_dF1"], "-o", label=cfg, color=color)
        axes[0].set_xscale("log")
        axes[0].axhline(0.00242, color="gray", linestyle="--", label="v1.0 baseline")
        axes[0].axhline(0.01, color="red", linestyle="--", label="H_C1_3 threshold")
        axes[0].set_xlabel("L2 C")
        axes[0].set_ylabel("Best ΔF1")
        axes[0].legend(fontsize=8)
        axes[0].set_title("ΔF1 vs L2 C (4 configurations)")
        axes[0].grid(alpha=0.3)

        for cfg, color in zip(["cfg_full", "cfg_drop_cov", "cfg_drop_nr", "cfg_ratio"],
                              ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]):
            sub = df_res[df_res["config"] == cfg]
            axes[1].plot(sub["C"], sub["max_abs_coef"], "-o", label=cfg, color=color)
        axes[1].set_xscale("log")
        axes[1].set_yscale("log")
        axes[1].set_xlabel("L2 C")
        axes[1].set_ylabel("max |coef| (log)")
        axes[1].legend(fontsize=8)
        axes[1].set_title("Coefficient magnitude vs L2 C")
        axes[1].grid(alpha=0.3, which="both")
        fig.tight_layout()
        fig.savefig(outdir / "figures" / "cycle1_step1_collinearity_comparison.png", dpi=150)
        plt.close(fig)
        log("Figure saved.", logf)


if __name__ == "__main__":
    main()
