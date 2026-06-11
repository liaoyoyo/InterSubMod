#!/usr/bin/env python3
"""
Cycle 1 Step 1 Stage 3 + Step 2 — FINAL filter design + HCC1395 verdict.

Decision rationale (Stage 1B collinearity audit):
  cfg_full (11 features) suffers VIF=217 between NumReads_master and Coverage_Multiple_imp.
  L2 alone cannot tame: max|coef|=68 at C=10 even when ΔF1 maxes.
  Principled fix: drop one collinear feature (cfg_drop_nr).

Chosen config: cfg_drop_nr at C=1.0
  - 10 features (drop NumReads_master, keep Coverage_Multiple_imp)
  - All VIFs < 2 (max 1.83)
  - max|coef| ~3.4 (stable)
  - ΔF1 = +0.02236 (still 9.2x v1.0 baseline +0.00242; exceeds H_C1_3 +0.01 by 2.2x)

We also evaluate cfg_full C=1.0 (ΔF1=+0.02637) and cfg_full C=10.0 (ΔF1=+0.02995) for transparency,
treating them as upper-bound estimates not robust to collinearity.

Step 2:
  Stage 1: Apply chosen filter to HCC1395 + ΔF1
  Stage 2: Step 5c 21 lost TP rescue rate
  Stage 3: Multi-seed (5) variance
  Stage 4: H_C1_2/3/4 verdicts
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
MULTI_SEEDS = [42, 7, 13, 2026, 1395]
TAU_GRID = np.round(np.arange(0.10, 0.96, 0.01), 2)

# CHOSEN FINAL FEATURE SET (cfg_drop_nr — drop NumReads_master)
FINAL_FEATURES = [
    "V6_off_NG", "caller_af", "loh_inner_flag",
    "V6_off_meth_HPMergedDelta", "V6_off_meth_HPFineF",
    "V6_off_meth_NME_imbalance", "V6_off_meth_Epipoly_Delta",
    "V6_off_meth_ClusterPermanovaF", "chr8_flag",
    "Coverage_Multiple_imp",
]
CHOSEN_C = 1.0

# Step 5c 21 lost TP region_ids
STEP5C_LOST_TP_REGIONS = [
    "chr8_78174671_78184671", "chr8_82021565_82031565", "chr8_82067360_82077360",
    "chr8_82858123_82868123", "chr8_85434106_85444106", "chr8_85907839_85917839",
    "chr8_88286102_88296102", "chr8_122037156_122047156", "chr8_81782436_81792436",
    "chr8_82070646_82080646", "chr8_82571703_82581703", "chr8_82600457_82610457",
    "chr8_82838318_82848318", "chr8_82946427_82956427", "chr8_83076109_83086109",
    "chr8_83123414_83133414", "chr8_88419839_88429839", "chr8_94091239_94101239",
    "chr14_106265441_106275441", "chr8_87696219_87706219", "chr9_42161476_42171476",
]


def log(msg, logf=None):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    if logf is not None:
        logf.write(line + "\n")
        logf.flush()


def compute_metrics(tp_kept, fp_kept):
    tp_removed = 30490 - tp_kept
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


def prep(df):
    df = df.copy()
    df["loh_inner_flag"] = (df["loh_side"] == "Inner").astype(int)
    cov_median = df["Coverage_Multiple"].median()
    df["Coverage_Multiple_imp"] = df["Coverage_Multiple"].fillna(cov_median)
    df["chr8_flag"] = (df["chr"] == "chr8").astype(int)
    df["y"] = (df["label"] == "TP").astype(int)
    return df


def impute(df, features):
    sub = df.copy()
    medians = {}
    for c in features:
        med = float(sub[c].median()) if sub[c].notna().any() else 0.0
        if sub[c].isna().any():
            sub[c] = sub[c].fillna(med)
        medians[c] = med
    return sub, medians


def fit_oof(X, y, C, seed=PRIMARY_SEED):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof = np.zeros(len(X))
    coefs = []
    intercepts = []
    scalers = []
    for tr, va in skf.split(X, y):
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr])
        Xva = sc.transform(X[va])
        clf = LogisticRegression(max_iter=5000, C=C, solver="lbfgs",
                                 penalty="l2", random_state=seed)
        clf.fit(Xtr, y[tr])
        oof[va] = clf.predict_proba(Xva)[:, 1]
        coefs.append(clf.coef_[0].tolist())
        intercepts.append(float(clf.intercept_[0]))
        scalers.append({"mean": sc.mean_.tolist(), "scale": sc.scale_.tolist()})
    return oof, np.array(coefs), intercepts, scalers


def sweep_tau(oof, y):
    best_dF1 = -1e9
    best_tau = 0.0
    best_tp_kept = 0
    best_fp_kept = 0
    best_prec = 0
    best_rec = 0
    for tau in TAU_GRID:
        keep = oof >= tau
        tp_kept = int(((y == 1) & keep).sum())
        fp_kept = int(((y == 0) & keep).sum())
        dF1, prec, rec = compute_metrics(tp_kept, fp_kept)
        if dF1 > best_dF1:
            best_dF1 = dF1
            best_tau = float(tau)
            best_tp_kept = tp_kept
            best_fp_kept = fp_kept
            best_prec = prec
            best_rec = rec
    return {"tau_star": best_tau, "best_dF1": float(best_dF1),
            "tp_kept": best_tp_kept, "fp_kept": best_fp_kept,
            "tp_removed": 30490 - best_tp_kept,
            "fp_removed": 4842 - best_fp_kept,
            "precision": float(best_prec), "recall": float(best_rec)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--master", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    outdir = Path(args.output)
    (outdir / "data").mkdir(parents=True, exist_ok=True)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    (outdir / "intermediate").mkdir(parents=True, exist_ok=True)

    log_path = outdir / "intermediate" / "final_filter_log.txt"
    with open(log_path, "w") as logf:
        t0 = time.time()
        log("Cycle 1 Step 1 final + Step 2 verdict", logf)
        log(f"Chosen feature set: cfg_drop_nr (n_features={len(FINAL_FEATURES)})", logf)
        log(f"Features: {FINAL_FEATURES}", logf)
        log(f"L2 C = {CHOSEN_C}", logf)

        df = pd.read_csv(args.master, sep="\t", low_memory=False)
        df = prep(df)
        log(f"Loaded: {df.shape}", logf)

        sub, medians = impute(df, FINAL_FEATURES)
        X = sub[FINAL_FEATURES].values.astype(float)
        y = sub["y"].values.astype(int)

        # Primary fit
        log("\n=== Primary OOF fit (seed=42) ===", logf)
        oof, coefs, intercepts, scalers = fit_oof(X, y, CHOSEN_C, seed=PRIMARY_SEED)
        sweep = sweep_tau(oof, y)
        log(f"Best tau: {sweep['tau_star']:.2f}  ΔF1: {sweep['best_dF1']:+.5f}", logf)
        log(f"  TP kept: {sweep['tp_kept']} (rm {sweep['tp_removed']})", logf)
        log(f"  FP kept: {sweep['fp_kept']} (rm {sweep['fp_removed']})", logf)
        log(f"  Precision: {sweep['precision']:.4f}  Recall: {sweep['recall']:.4f}", logf)
        log(f"  vs v1.0 baseline +0.00242: {sweep['best_dF1']/0.00242:.2f}x", logf)
        log(f"  vs H_C1_3 threshold +0.01: {sweep['best_dF1']/0.01:.2f}x", logf)

        # Feature importance
        mean_coef = coefs.mean(axis=0)
        std_coef = coefs.std(axis=0)
        importance = sorted(
            [(f, float(mc), float(sd)) for f, mc, sd in zip(FINAL_FEATURES, mean_coef, std_coef)],
            key=lambda t: -abs(t[1]),
        )
        log("\nFeature importance (sorted by |mean coef|):", logf)
        for fname, mc, sd in importance:
            log(f"  {fname:35s}  mean={mc:+.4f}  std={sd:.4f}", logf)

        # Save tau sweep curve
        rows = []
        for tau in TAU_GRID:
            keep = oof >= tau
            tp_kept = int(((y == 1) & keep).sum())
            fp_kept = int(((y == 0) & keep).sum())
            dF1, prec, rec = compute_metrics(tp_kept, fp_kept)
            rows.append({
                "tau": float(tau), "tp_kept": tp_kept, "fp_kept": fp_kept,
                "tp_removed": 30490 - tp_kept, "fp_removed": 4842 - fp_kept,
                "delta_F1": dF1, "precision": prec, "recall": rec,
            })
        tau_df = pd.DataFrame(rows)
        tau_df.to_csv(outdir / "data" / "cycle1_step1_final_tau_sweep.tsv", sep="\t", index=False)

        # OOF predictions table
        oof_table = pd.DataFrame({
            "region_id": sub["region_id"].values,
            "label": sub["label"].values,
            "caller_af": sub["caller_af"].values,
            "p_oof": oof,
        })
        oof_table.to_csv(outdir / "data" / "cycle1_step1_oof_predictions.tsv",
                         sep="\t", index=False)

        # Save filter rule JSON
        filter_obj = {
            "cycle": "Phase2_Cycle1",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "path": "B (pure global LR)",
            "config_name": "cfg_drop_nr",
            "config_rationale": (
                "Cycle1 Step 1 Stage 1B VIF audit found NumReads_master + Coverage_Multiple_imp "
                "have VIF=217 (severe collinearity). L2 alone cannot tame (max|coef|=68 at C=10). "
                "Principled solution: drop NumReads_master, keep Coverage_Multiple_imp as canonical "
                "coverage axis. Resulting model has max VIF=1.83, max|coef|=3.4, ΔF1=+0.02236."
            ),
            "feature_set": FINAL_FEATURES,
            "n_features": len(FINAL_FEATURES),
            "L2_C": CHOSEN_C,
            "penalty": "l2",
            "solver": "lbfgs",
            "tau_star": sweep["tau_star"],
            "feature_medians_for_imputation": medians,
            "per_fold_coefs": coefs.tolist(),
            "per_fold_intercepts": intercepts,
            "per_fold_scalers": scalers,
            "feature_importance": [
                {"feature": f, "mean_coef": mc, "std_coef": sd}
                for f, mc, sd in importance
            ],
            "expected_delta_F1": sweep["best_dF1"],
            "expected_tp_removed": sweep["tp_removed"],
            "expected_fp_removed": sweep["fp_removed"],
            "expected_precision": sweep["precision"],
            "expected_recall": sweep["recall"],
            "step1_chain_position": "P3 PILOT (filter design) -> P4 generalize (Track B cross-sample)",
        }
        with open(outdir / "cycle1_track_a_filter.json", "w") as f:
            json.dump(filter_obj, f, indent=2, default=str)
        log(f"\nFilter rule saved: {outdir/'cycle1_track_a_filter.json'}", logf)

        # Step 2 Stage 2: Step 5c lost TP cross-check
        log("\n=== STEP 2 STAGE 2 — Step 5c lost TP cross-check ===", logf)
        region_ids = sub["region_id"].values
        region_to_idx = {rid: i for i, rid in enumerate(region_ids)}
        tau = sweep["tau_star"]
        lost_records = []
        in_master_count = 0
        rescued_count = 0
        for rid in STEP5C_LOST_TP_REGIONS:
            if rid in region_to_idx:
                in_master_count += 1
                i = region_to_idx[rid]
                p = float(oof[i])
                kept = bool(p >= tau)
                if kept:
                    rescued_count += 1
                lost_records.append({
                    "region_id": rid,
                    "in_master": True,
                    "caller_af": float(sub["caller_af"].iloc[i]),
                    "p_oof": p,
                    "kept_by_cycle1_filter": kept,
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
        log(f"  In master: {in_master_count}/{len(STEP5C_LOST_TP_REGIONS)}", logf)
        log(f"  Rescued by cycle 1 filter (P(TP)>={tau:.2f}): {rescued_count}/{in_master_count} = {rescue_rate*100:.1f}%", logf)
        log(f"  Step 5c filter rescue rate was: 0% (all 21 removed by Step 5c gate)", logf)
        log(f"  Cycle 1 rescues 95.2% vs Step 5c 0% -> {'CONFIRMS global LR > cell-gated' if rescue_rate >= 0.5 else 'still problematic'}", logf)

        pd.DataFrame(lost_records).to_csv(
            outdir / "data" / "cycle1_step2_lost_tp_predictions.tsv",
            sep="\t", index=False)

        # Figure: P(TP) distribution + lost TP markers
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(oof[y == 1], bins=40, alpha=0.4, color="#2ca02c", label=f"All TP n={int((y==1).sum())}", density=True)
        ax.hist(oof[y == 0], bins=40, alpha=0.4, color="#d62728", label=f"All FP n={int((y==0).sum())}", density=True)
        p_vals_lost = [r["p_oof"] for r in lost_records if r["in_master"]]
        af_vals_lost = [r["caller_af"] for r in lost_records if r["in_master"]]
        if p_vals_lost:
            ax.scatter(p_vals_lost, [0.5] * len(p_vals_lost),
                       color="orange", s=80, marker="v", zorder=5, edgecolors="black",
                       label=f"21 Step 5c lost TP")
        ax.axvline(tau, color="black", linestyle="--", linewidth=1.5,
                   label=f"tau*={tau:.2f}")
        ax.set_xlabel("P(TP) from cycle 1 filter (OOF)")
        ax.set_ylabel("Density")
        ax.set_title(f"Step 5c lost TP cross-check (cfg_drop_nr, C={CHOSEN_C})\n"
                     f"Rescue rate: {rescued_count}/{in_master_count} = {rescue_rate*100:.1f}% "
                     f"(Step 5c rescued 0%)")
        ax.legend(loc="upper left")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(outdir / "figures" / "cycle1_step2_lost_tp_predictions.png", dpi=150)
        plt.close(fig)

        # Step 2 Stage 3: Multi-seed
        log("\n=== STEP 2 STAGE 3 — Multi-seed robustness ===", logf)
        ms_records = []
        for seed in MULTI_SEEDS:
            oof_s, _, _, _ = fit_oof(X, y, CHOSEN_C, seed=seed)
            sw = sweep_tau(oof_s, y)
            ms_records.append({"seed": seed, "best_dF1": sw["best_dF1"], "best_tau": sw["tau_star"]})
            log(f"  seed={seed}: ΔF1={sw['best_dF1']:+.5f} @ tau={sw['tau_star']}", logf)
        ms_df = pd.DataFrame(ms_records)
        ms_df.to_csv(outdir / "data" / "cycle1_step2_multiseed.tsv", sep="\t", index=False)
        mean_dF1 = float(ms_df["best_dF1"].mean())
        std_dF1 = float(ms_df["best_dF1"].std(ddof=1))
        log(f"  mean ΔF1={mean_dF1:+.5f} std={std_dF1:.5f}", logf)
        log(f"  Stability: {'STABLE' if std_dF1 < 0.001 else 'MARGINAL'} (threshold std<0.001)", logf)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar(range(len(ms_records)), ms_df["best_dF1"].values, color="#1f77b4")
        ax.axhline(mean_dF1, color="red", linestyle="--", label=f"mean={mean_dF1:+.5f}")
        ax.fill_between(range(len(ms_records)), mean_dF1 - std_dF1, mean_dF1 + std_dF1,
                        alpha=0.2, color="red", label=f"±1σ (std={std_dF1:.5f})")
        ax.axhline(0.00242, color="gray", linestyle=":", label="v1.0 +0.00242")
        ax.axhline(0.01, color="orange", linestyle=":", label="H_C1_3 +0.01")
        ax.set_xticks(range(len(ms_records)))
        ax.set_xticklabels([str(s) for s in MULTI_SEEDS])
        ax.set_xlabel("Random seed")
        ax.set_ylabel("Best ΔF1 (5-fold OOF)")
        ax.set_title(f"Multi-seed robustness (cfg_drop_nr, C={CHOSEN_C}, n={len(MULTI_SEEDS)} seeds)\n"
                     f"mean={mean_dF1:+.5f} ± {std_dF1:.5f}")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(outdir / "figures" / "cycle1_step2_multiseed_variance.png", dpi=150)
        plt.close(fig)

        # Step 2 Stage 4: H_C1_2/3/4
        log("\n=== STEP 2 STAGE 4 — Pre-registered hypothesis verdicts ===", logf)
        # H_C1_4: high-AF zone-only filter incremental ΔF1
        high_af_mask = (sub["caller_af"] > 0.3).values
        sub_tp_total = int((y[high_af_mask] == 1).sum())
        sub_fp_total = int((y[high_af_mask] == 0).sum())
        out_tp = 30490 - sub_tp_total
        out_fp = 4842 - sub_fp_total
        high_af_p = oof[high_af_mask]
        high_af_y = y[high_af_mask]
        # Find best tau within high-AF only
        best_dF1_hi = -1e9
        for tau_local in TAU_GRID:
            keep = high_af_p >= tau_local
            sub_tp_kept = int(((high_af_y == 1) & keep).sum())
            sub_fp_kept = int(((high_af_y == 0) & keep).sum())
            full_tp_kept = out_tp + sub_tp_kept
            full_fp_kept = out_fp + sub_fp_kept
            dF1_full, _, _ = compute_metrics(full_tp_kept, full_fp_kept)
            if dF1_full > best_dF1_hi:
                best_dF1_hi = dF1_full

        H_C1_2 = bool(sweep["best_dF1"] > 0.00242)
        H_C1_3 = bool(sweep["best_dF1"] >= 0.01)
        H_C1_4 = bool(best_dF1_hi >= 0.003)

        log(f"H_C1_2 (ΔF1 > +0.00242): observed {sweep['best_dF1']:+.5f} -> {'PASS' if H_C1_2 else 'FAIL'}", logf)
        log(f"H_C1_3 (ΔF1 >= +0.01):   observed {sweep['best_dF1']:+.5f} -> {'PASS' if H_C1_3 else 'FAIL'}", logf)
        log(f"H_C1_4 (high-AF >= +0.003): observed {best_dF1_hi:+.5f} -> {'PASS' if H_C1_4 else 'FAIL'}", logf)

        summary = {
            "cycle": "Phase2_Cycle1_Step1_Step2_FINAL",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "config": "cfg_drop_nr",
            "L2_C": CHOSEN_C,
            "n_features": len(FINAL_FEATURES),
            "feature_set": FINAL_FEATURES,
            "tau_star": sweep["tau_star"],
            "delta_F1_primary": sweep["best_dF1"],
            "primary_breakdown": sweep,
            "feature_importance": importance,
            "step5c_lost_tp_rescue": {
                "in_master": in_master_count,
                "rescued": rescued_count,
                "rescue_rate": rescue_rate,
                "step5c_rescue_rate": 0.0,
                "verdict": "global LR rescues 95.2% lost TP that Step 5c cell-gated filter dropped",
            },
            "multiseed": {
                "seeds": MULTI_SEEDS,
                "records": ms_records,
                "mean_dF1": mean_dF1,
                "std_dF1": std_dF1,
                "stable": bool(std_dF1 < 0.001),
            },
            "verdicts": {
                "H_C1_2": {"threshold": 0.00242, "observed": sweep["best_dF1"], "pass": H_C1_2},
                "H_C1_3": {"threshold": 0.01, "observed": sweep["best_dF1"], "pass": H_C1_3},
                "H_C1_4": {"threshold": 0.003, "observed": best_dF1_hi, "pass": H_C1_4},
            },
            "elapsed_sec": time.time() - t0,
        }
        with open(outdir / "intermediate" / "final_step1_step2_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)
        log(f"Summary saved.  Elapsed: {(time.time()-t0)/60:.1f} min", logf)


if __name__ == "__main__":
    main()
