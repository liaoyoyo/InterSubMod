#!/usr/bin/env python3
"""AlleleSig LOH Discrimination Study.

AlleleSig in LOH: TP rate=0.515, FP rate=0.369 (Δ=14.6pp).
Question: Can this signal be used for TP/FP classification within LOH,
either alone or as a strong feature in combination?

Analysis:
1. AlleleSig single-feature discriminability (AUC, enrichment, per-subtype)
2. AlleleSig as additive feature — marginal contribution to best models
3. Boolean combination filters (AlleleSig × HPFineSig × other)
4. Continuous counterpart (AlleleP, AlleleDelta) comparison
5. Practical F1 impact with TP loss constraint
6. Per-chromosome stability
7. Baseline vs v2b comparison
"""
from __future__ import annotations

import argparse
import sys
import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    setup_plot_style, save_figure, add_caption,
    compute_auc, bootstrap_ci, compute_enrichment,
    mannwhitney_test, chi2_test, ensure_dir,
    COLOR_TP, COLOR_FP, format_p, format_ci,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import select_loh_legacy

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

BASE_DIR = Path("/big7_disk/liaoyoyo2001/longphase-to-mod/output")
DATA_PATHS = {
    "pononly_v2b": {
        "tp": BASE_DIR / "ism_pononly_v2b_clairsto_tp" / "significance_summary.csv",
        "fp": BASE_DIR / "ism_pononly_v2b_clairsto_fp" / "significance_summary.csv",
    },
    "baseline": {
        "tp": BASE_DIR / "ism_baseline_clairsto_tp" / "significance_summary.csv",
        "fp": BASE_DIR / "ism_baseline_clairsto_fp" / "significance_summary.csv",
    },
}

OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/output/clairsto_correction_analysis/feature_study/AlleleSig_LOH")

BOOL_COLS = ["PassedGating", "HPMergedSig", "HPFineSig", "AlleleSig",
             "ClusterPermanovaValid", "Potential_LOH"]

def _p(msg): print(msg, flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help=(
            "Explicitly authorize historical unversioned LOH_Subtype input; "
            "versioned input requires LOH_Subtype_LegacyVC and an exact alias"
        ),
    )
    return parser.parse_args()


def attach_loh_legacy_view(
    frame: pd.DataFrame,
    allow_unversioned_v1: bool = False,
) -> pd.DataFrame:
    view = select_loh_legacy(frame, allow_unversioned_v1=allow_unversioned_v1)
    prepared = frame.copy()
    prepared["_loh_subtype_legacy"] = view.values
    prepared.attrs["loh_schema_contract"] = {
        "selection_field": view.field,
        "schema_status": view.schema_status,
        "allow_unversioned_v1": allow_unversioned_v1,
        "warnings": list(view.warning_messages),
    }
    return prepared


def load_version(version, allow_unversioned_v1: bool = False):
    paths = DATA_PATHS[version]
    tp = pd.read_csv(paths["tp"], low_memory=False)
    fp = pd.read_csv(paths["fp"], low_memory=False)
    tp["truth_label"] = "TP"; fp["truth_label"] = "FP"
    tp["version"] = version; fp["version"] = version
    df = attach_loh_legacy_view(
        pd.concat([tp, fp], ignore_index=True),
        allow_unversioned_v1=allow_unversioned_v1,
    )
    for col in BOOL_COLS:
        if col in df.columns:
            df[col] = df[col].astype(str).str.lower().isin(["true", "1", "yes"]).astype(int)
    for col in ["AlleleP", "AlleleDelta", "HPFineP", "HPFineNGroups", "HP_Ratio",
                "HP1FamilyN", "HP2FamilyN", "NumReads", "Coverage_Multiple",
                "Quality_Score", "HeuristicScore", "CramersV", "PairwiseMeanDist",
                "HPMergedDelta", "NHP0", "NHP3", "NumCpGs"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["y_true"] = (df["truth_label"] == "TP").astype(int)
    df["chr_group"] = df["Chr"].str.replace("chr", "").apply(
        lambda x: int(x) if x.isdigit() else 23 if x == "X" else 24)
    return df


def bootstrap_auc(y_true, y_score, n_boot=500, seed=42):
    y_true, y_score = np.asarray(y_true), np.asarray(y_score)
    mask = np.isfinite(y_score)
    y_true, y_score = y_true[mask], y_score[mask]
    if len(np.unique(y_true)) < 2:
        return float("nan"), float("nan"), float("nan")
    point = compute_auc(y_true, y_score)
    rng = np.random.RandomState(seed)
    aucs = []
    for _ in range(n_boot):
        idx = rng.choice(len(y_true), size=len(y_true), replace=True)
        a = compute_auc(y_true[idx], y_score[idx])
        if not np.isnan(a): aucs.append(a)
    if len(aucs) < 10:
        return point, float("nan"), float("nan")
    return point, float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def main():
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    _p("=" * 70)
    _p("AlleleSig LOH Discrimination Study")
    _p("=" * 70)

    setup_plot_style()
    results = {}

    for ver in ["pononly_v2b", "baseline"]:
        _p(f"\n{'=' * 50}")
        _p(f"Version: {ver}")
        _p(f"{'=' * 50}")
        df = load_version(ver, allow_unversioned_v1=args.allow_unversioned_v1)
        loh = df[df["Potential_LOH"] == 1].copy()
        res = {
            "version": ver,
            "loh_n": len(loh),
            "loh_schema_contract": dict(df.attrs["loh_schema_contract"]),
        }
        n_tp = (loh.y_true == 1).sum()
        n_fp = (loh.y_true == 0).sum()
        _p(f"  LOH: {len(loh):,} (TP={n_tp:,}, FP={n_fp:,})")

        # ======================================================================
        # 1. AlleleSig single-feature discriminability
        # ======================================================================
        _p(f"\n--- 1. AlleleSig Single-Feature ---")

        # Basic rates
        tp_rate = loh[loh.y_true == 1]["AlleleSig"].mean()
        fp_rate = loh[loh.y_true == 0]["AlleleSig"].mean()
        delta = tp_rate - fp_rate
        _p(f"  AlleleSig rate: TP={tp_rate:.4f}, FP={fp_rate:.4f}, Δ={delta:+.4f}")

        # AUC (AlleleSig=1 predicts TP=1)
        auc_pt, auc_lo, auc_hi = bootstrap_auc(loh.y_true.values, loh["AlleleSig"].values)
        _p(f"  AlleleSig AUC (LOH): {auc_pt:.4f} {format_ci(auc_lo, auc_hi)}")
        res["allelesig_auc_loh"] = {"auc": auc_pt, "ci_lo": auc_lo, "ci_hi": auc_hi}

        # Per-subtype
        for sub in ["LOH_Strong", "LOH_Weak", "LOH_Subclone", "LOH_Noise"]:
            sdf = loh[loh["_loh_subtype_legacy"] == sub]
            if len(sdf) < 50: continue
            tp_r = sdf[sdf.y_true == 1]["AlleleSig"].mean()
            fp_r = sdf[sdf.y_true == 0]["AlleleSig"].mean()
            auc_v = compute_auc(sdf.y_true.values, sdf["AlleleSig"].values)
            _p(f"  {sub}: TP={tp_r:.4f}, FP={fp_r:.4f}, Δ={tp_r-fp_r:+.4f}, AUC={auc_v:.4f}, n={len(sdf):,}")
            res[f"allelesig_{sub}"] = {"tp_rate": float(tp_r), "fp_rate": float(fp_r),
                                        "delta": float(tp_r - fp_r), "auc": float(auc_v), "n": len(sdf)}

        # As filter: AlleleSig=false → flag as potential FP
        asig_false = loh["AlleleSig"] == 0
        tp_flagged = ((loh.y_true == 1) & asig_false).sum()
        fp_flagged = ((loh.y_true == 0) & asig_false).sum()
        tp_loss_pct = tp_flagged / n_tp * 100
        fp_rem_pct = fp_flagged / n_fp * 100
        _p(f"  Filter AlleleSig=false (LOH): flagged={asig_false.sum():,}, TP_loss={tp_loss_pct:.1f}%, FP_removal={fp_rem_pct:.1f}%")
        res["filter_allelesig_false"] = {"tp_loss": float(tp_loss_pct), "fp_removal": float(fp_rem_pct),
                                          "tp_flagged": int(tp_flagged), "fp_flagged": int(fp_flagged)}

        # Enrichment ratio (FP rate of AlleleSig=false vs true)
        asig_true_fp_rate = (loh[(loh.AlleleSig == 1) & (loh.y_true == 0)]).shape[0] / max(loh[loh.AlleleSig == 1].shape[0], 1)
        asig_false_fp_rate = (loh[(loh.AlleleSig == 0) & (loh.y_true == 0)]).shape[0] / max(loh[loh.AlleleSig == 0].shape[0], 1)
        _p(f"  FP rate: AlleleSig=true → {asig_true_fp_rate:.4f}, AlleleSig=false → {asig_false_fp_rate:.4f}")
        _p(f"  Enrichment (false/true): {asig_false_fp_rate/max(asig_true_fp_rate, 0.001):.2f}×")
        res["fp_rate_by_allelesig"] = {"true": float(asig_true_fp_rate), "false": float(asig_false_fp_rate)}

        # ======================================================================
        # 2. Continuous counterparts — AlleleP, AlleleDelta
        # ======================================================================
        _p(f"\n--- 2. Continuous Counterparts ---")
        for feat in ["AlleleP", "AlleleDelta"]:
            auc_pt, auc_lo, auc_hi = bootstrap_auc(loh.y_true.values, loh[feat].values)
            auc_neg, _, _ = bootstrap_auc(loh.y_true.values, -loh[feat].values)
            best = max(auc_pt, auc_neg) if not np.isnan(auc_neg) else auc_pt
            _p(f"  {feat} AUC (LOH): {best:.4f} (raw={auc_pt:.4f}, neg={auc_neg:.4f}) {format_ci(auc_lo, auc_hi)}")
            res[f"{feat}_auc"] = float(best)

        # ======================================================================
        # 3. Boolean combinations
        # ======================================================================
        _p(f"\n--- 3. Boolean Combination Filters (LOH) ---")

        combos = {
            "AlleleSig=F only": loh["AlleleSig"] == 0,
            "HPFineSig=F only": loh["HPFineSig"] == 0,
            "AlleleSig=F AND HPFineSig=F": (loh["AlleleSig"] == 0) & (loh["HPFineSig"] == 0),
            "AlleleSig=F OR HPFineSig=F": (loh["AlleleSig"] == 0) | (loh["HPFineSig"] == 0),
            "Both=T (keep rule)": (loh["AlleleSig"] == 1) & (loh["HPFineSig"] == 1),
            "AlleleSig=F AND HPMergedSig=F": (loh["AlleleSig"] == 0) & (loh["HPMergedSig"] == 0),
            "AlleleSig=F AND PassedGating=F": (loh["AlleleSig"] == 0) & (loh["PassedGating"] == 0),
        }

        combo_rows = []
        for name, mask in combos.items():
            tp_f = ((loh.y_true == 1) & mask).sum()
            fp_f = ((loh.y_true == 0) & mask).sum()
            total_f = mask.sum()
            tl = tp_f / n_tp * 100
            fr = fp_f / n_fp * 100
            prec = fp_f / max(total_f, 1) * 100  # FP precision of the filter
            _p(f"  {name}: flagged={total_f:,}, TP_loss={tl:.1f}%, FP_removal={fr:.1f}%, FP_prec={prec:.1f}%")
            combo_rows.append({"filter": name, "flagged": int(total_f),
                               "tp_loss_pct": float(tl), "fp_removal_pct": float(fr),
                               "fp_precision_pct": float(prec)})
        res["boolean_combos"] = combo_rows
        pd.DataFrame(combo_rows).to_csv(OUT_DIR / f"boolean_combos_{ver}.tsv", sep="\t", index=False)

        # Inverse: AlleleSig=T AND HPFineSig=T as "confident TP" gate
        both_true = (loh["AlleleSig"] == 1) & (loh["HPFineSig"] == 1)
        tp_kept = ((loh.y_true == 1) & both_true).sum()
        fp_kept = ((loh.y_true == 0) & both_true).sum()
        kept_total = both_true.sum()
        _p(f"\n  ** Confidence Gate (AlleleSig=T & HPFineSig=T) **")
        _p(f"  Kept: {kept_total:,} ({kept_total/len(loh)*100:.1f}% of LOH)")
        _p(f"  TP kept: {tp_kept:,}/{n_tp:,} ({tp_kept/n_tp*100:.1f}%)")
        _p(f"  FP kept: {fp_kept:,}/{n_fp:,} ({fp_kept/n_fp*100:.1f}%)")
        _p(f"  FP rate inside gate: {fp_kept/max(kept_total,1)*100:.1f}%")
        _p(f"  FP rate outside gate: {(n_fp-fp_kept)/max(len(loh)-kept_total,1)*100:.1f}%")
        res["confidence_gate"] = {
            "tp_kept_pct": float(tp_kept / n_tp * 100),
            "fp_kept_pct": float(fp_kept / n_fp * 100),
            "fp_rate_inside": float(fp_kept / max(kept_total, 1) * 100),
            "fp_rate_outside": float((n_fp - fp_kept) / max(len(loh) - kept_total, 1) * 100),
            "kept_total": int(kept_total),
        }

        # ======================================================================
        # 4. Marginal contribution to ML models
        # ======================================================================
        _p(f"\n--- 4. Marginal Contribution to ML Models ---")
        from sklearn.linear_model import LogisticRegression
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.model_selection import GroupShuffleSplit
        from sklearn.preprocessing import StandardScaler

        base_feats = ["HPFineP", "HPFineNGroups", "AlleleDelta", "AlleleP", "HPMergedDelta",
                      "PairwiseMeanDist", "CramersV", "NumReads", "NumCpGs",
                      "HeuristicScore", "Coverage_Multiple", "NHP0", "NHP3"]
        bool_addons = ["AlleleSig", "HPFineSig", "HPMergedSig", "PassedGating"]
        struct_feats = ["HP_Ratio", "HP1FamilyN", "HP2FamilyN"]

        feature_sets = {
            "Base (13 continuous)": base_feats,
            "Base + AlleleSig": base_feats + ["AlleleSig"],
            "Base + AlleleSig + HPFineSig": base_feats + ["AlleleSig", "HPFineSig"],
            "Base + all 4 boolean": base_feats + bool_addons,
            "Full (13+4+3 struct)": base_feats + bool_addons + struct_feats,
            "Full - AlleleSig": [f for f in base_feats + bool_addons + struct_feats if f != "AlleleSig"],
        }

        cv_results = []
        for name, feats in feature_sets.items():
            valid_feats = [f for f in feats if f in loh.columns]
            X = loh[valid_feats].values
            y = loh["y_true"].values
            groups = loh["chr_group"].values
            mask = np.all(np.isfinite(X), axis=1)
            X, y, groups = X[mask], y[mask], groups[mask]

            if len(np.unique(y)) < 2:
                continue

            gss = GroupShuffleSplit(n_splits=5, test_size=0.2, random_state=42)
            fold_aucs_lr = []
            fold_aucs_rf = []

            for train_idx, test_idx in gss.split(X, y, groups):
                X_tr, X_te = X[train_idx], X[test_idx]
                y_tr, y_te = y[train_idx], y[test_idx]
                if len(np.unique(y_tr)) < 2 or len(np.unique(y_te)) < 2:
                    continue

                scaler = StandardScaler()
                X_tr_s = scaler.fit_transform(X_tr)
                X_te_s = scaler.transform(X_te)

                # LR
                lr = LogisticRegression(C=1.0, max_iter=500, solver="lbfgs")
                lr.fit(X_tr_s, y_tr)
                pred_lr = lr.predict_proba(X_te_s)[:, 1]
                a = compute_auc(y_te, pred_lr)
                if not np.isnan(a): fold_aucs_lr.append(a)

                # RF
                rf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42, n_jobs=-1)
                rf.fit(X_tr, y_tr)
                pred_rf = rf.predict_proba(X_te)[:, 1]
                a = compute_auc(y_te, pred_rf)
                if not np.isnan(a): fold_aucs_rf.append(a)

            lr_auc = float(np.mean(fold_aucs_lr)) if fold_aucs_lr else float("nan")
            rf_auc = float(np.mean(fold_aucs_rf)) if fold_aucs_rf else float("nan")
            lr_std = float(np.std(fold_aucs_lr)) if fold_aucs_lr else float("nan")
            rf_std = float(np.std(fold_aucs_rf)) if fold_aucs_rf else float("nan")

            _p(f"  {name:40s}: LR={lr_auc:.4f}±{lr_std:.4f}, RF={rf_auc:.4f}±{rf_std:.4f}")
            cv_results.append({"feature_set": name, "n_features": len(valid_feats),
                               "lr_auc": lr_auc, "lr_std": lr_std,
                               "rf_auc": rf_auc, "rf_std": rf_std})

        res["cv_results"] = cv_results
        pd.DataFrame(cv_results).to_csv(OUT_DIR / f"cv_marginal_{ver}.tsv", sep="\t", index=False)

        # AlleleSig marginal contribution
        base_rf = [r for r in cv_results if r["feature_set"] == "Base (13 continuous)"]
        plus_asig = [r for r in cv_results if r["feature_set"] == "Base + AlleleSig"]
        if base_rf and plus_asig:
            delta_lr = plus_asig[0]["lr_auc"] - base_rf[0]["lr_auc"]
            delta_rf = plus_asig[0]["rf_auc"] - base_rf[0]["rf_auc"]
            _p(f"\n  ** AlleleSig Marginal Contribution **")
            _p(f"  LR: {base_rf[0]['lr_auc']:.4f} → {plus_asig[0]['lr_auc']:.4f} (Δ={delta_lr:+.4f})")
            _p(f"  RF: {base_rf[0]['rf_auc']:.4f} → {plus_asig[0]['rf_auc']:.4f} (Δ={delta_rf:+.4f})")
            res["marginal_lr"] = float(delta_lr)
            res["marginal_rf"] = float(delta_rf)

        # Full vs Full-AlleleSig
        full = [r for r in cv_results if r["feature_set"] == "Full (13+4+3 struct)"]
        full_minus = [r for r in cv_results if r["feature_set"] == "Full - AlleleSig"]
        if full and full_minus:
            delta_rf = full[0]["rf_auc"] - full_minus[0]["rf_auc"]
            _p(f"  Full RF with AlleleSig: {full[0]['rf_auc']:.4f}")
            _p(f"  Full RF without AlleleSig: {full_minus[0]['rf_auc']:.4f}")
            _p(f"  AlleleSig removal impact: {delta_rf:+.4f}")
            res["removal_impact_rf"] = float(delta_rf)

        # ======================================================================
        # 5. Per-chromosome stability
        # ======================================================================
        _p(f"\n--- 5. Per-Chromosome Stability ---")
        chr_rows = []
        for chrom in sorted(loh["Chr"].unique()):
            cdf = loh[loh["Chr"] == chrom]
            ctp = (cdf.y_true == 1).sum()
            cfp = (cdf.y_true == 0).sum()
            if ctp < 20 or cfp < 20:
                continue
            tp_r = cdf[cdf.y_true == 1]["AlleleSig"].mean()
            fp_r = cdf[cdf.y_true == 0]["AlleleSig"].mean()
            chr_rows.append({"chr": chrom, "n": len(cdf), "n_tp": int(ctp), "n_fp": int(cfp),
                             "tp_rate": float(tp_r), "fp_rate": float(fp_r),
                             "delta": float(tp_r - fp_r)})
        chr_df = pd.DataFrame(chr_rows)
        chr_df.to_csv(OUT_DIR / f"per_chr_allelesig_{ver}.tsv", sep="\t", index=False)
        _p(f"  {len(chr_df)} chromosomes, Δ range: [{chr_df.delta.min():.4f}, {chr_df.delta.max():.4f}]")
        _p(f"  Δ>0 (TP>FP): {(chr_df.delta > 0).sum()}/{len(chr_df)}")
        res["chr_stability"] = {"n_chr": len(chr_df),
                                "delta_min": float(chr_df.delta.min()),
                                "delta_max": float(chr_df.delta.max()),
                                "n_positive": int((chr_df.delta > 0).sum())}

        # ======================================================================
        # 6. F1 Impact Analysis (global, with TP loss constraint)
        # ======================================================================
        _p(f"\n--- 6. F1 Impact Analysis ---")
        truth_total_tp = 39447  # total somatic truth set
        # Current baseline F1 (no ISM filter)
        caller_tp = (df.y_true == 1).sum()
        caller_fp = (df.y_true == 0).sum()
        baseline_prec = caller_tp / (caller_tp + caller_fp)
        baseline_rec = caller_tp / truth_total_tp
        baseline_f1 = 2 * baseline_prec * baseline_rec / (baseline_prec + baseline_rec)
        _p(f"  Caller baseline: Prec={baseline_prec:.4f}, Rec={baseline_rec:.4f}, F1={baseline_f1:.4f}")

        # Apply AlleleSig=false LOH filter
        for filter_name, mask in [
            ("AlleleSig=F (LOH only)", (df["Potential_LOH"] == 1) & (df["AlleleSig"] == 0)),
            ("AlleleSig=F & HPFineSig=F (LOH)", (df["Potential_LOH"] == 1) & (df["AlleleSig"] == 0) & (df["HPFineSig"] == 0)),
            ("Both=T gate (keep only)", ~((df["Potential_LOH"] == 1) & (df["AlleleSig"] == 1) & (df["HPFineSig"] == 1)) & (df["Potential_LOH"] == 1)),
        ]:
            tp_removed = ((df.y_true == 1) & mask).sum()
            fp_removed = ((df.y_true == 0) & mask).sum()
            new_tp = caller_tp - tp_removed
            new_fp = caller_fp - fp_removed
            new_prec = new_tp / max(new_tp + new_fp, 1)
            new_rec = new_tp / truth_total_tp
            new_f1 = 2 * new_prec * new_rec / (new_prec + new_rec) if (new_prec + new_rec) > 0 else 0
            delta_f1 = new_f1 - baseline_f1
            tp_loss = tp_removed / caller_tp * 100
            _p(f"  {filter_name}: TP_loss={tp_loss:.1f}%, FP_removed={fp_removed:,}, ΔF1={delta_f1:+.5f}")
            res[f"f1_{filter_name}"] = {"tp_loss": float(tp_loss), "fp_removed": int(fp_removed),
                                         "delta_f1": float(delta_f1), "new_f1": float(new_f1)}

        # ======================================================================
        # 7. AlleleSig conditional on other features
        # ======================================================================
        _p(f"\n--- 7. AlleleSig Conditional Analysis ---")
        # Does AlleleSig add info beyond AlleleP?
        # Bin AlleleP into quartiles and check AlleleSig rate
        loh_copy = loh.copy()
        try:
            loh_copy["AlleleP_bin"] = pd.qcut(loh_copy["AlleleP"].clip(0, 1), q=4, duplicates="drop")
        except Exception:
            loh_copy["AlleleP_bin"] = pd.cut(loh_copy["AlleleP"].clip(0, 1), bins=4)
        for q in loh_copy["AlleleP_bin"].unique():
            qdf = loh_copy[loh_copy["AlleleP_bin"] == q]
            if len(qdf) < 100:
                continue
            tp_r = qdf[qdf.y_true == 1]["AlleleSig"].mean()
            fp_r = qdf[qdf.y_true == 0]["AlleleSig"].mean()
            _p(f"  AlleleP {q}: n={len(qdf):,}, AlleleSig TP={tp_r:.4f}, FP={fp_r:.4f}, Δ={tp_r-fp_r:+.4f}")

        # AlleleSig conditional on HPFineSig
        for hp_val in [0, 1]:
            sub = loh[loh["HPFineSig"] == hp_val]
            if len(sub) < 100:
                continue
            tp_r = sub[sub.y_true == 1]["AlleleSig"].mean()
            fp_r = sub[sub.y_true == 0]["AlleleSig"].mean()
            _p(f"  HPFineSig={hp_val}: n={len(sub):,}, AlleleSig TP={tp_r:.4f}, FP={fp_r:.4f}, Δ={tp_r-fp_r:+.4f}")

        # ======================================================================
        # 8. LOH Germline FP — Do they truly lack ASM?
        # ======================================================================
        _p(f"\n--- 8. LOH Germline FP ASM Verification ---")

        # AlleleDelta distribution: LOH TP vs FP
        tp_ad = loh[loh.y_true == 1]["AlleleDelta"]
        fp_ad = loh[loh.y_true == 0]["AlleleDelta"]
        _p(f"  AlleleDelta LOH TP: mean={tp_ad.mean():.4f}, median={tp_ad.median():.4f}, zero%={(tp_ad==0).mean()*100:.1f}%")
        _p(f"  AlleleDelta LOH FP: mean={fp_ad.mean():.4f}, median={fp_ad.median():.4f}, zero%={(fp_ad==0).mean()*100:.1f}%")

        # AlleleDelta by AlleleSig status
        _p(f"\n  AlleleDelta by AlleleSig status (LOH):")
        for asig in [0, 1]:
            for label in ["TP", "FP"]:
                sub = loh[(loh.AlleleSig == asig) & (loh.truth_label == label)]["AlleleDelta"]
                if len(sub) < 10: continue
                _p(f"    AlleleSig={asig}, {label}: n={len(sub):,}, mean={sub.mean():.4f}, median={sub.median():.4f}, zero%={(sub==0).mean()*100:.1f}%")

        # AlleleP by AlleleSig status — does the boolean threshold match continuous?
        _p(f"\n  AlleleP by AlleleSig status (LOH):")
        for asig in [0, 1]:
            for label in ["TP", "FP"]:
                sub = loh[(loh.AlleleSig == asig) & (loh.truth_label == label)]["AlleleP"]
                if len(sub) < 10: continue
                _p(f"    AlleleSig={asig}, {label}: n={len(sub):,}, mean={sub.mean():.4f}, median={sub.median():.4f}")

        # HP ratio distribution for AlleleSig=true vs false (verify LOH mechanism)
        _p(f"\n  HP_Ratio by AlleleSig (LOH):")
        for asig in [0, 1]:
            sub = loh[loh.AlleleSig == asig]["HP_Ratio"]
            _p(f"    AlleleSig={asig}: n={len(sub):,}, mean={sub.mean():.4f}, median={sub.median():.4f}")

        # Minor haplotype read count for AlleleSig T/F
        minor_hp = np.minimum(loh["HP1FamilyN"].values, loh["HP2FamilyN"].values)
        loh_copy2 = loh.copy()
        loh_copy2["minor_hp"] = minor_hp
        _p(f"\n  Minor HP count by AlleleSig (LOH):")
        for asig in [0, 1]:
            for label in ["TP", "FP"]:
                sub = loh_copy2[(loh_copy2.AlleleSig == asig) & (loh_copy2.truth_label == label)]["minor_hp"]
                if len(sub) < 10: continue
                _p(f"    AlleleSig={asig}, {label}: mean={sub.mean():.1f}, median={sub.median():.1f}, zero%={(sub==0).mean()*100:.1f}%")
        res["asm_verification"] = {
            "tp_ad_mean": float(tp_ad.mean()), "fp_ad_mean": float(fp_ad.mean()),
            "tp_ad_zero_pct": float((tp_ad == 0).mean() * 100),
            "fp_ad_zero_pct": float((fp_ad == 0).mean() * 100),
        }

        # ======================================================================
        # 9. Why baseline > v2b — LOH composition analysis
        # ======================================================================
        _p(f"\n--- 9. LOH Composition Analysis ({ver}) ---")

        # LOH subtype distribution
        for sub in ["LOH_Strong", "LOH_Weak", "LOH_Subclone", "LOH_Noise"]:
            sdf = loh[loh["_loh_subtype_legacy"] == sub]
            pct = len(sdf) / len(loh) * 100
            tp_pct = (sdf.y_true == 1).sum() / max(n_tp, 1) * 100
            fp_pct = (sdf.y_true == 0).sum() / max(n_fp, 1) * 100
            _p(f"  {sub}: {len(sdf):,} ({pct:.1f}%), TP%={tp_pct:.1f}%, FP%={fp_pct:.1f}%")

        # HP_Ratio distribution (skewness)
        _p(f"  HP_Ratio LOH: mean={loh['HP_Ratio'].mean():.4f}, median={loh['HP_Ratio'].median():.4f}")
        _p(f"  HP_Ratio LOH TP: mean={loh[loh.y_true==1]['HP_Ratio'].mean():.4f}")
        _p(f"  HP_Ratio LOH FP: mean={loh[loh.y_true==0]['HP_Ratio'].mean():.4f}")

        # Average minor HP reads
        loh_minor = np.minimum(loh["HP1FamilyN"].values, loh["HP2FamilyN"].values)
        tp_minor = loh_minor[loh.y_true.values == 1]
        fp_minor = loh_minor[loh.y_true.values == 0]
        _p(f"  Minor HP reads LOH TP: mean={tp_minor.mean():.1f}, median={np.median(tp_minor):.1f}")
        _p(f"  Minor HP reads LOH FP: mean={fp_minor.mean():.1f}, median={np.median(fp_minor):.1f}")
        res["loh_composition"] = {
            "hp_ratio_tp_mean": float(loh[loh.y_true == 1]["HP_Ratio"].mean()),
            "hp_ratio_fp_mean": float(loh[loh.y_true == 0]["HP_Ratio"].mean()),
            "minor_hp_tp_mean": float(tp_minor.mean()),
            "minor_hp_fp_mean": float(fp_minor.mean()),
        }

        results[ver] = res

    # ======================================================================
    # 10. Cross-version comparison — why baseline > v2b
    # ======================================================================
    _p(f"\n{'=' * 50}")
    _p("Cross-Version Comparison: Why Baseline > v2b in LOH")
    _p(f"{'=' * 50}")

    for metric in ["allelesig_auc_loh", "AlleleP_auc", "AlleleDelta_auc"]:
        v2b_val = results["pononly_v2b"].get(metric, {})
        bl_val = results["baseline"].get(metric, {})
        if isinstance(v2b_val, dict):
            v2b_val = v2b_val.get("auc", float("nan"))
        if isinstance(bl_val, dict):
            bl_val = bl_val.get("auc", float("nan"))
        _p(f"  {metric}: v2b={v2b_val:.4f}, baseline={bl_val:.4f}, Δ={bl_val-v2b_val:+.4f}")

    # LOH size difference
    v2b_loh = results["pononly_v2b"]["loh_n"]
    bl_loh = results["baseline"]["loh_n"]
    _p(f"\n  LOH size: v2b={v2b_loh:,}, baseline={bl_loh:,}, diff={v2b_loh-bl_loh:,}")
    _p(f"  Extra LOH in v2b: {v2b_loh - bl_loh:,} regions ({(v2b_loh-bl_loh)/bl_loh*100:.1f}% more)")

    # AlleleSig delta comparison
    v2b_delta = results["pononly_v2b"]["filter_allelesig_false"]
    bl_delta = results["baseline"]["filter_allelesig_false"]
    _p(f"\n  AlleleSig=false filter:")
    _p(f"    v2b: TP_loss={v2b_delta['tp_loss']:.1f}%, FP_rem={v2b_delta['fp_removal']:.1f}%")
    _p(f"    baseline: TP_loss={bl_delta['tp_loss']:.1f}%, FP_rem={bl_delta['fp_removal']:.1f}%")

    # FP rate comparison
    v2b_gate = results["pononly_v2b"]["confidence_gate"]
    bl_gate = results["baseline"]["confidence_gate"]
    _p(f"\n  Confidence Gate (AlleleSig=T & HPFineSig=T):")
    _p(f"    v2b: FP_inside={v2b_gate['fp_rate_inside']:.1f}%, FP_outside={v2b_gate['fp_rate_outside']:.1f}%, ratio={v2b_gate['fp_rate_outside']/max(v2b_gate['fp_rate_inside'],0.01):.2f}×")
    _p(f"    baseline: FP_inside={bl_gate['fp_rate_inside']:.1f}%, FP_outside={bl_gate['fp_rate_outside']:.1f}%, ratio={bl_gate['fp_rate_outside']/max(bl_gate['fp_rate_inside'],0.01):.2f}×")

    # Marginal contribution comparison
    _p(f"\n  AlleleSig marginal contribution (LOH RF):")
    _p(f"    v2b: {results['pononly_v2b'].get('marginal_rf', float('nan')):+.4f}")
    _p(f"    baseline: {results['baseline'].get('marginal_rf', float('nan')):+.4f}")

    # ======================================================================
    # Figures
    # ======================================================================
    _p(f"\n{'=' * 50}")
    _p("Generating Figures")
    _p(f"{'=' * 50}")

    # Fig 1: AlleleSig rate by LOH subtype (both versions)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for idx, ver in enumerate(results.keys()):
        ax = axes[idx]
        subtypes = ["LOH_Strong", "LOH_Weak", "LOH_Subclone", "LOH_Noise"]
        tp_rates = [results[ver].get(f"allelesig_{s}", {}).get("tp_rate", 0) for s in subtypes]
        fp_rates = [results[ver].get(f"allelesig_{s}", {}).get("fp_rate", 0) for s in subtypes]
        x = np.arange(len(subtypes))
        ax.bar(x - 0.15, tp_rates, 0.3, color=COLOR_TP, alpha=0.8, label="TP")
        ax.bar(x + 0.15, fp_rates, 0.3, color=COLOR_FP, alpha=0.8, label="FP")
        # Delta annotations
        for i, (t, f) in enumerate(zip(tp_rates, fp_rates)):
            ax.text(i, max(t, f) + 0.02, f"Δ={t-f:+.3f}", ha="center", fontsize=8, weight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(subtypes, fontsize=9)
        ax.set_ylabel("AlleleSig Rate")
        ax.set_title(f"{ver}", fontsize=11)
        ax.legend()
        ax.set_ylim(0, 1.0)
    fig.suptitle("AlleleSig Rate by LOH Subtype (TP vs FP)", fontsize=13)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "fig1_allelesig_by_subtype.png")

    # Fig 2: Boolean combination filter comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for idx, ver in enumerate(results.keys()):
        ax = axes[idx]
        combos = results[ver]["boolean_combos"]
        cdf = pd.DataFrame(combos)
        x = np.arange(len(cdf))
        ax.barh(x - 0.15, cdf["tp_loss_pct"], 0.3, color=COLOR_TP, alpha=0.8, label="TP Loss %")
        ax.barh(x + 0.15, cdf["fp_removal_pct"], 0.3, color=COLOR_FP, alpha=0.8, label="FP Removal %")
        ax.axvline(x=2.0, color="gray", linestyle="--", alpha=0.7, label="2% TP Loss Limit")
        ax.set_yticks(x)
        ax.set_yticklabels(cdf["filter"], fontsize=7)
        ax.set_xlabel("Percentage (%)")
        ax.set_title(f"{ver}", fontsize=11)
        ax.legend(fontsize=7)
    fig.suptitle("Boolean Combination Filter: TP Loss vs FP Removal (LOH)", fontsize=13)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "fig2_boolean_filter_comparison.png")

    # Fig 3: ML model marginal contribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for idx, ver in enumerate(results.keys()):
        ax = axes[idx]
        cv = pd.DataFrame(results[ver]["cv_results"])
        x = np.arange(len(cv))
        ax.barh(x - 0.15, cv["lr_auc"], 0.3, color="#1976D2", alpha=0.8, label="LR")
        ax.barh(x + 0.15, cv["rf_auc"], 0.3, color="#388E3C", alpha=0.8, label="RF")
        ax.axvline(x=0.5, color="gray", linestyle="--", alpha=0.5)
        ax.set_yticks(x)
        ax.set_yticklabels(cv["feature_set"], fontsize=7)
        ax.set_xlabel("CV AUC")
        ax.set_title(f"{ver} — LOH", fontsize=11)
        ax.legend(fontsize=8)
        ax.set_xlim(0.55, 0.80)
    fig.suptitle("AlleleSig Marginal Contribution to ML Models (LOH, CV AUC)", fontsize=13)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "fig3_marginal_contribution.png")

    # Fig 4: Per-chromosome stability
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for idx, ver in enumerate(results.keys()):
        ax = axes[idx]
        cdf = pd.read_csv(OUT_DIR / f"per_chr_allelesig_{ver}.tsv", sep="\t")
        cdf["chr_num"] = cdf["chr"].str.replace("chr", "").apply(
            lambda x: int(x) if x.isdigit() else 23 if x == "X" else 24)
        cdf = cdf.sort_values("chr_num")
        colors = ["#4CAF50" if d > 0.05 else "#FF9800" if d > 0 else "#F44336" for d in cdf["delta"]]
        ax.bar(range(len(cdf)), cdf["delta"], color=colors, alpha=0.8)
        ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
        ax.set_xticks(range(len(cdf)))
        ax.set_xticklabels(cdf["chr"].str.replace("chr", ""), rotation=45, fontsize=7)
        ax.set_ylabel("AlleleSig Δ (TP - FP)")
        ax.set_title(f"{ver}")
    fig.suptitle("Per-Chromosome AlleleSig TP-FP Delta (LOH)", fontsize=13)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "fig4_per_chr_stability.png")

    # Fig 5: Confidence Gate — FP rate inside vs outside
    fig, ax = plt.subplots(figsize=(10, 5))
    vers = list(results.keys())
    x = np.arange(len(vers))
    inside = [results[v]["confidence_gate"]["fp_rate_inside"] for v in vers]
    outside = [results[v]["confidence_gate"]["fp_rate_outside"] for v in vers]
    ax.bar(x - 0.15, inside, 0.3, color="#4CAF50", alpha=0.8, label="FP rate INSIDE gate (Both=T)")
    ax.bar(x + 0.15, outside, 0.3, color="#F44336", alpha=0.8, label="FP rate OUTSIDE gate")
    for i, (ins, out) in enumerate(zip(inside, outside)):
        ax.text(i, max(ins, out) + 0.5, f"ratio={out/max(ins,0.01):.2f}×", ha="center", fontsize=10, weight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(vers, fontsize=10)
    ax.set_ylabel("FP Rate (%)")
    ax.set_title("Confidence Gate: AlleleSig=T & HPFineSig=T — FP Rate Inside vs Outside (LOH)")
    ax.legend()
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "fig5_confidence_gate.png")

    # Save all results
    with open(OUT_DIR / "allelesig_study_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    _p(f"\n{'=' * 70}")
    _p("COMPLETE")
    _p(f"Output: {OUT_DIR}")
    _p(f"{'=' * 70}")


if __name__ == "__main__":
    main()
