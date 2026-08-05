#!/usr/bin/env python3
"""TO ClairS-TO Feature Discrimination Deep Study — Q6: Feature Combinations.

Cross-validated multi-feature models comparing baseline vs PON-only v2b.
GroupShuffleSplit by chromosome, 5 splits, 80/20.

Tier 1 (safe): PairwiseMeanDist, CramersV, NumReads, NumCpGs, HeuristicScore, Coverage_Multiple, NHP0, NHP3
Tier 2 (within-stratum safe): HPFineP, HPFineNGroups, AlleleDelta, AlleleP, HPMergedDelta
Tier 3 (stratum-defining): HP_Ratio, HP1FamilyN, HP2FamilyN
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
    compute_auc, ensure_dir,
    COLOR_TP, COLOR_FP,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import select_loh_legacy

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
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

OUT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/output/clairsto_correction_analysis/feature_study")
Q6_DIR = OUT_ROOT / "Q6"

BOOLEAN_FEATURES = ["PassedGating", "HPMergedSig", "HPFineSig", "AlleleSig",
                     "ClusterPermanovaValid", "Potential_LOH"]

# Feature tiers
TIER1 = ["PairwiseMeanDist", "CramersV", "NumReads", "NumCpGs",
         "HeuristicScore", "Coverage_Multiple", "NHP0", "NHP3"]
TIER2 = ["HPFineP", "HPFineNGroups", "AlleleDelta", "AlleleP", "HPMergedDelta"]
TIER3 = ["HP_Ratio", "HP1FamilyN", "HP2FamilyN"]

# Named feature sets
STAT_TRIO = ["HPFineP", "AlleleDelta", "HPFineNGroups"]
STRUCT_TRIO = ["HP_Ratio", "HP2FamilyN", "NumReads"]

STRATA = ["All", "LOH", "Non-LOH", "LOH_Strong", "LOH_Weak"]

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


def load_version(
    version: str,
    allow_unversioned_v1: bool = False,
) -> pd.DataFrame:
    paths = DATA_PATHS[version]
    tp = pd.read_csv(paths["tp"], low_memory=False)
    fp = pd.read_csv(paths["fp"], low_memory=False)
    tp["truth_label"] = "TP"
    fp["truth_label"] = "FP"
    tp["version"] = version
    fp["version"] = version
    df = attach_loh_legacy_view(
        pd.concat([tp, fp], ignore_index=True),
        allow_unversioned_v1=allow_unversioned_v1,
    )
    for col in BOOLEAN_FEATURES + ["Potential_LOH"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.lower().isin(["true", "1", "yes"]).astype(int)
    for col in TIER1 + TIER2 + TIER3:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["y_true"] = (df["truth_label"] == "TP").astype(int)
    # Extract chromosome number for GroupShuffleSplit
    df["chr_group"] = df["Chr"].str.replace("chr", "").apply(
        lambda x: int(x) if x.isdigit() else 23 if x == "X" else 24 if x == "Y" else 25)
    return df


def stratify(df, stratum):
    if stratum == "All": return df
    elif stratum == "LOH": return df[df["Potential_LOH"] == 1]
    elif stratum == "Non-LOH": return df[df["Potential_LOH"] == 0]
    elif stratum == "LOH_Strong": return df[df["_loh_subtype_legacy"] == "LOH_Strong"]
    elif stratum == "LOH_Weak": return df[df["_loh_subtype_legacy"] == "LOH_Weak"]
    return df


def cv_evaluate(df, features, stratum, model_class, model_params, n_splits=5, seed=42):
    """Cross-validated AUC with GroupShuffleSplit by chromosome.

    Returns dict: mean_auc, std_auc, fold_aucs, feature_importances.
    """
    from sklearn.model_selection import GroupShuffleSplit
    from sklearn.preprocessing import StandardScaler

    sdf = stratify(df, stratum).copy()
    valid_feats = [f for f in features if f in sdf.columns]
    X = sdf[valid_feats].values
    y = sdf["y_true"].values
    groups = sdf["chr_group"].values

    # Drop rows with NaN
    mask = np.all(np.isfinite(X), axis=1)
    X, y, groups = X[mask], y[mask], groups[mask]

    if len(np.unique(y)) < 2 or len(X) < 50:
        return {"mean_auc": float("nan"), "std_auc": float("nan"), "fold_aucs": [],
                "n": len(X), "features": valid_feats}

    gss = GroupShuffleSplit(n_splits=n_splits, test_size=0.2, random_state=seed)
    fold_aucs = []
    importances_list = []

    for train_idx, test_idx in gss.split(X, y, groups):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
            continue

        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)

        model = model_class(**model_params)
        model.fit(X_train_s, y_train)

        if hasattr(model, "predict_proba"):
            y_pred = model.predict_proba(X_test_s)[:, 1]
        else:
            y_pred = model.decision_function(X_test_s)

        auc_v = compute_auc(y_test, y_pred)
        if not np.isnan(auc_v):
            fold_aucs.append(float(auc_v))

        # Feature importance
        if hasattr(model, "coef_"):
            importances_list.append(model.coef_[0])
        elif hasattr(model, "feature_importances_"):
            importances_list.append(model.feature_importances_)

    result = {
        "mean_auc": float(np.mean(fold_aucs)) if fold_aucs else float("nan"),
        "std_auc": float(np.std(fold_aucs)) if fold_aucs else float("nan"),
        "fold_aucs": fold_aucs,
        "n": len(X),
        "features": valid_feats,
    }

    if importances_list:
        avg_imp = np.mean(importances_list, axis=0)
        result["feature_importance"] = {f: float(v) for f, v in zip(valid_feats, avg_imp)}

    return result


def run_q6(dfs: dict[str, pd.DataFrame]):
    _p("\n" + "="*60)
    _p("Q6: Feature Combinations — Cross-Validated")
    _p("="*60)

    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.tree import DecisionTreeClassifier

    results = {}

    for ver, df in dfs.items():
        _p(f"\n{'='*40}")
        _p(f"Version: {ver}")
        _p(f"{'='*40}")
        res = {
            "version": ver,
            "loh_schema_contract": dict(df.attrs["loh_schema_contract"]),
        }

        # 1. Feature correlation matrix (for fig)
        all_feats = TIER1 + TIER2 + TIER3
        valid_feats = [f for f in all_feats if f in df.columns]
        corr = df[valid_feats].corr()
        corr.to_csv(Q6_DIR / f"q6_correlation_{ver}.tsv", sep="\t")

        # 2. Model comparison across strata
        models = {
            "LR_Tier12": (LogisticRegression, {"C": 1.0, "max_iter": 500, "solver": "lbfgs"}, TIER1 + TIER2),
            "RF_Tier123": (RandomForestClassifier, {"n_estimators": 100, "max_depth": 5, "random_state": 42, "n_jobs": -1}, TIER1 + TIER2 + TIER3),
            "LR_StatTrio": (LogisticRegression, {"C": 1.0, "max_iter": 500, "solver": "lbfgs"}, STAT_TRIO),
            "LR_StructTrio": (LogisticRegression, {"C": 1.0, "max_iter": 500, "solver": "lbfgs"}, STRUCT_TRIO),
        }

        model_results = {}
        for stratum in STRATA:
            _p(f"\n  Stratum: {stratum}")
            sdf = stratify(df, stratum)
            _p(f"    n={len(sdf):,}, TP={sum(sdf.y_true==1):,}, FP={sum(sdf.y_true==0):,}")

            stratum_res = {}

            # Single-feature baselines (top 5)
            single_aucs = {}
            for feat in valid_feats:
                auc_v = compute_auc(sdf.y_true.values, sdf[feat].values)
                auc_neg = compute_auc(sdf.y_true.values, -sdf[feat].values)
                best = max(auc_v, auc_neg) if not np.isnan(auc_neg) else auc_v
                single_aucs[feat] = float(best)
            top5 = sorted(single_aucs.items(), key=lambda x: -x[1])[:5]
            _p(f"    Top 5 single features:")
            for f, a in top5:
                _p(f"      {f}: AUC={a:.4f}")
            stratum_res["single_top5"] = [{"feature": f, "auc": a} for f, a in top5]

            # Top-5 feature subset LR
            top5_feats = [f for f, _ in top5]
            cv_res = cv_evaluate(df, top5_feats, stratum, LogisticRegression,
                                 {"C": 1.0, "max_iter": 500, "solver": "lbfgs"})
            _p(f"    LR_Top5: AUC={cv_res['mean_auc']:.4f}±{cv_res['std_auc']:.4f}")
            stratum_res["LR_Top5"] = cv_res

            # Named models
            for model_name, (model_cls, params, feats) in models.items():
                cv_res = cv_evaluate(df, feats, stratum, model_cls, params)
                _p(f"    {model_name}: AUC={cv_res['mean_auc']:.4f}±{cv_res['std_auc']:.4f}")
                stratum_res[model_name] = cv_res

            model_results[stratum] = stratum_res

        res["model_results"] = model_results

        # 3. LOH_Strong Decision Tree (interpretable rules)
        _p(f"\n  LOH_Strong Decision Tree:")
        loh_s = stratify(df, "LOH_Strong")
        if len(loh_s) > 100:
            dt_feats = TIER1 + TIER2 + TIER3
            dt_valid = [f for f in dt_feats if f in loh_s.columns]
            X = loh_s[dt_valid].values
            y = loh_s["y_true"].values
            mask = np.all(np.isfinite(X), axis=1)
            X, y = X[mask], y[mask]

            dt = DecisionTreeClassifier(max_depth=3, random_state=42)
            dt.fit(X, y)
            # Extract rules
            from sklearn.tree import export_text
            rules = export_text(dt, feature_names=dt_valid)
            _p(f"    Decision Tree rules:\n{rules}")
            res["loh_strong_dt_rules"] = rules
            res["loh_strong_dt_importance"] = {f: float(v) for f, v in zip(dt_valid, dt.feature_importances_) if v > 0}

            # DT AUC
            y_pred_dt = dt.predict_proba(X)[:, 1]
            dt_auc = compute_auc(y, y_pred_dt)
            _p(f"    DT (train) AUC: {dt_auc:.4f}")
            res["loh_strong_dt_auc"] = float(dt_auc)

        results[ver] = res

    # --- Q6 Figures ---
    setup_plot_style()

    # Q6-Fig1: Feature correlation matrix heatmap
    for ver, df in dfs.items():
        fig, ax = plt.subplots(figsize=(12, 10))
        corr = df[TIER1 + TIER2 + TIER3].corr()
        mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
        sns.heatmap(corr, mask=mask, cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                    annot=True, fmt=".2f", ax=ax, square=True,
                    annot_kws={"size": 6})
        ax.set_title(f"Q6: Feature Correlation Matrix — {ver}", fontsize=12)
        fig.tight_layout()
        save_figure(fig, Q6_DIR / f"Q6_fig1_correlation_{ver}.png")

    # Q6-Fig2: Model AUC comparison across strata (2×2 layout)
    fig, axes = plt.subplots(2, len(dfs), figsize=(7 * len(dfs), 10))
    if len(dfs) == 1:
        axes = axes.reshape(-1, 1)
    model_names = ["LR_Tier12", "RF_Tier123", "LR_StatTrio", "LR_StructTrio", "LR_Top5"]
    colors_model = ["#1976D2", "#388E3C", "#D32F2F", "#7B1FA2", "#F57C00"]

    for col_idx, ver in enumerate(dfs.keys()):
        # All strata comparison
        ax = axes[0, col_idx]
        x = np.arange(len(STRATA))
        width = 0.15
        for m_idx, model_name in enumerate(model_names):
            aucs = []
            for stratum in STRATA:
                mr = results[ver]["model_results"].get(stratum, {})
                r = mr.get(model_name, {})
                aucs.append(r.get("mean_auc", float("nan")))
            ax.bar(x + m_idx * width, aucs, width, label=model_name, color=colors_model[m_idx], alpha=0.8)
        ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5)
        ax.set_xticks(x + width * 2)
        ax.set_xticklabels(STRATA, rotation=30, fontsize=8)
        ax.set_ylabel("CV AUC")
        ax.set_title(f"{ver} — Model Comparison")
        ax.legend(fontsize=6, loc="upper right")
        ax.set_ylim(0.4, 0.8)

        # Single vs combo for LOH strata
        ax = axes[1, col_idx]
        loh_strata = ["LOH", "LOH_Strong", "LOH_Weak"]
        x = np.arange(len(loh_strata))
        for s_idx, stratum in enumerate(loh_strata):
            mr = results[ver]["model_results"].get(stratum, {})
            # Best single
            top5 = mr.get("single_top5", [])
            best_single = top5[0]["auc"] if top5 else float("nan")
            # Best combo
            combo_aucs = [mr.get(m, {}).get("mean_auc", float("nan")) for m in model_names]
            best_combo = max([a for a in combo_aucs if not np.isnan(a)], default=float("nan"))
            ax.bar(s_idx - 0.15, best_single, 0.3, color="#FF9800", alpha=0.8, label="Best Single" if s_idx == 0 else "")
            ax.bar(s_idx + 0.15, best_combo, 0.3, color="#2196F3", alpha=0.8, label="Best Combo" if s_idx == 0 else "")
            # Delta annotation
            delta = best_combo - best_single
            ax.text(s_idx, max(best_single, best_combo) + 0.01, f"Δ={delta:+.3f}",
                    ha="center", fontsize=8, color="black")
        ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(loh_strata, fontsize=9)
        ax.set_ylabel("AUC")
        ax.set_title(f"{ver} — Single vs Combo (LOH)")
        ax.legend(fontsize=8)
        ax.set_ylim(0.4, 0.8)

    fig.suptitle("Q6: Cross-Validated Feature Combination Performance", fontsize=14)
    fig.tight_layout()
    save_figure(fig, Q6_DIR / "Q6_fig2_model_comparison.png")

    # Q6-Fig3: LOH ROC 3-line comparison (best single, LR Tier12, RF Tier123)
    from sklearn.model_selection import GroupShuffleSplit
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import roc_curve

    for ver, df in dfs.items():
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        for s_idx, stratum in enumerate(["LOH", "LOH_Strong", "LOH_Weak"]):
            ax = axes[s_idx]
            sdf = stratify(df, stratum).copy()
            if len(sdf) < 100:
                continue

            # Best single feature ROC
            top5 = results[ver]["model_results"].get(stratum, {}).get("single_top5", [])
            if top5:
                best_feat = top5[0]["feature"]
                y = sdf.y_true.values
                s = sdf[best_feat].values
                mask = np.isfinite(s)
                fpr, tpr, _ = roc_curve(y[mask], s[mask])
                auc_v = compute_auc(y[mask], s[mask])
                ax.plot(fpr, tpr, color="#FF9800", linestyle="--", label=f"Best Single ({best_feat}): {auc_v:.3f}")

            # LR Tier12
            feats_lr = [f for f in TIER1 + TIER2 if f in sdf.columns]
            X = sdf[feats_lr].values
            y = sdf.y_true.values
            mask = np.all(np.isfinite(X), axis=1)
            X, y = X[mask], y[mask]
            if len(np.unique(y)) >= 2:
                scaler = StandardScaler()
                X_s = scaler.fit_transform(X)
                lr = LogisticRegression(C=1.0, max_iter=500, solver="lbfgs")
                lr.fit(X_s, y)
                y_pred = lr.predict_proba(X_s)[:, 1]
                fpr, tpr, _ = roc_curve(y, y_pred)
                auc_v = compute_auc(y, y_pred)
                ax.plot(fpr, tpr, color="#1976D2", label=f"LR Tier1+2: {auc_v:.3f}")

            # RF Tier123
            feats_rf = [f for f in TIER1 + TIER2 + TIER3 if f in sdf.columns]
            X = sdf[feats_rf].values
            y = sdf.y_true.values
            mask = np.all(np.isfinite(X), axis=1)
            X, y = X[mask], y[mask]
            if len(np.unique(y)) >= 2:
                rf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42, n_jobs=-1)
                rf.fit(X, y)
                y_pred = rf.predict_proba(X)[:, 1]
                fpr, tpr, _ = roc_curve(y, y_pred)
                auc_v = compute_auc(y, y_pred)
                ax.plot(fpr, tpr, color="#388E3C", label=f"RF Tier1+2+3: {auc_v:.3f}")

            ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
            ax.set_xlabel("FPR")
            ax.set_ylabel("TPR")
            ax.set_title(f"{stratum}", fontsize=11)
            ax.legend(fontsize=7)

        fig.suptitle(f"Q6: LOH ROC — {ver} (Note: train set ROC, CV AUC in text)", fontsize=13)
        fig.tight_layout()
        save_figure(fig, Q6_DIR / f"Q6_fig3_loh_roc_{ver}.png")

    # Q6-Fig4: LOH_Strong Decision Tree visualization
    from sklearn.tree import DecisionTreeClassifier, plot_tree
    for ver, df in dfs.items():
        loh_s = stratify(df, "LOH_Strong")
        if len(loh_s) < 100:
            continue
        dt_valid = [f for f in TIER1 + TIER2 + TIER3 if f in loh_s.columns]
        X = loh_s[dt_valid].values
        y = loh_s.y_true.values
        mask = np.all(np.isfinite(X), axis=1)
        X, y = X[mask], y[mask]
        dt = DecisionTreeClassifier(max_depth=3, random_state=42)
        dt.fit(X, y)
        fig, ax = plt.subplots(figsize=(20, 8))
        plot_tree(dt, feature_names=dt_valid, class_names=["FP", "TP"],
                  filled=True, rounded=True, ax=ax, fontsize=8)
        ax.set_title(f"Q6: LOH_Strong Decision Tree (depth=3) — {ver}", fontsize=12)
        save_figure(fig, Q6_DIR / f"Q6_fig4_dt_loh_strong_{ver}.png")

    # Q6-Fig5: Feature importance (RF)
    fig, axes = plt.subplots(1, len(dfs), figsize=(7 * len(dfs), 6))
    if len(dfs) == 1:
        axes = [axes]
    for idx, ver in enumerate(dfs.keys()):
        ax = axes[idx]
        mr = results[ver]["model_results"].get("LOH", {})
        rf_res = mr.get("RF_Tier123", {})
        imp = rf_res.get("feature_importance", {})
        if imp:
            imp_df = pd.DataFrame({"feature": list(imp.keys()), "importance": list(imp.values())})
            imp_df = imp_df.sort_values("importance", ascending=True).tail(15)
            ax.barh(imp_df["feature"], imp_df["importance"], color="#388E3C", alpha=0.8)
            ax.set_xlabel("Feature Importance")
            ax.set_title(f"{ver} — RF Feature Importance (LOH)")
    fig.suptitle("Q6: Random Forest Feature Importance", fontsize=13)
    fig.tight_layout()
    save_figure(fig, Q6_DIR / "Q6_fig5_feature_importance.png")

    # Q6-Fig6: AUC ceiling summary bar chart
    fig, ax = plt.subplots(figsize=(14, 6))
    bar_data = []
    for ver in dfs.keys():
        for stratum in STRATA:
            mr = results[ver]["model_results"].get(stratum, {})
            top5 = mr.get("single_top5", [])
            best_single = top5[0]["auc"] if top5 else float("nan")
            combo_aucs = {m: mr.get(m, {}).get("mean_auc", float("nan")) for m in model_names}
            best_combo_name = max(combo_aucs, key=lambda k: combo_aucs[k] if not np.isnan(combo_aucs[k]) else -1)
            best_combo = combo_aucs[best_combo_name]
            bar_data.append({
                "version": ver, "stratum": stratum,
                "best_single": float(best_single), "best_combo": float(best_combo),
                "best_combo_model": best_combo_name,
                "delta": float(best_combo - best_single),
            })

    bdf = pd.DataFrame(bar_data)
    bdf.to_csv(Q6_DIR / "q6_auc_ceiling_summary.tsv", sep="\t", index=False)

    x = np.arange(len(STRATA))
    width = 0.2
    for v_idx, ver in enumerate(dfs.keys()):
        vdf = bdf[bdf.version == ver]
        offset = (v_idx - 0.5) * width * 2
        ax.bar(x + offset - width/2, vdf["best_single"].values, width,
               label=f"{ver} Single", color="#FF9800" if v_idx == 0 else "#FFB74D", alpha=0.8)
        ax.bar(x + offset + width/2, vdf["best_combo"].values, width,
               label=f"{ver} Combo", color="#2196F3" if v_idx == 0 else "#64B5F6", alpha=0.8)
    ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(STRATA, fontsize=10)
    ax.set_ylabel("AUC")
    ax.set_title("Q6: AUC Ceiling — Best Single Feature vs Best Combination")
    ax.legend(fontsize=7, ncol=2)
    ax.set_ylim(0.4, 0.85)
    fig.tight_layout()
    save_figure(fig, Q6_DIR / "Q6_fig6_auc_ceiling.png")

    # Save results
    with open(Q6_DIR / "q6_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    _p("\n" + "="*60)
    _p("Q6 COMPLETE")
    _p("="*60)
    return results


def main() -> None:
    args = parse_args()
    Q6_DIR.mkdir(parents=True, exist_ok=True)
    # Load data
    available = {}
    for ver, paths in DATA_PATHS.items():
        if paths["tp"].exists() and paths["fp"].exists():
            available[ver] = True
            _p(f"  {ver}: available")
        else:
            _p(f"  {ver}: MISSING")

    dfs = {}
    for ver in available:
        _p(f"\nLoading {ver}...")
        df = load_version(ver, allow_unversioned_v1=args.allow_unversioned_v1)
        _p(f"  {len(df):,} rows")
        dfs[ver] = df

    import seaborn as sns
    run_q6(dfs)


if __name__ == "__main__":
    main()
