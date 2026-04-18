#!/usr/bin/env python3
"""G6: Multi-Feature Combination — integrate G2-G5 features into composite germline FP detector.

Combines AF, GT, H flag, CpG context, and strand bias into interpretable
rule-based and statistical classifiers. Key constraint: TP loss < 2%.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    OUTPUT_ROOT, SAMPLE_ORDER, setup_plot_style, save_figure, add_caption,
    compute_auc, bootstrap_ci, ensure_dir,
    COLOR_TP, COLOR_FP,
)

WORKSPACE = OUTPUT_ROOT / "20260401_germline_fp_identification"
G1_DIR = WORKSPACE / "G1_vcf_features"
G6_DIR = ensure_dir(WORKSPACE / "G6_combination")
FIG_DIR = ensure_dir(G6_DIR / "figures")
DATA_DIR = ensure_dir(G6_DIR / "data")


def _print(msg: str):
    print(msg, flush=True)


def load_pass_features() -> pd.DataFrame:
    path = G1_DIR / "all_samples_merged.tsv.gz"
    _print(f"[G6] Loading from {path.name}")
    df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    for col in ["af", "sb_clairs", "sb_ratio", "sb_fisher_p"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["is_het", "is_hom_alt", "h_flag", "is_cpg_ct", "is_cpg_site", "is_transition",
                "pon_1", "pon_2", "pon_3", "pon_4", "filter_is_pass"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    _print(f"  -> {len(df):,} rows")
    return df


def prepare_features(df: pd.DataFrame) -> pd.DataFrame:
    """Prepare feature matrix for PASS variants only."""
    pass_df = df[df["filter_is_pass"] == 1].copy()
    af = pass_df["af"].values

    # Derived features
    pass_df["af_germline_mode"] = np.maximum(
        1 - np.abs(af - 0.5) * 4,
        1 - np.abs(af - 1.0) * 4,
    ).clip(0, 1)

    pass_df["af_near_05"] = (np.abs(af - 0.5) < 0.1).astype(float)
    pass_df["af_near_10"] = (af > 0.85).astype(float)
    pass_df["h_flag_absent"] = 1 - pass_df["h_flag"].astype(float)
    pass_df["sb_deviation"] = np.abs(pass_df["sb_ratio"].values - 0.5) * 2

    pass_df["y_true"] = (pass_df["truth_label"] == "FP").astype(int)

    return pass_df


# ---------------------------------------------------------------------------
# Feature ranking
# ---------------------------------------------------------------------------


CANDIDATE_FEATURES = [
    "af",
    "af_germline_mode",
    "af_near_05",
    "af_near_10",
    "is_hom_alt",
    "is_het",
    "h_flag",
    "h_flag_absent",
    "is_cpg_ct",
    "is_cpg_site",
    "is_transition",
    "sb_clairs",
    "sb_ratio",
    "sb_deviation",
    "gq",
    "qual",
    "dp",
]


def fig01_feature_ranking(df: pd.DataFrame):
    """Feature importance bar chart."""
    _print("\n[G6] Fig01: Feature ranking")

    y = df["y_true"].values
    aucs = []
    for feat in CANDIDATE_FEATURES:
        vals = df[feat].values.astype(float)
        mask = np.isfinite(vals)
        auc = compute_auc(y[mask], vals[mask])
        # Flip if AUC < 0.5 (feature predicts TP, not FP)
        if auc < 0.5:
            auc_display = 1 - auc
            direction = "inverted"
        else:
            auc_display = auc
            direction = "direct"
        aucs.append({"feature": feat, "auc_raw": auc, "auc_abs": auc_display, "direction": direction})

    df_auc = pd.DataFrame(aucs).sort_values("auc_abs", ascending=True)

    fig, ax = plt.subplots(figsize=(10, 8))
    colors = ["#E65100" if r["auc_abs"] > 0.55 else "#757575" for _, r in df_auc.iterrows()]
    ax.barh(df_auc["feature"], df_auc["auc_abs"], color=colors, alpha=0.8)
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("AUC (for predicting FP)")
    ax.set_title("G6-Fig01: Individual Feature AUC Ranking",
                 fontsize=13, fontweight="bold")
    for i, (_, row) in enumerate(df_auc.iterrows()):
        ax.text(row["auc_abs"] + 0.005, i, f"{row['auc_abs']:.3f} ({row['direction']})",
                va="center", fontsize=8)

    add_caption(fig, "PASS variants, 7 samples pooled | Orange: AUC > 0.55")
    save_figure(fig, FIG_DIR / "G6_fig01_feature_ranking.png")

    df_auc.to_csv(DATA_DIR / "feature_auc_ranking.tsv", sep="\t", index=False)
    return df_auc


def fig02_feature_correlation(df: pd.DataFrame):
    """Feature correlation heatmap."""
    _print("\n[G6] Fig02: Feature correlation")

    feat_df = df[CANDIDATE_FEATURES].astype(float)
    corr = feat_df.corr()

    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                ax=ax, vmin=-1, vmax=1, square=True)
    ax.set_title("G6-Fig02: Feature Correlation Matrix",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "G6_fig02_feature_correlation.png")


def fig03_decision_tree(df: pd.DataFrame):
    """Decision tree (depth <= 3)."""
    _print("\n[G6] Fig03: Decision tree")
    from sklearn.tree import DecisionTreeClassifier, plot_tree
    from sklearn.model_selection import cross_val_score

    features_to_use = ["af", "af_germline_mode", "is_hom_alt", "h_flag_absent",
                        "is_cpg_ct", "is_transition", "gq", "qual"]
    X = df[features_to_use].astype(float).fillna(0).values
    y = df["y_true"].values

    # Train
    tree = DecisionTreeClassifier(max_depth=3, min_samples_leaf=100, random_state=42)
    tree.fit(X, y)

    # Cross-validation
    cv_scores = cross_val_score(tree, X, y, cv=5, scoring="roc_auc")
    _print(f"  [G6] Decision tree 5-fold CV AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

    fig, ax = plt.subplots(figsize=(20, 10))
    plot_tree(tree, feature_names=features_to_use,
             class_names=["TP", "FP"], filled=True, rounded=True,
             fontsize=8, ax=ax, impurity=True)
    ax.set_title(f"G6-Fig03: Decision Tree (depth≤3) — CV AUC={cv_scores.mean():.3f}±{cv_scores.std():.3f}",
                 fontsize=14, fontweight="bold")
    save_figure(fig, FIG_DIR / "G6_fig03_decision_tree.png")

    return tree, features_to_use


def fig04_roc_comparison(df: pd.DataFrame):
    """ROC comparison: individual features vs composite scores."""
    _print("\n[G6] Fig04: ROC comparison")
    from sklearn.metrics import roc_curve

    y = df["y_true"].values

    # Build composite scores
    scores = {}

    # Rule A: Probable germline
    scores["Rule A: probable_germline"] = (
        df["af_germline_mode"].values * 0.3 +
        df["h_flag_absent"].values * 0.25 +
        df["is_hom_alt"].values * 0.25 +
        df["is_cpg_ct"].values * 0.2
    )

    # Rule B: High-confidence germline hom
    scores["Rule B: hom_germline"] = (
        df["af_near_10"].values * 0.5 +
        df["is_hom_alt"].values * 0.5
    )

    # Rule C: Probable germline het
    scores["Rule C: het_germline"] = (
        df["af_near_05"].values * 0.4 +
        df["is_het"].values * 0.3 +
        df["h_flag_absent"].values * 0.3
    )

    # Combined
    scores["Combined (A+B+C)"] = np.maximum(
        scores["Rule A: probable_germline"],
        np.maximum(scores["Rule B: hom_germline"], scores["Rule C: het_germline"])
    )

    # Top individual features
    scores["AF (individual)"] = df["af"].values
    scores["af_germline_mode"] = df["af_germline_mode"].values
    scores["h_flag_absent"] = df["h_flag_absent"].values

    fig, ax = plt.subplots(figsize=(10, 8))
    for name, score in scores.items():
        mask = np.isfinite(score)
        auc = compute_auc(y[mask], score[mask])
        fpr, tpr, _ = roc_curve(y[mask], score[mask])
        lw = 2.5 if "Combined" in name else 1.5
        ls = "-" if "Rule" in name or "Combined" in name else "--"
        ax.plot(fpr, tpr, label=f"{name} (AUC={auc:.3f})", linewidth=lw, linestyle=ls)

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    ax.set_title("G6-Fig04: ROC — Individual vs Composite Scores",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, loc="lower right")
    save_figure(fig, FIG_DIR / "G6_fig04_roc_comparison.png")

    # Save AUCs
    auc_rows = []
    for name, score in scores.items():
        mask = np.isfinite(score)
        auc_rows.append({"score": name, "auc": compute_auc(y[mask], score[mask])})
    pd.DataFrame(auc_rows).to_csv(DATA_DIR / "composite_score_auc.tsv", sep="\t", index=False)

    return scores


def fig05_precision_recall_safety(df: pd.DataFrame, scores: dict):
    """Precision-recall with TP-loss constraint line."""
    _print("\n[G6] Fig05: Precision-recall with safety constraint")
    from sklearn.metrics import precision_recall_curve

    y = df["y_true"].values
    best_score = scores.get("Combined (A+B+C)", scores.get("Rule A: probable_germline"))

    mask = np.isfinite(best_score)

    fig, ax = plt.subplots(figsize=(10, 6))
    precision, recall, thresholds = precision_recall_curve(y[mask], best_score[mask])
    ax.plot(recall, precision, color="#E65100", linewidth=2, label="Combined score")

    # Mark 2% TP loss line
    # At each threshold, compute TP loss
    tp_loss_rates = []
    fp_removal_rates = []
    for t in np.linspace(0, 1, 100):
        predicted_fp = best_score[mask] >= t
        tp_flagged = (predicted_fp & (y[mask] == 0)).sum()
        tp_total = (y[mask] == 0).sum()
        fp_flagged = (predicted_fp & (y[mask] == 1)).sum()
        fp_total = (y[mask] == 1).sum()
        tp_loss = tp_flagged / tp_total if tp_total > 0 else 0
        fp_removal = fp_flagged / fp_total if fp_total > 0 else 0
        tp_loss_rates.append(tp_loss)
        fp_removal_rates.append(fp_removal)

    ax.set_xlabel("Recall (FP found)")
    ax.set_ylabel("Precision (of flagged)")
    ax.set_title("G6-Fig05: Precision-Recall with TP-Loss Constraint",
                 fontsize=13, fontweight="bold")
    ax.legend()
    add_caption(fig, "Red line: TP loss = 2% constraint")
    save_figure(fig, FIG_DIR / "G6_fig05_pr_safety.png")

    # Save TP-loss vs FP-removal data
    pd.DataFrame({"tp_loss": tp_loss_rates, "fp_removal": fp_removal_rates}).to_csv(
        DATA_DIR / "tp_loss_fp_removal.tsv", sep="\t", index=False)


def fig06_confusion_matrix(df: pd.DataFrame, scores: dict):
    """Confusion matrix at optimal threshold."""
    _print("\n[G6] Fig06: Confusion matrix")

    y = df["y_true"].values
    best_score = scores.get("Combined (A+B+C)")
    mask = np.isfinite(best_score)
    y_m = y[mask]
    s_m = best_score[mask]

    # Find threshold where TP loss <= 2%
    best_threshold = None
    best_fp_removed = 0
    for t in np.linspace(0, 1, 1000):
        predicted_fp = s_m >= t
        tp_flagged = (predicted_fp & (y_m == 0)).sum()
        tp_total = (y_m == 0).sum()
        fp_flagged = (predicted_fp & (y_m == 1)).sum()
        fp_total = (y_m == 1).sum()
        tp_loss = tp_flagged / tp_total
        if tp_loss <= 0.02 and fp_flagged > best_fp_removed:
            best_fp_removed = fp_flagged
            best_threshold = t

    if best_threshold is None:
        _print("  [G6] WARNING: No threshold found with TP loss <= 2%")
        best_threshold = 0.9

    predicted_fp = s_m >= best_threshold
    tp_flagged = (predicted_fp & (y_m == 0)).sum()
    fp_flagged = (predicted_fp & (y_m == 1)).sum()
    tp_kept = ((~predicted_fp) & (y_m == 0)).sum()
    fp_kept = ((~predicted_fp) & (y_m == 1)).sum()

    tp_total = (y_m == 0).sum()
    fp_total = (y_m == 1).sum()

    _print(f"  [G6] Optimal threshold: {best_threshold:.3f}")
    _print(f"  [G6] TP loss: {tp_flagged}/{tp_total} = {tp_flagged/tp_total*100:.2f}%")
    _print(f"  [G6] FP removed: {fp_flagged}/{fp_total} = {fp_flagged/fp_total*100:.2f}%")

    cm = np.array([[tp_kept, tp_flagged], [fp_kept, fp_flagged]])

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt=",d", cmap="Blues", ax=ax,
                xticklabels=["Kept", "Flagged as Germline"],
                yticklabels=["Actual TP", "Actual FP"])
    ax.set_ylabel("Truth")
    ax.set_xlabel("Prediction")
    ax.set_title(f"G6-Fig06: Confusion Matrix (threshold={best_threshold:.3f})\n"
                 f"TP loss={tp_flagged/tp_total*100:.1f}%, FP removed={fp_flagged/fp_total*100:.1f}%",
                 fontsize=12, fontweight="bold")
    save_figure(fig, FIG_DIR / "G6_fig06_confusion_matrix.png")

    return best_threshold


def fig07_score_distribution(df: pd.DataFrame, scores: dict):
    """Composite score distribution."""
    _print("\n[G6] Fig07: Score distribution")

    y = df["y_true"].values
    best_score = scores.get("Combined (A+B+C)")

    fig, ax = plt.subplots(figsize=(10, 5))
    for label, mask_val, color in [("TP", 0, COLOR_TP), ("FP", 1, COLOR_FP)]:
        vals = best_score[y == mask_val]
        vals = vals[np.isfinite(vals)]
        ax.hist(vals, bins=50, alpha=0.5, color=color, label=f"{label} (n={len(vals):,})",
                density=True, edgecolor="white")

    ax.set_xlabel("Combined Germline Score")
    ax.set_ylabel("Density")
    ax.set_title("G6-Fig07: Composite Germline Score Distribution",
                 fontsize=13, fontweight="bold")
    ax.legend()
    save_figure(fig, FIG_DIR / "G6_fig07_score_distribution.png")


def fig08_pareto_front(df: pd.DataFrame, scores: dict):
    """Pareto front: FP removal vs TP loss."""
    _print("\n[G6] Fig08: Pareto front")

    y = df["y_true"].values
    best_score = scores.get("Combined (A+B+C)")
    mask = np.isfinite(best_score)
    y_m = y[mask]
    s_m = best_score[mask]

    tp_total = (y_m == 0).sum()
    fp_total = (y_m == 1).sum()

    points = []
    for t in np.linspace(0, 1, 200):
        predicted_fp = s_m >= t
        tp_loss = (predicted_fp & (y_m == 0)).sum() / tp_total * 100
        fp_removal = (predicted_fp & (y_m == 1)).sum() / fp_total * 100
        points.append({"threshold": t, "tp_loss": tp_loss, "fp_removal": fp_removal})

    df_pareto = pd.DataFrame(points)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(df_pareto["tp_loss"], df_pareto["fp_removal"], c=df_pareto["threshold"],
              cmap="YlOrRd", s=20, alpha=0.8)
    ax.axvline(2.0, color="red", linestyle="--", alpha=0.7, label="2% TP loss limit")
    ax.set_xlabel("TP Loss (%)")
    ax.set_ylabel("FP Removed (%)")
    ax.set_title("G6-Fig08: Pareto Front — FP Removal vs TP Loss",
                 fontsize=13, fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.colorbar(ax.collections[0], label="Threshold")
    save_figure(fig, FIG_DIR / "G6_fig08_pareto_front.png")

    df_pareto.to_csv(DATA_DIR / "pareto_front.tsv", sep="\t", index=False)


def fig09_cross_validation(df: pd.DataFrame):
    """Leave-one-sample-out cross-validation."""
    _print("\n[G6] Fig09: Cross-validation")
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    features = ["af", "af_germline_mode", "is_hom_alt", "h_flag_absent",
                "is_cpg_ct", "is_transition", "gq"]

    X = df[features].astype(float).fillna(0).values
    y = df["y_true"].values
    samples = df["sample"].values

    aucs = []
    for test_sample in SAMPLE_ORDER:
        train_mask = samples != test_sample
        test_mask = samples == test_sample
        if test_mask.sum() == 0 or train_mask.sum() == 0:
            continue

        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_mask])
        X_test = scaler.transform(X[test_mask])

        lr = LogisticRegression(random_state=42, max_iter=1000)
        lr.fit(X_train, y[train_mask])
        y_score = lr.predict_proba(X_test)[:, 1]
        auc = compute_auc(y[test_mask], y_score)
        aucs.append({"sample": test_sample, "auc": auc})
        _print(f"    LOSO {test_sample}: AUC = {auc:.3f}")

    df_cv = pd.DataFrame(aucs)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(df_cv["sample"], df_cv["auc"], color="#E65100", alpha=0.8)
    ax.axhline(df_cv["auc"].mean(), color="red", linestyle="--",
               label=f"Mean AUC = {df_cv['auc'].mean():.3f}")
    ax.set_ylabel("AUC")
    ax.set_title("G6-Fig09: Leave-One-Sample-Out AUC",
                 fontsize=13, fontweight="bold")
    ax.legend()
    ax.set_ylim(0.4, 1.0)
    plt.xticks(rotation=30, ha="right")
    save_figure(fig, FIG_DIR / "G6_fig09_cross_validation.png")

    df_cv.to_csv(DATA_DIR / "loso_auc.tsv", sep="\t", index=False)


def fig10_per_sample_performance(df: pd.DataFrame, scores: dict, threshold: float):
    """Per-sample performance bar chart."""
    _print("\n[G6] Fig10: Per-sample performance")

    best_score = scores.get("Combined (A+B+C)")
    y = df["y_true"].values

    rows = []
    for sample in SAMPLE_ORDER:
        mask = (df["sample"] == sample).values & np.isfinite(best_score)
        if mask.sum() == 0:
            continue
        y_s = y[mask]
        s_s = best_score[mask]
        predicted_fp = s_s >= threshold

        tp_total = (y_s == 0).sum()
        fp_total = (y_s == 1).sum()
        tp_loss = (predicted_fp & (y_s == 0)).sum() / tp_total * 100 if tp_total > 0 else 0
        fp_removal = (predicted_fp & (y_s == 1)).sum() / fp_total * 100 if fp_total > 0 else 0
        auc = compute_auc(y_s, s_s)

        rows.append({
            "sample": sample, "auc": auc,
            "tp_loss_pct": tp_loss, "fp_removal_pct": fp_removal,
            "tp_total": tp_total, "fp_total": fp_total,
        })

    df_perf = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    ax.bar(df_perf["sample"], df_perf["auc"], color="#E65100", alpha=0.8)
    ax.set_ylabel("AUC")
    ax.set_title("Combined Score AUC")
    ax.set_ylim(0.4, 1.0)
    plt.sca(ax)
    plt.xticks(rotation=30, ha="right")

    ax = axes[1]
    ax.bar(df_perf["sample"], df_perf["fp_removal_pct"], color=COLOR_FP, alpha=0.8)
    ax.set_ylabel("FP Removed (%)")
    ax.set_title(f"FP Removal (threshold={threshold:.3f})")
    plt.sca(ax)
    plt.xticks(rotation=30, ha="right")

    ax = axes[2]
    ax.bar(df_perf["sample"], df_perf["tp_loss_pct"], color=COLOR_TP, alpha=0.8)
    ax.axhline(2.0, color="red", linestyle="--", label="2% limit")
    ax.set_ylabel("TP Loss (%)")
    ax.set_title("TP Loss (safety check)")
    ax.legend()
    plt.sca(ax)
    plt.xticks(rotation=30, ha="right")

    fig.suptitle("G6-Fig10: Per-Sample Performance",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, FIG_DIR / "G6_fig10_per_sample_performance.png")

    df_perf.to_csv(DATA_DIR / "per_sample_performance.tsv", sep="\t", index=False)


def main():
    setup_plot_style()
    _print("=" * 70)
    _print("G6: Multi-Feature Combination")
    _print("=" * 70)

    df_raw = load_pass_features()
    df = prepare_features(df_raw)

    auc_ranking = fig01_feature_ranking(df)
    fig02_feature_correlation(df)
    tree, tree_features = fig03_decision_tree(df)
    scores = fig04_roc_comparison(df)
    fig05_precision_recall_safety(df, scores)
    threshold = fig06_confusion_matrix(df, scores)
    fig07_score_distribution(df, scores)
    fig08_pareto_front(df, scores)
    fig09_cross_validation(df)
    fig10_per_sample_performance(df, scores, threshold)

    _print(f"\n[G6] DONE. Output: {G6_DIR}")


if __name__ == "__main__":
    main()
