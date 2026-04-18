#!/usr/bin/env python3
"""TO ClairS-TO Feature Discrimination Deep Study — Q1-Q5.

Compares longphase-to baseline vs PON-only v2b haplotagging,
both using correct ClairS-TO VCFs.

Q1: HPFineNGroups LOH data integrity
Q2: LOH_Weak HP Read Count Filter (AUC 0.72) — 3-level deconfounding
Q3: AUC heatmap red signals (FP > TP reversed features)
Q4: HPFineSig and AlleleSig in TO
Q5: AlleleDelta LOH inner/outer discrimination
"""
from __future__ import annotations

import sys
import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    setup_plot_style, save_figure, add_caption,
    compute_auc, bootstrap_ci, compute_enrichment,
    mannwhitney_test, chi2_test, ensure_dir,
    COLOR_TP, COLOR_FP, format_p, format_ci,
)

warnings.filterwarnings("ignore", category=FutureWarning)

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

OUT_ROOT = ensure_dir(Path("/big7_disk/liaoyoyo2001/InterSubMod/output/clairsto_correction_analysis/feature_study"))
DATA_OUT = ensure_dir(OUT_ROOT / "data")

# Feature lists
CONTINUOUS_FEATURES = [
    "NumReads", "NumCpGs", "GlobalP", "CramersV",
    "GlobalP_HPFamily", "CramersV_HPFamily", "GlobalP_HPFine", "CramersV_HPFine",
    "HPFine_NGroups_CF", "HeuristicScore", "PairwiseMeanDist", "PairwiseMedianDist",
    "ClusterPermanovaF", "ClusterPermanovaP",
    "HPMergedDelta", "HPMergedP", "HP1FamilyN", "HP2FamilyN",
    "HPFineF", "HPFineP", "HPFineNGroups",
    "AlleleDelta", "AlleleP",
    "LabelHPPermanovaF", "LabelHPPermanovaP",
    "LabelAllelePermanovaF", "LabelAllelePermanovaP",
    "UnassignedAffinity", "UnassignedAffinityP",
    "NHP3", "NHP0", "Stability",
    "HP_Ratio", "Coverage_Multiple", "Quality_Score",
]

BOOLEAN_FEATURES = [
    "PassedGating", "HPMergedSig", "HPFineSig", "AlleleSig",
    "ClusterPermanovaValid", "Potential_LOH",
]

STRATA = ["All", "LOH", "Non-LOH", "LOH_Strong", "LOH_Weak"]

def _p(msg): print(msg, flush=True)


# ---------------------------------------------------------------------------
# Data Loading
# ---------------------------------------------------------------------------
def load_version(version: str) -> pd.DataFrame:
    """Load TP + FP for one haplotagging version, add truth_label."""
    paths = DATA_PATHS[version]
    tp = pd.read_csv(paths["tp"], low_memory=False)
    fp = pd.read_csv(paths["fp"], low_memory=False)
    tp["truth_label"] = "TP"
    fp["truth_label"] = "FP"
    tp["version"] = version
    fp["version"] = version
    df = pd.concat([tp, fp], ignore_index=True)
    # Convert booleans
    for col in BOOLEAN_FEATURES + ["Potential_LOH"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.lower().isin(["true", "1", "yes"]).astype(int)
    # Ensure numeric
    for col in CONTINUOUS_FEATURES:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["y_true"] = (df["truth_label"] == "TP").astype(int)
    return df


def stratify(df: pd.DataFrame, stratum: str) -> pd.DataFrame:
    """Filter df by stratum."""
    if stratum == "All":
        return df
    elif stratum == "LOH":
        return df[df["Potential_LOH"] == 1]
    elif stratum == "Non-LOH":
        return df[df["Potential_LOH"] == 0]
    elif stratum == "LOH_Strong":
        return df[df["LOH_Subtype"] == "LOH_Strong"]
    elif stratum == "LOH_Weak":
        return df[df["LOH_Subtype"] == "LOH_Weak"]
    return df


def bootstrap_auc(y_true, y_score, n_boot=500, seed=42):
    """Bootstrap AUC with 95% CI."""
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
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
        if not np.isnan(a):
            aucs.append(a)
    if len(aucs) < 10:
        return point, float("nan"), float("nan")
    return point, float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


# ===========================================================================
# Q1: HPFineNGroups LOH Data Integrity
# ===========================================================================
def run_q1(dfs: dict[str, pd.DataFrame]):
    _p("\n" + "="*60)
    _p("Q1: HPFineNGroups LOH Data Integrity")
    _p("="*60)
    q1_dir = ensure_dir(OUT_ROOT / "Q1")

    results = {}
    for ver, df in dfs.items():
        _p(f"\n--- {ver} ---")
        res = {"version": ver}

        # 1. NGroups (0-4) × label contingency (LOH only)
        loh = stratify(df, "LOH")
        _p(f"  LOH: {len(loh):,} rows (TP={sum(loh.y_true==1):,}, FP={sum(loh.y_true==0):,})")

        ng_col = "HPFineNGroups"
        ct = pd.crosstab(loh[ng_col].clip(0, 4), loh["truth_label"])
        _p(f"  NGroups × label contingency:\n{ct.to_string()}")
        chi2_res = chi2_test(ct.values)
        _p(f"  Chi-squared: {chi2_res['chi2']:.2f}, p={format_p(chi2_res['p_value'])}, V={chi2_res['cramers_v']:.4f}")
        res["chi2"] = chi2_res

        # Per-bin Fisher exact
        per_bin = {}
        for ng in range(5):
            mask = loh[ng_col].clip(0, 4) == ng
            tp_c = int(((loh.y_true == 1) & mask).sum())
            fp_c = int(((loh.y_true == 0) & mask).sum())
            tp_tot = int((loh.y_true == 1).sum())
            fp_tot = int((loh.y_true == 0).sum())
            enr = compute_enrichment(tp_c, fp_c, tp_tot, fp_tot)
            per_bin[ng] = enr
            _p(f"  NGroups={ng}: TP={tp_c}({tp_c/tp_tot*100:.1f}%), FP={fp_c}({fp_c/fp_tot*100:.1f}%), enrichment={enr['enrichment']:.3f}")
        res["per_bin"] = per_bin

        # NGroups≤1 FP/TP ratio
        tp_le1 = ((loh.y_true == 1) & (loh[ng_col] <= 1)).sum() / max((loh.y_true == 1).sum(), 1)
        fp_le1 = ((loh.y_true == 0) & (loh[ng_col] <= 1)).sum() / max((loh.y_true == 0).sum(), 1)
        ratio_le1 = fp_le1 / tp_le1 if tp_le1 > 0 else float("nan")
        _p(f"  NGroups≤1: TP rate={tp_le1:.4f}, FP rate={fp_le1:.4f}, ratio={ratio_le1:.3f}")
        res["ngroups_le1_ratio"] = ratio_le1

        # 2. Per-chromosome NGroups distribution stability
        chr_stats = []
        for chrom in sorted(loh["Chr"].unique()):
            cdf = loh[loh["Chr"] == chrom]
            if len(cdf) < 20:
                continue
            tp_n = (cdf.y_true == 1).sum()
            fp_n = (cdf.y_true == 0).sum()
            if tp_n < 5 or fp_n < 5:
                continue
            tp_le1_c = ((cdf.y_true == 1) & (cdf[ng_col] <= 1)).sum() / max(tp_n, 1)
            fp_le1_c = ((cdf.y_true == 0) & (cdf[ng_col] <= 1)).sum() / max(fp_n, 1)
            r = fp_le1_c / tp_le1_c if tp_le1_c > 0 else float("nan")
            chr_stats.append({"chr": chrom, "n": len(cdf), "tp_n": int(tp_n), "fp_n": int(fp_n),
                              "tp_le1_rate": float(tp_le1_c), "fp_le1_rate": float(fp_le1_c), "ratio": float(r)})
        chr_df = pd.DataFrame(chr_stats)
        chr_df.to_csv(q1_dir / f"q1_per_chr_ngroups_{ver}.tsv", sep="\t", index=False)
        res["chr_stats"] = chr_stats
        _p(f"  Per-chromosome: {len(chr_stats)} chromosomes analyzed, ratio range [{chr_df.ratio.min():.2f}, {chr_df.ratio.max():.2f}]")

        # Check single chromosome dominance
        if len(chr_df) > 0:
            max_chr = chr_df.loc[chr_df["n"].idxmax()]
            max_pct = max_chr["n"] / len(loh) * 100
            _p(f"  Largest chromosome: {max_chr['chr']} ({max_pct:.1f}% of LOH data)")
            res["max_chr_pct"] = float(max_pct)

        # 3. HPFineNGroups AUC + bootstrap CI
        for stratum in ["LOH", "Non-LOH", "All"]:
            sdf = stratify(df, stratum)
            auc_pt, auc_lo, auc_hi = bootstrap_auc(sdf.y_true.values, sdf[ng_col].values * -1)  # lower NGroups → more FP → negate
            _p(f"  {ng_col} AUC ({stratum}): {auc_pt:.4f} {format_ci(auc_lo, auc_hi)}")
            res[f"auc_{stratum}"] = {"auc": auc_pt, "ci_lo": auc_lo, "ci_hi": auc_hi}

        # Also compute HPFine_NGroups_CF AUC
        cf_col = "HPFine_NGroups_CF"
        for stratum in ["LOH", "Non-LOH", "All"]:
            sdf = stratify(df, stratum)
            auc_pt, auc_lo, auc_hi = bootstrap_auc(sdf.y_true.values, sdf[cf_col].values * -1)
            _p(f"  {cf_col} AUC ({stratum}): {auc_pt:.4f} {format_ci(auc_lo, auc_hi)}")
            res[f"auc_cf_{stratum}"] = {"auc": auc_pt, "ci_lo": auc_lo, "ci_hi": auc_hi}

        results[ver] = res

    # --- Q1 Figures ---
    setup_plot_style()

    # Q1-Fig1: NGroups distribution (LOH vs Non-LOH, TP vs FP) — baseline and v2b side-by-side
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for col_idx, ver in enumerate(dfs.keys()):
        df = dfs[ver]
        for row_idx, stratum in enumerate(["LOH", "Non-LOH"]):
            ax = axes[row_idx, col_idx]
            sdf = stratify(df, stratum)
            for label, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
                vals = sdf[sdf.truth_label == label]["HPFineNGroups"].clip(0, 4)
                counts = vals.value_counts().sort_index()
                pcts = counts / len(vals) * 100
                ax.bar(pcts.index + (-0.15 if label == "TP" else 0.15),
                       pcts.values, width=0.3, color=color, alpha=0.8, label=label)
            ax.set_xlabel("HPFineNGroups")
            ax.set_ylabel("Percentage (%)")
            ax.set_title(f"{ver} — {stratum} (n={len(sdf):,})")
            ax.set_xticks(range(5))
            ax.legend()
    fig.suptitle("Q1: HPFineNGroups Distribution by LOH Status and Haplotagging Version", fontsize=14)
    fig.tight_layout()
    add_caption(fig, "Data: ClairS-TO VCF, HCC1395 100% purity")
    save_figure(fig, q1_dir / "Q1_fig1_ngroups_distribution.png")

    # Q1-Fig2: Per-chromosome NGroups=0 FP/TP ratio bar chart
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    for idx, ver in enumerate(dfs.keys()):
        ax = axes[idx]
        cdf = pd.DataFrame(results[ver]["chr_stats"])
        if len(cdf) == 0:
            continue
        # Sort by chromosome number
        cdf["chr_num"] = cdf["chr"].str.replace("chr", "").apply(
            lambda x: int(x) if x.isdigit() else 23 if x == "X" else 24)
        cdf = cdf.sort_values("chr_num")
        colors = ["#F44336" if r > 2.0 else "#FF9800" if r > 1.5 else "#4CAF50" for r in cdf["ratio"]]
        ax.bar(range(len(cdf)), cdf["ratio"].values, color=colors, alpha=0.8)
        ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)
        ax.axhline(y=results[ver]["ngroups_le1_ratio"], color="blue", linestyle=":", alpha=0.7, label=f"Global ratio={results[ver]['ngroups_le1_ratio']:.2f}")
        ax.set_xticks(range(len(cdf)))
        ax.set_xticklabels(cdf["chr"].str.replace("chr", ""), rotation=45, fontsize=7)
        ax.set_ylabel("NGroups≤1 FP/TP Ratio")
        ax.set_title(f"{ver}")
        ax.legend(fontsize=8)
    fig.suptitle("Q1: Per-Chromosome NGroups≤1 FP/TP Ratio (LOH)", fontsize=13)
    fig.tight_layout()
    save_figure(fig, q1_dir / "Q1_fig2_per_chr_ngroups_ratio.png")

    # Save results
    with open(q1_dir / "q1_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    _p("Q1 complete.")
    return results


# ===========================================================================
# Q2: LOH_Weak HP Read Count Filter — 3-Level Deconfounding
# ===========================================================================
def run_q2(dfs: dict[str, pd.DataFrame]):
    _p("\n" + "="*60)
    _p("Q2: LOH_Weak HP Read Count Filter — 3-Level Deconfounding")
    _p("="*60)
    q2_dir = ensure_dir(OUT_ROOT / "Q2")

    results = {}
    for ver, df in dfs.items():
        _p(f"\n--- {ver} ---")
        res = {"version": ver}

        # Three-level deconfounding AUC tests
        features_to_test = ["HP2FamilyN", "HP_Ratio", "HP1FamilyN"]
        for feat in features_to_test:
            _p(f"\n  Feature: {feat}")
            for test_name, stratum in [("A_LOH_Weak", "LOH_Weak"), ("B_All_LOH", "LOH"), ("C_All", "All")]:
                sdf = stratify(df, stratum)
                n_tp = (sdf.y_true == 1).sum()
                n_fp = (sdf.y_true == 0).sum()
                if feat == "HP_Ratio":
                    # HP_Ratio: LOH = very low or very high, so use distance from 0.5
                    score = np.abs(sdf[feat].values - 0.5)
                else:
                    score = sdf[feat].values
                auc_pt, auc_lo, auc_hi = bootstrap_auc(sdf.y_true.values, score)
                # Also try negated
                auc_neg, _, _ = bootstrap_auc(sdf.y_true.values, -score)
                best_auc = max(auc_pt, auc_neg) if not np.isnan(auc_neg) else auc_pt
                best_dir = "+" if auc_pt >= auc_neg else "-"
                _p(f"    Test {test_name} (n={len(sdf):,}, TP={n_tp}, FP={n_fp}): AUC={best_auc:.4f} (dir={best_dir}) {format_ci(auc_lo, auc_hi)}")
                res[f"{feat}_{test_name}"] = {"auc": best_auc, "ci_lo": auc_lo, "ci_hi": auc_hi,
                                               "n": len(sdf), "n_tp": int(n_tp), "n_fp": int(n_fp), "direction": best_dir}

        # HP2FamilyN threshold scan (LOH_Weak)
        loh_w = stratify(df, "LOH_Weak")
        if len(loh_w) > 0:
            _p(f"\n  HP2FamilyN threshold scan (LOH_Weak, n={len(loh_w):,}):")
            truth_total = (df.y_true == 1).sum()  # all TP as denominator for global F1
            thresholds = [1, 2, 3, 5, 8, 10, 15, 20, 30]
            scan_rows = []
            for thr in thresholds:
                # Filter: HP2FamilyN < thr → mark as FP
                flagged = loh_w["HP2FamilyN"] < thr
                tp_flagged = ((loh_w.y_true == 1) & flagged).sum()
                fp_flagged = ((loh_w.y_true == 0) & flagged).sum()
                tp_total = (loh_w.y_true == 1).sum()
                fp_total = (loh_w.y_true == 0).sum()
                tp_loss_pct = tp_flagged / max(tp_total, 1) * 100
                fp_removal_pct = fp_flagged / max(fp_total, 1) * 100
                precision = fp_flagged / max(tp_flagged + fp_flagged, 1)
                scan_rows.append({
                    "threshold": thr, "tp_flagged": int(tp_flagged), "fp_flagged": int(fp_flagged),
                    "tp_loss_pct": float(tp_loss_pct), "fp_removal_pct": float(fp_removal_pct),
                    "precision_of_filter": float(precision),
                })
                _p(f"    thr<{thr}: TP_loss={tp_loss_pct:.1f}%, FP_removal={fp_removal_pct:.1f}%, precision={precision:.3f}")
            scan_df = pd.DataFrame(scan_rows)
            scan_df.to_csv(q2_dir / f"q2_threshold_scan_{ver}.tsv", sep="\t", index=False)
            res["threshold_scan"] = scan_rows

        results[ver] = res

    # --- Q2 Figures ---
    setup_plot_style()

    # Q2-Fig1: Circular dependency DAG (conceptual)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlim(-0.5, 5.5)
    ax.set_ylim(-0.5, 3.5)
    boxes = {
        "HP1FamilyN\nHP2FamilyN": (0.5, 2.5),
        "HP_Ratio": (2, 2.5),
        "Potential_LOH": (3.5, 2.5),
        "HPMergedSig": (0.5, 0.5),
        "VerificationClass": (2, 0.5),
        "LOH_Subtype": (3.5, 0.5),
        "LOH_Weak": (5, 1.5),
    }
    for label, (x, y) in boxes.items():
        color = "#FFCDD2" if "LOH_Weak" in label else "#E3F2FD" if "HP" in label else "#FFF9C4"
        ax.add_patch(plt.Rectangle((x-0.45, y-0.3), 0.9, 0.6, facecolor=color, edgecolor="black", linewidth=1.5))
        ax.text(x, y, label, ha="center", va="center", fontsize=7, weight="bold")
    arrows = [
        ((0.95, 2.5), (1.55, 2.5)),
        ((2.45, 2.5), (3.05, 2.5)),
        ((0.95, 0.5), (1.55, 0.5)),
        ((2.45, 0.5), (3.05, 0.5)),
        ((3.5, 2.2), (3.5, 0.8)),
        ((2, 0.8), (2, 2.2)),
        ((3.95, 0.5), (4.55, 1.2)),
        ((3.95, 2.5), (4.55, 1.8)),
    ]
    for (x1, y1), (x2, y2) in arrows:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", color="black", lw=1.5))
    ax.set_title("Q2: Circular Dependency DAG — LOH_Weak Definition", fontsize=12)
    ax.axis("off")
    save_figure(fig, q2_dir / "Q2_fig1_circular_dependency_dag.png")

    # Q2-Fig2: HP2FamilyN distribution + ROC 3-line comparison
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for col_idx, ver in enumerate(dfs.keys()):
        df = dfs[ver]
        # Distribution
        ax = axes[0, col_idx]
        for stratum, ls in [("LOH_Weak", "-"), ("LOH", "--")]:
            sdf = stratify(df, stratum)
            for label, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
                vals = sdf[sdf.truth_label == label]["HP2FamilyN"].clip(0, 100)
                ax.hist(vals, bins=50, density=True, alpha=0.3, color=color, label=f"{label} ({stratum})", histtype="stepfilled", linestyle=ls)
        ax.set_xlabel("HP2FamilyN")
        ax.set_ylabel("Density")
        ax.set_title(f"{ver} — HP2FamilyN Distribution")
        ax.legend(fontsize=7)

        # ROC curves for 3 tests
        ax = axes[1, col_idx]
        from sklearn.metrics import roc_curve
        for test_name, stratum, color, ls in [
            ("A: LOH_Weak", "LOH_Weak", "#D32F2F", "-"),
            ("B: All LOH", "LOH", "#1976D2", "--"),
            ("C: All Data", "All", "#388E3C", ":"),
        ]:
            sdf = stratify(df, stratum)
            y = sdf.y_true.values
            s = sdf["HP2FamilyN"].values
            mask = np.isfinite(s)
            y, s = y[mask], s[mask]
            if len(np.unique(y)) >= 2:
                fpr, tpr, _ = roc_curve(y, s)
                auc_v = compute_auc(y, s)
                ax.plot(fpr, tpr, color=color, linestyle=ls, label=f"{test_name}: AUC={auc_v:.3f}")
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
        ax.set_xlabel("FPR")
        ax.set_ylabel("TPR")
        ax.set_title(f"{ver} — 3-Level Deconfounding ROC")
        ax.legend(fontsize=8)

    fig.suptitle("Q2: HP2FamilyN — Distribution & 3-Level Deconfounding", fontsize=14)
    fig.tight_layout()
    save_figure(fig, q2_dir / "Q2_fig2_hp2familyn_roc.png")

    # Q2-Fig3: Threshold scan — Precision-Recall-TP_loss curves
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for idx, ver in enumerate(dfs.keys()):
        ax = axes[idx]
        if "threshold_scan" not in results[ver]:
            continue
        scan = pd.DataFrame(results[ver]["threshold_scan"])
        ax.plot(scan["threshold"], scan["tp_loss_pct"], "o-", color=COLOR_TP, label="TP Loss %")
        ax.plot(scan["threshold"], scan["fp_removal_pct"], "s-", color=COLOR_FP, label="FP Removal %")
        ax.plot(scan["threshold"], scan["precision_of_filter"] * 100, "^-", color="#9C27B0", label="Precision %")
        ax.axhline(y=2.0, color="gray", linestyle="--", alpha=0.5, label="2% TP Loss Limit")
        ax.set_xlabel("HP2FamilyN Threshold (filter if < threshold)")
        ax.set_ylabel("Percentage (%)")
        ax.set_title(f"{ver} — LOH_Weak Filter Performance")
        ax.legend(fontsize=8)
        ax.set_ylim(0, 100)
    fig.suptitle("Q2: HP2FamilyN Threshold Scan — TP Loss vs FP Removal (LOH_Weak)", fontsize=13)
    fig.tight_layout()
    save_figure(fig, q2_dir / "Q2_fig3_threshold_scan.png")

    with open(q2_dir / "q2_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    _p("Q2 complete.")
    return results


# ===========================================================================
# Q3: AUC Heatmap Red Signals (FP > TP Reversed Features)
# ===========================================================================
def run_q3(dfs: dict[str, pd.DataFrame]):
    _p("\n" + "="*60)
    _p("Q3: AUC Heatmap — Reversed Features (AUC < 0.48)")
    _p("="*60)
    q3_dir = ensure_dir(OUT_ROOT / "Q3")

    results = {}
    for ver, df in dfs.items():
        _p(f"\n--- {ver} ---")
        res = {"version": ver}

        # Full feature × strata AUC heatmap
        auc_matrix = []
        for feat in CONTINUOUS_FEATURES:
            row = {"feature": feat}
            for stratum in STRATA:
                sdf = stratify(df, stratum)
                auc_v = compute_auc(sdf.y_true.values, sdf[feat].values)
                row[stratum] = auc_v
            auc_matrix.append(row)
        auc_df = pd.DataFrame(auc_matrix).set_index("feature")
        auc_df.to_csv(q3_dir / f"q3_auc_matrix_{ver}.tsv", sep="\t")
        res["auc_matrix"] = auc_df.to_dict()

        # Identify reversed features (AUC < 0.48)
        reversed_feats = []
        for feat in CONTINUOUS_FEATURES:
            for stratum in STRATA:
                auc_v = auc_df.loc[feat, stratum]
                if not np.isnan(auc_v) and auc_v < 0.48:
                    sdf = stratify(df, stratum)
                    # Bootstrap CI
                    auc_pt, auc_lo, auc_hi = bootstrap_auc(sdf.y_true.values, sdf[feat].values)
                    reversed_feats.append({
                        "feature": feat, "stratum": stratum,
                        "auc": float(auc_v), "ci_lo": auc_lo, "ci_hi": auc_hi,
                    })
        rev_df = pd.DataFrame(reversed_feats)
        if len(rev_df) > 0:
            rev_df = rev_df.sort_values("auc")
            rev_df.to_csv(q3_dir / f"q3_reversed_features_{ver}.tsv", sep="\t", index=False)
            _p(f"  Reversed features (AUC<0.48): {len(rev_df)}")
            for _, r in rev_df.head(10).iterrows():
                _p(f"    {r['feature']} @ {r['stratum']}: AUC={r['auc']:.4f} {format_ci(r['ci_lo'], r['ci_hi'])}")
        res["reversed"] = reversed_feats

        # AlleleSig Non-LOH deep-dive
        nloh = stratify(df, "Non-LOH")
        allele_sig_rate = {}
        for label in ["TP", "FP"]:
            sub = nloh[nloh.truth_label == label]
            rate = sub["AlleleSig"].mean() if len(sub) > 0 else float("nan")
            allele_sig_rate[label] = float(rate)
        _p(f"  AlleleSig Non-LOH: TP={allele_sig_rate['TP']:.4f}, FP={allele_sig_rate['FP']:.4f}")
        res["allelesig_nonloh"] = allele_sig_rate

        # Per-chromosome stability for top reversed features
        if len(rev_df) > 0:
            top_rev = rev_df.head(3)["feature"].unique()
            for feat in top_rev:
                chr_aucs = []
                for chrom in sorted(df["Chr"].unique()):
                    cdf = df[df["Chr"] == chrom]
                    if len(cdf) < 50 or len(cdf[cdf.y_true == 1]) < 10 or len(cdf[cdf.y_true == 0]) < 10:
                        continue
                    auc_v = compute_auc(cdf.y_true.values, cdf[feat].values)
                    chr_aucs.append({"chr": chrom, "auc": float(auc_v)})
                if chr_aucs:
                    _p(f"  {feat} per-chr AUC range: [{min(c['auc'] for c in chr_aucs):.3f}, {max(c['auc'] for c in chr_aucs):.3f}]")
                    res[f"chr_stability_{feat}"] = chr_aucs

        results[ver] = res

    # --- Q3 Figures ---
    setup_plot_style()

    # Q3-Fig1: Directional AUC heatmap (green=TP>FP, red=FP>TP)
    fig, axes = plt.subplots(1, 2, figsize=(12, max(10, len(CONTINUOUS_FEATURES) * 0.35)))
    for idx, ver in enumerate(dfs.keys()):
        ax = axes[idx]
        auc_df = pd.DataFrame(results[ver]["auc_matrix"])
        # Custom diverging colormap centered at 0.5
        from matplotlib.colors import TwoSlopeNorm
        norm = TwoSlopeNorm(vmin=0.35, vcenter=0.5, vmax=0.75)
        cmap = plt.cm.RdYlGn  # Red=low (FP>TP), Green=high (TP>FP)
        im = ax.imshow(auc_df.values, aspect="auto", cmap=cmap, norm=norm)
        ax.set_xticks(range(len(auc_df.columns)))
        ax.set_xticklabels(auc_df.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(auc_df.index)))
        ax.set_yticklabels(auc_df.index, fontsize=7)
        # Annotate cells
        for i in range(len(auc_df.index)):
            for j in range(len(auc_df.columns)):
                v = auc_df.values[i, j]
                if not np.isnan(v):
                    color = "white" if abs(v - 0.5) > 0.1 else "black"
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.5, color=color)
        ax.set_title(f"{ver}", fontsize=11)
        plt.colorbar(im, ax=ax, shrink=0.6, label="AUC")
    fig.suptitle("Q3: Feature AUC Heatmap (Green=TP favored, Red=FP favored)", fontsize=13)
    fig.tight_layout()
    save_figure(fig, q3_dir / "Q3_fig1_directional_auc_heatmap.png")

    # Q3-Fig2: Top 4 reversed features KDE
    for ver, df in dfs.items():
        rev_feats = [r["feature"] for r in results[ver]["reversed"][:4]]
        rev_feats = list(dict.fromkeys(rev_feats))[:4]  # unique, keep order
        if len(rev_feats) == 0:
            continue
        fig, axes = plt.subplots(1, min(4, len(rev_feats)), figsize=(4 * min(4, len(rev_feats)), 4))
        if len(rev_feats) == 1:
            axes = [axes]
        for i, feat in enumerate(rev_feats[:4]):
            ax = axes[i]
            for label, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
                vals = df[df.truth_label == label][feat].dropna()
                if len(vals) > 10:
                    ax.hist(vals, bins=50, density=True, alpha=0.4, color=color, label=label)
            auc_v = compute_auc(df.y_true.values, df[feat].values)
            ax.set_title(f"{feat}\nAUC={auc_v:.3f}", fontsize=9)
            ax.legend(fontsize=7)
        fig.suptitle(f"Q3: Top Reversed Features KDE — {ver}", fontsize=12)
        fig.tight_layout()
        save_figure(fig, q3_dir / f"Q3_fig2_reversed_kde_{ver}.png")

    with open(q3_dir / "q3_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    _p("Q3 complete.")
    return results


# ===========================================================================
# Q4: HPFineSig and AlleleSig in TO
# ===========================================================================
def run_q4(dfs: dict[str, pd.DataFrame]):
    _p("\n" + "="*60)
    _p("Q4: HPFineSig and AlleleSig in TO")
    _p("="*60)
    q4_dir = ensure_dir(OUT_ROOT / "Q4")

    results = {}
    for ver, df in dfs.items():
        _p(f"\n--- {ver} ---")
        res = {"version": ver}

        # 1. HPFineSig=false as LOH-only filter
        for stratum in ["LOH", "LOH_Strong", "LOH_Weak"]:
            sdf = stratify(df, stratum)
            if len(sdf) == 0:
                continue
            hp_sig_false = sdf["HPFineSig"] == 0
            tp_total = (sdf.y_true == 1).sum()
            fp_total = (sdf.y_true == 0).sum()
            tp_flagged = ((sdf.y_true == 1) & hp_sig_false).sum()
            fp_flagged = ((sdf.y_true == 0) & hp_sig_false).sum()
            tp_loss = tp_flagged / max(tp_total, 1) * 100
            fp_rem = fp_flagged / max(fp_total, 1) * 100
            _p(f"  HPFineSig=false ({stratum}): TP_loss={tp_loss:.1f}%, FP_removal={fp_rem:.1f}%")
            res[f"hpfinesig_filter_{stratum}"] = {
                "tp_loss_pct": float(tp_loss), "fp_removal_pct": float(fp_rem),
                "tp_flagged": int(tp_flagged), "fp_flagged": int(fp_flagged),
                "tp_total": int(tp_total), "fp_total": int(fp_total),
            }

        # 2. AlleleSig direction analysis (LOH vs Non-LOH)
        for stratum in ["LOH", "Non-LOH"]:
            sdf = stratify(df, stratum)
            tp_rate = sdf[sdf.y_true == 1]["AlleleSig"].mean()
            fp_rate = sdf[sdf.y_true == 0]["AlleleSig"].mean()
            direction = "TP>FP" if tp_rate > fp_rate else "FP>TP" if fp_rate > tp_rate else "equal"
            _p(f"  AlleleSig ({stratum}): TP={tp_rate:.4f}, FP={fp_rate:.4f}, dir={direction}")
            res[f"allelesig_{stratum}"] = {"tp_rate": float(tp_rate), "fp_rate": float(fp_rate), "direction": direction}

        # 3. Joint HPFineSig × AlleleSig 2×2 contingency (LOH)
        loh = stratify(df, "LOH")
        if len(loh) > 0:
            for label in ["TP", "FP"]:
                sub = loh[loh.truth_label == label]
                ct = pd.crosstab(sub["HPFineSig"], sub["AlleleSig"])
                _p(f"  Joint HPFineSig×AlleleSig ({label}):\n{ct.to_string()}")
                res[f"joint_2x2_{label}"] = ct.to_dict()

            # Joint filter: HPFineSig=false AND AlleleSig=false
            joint_false = (loh["HPFineSig"] == 0) & (loh["AlleleSig"] == 0)
            tp_j = ((loh.y_true == 1) & joint_false).sum()
            fp_j = ((loh.y_true == 0) & joint_false).sum()
            tp_tot = (loh.y_true == 1).sum()
            fp_tot = (loh.y_true == 0).sum()
            _p(f"  Joint HPFineSig=F & AlleleSig=F (LOH): TP_loss={tp_j/max(tp_tot,1)*100:.1f}%, FP_rem={fp_j/max(fp_tot,1)*100:.1f}%")
            res["joint_filter_loh"] = {
                "tp_loss_pct": float(tp_j / max(tp_tot, 1) * 100),
                "fp_removal_pct": float(fp_j / max(fp_tot, 1) * 100),
            }

        # 4. Per LOH_Subtype HPFineSig delta
        for sub in ["LOH_Strong", "LOH_Weak", "LOH_Subclone", "LOH_Noise"]:
            sdf = df[df["LOH_Subtype"] == sub]
            if len(sdf) < 20:
                continue
            tp_r = sdf[sdf.y_true == 1]["HPFineSig"].mean()
            fp_r = sdf[sdf.y_true == 0]["HPFineSig"].mean()
            delta = tp_r - fp_r
            _p(f"  HPFineSig rate ({sub}): TP={tp_r:.4f}, FP={fp_r:.4f}, Δ={delta:+.4f}")
            res[f"hpfinesig_rate_{sub}"] = {"tp_rate": float(tp_r), "fp_rate": float(fp_r), "delta": float(delta)}

        results[ver] = res

    # --- Q4 Figures ---
    setup_plot_style()

    # Q4-Fig1: HPFineSig/AlleleSig rate by stratum + joint heatmap
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for col_idx, ver in enumerate(dfs.keys()):
        df = dfs[ver]
        # HPFineSig rates
        ax = axes[0, col_idx]
        strata_plot = ["All", "LOH", "Non-LOH", "LOH_Strong", "LOH_Weak"]
        tp_rates = []
        fp_rates = []
        for s in strata_plot:
            sdf = stratify(df, s)
            tp_rates.append(sdf[sdf.y_true == 1]["HPFineSig"].mean())
            fp_rates.append(sdf[sdf.y_true == 0]["HPFineSig"].mean())
        x = np.arange(len(strata_plot))
        ax.bar(x - 0.15, tp_rates, 0.3, color=COLOR_TP, alpha=0.8, label="TP")
        ax.bar(x + 0.15, fp_rates, 0.3, color=COLOR_FP, alpha=0.8, label="FP")
        ax.set_xticks(x)
        ax.set_xticklabels(strata_plot, rotation=30, fontsize=8)
        ax.set_ylabel("HPFineSig Rate")
        ax.set_title(f"{ver} — HPFineSig")
        ax.legend()

        # AlleleSig rates
        ax = axes[1, col_idx]
        tp_rates = []
        fp_rates = []
        for s in strata_plot:
            sdf = stratify(df, s)
            tp_rates.append(sdf[sdf.y_true == 1]["AlleleSig"].mean())
            fp_rates.append(sdf[sdf.y_true == 0]["AlleleSig"].mean())
        ax.bar(x - 0.15, tp_rates, 0.3, color=COLOR_TP, alpha=0.8, label="TP")
        ax.bar(x + 0.15, fp_rates, 0.3, color=COLOR_FP, alpha=0.8, label="FP")
        ax.set_xticks(x)
        ax.set_xticklabels(strata_plot, rotation=30, fontsize=8)
        ax.set_ylabel("AlleleSig Rate")
        ax.set_title(f"{ver} — AlleleSig")
        ax.legend()

    fig.suptitle("Q4: Boolean Significance Features by Stratum", fontsize=14)
    fig.tight_layout()
    save_figure(fig, q4_dir / "Q4_fig1_boolean_rates.png")

    # Q4-Fig2: F1 impact analysis (filter effectiveness)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for idx, ver in enumerate(dfs.keys()):
        ax = axes[idx]
        filters = []
        for key in ["hpfinesig_filter_LOH", "hpfinesig_filter_LOH_Strong", "hpfinesig_filter_LOH_Weak"]:
            if key in results[ver]:
                r = results[ver][key]
                filters.append({"name": key.replace("hpfinesig_filter_", "HPFineSig=F\n"),
                                "tp_loss": r["tp_loss_pct"], "fp_removal": r["fp_removal_pct"]})
        if "joint_filter_loh" in results[ver]:
            r = results[ver]["joint_filter_loh"]
            filters.append({"name": "Joint\n(HPFine+Allele)", "tp_loss": r["tp_loss_pct"], "fp_removal": r["fp_removal_pct"]})
        if filters:
            fdf = pd.DataFrame(filters)
            x = np.arange(len(fdf))
            ax.bar(x - 0.15, fdf["tp_loss"], 0.3, color=COLOR_TP, alpha=0.8, label="TP Loss %")
            ax.bar(x + 0.15, fdf["fp_removal"], 0.3, color=COLOR_FP, alpha=0.8, label="FP Removal %")
            ax.axhline(y=2.0, color="gray", linestyle="--", alpha=0.5, label="2% TP Loss Limit")
            ax.set_xticks(x)
            ax.set_xticklabels(fdf["name"], fontsize=8)
            ax.set_ylabel("Percentage (%)")
            ax.set_title(f"{ver}")
            ax.legend(fontsize=8)
    fig.suptitle("Q4: Boolean Filter Performance (HPFineSig=false)", fontsize=13)
    fig.tight_layout()
    save_figure(fig, q4_dir / "Q4_fig2_filter_performance.png")

    with open(q4_dir / "q4_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    _p("Q4 complete.")
    return results


# ===========================================================================
# Q5: AlleleDelta LOH Inner/Outer Discrimination
# ===========================================================================
def run_q5(dfs: dict[str, pd.DataFrame]):
    _p("\n" + "="*60)
    _p("Q5: AlleleDelta LOH Inner/Outer Discrimination")
    _p("="*60)
    q5_dir = ensure_dir(OUT_ROOT / "Q5")

    results = {}
    for ver, df in dfs.items():
        _p(f"\n--- {ver} ---")
        res = {"version": ver}

        # 1. Zero-inflation analysis
        for stratum in ["LOH", "Non-LOH", "All"]:
            sdf = stratify(df, stratum)
            for label in ["TP", "FP"]:
                vals = sdf[sdf.truth_label == label]["AlleleDelta"]
                zero_pct = (vals == 0).sum() / max(len(vals), 1) * 100
                near_zero_pct = (vals < 0.001).sum() / max(len(vals), 1) * 100
                _p(f"  AlleleDelta zero% ({stratum}, {label}): exact_zero={zero_pct:.1f}%, <0.001={near_zero_pct:.1f}%")
                res[f"zero_pct_{stratum}_{label}"] = float(zero_pct)
                res[f"near_zero_pct_{stratum}_{label}"] = float(near_zero_pct)

        # 2. AlleleDelta stats by stratum
        for stratum in ["LOH", "Non-LOH"]:
            sdf = stratify(df, stratum)
            for label in ["TP", "FP"]:
                vals = sdf[sdf.truth_label == label]["AlleleDelta"]
                _p(f"  AlleleDelta ({stratum}, {label}): mean={vals.mean():.4f}, median={vals.median():.4f}, p75={vals.quantile(0.75):.4f}")

        # 3. AUC by stratum
        for stratum in ["LOH", "Non-LOH", "LOH_Strong", "LOH_Weak", "All"]:
            sdf = stratify(df, stratum)
            auc_pt, auc_lo, auc_hi = bootstrap_auc(sdf.y_true.values, sdf["AlleleDelta"].values)
            _p(f"  AlleleDelta AUC ({stratum}): {auc_pt:.4f} {format_ci(auc_lo, auc_hi)}")
            res[f"auc_{stratum}"] = {"auc": auc_pt, "ci_lo": auc_lo, "ci_hi": auc_hi}

        # 4. Minor HP count stratified AUC
        # min(HP1FamilyN, HP2FamilyN) as proxy for minor haplotype reads
        minor_hp = np.minimum(df["HP1FamilyN"].values, df["HP2FamilyN"].values)
        df_copy = df.copy()
        df_copy["minor_hp"] = minor_hp
        bins = [(0, 2), (3, 5), (6, 10), (11, 20), (21, 50), (51, 999)]
        bin_results = []
        for lo, hi in bins:
            mask = (df_copy["minor_hp"] >= lo) & (df_copy["minor_hp"] <= hi)
            sdf = df_copy[mask]
            if len(sdf) < 50 or len(sdf[sdf.y_true == 1]) < 10:
                continue
            auc_v = compute_auc(sdf.y_true.values, sdf["AlleleDelta"].values)
            bin_results.append({"bin": f"{lo}-{hi}", "n": len(sdf), "auc": float(auc_v),
                                "n_tp": int((sdf.y_true == 1).sum()), "n_fp": int((sdf.y_true == 0).sum())})
            _p(f"  AlleleDelta AUC (minor_hp={lo}-{hi}, n={len(sdf):,}): {auc_v:.4f}")
        res["minor_hp_auc"] = bin_results

        results[ver] = res

    # --- Q5 Figures ---
    setup_plot_style()

    # Q5-Fig1: AlleleDelta distribution 2×2 (LOH/Non-LOH × TP/FP)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for col_idx, ver in enumerate(dfs.keys()):
        df = dfs[ver]
        for row_idx, stratum in enumerate(["LOH", "Non-LOH"]):
            ax = axes[row_idx, col_idx]
            sdf = stratify(df, stratum)
            for label, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
                vals = sdf[sdf.truth_label == label]["AlleleDelta"].clip(0, 0.5)
                # Show zero spike as annotation
                zero_pct = (vals == 0).sum() / len(vals) * 100
                nonzero = vals[vals > 0]
                if len(nonzero) > 10:
                    ax.hist(nonzero, bins=50, density=True, alpha=0.4, color=color,
                            label=f"{label} (zero={zero_pct:.0f}%)")
            auc_v = compute_auc(sdf.y_true.values, sdf["AlleleDelta"].values)
            ax.set_title(f"{ver} — {stratum} (AUC={auc_v:.3f})", fontsize=10)
            ax.set_xlabel("AlleleDelta (non-zero only)")
            ax.set_ylabel("Density")
            ax.legend(fontsize=7)
    fig.suptitle("Q5: AlleleDelta Distribution (LOH vs Non-LOH, TP vs FP)", fontsize=13)
    fig.tight_layout()
    save_figure(fig, q5_dir / "Q5_fig1_alleledelta_distribution.png")

    # Q5-Fig2: Minor HP count stratified AUC + scatter
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for idx, ver in enumerate(dfs.keys()):
        ax = axes[idx]
        bins_data = results[ver].get("minor_hp_auc", [])
        if bins_data:
            bdf = pd.DataFrame(bins_data)
            x = range(len(bdf))
            colors = ["#4CAF50" if a > 0.55 else "#FF9800" if a > 0.50 else "#F44336" for a in bdf["auc"]]
            ax.bar(x, bdf["auc"].values, color=colors, alpha=0.8)
            ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5)
            ax.set_xticks(x)
            ax.set_xticklabels(bdf["bin"], rotation=30, fontsize=8)
            ax.set_xlabel("Minor HP Count Bin")
            ax.set_ylabel("AUC")
            ax.set_title(f"{ver} — AlleleDelta AUC by Minor HP")
            # Add n labels
            for i, row in bdf.iterrows():
                ax.text(i, row["auc"] + 0.005, f"n={row['n']:,}", ha="center", fontsize=6)
    fig.suptitle("Q5: AlleleDelta AUC Stratified by Minor Haplotype Read Count", fontsize=13)
    fig.tight_layout()
    save_figure(fig, q5_dir / "Q5_fig2_minor_hp_auc.png")

    with open(q5_dir / "q5_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    _p("Q5 complete.")
    return results


# ===========================================================================
# Main
# ===========================================================================
def main():
    _p("="*70)
    _p("TO ClairS-TO Feature Discrimination Deep Study — Q1-Q5")
    _p("="*70)

    # Check data availability
    available = {}
    for ver, paths in DATA_PATHS.items():
        tp_ok = paths["tp"].exists()
        fp_ok = paths["fp"].exists()
        _p(f"  {ver}: TP={'OK' if tp_ok else 'MISSING'}, FP={'OK' if fp_ok else 'MISSING'}")
        if tp_ok and fp_ok:
            available[ver] = True

    if not available:
        _p("ERROR: No data available. Exiting.")
        return

    # Load data
    dfs = {}
    for ver in available:
        _p(f"\nLoading {ver}...")
        df = load_version(ver)
        _p(f"  {ver}: {len(df):,} rows (TP={sum(df.y_true==1):,}, FP={sum(df.y_true==0):,})")
        _p(f"  LOH: {sum(df.Potential_LOH==1):,}, Non-LOH: {sum(df.Potential_LOH==0):,}")
        dfs[ver] = df

    # Phase 1: Quick comparison
    _p("\n" + "="*60)
    _p("Phase 1: Dual-Version Overview")
    _p("="*60)
    overview_rows = []
    for ver, df in dfs.items():
        row = {
            "version": ver,
            "total": len(df),
            "tp": int((df.y_true == 1).sum()),
            "fp": int((df.y_true == 0).sum()),
            "loh_pct": float((df.Potential_LOH == 1).sum() / len(df) * 100),
            "mean_num_reads": float(df["NumReads"].mean()),
            "mean_hp_ratio": float(df["HP_Ratio"].mean()),
            "mean_alleledelta": float(df["AlleleDelta"].mean()),
            "qs_auc": compute_auc(df.y_true.values, df["Quality_Score"].values),
        }
        overview_rows.append(row)
        _p(f"  {ver}: n={row['total']:,}, TP={row['tp']:,}, FP={row['fp']:,}, LOH={row['loh_pct']:.1f}%, QS_AUC={row['qs_auc']:.4f}")

    pd.DataFrame(overview_rows).to_csv(DATA_OUT / "phase1_overview.tsv", sep="\t", index=False)

    # Run Q1-Q5
    q1_res = run_q1(dfs)
    q2_res = run_q2(dfs)
    q3_res = run_q3(dfs)
    q4_res = run_q4(dfs)
    q5_res = run_q5(dfs)

    _p("\n" + "="*70)
    _p("Q1-Q5 COMPLETE")
    _p("="*70)
    _p(f"Output directory: {OUT_ROOT}")

    return {
        "overview": overview_rows,
        "q1": q1_res, "q2": q2_res, "q3": q3_res,
        "q4": q4_res, "q5": q5_res,
    }


if __name__ == "__main__":
    main()
