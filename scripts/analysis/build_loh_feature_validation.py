#!/usr/bin/env python3
"""
LOH 内 HPFineP / AlleleDelta TP-FP 区分力验证
===============================================
验证 O1 宣称: "TP 产生真正 HP fine 群组差异, FP 无系统性差异"
- HPFineP per-subtype distribution + AUC
- AlleleDelta per-subtype distribution + zero-inflation
- Per-chromosome stability (proxy for cross-sample robustness)
- Composition effect quantification
- Baseline vs v2b comparison

Output: output/clairsto_correction_analysis/feature_study/LOH_validation/
"""
import os, sys, json, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
plt.rcParams.update({"font.size": 9, "figure.dpi": 150})

# === Paths ===
BASE = "/big7_disk/liaoyoyo2001/longphase-to-mod/output"
VERSIONS = {
    "v2b": {
        "tp": f"{BASE}/ism_pononly_v2b_clairsto_tp/significance_summary.csv",
        "fp": f"{BASE}/ism_pononly_v2b_clairsto_fp/significance_summary.csv",
    },
    "baseline": {
        "tp": f"{BASE}/ism_baseline_clairsto_tp/significance_summary.csv",
        "fp": f"{BASE}/ism_baseline_clairsto_fp/significance_summary.csv",
    },
}
OUTDIR = "/big7_disk/liaoyoyo2001/InterSubMod/output/clairsto_correction_analysis/feature_study/LOH_validation"
os.makedirs(OUTDIR, exist_ok=True)

LOH_SUBTYPES = ["LOH_Strong", "LOH_Weak", "LOH_Noise"]
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]


def load_version(ver):
    """Load and merge TP/FP for a version."""
    tp = pd.read_csv(VERSIONS[ver]["tp"])
    fp = pd.read_csv(VERSIONS[ver]["fp"])
    tp["label"] = "TP"
    fp["label"] = "FP"
    df = pd.concat([tp, fp], ignore_index=True)
    # Parse chromosome
    df["chromosome"] = df["Chr"].astype(str)
    return df


def bootstrap_auc(y_true, y_score, n_boot=2000):
    """Bootstrap AUC with 95% CI."""
    n = len(y_true)
    if n < 10 or len(np.unique(y_true)) < 2:
        return np.nan, np.nan, np.nan
    aucs = []
    rng = np.random.RandomState(42)
    base_auc = roc_auc_score(y_true, y_score)
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        yt, ys = y_true[idx], y_score[idx]
        if len(np.unique(yt)) < 2:
            continue
        aucs.append(roc_auc_score(yt, ys))
    if len(aucs) < 100:
        return base_auc, np.nan, np.nan
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return base_auc, lo, hi


def compute_feature_stats(df_loh, feature, boolean=False):
    """Compute per-subtype stats for a feature."""
    results = []
    # Full LOH
    tp_vals = df_loh.loc[df_loh["label"] == "TP", feature].dropna()
    fp_vals = df_loh.loc[df_loh["label"] == "FP", feature].dropna()
    y = np.array([1]*len(tp_vals) + [0]*len(fp_vals))
    x = np.concatenate([tp_vals.values, fp_vals.values])

    if boolean:
        tp_rate = tp_vals.mean()
        fp_rate = fp_vals.mean()
        delta = tp_rate - fp_rate
        auc, lo, hi = bootstrap_auc(y, x.astype(float))
        results.append({
            "subtype": "All_LOH", "n_tp": len(tp_vals), "n_fp": len(fp_vals),
            "tp_rate": tp_rate, "fp_rate": fp_rate, "delta": delta,
            "auc": auc, "auc_lo": lo, "auc_hi": hi
        })
    else:
        auc, lo, hi = bootstrap_auc(y, x.astype(float))
        # Mann-Whitney U
        if len(tp_vals) > 0 and len(fp_vals) > 0:
            u_stat, u_p = stats.mannwhitneyu(tp_vals, fp_vals, alternative="two-sided")
        else:
            u_stat, u_p = np.nan, np.nan
        results.append({
            "subtype": "All_LOH", "n_tp": len(tp_vals), "n_fp": len(fp_vals),
            "tp_median": tp_vals.median(), "fp_median": fp_vals.median(),
            "tp_mean": tp_vals.mean(), "fp_mean": fp_vals.mean(),
            "auc": auc, "auc_lo": lo, "auc_hi": hi,
            "mann_whitney_p": u_p
        })

    # Per subtype
    for st in LOH_SUBTYPES:
        sub = df_loh[df_loh["LOH_Subtype"] == st]
        tp_v = sub.loc[sub["label"] == "TP", feature].dropna()
        fp_v = sub.loc[sub["label"] == "FP", feature].dropna()
        y_s = np.array([1]*len(tp_v) + [0]*len(fp_v))
        x_s = np.concatenate([tp_v.values, fp_v.values]) if len(tp_v) + len(fp_v) > 0 else np.array([])

        if boolean:
            tp_r = tp_v.mean() if len(tp_v) > 0 else np.nan
            fp_r = fp_v.mean() if len(fp_v) > 0 else np.nan
            d = tp_r - fp_r if not (np.isnan(tp_r) or np.isnan(fp_r)) else np.nan
            a, al, ah = bootstrap_auc(y_s, x_s.astype(float)) if len(x_s) > 10 else (np.nan, np.nan, np.nan)
            results.append({
                "subtype": st, "n_tp": len(tp_v), "n_fp": len(fp_v),
                "tp_rate": tp_r, "fp_rate": fp_r, "delta": d,
                "auc": a, "auc_lo": al, "auc_hi": ah
            })
        else:
            a, al, ah = bootstrap_auc(y_s, x_s.astype(float)) if len(x_s) > 10 else (np.nan, np.nan, np.nan)
            if len(tp_v) > 0 and len(fp_v) > 0:
                u_stat, u_p = stats.mannwhitneyu(tp_v, fp_v, alternative="two-sided")
            else:
                u_stat, u_p = np.nan, np.nan
            results.append({
                "subtype": st, "n_tp": len(tp_v), "n_fp": len(fp_v),
                "tp_median": tp_v.median() if len(tp_v) > 0 else np.nan,
                "fp_median": fp_v.median() if len(fp_v) > 0 else np.nan,
                "tp_mean": tp_v.mean() if len(tp_v) > 0 else np.nan,
                "fp_mean": fp_v.mean() if len(fp_v) > 0 else np.nan,
                "auc": a, "auc_lo": al, "auc_hi": ah,
                "mann_whitney_p": u_p
            })
    return pd.DataFrame(results)


def per_chromosome_auc(df_loh, feature, boolean=False):
    """Per-chromosome AUC for stability analysis."""
    rows = []
    for chrom in CHROMOSOMES:
        sub = df_loh[df_loh["chromosome"] == chrom]
        tp_v = sub.loc[sub["label"] == "TP", feature].dropna()
        fp_v = sub.loc[sub["label"] == "FP", feature].dropna()
        if len(tp_v) < 5 or len(fp_v) < 5:
            continue
        y = np.array([1]*len(tp_v) + [0]*len(fp_v))
        x = np.concatenate([tp_v.values, fp_v.values]).astype(float)
        try:
            auc = roc_auc_score(y, x)
        except:
            auc = np.nan
        rows.append({"chromosome": chrom, "n_tp": len(tp_v), "n_fp": len(fp_v), "auc": auc})
    return pd.DataFrame(rows)


def per_chromosome_subtype_auc(df_loh, feature, boolean=False):
    """Per-chromosome × per-subtype AUC matrix."""
    rows = []
    for chrom in CHROMOSOMES:
        for st in LOH_SUBTYPES:
            sub = df_loh[(df_loh["chromosome"] == chrom) & (df_loh["LOH_Subtype"] == st)]
            tp_v = sub.loc[sub["label"] == "TP", feature].dropna()
            fp_v = sub.loc[sub["label"] == "FP", feature].dropna()
            if len(tp_v) < 3 or len(fp_v) < 3:
                rows.append({"chromosome": chrom, "subtype": st, "auc": np.nan,
                             "n_tp": len(tp_v), "n_fp": len(fp_v)})
                continue
            y = np.array([1]*len(tp_v) + [0]*len(fp_v))
            x = np.concatenate([tp_v.values, fp_v.values]).astype(float)
            try:
                auc = roc_auc_score(y, x)
            except:
                auc = np.nan
            rows.append({"chromosome": chrom, "subtype": st, "auc": auc,
                         "n_tp": len(tp_v), "n_fp": len(fp_v)})
    return pd.DataFrame(rows)


# ============================================================
# Figure 1: HPFineP sig rate per LOH subtype (TP vs FP)
# 2x1 layout: v2b (left), baseline (right)
# ============================================================
def plot_fig1(data_dict):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for idx, (ver, df_loh) in enumerate(data_dict.items()):
        ax = axes[idx]
        # Compute HPFineSig boolean rate
        subtypes = ["All_LOH"] + LOH_SUBTYPES
        tp_rates, fp_rates, ns_tp, ns_fp = [], [], [], []
        for st in subtypes:
            sub = df_loh if st == "All_LOH" else df_loh[df_loh["LOH_Subtype"] == st]
            tp = sub[sub["label"] == "TP"]["HPFineSig"]
            fp = sub[sub["label"] == "FP"]["HPFineSig"]
            tp_rates.append(tp.mean() if len(tp) > 0 else 0)
            fp_rates.append(fp.mean() if len(fp) > 0 else 0)
            ns_tp.append(len(tp))
            ns_fp.append(len(fp))

        x = np.arange(len(subtypes))
        w = 0.35
        bars_tp = ax.bar(x - w/2, tp_rates, w, label="TP", color="#2196F3", alpha=0.8)
        bars_fp = ax.bar(x + w/2, fp_rates, w, label="FP", color="#F44336", alpha=0.8)

        # Add value labels
        for bar, rate, n in zip(bars_tp, tp_rates, ns_tp):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f"{rate:.1%}\n(n={n})", ha="center", va="bottom", fontsize=7)
        for bar, rate, n in zip(bars_fp, fp_rates, ns_fp):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f"{rate:.1%}\n(n={n})", ha="center", va="bottom", fontsize=7)

        ax.set_xticks(x)
        ax.set_xticklabels(subtypes, rotation=15)
        ax.set_ylabel("HPFineSig Rate")
        ax.set_title(f"{ver.upper()} — HPFineSig Rate by LOH Subtype")
        ax.legend()
        ax.set_ylim(0, 1.15)
        ax.axhline(0.5, color="gray", ls="--", lw=0.8, alpha=0.5)

    fig.suptitle("Fig 1: HPFineSig Rate — TP vs FP per LOH Subtype\n"
                 "O1 claim: 'FP lack systematic HP fine group differences'",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig1_hpfinesig_rate_by_subtype.png", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] Fig 1 saved")


# ============================================================
# Figure 2: HPFineP continuous distribution (violin + box)
# 2x2: rows = v2b/baseline, cols = LOH_Strong / LOH_Noise
# ============================================================
def plot_fig2(data_dict):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for row_idx, (ver, df_loh) in enumerate(data_dict.items()):
        for col_idx, st in enumerate(LOH_SUBTYPES):
            ax = axes[row_idx, col_idx]
            sub = df_loh[df_loh["LOH_Subtype"] == st]
            tp_p = sub.loc[sub["label"] == "TP", "HPFineP"].dropna()
            fp_p = sub.loc[sub["label"] == "FP", "HPFineP"].dropna()

            # Violin plot
            data_list = [tp_p.values, fp_p.values]
            labels = [f"TP (n={len(tp_p)})", f"FP (n={len(fp_p)})"]
            colors = ["#2196F3", "#F44336"]

            parts = ax.violinplot(data_list, positions=[0, 1], showmedians=True, showextrema=False)
            for pc, c in zip(parts["bodies"], colors):
                pc.set_facecolor(c)
                pc.set_alpha(0.5)
            parts["cmedians"].set_color("black")

            # Add box overlay
            bp = ax.boxplot(data_list, positions=[0, 1], widths=0.15,
                           patch_artist=True, showfliers=False)
            for patch, c in zip(bp["boxes"], colors):
                patch.set_facecolor(c)
                patch.set_alpha(0.3)

            # AUC annotation
            if len(tp_p) >= 5 and len(fp_p) >= 5:
                y = np.array([1]*len(tp_p) + [0]*len(fp_p))
                # For HPFineP, lower P = more significant = more TP-like
                # But we compute AUC as TP having LOWER values
                scores = np.concatenate([tp_p.values, fp_p.values])
                auc = roc_auc_score(y, -scores)  # negate: lower P → higher score for TP
                ax.text(0.5, 0.95, f"AUC={auc:.3f}\n(lower P → TP)",
                       transform=ax.transAxes, ha="center", va="top",
                       fontsize=8, bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))

            # Sig rate annotation
            tp_sig = (tp_p < 0.05).mean() if len(tp_p) > 0 else 0
            fp_sig = (fp_p < 0.05).mean() if len(fp_p) > 0 else 0
            ax.text(0.02, 0.02, f"TP sig={tp_sig:.1%}, FP sig={fp_sig:.1%}, Δ={tp_sig-fp_sig:+.1%}",
                   transform=ax.transAxes, fontsize=7, va="bottom",
                   bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

            ax.set_xticks([0, 1])
            ax.set_xticklabels(labels, fontsize=8)
            ax.set_ylabel("HPFineP (p-value)")
            ax.set_title(f"{ver} — {st}")
            ax.set_yscale("log")
            ax.set_ylim(1e-8, 2)
            ax.axhline(0.05, color="red", ls="--", lw=0.8, alpha=0.6, label="p=0.05")

    fig.suptitle("Fig 2: HPFineP Distribution per LOH Subtype — TP vs FP\n"
                 "Key question: Do FP also have low HPFineP (significant HP fine differences)?",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig2_hpfinep_distribution_by_subtype.png", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] Fig 2 saved")


# ============================================================
# Figure 3: AlleleDelta distribution + zero-inflation
# 2x3: rows = v2b/baseline, cols = LOH subtypes
# ============================================================
def plot_fig3(data_dict):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for row_idx, (ver, df_loh) in enumerate(data_dict.items()):
        for col_idx, st in enumerate(LOH_SUBTYPES):
            ax = axes[row_idx, col_idx]
            sub = df_loh[df_loh["LOH_Subtype"] == st]
            tp_d = sub.loc[sub["label"] == "TP", "AlleleDelta"].dropna()
            fp_d = sub.loc[sub["label"] == "FP", "AlleleDelta"].dropna()

            # Zero-inflation stats
            tp_zero = (tp_d == 0).mean() if len(tp_d) > 0 else np.nan
            fp_zero = (fp_d == 0).mean() if len(fp_d) > 0 else np.nan

            # Histogram (log scale on y)
            bins = np.linspace(0, 0.3, 40)
            ax.hist(tp_d.clip(0, 0.3), bins=bins, alpha=0.5, color="#2196F3",
                    label=f"TP (n={len(tp_d)})", density=True)
            ax.hist(fp_d.clip(0, 0.3), bins=bins, alpha=0.5, color="#F44336",
                    label=f"FP (n={len(fp_d)})", density=True)

            # Zero spike annotation
            ax.text(0.98, 0.95,
                    f"Zero%: TP={tp_zero:.1%}, FP={fp_zero:.1%}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=8,
                    bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

            # AUC
            if len(tp_d) >= 5 and len(fp_d) >= 5:
                y = np.array([1]*len(tp_d) + [0]*len(fp_d))
                scores = np.concatenate([tp_d.values, fp_d.values])
                auc = roc_auc_score(y, scores)
                ax.text(0.98, 0.78, f"AUC={auc:.3f}",
                       transform=ax.transAxes, ha="right", va="top", fontsize=9,
                       bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))

            ax.set_xlabel("AlleleDelta")
            ax.set_ylabel("Density")
            ax.set_title(f"{ver} — {st}")
            ax.legend(fontsize=7)

    fig.suptitle("Fig 3: AlleleDelta Distribution per LOH Subtype — TP vs FP\n"
                 "O1 claim: 'TP have nonzero AlleleDelta (somatic effect), FP have zero (germline no ASM)'",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig3_alleledelta_distribution_by_subtype.png", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] Fig 3 saved")


# ============================================================
# Figure 4: Per-chromosome AUC stability (HPFineP and AlleleDelta)
# 2x2: rows = feature, cols = v2b/baseline
# ============================================================
def plot_fig4(chr_auc_dict):
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    for col_idx, ver in enumerate(["v2b", "baseline"]):
        for row_idx, feat in enumerate(["HPFineSig", "AlleleDelta"]):
            ax = axes[row_idx, col_idx]
            df_chr = chr_auc_dict[(ver, feat)]
            if df_chr.empty:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
                continue

            colors = ["#4CAF50" if a > 0.5 else "#F44336" for a in df_chr["auc"]]
            bars = ax.bar(range(len(df_chr)), df_chr["auc"], color=colors, alpha=0.7)

            # Add chromosome labels
            ax.set_xticks(range(len(df_chr)))
            ax.set_xticklabels(df_chr["chromosome"].str.replace("chr", ""), rotation=45, fontsize=7)

            # Reference lines
            ax.axhline(0.5, color="red", ls="--", lw=1, label="Random (0.50)")
            mean_auc = df_chr["auc"].mean()
            ax.axhline(mean_auc, color="blue", ls=":", lw=1, label=f"Mean={mean_auc:.3f}")

            # Stability metrics
            std_auc = df_chr["auc"].std()
            n_above = (df_chr["auc"] > 0.5).sum()
            n_total = len(df_chr)
            ax.text(0.02, 0.95,
                    f"Mean={mean_auc:.3f}, SD={std_auc:.3f}\n"
                    f">{0.5}: {n_above}/{n_total} chr",
                    transform=ax.transAxes, va="top", fontsize=8,
                    bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

            ax.set_ylabel("AUC")
            ax.set_title(f"{ver} — {feat} — Per-chromosome AUC (LOH only)")
            ax.legend(fontsize=7)
            ax.set_ylim(0.3, 0.8)

    fig.suptitle("Fig 4: Per-Chromosome AUC Stability (LOH only)\n"
                 "Proxy for cross-sample robustness — consistent direction across chromosomes?",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig4_per_chromosome_auc_stability.png", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] Fig 4 saved")


# ============================================================
# Figure 5: Per-chromosome × per-subtype AUC heatmap
# 2x2: rows = HPFineSig/AlleleDelta, cols = v2b/baseline
# ============================================================
def plot_fig5(chr_sub_dict):
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    for col_idx, ver in enumerate(["v2b", "baseline"]):
        for row_idx, feat in enumerate(["HPFineSig", "AlleleDelta"]):
            ax = axes[row_idx, col_idx]
            df = chr_sub_dict[(ver, feat)]
            if df.empty:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
                continue

            # Pivot to matrix
            pivot = df.pivot_table(index="chromosome", columns="subtype", values="auc")
            # Reorder chromosomes
            chr_order = [f"chr{i}" for i in range(1, 23)]
            pivot = pivot.reindex([c for c in chr_order if c in pivot.index])

            # Heatmap
            im = ax.imshow(pivot.values, cmap="RdYlGn", vmin=0.3, vmax=0.7, aspect="auto")
            ax.set_xticks(range(len(pivot.columns)))
            ax.set_xticklabels(pivot.columns, fontsize=8)
            ax.set_yticks(range(len(pivot.index)))
            ax.set_yticklabels(pivot.index.str.replace("chr", ""), fontsize=7)

            # Annotate cells
            for i in range(len(pivot.index)):
                for j in range(len(pivot.columns)):
                    val = pivot.values[i, j]
                    if not np.isnan(val):
                        color = "white" if val < 0.4 or val > 0.6 else "black"
                        ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                                fontsize=6, color=color)

            ax.set_title(f"{ver} — {feat}")
            plt.colorbar(im, ax=ax, shrink=0.8, label="AUC")

    fig.suptitle("Fig 5: Per-Chromosome × Per-Subtype AUC Heatmap\n"
                 "Green (>0.5) = TP>FP direction, Red (<0.5) = reversed",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig5_chromosome_subtype_auc_heatmap.png", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] Fig 5 saved")


# ============================================================
# Figure 6: Composition effect quantification
# Bar chart: LOH-wide AUC vs within-subtype AUC
# ============================================================
def plot_fig6(stats_dict):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for idx, feat in enumerate(["HPFineSig", "AlleleDelta"]):
        ax = axes[idx]
        x_pos = np.arange(len(["v2b", "baseline"]))
        width = 0.15

        for ver_idx, ver in enumerate(["v2b", "baseline"]):
            st_df = stats_dict[(ver, feat)]
            all_auc = st_df.loc[st_df["subtype"] == "All_LOH", "auc"].values[0]
            sub_aucs = st_df.loc[st_df["subtype"].isin(LOH_SUBTYPES)]

            offset = ver_idx * 0.5
            # All LOH bar
            ax.bar(offset, all_auc, width, color="#1565C0", alpha=0.8, label="All LOH" if ver_idx == 0 else "")
            ax.text(offset, all_auc + 0.005, f"{all_auc:.3f}", ha="center", fontsize=8, fontweight="bold")

            # Per-subtype bars
            colors_st = {"LOH_Strong": "#4CAF50", "LOH_Weak": "#FFC107", "LOH_Noise": "#9E9E9E"}
            for s_idx, st in enumerate(LOH_SUBTYPES):
                row = sub_aucs[sub_aucs["subtype"] == st]
                if row.empty:
                    continue
                a = row["auc"].values[0]
                pos = offset + (s_idx + 1) * width
                ax.bar(pos, a, width, color=colors_st[st], alpha=0.8,
                       label=st if ver_idx == 0 else "")
                ax.text(pos, a + 0.005, f"{a:.3f}", ha="center", fontsize=7)

            # Composition effect arrow
            weighted_within = sub_aucs["auc"].mean()
            if not np.isnan(weighted_within) and not np.isnan(all_auc):
                comp_pct = (all_auc - weighted_within) / (all_auc - 0.5) * 100 if all_auc != 0.5 else 0
                ax.annotate(f"Composition\n{comp_pct:.0f}%",
                           xy=(offset, weighted_within),
                           xytext=(offset - 0.15, weighted_within + 0.04),
                           fontsize=7, color="red",
                           arrowprops=dict(arrowstyle="->", color="red", lw=1.5))

            ax.text(offset + width * 1.5, -0.03, ver.upper(), ha="center",
                    fontsize=9, fontweight="bold", transform=ax.get_xaxis_transform())

        ax.axhline(0.5, color="red", ls="--", lw=1, alpha=0.5)
        ax.set_ylabel("AUC")
        ax.set_title(f"{feat} — LOH-wide vs Within-Subtype AUC")
        ax.legend(fontsize=7, loc="upper right")
        ax.set_ylim(0.4, 0.75)
        ax.set_xticks([])

    fig.suptitle("Fig 6: Composition Effect Quantification\n"
                 "Gap between All_LOH and per-subtype AUC = composition effect",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig6_composition_effect_quantification.png", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] Fig 6 saved")


# ============================================================
# Figure 7: Minor HP reads analysis for AlleleDelta
# Show AlleleDelta vs min(HP1FamilyN, HP2FamilyN)
# ============================================================
def plot_fig7(data_dict):
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    for row_idx, (ver, df_loh) in enumerate(data_dict.items()):
        # Left: scatter AlleleDelta vs minor HP count
        ax = axes[row_idx, 0]
        df_loh = df_loh.copy()
        df_loh["minor_hp"] = df_loh[["HP1FamilyN", "HP2FamilyN"]].min(axis=1)

        tp = df_loh[df_loh["label"] == "TP"]
        fp = df_loh[df_loh["label"] == "FP"]

        ax.scatter(fp["minor_hp"], fp["AlleleDelta"], alpha=0.1, s=5, c="#F44336", label="FP", rasterized=True)
        ax.scatter(tp["minor_hp"], tp["AlleleDelta"], alpha=0.1, s=5, c="#2196F3", label="TP", rasterized=True)

        ax.set_xlabel("Minor HP Read Count (min(HP1, HP2))")
        ax.set_ylabel("AlleleDelta")
        ax.set_title(f"{ver} — AlleleDelta vs Minor HP Count (LOH)")
        ax.legend(fontsize=8, markerscale=5)
        ax.set_xlim(-1, 30)
        ax.set_ylim(-0.01, 0.3)

        # Right: AUC by minor HP bin
        ax2 = axes[row_idx, 1]
        bins = [0, 1, 3, 5, 10, 20, 999]
        bin_labels = ["0", "1-2", "3-4", "5-9", "10-19", "20+"]
        df_loh["minor_hp_bin"] = pd.cut(df_loh["minor_hp"], bins=bins, labels=bin_labels, right=False)

        aucs, ns = [], []
        for bl in bin_labels:
            sub = df_loh[df_loh["minor_hp_bin"] == bl]
            tp_v = sub.loc[sub["label"] == "TP", "AlleleDelta"].dropna()
            fp_v = sub.loc[sub["label"] == "FP", "AlleleDelta"].dropna()
            if len(tp_v) >= 5 and len(fp_v) >= 5:
                y = np.array([1]*len(tp_v) + [0]*len(fp_v))
                x = np.concatenate([tp_v.values, fp_v.values])
                auc = roc_auc_score(y, x)
                aucs.append(auc)
                ns.append(len(tp_v) + len(fp_v))
            else:
                aucs.append(np.nan)
                ns.append(len(tp_v) + len(fp_v))

        colors = ["#4CAF50" if (a and a > 0.5) else "#F44336" for a in aucs]
        bars = ax2.bar(range(len(bin_labels)), aucs, color=colors, alpha=0.7)
        for i, (a, n) in enumerate(zip(aucs, ns)):
            if not np.isnan(a):
                ax2.text(i, a + 0.005, f"{a:.3f}\n(n={n})", ha="center", fontsize=7)

        ax2.set_xticks(range(len(bin_labels)))
        ax2.set_xticklabels(bin_labels)
        ax2.set_xlabel("Minor HP Read Count Bin")
        ax2.set_ylabel("AUC")
        ax2.axhline(0.5, color="red", ls="--", lw=1)
        ax2.set_title(f"{ver} — AlleleDelta AUC by Minor HP Count (LOH)")
        ax2.set_ylim(0.35, 0.7)

    fig.suptitle("Fig 7: AlleleDelta Mechanism — Minor HP Reads Drive Discrimination\n"
                 "AlleleDelta=0 when minor HP reads=0 (physical inability, not biological signal)",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig7_alleledelta_minor_hp_mechanism.png", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] Fig 7 saved")


# ============================================================
# Figure 8: Combined summary — feature reliability scorecard
# ============================================================
def plot_fig8(all_stats):
    """Summary scorecard for all features across versions and subtypes."""
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.axis("off")

    # Build summary table
    headers = ["Feature", "Version", "All_LOH AUC", "Strong AUC", "Weak AUC",
               "Noise AUC", "Composition%", "Chr Consistency", "Verdict"]
    rows = []
    for ver in ["v2b", "baseline"]:
        for feat in ["HPFineSig", "AlleleDelta"]:
            st_df = all_stats[(ver, feat)]
            all_auc = st_df.loc[st_df["subtype"] == "All_LOH", "auc"].values[0]
            sub_aucs = {}
            for s in LOH_SUBTYPES:
                r = st_df.loc[st_df["subtype"] == s, "auc"]
                sub_aucs[s] = r.values[0] if len(r) > 0 else np.nan

            # Composition effect
            mean_within = np.nanmean(list(sub_aucs.values()))
            comp = (all_auc - mean_within) / (all_auc - 0.5) * 100 if all_auc != 0.5 else 0

            # Chr consistency from per-chr data
            chr_key = (ver, feat)
            # Will be added later

            # Verdict
            if all_auc < 0.55 and mean_within < 0.55:
                verdict = "NO SIGNAL"
            elif comp > 70:
                verdict = "COMPOSITION ONLY"
            elif comp > 40:
                verdict = "PARTIAL (weak within-subtype)"
            else:
                verdict = "REAL (within-subtype signal)"

            rows.append([feat, ver, f"{all_auc:.3f}", f"{sub_aucs.get('LOH_Strong', np.nan):.3f}",
                        f"{sub_aucs.get('LOH_Weak', np.nan):.3f}", f"{sub_aucs.get('LOH_Noise', np.nan):.3f}",
                        f"{comp:.0f}%", "—", verdict])

    table = ax.table(cellText=rows, colLabels=headers, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.8)

    # Color coding
    for i, row in enumerate(rows):
        verdict = row[-1]
        if "NO SIGNAL" in verdict:
            color = "#FFCDD2"  # red
        elif "COMPOSITION" in verdict:
            color = "#FFF9C4"  # yellow
        elif "PARTIAL" in verdict:
            color = "#FFE0B2"  # orange
        else:
            color = "#C8E6C9"  # green
        for j in range(len(headers)):
            table[i+1, j].set_facecolor(color)

    # Header color
    for j in range(len(headers)):
        table[0, j].set_facecolor("#E3F2FD")
        table[0, j].set_text_props(fontweight="bold")

    fig.suptitle("Fig 8: Feature Reliability Scorecard — LOH TP/FP Discrimination\n"
                 "Single sample (HCC1395 100% purity) — requires multi-sample validation",
                 fontweight="bold", fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig8_feature_reliability_scorecard.png", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] Fig 8 saved")


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 60)
    print("LOH Feature Validation: HPFineP + AlleleDelta")
    print("=" * 60)

    # Load data
    data_loh = {}
    all_stats = {}
    chr_auc = {}
    chr_sub_auc = {}

    for ver in ["v2b", "baseline"]:
        print(f"\n--- Loading {ver} ---")
        df = load_version(ver)
        print(f"  Total: {len(df)} ({(df['label']=='TP').sum()} TP, {(df['label']=='FP').sum()} FP)")

        df_loh = df[df["LOH_Subtype"].isin(LOH_SUBTYPES)].copy()
        print(f"  LOH: {len(df_loh)} ({(df_loh['label']=='TP').sum()} TP, {(df_loh['label']=='FP').sum()} FP)")
        for st in LOH_SUBTYPES:
            n = (df_loh["LOH_Subtype"] == st).sum()
            nt = ((df_loh["LOH_Subtype"] == st) & (df_loh["label"] == "TP")).sum()
            nf = ((df_loh["LOH_Subtype"] == st) & (df_loh["label"] == "FP")).sum()
            print(f"    {st}: {n} (TP={nt}, FP={nf})")

        data_loh[ver] = df_loh

        # Compute stats for both features
        for feat, is_bool in [("HPFineSig", True), ("AlleleDelta", False)]:
            print(f"\n  Computing {feat} stats...")
            st_df = compute_feature_stats(df_loh, feat, boolean=is_bool)
            all_stats[(ver, feat)] = st_df
            print(st_df.to_string(index=False))

            # Per-chromosome AUC
            chr_df = per_chromosome_auc(df_loh, feat, boolean=is_bool)
            chr_auc[(ver, feat)] = chr_df

            # Per-chromosome × per-subtype
            chr_sub_df = per_chromosome_subtype_auc(df_loh, feat, boolean=is_bool)
            chr_sub_auc[(ver, feat)] = chr_sub_df

    # Generate all figures
    print("\n\n=== Generating Figures ===")
    plot_fig1(data_loh)
    plot_fig2(data_loh)
    plot_fig3(data_loh)
    plot_fig4(chr_auc)
    plot_fig5(chr_sub_auc)
    plot_fig6(all_stats)
    plot_fig7(data_loh)
    plot_fig8(all_stats)

    # Save all stats to TSV
    print("\n=== Saving Stats ===")
    for (ver, feat), df in all_stats.items():
        path = f"{OUTDIR}/stats_{ver}_{feat}.tsv"
        df.to_csv(path, sep="\t", index=False)
        print(f"  [OK] {path}")

    for (ver, feat), df in chr_auc.items():
        path = f"{OUTDIR}/chr_auc_{ver}_{feat}.tsv"
        df.to_csv(path, sep="\t", index=False)
        print(f"  [OK] {path}")

    # Summary JSON
    summary = {
        "analysis": "LOH Feature Validation",
        "date": "2026-04-05",
        "sample": "HCC1395 (100% purity)",
        "versions": ["v2b (PON-only)", "baseline (longphase-to)"],
        "features_tested": ["HPFineSig/HPFineP", "AlleleDelta"],
        "key_findings": {},
        "single_sample_limitation": True,
        "cross_sample_validated": False,
    }

    for ver in ["v2b", "baseline"]:
        ver_findings = {}
        for feat in ["HPFineSig", "AlleleDelta"]:
            st_df = all_stats[(ver, feat)]
            all_row = st_df[st_df["subtype"] == "All_LOH"].iloc[0]
            sub_rows = st_df[st_df["subtype"].isin(LOH_SUBTYPES)]
            mean_within = sub_rows["auc"].mean()
            comp = (all_row["auc"] - mean_within) / (all_row["auc"] - 0.5) * 100 if all_row["auc"] != 0.5 else 0

            # Chr stability
            chr_df = chr_auc[(ver, feat)]
            chr_above = (chr_df["auc"] > 0.5).sum() if len(chr_df) > 0 else 0
            chr_total = len(chr_df)
            chr_mean = chr_df["auc"].mean() if len(chr_df) > 0 else np.nan
            chr_std = chr_df["auc"].std() if len(chr_df) > 0 else np.nan

            ver_findings[feat] = {
                "all_loh_auc": round(float(all_row["auc"]), 4),
                "mean_within_subtype_auc": round(float(mean_within), 4),
                "composition_effect_pct": round(float(comp), 1),
                "chr_consistency": f"{chr_above}/{chr_total}",
                "chr_mean_auc": round(float(chr_mean), 4),
                "chr_std_auc": round(float(chr_std), 4),
            }
        summary["key_findings"][ver] = ver_findings

    with open(f"{OUTDIR}/validation_summary.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"\n  [OK] {OUTDIR}/validation_summary.json")

    # Print final verdict
    print("\n" + "=" * 60)
    print("FINAL VERDICT")
    print("=" * 60)
    for ver in ["v2b", "baseline"]:
        print(f"\n--- {ver} ---")
        for feat in ["HPFineSig", "AlleleDelta"]:
            f = summary["key_findings"][ver][feat]
            comp = f["composition_effect_pct"]
            if comp > 70:
                v = "COMPOSITION ONLY"
            elif comp > 40:
                v = "PARTIAL (weak within-subtype)"
            else:
                v = "REAL"
            print(f"  {feat}: LOH AUC={f['all_loh_auc']}, Within={f['mean_within_subtype_auc']}, "
                  f"Comp={comp:.0f}%, Chr={f['chr_consistency']}, → {v}")

    print(f"\n⚠ ALL RESULTS FROM SINGLE SAMPLE (HCC1395 100% purity)")
    print(f"  No other TO sample ISM data available for cross-validation")
    print(f"  Per-chromosome stability is the ONLY proxy for robustness")


if __name__ == "__main__":
    main()
