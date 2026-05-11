#!/usr/bin/env python3
"""R1 cross-mode concordance analysis.

Goal: 在 TO arm 和 paired arm 共同 region 上，檢查：
  1. Normal BAM 結果的 sanity（coverage, SampleASM 分布合理？）
  2. TO SampleASM_Delta 與 paired Fisher_Frac_Sig 是否捕捉同一訊號？
  3. 兩模式 TP/FP 分類一致性
  4. Cross-mode feature correlation matrix

Outputs:
  figures/fig6_normal_bam_sanity.png
  figures/fig7_crossmode_signal_scatter.png
  figures/fig8_crossmode_feature_correlation.png
  figures/fig9_crossmode_tp_agreement.png
  observations/step7_crossmode_concordance.md 資料補充
"""
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
PILOT = ROOT / "output/hcc1395_normal_pilot"
DATA = ROOT / "research/F_hpfinengroups_deepening/data"
FIG = ROOT / "research/F_hpfinengroups_deepening/figures"

sns.set_theme(style="whitegrid", context="talk")


def load_arm(mode):
    sig = pd.read_csv(PILOT / mode / "significance_summary.csv", low_memory=False)
    truth = pd.read_csv(DATA / f"r1_hcc1395_filter_snvs_{'to' if mode == 'TO' else mode}.tsv", sep="\t")
    sig["key"] = sig["Chr"].astype(str) + ":" + sig["Pos"].astype(str)
    truth["key"] = truth["Chr"].astype(str) + ":" + truth["Pos"].astype(str)
    merged = sig.merge(truth[["key", "truth_label"]], on="key", how="left")
    merged["is_tp"] = (merged["truth_label"].astype(str) == "TP").astype(int)
    return merged


def fig_normal_sanity(df_to, df_pa, out):
    """Sanity check: Normal BAM coverage, Sample ASM distributions."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))

    # Row 1: TO arm
    ax = axes[0, 0]
    ax.hist(df_to["NormalBaseline_Coverage"].dropna(), bins=40, color="#2c7fb8", edgecolor="k")
    ax.set_title(f"TO arm — NormalBaseline_Coverage\nmean={df_to['NormalBaseline_Coverage'].mean():.1f} median={df_to['NormalBaseline_Coverage'].median():.1f}")
    ax.set_xlabel("Normal coverage (read count)")
    ax.set_ylabel("region count")

    ax = axes[0, 1]
    ax.hist(df_to["NormalBaseline_Mean"].dropna(), bins=40, color="#2c7fb8", edgecolor="k")
    ax.set_title(f"TO arm — NormalBaseline_Mean\nmean={df_to['NormalBaseline_Mean'].mean():.3f}")
    ax.set_xlabel("Normal methylation baseline (mean)")

    ax = axes[0, 2]
    ax.hist(df_to["SampleASM_Delta"].dropna(), bins=40, color="#fdb462", edgecolor="k")
    ax.axvline(0, color="red", ls="--")
    ax.set_title(f"TO arm — SampleASM_Delta\nmean={df_to['SampleASM_Delta'].mean():.4f}")
    ax.set_xlabel("Sample ASM Δ (tumor − normal)")

    # Row 2: paired arm
    ax = axes[1, 0]
    ax.hist(df_pa["NormalBaseline_Coverage"].dropna(), bins=40, color="#2c7fb8", edgecolor="k")
    ax.set_title(f"paired arm — NormalBaseline_Coverage\nmean={df_pa['NormalBaseline_Coverage'].mean():.1f} median={df_pa['NormalBaseline_Coverage'].median():.1f}")
    ax.set_xlabel("Normal coverage (read count)")
    ax.set_ylabel("region count")

    ax = axes[1, 1]
    ax.hist(df_pa["NormalBaseline_Mean"].dropna(), bins=40, color="#2c7fb8", edgecolor="k")
    ax.set_title(f"paired arm — NormalBaseline_Mean\nmean={df_pa['NormalBaseline_Mean'].mean():.3f}")
    ax.set_xlabel("Normal methylation baseline (mean)")

    ax = axes[1, 2]
    ax.hist(df_pa["SampleASM_Delta"].dropna(), bins=40, color="#fdb462", edgecolor="k")
    ax.axvline(0, color="red", ls="--")
    ax.set_title(f"paired arm — SampleASM_Delta\nmean={df_pa['SampleASM_Delta'].mean():.4f}")
    ax.set_xlabel("Sample ASM Δ (tumor − normal)")

    plt.suptitle("HCC1395 Normal BAM sanity check (region-subset 1.8GB)", fontsize=16, y=1.01)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()

    # Return summary stats
    stats_rows = []
    for name, df in [("TO", df_to), ("paired", df_pa)]:
        for col in ["NormalBaseline_Coverage", "NormalBaseline_Mean", "SampleASM_Delta"]:
            s = df[col].dropna()
            stats_rows.append({
                "arm": name, "feature": col,
                "n": len(s), "mean": s.mean(), "median": s.median(),
                "std": s.std(), "min": s.min(), "max": s.max(),
                "n_zero": int((s == 0).sum()),
                "n_nan": int(df[col].isna().sum()),
            })
    return pd.DataFrame(stats_rows)


def match_regions(df_to, df_pa):
    """Inner join on Chr:Pos."""
    common = df_to.merge(
        df_pa, on="key", how="inner", suffixes=("_TO", "_paired")
    )
    return common


def fig_crossmode_scatter(common, out):
    """TO SampleASM_Delta vs paired Fisher_Frac_Sig on matched regions."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))

    # Panel 1: TO SampleASM vs paired Fisher
    ax = axes[0]
    sub = common.dropna(subset=["SampleASM_Delta_TO", "Fisher_Frac_Sig_paired"])
    for label, color in [(1, "#2c7fb8"), (0, "#d7191c")]:
        s = sub[sub["is_tp_TO"] == label]
        ax.scatter(s["SampleASM_Delta_TO"], s["Fisher_Frac_Sig_paired"],
                   alpha=0.5, s=20, c=color, label=f"{'TP' if label else 'FP'} (n={len(s)})")
    rho, p = stats.spearmanr(sub["SampleASM_Delta_TO"], sub["Fisher_Frac_Sig_paired"])
    ax.set_xlabel("TO arm — SampleASM_Delta")
    ax.set_ylabel("paired arm — Fisher_Frac_Sig")
    ax.set_title(f"Cross-mode signal scatter\nSpearman ρ={rho:.3f}, p={p:.2e}")
    ax.legend()

    # Panel 2: TO NormalBaseline_Coverage vs paired NormalBaseline_Coverage (should correlate, same normal BAM)
    ax = axes[1]
    sub = common.dropna(subset=["NormalBaseline_Coverage_TO", "NormalBaseline_Coverage_paired"])
    ax.scatter(sub["NormalBaseline_Coverage_TO"], sub["NormalBaseline_Coverage_paired"],
               alpha=0.4, s=15, c="gray")
    max_v = max(sub["NormalBaseline_Coverage_TO"].max(), sub["NormalBaseline_Coverage_paired"].max())
    ax.plot([0, max_v], [0, max_v], "r--", alpha=0.5, label="y=x")
    rho, _ = stats.spearmanr(sub["NormalBaseline_Coverage_TO"], sub["NormalBaseline_Coverage_paired"])
    ax.set_xlabel("TO arm — NormalBaseline_Coverage")
    ax.set_ylabel("paired arm — NormalBaseline_Coverage")
    ax.set_title(f"Normal coverage sanity (same Normal BAM)\nSpearman ρ={rho:.3f}")
    ax.legend()

    # Panel 3: TO is_tp vs paired is_tp agreement
    ax = axes[2]
    sub = common.dropna(subset=["is_tp_TO", "is_tp_paired"])
    # 2x2 contingency
    ct = pd.crosstab(sub["is_tp_TO"], sub["is_tp_paired"])
    sns.heatmap(ct, annot=True, fmt="d", cmap="Blues", ax=ax, cbar=False, linewidths=0.5, linecolor="white")
    ax.set_xlabel("paired truth_label (1=TP)")
    ax.set_ylabel("TO truth_label (1=TP)")
    ax.set_title(f"TO×paired truth agreement\n(common regions n={len(sub)})")

    plt.suptitle("Cross-mode concordance analysis (matched regions)", fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()


def fig_feature_correlation(common, out):
    """Cross-mode feature correlation heatmap."""
    to_feats = ["SampleASM_Delta_TO", "NormalBaseline_Coverage_TO", "HP_Residual_Delta_TO",
                "HPFineF_TO", "Tumor_HP_Delta_TO", "Epipoly_Delta_TO"]
    pa_feats = ["Fisher_Frac_Sig_paired", "Epipoly_Delta_paired", "Tumor_HP_Delta_paired",
                "HPFineF_paired", "SampleASM_Delta_paired", "NormalBaseline_Coverage_paired"]
    cols = [c for c in to_feats + pa_feats if c in common.columns]
    sub = common[cols].dropna()
    if len(sub) < 20:
        print(f"[!] Not enough matched data: {len(sub)}")
        return
    corr = sub.corr(method="spearman")

    fig, ax = plt.subplots(figsize=(14, 11))
    mask = np.zeros_like(corr, dtype=bool)
    # highlight cross-mode (TO × paired) quadrant
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="RdBu_r", center=0, vmin=-0.6, vmax=0.6,
                ax=ax, cbar_kws={"label": "Spearman ρ"}, linewidths=0.5, linecolor="white")
    ax.set_title(f"Cross-mode feature correlation (Spearman, n={len(sub)} matched regions)\n"
                 "Upper-left = TO×TO; Lower-right = paired×paired; Off-diagonal = TO×paired independence")
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()


def fig_tp_agreement(common, out):
    """Per-region TP/FP agreement + signal direction."""
    sub = common.dropna(subset=["is_tp_TO", "is_tp_paired"])
    sub["agreement"] = np.where(
        sub["is_tp_TO"] == sub["is_tp_paired"],
        np.where(sub["is_tp_TO"] == 1, "both_TP", "both_FP"),
        np.where(sub["is_tp_TO"] == 1, "TO_TP_paired_FP", "TO_FP_paired_TP")
    )

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Left: agreement pie
    counts = sub["agreement"].value_counts()
    colors = {"both_TP": "#2c7fb8", "both_FP": "#d7191c",
              "TO_TP_paired_FP": "#fdb462", "TO_FP_paired_TP": "#80b1d3"}
    pie_colors = [colors[k] for k in counts.index]
    axes[0].pie(counts.values, labels=[f"{k}\n(n={v})" for k, v in counts.items()],
                colors=pie_colors, autopct="%.1f%%", startangle=90)
    axes[0].set_title(f"TP/FP agreement across modes (n={len(sub)})")

    # Right: SampleASM_Delta_TO by agreement group
    ax = axes[1]
    for group, color in colors.items():
        vals = sub[sub["agreement"] == group]["SampleASM_Delta_TO"].dropna()
        if len(vals) >= 3:
            ax.scatter([group] * len(vals), vals, alpha=0.5, s=30, color=color)
            ax.scatter([group], [vals.median()], marker="_", s=800, color="black", zorder=10)
    ax.set_ylabel("TO arm — SampleASM_Delta")
    ax.set_xlabel("Agreement group")
    ax.set_title("SampleASM_Delta by agreement group")
    ax.tick_params(axis="x", rotation=20)
    ax.axhline(0, color="red", ls="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()

    return sub["agreement"].value_counts().to_dict()


def main():
    df_to = load_arm("TO")
    df_pa = load_arm("paired")
    print(f"TO rows: {len(df_to)}, paired rows: {len(df_pa)}")

    # Sanity figures
    print("\n--- Figure 6: Normal BAM sanity ---")
    stats_df = fig_normal_sanity(df_to, df_pa, FIG / "fig6_normal_bam_sanity.png")
    print(stats_df.to_string(index=False))
    stats_df.to_csv(FIG / "fig6_normal_bam_sanity.tsv", sep="\t", index=False)

    # Cross-mode match
    common = match_regions(df_to, df_pa)
    print(f"\n--- Matched common regions: {len(common)} ---")
    print(f"TO has {len(df_to)} regions; paired has {len(df_pa)} regions; intersection {len(common)}")

    print("\n--- Figure 7: Cross-mode signal scatter ---")
    fig_crossmode_scatter(common, FIG / "fig7_crossmode_signal_scatter.png")

    print("--- Figure 8: Feature correlation ---")
    fig_feature_correlation(common, FIG / "fig8_crossmode_feature_correlation.png")

    print("--- Figure 9: TP agreement ---")
    agg = fig_tp_agreement(common, FIG / "fig9_crossmode_tp_agreement.png")
    print(f"Agreement distribution: {agg}")

    # Independence tests
    print("\n--- Cross-mode independence ---")
    sub = common.dropna(subset=["SampleASM_Delta_TO", "Fisher_Frac_Sig_paired"])
    if len(sub) >= 20:
        rho, p = stats.spearmanr(sub["SampleASM_Delta_TO"], sub["Fisher_Frac_Sig_paired"])
        print(f"TO SampleASM_Delta vs paired Fisher_Frac_Sig: ρ={rho:.4f}, p={p:.3e}, n={len(sub)}")

    sub = common.dropna(subset=["NormalBaseline_Coverage_TO", "NormalBaseline_Coverage_paired"])
    if len(sub) >= 20:
        rho, p = stats.spearmanr(sub["NormalBaseline_Coverage_TO"], sub["NormalBaseline_Coverage_paired"])
        print(f"Normal coverage TO vs paired (same BAM, sanity): ρ={rho:.4f}, p={p:.3e}, n={len(sub)}")

    # Write supplementary concordance data
    agreement_tsv = FIG / "fig9_crossmode_tp_agreement.tsv"
    pd.DataFrame([agg]).T.rename(columns={0: "n"}).to_csv(agreement_tsv, sep="\t")
    print(f"\n[+] All figures: {FIG}/fig{{6,7,8,9}}*.png")


if __name__ == "__main__":
    main()
