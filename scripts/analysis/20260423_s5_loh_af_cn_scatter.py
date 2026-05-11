#!/usr/bin/env python3
"""
S5 圖：各樣本 × (LOH Inner / Outer) × TP/FP 散點圖
X: Coverage_Multiple (CN proxy)   Y: AF (caller)
AF 三分類背景 + 虛線分隔
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

SRC = "/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"
# HCC1395 的 phase1_new run LOH annotation 不完整，改用 stale master 的 archive 資料
# （NumReads/54 作 KDE-equivalent CovM；54 為 SEQC2 neutral median，HCC1395 pilot bx=53.0×）
SRC_HCC1395_ARCHIVE = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz"
HCC1395_KDE_BASELINE = 54.0  # SEQC2 neutral median

OUTDIR = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260423_s5_loh_af_cn_scatter"
os.makedirs(OUTDIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# Palette
C_TP = "#0072B2"
C_FP = "#D55E00"
C_EXTREME = "#F0F0F0"
C_NEARHALF = "#E0F2E5"
C_INTER = "#FFF0DC"
C_DASH = "#888888"

# 樣本順序（6 TO 樣本 · 按 LOH 佔比 / purity 排序）
SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H2009", "H1437"]

# 各樣本 KDE baseline (bx)：HCC1395 54 (SEQC2 neutral)，其他來自 merged Diploid_Coverage_Used
SAMPLE_KDE_BX = {
    "HCC1395": 54.0,
    "HCC1395_DORADO": 53.0,
    "HCC1937": 91.0,
    "HCC1954": 61.0,
    "H2009": 79.0,
    "H1437": 71.0,
}


def load_all_from_archive():
    """All 6 TO samples from stale master archive, unified schema.
    用 caller_af 作 Y (而非 merged 檔的 AF — 後者非真正 VAF)；
    用 NumReads / per-sample KDE bx 作 X (KDE-equivalent CovM)。
    """
    arch = pd.read_csv(SRC_HCC1395_ARCHIVE, sep="\t", low_memory=False)
    # 只取 TO
    arch_to = arch[arch["mode"] == "to"].copy()
    # 只留我們要的 6 樣本
    arch_to = arch_to[arch_to["sample"].isin(SAMPLE_KDE_BX.keys())].copy()
    # 映射 bx
    arch_to["bx"] = arch_to["sample"].map(SAMPLE_KDE_BX)
    # KDE-equivalent CovM
    arch_to["Coverage_Multiple"] = arch_to["NumReads"] / arch_to["bx"]
    arch_to["Diploid_Coverage_Used"] = arch_to["bx"]
    # 使用真正 caller AF
    arch_to["AF"] = arch_to["caller_af"]
    arch_to["loh_side"] = arch_to["Potential_LOH"].map({True: "Inner", False: "Outer"})
    arch_to["tp_label"] = (arch_to["truth_label"] == "TP").astype(int)
    arch_to["mode"] = "to_pileup"
    arch_to["kde_status"] = "archive_unified_kde_equiv"
    return arch_to[["sample", "mode", "AF", "Coverage_Multiple", "loh_side", "tp_label",
                    "Diploid_Coverage_Used", "kde_status", "LOH_Subtype"]]


def main():
    # 6 個 TO 樣本全部從 stale master archive 讀，統一使用 caller_af 與 per-sample KDE bx
    df = load_all_from_archive()

    df = df.dropna(subset=["AF", "Coverage_Multiple", "loh_side", "tp_label"])
    df["tpfp"] = df["tp_label"].apply(
        lambda x: "TP" if (x in (1, 1.0, True) or str(x).strip().lower() in ("1", "true", "tp"))
        else "FP")
    df["CovM_clip"] = df["Coverage_Multiple"].clip(0.05, 3.5)
    sample_baseline = df.groupby("sample")["Diploid_Coverage_Used"].median().to_dict()

    # Per-sample baseline TP rate（用於 panel title 的 Δpp 對照）
    sample_base_tp = df.groupby("sample")["tp_label"].mean().to_dict()

    # 先計算真實 TP/FP 數量 (subsample 前)
    real_counts = df.groupby(["sample", "loh_side", "tpfp"]).size().unstack(fill_value=0)
    real_counts = real_counts.reindex(columns=["TP", "FP"], fill_value=0)

    # 子取樣：每 (sample, loh_side) 固定最多 6000 點（統一密度）
    def sub(g):
        n = min(len(g), 6000)
        return g.sample(n=n, random_state=42) if len(g) > n else g
    df_plot = df.groupby(["sample", "loh_side"], group_keys=False).apply(sub)

    # caller FP 背景 (Outer Extreme AF 全體 FP rate, 用於標記 HCC1954/HCC1937 異常)
    outer_extreme_fp = {}
    for s in df["sample"].unique():
        sub_oe = df[(df["sample"] == s) & (df["loh_side"] == "Outer") &
                    ((df["AF"] < 0.1) | (df["AF"] > 0.9))]
        if len(sub_oe):
            fp_pct = (sub_oe["tpfp"] == "FP").sum() / len(sub_oe) * 100
            outer_extreme_fp[s] = fp_pct

    n_samples = len(SAMPLE_ORDER)
    fig, axes = plt.subplots(n_samples, 2, figsize=(12, 2.2 * n_samples), sharex=True, sharey=True)

    for row, sample in enumerate(SAMPLE_ORDER):
        for col, side in enumerate(["Inner", "Outer"]):
            ax = axes[row, col] if n_samples > 1 else axes[col]
            sub_df = df_plot[(df_plot["sample"] == sample) & (df_plot["loh_side"] == side)]
            # 真實數量
            try:
                n_tp = int(real_counts.loc[(sample, side), "TP"])
                n_fp = int(real_counts.loc[(sample, side), "FP"])
            except KeyError:
                n_tp = n_fp = 0
            n_total = n_tp + n_fp
            tp_rate = n_tp / n_total if n_total else 0

            # AF 背景分區
            ax.axhspan(0.0, 0.1, color=C_EXTREME, alpha=0.75, zorder=0)
            ax.axhspan(0.9, 1.0, color=C_EXTREME, alpha=0.75, zorder=0)
            ax.axhspan(0.4, 0.6, color=C_NEARHALF, alpha=0.75, zorder=0)
            # Intermediate = 剩餘兩塊（0.1-0.4 與 0.6-0.9）
            ax.axhspan(0.1, 0.4, color=C_INTER, alpha=0.4, zorder=0)
            ax.axhspan(0.6, 0.9, color=C_INTER, alpha=0.4, zorder=0)

            # AF 三分類虛線
            for y in (0.1, 0.4, 0.6, 0.9):
                ax.axhline(y, color=C_DASH, linestyle="--", linewidth=0.6, zorder=1)
            # CN tier 虛線 (SEQC2 empirical 0.65 / 0.99 / 1.33 / 1.82)
            for x in (0.65, 0.99, 1.33, 1.82):
                ax.axvline(x, color=C_DASH, linestyle=":", linewidth=0.5, alpha=0.6, zorder=1)

            # 畫 TP 先（底層、稍大），再 FP（上層、小點）— 讓 FP 可見但不蓋過 TP 密度
            fp = sub_df[sub_df["tpfp"] == "FP"]
            tp = sub_df[sub_df["tpfp"] == "TP"]
            ax.scatter(tp["CovM_clip"], tp["AF"], s=5, c=C_TP, alpha=0.4, zorder=2)
            ax.scatter(fp["CovM_clip"], fp["AF"], s=2.2, c=C_FP, alpha=0.6, zorder=3)

            ax.set_xlim(0, 3.5)
            ax.set_ylim(-0.02, 1.02)

            # 標題：加 baseline TP rate + Δ (相對 baseline 改善)，讓讀者看到 Inner/Outer 對 baseline 的真實增益
            tag = "LOH Inner" if side == "Inner" else "LOH Outer"
            bx = sample_baseline.get(sample, None)
            base_tp = sample_base_tp.get(sample, tp_rate)
            delta = (tp_rate - base_tp) * 100
            delta_color = "#2E7D32" if delta > 0 else "#C62828"
            bx_str = f"  bx={bx:.0f}×" if bx else ""
            ax.set_title(
                f"{sample} · {tag}  TP={tp_rate:.2f} (base={base_tp:.2f}, Δ{delta:+.1f}pp)  n={n_total:,}{bx_str}",
                fontsize=9, pad=3, color=delta_color)

            if row == n_samples - 1:
                ax.set_xlabel("Coverage_Multiple (KDE-calibrated CN proxy)  —  tier: 0.65 / 0.99 / 1.33 / 1.82")
            if col == 0:
                ax.set_ylabel("caller AF")

            # AF 分區標籤（僅第一列）
            if row == 0 and col == 1:
                ax.text(3.47, 0.05, "Extreme", ha="right", va="center", fontsize=7, color="#555", style="italic")
                ax.text(3.47, 0.50, "Near-half", ha="right", va="center", fontsize=7, color="#2E7D32", style="italic")
                ax.text(3.47, 0.25, "Inter.", ha="right", va="center", fontsize=7, color="#BF6000", style="italic")
                ax.text(3.47, 0.75, "Inter.", ha="right", va="center", fontsize=7, color="#BF6000", style="italic")
                ax.text(3.47, 0.95, "Extreme", ha="right", va="center", fontsize=7, color="#555", style="italic")

    # 全圖標題（全英文避免 glyph 問題）
    fig.suptitle("TP/FP scatter per sample × (LOH Inner / Outer)  —  "
                 "X: NumReads / KDE bx (per-sample CN proxy)   Y: caller_af (真 VAF)\n"
                 "AF tri-class (Extreme <0.1 或 >0.9 / Near-half 0.4–0.6 / Intermediate 其他) · Title shows Δpp vs sample baseline\n"
                 "All 6 TO samples unified source: stale master archive (complete LOH) · caller_af · per-sample KDE bx",
                 fontsize=10, y=0.998, fontweight="bold", color="#1E2A44")

    handles = [
        mpatches.Patch(color=C_EXTREME, label="Extreme (AF<0.1 or >0.9)"),
        mpatches.Patch(color=C_NEARHALF, label="Near-half (0.4-0.6)"),
        mpatches.Patch(color=C_INTER, alpha=0.4, label="Intermediate (0.1-0.4 or 0.6-0.9)"),
        mpatches.Patch(color=C_TP, label="TP (true positive)"),
        mpatches.Patch(color=C_FP, label="FP (false positive)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=9,
               bbox_to_anchor=(0.5, 0.0), frameon=False)

    plt.tight_layout(rect=[0, 0.02, 1, 0.985])
    out = os.path.join(OUTDIR, "fig_s5_loh_inner_outer_af_cn_per_sample.png")
    plt.savefig(out, bbox_inches="tight", dpi=130)
    plt.close()
    print(f"[OK] {out}")
    print(f"     size: {os.path.getsize(out)/1024:.1f} KB")


if __name__ == "__main__":
    main()
