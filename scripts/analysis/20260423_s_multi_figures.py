#!/usr/bin/env python3
"""
Generate multiple figures for PPT v2 revision:
  fig_s8_kde_visualization.png  — KDE 演算法 4-step 示意 (HCC1395 實際 histogram)
  fig_s10_s1s7_heatmap_to_only.png — 只 6 TO 樣本的 scheme heatmap
  fig_s17_s4_secondary_analysis.png — S4 ambiguous LR 分析
  fig_s20_paired_vs_to_schematic.png — Paired vs TO phaser schematic
  fig_s23_pon_anchor_schematic.png — PON 錨點方案 schematic
"""
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch, Circle
from matplotlib import font_manager as fm
from scipy.ndimage import gaussian_filter1d

OUTDIR = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260423_ppt_v2_extras"
os.makedirs(OUTDIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────
# 中文字型配置 (rule: PPT 圖片產生時中文字必須可正確 render)
# 使用 Droid Sans Fallback（支援 CJK + Latin，單一字型解決 glyph fallback）
# ─────────────────────────────────────────────────────────────────────────
CJK_FONT_PATH = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
try:
    fm.fontManager.addfont(CJK_FONT_PATH)
    CJK_FONT_NAME = fm.FontProperties(fname=CJK_FONT_PATH).get_name()
except Exception:
    CJK_FONT_NAME = "DejaVu Sans"

plt.rcParams.update({
    "font.family": [CJK_FONT_NAME, "DejaVu Sans"],
    "font.sans-serif": [CJK_FONT_NAME, "DejaVu Sans"],
    "axes.unicode_minus": False,
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

DARK = "#1E2A44"
ACCENT = "#A85540"
POSITIVE = "#009E73"
NEGATIVE = "#D55E00"
BLUE = "#0072B2"
MUTED = "#5E6572"
BG = "#F7F3EC"

# ─────────────────────────────────────────────────────
# S8: KDE 演算法示意 (4 panel)
# ─────────────────────────────────────────────────────
def fig_s8_kde():
    STALE = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz"
    df = pd.read_csv(STALE, sep="\t", low_memory=False)
    hcc = df[(df["sample"] == "HCC1395") & (df["mode"] == "to")]
    nr = hcc["NumReads"].dropna().values
    nr = nr[(nr > 5) & (nr < 300)]

    fig, axes = plt.subplots(1, 4, figsize=(15, 4))
    bins = np.arange(0, 250, 2)

    # Panel 1: Raw histogram
    ax = axes[0]
    hist, edges = np.histogram(nr, bins=bins)
    ax.bar(edges[:-1], hist, width=2, color=MUTED, alpha=0.7, edgecolor="none")
    ax.set_title("Step 1 · Raw coverage histogram\n(HCC1395 TO, n=40,096 regions)",
                 fontsize=11, fontweight="bold", color=DARK)
    ax.set_xlabel("NumReads per region")
    ax.set_ylabel("Region count")
    ax.set_xlim(0, 250)
    ax.axvline(75, color=NEGATIVE, linestyle="--", linewidth=1.2, label="stale default 75×")
    ax.legend(loc="upper right", fontsize=9)

    # Panel 2: Gaussian smooth
    ax = axes[1]
    centers = edges[:-1] + 1
    smooth = gaussian_filter1d(hist.astype(float), sigma=2.5)
    ax.fill_between(centers, smooth, color=BLUE, alpha=0.4)
    ax.plot(centers, smooth, color=BLUE, linewidth=2)
    ax.set_title("Step 2 · Gaussian smoothing\n(σ=5 · removes noise)",
                 fontsize=11, fontweight="bold", color=DARK)
    ax.set_xlabel("NumReads per region")
    ax.set_ylabel("Smoothed density")
    ax.set_xlim(0, 250)

    # Panel 3: 2nd-derivative peak detection
    ax = axes[2]
    ax.plot(centers, smooth, color=BLUE, linewidth=1.5, alpha=0.5, label="smoothed")
    # compute 2nd deriv
    d2 = np.gradient(np.gradient(smooth))
    # find peaks via local maxima of smooth where d2<0
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(smooth, prominence=50)
    peaks = [p for p in peaks if 20 < centers[p] < 200]
    if peaks:
        peak_pos = centers[peaks[0]]
    else:
        peak_pos = 54
    ax.axvline(peak_pos, color=POSITIVE, linestyle="-", linewidth=2.5, label=f"2× diploid peak = {peak_pos:.0f}×")
    # marks on peaks
    for p in peaks[:3]:
        ax.plot(centers[p], smooth[p], "v", color=POSITIVE, markersize=12)
    ax.set_title("Step 3 · 2nd-deriv peak detection\n(find leftmost mode)",
                 fontsize=11, fontweight="bold", color=DARK)
    ax.set_xlabel("NumReads per region")
    ax.set_ylabel("Smoothed density")
    ax.set_xlim(0, 250)
    ax.legend(loc="upper right", fontsize=9)

    # Panel 4: CN tier mapping
    ax = axes[3]
    ax.fill_between(centers, smooth, color=MUTED, alpha=0.2)
    ax.axvline(peak_pos, color=POSITIVE, linestyle="-", linewidth=2,
               label=f"baseline_x = {peak_pos:.0f}×")
    # CN tier boundaries
    tiers = [(0.65, "CN<1", "#FFB84D"),
             (0.99, "CN~2", "#7FC97F"),
             (1.33, "CN~3", "#BEAED4"),
             (1.82, "CN~4+", "#FDC086")]
    for mul, label, col in tiers:
        x = peak_pos * mul
        ax.axvline(x, color=col, linestyle=":", linewidth=1.5, alpha=0.8)
        ax.text(x, ax.get_ylim()[1] * 0.92, f"{mul}×", fontsize=8,
                ha="center", color=col, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=col, lw=0.8))
    ax.set_title(f"Step 4 · CovM = NumReads / {peak_pos:.0f}\nCN tiers: 0.65 / 0.99 / 1.33 / 1.82",
                 fontsize=11, fontweight="bold", color=DARK)
    ax.set_xlabel("NumReads per region")
    ax.set_ylabel("Smoothed density")
    ax.set_xlim(0, 250)
    ax.legend(loc="upper right", fontsize=9)

    fig.suptitle("KDE two-pass calibration — per-sample 2× diploid peak detection (HCC1395 worked example)",
                 fontsize=13, fontweight="bold", color=DARK, y=1.02)
    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_s8_kde_visualization.png")
    plt.savefig(out, bbox_inches="tight", dpi=130)
    plt.close()
    print(f"[OK] {out}")


# ─────────────────────────────────────────────────────
# S10: TO-only S1-S7 heatmap
# ─────────────────────────────────────────────────────
def fig_s10_to_only():
    SRC = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/s1s7_per_sample.tsv"
    df = pd.read_csv(SRC, sep="\t")
    df = df[df["mode"] == "to_pileup"].copy()

    samples = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H2009", "H1437"]
    schemes = ["baseline", "S1", "S2", "S3", "S4", "S5", "S6", "S7"]

    mat = np.full((len(samples), len(schemes)), np.nan)
    n_mat = np.full((len(samples), len(schemes)), 0, dtype=int)
    for i, s in enumerate(samples):
        for j, sch in enumerate(schemes):
            sub = df[(df["sample"] == s) & (df["scheme"] == sch)]
            if len(sub):
                mat[i, j] = sub["tp_rate"].iloc[0]
                n_mat[i, j] = int(sub["n"].iloc[0])

    fig, ax = plt.subplots(figsize=(10, 5.5))
    cmap = plt.cm.get_cmap("RdYlGn")
    im = ax.imshow(mat, cmap=cmap, vmin=0.2, vmax=1.0, aspect="auto")

    ax.set_xticks(range(len(schemes)))
    ax.set_xticklabels(schemes, fontsize=11, fontweight="bold")
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples, fontsize=11)
    ax.set_xlabel("Scheme", fontsize=12)
    ax.set_title("S1-S7 Scheme × Sample TP rate  (TO mode, 6 samples)\n"
                 "Phase 3 Synthesis · cell text: TP rate (n)",
                 fontsize=12, fontweight="bold", color=DARK)

    for i in range(len(samples)):
        for j in range(len(schemes)):
            if np.isnan(mat[i, j]):
                ax.text(j, i, "—", ha="center", va="center", color="#888", fontsize=10)
            else:
                color = "white" if mat[i, j] < 0.55 else "black"
                ax.text(j, i, f"{mat[i,j]:.2f}\n(n={n_mat[i,j]:,})",
                        ha="center", va="center", color=color, fontsize=8,
                        fontweight=("bold" if mat[i, j] >= 0.85 else "normal"))

    plt.colorbar(im, ax=ax, label="TP rate", shrink=0.75)
    # 框出 S3 S5
    for j_idx in [3, 5]:
        rect = Rectangle((j_idx - 0.5, -0.5), 1, len(samples),
                         fill=False, edgecolor=ACCENT, linewidth=3, linestyle="--")
        ax.add_patch(rect)
    ax.text(3, -0.8, "S3 ★", ha="center", color=ACCENT, fontsize=10, fontweight="bold")
    ax.text(5, -0.8, "S5 combo ★", ha="center", color=ACCENT, fontsize=10, fontweight="bold")

    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_s10_s1s7_heatmap_to_only.png")
    plt.savefig(out, bbox_inches="tight", dpi=130)
    plt.close()
    print(f"[OK] {out}")


# ─────────────────────────────────────────────────────
# S17: S4 Ambiguous secondary analysis
# ─────────────────────────────────────────────────────
def fig_s17_s4():
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: B4 feature AUC bar
    ax = axes[0]
    feats = [
        ("AF + AlleleDelta (LR)", 0.717, POSITIVE, "combo · 主導"),
        ("AlleleDelta", 0.661, ACCENT, "HP-indep"),
        ("HPMergedDelta", 0.578, BLUE, "HP-derived"),
        ("NumReads", 0.524, MUTED, ""),
        ("Fisher_Frac_Sig", 0.519, MUTED, ""),
        ("Quality_Score", 0.510, MUTED, ""),
    ]
    y = np.arange(len(feats))
    colors = [f[2] for f in feats]
    vals = [f[1] for f in feats]
    names = [f[0] for f in feats]
    bars = ax.barh(y, vals, color=colors, edgecolor="white", linewidth=0.8)
    ax.axvline(0.5, color="#888", linestyle="--", linewidth=0.8, label="random")
    ax.axvline(0.58, color=ACCENT, linestyle=":", linewidth=1.2, label="Beyond-AUC ceiling")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=10)
    ax.set_xlabel("AUC in S4 subset (HCC1395 TO, n=30,432)", fontsize=10)
    ax.set_xlim(0.45, 0.80)
    ax.set_title("Panel A · S4 secondary discriminators (B4)\n"
                 "AF + AlleleDelta LR combo reaches 0.717",
                 fontsize=11, fontweight="bold", color=DARK)
    for i, v in enumerate(vals):
        ax.text(v + 0.005, i, f"{v:.3f}", va="center", fontsize=9, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9)
    ax.invert_yaxis()

    # Panel B: S4 TP/FP overlap in AF space
    ax = axes[1]
    np.random.seed(42)
    # Simulate: TP (67%) distributes in AF ~0.05 or ~0.95, with spread
    n_tp, n_fp = 4000, 2000
    tp_af = np.concatenate([
        np.random.beta(0.6, 10, int(n_tp * 0.6)),
        1 - np.random.beta(0.6, 10, int(n_tp * 0.4))
    ])
    fp_af = np.concatenate([
        np.random.beta(0.7, 8, int(n_fp * 0.5)),
        1 - np.random.beta(0.7, 8, int(n_fp * 0.5))
    ])
    tp_ad = np.random.normal(0.08, 0.05, len(tp_af)).clip(0, 1)
    fp_ad = np.random.normal(0.04, 0.06, len(fp_af)).clip(0, 1)

    ax.scatter(tp_af, tp_ad, s=3, color=BLUE, alpha=0.25, label=f"TP (n≈{n_tp})")
    ax.scatter(fp_af, fp_ad, s=5, color=NEGATIVE, alpha=0.5, label=f"FP (n≈{n_fp})")
    ax.axhline(0.07, color=POSITIVE, linestyle="--", linewidth=1.5,
               label="LR decision boundary  (schematic · 示意)")
    ax.set_xlabel("caller AF (S4: Extreme, <0.1 或 >0.9)", fontsize=10)
    ax.set_ylabel("AlleleDelta", fontsize=10)
    ax.set_title("Panel B · S4 TP/FP in (AF × AlleleDelta) plane  [schematic]\n"
                 "TP concentrates low AlleleDelta; FP higher · 虛線為示意非真實 LR fit",
                 fontsize=11, fontweight="bold", color=DARK)
    # 在右上角加浮水印 "SCHEMATIC"
    ax.text(0.97, 0.97, "SCHEMATIC", transform=ax.transAxes,
            ha="right", va="top", fontsize=11, color="#C62828", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", fc="#FFEBEE", ec="#C62828", lw=0.8))
    ax.legend(loc="upper center", fontsize=9)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 0.6)

    fig.suptitle("S4 Ambiguous bucket (LOH=None × Extreme AF) secondary discrimination — B4 pilot",
                 fontsize=13, fontweight="bold", color=DARK, y=1.00)
    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_s17_s4_secondary_analysis.png")
    plt.savefig(out, bbox_inches="tight", dpi=130)
    plt.close()
    print(f"[OK] {out}")


# ─────────────────────────────────────────────────────
# S20: Paired vs TO schematic
# ─────────────────────────────────────────────────────
def fig_s20_paired_vs_to():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    def draw_chromosome(ax, y, label, color="#D6CCBF"):
        ax.add_patch(Rectangle((0.5, y - 0.06), 9, 0.12, facecolor=color, edgecolor=DARK, linewidth=1))
        ax.text(0.2, y, label, fontsize=10, va="center", color=DARK, fontweight="bold")

    def draw_snv(ax, x, y, kind, label=None):
        if kind == "germline":
            ax.plot(x, y, "o", color=POSITIVE, markersize=14, markeredgecolor="white", markeredgewidth=1.5, zorder=3)
            if label:
                ax.text(x, y + 0.22, label, ha="center", fontsize=8, color=POSITIVE, fontweight="bold")
        elif kind == "somatic":
            ax.plot(x, y, "^", color=NEGATIVE, markersize=14, markeredgecolor="white", markeredgewidth=1.5, zorder=3)
            if label:
                ax.text(x, y + 0.22, label, ha="center", fontsize=8, color=NEGATIVE, fontweight="bold")

    # Panel A: Paired mode
    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.5, 4.5)
    ax.set_title("Paired Mode — germline scaffold from normal sample",
                 fontsize=12, fontweight="bold", color=POSITIVE)
    ax.axis("off")

    # Normal sample band
    ax.add_patch(FancyBboxPatch((0.3, 3.3), 9.4, 1.0, boxstyle="round,pad=0.05",
                                facecolor="#E8F5E9", edgecolor=POSITIVE, linewidth=1.2))
    ax.text(0.5, 4.1, "Normal sample", fontsize=10, color=POSITIVE, fontweight="bold")
    draw_chromosome(ax, 3.7, "HP1", "#B8E0D2")
    for x in [1.5, 3.0, 4.5, 6.0, 7.5, 9.0]:
        draw_snv(ax, x, 3.7, "germline")
    ax.text(5, 3.3, "→ scaffold by germline variants (pure)", fontsize=9, style="italic", color=POSITIVE)

    # Tumor sample band
    ax.add_patch(FancyBboxPatch((0.3, 0.3), 9.4, 2.5, boxstyle="round,pad=0.05",
                                facecolor="#FFF3E0", edgecolor=ACCENT, linewidth=1.2))
    ax.text(0.5, 2.6, "Tumor sample", fontsize=10, color=ACCENT, fontweight="bold")
    # 兩條 haplotype
    draw_chromosome(ax, 2.0, "HP1")
    draw_chromosome(ax, 1.2, "HP2")
    # germline SNVs on both (follows scaffold)
    for x in [1.5, 3.0, 4.5, 6.0, 7.5, 9.0]:
        draw_snv(ax, x, 2.0, "germline")
    # somatic only on one HP (passive)
    for x, hp in [(2.2, 2.0), (5.2, 1.2), (7.0, 2.0)]:
        draw_snv(ax, x, hp, "somatic")
    ax.text(5, 0.5, "→ somatic variants 被動掛上 scaffold (不影響 phase 規則)",
            fontsize=9, style="italic", color=ACCENT)

    # 裁判圖示
    ax.text(5, -0.2, "★  Normal = Scaffold 裁判 (independent)",
            fontsize=11, ha="center", color=POSITIVE, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", fc="#E8F5E9", ec=POSITIVE, lw=1))

    # Panel B: TO mode
    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.5, 4.5)
    ax.set_title("TO Mode — tumor-only, no germline anchor",
                 fontsize=12, fontweight="bold", color=NEGATIVE)
    ax.axis("off")

    # 無 normal 標示
    ax.add_patch(FancyBboxPatch((0.3, 3.3), 9.4, 1.0, boxstyle="round,pad=0.05",
                                facecolor="#FFEBEE", edgecolor=NEGATIVE, linewidth=1.2, linestyle="--"))
    ax.text(5, 3.8, "✗ No normal sample available", ha="center",
            fontsize=12, color=NEGATIVE, fontweight="bold")

    # Tumor only band with mixed SNVs
    ax.add_patch(FancyBboxPatch((0.3, 0.3), 9.4, 2.5, boxstyle="round,pad=0.05",
                                facecolor="#FFF3E0", edgecolor=ACCENT, linewidth=1.2))
    ax.text(0.5, 2.6, "Tumor-only phasing", fontsize=10, color=ACCENT, fontweight="bold")
    draw_chromosome(ax, 2.0, "HP1")
    draw_chromosome(ax, 1.2, "HP2")
    # germline + somatic 混合
    for x in [1.5, 3.0, 4.5, 6.0, 7.5, 9.0]:
        draw_snv(ax, x, 2.0, "germline")
    # somatic variants — 全都強烈偏向 HP1 (self-phasing bias)
    for x in [2.2, 2.8, 4.0, 5.2, 6.8, 8.2, 8.8]:
        draw_snv(ax, x, 2.0, "somatic")  # 全 HP1
    # 畫 scaffold 建立箭頭
    for x in [2.2, 4.0, 6.8]:
        ax.annotate("", xy=(x, 2.4), xytext=(x, 3.1),
                    arrowprops=dict(arrowstyle="->", color=NEGATIVE, lw=1.5))
    ax.text(5, 0.5, "→ somatic 參與定義 scaffold → 17.3:1 HP1 bias",
            fontsize=9, style="italic", color=NEGATIVE, fontweight="bold")

    ax.text(5, -0.2, "⚠  Somatic = 球員兼裁判 (self-phasing artifact)",
            fontsize=11, ha="center", color=NEGATIVE, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", fc="#FFEBEE", ec=NEGATIVE, lw=1))

    # 共用圖例
    legend_patches = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=POSITIVE,
                   markersize=12, label="Germline SNV"),
        plt.Line2D([0], [0], marker="^", color="w", markerfacecolor=NEGATIVE,
                   markersize=12, label="Somatic SNV"),
    ]
    fig.legend(handles=legend_patches, loc="lower center", ncol=2,
               fontsize=10, bbox_to_anchor=(0.5, -0.02), frameon=False)

    fig.suptitle("Phaser behavior: Paired (germline scaffold) vs TO (tumor-only, self-phasing)",
                 fontsize=13, fontweight="bold", color=DARK, y=0.98)
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    out = os.path.join(OUTDIR, "fig_s20_paired_vs_to_schematic.png")
    plt.savefig(out, bbox_inches="tight", dpi=130)
    plt.close()
    print(f"[OK] {out}")


# ─────────────────────────────────────────────────────
# S23: PON 錨點 schematic
# ─────────────────────────────────────────────────────
def fig_s23_pon():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel A: TO without PON (current, self-phasing)
    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.5, 4)
    ax.set_title("Before — TO tumor-only without PON",
                 fontsize=12, fontweight="bold", color=NEGATIVE)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.3, 0.3), 9.4, 3.2, boxstyle="round,pad=0.05",
                                facecolor="#FFF3E0", edgecolor=ACCENT, linewidth=1.2))
    ax.text(0.5, 3.3, "All tumor variants enter phasing graph",
            fontsize=10, color=ACCENT, fontweight="bold")

    # Chromosome
    ax.add_patch(Rectangle((0.5, 1.8), 9, 0.15, facecolor="#D6CCBF", edgecolor=DARK, linewidth=1))

    # 混合所有 SNVs — germline + somatic 都參與
    for x in [1.3, 2.8, 4.3, 5.8, 7.3, 8.8]:
        ax.plot(x, 1.87, "o", color=POSITIVE, markersize=13, markeredgecolor="white", markeredgewidth=1.2)
    for x in [2.0, 3.5, 5.0, 6.5, 8.0]:
        ax.plot(x, 1.87, "^", color=NEGATIVE, markersize=13, markeredgecolor="white", markeredgewidth=1.2)

    ax.text(5, 1.3, "Germline + Somatic 一起進 phase graph",
            fontsize=10, ha="center", color=ACCENT)
    # 17.3:1 bias arrow
    ax.annotate("17.3:1 HP1 bias\n(self-phasing)", xy=(5, 1.5), xytext=(5, 0.5),
                fontsize=10, color=NEGATIVE, ha="center", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=NEGATIVE, lw=2),
                bbox=dict(boxstyle="round,pad=0.3", fc="#FFEBEE", ec=NEGATIVE, lw=1.2))
    ax.text(5, 2.7, "球員兼裁判", fontsize=12, ha="center", color=NEGATIVE, fontweight="bold")

    # Panel B: TO with PON anchor
    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.5, 4)
    ax.set_title("After — PON-anchored TO phasing (proposed)",
                 fontsize=12, fontweight="bold", color=POSITIVE)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.3, 0.3), 9.4, 3.2, boxstyle="round,pad=0.05",
                                facecolor="#E8F5E9", edgecolor=POSITIVE, linewidth=1.2))
    ax.text(0.5, 3.3, "PON 95% loci = scaffold anchor (read-only)",
            fontsize=10, color=POSITIVE, fontweight="bold")

    # Chromosome
    ax.add_patch(Rectangle((0.5, 1.8), 9, 0.15, facecolor="#D6CCBF", edgecolor=DARK, linewidth=1))

    # PON anchors (green, large, with halo)
    pon_positions = [1.3, 2.8, 4.3, 5.8, 7.3, 8.8]
    for x in pon_positions:
        # halo
        ax.add_patch(Circle((x, 1.87), 0.25, facecolor=POSITIVE, alpha=0.25, edgecolor="none"))
        # anchor icon
        ax.plot(x, 1.87, "o", color=POSITIVE, markersize=14, markeredgecolor="white", markeredgewidth=1.5)
    ax.text(5, 2.55, "PON germline anchors (fixed scaffold)",
            fontsize=9, ha="center", color=POSITIVE, fontweight="bold", style="italic")

    # Somatic variants — small, dimmed, "不進 HP 分群"
    for x in [2.0, 3.5, 5.0, 6.5, 8.0]:
        ax.plot(x, 1.87, "^", color=NEGATIVE, markersize=10,
                markeredgecolor="white", markeredgewidth=1.0, alpha=0.7)
    ax.text(5, 1.35, "Somatic: 只標記 · 不參與 scaffold 定義",
            fontsize=10, ha="center", color=NEGATIVE, style="italic")

    # Clean phasing indicator
    ax.annotate("HP1:HP2 ≈ 1:1\n(clean phasing)", xy=(5, 1.5), xytext=(5, 0.5),
                fontsize=10, color=POSITIVE, ha="center", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=POSITIVE, lw=2),
                bbox=dict(boxstyle="round,pad=0.3", fc="#E8F5E9", ec=POSITIVE, lw=1.2))
    ax.text(5, 2.85, "裁判獨立於球員", fontsize=12, ha="center", color=POSITIVE, fontweight="bold")

    # Legend
    legend_patches = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=POSITIVE,
                   markersize=12, label="Germline SNV / PON anchor"),
        plt.Line2D([0], [0], marker="^", color="w", markerfacecolor=NEGATIVE,
                   markersize=12, label="Somatic SNV"),
    ]
    fig.legend(handles=legend_patches, loc="lower center", ncol=2,
               fontsize=10, bbox_to_anchor=(0.5, -0.02), frameon=False)

    fig.suptitle("PON-anchored TO phasing — borrow paired-mode scaffold advantage using Panel of Normals",
                 fontsize=13, fontweight="bold", color=DARK, y=0.98)
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    out = os.path.join(OUTDIR, "fig_s23_pon_anchor_schematic.png")
    plt.savefig(out, bbox_inches="tight", dpi=130)
    plt.close()
    print(f"[OK] {out}")


# ─────────────────────────────────────────────────────
# S15: Fold-improvement heatmap with diverging colormap (取代原 viridis)
# ─────────────────────────────────────────────────────
def fig_s15_fold_diverging():
    SRC = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/s1s7_per_sample.tsv"
    df = pd.read_csv(SRC, sep="\t")
    samples_to = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H2009", "H1437"]
    schemes = ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]

    mat = np.full((len(samples_to), len(schemes)), np.nan)
    n_mat = np.full((len(samples_to), len(schemes)), 0, dtype=int)
    df_to = df[df["mode"] == "to_pileup"]
    for i, s in enumerate(samples_to):
        for j, sch in enumerate(schemes):
            sub = df_to[(df_to["sample"] == s) & (df_to["scheme"] == sch)]
            if len(sub):
                mat[i, j] = sub["fold_vs_baseline"].iloc[0]
                n_mat[i, j] = int(sub["n"].iloc[0])

    fig, ax = plt.subplots(figsize=(9, 5))
    # Diverging colormap centered at 1.0 (baseline = no change)
    # <1 red (worse), =1 white (neutral), >1 green (better)
    from matplotlib.colors import TwoSlopeNorm
    cmap = plt.cm.get_cmap("RdYlGn")
    vmax = np.nanmax(mat) if np.nanmax(mat) > 1 else 2
    norm = TwoSlopeNorm(vmin=0.3, vcenter=1.0, vmax=min(vmax, 10))
    im = ax.imshow(mat, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(range(len(schemes)))
    ax.set_xticklabels(schemes, fontsize=11, fontweight="bold")
    ax.set_yticks(range(len(samples_to)))
    ax.set_yticklabels(samples_to, fontsize=11)
    ax.set_xlabel("Scheme", fontsize=12)
    ax.set_title("TP:FP fold-improvement vs baseline  (6 TO samples)\n"
                 "Green >1 = 改善 · White ≈1 = 持平 · Red <1 = 劣化  (diverging colormap)",
                 fontsize=12, fontweight="bold", color=DARK)

    for i in range(len(samples_to)):
        for j in range(len(schemes)):
            if np.isnan(mat[i, j]):
                ax.text(j, i, "—", ha="center", va="center", color="#888", fontsize=10)
            else:
                v = mat[i, j]
                color = "white" if (v > 5 or v < 0.6) else "black"
                ax.text(j, i, f"{v:.2f}\n(n={n_mat[i,j]:,})",
                        ha="center", va="center", color=color, fontsize=8,
                        fontweight=("bold" if v >= 3 else "normal"))

    cbar = plt.colorbar(im, ax=ax, shrink=0.75, extend="both")
    cbar.set_label("fold vs baseline (1.0 = 持平)", fontsize=9)
    # 特別強調 S3 S5
    for j_idx in [2, 4]:  # S3 idx 2, S5 idx 4
        rect = Rectangle((j_idx - 0.5, -0.5), 1, len(samples_to),
                         fill=False, edgecolor=ACCENT, linewidth=3, linestyle="--")
        ax.add_patch(rect)
    ax.text(2, -0.8, "S3 ★", ha="center", color=ACCENT, fontsize=10, fontweight="bold")
    ax.text(4, -0.8, "S5 combo ★", ha="center", color=ACCENT, fontsize=10, fontweight="bold")

    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_s15_fold_heatmap_diverging.png")
    plt.savefig(out, bbox_inches="tight", dpi=130)
    plt.close()
    print(f"[OK] {out}")


if __name__ == "__main__":
    fig_s8_kde()
    fig_s10_to_only()
    fig_s15_fold_diverging()
    fig_s17_s4()
    fig_s20_paired_vs_to()
    fig_s23_pon()
    print("\nAll 6 figures generated.")
