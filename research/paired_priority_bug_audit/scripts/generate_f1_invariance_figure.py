#!/usr/bin/env python3
"""Generate PPT-ready figure showing caller F1 invariance across baseline / V3F / V5 / V6.

Data from: docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md §3a
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Rectangle

# CJK font injection
for cand in [
    "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf",
    "/usr/share/fonts/truetype/wqy/wqy-microhei.ttc",
]:
    if Path(cand).exists():
        font_manager.fontManager.addfont(cand)
plt.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "WenQuanYi Micro Hei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# Data — 4/30 actual F1 results vs SEQC2 truth
versions = ["ClairS-TO\n(caller raw)", "baseline\nLongPhase-TO", "V3F\nbinary", "V5\nbinary", "V6\n(reuse V5 VCF)"]
f1_093 = [0.7166, 0.7166, 0.7166, 0.7166, 0.7166]
f1_06 = [0.6273, 0.6273, 0.6273, 0.6273, 0.6273]
tp_093 = [28509, 28509, 28509, 28509, 28509]
fp_093 = [11606, 11606, 11606, 11606, 11606]
fn_093 = [10938, 10938, 10938, 10938, 10938]

OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/paired_priority_bug_audit/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def panel1_f1_bars(ax):
    x = np.arange(len(versions))
    width = 0.35
    bars1 = ax.bar(x - width / 2, f1_093, width, label="purity 0.93", color="#2E7D32")
    bars2 = ax.bar(x + width / 2, f1_06, width, label="purity 0.6", color="#1565C0")

    for b, v in zip(bars1, f1_093):
        ax.annotate(f"{v:.4f}", (b.get_x() + b.get_width() / 2, v), ha="center", va="bottom", fontsize=9)
    for b, v in zip(bars2, f1_06):
        ax.annotate(f"{v:.4f}", (b.get_x() + b.get_width() / 2, v), ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(versions, fontsize=9)
    ax.set_ylabel("F1 vs SEQC2 truth", fontsize=11)
    ax.set_title("caller F1 跨 5 階段完全不變（每個小數位）", fontsize=12)
    ax.legend(loc="lower right", fontsize=10)
    ax.set_ylim(0.5, 0.8)
    ax.grid(axis="y", alpha=0.3)

    # Highlight invariance
    ax.axhline(y=0.7166, color="#2E7D32", linestyle="--", linewidth=1, alpha=0.4)
    ax.axhline(y=0.6273, color="#1565C0", linestyle="--", linewidth=1, alpha=0.4)


def panel2_tpfpfn_table(ax):
    ax.axis("off")
    headers = ["Version", "TP", "FP", "FN", "Precision", "Recall", "F1"]
    rows_093 = [
        ["ClairS-TO raw", "28,509", "11,606", "10,938", "0.7107", "0.7227", "0.7166"],
        ["A1 baseline", "28,509", "11,606", "10,938", "0.7107", "0.7227", "0.7166"],
        ["A3 V3F (no PON)", "28,509", "11,606", "10,938", "0.7107", "0.7227", "0.7166"],
        ["A5 V5", "28,509", "11,606", "10,938", "0.7107", "0.7227", "0.7166"],
        ["V6 (=V5 VCF)", "28,509", "11,606", "10,938", "0.7107", "0.7227", "0.7166"],
    ]
    rows_06 = [
        ["ClairS-TO raw", "24,190", "13,487", "15,257", "0.6420", "0.6132", "0.6273"],
        ["B1 baseline", "24,190", "13,487", "15,257", "0.6420", "0.6132", "0.6273"],
        ["B3 V3F (no PON)", "24,190", "13,487", "15,257", "0.6420", "0.6132", "0.6273"],
        ["B5 V5", "24,190", "13,487", "15,257", "0.6420", "0.6132", "0.6273"],
        ["V6 (=V5 VCF)", "24,190", "13,487", "15,257", "0.6420", "0.6132", "0.6273"],
    ]

    table_data = (
        [["@ purity 0.93", "", "", "", "", "", ""]]
        + rows_093
        + [["@ purity 0.6", "", "", "", "", "", ""]]
        + rows_06
    )

    tbl = ax.table(
        cellText=table_data,
        colLabels=headers,
        loc="center",
        cellLoc="center",
        colColours=["#E0E0E0"] * len(headers),
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.5)

    # Style purity dividers
    for i, row in enumerate(table_data):
        if row[0].startswith("@ purity"):
            for j in range(len(headers)):
                tbl[i + 1, j].set_facecolor("#FFE082")
                tbl[i + 1, j].set_text_props(weight="bold")

    ax.set_title("F1 完整 TP/FP/FN 數據（4/30 已實證）", fontsize=12)


def panel3_pipeline_diagram(ax):
    ax.axis("off")
    stages = [
        ("ClairS-TO\ncaller", "F1=0.7166", "#FFCDD2"),
        ("baseline\nphased VCF", "F1=0.7166", "#FFE0B2"),
        ("V3F\nphased VCF", "F1=0.7166", "#FFF9C4"),
        ("V5\nphased VCF", "F1=0.7166", "#C8E6C9"),
        ("V6\nuse V5 VCF", "F1=0.7166", "#B2DFDB"),
    ]
    y = 0.5
    n = len(stages)
    for i, (name, f1, color) in enumerate(stages):
        x = 0.05 + (i / (n - 1)) * 0.85
        rect = Rectangle((x - 0.07, y - 0.18), 0.14, 0.36, facecolor=color, edgecolor="black", lw=1.5)
        ax.add_patch(rect)
        ax.annotate(name, (x, y + 0.05), ha="center", va="center", fontsize=10, fontweight="bold")
        ax.annotate(f1, (x, y - 0.10), ha="center", va="center", fontsize=10, color="#1B5E20")
        if i < n - 1:
            ax.annotate(
                "",
                xy=(x + 0.10, y),
                xytext=(x + 0.07, y),
                arrowprops=dict(arrowstyle="->", color="black", lw=2),
            )
            # mechanism label
            mechanisms = [
                "phase\n(不動 FILTER)",
                "V3F\nhaplotag\n(不碰 VCF)",
                "V5 phase\n(只動 PS/GT)",
                "V6\n(重用 V5 VCF)",
            ]
            ax.annotate(
                mechanisms[i],
                (x + 0.085, y + 0.25),
                ha="center",
                fontsize=7,
                color="#424242",
                fontstyle="italic",
            )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("F1 演進證據鏈：5 階段不變", fontsize=12)


def panel4_invariance_proof(ax):
    ax.axis("off")
    text = """三層獨立證據 確認 F1 不變：

【證據 ①】直接實證 (4/30 報告)
   • 6 個 phased VCF × SEQC2 truth
   • 每個小數位 TP/FP/FN/F1 完全相同
   • 跨 purity 0.93 + 0.6 雙重驗證

【證據 ②】檔案層級 FILTER 計數
   • baseline / V2b / V5 phased VCF
   • PASS = 47,798 (3 版本相同)
   • total = 3,187,275 (3 版本相同)
   • FILTER 子類別逐項相同

【證據 ③】V6 檔案 identity 數學保證
   • V6 重用 V5 phased VCF (same file)
   • V6 haplotag 只動 BAM
   • 完全不碰 VCF
   → V6 F1 = V5 F1 by file identity

【機制證明】
   • longphase-to phase 只動 GT/PS/GT2/GT3
   • 不改 FILTER 欄位
     (PASS/LowQual/NonSomatic)
   • haplotag 完全不動 VCF
   • F1 只看 PASS variants
   → 數學上不可能變"""
    ax.annotate(
        text,
        (0.05, 0.95),
        ha="left",
        va="top",
        fontsize=10,
        family="sans-serif",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#F5F5F5", edgecolor="#9E9E9E"),
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("F1 不變性 — 三層證據 + 機制", fontsize=12)


# Build figure
fig = plt.figure(figsize=(18, 12))
gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.20)

ax1 = fig.add_subplot(gs[0, 0])
panel1_f1_bars(ax1)

ax2 = fig.add_subplot(gs[0, 1])
panel2_tpfpfn_table(ax2)

ax3 = fig.add_subplot(gs[1, 0])
panel3_pipeline_diagram(ax3)

ax4 = fig.add_subplot(gs[1, 1])
panel4_invariance_proof(ax4)

fig.suptitle(
    "caller F1 vs SEQC2 truth — ClairS-TO → baseline → V3F → V5 → V6 完全不變性驗證",
    fontsize=14,
    fontweight="bold",
)

out_path = OUT_DIR / "f1_invariance_5stage.png"
plt.savefig(out_path, dpi=160, bbox_inches="tight")
print(f"✓ figure: {out_path}")

# Also a slide-ready single panel
fig2, ax = plt.subplots(figsize=(12, 7))
panel2_tpfpfn_table(ax)
out2 = OUT_DIR / "f1_table_only.png"
plt.savefig(out2, dpi=160, bbox_inches="tight")
print(f"✓ table-only: {out2}")

fig3, ax = plt.subplots(figsize=(14, 6))
panel3_pipeline_diagram(ax)
out3 = OUT_DIR / "f1_pipeline_chain.png"
plt.savefig(out3, dpi=160, bbox_inches="tight")
print(f"✓ pipeline-chain: {out3}")
