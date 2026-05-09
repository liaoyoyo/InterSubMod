#!/usr/bin/env python3
"""Generate 5 figures for the Self-Phasing 完整觀察整合報告（2026-05-08）.

Output → InterSubMod/docs/reports/validated/2026/05/figures/20260508_self_phasing_整合/
  F1_priority_bug_mechanism.png            — getVote() vector 順序檢查 + break early before/after
  F2_priority_bug_per_chr_enrichment.png   — per-chr priority bug victim 長條圖（N + ‰ 雙 panel）
  F3_binary_commit_timeline.png            — 5 commits × 2 layer 時間軸
  F4_chr19_752_victims_scatter.png         — chr19 1Mb window × victim count 散點 + 30M/27M 標註
  F5_layer15_zero_sum_4quadrant.png        — Layer 1.5 trigger 4-quadrant (germline 0/>0 × V3F/V5)

Inputs:
  - vote_dump_genome.tsv.gz × 3（baseline / V3F / V5）— 用於 F2 per-chr enrichment
  - 已知數值（PI 報告 / V5 audit）— F1/F3/F4/F5 直接 hardcode（這些不是 ad-hoc 計算結果）

Run:
  python3 InterSubMod/scripts/analysis/self_phasing_report_figures.py
"""
from __future__ import annotations

import gzip
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
import numpy as np
import pandas as pd

# CJK font injection
rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
rcParams["axes.unicode_minus"] = False

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
DUMP_DIR = REPO / "research/v5_provenance_followup/T1_2_read_level_audit"
OUT_DIR = REPO / "docs/reports/validated/2026/05/figures/20260508_self_phasing_整合"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def fig1_priority_bug_mechanism():
    """getVote() vector 順序檢查 + break early before/after."""
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(14, 6))

    # Left: baseline (priority bug)
    ax_left.set_xlim(0, 10); ax_left.set_ylim(0, 10); ax_left.axis("off")
    ax_left.set_title("baseline `getVote()` — vector 順序檢查 + break early ❌",
                      fontsize=12, color="#c0392b", pad=10)
    boxes_left = [
        ((0.5, 8.0), 9, 0.7, "vector<pair> = {", "#fff", "#333"),
        ((1.0, 7.0), 8, 0.7, "①  {HP1_1, HP2_1}  ← FIRST  (somatic pair)", "#ffd5d5", "#c0392b"),
        ((1.0, 6.2), 8, 0.7, "②  {HP3, HP2_1}  (mixed pair)", "#fff", "#666"),
        ((1.0, 5.4), 8, 0.7, "③  {HP1, HP2}    ← LAST  (germline pair)", "#fff", "#333"),
        ((0.5, 4.6), 9, 0.7, "}", "#fff", "#333"),
        ((0.5, 3.4), 9, 0.7, "for (pair : vector) {", "#fff", "#333"),
        ((1.0, 2.6), 8, 0.7, "if (votes>0) { hp = ...; break; ← 1 票就停", "#ffd5d5", "#c0392b"),
        ((0.5, 1.8), 9, 0.7, "}", "#fff", "#333"),
    ]
    for (x, y), w, h, text, fc, ec in boxes_left:
        ax_left.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
                                                    fc=fc, ec=ec, lw=1.0))
        ax_left.text(x + 0.2, y + h/2, text, fontsize=9, va="center", family="monospace")

    # 範例
    ax_left.text(5, 0.8, "範例：germline_HP2=5, somatic_HP1_1=1\n→ ① 命中 → hp=HP:i:11 → break ❌\n  germline 5 票被忽略！正解 hp=HP:i:21",
                 fontsize=9, ha="center", color="#c0392b", style="italic",
                 bbox=dict(boxstyle="round,pad=0.4", fc="#ffeaea", ec="#c0392b"))

    # Right: V3F + V5 (fix)
    ax_right.set_xlim(0, 10); ax_right.set_ylim(0, 10); ax_right.axis("off")
    ax_right.set_title("V3F + V5 `getVote()` — explicit Layer 1 / 1.5 / 2 ✅",
                       fontsize=12, color="#27ae60", pad=10)
    boxes_right = [
        ((0.5, 8.5), 9, 0.7, "Layer 1: germline only (HP1 vs HP2)", "#d5f5d5", "#27ae60"),
        ((1.0, 7.7), 8, 0.7, "if (HP1>0 || HP2>0):  germline_maj = arg max", "#fff", "#333"),
        ((0.5, 6.7), 9, 0.7, "Layer 1.5 (V5 only): germline 缺席時 fallback", "#d5e8f5", "#2980b9"),
        ((1.0, 5.9), 8, 0.7, "elif (HP1_1>0 || HP2_1>0):  germline_maj = somatic_maj", "#fff", "#333"),
        ((0.5, 4.9), 9, 0.7, "Layer 2: somatic annotation", "#fff5d5", "#d68910"),
        ((1.0, 4.1), 8, 0.7, "if (somatic_total>0):", "#fff", "#333"),
        ((1.5, 3.3), 7.5, 0.7, "hp = 11 / 21 / 33  (基於 germline_maj)", "#fff", "#333"),
        ((1.0, 2.5), 8, 0.7, "else:  hp = germline_maj", "#fff", "#333"),
    ]
    for (x, y), w, h, text, fc, ec in boxes_right:
        ax_right.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
                                                    fc=fc, ec=ec, lw=1.0))
        ax_right.text(x + 0.2, y + h/2, text, fontsize=9, va="center", family="monospace")

    ax_right.text(5, 0.8, "同樣範例：germline_HP2=5, somatic_HP1_1=1\n→ Layer 1: HP2=5 主導 → germline_maj=2\n→ Layer 2: somatic 有票 → hp=HP:i:21 ✅",
                  fontsize=9, ha="center", color="#27ae60", style="italic",
                  bbox=dict(boxstyle="round,pad=0.4", fc="#eaffea", ec="#27ae60"))

    plt.suptitle("Figure F1 — `getVote()` priority bug 機制 before / after",
                 fontsize=14, fontweight="bold", y=1.0)
    plt.tight_layout()
    out = OUT_DIR / "F1_priority_bug_mechanism.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  F1 → {out}")


def fig2_per_chr_enrichment():
    """Per-chr priority bug enrichment 雙 panel（N 排序 + ‰ 排序）."""
    # 直接從 T1.2-F1 audit 報告的數據
    data = [
        ("chr1", 2674, 4942520), ("chr2", 2792, 4525538), ("chr3", 1595, 4102109),
        ("chr4", 2451, 3737586), ("chr5", 2194, 3506133), ("chr6", 1824, 3319118),
        ("chr7", 3508, 4852872), ("chr8", 666, 3324020),  ("chr9", 1696, 2392136),
        ("chr10", 1010, 2452467), ("chr11", 776, 2213598), ("chr12", 1024, 2406447),
        ("chr13", 819, 1497933),  ("chr14", 1653, 2266658), ("chr15", 1205, 1791787),
        ("chr16", 2584, 2267135), ("chr17", 285, 1602470), ("chr18", 1084, 1700337),
        ("chr19", 752, 1069832),  ("chr20", 2101, 1609083), ("chr21", 792, 619219),
        ("chr22", 673, 1098683),  ("chrX", 630, 1719017),   ("chrY", 67, 45137),
    ]
    df = pd.DataFrame(data, columns=["chr", "victims", "total_reads"])
    df["enrichment_permille"] = df["victims"] / df["total_reads"] * 1000

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Top panel: by victim N
    df_n = df.sort_values("victims", ascending=False)
    colors_n = ["#c0392b" if c == "chr19" else ("#3498db" if c == "chr8" else "#7f8c8d") for c in df_n["chr"]]
    bars = ax1.barh(df_n["chr"], df_n["victims"], color=colors_n)
    ax1.set_xlabel("Priority bug victim count")
    ax1.set_title("Per-chr by victim N (chr7/chr2/chr1 主要分佈，chr19 排 19，chr8 冷區)",
                   fontsize=11)
    ax1.invert_yaxis()
    for bar, n in zip(bars, df_n["victims"]):
        ax1.text(n + 50, bar.get_y() + bar.get_height()/2, str(n), va="center", fontsize=8)

    # Bottom panel: by enrichment ‰
    df_p = df.sort_values("enrichment_permille", ascending=False)
    colors_p = ["#c0392b" if c == "chr19" else ("#3498db" if c == "chr8" else "#7f8c8d") for c in df_p["chr"]]
    bars2 = ax2.barh(df_p["chr"], df_p["enrichment_permille"], color=colors_p)
    ax2.set_xlabel("Priority bug enrichment ‰  (victims / total reads × 1000)")
    ax2.set_title("Per-chr by enrichment ‰ (chrY/chr20/chr21 高，chr8 冷區 0.34× avg)",
                   fontsize=11)
    ax2.invert_yaxis()
    for bar, p in zip(bars2, df_p["enrichment_permille"]):
        ax2.text(p + 0.02, bar.get_y() + bar.get_height()/2, f"{p:.2f}", va="center", fontsize=8)

    # Genome avg reference
    avg = df["victims"].sum() / df["total_reads"].sum() * 1000
    ax2.axvline(avg, color="orange", linestyle="--", alpha=0.7, label=f"genome avg = {avg:.2f}‰")
    ax2.legend(loc="lower right", fontsize=9)

    plt.suptitle("Figure F2 — Priority bug per-chromosome 分佈（HCC1395 5kHz, 全基因組 34,855 victims）",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    out = OUT_DIR / "F2_priority_bug_per_chr_enrichment.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  F2 → {out}")


def fig3_binary_commit_timeline():
    """5 commits × 2 layer 時間軸."""
    commits = [
        ("8b8c1fd", "2026-04-09", "PON-only flag", "phasing", "球員兼裁判 fix"),
        ("41ff147", "2026-04-10", "two-layer getVote", "tagging", "priority bug fix ★"),
        ("380e8d2", "2026-04-25", "INDEL guard", "tagging", "OOB UB fix"),
        ("d0bcd8c", "2026-04-30", "ploidy fix + Layer 1.5", "phasing+tagging", "Pass 2 觸發 + germline 缺席 fallback"),
        ("938f0df", "2026-04-30", "threshold 0.95→0.9", "phasing", "Pass 2 觸發門檻放寬"),
    ]

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.set_xlim(-0.5, len(commits) - 0.5); ax.set_ylim(-2, 4); ax.axis("off")

    layer_color = {"phasing": "#3498db", "tagging": "#27ae60", "phasing+tagging": "#9b59b6"}

    # Timeline arrow
    ax.annotate("", xy=(len(commits) - 0.5, 0), xytext=(-0.5, 0),
                arrowprops=dict(arrowstyle="-|>", lw=2, color="#333"))

    # Stack labels at top
    ax.text(0.5, 3.6, "baseline", fontsize=11, ha="center", color="#7f8c8d", fontweight="bold")
    ax.annotate("", xy=(0.5, 3.3), xytext=(0.5, 3.5), arrowprops=dict(arrowstyle="-", lw=1, color="#7f8c8d"))
    ax.plot([-0.3, 0.3], [3.3, 3.3], color="#7f8c8d", lw=2)

    ax.text(2, 3.6, "V3-Fixed", fontsize=11, ha="center", color="#27ae60", fontweight="bold")
    ax.plot([0.7, 2.3], [3.3, 3.3], color="#27ae60", lw=2)

    ax.text(4, 3.6, "V5", fontsize=11, ha="center", color="#9b59b6", fontweight="bold")
    ax.plot([2.7, 4.3], [3.3, 3.3], color="#9b59b6", lw=2)

    for i, (h, d, title, layer, note) in enumerate(commits):
        color = layer_color[layer]
        ax.scatter(i, 0, s=400, color=color, zorder=5, edgecolor="black", lw=1.5)
        # Label above
        ax.text(i, 1.0, h, fontsize=9, ha="center", fontweight="bold", family="monospace")
        ax.text(i, 1.5, title, fontsize=10, ha="center", color=color)
        ax.text(i, 2.0, d, fontsize=8, ha="center", color="#666", style="italic")
        # Note below
        ax.text(i, -0.6, note, fontsize=9, ha="center", color="#333", wrap=True,
                bbox=dict(boxstyle="round,pad=0.3", fc=color + "33", ec=color, alpha=0.8))
        ax.text(i, -1.4, f"layer: {layer}", fontsize=8, ha="center", color=color, style="italic")

    plt.title("Figure F3 — Self-phasing 修補 5 commits 時間軸 + 兩層對應",
              fontsize=14, fontweight="bold", pad=15)
    plt.tight_layout()
    out = OUT_DIR / "F3_binary_commit_timeline.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  F3 → {out}")


def fig4_chr19_victims_scatter():
    """chr19 1Mb window × victim count 散點."""
    # 從 T1.2 chr19 audit 數據
    chr19_top10 = [
        (30, 215, 24732), (27, 133, 23639), (16, 41, 22078), (14, 37, 23845),
        (38, 23, 18057), (21, 27, 22114), (31, 27, 24871), (56, 20, 18814),
        (18, 22, 23371), (29, 21, 24773),
    ]
    df = pd.DataFrame(chr19_top10, columns=["window_Mb", "victims", "total"])
    df["enrichment_pct"] = df["victims"] / df["total"] * 100

    fig, ax = plt.subplots(figsize=(14, 6))
    sc = ax.scatter(df["window_Mb"], df["victims"],
                     s=df["enrichment_pct"] * 200, c=df["enrichment_pct"], cmap="Reds",
                     edgecolors="black", linewidths=1, alpha=0.85)

    # 標出 SP1/SP2/SP3 位置 (chr19:17.5M / 12.5M / 12.5M)
    sp_data = [(17.5, "SP1\n113:0", "#c0392b"), (12.5, "SP2/SP3\n109:1, 108:0", "#c0392b")]
    for sp_x, label, color in sp_data:
        ax.axvline(sp_x, color=color, linestyle="--", alpha=0.5)
        ax.text(sp_x, df["victims"].max() * 0.9, label, fontsize=10, color=color,
                ha="left", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", fc="#fff", ec=color))

    # 標 30M / 27M hotspot
    for _, r in df.iterrows():
        ax.text(r["window_Mb"], r["victims"] + 8, f"{int(r['window_Mb'])}M",
                fontsize=8, ha="center", color="#666")

    ax.set_xlabel("chr19 position (1 Mb window)")
    ax.set_ylabel("Priority bug victim count")
    ax.set_title("Figure F4 — chr19 priority bug top-10 1 Mb windows\n(30M=215 / 27M=133 累計 46% chr19 victims；SP1/2/3 在 12-17 M 區段)",
                 fontsize=12, fontweight="bold")
    cbar = plt.colorbar(sc, ax=ax, label="Enrichment %")
    plt.tight_layout()
    out = OUT_DIR / "F4_chr19_752_victims_scatter.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  F4 → {out}")


def fig5_layer15_zero_sum():
    """Layer 1.5 trigger 4-quadrant: germline (0/>0) × V3F/V5 tag count."""
    # Data: zero-sum 重分配
    data = {
        "germline=0":   {"V3F": 0,          "V5": 560_881},   # Layer 1.5 trigger
        "germline>0":   {"V3F": 18_895_432, "V5": 18_334_551},
    }
    diff = {
        "germline=0":   data["germline=0"]["V5"]   - data["germline=0"]["V3F"],
        "germline>0":   data["germline>0"]["V5"]   - data["germline>0"]["V3F"],
    }

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: stacked bar (V3F vs V5 在兩個 region)
    ax = axes[0]
    regions = ["germline_vote = 0\n(21.7 M reads, 36.8%)", "germline_vote > 0\n(37.3 M reads, 63.2%)"]
    v3f_vals = [data["germline=0"]["V3F"], data["germline>0"]["V3F"]]
    v5_vals = [data["germline=0"]["V5"], data["germline>0"]["V5"]]
    x = np.arange(len(regions))
    width = 0.35

    b1 = ax.bar(x - width/2, np.array(v3f_vals) / 1e6, width, label="V3F", color="#27ae60")
    b2 = ax.bar(x + width/2, np.array(v5_vals) / 1e6, width, label="V5", color="#9b59b6")

    # Annotate values
    for bars, vals in [(b1, v3f_vals), (b2, v5_vals)]:
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                    f"{v:,}", ha="center", fontsize=8, rotation=0)

    ax.set_xticks(x); ax.set_xticklabels(regions, fontsize=10)
    ax.set_ylabel("Tagged reads (×10⁶)")
    ax.set_title("V3F vs V5 各 region tag count", fontsize=11)
    ax.legend()

    # Right: zero-sum bar (V5 - V3F)
    ax2 = axes[1]
    diff_vals = [diff["germline=0"], diff["germline>0"]]
    colors = ["#27ae60" if v > 0 else "#c0392b" for v in diff_vals]
    bars3 = ax2.bar(["germline=0\n(Layer 1.5 多 tag)", "germline>0\n(V5 重新分類少 tag)"], diff_vals, color=colors)
    for bar, v in zip(bars3, diff_vals):
        ax2.text(bar.get_x() + bar.get_width()/2, v + (15000 if v > 0 else -15000),
                 f"{v:+,}", ha="center", fontsize=11, fontweight="bold",
                 va="bottom" if v > 0 else "top")
    ax2.axhline(0, color="black", lw=0.8)
    ax2.set_ylabel("Δ tag count (V5 − V3F)")
    ax2.set_title("Zero-sum 重分配（總和 = 0）\nV5 在 germline=0 區 +560 K，germline>0 區 -560 K", fontsize=11)

    plt.suptitle("Figure F5 — Layer 1.5 zero-sum 重分配（HCC1395 5kHz 全基因組）",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    out = OUT_DIR / "F5_layer15_zero_sum_4quadrant.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  F5 → {out}")


if __name__ == "__main__":
    print(f"[F1-F5] generating figures to {OUT_DIR}")
    fig1_priority_bug_mechanism()
    fig2_per_chr_enrichment()
    fig3_binary_commit_timeline()
    fig4_chr19_victims_scatter()
    fig5_layer15_zero_sum()
    print("done.")
