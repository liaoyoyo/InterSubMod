"""
gen_concept_figures.py — HKU 合作 ISM 概念示意圖產生器

Outputs 6 PNG figures into ../figures/:
    F1_pipeline_flow.png        — 6-node ISM pipeline 全流程
    F2_two_locus_6state.png     — 兩位點 R/A/X 3x3 matrix
    F3_missing_value_rules.png  — X 缺失值補完 4 panels
    F4_hp_subdivision_tree.png  — HP1/HP2 細分樹（含 HCC1395 baseline 範例）
    F5_ism_readcpg_heatmap.png  — read × CpG 矩陣 + dendrogram
    F6_caller_window_vs_phasing.png — caller ±128bp 視窗 vs 跨視窗 SNV

Design constraints:
    - matplotlib + seaborn，CJK 字型 (Droid Sans Fallback / Noto CJK)
    - DPI=150, PNG, 16:9 friendly
    - 每張圖獨立可讀 (title / subtitle / legend / source note)
    - 示意數據（非真實 HCC1395 數字部分皆標 "示意"）

Repro:
    cd /big7_disk/liaoyoyo2001/InterSubMod
    python research/hku_collaboration/scripts/gen_concept_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure project root is importable for plot_setup
PROJECT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.lib.plot_setup import setup_plot_style  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle  # noqa: E402
import numpy as np  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy.cluster.hierarchy import linkage, dendrogram  # noqa: E402

# Apply central font chain (CJK fallback)
setup_plot_style(base_size=11, dpi=150)

# Override: matplotlib 3.6.2 sometimes fails to traverse the chain when DejaVu Sans
# is first. Force a CJK-capable font to lead so Chinese glyphs render correctly.
# Both Noto Sans CJK JP and Droid Sans Fallback contain full CJK + Latin glyphs.
import matplotlib  # noqa: E402

# Minimal 2-font chain (verified per-glyph fallback works on matplotlib 3.6.2):
#   - DejaVu Sans for Latin (numbers/letters/punctuation)
#   - Droid Sans Fallback for CJK (full Traditional + Simplified coverage)
# Adding Noto Sans CJK JP in between breaks fallback (covers most but missing a few TC chars).
_cjk_first_chain = [
    "DejaVu Sans",
    "Droid Sans Fallback",
]
matplotlib.rcParams["font.family"] = _cjk_first_chain
matplotlib.rcParams["font.sans-serif"] = _cjk_first_chain

FIG_DIR = PROJECT_ROOT / "research" / "hku_collaboration" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Color palette (consistent across figures)
COLOR_INPUT = "#74C0FC"      # light blue
COLOR_PROCESS = "#FFD43B"    # yellow
COLOR_OUTPUT = "#69DB7C"     # light green
COLOR_REF = "#4DABF7"        # blue (Ref allele)
COLOR_ALT = "#FF6B6B"        # red (Alt allele)
COLOR_MISSING = "#CED4DA"    # gray (missing/X)
COLOR_HP1 = "#5C7CFA"        # HP1
COLOR_HP2 = "#FA5252"        # HP2
COLOR_HP3 = "#868E96"        # HP3 (low confidence)

SOURCE_NOTE = "Source: HKU collaboration concept figures · gen_concept_figures.py"


# ============================================================
# F1. Pipeline 流程圖
# ============================================================
def fig1_pipeline_flow():
    fig, ax = plt.subplots(figsize=(13, 7.3), dpi=150)
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 7.3)
    ax.axis("off")

    # 6 nodes: (x, y, label, sub_label, color, role)
    nodes = [
        (1.7, 6.4, "mod.bam", "含 5mC 甲基訊號 (MM/ML tag)", COLOR_INPUT, "Input"),
        (1.7, 5.0, "tp.vcf", "somatic SNV 候選 (paired ClairS pileup)", COLOR_PROCESS, "Caller"),
        (1.7, 3.6, "tagged BAM", "HP1 / HP2 / HP1-1 / HP2-1 / HP3", COLOR_PROCESS, "Phasing"),
        (1.7, 2.2, "read × CpG matrix", "ISM 距離矩陣 (Manhattan / Euclidean)", COLOR_PROCESS, "ISM"),
        (1.7, 0.8, "甲基聚類結果", "hierarchical clustering · 補 longphase 缺失值", COLOR_OUTPUT, "Output"),
    ]

    # Right-side metric annotations
    metrics = [
        (8.5, 6.4, "輸入", "long-read ONT BAM with methylation"),
        (8.5, 5.0, "工具", "ClairS (paired mode, pileup)"),
        (8.5, 3.6, "工具", "longphase-S (用 normal germline SNP phase)"),
        (8.5, 2.2, "metric", "read-read distance · CpG 重疊 ≥ τ"),
        (8.5, 0.8, "metric", "cluster k=2 · 補完 HP3 / X 缺失"),
    ]

    step_labels = [
        "→ paired ClairS pileup",
        "→ longphase-S (germline SNP)",
        "→ ISM read × CpG distance",
        "→ hierarchical clustering",
    ]

    # Draw nodes
    for (x, y, label, sub, color, role) in nodes:
        box = FancyBboxPatch(
            (x - 1.4, y - 0.45), 2.8, 0.9,
            boxstyle="round,pad=0.05,rounding_size=0.12",
            linewidth=1.5, edgecolor="#343A40", facecolor=color, alpha=0.85,
        )
        ax.add_patch(box)
        ax.text(x, y + 0.13, label, ha="center", va="center",
                fontsize=12, fontweight="bold", color="#212529")
        ax.text(x, y - 0.22, sub, ha="center", va="center",
                fontsize=8, color="#495057")
        # role badge top-left
        ax.text(x - 1.35, y + 0.55, f"[{role}]", fontsize=7,
                color="#868E96", style="italic")

    # Draw arrows between nodes + step label
    for i in range(len(nodes) - 1):
        y_top = nodes[i][1] - 0.45
        y_bot = nodes[i + 1][1] + 0.45
        arrow = FancyArrowPatch(
            (nodes[i][0], y_top), (nodes[i + 1][0], y_bot),
            arrowstyle="->", mutation_scale=18,
            color="#495057", linewidth=1.8,
        )
        ax.add_patch(arrow)
        ax.text(nodes[i][0] + 0.3, (y_top + y_bot) / 2,
                step_labels[i], fontsize=9, color="#495057",
                ha="left", va="center")

    # Right metric annotations
    for (x, y, tag, text) in metrics:
        ax.text(x, y + 0.13, f"[{tag}]", fontsize=8, color="#868E96",
                ha="left", va="center", style="italic")
        ax.text(x, y - 0.13, text, fontsize=9.5, color="#212529",
                ha="left", va="center")

    # Vertical separator
    ax.axvline(7.4, ymin=0.05, ymax=0.95, color="#DEE2E6", linewidth=1, linestyle="--")

    # Legend
    legend_elems = [
        mpatches.Patch(facecolor=COLOR_INPUT, edgecolor="#343A40", label="Input"),
        mpatches.Patch(facecolor=COLOR_PROCESS, edgecolor="#343A40", label="Process"),
        mpatches.Patch(facecolor=COLOR_OUTPUT, edgecolor="#343A40", label="Output"),
    ]
    ax.legend(handles=legend_elems, loc="lower right",
              fontsize=9, frameon=True, ncol=3,
              bbox_to_anchor=(0.98, -0.02))

    fig.suptitle("F1 · HKU 合作 ISM Pipeline 全流程",
                 fontsize=15, fontweight="bold", y=0.97)
    ax.text(6.5, 7.0, "從 mod.bam 甲基訊號到 ISM 補完 longphase 缺失的 5 階段資料流",
            fontsize=10, color="#495057", ha="center", style="italic")
    ax.text(13, 0.05, SOURCE_NOTE, fontsize=7, color="#ADB5BD",
            ha="right", va="bottom")

    out = FIG_DIR / "F1_pipeline_flow.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[F1] {out}")


# ============================================================
# F2. 兩位點 R/A/X 6 狀態 matrix
# ============================================================
def fig2_two_locus_matrix():
    fig = plt.figure(figsize=(13, 7.5), dpi=150)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.25)
    ax = fig.add_subplot(gs[0])
    ax_text = fig.add_subplot(gs[1])
    ax_text.axis("off")

    # Build 3x3 intensity matrix (encodes Alt signal strength: R=0, X=0.5 ambiguous, A=1)
    # rows = first locus, cols = second locus, ordering: X, Ref, Alt
    states = ["X", "Ref", "Alt"]
    # Alt signal heuristic: count of Alt across 2 loci / 2 (X counted 0.5)
    alt_signal = np.array([
        [0.25, 0.0, 0.5],   # 1st X
        [0.0, 0.0, 0.5],    # 1st Ref
        [0.5, 0.5, 1.0],    # 1st Alt
    ])
    # State cell labels
    labels_state = np.array([
        ["XX", "XR", "XA"],
        ["RX", "RR", "RA"],
        ["AX", "AR", "AA"],
    ])
    # Original HP rule mapping (示意命名，待主 agent 確認)
    labels_rule = np.array([
        ["待補(X×X)", "待補(X×Ref)", "待補(X×Alt)"],
        ["待補(Ref×X)", "HP1-1\n(原規則)", "HP1-1-1\n(同甲基群)"],
        ["待補(Alt×X)", "HP1-1-2\n(異甲基群)", "HP1-2-1\n(雙 Alt)"],
    ])

    # Custom diverging colormap (blue->white->red)
    cmap = sns.diverging_palette(220, 20, as_cmap=True)

    sns.heatmap(
        alt_signal,
        annot=False, cmap=cmap, vmin=0, vmax=1,
        cbar_kws={"label": "Alt 信號強度 (Ref→Alt)", "shrink": 0.7},
        linewidths=1.2, linecolor="white",
        xticklabels=[f"第二點 {s}" for s in states],
        yticklabels=[f"第一點 {s}" for s in states],
        ax=ax, square=True,
    )

    # Two-line annotation: top=state, bottom=rule
    for i in range(3):
        for j in range(3):
            # state label top
            ax.text(j + 0.5, i + 0.30, labels_state[i, j],
                    ha="center", va="center",
                    fontsize=14, fontweight="bold",
                    color="white" if alt_signal[i, j] > 0.6 or alt_signal[i, j] < 0.05 else "#212529")
            # rule label bottom (smaller)
            ax.text(j + 0.5, i + 0.68, labels_rule[i, j],
                    ha="center", va="center",
                    fontsize=8,
                    color="white" if alt_signal[i, j] > 0.6 else "#495057",
                    style="italic")

    ax.set_title("3 × 3 = 9 cell · 對角為兩位點同型, 非對角為混合型",
                 fontsize=10, color="#495057", pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=10)

    # Right inset rule explanation
    # NOTE: avoid "─" (U+2500) — matplotlib + Droid Sans Fallback rendering
    # produces text duplication artifact. Use plain "-" instead.
    rules_text = (
        "改進規則 (placeholder — 待主 agent / 用戶確認)\n"
        + "-" * 50 + "\n\n"
        "RR  → HP1-1   (Ref-Ref · 同單體型)\n"
        "RA  → HP1-1-1 (Ref-Alt · 第二點 somatic)\n"
        "AR  → HP1-1-2 (Alt-Ref · 第一點 somatic)\n"
        "AA  → HP1-2-1 (雙 Alt · 同 somatic 群)\n\n"
        "X 狀態 (longphase 無法判定):\n"
        "  RX / XR / AX / XA / XX\n"
        "  → 進入 F3 ISM 甲基聚類補完\n\n"
        "顏色編碼: 紅 = Alt 信號強, 藍 = Ref 信號強\n"
        "         白 = X 缺失或 Ref/Ref 中性\n\n"
        "Note: 示意命名, 非最終 schema."
    )
    ax_text.text(
        0.02, 0.98, rules_text,
        transform=ax_text.transAxes,
        fontsize=10, va="top", ha="left",
        # NOTE: cannot use family="monospace" — breaks CJK fallback chain.
        # Default sans-serif chain (Droid Sans Fallback) covers all CJK glyphs.
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#F8F9FA",
                  edgecolor="#DEE2E6", linewidth=1.2),
    )

    fig.suptitle("F2 · 兩位點 R/A/X 六狀態組合 matrix",
                 fontsize=15, fontweight="bold", y=0.99)
    fig.text(0.5, 0.93,
             "每條 read 在相鄰兩個 SNV 位點的等位基因組合 → 9 種狀態 → 對應 HP 細分規則",
             ha="center", fontsize=10, color="#495057", style="italic")
    fig.text(0.99, 0.01, SOURCE_NOTE, fontsize=7, color="#ADB5BD",
             ha="right", va="bottom")

    out = FIG_DIR / "F2_two_locus_6state.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[F2] {out}")


# ============================================================
# F3. 缺失值 X 補完規則圖 (2x2 panels)
# ============================================================
def fig3_missing_value_rules():
    fig, axes = plt.subplots(2, 2, figsize=(14, 8.5), dpi=150)
    fig.subplots_adjust(hspace=0.4, wspace=0.25)

    panels = [
        ("RX", "RR or RA", "第二點 X → 用甲基聚類分到 RR 或 RA",
         [COLOR_REF, COLOR_MISSING], [COLOR_REF, COLOR_REF, COLOR_REF, COLOR_ALT]),
        ("AX", "AA or AR", "第二點 X → 用甲基聚類分到 AA 或 AR",
         [COLOR_ALT, COLOR_MISSING], [COLOR_ALT, COLOR_ALT, COLOR_ALT, COLOR_REF]),
        ("XA", "AA or RA", "第一點 X → 用甲基聚類分到 AA 或 RA",
         [COLOR_MISSING, COLOR_ALT], [COLOR_ALT, COLOR_ALT, COLOR_REF, COLOR_ALT]),
        ("XR", "RR or AR", "第一點 X → 用甲基聚類分到 RR 或 AR",
         [COLOR_MISSING, COLOR_REF], [COLOR_REF, COLOR_REF, COLOR_ALT, COLOR_REF]),
    ]

    for ax, (start, end, title, before_colors, after_colors) in zip(axes.flat, panels):
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 5)
        ax.axis("off")

        # Title
        ax.text(5, 4.65, f"{start}  →  {end}", ha="center", va="center",
                fontsize=15, fontweight="bold", color="#212529")
        ax.text(5, 4.2, title, ha="center", va="center",
                fontsize=9.5, color="#495057", style="italic")

        # Left: "未補" state (single read with X)
        ax.text(1.5, 3.55, "未補狀態", ha="center", fontsize=10,
                fontweight="bold", color="#495057")
        # Draw 2 SNV positions as circles
        for i, c in enumerate(before_colors):
            circle = plt.Circle((1.0 + i * 1.0, 2.5), 0.25,
                                color=c, ec="#343A40", linewidth=1.5)
            ax.add_patch(circle)
            label = ["X", "R", "A"][[COLOR_MISSING, COLOR_REF, COLOR_ALT].index(c)]
            ax.text(1.0 + i * 1.0, 2.5, label, ha="center", va="center",
                    fontsize=10, fontweight="bold", color="white")
        # read line
        ax.plot([0.5, 2.5], [2.5, 2.5], color="#495057", linewidth=2, zorder=0)
        ax.text(1.5, 1.85, "一條 read\n(longphase 無法判定 X 位點)",
                ha="center", va="center", fontsize=8, color="#868E96")

        # Center: ISM cluster icon
        ax.text(5.0, 3.55, "ISM 甲基聚類", ha="center", fontsize=10,
                fontweight="bold", color="#495057")
        # Draw mini heatmap icon
        mini = np.array([[0.1, 0.9, 0.2, 0.85, 0.15],
                         [0.85, 0.15, 0.9, 0.1, 0.95],
                         [0.2, 0.8, 0.15, 0.9, 0.2]])
        for i in range(3):
            for j in range(5):
                rect = Rectangle((4.2 + j * 0.18, 2.85 - i * 0.18),
                                 0.16, 0.16,
                                 facecolor=plt.cm.RdYlBu_r(mini[i, j]),
                                 edgecolor="#DEE2E6", linewidth=0.4)
                ax.add_patch(rect)
        ax.text(5.0, 2.2, "read × CpG\nmethyl. cluster",
                ha="center", fontsize=8, color="#495057")
        # arrow center
        arrow_in = FancyArrowPatch((2.8, 2.5), (4.1, 2.5),
                                    arrowstyle="->", mutation_scale=14,
                                    color="#868E96", linewidth=1.5)
        ax.add_patch(arrow_in)
        arrow_out = FancyArrowPatch((5.95, 2.5), (7.3, 2.5),
                                     arrowstyle="->", mutation_scale=14,
                                     color="#868E96", linewidth=1.5)
        ax.add_patch(arrow_out)

        # Right: 補完後 (2 reads, one to each cluster)
        ax.text(8.5, 3.55, f"補完後 → {end}", ha="center",
                fontsize=10, fontweight="bold", color="#495057")
        # Top read (first option)
        for i, c in enumerate(after_colors[:2]):
            circle = plt.Circle((7.7 + i * 1.0, 3.0), 0.2,
                                color=c, ec="#343A40", linewidth=1.5)
            ax.add_patch(circle)
            label = ["X", "R", "A"][[COLOR_MISSING, COLOR_REF, COLOR_ALT].index(c)]
            ax.text(7.7 + i * 1.0, 3.0, label, ha="center", va="center",
                    fontsize=9, fontweight="bold", color="white")
        ax.plot([7.4, 9.0], [3.0, 3.0], color="#495057", linewidth=2, zorder=0)
        # Bottom read (second option)
        for i, c in enumerate(after_colors[2:]):
            circle = plt.Circle((7.7 + i * 1.0, 2.0), 0.2,
                                color=c, ec="#343A40", linewidth=1.5)
            ax.add_patch(circle)
            label = ["X", "R", "A"][[COLOR_MISSING, COLOR_REF, COLOR_ALT].index(c)]
            ax.text(7.7 + i * 1.0, 2.0, label, ha="center", va="center",
                    fontsize=9, fontweight="bold", color="white")
        ax.plot([7.4, 9.0], [2.0, 2.0], color="#495057", linewidth=2, zorder=0)
        # cluster labels
        opt = end.split(" or ")
        ax.text(8.2, 3.42, f"cluster A → {opt[0]}", ha="center",
                fontsize=7.5, color="#1864AB")
        ax.text(8.2, 1.58, f"cluster B → {opt[1]}", ha="center",
                fontsize=7.5, color="#C92A2A")

    # Legend (bottom)
    legend_elems = [
        mpatches.Patch(facecolor=COLOR_REF, edgecolor="#343A40", label="Ref allele"),
        mpatches.Patch(facecolor=COLOR_ALT, edgecolor="#343A40", label="Alt allele"),
        mpatches.Patch(facecolor=COLOR_MISSING, edgecolor="#343A40", label="X (missing / longphase 無法判定)"),
    ]
    fig.legend(handles=legend_elems, loc="lower center",
               fontsize=10, frameon=True, ncol=3,
               bbox_to_anchor=(0.5, -0.01))

    fig.suptitle("F3 · X 缺失值 ISM 甲基聚類補完規則 (4 種 X 模式)",
                 fontsize=15, fontweight="bold", y=0.98)
    fig.text(0.5, 0.94,
             "longphase 無法判定的位點 (X) 透過 read 周邊 CpG 甲基化模式分群, 補回 Ref 或 Alt",
             ha="center", fontsize=10, color="#495057", style="italic")
    fig.text(0.99, 0.005, SOURCE_NOTE, fontsize=7, color="#ADB5BD",
             ha="right", va="bottom")

    out = FIG_DIR / "F3_missing_value_rules.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[F3] {out}")


# ============================================================
# F4. HP 細分樹
# ============================================================
def fig4_hp_subdivision_tree():
    fig, ax = plt.subplots(figsize=(14, 8), dpi=150)
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.axis("off")

    # Tree nodes: (x, y, label, count_label, color)
    # HCC1395 baseline 真實數字: HP1=863K, HP2=266K, HP1-1=79K, HP2-1=23K, HP3=630
    nodes = {
        "root":    (7.0, 7.2, "Normal = Ref\n(germline het)", "起始", "#ADB5BD"),
        "HP1":     (3.0, 5.7, "HP1", "n ≈ 863K reads\n(HCC1395 baseline)", COLOR_HP1),
        "HP2":     (11.0, 5.7, "HP2", "n ≈ 266K reads", COLOR_HP2),
        "HP1-1":   (1.8, 4.0, "HP1-1", "n ≈ 79K\n(1 somatic on HP1)", "#748FFC"),
        "HP1-1-1": (0.8, 2.2, "HP1-1-1", "≥2 somatic\n同甲基群", "#91A7FF"),
        "HP1-1-2": (2.8, 2.2, "HP1-1-2", "≥2 somatic\n異甲基群", "#91A7FF"),
        "HP2-1":   (9.2, 4.0, "HP2-1", "n ≈ 23K\n(1 somatic on HP2)", "#FF8787"),
        "HP2-1-1": (8.2, 2.2, "HP2-1-1", "≥2 somatic\n同甲基群", "#FFA8A8"),
        "HP2-1-2": (10.2, 2.2, "HP2-1-2", "≥2 somatic\n異甲基群", "#FFA8A8"),
        "HP3":     (12.5, 4.0, "HP3", "n ≈ 630 reads\n(low confidence)", COLOR_HP3),
    }

    edges = [
        ("root", "HP1",     "1 條單純 (no somatic)"),
        ("root", "HP2",     "1 條單純 (no somatic)"),
        ("HP1", "HP1-1",    "經過 1 somatic"),
        ("HP1-1", "HP1-1-1", "+1 (同群)"),
        ("HP1-1", "HP1-1-2", "+1 (異群)"),
        ("HP2", "HP2-1",    "經過 1 somatic"),
        ("HP2-1", "HP2-1-1", "+1 (同群)"),
        ("HP2-1", "HP2-1-2", "+1 (異群)"),
        ("root", "HP3",     "low conf · drop"),
    ]

    # Draw edges first
    for parent, child, label in edges:
        x1, y1, *_ = nodes[parent]
        x2, y2, *_ = nodes[child]
        ax.plot([x1, x2], [y1 - 0.3, y2 + 0.4],
                color="#495057", linewidth=1.5, zorder=0)
        # mid label
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text(mx + 0.15, my, label, fontsize=7.5,
                color="#495057", style="italic",
                ha="left", va="center")

    # Draw nodes (rounded box)
    for key, (x, y, label, count, color) in nodes.items():
        w, h = 1.55, 0.85
        box = FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h,
            boxstyle="round,pad=0.05,rounding_size=0.12",
            linewidth=1.5, edgecolor="#212529",
            facecolor=color, alpha=0.85,
        )
        ax.add_patch(box)
        ax.text(x, y + 0.13, label, ha="center", va="center",
                fontsize=10.5, fontweight="bold", color="#212529")
        ax.text(x, y - 0.22, count, ha="center", va="center",
                fontsize=7.5, color="#212529")

    # Legend
    legend_elems = [
        mpatches.Patch(facecolor=COLOR_HP1, label="HP1 lineage"),
        mpatches.Patch(facecolor=COLOR_HP2, label="HP2 lineage"),
        mpatches.Patch(facecolor=COLOR_HP3, label="HP3 (low confidence, 通常 drop)"),
        mpatches.Patch(facecolor="#ADB5BD", label="root (germline het Ref)"),
    ]
    ax.legend(handles=legend_elems, loc="lower left", fontsize=9,
              frameon=True, bbox_to_anchor=(0.01, 0.01))

    # Annotation: counts source
    ax.text(13.95, 7.5,
            "Count 範例: HCC1395 baseline (paired pileup mode)\n"
            "    HP1 ≈ 863K · HP2 ≈ 266K · HP3 ≈ 630 reads\n"
            "    HP1-1 ≈ 79K · HP2-1 ≈ 23K (post-haplotag-fix)",
            ha="right", va="top", fontsize=8, color="#495057",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#F8F9FA",
                      edgecolor="#DEE2E6"))

    fig.suptitle("F4 · HP 細分樹 — somatic 累積 + 甲基同群/異群分支",
                 fontsize=15, fontweight="bold", y=0.97)
    ax.text(7.0, 7.85,
            "從 germline het Ref 開始, 沿 HP1/HP2 累積 somatic 變異, 再依甲基聚類異同細分至 HP*-1-1/HP*-1-2",
            ha="center", fontsize=10, color="#495057", style="italic")
    ax.text(13.95, 0.05, SOURCE_NOTE, fontsize=7, color="#ADB5BD",
            ha="right", va="bottom")

    out = FIG_DIR / "F4_hp_subdivision_tree.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[F4] {out}")


# ============================================================
# F5. ISM read × CpG matrix heatmap + dendrogram + HP annotation
# ============================================================
def fig5_ism_readcpg_heatmap():
    rng = np.random.default_rng(seed=42)
    n_reads, n_cpg = 10, 20
    # Two clusters: cluster A (HP1 dominant) low methyl in left half, cluster B (HP2) high methyl
    cluster_a_size = 6
    data = np.zeros((n_reads, n_cpg))
    # Cluster A: reads 0-5, mostly low methyl with some high
    for i in range(cluster_a_size):
        data[i, :10] = rng.beta(1.5, 6, size=10)   # low methyl on left
        data[i, 10:] = rng.beta(2.5, 3, size=10)   # mid on right
    for i in range(cluster_a_size, n_reads):
        data[i, :10] = rng.beta(6, 1.5, size=10)   # high methyl left
        data[i, 10:] = rng.beta(3, 2.5, size=10)   # mid-high right
    data = np.clip(data, 0, 1)

    # Hierarchical clustering on reads (rows)
    Z = linkage(data, method="ward")

    fig = plt.figure(figsize=(13, 7), dpi=150)
    gs = fig.add_gridspec(
        1, 4, width_ratios=[0.6, 0.3, 5, 0.3], wspace=0.05,
        left=0.08, right=0.96, top=0.86, bottom=0.12,
    )

    # 1. Dendrogram (left)
    ax_dend = fig.add_subplot(gs[0])
    ddata = dendrogram(
        Z, orientation="left", ax=ax_dend,
        color_threshold=Z[-1, 2] * 0.6,
        above_threshold_color="#868E96",
    )
    ax_dend.axis("off")
    ax_dend.invert_yaxis()
    read_order = ddata["leaves"]

    # 2. HP annotation bar
    ax_hp = fig.add_subplot(gs[1])
    # Assign HP1/HP2 by original ordering: 0-5 HP1, 6-9 HP2
    hp_orig = ["HP1"] * cluster_a_size + ["HP2"] * (n_reads - cluster_a_size)
    hp_reordered = [hp_orig[i] for i in read_order]
    hp_color_map = {"HP1": COLOR_HP1, "HP2": COLOR_HP2}
    for i, hp in enumerate(hp_reordered):
        rect = Rectangle((0, i), 1, 1,
                         facecolor=hp_color_map[hp],
                         edgecolor="white", linewidth=0.5)
        ax_hp.add_patch(rect)
        ax_hp.text(0.5, i + 0.5, hp, ha="center", va="center",
                   fontsize=7, color="white", fontweight="bold")
    ax_hp.set_xlim(0, 1)
    ax_hp.set_ylim(0, n_reads)
    ax_hp.invert_yaxis()
    ax_hp.set_xticks([])
    ax_hp.set_yticks([])
    ax_hp.set_title("HP", fontsize=9, pad=4)
    for spine in ax_hp.spines.values():
        spine.set_visible(False)

    # 3. Main heatmap
    ax_main = fig.add_subplot(gs[2])
    data_reordered = data[read_order, :]
    sns.heatmap(
        data_reordered, ax=ax_main, cmap="RdYlBu_r",
        vmin=0, vmax=1, cbar=False,
        xticklabels=[f"CpG_{i+1}" for i in range(n_cpg)],
        yticklabels=[f"read_{i+1}" for i in read_order],
        linewidths=0.3, linecolor="#F8F9FA",
    )
    ax_main.set_xlabel("CpG site (5' → 3')", fontsize=10)
    ax_main.set_ylabel("")
    ax_main.tick_params(axis="x", rotation=45, labelsize=7.5)
    ax_main.tick_params(axis="y", labelsize=8)
    for label in ax_main.get_xticklabels():
        label.set_horizontalalignment("right")

    # 4. Colorbar (right)
    ax_cb = fig.add_subplot(gs[3])
    sm = plt.cm.ScalarMappable(cmap="RdYlBu_r",
                                norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=ax_cb)
    cb.set_label("甲基化機率 P(5mC)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    # Cluster dividing line (between cluster A and B in reordered space)
    # Find first row in reordered ordering that switches HP majority
    ax_main.axhline(cluster_a_size, color="#212529",
                    linewidth=2, linestyle="--", alpha=0.6)
    ax_main.text(n_cpg + 0.5, cluster_a_size - 0.3,
                 "cluster cut", fontsize=8,
                 color="#212529", ha="left", va="center", style="italic")

    fig.suptitle("F5 · ISM read × CpG matrix · hierarchical clustering 補完 HP",
                 fontsize=14, fontweight="bold", y=0.97)
    fig.text(0.5, 0.92,
             "10 條 read × 20 個 CpG · cluster k=2 對應 HP1/HP2; 示意數據 (非真實 HCC1395)",
             ha="center", fontsize=9.5, color="#495057", style="italic")
    fig.text(0.99, 0.01, SOURCE_NOTE, fontsize=7, color="#ADB5BD",
             ha="right", va="bottom")

    out = FIG_DIR / "F5_ism_readcpg_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[F5] {out}")


# ============================================================
# F6. caller ±128bp 視窗 vs 跨視窗 SNV pair
# ============================================================
def fig6_caller_window():
    fig, ax = plt.subplots(figsize=(14, 7), dpi=150)
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 7)
    ax.axis("off")

    # Genomic axis
    axis_y = 3.0
    ax.plot([0.5, 13.5], [axis_y, axis_y], color="#343A40",
            linewidth=2.2, zorder=1)
    # Tick marks every ~2 units
    for x in [0.5, 3.0, 5.5, 7.0, 8.5, 11.0, 13.5]:
        ax.plot([x, x], [axis_y - 0.08, axis_y + 0.08],
                color="#343A40", linewidth=1.5, zorder=1)

    # caller center position
    center_x = 7.0
    ax.plot([center_x], [axis_y], marker="v", markersize=18,
            color="#E03131", zorder=3)
    ax.text(center_x, axis_y + 0.5, "caller target SNV\nchr2:18,086,020",
            ha="center", va="bottom", fontsize=10,
            fontweight="bold", color="#E03131")

    # ±128bp window (red box)
    win_left, win_right = 5.5, 8.5
    rect = Rectangle((win_left, axis_y - 0.8), win_right - win_left, 1.6,
                     facecolor="#FFE3E3", edgecolor="#E03131",
                     linewidth=2, linestyle="--", alpha=0.6, zorder=0)
    ax.add_patch(rect)
    ax.text((win_left + win_right) / 2, axis_y - 1.05,
            "caller 視窗 ±128 bp",
            ha="center", fontsize=10, color="#C92A2A",
            fontweight="bold", style="italic")

    # In-window SNVs (3 nearby — caller sees these)
    in_window = [(6.0, "snv_a"), (6.5, "snv_b"), (8.0, "snv_c")]
    for (x, name) in in_window:
        ax.plot([x], [axis_y], marker="o", markersize=10,
                markerfacecolor="#FFA94D", markeredgecolor="#212529",
                markeredgewidth=1.2, zorder=2)
        ax.text(x, axis_y - 0.32, name, ha="center", fontsize=8,
                color="#5C3D00")

    # Out-of-window SNVs (left side, 3)
    out_left = [(1.2, "snv_L1"), (2.3, "snv_L2"), (3.5, "snv_L3")]
    for (x, name) in out_left:
        ax.plot([x], [axis_y], marker="o", markersize=10,
                markerfacecolor="#74C0FC", markeredgecolor="#212529",
                markeredgewidth=1.2, zorder=2)
        ax.text(x, axis_y - 0.32, name, ha="center", fontsize=8,
                color="#1864AB")

    # Out-of-window SNVs (right side, 3)
    out_right = [(10.0, "snv_R1"), (11.5, "snv_R2"), (12.8, "snv_R3")]
    for (x, name) in out_right:
        ax.plot([x], [axis_y], marker="o", markersize=10,
                markerfacecolor="#74C0FC", markeredgecolor="#212529",
                markeredgewidth=1.2, zorder=2)
        ax.text(x, axis_y - 0.32, name, ha="center", fontsize=8,
                color="#1864AB")

    # Tag arrows: longphase tag links out-of-window into target
    # Draw curved arrows from each out-of-window SNV to target
    for (x, _) in out_left + out_right:
        arrow = FancyArrowPatch(
            (x, axis_y + 0.18), (center_x, axis_y + 1.7),
            arrowstyle="->", mutation_scale=12,
            color="#5C7CFA", linewidth=1.0, alpha=0.55,
            connectionstyle="arc3,rad=-0.2" if x < center_x else "arc3,rad=0.2",
            zorder=2,
        )
        ax.add_patch(arrow)

    # PS (phase set) annotation
    ax.text(center_x, axis_y + 1.95,
            "longphase tag (HP/PS) 串入跨視窗資訊\n(同一 PS 內可拼接 SNV)",
            ha="center", fontsize=10, color="#3B5BDB",
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#E7F5FF",
                      edgecolor="#5C7CFA", linewidth=1.2))

    # PS span bar (under axis)
    ps_y = axis_y - 1.55
    ax.plot([0.8, 13.2], [ps_y, ps_y], color="#5C7CFA",
            linewidth=4, alpha=0.5, zorder=0)
    ax.text(7.0, ps_y - 0.25, "phase set (PS) span · 通常 10-100 kb",
            ha="center", fontsize=9, color="#3B5BDB", style="italic")

    # Annotation boxes
    # Top-left: caller window meaning
    ax.text(0.3, 6.4,
            "caller (ClairS pileup) 視角:\n"
            "  ◆ 只看 target ±128 bp 內\n"
            "  ◆ 3 個 in-window SNV 可直接組合\n"
            "  ◆ 看不到視窗外 6 個 SNV\n"
            "  ◆ 限制: 視窗內訊息不足以判 phase",
            fontsize=9, color="#C92A2A", ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#FFF5F5",
                      edgecolor="#FFA8A8", linewidth=1.2))

    # Top-right: longphase tag meaning
    ax.text(13.7, 6.4,
            "longphase tag 視角:\n"
            "  ◆ 同 PS 內所有 SNV 可串接\n"
            "  ◆ 補入 6 個視窗外 SNV 為 context\n"
            "  ◆ HP1/HP2 標籤跨越 ±128 限制\n"
            "  ◆ ISM 可拿到全 PS 的 read 集合",
            fontsize=9, color="#3B5BDB", ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#E7F5FF",
                      edgecolor="#74C0FC", linewidth=1.2))

    # Legend
    legend_elems = [
        plt.Line2D([0], [0], marker="v", color="w",
                   markerfacecolor="#E03131", markersize=11, label="caller target SNV"),
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="#FFA94D", markeredgecolor="#212529",
                   markersize=9, label="in-window SNV (caller 可見)"),
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="#74C0FC", markeredgecolor="#212529",
                   markersize=9, label="out-of-window SNV (僅 phasing 可見)"),
    ]
    ax.legend(handles=legend_elems, loc="lower center",
              fontsize=9, frameon=True, ncol=3,
              bbox_to_anchor=(0.5, -0.02))

    fig.suptitle("F6 · caller ±128bp 視窗 vs longphase 跨視窗 SNV 串接",
                 fontsize=15, fontweight="bold", y=0.96)
    ax.text(7.0, 6.85,
            "caller 受 ±128bp 限制無法判 phase; longphase tag 在同 PS 內補入跨視窗 SNV 資訊",
            ha="center", fontsize=10, color="#495057", style="italic")
    ax.text(13.95, 0.05, SOURCE_NOTE, fontsize=7, color="#ADB5BD",
            ha="right", va="bottom")

    out = FIG_DIR / "F6_caller_window_vs_phasing.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[F6] {out}")


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    print(f"Output dir: {FIG_DIR}")
    fig1_pipeline_flow()
    fig2_two_locus_matrix()
    fig3_missing_value_rules()
    fig4_hp_subdivision_tree()
    fig5_ism_readcpg_heatmap()
    fig6_caller_window()
    print("All 6 figures generated.")
