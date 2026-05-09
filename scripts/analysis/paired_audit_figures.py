#!/usr/bin/env python3
"""Generate F6 + F7 figures for paired audit / §8.6 strengthening.

F6: paired vs TO HP distribution comparison (chr19)
F7: germline-absent cross-tab — V3F hp=33 vs V5 4.19:1 偏移

Output → InterSubMod/research/paired_priority_bug_audit/figures/
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np

rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
rcParams["axes.unicode_minus"] = False

OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/paired_priority_bug_audit/figures")
OUT.mkdir(parents=True, exist_ok=True)


def fig6_hp_distribution_compare():
    """Paired vs TO HP tag distribution chr19."""
    # Paired chr19 (HP:Z: 字串)
    paired = {"HP:1": 143760, "HP:2": 183309, "HP:1-1": 12401, "HP:2-1": 14504, "HP:3": 1145}
    # TO baseline chr19 從 dump 整體分布 (近似估算 from cross-tab data)
    # baseline 全 events: hp=1: 156481, hp=2: 144175, hp=11: 26190, hp=21: 13045, hp=33: 716, hp=0: 47163
    # 從 dump 計算: 462,435 reads x 1.18 events = 549,206 events
    to_baseline = {"hp=0": 47163, "hp=1": 156481, "hp=2": 144175, "hp=11": 26190, "hp=21": 13045, "hp=33": 716}

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: paired
    ax = axes[0]
    cats = ["HP:1", "HP:2", "HP:1-1", "HP:2-1", "HP:3"]
    vals = [paired[c] for c in cats]
    colors = ["#3498db", "#e74c3c", "#5dade2", "#ec7063", "#9b59b6"]
    bars = ax.bar(cats, vals, color=colors)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, v + 3000, f"{v:,}", ha="center", fontsize=9)
    ax.set_ylabel("Reads")
    ax.set_title("Paired (longphase-s) chr19 HP:Z: tag distribution\n"
                 f"germline HP1:HP2 = 1:{paired['HP:2']/paired['HP:1']:.3f} (近 1:1, 無偏移)\n"
                 f"somatic HP1-1:HP2-1 = 1:{paired['HP:2-1']/paired['HP:1-1']:.3f} (近 1:1)",
                 fontsize=10)
    ax.set_ylim(0, max(vals) * 1.15)

    # Right: TO baseline
    ax2 = axes[1]
    cats2 = ["hp=0", "hp=1", "hp=2", "hp=11", "hp=21", "hp=33"]
    vals2 = [to_baseline[c] for c in cats2]
    colors2 = ["#7f8c8d", "#3498db", "#e74c3c", "#5dade2", "#ec7063", "#9b59b6"]
    bars2 = ax2.bar(cats2, vals2, color=colors2)
    for bar, v in zip(bars2, vals2):
        ax2.text(bar.get_x() + bar.get_width()/2, v + 3000, f"{v:,}", ha="center", fontsize=9)
    ax2.set_ylabel("Events")
    ratio_germ = to_baseline['hp=1'] / to_baseline['hp=2']
    ratio_som = to_baseline['hp=11'] / to_baseline['hp=21']
    ax2.set_title("TO baseline (longphase-to) chr19 HP:i: hpResult distribution\n"
                  f"germline hp=1:hp=2 = {ratio_germ:.3f}:1\n"
                  f"somatic hp=11:hp=21 = {ratio_som:.3f}:1 (偏 HP1)",
                  fontsize=10)
    ax2.set_ylim(0, max(vals2) * 1.15)

    plt.suptitle("Figure F6 — Paired vs TO HP tag distribution chr19 (HCC1395 5kHz)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    out = OUT / "F6_paired_vs_TO_HP_distribution.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  F6 -> {out}")


def fig7_germline_absent_crosstab():
    """Germline-absent cross-tab — V3F vs V5 vs baseline 在 5,789 events 上的表現."""
    # Data from cross_tab_summary.txt
    # baseline / V3F / V5 在 germline-absent 的分布
    versions = ["baseline", "V3F", "V5"]
    hp11 = [3312, 0, 3313]      # HP1 系列
    hp21 = [791, 0, 790]         # HP2 系列
    hp33 = [80, 5789, 53]        # 純 somatic ambiguous

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: stacked bar — 三版的處理方式
    ax = axes[0]
    x = np.arange(len(versions))
    width = 0.6

    bars1 = ax.bar(x, hp11, width, label="hp=11 (HP1 系列)", color="#3498db")
    bars2 = ax.bar(x, hp21, width, bottom=hp11, label="hp=21 (HP2 系列)", color="#e74c3c")
    hp1121 = [a+b for a,b in zip(hp11, hp21)]
    bars3 = ax.bar(x, hp33, width, bottom=hp1121, label="hp=33 (純 somatic, 不選邊)", color="#9b59b6")

    for i, v in enumerate(versions):
        total = hp11[i] + hp21[i] + hp33[i]
        if hp11[i] > 100:
            ax.text(i, hp11[i]/2, f"{hp11[i]:,}", ha="center", color="white", fontweight="bold", fontsize=10)
        if hp21[i] > 100:
            ax.text(i, hp11[i] + hp21[i]/2, f"{hp21[i]:,}", ha="center", color="white", fontweight="bold", fontsize=10)
        if hp33[i] > 100:
            ax.text(i, hp1121[i] + hp33[i]/2, f"{hp33[i]:,}", ha="center", color="white", fontweight="bold", fontsize=10)
        # Annotate ratio
        if hp11[i] > 0 and hp21[i] > 0:
            r = hp11[i] / hp21[i]
            ax.text(i, total + 200, f"HP1:HP2 = {r:.2f}:1\n{'(priority bug 偏 HP1!)' if r > 2 else ''}",
                    ha="center", fontsize=9, color="#c0392b" if r > 2 else "#27ae60")
        else:
            ax.text(i, total + 200, "全 hp=33\n(保守不選邊 ✓)",
                    ha="center", fontsize=9, color="#27ae60")

    ax.set_xticks(x)
    ax.set_xticklabels(versions, fontsize=11, fontweight="bold")
    ax.set_ylabel("germline-absent events")
    ax.set_title("Germline-absent events (5,789 events) 三版處理對比\n"
                 "baseline + V5 都 4.19:1 偏 HP1; V3F 全 hp=33 保守處理",
                 fontsize=11)
    ax.legend(loc="upper right")
    ax.set_ylim(0, 7500)

    # Right: V5 與 baseline 完全相同的證據
    ax2 = axes[1]
    paired_hps = ["HP:Z:1-1", "HP:Z:2-1", "HP:Z:3"]
    paired_counts = [2040, 1588, 530]
    baseline_hp11 = [1679, 1291, 342]
    baseline_hp21 = [318, 295, 178]
    v5_hp11 = [1679, 1291, 343]
    v5_hp21 = [318, 295, 177]

    x2 = np.arange(len(paired_hps))
    w = 0.2

    ax2.bar(x2 - 1.5*w, baseline_hp11, w, label="baseline hp=11", color="#3498db", alpha=0.8)
    ax2.bar(x2 - 0.5*w, baseline_hp21, w, label="baseline hp=21", color="#e74c3c", alpha=0.8)
    ax2.bar(x2 + 0.5*w, v5_hp11, w, label="V5 hp=11", color="#3498db", hatch="///", edgecolor="black")
    ax2.bar(x2 + 1.5*w, v5_hp21, w, label="V5 hp=21", color="#e74c3c", hatch="///", edgecolor="black")

    ax2.set_xticks(x2)
    ax2.set_xticklabels(paired_hps, fontsize=10)
    ax2.set_ylabel("events")
    ax2.set_title("baseline vs V5 cross-tab 完全相同\n(V5 Layer 1.5 = priority bug feature 化, 非修補)",
                  fontsize=11)
    ax2.legend(loc="upper right", fontsize=8)

    plt.suptitle("Figure F7 — Germline-absent events: V5 Layer 1.5 設計缺陷揭露 (5/9 Step D)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    out = OUT / "F7_germline_absent_crosstab.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  F7 -> {out}")


if __name__ == "__main__":
    print(f"Generating F6+F7 figures to {OUT}")
    fig6_hp_distribution_compare()
    fig7_germline_absent_crosstab()
    print("done.")
