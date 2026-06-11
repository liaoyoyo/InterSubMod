#!/usr/bin/env python3
"""
Schematic: WHY long-read enables allele-specific methylation that short-read cannot.
Two-panel comparison + full pipeline flow diagram.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow, FancyBboxPatch
import numpy as np

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
import matplotlib.font_manager as fm
fm.fontManager.addfont("/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf")
plt.rcParams["font.family"] = ["Droid Sans Fallback", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

fig, (axS, axL) = plt.subplots(2, 1, figsize=(14, 9))

# Genomic landmarks (schematic coordinates 0-1000)
VAR_X = 500          # somatic variant position
CPG_XS = [120, 210, 290, 380, 620, 700, 790, 880]  # CpG sites

def draw_genome_axis(ax, label):
    ax.set_xlim(-50, 1050)
    ax.set_ylim(-1, 11)
    ax.axis("off")
    # Genome line
    ax.plot([0, 1000], [10, 10], color="#334155", lw=2)
    # Variant marker
    ax.plot(VAR_X, 10, marker="v", markersize=14, color="purple", zorder=5)
    ax.text(VAR_X, 10.6, "somatic SNV\n(身分證)", ha="center", fontsize=9, color="purple", fontweight="bold")
    # CpG markers
    for cx in CPG_XS:
        ax.plot(cx, 10, marker="o", markersize=7, markerfacecolor="#16a34a", markeredgecolor="black", zorder=4)
    ax.text(CPG_XS[0]-40, 10.6, "CpG (甲基化點)", ha="left", fontsize=8, color="#16a34a")
    ax.text(-50, 9.0, label, ha="left", fontsize=12, fontweight="bold")

# ====== PANEL 1: SHORT-READ ======
draw_genome_axis(axS, "① 短讀 Short-read (Illumina / bisulfite array) — 以前的方法")
# Short reads: ~150bp = ~120 schematic units, scattered, each covers EITHER variant OR some CpG, never both well
short_reads = [
    (60, 180, "covers CpG only", "#16a34a"),
    (250, 370, "covers CpG only", "#16a34a"),
    (440, 560, "covers SNV only", "#a855f7"),
    (450, 570, "covers SNV only", "#a855f7"),
    (600, 720, "covers CpG only", "#16a34a"),
    (770, 890, "covers CpG only", "#16a34a"),
]
y = 7.5
for x0, x1, lab, c in short_reads:
    axS.add_patch(FancyBboxPatch((x0, y-0.25), x1-x0, 0.5, boxstyle="round,pad=0.02",
                                  facecolor=c, alpha=0.4, edgecolor=c))
    # show what's on the read
    for cx in CPG_XS:
        if x0 <= cx <= x1:
            axS.plot(cx, y, marker="o", markersize=5, markerfacecolor="#16a34a", markeredgecolor="k", zorder=5)
    if x0 <= VAR_X <= x1:
        axS.plot(VAR_X, y, marker="v", markersize=9, color="purple", zorder=5)
    y -= 0.9
axS.text(500, 1.2, "✗ 每條 read 太短：要嘛看到 SNV(身分證)、要嘛看到 CpG(甲基) — 無法同時\n"
                   "✗ 無法知道甲基化來自哪條染色體 → 只能算『兩條 allele 的平均』",
         ha="center", fontsize=10, color="#b91c1c",
         bbox=dict(boxstyle="round", facecolor="#fee2e2", edgecolor="#b91c1c"))

# ====== PANEL 2: LONG-READ ======
draw_genome_axis(axL, "② 長讀 Long-read (ONT + 5mC + longphase-S) — 我們用的方法")
# Long reads: each spans variant + multiple CpG + long enough to phase. 2 haplotypes.
# HP1 (germline) - methylated; HP1-1 (somatic) - hypomethylated
long_reads = [
    (30, 970, "HP1 germline", "#2563eb", True),    # methylated
    (40, 960, "HP1 germline", "#2563eb", True),
    (20, 980, "HP1-1 somatic", "#dc2626", False),  # hypomethylated, carries variant
    (35, 965, "HP1-1 somatic", "#dc2626", False),
]
y = 7.8
for x0, x1, lab, c, methylated in long_reads:
    axL.add_patch(FancyBboxPatch((x0, y-0.28), x1-x0, 0.56, boxstyle="round,pad=0.02",
                                  facecolor=c, alpha=0.25, edgecolor=c, lw=1.5))
    # CpG on read: filled black if methylated, white if not
    for cx in CPG_XS:
        if x0 <= cx <= x1:
            fc = "black" if methylated else "white"
            axL.plot(cx, y, marker="o", markersize=6, markerfacecolor=fc, markeredgecolor="k", zorder=5)
    # variant: somatic reads carry ALT (red), germline carry REF
    if "somatic" in lab:
        axL.plot(VAR_X, y, marker="v", markersize=11, color="purple", zorder=6)
    axL.text(x1+12, y, lab, va="center", fontsize=8, color=c, fontweight="bold")
    y -= 0.95
axL.text(500, 1.4, "✓ 每條 read 夠長：同時帶 SNV(身分證) + 沿路所有 CpG(甲基) + 夠長能 phase\n"
                   "✓ 能分『腫瘤 allele vs 正常 allele 各自的甲基化』→ allele-specific resolution\n"
                   "黑點=甲基化  白點=去甲基  → HP1-1(腫瘤,紅)明顯比 HP1(正常,藍)少甲基",
         ha="center", fontsize=10, color="#15803d",
         bbox=dict(boxstyle="round", facecolor="#dcfce7", edgecolor="#15803d"))

plt.suptitle("為什麼 long-read 能做到 allele-specific methylation，short-read 不行",
             fontsize=14, fontweight="bold", y=0.99)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(f"{BASE}/figures/method_shortread_vs_longread.png", dpi=130, bbox_inches="tight")
print("Saved figures/method_shortread_vs_longread.png")
