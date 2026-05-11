#!/usr/bin/env python3
"""
Generate Thread C quantitative figures for weekly report 2026-04-23.

Source data: hardcoded from docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_HCC1395_01.md
(Phase 1 HCC1395 TO full validation)

Outputs: docs/experiments/in_progress/2026/04/figures/20260423_weekly_thread_c/
  fig_c1_auc_before_after.png          Bar chart: 18 HP-related features AUC flag=off vs on
  fig_c2_ngroups_distribution.png      Stacked: NG=0..4 counts TP/FP x off/on
  fig_c3_hpratio_shift.png             Boxplot: HP_Ratio by LOH status x flag
  fig_c4_audit_independence.png        Bar chart: NHP_Somatic11/21/33 sum TP/FP x off/on
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUTDIR = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260423_weekly_thread_c"
os.makedirs(OUTDIR, exist_ok=True)

# Style
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "figure.dpi": 110,
})

C_OFF = "#5B8DEF"  # blue
C_ON = "#E07B60"   # coral
C_DARK = "#1E2A44"
C_ACCENT = "#A85540"
C_GREEN = "#009E73"
C_RED = "#D55E00"

# -----------------------------------------------------------------------------
# Figure C1 : AUC before/after bar chart (18 HP-related features)
# -----------------------------------------------------------------------------

auc_rows = [
    # (feature, auc_off, auc_on, category)
    ("HP_Ratio", 0.5257, 0.5137, "HP-derived"),
    ("HPFineNGroups", 0.5359, 0.5099, "HP-derived"),
    ("HeuristicScore", 0.5066, 0.5088, "Composite"),
    ("Quality_Score", 0.5251, 0.5148, "Composite"),
    ("Coverage_Multiple", 0.5126, 0.5126, "HP-independent"),
    ("CramersV", 0.5004, 0.5010, "HP-derived"),
    ("CramersV_HPFine", 0.5019, 0.5012, "HP-derived"),
    ("HPFineF", 0.5117, 0.5045, "HP-derived"),
    ("HP1FamilyN", 0.5086, 0.5171, "HP-derived"),
    ("HP2FamilyN", 0.5420, 0.5295, "HP-derived"),
    ("NHP3", 0.5350, 0.5000, "HP-derived"),
    ("NHP0", 0.5563, 0.5528, "HP-derived"),
    ("NumReads", 0.5085, 0.5085, "HP-independent"),
    ("Fisher_Frac_Sig", 0.5047, 0.5127, "Methylation"),
    ("Stability", 0.5000, 0.5000, "Composite"),
    ("HPMergedDelta", 0.5406, 0.5154, "HP-derived"),
    ("AlleleDelta", 0.6294, 0.6294, "HP-independent"),
    ("HPFine_NGroups_CF", 0.5359, 0.5099, "HP-derived"),
]
df_auc = pd.DataFrame(auc_rows, columns=["feature", "auc_off", "auc_on", "category"])
df_auc["delta"] = df_auc["auc_on"] - df_auc["auc_off"]
df_auc = df_auc.sort_values("delta", ascending=True).reset_index(drop=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 7), gridspec_kw={"width_ratios": [1.3, 1]})

# Panel A: side-by-side bars
ax = axes[0]
y = np.arange(len(df_auc))
bar_h = 0.38
ax.barh(y - bar_h / 2, df_auc["auc_off"], bar_h, color=C_OFF, label="flag=off (baseline)", edgecolor="white")
ax.barh(y + bar_h / 2, df_auc["auc_on"], bar_h, color=C_ON, label="flag=on (germline-hp-only)", edgecolor="white")
ax.axvline(0.5, color="#999", linestyle="--", linewidth=0.8, label="random (AUC=0.5)")
ax.axvline(0.58, color=C_ACCENT, linestyle=":", linewidth=1.2, label="Beyond-AUC ceiling 0.58")
ax.set_yticks(y)
ax.set_yticklabels(df_auc["feature"])
ax.set_xlim(0.48, 0.66)
ax.set_xlabel("max(AUC, 1-AUC)  — HCC1395 TO, n=40,115")
ax.set_title("Panel A · AUC before/after --germline-hp-only\n(HP-derived features dominant)")
ax.legend(loc="lower right", fontsize=9)

# Panel B: delta bar (signed)
ax = axes[1]
colors = []
for _, r in df_auc.iterrows():
    if r["delta"] >= 0.02:
        colors.append(C_GREEN)
    elif r["delta"] <= -0.02:
        colors.append(C_RED)
    else:
        colors.append("#999")
ax.barh(y, df_auc["delta"], 0.65, color=colors, edgecolor="white")
ax.axvline(0, color="black", linewidth=0.7)
ax.axvline(0.02, color=C_GREEN, linestyle="--", linewidth=1, alpha=0.7, label="Gate ≥ +0.02")
ax.axvline(-0.02, color=C_RED, linestyle="--", linewidth=1, alpha=0.7)
ax.set_yticks(y)
ax.set_yticklabels(df_auc["feature"])
ax.set_xlim(-0.04, 0.04)
ax.set_xlabel("ΔAUC (on − off)")
ax.set_title("Panel B · Signed ΔAUC\nNone reach +0.02 Gate; 4 features drop ≤ −0.025")
ax.legend(loc="lower right", fontsize=9)
for i, r in df_auc.iterrows():
    x = r["delta"]
    ha = "left" if x >= 0 else "right"
    dx = 0.001 if x >= 0 else -0.001
    ax.text(x + dx, i, f"{x:+.4f}", va="center", ha=ha, fontsize=8, color=C_DARK)

fig.suptitle("Thread C · Phase 1 AUC Gate — all 18 HP-related features fail the +0.02 Gate",
             fontsize=14, fontweight="bold", color=C_DARK)
plt.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(OUTDIR, "fig_c1_auc_before_after.png")
plt.savefig(out, bbox_inches="tight", dpi=150)
plt.close()
print(f"  [OK] {out}")

# -----------------------------------------------------------------------------
# Figure C2 : HPFineNGroups distribution TP/FP x flag
# -----------------------------------------------------------------------------
# From Phase 1 report §驗證 5 (HPFineNGroups distribution convergence table)
ng_data = {
    "NGroups": [0, 1, 2, 3, 4],
    "TP_off": [11, 140, 5626, 18729, 4003],
    "TP_on":  [38, 742, 27729, 0, 0],
    "FP_off": [1, 34, 3423, 6385, 1763],
    "FP_on":  [8, 80, 11518, 0, 0],
}
df_ng = pd.DataFrame(ng_data)
df_ng["TP_off_pct"] = df_ng["TP_off"] / df_ng["TP_off"].sum() * 100
df_ng["TP_on_pct"] = df_ng["TP_on"] / df_ng["TP_on"].sum() * 100
df_ng["FP_off_pct"] = df_ng["FP_off"] / df_ng["FP_off"].sum() * 100
df_ng["FP_on_pct"] = df_ng["FP_on"] / df_ng["FP_on"].sum() * 100

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
cats = ["flag=off", "flag=on"]
ng_colors = ["#D8D8D8", "#FFC896", "#E07B60", "#8A5A6E", "#3E4E78"]
ng_labels = ["NG=0", "NG=1", "NG=2", "NG=3 ⚠", "NG=4 ⚠"]

# TP
ax = axes[0]
x = np.arange(2)
bottom_tp_off = 0
bottom_tp_on = 0
for i, (lbl, col) in enumerate(zip(ng_labels, ng_colors)):
    off_v = df_ng.iloc[i]["TP_off_pct"]
    on_v = df_ng.iloc[i]["TP_on_pct"]
    ax.bar(0, off_v, 0.6, bottom=bottom_tp_off, color=col, edgecolor="white",
           label=lbl if off_v + on_v > 0.5 else None)
    ax.bar(1, on_v, 0.6, bottom=bottom_tp_on, color=col, edgecolor="white")
    if off_v > 3:
        ax.text(0, bottom_tp_off + off_v / 2, f"{off_v:.1f}%", ha="center", va="center",
                color="white" if i >= 3 else C_DARK, fontsize=9, fontweight="bold")
    if on_v > 3:
        ax.text(1, bottom_tp_on + on_v / 2, f"{on_v:.1f}%", ha="center", va="center",
                color="white" if i >= 3 else C_DARK, fontsize=9, fontweight="bold")
    bottom_tp_off += off_v
    bottom_tp_on += on_v
ax.set_xticks(x)
ax.set_xticklabels(cats)
ax.set_ylim(0, 105)
ax.set_ylabel("Proportion of regions (%)")
ax.set_title("TP regions (n=28,509 flag=off / 28,509 flag=on)\n97.3% collapse to NG=2 when flag=on")
ax.legend(loc="upper right", fontsize=9)

# FP
ax = axes[1]
bottom_off = 0
bottom_on = 0
for i, (lbl, col) in enumerate(zip(ng_labels, ng_colors)):
    off_v = df_ng.iloc[i]["FP_off_pct"]
    on_v = df_ng.iloc[i]["FP_on_pct"]
    ax.bar(0, off_v, 0.6, bottom=bottom_off, color=col, edgecolor="white")
    ax.bar(1, on_v, 0.6, bottom=bottom_on, color=col, edgecolor="white")
    if off_v > 3:
        ax.text(0, bottom_off + off_v / 2, f"{off_v:.1f}%", ha="center", va="center",
                color="white" if i >= 3 else C_DARK, fontsize=9, fontweight="bold")
    if on_v > 3:
        ax.text(1, bottom_on + on_v / 2, f"{on_v:.1f}%", ha="center", va="center",
                color="white" if i >= 3 else C_DARK, fontsize=9, fontweight="bold")
    bottom_off += off_v
    bottom_on += on_v
ax.set_xticks(x)
ax.set_xticklabels(cats)
ax.set_title("FP regions (n=11,606 flag=off / 11,606 flag=on)\nSame NG=2 collapse pattern — mechanism, not TP-specific")

fig.suptitle("Thread C · HPFineNGroups distribution collapses to NG=2 under flag=on\n(NG≥3 categories → 0 in both TP and FP; 22,732+8,148 regions reassigned)",
             fontsize=13, fontweight="bold", color=C_DARK)
plt.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.join(OUTDIR, "fig_c2_ngroups_distribution.png")
plt.savefig(out, bbox_inches="tight", dpi=150)
plt.close()
print(f"  [OK] {out}")

# -----------------------------------------------------------------------------
# Figure C3 : HP_Ratio median shift by stratum
# -----------------------------------------------------------------------------
# From Phase 1 §驗證 2 分層表
hprow = [
    ("TP · All",                  0.549, 0.529, 28509),
    ("TP · Potential_LOH=False",  0.554, 0.531, 26439),
    ("TP · Potential_LOH=True",   0.091, 0.093, 2070),
    ("FP · All",                  0.516, 0.500, 11606),
]
df_hp = pd.DataFrame(hprow, columns=["stratum", "off_med", "on_med", "n"])
df_hp["delta"] = df_hp["on_med"] - df_hp["off_med"]

fig, axes = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [1.5, 1]})

# Panel A: dumbbell
ax = axes[0]
y = np.arange(len(df_hp))
for i, r in df_hp.iterrows():
    ax.plot([r["off_med"], r["on_med"]], [i, i], color="#BBB", linewidth=2, zorder=1)
    ax.scatter(r["off_med"], i, s=140, color=C_OFF, edgecolor="white", zorder=2, label="flag=off" if i == 0 else None)
    ax.scatter(r["on_med"], i, s=140, color=C_ON, edgecolor="white", zorder=3, label="flag=on" if i == 0 else None)
    ax.text(r["off_med"], i + 0.22, f"{r['off_med']:.3f}", ha="center", fontsize=9, color=C_OFF)
    ax.text(r["on_med"], i - 0.30, f"{r['on_med']:.3f}", ha="center", fontsize=9, color=C_ON)
ax.axvline(0.5, color="black", linestyle="--", linewidth=0.8, alpha=0.5, label="balanced (0.5)")
ax.set_yticks(y)
ax.set_yticklabels([f"{r['stratum']}\n(n={r['n']:,})" for _, r in df_hp.iterrows()])
ax.set_xlim(0.0, 0.8)
ax.set_xlabel("HP_Ratio median (per-region)")
ax.set_title("Panel A · HP_Ratio median shift\nNon-LOH direction ✓ but magnitude limited; LOH stratum already extreme")
ax.legend(loc="upper right", fontsize=9)

# Panel B: delta
ax = axes[1]
colors = [C_RED if abs(d) >= 0.10 else (C_GREEN if abs(d) >= 0.05 else "#999") for d in df_hp["delta"]]
bars = ax.barh(y, df_hp["delta"], 0.55, color=colors, edgecolor="white")
ax.axvline(0, color="black", linewidth=0.7)
ax.set_yticks(y)
ax.set_yticklabels([])
ax.set_xlim(-0.15, 0.05)
ax.set_xlabel("Δ median (on − off)")
ax.set_title("Panel B · Shift magnitude\nPlan expected 0.836 → 0.55-0.65; actual <0.03")
for i, d in enumerate(df_hp["delta"]):
    ax.text(d + (-0.003 if d < 0 else 0.003), i, f"{d:+.3f}", va="center",
            ha="right" if d < 0 else "left", fontsize=9)

fig.suptitle("Thread C · HP_Ratio shift minor — plan R3 trigger condition NOT met\n(bug not in upstream hp_tag parsing)",
             fontsize=13, fontweight="bold", color=C_DARK)
plt.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.join(OUTDIR, "fig_c3_hpratio_shift.png")
plt.savefig(out, bbox_inches="tight", dpi=150)
plt.close()
print(f"  [OK] {out}")

# -----------------------------------------------------------------------------
# Figure C4 : Audit independence (NHP_Somatic11/21/33 sum)
# -----------------------------------------------------------------------------
# From §驗證 1 Audit table
audit = [
    ("NHP_Somatic11", 347213, 347213, 122307, 122307),
    ("NHP_Somatic21", 308762, 308762, 149308, 149308),
    ("NHP_Somatic33", 124096, 124096, 50081, 50081),
]
df_au = pd.DataFrame(audit, columns=["tag", "tp_off", "tp_on", "fp_off", "fp_on"])

fig, ax = plt.subplots(figsize=(10, 5.5))
x = np.arange(len(df_au))
bar_w = 0.2
ax.bar(x - 1.5 * bar_w, df_au["tp_off"], bar_w, color=C_OFF, alpha=0.95, label="TP flag=off")
ax.bar(x - 0.5 * bar_w, df_au["tp_on"], bar_w, color=C_OFF, alpha=0.55, label="TP flag=on", edgecolor=C_DARK, linewidth=0.8)
ax.bar(x + 0.5 * bar_w, df_au["fp_off"], bar_w, color=C_ON, alpha=0.95, label="FP flag=off")
ax.bar(x + 1.5 * bar_w, df_au["fp_on"], bar_w, color=C_ON, alpha=0.55, label="FP flag=on", edgecolor=C_DARK, linewidth=0.8)

for i, r in df_au.iterrows():
    ax.text(i - 1.5 * bar_w, r["tp_off"] + 5000, f"{r['tp_off']:,}", ha="center", fontsize=8)
    ax.text(i + 0.5 * bar_w, r["fp_off"] + 5000, f"{r['fp_off']:,}", ha="center", fontsize=8)
    ax.text(i, r["tp_off"] * 0.55, "✓ identical", ha="center", fontsize=9,
            color=C_GREEN, fontweight="bold")

ax.set_xticks(x)
ax.set_xticklabels(df_au["tag"])
ax.set_ylabel("Genome-wide sum of somatic tag counts")
ax.set_title("Thread C · Audit independence ✓\nSomatic tag counts identical in flag=off / flag=on (TP+FP merged)\n"
             "→ Mechanism valid: flag demotes rather than removes somatic reads",
             fontsize=13, fontweight="bold", color=C_DARK)
ax.legend(loc="upper right", ncol=2, fontsize=9)
ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.tight_layout()
out = os.path.join(OUTDIR, "fig_c4_audit_independence.png")
plt.savefig(out, bbox_inches="tight", dpi=150)
plt.close()
print(f"  [OK] {out}")

print("\nAll Thread C figures generated.")
