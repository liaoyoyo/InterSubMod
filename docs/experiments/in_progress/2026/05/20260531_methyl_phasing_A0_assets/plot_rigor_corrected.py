#!/usr/bin/env python3
"""修正後決定性圖：T3 matched 配對檢定 + 3-target verdict。數字自 rigor_t3_chr*.json 算（與獨立驗證一致）。"""
import json, glob, os
import numpy as np
from scipy.stats import wilcoxon
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK TC", "Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
AD = os.path.dirname(os.path.abspath(__file__))

sub_noise, ctrl_noise = [], []
for f in glob.glob(AD + "/rigor_t3_chr*.json"):
    for r in json.load(open(f))["rows"]:
        if "t3_h1" in r and "noise_floor_1" in r: sub_noise.append(r["t3_h1"]["dbeta"] - r["noise_floor_1"]["dbeta"])
        if "ctrl_h1_vs_h2" in r and "noise_floor_1" in r: ctrl_noise.append(r["ctrl_h1_vs_h2"]["dbeta"] - r["noise_floor_1"]["dbeta"])
sn, cn = np.array(sub_noise), np.array(ctrl_noise)
p_sn, p_cn = wilcoxon(sn).pvalue, wilcoxon(cn).pvalue

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5))
# 左：matched 配對差 boxplot
bp = ax1.boxplot([sn, cn], labels=[f"T3 亞群−噪音\n(同窗配對)\nmedian={np.median(sn):.4f}\np={p_sn:.3f}",
                                    f"正控 H1vsH2−噪音\n(同窗配對)\nmedian={np.median(cn):.4f}\np={p_cn:.0e}"],
                 patch_artist=True, showfliers=False, widths=0.5)
bp["boxes"][0].set_facecolor("#C2410C"); bp["boxes"][0].set_alpha(0.6)
bp["boxes"][1].set_facecolor("#15803D"); bp["boxes"][1].set_alpha(0.6)
ax1.axhline(0, ls="--", c="gray", lw=1)
ax1.set_ylabel("同窗 Δβ 配對差 (target − 噪音地板)")
ax1.set_title("T3 修正：matched 配對檢定\n亞群−噪音≈0 (p=0.92, 不可區分) | 正控存活 (p≈0)\n→ 亞群無獨特甲基, 用戶必要條件不滿足")

# 右：3-target verdict
labs = ["T1 unphase→H1/H2", "T2 H3→分支", "T3 拆亞群"]
# 用「驗證強度」概念條：T1 genome-wide 0.885; T2 只在有真值 0.90/H3未證; T3 negative
vals = [0.885, 0.90, 0.50]
colors = ["#15803D", "#A16207", "#C2410C"]
notes = ["[SUPPORTED+caveats]\nheadline=V6全基因組0.885\n0.926僅可分子集\nLOH未測",
         "[OVERSTATED]\n0.90僅在有真值1-1/2-1\nH3外推未證:無AUC/null\n僅15-18%可指派",
         "[NEGATIVE]\nmatched p=0.924\n亞群=噪音\n必要條件不滿足"]
ax2.bar(range(3), vals, color=colors, alpha=0.85)
ax2.axhline(0.5, ls=":", c="gray", lw=1)
for i, (v, n) in enumerate(zip(vals, notes)):
    ax2.text(i, 0.06, n, ha="center", fontsize=8, va="bottom")
ax2.set_xticks(range(3)); ax2.set_xticklabels(labs, fontsize=9.5); ax2.set_ylim(0, 1.0)
ax2.set_ylabel("驗證強度 (AUC/acc; T3=隨機線)")
ax2.set_title("4-agent 對抗審查後最終 verdict\n(2 blocking 全在 T2; 全線單樣本 -> tier 3 上限)")
fig.tight_layout(); fig.savefig(AD + "/fig_rigor_corrected.png", dpi=130, bbox_inches="tight"); plt.close()
print(f"saved fig_rigor_corrected.png")
print(f"T3 matched: 亞群−噪音 median={np.median(sn):.4f} p={p_sn:.4f} (n={len(sn)}); 正控−噪音 median={np.median(cn):.4f} p={p_cn:.3g} (n={len(cn)})")
