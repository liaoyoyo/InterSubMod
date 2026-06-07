#!/usr/bin/env python3
"""V10 圖：tumor vs normal allele-ASM AUC 配對對照（證明非 copy）。數字全讀 allele_asm_aggregate.json。"""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK TC", "Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
AD = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(AD + "/allele_asm_aggregate.json"))

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.6))
T_C, N_C = "#C2410C", "#1E3A8A"

# 圖1: per-chr tumor vs normal
bc = d["by_chrom"]
chrs = sorted(bc.keys(), key=lambda c: int(c.replace("chr", "")))
tv = [bc[c]["tumor_auc_median"] for c in chrs]
nv = [bc[c]["normal_auc_median"] for c in chrs]
x = np.arange(len(chrs))
ax1.bar(x - 0.2, tv, 0.4, color=T_C, alpha=0.85, label="tumor")
ax1.bar(x + 0.2, nv, 0.4, color=N_C, alpha=0.85, label="normal (copy-clean 二倍體)")
ax1.set_xticks(x); ax1.set_xticklabels(chrs, fontsize=9)
ax1.set_ylabel("allele-ASM AUC (median)"); ax1.set_ylim(0, 1.08)
ax1.axhline(0.5, ls=":", c="gray", lw=1)
ax1.set_title("圖1: 6/6 染色體 normal > tumor\n→ copy 非分離成因（若是 copy，normal 應低）")
ax1.legend(fontsize=8, loc="lower right")

# 圖2: by tumor CN status — tumor vs normal@same loci
bs = d["by_tumor_status"]
order = [("neutral", "CN=2\nneutral"), ("gain", "gain\n(tumor 多拷貝)"), ("loh", "LOH\n(tumor 失一條)")]
labels = [o[1] for o in order]
tvs = [bs[o[0]]["tumor"]["auc_median"] for o in order]
nvs = [bs[o[0]]["normal_same_loci"]["auc_median"] for o in order]
x2 = np.arange(len(order))
ax2.bar(x2 - 0.2, tvs, 0.4, color=T_C, alpha=0.85, label="tumor")
ax2.bar(x2 + 0.2, nvs, 0.4, color=N_C, alpha=0.85, label="normal @ 同位點")
ax2.set_xticks(x2); ax2.set_xticklabels(labels, fontsize=8.5)
ax2.set_ylabel("AUC (median)"); ax2.set_ylim(0, 1.08)
ax2.set_title("圖2: tumor 的 copy 事件「降低」分離\nLOH 最低 0.78（copy artifact 假說的反面）")
ax2.legend(fontsize=8, loc="lower left")

# 圖3: 誠實 effect — real vs shuffle null (overall)
ov = d["overall"]
grp = ["tumor", "normal"]
real = [ov["tumor"]["auc_median"], ov["normal"]["auc_median"]]
shuf = [ov["tumor"]["shuffle_p95_median"], ov["normal"]["shuffle_p95_median"]]
dm = [ov["tumor"]["auc_depthmatched_median"], ov["normal"]["auc_depthmatched_median"]]
x3 = np.arange(len(grp))
ax3.bar(x3 - 0.25, real, 0.25, color="#15803D", alpha=0.85, label="real AUC")
ax3.bar(x3, dm, 0.25, color="#6B6B66", alpha=0.7, label="depth-matched")
ax3.bar(x3 + 0.25, shuf, 0.25, color="#D97757", alpha=0.7, label="shuffle null p95")
ax3.set_xticks(x3); ax3.set_xticklabels(grp, fontsize=9)
ax3.set_ylabel("AUC (median)"); ax3.set_ylim(0, 1.08)
ax3.set_title("圖3: 誠實 effect = real 超 null 的量\n(絕對 AUC 有方法樂觀；depth-matched≈real 否證 P-06)")
ax3.legend(fontsize=8, loc="lower center")

fig.suptitle("V10 決定性對照：甲基 allele 差異「不是 copy」— matched normal (copy-clean) 分離度 ≥ tumor", fontsize=12.5, y=1.04)
fig.tight_layout()
fig.savefig(AD + "/fig_allele_asm_tumor_vs_normal.png", dpi=130, bbox_inches="tight")
plt.close()
print("saved fig_allele_asm_tumor_vs_normal.png")
print(f"check: overall tumor={real[0]} normal={real[1]}; chr normal>tumor = {sum(1 for t,n in zip(tv,nv) if n>t)}/{len(chrs)}")
