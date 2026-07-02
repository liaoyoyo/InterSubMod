#!/usr/bin/env python3
"""切群方法流程圖(讓使用者確認 workflow)。純手繪 matplotlib boxes+arrows。"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams["font.family"]=["Droid Sans Fallback","DejaVu Sans"]; plt.rcParams["axes.unicode_minus"]=False
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"

fig,ax=plt.subplots(figsize=(11,12)); ax.set_xlim(0,10); ax.set_ylim(0,24); ax.axis("off")
def box(x,y,w,h,txt,fc,ec="#333",fs=9.5,tc="#111"):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.1",fc=fc,ec=ec,lw=1.3))
    ax.text(x,y,txt,ha="center",va="center",fontsize=fs,color=tc,wrap=True)
def arr(x1,y1,x2,y2,txt=""):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=16,lw=1.4,color="#555"))
    if txt: ax.text((x1+x2)/2+0.3,(y1+y2)/2,txt,fontsize=8,color="#777",ha="left")

box(5,23,8.2,1.2,"INPUT｜ISM binary 每位點: read×read BERNOULLI 距離 + read×CpG 甲基 + reads.tsv(hp/allele 標籤)","#e0f2fe",fs=9)
arr(5,22.4,5,21.6,"過濾 tumor + 有效 hp")
box(5,21,7.5,1.0,"PEEL｜移除無效距離(NaN/-1)的 read","#f0ebe2")
arr(5,20.5,5,19.7)
box(5,19,8.2,1.3,"UPGMA(average linkage) → 在『顯著 gap(>2×中位)』切 + QUARANTINE(群<3 read=離群)\n→ 列舉候選解析度 k=2..6","#fef3c7",fs=9)
arr(5,18.35,5,17.5,"每解析度算 4 個量")
# 4 gates
box(2.3,16.3,3.6,1.3,"① balance\n核心群 ≥3 read\n(業界 caller 最低)","#dde6ff",fs=8.5)
box(5.0,16.3,3.6,1.3,"② null-excess ≥0.10\nclusterboot Jaccard −\nwithin-1-group null(95%)\n= 真實/明顯","#dde6ff",fs=8.5)
box(7.7,16.3,3.6,1.3,"③ gap-ratio ≥8×\n分支跳躍\n= 結構支持","#dde6ff",fs=8.5)
box(5.0,14.5,7.5,1.0,"④ alignment｜CramérV vs hp/carrier/allele (V≥0.3, p<.05, Cochran e≥5) = 歸因","#ddeeff",fs=8.5)
arr(2.3,15.65,3.5,15.0); arr(5,15.65,5,15.0); arr(7.7,15.65,6.5,15.0)
arr(5,14.0,5,13.2,"分類(per 解析度)")
# classification
box(1.9,12.3,3.3,1.5,"CONFIRMED\n真實 + 可靠對齊 germline\n(cis 骨幹)","#cdeee6",fs=8.5)
box(5.0,12.3,3.3,1.5,"REAL_NOVEL\n真實 + 大跳, 不對齊\n= subclone 候選 ⭐","#e9d8fb",fs=8.5)
box(8.1,12.3,3.3,1.5,"REAL_DIFFUSE\n真實但無大跳+無對齊\n(低信心候選)","#fde2ef",fs=8.5)
box(5.0,10.6,4.0,0.9,"NO_CLEAR｜全無真實(excess<0.10) = 真的無法切","#e5e7eb",fs=8.5)
arr(5,11.55,5,11.05)
arr(3,11.55,3,11.0); arr(7,11.55,7,11.0)
arr(5,10.15,5,9.4,"輸出")
# output
box(2.4,8.3,3.8,1.3,"coarse\n最粗 CONFIRMED 骨幹\n(germline 層)","#d1fae5",fs=8.5)
box(7.0,8.3,3.8,1.3,"fine\n最細 supported(大跳 or 對齊)\n全真實結構","#d1fae5",fs=8.5)
box(5.0,6.6,8.0,1.1,"+ edge 群(≥3 離群成一組) + per-read 群歸屬(core_confirmed/core_novel/edge/outlier)","#d1fae5",fs=9)
arr(2.4,7.65,3.5,7.2); arr(7.0,7.65,6.0,7.2)
# 精細度天花板註
box(5.0,4.7,8.6,1.4,"🔑 精細度天花板 = ② null-excess(真實?) + ③ gap-ratio 或 ④ alignment(結構支持?)\n"
    "excess 會隨 k 無限上升 → 不能單用; gap-ratio/對齊才定『切到多細』\n"
    "「無法切」= 全無真實才算; 真實但 diffuse → REAL_DIFFUSE 不丟棄","#fff7ed",ec="#d97706",fs=8.8)
box(5.0,2.7,8.6,1.2,"⚠ 單樣本 ⭐2-3: REAL_NOVEL/DIFFUSE = subclone 候選, 需 normal cis-control 排 cis-ASM 才能確認\n"
    "⚠ GAP_TH=8× 對大樹(n大)偏鬆(chr4 仍 k6) → gap-ratio scale-invariance 待修","#fef2f2",ec="#dc2626",fs=8.8)
ax.text(5,1.0,"切群方法流程 — 三閘(balance+null-excess+alignment)+ gap 天花板 + coarse/fine + REAL_DIFFUSE",
        ha="center",fontsize=11,fontweight="bold")
fig.savefig(f"{A}/figs_dendro/method_flowchart.png",dpi=120,bbox_inches="tight"); print("flowchart → figs_dendro/method_flowchart.png")
