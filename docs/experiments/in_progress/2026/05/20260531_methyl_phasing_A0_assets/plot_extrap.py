import json, glob, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"]=["Noto Sans CJK TC","Droid Sans Fallback","DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"]=False
AD=os.path.dirname(os.path.abspath(__file__))
allb=[]; perchr=[]
for fp in sorted(glob.glob(AD+"/extrapolation_chr*.json")):
    d=json.load(open(fp))
    if d['n_blocks_tested']==0: continue
    perchr.append((d['chrom'],d['accuracy_median'],d['null_acc_median']))
    allb.extend(d['blocks'])
accs=[b['accuracy'] for b in allb]
nulls=[b['null_acc_mean'] for b in allb if b['null_acc_mean']]

# 圖1: accuracy vs null 分布 (直方圖對照)
fig,ax=plt.subplots(figsize=(7.5,4.5))
ax.hist(accs,bins=20,range=(0.3,1),alpha=0.65,color="#15803D",label=f"救援 accuracy (median {np.median(accs):.3f})")
ax.hist(nulls,bins=20,range=(0.3,1),alpha=0.55,color="#9CA3AF",label=f"shuffle null (median {np.median(nulls):.3f})")
ax.axvline(0.5,ls="--",c="gray",lw=1)
ax.set_xlabel("救援 accuracy (甲基預測 unphase read 真 HP)");ax.set_ylabel("窗數")
ax.set_title("圖1: 外推驗證 — 甲基救援 vs 隨機 (183 窗, 20 染色體)\n救援 median 0.885 遠超 null 0.524 → 甲基能救 unphase (88.5%正確)")
ax.legend(fontsize=9)
fig.tight_layout();fig.savefig(AD+"/fig_extrap_hist.png",dpi=130);plt.close()

# 圖2: per-chr accuracy (bar, 排序)
fig,ax=plt.subplots(figsize=(9,4.5))
perchr_s=sorted(perchr,key=lambda x:-x[1])
chrs=[p[0] for p in perchr_s]; ac=[p[1] for p in perchr_s]; nu=[p[2] for p in perchr_s]
x=np.arange(len(chrs))
ax.bar(x-0.2,ac,0.4,color="#15803D",alpha=0.7,label="救援 accuracy")
ax.bar(x+0.2,nu,0.4,color="#9CA3AF",alpha=0.7,label="shuffle null")
ax.axhline(0.5,ls="--",c="gray",lw=1)
ax.set_xticks(x);ax.set_xticklabels(chrs,rotation=45,ha="right",fontsize=8)
ax.set_ylabel("accuracy median");ax.set_ylim(0,1.05)
ax.set_title("圖2: 各染色體甲基救援正確率 (vs null)\n多數 0.76-0.97；chr3(n=1)/chr4 偏低")
ax.legend(fontsize=9)
fig.tight_layout();fig.savefig(AD+"/fig_extrap_perchr.png",dpi=130);plt.close()
print("2 figs:")
for f in ["fig_extrap_hist.png","fig_extrap_perchr.png"]:
    p=AD+"/"+f;print(f"  {f}: {'OK' if os.path.exists(p) else 'MISS'} {os.path.getsize(p)//1024 if os.path.exists(p) else 0}KB")
