import json,glob,os,numpy as np
import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"]=["Noto Sans CJK TC","Droid Sans Fallback","DejaVu Sans"];matplotlib.rcParams["axes.unicode_minus"]=False
AD=os.path.dirname(os.path.abspath(__file__))
perchr=[]
allw=[]
for fp in glob.glob(AD+"/admr_chr*.json"):
    d=json.load(open(fp))
    if d.get('n_windows',0)==0:continue
    perchr.append((d['chrom'],d['admr_frac_in_cnvloh'],d['bg_frac_in_cnvloh']))
    allw.extend(d['windows'])
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,4.6))
# 1 per-chr aDMR落LOH vs 背景 (證明 ≈ 背景, 無富集)
perchr.sort(key=lambda x:int(''.join(filter(str.isdigit,x[0].split('chr')[1]))))
chrs=[p[0] for p in perchr];af=[p[1] for p in perchr];bg=[p[2] for p in perchr]
x=np.arange(len(chrs))
ax1.bar(x-0.2,af,0.4,color="#C2410C",alpha=0.7,label="aDMR 落 CNV/LOH")
ax1.bar(x+0.2,bg,0.4,color="#9CA3AF",alpha=0.7,label="背景率(全窗)")
ax1.set_xticks(x);ax1.set_xticklabels(chrs,rotation=45,ha="right",fontsize=8)
ax1.set_ylabel("比例");ax1.set_ylim(0,1.1);ax1.axhline(0.79,ls="--",c="blue",lw=1,label="文獻 79%")
ax1.set_title("圖1: aDMR 落 CNV/LOH ≈ 背景率 (各染色體)\n→ 無富集(OR=1.46 p=0.38); 高比例是 HCC1395 背景使然")
ax1.legend(fontsize=8)
# 2 aDMR vs 非 maxdelta (證明 aDMR 是真訊號)
admr=[w['max_delta'] for w in allw if w['is_admr']];non=[w['max_delta'] for w in allw if not w['is_admr']]
ax2.hist(non,bins=25,range=(0,1),alpha=0.6,color="#9CA3AF",label=f"非aDMR (median {np.median(non):.2f})")
ax2.hist(admr,bins=25,range=(0,1),alpha=0.7,color="#15803D",label=f"aDMR (median {np.median(admr):.2f})")
ax2.axvline(0.25,ls="--",c="gray",lw=1,label="Δβ=0.25 門檻")
ax2.set_xlabel("窗內 max |Δβ| (HP1 vs HP2 甲基差)");ax2.set_ylabel("窗數")
ax2.set_title("圖2: aDMR 是真訊號 (maxΔ 0.80 vs 非 0.27)\n但不偏好 LOH 區(圖1)")
ax2.legend(fontsize=8)
fig.suptitle("aDMR × LOH 富集：對照文獻 79% — 數字相符但無真富集(背景 confound)",fontsize=12,y=1.04)
fig.tight_layout();fig.savefig(AD+"/fig_admr_enrichment.png",dpi=130,bbox_inches="tight");plt.close()
print("saved fig_admr_enrichment.png")
