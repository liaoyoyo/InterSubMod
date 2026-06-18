#!/usr/bin/env python3
"""§8 數據: paired(pipeline Fisher/CramersV) vs tumor-only 對照 + CramersV 分布圖。鎖 section8.json。"""
import json, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
ASSET="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
SUM="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_wg_bdcprime_verify/significance_summary.csv"
def fl(x):
    try: return float(x)
    except: return None
def b(x): return str(x).strip().lower() in ("true","1","yes")
def iv(x):
    try: return int(float(x))
    except: return 0
rows=list(csv.DictReader(open(SUM))); N=len(rows)

# paired
pal_p=[fl(r["GlobalP"]) for r in rows]; pal_v=[fl(r["CramersV"]) for r in rows]
phpf_p=[fl(r["GlobalP_HPFine"]) for r in rows]; phpf_v=[fl(r["CramersV_HPFine"]) for r in rows]
# tumor-only from contingency
recs=json.load(open(f"{ASSET}/records_wg2.json"))
to_p=[]; to_v=[]
for r in recs:
    cont=r["all"]["contingency"]
    if not cont: to_p.append(None); to_v.append(None); continue
    labs=sorted(cont); clus=sorted({c for d in cont.values() for c in d})
    M=np.array([[cont[L].get(c,0) for c in clus] for L in labs],dtype=float)
    if M.shape[0]<2 or M.shape[1]<2 or M.sum()<6: to_p.append(None); to_v.append(None); continue
    try:
        chi2,p,_,_=chi2_contingency(M); v=np.sqrt(chi2/(M.sum()*(min(M.shape)-1)))
        to_p.append(float(p)); to_v.append(round(float(v),3))
    except: to_p.append(None); to_v.append(None)

def sig(ps,th=0.05): return sum(1 for x in ps if x is not None and x<th)
def vbins(vs):
    v=[x for x in vs if x is not None]
    out={"zero":sum(1 for x in v if x==0)}
    for lo,hi,k in [(0.0001,0.1,"a"),(0.1,0.3,"b"),(0.3,0.5,"c"),(0.5,0.7,"d"),(0.7,1.01,"e")]:
        out[k]=sum(1 for x in v if lo<=x<hi)
    pos=[x for x in v if x>0]; out["nonzero"]=len(pos); out["med"]=round(float(np.median(pos)),3) if pos else 0
    return out

D={"N":N,
   "paired_allele":{"sig05":sig(pal_p),"sig01":sig(pal_p,0.01),"sig001":sig(pal_p,0.001),"v":vbins(pal_v)},
   "paired_hpfine":{"sig05":sig(phpf_p),"v":vbins(phpf_v)},
   "tumor_only":{"computable":sum(1 for x in to_p if x is not None),"sig05":sig(to_p),"v":vbins(to_v)},
   "passed_gating":sum(1 for r in rows if b(r["PassedGating"])),
   "percpg_tested":sum(1 for r in rows if iv(r["Fisher_N_Tested"])>0),
   "percpg_sig":sum(1 for r in rows if iv(r["Fisher_N_Sig"])>0)}
json.dump(D,open(f"{ASSET}/section8.json","w"),ensure_ascii=False,indent=1)
print("section8.json:",json.dumps({k:(v if not isinstance(v,dict) else "...") for k,v in D.items()},ensure_ascii=False))

# figure: CramersV>0 overlay paired vs tumor-only
plt.rcParams["font.sans-serif"]=["Noto Sans CJK TC","DejaVu Sans"]; plt.rcParams["axes.unicode_minus"]=False
fig,ax=plt.subplots(figsize=(6.6,3.4))
pv=[x for x in pal_v if x and x>0]; tv=[x for x in to_v if x and x>0]
bins=np.linspace(0,1,26)
ax.hist(pv,bins=bins,density=True,alpha=0.55,color="#4A77A8",label=f"PAIRED allele (n={len(pv)}, med 0.561)")
ax.hist(tv,bins=bins,density=True,alpha=0.55,color="#D97757",label=f"TUMOR-only (n={len(tv)}, med 0.750)")
ax.set_xlabel("Cramér's V (>0 者)",fontsize=9); ax.set_ylabel("density",fontsize=9)
ax.set_title("CramersV 分布: tumor-only 切得出群者效應更強(右移)\npaired 含 tumor-vs-normal 普遍但效應較分散",fontsize=9.5)
ax.legend(fontsize=8.5); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{ASSET}/figs/cramersv_paired_vs_tumor.png",dpi=100); plt.close(fig)
print("WROTE figs/cramersv_paired_vs_tumor.png")
