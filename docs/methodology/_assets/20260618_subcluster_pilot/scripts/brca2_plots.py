import csv, json, sys
import numpy as np
sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts")
import ism_heatmap_std as H
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"]=["Droid Sans Fallback","DejaVu Sans"]; plt.rcParams["axes.unicode_minus"]=False
RD="output/_brca2_region/mini/chr13/chr13_32315128/chr13_32310128_32320128"
AS="docs/methodology/_assets/20260618_subcluster_pilot"; FIGS=f"{AS}/figs_brca2"
VAR=32315128
mo=json.load(open('/tmp/_brca2_multiaxis_out.json'))
cpgs=mo['cpgs']; pc=mo['percpg']; sig={k:set(v) for k,v in mo['sig'].items()}
# --- track: carrier Δβ per CpG vs dist ---
fig,ax=plt.subplots(figsize=(11,3.4))
car=[r for r in pc['carrier'] if r]
xs=[r['dist'] for r in car]; ys=[r['dbeta'] for r in car]
issig=[r['pos'] in sig['carrier'] for r in car]
multi=[r['pos'] in sig['allele'] for r in car]  # 也對 allele
ax.axhline(0,color='#888',lw=0.8); ax.axvline(0,color='#c00',lw=1,ls='--')
ax.text(0,0.9,'sSNV',color='#c00',fontsize=9,ha='center')
for x,y,s,m in zip(xs,ys,issig,multi):
    if s and m: ax.scatter(x,y,c='#7c3aed',s=42,zorder=3,marker='D',label='_')  # carrier+allele 多軸
    elif s: ax.scatter(x,y,c='#dc2626',s=30,zorder=2)
    else: ax.scatter(x,y,c='#cbd5e1',s=10,zorder=1)
ax.axvspan(-600,-169,color='#fde68a',alpha=0.3,zorder=0)
ax.text(-385,-0.95,'焦點窗 −600~−169 (18 carrier-sig)',fontsize=8.5,ha='center',color='#92400e')
ax.set_xlabel('CpG 距 sSNV (bp)'); ax.set_ylabel('Δβ (somatic − germline)')
ax.set_title('BRCA2 chr13:32315128 — carrier 軸 per-CpG Δβ 定位 track（紅=carrier顯著,紫菱=兼allele多軸,負=somatic低甲基）',fontsize=10)
ax.set_ylim(-1.05,0.55)
plt.tight_layout(); plt.savefig(f"{FIGS}/track_carrier_dbeta.png",dpi=120,bbox_inches='tight'); plt.close()
# --- dual-panel: read×CpG meth + read×read dist ---
reads={r['read_id']:r for r in csv.DictReader(open(f"{RD}/reads/reads.tsv"),delimiter='\t')}
lines=open(f"{RD}/methylation/methylation.csv").read().strip().split('\n')
mids=[]; Me={}
for ln in lines[1:]:
    q=ln.split(','); mids.append(q[0]); Me[q[0]]=[np.nan if v in('','NA','nan','NaN') else float(v) for v in q[1:]]
dl=open(f"{RD}/distance/BERNOULLI/matrix.csv").read().strip().split('\n')
dids=dl[0].split(',')[1:] if False else [x.split(',')[0] for x in dl[1:]]
D={};
for ln in dl[1:]:
    q=ln.split(','); D[q[0]]=[np.nan if v in('','NA','nan','NaN') else float(v) for v in q[1:]]
didx={r:i for i,r in enumerate(dids)}
# 用所有 read，排序 carrier→allele→mean meth
def carr(r): return 0 if reads[r]['hp'] in('1','2') else 1
def alle(r): return {'REF':0,'ALT':1}.get(reads[r]['alt_support'],2)
use=[r for r in mids if r in reads and r in didx]
rowm={r:(np.nanmean(Me[r]) if np.any(np.isfinite(Me[r])) else 0) for r in use}
order=sorted(use,key=lambda r:(carr(r),alle(r),rowm[r]))
meth=np.array([Me[r] for r in order])
Dm=np.array([[D[a][didx[b]] for b in order] for a in order]); np.fill_diagonal(Dm,0); Dm=np.maximum(Dm,Dm.T)
# sidebars: carrier / HP / ALT / T-N / strand
CARC={0:'#0ea5e9',1:'#ef4444'}  # germ blue / somatic red
carrier_col=[CARC[carr(r)] for r in order]
sb=[("carrier(germ藍/som紅)",carrier_col)]
sb+=H.sidebar_specs({r:reads[r] for r in order},order,cluster_of=None,include_tn=True,include_strand=True)
title=f"BRCA2 chr13:32315128  read×CpG 甲基 + read×read 距離  (排序 carrier→allele→β; n={len(order)} reads, {len(cpgs)} CpG)"
H.mpl_dual_panel(meth,Dm,sb,cpgs,VAR,title,f"{FIGS}/dualpanel_brca2.png",n_cluster=2,dpi=95)
print("WROTE", FIGS+"/track_carrier_dbeta.png", "+", FIGS+"/dualpanel_brca2.png")
print(f"track: {len(car)} CpG plotted, carrier-sig {sum(issig)}")
