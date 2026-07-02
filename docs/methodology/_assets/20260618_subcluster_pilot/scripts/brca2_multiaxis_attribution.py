import csv, json
import numpy as np
from scipy.stats import fisher_exact
RD="output/_brca2_region/mini/chr13/chr13_32315128/chr13_32310128_32320128"
VAR=32315128
def bh(ps):
    p=np.array(ps); m=len(p)
    if m==0: return p
    o=np.argsort(p); r=p[o]*m/np.arange(1,m+1); q=np.minimum.accumulate(r[::-1])[::-1]
    out=np.empty(m); out[o]=np.clip(q,0,1); return out
reads={r['read_id']:r for r in csv.DictReader(open(f"{RD}/reads/reads.tsv"),delimiter='\t')}
# methylation.csv: read × CpG (0/1/NaN)
lines=open(f"{RD}/methylation/methylation.csv").read().strip().split('\n')
cpgs=[int(x) for x in lines[0].split(',')[1:]]
M={}; 
for ln in lines[1:]:
    q=ln.split(','); rid=q[0]
    M[rid]=[np.nan if v in('','NA','nan','NaN') else float(v) for v in q[1:]]
# 只用 tumor read (somatic 軸)
tum=[rid for rid in M if rid in reads and reads[rid]['is_tumor'] in('1','true','True')]
def lab(rid,axis):
    h=reads[rid]['hp']; a=reads[rid]['alt_support']
    if axis=='hp': return 'HP1fam' if h in('1','1-1') else 'HP2fam'   # 物理單倍型
    if axis=='carrier': return 'germ' if h in('1','2') else 'somatic'  # germline vs somatic-tag
    if axis=='allele': return a if a in('REF','ALT') else None
AX=['hp','carrier','allele']
# per-CpG per-axis Fisher
percpg={ax:[] for ax in AX}
for ci,cp in enumerate(cpgs):
    for ax in AX:
        g={}
        for rid in tum:
            l=lab(rid,ax)
            if l is None: continue
            v=M[rid][ci]
            if np.isnan(v): continue
            g.setdefault(l,[]).append(1 if v>=0.5 else 0)
        ks=sorted(g)
        if len(ks)<2 or any(len(g[k])<5 for k in ks): 
            percpg[ax].append(None); continue
        a0=g[ks[0]]; a1=g[ks[1]]
        m0=sum(a0); m1=sum(a1)
        tab=[[m1,len(a1)-m1],[m0,len(a0)-m0]]
        try: _,p=fisher_exact(tab)
        except: percpg[ax].append(None); continue
        db=np.mean(a1)-np.mean(a0)  # ks[1]-ks[0]
        percpg[ax].append({'pos':cp,'dist':cp-VAR,'p':p,'dbeta':db,'g0':ks[0],'g1':ks[1],'n0':len(a0),'n1':len(a1)})
# BH per 軸
sigmask={}
for ax in AX:
    idx=[i for i,r in enumerate(percpg[ax]) if r]
    ps=[percpg[ax][i]['p'] for i in idx]; qs=bh(ps)
    for j,i in enumerate(idx): percpg[ax][i]['q']=qs[j]
    sigmask[ax]=set(percpg[ax][i]['pos'] for k,i in enumerate(idx) if qs[k]<0.05)
print("=== BRCA2 多軸 per-CpG 歸屬 (tumor reads, BH-FDR q<0.05) ===")
for ax in AX:
    n=len([r for r in percpg[ax] if r]); ns=len(sigmask[ax])
    print(f"  {ax:8s}: 可測 {n:3d} CpG, 顯著 {ns:3d}")
# 多軸交集
allsig=set().union(*sigmask.values())
print(f"\n各 CpG 對齊軸組合 (顯著聯集 {len(allsig)} CpG):")
from collections import Counter
combo=Counter()
for cp in allsig:
    tag=tuple(ax for ax in AX if cp in sigmask[ax])
    combo['+'.join(tag)]+=1
for k,v in sorted(combo.items(),key=lambda x:-x[1]): print(f"  {k:25s}: {v}")
# carrier × allele 共線診斷
print("\n=== carrier × allele 共線診斷 (tumor reads) ===")
ct=Counter((lab(r,'carrier'),lab(r,'allele')) for r in tum if lab(r,'allele'))
print("  cross-tab (carrier,allele):",dict(ct))
json.dump({'cpgs':cpgs,'percpg':percpg,'sig':{k:sorted(v) for k,v in sigmask.items()}},open('/tmp/_brca2_multiaxis_out.json','w'))
