#!/usr/bin/env python3
"""
Imprinting positive control on NORMAL (intact imprinting) vs TUMOR (cancer CNV/LOH/LOI).
Clean positive control: NORMAL imprinted DMRs should show STRONG HP1 vs HP2 cluster
(proves ruler detects REAL ASM). TUMOR weaker = biologically meaningful (LOI/CNV).
Also runs germline-het sites in normal as Layer-A reference.
"""
import numpy as np, random, sys
import pysam
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score
from collections import Counter
random.seed(42); np.random.seed(42)
TUMOR = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
NORMAL = "/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam"
WINDOW = 750; ML_HIGH, ML_LOW = 200, 50; MIN_CPG = 4; MIN_SHARED = 3; MIN_GROUP = 8
# imprinted DMRs hg38 (refined)
IMPRINTED = [
    ("chr11", 2001800, "H19/IGF2 IG-DMR"), ("chr11", 2700000, "KvDMR1"),
    ("chr15", 24955000, "SNRPN/PWS-IC"), ("chr20", 58851000, "GNAS-XL"),
    ("chr20", 58889000, "GNAS A/B"), ("chr14", 100809000, "MEG3 IG-DMR"),
    ("chr19", 56840000, "PEG3"), ("chr7", 130487000, "MEST/PEG1"),
    ("chr6", 144006000, "PLAGL1"), ("chr7", 94656000, "PEG10"),
    ("chr6", 3849000, "FAM50B"), ("chr13", 48319000, "RB1 DMR"),
]
def hp(r):
    for t,v in r.tags:
        if t=="HP": return str(v)
    return None
def mcalls(r):
    o={}
    try: mod=r.modified_bases
    except: return o
    if not mod: return o
    r2={a:b for a,b in r.get_aligned_pairs(matches_only=False) if a is not None and b is not None}
    for k,calls in mod.items():
        if k[2]!='m': continue
        for rp,ml in calls:
            rf=r2.get(rp)
            if rf is not None: o[rf]=1 if ml>=ML_HIGH else (0 if ml<=ML_LOW else -1)
    return o
def observed_only(reads):
    n=len(reads); cov=Counter()
    for r in reads:
        for c in r: cov[c]+=1
    core={c for c,k in cov.items() if k>=max(3,0.25*n)}
    reads=[{c:r[c] for c in r if c in core} for r in reads]
    keep=[i for i,r in enumerate(reads) if len(r)>=MIN_CPG]; reads=[reads[i] for i in keep]; n=len(reads)
    if n<2*MIN_GROUP: return None,keep
    for _ in range(n):
        if n<=2*MIN_GROUP: break
        bad=np.array([sum(1 for j in range(n) if i!=j and len(set(reads[i])&set(reads[j]))<MIN_SHARED)/(n-1) for i in range(n)])
        if bad.max()<=0.45: break
        d=int(bad.argmax()); reads.pop(d); keep.pop(d); n-=1
    D=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            sh=set(reads[i])&set(reads[j])
            D[i,j]=D[j,i]=np.mean([reads[i][c]!=reads[j][c] for c in sh]) if len(sh)>=MIN_SHARED else np.nan
    if np.isnan(D).any(): D[np.isnan(D)]=np.nanmax(D) if not np.isnan(D).all() else 0.5
    return D,keep
def blind_ari(D,lab):
    preds=[AgglomerativeClustering(n_clusters=2,metric='precomputed',linkage=l).fit_predict(D) for l in ('average','complete')]
    try: preds.append(SpectralClustering(n_clusters=2,affinity='precomputed',random_state=42).fit_predict(np.exp(-D/(D.std()+1e-9))))
    except: pass
    return max(adjusted_rand_score(lab,p) for p in preds)
def run_locus(bamf, chrom, pos):
    bam=pysam.AlignmentFile(bamf,"rb")
    var0=pos-1; s,e=var0-WINDOW,var0+WINDOW; g={'1':[],'2':[]}
    for r in bam.fetch(chrom,max(0,s),e):
        if r.flag&0x900 or r.flag&0x400: continue
        h=hp(r)
        if h not in ('1','2'): continue
        m={p:st for p,st in mcalls(r).items() if s<=p<=e and st>=0}
        if len(m)>=MIN_CPG: g[h].append(m)
    bam.close()
    n1,n2=len(g['1']),len(g['2'])
    if n1<MIN_GROUP or n2<MIN_GROUP: return None,(n1,n2)
    reads=g['1']+g['2']; lab0=[0]*n1+[1]*n2
    D,keep=observed_only(reads)
    if D is None: return None,(n1,n2)
    lab=np.array([lab0[i] for i in keep])
    if lab.sum()<MIN_GROUP or (len(lab)-lab.sum())<MIN_GROUP: return None,(n1,n2)
    return blind_ari(D,lab),(int((lab==0).sum()),int((lab==1).sum()))

print("=== IMPRINTING: NORMAL (intact) vs TUMOR (HCC1395 cancer) ===")
print(f"{'DMR':<18}{'NORMAL ARI':>14}{'TUMOR ARI':>14}   interpretation")
nor_ari, tum_ari = [], []
for chrom,pos,name in IMPRINTED:
    na,nn=run_locus(NORMAL,chrom,pos); ta,tn=run_locus(TUMOR,chrom,pos)
    ns=f"{na:.3f}({nn[0]}/{nn[1]})" if na is not None else f"NA({nn[0]}/{nn[1]})"
    ts=f"{ta:.3f}({tn[0]}/{tn[1]})" if ta is not None else f"NA({tn[0]}/{tn[1]})"
    if na is not None: nor_ari.append(na)
    if ta is not None: tum_ari.append(ta)
    interp=""
    if na is not None and ta is not None:
        if na>=0.30 and ta<0.15: interp="LOI in tumor (normal ASM lost)"
        elif na>=0.30 and ta>=0.30: interp="imprinting preserved"
        elif na<0.15: interp="weak/coord-off"
    print(f"{name:<18}{ns:>14}{ts:>14}   {interp}")
print()
if nor_ari:
    print(f"NORMAL imprinted: median blind-ARI={np.median(nor_ari):.3f} (n={len(nor_ari)}), high(>=0.3)={sum(1 for a in nor_ari if a>=0.3)}/{len(nor_ari)}")
if tum_ari:
    print(f"TUMOR  imprinted: median blind-ARI={np.median(tum_ari):.3f} (n={len(tum_ari)})")
print(">>> 對照 somatic TP=0.135 / het-NULL=0.177 (廣掃)")
print(">>> 若 NORMAL imprinted median >> somatic -> 尺證明有效 (偵測真 ASM); tumor < normal -> LOI/CNV")
