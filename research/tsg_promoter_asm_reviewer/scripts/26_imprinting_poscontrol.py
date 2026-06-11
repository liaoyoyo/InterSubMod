#!/usr/bin/env python3
"""
Positive control on REAL biology: validated ruler on known imprinted DMRs.
Imprinted loci: one germline haplotype methylated, other not -> STRONG HP1 vs HP2 cluster.
If ruler gives high blind-ARI here -> proves it detects REAL ASM (not just simulation),
and the somatic NEGATIVE is real biology. Also establishes Layer A baseline reference.
Reuses observed-only + blind-ARI from the validated framework.
"""
import numpy as np, random
import pysam
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score
from collections import Counter
random.seed(42); np.random.seed(42)
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
WINDOW = 750; ML_HIGH, ML_LOW = 200, 50; MIN_CPG = 4; MIN_SHARED = 3; MIN_GROUP = 8
bam = pysam.AlignmentFile(BAM, "rb")

# Known imprinted DMRs (hg38). Split reads by HP1 vs HP2 (germline haplotypes).
IMPRINTED = [
    ("chr11", 2018000, "H19/IGF2 ICR"),
    ("chr11", 2699900, "KvDMR1/KCNQ1OT1"),
    ("chr15", 24954900, "SNRPN/PWS-IC"),
    ("chr20", 58851000, "GNAS A/B DMR"),
    ("chr14", 100809000, "MEG3 IG-DMR"),
    ("chr19", 56840000, "PEG3 DMR"),
    ("chr7", 130490000, "MEST/PEG1"),
    ("chr6", 144006000, "PLAGL1/ZAC1"),
    ("chr20", 37521500, "NNAT DMR"),
    ("chr7", 94656000, "PEG10"),
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
    return max(adjusted_rand_score(lab,p) for p in preds), max(adjusted_mutual_info_score(lab,p) for p in preds)

def collect_hp(chrom,pos):
    var0=pos-1; s,e=var0-WINDOW,var0+WINDOW; g={'1':[],'2':[]}
    for r in bam.fetch(chrom,max(0,s),e):
        if r.flag&0x900 or r.flag&0x400: continue
        h=hp(r)
        if h not in ('1','2'): continue
        m={p:st for p,st in mcalls(r).items() if s<=p<=e and st>=0}
        if len(m)>=MIN_CPG: g[h].append(m)
    return g

print("=== IMPRINTING POSITIVE CONTROL (validated ruler on REAL germline ASM) ===")
print(f"{'DMR':<22}{'chrom:pos':<18}{'n1/n2':>9}{'blind_ARI':>10}{'sil':>7}  expect=HIGH")
aris=[]
for chrom,pos,name in IMPRINTED:
    g=collect_hp(chrom,pos)
    if len(g['1'])<MIN_GROUP or len(g['2'])<MIN_GROUP: print(f"{name:<22}{chrom+':'+str(pos):<18} insufficient ({len(g['1'])}/{len(g['2'])})"); continue
    reads=g['1']+g['2']; lab0=[0]*len(g['1'])+[1]*len(g['2'])
    D,keep=observed_only(reads)
    if D is None: print(f"{name:<22} drop too many"); continue
    lab=np.array([lab0[i] for i in keep])
    if lab.sum()<MIN_GROUP or (len(lab)-lab.sum())<MIN_GROUP: print(f"{name:<22} unbalanced after drop"); continue
    ari,ami=blind_ari(D,lab); sil=silhouette_score(D,lab,metric='precomputed')
    aris.append(ari)
    flag="STRONG ASM ✓" if ari>=0.30 else ("weak" if ari>=0.15 else "no cluster")
    print(f"{name:<22}{chrom+':'+str(pos):<18}{str(int((lab==0).sum()))+'/'+str(int((lab==1).sum())):>9}{ari:>10.3f}{sil:>7.3f}  {flag}")
print()
if aris:
    print(f"imprinted median blind-ARI = {np.median(aris):.3f} (n={len(aris)})  | high(>=0.30) = {sum(1 for a in aris if a>=0.30)}/{len(aris)}")
    print(f">>> 對照: somatic TP median ARI=0.135, het-NULL=0.177 (廣掃)")
    print(f">>> {'尺在真 imprinted ASM 偵測到強 cluster -> 證明尺有效 + somatic NEGATIVE 是真 biology' if np.median(aris)>0.30 else '若 imprinted 也低 -> 需檢討尺敏感度'}")
