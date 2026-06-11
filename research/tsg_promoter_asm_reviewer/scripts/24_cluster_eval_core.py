#!/usr/bin/env python3
"""
Cluster-quality eval CORE (priority subset, framework §5 + 5 user corrections).
Decisive gates: M0 observed-only distance (NO impute) + M1 blind-ARI recovery + M8 placebo-anchor collider.
PLUS meta-validation (驗尺): PC1 simulated 2-cluster + NC1 label-shuffle + NC-collider.
"Validate the ruler FIRST, then use it" (correction #5).
"""
import numpy as np, random, sys
import pysam
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score
from skbio.stats.distance import DistanceMatrix, permanova
random.seed(42); np.random.seed(42)

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
WINDOW = 600; ML_HIGH, ML_LOW = 200, 50
MIN_CORE_CPG = 4; MIN_SHARED = 3; MIN_GROUP = 8
bam = pysam.AlignmentFile(BAM, "rb")

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

def collect(chrom,pos,groups):
    var0=pos-1; s,e=var0-WINDOW,var0+WINDOW
    g={k:[] for k in groups}; meta={k:[] for k in groups}
    for r in bam.fetch(chrom,max(0,s),e):
        if r.flag & 0x900 or r.flag & 0x400: continue
        h=hp(r)
        if h not in groups: continue
        m={p:st for p,st in mcalls(r).items() if s<=p<=e and st>=0}
        if len(m)>=MIN_CORE_CPG:
            g[h].append(m); meta[h].append((r.reference_start,r.reference_end,r.mapping_quality,
                                            r.query_length or (r.reference_end-r.reference_start)))
    return g,meta

def observed_only_distance(reads):
    """M0: Hamming over shared CpG, NO impute. Drop poorly-connected reads iteratively."""
    n=len(reads)
    # core CpGs covered by >=30% reads
    from collections import Counter
    cov=Counter()
    for r in reads:
        for c in r: cov[c]+=1
    core=[c for c,k in cov.items() if k>=max(3,0.3*n)]
    reads=[{c:r[c] for c in r if c in core} for r in reads]
    keep=[i for i,r in enumerate(reads) if len(r)>=MIN_CORE_CPG]
    reads=[reads[i] for i in keep]
    n=len(reads)
    if n<2*MIN_GROUP: return None,None,keep
    # iterative drop: read whose >40% pairs lack MIN_SHARED
    while True:
        bad=np.zeros(n)
        for i in range(n):
            inc=0
            for j in range(n):
                if i==j: continue
                if len(set(reads[i])&set(reads[j]))<MIN_SHARED: inc+=1
            bad[i]=inc/(n-1)
        if bad.max()<=0.4 or n<=2*MIN_GROUP: break
        drop=int(bad.argmax())
        reads.pop(drop); keep.pop(drop); n-=1
    D=np.zeros((n,n)); ok=True
    for i in range(n):
        for j in range(i+1,n):
            sh=set(reads[i])&set(reads[j])
            if len(sh)<MIN_SHARED: ok=False; D[i,j]=D[j,i]=np.nan
            else:
                d=np.mean([reads[i][c]!=reads[j][c] for c in sh]); D[i,j]=D[j,i]=d
    if np.isnan(D).any():  # still incomplete -> fill with max observed (conservative, not median)
        D[np.isnan(D)]=np.nanmax(D)
    return D, reads, keep

def blind_ari(D,true_labels):
    """M1: blind k=2 cluster, max ARI vs true labels."""
    preds=[]
    for link in ('average','complete'):
        preds.append(AgglomerativeClustering(n_clusters=2,metric='precomputed',linkage=link).fit_predict(D))
    try:
        aff=np.exp(-D/(D.std()+1e-9))
        preds.append(SpectralClustering(n_clusters=2,affinity='precomputed',random_state=42).fit_predict(aff))
    except Exception: pass
    aris=[adjusted_rand_score(true_labels,p) for p in preds]
    amis=[adjusted_mutual_info_score(true_labels,p) for p in preds]
    return max(aris), max(amis)

def permR2(D,labels):
    dm=DistanceMatrix(D,[str(i) for i in range(len(labels))])
    res=permanova(dm,labels,permutations=499)
    F=res['test statistic']; R2=F/(F+(len(labels)-2))
    return R2,res['p-value']

def evaluate(reads_g0, reads_g1, name):
    reads=reads_g0+reads_g1
    labels0=[0]*len(reads_g0)+[1]*len(reads_g1)
    D,kreads,keep=observed_only_distance(reads)
    if D is None: return None
    labels=np.array([labels0[i] for i in keep])
    if labels.sum()<MIN_GROUP or (len(labels)-labels.sum())<MIN_GROUP: return None
    ari,ami=blind_ari(D,labels)
    R2,p=permR2(D,labels.tolist())
    sil=silhouette_score(D,labels,metric='precomputed')
    return dict(name=name,n0=int((labels==0).sum()),n1=int((labels==1).sum()),
                ari=ari,ami=ami,R2=R2,p=p,sil=sil,n_drop=len(reads)-len(keep))

# ================= META-VALIDATION (驗尺) =================
print("="*70); print("META-VALIDATION (先驗尺再用尺)"); print("="*70)
def sim_reads(n,p_meth,ncpg=20,miss=0.4):
    out=[]
    for _ in range(n):
        r={}
        for c in range(ncpg):
            if random.random()>miss:
                r[c]=1 if random.random()<p_meth else 0
        if len(r)>=MIN_CORE_CPG: out.append(r)
    return out
# PC1: two real clusters Δβ=0.3 -> ruler MUST detect (ARI high)
pc1=evaluate(sim_reads(40,0.65),sim_reads(40,0.35),"PC1 simulated Δβ=0.3")
# NC1: one population split randomly -> ruler MUST reject (ARI~0)
one=sim_reads(80,0.5); nc1=evaluate(one[:40],one[40:],"NC1 random-split (no real cluster)")
for r in [pc1,nc1]:
    if r: print(f"  {r['name']:<30} ARI={r['ari']:.3f} AMI={r['ami']:.3f} R2={r['R2']:.3f} p={r['p']:.3f} sil={r['sil']:.3f}")
ruler_ok = pc1 and nc1 and pc1['ari']>0.5 and nc1['ari']<0.15
print(f"  >>> 尺驗證: PC1 ARI>0.5={pc1['ari']>0.5 if pc1 else '?'} AND NC1 ARI<0.15={nc1['ari']<0.15 if nc1 else '?'} -> {'VALIDATED' if ruler_ok else 'RULER BROKEN'}")

# ================= APPLY to BRCA2 + top loci (M0+M1+M8) =================
print("\n"+"="*70); print("APPLY: somatic loci HP1 vs HP1-1 (M0 observed-only + M1 blind-ARI)"); print("="*70)
LOCI=[("chr13",32315128,"BRCA2","1","1-1"),("chr20",50267392,"chr20","1","1-1"),
      ("chr12",31601630,"chr12","2","2-1"),("chr9",42376881,"chr9","2","2-1"),
      ("chr16",17774746,"chr16","2","2-1")]
print(f"{'locus':<16}{'n0/n1':>9}{'blind_ARI':>10}{'AMI':>7}{'R2':>7}{'p':>7}{'sil':>7}  verdict")
results=[]
for chrom,pos,lab,gg,gs in LOCI:
    g,meta=collect(chrom,pos,[gg,gs])
    if len(g[gg])<MIN_GROUP or len(g[gs])<MIN_GROUP: print(f"{lab:<16} insufficient"); continue
    r=evaluate(g[gg],g[gs],lab)
    if not r: print(f"{lab:<16} eval failed (observed-only drop too many)"); continue
    # M8 placebo: split HP1(germline gg) reads by LENGTH (mimic ALT-read spanning selection), blind-ARI
    gm=meta[gg]; lengths=np.array([m[3] for m in gm])
    med=np.median(lengths); pseudo=[i for i in range(len(g[gg])) if lengths[i]>=med]; rest=[i for i in range(len(g[gg])) if lengths[i]<med]
    placebo_ari=np.nan
    if len(pseudo)>=MIN_GROUP and len(rest)>=MIN_GROUP:
        pr=evaluate([g[gg][i] for i in rest],[g[gg][i] for i in pseudo],lab+"_placebo")
        placebo_ari=pr['ari'] if pr else np.nan
    verdict="Tier-A候選" if (r['ari']>=0.30 and (np.isnan(placebo_ari) or placebo_ari<0.10)) else \
            ("PLACEBO重現(collider?)" if placebo_ari>=0.30 else "blind-ARI弱")
    results.append((r,placebo_ari))
    print(f"{lab:<16}{str(r['n0'])+'/'+str(r['n1']):>9}{r['ari']:>10.3f}{r['ami']:>7.3f}{r['R2']:>7.3f}{r['p']:>7.3f}{r['sil']:>7.3f}  {verdict} | placebo_ARI={placebo_ari:.3f}")

print("\n"+"="*70)
print("關鍵讀法: blind_ARI = 盲分群能否恢復 somatic/germline split (M1, 最難造假)")
print("  ARI>=0.30 = 真有可分 cluster; ARI~0 = 需 label 才看得到(假); placebo_ARI 高 = collider artifact")
print("  R2/sil 是 supervised (餵 label), 僅參考; blind_ARI 才是主判準")
