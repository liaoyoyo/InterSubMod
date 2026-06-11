#!/usr/bin/env python3
"""
36_modtype_readrelationship_poc.py — PoC: is the 5mC-vs-5hmC read-relationship test viable?

Locus: chr4:10739757 (TP_LOH, REF=G ALT=A) — has 5hmC signal (dbeta_5hmC=+0.022, 5mC=-0.137).
LOH -> use ALLELE-axis (ALT=somatic vs REF=germline).

Tests:
 (A) 5hmC abundance — is it even dense enough to cluster?
 (B) per-CpG Δ5mC vs Δ5hmC correlation (ALT−REF) — conversion signature (5mC↓ where 5hmC↑?)
 (C) cluster-ARI: cluster reads by 5mC matrix vs by 5hmC matrix — same read partition?
 (D) per-read mean-5mC vs mean-5hmC correlation (note mechanical mutual-exclusivity floor)

Window-restricted to ±1kb (FIX from 35). Output: console + figure.
"""
import numpy as np, pysam
from collections import defaultdict
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
try:
    import sys; sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod")
    from scripts.lib.plot_setup import setup_plot_style; setup_plot_style(base_size=10,dpi=150)
except Exception: pass

BAM="/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
CHROM,SPOS,REF,ALT="chr4",10739757,"G","A"
WIN=1000; THR=0.5; MIN_CPG=4
bam=pysam.AlignmentFile(BAM)

# per read: {cpg: {'m':prob,'h':prob}} (tumor reads only) + allele
reads={}; allele={}
for r in bam.fetch(CHROM,SPOS-WIN,SPOS+WIN):
    if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
    q2r=dict((q,rp) for q,rp in r.get_aligned_pairs(matches_only=True))
    r2q={v:k for k,v in q2r.items()}
    al='na'
    if (SPOS-1) in r2q:
        b=r.query_sequence[r2q[SPOS-1]].upper(); al='alt' if b==ALT else ('ref' if b==REF else 'oth')
    mb=r.modified_bases
    if not mb: continue
    rid=r.query_name; d=reads.setdefault(rid,{})
    for (canon,strand,mod),calls in mb.items():
        mc='m' if str(mod) in('m','0') else('h' if str(mod) in('h','17') else None)
        if mc is None: continue
        for qp,ql in calls:
            rp=q2r.get(qp)
            if rp is None or abs(rp-(SPOS-1))>WIN: continue
            d.setdefault(rp,{})[mc]=ql/255.0
    allele[rid]=al
reads={k:v for k,v in reads.items() if len(v)>=MIN_CPG}
print(f"reads (>= {MIN_CPG} CpG in ±{WIN}bp): {len(reads)}  allele: alt={sum(1 for a in allele.values() if a=='alt')} ref={sum(1 for a in allele.values() if a=='ref')}")

# (A) 5hmC abundance
all_m=[rec.get('m',np.nan) for d in reads.values() for rec in d.values()]
all_h=[rec.get('h',np.nan) for d in reads.values() for rec in d.values()]
frac_m=np.nanmean([1 if x>=THR else 0 for x in all_m if not np.isnan(x)])
frac_h=np.nanmean([1 if x>=THR else 0 for x in all_h if not np.isnan(x)])
print(f"\n(A) abundance: 5mC frac>={THR} = {frac_m:.3f}  5hmC frac>={THR} = {frac_h:.3f}  (5hmC/5mC = {frac_h/max(frac_m,1e-9):.2f})")

# (B) per-CpG Δ5mC vs Δ5hmC (ALT−REF)
cpgs=sorted(set(rp for d in reads.values() for rp in d))
def per_cpg_beta(mc, al_set):
    out={}
    for cpg in cpgs:
        vals=[reads[rid][cpg].get(mc) for rid in reads if allele.get(rid) in al_set and cpg in reads[rid] and mc in reads[rid][cpg]]
        if len(vals)>=3: out[cpg]=np.mean([1 if v>=THR else 0 for v in vals])
    return out
m_alt=per_cpg_beta('m',{'alt'}); m_ref=per_cpg_beta('m',{'ref'})
h_alt=per_cpg_beta('h',{'alt'}); h_ref=per_cpg_beta('h',{'ref'})
shared=[c for c in cpgs if c in m_alt and c in m_ref and c in h_alt and c in h_ref]
dm=np.array([m_alt[c]-m_ref[c] for c in shared]); dh=np.array([h_alt[c]-h_ref[c] for c in shared])
print(f"\n(B) per-CpG (ALT−REF) on {len(shared)} shared CpG: mean Δ5mC={dm.mean():.3f} Δ5hmC={dh.mean():.3f}")
if len(shared)>=5:
    rho,p=stats.spearmanr(dm,dh)
    print(f"    Spearman(Δ5mC, Δ5hmC) = {rho:.3f} (p={p:.3f})  {'← 負相關=5mC↓處5hmC↑=轉換signature' if rho<-0.2 else '← 無明顯轉換耦合'}")

# (C) cluster-ARI: 5mC matrix vs 5hmC matrix
rid_list=list(reads);
def dist_matrix(mc):
    n=len(rid_list); D=np.full((n,n),np.nan)
    for i in range(n): D[i,i]=0
    for i in range(n):
        di=reads[rid_list[i]]
        for j in range(i+1,n):
            dj=reads[rid_list[j]]
            sh=[c for c in set(di)&set(dj) if mc in di[c] and mc in dj[c]]
            if len(sh)>=MIN_CPG:
                D[i,j]=D[j,i]=np.mean([abs(di[c][mc]-dj[c][mc]) for c in sh])
    return D
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score
def cluster(D):
    Df=np.where(np.isnan(D),np.nanmax(D[~np.isnan(D)]),D)
    return AgglomerativeClustering(n_clusters=2,metric="precomputed",linkage="average").fit_predict(Df)
Dm=dist_matrix('m'); Dh=dist_matrix('h')
lab_m=cluster(Dm); lab_h=cluster(Dh)
ari=adjusted_rand_score(lab_m,lab_h)
# also each vs allele truth
al_truth=[0 if allele.get(rid)=='ref' else 1 for rid in rid_list]
ari_m_allele=adjusted_rand_score(al_truth,lab_m); ari_h_allele=adjusted_rand_score(al_truth,lab_h)
print(f"\n(C) cluster-ARI(5mC-labels, 5hmC-labels) = {ari:.3f}  {'← 同一批read分群=耦合' if ari>0.3 else '← 不同read分群=獨立/或5hmC太稀疏'}")
print(f"    5mC-cluster vs allele ARI={ari_m_allele:.3f}  |  5hmC-cluster vs allele ARI={ari_h_allele:.3f}")

# (D) per-read mean 5mC vs mean 5hmC
rm=[]; rh=[]
for rid in rid_list:
    ms=[rec['m'] for rec in reads[rid].values() if 'm' in rec]; hs=[rec['h'] for rec in reads[rid].values() if 'h' in rec]
    if ms and hs:
        rm.append(np.mean([1 if x>=THR else 0 for x in ms])); rh.append(np.mean([1 if x>=THR else 0 for x in hs]))
rm=np.array(rm); rh=np.array(rh)
rho2,p2=stats.spearmanr(rm,rh) if len(rm)>=5 else (np.nan,np.nan)
print(f"\n(D) per-read Spearman(mean5mC, mean5hmC) = {rho2:.3f} (p={p2:.3f}) over {len(rm)} reads")
print(f"    ⚠ caveat: 同一C只能是5mC或5hmC其一 → 機械上有負向floor；需與per-CpG(B)合看")

# figure
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(11,4.6))
if len(shared)>=5:
    ax1.scatter(dm,dh,s=40,c="#D97757",alpha=0.7,edgecolors="k",linewidths=0.4)
    ax1.axhline(0,color="#999",lw=0.8); ax1.axvline(0,color="#999",lw=0.8)
    ax1.set_xlabel("Δ5mC (ALT−REF) per CpG"); ax1.set_ylabel("Δ5hmC (ALT−REF) per CpG")
    ax1.set_title(f"(B) per-CpG Δ5mC vs Δ5hmC  Spearman={rho:.2f}\n負相關=去甲基化轉換 signature")
ax2.scatter(rm,rh,s=40,c="#1E3A8A",alpha=0.6,edgecolors="k",linewidths=0.4)
ax2.set_xlabel("per-read mean 5mC"); ax2.set_ylabel("per-read mean 5hmC")
ax2.set_title(f"(D) per-read 5mC vs 5hmC  Spearman={rho2:.2f}\n(⚠機械互斥floor)  cluster-ARI={ari:.2f}")
plt.tight_layout(); plt.savefig("figures/modtype_readrel_poc_chr4.png",dpi=150,bbox_inches="tight")
print("\nwrote figures/modtype_readrel_poc_chr4.png")
print(f"\n=== 方法可行性 verdict ===")
print(f"5hmC 豐度 {frac_h:.3f} | per-CpG shared {len(shared)} | per-read {len(rm)} | cluster-ARI {ari:.3f}")
