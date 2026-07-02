import json, os, csv
import numpy as np
from collections import Counter
import random
random.seed(7); np.random.seed(7)
WT='/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra'
cands=json.load(open('cis_candidates_resolved.json'))
bykey={f"{c['chrom']}:{c['pos']}":c for c in cands}
res=json.load(open('nact_results.json'))
cand=[r for r in res if r['nact_verdict']=='candidate_subclone']
DTHR=0.2
def load(rd):
    mp=f'{rd}/methylation/methylation.csv';rt=f'{rd}/reads/reads.tsv';pg=f'{rd}/clustering/phylo_groups.tsv'
    mr=open(mp).read().strip().split('\n');cpgs=[int(c) for c in mr[0].split(',')[1:]]
    M={}
    for ln in mr[1:]:
        q=ln.split(',');M[q[0]]=np.array([np.nan if v in ('','NA','nan','NaN') else float(v) for v in q[1:]])
    meta={};rows=open(rt).read().splitlines();hdr=rows[0].split('\t')
    ix={k:hdr.index(k) for k in ('is_tumor','hp') if k in hdr}
    for r in rows[1:]:
        c=r.split('\t');meta[c[0]]={'is_tumor':c[ix['is_tumor']],'hp':c[ix['hp']]}
    phylo={}
    if os.path.exists(pg):
        for row in csv.DictReader(open(pg),delimiter='\t'):phylo[row['read_id']]={'coarse':row.get('coarse_label',''),'out':row.get('is_outlier','0')}
    return M,cpgs,meta,phylo
def groups(M,meta,phylo):
    T=[r for r in M if meta.get(r,{}).get('is_tumor')=='1']
    lab={r:phylo[r]['coarse'] for r in T if r in phylo and phylo[r]['coarse'] not in ('','outlier','other') and phylo[r].get('out','0') not in ('1','True','true')}
    leaves=sorted(set(lab.values()))
    if len(leaves)<2:return None
    Ln=2
    for L in range(2,8):
        if len(set('-'.join(x.split('-')[:L]) for x in leaves))>=2:Ln=L;break
    def pref(x):return '-'.join(x.split('-')[:Ln])
    top2=[k for k,_ in Counter(pref(v) for v in lab.values()).most_common(2)]
    if len(top2)<2:return None
    gA=[r for r in lab if pref(lab[r])==top2[0]];gB=[r for r in lab if pref(lab[r])==top2[1]]
    if len(gA)<4 or len(gB)<4:return None
    return gA,gB
# Vectorized mean-diff signature size: count CpGs with |mean(B)-mean(A)|>=DTHR & both sides >=3 nonNaN
def sig_count(Mat, idxA, idxB):
    A=Mat[idxA]; B=Mat[idxB]
    nA=np.sum(~np.isnan(A),axis=0); nB=np.sum(~np.isnan(B),axis=0)
    mA=np.nanmean(np.where(np.isnan(A),np.nan,A),axis=0)
    mB=np.nanmean(np.where(np.isnan(B),np.nan,B),axis=0)
    ok=(nA>=3)&(nB>=3)
    diff=np.abs(mB-mA)
    return int(np.sum(ok & (diff>=DTHR)))
NPERM=200
sample=cand[:20]
out=[]
for r in sample:
    key=f"{r['chrom']}:{r['pos']}";rd=bykey[key]['region_dir']
    try:M,cpgs,meta,phylo=load(rd)
    except:continue
    g=groups(M,meta,phylo)
    if not g:continue
    gA,gB=g
    allT=gA+gB
    Mat=np.array([M[x] for x in allT])  # reads x CpG
    nA=len(gA); idx=np.arange(len(allT))
    obs=sig_count(Mat, idx[:nA], idx[nA:])
    if obs==0:continue
    ge=0
    for _ in range(NPERM):
        perm=idx.copy(); np.random.shuffle(perm)
        if sig_count(Mat, perm[:nA], perm[nA:])>=obs: ge+=1
    pperm=(ge+1)/(NPERM+1)
    out.append((key,len(gA),len(gB),obs,round(pperm,3)))
print('Tumor-side label-permutation null on cluster_split (mean-diff signature size):')
print('key, n_gA, n_gB, obs_sigCpG, p_perm')
for x in out:print('  ',x)
weak=[x for x in out if x[4]>=0.05]
print('survivors where random split reaches same signature (p_perm>=0.05, split not special):', len(weak),'/',len(out))
print('p_perm: min',min(x[4] for x in out),'median',round(float(np.median([x[4] for x in out])),3),'max',max(x[4] for x in out))
