#!/usr/bin/env python3
"""全資料比較: per-CpG-by-label (modkit-equiv marginal) vs ISM joint PERMANOVA / cluster。
merge percpg_records (per-CpG n_sig/axis) + decisionflow_records_tumor (joint perm/axis + byT cluster) by (chrom,pos)。
per 軸 4-cell Venn (joint-sig × percpg-sig) + 按 decisionflow state 的 per-CpG 分佈。
ISM-only = joint-sig 但 per-CpG 貧乏 (ISM 淨優勢) / modkit-only = per-CpG-sig 但 joint 不顯著。
輸出 percpg_compare_summary.json。"""
import json, os
import numpy as np
A=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AX=('hp','carrier','allele')
def bh(ps):
    p=np.asarray(ps,float); m=len(p)
    if m==0: return p
    o=np.argsort(p); r=p[o]*m/np.arange(1,m+1); q=np.minimum.accumulate(r[::-1])[::-1]
    out=np.empty(m); out[o]=np.clip(q,0,1); return out
def pct(a,b): return round(100*a/b,2) if b else 0.0
pc={(r['chrom'],r['pos']):r for r in json.load(open(f"{A}/percpg_records.json"))}
df={(r['chrom'],r['pos']):r for r in json.load(open(f"{A}/decisionflow_records_tumor.json"))}
keys=[k for k in df if k in pc]
print(f"merge: percpg={len(pc)} decisionflow={len(df)} 共同={len(keys)}")
# joint BH-FDR across loci per axis (用 decisionflow perm[axis].p)
jointq={ax:{} for ax in AX}
for ax in AX:
    idx=[k for k in keys if df[k]['perm'].get(ax)]
    ps=[df[k]['perm'][ax]['p'] for k in idx]; qs=bh(ps)
    for k,q in zip(idx,qs): jointq[ax][k]=q
def joint_sig(k,ax): return k in jointq[ax] and jointq[ax][k]<0.05
def percpg_sig(k,ax,K):
    a=pc[k].get('axes')
    if not a or not a.get(ax): return None
    return a[ax]['n_sig']>=K
res={'n_merged':len(keys)}
# per-axis 4-cell Venn (K=1 與 K=3)
for K in (1,3):
    res[f'venn_K{K}']={}
    for ax in AX:
        both=ism_only=mod_only=neither=ntest=0
        for k in keys:
            js=joint_sig(k,ax); ps=percpg_sig(k,ax,K)
            if ps is None or k not in jointq[ax]: continue  # 兩邊都可測才比
            ntest+=1
            if js and ps: both+=1
            elif js and not ps: ism_only+=1
            elif ps and not js: mod_only+=1
            else: neither+=1
        res[f'venn_K{K}'][ax]=dict(n_test=ntest,both=both,ism_only=ism_only,mod_only=mod_only,neither=neither,
            ism_only_pct=pct(ism_only,ntest),mod_only_pct=pct(mod_only,ntest))
# overall any-axis (K=1): ISM 偵測(joint any-sig) × modkit 偵測(percpg any-sig)
for K in (1,3):
    both=ism_only=mod_only=neither=0
    for k in keys:
        ja=any(joint_sig(k,ax) for ax in AX if k in jointq[ax])
        pa=any((percpg_sig(k,ax,K) or False) for ax in AX)
        if ja and pa: both+=1
        elif ja and not pa: ism_only+=1
        elif pa and not ja: mod_only+=1
        else: neither+=1
    res[f'overall_K{K}']=dict(both=both,ism_only=ism_only,mod_only=mod_only,neither=neither,
        ism_only_pct=pct(ism_only,len(keys)),mod_only_pct=pct(mod_only,len(keys)))
# 按 decisionflow state 的 per-CpG n_sig 分佈 (max over axes)
def state(k):
    r=df[k]; nc=r.get('n_complete')
    if nc is None or nc<6: return 'S1'
    bt=r['byT'].get('4',{})
    if bt.get('split'): return 'S5' if bt.get('aligned') else 'S4'
    sig=[ax for ax in AX if df[k]['perm'].get(ax) and jointq[ax].get(k,1)<0.05]
    return 'S3' if sig else 'S2'
bystate={}
for k in keys:
    s=state(k); a=pc[k].get('axes')
    mx=max((a[ax]['n_sig'] for ax in AX if a and a.get(ax)),default=0)
    bystate.setdefault(s,[]).append(mx)
res['by_state']={s:dict(n=len(v),percpg_zero=int(sum(1 for x in v if x==0)),
    percpg_zero_pct=pct(sum(1 for x in v if x==0),len(v)),
    median_nsig=float(np.median(v)),mean_nsig=round(float(np.mean(v)),2)) for s,v in bystate.items()}
json.dump(res,open(f"{A}/percpg_compare_summary.json","w"),indent=1)
print("\n=== per-axis Venn (K=1: ≥1 sig CpG = modkit 偵測到) ===")
for ax in AX:
    v=res['venn_K1'][ax]; print(f"  {ax:8s} (n={v['n_test']}): both={v['both']} ISM-only={v['ism_only']}({v['ism_only_pct']}%) modkit-only={v['mod_only']}({v['mod_only_pct']}%) neither={v['neither']}")
print(f"\n=== overall any-axis K=1 ===\n  both={res['overall_K1']['both']} ISM-only={res['overall_K1']['ism_only']}({res['overall_K1']['ism_only_pct']}%) modkit-only={res['overall_K1']['mod_only']}({res['overall_K1']['mod_only_pct']}%) neither={res['overall_K1']['neither']}")
print("\n=== per-CpG n_sig by decisionflow state (max over axes) ===")
for s in ('S5','S4','S3','S2','S1'):
    if s in res['by_state']: d=res['by_state'][s]; print(f"  {s}: n={d['n']} median_nsig={d['median_nsig']} percpg=0: {d['percpg_zero']}({d['percpg_zero_pct']}%)")
print("\nWROTE percpg_compare_summary.json")
