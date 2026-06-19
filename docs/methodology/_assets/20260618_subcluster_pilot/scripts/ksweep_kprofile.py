#!/usr/bin/env python3
"""per-locus k-profile: 切法品質(夠好 vs 次好) + 唯一性 + 多解析度分類。
讀 ksweep_wg_records.json(tumor) + ksweep_records_merged.json(merged)。零 compute。
margin = best-k 與 2nd-best-k 分數差(silhouette + alignment-V 兩種)。
class: single-k-forced(len=1,被迫k=2) / [k-choice 子集:] multi-resolution / confident-unique / ambiguous-near-tie。
meaningful k gate: 任一軸 V>=0.3 & Cochran e>=5 & chi2 p<0.05。
輸出 kprofile_summary.json + kprofile_loci_{tumor,merged}.json。"""
import json
import numpy as np
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
VTH=0.30; ETH=5.0; PTH=0.05; UNIQ=0.15; AMB=0.05
AXES=("hp","carrier","allele")

def meaningful(d):
    """此 k 是否顯著對齊任一 a-priori 軸。回 (bool, best_axis, best_V)。"""
    best=None
    for ax in AXES:
        V,p,e=d[f"V_{ax}"],d[f"p_{ax}"],d[f"e_{ax}"]
        if V is not None and V>=VTH and e is not None and e>=ETH and p is not None and p<PTH:
            if best is None or V>best[2]: best=(True,ax,V)
    return best if best else (False,None,None)

def primary_axis(pk):
    """全 k 中 max-V 的軸。"""
    bestV=-1; bestax=None
    for d in pk:
        for ax in AXES:
            V=d[f"V_{ax}"]
            if V is not None and V>bestV: bestV=V; bestax=ax
    return bestax

def peaks(vals):
    if len(vals)<3: return 0
    p=sum(1 for i in range(1,len(vals)-1) if vals[i]>vals[i-1] and vals[i]>vals[i+1])
    p+=(1 if vals[0]>vals[1] else 0)+(1 if vals[-1]>vals[-2] else 0)
    return p

def classify(r):
    pk=r["per_k"]; nk=len(pk)
    sils=[d["sil"] for d in pk]
    # meaningful k
    mk=[]; mk_axes={}
    for d in pk:
        ok,ax,V=meaningful(d)
        if ok: mk.append(d["k"]); mk_axes[d["k"]]=(ax,round(V,3))
    n_mean=len(mk)
    # margins (需 >=2 k)
    sil_margin=align_margin=None
    pax=primary_axis(pk)
    if nk>=2:
        ss=sorted(sils,reverse=True); sil_margin=round(ss[0]-ss[1],3)
        if pax is not None:
            vs=sorted([d[f"V_{pax}"] for d in pk if d[f"V_{pax}"] is not None],reverse=True)
            if len(vs)>=2: align_margin=round(vs[0]-vs[1],3)
    best_margin=max([m for m in (sil_margin,align_margin) if m is not None],default=None)
    npk=peaks(sils)
    # classify
    if nk==1:
        cls="single-k-forced"; sub="aligned" if n_mean>=1 else "no-axis"
    elif n_mean>=2:
        cls="multi-resolution"; sub=f"{n_mean}meaningful-k"
    elif n_mean==1 and best_margin is not None and best_margin>=UNIQ:
        cls="confident-unique"; sub="margin>=0.15"
    else:
        cls="ambiguous-near-tie"
        sub="no-axis" if n_mean==0 else ("near-tie" if (best_margin is not None and best_margin<AMB) else "weak")
    return dict(chrom=r["chrom"],pos=r["pos"],n=r["n"],nk=nk,cls=cls,sub=sub,
        best_k=(max(pk,key=lambda d:d["sil"])["k"]),meaningful_ks=mk,mk_axes=mk_axes,
        sil_margin=sil_margin,align_margin=align_margin,peaks=npk,primary_axis=primary_axis(pk))

def run(path):
    recs=json.load(open(path)); loci=[classify(r) for r in recs]; N=len(loci)
    from collections import Counter
    cc=Counter(l["cls"] for l in loci)
    # margin 分佈 (>=2 k)
    sm=[l["sil_margin"] for l in loci if l["sil_margin"] is not None]
    am=[l["align_margin"] for l in loci if l["align_margin"] is not None]
    def band(xs):
        x=np.array(xs); n=len(x)
        return dict(n=n,mean=round(float(x.mean()),3),median=round(float(np.median(x)),3),
            uniq=int((x>=UNIQ).sum()),mid=int(((x>=AMB)&(x<UNIQ)).sum()),amb=int((x<AMB).sum()))
    # k-choice 子集(nk>=2) 的三態
    choice=[l for l in loci if l["nk"]>=2]
    cc_choice=Counter(l["cls"] for l in choice)
    summ=dict(N=N,classes=dict(cc),
        single_k_forced=cc["single-k-forced"],
        k_choice_n=len(choice),
        three_state={k:cc_choice.get(k,0) for k in ("multi-resolution","confident-unique","ambiguous-near-tie")},
        sil_margin=band(sm),align_margin=band(am),
        align_gt_sil=int(sum(1 for l in loci if l["sil_margin"] is not None and l["align_margin"] is not None and l["align_margin"]>l["sil_margin"])),
        both_margin_n=int(sum(1 for l in loci if l["sil_margin"] is not None and l["align_margin"] is not None)))
    return summ,loci

st_t,loci_t=run(f"{A}/ksweep_wg_records.json")
st_m,loci_m=run(f"{A}/ksweep_records_merged.json")
json.dump(loci_t,open(f"{A}/kprofile_loci_tumor.json","w"))
json.dump(loci_m,open(f"{A}/kprofile_loci_merged.json","w"))
json.dump({"tumor":st_t,"merged":st_m,"gate":dict(V=VTH,e=ETH,p=PTH,uniq=UNIQ,amb=AMB)},
          open(f"{A}/kprofile_summary.json","w"),indent=1)
# print
for nm,st in (("TUMOR",st_t),("MERGED",st_m)):
    print(f"\n===== {nm} (可分群 N={st['N']}) =====")
    print(f"  single-k-forced(len=1,被迫k=2): {st['single_k_forced']} ({100*st['single_k_forced']/st['N']:.1f}%)")
    print(f"  --- k-choice 子集(>=2 k) n={st['k_choice_n']} 的三態 ---")
    for k,v in st['three_state'].items():
        print(f"    {k:<20}: {v} ({100*v/st['k_choice_n']:.1f}% of choice)")
    s=st['sil_margin']; a=st['align_margin']
    print(f"  sil-margin   : mean={s['mean']} 唯一(>=.15)={s['uniq']} 中={s['mid']} 模糊(<.05)={s['amb']} (n={s['n']})")
    print(f"  align-margin : mean={a['mean']} 唯一(>=.15)={a['uniq']} 中={a['mid']} 模糊(<.05)={a['amb']} (n={a['n']})")
    print(f"  align-margin > sil-margin: {st['align_gt_sil']}/{st['both_margin_n']} ({100*st['align_gt_sil']/st['both_margin_n']:.1f}%)")
    print(f"  sum-check: single + choice = {st['single_k_forced']+st['k_choice_n']} == N {st['N']} {'OK' if st['single_k_forced']+st['k_choice_n']==st['N'] else 'FAIL'}")
