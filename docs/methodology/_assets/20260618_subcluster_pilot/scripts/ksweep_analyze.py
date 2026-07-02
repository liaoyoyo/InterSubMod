#!/usr/bin/env python3
"""k-sweep 後分析: 對齊最佳 k vs silhouette 最佳 k + 顯著性 gate + 雙層多重比較校正。
讀 ksweep_wg_records.json → 量化「多少 loci 的 silhouette-k 其實不是最對齊的 k」。
校正: (1) within-locus Bonferroni × K_test(該軸掃的 k 數); (2) genome-wide BH-FDR over comparable loci。
gate: 改善 ΔV>DELTA + Cochran 可靠(min expected cell>=5)。便宜可重算。"""
import json, sys
import numpy as np
ASSET="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
DELTA=float(sys.argv[1]) if len(sys.argv)>1 else 0.10
FDR=float(sys.argv[2]) if len(sys.argv)>2 else 0.05
INP=sys.argv[3] if len(sys.argv)>3 else f"{ASSET}/ksweep_wg_records.json"
OUTP=sys.argv[4] if len(sys.argv)>4 else f"{ASSET}/ksweep_wg_summary.json"
recs=json.load(open(INP))

def bh_reject(pvals, q=0.05):
    """BH-FDR; return boolean reject array aligned to input order + threshold."""
    m=len(pvals);
    if m==0: return [], None
    order=sorted(range(m), key=lambda i:pvals[i])
    thr=0.0; cut=-1
    for rank,i in enumerate(order,1):
        if pvals[i] <= q*rank/m: cut=rank; thr=q*rank/m
    rej=[False]*m
    if cut>0:
        for rank,i in enumerate(order,1):
            if rank<=cut: rej[i]=True
    return rej, thr

def analyze(axis):
    Vk=f"V_{axis}"; pk=f"p_{axis}"; ek=f"e_{axis}"
    n_eval=n_comp=n_eq=n_ne=n_ne_hi=n_ne_lo=0
    p_locus=[]; meta=[]; examples=[]
    for r in recs:
        per=r["per_k"]
        valid=[d for d in per if d[Vk] is not None]
        if not valid: continue
        n_eval+=1
        sil_opt=max(per,key=lambda d:d["sil"])
        V_sil=sil_opt[Vk]
        align=max(valid,key=lambda d:d[Vk])
        K_test=len(valid)
        if V_sil is None:   # silhouette 的 k 該軸無法成表 → 不列入可比較分母
            p_locus.append(1.0); meta.append(None); continue
        n_comp+=1
        if align["k"]==sil_opt["k"]: n_eq+=1
        else:
            n_ne+=1
            if align["k"]>sil_opt["k"]: n_ne_hi+=1
            else: n_ne_lo+=1
        # improving k (≠sil_k, ΔV>DELTA, Cochran reliable)
        imp=[d for d in valid if d["k"]!=sil_opt["k"] and (d[Vk]-V_sil)>DELTA and (d[ek] is not None and d[ek]>=5)]
        if imp:
            best=min(imp,key=lambda d:d[pk])
            p_bonf=min(1.0, best[pk]*K_test)
            p_locus.append(p_bonf)
            meta.append({"pos":f"{r['chrom']}:{r['pos']}","n":r["n"],"sil_k":sil_opt["k"],"V_sil":round(V_sil,3),
                "align_k":best["k"],"V_align":round(best[Vk],3),"p_raw":best[pk],"K_test":K_test,"p_bonf":p_bonf})
        else:
            p_locus.append(1.0); meta.append(None)
    rej,thr=bh_reject(p_locus,FDR)
    n_sig_bonf=sum(1 for p in p_locus if p<0.05)
    n_fdr=sum(1 for i,m in enumerate(meta) if m and rej[i])
    # examples: strongest FDR-significant
    sig_meta=[meta[i] for i in range(len(meta)) if meta[i] and rej[i]]
    sig_meta.sort(key=lambda m:m["p_bonf"])
    examples=sig_meta[:8]
    return {"axis":axis,"n_axis_evaluable":n_eval,"n_comparable":n_comp,
        "n_silK_eq_alignK":n_eq,"n_silK_ne_alignK_raw":n_ne,
        "  of_ne_higher_k":n_ne_hi,"  of_ne_lower_k":n_ne_lo,
        "n_sil_not_most_aligned_BonfSig":n_sig_bonf,
        "n_sil_not_most_aligned_FDR":n_fdr,
        "pct_FDR_of_comparable":round(100*n_fdr/n_comp,2) if n_comp else None,
        "selfcheck_eq_plus_ne":n_eq+n_ne, "selfcheck_eq_ne_eq_comparable":(n_eq+n_ne==n_comp),
        "examples":examples}

N=len(recs)
out={"sample":"HCC1395 WG (all 22 chr, tumor reads, BERNOULLI ±5000)","DELTA_V":DELTA,"FDR_q":FDR,
     "n_clusterable_regions":N,
     "HP_axis":analyze("hp"),"CARRIER_axis":analyze("carrier"),"ALLELE_axis":analyze("allele")}
json.dump(out,open(OUTP,"w"),indent=1)
print(json.dumps(out,ensure_ascii=False,indent=1))
