#!/usr/bin/env python3
"""判別流程 5-態分類 + 比例 + 門檻敏感度 (tumor / merged 分開)。
讀 decisionflow_records_{tumor,merged}.json → 套判別流程 → decisionflow_summary.json。

5 態:
 S1 ①不可驗證       : n_complete<6 (read 太少/去 NaN 後無法可靠聚類)
 S2 ②1群/無訊號     : 切不出 & PERMANOVA-FDR 不顯著 (a-priori 也無均值差)
 S3 ③監督可分        : 切不出 但 PERMANOVA-FDR 顯著 ("切不出≠沒訊號")
     S3_loc location-clean (∃軸 PERMANOVA-sig 且 PERMDISP dispP>=.05) / S3_disp dispersion-confounded
 S4 ④可分未對齊      : 離散切得出 但不對 HP/carrier/allele (epiallele? 雜訊?)
 S5 ⑤確認真結構      : 離散切得出 且對齊生物軸 (CramerV gate)

門檻: valid 群 size>=T, T in {3,4,5} 各算; headline=T4 (≤3 當 outlier)。
FDR: BH per axis 套在 "該 T 下切不出 + precondition-pass" 的位點池 (家庭=判 ②vs③ 的測試)。"""
import json, os, sys
import numpy as np
A=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AXES=("hp","carrier","allele")
def bh_q(pvals):
    p=np.asarray(pvals,float); m=len(p)
    if m==0: return np.array([])
    order=np.argsort(p); ranked=p[order]*m/np.arange(1,m+1)
    qmono=np.minimum.accumulate(ranked[::-1])[::-1]
    q=np.empty(m); q[order]=np.clip(qmono,0,1); return q
def classify(records,T):
    Tk=str(T); n=len(records); cls=[None]*n
    split_minor=[]; split_out=[]; split_minor_align=[]; split_minor_unalign=[]
    # pass1: S1 / split(S4/S5) / non-split(defer)
    nonsplit=[]  # idx of precond-pass, can't split at T
    for i,r in enumerate(records):
        nc=r.get("n_complete")
        if nc is None or nc<6: cls[i]="S1"; continue
        bt=r.get("byT",{}).get(Tk)
        if bt and bt.get("split"):
            mn=bt.get("minority"); split_minor.append(mn); split_out.append(bt.get("n_out",0))
            if bt.get("aligned"): cls[i]="S5"; split_minor_align.append(mn)
            else: cls[i]="S4"; split_minor_unalign.append(mn)
        else:
            nonsplit.append(i)
    # pass2: BH per axis over non-split pool → S2/S3
    qmap={ax:{} for ax in AXES}
    for ax in AXES:
        idxs=[i for i in nonsplit if records[i]["perm"].get(ax)]
        if not idxs: continue
        ps=[records[i]["perm"][ax]["p"] for i in idxs]
        qs=bh_q(ps)
        for j,i in enumerate(idxs): qmap[ax][i]=qs[j]
    s3_loc=s3_disp=s2=s2_untest=0
    for i in nonsplit:
        sig_axes=[ax for ax in AXES if i in qmap[ax] and qmap[ax][i]<0.05]
        testable=any(records[i]["perm"].get(ax) for ax in AXES)
        if sig_axes:
            loc=any((records[i]["perm"][ax].get("dispP") is not None and records[i]["perm"][ax]["dispP"]>=0.05) for ax in sig_axes)
            if loc: cls[i]="S3_loc"; s3_loc+=1
            else: cls[i]="S3_disp"; s3_disp+=1
        else:
            cls[i]="S2"; s2+=1
            if not testable: s2_untest+=1
    cnt={k:0 for k in ("S1","S2","S3_loc","S3_disp","S4","S5")}
    for c in cls: cnt[c]=cnt.get(c,0)+1
    return dict(T=T,n_total=n,counts=cnt,
        S1=cnt["S1"],S2=cnt["S2"],S3=cnt["S3_loc"]+cnt["S3_disp"],S3_loc=cnt["S3_loc"],S3_disp=cnt["S3_disp"],
        S4=cnt["S4"],S5=cnt["S5"],s2_untestable=s2_untest,
        n_split=cnt["S4"]+cnt["S5"],n_nonsplit_precond=len(nonsplit),
        split_minor=split_minor,split_out=split_out,
        split_minor_align=split_minor_align,split_minor_unalign=split_minor_unalign,cls=cls)
def pct(a,b): return round(100*a/b,2) if b else 0.0
def summarize(records,name):
    n=len(records)
    precond=[r for r in records if r.get("n_complete") and r["n_complete"]>=6]
    npre=len(precond)
    byT={}
    for T in (3,4,5): byT[T]=classify(records,T)
    h=byT[4]  # headline
    # minority-size distribution (T4 split loci)
    mn=np.array(h["split_minor"],float) if h["split_minor"] else np.array([])
    mna=np.array(h["split_minor_align"],float) if h["split_minor_align"] else np.array([])
    mnu=np.array(h["split_minor_unalign"],float) if h["split_minor_unalign"] else np.array([])
    def qd(a):
        return {} if a.size==0 else dict(n=int(a.size),min=int(a.min()),p25=float(np.percentile(a,25)),
            median=float(np.median(a)),p75=float(np.percentile(a,75)),max=int(a.max()),mean=round(float(a.mean()),2))
    # minority size histogram (1..12+)
    def hist(a):
        h={}
        for v in a: k=str(int(v)) if v<12 else "12+"; h[k]=h.get(k,0)+1
        return h
    out=dict(name=name,n_total=n,n_precond=npre,pct_precond=pct(npre,n),
        headline_T=4,
        states_T4=dict(S1=h["S1"],S2=h["S2"],S3=h["S3"],S3_loc=h["S3_loc"],S3_disp=h["S3_disp"],S4=h["S4"],S5=h["S5"]),
        pct_T4={k:pct(h[k] if k in("S1","S2","S4","S5") else (h["S3"] if k=="S3" else 0),n) for k in("S1","S2","S3","S4","S5")},
        # effectiveness metrics
        split_align_rate=pct(h["S5"],h["n_split"]),          # ⑤/(④+⑤) 切群是否找到生物結構
        nonsplit_meanshift_rate=pct(h["S3"],h["n_nonsplit_precond"]),  # ③/(②+③) 切不出≠沒訊號
        nonsplit_loc_rate=pct(h["S3_loc"],h["n_nonsplit_precond"]),
        s3_loc_of_s3=pct(h["S3_loc"],h["S3"]) if h["S3"] else 0.0,
        minority_dist=dict(all=qd(mn),aligned=qd(mna),unaligned=qd(mnu),
            hist_all=hist(mn),hist_aligned=hist(mna),hist_unaligned=hist(mnu)),
        threshold_sensitivity={T:dict(S1=byT[T]["S1"],S2=byT[T]["S2"],S3=byT[T]["S3"],
            S4=byT[T]["S4"],S5=byT[T]["S5"],n_split=byT[T]["n_split"],
            split_align_rate=pct(byT[T]["S5"],byT[T]["n_split"])) for T in (3,4,5)},
        cls_T4=byT[4]["cls"])
    return out
def main():
    res={}
    for nm in ("tumor","merged"):
        p=f"{A}/decisionflow_records_{nm}.json"
        if not os.path.exists(p): print(f"MISSING {p}",file=sys.stderr); continue
        recs=json.load(open(p))
        res[nm]=summarize(recs,nm)
        s=res[nm]; h=s["states_T4"]
        print(f"\n===== {nm} (N={s['n_total']}, precond-pass={s['n_precond']} {s['pct_precond']}%) =====")
        print(f"  T4 states: S1={h['S1']} S2={h['S2']} S3={h['S3']}(loc{h['S3_loc']}/disp{h['S3_disp']}) S4={h['S4']} S5={h['S5']}")
        print(f"  pct(of all): {s['pct_T4']}")
        print(f"  切群對齊率 ⑤/(④+⑤)={s['split_align_rate']}%  切不出有訊號 ③/(②+③)={s['nonsplit_meanshift_rate']}% (location-clean {s['nonsplit_loc_rate']}%)")
        print(f"  minority(all) median={s['minority_dist']['all'].get('median')} aligned median={s['minority_dist']['aligned'].get('median')} unaligned median={s['minority_dist']['unaligned'].get('median')}")
        print(f"  threshold sens: "+" ".join(f"T{T}:split{v['n_split']}(align{v['split_align_rate']}%)" for T,v in s['threshold_sensitivity'].items()))
    # strip cls arrays from json (large) into separate file
    full={nm:{k:v for k,v in s.items() if k!="cls_T4"} for nm,s in res.items()}
    json.dump(full,open(f"{A}/decisionflow_summary.json","w"),indent=1)
    json.dump({nm:s["cls_T4"] for nm,s in res.items()},open(f"{A}/decisionflow_cls_T4.json","w"))
    print(f"\nWROTE decisionflow_summary.json + decisionflow_cls_T4.json")
if __name__=="__main__": main()
