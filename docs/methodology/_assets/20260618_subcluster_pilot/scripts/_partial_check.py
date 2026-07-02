#!/usr/bin/env python3
"""進度快照分類確認: 對已完成 chr 的 records 跑 5-態分類 + 各態實例。標 PARTIAL。"""
import sys, json, os
from collections import Counter
A=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,os.path.join(A,"scripts"))
import decisionflow_analyze as DA
STLAB={"S1":"①不可驗證","S2":"②1群/無訊號","S3_loc":"③mean-shift(location淨)","S3_disp":"③mean-shift(dispersion)","S4":"④可分未對齊","S5":"⑤確認真結構"}
for nm in ("tumor","merged"):
    snap=f"{A}/_snap_{nm}.json"
    if not os.path.exists(snap): print(f"[{nm}] no snapshot"); continue
    recs=json.load(open(snap))
    chroms=Counter(r["chrom"] for r in recs)
    s=DA.summarize(recs,nm)
    cls=DA.classify(recs,4)["cls"]
    h=s["states_T4"]
    print(f"\n========== {nm}  (PARTIAL: {len(recs)} 位點, chr={sorted(chroms, key=lambda x:int(x[3:]))}) ==========")
    print(f"  precond-pass={s['n_precond']} ({s['pct_precond']}%)")
    print(f"  ① 不可驗證        : {h['S1']:>6}  ({100*h['S1']/len(recs):.1f}%)")
    print(f"  ② 1群/無訊號      : {h['S2']:>6}  ({100*h['S2']/len(recs):.1f}%)")
    print(f"  ③ 監督可分mean-shift: {h['S3']:>6}  ({100*h['S3']/len(recs):.1f}%)  [location淨 {h['S3_loc']} / dispersion {h['S3_disp']}]")
    print(f"  ④ 可分未對齊      : {h['S4']:>6}  ({100*h['S4']/len(recs):.1f}%)")
    print(f"  ⑤ 確認真結構      : {h['S5']:>6}  ({100*h['S5']/len(recs):.1f}%)")
    print(f"  → 切群對齊率 ⑤/(④+⑤) = {s['split_align_rate']}%   切不出有訊號 ③/(②+③) = {s['nonsplit_meanshift_rate']}% (loc {s['nonsplit_loc_rate']}%)")
    # 各態 2 實例
    by={}
    for i,c in enumerate(cls): by.setdefault(c,[]).append(i)
    print("  -- 各態實例 --")
    for st in ("S1","S2","S3_loc","S3_disp","S4","S5"):
        idxs=by.get(st,[])[:2]
        for i in idxs:
            r=recs[i]; bt=r.get("byT",{}).get("4",{}); pm=r.get("perm",{})
            ps=";".join(f"{ax}:F{pm[ax]['F']}/p{pm[ax]['p']}/dP{pm[ax]['dispP']}" for ax in("hp","carrier","allele") if pm.get(ax))
            sp=f"split n_valid={bt.get('n_valid')} minority={bt.get('minority')} out={bt.get('n_out')} aligned={bt.get('aligned')}({bt.get('ax')})" if bt.get("split") else "切不出"
            print(f"    {STLAB[st]:<22} {r['chrom']}:{r['pos']} n={r['n']} comp={r['n_complete']} | {sp}")
            if ps: print(f"        perm: {ps}")
