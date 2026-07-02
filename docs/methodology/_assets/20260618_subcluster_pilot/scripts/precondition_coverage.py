#!/usr/bin/env python3
"""驗證前提覆蓋率: 各軸需 >=2 種不同標籤(各>=3 reads)才能測。
從 records_wg2.json 的 tumor 標籤計數(1/1-1/2/2-1)直接算 + significance allele valid。
回答「會不會很多全是1個標籤」「多少位點不同軸可驗證」。"""
import json, csv
import numpy as np
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
recs=json.load(open(f"{A}/records_wg2.json")); N=len(recs)
MIN=3
from collections import Counter
ndist=Counter(); hp_both=hp_both3=carr_both=carr_both3=0
hpfine_ge2=hpfine3=hpfine4=0
keyset={f"{r['chrom']}_{r['pos']}":r for r in recs}
for r in recs:
    L=r['all']['labels']
    c1=L.get('1',0); c11=L.get('1-1',0); c2=L.get('2',0); c21=L.get('2-1',0)
    present=[x for x in (c1,c11,c2,c21) if x>0]
    ndist[len(present)]+=1
    HP1=c1+c11; HP2=c2+c21; germ=c1+c2; carr=c11+c21
    if HP1>0 and HP2>0: hp_both+=1
    if HP1>=MIN and HP2>=MIN: hp_both3+=1
    if germ>0 and carr>0: carr_both+=1
    if germ>=MIN and carr>=MIN: carr_both3+=1
    ge3=sum(1 for x in (c1,c11,c2,c21) if x>=MIN)
    if ge3>=2: hpfine_ge2+=1
    if ge3>=3: hpfine3+=1
    if ge3==4: hpfine4+=1
# allele 前提: significance LabelAllelePermanovaValid
SIG=f"{WT}/output/_wg_bdcprime_verify/significance_summary.csv"
def tru(v): return str(v).strip().lower() in("true","1","yes")
al_valid=hp_valid_csv=0
for row in csv.DictReader(open(SIG)):
    if f"{row['Chr']}_{row['Pos']}" not in keyset: continue
    if tru(row.get('LabelAllelePermanovaValid')): al_valid+=1
    if tru(row.get('LabelHPPermanovaValid')): hp_valid_csv+=1
def pc(x): return round(100*x/N,1)
out={"N":N,"MIN_reads":MIN,
 "n_distinct_labels":{str(k):ndist.get(k,0) for k in (1,2,3,4)},
 "n_distinct_pct":{str(k):pc(ndist.get(k,0)) for k in (1,2,3,4)},
 "single_label":ndist.get(1,0),"single_pct":pc(ndist.get(1,0)),
 "multi_label":N-ndist.get(1,0),"multi_pct":pc(N-ndist.get(1,0)),
 "HP_axis":{"both_present":hp_both,"both_present_pct":pc(hp_both),"both_ge3":hp_both3,"both_ge3_pct":pc(hp_both3),"csv_valid":hp_valid_csv,"csv_valid_pct":pc(hp_valid_csv)},
 "carrier_axis":{"both_present":carr_both,"both_present_pct":pc(carr_both),"both_ge3":carr_both3,"both_ge3_pct":pc(carr_both3)},
 "allele_axis":{"csv_valid":al_valid,"csv_valid_pct":pc(al_valid)},
 "hpfine":{"ge2_of4":hpfine_ge2,"ge2_pct":pc(hpfine_ge2),"ge3_of4":hpfine3,"ge3_pct":pc(hpfine3),"all4":hpfine4,"all4_pct":pc(hpfine4)}}
json.dump(out,open(f"{A}/precondition_coverage.json","w"),indent=1)
print(f"N={N} (tumor 標籤)")
print(f"\n單一標籤(只1種HP tag): {out['single_label']} ({out['single_pct']}%)  ← 全1標籤、無法做標籤檢定")
print(f"多標籤(>=2種): {out['multi_label']} ({out['multi_pct']}%)")
print(f"標籤種類分佈(1/2/3/4種): {out['n_distinct_labels']}  %={out['n_distinct_pct']}")
print(f"\n各軸前提(>=2 組各>=3 reads 才能測):")
print(f"  HP 軸(HP1-fam vs HP2-fam): both>=3 = {hp_both3} ({pc(hp_both3)}%)  [both present {hp_both} ({pc(hp_both)}%)]")
print(f"  CARRIER 軸(germline vs carrier): both>=3 = {carr_both3} ({pc(carr_both3)}%)  ← subclone 前提,最嚴")
print(f"  ALLELE 軸(REF vs ALT, csv valid): {al_valid} ({pc(al_valid)}%)")
print(f"\n  HP-fine 4組: >=2組各>=3: {hpfine_ge2} ({pc(hpfine_ge2)}%); >=3組: {hpfine3} ({pc(hpfine3)}%); 全4組: {hpfine4} ({pc(hpfine4)}%)")
