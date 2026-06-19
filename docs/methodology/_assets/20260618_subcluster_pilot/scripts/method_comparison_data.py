#!/usr/bin/env python3
"""統一 Fisher+V vs PERMANOVA 方法對照數字 → method_comparison.json (供 HTML §13-A 注入)。
Q1 軸/多組 / Q2 dispersion / Q3 raw+clean Venn / Q4 Δβ-by-class。讀 records + significance_csv + permanova_clean_4group.json。"""
import json, csv
import numpy as np
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
SIG=f"{WT}/output/_wg_bdcprime_verify/significance_summary.csv"
recs={f"{r['chrom']}_{r['pos']}":r for r in json.load(open(f"{A}/records_wg2.json"))}
N=len(recs)
def tru(v): return str(v).strip().lower() in("true","1","yes")
def ff(v):
    try: return float(v)
    except: return None
cansplit={k for k,r in recs.items() if r['all']['best_k'] is not None and r['all']['best_k']>=2}
hp_sig=hp_disp=hp_clean=al_sig=al_disp=0
permraw=set(); clean=set(); dclean=[]; ddisp=[]; dnon=[]
for row in csv.DictReader(open(SIG)):
    key=f"{row['Chr']}_{row['Pos']}"
    if key not in recs: continue
    hv=tru(row.get('LabelHPPermanovaValid')); hp=ff(row.get('LabelHPPermanovaP')); hw=tru(row.get('LabelHPDispersionWarn'))
    av=tru(row.get('LabelAllelePermanovaValid')); ap=ff(row.get('LabelAllelePermanovaP')); aw=tru(row.get('LabelAlleleDispersionWarn'))
    gm=ff(row.get('GermlineAsmDbeta'))
    hsig=hv and hp is not None and hp<0.05; asig=av and ap is not None and ap<0.05
    if hsig:
        hp_sig+=1
        if hw: hp_disp+=1
        else: hp_clean+=1
    if asig:
        al_sig+=1
        if aw: al_disp+=1
    if hsig or asig: permraw.add(key)
    if (hsig and not hw) or (asig and not aw): clean.add(key)
    if gm is not None and hv:
        g=abs(gm)
        if hsig and not hw: dclean.append(g)
        elif hsig and hw: ddisp.append(g)
        elif not hsig: dnon.append(g)
def pc(x): return round(100*x/N,1)
A_,B,C=cansplit,permraw,clean
q2=json.load(open(f"{A}/permanova_clean_4group.json"))["Q2"]
out={
 "N":N,
 "fisher_v":{"cansplit":len(A_),"pct":pc(len(A_))},
 "permanova_raw":{"hp":hp_sig,"hp_pct":pc(hp_sig),"al":al_sig,"al_pct":pc(al_sig),"any":len(B),"any_pct":pc(len(B))},
 "dispersion":{"hp_sig":hp_sig,"hp_disp":hp_disp,"hp_disp_pct":round(100*hp_disp/hp_sig,1),
   "hp_clean":hp_clean,"hp_clean_pct":round(100*hp_clean/hp_sig,1),"hp_clean_of_all":pc(hp_clean),
   "al_disp":al_disp,"al_disp_pct":round(100*al_disp/al_sig,1)},
 "venn_raw":{"AnB":len(A_&B),"A_only":len(A_-B),"B_only":len(B-A_),"neither":N-len(A_|B),
   "AnB_pct":pc(len(A_&B)),"B_only_pct":pc(len(B-A_)),"neither_pct":pc(N-len(A_|B))},
 "venn_clean":{"C":len(C),"C_pct":pc(len(C)),"AnC":len(A_&C),"AnC_pct":pc(len(A_&C)),
   "A_only":len(A_-C),"A_only_pct":pc(len(A_-C)),"C_only":len(C-A_),"C_only_pct":pc(len(C-A_)),
   "jaccard":round(len(A_&C)/len(A_|C),3)},
 "delta_by_class":{"clean_med":round(float(np.median(dclean)),3),"clean_n":len(dclean),
   "disp_med":round(float(np.median(ddisp)),3),"disp_n":len(ddisp),
   "nonsig_med":round(float(np.median(dnon)),3),"nonsig_n":len(dnon)},
 "fourgroup":q2,
}
json.dump(out,open(f"{A}/method_comparison.json","w"),indent=1)
print(json.dumps({k:out[k] for k in("fisher_v","permanova_raw","dispersion","venn_raw","venn_clean","delta_by_class")},ensure_ascii=False,indent=1))
print("4-group:",out["fourgroup"])
