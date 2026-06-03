#!/usr/bin/env python3
"""彙總 admr_chr*.json：算 aDMR 落 CNV/LOH 的真富集（OR + Fisher），對照文獻 79% 與背景率。"""
import json,glob,os,numpy as np
from scipy.stats import fisher_exact
AD=os.path.dirname(os.path.abspath(__file__))
wins=[]
for fp in glob.glob(AD+"/admr_chr*.json"):
    d=json.load(open(fp))
    if d.get("n_windows",0)>0: wins.extend(d["windows"])
admr=[w for w in wins if w["is_admr"]]
non=[w for w in wins if not w["is_admr"]]
# 2x2: aDMR vs non, in_cnvloh vs not
a=sum(1 for w in admr if w["in_cnvloh"]); b=len(admr)-a
c=sum(1 for w in non if w["in_cnvloh"]); d2=len(non)-c
OR=None;pf=None
if b>0 and c>0 and d2>0:
    try:
        OR,pf=fisher_exact([[a,b],[c,d2]])
    except Exception:pass
elif b==0 or d2==0:
    OR=float('inf') if (b==0 and a>0) else OR
def frac(g,k): return round(sum(1 for w in g if w[k])/len(g),3) if g else None
out={
  "n_windows":len(wins),"n_admr":len(admr),"frac_admr":round(len(admr)/len(wins),3) if wins else None,
  "admr_frac_in_cnvloh":frac(admr,"in_cnvloh"),
  "admr_frac_in_loh":frac(admr,"in_loh"),
  "nonadmr_frac_in_cnvloh":frac(non,"in_cnvloh"),
  "background_frac_in_cnvloh":frac(wins,"in_cnvloh"),
  "enrichment_OR_admr_vs_non":round(float(OR),3) if OR not in (None,float('inf')) else str(OR),
  "fisher_p":round(float(pf),5) if pf is not None else None,
  "literature_ref":"79% aDMR in CNV/LOH (advanced cancer ONT cohort)",
  "interpretation_note":"HCC1395 背景 CNV/LOH 覆蓋極高(~90%); 須看 OR 是否>1(真富集) 而非絕對%",
  # delta 分層：aDMR 強度 vs non
  "admr_maxdelta_median":round(float(np.median([w["max_delta"] for w in admr])),3) if admr else None,
  "nonadmr_maxdelta_median":round(float(np.median([w["max_delta"] for w in non])),3) if non else None,
}
json.dump(out,open(AD+"/admr_aggregate.json","w"),ensure_ascii=False,indent=2)
print(json.dumps(out,ensure_ascii=False,indent=2))
