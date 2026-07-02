#!/usr/bin/env python3
"""多角度驗證「切不出 ≠ 沒訊號」: 效應量(Δβ) + 置換null對照(PERMANOVA是999置換檢定) + per-CpG + 軸分解。
切不出(best_k None,n>=6) vs 切得出(best_k>=2) vs insuff(n<6)。零 compute(讀 records + significance_csv)。
存 cantsplit_validation.json (含直方圖 bins 供畫圖)。"""
import json, csv
import numpy as np
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
SIG=f"{WT}/output/_wg_bdcprime_verify/significance_summary.csv"
recs=json.load(open(f"{A}/records_wg2.json"))
def lab_n(r): return len([l for l,c in r["all"]["labels"].items() if c>0])
grp={}
for r in recs:
    bk=r["all"]["best_k"]; n=r["all"]["n"]
    if n<6: g="insuff"
    elif bk is None: g="切不出"+("·multi" if lab_n(r)>=2 else "·single")
    elif bk>=2: g="切得出"
    else: g="other"
    grp[f"{r['chrom']}_{r['pos']}"]=g
def ff(v):
    try: return float(v)
    except: return None
def truthy(v): return str(v).strip().lower() in("true","1","yes")
data={k:{"hpP":[],"hpValid":0,"germ":[],"sub":[],"frac":[],"epi":[],"alleP":[]} for k in set(grp.values())}
for row in csv.DictReader(open(SIG)):
    key=f"{row['Chr']}_{row['Pos']}"; g=grp.get(key)
    if not g: continue
    d=data[g]
    if truthy(row.get("LabelHPPermanovaValid")):
        d["hpValid"]+=1; p=ff(row.get("LabelHPPermanovaP"))
        if p is not None: d["hpP"].append(p)
    gm=ff(row.get("GermlineAsmDbeta"))
    if gm is not None: d["germ"].append(abs(gm))
    s1=ff(row.get("SubcloneDbeta_HP1")); s2=ff(row.get("SubcloneDbeta_HP2"))
    sv=[abs(x) for x in (s1,s2) if x is not None]
    if sv: d["sub"].append(max(sv))
    fr=ff(row.get("Fisher_Frac_Sig"))
    if fr is not None: d["frac"].append(fr)
    ep=ff(row.get("Epipoly_HP1"))
    if ep is not None: d["epi"].append(ep)
    ap=ff(row.get("LabelAllelePermanovaP"))
    if ap is not None: d["alleP"].append(ap)

def summ(d,n_total):
    hp=np.array(d["hpP"]);
    def med(x): return round(float(np.median(x)),3) if len(x) else None
    return dict(
        matched=n_total,
        hp_valid=d["hpValid"],
        hp_p05=int((hp<0.05).sum()) if len(hp) else 0,
        hp_p05_pct=round(100*(hp<0.05).mean(),1) if len(hp) else None,
        hp_p001_pct=round(100*(hp<=0.0011).mean(),1) if len(hp) else None,  # 999-perm floor
        germ_med=med(d["germ"]), germ_n=len(d["germ"]),
        germ_gt01=int((np.array(d["germ"])>0.1).sum()) if d["germ"] else 0,
        sub_med=med(d["sub"]),
        frac_med=med(d["frac"]), frac_gt0=int((np.array(d["frac"])>0).sum()) if d["frac"] else 0,
        epi_med=med(d["epi"]))
import collections
cnt=collections.Counter(grp.values())
out={"groups":{},"hist":{}}
for g,d in data.items():
    out["groups"][g]=summ(d,cnt[g])
# 直方圖: germ |Δβ| 切不出·multi vs 切得出 (0..0.6, 12 bins)
bins=np.linspace(0,0.6,13)
for g in ("切不出·multi","切得出"):
    h,_=np.histogram(np.clip(data[g]["germ"],0,0.6),bins=bins)
    out["hist"][f"germ_{g}"]=h.tolist()
out["hist"]["bins"]=[round(x,3) for x in bins.tolist()]
# 置換 null 期望: random 下 p<0.05 應 ~5%, p<=.001 應 ~0.1%
out["null_expect"]={"p05":5.0,"p001":0.1}
json.dump(out,open(f"{A}/cantsplit_validation.json","w"),indent=1)
# print
print("group                    matched  HP-valid  HP p<.05    HP p<=.001  |germΔβ|med  germ>0.1  Fisher_frac>0  subΔβ-med")
for g in ("切不出·multi","切不出·single","切得出","insuff"):
    s=out["groups"].get(g)
    if not s: continue
    print(f"{g:<22} {s['matched']:>7}  {s['hp_valid']:>7}  {str(s['hp_p05'])+' ('+str(s['hp_p05_pct'])+'%)':<12} {str(s['hp_p001_pct'])+'%':<10} {str(s['germ_med']):<11} {s['germ_gt01']:>7}  {str(s['frac_med'])+' ('+str(s['frac_gt0'])+')':<14} {s['sub_med']}")
print(f"\n置換 null 期望(隨機): p<.05 應~5%, p<=.001 應~0.1%  → 切不出 遠超 = 訊號真(非 chance)")
