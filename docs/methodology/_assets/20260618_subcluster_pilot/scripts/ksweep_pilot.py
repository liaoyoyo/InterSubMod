#!/usr/bin/env python3
"""k-sweep pilot: 比較『silhouette 最佳 k(現行)』vs『對齊 a-priori 軸 CramerV 最佳 k(supervised)』。
對 chr22 每個 region 掃 k=2..min(5,n//3)，每 k 算 silhouette + 各 a-priori 軸(HP/CARRIER/ALLELE)的 CramerV+chi2 p。
回答: 對齊最佳 k 是否 != silhouette 最佳 k? 其他 k 是否更顯著(對齊更好)?
PILOT — 單樣本 HCC1395 chr22 subset。chi2 p(非 Fisher-FH)+ 未做 k-sweep 多重比較校正(故意, 看上限)。"""
import os, csv, glob, json
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.metrics import silhouette_score
from collections import Counter, defaultdict

OUT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_tmp_ksweep_chr22"
ASSET="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
MIN_SZ=3; LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}

def _f(x):
    x=x.strip(); return np.nan if x in("","NA","nan","NaN") else float(x)
def load_matrix(path):
    rows=open(path).read().strip().split("\n"); ids=[]; M=[]
    for ln in rows[1:]:
        p=ln.split(","); ids.append(p[0]); M.append([_f(x) for x in p[1:]])
    return ids,np.array(M)
def peel(sub):
    idx=list(range(sub.shape[0]))
    while True:
        S=sub[np.ix_(idx,idx)]; bad=(S<0)|np.isnan(S); np.fill_diagonal(bad,False)
        if not bad.any(): return idx
        idx.remove(idx[int(np.argmax(bad.sum(1)))])
        if len(idx)<MIN_SZ*2: return idx
def cramerv(groups, clusters):
    # groups: list of 2-level axis label; clusters: cluster id per read (same order)
    g=sorted(set(groups)); c=sorted(set(clusters))
    if len(g)<2 or len(c)<2: return None,None
    tab=np.zeros((len(g),len(c)))
    gi={x:i for i,x in enumerate(g)}; ci={x:i for i,x in enumerate(c)}
    for a,b in zip(groups,clusters): tab[gi[a],ci[b]]+=1
    if (tab.sum(0)==0).any() or (tab.sum(1)==0).any():
        tab=tab[tab.sum(1)>0][:,tab.sum(0)>0]
    if tab.shape[0]<2 or tab.shape[1]<2: return None,None
    try: chi2,p,_,_=chi2_contingency(tab)
    except Exception: return None,None
    n=tab.sum(); v=np.sqrt(chi2/(n*(min(tab.shape)-1))) if n>0 else None
    return (float(v) if v is not None else None), float(p)

mats=glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True)
recs=[]
for mp in mats:
    rdir=os.path.dirname(os.path.dirname(os.path.dirname(mp)))
    rt=os.path.join(rdir,"reads","reads.tsv")
    if not os.path.exists(rt): continue
    pos=None
    for part in rdir.split("/"):
        if part.startswith("chr22_") and part.count("_")==1: pos=part.split("_")[1]; break
    reads={r["read_id"]:r for r in csv.DictReader(open(rt),delimiter="\t")}
    ids,M=load_matrix(mp); id2row={rid:i for i,rid in enumerate(ids)}
    def istum(t): return str(t) in ("1","true","True")
    rows=[]; lab=[]; allele=[]
    for rid in ids:
        r=reads.get(rid)
        if r and istum(r["is_tumor"]) and r["hp"] in LABMAP:
            rows.append(id2row[rid]); lab.append(LABMAP[r["hp"]]); allele.append(r["alt_support"])
    if len(rows)<MIN_SZ*2: continue
    sub=M[np.ix_(rows,rows)]; keep=peel(sub)
    if len(keep)<MIN_SZ*2: continue
    lab=[lab[i] for i in keep]; allele=[allele[i] for i in keep]
    sub=sub[np.ix_(keep,keep)].copy(); np.fill_diagonal(sub,0.0); sub=np.maximum(sub,sub.T)
    n=len(keep)
    try: Z=linkage(squareform(sub,checks=False),method="average")
    except Exception: continue
    # axes: 2-level groupings
    hp=["HP1" if l in ("1","1-1") else "HP2" for l in lab]
    carrier=["G" if l in ("1","2") else "C" for l in lab]   # germline vs carrier
    per_k=[]
    for k in range(2,min(5,n//MIN_SZ)+1):
        cl=fcluster(Z,k,criterion="maxclust"); sz=Counter(cl)
        if len(sz)<k or min(sz.values())<MIN_SZ: continue
        try: sil=float(silhouette_score(sub,cl,metric="precomputed"))
        except Exception: continue
        vhp,php=cramerv(hp,list(cl)); vc,pc=cramerv(carrier,list(cl)); va,pa=cramerv(allele,list(cl))
        per_k.append({"k":k,"sil":round(sil,3),
                      "V_hp":vhp,"p_hp":php,"V_carrier":vc,"p_carrier":pc,"V_allele":va,"p_allele":pa})
    if not per_k: continue
    sil_opt=max(per_k,key=lambda d:d["sil"])
    def align_opt(axis):
        cand=[d for d in per_k if d[f"V_{axis}"] is not None]
        return max(cand,key=lambda d:d[f"V_{axis}"]) if cand else None
    recs.append({"pos":pos,"n":n,"ks":[d["k"] for d in per_k],
        "sil_opt_k":sil_opt["k"],"sil_opt_sil":sil_opt["sil"],
        "sil_opt_V_hp":sil_opt["V_hp"],"sil_opt_V_carrier":sil_opt["V_carrier"],"sil_opt_V_allele":sil_opt["V_allele"],
        "alignopt_hp":align_opt("hp"),"alignopt_carrier":align_opt("carrier"),"alignopt_allele":align_opt("allele"),
        "per_k":per_k})

# ---- aggregate ----
def summarize(axis):
    out={"n_eval":0,"align_eq_sil":0,"align_diff_sil":0,"higher_k_more_sig":0,"dV_gt_0.1_at_other_k":0,
         "examples":[]}
    for r in recs:
        ao=r[f"alignopt_{axis}"]; sv=r[f"sil_opt_V_{axis}"]
        if ao is None or sv is None: continue
        out["n_eval"]+=1
        if ao["k"]==r["sil_opt_k"]: out["align_eq_sil"]+=1
        else:
            out["align_diff_sil"]+=1
            if ao["k"]>r["sil_opt_k"]: out["higher_k_more_sig"]+=1
        if ao[f"V_{axis}"]-sv>0.1:
            out["dV_gt_0.1_at_other_k"]+=1
            if len(out["examples"])<6:
                out["examples"].append({"pos":r["pos"],"n":r["n"],"sil_opt_k":r["sil_opt_k"],
                    f"sil_V_{axis}":round(sv,3),"align_k":ao["k"],f"align_V_{axis}":round(ao[f"V_{axis}"],3),
                    f"align_p_{axis}":ao[f"p_{axis}"]})
    return out

agg={"sample":"HCC1395 chr22 (PILOT subset 1/22 chr)","n_regions_clusterable":len(recs),
     "HP_axis":summarize("hp"),"CARRIER_axis":summarize("carrier"),"ALLELE_axis":summarize("allele")}
json.dump({"agg":agg,"records":recs},open(f"{ASSET}/ksweep_pilot_chr22.json","w"),indent=1)
print(json.dumps(agg,ensure_ascii=False,indent=1))
