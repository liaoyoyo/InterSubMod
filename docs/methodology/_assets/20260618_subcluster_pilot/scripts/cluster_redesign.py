#!/usr/bin/env python3
"""新切群模組（純 Python，operate on 既有 BERNOULLI 距離矩陣）— 對應三個觀察：
  FM1 隔離離群：<min_sz 的小群當「離群 read」單獨記錄，不占群名額、不害位點被丟。
  FM2 自然 gap 切：依樹高間隙(distance 門檻)切，取代「強制 k 群」maxclust。
  FM3 列舉解析度：不折成一個 k，列出所有乾淨分離解析度 + 每 read 群歸屬 + 各軸對齊。
匯入用（函式）或直跑（對 29 example 出 cluster_redesign.json）。切群法基底 = UPGMA(average)。
"""
import os, csv, glob, json
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.metrics import silhouette_score
from collections import Counter

MIN_SZ=3; MAX_OUT_FRAC=0.20
KMAX=6; SIL_MIN=0.15; GAP_MULT=2.0   # 分離度閘：只在顯著 gap 切 + silhouette 門檻 + 群數上限（防過度碎裂）
# 群下限 = MIN_SZ=3（業界 caller 最低支持 read 標準；「多數足以確認非雜訊」）。
# 「明顯」不靠群大小門檻 → 交給 null-excess(扣 within-1-group null≥0.10) + 分離度判定。
def grp_min(n): return MIN_SZ
LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}

def _f(x):
    x=str(x).strip(); return np.nan if x in("","NA","nan","NaN") else float(x)
def loadm(p):
    r=open(p).read().strip().split("\n"); ids=[];M=[]
    for ln in r[1:]:
        q=ln.split(","); ids.append(q[0]); M.append([_f(x) for x in q[1:]])
    return ids,np.array(M)
def peel(s):
    idx=list(range(s.shape[0]))
    while True:
        S=s[np.ix_(idx,idx)];b=(S<0)|np.isnan(S);np.fill_diagonal(b,False)
        if not b.any(): return idx
        idx.remove(idx[int(np.argmax(b.sum(1)))])
        if len(idx)<MIN_SZ*2: return idx
def linkZ(sub):
    s=sub.copy(); np.fill_diagonal(s,0.0); s=np.maximum(s,s.T)
    return linkage(squareform(s,checks=False),method="average"), s

# ---------- within-1-group null（定「真實」: 扣掉 read-內相關造的假穩）----------
def bernoulli_dist(P):
    """向量化重算 BERNOULLI 距離(與 C++ calculate_bernoulli 一致)。NaN=missing。"""
    n,C=P.shape; Mk=(~np.isnan(P)).astype(float); Pz=np.where(np.isnan(P),0.0,P)
    Aw=2.0*np.abs(Pz-0.5)*Mk; Bw=Aw*Pz
    sw=Aw@Aw.T; swd=Bw@Aw.T+Aw@Bw.T-2.0*(Bw@Bw.T); common=Mk@Mk.T
    with np.errstate(divide="ignore",invalid="ignore"): Dm=swd/sw
    Dm[(common<MIN_SZ)|(sw<1e-9)]=-1.0; np.fill_diagonal(Dm,0.0); return Dm
def _clusterboot_meanj(sub, cl0, k, rng, B=200, frac=0.8):
    n=sub.shape[0]; m=max(MIN_SZ*2,int(round(frac*n)))
    orig=[set(np.where(cl0==c)[0]) for c in sorted(set(cl0))]; acc=[[] for _ in orig]
    for _ in range(B):
        samp=np.sort(rng.choice(n,size=m,replace=False)); ss=sub[np.ix_(samp,samp)]
        if ((ss<0)|np.isnan(ss)).any(): continue
        try: cl=fcluster(linkage(squareform(np.maximum(ss,ss.T),checks=False),"average"),k,"maxclust")
        except Exception: continue
        bc=[set(samp[cl==c].tolist()) for c in sorted(set(cl))]; sset=set(samp.tolist())
        for oi,oc in enumerate(orig):
            ocs=oc&sset
            if len(ocs)<2: continue
            acc[oi].append(max((len(ocs&b)/len(ocs|b) for b in bc),default=0.0))
    pc=[float(np.mean(a)) if a else 0.0 for a in acc]; return float(np.mean(pc)) if pc else 0.0
def stab_excess(sub, P, k, rng, Rnull=25):
    """real clusterboot Jaccard − within-1-group null(95pct)。"""
    try: cl0=fcluster(linkZ(sub)[0],k,"maxclust")
    except Exception: return None
    if len(set(cl0))<k or min(Counter(cl0).values())<MIN_SZ: return None
    real=_clusterboot_meanj(sub,cl0,k,rng)
    nulls=[]
    n,C=P.shape
    for _ in range(Rnull):
        Pn=P.copy()
        for c in range(C):
            col=Pn[:,c]; vi=np.where(~np.isnan(col))[0]
            if vi.size>1: Pn[vi,c]=col[rng.permutation(vi)]
        Dn=bernoulli_dist(Pn); sub2,keep=clean_sub(Dn)
        if sub2 is None or sub2.shape[0]<k*MIN_SZ: nulls.append(0.0); continue
        try: cln=fcluster(linkZ(sub2)[0],k,"maxclust")
        except Exception: nulls.append(0.0); continue
        if len(set(cln))<k or min(Counter(cln).values())<MIN_SZ: nulls.append(0.0); continue
        nulls.append(_clusterboot_meanj(sub2,cln,k,rng,B=60))
    return round(real-float(np.percentile(nulls,95)),3)
def clean_sub(D):
    keep=peel(D)
    if len(keep)<MIN_SZ*2: return None,None
    sub=D[np.ix_(keep,keep)].copy()
    if ((sub<0)|np.isnan(sub)).any(): return None,None
    np.fill_diagonal(sub,0.0); sub=np.maximum(sub,sub.T); return sub,keep
def cramerv(groups,clusters):
    g=sorted(set(groups)); c=sorted(set(clusters))
    if len(g)<2 or len(c)<2: return None,None,None
    tab=np.zeros((len(g),len(c))); gi={x:i for i,x in enumerate(g)}; ci={x:i for i,x in enumerate(c)}
    for a,b in zip(groups,clusters): tab[gi[a],ci[b]]+=1
    tab=tab[tab.sum(1)>0][:,tab.sum(0)>0]
    if tab.shape[0]<2 or tab.shape[1]<2: return None,None,None
    try: chi2,p,dof,exp=chi2_contingency(tab)
    except Exception: return None,None,None
    nn=tab.sum(); v=float(np.sqrt(chi2/(nn*(min(tab.shape)-1)))) if nn>0 else None
    return v,float(p),float(exp.min())

# ---------- OLD: maxclust + min3-reject（現行）----------
def old_cut(sub, min_sz=MIN_SZ):
    n=sub.shape[0]; Z,s=linkZ(sub); accepted=[]
    for k in range(2,min(6,n//min_sz)+1):
        if k>=n: break
        cl=fcluster(Z,k,"maxclust"); sz=Counter(cl)
        if len(sz)>=k and min(sz.values())>=min_sz: accepted.append((k,cl.tolist(),sorted(sz.values(),reverse=True)))
    if not accepted:
        return {"status":"CANT_SPLIT","k_groups":1,"accepted_ks":[],"labels":None,"sizes":[n]}
    k,cl,sz=accepted[-1]  # richest accepted maxclust k
    return {"status":"SPLIT","k_groups":k,"accepted_ks":[a[0] for a in accepted],"labels":cl,"sizes":sz}

# ---------- NEW: quarantine + gap-cut + enumerate ----------
def new_cut(sub, min_sz=MIN_SZ, max_out_frac=MAX_OUT_FRAC, kmax=KMAX, sil_min=SIL_MIN, gap_mult=GAP_MULT):
    """隔離離群 + 只在『顯著 gap』切 + silhouette 分離度閘 + k_core 上限（防過度碎裂）。"""
    n=sub.shape[0]; Z,s=linkZ(sub); gm=grp_min(n)
    heights=np.sort(Z[:,2]); gaps=np.diff(heights)
    pos=gaps[gaps>0]; med=float(np.median(pos)) if pos.size else 0.0
    res={}  # k_core -> best (by core silhouette)
    for i in range(len(gaps)):
        if med>0 and gaps[i] < gap_mult*med: continue          # 只在顯著 gap 切
        t=(heights[i]+heights[i+1])/2.0
        cl=fcluster(Z,t,"distance"); cnt=Counter(cl)
        core_lbls=[c for c,sz in cnt.items() if sz>=gm]         # balance 閘：核心群 ≥ grp_min
        out_idx=[j for j in range(n) if cnt[cl[j]]<gm]          # 小群(含 3-read 尾巴)→離群
        kc=len(core_lbls)
        if not(2<=kc<=kmax): continue                          # k_core 上限 + 至少 2 群
        if len(out_idx)/n>max_out_frac: continue               # 離群比例上限
        cm=np.array([cnt[cl[j]]>=gm for j in range(n)]); ci=np.where(cm)[0]
        if len(ci)<kc*MIN_SZ: continue
        try: sep=float(silhouette_score(s[np.ix_(ci,ci)],cl[cm],metric="precomputed"))
        except Exception: continue
        if sep<sil_min: continue                               # 分離度閘
        prev=res.get(kc)
        if prev is None or sep>prev["sep"]:
            res[kc]={"k_core":kc,"t":round(t,6),"labels":cl.tolist(),"out_idx":out_idx,
                     "core_sizes":sorted([cnt[c] for c in core_lbls],reverse=True),
                     "n_out":len(out_idx),"sep":sep}
    return Z,res

def align_resolution(labels, out_idx, axis_lab):
    """對核心 read 算 CramerV(cut vs 軸)。回 (V,p,e)。"""
    core=[i for i in range(len(labels)) if i not in set(out_idx)]
    g=[axis_lab[i] for i in core if axis_lab[i] is not None]
    c=[labels[i] for i in core if axis_lab[i] is not None]
    if len(set(g))<2 or len(set(c))<2: return None,None,None
    return cramerv(g,c)

EXC_MIN=0.10   # null-excess 門檻（真實結構 = 扣 within-1-group null 後）

def analyze_locus(sub, P, hp_lab, alle_lab, rng):
    """三閘切群: balance(new_cut 內) + null-excess(真實?) + alignment(歸因?)。
    每解析度分類: CONFIRMED(真實&可靠對齊) / REAL_NOVEL(真實但不對齊=subclone 候選) / NOT_REAL(扣null後弱→不切)。
    邊緣群: 離群 read 夠多(≥grp_min) → 標 edge_group。輸出 per-read 群歸屬。"""
    n=sub.shape[0]; gm=grp_min(n)
    old=old_cut(sub); Z,res=new_cut(sub)
    axes={"hp":["HP1" if l in("1","1-1") else "HP2" for l in hp_lab],
          "carrier":["G" if l in("1","2") else "C" for l in hp_lab],
          "allele":[a if a in("REF","ALT") else None for a in alle_lab]}
    resolutions=[]
    for kc in sorted(res):
        r=res[kc]; al={}; reliable_axes=[]
        for ax,lab in axes.items():
            V,p,e=align_resolution(r["labels"],r["out_idx"],lab)
            rel=bool(V is not None and V>=0.3 and (p if p is not None else 1)<0.05 and (e or 0)>=5)
            al[ax]={"V":round(V,3) if V is not None else None,"p":round(p,4) if p is not None else None,
                    "e":round(e,1) if e is not None else None,"reliable":rel}
            if rel: reliable_axes.append(ax)
        exc=stab_excess(sub,P,kc,rng)
        real=(exc is not None and exc>=EXC_MIN)
        conf=("CONFIRMED" if (real and reliable_axes) else
              "REAL_NOVEL" if real else "NOT_REAL")
        n_out=r["n_out"]; edge=bool(n_out>=gm)   # 邊緣群: 離群夠多→成一組
        resolutions.append({"k_core":kc,"core_sizes":r["core_sizes"],"n_outliers":n_out,
                            "edge_group":edge,"separation":round(r["sep"],3) if r["sep"] is not None else None,
                            "stab_excess":exc,"real":real,"confidence":conf,
                            "aligned_axes":reliable_axes,"align":al,"labels":r["labels"],"out_idx":r["out_idx"]})
    confirmed=[r for r in resolutions if r["confidence"]=="CONFIRMED"]
    novel=[r for r in resolutions if r["confidence"]=="REAL_NOVEL"]
    real_res=confirmed+novel
    def perread_of(res):
        if not res: return None
        lbl=res["labels"]; out=set(res["out_idx"]); gtype=res["confidence"]; cmap={}; cid=0
        for i in range(n):
            if i in out: continue
            if lbl[i] not in cmap: cmap[lbl[i]]=cid; cid+=1
        pr=[]
        for i in range(n):
            if i in out: pr.append({"grp":"edge" if res["edge_group"] else "outlier","type":"edge" if res["edge_group"] else "outlier"})
            else: pr.append({"grp":cmap[lbl[i]],"type":"core_confirmed" if gtype=="CONFIRMED" else "core_novel"})
        return pr
    # 雙輸出: coarse=最粗 CONFIRMED 骨幹; fine=最細全真實結構(confirmed∪novel)
    coarse=min(confirmed,key=lambda r:r["k_core"]) if confirmed else None
    fine=max(real_res,key=lambda r:r["k_core"]) if real_res else None
    old_k=old["k_groups"] if old["status"]=="SPLIT" else 0
    delta=[]
    if old["status"]=="CANT_SPLIT" and real_res: delta.append("RESCUED_cant_split→split")
    if len(confirmed)>=2: delta.append(f"MULTIRES_confirmed={[r['k_core'] for r in confirmed]}")
    if novel: delta.append(f"NOVEL_k={[r['k_core'] for r in novel]}(subclone候選)")
    return {"old":old,
            "coarse_k":(coarse["k_core"] if coarse else 0),"coarse_confidence":(coarse["confidence"] if coarse else "NONE"),
            "fine_k":(fine["k_core"] if fine else 0),"fine_confidence":(fine["confidence"] if fine else "NO_CLEAR_SPLIT"),
            "primary_k":(fine["k_core"] if fine else 0),"primary_confidence":(fine["confidence"] if fine else "NO_CLEAR_SPLIT"),
            "n_confirmed":len(confirmed),"n_novel":len(novel),
            "confirmed_ks":[r["k_core"] for r in confirmed],"novel_ks":[r["k_core"] for r in novel],
            "new_resolutions":resolutions,"perread_coarse":perread_of(coarse),"perread_fine":perread_of(fine),
            "perread":perread_of(fine),"delta":delta}

# ================= 直跑：對 29 example =================
if __name__=="__main__":
    WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
    A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUT=f"{WT}/output/_kprofile_heatmap"
    dirmap={}
    for mp in glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True):
        rd=mp.rsplit("/distance/",1)[0]
        for part in rd.split("/"):
            if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd
    items=json.load(open(f"{A}/kprofile_examples.json"))["items"]
    rng=np.random.default_rng(20260622); out=[]
    for L in items:
        key=f"{L['chrom']}_{L['pos']}"; rd=dirmap.get(key)
        if not rd: continue
        reads={r["read_id"]:r for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
        dids,D=loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={r:i for i,r in enumerate(dids)}
        mids,Me=loadm(f"{rd}/methylation/methylation.csv"); mi={r:i for i,r in enumerate(mids)}
        def istum(t): return str(t) in("1","true","True")
        kids=[r for r in dids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP and r in mi]
        if len(kids)<MIN_SZ*2: continue
        sub=D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]; kp=peel(sub)
        if len(kp)<MIN_SZ*2: continue
        kids=[kids[i] for i in kp]; sub=sub[np.ix_(kp,kp)]
        P=np.array([Me[mi[r]] for r in kids])
        hp_lab=[LABMAP[reads[r]["hp"]] for r in kids]; alle_lab=[reads[r]["alt_support"] for r in kids]
        a=analyze_locus(sub,P,hp_lab,alle_lab,rng)
        a.update({"key":key,"group":L["group"],"n":len(kids),"primary_axis":L.get("primary_axis"),
                  "read_ids":kids,"hp_lab":hp_lab,"alle_lab":alle_lab})  # for plotting reorder
        out.append(a)
        print(f"  {L['group']:18s} {key:16s} n={a['n']:3d} OLD k={a['old']['k_groups']} "
              f"→ COARSE k={a['coarse_k']}({a['coarse_confidence']}) | FINE k={a['fine_k']}({a['fine_confidence']}) "
              f"confirmed={a['confirmed_ks']} novel={a['novel_ks']}",flush=True)
    # summary
    multires=[r for r in out if r["n_confirmed"]>=2]
    has_novel=[r for r in out if r["n_novel"]>=1]
    has_edge=[r for r in out if any(res["edge_group"] for res in r["new_resolutions"])]
    from collections import Counter as _C
    summ={"n_loci":len(out),"params":{"min_sz":MIN_SZ,"grp_min":"=MIN_SZ=3(業界caller最低)","max_out_frac":MAX_OUT_FRAC,"kmax":KMAX,
            "gap_mult":GAP_MULT,"exc_min":EXC_MIN,"seed":20260622,
            "三閘":"balance(群≥3) + null-excess(扣 within-1-group null≥0.10=真實/明顯) + alignment(歸因)",
            "confidence":"CONFIRMED(真實&可靠對齊germline) / REAL_NOVEL(真實但不對齊=subclone候選) / NOT_REAL(不切)",
            "edge_group":"離群 read ≥3 → 歸成邊緣一組","dual_output":"coarse=最粗CONFIRMED骨幹 / fine=最細全真實結構(含novel)"},
          "FM3_multiresolution_confirmed(≥2)":len(multires),"has_REAL_NOVEL(subclone候選)":len(has_novel),
          "has_edge_group":len(has_edge),
          "primary_confidence_dist":dict(_C(r["primary_confidence"] for r in out)),
          "loci":out}
    json.dump(summ,open(f"{A}/cluster_redesign.json","w"),indent=1)
    print(f"\nDONE n={len(out)}  multi-res(≥2 confirmed)={len(multires)}  has_REAL_NOVEL={len(has_novel)}  has_edge={len(has_edge)}")
    print("primary confidence:",dict(_C(r['primary_confidence'] for r in out)))
