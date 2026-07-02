#!/usr/bin/env python3
"""切割穩定性 compute (per-locus): 對 kprofile_examples.json 中已有 binary 落檔的代表位點，
量化「同一刀在擾動下是否重現」+ chance-corrected null，避開「穩定≠正確」陷阱。

三軸（CLAUDE 主軸 §3）:
  ① 存在性  = 取自 ksweep 切群成功(本檔位點皆 splittable) + per-k silhouette
  ② 穩定性  = clusterboot Jaccard(80% subsample, Hennig) + consensus 矩陣 PAC
               對「within-1-group null」(逐 CpG 欄內非 NaN 重排 → 重算 BERNOULLI) 校正
  ③ 非循環  = CramerV(cut vs primary a-priori 軸) 對「label-shuffle null」校正

null（雙 null，user 06-22 指定由資料特性選 → 選雙 null 因 read-內甲基相關是 81% A-path 失敗根因）:
  - stability null: within-1-group(打散 read 間共結構、保 per-CpG 邊際+per-read coverage) → 假穩 ceiling
  - alignment null: label-shuffle → 假對齊 ceiling

切群法與 ksweep_wg_v2.py 完全一致: peel → UPGMA(average) → fcluster(maxclust), MIN_SZ=3。
BERNOULLI 重算與 src/core/DistanceMatrix.cpp:254 calculate_bernoulli 一致(C_min=3)，
腳本內含對既有 matrix.csv 的自檢(max abs diff)。

輸出: kprofile_stability.json（每位點 per-k 穩定性 + null + 對齊 + 3 軸 verdict）。
"""
import os, csv, glob, json, sys
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.metrics import silhouette_score
from collections import Counter

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = f"{WT}/output/_kprofile_heatmap"
MIN_SZ = 3
CMIN = 3                      # min_common_coverage (Config.hpp:37)
SUBFRAC = 0.80               # clusterboot subsample fraction
EXCESS_MIN = 0.10            # meaningful excess-over-null margin (convention)
B = 300                      # subsamples per locus per k
RNULL = 40                   # within-1-group null draws
PSHUF = 500                  # label-shuffle permutations
SEED = 20260622
LABMAP = {"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}
rng = np.random.default_rng(SEED)


def _f(x):
    x = str(x).strip()
    return np.nan if x in ("","NA","nan","NaN") else float(x)

def loadm(p, want_cols=False):
    r = open(p).read().strip().split("\n"); ids=[]; M=[]; cols=r[0].split(",")[1:]
    for ln in r[1:]:
        q = ln.split(","); ids.append(q[0]); M.append([_f(x) for x in q[1:]])
    return (ids, np.array(M), cols) if want_cols else (ids, np.array(M))

def peel(s):
    idx = list(range(s.shape[0]))
    while True:
        S = s[np.ix_(idx,idx)]; b=(S<0)|np.isnan(S); np.fill_diagonal(b,False)
        if not b.any(): return idx
        idx.remove(idx[int(np.argmax(b.sum(1)))])
        if len(idx) < MIN_SZ*2: return idx

def cramerv(groups, clusters):
    g=sorted(set(groups)); c=sorted(set(clusters))
    if len(g)<2 or len(c)<2: return None,None,None
    tab=np.zeros((len(g),len(c))); gi={x:i for i,x in enumerate(g)}; ci={x:i for i,x in enumerate(c)}
    for a,b in zip(groups,clusters): tab[gi[a],ci[b]]+=1
    tab=tab[tab.sum(1)>0][:,tab.sum(0)>0]
    if tab.shape[0]<2 or tab.shape[1]<2: return None,None,None
    try: chi2,p,dof,exp=chi2_contingency(tab)
    except Exception: return None,None,None
    nn=tab.sum(); v=float(np.sqrt(chi2/(nn*(min(tab.shape)-1)))) if nn>0 else None
    return v, float(p), float(exp.min())

def bernoulli_dist(P):
    """向量化重算 BERNOULLI 距離矩陣 (n×n)，NaN=missing。
    與 calculate_bernoulli 等價: w=2|p-.5|, delta=p+q-2pq, D=sum(w_i w_j delta)/sum(w_i w_j)。
    common<CMIN 或 sum_w<1e-9 → -1（invalid，對應 SKIP）。"""
    n, C = P.shape
    Mk = (~np.isnan(P)).astype(float)            # valid mask
    Pz = np.where(np.isnan(P), 0.0, P)           # p, NaN→0
    Aw = 2.0*np.abs(Pz-0.5)*Mk                   # w (0 where invalid)
    Bw = Aw*Pz                                   # w*p
    sw  = Aw @ Aw.T                              # sum_c w_i w_j
    swd = Bw @ Aw.T + Aw @ Bw.T - 2.0*(Bw @ Bw.T)  # sum_c w_i w_j (p_i+p_j-2 p_i p_j)
    common = Mk @ Mk.T
    with np.errstate(divide="ignore", invalid="ignore"):
        D = swd/sw
    D[(common < CMIN) | (sw < 1e-9)] = -1.0
    np.fill_diagonal(D, 0.0)
    return D

def clean_sub(D):
    """peel + 清理成可 squareform 的對稱矩陣，回 (sub, keep_idx) 或 (None,None)。"""
    keep = peel(D)
    if len(keep) < MIN_SZ*2: return None, None
    sub = D[np.ix_(keep,keep)].copy()
    if ((sub<0)|np.isnan(sub)).any(): return None, None
    np.fill_diagonal(sub,0.0); sub=np.maximum(sub,sub.T)
    return sub, keep

def cut(sub, k):
    Z = linkage(squareform(sub, checks=False), method="average")
    cl = fcluster(Z, k, criterion="maxclust")
    return cl

def kmax_of(n): return min(5, n//MIN_SZ)

def clusterboot(sub, cl0, k, B, rng):
    """Hennig cluster-wise Jaccard via 80% subsample。回 (per_cluster_jaccard, mean, consensus_PAC)。"""
    n = sub.shape[0]
    m = max(MIN_SZ*2, int(round(SUBFRAC*n)))
    orig_clusters = [set(np.where(cl0==c)[0]) for c in sorted(set(cl0))]
    jac_acc = [[] for _ in orig_clusters]
    co = np.zeros((n,n)); both = np.zeros((n,n))
    for _ in range(B):
        samp = np.sort(rng.choice(n, size=m, replace=False))
        ss = sub[np.ix_(samp,samp)]
        if ((ss<0)|np.isnan(ss)).any(): continue
        try: cl = cut(ss, k)
        except Exception: continue
        sampset = set(samp.tolist())
        # consensus co-assignment
        both[np.ix_(samp,samp)] += 1
        for ci in sorted(set(cl)):
            members = samp[cl==ci]
            co[np.ix_(members,members)] += 1
        # cluster-wise jaccard: each original cluster vs best bootstrap cluster (restricted to sampled)
        boot_clusters = [set(samp[cl==ci].tolist()) for ci in sorted(set(cl))]
        for oi, oc in enumerate(orig_clusters):
            oc_s = oc & sampset
            if len(oc_s) < 2: continue
            best = max((len(oc_s & bc)/len(oc_s | bc) for bc in boot_clusters), default=0.0)
            jac_acc[oi].append(best)
    per_cluster = [float(np.mean(j)) if j else 0.0 for j in jac_acc]
    mean_j = float(np.mean(per_cluster)) if per_cluster else 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        M = np.where(both>0, co/both, np.nan)
    offdiag = M[np.triu_indices(n,1)]
    valid = offdiag[~np.isnan(offdiag)]
    PAC = float(np.mean((valid>0.1)&(valid<0.9))) if valid.size else 1.0
    return per_cluster, mean_j, PAC, M

def null_within1group(P, k, RNULL, B, rng):
    """within-1-group null: 逐 CpG 欄內非 NaN 重排 → 重算 BERNOULLI → 同 clusterboot。
    回 null mean-jaccard 分布(每個 null 一個 mean)。"""
    n, C = P.shape
    null_means = []
    for _ in range(RNULL):
        Pn = P.copy()
        for c in range(C):
            col = Pn[:,c]; vi = np.where(~np.isnan(col))[0]
            if vi.size > 1:
                perm = rng.permutation(vi)
                Pn[vi,c] = col[perm]
        Dn = bernoulli_dist(Pn)
        sub, keep = clean_sub(Dn)
        if sub is None or sub.shape[0] < k*MIN_SZ:
            null_means.append(0.0); continue
        try: cl0 = cut(sub, k)
        except Exception:
            null_means.append(0.0); continue
        if len(set(cl0)) < k or min(Counter(cl0).values()) < MIN_SZ:
            null_means.append(0.0); continue
        _, mean_j, _, _ = clusterboot(sub, cl0, k, max(60,B//4), rng)
        null_means.append(mean_j)
    return null_means

def align_test(cl0, labels, PSHUF, rng):
    """CramerV(cut, primary 軸) + label-shuffle null p。"""
    V, p_chi, e = cramerv(labels, list(cl0))
    if V is None: return None, None, None, None
    null=[]
    lab=np.array(labels)
    for _ in range(PSHUF):
        Vs,_,_ = cramerv(list(rng.permutation(lab)), list(cl0))
        null.append(Vs if Vs is not None else 0.0)
    null=np.array(null); shuf_p=float((null>=V).mean())
    return float(V), float(e if e else 0), shuf_p, float(np.percentile(null,95))


# ---- locate binary outputs (same dirmap as plot_kprofile_heatmaps) ----
dirmap={}
for mp in glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd

items = json.load(open(f"{A}/kprofile_examples.json"))["items"]
results=[]; selfcheck_max=0.0
for L in items:
    key=f"{L['chrom']}_{L['pos']}"; rd=dirmap.get(key)
    if not rd: continue
    reads={r["read_id"]:r for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D = loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={r:i for i,r in enumerate(dids)}
    mids,Me,cpgs = loadm(f"{rd}/methylation/methylation.csv", want_cols=True); mi={r:i for i,r in enumerate(mids)}
    def istum(t): return str(t) in ("1","true","True")
    kids=[r for r in dids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP and r in mi]
    if len(kids) < MIN_SZ*2: continue
    # tumor sub-distance + peel
    subD = D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]
    sub, keep = clean_sub(subD)
    if sub is None: continue
    kids=[kids[i] for i in keep]; n=len(kids)
    P = np.array([Me[mi[r]] for r in kids])           # read×CpG raw p (tumor, peeled)
    # --- self-check: recompute BERNOULLI from P, compare to stored sub ---
    Dchk = bernoulli_dist(P)
    valid = (sub>=0)&(Dchk>=0);
    if valid.any():
        selfcheck_max=max(selfcheck_max, float(np.abs(sub[valid]-Dchk[valid]).max()))
    # labels
    lab=[LABMAP[reads[r]["hp"]] for r in kids]
    axmap={"hp":["HP1" if l in("1","1-1") else "HP2" for l in lab],
           "carrier":["G" if l in("1","2") else "C" for l in lab],
           "allele":[reads[r]["alt_support"] for r in kids]}
    prim = L.get("primary_axis") or "hp"
    if prim=="allele":
        keepA=[i for i in range(n) if axmap["allele"][i] in("REF","ALT")]
        prim_lab=[axmap["allele"][i] for i in keepA]; prim_idx=keepA
    else:
        prim_lab=axmap[prim]; prim_idx=list(range(n))
    # per-k stability (real + null) + alignment
    per_k=[]
    for k in range(2, kmax_of(n)+1):
        try: cl0=cut(sub,k)
        except Exception: continue
        if len(set(cl0))<k or min(Counter(cl0).values())<MIN_SZ: continue
        try: sil=float(silhouette_score(sub,cl0,metric="precomputed"))
        except Exception: sil=None
        pc, meanj, PAC, _ = clusterboot(sub, cl0, k, B, rng)
        nulls = null_within1group(P, k, RNULL, B, rng)
        null95 = float(np.percentile(nulls,95)); null_mean=float(np.mean(nulls))
        # alignment at this k
        if prim=="allele":
            cl_p=[cl0[i] for i in prim_idx]
            Va,ea,shp,_ = align_test(np.array(cl_p), prim_lab, PSHUF, rng)
        else:
            Va,ea,shp,_ = align_test(cl0, prim_lab, PSHUF, rng)
        ps = bool(meanj>null95 and (meanj-null95)>=EXCESS_MIN and meanj>=0.6)
        pa = bool(Va is not None and Va>=0.3 and shp is not None and shp<0.05 and (ea or 0)>=5)
        per_k.append({"k":k,"sil":sil,"n_clusters":len(set(cl0)),
            "jac_real_mean":round(meanj,3),"jac_real_per_cluster":[round(x,3) for x in pc],
            "jac_null_mean":round(null_mean,3),"jac_null_p95":round(null95,3),
            "stab_excess":round(meanj-null95,3),"pass_stability":ps,
            "PAC":round(PAC,3),
            "align_V":round(Va,3) if Va is not None else None,"align_e":round(ea,1) if ea is not None else None,
            "align_shuffle_p":round(shp,4) if shp is not None else None,
            "pass_align":pa,"pass_both":bool(ps and pa)})
    if not per_k: continue
    # coherent-best k = 同一 k 同時過 stab+align，取 excess 最大者；無則取整體 max-excess
    coherent = [d for d in per_k if d["pass_both"]]
    head_k = (max(coherent, key=lambda d:d["stab_excess"]) if coherent
              else max(per_k, key=lambda d:d["stab_excess"]))
    excess_k = max(per_k, key=lambda d:d["stab_excess"])
    any_pass_stab = any(d["pass_stability"] for d in per_k)
    any_pass_align = any(d["pass_align"] for d in per_k)
    any_both = bool(coherent)
    verdict = ("GOOD_CUT" if any_both else
               "STABLE_BUT_UNALIGNED" if any_pass_stab and not any_pass_align else
               "ALIGNED_BUT_UNSTABLE" if any_pass_align and not any_pass_stab else
               "STABLE+ALIGNED_DIFF_K" if any_pass_stab and any_pass_align else
               "NEITHER")
    results.append({"chrom":L["chrom"],"pos":L["pos"],"group":L["group"],"n":n,
        "primary_axis":prim,"meaningful_ks":L.get("meaningful_ks"),
        "best_k_silhouette":L.get("best_k"),"headline_k":head_k["k"],"max_excess_k":excess_k["k"],
        "n_reliable_align_k":sum(1 for d in per_k if (d["align_e"] or 0)>=5),
        "per_k":per_k,"verdict":verdict,
        "headline":{"k":head_k["k"],"jac_real":head_k["jac_real_mean"],"jac_null_p95":head_k["jac_null_p95"],
                    "stab_excess":head_k["stab_excess"],
                    "align_V":head_k["align_V"],"align_e":head_k["align_e"],"align_shuffle_p":head_k["align_shuffle_p"]}})
    print(f"  {L['group']:18s} {key:18s} n={n} verdict={verdict:21s} "
          f"k={head_k['k']} jac={head_k['jac_real_mean']} null95={head_k['jac_null_p95']} "
          f"exc={head_k['stab_excess']:+.3f} V={head_k['align_V']} e={head_k['align_e']}", flush=True)

summary={"n_loci":len(results),
    "params":{"subsample_frac":SUBFRAC,"B_subsamples":B,"R_null":RNULL,"P_shuffle":PSHUF,
              "min_sz":MIN_SZ,"cmin":CMIN,"excess_min":EXCESS_MIN,"seed":SEED,
              "null_stability":"within-1-group (per-CpG within-column non-NaN permute → recompute BERNOULLI)",
              "null_alignment":"label-shuffle","clustering":"peel→UPGMA(average)→fcluster maxclust",
              "stability_pass":"jac_real_mean > null_p95 AND excess>=0.10 AND jac_real_mean>=0.6 (Hennig)",
              "align_pass":"CramerV>=0.3 AND shuffle_p<0.05 AND Cochran_e>=5 (e>=5 = 可靠性閘)",
              "verdict_rule":"GOOD_CUT 需同一 k 同時過 stab+align(coherent); headline_k=coherent 中 max-excess(無則整體 max-excess)"},
    "bernoulli_selfcheck_max_absdiff":round(selfcheck_max,6),
    "verdict_counts":dict(Counter(r["verdict"] for r in results)),
    "by_group":{g:dict(Counter(r["verdict"] for r in results if r["group"]==g))
                for g in sorted(set(r["group"] for r in results))},
    "loci":results}
json.dump(summary, open(f"{A}/kprofile_stability.json","w"), indent=1)
print(f"\nDONE n={len(results)}  BERNOULLI self-check max|Δ|={selfcheck_max:.2e}  "
      f"(0≈完全一致 → 重算正確)")
print("verdicts:", summary["verdict_counts"])
print("by group:", json.dumps(summary["by_group"], ensure_ascii=False))
