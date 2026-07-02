#!/usr/bin/env python3
"""每 k-profile 代表位點畫 dual-panel: 左甲基 read×CpG + 右 read×read 距離 (SoT scripts/ism_heatmap_std.py)。
側欄 = 多-k cut(k2/k3/k4) + 標準 HP|ALT|Strand (HP 4 色已含 germline/carrier)。
🔴 必含距離熱圖一起觀察 (feedback_ism_case_heatmap_standard_sidebars line 31)。"""
import os, glob, csv, json, sys
import numpy as np
sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts")
import ism_heatmap_std as H
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from collections import Counter
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
OUT=f"{WT}/output/_kprofile_heatmap"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
FIGS=f"{A}/figs_kprofile"; os.makedirs(FIGS,exist_ok=True)
MIN_SZ=3; LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}
def _f(x):
    x=str(x).strip(); return np.nan if x in("","NA","nan","NaN") else float(x)
def loadm(p,want_cols=False):
    r=open(p).read().strip().split("\n"); ids=[];M=[]; cols=r[0].split(",")[1:]
    for ln in r[1:]:
        q=ln.split(","); ids.append(q[0]); M.append([_f(x) for x in q[1:]])
    return (ids,np.array(M),cols) if want_cols else (ids,np.array(M))
def peel(s):
    idx=list(range(s.shape[0]))
    while True:
        S=s[np.ix_(idx,idx)];b=(S<0)|np.isnan(S);np.fill_diagonal(b,False)
        if not b.any(): return idx
        idx.remove(idx[int(np.argmax(b.sum(1)))])
        if len(idx)<MIN_SZ*2: return idx

dirmap={}
for mp in glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True):
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd
items=json.load(open(f"{A}/kprofile_examples.json"))["items"]; done=[]
for L in items:
    key=f"{L['chrom']}_{L['pos']}"; rd=dirmap.get(key)
    if not rd: continue
    reads={r["read_id"]:r for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={r:i for i,r in enumerate(dids)}
    mids,Me,cpgs=loadm(f"{rd}/methylation/methylation.csv",want_cols=True); mi={r:i for i,r in enumerate(mids)}
    def istum(t): return str(t) in("1","true","True")
    kids=[r for r in dids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP and r in mi]
    if len(kids)<MIN_SZ*2: continue
    sub=D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]; kp=peel(sub)
    if len(kp)<MIN_SZ*2: continue
    kids=[kids[i] for i in kp]; sub=sub[np.ix_(kp,kp)].copy()
    np.fill_diagonal(sub,0);sub=np.maximum(sub,sub.T)
    Z=linkage(squareform(sub,checks=False),method="average"); n=len(kids)
    kmax=min(4,n//MIN_SZ)
    cuts={k:list(fcluster(Z,k,'maxclust')) for k in (2,3,4) if k<=kmax}
    sortk=max(cuts)
    lab=[LABMAP[reads[r]["hp"]] for r in kids]
    prim={"hp":["HP1" if l in("1","1-1") else "HP2" for l in lab],
          "carrier":["G" if l in("1","2") else "C" for l in lab],
          "allele":[reads[r]["alt_support"] for r in kids]}.get(L["primary_axis"],lab)
    order=sorted(range(n),key=lambda i:(cuts[sortk][i],prim[i]))
    # 排序後重排
    kids=[kids[i] for i in order]
    sub=sub[np.ix_(order,order)]
    meth=np.array([Me[mi[r]] for r in kids])
    cutsO={k:[cuts[k][i] for i in order] for k in cuts}
    # sidebars: 主 cut(sortk, sidebars[0] 供分隔線) + 標準 HP|ALT|Strand + 其餘 k cut
    def clhex(seq):
        cl=sorted(set(seq)); cm={c:H.CLUSTER_COL[i%len(H.CLUSTER_COL)] for i,c in enumerate(cl)}
        return [cm[v] for v in seq]
    sidebars=[(f"k{sortk}",clhex(cutsO[sortk]))]
    rd_dict={r:reads[r] for r in kids}
    sidebars+=H.sidebar_specs(rd_dict,kids,cluster_of=None,include_tn=False,include_strand=True)  # HP|ALT|Strand
    sidebars+=[(f"k{k}",clhex(cutsO[k])) for k in (2,3,4) if k in cutsO and k!=sortk]
    cpg_int=[int(c) for c in cpgs]
    mkax="; ".join(f"k{k}→{L['mk_axes'][str(k)][0]}({L['mk_axes'][str(k)][1]})" for k in sorted(map(int,L['mk_axes'].keys()))) if L['mk_axes'] else "—"
    title=(f"[{L['group']}] {L['chrom']}:{L['pos']} n={n} | meaningful_ks={L['meaningful_ks']} {mkax} | "
           f"sil_m={L['sil_margin']} align_m={L['align_margin']}")
    fn=f"{FIGS}/kp_{L['group']}_{key}.png"
    H.mpl_dual_panel(meth,sub,sidebars,cpg_int,int(L["pos"]),title,fn,n_cluster=len(set(cutsO[sortk])),dpi=92)
    done.append({"key":key,"group":L["group"],"chrom":L["chrom"],"pos":L["pos"],"n":n,
        "meaningful_ks":L["meaningful_ks"],"mk_axes":L["mk_axes"],"sil_margin":L["sil_margin"],
        "align_margin":L["align_margin"],"primary_axis":L["primary_axis"],"png":f"figs_kprofile/kp_{L['group']}_{key}.png"})
json.dump({"n":len(done),"items":done},open(f"{A}/kprofile_heatmap_index.json","w"),indent=1)
print(f"plotted {len(done)}/{len(items)} dual-panel(甲基+距離) -> {FIGS}")
print("by group:",dict(Counter(d['group'] for d in done)))
