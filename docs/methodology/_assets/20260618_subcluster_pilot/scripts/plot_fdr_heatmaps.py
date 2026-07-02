#!/usr/bin/env python3
"""每 FDR-顯著位點畫 dual-panel: 左甲基 read×CpG + 右 read×read 距離 (SoT ism_heatmap_std)。
側欄 = cut(align-k, k2) + 標準 HP|ALT|Strand。🔴 必含距離熱圖一起觀察。
取代 06-18 單-panel 版（漏距離 panel + 自定色盤）。輸出 fdr_heatmap_index.json(同格式供 workstation)。"""
import os, glob, csv, json, sys
import numpy as np
sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts")
import ism_heatmap_std as H
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
MIN_SZ=3; LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
OUT=f"{WT}/output/_fdr_loci_heatmap"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
FIGS=f"{A}/figs_fdr"; os.makedirs(FIGS,exist_ok=True)
def _f(x):
    x=str(x).strip(); return np.nan if x in("","NA","nan","NaN") else float(x)
def loadm(p,cols=False):
    r=open(p).read().strip().split("\n"); ids=[];M=[];c=r[0].split(",")[1:]
    for ln in r[1:]:
        q=ln.split(","); ids.append(q[0]); M.append([_f(x) for x in q[1:]])
    return (ids,np.array(M),c) if cols else (ids,np.array(M))
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
loci=json.load(open(f"{A}/ksweep_fdr_loci.json"))["loci"]; done=[]
for L in loci:
    key=f"{L['chrom']}_{L['pos']}"; rd=dirmap.get(key)
    if not rd: continue
    reads={r["read_id"]:r for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={r:i for i,r in enumerate(dids)}
    mids,Me,cpgs=loadm(f"{rd}/methylation/methylation.csv",cols=True); mi={r:i for i,r in enumerate(mids)}
    def istum(t): return str(t) in("1","true","True")
    kids=[r for r in dids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP and r in mi]
    if len(kids)<MIN_SZ*2: continue
    sub=D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]; kp=peel(sub)
    if len(kp)<MIN_SZ*2: continue
    kids=[kids[i] for i in kp]; sub=sub[np.ix_(kp,kp)].copy(); np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T)
    Z=linkage(squareform(sub,checks=False),method="average"); n=len(kids)
    prim=L["primary"]; ak=L["axes"][prim]["align_k"]
    cuts={2:list(fcluster(Z,2,'maxclust'))}
    if ak<=min(5,n//MIN_SZ): cuts[ak]=list(fcluster(Z,ak,'maxclust'))
    sortk=ak if ak in cuts else 2
    pa={"hp":["HP1" if LABMAP[reads[r]["hp"]] in("1","1-1") else "HP2" for r in kids],
        "carrier":["G" if LABMAP[reads[r]["hp"]] in("1","2") else "C" for r in kids],
        "allele":[reads[r]["alt_support"] for r in kids]}.get(prim)
    order=sorted(range(n),key=lambda i:(cuts[sortk][i],pa[i]))
    kids=[kids[i] for i in order]; sub=sub[np.ix_(order,order)]
    meth=np.array([Me[mi[r]] for r in kids])
    cutsO={k:[cuts[k][i] for i in order] for k in cuts}
    def clhex(seq):
        cl=sorted(set(seq)); cm={c:H.CLUSTER_COL[i%len(H.CLUSTER_COL)] for i,c in enumerate(cl)}
        return [cm[v] for v in seq]
    sidebars=[(f"k{sortk}*",clhex(cutsO[sortk]))]
    sidebars+=H.sidebar_specs({r:reads[r] for r in kids},kids,cluster_of=None,include_tn=False,include_strand=True)
    if 2!=sortk and 2 in cutsO: sidebars+=[("k2",clhex(cutsO[2]))]
    a=L["axes"][prim]
    title=(f"[FDR {prim.upper()}] {L['chrom']}:{L['pos']} n={n} | sil k=2(V={a['V_sil']}) → align k={ak}(V={a['V_align']}) "
           f"| dV={a['dV']} p_bonf={a['p_bonf']:.0e}")
    fn=f"{FIGS}/fdr_{key}.png"
    H.mpl_dual_panel(meth,sub,sidebars,[int(c) for c in cpgs],int(L["pos"]),title,fn,n_cluster=len(set(cutsO[sortk])),dpi=90)
    done.append({"key":key,"chrom":L["chrom"],"pos":L["pos"],"n":n,"primary":prim,
        "axes":L["axes"],"png":f"figs_fdr/fdr_{key}.png"})
json.dump({"n":len(done),"items":done},open(f"{A}/fdr_heatmap_index.json","w"),indent=1)
print(f"plotted {len(done)}/{len(loci)} dual-panel(甲基+距離) -> {FIGS}")
