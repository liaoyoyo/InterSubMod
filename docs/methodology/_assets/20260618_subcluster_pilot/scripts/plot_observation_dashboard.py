#!/usr/bin/env python3
"""觀察儀表板格式: 標準 ISM dual-panel(原始甲基 read×CpG + 距離 read×read)+ 新切群側欄。
側欄(左→右,熱圖在右): fine群 | coarse群 | HP | ALT | Strand (tumor-only,無 T/N)。
重用 scripts/ism_heatmap_std.mpl_dual_panel(顏色/佈局=既有觀察標準)。讀 cluster_redesign.json 的 per-read 群歸屬。"""
import os, csv, glob, json, sys
import numpy as np
sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0,os.path.dirname(__file__)); import cluster_redesign as CR
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUT=f"{WT}/output/_kprofile_heatmap"
FIGS=f"{A}/figs_dashboard"; os.makedirs(FIGS,exist_ok=True)
EDGE="#111111"; OUTL="#c9c9c9"
dirmap={}
for mp in glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True):
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd
red={r["key"]:r for r in json.load(open(f"{A}/cluster_redesign.json"))["loci"]}

def grpcol(g):  # group id(int)→cluster color; edge/outlier→特殊
    if g=="edge": return EDGE
    if g=="outlier": return OUTL
    return H.CLUSTER_COL[int(g)%len(H.CLUSTER_COL)]

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
def plot(key):
    r=red[key]; rd=dirmap[key]
    reads={x["read_id"]:x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={x:i for i,x in enumerate(dids)}
    rows=open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); cpgs=[int(c) for c in rows[0].split(",")[1:]]
    mi={}; M=[]
    for j,ln in enumerate(rows[1:]):
        q=ln.split(","); mi[q[0]]=j; M.append([np.nan if v in("","NA","nan","NaN") else float(v) for v in q[1:]])
    M=np.array(M)
    ids=[x for x in r["read_ids"] if x in di and x in mi]
    coarse={x:g for x,g in zip(r["read_ids"],[p["grp"] for p in (r["perread_coarse"] or [{"grp":"outlier"}]*len(r["read_ids"]))])}
    fine={x:g for x,g in zip(r["read_ids"],[p["grp"] for p in (r["perread_fine"] or [{"grp":"outlier"}]*len(r["read_ids"]))])}
    sub=D[np.ix_([di[x] for x in ids],[di[x] for x in ids])].copy()
    Z,_=CR.linkZ(sub)
    h=np.sort(Z[:,2]); fk=r.get("fine_k",2) or 2; tf=(h[-(fk-1)]+h[-fk])/2 if 2<=fk<=len(h) else None
    dn=dendrogram(Z,orientation="left",no_plot=True)
    order=dn["leaves"][::-1]                      # 上→下(imshow row0=top)
    ids_o=[ids[i] for i in order]
    meth=np.array([M[mi[x]] for x in ids_o]); dist=sub[np.ix_(order,order)].copy()
    np.fill_diagonal(dist,0); dist[dist<0]=np.nan
    sb=[("fine",[grpcol(fine.get(x,"outlier")) for x in ids_o]),
        ("coarse",[grpcol(coarse.get(x,"outlier")) for x in ids_o])]
    sb+=H.sidebar_specs({x:reads[x] for x in ids_o}, ids_o, cluster_of=None, include_tn=False, include_strand=True)
    nclust=len({fine.get(x) for x in ids_o if fine.get(x) not in("edge","outlier")})
    # 組合面板: dendro | sb | meth | gap | sb | dist | (legend row)
    mc,dc=H.mpl_cmaps(); nsb=len(sb); n,ncol=meth.shape
    wr=[1.4]+[0.05]*nsb+[1.0,0.16]+[0.05]*nsb+[1.0]
    fig=plt.figure(figsize=(11.5,5.2)); gs=fig.add_gridspec(2,len(wr),width_ratios=wr,height_ratios=[1,0.4],wspace=0.04,hspace=0.02)
    cl_seq=sb[0][1]; c=0
    axdn=fig.add_subplot(gs[0,c]); c+=1
    dendrogram(Z,orientation="left",color_threshold=tf,above_threshold_color="#999",ax=axdn,no_labels=True)
    axdn.set_xticks([]); axdn.set_yticks([]); axdn.set_title("UPGMA 樹",fontsize=8); [axdn.spines[s].set_visible(False) for s in axdn.spines]
    for lab,hexes in sb: H._sb(fig.add_subplot(gs[0,c]),hexes,lab); c+=1
    axm=fig.add_subplot(gs[0,c]); c+=1
    axm.imshow(meth,aspect="auto",cmap=mc,vmin=0,vmax=1,interpolation="nearest")
    sx=H.snv_fractional_x(cpgs,int(key.split('_')[1]))
    if sx is not None: axm.axvline(sx,color=H.SNV_COL,lw=2,zorder=4)
    prev=None
    for j in range(n):
        if prev is not None and cl_seq[j]!=prev: axm.axhline(j-0.5,color="black",lw=0.8)
        prev=cl_seq[j]
    axm.set_xticks([]); axm.set_yticks([]); axm.set_xlabel(f"{ncol} CpG·SNV標橙",fontsize=7); axm.set_title("甲基 read×CpG",fontsize=8.5)
    fig.add_subplot(gs[0,c]).axis("off"); c+=1
    for lab,hexes in sb: H._sb(fig.add_subplot(gs[0,c]),hexes,lab); c+=1
    axd=fig.add_subplot(gs[0,c])
    vmax=max(0.5,float(np.nanmax(dist)) if np.isfinite(np.nanmax(dist)) else 0.5)
    axd.imshow(dist,aspect="auto",cmap=dc,vmin=0,vmax=vmax,interpolation="nearest")
    prev=None
    for j in range(n):
        if prev is not None and cl_seq[j]!=prev: axd.axhline(j-0.5,color="white",lw=0.6,alpha=0.6)
        prev=cl_seq[j]
    axd.set_xticks([]); axd.set_yticks([]); axd.set_xlabel("read×read 距離·暗=近",fontsize=7); axd.set_title("距離(對角塊=群)",fontsize=8.5)
    H.grouped_legend(fig.add_subplot(gs[1,:]),[s[0] for s in sb],nclust)
    fig.suptitle(f"{key} [{r['group']}] n={n} tumor-only | coarse k={r['coarse_k']}({r['coarse_confidence']}) / "
                 f"fine k={r['fine_k']}({r['fine_confidence']}) · confirmed={r['confirmed_ks']} novel={r['novel_ks']}",fontsize=9.5,y=1.0)
    fn=f"{FIGS}/dash_{key}.png"; fig.savefig(fn,dpi=115,bbox_inches="tight"); plt.close(fig)
    print(f"  {key}: fine k={r['fine_k']} → dash_{key}.png")
    return f"figs_dashboard/dash_{key}.png"

KEYS=list(red.keys())
done=[plot(k) for k in KEYS if k in dirmap]
print(f"dashboard+dendro POC → {FIGS}")
