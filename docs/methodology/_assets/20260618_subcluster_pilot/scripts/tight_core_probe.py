#!/usr/bin/env python3
"""[A pilot] bottom-up 小緊密核偵測 — 找 fine(top-down)併掉/丟掉的『小而緊密、清楚分離』的 read 核。
判準: 子樹群(3≤size≤15)且 分離比=群外均距/群內均距≥SEP_MIN 且 群內均距<整體中位(緊密)。
與 fine 比較: 核的 read 在 fine 是否被併進同一大群/標離群(=fine 漏掉)。渲染儀表板(+tight_core 側欄)。
重用 _ws_render 快取矩陣(不重跑 binary)。"""
import os, csv, glob, json, sys
import numpy as np
sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0,os.path.dirname(__file__)); import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD=f"{WT}/output/_ws_render"
FIGS=f"{A}/figs_tightcore"; os.makedirs(FIGS,exist_ok=True)
SEP_MIN=1.6; MINSZ=3; MAXSZ=15; CORECOL=['#e11d48','#0891b2','#65a30d','#7c3aed','#d97706','#0d9488','#be123c','#4338ca']

def _sep(s,S,rest):
    ss=s[np.ix_(S,rest)]; ss=ss[ss>=0]; sw=s[np.ix_(S,S)][np.triu_indices(len(S),1)]; sw=sw[sw>=0]
    if sw.size==0 or ss.size==0: return None,None
    wm=float(sw.mean()); return (float(ss.mean())/wm if wm>1e-6 else 99), wm
def find_tight_cores(sub,P,rng,RNULL=12):
    """子樹小群(3-15)+ 分離比≥SEP_MIN + 群內<中位 + 分離比顯著高於 within-1-group null(95pct)= 真小核。"""
    n=sub.shape[0]; Z,s=CR.linkZ(sub)
    med=np.median(s[np.triu_indices(n,1)][s[np.triu_indices(n,1)]>=0]) if n>1 else 0
    desc={i:[i] for i in range(n)}
    for i in range(len(Z)):
        a,b=int(Z[i,0]),int(Z[i,1]); desc[n+i]=desc[a]+desc[b]
    # 候選(緊密+分離)
    cand=[]
    for node,lv in desc.items():
        if not(MINSZ<=len(lv)<=MAXSZ) or len(lv)>=n-MINSZ: continue
        S=np.array(lv); mask=np.ones(n,bool); mask[S]=False; rest=np.where(mask)[0]
        if rest.size<MINSZ: continue
        sep,wm=_sep(s,S,rest)
        if sep is not None and sep>=SEP_MIN and wm<med: cand.append((set(lv),S,rest,sep,wm))
    if not cand: return []
    # within-1-group null: 打散 read 間結構(逐 CpG 欄內重排)→ 重算距離 → 同 S 的分離比分布
    C=P.shape[1]; null_dists=[]
    for _ in range(RNULL):
        Pn=P.copy()
        for cc in range(C):
            col=Pn[:,cc]; vi=np.where(~np.isnan(col))[0]
            if vi.size>1: Pn[vi,cc]=col[rng.permutation(vi)]
        Dn=CR.bernoulli_dist(Pn); np.fill_diagonal(Dn,0); Dn=np.maximum(Dn,Dn.T); null_dists.append(Dn)
    keep=[]
    for Sset,S,rest,sep,wm in cand:
        ns=[_sep(Dn,S,rest)[0] for Dn in null_dists]; ns=[x for x in ns if x is not None]
        null95=float(np.percentile(ns,95)) if ns else 0
        if sep>null95:   # 真實: 分離比顯著超過「無結構」基線
            keep.append({"S":Sset,"sep":round(sep,2),"win":round(wm,3),"sz":len(S),"null95":round(null95,2)})
    keep.sort(key=lambda c:(-c["sz"],-c["sep"]))
    out2=[]
    for c in keep:
        if not any(c["S"]<k["S"] for k in out2): out2.append(c)   # 取最大真核(非他核子集)
    return out2

dirmap={}
for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv",recursive=True):
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd
items=json.load(open(f"{A}/ws_items.json")); rng=np.random.default_rng(20260622); out=[]
def grpcol(g): return "#111111" if g=="edge" else "#c9c9c9" if g=="outlier" else H.CLUSTER_COL[int(g)%len(H.CLUSTER_COL)]
for it in items:
    key=it["key"]; rd=dirmap.get(key)
    if not rd: continue
    reads={x["read_id"]:x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={x:i for i,x in enumerate(dids)}
    rows=open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); cpgs=[int(c) for c in rows[0].split(",")[1:]]
    mi={};M=[]
    for j,ln in enumerate(rows[1:]):
        q=ln.split(","); mi[q[0]]=j; M.append([np.nan if v in("","NA","nan","NaN") else float(v) for v in q[1:]])
    M=np.array(M); itf=lambda t:str(t) in("1","true","True")
    ids=[x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    if len(ids)<MINSZ*2: continue
    sub=D[np.ix_([di[x] for x in ids],[di[x] for x in ids])]; kp=CR.peel(sub)
    ids=[ids[i] for i in kp]; sub=sub[np.ix_(kp,kp)]; n=len(ids)
    cores=find_tight_cores(sub)
    a=CR.analyze_locus(sub,np.array([M[mi[x]] for x in ids]),[CR.LABMAP[reads[x]["hp"]] for x in ids],[reads[x]["alt_support"] for x in ids],rng)
    fine=[p["grp"] for p in (a["perread_fine"] or [{"grp":"outlier"}]*n)]
    # 核 → read index 著色 + 判斷 fine 是否漏(核 reads 在 fine 同一群 或 全離群)
    core_lab=[-1]*n; missed=[]
    for ci,c in enumerate(cores):
        for i in c["S"]: core_lab[i]=ci
        fg={fine[i] for i in c["S"]}
        is_outl=all(isinstance(fine[i],str) for i in c["S"])  # 全離群/edge
        merged=(len([g for g in fg if isinstance(g,int)])<=1)  # fine 內 ≤1 群(被併)
        if is_outl or merged: missed.append({"core":ci,"sz":c["sz"],"sep":c["sep"],"reason":"離群" if is_outl else "併入大群"})
    # render
    Z,_=CR.linkZ(sub); dn=dendrogram(Z,orientation="left",no_plot=True); order=dn["leaves"][::-1]
    ids_o=[ids[i] for i in order]; meth=np.array([M[mi[x]] for x in ids_o]); dist=sub[np.ix_(order,order)].copy()
    np.fill_diagonal(dist,0); dist[dist<0]=np.nan
    cl_o=[core_lab[i] for i in order]; fine_o=[fine[i] for i in order]
    sb=[("tight核",["#eeeeee" if c<0 else CORECOL[c%len(CORECOL)] for c in cl_o]),
        ("fine",[grpcol(g) for g in fine_o])]
    sb+=H.sidebar_specs({x:reads[x] for x in ids_o},ids_o,cluster_of=None,include_tn=False,include_strand=True)
    mc,dc=H.mpl_cmaps(); nsb=len(sb); ncol=meth.shape[1]
    wr=[1.2]+[0.05]*nsb+[1.0,0.16]+[0.05]*nsb+[1.0]
    fig=plt.figure(figsize=(11,5.0)); gs=fig.add_gridspec(2,len(wr),width_ratios=wr,height_ratios=[1,0.4],wspace=0.04,hspace=0.02)
    c=0; axdn=fig.add_subplot(gs[0,c]); c+=1
    dendrogram(Z,orientation="left",above_threshold_color="#bbb",ax=axdn,no_labels=True)
    axdn.set_xticks([]); axdn.set_yticks([]); axdn.set_title("UPGMA 樹",fontsize=8); [axdn.spines[s2].set_visible(False) for s2 in axdn.spines]
    for lab,hx in sb: H._sb(fig.add_subplot(gs[0,c]),hx,lab); c+=1
    axm=fig.add_subplot(gs[0,c]); c+=1; axm.imshow(meth,aspect="auto",cmap=mc,vmin=0,vmax=1,interpolation="nearest")
    sx=H.snv_fractional_x(cpgs,int(key.split('_')[1]))
    if sx is not None: axm.axvline(sx,color=H.SNV_COL,lw=2)
    axm.set_xticks([]); axm.set_yticks([]); axm.set_title("甲基 read×CpG",fontsize=8.5)
    fig.add_subplot(gs[0,c]).axis("off"); c+=1
    for lab,hx in sb: H._sb(fig.add_subplot(gs[0,c]),hx,lab); c+=1
    axd=fig.add_subplot(gs[0,c]); vmax=max(0.5,float(np.nanmax(dist)) if np.isfinite(np.nanmax(dist)) else 0.5)
    axd.imshow(dist,aspect="auto",cmap=dc,vmin=0,vmax=vmax,interpolation="nearest")
    axd.set_xticks([]); axd.set_yticks([]); axd.set_title("距離(對角塊=群)",fontsize=8.5)
    H.grouped_legend(fig.add_subplot(gs[1,:]),[s2[0] for s2 in sb],0)
    fig.suptitle(f"{key} n={n} fine={a['fine_confidence']}(k{a['fine_k']}) | tight核={len(cores)} 個, fine 漏掉={len(missed)} 個",fontsize=9.5,y=1.0)
    fn=f"{FIGS}/tc_{key}.png"; fig.savefig(fn,dpi=108,bbox_inches="tight"); plt.close(fig)
    out.append({"key":key,"n":n,"fine_conf":a["fine_confidence"],"fine_k":a["fine_k"],
        "n_cores":len(cores),"cores":[{"sz":c["sz"],"sep":c["sep"],"win":c["win"]} for c in cores],
        "n_missed":len(missed),"missed":missed,"png":f"figs_tightcore/tc_{key}.png"})
    print(f"  {key} n={n} fine={a['fine_confidence']} tight核={len(cores)} fine漏={len(missed)} {[(m['sz'],m['reason']) for m in missed]}",flush=True)
json.dump(out,open(f"{A}/tight_core_probe.json","w"),indent=1)
tot_missed=sum(o["n_missed"] for o in out); withmiss=sum(1 for o in out if o["n_missed"]>0)
print(f"\nDONE {len(out)} 位點 · 偵測到 tight 核共 {sum(o['n_cores'] for o in out)} · fine 漏掉小核 {tot_missed}(分佈在 {withmiss} 位點)")
