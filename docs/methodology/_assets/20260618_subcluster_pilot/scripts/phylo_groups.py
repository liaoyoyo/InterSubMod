#!/usr/bin/env python3
"""[A pilot] 演化群階層標籤 — 遞迴遍歷 UPGMA 樹,每分裂測『兩子群是否都穩定(分離比顯著超 within-1-group null)』。
穩定→分,**較多 read 子群編號較小**(主群),階層標籤 1 / 1-1 / 1-2 / 2 / 2-1(1-1,1-2=『1』的相近子群)。
不穩定/太小→終止為該層標籤;非 cohesive 終端→離群。輸出最多穩定群 + 階層關係 + 各群對齊。渲染儀表板。
重用 _ws_render 快取矩陣(免 binary)。"""
import os, csv, glob, json, sys
import numpy as np
sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0,os.path.dirname(__file__)); import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD=f"{WT}/output/_ws_render"
FIGS=f"{A}/figs_phylo"; os.makedirs(FIGS,exist_ok=True)
MINSZ=3; SEP_MIN=1.3; RNULL=12
# 標籤色: 1-家族=暖色系深淺, 2-家族=冷色系, 3=綠系 ...
FAM={"1":["#7f1d1d","#dc2626","#f87171","#fca5a5"],"2":["#1e3a8a","#2563eb","#60a5fa","#93c5fd"],
     "3":["#14532d","#16a34a","#4ade80"],"4":["#713f12","#d97706","#fbbf24"],"5":["#581c87","#9333ea","#c084fc"]}
def labcol(lab):
    if lab is None or lab in("outlier","-"): return "#cfcfcf"
    fam=lab.split("-")[0]; depth=lab.count("-"); pal=FAM.get(fam,["#555","#888","#aaa"]); return pal[min(depth,len(pal)-1)]

def _bw(s,S1,S2):
    """between(C1,C2) / within(C1)+within(C2) 平均比 — 兄弟子群分離度。"""
    bet=s[np.ix_(S1,S2)]; bet=bet[bet>=0]
    w1=s[np.ix_(S1,S1)][np.triu_indices(len(S1),1)]; w1=w1[w1>=0]
    w2=s[np.ix_(S2,S2)][np.triu_indices(len(S2),1)]; w2=w2[w2>=0]
    wm=np.concatenate([w1,w2]) if (w1.size or w2.size) else np.array([])
    if bet.size==0 or wm.size==0 or wm.mean()<=1e-6: return None
    return float(bet.mean())/float(wm.mean())

def phylo_label(sub,P,rng):
    n=sub.shape[0]; Z,s=CR.linkZ(sub)
    desc={i:[i] for i in range(n)}; ch={}
    for i in range(len(Z)):
        a,b=int(Z[i,0]),int(Z[i,1]); desc[n+i]=desc[a]+desc[b]; ch[n+i]=(a,b)
    C=P.shape[1]; nulls=[]
    for _ in range(RNULL):
        Pn=P.copy()
        for cc in range(C):
            col=Pn[:,cc]; vi=np.where(~np.isnan(col))[0]
            if vi.size>1: Pn[vi,cc]=col[rng.permutation(vi)]
        Dn=CR.bernoulli_dist(Pn); np.fill_diagonal(Dn,0); nulls.append(np.maximum(Dn,Dn.T))
    def split_real(c1,c2):
        """此分裂(兄弟 C1 vs C2)是否真實: between/within ≥SEP_MIN 且 顯著超 within-1-group null。"""
        S1,S2=np.array(desc[c1]),np.array(desc[c2])
        if min(len(S1),len(S2))<MINSZ: return False
        r=_bw(s,S1,S2)
        if r is None or r<SEP_MIN: return False
        ns=[_bw(Dn,S1,S2) for Dn in nulls]; ns=[x for x in ns if x is not None]
        return r>(np.percentile(ns,95) if ns else 0)
    lab=[None]*n
    def rec(node,label):
        leaves=desc[node]
        if node not in ch or len(leaves)<2*MINSZ:
            for i in leaves: lab[i]=label; return
        c1,c2=ch[node]
        if split_real(c1,c2):
            big,small=(c1,c2) if len(desc[c1])>=len(desc[c2]) else (c2,c1)
            rec(big,label+"-1"); rec(small,label+"-2")
        else:
            for i in leaves: lab[i]=label   # 不再分裂 = 此層終端群
    root=2*n-2; c1,c2=ch[root]
    if split_real(c1,c2):
        big,small=(c1,c2) if len(desc[c1])>=len(desc[c2]) else (c2,c1); rec(big,"1"); rec(small,"2")
    else:
        for i in range(n): lab[i]="1"   # 無真實分裂 = 單一未分化群(≈NO_CLEAR)
    return lab,Z

dirmap={}
for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv",recursive=True):
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd
def grpcol(g): return "#111111" if g=="edge" else "#c9c9c9" if g=="outlier" else H.CLUSTER_COL[int(g)%len(H.CLUSTER_COL)]
items=json.load(open(f"{A}/ws_items.json")); rng=np.random.default_rng(20260622); out=[]
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
    ids=[ids[i] for i in kp]; sub=sub[np.ix_(kp,kp)]; n=len(ids); P=np.array([M[mi[x]] for x in ids])
    lab,Z=phylo_label(sub,P,rng); lab=[l if l else "outlier" for l in lab]
    from collections import Counter as _C
    _sm={L for L,c in _C(l for l in lab if l!="outlier").items() if c<MINSZ}
    lab=[("outlier" if l in _sm else l) for l in lab]
    from collections import Counter
    labct=Counter(l for l in lab if l!="outlier"); nterm=len(labct); nout=sum(1 for l in lab if l=="outlier")
    # 各群對齊(dominant germline 軸)
    def dom(idxs,field):
        c=Counter(reads[ids[i]][field] for i in idxs); return c.most_common(1)[0] if c else ("-",0)
    galign={}
    for L in labct:
        idxs=[i for i in range(n) if lab[i]==L]
        hpd=dom(idxs,"hp"); ald=dom(idxs,"alt_support")
        galign[L]={"n":len(idxs),"hp":f"{hpd[0]}({100*hpd[1]//len(idxs)}%)","allele":f"{ald[0]}({100*ald[1]//len(idxs)}%)"}
    # render: 樹依終端標籤著色 + 標籤側欄
    dn=dendrogram(Z,orientation="left",no_plot=True); order=dn["leaves"][::-1]
    ids_o=[ids[i] for i in order]; meth=np.array([P[i] for i in order]); dist=sub[np.ix_(order,order)].copy()
    np.fill_diagonal(dist,0); dist[dist<0]=np.nan; lab_o=[lab[i] for i in order]
    # tree link color by terminal label
    co={i:lab[i] for i in range(n)}
    desc={i:[i] for i in range(n)}; ch={}
    for i in range(len(Z)):
        a,b=int(Z[i,0]),int(Z[i,1]); desc[n+i]=desc[a]+desc[b]; ch[n+i]=(a,b)
    nodelab={}
    for nd in range(2*n-1):
        ll={lab[i] for i in desc[nd]}; nodelab[nd]=(list(ll)[0] if len(ll)==1 else None)
    def lcf(k): return labcol(nodelab.get(k)) if nodelab.get(k) else "#cccccc"
    sb=[("演化群",[labcol(l) for l in lab_o])]
    sb+=H.sidebar_specs({x:reads[x] for x in ids_o},ids_o,cluster_of=None,include_tn=False,include_strand=True)
    mc,dc=H.mpl_cmaps(); nsb=len(sb); ncol=meth.shape[1]
    wr=[1.3]+[0.06]*nsb+[1.0,0.16]+[0.06]*nsb+[1.0]
    fig=plt.figure(figsize=(11,5.2)); gs=fig.add_gridspec(2,len(wr),width_ratios=wr,height_ratios=[1,0.42],wspace=0.04,hspace=0.02)
    c=0; axdn=fig.add_subplot(gs[0,c]); c+=1
    dendrogram(Z,orientation="left",link_color_func=lcf,ax=axdn,no_labels=True)
    axdn.set_xticks([]); axdn.set_yticks([]); axdn.set_title("UPGMA 樹(分支=演化群)",fontsize=8); [axdn.spines[s2].set_visible(False) for s2 in axdn.spines]
    for lb,hx in sb: H._sb(fig.add_subplot(gs[0,c]),hx,lb); c+=1
    axm=fig.add_subplot(gs[0,c]); c+=1; axm.imshow(meth,aspect="auto",cmap=mc,vmin=0,vmax=1,interpolation="nearest")
    sx=H.snv_fractional_x(cpgs,int(key.split('_')[1]))
    if sx is not None: axm.axvline(sx,color=H.SNV_COL,lw=2)
    axm.set_xticks([]); axm.set_yticks([]); axm.set_title("甲基 read×CpG",fontsize=8.5)
    fig.add_subplot(gs[0,c]).axis("off"); c+=1
    for lb,hx in sb: H._sb(fig.add_subplot(gs[0,c]),hx,lb); c+=1
    axd=fig.add_subplot(gs[0,c]); vmax=max(0.5,float(np.nanmax(dist)) if np.isfinite(np.nanmax(dist)) else 0.5)
    axd.imshow(dist,aspect="auto",cmap=dc,vmin=0,vmax=vmax,interpolation="nearest")
    axd.set_xticks([]); axd.set_yticks([]); axd.set_title("距離(對角塊=群)",fontsize=8.5)
    # 群+對齊 圖例
    axl=fig.add_subplot(gs[1,:]); axl.axis("off")
    txt=" · ".join(f"{L}:n{galign[L]['n']} hp={galign[L]['hp']} alle={galign[L]['allele']}" for L in sorted(labct))
    axl.text(0.0,0.7,f"穩定群 {nterm} 個 + 離群 {nout}:",fontsize=8.5,fontweight="bold",transform=axl.transAxes)
    axl.text(0.0,0.3,txt,fontsize=7.5,transform=axl.transAxes,wrap=True)
    fig.suptitle(f"{key} n={n} fine={it['fine_conf']} → 演化群 {nterm} 個(標籤 {sorted(labct)})離群 {nout}",fontsize=9.5,y=1.0)
    fn=f"{FIGS}/phy_{key}.png"; fig.savefig(fn,dpi=108,bbox_inches="tight"); plt.close(fig)
    out.append({"key":key,"n":n,"fine_conf":it["fine_conf"],"n_groups":nterm,"n_outlier":nout,
        "labels":dict(labct),"align":galign,"png":f"figs_phylo/phy_{key}.png"})
    print(f"  {key} n={n} fine={it['fine_conf']} → 群{nterm} {sorted(labct)} 離群{nout}",flush=True)
json.dump(out,open(f"{A}/phylo_groups.json","w"),indent=1)
print(f"\nDONE {len(out)} 位點 · 平均穩定群 {np.mean([o['n_groups'] for o in out]):.1f}")
