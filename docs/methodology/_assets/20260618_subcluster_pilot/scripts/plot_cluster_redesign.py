#!/usr/bin/env python3
"""old vs new 切群對照圖（eyeball 驗證）。距離熱圖依 NEW primary 重排，左側欄:
OLD maxclust | NEW primary(核心群+離群) | hp | carrier | allele。多解析度位點加一條 confirmed 次解析度。
純讀 cluster_redesign.json + 重載既有距離矩陣（不重跑 binary）。"""
import os, csv, glob, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams["font.family"]=["Droid Sans Fallback","DejaVu Sans"]; plt.rcParams["axes.unicode_minus"]=False

WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUT=f"{WT}/output/_kprofile_heatmap"
FIGS=f"{A}/figs_redesign"; os.makedirs(FIGS,exist_ok=True)
CLUSTER_COL=['#db2777','#0d9488','#0891b2','#65a30d','#7c3aed','#d97706','#475569','#be123c','#0369a1','#4d7c0f']
OUTCOL="#111111"
def _f(x):
    x=str(x).strip(); return np.nan if x in("","NA","nan","NaN") else float(x)
def loadm(p):
    r=open(p).read().strip().split("\n"); ids=[];M=[]
    for ln in r[1:]:
        q=ln.split(","); ids.append(q[0]); M.append([_f(x) for x in q[1:]])
    return ids,np.array(M)
dirmap={}
for mp in glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True):
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd

D=json.load(open(f"{A}/cluster_redesign.json"))
loci={r["key"]:r for r in D["loci"]}
# join stability-excess（real? 證據）
try:
    ST={(r["chrom"]+"_"+r["pos"]):{p["k"]:p["stab_excess"] for p in r["per_k"]}
        for r in json.load(open(f"{A}/kprofile_stability.json"))["loci"]}
except Exception: ST={}

def strip(ax, colidx, palette, x0, label):
    """畫一條 n×1 色帶。colidx: list of int(-1=outlier→黑)。"""
    n=len(colidx); img=np.zeros((n,1,3))
    for i,c in enumerate(colidx):
        hexc=OUTCOL if c<0 else palette[c%len(palette)]
        img[i,0]=[int(hexc[j:j+2],16)/255 for j in (1,3,5)]
    ax.imshow(img,aspect="auto",extent=[x0,x0+0.85,n,0])
    ax.text(x0+0.42,-n*0.012,label,ha="left",va="bottom",fontsize=7.5,rotation=90)

def plot_locus(key, extra_res_k=None):
    r=loci[key]; rd=dirmap[key]
    dids,Dm=loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={x:i for i,x in enumerate(dids)}
    ids=r["read_ids"]; sub=Dm[np.ix_([di[x] for x in ids],[di[x] for x in ids])].copy()
    np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T); sub[sub<0]=np.nan
    n=len(ids)
    # primary labels
    prim=None
    for res in r["new_resolutions"]:
        if res["k_core"]==r["primary_k"]: prim=res
    if prim is None: prim=r["new_resolutions"][-1] if r["new_resolutions"] else None
    new_lab=prim["labels"] if prim else [1]*n; new_out=set(prim["out_idx"]) if prim else set()
    old_lab=r["old"]["labels"] if r["old"]["labels"] else [1]*n
    hp=r["hp_lab"]; carrier=["G" if l in("1","2") else "C" for l in hp]; allele=r["alle_lab"]
    hpc=[0 if l in("1","1-1") else 1 for l in hp]; carc=[0 if c=="G" else 1 for c in carrier]
    allc=[0 if a=="REF" else (1 if a=="ALT" else -1) for a in allele]
    # remap cluster labels to 0..k; outliers -1
    def remap(lab,out=set()):
        u=sorted(set(lab[i] for i in range(len(lab)) if i not in out)); m={v:i for i,v in enumerate(u)}
        return [(-1 if i in out else m[lab[i]]) for i in range(len(lab))]
    new_c=remap(new_lab,new_out); old_c=remap(old_lab)
    # extra resolution
    extra_c=None; extra_lbl=""
    if extra_res_k:
        for res in r["new_resolutions"]:
            if res["k_core"]==extra_res_k:
                extra_c=remap(res["labels"],set(res["out_idx"]))
                ax_a=[a for a in res["aligned_axes"]]; extra_lbl=f"k{extra_res_k}\n{'/'.join(ax_a) or '?'}"
    # order: primary core group, then carrier, then allele; outliers last
    order=sorted(range(n),key=lambda i:(99 if i in new_out else new_c[i], carc[i], allc[i]))
    sub=sub[np.ix_(order,order)]
    reidx=lambda L:[L[i] for i in order]
    new_c,old_c,hpc,carc,allc=map(reidx,[new_c,old_c,hpc,carc,allc])
    if extra_c is not None: extra_c=reidx(extra_c)
    # figure: sidebars + heatmap
    bars=[("OLD maxclust",old_c,CLUSTER_COL),("NEW primary",new_c,CLUSTER_COL)]
    if extra_c is not None: bars.append((extra_lbl,extra_c,CLUSTER_COL))
    bars+=[("hp",hpc,['#2563eb','#f59e0b']),("carrier",carc,['#16a34a','#dc2626']),("allele",allc,['#0891b2','#db2777'])]
    nb=len(bars)
    fig=plt.figure(figsize=(9.2,7.2))
    gs=fig.add_gridspec(1,2,width_ratios=[nb*0.95,10],wspace=0.04)
    axs=fig.add_subplot(gs[0]); axh=fig.add_subplot(gs[1])
    for j,(lbl,col,pal) in enumerate(bars): strip(axs,col,pal,j,lbl)
    axs.set_xlim(0,nb); axs.set_ylim(n,0); axs.axis("off")
    im=axh.imshow(sub,aspect="auto",cmap="magma_r",vmin=0,vmax=np.nanmax(sub)); axh.set_xticks([]); axh.set_yticks([])
    axh.set_title("read×read 距離(依 NEW primary 重排;亮黃=近=同群)",fontsize=9)
    cks=r["confirmed_ks"]; nks=r.get("novel_ks",[])
    exc=prim.get("stab_excess") if prim else None
    excs=f"  excess(vs null)={exc:+.2f}{'(真實)' if (exc or 0)>=0.10 else '(弱)'}" if exc is not None else ""
    edge="  +邊緣群" if (prim and prim.get("edge_group")) else ""
    fig.suptitle(f"{key} [{r['group']}] n={n}   OLD maxclust k={r['old']['k_groups']}  →  "
                 f"NEW primary k={r['primary_k']}({r['primary_confidence']}){excs}{edge}\n"
                 f"confirmed_k={cks}  novel_k(subclone候選)={nks}  離群(黑)={len(new_out)}   "
                 f"側欄: hp 藍/橙 · carrier 綠/紅 · allele 青/桃",fontsize=9.5,y=1.0)
    fig.colorbar(im,ax=axh,fraction=0.025,pad=0.01)
    fn=f"{FIGS}/redesign_{key}.png"; fig.tight_layout(rect=[0,0,1,0.94]); fig.savefig(fn,dpi=110,bbox_inches="tight"); plt.close(fig)
    return f"figs_redesign/redesign_{key}.png"

# 代表位點（三閘收斂後跨 confidence 光譜）
sel=[("chr8_132099309",2),    # 多解析度 CONFIRMED: k=5(carrier)+k=2(allele)
     ("chr4_190112507",None), # null 修掉過切→回 k=2 CONFIRMED(高覆蓋乾淨)
     ("chr8_138619384",3),    # 混合: k=2 CONFIRMED(germline) + k=3-6 NOVEL(subclone候選)
     ("chr8_64387206",None),  # REAL_NOVEL k=4(confident-unique 其實藏 subclone 候選)
     ("chr17_79567788",None), # REAL_NOVEL k=2(不對齊 germline)
     ("chr1_218115208",None), # balance+null 後 CONFIRMED k=2(太鬆已修)
     ("chr10_55938697",None), # CONFIRMED 但鬆散
     ("chr2_44798355",None)]  # NO_CLEAR_SPLIT
done=[]
for key,extra in sel:
    if key in loci and key in dirmap:
        png=plot_locus(key,extra); r=loci[key]
        done.append({"key":key,"png":png,"primary_k":r["primary_k"],"confidence":r["primary_confidence"],
                     "confirmed_ks":r["confirmed_ks"],"old_k":r["old"]["k_groups"]})
        print(f"  {key}: OLD k={r['old']['k_groups']} → primary k={r['primary_k']}({r['primary_confidence']}) confirmed={r['confirmed_ks']}")
json.dump({"figs":done},open(f"{A}/cluster_redesign_figindex.json","w"),indent=1)
print(f"plotted {len(done)} → {FIGS}")
