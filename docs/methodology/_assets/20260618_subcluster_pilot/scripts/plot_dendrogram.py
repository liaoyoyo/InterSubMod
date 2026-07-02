#!/usr/bin/env python3
"""UPGMA 距離樹圖(dendrogram)+ 標色分支 + a-priori 葉標 + 距離熱圖。
讓使用者看樹結構判斷分類(如 chr4 大群底下對齊 ALT 的子群是否該切)。
標色: 分支依 fine-k cut 著色;另畫 k=2/3/4 cut 高度虛線。葉旁色帶 = allele/hp/carrier。"""
import os, csv, glob, json, sys
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]=["Droid Sans Fallback","DejaVu Sans"]; plt.rcParams["axes.unicode_minus"]=False

WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUT=f"{WT}/output/_kprofile_heatmap"
FIGS=f"{A}/figs_dendro"; os.makedirs(FIGS,exist_ok=True)
CL=['#db2777','#0d9488','#0891b2','#65a30d','#7c3aed','#d97706','#be123c','#0369a1']
dirmap={}
for mp in glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True):
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd
red={r["key"]:r for r in json.load(open(f"{A}/cluster_redesign.json"))["loci"]}

def load(key):
    rd=dirmap[key]
    reads={r["read_id"]:r for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={r:i for i,r in enumerate(dids)}
    it=lambda t:str(t) in("1","true","True")
    kids=[r for r in dids if r in reads and it(reads[r]["is_tumor"]) and reads[r]["hp"] in CR.LABMAP]
    sub=D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]; kp=CR.peel(sub)
    kids=[kids[i] for i in kp]; sub=sub[np.ix_(kp,kp)]
    hp=[CR.LABMAP[reads[r]["hp"]] for r in kids]; alle=[reads[r]["alt_support"] for r in kids]
    return sub,hp,alle

def cut_height(Z,k):
    h=np.sort(Z[:,2]); return (h[-(k-1)]+h[-k])/2 if k>=2 and k<=len(h) else None

def strip(ax,x0,colidx,pal):
    n=len(colidx); img=np.zeros((n,1,3))
    for i,c in enumerate(colidx):
        hexc="#bbbbbb" if c<0 else pal[c%len(pal)]; img[i,0]=[int(hexc[j:j+2],16)/255 for j in (1,3,5)]
    ax.imshow(img,aspect="auto",extent=[x0,x0+0.9,n,0])

def plot(key):
    sub,hp,alle=load(key); n=sub.shape[0]; Z,s=CR.linkZ(sub); r=red.get(key,{})
    fk=r.get("fine_k",2) or 2; tf=cut_height(Z,fk)
    fig=plt.figure(figsize=(12,max(5,n*0.06)))
    gs=fig.add_gridspec(1,5,width_ratios=[4,0.25,0.25,0.25,5],wspace=0.03)
    axd=fig.add_subplot(gs[0])
    dn=dendrogram(Z,orientation="left",color_threshold=tf,above_threshold_color="#999",
                  link_color_func=None,ax=axd,no_labels=True)
    order=dn["leaves"][::-1]   # 上→下
    # k=2/3/4 cut 高度虛線
    for k in (2,3,4):
        t=cut_height(Z,k)
        if t: axd.axvline(t,ls="--",c="#d97706",lw=0.8); axd.text(t,n*1.01,f"k{k}",fontsize=7,color="#d97706",ha="center")
    axd.set_title(f"UPGMA 樹(分支色=fine k={fk} cut)",fontsize=9); axd.set_xlabel("距離"); axd.set_yticks([])
    hp_c=[0 if l in("1","1-1") else 1 for l in hp]; carc=[0 if l in("1","2") else 1 for l in hp]
    allc=[0 if a=="REF" else (1 if a=="ALT" else -1) for a in alle]
    reb=lambda L:[L[i] for i in order]
    for j,(lbl,col,pal) in enumerate([("allele",reb(allc),['#0891b2','#db2777']),
                                       ("hp",reb(hp_c),['#2563eb','#f59e0b']),
                                       ("carrier",reb(carc),['#16a34a','#dc2626'])]):
        ax=fig.add_subplot(gs[1+j]); strip(ax,0,col,pal); ax.set_xlim(0,1); ax.set_ylim(n,0); ax.axis("off")
        ax.set_title(lbl,fontsize=7.5)
    axh=fig.add_subplot(gs[4]); subo=sub[np.ix_(order,order)]; np.fill_diagonal(subo,0); subo[subo<0]=np.nan
    im=axh.imshow(subo,aspect="auto",cmap="magma_r",vmin=0,vmax=np.nanmax(subo)); axh.set_xticks([]); axh.set_yticks([])
    axh.set_title("距離熱圖(依樹葉序;亮黃=近)",fontsize=9)
    fig.colorbar(im,ax=axh,fraction=0.025,pad=0.01)
    fig.suptitle(f"{key} [{r.get('group','')}] n={n}  coarse k={r.get('coarse_k')}({r.get('coarse_confidence')}) | "
                 f"fine k={r.get('fine_k')}({r.get('fine_confidence')})  confirmed={r.get('confirmed_ks')} novel={r.get('novel_ks')}",
                 fontsize=10.5,y=1.0)
    fn=f"{FIGS}/dendro_{key}.png"; fig.tight_layout(rect=[0,0,1,0.96]); fig.savefig(fn,dpi=110,bbox_inches="tight"); plt.close(fig)
    # gap 階梯: 各 k 的合併高度跳躍 = h[-(k-1)]-h[-k]，相對中位 gap 的倍數
    h=np.sort(Z[:,2]); med=np.median(np.diff(h[h>0])) if (np.diff(h)>0).any() else 1e-9
    gapr={k:round((h[-(k-1)]-h[-k])/med,1) for k in range(2,7) if k<=len(h)}
    print(f"  {key:16s} n={n} fine_k={r.get('fine_k')} gap_ratio(×median): "
          +" ".join(f"k{k}={g}" for k,g in gapr.items()))
    GAPS[key]={"n":n,"group":r.get("group"),"coarse_k":r.get("coarse_k"),"coarse_confidence":r.get("coarse_confidence"),
               "fine_k":r.get("fine_k"),"fine_confidence":r.get("fine_confidence"),
               "confirmed_ks":r.get("confirmed_ks"),"novel_ks":r.get("novel_ks"),"gap_ratio":gapr,
               "png":f"figs_dendro/dendro_{key}.png"}
    return f"figs_dendro/dendro_{key}.png"

KEYS=["chr4_190112507",                                   # LOH ALT 子群
      "chr8_132099309","chr14_97275848","chr2_90072616",  # germline-cis CONFIRMED
      "chr8_138619384","chr15_94640645","chr8_64387206",  # 骨幹+novel
      "chr1_218115208","chr17_79567788",                  # 純novel
      "chr1_119446237","chr3_1026519","chr3_131113525","chr2_44798355"]  # 4 個 REAL_DIFFUSE(原誤判 NO_CLEAR)
GAPS={}
done=[plot(k) for k in KEYS if k in dirmap]
json.dump({"figs":done,"gaps":GAPS},open(f"{A}/dendro_figindex.json","w"),indent=1)
print(f"dendrograms → {FIGS}")
