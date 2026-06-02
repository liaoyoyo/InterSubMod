#!/usr/bin/env python3
"""精選 loci 展示熱圖：per-read 甲基 + 救援前後 HP 側欄（可驗證）。
數據從 BAM 真實抽，不產生統計數字（只畫圖）。"""
import sys, argparse
import numpy as np
import pysam
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
matplotlib.rcParams["font.sans-serif"]=["Noto Sans CJK TC","Droid Sans Fallback","DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"]=False

BAM="/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
AD="/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets"

def get_hp(r):
    try: return str(r.get_tag("HP"))
    except KeyError: return "."

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--region",required=True); ap.add_argument("--title",default="")
    ap.add_argument("--out",required=True); ap.add_argument("--win",type=int,default=2000)
    a=ap.parse_args()
    chrom,pos=a.region.split(":"); pos=int(pos.replace(",",""))
    ws,we=pos-a.win,pos+a.win
    bam=pysam.AlignmentFile(BAM,"rb")
    rows=[];seen=set();allpos=set()
    for read in bam.fetch(chrom,ws,we):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:continue
        if read.query_name in seen:continue
        hp=get_hp(read)
        mods=read.modified_bases or {}
        mc=None
        for k,c in mods.items():
            if k[0] in ("C",b"C") and k[2] in ("m",27551):mc=c;break
        if not mc:continue
        q2r={q:r for q,r in read.get_aligned_pairs(matches_only=True)}
        meth={}
        for qp,qual in mc:
            rp=q2r.get(qp)
            if rp is None or rp<ws or rp>we:continue
            meth[rp]=qual/255.0
        if len(meth)<5:continue
        seen.add(read.query_name); allpos.update(meth.keys())
        rows.append((hp,meth))
    bam.close()
    if len(rows)<10:
        print(f"  too few reads {a.region}");return
    positions=sorted(allpos)
    cov={p:sum(1 for _,m in rows if p in m) for p in positions}
    positions=[p for p in positions if cov[p]>=0.3*len(rows)]
    # 救援：unphase read 用 anchor 質心預測
    def vec(m):return np.array([m.get(p,np.nan) for p in positions])
    anc1=[m for h,m in rows if h=="1"];anc2=[m for h,m in rows if h=="2"]
    pred={}
    if len(anc1)>=4 and len(anc2)>=4:
        c1=np.nanmean([vec(m) for m in anc1],axis=0);c2=np.nanmean([vec(m) for m in anc2],axis=0)
        for i,(h,m) in enumerate(rows):
            if h==".":
                v=vec(m);mask=~np.isnan(v)&~np.isnan(c1)&~np.isnan(c2)
                if mask.sum()>=3:
                    d1=np.sqrt(np.mean((v[mask]-c1[mask])**2));d2=np.sqrt(np.mean((v[mask]-c2[mask])**2))
                    pred[i]="1pred" if d1<d2 else "2pred"
    # 排序：HP1, HP1pred, unphase, HP2pred, HP2
    def keyf(i):
        h=rows[i][0]
        if h=="1":return 0
        if pred.get(i)=="1pred":return 1
        if h==".":return 2
        if pred.get(i)=="2pred":return 3
        if h=="2":return 4
        return 5
    order=sorted(range(len(rows)),key=keyf)
    M=np.full((len(rows),len(positions)),np.nan)
    for r_i,i in enumerate(order):
        for j,p in enumerate(positions):
            if p in rows[i][1]:M[r_i,j]=rows[i][1][p]
    # 側欄顏色：原 HP + 救援
    Cmap={"1":"#2563EB","2":"#DC2626","1-1":"#60A5FA","2-1":"#F87171","3":"#A16207",".":"#9CA3AF"}
    orig_col=[Cmap.get(rows[i][0],"#000") for i in order]
    resc_col=[("#93C5FD" if pred.get(i)=="1pred" else "#FCA5A5" if pred.get(i)=="2pred" else Cmap.get(rows[i][0],"#9CA3AF")) for i in order]
    fig,(axo,axr,axh)=plt.subplots(1,3,figsize=(13,8),gridspec_kw={"width_ratios":[1,1,40],"wspace":0.02})
    for ax,cols,t in [(axo,orig_col,"原HP"),(axr,resc_col,"救援後")]:
        ax.imshow([[0]],aspect="auto",extent=[0,1,0,len(order)],alpha=0)
        for k,c in enumerate(cols):ax.add_patch(plt.Rectangle((0,len(order)-1-k),1,1,color=c))
        ax.set_xlim(0,1);ax.set_ylim(0,len(order));ax.set_xticks([]);ax.set_yticks([]);ax.set_title(t,fontsize=9)
    Mfill=np.where(np.isnan(M),np.nanmean(M),M)
    axh.imshow(Mfill,aspect="auto",cmap="RdBu_r",vmin=0,vmax=1,interpolation="nearest")
    axh.set_xlabel(f"CpG 位點 (n={len(positions)})");axh.set_yticks([])
    axh.set_title(f"{a.title}\n{a.region} | read n={len(rows)} | 救援 {len(pred)} unphase read",fontsize=10)
    leg=[Patch(facecolor="#2563EB",label="HP1"),Patch(facecolor="#DC2626",label="HP2"),
         Patch(facecolor="#9CA3AF",label="unphase"),Patch(facecolor="#93C5FD",label="→救HP1"),Patch(facecolor="#FCA5A5",label="→救HP2")]
    axh.legend(handles=leg,bbox_to_anchor=(1.01,1),loc="upper left",fontsize=8)
    fig.savefig(a.out,dpi=120,bbox_inches="tight");plt.close()
    print(f"  saved {a.out} (reads={len(rows)}, rescued={len(pred)})")

if __name__=="__main__":main()
