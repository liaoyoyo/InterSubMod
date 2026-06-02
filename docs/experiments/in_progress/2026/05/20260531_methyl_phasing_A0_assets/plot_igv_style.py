#!/usr/bin/env python3
"""IGV-style read pileup：read 依基因組位置排列，HP 著色，CpG 甲基標點。
取代無法 headless 的 desktop IGV。數據從 BAM 真抽。"""
import argparse
import numpy as np
import pysam
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
matplotlib.rcParams["font.sans-serif"]=["Noto Sans CJK TC","Droid Sans Fallback","DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"]=False
BAM="/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"

def get_hp(r):
    try:return str(r.get_tag("HP"))
    except KeyError:return "."

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--region",required=True);ap.add_argument("--title",default="")
    ap.add_argument("--out",required=True);ap.add_argument("--win",type=int,default=2000)
    ap.add_argument("--max-reads",type=int,default=120)
    a=ap.parse_args()
    chrom,pos=a.region.split(":");pos=int(pos.replace(",",""))
    ws,we=pos-a.win,pos+a.win
    bam=pysam.AlignmentFile(BAM,"rb")
    reads=[]
    for r in bam.fetch(chrom,ws,we):
        if r.is_secondary or r.is_supplementary or r.is_unmapped:continue
        hp=get_hp(r)
        mods=r.modified_bases or {}
        mc=None
        for k,c in mods.items():
            if k[0] in ("C",b"C") and k[2] in ("m",27551):mc=c;break
        cpg=[]
        if mc:
            q2r={q:rr for q,rr in r.get_aligned_pairs(matches_only=True)}
            for qp,qual in mc:
                rp=q2r.get(qp)
                if rp and ws<=rp<=we:cpg.append((rp,qual/255.0))
        reads.append((r.reference_start,r.reference_end,hp,cpg))
    bam.close()
    # 依 HP 分組排序（HP1 上, unphase 中, HP2 下）
    grp={"1":0,"1-1":1,".":2,"3":2,"2-1":3,"2":4}
    reads.sort(key=lambda x:(grp.get(x[2],5),x[0]))
    reads=reads[:a.max_reads]
    Cmap={"1":"#2563EB","2":"#DC2626","1-1":"#60A5FA","2-1":"#F87171","3":"#A16207",".":"#9CA3AF"}
    fig,ax=plt.subplots(figsize=(13,8))
    for i,(s,e,hp,cpg) in enumerate(reads):
        y=len(reads)-1-i
        ax.add_patch(Rectangle((s,y+0.1),e-s,0.8,color=Cmap.get(hp,"#000"),alpha=0.35,lw=0))
        # 甲基點：紅=甲基化(>0.5) 藍=未(<0.5)
        for rp,pr in cpg:
            ax.plot(rp,y+0.5,marker="o",ms=2.2,color=("#7F1D1D" if pr>0.5 else "#1E3A8A"),alpha=0.85)
    ax.axvline(pos,color="black",ls="--",lw=1,alpha=0.6)
    ax.set_xlim(ws,we);ax.set_ylim(0,len(reads))
    ax.set_xlabel(f"{chrom} 位置 (bp)  — 虛線=中心 {pos:,}");ax.set_yticks([])
    ax.set_ylabel(f"read (n={len(reads)}, 依 HP 排序)")
    ax.set_title(f"{a.title}\n{a.region} | IGV-style: 橫條=read(HP著色), 點=CpG(紅甲基/藍未)",fontsize=10)
    leg=[Patch(facecolor="#2563EB",alpha=.4,label="HP1 read"),Patch(facecolor="#DC2626",alpha=.4,label="HP2 read"),
         Patch(facecolor="#9CA3AF",alpha=.4,label="unphase read"),
         plt.Line2D([],[],marker="o",ls="",color="#7F1D1D",label="CpG 甲基化"),
         plt.Line2D([],[],marker="o",ls="",color="#1E3A8A",label="CpG 未甲基")]
    ax.legend(handles=leg,bbox_to_anchor=(1.01,1),loc="upper left",fontsize=8)
    fig.savefig(a.out,dpi=120,bbox_inches="tight");plt.close()
    from collections import Counter
    print(f"  saved {a.out} reads={len(reads)} HP={dict(Counter(r[2] for r in reads))}")

if __name__=="__main__":main()
