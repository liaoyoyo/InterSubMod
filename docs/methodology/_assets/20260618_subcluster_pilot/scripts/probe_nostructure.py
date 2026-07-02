#!/usr/bin/env python3
"""極端「真無結構」案例(decisionflow state② S2)— big7 本機跑 binary → 套新四閘 → 儀表板出圖供判斷。
驗證:方法是否正確把『多 read 但無結構』判成 NO_CLEAR(真無法切),而非誤切。"""
import os, csv, glob, json, sys, subprocess
import numpy as np
sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0,os.path.dirname(__file__)); import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; FIGS=f"{A}/figs_dashboard"
BIN=f"{WT}/build/bin/inter_sub_mod"
# big7 本機資料(非 big8 NFS)
TUMOR="/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NORMAL="/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF="/big7_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"
VCFDIR="/big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
os.environ["TMPDIR"]="/big7_disk/liaoyoyo2001/tmp"; os.makedirs("/big7_disk/liaoyoyo2001/tmp",exist_ok=True)
OUTD=f"{WT}/output/_nostructure_probe"; EDGE="#111111"; OUTL="#c9c9c9"
S2=[("chr1","23041213"),("chr1","196172653"),("chr1","100111922"),("chr1","62057441")]  # 高覆蓋 S2

# --- mini-VCF + binary (big7) ---
os.makedirs(OUTD,exist_ok=True)
poss={p for _,p in S2}
raw=subprocess.run(["zcat",f"{VCFDIR}/filtered_snv_tp_chr1.vcf.gz"],capture_output=True,text=True).stdout.splitlines()
hdr=[l for l in raw if l.startswith("#")]; body=[l for l in raw if not l.startswith("#") and l.split("\t")[1] in poss]
mini=f"{OUTD}/mini.vcf"; open(mini,"w").write("\n".join(hdr)+"\n"+"\n".join(body)+"\n")
subprocess.run(["bgzip","-f",mini]); subprocess.run(["tabix","-f","-p","vcf",mini+".gz"])
print(f"mini-VCF records={len(body)} → 跑 binary(big7)...",flush=True)
import shutil; shutil.rmtree(f"{OUTD}/run",ignore_errors=True); os.makedirs(f"{OUTD}/run")
subprocess.run([BIN,"-t",TUMOR,"-n",NORMAL,"-r",REF,"-v",mini+".gz","-w","5000","-j","16",
    "--distance-metric","BERNOULLI","--nan-distance-strategy","SKIP","-o",f"{OUTD}/run"],
    stdout=open(f"{OUTD}/run.log","w"),stderr=subprocess.STDOUT)
print("binary done",flush=True)

def grpcol(g):
    return EDGE if g=="edge" else OUTL if g=="outlier" else H.CLUSTER_COL[int(g)%len(H.CLUSTER_COL)]
dirmap={}
for mp in glob.glob(f"{OUTD}/run/**/distance/BERNOULLI/matrix.csv",recursive=True):
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd

rng=np.random.default_rng(20260622); done=[]
for chrom,pos in S2:
    key=f"{chrom}_{pos}"; rd=dirmap.get(key)
    if not rd: print(f"  {key}: 無輸出(可能 read 太少)"); continue
    reads={x["read_id"]:x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={x:i for i,x in enumerate(dids)}
    rows=open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); cpgs=[int(c) for c in rows[0].split(",")[1:]]
    mi={}; M=[]
    for j,ln in enumerate(rows[1:]):
        q=ln.split(","); mi[q[0]]=j; M.append([np.nan if v in("","NA","nan","NaN") else float(v) for v in q[1:]])
    M=np.array(M)
    it=lambda t:str(t) in("1","true","True")
    ids=[x for x in dids if x in reads and it(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    if len(ids)<CR.MIN_SZ*2: print(f"  {key}: tumor read 太少"); continue
    sub=D[np.ix_([di[x] for x in ids],[di[x] for x in ids])]; kp=CR.peel(sub)
    ids=[ids[i] for i in kp]; sub=sub[np.ix_(kp,kp)]
    P=np.array([M[mi[x]] for x in ids]); hp=[CR.LABMAP[reads[x]["hp"]] for x in ids]; al=[reads[x]["alt_support"] for x in ids]
    a=CR.analyze_locus(sub,P,hp,al,rng)
    fine={x:g for x,g in zip(ids,[p["grp"] for p in (a["perread_fine"] or [{"grp":"outlier"}]*len(ids))])}
    coarse={x:g for x,g in zip(ids,[p["grp"] for p in (a["perread_coarse"] or [{"grp":"outlier"}]*len(ids))])}
    Z,_=CR.linkZ(sub); dn=dendrogram(Z,orientation="left",no_plot=True); order=dn["leaves"][::-1]
    ids_o=[ids[i] for i in order]; meth=np.array([P[order[k2]] for k2 in range(len(order))])
    dist=sub[np.ix_(order,order)].copy(); np.fill_diagonal(dist,0); dist[dist<0]=np.nan
    sb=[("fine",[grpcol(fine.get(x,"outlier")) for x in ids_o]),("coarse",[grpcol(coarse.get(x,"outlier")) for x in ids_o])]
    sb+=H.sidebar_specs({x:reads[x] for x in ids_o},ids_o,cluster_of=None,include_tn=False,include_strand=True)
    mc,dc=H.mpl_cmaps(); nsb=len(sb); n,ncol=meth.shape
    wr=[1.4]+[0.05]*nsb+[1.0,0.16]+[0.05]*nsb+[1.0]
    fig=plt.figure(figsize=(11.5,5.2)); gs=fig.add_gridspec(2,len(wr),width_ratios=wr,height_ratios=[1,0.4],wspace=0.04,hspace=0.02)
    cl_seq=sb[0][1]; c=0; axdn=fig.add_subplot(gs[0,c]); c+=1
    dendrogram(Z,orientation="left",above_threshold_color="#999",ax=axdn,no_labels=True)
    axdn.set_xticks([]); axdn.set_yticks([]); axdn.set_title("UPGMA 樹",fontsize=8); [axdn.spines[s].set_visible(False) for s in axdn.spines]
    for lab,hexes in sb: H._sb(fig.add_subplot(gs[0,c]),hexes,lab); c+=1
    axm=fig.add_subplot(gs[0,c]); c+=1; axm.imshow(meth,aspect="auto",cmap=mc,vmin=0,vmax=1,interpolation="nearest")
    sx=H.snv_fractional_x(cpgs,int(pos));  axm.axvline(sx,color=H.SNV_COL,lw=2) if sx is not None else None
    axm.set_xticks([]); axm.set_yticks([]); axm.set_xlabel(f"{ncol} CpG·SNV標橙",fontsize=7); axm.set_title("甲基 read×CpG",fontsize=8.5)
    fig.add_subplot(gs[0,c]).axis("off"); c+=1
    for lab,hexes in sb: H._sb(fig.add_subplot(gs[0,c]),hexes,lab); c+=1
    axd=fig.add_subplot(gs[0,c]); vmax=max(0.5,float(np.nanmax(dist)) if np.isfinite(np.nanmax(dist)) else 0.5)
    axd.imshow(dist,aspect="auto",cmap=dc,vmin=0,vmax=vmax,interpolation="nearest")
    axd.set_xticks([]); axd.set_yticks([]); axd.set_xlabel("read×read 距離·暗=近",fontsize=7); axd.set_title("距離(對角無塊=無結構)",fontsize=8.5)
    H.grouped_legend(fig.add_subplot(gs[1,:]),[s[0] for s in sb],0)
    fig.suptitle(f"[S2 無訊號] {key} n={n} tumor-only → fine={a['fine_confidence']} (coarse k={a['coarse_k']}/{a['coarse_confidence']})",fontsize=9.5,y=1.0)
    fn=f"{FIGS}/nostruct_{key}.png"; fig.savefig(fn,dpi=115,bbox_inches="tight"); plt.close(fig)
    done.append({"key":key,"n":n,"verdict":a["fine_confidence"],"png":f"figs_dashboard/nostruct_{key}.png"})
    print(f"  {key} n={n} → verdict={a['fine_confidence']}",flush=True)
json.dump(done,open(f"{A}/nostructure_probe.json","w"),indent=1)
print("DONE",[d['verdict'] for d in done])
