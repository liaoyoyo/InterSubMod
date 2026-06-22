#!/usr/bin/env python3
"""觀察工作站資料層: 每類抽 N 代表(chr20-22, big7 binary)→ analyze_locus + 渲染儀表板圖 + 擷取卡片資料。
輸出 ws_items.json(每位點: 分類/excess/對齊/群大小/per-read + png)供 build_observation_workstation.py。"""
import os, csv, glob, json, sys, subprocess, shutil
os.environ.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
import numpy as np
sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0,os.path.dirname(__file__)); import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; FIGS=f"{A}/figs_ws"; os.makedirs(FIGS,exist_ok=True)
BIN=f"{WT}/build/bin/inter_sub_mod"
TUMOR="/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"; NORMAL="/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF="/big7_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"; VCFDIR="/big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
os.environ["TMPDIR"]="/big7_disk/liaoyoyo2001/tmp"; OUTD=f"{WT}/output/_ws_render"; EDGE="#111111"; OUTL="#c9c9c9"
NPER=6; CLASSES=["CONFIRMED","NEAR_CONFIRMED","REAL_NOVEL","REAL_DIFFUSE","NO_CLEAR_SPLIT"]

rec=[r for r in json.load(open(f"{A}/cluster_redesign_wg_records.json")) if "err" not in r and r["set"]=="TP"]
# 抽樣: chr20-22 優先, 每類取 n 分散(小/中/大)
sel={}
for cl in CLASSES:
    pool=[r for r in rec if r["fine_conf"]==cl and r["chrom"] in("chr20","chr21","chr22")]
    if len(pool)<NPER: pool=[r for r in rec if r["fine_conf"]==cl]  # NEAR 等少數 → 全基因組
    pool=sorted(pool,key=lambda r:r["n"])
    idx=np.linspace(0,len(pool)-1,min(NPER,len(pool))).astype(int) if pool else []
    sel[cl]=[pool[i] for i in idx]
picks=[r for cl in CLASSES for r in sel[cl]]
chrs=sorted({r["chrom"] for r in picks})
print("抽樣:",{cl:len(sel[cl]) for cl in CLASSES},"chr:",chrs,flush=True)

# mini-VCF + binary(big7), per chr
shutil.rmtree(OUTD,ignore_errors=True); os.makedirs(OUTD)
dirmap={}
for chrom in chrs:
    poss={r["pos"] for r in picks if r["chrom"]==chrom}
    raw=subprocess.run(["zcat",f"{VCFDIR}/filtered_snv_tp_{chrom}.vcf.gz"],capture_output=True,text=True).stdout.splitlines()
    hdr=[l for l in raw if l.startswith("#")]; body=[l for l in raw if not l.startswith("#") and l.split("\t")[1] in poss]
    mini=f"{OUTD}/{chrom}.vcf"; open(mini,"w").write("\n".join(hdr)+"\n"+"\n".join(body)+"\n")
    subprocess.run(["bgzip","-f",mini]); subprocess.run(["tabix","-f","-p","vcf",mini+".gz"])
    od=f"{OUTD}/run_{chrom}"; os.makedirs(od)
    subprocess.run([BIN,"-t",TUMOR,"-n",NORMAL,"-r",REF,"-v",mini+".gz","-w","5000","-j","16",
        "--distance-metric","BERNOULLI","--nan-distance-strategy","SKIP","-o",od],stdout=open(f"{od}/log","w"),stderr=subprocess.STDOUT)
    for mp in glob.glob(f"{od}/**/distance/BERNOULLI/matrix.csv",recursive=True):
        rd=mp.rsplit("/distance/",1)[0]
        for part in rd.split("/"):
            if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd
    print(f"  {chrom}: binary done, {len([1 for p in picks if p['chrom']==chrom])} picks",flush=True)

def grpcol(g): return EDGE if g=="edge" else OUTL if g=="outlier" else H.CLUSTER_COL[int(g)%len(H.CLUSTER_COL)]
def render(key,sub,P,reads,ids,a,cpgs,pos,cl):
    fine={x:g for x,g in zip(ids,[p["grp"] for p in (a["perread_fine"] or [{"grp":"outlier"}]*len(ids))])}
    coarse={x:g for x,g in zip(ids,[p["grp"] for p in (a["perread_coarse"] or [{"grp":"outlier"}]*len(ids))])}
    Z,_=CR.linkZ(sub); dn=dendrogram(Z,orientation="left",no_plot=True); order=dn["leaves"][::-1]
    ids_o=[ids[i] for i in order]; meth=np.array([P[i] for i in order]); dist=sub[np.ix_(order,order)].copy()
    np.fill_diagonal(dist,0); dist[dist<0]=np.nan
    sb=[("fine",[grpcol(fine.get(x,"outlier")) for x in ids_o]),("coarse",[grpcol(coarse.get(x,"outlier")) for x in ids_o])]
    sb+=H.sidebar_specs({x:reads[x] for x in ids_o},ids_o,cluster_of=None,include_tn=False,include_strand=True)
    mc,dc=H.mpl_cmaps(); nsb=len(sb); n,ncol=meth.shape
    wr=[1.3]+[0.05]*nsb+[1.0,0.16]+[0.05]*nsb+[1.0]
    fig=plt.figure(figsize=(11,5.0)); gs=fig.add_gridspec(2,len(wr),width_ratios=wr,height_ratios=[1,0.4],wspace=0.04,hspace=0.02)
    c=0; axdn=fig.add_subplot(gs[0,c]); c+=1
    dendrogram(Z,orientation="left",above_threshold_color="#999",ax=axdn,no_labels=True)
    axdn.set_xticks([]); axdn.set_yticks([]); axdn.set_title("UPGMA 樹",fontsize=8); [axdn.spines[s].set_visible(False) for s in axdn.spines]
    for lab,hx in sb: H._sb(fig.add_subplot(gs[0,c]),hx,lab); c+=1
    axm=fig.add_subplot(gs[0,c]); c+=1; axm.imshow(meth,aspect="auto",cmap=mc,vmin=0,vmax=1,interpolation="nearest")
    sx=H.snv_fractional_x(cpgs,int(pos))
    if sx is not None: axm.axvline(sx,color=H.SNV_COL,lw=2)
    axm.set_xticks([]); axm.set_yticks([]); axm.set_title("甲基 read×CpG",fontsize=8.5); axm.set_xlabel(f"{ncol}CpG·SNV橙",fontsize=7)
    fig.add_subplot(gs[0,c]).axis("off"); c+=1
    for lab,hx in sb: H._sb(fig.add_subplot(gs[0,c]),hx,lab); c+=1
    axd=fig.add_subplot(gs[0,c]); vmax=max(0.5,float(np.nanmax(dist)) if np.isfinite(np.nanmax(dist)) else 0.5)
    axd.imshow(dist,aspect="auto",cmap=dc,vmin=0,vmax=vmax,interpolation="nearest")
    axd.set_xticks([]); axd.set_yticks([]); axd.set_title("距離(對角塊=群)",fontsize=8.5); axd.set_xlabel("暗=近",fontsize=7)
    H.grouped_legend(fig.add_subplot(gs[1,:]),[s[0] for s in sb],0)
    fig.suptitle(f"[{cl}] {key} n={n} | coarse k={a['coarse_k']}({a['coarse_confidence']}) fine k={a['fine_k']}({a['fine_confidence']})",fontsize=9.5,y=1.0)
    fn=f"{FIGS}/ws_{key}.png"; fig.savefig(fn,dpi=108,bbox_inches="tight"); plt.close(fig); return f"figs_ws/ws_{key}.png"

rng=np.random.default_rng(20260622); items=[]
for cl in CLASSES:
    for r in sel[cl]:
        key=f"{r['chrom']}_{r['pos']}"; rd=dirmap.get(key)
        if not rd: continue
        reads={x["read_id"]:x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
        dids,D=CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={x:i for i,x in enumerate(dids)}
        rows=open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); cpgs=[int(c) for c in rows[0].split(",")[1:]]
        mi={};M=[]
        for j,ln in enumerate(rows[1:]):
            q=ln.split(","); mi[q[0]]=j; M.append([np.nan if v in("","NA","nan","NaN") else float(v) for v in q[1:]])
        M=np.array(M); it=lambda t:str(t) in("1","true","True")
        ids=[x for x in dids if x in reads and it(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
        if len(ids)<CR.MIN_SZ*2: continue
        sub=D[np.ix_([di[x] for x in ids],[di[x] for x in ids])]; kp=CR.peel(sub)
        if len(kp)<CR.MIN_SZ*2: continue
        ids=[ids[i] for i in kp]; sub=sub[np.ix_(kp,kp)]
        P=np.array([M[mi[x]] for x in ids]); hp=[CR.LABMAP[reads[x]["hp"]] for x in ids]; al=[reads[x]["alt_support"] for x in ids]
        a=CR.analyze_locus(sub,P,hp,al,rng)
        # 取 fine 解析度細節
        fk=a["fine_k"]; fres=next((x for x in a["new_resolutions"] if x["k_core"]==fk),None)
        bestax=("-",None,None)
        if fres:
            for ax,d in fres["align"].items():
                if d["V"] is not None and (bestax[1] is None or d["V"]>bestax[1]): bestax=(ax,d["V"],d["e"])
        png=render(key,sub,P,reads,ids,a,cpgs,r["pos"],cl)
        from collections import Counter
        prt=Counter(p["type"] for p in (a["perread_fine"] or []))
        items.append({"key":key,"class":cl,"n":len(ids),"coarse_k":a["coarse_k"],"coarse_conf":a["coarse_confidence"],
            "fine_k":fk,"fine_conf":a["fine_confidence"],"excess":(fres["stab_excess"] if fres else None),
            "gap_ratio":(fres["gap_ratio"] if fres else None),"big_gap":(fres["big_gap"] if fres else None),
            "align_axis":bestax[0],"align_V":bestax[1],"align_e":bestax[2],
            "core_sizes":(fres["core_sizes"] if fres else []),"n_outliers":(fres["n_outliers"] if fres else 0),
            "perread":{"core":prt.get("core_confirmed",0)+prt.get("core_novel",0),"edge":prt.get("edge",0),"outlier":prt.get("outlier",0)},
            "png":png})
        print(f"  {cl:16s} {key} n={len(ids)} fine={a['fine_confidence']} exc={fres['stab_excess'] if fres else None} V={bestax[1]} e={bestax[2]}",flush=True)
json.dump(items,open(f"{A}/ws_items.json","w"),indent=1)
shutil.rmtree(OUTD,ignore_errors=True)
print(f"DONE {len(items)} items → ws_items.json")
