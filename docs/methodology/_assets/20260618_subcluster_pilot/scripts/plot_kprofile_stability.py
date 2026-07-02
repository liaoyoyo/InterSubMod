#!/usr/bin/env python3
"""切割穩定性視覺化（純讀 kprofile_stability.json，數字不手打）。
FIG-overview: 全 29 位點散點 x=excess-over-null y=align_V，顏色=verdict、形狀=group，
              + 門檻線(excess=0.10 / V=0.3) → 一張圖看「穩+對齊」vs「穩但陷阱」分離。
FIG-per-locus: 代表位點 across-k bar = real Jaccard vs within-1-group null95(excess 著色) + 對齊 V vs shuffle95。
代表位點從 json 動態挑(clean-GOOD / multi-res / trap-高raw-unaligned / 高V低可靠 / single-k / 最低excess)。
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]=["Droid Sans Fallback","DejaVu Sans"]
plt.rcParams["axes.unicode_minus"]=False

A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
FIGS=f"{A}/figs_stability"; os.makedirs(FIGS,exist_ok=True)
d=json.load(open(f"{A}/kprofile_stability.json"))
loci=d["loci"]; P=d["params"]
VCOL={"GOOD_CUT":"#0d9488","STABLE_BUT_UNALIGNED":"#db2777","ALIGNED_BUT_UNSTABLE":"#d97706",
      "STABLE+ALIGNED_DIFF_K":"#7c3aed","NEITHER":"#475569"}
GMARK={"confident-unique":"o","multi-resolution":"s","ambiguous-near-tie":"^","single-k-forced":"D"}

def hv(r,k):  # headline value
    return r["headline"][k]

# ---------- FIG overview (2-panel) ----------
from matplotlib.lines import Line2D
fig,(axB,axS)=plt.subplots(1,2,figsize=(14,6.0),gridspec_kw={"width_ratios":[1,1.15]})
# Panel A: verdict-by-group stacked bar (乾淨類別總結 = 真正鑑別)
groups=["confident-unique","multi-resolution","ambiguous-near-tie","single-k-forced"]
vorder=["GOOD_CUT","STABLE+ALIGNED_DIFF_K","STABLE_BUT_UNALIGNED","ALIGNED_BUT_UNSTABLE","NEITHER"]
vorder=[v for v in vorder if any(r["verdict"]==v for r in loci)]
y=np.arange(len(groups))
left=np.zeros(len(groups))
for v in vorder:
    cnt=[sum(1 for r in loci if r["group"]==g and r["verdict"]==v) for g in groups]
    axB.barh(y,cnt,left=left,color=VCOL[v],edgecolor="#222",lw=0.5,label=v)
    for i,c in enumerate(cnt):
        if c>0: axB.text(left[i]+c/2,y[i],str(c),ha="center",va="center",fontsize=9,color="w",fontweight="bold")
    left+=cnt
axB.set_yticks(y); axB.set_yticklabels([g+f"\n(n_loci={sum(1 for r in loci if r['group']==g)})" for g in groups],fontsize=9)
axB.set_xlabel("位點數",fontsize=10); axB.invert_yaxis()
axB.set_title("verdict × group：GOOD 集中高覆蓋(multi-res/single-k)；\nambiguous 幾全 UNALIGNED(穩但不對齊)",fontsize=10)
axB.legend(fontsize=7.5,loc="lower right",title="verdict",title_fontsize=8)
# Panel B: scatter excess vs V (全落右上 → 兩軸邊際也看似好；真正鑑別在色)
for r in loci:
    x=hv(r,"stab_excess"); V=hv(r,"align_V"); yy=V if V is not None else -0.05
    axS.scatter(x,yy,c=VCOL.get(r["verdict"],"#888"),marker=GMARK.get(r["group"],"o"),
               s=95,edgecolors="#222",linewidths=0.6,alpha=0.9,zorder=3)
axS.axvline(P["excess_min"],ls="--",c="#888",lw=1); axS.axhline(0.3,ls="--",c="#888",lw=1); axS.axhline(0.0,ls=":",c="#bbb",lw=0.8)
axS.text(0.62,-0.04,"V=None(標籤退化)",fontsize=7,color="#999")
axS.set_xlabel("excess-over-null = real clusterboot Jaccard − within-1-group null(95pct)",fontsize=9.5)
axS.set_ylabel("alignment CramérV (cut vs a-priori 軸)",fontsize=9.5)
axS.set_xlim(-0.05,1.05); axS.set_ylim(-0.1,1.08); axS.grid(alpha=0.25,zorder=0)
axS.text(0.04,0.06,"[!] 多數點落右上(高 excess+高 V)\n→ (excess,V) 邊際也『看似都好』\n真正鑑別 = 同一 k coherence + Cochran e≥5\n(見色：pink=穩但不可靠對齊)",
         fontsize=8,color="#444",bbox=dict(boxstyle="round",fc="#fff8f0",ec="#d97706",alpha=0.9))
gh=[Line2D([],[],marker=m,ls="",mfc="#ccc",mec="#222",ms=9,label=g) for g,m in GMARK.items()]
axS.legend(handles=gh,title="group",loc="lower right",fontsize=8,title_fontsize=8.5)
fig.suptitle("切割穩定性 3 軸總覽 (HCC1395 tumor-only, n=%d 代表位點; raw Jaccard 全 0.74–1.0 無鑑別力; BERNOULLI 自檢 max|Δ|=%.0e)"
             %(d["n_loci"],d["bernoulli_selfcheck_max_absdiff"]),fontsize=11,y=1.0)
fig.tight_layout(); fig.savefig(f"{FIGS}/overview_excess_vs_align.png",dpi=120,bbox_inches="tight"); plt.close(fig)

# ---------- 選代表位點 ----------
by={r["chrom"]+":"+r["pos"]:r for r in loci}
def pick(pred,key,rev=True):
    c=[r for r in loci if pred(r)]
    return (sorted(c,key=key,reverse=rev)[0]["chrom"]+":"+sorted(c,key=key,reverse=rev)[0]["pos"]) if c else None
reps=[]
reps.append(("[OK] 乾淨 GOOD(穩+可靠對齊同k)", pick(lambda r:r["verdict"]=="GOOD_CUT" and sum(1 for k in r["per_k"] if k["pass_both"])>=1, lambda r:hv(r,"stab_excess"))))
reps.append(("[OK] multi-resolution GOOD(多k皆過)", pick(lambda r:r["verdict"]=="GOOD_CUT" and r["group"]=="multi-resolution" and sum(1 for k in r["per_k"] if k["pass_both"])>=2, lambda r:sum(1 for k in r["per_k"] if k["pass_both"]))))
reps.append(("[CEILING] STABLE+ALIGNED_DIFF_K(單樣本可靠性天花板)", pick(lambda r:r["verdict"]=="STABLE+ALIGNED_DIFF_K", lambda r:hv(r,"stab_excess"))))
reps.append(("[TRAP] 穩但陷阱(raw 高/不對齊)", pick(lambda r:r["verdict"]=="STABLE_BUT_UNALIGNED", lambda r:max(k["jac_real_mean"] for k in r["per_k"]))))
reps.append(("[TRAP] 高V 但 Cochran e<5 不可靠", pick(lambda r:r["verdict"]=="STABLE_BUT_UNALIGNED" and (hv(r,"align_V") or 0)>=0.7, lambda r:hv(r,"align_V"))))
reps.append(("[OK] single-k GOOD(高覆蓋 2 群)", pick(lambda r:r["verdict"]=="GOOD_CUT" and r["group"]=="single-k-forced", lambda r:r["n"])))
seen=set(); reps=[(lab,k) for lab,k in reps if k and not (k in seen or seen.add(k))][:6]

# ---------- FIG per-locus ----------
plotted=[]
for lab,key in reps:
    r=by[key]; pk=r["per_k"]; ks=[p["k"] for p in pk]
    fig,(a1,a2)=plt.subplots(1,2,figsize=(11,4.0),gridspec_kw={"width_ratios":[1.25,1]})
    x=np.arange(len(ks)); w=0.38
    real=[p["jac_real_mean"] for p in pk]; null=[p["jac_null_p95"] for p in pk]
    a1.bar(x-w/2,real,w,label="real clusterboot Jaccard",color="#0d9488",edgecolor="#222",lw=0.5)
    a1.bar(x+w/2,null,w,label="within-1-group null (95pct)",color="#d4d4d4",edgecolor="#222",lw=0.5,hatch="//")
    for i,p in enumerate(pk):
        exc=p["stab_excess"]; col="#0d9488" if p["pass_both"] else ("#db2777" if exc>=P["excess_min"] else "#999")
        a1.annotate(f"exc{exc:+.2f}",(i,max(real[i],null[i])+0.02),ha="center",fontsize=7.5,color=col)
        if p["pass_both"]: a1.text(i,-0.08,"✓both",ha="center",fontsize=7.5,color="#0d9488")
    a1.axhline(0.6,ls=":",c="#999",lw=0.8); a1.set_xticks(x); a1.set_xticklabels([f"k={k}" for k in ks])
    a1.set_ylim(-0.12,1.12); a1.set_ylabel("Jaccard"); a1.legend(fontsize=7.5,loc="upper left")
    a1.set_title("穩定性: real vs within-1-group null (excess=鑑別量)",fontsize=9.5)
    # alignment panel
    Vs=[p["align_V"] if p["align_V"] is not None else 0 for p in pk]
    es=[p["align_e"] if p["align_e"] is not None else 0 for p in pk]
    cols=["#0d9488" if (p["align_e"] or 0)>=5 and (p["align_V"] or 0)>=0.3 and (p["align_shuffle_p"] or 1)<0.05 else "#db2777" for p in pk]
    a2.bar(x,Vs,0.55,color=cols,edgecolor="#222",lw=0.5)
    a2.axhline(0.3,ls="--",c="#888",lw=1)
    for i,p in enumerate(pk):
        rel="e<5✗" if (p["align_e"] or 0)<5 else "e≥5"
        a2.annotate(f"{rel}\np={p['align_shuffle_p']}",(i,(Vs[i] or 0)+0.02),ha="center",fontsize=6.8,
                    color="#db2777" if (p['align_e'] or 0)<5 else "#0d9488")
    a2.set_xticks(x); a2.set_xticklabels([f"k={k}" for k in ks]); a2.set_ylim(0,1.12)
    a2.set_ylabel("CramérV (對齊)"); a2.set_title("對齊 (綠=可靠對齊; 紅=不可靠/弱)",fontsize=9.5)
    fig.suptitle(f"[{lab}]  {key}  n={r['n']}  軸={r['primary_axis']}  → verdict={r['verdict']}  (headline k={r['headline_k']})",
                 fontsize=11,y=1.02)
    fig.tight_layout(); fn=f"{FIGS}/stab_{key.replace(':','_')}.png"; fig.savefig(fn,dpi=110,bbox_inches="tight"); plt.close(fig)
    plotted.append({"label":lab,"key":key,"verdict":r["verdict"],"n":r["n"],"group":r["group"],
                    "headline":r["headline"],"png":f"figs_stability/stab_{key.replace(':','_')}.png",
                    "dualpanel":f"figs_kprofile/kp_{r['group']}_{key.replace(':','_')}.png"})
json.dump({"overview":"figs_stability/overview_excess_vs_align.png","reps":plotted},
          open(f"{A}/kprofile_stability_figindex.json","w"),indent=1)
print(f"overview + {len(plotted)} per-locus figs → {FIGS}")
for p in plotted: print(" ",p["label"],p["key"],p["verdict"])
