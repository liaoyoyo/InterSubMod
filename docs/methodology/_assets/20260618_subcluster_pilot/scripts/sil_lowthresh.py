#!/usr/bin/env python3
"""補 0.1 以上門檻 + 證明『有 sil 值 ≠ 有訊號』: best_sil 分布 by VC(Noise vs 結構)。
若 Noise 位點也有 silhouette mass 在低-中段 → 低 sil 與噪音重疊 → sil 值本身非訊號判準。"""
import json, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
ASSET="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
FIGS=f"{ASSET}/figs"
SIG="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_wg_bdcprime_verify/significance_summary.csv"
recs=json.load(open(f"{ASSET}/records_wg2.json")); N=len(recs)
NOISE={"Noise_Uniform","Noise_Uncorrelated","Noise_Chaotic"}
vc_by_pos={r["Pos"]:r["VerificationClass"] for r in csv.DictReader(open(SIG))}
MIN_CELL=3; TH_WITHIN=0.5
def has_within(w): return any(w[l]["best_sil"] is not None and w[l]["best_sil"]>=TH_WITHIN for l in ("1-1","2-1"))
def corr(a,th):
    hs=a["best_k"] is not None and a["best_k"]>=2 and a["best_sil"] is not None and a["best_sil"]>=th
    if not(hs and a["contingency"]): return hs,False,False
    clu=defaultdict(set); lab=defaultdict(set)
    for L,d in a["contingency"].items():
        for cl,c in d.items():
            if c>=MIN_CELL: clu[cl].add(L); lab[L].add(cl)
    return hs,(not any(len(s)>=2 for s in clu.values()) and not any(len(s)>=2 for s in lab.values())),any(len(s)>=2 for s in lab.values())

# ---- extended sensitivity 0.1..0.5 ----
THS=[0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50]
rows=[]
for th in THS:
    hs=o2o=m1=0
    for r in recs:
        a,b,c=corr(r["all"],th); hs+=a; o2o+=b; m1+=c
    nosig=sum(1 for r in recs if not(r["all"]["best_k"] is not None and r["all"]["best_k"]>=2 and r["all"]["best_sil"] is not None and r["all"]["best_sil"]>=th) and not has_within(r["within"]))
    rows.append({"th":th,"has_struct":hs,"one2one":o2o,"many2one":m1,"nosig":nosig})
json.dump({"N":N,"rows":rows},open(f"{ASSET}/sensitivity_ext.json","w"),indent=1)
print(f"=== 擴充門檻 0.1→0.5 (% of {N}) ===")
print(f"{'TH':>5}{'有結構':>9}{'1對1':>8}{'多對1':>8}{'沒訊號':>9}")
for r in rows: print(f"{r['th']:>5}{100*r['has_struct']/N:>8.1f}%{100*r['one2one']/N:>7.1f}%{100*r['many2one']/N:>7.1f}%{100*r['nosig']/N:>8.1f}%")

# ---- best_sil by VC (Noise vs 結構) — 證明重疊 ----
sil_noise=[r["all"]["best_sil"] for r in recs if r["all"]["best_sil"] is not None and vc_by_pos.get(r["pos"]) in NOISE]
sil_struct=[r["all"]["best_sil"] for r in recs if r["all"]["best_sil"] is not None and vc_by_pos.get(r["pos"]) not in NOISE and vc_by_pos.get(r["pos"]) not in("?",None)]
def q(x):
    if not x: return "n=0"
    x=sorted(x); n=len(x); return f"n={n} med={x[n//2]:.3f} q25={x[n//4]:.3f} q75={x[3*n//4]:.3f}"
print(f"\n=== best_sil 分布 by pipeline VC (有切出群者) ===")
print(f"  VC=Noise:   {q(sil_noise)}")
print(f"  VC=結構:    {q(sil_struct)}")
for th in (0.3,0.4,0.5):
    nn=sum(1 for s in sil_noise if s>=th); ss=sum(1 for s in sil_struct if s>=th)
    print(f"  sil>={th}: Noise {nn}/{len(sil_noise)} ({round(100*nn/len(sil_noise),0):.0f}%) vs 結構 {ss}/{len(sil_struct)} ({round(100*ss/len(sil_struct),0):.0f}%)")

fig,ax=plt.subplots(figsize=(6.6,3.5))
plt.rcParams["font.sans-serif"]=["Noto Sans CJK TC","DejaVu Sans"]; plt.rcParams["axes.unicode_minus"]=False
bins=np.linspace(0,1,33)
ax.hist(sil_struct,bins=bins,density=True,alpha=0.55,color="#4A77A8",label=f"結構-VC (n={len(sil_struct)})")
ax.hist(sil_noise,bins=bins,density=True,alpha=0.55,color="#C0432F",label=f"Noise-VC (n={len(sil_noise)})")
for th,c in [(0.4,"#C0432F")]: ax.axvline(th,color=c,lw=1.5,ls="--")
ax.text(0.4,ax.get_ylim()[1]*0.9,"0.4",color="#C0432F",fontsize=9,ha="center")
ax.set_xlabel("all-read best silhouette",fontsize=9); ax.set_ylabel("density",fontsize=9)
ax.set_title("關鍵: Noise-VC 位點也有相當 silhouette mass\n→ 低-中 sil 與噪音重疊 = 有 sil 值本身不代表有訊號",fontsize=9.5)
ax.legend(fontsize=8.5); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{FIGS}/sil_noise_overlap.png",dpi=100); plt.close(fig)
print("\nWROTE figs/sil_noise_overlap.png + sensitivity_ext.json")
# answer metric
ov=sum(1 for s in sil_noise if 0.2<=s<0.5)
print(f"\n答: VC=Noise 中 best_sil 落 0.2-0.5 的有 {ov}/{len(sil_noise)} ({round(100*ov/len(sil_noise),0):.0f}%) → 低-中 sil 不可靠")
