#!/usr/bin/env python3
"""paired(pipeline GlobalTest) vs tumor-only(我的 cluster×label pass) 的 Fisher/CramersV 對照。
paired: significance_summary.csv GlobalP/CramersV(allele軸) + HPFine軸. tumor-only: records_wg2 列聯表算 CramersV+chi2p."""
import json, csv
import numpy as np
from scipy.stats import chi2_contingency
ASSET="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
SUM="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_wg_bdcprime_verify/significance_summary.csv"
def fl(x):
    try: return float(x)
    except: return None
rows=list(csv.DictReader(open(SUM))); N=len(rows)
print(f"N = {N}\n")

def dist(vals,name,is_p=False):
    v=[x for x in vals if x is not None]
    print(f"  {name}: 有值 {len(v)}")
    if is_p:
        for th in (0.05,0.01,0.001):
            s=sum(1 for x in v if x<th); print(f"    p<{th}: {s} ({round(100*s/N,1)}% of WG, {round(100*s/len(v),1)}% of 有值)")
    else:
        a=np.array([x for x in v if x>0]); z=sum(1 for x in v if x==0)
        print(f"    =0: {z} ({round(100*z/N,1)}%); >0: {len(a)} ({round(100*len(a)/N,1)}%)" + (f" med={np.median(a):.3f} q75={np.percentile(a,75):.3f}" if len(a) else ""))

# ===== PAIRED (pipeline) =====
print("=== PAIRED (pipeline GlobalTest, tumor+normal reads) ===")
print(" [allele 軸 ALT/REF — 預設 GlobalP/CramersV]")
dist([fl(r["GlobalP"]) for r in rows],"GlobalP(Fisher-FFH)",is_p=True)
dist([fl(r["CramersV"]) for r in rows],"CramersV")
print(" [hp_fine 軸 1-1/2-1...]")
dist([fl(r["GlobalP_HPFine"]) for r in rows],"GlobalP_HPFine",is_p=True)
dist([fl(r["CramersV_HPFine"]) for r in rows],"CramersV_HPFine")

# ===== TUMOR-ONLY (my pass) =====
recs=json.load(open(f"{ASSET}/records_wg2.json"))
to_v=[]; to_p=[]; to_has=0
for r in recs:
    cont=r["all"]["contingency"]
    if not cont:
        to_v.append(None); to_p.append(None); continue
    labels=sorted(cont.keys()); clusters=sorted({c for d in cont.values() for c in d})
    M=np.array([[cont[L].get(c,0) for c in clusters] for L in labels],dtype=float)
    if M.shape[0]<2 or M.shape[1]<2 or M.sum()<6:
        to_v.append(None); to_p.append(None); continue
    try:
        chi2,p,_,_=chi2_contingency(M)
        n=M.sum(); v=np.sqrt(chi2/(n*(min(M.shape)-1)))
        to_v.append(round(float(v),3)); to_p.append(float(p)); to_has+=1
    except Exception:
        to_v.append(None); to_p.append(None)
print(f"\n=== TUMOR-ONLY (我的 cluster×label pass, tumor reads, haplotag 1/1-1/2/2-1) ===")
print(f"  可算列聯表(tumor-only 有切出群): {to_has} ({round(100*to_has/N,1)}% of WG)")
dist(to_p,"chi2-p(近似 Fisher)",is_p=True)
dist(to_v,"CramersV(tumor-only)")

print(f"\n=== 對照摘要 ===")
pg_sig=sum(1 for r in rows if fl(r['GlobalP']) is not None and fl(r['GlobalP'])<0.05)
hpf_sig=sum(1 for r in rows if fl(r['GlobalP_HPFine']) is not None and fl(r['GlobalP_HPFine'])<0.05)
to_sig=sum(1 for x in to_p if x is not None and x<0.05)
print(f"  顯著(p<0.05): PAIRED-allele {round(100*pg_sig/N,1)}% | PAIRED-hpfine {round(100*hpf_sig/N,1)}% | TUMOR-only {round(100*to_sig/N,1)}%")
toa=np.array([x for x in to_v if x and x>0])
print(f"  CramersV>0 med: PAIRED {np.median([fl(r['CramersV']) for r in rows if fl(r['CramersV']) and fl(r['CramersV'])>0]):.3f} | TUMOR-only {np.median(toa) if len(toa) else 0:.3f}")
json.dump({"N":N,"paired_allele_sig":pg_sig,"paired_hpfine_sig":hpf_sig,"tumor_only_sig":to_sig,
  "tumor_only_computable":to_has},open(f"{ASSET}/tumor_vs_paired.json","w"),indent=1)
print("\n[-> tumor_vs_paired.json]")
