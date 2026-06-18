#!/usr/bin/env python3
"""拆解『沒訊號』位點的原因並分類。沒訊號 = 無 all-read 結構(sil<0.4) 且 無 within-carrier split。
原因(互斥,優先序): 讀數不足 / 矩陣不完整 / 無平衡分群 / 弱結構sub-threshold / 很弱接近均勻。
+ 每類 cross pipeline VC (確認是否仍有 label-PERMANOVA 訊號)。"""
import json, csv
from collections import Counter, defaultdict
ASSET="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
SIG="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_wg_bdcprime_verify/significance_summary.csv"
TH_STRUCT=0.4; TH_WITHIN=0.5; MIN_SZ=3
NOISE={"Noise_Uniform","Noise_Uncorrelated","Noise_Chaotic"}
recs=json.load(open(f"{ASSET}/records_wg2.json")); N=len(recs)
vc_by_pos={r["Pos"]:r["VerificationClass"] for r in csv.DictReader(open(SIG))}

def has_struct(a): return a["best_k"] is not None and a["best_k"]>=2 and a["best_sil"] is not None and a["best_sil"]>=TH_STRUCT
def has_within(w): return any(w[l]["best_sil"] is not None and w[l]["best_sil"]>=TH_WITHIN for l in ("1-1","2-1"))

def reason(r):
    a=r["all"]
    if a["n"] < MIN_SZ*2:        return "1_讀數不足(n<6)"
    if a["n_complete"] is None or a["n_complete"] < MIN_SZ*2: return "2_矩陣不完整(NaN/覆蓋)"
    if a["best_k"] is None:      return "3_無平衡分群(切不出≥2群各≥3)"
    s=a["best_sil"]
    if s is None:                return "3_無平衡分群(切不出≥2群各≥3)"
    if s>=0.3:                   return "4_弱結構sub-threshold(0.3≤sil<0.4)"
    return "5_很弱/接近均勻(sil<0.3)"

nosig=[r for r in recs if not has_struct(r["all"]) and not has_within(r["within"])]
cats=Counter(reason(r) for r in nosig)
NS=len(nosig)
print(f"全位點 N={N}; 沒訊號 = {NS} ({round(100*NS/N,1)}%)")
print(f"\n=== 沒訊號的原因分類 (互斥) ===")
print(f"{'原因':<34}{'位點':>7}{'佔沒訊號':>9}{'佔全WG':>8}")
ABSTAIN=0; WEAK=0
for k in sorted(cats):
    v=cats[k]
    print(f"{k:<34}{v:>7}{round(100*v/NS,1):>8}%{round(100*v/N,1):>7}%")
    if k[0] in "123": ABSTAIN+=v
    else: WEAK+=v
print(f"\n  ├ 『根本問不了』(讀數/覆蓋/切不出, 類1-3): {ABSTAIN} ({round(100*ABSTAIN/NS,1)}% of 沒訊號, {round(100*ABSTAIN/N,1)}% of WG)")
print(f"  └ 『問了但結構弱』(類4-5):                {WEAK} ({round(100*WEAK/NS,1)}% of 沒訊號, {round(100*WEAK/N,1)}% of WG)")

print(f"\n=== 每類 × pipeline VC (確認是否仍有 label-PERMANOVA 訊號) ===")
print(f"{'原因':<34}{'VC=Noise':>10}{'非-Noise(有label結構)':>22}")
for k in sorted(cats):
    sub=[r for r in nosig if reason(r)==k]
    nv=sum(1 for r in sub if vc_by_pos.get(r["pos"]) in NOISE)
    nn=len(sub)-nv-sum(1 for r in sub if vc_by_pos.get(r["pos"]) in ("?",None))
    print(f"{k:<34}{nv:>10}{nn:>16} ({round(100*nn/len(sub),0):.0f}%)")
totnoise=sum(1 for r in nosig if vc_by_pos.get(r["pos"]) in NOISE)
print(f"\n  沒訊號中 pipeline VC=Noise(真均勻): {totnoise} ({round(100*totnoise/NS,1)}%)")
print(f"  沒訊號中 pipeline 非-Noise(有 label-PERMANOVA 訊號但無乾淨無監督群): {NS-totnoise} ({round(100*(NS-totnoise)/NS,1)}%)")

# save
json.dump({"N":N,"nosig":NS,"reasons":dict(cats),"abstain":ABSTAIN,"weak":WEAK,
  "nosig_vc_noise":totnoise,"nosig_vc_nonnoise":NS-totnoise},
  open(f"{ASSET}/nosignal_breakdown.json","w"),ensure_ascii=False,indent=1)
print(f"\n[-> nosignal_breakdown.json]")
