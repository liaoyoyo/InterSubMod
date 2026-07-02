#!/usr/bin/env python3
"""從 records_wg2.json 分類每位點 1對1/1對多/多對1 + 有訊號/沒訊號。門檻可調,便宜重算。
correspondence 只在『有結構』(all-read sil>=TH_struct)的聯列表上算(consistent Venn);
within-label carrier split = subclone-direction refinement (獨立)。
用法: classify_contingency.py [TH_within=0.5] [TH_struct=0.4] [MIN_CELL=3]"""
import sys, json, csv
from collections import Counter, defaultdict
ASSET="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
SIG="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_wg_bdcprime_verify/significance_summary.csv"
TH_WITHIN=float(sys.argv[1]) if len(sys.argv)>1 else 0.5
TH_STRUCT=float(sys.argv[2]) if len(sys.argv)>2 else 0.4
MIN_CELL=int(sys.argv[3]) if len(sys.argv)>3 else 3
NOISE={"Noise_Uniform","Noise_Uncorrelated","Noise_Chaotic"}

recs=json.load(open(f"{ASSET}/records_wg2.json"))
vc_by_pos={}
try:
    for r in csv.DictReader(open(SIG)): vc_by_pos[r["Pos"]]=r["VerificationClass"]
except FileNotFoundError: pass

def classify(r):
    a=r["all"]; w=r["within"]
    has_struct = a["best_k"] is not None and a["best_k"]>=2 and a["best_sil"] is not None and a["best_sil"]>=TH_STRUCT
    has_within_carrier = any(w[l]["best_sil"] is not None and w[l]["best_sil"]>=TH_WITHIN for l in ("1-1","2-1"))
    has_1many=has_many1=one2one=False
    if has_struct and a["contingency"]:
        clu_labels=defaultdict(set); lab_clusters=defaultdict(set)
        for lab,d in a["contingency"].items():
            for cl,cnt in d.items():
                if cnt>=MIN_CELL: clu_labels[cl].add(lab); lab_clusters[lab].add(cl)
        has_1many=any(len(s)>=2 for s in clu_labels.values())     # 一群多標籤
        has_many1=any(len(s)>=2 for s in lab_clusters.values())   # 一標籤多群(joint)
        one2one = not has_1many and not has_many1                 # 乾淨對齊
    sig_strict = one2one
    sig_loose  = one2one or has_many1 or has_within_carrier
    sig_any    = has_struct or has_within_carrier
    return dict(has_struct=has_struct,has_many1=has_many1,has_within_carrier=has_within_carrier,
                has_1many=has_1many,one2one=one2one,sig_strict=sig_strict,sig_loose=sig_loose,sig_any=sig_any)

cls=[classify(r) for r in recs]; N=len(recs)
def pct(x): return f"{x} ({round(100*x/N,1)}%)"
g=lambda key:sum(c[key] for c in cls)
print(f"N={N}  TH_within={TH_WITHIN} TH_struct={TH_STRUCT} MIN_CELL={MIN_CELL}")
print("\n=== cluster×label 對應 (位點數; 只在有結構聯列表上算) ===")
print(f"  有 all-read 結構(k>=2,sil>={TH_STRUCT}): {pct(g('has_struct'))}")
print(f"  1對1 (乾淨對齊):          {pct(g('one2one'))}")
print(f"  1對多 (一群多標籤):       {pct(g('has_1many'))}")
print(f"  多對1 (一標籤多群,joint): {pct(g('has_many1'))}")
print(f"  within-carrier split(subclone方向,獨立): {pct(g('has_within_carrier'))}")
both=sum(c['has_1many'] and c['has_many1'] for c in cls)
print(f"  1對多 ∩ 多對1:            {pct(both)}")
print(f"\n=== 文氏圖 (有結構 {g('has_struct')} 內: A=多對1 B=1對多) ===")
A=set(i for i,c in enumerate(cls) if c['has_many1']); B=set(i for i,c in enumerate(cls) if c['has_1many'])
print(f"  1對1(乾淨): {pct(g('one2one'))}")
print(f"  只 多對1: {pct(len(A-B))}   只 1對多: {pct(len(B-A))}   交集: {pct(len(A&B))}")
print(f"  無結構: {pct(N-g('has_struct'))}")
print(f"\n=== 有訊號 vs 沒訊號 (以位點為單位) ===")
print(f"  有訊號(嚴格=1對1可判別):      {pct(g('sig_strict'))}")
print(f"  有訊號(+多對1/within-carrier): {pct(g('sig_loose'))}")
print(f"  有訊號(任何結構):             {pct(g('sig_any'))}")
print(f"  沒訊號:                       {pct(N-g('sig_any'))}")
if vc_by_pos:
    nosig=[recs[i] for i,c in enumerate(cls) if not c['sig_any']]
    nv=sum(1 for r in nosig if vc_by_pos.get(r['pos']) in NOISE)
    pstruct=sum(1 for r in recs if vc_by_pos.get(r['pos']) not in NOISE and vc_by_pos.get(r['pos']) not in ('?',None))
    print(f"\n=== 對照 pipeline VC ===")
    print(f"  我『沒訊號』{len(nosig)} 中 pipeline VC=Noise: {nv} ({round(100*nv/len(nosig),1) if nosig else 0}%)")
    print(f"  pipeline 非-Noise(有label-PERMANOVA結構): {pstruct} ({round(100*pstruct/N,1)}%) — 對比我無監督 8% (定義不同: label-PERMANOVA vs 無監督群)")
json.dump({"N":N,"th_within":TH_WITHIN,"th_struct":TH_STRUCT,"min_cell":MIN_CELL,
  "has_struct":g('has_struct'),"one2one":g('one2one'),"one2many":g('has_1many'),"many2one":g('has_many1'),
  "within_carrier":g('has_within_carrier'),"both":both,"onlyA":len(A-B),"onlyB":len(B-A),"AB":len(A&B),
  "nostruct":N-g('has_struct'),"sig_strict":g('sig_strict'),"sig_loose":g('sig_loose'),"sig_any":g('sig_any'),
  "nosig":N-g('sig_any')},open(f"{ASSET}/contingency_summary.json","w"),indent=1)
print(f"\n[summary -> contingency_summary.json]")
