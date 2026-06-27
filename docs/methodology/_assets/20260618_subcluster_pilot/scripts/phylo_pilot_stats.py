#!/usr/bin/env python3
"""[統計] 30 pilot 位點 v3 完整分類統計 + 邊界案例。
有結構/無結構比例 · 群數分佈 · 分層深度 · 對齊型態(allele軸/hp軸/within-germline) · 4-gate vs v3 交叉 · 邊界案例。
純讀 json(L1)。"""
import json
from collections import Counter
A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
v3 = json.load(open(f"{A}/phylo_v3_render.json"))  # v2_ngroups 欄=v3 群數
ws = {x["key"]: x for x in json.load(open(f"{A}/ws_items.json"))}
N = len(v3)


def align_type(o):
    al = o["align"]
    if o["v2_ngroups"] < 2: return "no_structure"
    alleles = {al[L]["allele"].split("(")[0] for L in al}
    hps = {al[L]["hp"].split("(")[0] for L in al}
    if len(alleles) > 1: return "allele_axis"      # REF/ALT 分裂 = cis-ASM
    if len(hps) > 1: return "hp_axis"              # germline 單倍型軸
    return "within_germline"                       # 同 allele+hp = subclone 候選


for o in v3: o["atype"] = align_type(o); o["ng"] = o["v2_ngroups"]
multi = [o for o in v3 if o["ng"] >= 2]; single = [o for o in v3 if o["ng"] < 2]

print(f"=== A. 結構/無結構（N={N} pilot 位點，curated 非隨機）===")
print(f"有結構(≥2群): {len(multi)} ({100*len(multi)/N:.1f}%) | 無結構(1群): {len(single)} ({100*len(single)/N:.1f}%)")
print(f"\n=== B. 群數分佈 ===")
gd = Counter(o["ng"] for o in v3)
for k in sorted(gd): print(f"  {k} 群: {gd[k]} ({100*gd[k]/N:.1f}%)")
print(f"\n=== C. 分層深度（多群位點，標籤最深 '-' 數）===")
dd = Counter(max((L.count('-') for L in o['labels']), default=0) for o in multi)
for k in sorted(dd): print(f"  深度 {k}(如 {'1' if k==0 else '1-'+'1-'*(k-1)+'1'}): {dd[k]} 位點")
print(f"\n=== D. 對齊型態分類（多群 {len(multi)} 位點）===")
at = Counter(o["atype"] for o in multi)
labels = {"allele_axis": "allele 軸 REF/ALT (cis-ASM)", "hp_axis": "hp 軸 germline 單倍型", "within_germline": "🔴 同 allele+hp (subclone 候選)"}
for t in ["allele_axis", "hp_axis", "within_germline"]:
    n = at.get(t, 0); print(f"  {labels[t]}: {n} ({100*n/len(multi):.0f}% of multi, {100*n/N:.1f}% of all)")
print(f"\n=== E. 4-gate fine_class vs v3 群數 交叉 ===")
print(f"{'4-gate class':>16} | {'位點':>4} | v3 多群數 | v3 單群數")
for cl in ["CONFIRMED", "NEAR_CONFIRMED", "REAL_NOVEL", "REAL_DIFFUSE", "NO_CLEAR_SPLIT"]:
    g = [o for o in v3 if o["fine_class"] == cl]
    if not g: continue
    mm = sum(1 for o in g if o["ng"] >= 2); ss = len(g) - mm
    print(f"{cl:>16} | {len(g):>4} | {mm:>8} | {ss:>8}")
print(f"\n=== F. 邊界案例（需特別觀察）===")
# F1 unstable: chr20_42981498 (sensitivity 證跨 seed 翻)
print("  [F1 verdict 不穩] chr20_42981498 — RNULL=80 仍 1/2/3 翻動（資料在決策邊界）→ 該標 ambiguous")
# F2 within-germline (subclone 候選)
wg = [o for o in multi if o["atype"] == "within_germline"]
print(f"  [F2 within-germline] {len(wg)} 位點同 allele+hp 多群（subclone 候選 vs cis-ASM 子結構需細看）: {[o['key'] for o in wg]}")
# F3 deep hierarchy (depth≥2)
deep = [o for o in multi if max((L.count('-') for L in o['labels']), default=0) >= 2]
print(f"  [F3 深層階層 depth≥2] {len(deep)} 位點: {[(o['key'], sorted(o['labels'])) for o in deep]}")
# F4 low CpG (FP 風險區)
lowc = [o for o in v3 if o["C"] <= 40]
print(f"  [F4 低 CpG(≤40, FP 略升區)] {len(lowc)} 位點: {[(o['key'], o['C'], o['ng']) for o in lowc]}")
# F5 大 n 高信心 multi
bign = [o for o in multi if o["n"] >= 60]
print(f"  [F5 大n(≥60)高信心多群] {len(bign)} 位點: {[(o['key'], o['n'], o['ng']) for o in bign]}")

stats = {"N": N, "structure": len(multi), "no_structure": len(single),
         "ngroups_dist": dict(gd), "depth_dist": dict(dd), "atype_dist": dict(at),
         "edge_unstable": ["chr20_42981498"], "edge_within_germline": [o["key"] for o in wg],
         "edge_deep": [o["key"] for o in deep], "edge_lowcpg": [o["key"] for o in lowc], "edge_bign": [o["key"] for o in bign]}
json.dump(stats, open(f"{A}/phylo_pilot_stats.json", "w"), ensure_ascii=False, indent=1)
print(f"\nDONE → phylo_pilot_stats.json")
