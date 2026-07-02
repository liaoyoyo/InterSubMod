#!/usr/bin/env python3
"""[邊界案例] 全基因組 v3 records 的邊界案例盤點 — 供確認是否合理。
0群退化 / 高群數(≥5) / 深層階層(depth≥4) / FP帶結構(unaligned) / 低CpG結構 / 大n。
列各類 n·C·aligned 分佈 + 抽代表位點(供 mini-VCF 渲染)。純讀 records(L1)。"""
import json
import numpy as np
from collections import Counter
A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
rec = json.load(open(f"{A}/phylo_v3_wg_records.json"))
TP = [r for r in rec if r["set"] == "TP"]; FP = [r for r in rec if r["set"] == "FP"]


def desc(g, name):
    if not g: print(f"  {name}: 0"); return
    ns = [r["n"] for r in g]; cs = [r["C"] for r in g]
    print(f"  {name}: {len(g)} 個 | n中位={int(np.median(ns))}(範圍{min(ns)}-{max(ns)}) | C中位={int(np.median(cs))}(範圍{min(cs)}-{max(cs)})")


print(f"=== 全基因組 v3: TP {len(TP)} / FP {len(FP)} ===\n")

print("【E1 0群退化】(ngroups=0 = 全部 read 變離群)")
e1tp = [r for r in TP if r["ngroups"] == 0]; e1fp = [r for r in FP if r["ngroups"] == 0]
desc(e1tp, "TP"); desc(e1fp, "FP")
print(f"  → n_outlier 中位={int(np.median([r['n_outlier'] for r in e1tp])) if e1tp else 0}; 多為小 n 全離群? n≤9 佔 {sum(1 for r in e1tp if r['n']<=9)}/{len(e1tp)}")

print("\n【E2 高群數 ≥5】")
e2tp = [r for r in TP if r["ngroups"] >= 5]; e2fp = [r for r in FP if r["ngroups"] >= 5]
desc(e2tp, "TP"); desc(e2fp, "FP")
print(f"  TP 群數分佈: {dict(sorted(Counter(r['ngroups'] for r in e2tp).items()))}")
print(f"  aligned 比例: TP {sum(1 for r in e2tp if r['aligned'])}/{len(e2tp)} | FP {sum(1 for r in e2fp if r['aligned'])}/{len(e2fp)}")

print("\n【E3 深層階層 depth≥4】")
e3tp = [r for r in TP if r["maxdepth"] >= 4]; e3fp = [r for r in FP if r["maxdepth"] >= 4]
desc(e3tp, "TP"); desc(e3fp, "FP")
print(f"  TP depth 分佈: {dict(sorted(Counter(r['maxdepth'] for r in e3tp).items()))}")

print("\n【E4 FP 帶結構 unaligned】(FP 不該有真 somatic 結構 → 應為噪音/cis 假象)")
e4 = [r for r in FP if r["ngroups"] >= 2 and not r["aligned"]]
desc(e4, "FP unaligned multi")
print(f"  佔 FP {100*len(e4)/len(FP):.1f}% | n中位={int(np.median([r['n'] for r in e4])) if e4 else 0} | V_allele 中位={np.median([r['V_allele'] for r in e4]) if e4 else 0:.2f}")
e4tp = [r for r in TP if r["ngroups"] >= 2 and not r["aligned"]]
print(f"  對照 TP unaligned multi: {len(e4tp)} ({100*len(e4tp)/len(TP):.1f}%) → FP/TP 比 = {len(e4)/len(FP)/(len(e4tp)/len(TP)):.2f}× (>1=FP偏多=非subclone)")

print("\n【E5 低 CpG(≤20) 帶結構】(校準: C≤20 FP 升 16-24%)")
e5 = [r for r in TP if r["C"] <= 20 and r["ngroups"] >= 2]
allc = [r for r in TP if r["C"] <= 20]
print(f"  C≤20 的 TP: {len(allc)} 個, 其中帶結構 {len(e5)} ({100*len(e5)/len(allc) if allc else 0:.1f}%) → 這些 structure 可信度低(FP風險區)")
print(f"  C≤40 帶結構: {sum(1 for r in TP if r['C']<=40 and r['ngroups']>=2)} / C≤40 共 {sum(1 for r in TP if r['C']<=40)}")

print("\n【E6 大 n(≥150) 高信心】")
e6 = [r for r in TP if r["n"] >= 150 and r["ngroups"] >= 2]
desc(e6, "TP 大n多群")
print(f"  aligned: {sum(1 for r in e6 if r['aligned'])}/{len(e6)}")

# 抽代表位點供渲染(每類 3 個, 有 pos/chrom)
def pick(g, k=3):
    g2 = [r for r in g if r.get("chrom") and r.get("pos")]
    g2.sort(key=lambda r: -r["n"]); return [(r["chrom"], r["pos"], r["set"], r["n"], r["C"], r["ngroups"], r["maxdepth"], r["aligned"]) for r in g2[:k]]


samples = {"E1_degenerate": pick(e1tp), "E2_highk": pick(e2tp), "E3_deep": pick(e3tp),
           "E4_fp_unaligned": pick(e4), "E5_lowcpg": pick(e5), "E6_bign": pick(e6)}
json.dump(samples, open(f"{A}/phylo_edge_samples.json", "w"), ensure_ascii=False, indent=1)
print("\n=== 代表位點(供 mini-VCF 渲染) → phylo_edge_samples.json ===")
for cat, s in samples.items():
    print(f"  {cat}: {[(c, p, f'n{n}', f'C{C}', f'{ng}群', f'd{d}', 'aln' if al else 'unaln') for c, p, st, n, C, ng, d, al in s]}")
