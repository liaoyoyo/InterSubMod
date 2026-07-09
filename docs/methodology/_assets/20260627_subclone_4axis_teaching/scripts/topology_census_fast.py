#!/usr/bin/env python3
"""
topology_census_fast.py — 拓撲普查(讀 region-view,不重跑 solver;2026-07-09)。
全樣本有哪些樹拓撲、各佔多少、生物意義;HCC1395 5khz vs DORADO 同區一致性。
🔴 多數 ambiguous 區「形狀」確定(n_distinct_shapes=1)→ 拓撲可確認即使精確樹欠定。
拓撲→生物: single=單clone / linear(A⊂B⊂C)=巢狀subclone後代 / branched_internal=姊妹subclone分歧 / star_root=多獨立lineage。
用法: python3 topology_census_fast.py
"""
import os, json, glob
from collections import Counter, defaultdict

MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]

def rv_path(s):
    return f"{PILOT}/layered_region_view_HCC1395.json" if s == "HCC1395" else f"{MSROOT}/{s}/layered_region_view_{s}.json"

def tree_shape(edges):
    """由樹邊分類形狀。回 shape。"""
    if not edges:
        return "single_or_germline"
    ch = defaultdict(list); nodes = set()
    for p, c in edges:
        ch[p].append(c); nodes.add(p); nodes.add(c)
    if len([n for n in nodes if n != "ROOT"]) <= 1:
        return "single"
    childcnt = {p: len(cs) for p, cs in ch.items()}
    if any(v >= 2 for p, v in childcnt.items() if p != "ROOT"):
        return "branched_internal"   # 姊妹subclone(clone內分歧)
    if childcnt.get("ROOT", 0) >= 2:
        return "star_root"           # germline→多平行lineage
    return "linear"                  # 巢狀後代 A⊂B⊂C

BIO = {"single": "單clone(1突變集)", "single_or_germline": "單clone/germline-only",
       "linear": "巢狀subclone後代(序列累積)", "branched_internal": "姊妹subclone(clone內分歧)",
       "star_root": "多獨立lineage(germline→平行)", "capped": "太密(未定)", "mixed_shape": "形狀未定(多形)"}

def unit_shape(L):
    """回 (shape, determined_shape)。ambiguous 但單一形狀 → 形狀確定。"""
    if L.get("capped"):
        return "capped", False
    trees = L.get("trees") or []
    if not trees:
        return "single_or_germline", True   # determined 但 germline-only/單群(edges=[])
    shapes = {tree_shape(t["edges"]) for t in trees}
    if len(shapes) == 1:
        return next(iter(shapes)), True
    return "mixed_shape", False

def run(s):
    d = json.load(open(rv_path(s)))["regions"]
    sc = Counter(); det = Counter(); region_shape = {}
    for r in d:
        for L in r["lineages"]:
            if L["family"] not in ("1", "2"):
                continue
            sh, dd = unit_shape(L)
            sc[sh] += 1
            det["形狀確定" if dd else "形狀未定"] += 1
            if sh not in ("capped", "mixed_shape"):
                region_shape[(r["chrom"], r["start"], L["family"])] = sh
    return sc, det, region_shape

print("=" * 92)
print("拓撲普查(新骨幹;region-view):全樣本有哪些樹拓撲、各佔多少、生物意義")
print("拓撲→生物: single=單clone / linear=巢狀subclone後代 / branched_internal=姊妹subclone分歧 / star_root=多獨立lineage")
print("=" * 92)
allrs = {}
ORDER = ["single", "single_or_germline", "linear", "branched_internal", "star_root", "mixed_shape", "capped"]
for s in SAMPLES:
    sc, det, rs = run(s); allrs[s] = rs
    tot = sum(sc.values())
    # 生物聚合:有subclone結構(linear+branched+star) vs 單clone
    substruct = sc["linear"] + sc["branched_internal"] + sc["star_root"]
    sister = sc["branched_internal"] + sc["star_root"]
    print(f"\n### {s}  (germline lineage={tot}；形狀確定 {det.get('形狀確定',0)}={100*det.get('形狀確定',0)/tot:.0f}%)")
    for k in ORDER:
        if sc.get(k):
            print(f"    {sc[k]:6d} ({100*sc[k]/tot:4.1f}%)  {k}  — {BIO.get(k,'')}")
    print(f"    生物: 有subclone結構(linear+branched) {substruct}({100*substruct/tot:.0f}%) | 其中姊妹分支(branched+star) {sister}({100*sister/tot:.0f}%) | linear巢狀 {sc['linear']}({100*sc['linear']/tot:.0f}%)")

print("\n" + "=" * 92)
print("HCC1395(5khz) vs DORADO 同細胞株同區拓撲一致性")
print("=" * 92)
r5, rd = allrs["HCC1395"], allrs["HCC1395_DORADO"]
common = set(r5) & set(rd)
same = sum(1 for k in common if r5[k] == rd[k])
print(f"  共同 (chrom,start,fam): {len(common)};拓撲一致 {same} ({100*same/len(common):.1f}%)" if common else "無共同")
# 生物層一致(subclone有無 = single vs 有結構)
def biocat(sh): return "單clone" if sh in ("single", "single_or_germline") else ("姊妹" if sh in ("branched_internal", "star_root") else "巢狀")
biosame = sum(1 for k in common if biocat(r5[k]) == biocat(rd[k]))
print(f"  生物類別一致(單clone/巢狀/姊妹) {biosame} ({100*biosame/len(common):.1f}%)" if common else "")
mism = Counter((r5[k], rd[k]) for k in common if r5[k] != rd[k])
if mism:
    print("  不一致(top):")
    for (a, b), v in mism.most_common(5):
        print(f"    {v:4d}  5khz={a} / DORADO={b}")
