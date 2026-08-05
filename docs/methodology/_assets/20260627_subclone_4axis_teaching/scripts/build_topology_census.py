#!/usr/bin/env python3
"""
build_topology_census.py — 拓撲普查 + 生物意義 + HCC1395/DORADO 一致性(2026-07-09 使用者要求)。
問題:全樣本有哪些拓撲、各佔多少;哪些=subclone(巢狀後代)、哪些=姊妹clone(分支);同細胞株同區是否一致。
🔴 關鍵:多數 ambiguous 區的「形狀」確定(n_distinct_shapes=1,樹僅 label/順序變體)→ 拓撲類型可確認,即使精確樹欠定。
拓撲→生物:single=單clone / linear(A⊂B⊂C)=巢狀subclone後代 / branched-internal=姊妹subclone分歧 / branched-root=多獨立lineage。
用法: python3 build_topology_census.py
"""
import os, sys, json, glob
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tree_enumeration_solver as S

MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]

def mlhp_paths(s):
    wd = PILOT if s == "HCC1395" else f"{MSROOT}/{s}"
    return sorted(glob.glob(f"{wd}/mlhp_part_*.json"))

def tree_shape(edges):
    """由樹邊分類形狀 + 生物意義。回 (shape, bio)。"""
    if not edges:
        return "single", "單clone/germline-only"
    ch = defaultdict(list)
    nodes = set()
    for p, c in edges:
        ch[p].append(c); nodes.add(p); nodes.add(c)
    n_nonroot = len([n for n in nodes if n != "ROOT"])
    if n_nonroot <= 1:
        return "single", "單clone(1突變集)"
    childcnt = {p: len(cs) for p, cs in ch.items()}
    internal_branch = any(v >= 2 for p, v in childcnt.items() if p != "ROOT")
    root_branch = childcnt.get("ROOT", 0) >= 2
    if internal_branch:
        return "branched_internal", "姊妹subclone(clone 內分歧→平行 subclone)"
    if root_branch:
        return "star_root", "多獨立lineage(germline→多平行,含姊妹clone)"
    return "linear", "巢狀subclone後代(A⊂B⊂C 序列累積)"

def unit_shape(full, part, k):
    """回 (shape, bio, determined_shape:bool, n_distinct, base_class)。ambiguous 但單一形狀 → 形狀確定。"""
    res = S.enumerate_min_trees(full, part, k)
    base = S.classify(res)
    if res["capped"]:
        return "capped", "太密(枚舉未完)", False, None, base
    if not res["trees"]:
        return "underdetermined", "無樹", False, 0, base
    shapes = set()
    for t in res["trees"]:
        shapes.add(tree_shape(t["edges"])[0])
    if len(shapes) == 1:
        sh, bio = tree_shape(res["trees"][0]["edges"])
        return sh, bio, True, 1, base            # 形狀確定(即使精確樹 ambiguous)
    return "mixed_shape", f"形狀未定({'/'.join(sorted(shapes))})", False, len(shapes), base

def run_sample(s):
    groups = []
    for f in mlhp_paths(s):
        groups += json.load(open(f))["groups"]
    shape_c = Counter(); bio_c = Counter(); shape_det = Counter()
    region_shape = {}   # (chrom,start,fam) -> shape (供跨樣本比對)
    for g in groups:
        if g.get("n_sSNV", 0) < 2:
            continue
        pbh = g.get("populations_by_hp", {}) or {}; sbh = g.get("subread_groups_by_hp", {}) or {}
        for fam in ("1", "2"):
            full = pbh.get(fam, {}) or {}; part = list((sbh.get(fam, {}) or {}).keys())
            if not full and not part:
                continue
            k = len(next(iter(full))) if full else len(part[0])
            sh, bio, det, nd, base = unit_shape(full, part, k)
            shape_c[sh] += 1
            if sh not in ("capped", "underdetermined", "mixed_shape"):
                bio_c[bio] += 1
                shape_det["形狀確定" if det else "形狀未定"] += 1
                region_shape[(g["chrom"], g["start"], fam)] = sh
            else:
                shape_det["形狀未定"] += 1
    return shape_c, bio_c, shape_det, region_shape

# === 跑全 7 樣本 ===
print("=" * 100)
print("拓撲普查:全樣本有哪些樹拓撲、各佔多少、生物意義")
print("=" * 100)
all_rs = {}
for s in SAMPLES:
    sc, bc, sd, rs = run_sample(s); all_rs[s] = rs
    tot = sum(sc.values())
    print(f"\n### {s}  (germline lineage units={tot}；形狀確定 {sd.get('形狀確定',0)} / 未定 {sd.get('形狀未定',0)})")
    for k, v in sc.most_common():
        print(f"    {v:6d} ({100*v/tot:4.1f}%)  {k}")

# === HCC1395 5khz vs DORADO 同區一致性 ===
print("\n" + "=" * 100)
print("HCC1395(5khz) vs DORADO 同細胞株同區拓撲一致性")
print("=" * 100)
r5 = all_rs["HCC1395"]; rd = all_rs["HCC1395_DORADO"]
common = set(r5) & set(rd)
same = sum(1 for k in common if r5[k] == rd[k])
print(f"  兩者共同 (chrom,start,fam) 單位: {len(common)}")
print(f"  拓撲一致: {same} ({100*same/len(common):.1f}%)" if common else "  無共同")
mism = Counter()
for k in common:
    if r5[k] != rd[k]:
        mism[(r5[k], rd[k])] += 1
if mism:
    print("  不一致組合(top):")
    for (a, b), v in mism.most_common(5):
        print(f"    {v:4d}  5khz={a} / DORADO={b}")
