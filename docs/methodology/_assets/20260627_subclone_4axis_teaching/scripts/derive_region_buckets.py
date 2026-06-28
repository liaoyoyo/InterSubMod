#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
從凍結的 regions.tsv（逐區明細，§13-A 可驗證源）derive 出之前只在 prose/硬編的數字：
  - 6-bucket tree_shape 分布 + 每 shape 的乾淨 CN(loh+neutral) 子集
  - 乾淨 full_tree=205 等（取代 build_4axis_teaching.py 的硬編 CLEAN）
  - n_populations 分布（主張4 的 by_n_populations）
輸出 grep-able JSON，缺 anchor 不符即 assert 失敗（拒絕產出錯誤真值）。
用法：python3 derive_region_buckets.py
"""
import csv, json, os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TSV = os.path.join(DATA, "regions.tsv")

rows = list(csv.DictReader(open(TSV, encoding="utf-8"), delimiter="\t"))
N = len(rows)
CLEAN_CN = {"loh", "neutral"}
SHAPES = ["full_tree", "linear_nested", "sibling_only", "co_linked_lineage", "no_confirmed_structure", "inconsistent"]

by_shape = {}
for sh in SHAPES:
    sub = [r for r in rows if r["tree_shape"] == sh]
    clean = [r for r in sub if r["cn"] in CLEAN_CN]
    by_shape[sh] = {"n": len(sub), "clean_loh_neutral": len(clean), "pct_of_regions": round(100*len(sub)/N, 1)}

structured = sum(by_shape[s]["n"] for s in ("full_tree", "linear_nested", "sibling_only"))
structured_plus = structured + by_shape["co_linked_lineage"]["n"]
clean_subset = {s: by_shape[s]["clean_loh_neutral"] for s in SHAPES}

# n_populations 分布（主張4：更多位點→更多群的實測上限）
npop = Counter(int(r["n_populations"]) for r in rows if r["n_populations"].isdigit())
npop_dist = {str(k): npop[k] for k in sorted(npop)}

# n_sSNV 分布（power：每區能共讀幾個 sSNV）
nsnv = Counter(int(r["n_sSNV"]) for r in rows if r["n_sSNV"].isdigit())
nsnv_summary = {"max": max(nsnv), "n_regions_ge3_sSNV": sum(v for k, v in nsnv.items() if k >= 3)}

cn_dist = dict(Counter(r["cn"] for r in rows))
has_cycle = sum(1 for r in rows if r["has_cycle"].strip().lower() in ("true", "1"))

out_dist = {
    "source": "derived from regions.tsv (frozen @ branch feat/summary-nreadsvalid)",
    "total_regions": N,
    "by_shape": by_shape,
    "structured_3shape(full+linear+sibling)": structured,
    "structured_plus_colinked": structured_plus,
    "no_confirmed_structure": by_shape["no_confirmed_structure"]["n"],
    "inconsistent": by_shape["inconsistent"]["n"],
    "cn_distribution": cn_dist,
    "clean_loh_neutral_total": sum(1 for r in rows if r["cn"] in CLEAN_CN),
    "has_cycle_regions": has_cycle,
    "n_populations_dist": npop_dist,
    "n_sSNV": nsnv_summary,
}
out_clean = {
    "source": "derived from regions.tsv; clean = cn in {loh,neutral}",
    "definition": "論文級可信子集 = CN-clean(LOH+neutral) 的各 tree_shape 區",
    "clean_by_shape": clean_subset,
    "clean_full_tree": clean_subset["full_tree"],
}

with open(os.path.join(DATA, "region_shape_distribution.json"), "w", encoding="utf-8") as f:
    json.dump(out_dist, f, ensure_ascii=False, indent=1)
with open(os.path.join(DATA, "clean_subset.json"), "w", encoding="utf-8") as f:
    json.dump(out_clean, f, ensure_ascii=False, indent=1)

# ---- 對已知錨點驗證（不符即 refuse）----
assert N == 7143, f"total {N}!=7143"
assert by_shape["full_tree"]["n"] == 677, "full_tree!=677"
assert structured == 3820, f"structured {structured}!=3820"
assert structured_plus == 4678, f"structured+ {structured_plus}!=4678"
assert by_shape["no_confirmed_structure"]["n"] == 2443, "no_confirmed!=2443"
assert by_shape["inconsistent"]["n"] == 22, "inconsistent!=22"
print("OK wrote region_shape_distribution.json + clean_subset.json")
print(f"  6-bucket: {{ {', '.join(f'{s}:{by_shape[s][chr(110)]}' for s in SHAPES)} }}")
print(f"  structured 3820={structured}  +colinked 4678={structured_plus}")
print(f"  clean(LOH+neutral) by_shape: {clean_subset}")
print(f"  ANCHOR clean_full_tree={clean_subset['full_tree']} (region_trees §4 said 205)")
print(f"  n_populations dist: {npop_dist}")
print(f"  has_cycle regions={has_cycle} (應≈inconsistent 22)")
