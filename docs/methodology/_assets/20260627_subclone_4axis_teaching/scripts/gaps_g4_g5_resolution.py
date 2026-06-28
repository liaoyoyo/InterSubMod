#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G4 isolated 計數收斂 + G5 覆蓋功率計算(純 JSON/census,不需 BAM)。
G4: 四套「isolated/single」數字測不同東西 → 定 canonical + 各自定義。
G5: underpowered(有 partner 無 co-read link)→ 用 max_coread × 覆蓋倍數估 2x/4x 可救回數。
輸出 gaps_g4_g5_resolution.json。compute batch(§13.0)。
"""
import json, os, csv
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

# ---- G4: isolated 計數收斂 ----
led = json.load(open(f"{DATA}/sm_completeness_ledger.json"))
ssa = json.load(open(f"{DATA}/single_snv_accounting.json"))
lms = json.load(open(f"{DATA}/sm_locus_master_summary.json"))
g4 = {
    "definitions": {
        "isolated_singleton_8320": {"value": led["buckets"]["isolated_singleton"],
            "basis": "per-sSNV;n_partners_le50k==0(read-span≤50kb 內無任何 partner)",
            "denominator": led["universe_total"], "role": "🔵 CANONICAL『isolated』(骨幹連鎖本義:無法被 read 連)"},
        "tree_role_all_isolated_11520": {"value": lms.get("tree_role_all", {}).get("isolated"),
            "basis": "全 loci 的 tree_role 分類(含 germline/未過濾,基底更廣)", "role": "不同分母,非 isolation 本義"},
        "tree_role_somatic_isolated_266": {"value": lms.get("tree_role_somatic", {}).get("isolated"),
            "basis": "只 somatic-confirmed loci 的 tree_role isolated", "role": "somatic 子集 tree-role"},
        "single_sSNV_region_13778": {"value": ssa["single_sSNV_total"],
            "basis": "區層級:只含 1 個 sSNV 的『區』(非 per-sSNV isolation)", "denominator": ssa["universe_total"],
            "role": "區計數,與 isolation 正交"},
    },
    "canonical": "isolated = 8,320 sSNV(n_partners_le50k==0);TP|isolated 7555 / FP|isolated 765",
    "why_not_conflict": "四數測不同東西:8320=per-sSNV 無 partner(本義);11520=tree_role 全 loci(廣基底);266=somatic tree_role;13778=單-sSNV 區。非矛盾,是混用未標基底。",
    "presentation_gap": "isolated 8320 loci 無 tree(不在 topology detail 3885)→ 最終輸出需獨立『isolated 表』:每 sSNV 有 caller VAF(可放 clonal 譜)+ 可能 same-PS partner(Tier-PS 放單倍型,非巢狀)。",
}

# ---- G5: 覆蓋功率(underpowered 救回估計)----
# per_sSNV_census: underpowered = n_partners>0 AND n_realized_links==0(有 partner 但 coread<6 無 link)
MIN_COREAD = 6
rows = []
with open(f"{DATA}/per_sSNV_census.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try:
            npart = int(r["n_partners_le50k"]); nlink = int(r["n_realized_links"]); mc = int(r["max_coread"])
        except Exception:
            continue
        rows.append((npart, nlink, mc))
underpowered = [(npart, nlink, mc) for npart, nlink, mc in rows if npart > 0 and nlink == 0]
# 估:coread 隨覆蓋線性 → m 倍後 max_coread*m;達 MIN_COREAD 即可建至少 1 link
def rescued_at(m):
    return sum(1 for npart, nlink, mc in underpowered if mc * m >= MIN_COREAD)
n_up = len(underpowered)
g5 = {
    "underpowered_definition": "n_partners_le50k>0 AND n_realized_links==0(有 partner 但 max_coread<6 無 link)",
    "n_underpowered": n_up,
    "min_coread_for_link": MIN_COREAD,
    "current_max_coread_dist": dict(sorted(Counter(mc for _, _, mc in underpowered).items())[:8]),
    "rescued_estimate": {
        "at_2x": rescued_at(2), "at_3x": rescued_at(3), "at_4x": rescued_at(4),
        "at_2x_pct": round(100 * rescued_at(2) / n_up, 1) if n_up else 0,
        "at_4x_pct": round(100 * rescued_at(4) / n_up, 1) if n_up else 0},
    "model_note": "線性近似:coread ∝ 覆蓋;m 倍後 max_coread×m≥6 即可建 ≥1 link(救回為樹可用)。"
                  "保守(忽略 read 長/取樣分散);僅給數量級,非精確功率。C_underdetermined(550 單群)亦受此限。",
    "caveat": "此為『可建 link』的上界估計;實際還需該 link 落在能定樹的位置 + 非 CN-gain multiplicity 假象。",
}

out = {"G4_isolated_reconciliation": g4, "G5_coverage_power": g5}
json.dump(out, open(f"{DATA}/gaps_g4_g5_resolution.json", "w"), ensure_ascii=False, indent=1)
print("G4/G5 DONE")
print("G4 canonical:", g4["canonical"])
print(f"G5 underpowered={n_up};救回估 2x={g5['rescued_estimate']['at_2x']}({g5['rescued_estimate']['at_2x_pct']}%) 4x={g5['rescued_estimate']['at_4x']}({g5['rescued_estimate']['at_4x_pct']}%)")
print("G5 max_coread 分布(前8):", g5["current_max_coread_dist"])
