#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
完整量化普查(2026-07-10,使用者階層):S → W(k=1/k>1) → C 群組合 → HP 分類 → 拓撲確認 → 拓撲類型。
總數 + 相對比例。資料 = 最新 layered(per-HP-家族)+ mlhp(region-level C)+ completeness ledger(S)。
只讀不寫源碼。partial flag: HCC1395。用法: python3 full_census_hierarchy.py
"""
import json, os, glob
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.normpath(os.path.join(HERE, "..", "data"))
PILOT = os.path.normpath(os.path.join(HERE, "..", "..", "20260618_subcluster_pilot"))

led = json.load(open(f"{D}/sm_completeness_ledger.json"))
lay = json.load(open(f"{PILOT}/layered_reconstruction_HCC1395.json"))
detail = lay["detail"]
# mlhp region-level(n_populations_with_ALT = C)
mreg = {}
for f in sorted(glob.glob(f"{PILOT}/mlhp_part_*.json")):
    for g in json.load(open(f)).get("groups", []):
        mreg[f"{g['chrom']}:{g['start']}-{g['end']}"] = g

def pct(n, d): return round(n / d * 100, 1) if d else 0.0

out = {"partial_flag": "HCC1395 / 最新 layered per-HP-家族框架"}

# ── A. sSNV 總數 S + k=1/k>1(sSNV 層 bucket)──
S = led["universe_total"]; b = led["buckets"]
out["A_sSNV"] = {
    "S_total_sSNV": S, "TP": led["universe_tp"], "FP": led["universe_fp"],
    "isolated_singleton(k=1·單點·無樹)": b["isolated_singleton"], "iso_pct": pct(b["isolated_singleton"], S),
    "underpowered(有partner無足夠共讀)": b["underpowered"], "under_pct": pct(b["underpowered"], S),
    "linked(有共現·可建樹)": b["linked"], "linked_pct": pct(b["linked"], S),
}

# ── B. 區域 W:k=1 vs k>1(用 layered 單位的 region 集 + n_sSNV)──
region_fams = defaultdict(set); region_nss = {}
for u in detail:
    region_fams[u["region"]].add(u["family"])
    region_nss[u["region"]] = u.get("n_sSNV", 0)
W = len(region_fams)
k1 = sum(1 for r, n in region_nss.items() if n == 1)
k2 = sum(1 for r, n in region_nss.items() if n >= 2)
out["B_regions"] = {"W_total_regions": W, "k=1(單點)": k1, "k1_pct": pct(k1, W),
                    "k>1(兩點以上·可建樹)": k2, "k2_pct": pct(k2, W)}

# ── C. C 群組合(每 lineage 單位的 full-pop 群數)× 確認樹結構 ──
#   C = n_full_pops(該家族單位的全覆蓋基因型群數);確認=determined、不確認=ambiguous/capped/recurrence
lineage = [u for u in detail if u.get("is_lineage")]
cbucket = defaultdict(lambda: Counter())
for u in lineage:
    C = u.get("n_full_pops", 0)
    ck = str(C) if C <= 6 else "7+"
    cls = u["L1_class"]
    kind = "確認樹(determined)" if cls == "determined" else ("多可能性(ambiguous)" if "ambiguous" in cls else ("capped" if "capped" in cls else "recurrence"))
    cbucket[ck][kind] += 1
out["C_groups"] = {ck: dict(v) for ck, v in sorted(cbucket.items(), key=lambda x: (len(x[0]), x[0]))}
out["C_groups_note"] = "C = 每 lineage 單位的全覆蓋基因型群數(n_full_pops);0=純 partial 無 full-pop"

# ── D. HP 分類(region 層):HP1 xor/and HP2 × HP3 ──
d = {"HP1_xor_HP2": {"not_HP3": 0, "HP3": 0}, "HP1_and_HP2": {"not_HP3": 0, "HP3": 0},
     "no_germline(only 3/none)": 0}
for r, fams in region_fams.items():
    germ = fams & {"1", "2"}; has3 = "3" in fams
    if len(germ) == 1:
        d["HP1_xor_HP2"]["HP3" if has3 else "not_HP3"] += 1
    elif len(germ) == 2:
        d["HP1_and_HP2"]["HP3" if has3 else "not_HP3"] += 1
    else:
        d["no_germline(only 3/none)"] += 1
xor = sum(d["HP1_xor_HP2"].values()); andd = sum(d["HP1_and_HP2"].values())
out["D_HP"] = {**d, "HP1_xor_HP2_total": xor, "xor_pct": pct(xor, W),
               "HP1_and_HP2_total": andd, "and_pct": pct(andd, W)}

# ── E. 拓撲確認狀況(lineage 單位層 determinacy)+ 可能性數量(n_trees)──
nl = len(lineage)
det = Counter(u["L1_class"] for u in lineage)
out["E_determinacy"] = {"n_lineage_units": nl,
                        "分類": {k: {"n": v, "pct": pct(v, nl)} for k, v in det.most_common()}}
# 可能性數量(n_trees)分佈 — 只對 ambiguous(多候選)
amb = [u for u in lineage if "ambiguous" in u["L1_class"]]
ntree_dist = Counter(min(u.get("n_trees", 0), 10) if u.get("n_trees", 0) <= 10 else "10+" for u in amb)
out["E_possibility_count"] = {"determined(唯一樹=1可能)": det["determined"],
                              "ambiguous 的 n_trees 分佈": {str(k): v for k, v in sorted(ntree_dist.items(), key=lambda x: str(x[0]))},
                              "note": "n_trees = 該單位相容的候選樹數(可能性數量)"}

# ── F. 拓撲類型(用 n_distinct_shapes / n_hidden 概觀)──
shape = Counter()
for u in lineage:
    ds = u.get("n_distinct_shapes")
    nh = u.get("n_hidden", 0)
    if u["L1_class"] == "determined":
        shape["determined·單一形狀" + ("(含隱藏祖先)" if nh > 0 else "(純觀測)")] += 1
    elif "ambiguous" in u["L1_class"]:
        shape[f"ambiguous·{ds if ds else '?'}種形狀"] += 1
out["F_topology_shapes"] = dict(shape.most_common(12))

# ── G. 完整樹? ──
out["G_complete_tree"] = {
    "能否確認樣本完整癌症樹": "否",
    "理由": f"單一 bulk:{pct(det['ambiguous_structure(多完成/多結構)'], nl)}% 單位多候選、跨區整樹非可辨識(資訊論極限);需 single-cell/multi-region",
    "封頂": "⭐3(characterize 非 confirm)"}

json.dump({"summary": out}, open(f"{D}/full_census_hierarchy.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(out, ensure_ascii=False, indent=1))
