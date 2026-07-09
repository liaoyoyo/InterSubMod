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

# ── A. sSNV 總數 S(=ClairS PASS = somatic_pass.vcf.gz;07-09 layered 實際 universe)──
#   🔴 2026-07-10 修:07-09 pipeline 骨幹=ClairS PASS(somatic_pass),移除舊 is_somatic 檢+SEQC2 TP/FP 子集。
#   S=113997 由 zcat somatic_pass.vcf.gz PASS SNV 實測(canonical HCC1395 longphase_s)。
#   舊 ledger 的 35332(TP 30490/FP 4842)= SEQC2-可評估子集(僅 31%),只作觀察標籤。
S = 113997
SEQC2_eval = led["universe_total"]  # 35332 = SEQC2-可評估子集(TP+FP)
# 多位點 sSNV = mlhp 區的 n_sSNV 和(區不重疊);isolated/單點 = S − 多位點
multi_ssnv = sum(g.get("n_sSNV", 0) for g in mreg.values())
iso = S - multi_ssnv
out["A_sSNV"] = {
    "S_total_ClairS_PASS": S, "note": "= somatic_pass.vcf.gz PASS SNV(LongPhase-S/ClairS);移除舊 is_somatic 檢",
    "SEQC2_evaluable_subset(TP+FP·僅觀察)": SEQC2_eval, "SEQC2_TP": led["universe_tp"], "SEQC2_FP": led["universe_fp"],
    "SEQC2_pct_of_S": pct(SEQC2_eval, S),
    "k=1_isolated(單點·無樹)": iso, "iso_pct": pct(iso, S),
    "k>1_multilocus_sSNV(可建樹)": multi_ssnv, "multi_pct": pct(multi_ssnv, S),
}

# ── B. 區域 W:k=1 vs k>1 ──
region_fams = defaultdict(set)
for u in detail:
    region_fams[u["region"]].add(u["family"])
W = len(region_fams)
out["B_regions"] = {"S_ClairS_PASS": S, "k=1_單點_sSNV(無樹)": iso,
                    "k>1_多位點_區數": len(mreg), "k>1_多位點_sSNV": multi_ssnv,
                    "regions_with_lineage_units": W, "note": "k>1 多位點區→依 HP 家族拆 lineage 單位(下方 C-F 母體)"}

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
