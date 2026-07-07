#!/usr/bin/env python3
"""
build_region_view.py — layered_reconstruction_{S}.json → layered_region_view_{S}.json (2026-07-07)
region-centric「最終數據格式」(資料模型 spec §5;主分母=region)。含 census(HP-multiplicity + region-determinacy
+ U1 somatic sSNV + U5 樹 + U6 隱藏祖先)+ per-region 逐 lineage(枚舉樹 + n_distinct_shapes + 判斷軌跡)。
env: SM_LAYERED=<layered_reconstruction.json> SM_OUT=<region_view.json>
     SM_CENSUS=<sm_linkage_genomewide.json>(取 somatic sSNV total;可省) SM_INTEGRATION=<sm_region_integration.json>(取 linkage 全 span;可省) SM_SAMPLE=<name>
"""
import os, sys, json
from collections import Counter, defaultdict

LAY = os.environ.get("SM_LAYERED")
OUT = os.environ.get("SM_OUT")
SAMPLE = os.environ.get("SM_SAMPLE", "?")
if not LAY or not OUT:
    sys.exit("需 SM_LAYERED + SM_OUT")
L = json.load(open(LAY, encoding="utf-8"))
detail = L["detail"]

def region_det(units):
    gl = [u for u in units if u["family"] in ("1", "2")]
    if not gl:
        return "no_germline_lineage"
    def sh(c):
        return 'D' if c == 'determined' else ('R' if c.startswith('recurrence') else ('C' if c.startswith('capped') else 'A'))
    s = [sh(u["L1_class"]) for u in gl]
    if all(x == 'D' for x in s): return "all_determined"
    if any(x == 'R' for x in s): return "has_recurrence"
    if any(x == 'C' for x in s): return "has_capped"
    return "has_ambiguous"

byreg = defaultdict(list)
for u in detail:
    byreg[u["region"]].append(u)

regions = []
for rid, units in byreg.items():
    u0 = units[0]
    gl = [u for u in units if u["family"] in ("1", "2")]
    regions.append({
        "region": rid, "chrom": u0["chrom"], "start": u0["start"], "end": u0["end"],
        "n_sSNV": u0["n_sSNV"], "cn": u0["cn"],
        "hp_multiplicity": len(gl), "is_multiHP": len(gl) >= 2,
        "region_determinacy": region_det(units),
        "lineages": sorted([{
            "family": u["family"], "fam_label": u["fam_label"], "is_germline": u["family"] in ("1", "2"),
            "L1_class": u["L1_class"], "n_trees": u["n_trees"], "n_trees_stored": u.get("n_trees_stored"),
            "n_distinct_shapes": u.get("n_distinct_shapes"), "n_hidden": u["n_hidden"], "capped": u["capped"],
            "n_reads": u["n_reads"], "n_full_pops": u["n_full_pops"], "n_partial": u["n_partial"],
            "trees": u["trees"], "L2_cn_verdict": u["L2_cn_verdict"], "verify_pass": u["verify_pass"], "trace": u["trace"],
        } for u in units], key=lambda x: (x["family"] == "none", x["family"] == "3", x["family"])),
    })
regions.sort(key=lambda r: (-r["hp_multiplicity"], -r["n_sSNV"]))

# census
lin = [u for u in detail if u["is_lineage"]]
noncap = [u for u in lin if not u["capped"]]
sum_trees = sum(u["n_trees"] for u in noncap)
ntdist = Counter()
for u in lin:
    n = u["n_trees"]
    ntdist['1' if n == 1 else ('2-5' if n <= 5 else ('6-20' if n <= 20 else '>20'))] += 1
nh = Counter(u["n_hidden"] for u in lin)
census = {"n_regions": len(regions),
          "hp_multiplicity": dict(Counter(r["hp_multiplicity"] for r in regions)),
          "region_determinacy": dict(Counter(r["region_determinacy"] for r in regions)),
          "L0": L["L0_hp_family"], "L1": L["L1_ssnv_algorithm"], "L2": L["L2_cn"],
          "U5_trees": {"sum_ntrees_noncapped": sum_trees, "ntrees_dist": dict(ntdist)},
          "U6_hidden": {"sum_hidden": sum(u["n_hidden"] for u in lin), "dist": dict(sorted(nh.items()))}}
# U1 somatic sSNV total(可選)
cen_path = os.environ.get("SM_CENSUS")
if cen_path and os.path.exists(cen_path):
    cen = json.load(open(cen_path, encoding="utf-8")).get("census", {})
    census["U1_sSNV_somatic_total"] = sum(1 for v in cen.values() if v.get("somatic") is True)
    census["U1_census_total_positions"] = len(cen)
intg = os.environ.get("SM_INTEGRATION")
if intg and os.path.exists(intg):
    census["U3_linkage_regions_full_span"] = json.load(open(intg, encoding="utf-8"))["aggregate"]["n_regions"]

out = {"sample": SAMPLE, "data_model_spec": "docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md",
       "census": census, "regions": regions}
json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False)
print(f"OK region-view [{SAMPLE}] n_regions={len(regions)} region_det={dict(Counter(r['region_determinacy'] for r in regions))} -> {OUT}")
