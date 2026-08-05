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

# 原始觀測(2026-07-10):join mlhp groups → 每 lineage 帶 positions/col_coverage/populations/subreads
# (region_view 原本只存重建後的樹,丟失原始 read 證據 → 補回供工作站顯示 read 矩陣/pairwise/locus/位點證據)
import glob as _glob
raw_by_region = {}
_ml_glob = os.environ.get("SM_ML_GLOB")
if _ml_glob:
    for _f in _glob.glob(_ml_glob):
        try:
            _m = json.load(open(_f, encoding="utf-8"))
        except Exception:
            continue
        for g in _m.get("groups", []):
            if g.get("n_sSNV", 0) < 2:
                continue
            rid = g.get("region_id") or "{}:{}-{}".format(g.get("chrom"), g.get("start"), g.get("end"))
            if rid in raw_by_region:
                raise RuntimeError(f"duplicate raw region identity: {rid}")
            raw_by_region[rid] = g

def region_det(units):
    primary = [u for u in units if u.get("is_primary_lineage", u.get("is_lineage", False))]
    if not primary:
        return "no_primary_lineage"
    def sh(c):
        return 'D' if c == 'determined' else ('R' if c.startswith('recurrence') else ('C' if c.startswith('capped') else 'A'))
    s = [sh(u["L1_class"]) for u in primary]
    if all(x == 'D' for x in s): return "all_determined"
    if any(x == 'R' for x in s): return "has_recurrence"
    if any(x == 'C' for x in s): return "has_capped"
    return "has_ambiguous"

byreg = defaultdict(list)
for u in detail:
    byreg[u["region"]].append(u)

regions = []
region_ids = list(byreg)
region_ids.extend(sorted(rid for rid in raw_by_region if rid not in byreg))
for rid in region_ids:
    units = byreg.get(rid, [])
    _raw = raw_by_region.get(rid, {})
    if not units and not _raw:
        continue
    u0 = units[0] if units else _raw
    primary = [u for u in units if u.get("is_primary_lineage", u.get("is_lineage", False))]
    primary_families = {u["family"] for u in primary}
    reference_controls = [u for u in units if u.get("reference_only", False)]
    h3_aux = [u for u in units if u.get("is_h3_auxiliary", False) and u.get("mutation_bearing", False)]
    h4_aux = [u for u in units if u.get("is_h4_auxiliary", False) and u.get("mutation_bearing", False)]
    _cov = _raw.get("col_coverage_by_hp") or {}
    _pop = _raw.get("populations_by_hp") or {}
    _sub = _raw.get("subread_groups_by_hp") or {}
    regions.append({
        "region": rid, "chrom": u0["chrom"], "start": u0["start"], "end": u0["end"],
        "n_sSNV": u0["n_sSNV"], "cn": u0["cn"],
        "analysis_unit": u0.get("analysis_unit"),
        "phase_set": u0.get("phase_set"), "phase_set_status": u0.get("phase_set_status"),
        "source_unit_id": u0.get("source_unit_id", _raw.get("unit_id")),
        "source_component_id": u0.get("source_component_id", _raw.get("component_id")),
        "source_block_id": u0.get("source_block_id", _raw.get("block_id")),
        "hp_multiplicity": len(primary_families), "is_multiHP": len(primary_families) >= 2,
        "hp_multiplicity_definition": "mutation-bearing HP1/HP2 primary lineages",
        "n_primary_lineages": len(primary),
        "n_reference_only_controls": len(reference_controls),
        "n_H3_auxiliary": len(h3_aux),
        "n_H4_auxiliary": len(h4_aux),
        "region_determinacy": region_det(units),
        "positions": _raw.get("positions"),               # [pos1, pos2, ...] 每 sSNV 實際座標
        "n_full_cov_reads": _raw.get("n_full_cov_reads"),  # 跨全部位點的 read 數(co-phase 能力)
        "lineages": sorted([{
            "family": u["family"], "fam_label": u["fam_label"], "is_germline": u["family"] in ("1", "2"),
            "unit_role": u.get("unit_role"), "mutation_bearing": u.get("mutation_bearing"),
            "reference_only": u.get("reference_only"), "is_primary_lineage": u.get("is_primary_lineage"),
            "is_H3_auxiliary": u.get("is_h3_auxiliary"),
            "is_H4_auxiliary": u.get("is_h4_auxiliary"),
            "L1_class": u["L1_class"], "n_trees": u["n_trees"], "n_trees_stored": u.get("n_trees_stored"),
            "display_trees_complete": u.get("display_trees_complete"),
            "analysis_candidate_set_complete": u.get("analysis_candidate_set_complete"),
            "analysis_tree_digest_sha256": u.get("analysis_tree_digest_sha256"),
            "n_distinct_shapes": u.get("n_distinct_shapes"),
            "n_distinct_shapes_exact": u.get("n_distinct_shapes_exact"),
            "n_distinct_shapes_stored": u.get("n_distinct_shapes_stored"),
            "n_hidden": u["n_hidden"], "capped": u["capped"],
            "n_reads": u["n_reads"], "n_full_pops": u["n_full_pops"], "n_partial": u["n_partial"],
            "trees": u["trees"], "L2_cn_verdict": u["L2_cn_verdict"],
            "verification_status": u.get("verification_status"),
            "verification_complete": u.get("verification_complete"),
            "verification_skipped": u.get("verification_skipped"),
            "verify_pass": u["verify_pass"], "trace": u["trace"],
            "phase_set": u.get("phase_set"), "source_unit_id": u.get("source_unit_id"),
            "source_component_id": u.get("source_component_id"),
            "source_block_id": u.get("source_block_id"),
            "vaf_eligible": u.get("vaf_eligible"),
            "coverage_interpretation": u.get("coverage_interpretation"),
            # 原始觀測(此 germline-HP 家族):每位點 [REF,ALT] read 數 / 全跨 populations / partial subreads
            "obs_col_coverage": _cov.get(u["family"]), "obs_populations": _pop.get(u["family"]),
            "obs_subreads": _sub.get(u["family"]),
        } for u in units], key=lambda x: (x["family"] == "none", x["family"] in ("3", "4"), x["family"])),
    })
regions.sort(key=lambda r: (-r["hp_multiplicity"], -r["n_sSNV"]))

# census
lin = [u for u in detail if u.get("is_primary_lineage", u.get("is_lineage", False))]
noncap = [u for u in lin if not u["capped"]]
sum_trees = sum(u["n_trees"] for u in noncap)
ntdist = Counter()
for u in lin:
    n = u["n_trees"]
    ntdist['1' if n == 1 else ('2-5' if n <= 5 else ('6-20' if n <= 20 else '>20'))] += 1
nh = Counter(u["n_hidden"] for u in lin)
census = {"n_regions": len(regions),
          "n_regions_with_primary_lineage": sum(1 for r in regions if r["n_primary_lineages"] > 0),
          "n_regions_without_primary_lineage": sum(1 for r in regions if r["n_primary_lineages"] == 0),
          "n_regions_with_reference_only_control": sum(1 for r in regions if r["n_reference_only_controls"] > 0),
          "n_regions_with_H3_auxiliary": sum(1 for r in regions if r["n_H3_auxiliary"] > 0),
          "n_regions_with_H4_auxiliary": sum(1 for r in regions if r["n_H4_auxiliary"] > 0),
          "hp_multiplicity": dict(Counter(r["hp_multiplicity"] for r in regions)),
          "hp_multiplicity_definition": "mutation-bearing HP1/HP2 primary lineages",
          "region_determinacy": dict(Counter(r["region_determinacy"] for r in regions)),
          "L0": L["L0_hp_family"], "L1": L["L1_ssnv_algorithm"], "L2": L["L2_cn"],
          "U5_trees": {"sum_ntrees_noncapped": sum_trees, "ntrees_dist": dict(ntdist)},
          "U5_shapes": {"exact_for_complete_noncapped": True,
                        "n_units_incomplete": sum(1 for u in lin if not u.get("analysis_candidate_set_complete", False))},
          "U6_hidden": {"sum_hidden": sum(u["n_hidden"] for u in lin), "dist": dict(sorted(nh.items()))},
          "analysis_scope": "chr1-22 primary; chrX/chrY out-of-scope census only",
          "claim_scope": "regional mutation-state trees; not confirmed cell clones"}
census["read_tag_census"] = L.get("read_tag_census", {})
# U1 tree sSNV total: production backbone = LongPhase-S _sc.vcf FILTER=PASS; legacy = census somatic==True.
somatic_vcf = os.environ.get("SM_SOMATIC_VCF")
if somatic_vcf and os.path.exists(somatic_vcf):
    import gzip as _gz
    _op = _gz.open if somatic_vcf.endswith(".gz") else open
    n = 0
    n_auto = 0
    n_out = 0
    auto_chroms = {f"chr{i}" for i in range(1, 23)}
    with _op(somatic_vcf, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            c = ln.split("\t")
            if len(c) > 4 and len(c[3]) == 1 and len(c[4].strip()) == 1:
                n += 1
                if c[0] in auto_chroms:
                    n_auto += 1
                else:
                    n_out += 1
    census["U1_sSNV_somatic_total"] = n
    census["U1_sSNV_scope_chr1_22"] = n_auto
    census["U1_sSNV_out_of_scope"] = n_out
    census["U1_backbone_source"] = os.environ.get(
        "SM_BACKBONE_SOURCE", "LongPhase-S recalibrated FILTER=PASS tree input")
    census["U1_input_contract"] = "LongPhase-S _sc.vcf FILTER=PASS biallelic sSNV"
    f = dict(L.get("input_funnel", {}) or {})
    if f.get("funnel_contract") in {
        "exact_ps_membership_v1",
        "strict_endpoint_exact_ps_membership_v2",
    }:
        f.update({"L1_all_pass_universe": n,
                  "L1_tree_input_PASS_universe": n,
                  "L2_out_of_scope_non_autosomal": n_out,
                  "autosomal_chr1_22": n_auto})
        f["check_autosomal_matches_upstream"] = n_auto == f.get("candidate_loci_S")
        f["check_six_branch_conservation"] = None
        f["six_branch_conservation_status"] = "not_applicable_exact_ps_membership_funnel"
    else:
        f.update({"L1_all_pass_universe": n,
              "L1_tree_input_PASS_universe": n,
              "L2_out_of_scope_non_autosomal": n_out,
              "autosomal_chr1_22": n_auto,
              "L3_positional_singleton": f.get("n_positional_singleton", 0),
              "L4_multilocus_pre_cap_sSNV": f.get("n_multilocus_pre_cap_sSNV", 0),
              "L5_cap_excluded_sSNV": f.get("n_sSNV_cap_excluded", 0),
              "L5_read_unsupported_sSNV": f.get("n_sSNV_read_unsupported", 0),
              "L6_retained_sSNV": f.get("n_sSNV_retained", 0)})
        f["check_autosomal_matches_upstream"] = n_auto == f.get("n_sSNV_scope_input")
        f["check_six_branch_conservation"] = n_auto == (
            f["L3_positional_singleton"] + f["L5_cap_excluded_sSNV"]
            + f["L5_read_unsupported_sSNV"] + f["L6_retained_sSNV"])
    census["funnel"] = f
else:
    cen_path = os.environ.get("SM_CENSUS")
    if cen_path and os.path.exists(cen_path):
        cen = json.load(open(cen_path, encoding="utf-8")).get("census", {})
        census["U1_sSNV_somatic_total"] = sum(1 for v in cen.values() if v.get("somatic") is True)
        census["U1_census_total_positions"] = len(cen)
intg = os.environ.get("SM_INTEGRATION")
if intg and os.path.exists(intg):
    census["U3_linkage_regions_full_span"] = json.load(open(intg, encoding="utf-8"))["aggregate"]["n_regions"]

out = {"schema_version": L.get("schema_version", "2.0"), "sample": SAMPLE,
       "data_model_spec": "docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md",
       "analysis_contract": L.get("analysis_contract"), "L3_methyl": L.get("L3_methyl"),
       "census": census, "regions": regions}
with open(OUT, "w", encoding="utf-8") as handle:
    json.dump(out, handle, ensure_ascii=False)
print(f"OK region-view [{SAMPLE}] n_regions={len(regions)} region_det={dict(Counter(r['region_determinacy'] for r in regions))} -> {OUT}")
