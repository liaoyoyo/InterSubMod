#!/usr/bin/env python3
"""
layered_tree_reconstruction.py — 分層樹重建 driver (2026-07-06 使用者定案)

使用者定案的嚴格分層(HP 家族優先於算法):
  L0 HP 家族  : germline HP 家族(1/2)分群 read;不同家族=不同親代染色體=分開樹;'3'/'4'/'none' 分開報告。
  L1 sSNV+算法: 每 region×germline 家族內,用 tree_enumeration_solver 枚舉全最小樹(§7.5)。= 「只用 sSNV+算法」的範圍與數量。
  L2 CN       : 建完樹「才」用 CN m-通道拆 recurrence(artifact/candidate/LOH)。不覆蓋 L1 結構,只標記。
  L3 甲基     : 事後 bounded-auxiliary 標記(本 driver 佔位;非循環,絕不 rank 樹集)。

每層顯示數字/比例 + 每 region×family 的判斷軌跡(trace)。V1-V7 逐 unit 驗證。

輸入: ml_part JSON(sm_multilocus,含 populations_by_hp / subread_groups_by_hp / col_coverage_by_hp)。
  env SM_ML=<ml.json>(單檔) 或 SM_ML_GLOB=<dir/ml_part_*.json>(多片合併);SM_OUT=<out.json>。
  CN 整數 bed: SM_CN_INT_GAIN / SM_CN_INT_LOSS(預設 HCC1395 SEQC2)。
"""
import sys, os, json, glob, hashlib
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tree_enumeration_solver as S

# ---- L2 CN 整數通道(copy 自 topology_analysis.py gap#3,保持 standalone) ----
CN_INT_GAIN = os.environ.get("SM_CN_INT_GAIN", "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_gain_cn.bed")
CN_INT_LOSS = os.environ.get("SM_CN_INT_LOSS", "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_loss_cn.bed")

def _load_cn_int(path):
    segs = defaultdict(list)
    if path and os.path.exists(path):
        for ln in open(path):
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 4:
                try: segs[p[0]].append((int(p[1]), int(p[2]), float(p[3])))
                except ValueError: pass
    return segs

_CN_GAIN = _load_cn_int(CN_INT_GAIN); _CN_LOSS = _load_cn_int(CN_INT_LOSS)

def _cn_int_lookup(segs, chrom, pos):
    for s, e, v in segs.get(chrom, []):
        if s <= pos <= e: return v
    return None

def m_channel_split(chrom, mid, cn, max_vaf):
    """Interpret recurrence after tree construction; missing CN is never neutral."""
    if cn in (None, "unknown", "unavailable"):
        return "m通道不可用", {"verdict": "unavailable", "cn_int": None}
    gcn = _cn_int_lookup(_CN_GAIN, chrom, mid); lcn = _cn_int_lookup(_CN_LOSS, chrom, mid)
    if gcn is not None and gcn >= 3:
        return "artifact(m>1;CN-amp)", {"verdict": "artifact_drop", "cn_int": gcn}
    if lcn is not None and lcn == 0:
        return "artifact(CN=0)", {"verdict": "artifact_drop", "cn_int": 0.0}
    if cn == "gain":
        return "CN-gain_confounded", {"verdict": "gain_confounded", "cn_int": gcn}
    if cn == "neutral":
        return "candidate(m=1)", {"verdict": "candidate_keep", "cn_int": 2}
    if cn == "loss":
        return "CN-loss_confounded", {"verdict": "loss_confounded", "cn_int": lcn}
    if cn == "loh":
        vf = "likely_artifact(高VAF)" if (max_vaf is not None and max_vaf > 0.7) else "likely_recurrence(低VAF)"
        return "LOH_unresolved", {"verdict": "LOH_unresolved", "cn_int": 2,
                                  "max_vaf": round(max_vaf, 3) if max_vaf else None, "vaf_L3": vf}
    return "m通道不可用", {"verdict": "unavailable", "cn_int": None}

# ---- 分層處理 ----
GERMLINE_FAMS = ("1", "2")
FAMILY_FAMS = ("1", "2", "3", "4")

def _canon_shape(edges):
    """樹的 canonical 形狀字串(忽略節點 label/姊妹順序)→ 用於數「不同拓撲形狀」。
    n 棵等機率樹可能是同一形狀的 label/順序變體;此函式讓 UI 顯示『N 棵樹 = M 種不同形狀』。"""
    ch = {}
    for p, c in edges:
        ch.setdefault(p, []).append(c)
    roots = [n for n in ch if all(n != cc for cs in ch.values() for cc in cs)]
    def rec(n):
        return "(" + "".join(sorted(rec(x) for x in ch.get(n, []))) + ")"
    return "|".join(sorted(rec(r) for r in roots)) if roots else "()"
ANALYSIS_TREE_CAP = int(os.environ.get("SM_ANALYSIS_TREE_CAP", "0"))
DISPLAY_TREE_CAP = int(os.environ.get("SM_DISPLAY_TREE_CAP", "32"))
VERIFY_EVERY = int(os.environ.get("SM_VERIFY_EVERY", "1"))  # 每 N 個 unit 跑完整 V4/V5(昂貴);其餘 light(V1/2/3/6/7)。1=全跑
_unit_idx = [0]              # 全域 unit 計數(抽樣 V4/V5)


def _mutation_bearing(full, partial):
    """A model lineage needs a MINREAD-supported genotype containing ALT."""
    return any("A" in g for g in full) or any("A" in g for g in partial)


def _tree_digest(trees):
    payload = [sorted([list(e) for e in t["edges"]]) for t in trees]
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

def coarse_class(result):
    """L1 粗分類:capped 優先(太密);ambiguous 細分順序未定 vs 多結構;否則 determined/recurrence/underdetermined。
    ambiguous 細分(2026-07-06):feasible_N 只 1 個(同節點集,只 parent/順序不同)=順序未定(突變集合已知);
    >1 個(不同節點集/IDP 多完成)=多結構。"""
    if result.get("capped"):
        return "capped(太密;枚舉未完)"
    cls = S.classify(result)
    if cls == "ambiguous":
        n_node_sets = len(result.get("_feasible_N", []))
        return "ambiguous_order(順序未定;同節點集)" if n_node_sets == 1 else "ambiguous_structure(多完成/多結構)"
    return cls

def max_vaf_of(colcov):
    """col_coverage {pos:[nREF,nALT]} → max ALT 比例。"""
    best = None
    for p, (nr, na) in (colcov or {}).items():
        tot = nr + na
        if tot > 0:
            v = na / tot
            best = v if best is None else max(best, v)
    return best

def process_group(r):
    """一個 region → 逐 germline 家族(+H3/H4/none)產 unit 記錄。回 [unit,...]。"""
    chrom = r["chrom"]; positions = r.get("positions", [])
    start = r.get("start"); end = r.get("end"); cn = r.get("cn")
    region_id = r.get("region_id") or f"{chrom}:{start}-{end}"
    phase_set = r.get("phase_set")
    mid = (positions[0] + positions[-1]) // 2 if positions else start
    pbh = r.get("populations_by_hp", {}) or {}
    sbh = r.get("subread_groups_by_hp", {}) or {}
    cbh = r.get("col_coverage_by_hp", {}) or {}
    reads_by_hp = r.get("reads_by_hp", {}) or {}
    fams = sorted(set(pbh) | set(sbh))
    units = []
    for fam in fams:
        full = pbh.get(fam, {}) or {}
        part = list((sbh.get(fam, {}) or {}).keys())
        if not full and not part:
            continue
        k = len(next(iter(full))) if full else len(part[0])
        res = S.enumerate_min_trees(full, part, k, tree_cap=ANALYSIS_TREE_CAP)
        cls = coarse_class(res)
        _unit_idx[0] += 1
        _light = (VERIFY_EVERY > 1) and (_unit_idx[0] % VERIFY_EVERY != 0)
        ver = S.verify_all(res, full, part, k, light=_light)
        is_germ = fam in GERMLINE_FAMS
        mutation_bearing = _mutation_bearing(full, part)
        reference_only = fam in FAMILY_FAMS and not mutation_bearing
        is_h3_aux = fam == "3"
        is_h4_aux = fam == "4"
        is_primary_lineage = is_germ and mutation_bearing
        is_lineage = is_primary_lineage
        legacy_is_lineage = fam in FAMILY_FAMS
        if reference_only:
            unit_role = "reference_only_control"
        elif is_primary_lineage:
            unit_role = "primary_mutation_lineage"
        elif is_h3_aux:
            unit_role = "unresolved_H3_auxiliary"
        elif is_h4_aux:
            unit_role = "shared_H4_auxiliary"
        else:
            unit_role = "unphased_auxiliary"
        fam_label = {"1": "germline-hap1", "2": "germline-hap2",
                     "3": "H3?(somatic-integrated unresolved)",
                     "4": "H4?(somatic ALT shared by germline HP1/HP2)"}.get(fam, "unphased(none)")
        colcov = cbh.get(fam, {})
        mvaf = max_vaf_of(colcov) if r.get("vaf_eligible", True) else None
        # L2 CN(只對 recurrence_required 非 capped)
        l2_label, l2_dict = "n/a(非recurrence)", None
        base_cls = S.classify(res)
        if mutation_bearing and base_cls == "recurrence_required" and not res.get("capped"):
            l2_label, l2_dict = m_channel_split(chrom, mid, cn, mvaf)
        elif res.get("capped") and "recurrence" in cls:
            l2_label = "n/a(capped;先不送CN)"
        executed = [name for name in ("V1", "V2", "V3", "V4", "V5", "V6", "V7") if ver[name][0] is not None]
        skipped = [name for name in ("V1", "V2", "V3", "V4", "V5", "V6", "V7") if ver[name][0] is None]
        failed = [name for name in executed if ver[name][0] is False]
        verification_complete = not res.get("capped") and not skipped and not failed
        if failed:
            verification_status = "fail"
        elif res.get("capped"):
            verification_status = "not_applicable_capped"
        elif verification_complete:
            verification_status = "full_pass"
        else:
            verification_status = "partial_pass"
        analysis_complete = not res.get("capped") and res.get("trees_complete", False)
        all_trees = res["trees"]
        display_trees = all_trees if DISPLAY_TREE_CAP <= 0 else all_trees[:DISPLAY_TREE_CAP]
        exact_shapes = len({_canon_shape(t["edges"]) for t in all_trees}) if analysis_complete else None
        stored_shapes = len({_canon_shape(t["edges"]) for t in display_trees})
        trace = [
            f"L0: unit={region_id}, exact_PS={phase_set if phase_set is not None else 'legacy-unbounded'}, HP家族={fam}({fam_label}), role={unit_role}, {reads_by_hp.get(fam,'?')} reads, {len(full)} full-pop / {len(part)} partial",
            f"L1: sSNV+算法 → {cls}  (n_trees={res['n_trees']}, n_hidden={res['n_hidden']}, capped={res['capped']})",
            f"L2: cn={cn} → {l2_label}",
            "L3: not_evaluated(bounded-auxiliary;禁止 rank/confirm 樹集)",
        ]
        units.append({
            "region": region_id, "chrom": chrom, "start": start, "end": end,
            "analysis_unit": r.get("analysis_unit", "legacy_geometric_region"),
            "phase_set": phase_set, "phase_set_status": r.get("phase_set_status"),
            "source_unit_id": r.get("unit_id"), "source_component_id": r.get("component_id"),
            "source_block_id": r.get("block_id"), "source_block_index": r.get("block_index"),
            "linkage_basis": r.get("linkage_basis"),
            "family": fam, "is_germline_family": is_germ, "is_family_unit": fam in FAMILY_FAMS,
            "mutation_bearing": mutation_bearing, "reference_only": reference_only,
            "is_h3_auxiliary": is_h3_aux, "is_h4_auxiliary": is_h4_aux,
            "is_primary_lineage": is_primary_lineage,
            "is_lineage": is_lineage, "legacy_is_lineage": legacy_is_lineage,
            "unit_role": unit_role, "fam_label": fam_label,
            "n_sSNV": r.get("n_sSNV"), "cn": cn,
            "n_reads": reads_by_hp.get(fam), "n_full_pops": len(full), "n_partial": len(part),
            "L1_class": cls, "L1_base_class": base_cls, "n_trees": res["n_trees"],
            "n_hidden": res["n_hidden"], "capped": res["capped"], "cap_reason": res.get("cap_reason"),
            "trees": display_trees, "n_trees_stored": len(display_trees),
            "display_trees_complete": len(display_trees) == res["n_trees"],
            "analysis_trees_generated": len(all_trees), "analysis_candidate_set_complete": analysis_complete,
            "analysis_tree_digest_sha256": _tree_digest(all_trees) if analysis_complete else None,
            "n_distinct_shapes": exact_shapes, "n_distinct_shapes_exact": exact_shapes,
            "n_distinct_shapes_stored": stored_shapes,
            "L2_cn_verdict": l2_label, "L2_m_channel": l2_dict,
            "verify": {kk: vv[0] for kk, vv in ver.items() if kk != "overall"},
            "verification_status": verification_status, "verification_complete": verification_complete,
            "verification_executed": executed, "verification_skipped": skipped,
            "verification_failed": failed, "verify_pass": verification_status == "full_pass",
            "max_vaf": round(mvaf, 3) if mvaf is not None else None,
            "vaf_eligible": bool(r.get("vaf_eligible", True)),
            "coverage_interpretation": r.get("coverage_interpretation"),
            "trace": trace,
        })
    return units


def _merge_funnels(docs):
    numeric = Counter()
    rules = set()
    all_checks = True
    contracts = set()
    for doc in docs:
        f = doc.get("input_funnel", {}) or {}
        if f.get("funnel_contract"):
            contracts.add(f["funnel_contract"])
        for key, value in f.items():
            if isinstance(value, int) and not isinstance(value, bool):
                numeric[key] += value
        if f.get("grouping_rule"):
            rules.add(f["grouping_rule"])
        if f.get("funnel_contract") != "exact_ps_membership_v1":
            all_checks = all_checks and bool(f.get("check_scope_conservation", False))
    out = dict(numeric)
    if contracts == {"exact_ps_membership_v1"}:
        out["funnel_contract"] = "exact_ps_membership_v1"
        out["check_scope_conservation"] = all(
            bool((doc.get("input_funnel") or {}).get("check_constraint_weight_conservation"))
            for doc in docs
        )
    else:
        if contracts:
            out["funnel_contract"] = sorted(contracts)
        out["check_scope_conservation"] = (
            all_checks and out.get("n_sSNV_scope_input", 0) == out.get("n_sSNV_accounted", -1))
    out["grouping_rule"] = sorted(rules)
    return out

def _merge_tag_census(docs):
    numeric = Counter()
    nested = defaultdict(Counter)
    frozen_fields = (
        "identity_schema", "sidecar_sha256", "sidecar_index_sha256",
        "alignment_payload_identity_sha256", "sidecar_duplicates",
        "duplicate_identity_policy", "sidecar_extra", "sidecar_malformed",
        "phase_set_policy",
    )
    frozen = {key: set() for key in frozen_fields}
    for doc in docs:
        census = doc.get("read_tag_census", {}) or {}
        for key, value in census.items():
            if key in frozen:
                frozen[key].add(value)
            elif isinstance(value, int) and not isinstance(value, bool):
                numeric[key] += value
            elif isinstance(value, dict):
                nested[key].update(value)
    out = dict(numeric)
    out.update({key: dict(value) for key, value in nested.items()})
    exact_ps_mode = bool(docs) and all(
        (doc.get("params") or {}).get("input_mode") == "exact_ps_partition_adapter"
        for doc in docs
    )
    if exact_ps_mode:
        phase_set_policies = {
            (doc.get("read_tag_census") or {}).get("phase_set_policy") for doc in docs
        }
        phase_set_policies.discard(None)
        out["phase_set_policy"] = (
            next(iter(phase_set_policies)) if len(phase_set_policies) == 1 else None
        )
        out["evidence_binding"] = (
            "validated adapter receipt + partition run receipt + per-chromosome stage receipts"
        )
        out["check_partition_receipts_all_pass"] = all(
            all(
                bool((doc.get("input_funnel") or {}).get(key))
                for key in (
                    "check_constraint_weight_conservation",
                    "check_cpp_parity_required",
                    "check_cross_hp_zero",
                    "check_cross_ps_zero",
                    "check_stable_region_ids_unique",
                )
            )
            for doc in docs
        )
        out["check_exact_sidecar_join"] = None
        out["check_exact_sidecar_join_status"] = (
            "not_applicable; exact-PS adapter consumes already validated partition evidence"
        )
        return out
    frozen_consistent = all(len(values) == 1 for values in frozen.values())
    for key, values in frozen.items():
        out[key] = next(iter(values)) if len(values) == 1 else None
    out["check_frozen_sidecar_binding_consistent"] = frozen_consistent
    out["check_exact_sidecar_join"] = (
        out.get("sidecar_missing", 0) == 0 and out.get("sidecar_conflicts", 0) == 0
        and out.get("alignment_identity_allele_conflicts", 0) == 0
        and out.get("sidecar_exact_matches", -1) == out.get("alignment_group_exposures", -2)
        and out.get("duplicate_identity_policy") == "collapse_redundant_rows_with_identical_HP_PS"
        and frozen_consistent)
    return out

def main():
    ml_path = os.environ.get("SM_ML")
    ml_glob = os.environ.get("SM_ML_GLOB")
    out_path = os.environ.get("SM_OUT", "layered_reconstruction.json")
    groups = []
    ml_docs = []
    if ml_path:
        with open(ml_path, encoding="utf-8") as handle:
            doc = json.load(handle)
        ml_docs.append(doc)
        groups = doc["groups"]
    elif ml_glob:
        for f in sorted(glob.glob(ml_glob)):
            with open(f, encoding="utf-8") as handle:
                doc = json.load(handle)
            ml_docs.append(doc)
            groups += doc["groups"]
    else:
        sys.exit("需 SM_ML 或 SM_ML_GLOB")
    regs = [g for g in groups if g.get("n_sSNV", 0) >= 2]

    # L0 家族分布
    L0 = {"regions_with_germ_family": Counter(), "regions_mixed_germline": 0,
          "reads_by_family_total": Counter(), "regions_total": len(regs)}
    for r in regs:
        fams = set(r.get("populations_by_hp", {})) | set(r.get("subread_groups_by_hp", {}))
        germ = [f for f in fams if f in GERMLINE_FAMS]
        for f in fams:
            L0["regions_with_germ_family"][f] += 1
        if len(germ) >= 2:
            L0["regions_mixed_germline"] += 1
        for f, n in (r.get("reads_by_hp", {}) or {}).items():
            L0["reads_by_family_total"][f] += n

    # L1/L2 per unit
    detail = []
    for r in regs:
        detail.extend(process_group(r))

    family_units = [u for u in detail if u["is_family_unit"]]
    primary_units = [u for u in detail if u["is_primary_lineage"]]
    reference_units = [u for u in detail if u["reference_only"]]
    h3_units = [u for u in detail if u["is_h3_auxiliary"] and u["mutation_bearing"]]
    h4_units = [u for u in detail if u["is_h4_auxiliary"] and u["mutation_bearing"]]
    germ_units_all = [u for u in detail if u["is_germline_family"]]
    none_units = [u for u in detail if u["family"] == "none"]
    eligible = [u for u in detail if not u["capped"]]
    primary_eligible = [u for u in primary_units if not u["capped"]]
    status_counts = Counter(u["verification_status"] for u in detail)
    eligible_skipped = [u for u in eligible if u["verification_skipped"]]
    L1 = {"n_family_units_total": len(family_units),
          "n_units_total_including_unphased": len(detail),
          "n_primary_lineage_units": len(primary_units),
          "n_lineage_units": len(primary_units),
          "n_reference_only_controls": len(reference_units),
          "n_unresolved_H3_auxiliary": len(h3_units),
          "n_shared_H4_auxiliary": len(h4_units),
          "n_germline_family_units_all": len(germ_units_all),
          "n_none_units": len(none_units),
          "determinacy_primary_lineage": Counter(u["L1_class"] for u in primary_units),
          "determinacy_lineage": Counter(u["L1_class"] for u in primary_units),
          "determinacy_germline_1_2_all": Counter(u["L1_class"] for u in germ_units_all),
          "determinacy_H3_auxiliary": Counter(u["L1_class"] for u in h3_units),
          "determinacy_H4_auxiliary": Counter(u["L1_class"] for u in h4_units),
          "determinacy_reference_controls": Counter(u["L1_class"] for u in reference_units),
          "determinacy_unphased_none": Counter(u["L1_class"] for u in none_units),
          "verification_status": status_counts,
          "n_verification_fail": status_counts.get("fail", 0),
          "n_verification_not_applicable_capped": status_counts.get("not_applicable_capped", 0),
          "n_eligible_units": len(eligible),
          "n_eligible_skipped_V4V5": len(eligible_skipped),
          "all_eligible_V1V7_pass": bool(eligible) and all(u["verification_status"] == "full_pass" for u in eligible),
          "all_primary_noncapped_V1V7_pass": bool(primary_eligible) and all(u["verification_status"] == "full_pass" for u in primary_eligible),
          "verification_scope": "all non-capped units; capped=not_applicable_capped"}
    L1["all_V1V7_pass"] = L1["all_eligible_V1V7_pass"]
    L1["proportion_primary_lineage"] = {
        k: round(v / max(1, len(primary_units)), 4) for k, v in L1["determinacy_primary_lineage"].items()}
    L1["proportion_lineage"] = dict(L1["proportion_primary_lineage"])

    primary_by_region = defaultdict(set)
    reference_by_region = defaultdict(set)
    h3_by_region = defaultdict(set)
    h4_by_region = defaultdict(set)
    for u in detail:
        if u["is_primary_lineage"]:
            primary_by_region[u["region"]].add(u["family"])
        if u["reference_only"]:
            reference_by_region[u["region"]].add(u["family"])
        if u["is_h3_auxiliary"] and u["mutation_bearing"]:
            h3_by_region[u["region"]].add(u["family"])
        if u["is_h4_auxiliary"] and u["mutation_bearing"]:
            h4_by_region[u["region"]].add(u["family"])
    L0["regions_with_primary_lineage"] = len(primary_by_region)
    L0["regions_mixed_primary_HP1_HP2"] = sum(1 for fams in primary_by_region.values() if len(fams) >= 2)
    L0["regions_with_reference_only_control"] = len(reference_by_region)
    L0["regions_with_H3_auxiliary"] = len(h3_by_region)
    L0["regions_with_H4_auxiliary"] = len(h4_by_region)

    rec_units = [u for u in primary_units if u["L1_base_class"] == "recurrence_required" and not u["capped"]]
    rec_h3 = [u for u in h3_units if u["L1_base_class"] == "recurrence_required" and not u["capped"]]
    L2 = {"n_primary_recurrence_sent_to_cn": len(rec_units),
          "n_recurrence_sent_to_cn": len(rec_units),
          "primary_cn_split": Counter(u["L2_cn_verdict"] for u in rec_units),
          "cn_split": Counter(u["L2_cn_verdict"] for u in rec_units),
          "n_H3_auxiliary_recurrence_annotated": len(rec_h3),
          "H3_auxiliary_cn_split": Counter(u["L2_cn_verdict"] for u in rec_h3)}

    exact_ps_mode = bool(ml_docs) and all(
        (doc.get("params") or {}).get("input_mode") == "exact_ps_partition_adapter"
        for doc in ml_docs
    )
    strict_endpoint_mode = exact_ps_mode and all(
        (doc.get("params") or {}).get("strict_endpoint_receipt_required") is True
        for doc in ml_docs
    )
    input_samples = {doc.get("sample") for doc in ml_docs if doc.get("sample")}
    if len(input_samples) > 1:
        raise RuntimeError(f"multiple samples in layered input: {sorted(input_samples)}")
    input_sample = next(iter(input_samples)) if input_samples else None
    out = {"schema_version": "2.1" if exact_ps_mode else "2.0",
           "sample": input_sample,
           "layer_note": "L0 family units → L1 primary mutation-bearing HP1/HP2 trees → L2 post-tree CN annotation → L3 bounded auxiliary",
           "analysis_contract": {"primary_lineage": "mutation-bearing HP1/HP2 only",
                                 "reference_only": "control, excluded from lineage/determinacy/multi-HP denominators",
                                 "HP3": "H3? unresolved auxiliary, excluded from primary denominators",
                                 "HP4": "H4? shared-somatic auxiliary, excluded from primary denominators",
                                 "PS": ("exact non-missing primary evidence boundary; cross-PS pooling forbidden"
                                        if exact_ps_mode else
                                        "preserved for phase-block QC; not used as a topology edge or lineage label"),
                                 "analysis_unit": ("exact PS x HP x read-linked component x bounded block"
                                                   if exact_ps_mode else "legacy geometric region x HP family"),
                                 "read_weight_role": ("MINREAD and segmentation audit only; solver uses pattern presence"
                                                      if exact_ps_mode else "MINREAD filtering"),
                                 "claim_scope": "regional mutation-state trees, not confirmed cell clones",
                                 "region_rule": (
                                     ("exact PS x HP strict fixed-R/A endpoint-pair connected component; "
                                      "every retained graph edge has direct distinct-molecule support; "
                                      "k>12 scientific components are split only into bounded k<=12 "
                                      "computational blocks; no cross-PS join")
                                     if strict_endpoint_mode else
                                     ("exact PS x HP upstream read-linked component; k>12 units are split into "
                                      "bounded k<=12 blocks by retained-read-weight DP; no cross-PS join")
                                     if exact_ps_mode else
                                     "adjacent-gap connected component; total span may exceed threshold")},
           "params": {"VERIFY_EVERY": VERIFY_EVERY, "ANALYSIS_TREE_CAP": ANALYSIS_TREE_CAP,
                      "DISPLAY_TREE_CAP": DISPLAY_TREE_CAP,
                      "input_params": [d.get("params", {}) for d in ml_docs]},
           "input_funnel": _merge_funnels(ml_docs),
           "read_tag_census": _merge_tag_census(ml_docs),
           "L0_hp_family": {k: (dict(v) if isinstance(v, Counter) else v) for k, v in L0.items()},
           "L1_ssnv_algorithm": {k: (dict(v) if isinstance(v, Counter) else v) for k, v in L1.items()},
           "L2_cn": {k: (dict(v) if isinstance(v, Counter) else v) for k, v in L2.items()},
           "L3_methyl": {"status": "not_evaluated", "role": "bounded_auxiliary",
                         "allowed_uses": ["negative_screen", "residual_flag"],
                         "prohibited_uses": ["tree_ranking", "lineage_confirmation", "clone_confirmation"],
                         "reason": "No orthogonal per-unit methyl validation is joined in this driver"},
           "n_detail_units": len(detail), "detail": detail}
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(out, handle, ensure_ascii=False)

    print("=" * 84)
    print("分層樹重建 census")
    print("=" * 84)
    print(f"L0 家族: regions={L0['regions_total']}  含各家族區={dict(L0['regions_with_germ_family'])}  混合germline(1&2)={L0['regions_mixed_germline']}")
    print(f"    reads/family total={dict(L0['reads_by_family_total'])}")
    print(f"L1 primary mutation lineages={len(primary_units)}; reference controls={len(reference_units)}; "
          f"H3? auxiliary={len(h3_units)}; H4? auxiliary={len(h4_units)}")
    for k, v in sorted(L1["determinacy_primary_lineage"].items(), key=lambda x: -x[1]):
        print(f"    {v:6d} ({L1['proportion_primary_lineage'][k]*100:5.1f}%)  {k}")
    print(f"    [germline1/2 all family units] {dict(L1['determinacy_germline_1_2_all'])}")
    print(f"    [H3? auxiliary] {dict(L1['determinacy_H3_auxiliary'])}")
    print(f"    [H4? auxiliary] {dict(L1['determinacy_H4_auxiliary'])}")
    print(f"    [reference-only controls] {dict(L1['determinacy_reference_controls'])}")
    print(f"    [unphased none, 分開] {dict(L1['determinacy_unphased_none'])}")
    v7line = ("ALL ELIGIBLE PASS" if L1["all_eligible_V1V7_pass"]
              else f"fail={L1['n_verification_fail']} eligible_skipped={L1['n_eligible_skipped_V4V5']}")
    print(f"    V1-V7(non-capped eligible): {v7line}; capped n/a={L1['n_verification_not_applicable_capped']}")
    print(f"L2 CN(primary recurrence→CN, n={len(rec_units)}):  {dict(L2['primary_cn_split'])}")
    print(f"→ {out_path}  ({len(detail)} units)")

if __name__ == "__main__":
    main()
