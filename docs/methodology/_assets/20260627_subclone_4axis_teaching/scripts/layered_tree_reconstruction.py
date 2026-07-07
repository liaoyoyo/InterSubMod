#!/usr/bin/env python3
"""
layered_tree_reconstruction.py — 分層樹重建 driver (2026-07-06 使用者定案)

使用者定案的嚴格分層(HP 家族優先於算法):
  L0 HP 家族  : germline HP 家族(1/2)分群 read;不同家族=不同親代染色體=分開樹;'3'/'none' 分開報告。
  L1 sSNV+算法: 每 region×germline 家族內,用 tree_enumeration_solver 枚舉全最小樹(§7.5)。= 「只用 sSNV+算法」的範圍與數量。
  L2 CN       : 建完樹「才」用 CN m-通道拆 recurrence(artifact/candidate/LOH)。不覆蓋 L1 結構,只標記。
  L3 甲基     : 事後 bounded-auxiliary 標記(本 driver 佔位;非循環,絕不 rank 樹集)。

每層顯示數字/比例 + 每 region×family 的判斷軌跡(trace)。V1-V7 逐 unit 驗證。

輸入: ml_part JSON(sm_multilocus,含 populations_by_hp / subread_groups_by_hp / col_coverage_by_hp)。
  env SM_ML=<ml.json>(單檔) 或 SM_ML_GLOB=<dir/ml_part_*.json>(多片合併);SM_OUT=<out.json>。
  CN 整數 bed: SM_CN_INT_GAIN / SM_CN_INT_LOSS(預設 HCC1395 SEQC2)。
"""
import sys, os, json, glob, math
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
    """L2:回 (verdict_label, m_dict)。整數 CN≥3(gain)/CN=0 → artifact;neutral/loss≥1 → candidate;LOH → VAF 軟旗標。"""
    if cn == "unknown":
        return "m通道不可用", {"verdict": "unavailable", "cn_int": None}
    gcn = _cn_int_lookup(_CN_GAIN, chrom, mid); lcn = _cn_int_lookup(_CN_LOSS, chrom, mid)
    if gcn is not None and gcn >= 3:
        return "artifact(m>1;CN-amp)", {"verdict": "artifact_drop", "cn_int": gcn}
    if lcn is not None and lcn == 0:
        return "artifact(CN=0)", {"verdict": "artifact_drop", "cn_int": 0.0}
    if cn == "neutral" or (cn == "loss" and lcn is not None and lcn >= 1):
        return "candidate(m=1)", {"verdict": "candidate_keep", "cn_int": (lcn if cn == "loss" else 2)}
    if cn == "loh":
        vf = "likely_artifact(高VAF)" if (max_vaf is not None and max_vaf > 0.7) else "likely_recurrence(低VAF)"
        return "LOH_unresolved", {"verdict": "LOH_unresolved", "cn_int": 2,
                                  "max_vaf": round(max_vaf, 3) if max_vaf else None, "vaf_L3": vf}
    return "m通道不可用", {"verdict": "unavailable", "cn_int": None}

# ---- 分層處理 ----
GERMLINE_FAMS = ("1", "2")        # germline HP 家族(有 germline 根錨)
LINEAGE_FAMS = ("1", "2", "3")    # lineage(樹分開建):1/2=germline、3=somatic-integrated(使用者定案 2026-07-06 當第三 lineage)
#   none(unphased)=無家族,分開報告不建 lineage 樹

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
MAX_STORE_TREES = 32         # 每 unit 存樹上限(供 UI 切換;n_trees 記全數;>此數註明「顯示前 N of M」)
VERIFY_EVERY = int(os.environ.get("SM_VERIFY_EVERY", "1"))  # 每 N 個 unit 跑完整 V4/V5(昂貴);其餘 light(V1/2/3/6/7)。1=全跑
_unit_idx = [0]              # 全域 unit 計數(抽樣 V4/V5)

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
    """一個 region → 逐 germline 家族(+3/none)產 unit 記錄。回 [unit,...]。"""
    chrom = r["chrom"]; positions = r.get("positions", [])
    start = r.get("start"); end = r.get("end"); cn = r.get("cn")
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
        res = S.enumerate_min_trees(full, part, k)
        cls = coarse_class(res)
        _unit_idx[0] += 1
        _light = (VERIFY_EVERY > 1) and (_unit_idx[0] % VERIFY_EVERY != 0)
        ver = S.verify_all(res, full, part, k, light=_light)
        is_germ = fam in GERMLINE_FAMS
        is_lineage = fam in LINEAGE_FAMS   # 1/2/3 建 lineage 樹;none 不算
        fam_label = {"1": "germline-hap1", "2": "germline-hap2", "3": "somatic-integrated(第三lineage)"}.get(fam, "unphased(none)")
        colcov = cbh.get(fam, {})
        mvaf = max_vaf_of(colcov)
        # L2 CN(只對 recurrence_required 非 capped)
        l2_label, l2_dict = "n/a(非recurrence)", None
        base_cls = S.classify(res)
        if base_cls == "recurrence_required" and not res.get("capped"):
            l2_label, l2_dict = m_channel_split(chrom, mid, cn, mvaf)
        elif res.get("capped") and "recurrence" in cls:
            l2_label = "n/a(capped;先不送CN)"
        # trace(每區判斷軌跡)
        trace = [
            f"L0: HP家族={fam}({fam_label}), {reads_by_hp.get(fam,'?')} reads, {len(full)} full-pop / {len(part)} partial",
            f"L1: sSNV+算法 → {cls}  (n_trees={res['n_trees']}, n_hidden={res['n_hidden']}, capped={res['capped']})",
            f"L2: cn={cn} → {l2_label}",
            "L3: 甲基 pending(bounded-auxiliary;不 rank 樹集)",
        ]
        units.append({
            "region": f"{chrom}:{start}-{end}", "chrom": chrom, "start": start, "end": end,
            "family": fam, "is_germline_family": is_germ, "is_lineage": is_lineage, "fam_label": fam_label,
            "n_sSNV": r.get("n_sSNV"), "cn": cn,
            "n_reads": reads_by_hp.get(fam), "n_full_pops": len(full), "n_partial": len(part),
            "L1_class": cls, "L1_base_class": base_cls, "n_trees": res["n_trees"],
            "n_hidden": res["n_hidden"], "capped": res["capped"], "cap_reason": res.get("cap_reason"),
            "trees": res["trees"][:MAX_STORE_TREES], "n_trees_stored": min(len(res["trees"]), MAX_STORE_TREES),
            "n_distinct_shapes": len({_canon_shape(t["edges"]) for t in res["trees"][:MAX_STORE_TREES]}),
            "L2_cn_verdict": l2_label, "L2_m_channel": l2_dict,
            "verify": {kk: vv[0] for kk, vv in ver.items() if kk != "overall"},
            "verify_pass": ver["overall"], "max_vaf": round(mvaf, 3) if mvaf is not None else None,
            "trace": trace,
        })
    return units

def main():
    ml_path = os.environ.get("SM_ML")
    ml_glob = os.environ.get("SM_ML_GLOB")
    out_path = os.environ.get("SM_OUT", "layered_reconstruction.json")
    groups = []
    if ml_path:
        groups = json.load(open(ml_path))["groups"]
    elif ml_glob:
        for f in sorted(glob.glob(ml_glob)):
            groups += json.load(open(f))["groups"]
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

    lineage_units = [u for u in detail if u["is_lineage"]]              # 1/2/3 建 lineage 樹
    germ_units = [u for u in detail if u["is_germline_family"]]         # 1/2(有 germline 根)
    som_units = [u for u in detail if u["family"] == "3"]               # 3(somatic-integrated 第三 lineage)
    none_units = [u for u in detail if u["family"] == "none"]           # unphased(分開報告)
    L1 = {"n_family_units_total": len(detail),
          "n_lineage_units": len(lineage_units),
          "determinacy_lineage": Counter(u["L1_class"] for u in lineage_units),
          "determinacy_germline_1_2": Counter(u["L1_class"] for u in germ_units),
          "determinacy_somatic_3": Counter(u["L1_class"] for u in som_units),
          "determinacy_unphased_none": Counter(u["L1_class"] for u in none_units),
          "n_germline_units": len(germ_units), "n_somatic3_units": len(som_units), "n_none_units": len(none_units),
          "n_verify_fail": sum(1 for u in detail if not u["verify_pass"]),
          "all_V1V7_pass": all(u["verify_pass"] for u in detail)}
    L1["proportion_lineage"] = {k: round(v / max(1, len(lineage_units)), 4) for k, v in L1["determinacy_lineage"].items()}

    rec_units = [u for u in lineage_units if u["L1_base_class"] == "recurrence_required" and not u["capped"]]
    L2 = {"n_recurrence_sent_to_cn": len(rec_units),
          "cn_split": Counter(u["L2_cn_verdict"] for u in rec_units)}

    out = {"layer_note": "L0 HP家族 → L1 sSNV+算法(枚舉全最小樹) → L2 CN → L3 甲基;家族優先於算法",
           "L0_hp_family": {k: (dict(v) if isinstance(v, Counter) else v) for k, v in L0.items()},
           "L1_ssnv_algorithm": {k: (dict(v) if isinstance(v, Counter) else v) for k, v in L1.items()},
           "L2_cn": {k: (dict(v) if isinstance(v, Counter) else v) for k, v in L2.items()},
           "L3_methyl": {"status": "pending(bounded-auxiliary;非循環;不 rank 樹集)"},
           "n_detail_units": len(detail), "detail": detail}
    json.dump(out, open(out_path, "w"), ensure_ascii=False)

    print("=" * 84)
    print("分層樹重建 census")
    print("=" * 84)
    print(f"L0 家族: regions={L0['regions_total']}  含各家族區={dict(L0['regions_with_germ_family'])}  混合germline(1&2)={L0['regions_mixed_germline']}")
    print(f"    reads/family total={dict(L0['reads_by_family_total'])}")
    print(f"L1 sSNV+算法(lineage units={len(lineage_units)} = germline1/2:{len(germ_units)} + somatic3:{len(som_units)}):")
    for k, v in sorted(L1["determinacy_lineage"].items(), key=lambda x: -x[1]):
        print(f"    {v:6d} ({L1['proportion_lineage'][k]*100:5.1f}%)  {k}")
    print(f"    [germline1/2] {dict(L1['determinacy_germline_1_2'])}")
    print(f"    [somatic3]    {dict(L1['determinacy_somatic_3'])}")
    print(f"    [unphased none, 分開] {dict(L1['determinacy_unphased_none'])}")
    v7line = "ALL PASS ✓" if L1["all_V1V7_pass"] else f"{L1['n_verify_fail']} FAIL ✗"
    print(f"    V1-V7: {v7line}")
    print(f"L2 CN(recurrence→CN, n={len(rec_units)}):  {dict(L2['cn_split'])}")
    print(f"→ {out_path}  ({len(detail)} units)")

if __name__ == "__main__":
    main()
