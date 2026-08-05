"""[多位點組合 census] 從 2-locus pairwise 升級到 multi-locus (≥3 sSNV) 的 per-read 基因型組合。
每個 group 取 somatic sSNV (census somatic==True)，tumor per-read 抽各 sSNV 等位 → 對「覆蓋全部 sSNV」的 read
建基因型向量字串 (R/A) → 不同向量 = 局部細胞 population (如 chr17: RRRR=ancestral / γ / α-only / α+β)。
每 group 記 n_sSNV / n_full_cov_reads / n_populations / 各 combo:count / CN state / sameHP。
genome-wide 聚合: 每 group population 數分布 + 比例 (≥2/≥3/≥4) + CN 分層 (clean LOH/neutral vs gain)。
argv: chroms out_path。輸出含 per-group list。
"""
import sys
import json
from collections import Counter, defaultdict
import os
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # import 同夾的(env-var 化) sm_linkage
import sm_linkage_genomewide as M

# CNBED env-var; an absent file is missing data, never a neutral measurement.
CNBED = os.environ.get("SM_CNBED", "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed")
CN_SOURCE = os.environ.get("SM_CN_SOURCE", "unspecified")
MINREAD = int(os.environ.get("SM_MINREAD", "3"))
MAX_SNV = int(os.environ.get("SM_MAX_SNV", "8"))
# 骨幹來源(2026-07-11 使用者定案:改用 LongPhase-S recalibrated PASS,移除 is_somatic 粗重檢)。
#   SM_SOMATIC_VCF 設 = LongPhase-S _sc.vcf FILTER=PASS（ClairS PASS 是上游 LongPhase-S input）→
#   取代舊 filtered_snv_tp/fp(SEQC2 區限制)+census is_somatic(粗閾值誤殺真 somatic,見 20260706 資料模型 spec)。
#   未設 → 舊行為(load_union tp/fp + is_somatic),供對照/重現。
SOMATIC_VCF = os.environ.get("SM_SOMATIC_VCF", "")
FUNNEL_FIELDS = (
    "n_sSNV_scope_input",
    "n_positional_singleton",
    "n_multilocus_pre_cap_groups",
    "n_multilocus_pre_cap_sSNV",
    "n_groups_retained",
    "n_groups_read_unsupported",
    "n_sSNV_retained",
    "n_sSNV_read_unsupported",
    "n_groups_capped_by_MAX_SNV",
    "n_sSNV_cap_excluded",
    "n_sSNV_accounted",
)


def new_funnel():
    return Counter({field: 0 for field in FUNNEL_FIELDS})


def load_somatic_universe(chrom):
    """從 SM_SOMATIC_VCF(LongPhase-S recalibrated PASS tree input)取該 chrom 全 SNV 位點。"""
    idx = (SOMATIC_VCF + ".tbi") if os.path.exists(SOMATIC_VCF + ".tbi") else (
        (SOMATIC_VCF + ".csi") if os.path.exists(SOMATIC_VCF + ".csi") else None)
    out = []
    try:
        with (pysam.TabixFile(SOMATIC_VCF, index=idx) if idx else pysam.TabixFile(SOMATIC_VCF)) as tb:
            for ln in tb.fetch(chrom):
                p = ln.split("\t")
                ref, alt = p[3].upper(), p[4].strip().upper()
                if len(ref) == 1 and len(alt) == 1:  # Biallelic SNV-only operational universe.
                    out.append((int(p[1]), ref, alt, "PASS"))
    except (ValueError, OSError):
        return []
    return sorted(out)


def load_cn():
    if not CNBED or not os.path.exists(CNBED):
        return None
    iv = defaultdict(list)
    for ln in open(CNBED):
        p = ln.split()
        if len(p) < 4 or not p[1].isdigit():
            continue
        iv[p[0]].append((int(p[1]), int(p[2]), p[3].lower()))
    for c in iv:
        iv[c].sort()
    return iv


def cn_state(cn, chrom, pos):
    if cn is None:
        return "unavailable"
    for s, e, lab in cn.get(chrom, []):
        if s <= pos <= e:
            return lab
    return "neutral"


def germ_family(hptag):
    """L0 germline HP 家族(2026-07-06 使用者定案:家族優先於算法)。
    HP tag 完整語彙: ./1/2/3/4/1-1/2-1/1-2/2-2。
    germline 家族 = 第一碼:1/1-1/1-2→'1'、2/2-1/2-2→'2';3/4 各自為 unresolved/shared auxiliary。
    🔴 只有同 germline 家族的 read 才依 sSNV 連接建樹(不同家族=不同親代染色體=分開樹)。"""
    if hptag is None:
        return "none"
    s = str(hptag)
    if s.startswith("1"):
        return "1"
    if s.startswith("2"):
        return "2"
    if s == "3":
        return "3"
    if s == "4":
        return "4"
    return "none"


def main(chroms, out_path):
    # SM_SOMATIC_VCF mode uses the explicit production tree backbone; legacy mode reuses census/is_somatic.
    cen = {} if SOMATIC_VCF else json.load(open(f"{M.A}/sm_linkage_genomewide.json"))["census"]
    cn = load_cn()
    tb = pysam.AlignmentFile(M.TBAM, "rb")
    groups_out = []
    agg = Counter()
    agg_clean = Counter()
    agg_neutral = Counter()
    tag_diagnostics = Counter()
    raw_hp_total = Counter()
    raw_hp_ps_total = Counter()
    alignment_class_total = Counter()
    phase_set_region_counts = Counter()
    funnel = new_funnel()
    for chrom in chroms:
        snvs = load_somatic_universe(chrom) if SOMATIC_VCF else M.load_union(chrom)
        for g in M.make_groups(snvs):
            # Production tree VCF is already caller + LongPhase-S filtered; legacy mode uses census somatic.
            som = list(g) if SOMATIC_VCF else [s for s in g if cen.get(f"{chrom}:{s[0]}", {}).get("somatic") is True]
            funnel["n_sSNV_scope_input"] += len(som)
            if len(som) < 2:
                funnel["n_positional_singleton"] += len(som)
                continue
            pre_cap_n = len(som)
            funnel["n_multilocus_pre_cap_groups"] += 1
            funnel["n_multilocus_pre_cap_sSNV"] += pre_cap_n
            # 若 >MAX_SNV, 取最密集連續 MAX_SNV 個 (span 最小) 以利 full-coverage
            if len(som) > MAX_SNV:
                best = min(range(len(som) - MAX_SNV + 1),
                           key=lambda i: som[i + MAX_SNV - 1][0] - som[i][0])
                som = som[best:best + MAX_SNV]
                funnel["n_groups_capped_by_MAX_SNV"] += 1
                funnel["n_sSNV_cap_excluded"] += pre_cap_n - len(som)
            calls, hp, ps = M.per_read_calls(tb, chrom, som)
            tag_diag = M.last_tag_diagnostics()
            for key in (
                "alignment_group_exposures", "sidecar_exact_matches", "sidecar_missing",
                "sidecar_conflicts", "alignment_identity_allele_conflicts",
            ):
                tag_diagnostics[key] += int(tag_diag.get(key, 0))
            raw_hp_total.update(tag_diag.get("raw_HP_counts", {}))
            raw_hp_ps_total.update(tag_diag.get("raw_HP_with_PS_counts", {}))
            alignment_class_total.update(tag_diag.get("alignment_class_counts", {}))
            positions = [s[0] for s in som]
            combo = Counter()
            hpset = Counter()
            # gap#1(2026-07-04):同時記錄 partial-read subcube-groups('X'=未覆蓋位點=cube 自由維度)+
            #   per-column coverage(含 partial),供 group-Steiner 覆蓋條件。full-cov 路徑(combo)不變 → populations 向後相容。
            subread = Counter()
            colcov = {p: [0, 0] for p in positions}  # pos -> [nREF, nALT]
            # L0 per-HP-family split(2026-07-06 使用者定案):樹 per germline 家族分開建。
            #   pooled(combo/subread)保留供對照;下游 layered driver 用 *_by_hp 建 per-family 樹。
            combo_fam = defaultdict(Counter)
            subread_fam = defaultdict(Counter)
            colcov_fam = defaultdict(lambda: {p: [0, 0] for p in positions})
            reads_fam = Counter()
            raw_hp_group = Counter()
            raw_hp_ps_group = Counter()
            phase_sets_group = set()
            phase_set_counts = Counter()
            phase_set_hp_counts = Counter()
            for rn, c in calls.items():
                vec = "".join("A" if c.get(p) == "ALT" else ("R" if c.get(p) == "REF" else "X") for p in positions)
                ncov = len(positions) - vec.count("X")
                if ncov == 0:
                    continue
                fam = germ_family(hp.get(rn))  # L0: read 的 germline 家族
                raw_tag = str(hp.get(rn, "."))
                raw_hp_group[raw_tag] += 1
                if rn in ps:
                    raw_hp_ps_group[raw_tag] += 1
                    phase_set = str(ps[rn])
                    phase_sets_group.add(phase_set)
                    phase_set_counts[phase_set] += 1
                    phase_set_hp_counts[f"{phase_set}|{raw_tag}"] += 1
                reads_fam[fam] += 1
                subread[vec] += 1
                subread_fam[fam][vec] += 1
                for i, p in enumerate(positions):
                    if vec[i] == "R":
                        colcov[p][0] += 1; colcov_fam[fam][p][0] += 1
                    elif vec[i] == "A":
                        colcov[p][1] += 1; colcov_fam[fam][p][1] += 1
                if ncov == len(positions):  # full-cov 路徑不變
                    combo[vec] += 1
                    combo_fam[fam][vec] += 1
                    if rn in hp and "A" in vec:
                        hpset[hp[rn]] += 1
            nfull = sum(combo.values())
            subread_groups = {v: n for v, n in subread.items() if n >= MINREAD}
            # gap#1: 只要有 full-cov population 或 partial subcube-group(≥MINREAD)就保留(原:無 full-cov 即整組丟 → 40.4% 空)
            if nfull < MINREAD and not subread_groups:
                funnel["n_groups_read_unsupported"] += 1
                funnel["n_sSNV_read_unsupported"] += len(som)
                continue
            pops = {k: v for k, v in combo.items() if v >= MINREAD}
            # L0 per-family(≥MINREAD 過濾同 pooled;去空家族)
            pops_by_hp = {f: {k: v for k, v in cc.items() if v >= MINREAD} for f, cc in combo_fam.items()}
            pops_by_hp = {f: d for f, d in pops_by_hp.items() if d}
            subread_by_hp = {f: {k: v for k, v in cc.items() if v >= MINREAD} for f, cc in subread_fam.items()}
            subread_by_hp = {f: d for f, d in subread_by_hp.items() if d}
            colcov_by_hp = {f: {str(p): colcov_fam[f][p] for p in positions} for f in set(pops_by_hp) | set(subread_by_hp)}
            npop = len(pops)
            # populations with >=1 ALT (non-ancestral) = mutation-bearing 群
            npop_alt = sum(1 for k in pops if "A" in k)
            cst = cn_state(cn, chrom, positions[len(positions) // 2])
            same_hp = (len(hpset) == 1)
            rec = {"chrom": chrom, "start": positions[0], "end": positions[-1],
                   "span": positions[-1] - positions[0], "n_sSNV": len(som),
                   "pre_cap_n_sSNV": pre_cap_n, "n_sSNV_cap_excluded": pre_cap_n - len(som),
                   "positions": positions, "n_full_cov_reads": nfull,
                   "n_populations": npop, "n_populations_with_ALT": npop_alt,
                   "populations": dict(sorted(pops.items(), key=lambda x: -x[1])),
                   "subread_groups": dict(sorted(subread_groups.items(), key=lambda x: -x[1])),  # gap#1: partial subcube-groups
                   "col_coverage": {str(p): colcov[p] for p in positions},                        # gap#1: per-locus REF/ALT 覆蓋(含 partial)
                   # L0 per-HP-family(2026-07-06):樹 per germline 家族分開建的底料
                   "populations_by_hp": {f: dict(sorted(d.items(), key=lambda x: -x[1])) for f, d in pops_by_hp.items()},
                   "subread_groups_by_hp": {f: dict(sorted(d.items(), key=lambda x: -x[1])) for f, d in subread_by_hp.items()},
                   "col_coverage_by_hp": colcov_by_hp,
                   "reads_by_hp": dict(reads_fam),
                   "raw_HP_counts": dict(raw_hp_group),
                   "raw_HP_with_PS_counts": dict(raw_hp_ps_group),
                   "n_unique_phase_sets": len(phase_sets_group),
                   "phase_set_counts": dict(phase_set_counts),
                   "phase_set_HP_counts": dict(phase_set_hp_counts),
                   "read_tag_diagnostics": tag_diag,
                   "cn": cst, "cn_data_available": cn is not None,
                   "dominant_hp_count": len(hpset), "same_hp": same_hp}
            groups_out.append(rec)
            phase_set_region_counts[
                "none" if not phase_sets_group else ("one" if len(phase_sets_group) == 1 else "multiple")
            ] += 1
            funnel["n_groups_retained"] += 1
            funnel["n_sSNV_retained"] += len(som)
            agg[f"npop_{min(npop, 5)}"] += 1
            agg[f"nsnv_{min(len(som), 6)}"] += 1
            if npop >= 2:
                agg["groups_ge2_pop"] += 1
            if npop >= 3:
                agg["groups_ge3_pop"] += 1
            if cst in ("loh", "neutral"):
                agg_clean[f"npop_{min(npop, 5)}"] += 1
                if npop >= 2:
                    agg_clean["groups_ge2_pop"] += 1
                if npop >= 3:
                    agg_clean["groups_ge3_pop"] += 1
            if cst == "neutral":
                agg_neutral[f"npop_{min(npop, 5)}"] += 1
                if npop >= 2:
                    agg_neutral["groups_ge2_pop"] += 1
                if npop >= 3:
                    agg_neutral["groups_ge3_pop"] += 1
    tb.close()
    funnel["n_sSNV_accounted"] = (funnel["n_positional_singleton"] + funnel["n_sSNV_cap_excluded"]
                                    + funnel["n_sSNV_read_unsupported"] + funnel["n_sSNV_retained"])
    funnel_out = dict(funnel)
    funnel_out["check_scope_conservation"] = (
        funnel["n_sSNV_scope_input"] == funnel["n_sSNV_accounted"])
    funnel_out["grouping_rule"] = f"adjacent_gap<={M.TIER_R}; total_span_not_bounded"
    tag_diagnostics_out = {
        key: int(tag_diagnostics.get(key, 0))
        for key in (
            "alignment_group_exposures", "sidecar_exact_matches", "sidecar_missing",
            "sidecar_conflicts", "alignment_identity_allele_conflicts",
        )
    }
    out = {"schema_version": "2.0",
           "params": {"MINREAD": MINREAD, "MAX_SNV": MAX_SNV, "TIER_R": M.TIER_R,
                      "MAPQ_MIN": M.MAPQ_MIN, "BASEQ_MIN": M.BASEQ_MIN,
                      "read_tag_source": M.READ_TAG_SIDECAR or "BAM_HP_PS",
                      "read_identity": "QNAME+chrom+start+end+FLAG+BLAKE2b8(CIGAR)",
                      "phase_set_role": "preserved_for_QC_not_used_as_topology_edge",
                      "somatic_backbone": SOMATIC_VCF or "legacy_tp_fp_census",
                      "cn_data_available": cn is not None,
                      "cn_source": CN_SOURCE if cn is not None else "unavailable",
                      "cn_interpretation": ("segments plus unlisted=neutral" if cn is not None else "missing=unavailable")},
           "input_funnel": funnel_out,
           "read_tag_census": {"raw_HP_counts": dict(raw_hp_total),
                               "raw_HP_with_PS_counts": dict(raw_hp_ps_total),
                               "alignment_class_counts": dict(alignment_class_total),
                               "phase_set_region_counts": dict(phase_set_region_counts),
                               "n_regions_with_phase_set_census": len(groups_out),
                               **tag_diagnostics_out},
           "n_groups_analyzed": len(groups_out),
           "aggregate_all": dict(agg), "aggregate_clean_loh_neutral": dict(agg_clean),
           "aggregate_strict_neutral": dict(agg_neutral),
           "groups": groups_out}
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(out, handle, ensure_ascii=False)
    print(f"DONE chroms={chroms} groups={len(groups_out)} ge2pop={agg['groups_ge2_pop']} "
          f"ge3pop={agg['groups_ge3_pop']} -> {out_path}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","), sys.argv[2])
