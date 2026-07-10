#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""權威六層 funnel 分母字典 + tree-level 誠實重分類（HCC1395，2026-07-10 P0 修正）。

背景：2026-07-10 peer review 揪出 3 個 P0（獨立重算精確吻合）：
  P0-1 舊 census 的「isolated 88,358」把 out-of-scope(chrX) + densest-8 截斷丟棄 當「生物單點」；
       真 autosomal positional singleton 只有 8,447。
  P0-2 「determined 7,123」含 1,739 root-only(純 REF 家族、樹只有 ROOT 無 somatic edge)；
       mutation-bearing determined = 5,384。
  P0-3 「multi-HP 4,443(56%)」被 REF-only germline 家族膨脹；兩 germline 都 mutation-bearing = 2,765(34.9%)。

本腳本 = 單一真值來源：從 somatic_pass.vcf.gz（六層 funnel）+ layered JSON（tree-level 重分類）
重算全部 → 落 funnel_census_HCC1395.json。只讀不寫源碼。HTML/PI 由此 JSON 注入（不手打）。
scope = chr1–22（使用者定案 2026-07-10）；chrX 標 out-of-scope。
用法：python3 funnel_census_HCC1395.py
"""
import json, os, re, gzip, glob
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.normpath(os.path.join(HERE, "..", "data"))
PILOT = os.path.normpath(os.path.join(HERE, "..", "..", "20260618_subcluster_pilot"))
SP = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/somatic_pass.vcf.gz"

TIER_R = 50000   # linkage：相鄰 gap>50kb 切單連通群（sm_linkage_genomewide.py:32）
MAX_SNV = 8      # densest-8 截斷（sm_multilocus_combinations.py:19）
AUTO = re.compile(r"^chr([1-9]|1[0-9]|2[0-2])$")


def pct(n, d): return round(n / d * 100, 1) if d else 0.0


# ── 六層 funnel（sSNV 層，純從 VCF 位置 + mlhp retained 重建）──
def compute_funnel():
    by_chrom = defaultdict(list); total = 0; xy = 0
    with gzip.open(SP, "rt") as f:
        for ln in f:
            if ln.startswith("#"): continue
            p = ln.split("\t")
            if len(p[3]) != 1 or len(p[4].strip()) != 1: continue  # SNV only
            total += 1
            c = p[0]
            if c in ("chrX", "chrY"): xy += 1
            elif AUTO.match(c): by_chrom[c].append(int(p[1]))
    auto = sum(len(v) for v in by_chrom.values())
    singleton = 0; in_group = 0; capped_positional = 0; ngroups = 0
    for c, ps in by_chrom.items():
        ps.sort(); grp = [ps[0]]
        def flush(g):
            nonlocal singleton, in_group, capped_positional, ngroups
            if len(g) == 1:
                singleton += 1
            else:
                ngroups += 1; in_group += len(g); capped_positional += min(len(g), MAX_SNV)
        for x in ps[1:]:
            if x - grp[-1] > TIER_R: flush(grp); grp = [x]
            else: grp.append(x)
        flush(grp)
    # mlhp 實際 retained（讀已分析區的 n_sSNV 和，已含 densest-8 截斷 + read-support 過濾）
    retained = 0; nreg = 0
    for fp in sorted(glob.glob(f"{PILOT}/mlhp_part_*.json")):
        for g in json.load(open(fp)).get("groups", []):
            retained += g.get("n_sSNV", 0); nreg += 1
    cap_excluded = in_group - capped_positional          # densest-8 截斷丟棄
    read_unsupported = capped_positional - retained       # 通過位置+cap 但 read 不足被丟
    return {
        "L1_all_pass_universe": total, "note_universe": "ClairS paired PASS SNV 全基因組=操作骨幹，非生物真值",
        "L2_out_of_scope_chrXY": xy, "note_scope": "chrX/chrY 未做 read 分析（使用者定案 scope=chr1–22）",
        "autosomal_chr1_22": auto,
        "L3_positional_singleton": singleton, "note_singleton": "相鄰 gap>50kb 兩側=真孤立單點（無共現可建樹）",
        "in_multilocus_group": in_group, "n_positional_groups": ngroups,
        "L4_cap_excluded_densest8": cap_excluded, "note_cap": "大群只取最密 8 個，其餘 sSNV 未進 solver（精確窮舉工程界）",
        "L5_read_unsupported": read_unsupported, "note_unsupp": "位置成群但 full-cov/subread read<MINREAD=3 被丟",
        "L6_retained_sSNV": retained, "retained_regions": nreg,
        "note_retained": "chr1–22、經 densest-8 與 read-support 後真正送入 group-Steiner solver 的 sSNV",
        "funnel_check": {
            "sum_layers": xy + singleton + cap_excluded + read_unsupported + retained,
            "equals_universe": (xy + singleton + cap_excluded + read_unsupported + retained) == total,
        },
    }


# ── tree-level 誠實重分類（layered JSON）──
def has_somatic_edge(u):
    for t in (u.get("trees") or []):
        if t.get("edges"): return True
    return False


def compute_tree_level():
    lay = json.load(open(f"{PILOT}/layered_reconstruction_HCC1395.json"))
    detail = lay["detail"]
    lineage = [u for u in detail if u.get("is_lineage")]
    nl = len(lineage)
    det = [u for u in lineage if u["L1_class"] == "determined"]
    root_only = [u for u in det if not has_somatic_edge(u)]       # 純 REF 家族=樹只有 ROOT
    det_mut = len(det) - len(root_only)
    amb = sum(1 for u in lineage if "ambiguous" in u["L1_class"])
    cap = sum(1 for u in lineage if "capped" in u["L1_class"])
    rec = sum(1 for u in lineage if u["L1_class"] == "recurrence_required")
    nonroot = nl - len(root_only)
    # multi-HP：region 兩 germline 家族 both mutation-bearing vs 含 REF-only
    reg_allfam = defaultdict(set); reg_mutfam = defaultdict(set)
    for u in detail:
        reg_allfam[u["region"]].add(u["family"])
        if has_somatic_edge(u): reg_mutfam[u["region"]].add(u["family"])
    W = len(reg_allfam)
    span_two_germ = sum(1 for r, fs in reg_allfam.items() if len(fs & {"1", "2"}) == 2)
    both_germ_mut = sum(1 for r in reg_allfam if len(reg_mutfam.get(r, set()) & {"1", "2"}) == 2)
    hp3_regions = sum(1 for r, fs in reg_allfam.items() if "3" in fs)
    return {
        "n_lineage_units": nl,
        "determined_all": len(det),
        "L2a_root_only_reference_family": len(root_only),
        "note_root_only": "純 REF 家族：樹只有 ROOT、無 somatic edge → 不算重建出 subclone 樹（trivially determined）",
        "determined_mutation_bearing": det_mut,
        "det_all_pct_of_units": pct(len(det), nl),
        "det_mut_pct_of_nonroot": pct(det_mut, nonroot),
        "note_det_pct": f"誠實 determined = {det_mut}/{nonroot}(排除 root-only 家族分母) = {pct(det_mut, nonroot)}%（舊 {pct(len(det), nl)}% 含 root-only）",
        "ambiguous": amb, "capped": cap, "recurrence": rec,
        "W_regions": W,
        "multi_hp_span_two_germline_incl_refonly": span_two_germ, "multi_hp_span_pct": pct(span_two_germ, W),
        "multi_hp_both_germline_mutation_bearing": both_germ_mut, "multi_hp_mut_pct": pct(both_germ_mut, W),
        "note_multi_hp": f"雙 germline lineage 誠實值 = {both_germ_mut}/{W}={pct(both_germ_mut, W)}%（舊 {span_two_germ}={pct(span_two_germ, W)}% 被 REF-only germline 家族膨脹）",
        "hp3_regions": hp3_regions,
        "note_hp3": "HP3=整合 somatic 的 ambiguous haplotype（KB）；標 H3?、維持計數但非確認第三獨立 germline lineage（使用者定案：維持但加註）",
    }


out = {
    "meta": {"sample": "HCC1395", "date": "2026-07-10", "scope": "chr1–22（chrX out-of-scope）",
             "backbone": "ClairS paired PASS = somatic_pass.vcf.gz（操作骨幹非生物真值）",
             "claim": "regional mutation-state tree 可驗證；生物 subclone 可如其他軟體般推論（inferable，非 confirmed）；單樣本 ⭐3",
             "provenance": "六層 funnel 由 somatic_pass.vcf.gz 位置重建；tree-level 由 layered_reconstruction_HCC1395.json 重分類"},
    "funnel_sSNV": compute_funnel(),
    "tree_level": compute_tree_level(),
}
json.dump(out, open(f"{D}/funnel_census_HCC1395.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(out, ensure_ascii=False, indent=1))
