#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""7 樣本六層 funnel + tree-level 誠實重分類（2026-07-10 P0 修正，同 HCC1395 口徑擴 7 樣本）。

不需重讀 BAM：全 7 樣本 layered_reconstruction_{S}.json 已於 2026-07-09 用新骨幹（ClairS PASS）跑完。
每樣本：
  ① 六層 funnel（somatic_pass.vcf.gz 位置重建 50kb 單連通群 + densest-8 截斷；retained=layered 區 n_sSNV 和）
  ② tree-level 重分類（determined 拆 root-only 純 REF vs mutation-bearing；multi-HP both-germline-mutation-bearing）
  ③ P0-4：無 CN bed 的 6 樣本，cn='neutral' 是預設非量測 → 報告層標 'unavailable'（recurrence CN 裁決同）。
scope=chr1–22（chrX out-of-scope）。只讀不寫源碼。輸出 funnel_census_7samples.json。
"""
import json, os, re, gzip
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.normpath(os.path.join(HERE, "..", "data"))
C = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical"
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
TIER_R = 50000; MAX_SNV = 8
AUTO = re.compile(r"^chr([1-9]|1[0-9]|2[0-2])$")

# 只有 HCC1395 有 SEQC2 CN bed；其餘 6 樣本 cn=unknown（run_layered_7samples_newbb.sh CNBED=""）
# 每樣本：(name, layered JSON, work-dir〔含 mlhp_part_*.json〕, has_cn)
SAMPLES = [
    ("HCC1395",        f"{PILOT}/layered_reconstruction_HCC1395.json",              PILOT,                     True),
    ("COLO829",        f"{MSROOT}/COLO829/layered_reconstruction_COLO829.json",     f"{MSROOT}/COLO829",       False),
    ("H1437",          f"{MSROOT}/H1437/layered_reconstruction_H1437.json",         f"{MSROOT}/H1437",         False),
    ("H2009",          f"{MSROOT}/H2009/layered_reconstruction_H2009.json",         f"{MSROOT}/H2009",         False),
    ("HCC1395_DORADO", f"{MSROOT}/HCC1395_DORADO/layered_reconstruction_HCC1395_DORADO.json", f"{MSROOT}/HCC1395_DORADO", False),
    ("HCC1937",        f"{MSROOT}/HCC1937/layered_reconstruction_HCC1937.json",     f"{MSROOT}/HCC1937",       False),
    ("HCC1954",        f"{MSROOT}/HCC1954/layered_reconstruction_HCC1954.json",     f"{MSROOT}/HCC1954",       False),
]


def pct(n, d): return round(n / d * 100, 1) if d else 0.0


def somatic_pass_path(s):
    import glob
    lp = glob.glob(f"{C}/{s}/paired_full/*complete_matrix/longphase_s")
    return f"{lp[0]}/somatic_pass.vcf.gz" if lp else None


def has_somatic_edge(u):
    for t in (u.get("trees") or []):
        if t.get("edges"): return True
    return False


def mlhp_retained(wd):
    """retained sSNV = mlhp 各區 n_sSNV 和（已含 densest-8 截斷 + read-support；與單樣本 census 同源）。"""
    import glob
    retained = 0; nreg = 0
    for fp in sorted(glob.glob(f"{wd}/mlhp_part_*.json")):
        for g in json.load(open(fp)).get("groups", []):
            retained += g.get("n_sSNV", 0); nreg += 1
    return retained, nreg


def compute_funnel(s, wd):
    sp = somatic_pass_path(s)
    by_chrom = defaultdict(list); total = 0; xy = 0
    with gzip.open(sp, "rt") as f:
        for ln in f:
            if ln.startswith("#"): continue
            p = ln.split("\t")
            if len(p[3]) != 1 or len(p[4].strip()) != 1: continue
            total += 1; c = p[0]
            if c in ("chrX", "chrY"): xy += 1
            elif AUTO.match(c): by_chrom[c].append(int(p[1]))
    auto = sum(len(v) for v in by_chrom.values())
    singleton = 0; in_group = 0; capped_positional = 0; ngroups = 0
    for c, ps in by_chrom.items():
        ps.sort(); grp = [ps[0]]
        def flush(g):
            nonlocal singleton, in_group, capped_positional, ngroups
            if len(g) == 1: singleton += 1
            else:
                ngroups += 1; in_group += len(g); capped_positional += min(len(g), MAX_SNV)
        for x in ps[1:]:
            if x - grp[-1] > TIER_R: flush(grp); grp = [x]
            else: grp.append(x)
        flush(grp)
    retained, nreg = mlhp_retained(wd)
    cap_excluded = in_group - capped_positional
    read_unsupported = capped_positional - retained
    return {
        "L1_all_pass_universe": total, "L2_out_of_scope_chrXY": xy, "autosomal_chr1_22": auto,
        "L3_positional_singleton": singleton, "in_multilocus_group": in_group, "n_positional_groups": ngroups,
        "L4_cap_excluded_densest8": cap_excluded, "L5_read_unsupported": read_unsupported,
        "L6_retained_sSNV": retained, "retained_regions": nreg,
        "singleton_pct_of_universe": pct(singleton, total),
        "check_sum": xy + singleton + cap_excluded + read_unsupported + retained,
        "check_ok": (xy + singleton + cap_excluded + read_unsupported + retained) == total,
    }


def compute_tree_level(detail, has_cn):
    lineage = [u for u in detail if u.get("is_lineage")]
    nl = len(lineage)
    det = [u for u in lineage if u["L1_class"] == "determined"]
    root_only = [u for u in det if not has_somatic_edge(u)]
    det_mut = len(det) - len(root_only)
    amb = sum(1 for u in lineage if "ambiguous" in u["L1_class"])
    cap = sum(1 for u in lineage if "capped" in u["L1_class"])
    rec = sum(1 for u in lineage if u["L1_class"] == "recurrence_required")
    nonroot = nl - len(root_only)
    reg_allfam = defaultdict(set); reg_mutfam = defaultdict(set)
    for u in detail:
        reg_allfam[u["region"]].add(u["family"])
        if has_somatic_edge(u): reg_mutfam[u["region"]].add(u["family"])
    W = len(reg_allfam)
    span_two = sum(1 for r, fs in reg_allfam.items() if len(fs & {"1", "2"}) == 2)
    both_mut = sum(1 for r in reg_allfam if len(reg_mutfam.get(r, set()) & {"1", "2"}) == 2)
    hp3 = sum(1 for r, fs in reg_allfam.items() if "3" in fs)
    return {
        "n_lineage_units": nl, "W_regions": W,
        "determined_all": len(det), "root_only_reference_family": len(root_only),
        "determined_mutation_bearing": det_mut,
        "det_all_pct": pct(len(det), nl), "det_mut_pct_of_nonroot": pct(det_mut, nonroot),
        "ambiguous": amb, "capped": cap, "recurrence": rec,
        "multi_hp_span_two_germline": span_two, "multi_hp_span_pct": pct(span_two, W),
        "multi_hp_both_mutation_bearing": both_mut, "multi_hp_mut_pct": pct(both_mut, W),
        "hp3_regions": hp3,
        # P0-4：無 CN → recurrence CN 裁決 unavailable（非 neutral 預設）
        "recurrence_cn_verdict": ("有 CN 可判（artifact/LOH）" if has_cn else "unavailable（無 CN bed·非 neutral 預設）"),
        "cn_status": ("SEQC2 CN bed" if has_cn else "unavailable（cn=unknown）"),
    }


out = {"meta": {"date": "2026-07-10", "scope": "chr1–22（chrX out-of-scope）",
                "note": "P0 修正同口徑擴 7 樣本；不重讀 BAM（layered 07-09 已跑）；P0-4 無CN→unavailable",
                "cn_available": ["HCC1395"], "cn_unavailable": [s for s, _, _, h in SAMPLES if not h]},
       "samples": {}}
for s, lay_path, wd, has_cn in SAMPLES:
    detail = json.load(open(lay_path))["detail"]
    out["samples"][s] = {"has_cn": has_cn, "funnel": compute_funnel(s, wd),
                         "tree_level": compute_tree_level(detail, has_cn)}
json.dump(out, open(f"{D}/funnel_census_7samples.json", "w"), ensure_ascii=False, indent=1)

# 精簡列印
print(f"{'sample':16} {'universe':>9} {'chrX':>7} {'singl':>7} {'retain':>7} {'units':>6} {'detAll':>7} {'detMut':>7} {'mut%':>5} {'bothHP':>7} {'CN':>12}")
for s, _, _, has_cn in SAMPLES:
    f = out["samples"][s]["funnel"]; t = out["samples"][s]["tree_level"]
    print(f"{s:16} {f['L1_all_pass_universe']:>9} {f['L2_out_of_scope_chrXY']:>7} {f['L3_positional_singleton']:>7} "
          f"{f['L6_retained_sSNV']:>7} {t['n_lineage_units']:>6} {t['determined_all']:>7} {t['determined_mutation_bearing']:>7} "
          f"{t['det_mut_pct_of_nonroot']:>5} {t['multi_hp_both_mutation_bearing']:>7} {'SEQC2' if has_cn else 'unavail':>12}  funnel_ok={f['check_ok']}")
