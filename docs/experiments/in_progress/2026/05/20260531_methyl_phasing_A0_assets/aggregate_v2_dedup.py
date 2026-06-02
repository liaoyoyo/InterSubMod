#!/usr/bin/env python3
"""v2 去重疊-window 彙總（修取樣 bug：±2kb window 相鄰<4kb 共享同批 read = 重複計數）。
全統計落 seqc2_v2_dedup_aggregate.json，報告 Read 後才引用。"""
import json, glob, os
import numpy as np
from scipy.stats import mannwhitneyu, wilcoxon
from collections import Counter

AD = os.path.dirname(os.path.abspath(__file__))
OUT = AD + "/seqc2_v2_dedup_aggregate.json"

def med(xs):
    xs = [x for x in xs if x is not None]
    return round(float(np.median(xs)), 4) if xs else None
def mean(xs):
    xs = [x for x in xs if x is not None]
    return round(float(np.mean(xs)), 4) if xs else None

rows = []
chroms = []
for fp in sorted(glob.glob(AD + "/seqc2_cn_methyl_v2_chr*.json")):
    d = json.load(open(fp)); chroms.append(d["chrom"]); rows.extend(d["regions"])

def dedup(rs, min_gap=4000):
    by = {}
    for r in rs: by.setdefault(r["chrom"], []).append(r)
    kept = []
    for c, lst in by.items():
        lst = sorted(lst, key=lambda r: r["pos"]); last = -1e18
        for r in lst:
            if r["pos"] - last >= min_gap:
                kept.append(r); last = r["pos"]
    return kept

ded = dedup(rows)
out = {"n_chroms": len(chroms), "chroms": sorted(chroms),
       "n_regions_raw": len(rows), "n_regions_dedup": len(ded),
       "dedup_note": "±2kb window 相鄰<4kb 視為重疊共享 read，貪婪保留間距>=4kb"}

def sep(r):
    return (r["anchor_auc_raw"] is not None and r["shuffle_p95_raw"] is not None
            and r["anchor_auc_raw"] - r["shuffle_p95_raw"] > 0.05)

# by status (dedup)
by_status = {}
for s in ("gain", "loss", "loh", "neutral"):
    sub = [r for r in ded if r["seqc2_status"] == s]
    sepsub = [r for r in sub if sep(r)]
    by_status[s] = {
        "n": len(sub),
        "auc_raw_median": med([r["anchor_auc_raw"] for r in sub]),
        "auc_excl_median": med([r["anchor_auc_cpgsnp_excl"] for r in sub]),
        "n_separable": len(sepsub),
        "frac_separable": round(len(sepsub)/len(sub), 3) if sub else None,
        "by_chr_n": dict(Counter(r["chrom"] for r in sub)),
        "by_chr_sep": dict(Counter(r["chrom"] for r in sepsub)),
    }
out["by_status_dedup"] = by_status

# (a) CpG-SNP 排除（dedup）
drops = [r["auc_drop_after_excl"] for r in ded if r["auc_drop_after_excl"] is not None]
paired = [(r["anchor_auc_raw"], r["anchor_auc_cpgsnp_excl"]) for r in ded
          if r["anchor_auc_raw"] is not None and r["anchor_auc_cpgsnp_excl"] is not None]
a = {
    "n": len(drops), "drop_median": med(drops), "drop_mean": mean(drops),
    "drop_max": round(float(np.max(drops)), 4) if drops else None,
    "frac_drop_gt_0.05": round(float(np.mean([d > 0.05 for d in drops])), 3) if drops else None,
    "tainted_cpg_median": med([r["n_cpg_tainted"] for r in ded]),
}
if len(paired) >= 8:
    try:
        _, p = wilcoxon([x[0] for x in paired], [x[1] for x in paired])
        a["wilcoxon_raw_vs_excl_p"] = round(float(p), 5)
    except Exception:
        a["wilcoxon_raw_vs_excl_p"] = None
out["a_cpgsnp_exclusion_dedup"] = a

# (c) LOH 個案（dedup）：分得開的 LOH + 異常檢查
loh = [r for r in ded if r["seqc2_status"] == "loh"]
loh_sep = [r for r in loh if sep(r)]
loh_not = [r for r in loh if not sep(r)]
loh_sorted = sorted(loh, key=lambda r: -(r["anchor_auc_raw"] or 0))
# LOH vs neutral Mann-Whitney (dedup)
neu = [r["anchor_auc_raw"] for r in ded if r["seqc2_status"] == "neutral"]
lohp = None
if len(loh) >= 5 and len(neu) >= 5:
    try:
        _, lohp = mannwhitneyu([r["anchor_auc_raw"] for r in loh], neu, alternative="two-sided")
        lohp = round(float(lohp), 5)
    except Exception:
        pass
out["c_loh_dedup"] = {
    "n_loh": len(loh), "n_separable": len(loh_sep),
    "frac_separable": round(len(loh_sep)/len(loh), 3) if loh else None,
    "loh_auc_median": med([r["anchor_auc_raw"] for r in loh]),
    "loh_vs_neutral_mannwhitney_p": lohp,
    "by_chr_n": dict(Counter(r["chrom"] for r in loh)),
    "by_chr_sep": dict(Counter(r["chrom"] for r in loh_sep)),
    # 異常：分得開 vs 分不開 LOH 的 depth/bic
    "sep_depth_median": med([r["mean_depth"] for r in loh_sep]),
    "notsep_depth_median": med([r["mean_depth"] for r in loh_not]),
    "sep_bic_median": med([r["gmm_bic_diff"] for r in loh_sep]),
    "notsep_bic_median": med([r["gmm_bic_diff"] for r in loh_not]),
    "separable_top10": [
        {"region": f"{r['chrom']}:{r['pos']}", "auc_raw": r["anchor_auc_raw"],
         "auc_excl": r["anchor_auc_cpgsnp_excl"], "shuffle_p95": r["shuffle_p95_raw"],
         "depth": r["mean_depth"], "n_cpg": r["n_cpg_raw"], "bic": r["gmm_bic_diff"],
         "n_hp1": r["n_hp1"], "n_hp2": r["n_hp2"]}
        for r in loh_sorted[:10]
    ],
}
# v1 對照（chr8 兩次重現性）
v1 = []
for fp in glob.glob(AD + "/seqc2_cn_methyl_chr*.json"):
    if "v2" in fp: continue
    v1.extend(json.load(open(fp))["regions"])
v1_chr8_loh = [r["anchor_auc"] for r in v1 if r["chrom"] == "chr8" and r["seqc2_status"] == "loh"]
v2_chr8_loh = [r["anchor_auc_raw"] for r in ded if r["chrom"] == "chr8" and r["seqc2_status"] == "loh"]
out["chr8_loh_reproducibility"] = {
    "v1_n": len(v1_chr8_loh), "v1_median": med(v1_chr8_loh),
    "v2_n": len(v2_chr8_loh), "v2_median": med(v2_chr8_loh),
}

json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=2)
print(json.dumps({k: v for k, v in out.items() if "top10" not in str(k)}, ensure_ascii=False, indent=2)[:2000])
print(f"\n[dedup-agg] -> {OUT}")
