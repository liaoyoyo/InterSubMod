#!/usr/bin/env python3
"""彙總 v2（含 CpG-SNP 排除 a + LOH 個案 c）。全統計落 JSON，報告 Read 後才引用。"""
import json, glob, os
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu, wilcoxon

AD = os.path.dirname(os.path.abspath(__file__))
OUT = AD + "/seqc2_v2_aggregate.json"

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

out = {"n_chroms": len(chroms), "chroms": sorted(chroms), "n_regions": len(rows)}

# ===== (a) CpG-SNP 排除效果 =====
paired = [(r["anchor_auc_raw"], r["anchor_auc_cpgsnp_excl"]) for r in rows
          if r["anchor_auc_raw"] is not None and r["anchor_auc_cpgsnp_excl"] is not None]
drops = [r["auc_drop_after_excl"] for r in rows if r["auc_drop_after_excl"] is not None]
taint = [r["n_cpg_tainted"] for r in rows]
cpg_raw = [r["n_cpg_raw"] for r in rows]
a = {
    "n_with_both": len(paired),
    "auc_raw_median": med([x[0] for x in paired]),
    "auc_excl_median": med([x[1] for x in paired]),
    "drop_median": med(drops), "drop_mean": mean(drops),
    "drop_max": round(float(np.max(drops)), 4) if drops else None,
    "drop_p95": round(float(np.percentile(drops, 95)), 4) if drops else None,
    "frac_drop_gt_0.05": round(float(np.mean([d > 0.05 for d in drops])), 3) if drops else None,
    "tainted_cpg_median": med(taint),
    "tainted_frac_of_cpg_median": round(float(np.median([t/c for t, c in zip(taint, cpg_raw) if c > 0])), 4) if cpg_raw else None,
}
# 配對 Wilcoxon: raw vs excl 是否顯著不同
if len(paired) >= 8:
    try:
        w, p = wilcoxon([x[0] for x in paired], [x[1] for x in paired])
        a["wilcoxon_raw_vs_excl_p"] = round(float(p), 5)
    except Exception:
        a["wilcoxon_raw_vs_excl_p"] = None
out["a_cpgsnp_exclusion"] = a

# ===== Q1/Q3 重算（v2 用 raw auc）=====
by_status = {}
for s in ("gain", "loss", "loh", "neutral"):
    sub = [r for r in rows if r["seqc2_status"] == s]
    by_status[s] = {
        "n": len(sub),
        "auc_raw_median": med([r["anchor_auc_raw"] for r in sub]),
        "auc_excl_median": med([r["anchor_auc_cpgsnp_excl"] for r in sub]),
    }
out["by_status"] = by_status

# ===== (c) LOH 個案：哪些 LOH 反而分得開 + 異常 =====
loh = [r for r in rows if r["seqc2_status"] == "loh"]
# 排序：anchor_auc_raw 高 = 反而分得開
loh_sorted = sorted(loh, key=lambda r: -(r["anchor_auc_raw"] or 0))
# 定「分得開」= auc_raw > shuffle_p95 + 0.05 (顯著超 null)
def separable(r):
    return (r["anchor_auc_raw"] is not None and r["shuffle_p95_raw"] is not None
            and r["anchor_auc_raw"] - r["shuffle_p95_raw"] > 0.05)
loh_sep = [r for r in loh if separable(r)]
loh_notsep = [r for r in loh if not separable(r)]
out["c_loh"] = {
    "n_loh": len(loh),
    "n_separable": len(loh_sep),
    "frac_separable": round(len(loh_sep)/len(loh), 3) if loh else None,
    "loh_auc_median": med([r["anchor_auc_raw"] for r in loh]),
    # 異常檢查：分得開的 LOH 是否有高 depth(疑多拷貝混LOH) / 高 tainted / 多 CpG
    "separable_loh_top": [
        {"region": f"{r['chrom']}:{r['pos']}", "auc_raw": r["anchor_auc_raw"],
         "auc_excl": r["anchor_auc_cpgsnp_excl"], "shuffle_p95": r["shuffle_p95_raw"],
         "depth": r["mean_depth"], "n_cpg": r["n_cpg_raw"], "tainted": r["n_cpg_tainted"],
         "bic": r["gmm_bic_diff"], "n_hp1": r["n_hp1"], "n_hp2": r["n_hp2"]}
        for r in loh_sorted[:15]
    ],
    # 異常統計：分得開 vs 分不開 LOH 的 depth/cpg 對比
    "separable_depth_median": med([r["mean_depth"] for r in loh_sep]),
    "notsep_depth_median": med([r["mean_depth"] for r in loh_notsep]),
    "separable_bic_median": med([r["gmm_bic_diff"] for r in loh_sep]),
    "notsep_bic_median": med([r["gmm_bic_diff"] for r in loh_notsep]),
}

json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=2)
print(json.dumps({k: v for k, v in out.items() if k != "c_loh"}, ensure_ascii=False, indent=2))
print("\n--- c_loh ---")
print(json.dumps(out["c_loh"], ensure_ascii=False, indent=2)[:1500])
print(f"\n[aggregate-v2] -> {OUT}")
