#!/usr/bin/env python3
"""
79 - (A) Why high-Δβ loci cluster as FP: condition→FP statistics.
     (B) Consensus high-methylation-difference loci (high by ALL methods).
HCC1395 (only sample with usable FP=627 + ISM built-in LOH columns).
§13 layer-A. Output: genome_survey_v2/cn_confound/cross_sample/condition_fp_consensus.json
"""
import os, csv, json
import numpy as np
from scipy.stats import fisher_exact

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
CS = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample")
OUT = f"{CS}/condition_fp_consensus.json"
TESTS = ["ClusterPermanovaP", "LabelHPPermanovaP", "LabelAllelePermanovaP", "HPFineP"]


def num(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def load(s, c):
    rows = list(csv.DictReader(open(f"{EX}/{s}_{c}/significance_summary.csv")))
    for r in rows:
        r["_db"] = abs(num(r, "HPMergedDelta") or 0.0)
        r["_cv"] = num(r, "CramersV") or 0.0
        r["_nr"] = num(r, "NumReads") or 0.0
        r["_ntests"] = sum(1 for t in TESTS if (num(r, t) is not None and num(r, t) <= 0.05))
        r["_loh"] = tb(r, "Potential_LOH") or tb(r, "LOH_Bed_Overlap")
        r["_n1s"] = num(r, "HPFineN_HP1S") or 0
        r["_n2s"] = num(r, "HPFineN_HP2S") or 0
        r["_subhap_small"] = (max(r["_n1s"], r["_n2s"]) < 10)
        r["_nbase"] = num(r, "NormalBaseline_Mean")
        r["_extreme_base"] = (r["_nbase"] is not None and (r["_nbase"] <= 0.15 or r["_nbase"] >= 0.85))
        r["_sig"] = tb(r, "Significant")
    return rows


tp = load("HCC1395", "tp")
fp = load("HCC1395", "fp")


def fisher_cond(name, hi_tp, hi_fp):
    """hi_tp/hi_fp = lists; count condition True vs False, OR for FP enrichment."""
    a = sum(1 for r in hi_tp if r["_cond"])
    b = len(hi_tp) - a
    c = sum(1 for r in hi_fp if r["_cond"])
    d = len(hi_fp) - c
    # OR for FP given condition (rows=condition, cols=FP/TP)
    try:
        orv, p = fisher_exact([[c, d], [a, b]])  # FP-vs-TP for cond
    except Exception:
        orv, p = float("nan"), float("nan")
    tp_rate = a / len(hi_tp) if hi_tp else None
    fp_rate = c / len(hi_fp) if hi_fp else None
    return dict(condition=name, tp_with=a, tp_n=len(hi_tp), fp_with=c, fp_n=len(hi_fp),
                tp_rate=round(tp_rate, 4) if tp_rate is not None else None,
                fp_rate=round(fp_rate, 4) if fp_rate is not None else None,
                fp_vs_tp_OR=round(orv, 3), p=round(p, 6))


# ---- (A) among HIGH-Δβ loci (|Δβ|>=0.1), which conditions distinguish FP? ----
hi_tp = [r for r in tp if r["_db"] >= 0.10]
hi_fp = [r for r in fp if r["_db"] >= 0.10]
condA = []
for name, fn in [
    ("no_clustering (CramersV<0.1)", lambda r: r["_cv"] < 0.1),
    ("LOH region", lambda r: r["_loh"]),
    ("small somatic subhap (<10 reads)", lambda r: r["_subhap_small"]),
    ("extreme normal baseline (<=0.15 or >=0.85)", lambda r: r["_extreme_base"]),
    ("low coverage (<30 reads)", lambda r: r["_nr"] < 30),
]:
    for r in hi_tp + hi_fp:
        r["_cond"] = fn(r)
    condA.append(fisher_cond(name, hi_tp, hi_fp))

# combined: high-Δβ AND no-clustering AND LOH -> FP rate
def combo(rows):
    n = sum(1 for r in rows if r["_db"] >= 0.10 and r["_cv"] < 0.1 and r["_loh"])
    return n
combo_tp = combo(tp); combo_fp = combo(fp)

# ---- (B) consensus high-diff loci: BOTH ISM-clustering AND Δβ strong AND multi-test ----
consensus = []
for r in tp:
    if r["_cv"] >= 0.1 and r["_db"] >= 0.10 and r["_ntests"] >= 3 and r["_sig"]:
        consensus.append(dict(locus=f"{r.get('Chr')}:{r.get('Pos')}",
                              cramers=round(r["_cv"], 3), dbeta=round(num(r, "HPMergedDelta"), 3),
                              n_tests=r["_ntests"], num_reads=int(r["_nr"]),
                              loh=r["_loh"], significant=r["_sig"]))
consensus.sort(key=lambda x: -(abs(x["dbeta"]) + x["cramers"]))

out = dict(
    meta=dict(script="79_condition_fp_and_consensus.py", sample="HCC1395",
              note="(A) among high-Δβ (|Δβ|>=0.1) loci, Fisher FP-vs-TP OR per condition; "
                   "(B) consensus = CramersV>=0.1 AND |Δβ|>=0.1 AND >=3/4 tests AND Significant.",
              hi_db_tp=len(hi_tp), hi_db_fp=len(hi_fp)),
    A_condition_fp=condA,
    A_combo_highDb_noClust_LOH=dict(tp=combo_tp, fp=combo_fp,
        tp_rate=round(combo_tp / len(tp), 4), fp_rate=round(combo_fp / len(fp), 4),
        interpretation="高Δβ + 無分群 + LOH 三條件齊備時的 TP/FP 率"),
    B_consensus=dict(n=len(consensus), top=consensus[:30]),
)
with open(OUT, "w") as f:
    json.dump(out, f, indent=2, default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
print(f"[79] wrote {OUT}")
print(f"\n=== (A) high-Δβ (TP n={len(hi_tp)} / FP n={len(hi_fp)}) 條件→FP ===")
for c in condA:
    print(f"  {c['condition']:42s}: TP {c['tp_rate']} vs FP {c['fp_rate']}  FP/TP OR={c['fp_vs_tp_OR']} p={c['p']}")
print(f"\n  combo 高Δβ+無分群+LOH: TP {combo_tp}/{len(tp)} ({combo_tp/len(tp)*100:.2f}%)  FP {combo_fp}/{len(fp)} ({combo_fp/len(fp)*100:.2f}%)")
print(f"\n=== (B) consensus 高甲基差異(高 CramersV+高Δβ+多檢定) n={len(consensus)} top 8 ===")
for x in consensus[:8]:
    print(f"  {x['locus']:18s} CramersV={x['cramers']} Δβ={x['dbeta']:+.3f} tests={x['n_tests']} reads={x['num_reads']} LOH={x['loh']}")
