#!/usr/bin/env python3
"""[弱結構保留分類] 把 phylo-v4.1 genome-wide 每位點記錄分成「互斥 6 態」決策流,
目的: 不把弱結構丟進 no_structure 黑洞,而是顯式分桶讓人可觀察確認。
讀 phylo_cpp_wg_records.json(C++ native 真值,34736) → 每位點走 ordered first-match 決策流 → 互斥分類。
輸出 counts(TP/FP, sum-check) + 每位點分類(dashboard 用)。純讀已落檔記錄,不重算 binary。"""
import json
import os
from collections import Counter, defaultdict

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
REC = f"{A}/phylo_cpp_wg_records.json"

# 6-態互斥決策流(ordered, first-match-wins)。強→弱→無。
#   S1 CONFIRMED_MULTI_ALIGNED : coarse>=2 & stable & aligned   — 最高信心多群,顯著+穩定+對齊a-priori
#   S2 CONFIRMED_MULTI_UNALIGNED: coarse>=2 & stable & !aligned  — 顯著穩定多群但不對齊(cis-ASM候選/無標籤對應)
#   S3 UNSTABLE_MULTI          : coarse>=2 & unstable           — coarse報多群但seed分歧(modal<0.7)=borderline,留候選
#   S4 WEAK_FINE_ONLY          : coarse<2 & fine>=2  — ⚠ 僅寬閘(null90,alpha0.10)過、嚴閘(null95,alpha0.05)擋。
#                                 TP 5.97% vs FP 5.19% = ratio 1.15 (chance-level,與噪音統計不可分) → 非「保留的真弱結構」,
#                                 是 sub-threshold candidates pending validation。(對抗驗證 wmgsefb9d Lens1 NEEDS_CORRECTION)
#   S5 WEAK_RESIDUAL           : ⚠ 結構性不可達(永遠=0)。n_other 與 hidden_het 皆為「coarse 已分裂」的產物
#                                 (n_other>0 → 100% coarse>=2),與 coarse<2 前提矛盾 = unsatisfiable。弱結構在 coarse<2 全由 S4 承載。
#   S6 NO_STRUCTURE_CLEAN      : coarse<2 & fine<2 & n_other==0 & !hidden_het — 乾淨單群桶(但 ~21852 個 n>=12 的 S6
#                                 在 scipy 0.7×max 幾何平切下會切出≥2群=幾何-顯著性分歧,目前無 geometry_ng 欄位故不可見)


def classify(r):
    cng, fng = r["coarse_ng"], r["fine_ng"]
    nother, unstable, hidden = r["n_other"], r["unstable"], r["hidden_het"]
    aligned = r["aligned"]
    if cng >= 2 and not unstable and aligned:
        return "S1_CONFIRMED_MULTI_ALIGNED"
    if cng >= 2 and not unstable and not aligned:
        return "S2_CONFIRMED_MULTI_UNALIGNED"
    if cng >= 2 and unstable:
        return "S3_UNSTABLE_MULTI"
    if cng < 2 and fng >= 2:
        return "S4_WEAK_FINE_ONLY"
    if cng < 2 and fng < 2 and (nother > 0 or hidden):
        return "S5_WEAK_RESIDUAL"
    return "S6_NO_STRUCTURE_CLEAN"


STATES = ["S1_CONFIRMED_MULTI_ALIGNED", "S2_CONFIRMED_MULTI_UNALIGNED", "S3_UNSTABLE_MULTI",
          "S4_WEAK_FINE_ONLY", "S5_WEAK_RESIDUAL", "S6_NO_STRUCTURE_CLEAN"]
# 弱結構「候選」桶 = S3(不穩定多群)+S4(僅 null90 過). ⚠ 非已驗證真結構:S4 chance-level、S3 seed-borderline.
# 對抗驗證後正名 candidates pending validation(非 "retained real structure"). S5 結構性=0 故不計。
WEAK_BUCKETS = ["S3_UNSTABLE_MULTI", "S4_WEAK_FINE_ONLY", "S5_WEAK_RESIDUAL"]


def main():
    recs = json.load(open(REC))
    per = []
    by_set = defaultdict(Counter)
    for r in recs:
        st = classify(r)
        r2 = dict(r)
        r2["wstate"] = st
        per.append(r2)
        by_set[r["set"]][st] += 1
    total = Counter(p["wstate"] for p in per)

    def block(setlab):
        c = by_set[setlab]
        n = sum(c.values())
        rows = {s: {"n": c[s], "pct": round(100 * c[s] / n, 2)} for s in STATES}
        # 關鍵衍生: 舊 no_structure(coarse<2)=S4+S5+S6;其中弱結構(S4+S5)被舊二分黑洞吃掉的量
        coarse_lt2 = c["S4_WEAK_FINE_ONLY"] + c["S5_WEAK_RESIDUAL"] + c["S6_NO_STRUCTURE_CLEAN"]
        weak_hidden_in_nostruct = c["S4_WEAK_FINE_ONLY"] + c["S5_WEAK_RESIDUAL"]
        return {
            "n": n,
            "states": rows,
            "multi_confident_pct": round(100 * (c["S1_CONFIRMED_MULTI_ALIGNED"] + c["S2_CONFIRMED_MULTI_UNALIGNED"]) / n, 2),
            "weak_candidate_pct": round(100 * sum(c[s] for s in WEAK_BUCKETS) / n, 2),  # S3+S4 = candidates pending validation (NOT confirmed weak structure)
            "truly_clean_pct": round(100 * c["S6_NO_STRUCTURE_CLEAN"] / n, 2),
            "old_no_structure_n": coarse_lt2,
            "old_no_structure_pct": round(100 * coarse_lt2 / n, 2),
            "weak_hidden_inside_old_no_structure_n": weak_hidden_in_nostruct,
            "weak_hidden_inside_old_no_structure_pct": round(100 * weak_hidden_in_nostruct / n, 2),
            "sum_check": sum(rows[s]["n"] for s in STATES) == n,
        }

    out = {
        "source": "phylo_cpp_wg_records.json (C++ native phylo-v4.1, genome-wide HCC1395)",
        "n_total": len(recs),
        "states_def": {
            "S1_CONFIRMED_MULTI_ALIGNED": "coarse>=2 & stable & aligned",
            "S2_CONFIRMED_MULTI_UNALIGNED": "coarse>=2 & stable & !aligned",
            "S3_UNSTABLE_MULTI": "coarse>=2 & unstable(modal<0.7)",
            "S4_WEAK_FINE_ONLY": "coarse<2 & fine>=2 (ONLY looser null90 passed; TP/FP ratio 1.15 = chance-level, statistically indistinguishable from noise -> sub-threshold candidates pending validation, NOT retained real structure)",
            "S5_WEAK_RESIDUAL": "coarse<2 & fine<2 & (n_other>0 or hidden_het) -- STRUCTURALLY UNREACHABLE: both triggers require a passed coarse split (coarse>=2), contradicting coarse<2; count always 0 by construction, not a bug",
            "S6_NO_STRUCTURE_CLEAN": "coarse<2 & fine<2 & no residual = truly single group (but ~21852 with n>=12 would flat-cut >=2 groups = geometry-vs-null divergence, no geometry_ng column = invisible)",
        },
        "caveats": {
            "S4_not_real_structure": "S4 (null90-only) is NOT validated weak structure; TP/FP enrichment 1.15 (chance). Reframe as sub-threshold candidates pending validation.",
            "S5_unreachable": "S5 structurally empty (n_other/hidden_het require coarse>=2). Effectively a 5-state partition presented as 6.",
            "aligned_NA_for_coarse_lt2": "'aligned'/V_hp/V_allele only computed when coarse>=2 (8104/8104 aligned have coarse>=2; 0/24180 coarse==1). For S4/S6 aligned=False is NOT empirical -> treat as NA/untested.",
            "confident_multi_anti_discriminative": "S1+S2 confident-multi is FP-enriched (FP 35.65% > TP 24.66%); 'aligned'=cis-ASM fingerprint (max(V_hp,V_allele)>=0.3), NOT a subclone test. Observation/characterization layer only, NOT a TP/FP filter.",
            "multiple_testing": "No FDR across 34736 loci nor recursive within-locus split tests; RNULL=40 gives high-SE percentile -> null90 vs null95 distinction is coarse.",
            "geometry_gap": "No geometry_ng (flat-cut) column exists; ~21852 S6 loci with n>=12 would flat-cut into >=2 groups (geometry-vs-null divergence) and are invisible. Needs subset rerun; cannot be backfilled (matrices deleted).",
            "scope": "single-sample HCC1395, exploratory tier 2-3, no external truth; resolving cis-ASM vs subclone needs normal cis-control / multi-sSNV CCF / single-cell.",
        },
        "weak_buckets": WEAK_BUCKETS,
        "TP": block("TP"),
        "FP": block("FP"),
        "ALL": {s: {"n": total[s], "pct": round(100 * total[s] / len(recs), 2)} for s in STATES},
        "ALL_sum_check": sum(total[s] for s in STATES) == len(recs),
    }
    json.dump(out, open(f"{A}/phylo_weak_structure_counts.json", "w"), ensure_ascii=False, indent=2)
    json.dump(per, open(f"{A}/phylo_weak_structure_classified.json", "w"), ensure_ascii=False)
    # console echo (供 Read 回真值)
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
