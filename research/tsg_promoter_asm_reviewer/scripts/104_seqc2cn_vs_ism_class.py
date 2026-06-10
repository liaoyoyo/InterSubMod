#!/usr/bin/env python3
"""
104 - Relationship between SEQC2 CN ground-truth and ISM ASM classification.

Question: does the ISM ASM-structure call (PASS / MISSED=over-strict / FILTERED, Tier,
CramersV, Δβ, PERMANOVA, over-strict, independent-validation) ASSOCIATE with SEQC2 copy
number state (gain/loss/loh/neutral)? If structure concentrates in LOH/gain, copy number
may confound the ASM signal (LOH -> one haplotype lost -> imbalanced HP -> sparse table;
gain -> more reads -> more clustering power).

Descriptive cross-tabs + enrichment + per-class metric distributions over all 30,350 loci
(TP+FP). Read-only, all numbers from frozen existence summary + seqc2_cn.json.

Output: display_v2/seqc2_vs_ism.json + printed report.
"""
import csv, json, math
from collections import Counter, defaultdict

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")
SQ = ["neutral", "gain", "loss", "loh"]


def num(r, k):
    try:
        v = float(r.get(k, "")); return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def cvmax(r):
    return max([num(r, k) or 0 for k in ("CramersV", "CramersV_HPFamily", "CramersV_HPFine")])


def cluster_real(r):
    return tb(r, "ClusterPermanovaValid") and (num(r, "ClusterPermanovaP") or 1) <= 0.05 and not tb(r, "ClusterDispersionWarn")


def hp_consistent(r):
    return tb(r, "LabelHPPermanovaValid") and (num(r, "LabelHPPermanovaP") or 1) <= 0.05 and not tb(r, "LabelHPDispersionWarn")


def cramers_v_from_table(tab):
    """tab: dict[(row,col)]=count -> Cramér's V."""
    rows = sorted({r for r, c in tab}); cols = sorted({c for r, c in tab})
    n = sum(tab.values())
    if n == 0:
        return 0.0
    rt = {r: sum(tab.get((r, c), 0) for c in cols) for r in rows}
    ct = {c: sum(tab.get((r, c), 0) for r in rows) for c in cols}
    chi2 = 0.0
    for r in rows:
        for c in cols:
            e = rt[r] * ct[c] / n
            if e > 0:
                chi2 += (tab.get((r, c), 0) - e) ** 2 / e
    k = min(len(rows), len(cols))
    return math.sqrt(chi2 / (n * (k - 1))) if k > 1 else 0.0


def main():
    seqc2 = json.load(open(f"{DV}/seqc2_cn.json"))
    tier = json.load(open(f"{DV}/tier_assignment.json"))
    vloci = json.load(open(f"{DV}/validation3.json"))["loci"]

    L = []
    for cls in ("tp", "fp"):
        for r in csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")):
            key = f'{r["Chr"]}_{r["Pos"]}'
            cv = cvmax(r); db = abs(num(r, "HPMergedDelta") or 0); reads = num(r, "NumReads") or 0
            sig = tb(r, "Significant")
            overstrict = (reads >= 20 and tb(r, "PassedGating") and cluster_real(r) and hp_consistent(r) and cv < 0.1)
            cat = 0 if sig else (1 if overstrict else 2)        # PASS / MISSED / FILTERED
            sqc = seqc2.get(key, [0, 2.0])[0]
            t = {"A": 1, "B": 2}.get(tier.get(key, ""), 0)
            v = vloci.get(key, {})
            val = 1 if v.get("validated") else (0 if v.get("ok") else -1)
            mhn = min(num(r, "HP1FamilyN") or 0, num(r, "HP2FamilyN") or 0)
            L.append(dict(cls=cls, cat=cat, sqc=sqc, tier=t, val=val, cv=cv, db=db,
                          permf=num(r, "ClusterPermanovaF") or 0, reads=reads, over=overstrict, mhn=mhn,
                          ismloh=tb(r, "Potential_LOH")))

    n = len(L)
    out = {"n": n, "seqc2_dist": dict(Counter(SQ[x["sqc"]] for x in L))}

    # A. SEQC2 class x ISM cat contingency + row%
    tab = Counter((x["sqc"], x["cat"]) for x in L)
    out["A_contingency"] = {SQ[s]: {["PASS", "MISSED", "FILTERED"][c]: tab.get((s, c), 0) for c in range(3)} for s in range(4)}
    out["A_cramersV"] = round(cramers_v_from_table(tab), 3)
    # B. per-class rates
    by = defaultdict(list)
    for x in L:
        by[x["sqc"]].append(x)
    rates = {}
    for s in range(4):
        g = by[s]; ng = len(g) or 1
        rates[SQ[s]] = dict(
            n=len(g),
            pass_pct=round(100 * sum(1 for x in g if x["cat"] == 0) / ng, 1),
            missed_pct=round(100 * sum(1 for x in g if x["cat"] == 1) / ng, 1),
            overstrict_pct=round(100 * sum(1 for x in g if x["over"]) / ng, 1),
            tierA_pct=round(100 * sum(1 for x in g if x["tier"] == 1) / ng, 1),
            ismLOH_pct=round(100 * sum(1 for x in g if x["ismloh"]) / ng, 1),
            median_reads=sorted(x["reads"] for x in g)[len(g) // 2] if g else 0,
            mean_cv=round(sum(x["cv"] for x in g) / ng, 3),
            mean_db=round(sum(x["db"] for x in g) / ng, 3),
            mean_permf=round(sum(x["permf"] for x in g) / ng, 1),
            val_pass_rate=round(100 * sum(1 for x in g if x["val"] == 1) / max(1, sum(1 for x in g if x["val"] >= 0)), 1),
        )
    out["B_per_class_rates"] = rates
    # C. over-strict vs ISM-LOH and SEQC2-loh (the sparse-table mechanism)
    over = [x for x in L if x["over"]]
    out["C_overstrict_profile"] = dict(
        n_overstrict=len(over),
        pct_in_seqc2_loh=round(100 * sum(1 for x in over if x["sqc"] == 3) / max(1, len(over)), 1),
        pct_in_seqc2_gain=round(100 * sum(1 for x in over if x["sqc"] == 1) / max(1, len(over)), 1),
        pct_ismLOH=round(100 * sum(1 for x in over if x["ismloh"]) / max(1, len(over)), 1),
        median_minHP=sorted(x["mhn"] for x in over)[len(over) // 2] if over else 0,
        baseline_seqc2_loh_pct=round(100 * sum(1 for x in L if x["sqc"] == 3) / n, 1),
    )
    # D. PASS enrichment vs base per SEQC2 class (odds ratio-ish)
    base_pass = sum(1 for x in L if x["cat"] == 0) / n
    out["D_pass_enrichment"] = {SQ[s]: round((rates[SQ[s]]["pass_pct"] / 100) / base_pass, 2) for s in range(4)}

    json.dump(out, open(f"{DV}/seqc2_vs_ism.json", "w"), indent=1)

    # print
    print(f"n={n}  SEQC2 dist: {out['seqc2_dist']}")
    print(f"\n=== A. SEQC2 class × ISM 分類（association Cramér's V = {out['A_cramersV']}）===")
    print(f"{'SEQC2':9s} {'n':>7s} {'PASS%':>7s} {'MISSED%':>8s} {'overstrict%':>12s} {'TierA%':>7s} {'ISM-LOH%':>9s} {'medReads':>9s} {'meanCV':>7s} {'meanΔβ':>7s} {'val通過%':>8s}")
    for s in range(4):
        r = rates[SQ[s]]
        print(f"{SQ[s]:9s} {r['n']:>7d} {r['pass_pct']:>7} {r['missed_pct']:>8} {r['overstrict_pct']:>12} {r['tierA_pct']:>7} {r['ismLOH_pct']:>9} {r['median_reads']:>9} {r['mean_cv']:>7} {r['mean_db']:>7} {r['val_pass_rate']:>8}")
    c = out["C_overstrict_profile"]
    print(f"\n=== C. over-strict(稀疏歸零) 的 CN profile ===")
    print(f"  over-strict n={c['n_overstrict']}; 在 SEQC2-loh={c['pct_in_seqc2_loh']}% (基線 {c['baseline_seqc2_loh_pct']}%); 在 gain={c['pct_in_seqc2_gain']}%; ISM-LOH={c['pct_ismLOH']}%; median minHP={c['median_minHP']}")
    print(f"\n=== D. PASS 富集倍數 (vs 全體 base) ===  {out['D_pass_enrichment']}")
    print(f"\n[104] wrote {DV}/seqc2_vs_ism.json")


if __name__ == "__main__":
    main()
