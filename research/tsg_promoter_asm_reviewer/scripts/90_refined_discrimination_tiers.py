#!/usr/bin/env python3
"""
90 - (a) Confirm the TWO causes of "clean HP-aligned structure but gated CramersV=0"
        (over-strict), and (b) design a tiered high-confidence HP-tag-consistent ASM
        criterion that minimises misjudgment (false structure + over-strict miss).

(a) Cause split (uses raw CramersV from the curated re-run manifest where available):
    - reliability-gated  : raw CramersV high (>=0.3) but gated->0 (sparse contingency,
                           Cochran min-expected<5). chi2 IS strong but untrustworthy.
    - partition-mismatch : even raw CramersV low (<0.1) yet PERMANOVA HP-aligned.
                           The dendrogram cut (optimal_k, discrete) doesn't match HP,
                           but reads DO separate by HP in continuous distance space.

(b) Tiers over full TP/FP (frozen summary). "與 tag 高度一致 + 減少誤判":
    cluster_real  = ClusterPermanova valid&p<=0.05 & no dispersion warn
    hp_consistent = LabelHPPermanova valid&p<=0.05 & no dispersion warn (cluster~HP tag)
    adequate      = min(HP1N,HP2N)>=MINN & NumReads>=20
    reliable_cv   = gated CramersV_max>=0.1   (reliable chi2 agreement)
    Tier A (highest conf): reliable_cv & cluster_real & hp_consistent & adequate
    Tier B (recovered)   : !reliable_cv & cluster_real & hp_consistent & adequate
    reject               : dispersion warn / inadequate / not hp_consistent
    Report count + TP:FP (=discrimination/misjudgment proxy) per tier, sweep MINN.

Output: display_v2/diag90.json
"""
import csv, json

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")


def num(r, k):
    try:
        v = float(r.get(k, "")); return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def load(cls):
    return list(csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")))


def cvmax(r):
    return max([num(r, k) or 0 for k in ("CramersV", "CramersV_HPFamily", "CramersV_HPFine")])


def minhpn(r):
    return min(num(r, "HP1FamilyN") or 0, num(r, "HP2FamilyN") or 0)


def cluster_real(r):
    return tb(r, "ClusterPermanovaValid") and (num(r, "ClusterPermanovaP") or 1) <= 0.05 and not tb(r, "ClusterDispersionWarn")


def hp_consistent(r):
    return tb(r, "LabelHPPermanovaValid") and (num(r, "LabelHPPermanovaP") or 1) <= 0.05 and not tb(r, "LabelHPDispersionWarn")


def main():
    tp, fp = load("tp"), load("fp")
    base = len(tp) / len(fp)
    out = {"base_ratio": round(base, 1), "n_tp": len(tp), "n_fp": len(fp)}

    # ---- (a) cause split using raw CramersV from curated manifest ----
    man = {f"{m['chr']}_{m['pos']}": m for m in json.load(open(f"{DV}/manifest.json"))}
    overstrict_tp = [r for r in tp if (num(r, "NumReads") or 0) >= 20 and tb(r, "PassedGating")
                     and cluster_real(r) and hp_consistent(r) and cvmax(r) < 0.1]
    have_raw = reliability = mismatch = mid = 0
    for r in overstrict_tp:
        m = man.get(f"{r['Chr']}_{r['Pos']}")
        if not m or "raw_max" not in m:
            continue
        have_raw += 1
        rm = m["raw_max"]
        if rm >= 0.3:
            reliability += 1     # strong raw chi2 but gated (sparse/unreliable)
        elif rm < 0.1:
            mismatch += 1        # even raw chi2 low -> partition vs continuous mismatch
        else:
            mid += 1
    out["cause_split"] = dict(
        overstrict_tp=len(overstrict_tp), with_raw_available=have_raw,
        reliability_gated_rawCV_ge0_3=reliability, partition_mismatch_rawCV_lt0_1=mismatch, mid_0_1to0_3=mid,
        note="subset with raw CramersV = curated high-|Δβ| loci; full-genome split needs full re-run")

    # ---- (b) tiered criterion sweep ----
    tiers = []
    for minn in (5, 8, 10, 15):
        def classify(rows):
            A = B = 0
            for r in rows:
                if (num(r, "NumReads") or 0) < 20 or not cluster_real(r) or not hp_consistent(r) or minhpn(r) < minn:
                    continue
                if cvmax(r) >= 0.1:
                    A += 1
                else:
                    B += 1
            return A, B
        a_tp, b_tp = classify(tp)
        a_fp, b_fp = classify(fp)
        ab_tp, ab_fp = a_tp + b_tp, a_fp + b_fp
        tiers.append(dict(
            minN=minn,
            tierA_tp=a_tp, tierA_fp=a_fp, tierA_ratio=round(a_tp / max(a_fp, 1), 1),
            tierA_enrich=round((a_tp / max(a_fp, 1)) / base, 2),
            tierB_tp=b_tp, tierB_fp=b_fp, tierB_ratio=round(b_tp / max(b_fp, 1), 1),
            tierB_enrich=round((b_tp / max(b_fp, 1)) / base, 2),
            AB_tp=ab_tp, AB_fp=ab_fp, AB_ratio=round(ab_tp / max(ab_fp, 1), 1),
            AB_enrich=round((ab_tp / max(ab_fp, 1)) / base, 2)))
    out["tiers"] = tiers
    out["orig_significant"] = dict(tp=sum(1 for r in tp if tb(r, "Significant")),
                                   fp=sum(1 for r in fp if tb(r, "Significant")))
    o = out["orig_significant"]
    out["orig_significant"]["ratio"] = round(o["tp"] / max(o["fp"], 1), 1)
    out["orig_significant"]["enrich"] = round((o["tp"] / max(o["fp"], 1)) / base, 2)

    with open(f"{DV}/diag90.json", "w") as f:
        json.dump(out, f, indent=2)

    cs = out["cause_split"]
    print(f"base TP:FP={base:.1f}:1\n")
    print("== (a) cause of over-strict (CramersV=0 despite clean HP structure) ==")
    print(f"   over-strict TP={cs['overstrict_tp']}; raw available (curated)={cs['with_raw_available']}")
    print(f"     reliability-gated (rawCV>=0.3, sparse)   = {cs['reliability_gated_rawCV_ge0_3']}")
    print(f"     partition-mismatch (rawCV<0.1)           = {cs['partition_mismatch_rawCV_lt0_1']}")
    print(f"     mid (0.1-0.3)                            = {cs['mid_0_1to0_3']}")
    print(f"\n== (b) tiered HP-consistent criterion (TP:FP = lower misjudgment) ==")
    print(f"   orig Significant: TP={o['tp']} FP={o['fp']} ratio={o['ratio']}:1 ({o['enrich']}x base)")
    print(f"   {'minN':>4} {'TierA(reliable)':>22} {'TierB(recovered)':>22} {'A+B':>18}")
    for t in tiers:
        print(f"   {t['minN']:>4}  A:{t['tierA_tp']:>4}/{t['tierA_fp']:<3}={t['tierA_ratio']:>5}:1({t['tierA_enrich']}x)  "
              f"B:{t['tierB_tp']:>4}/{t['tierB_fp']:<3}={t['tierB_ratio']:>5}:1({t['tierB_enrich']}x)  "
              f"AB:{t['AB_tp']:>4}/{t['AB_fp']:<3}={t['AB_ratio']:>5}:1")
    print(f"\n[90] wrote {DV}/diag90.json")


if __name__ == "__main__":
    main()
