#!/usr/bin/env python3
"""
91 - From the frozen existence-scan summary (gated CramersV + PERMANOVA), classify
every TP/FP locus into the confidence tier and export:
  - tierA_high_confidence.tsv : the 898 highest-confidence HP-consistent ASM loci
    (reliable CramersV >=0.1 AND PERMANOVA cluster+HP clean AND min-HP-N>=10)
  - overstrict_tp_pos.txt / overstrict_fp_pos.txt : positions for the genome-wide raw
    CramersV re-run (script 92) confirming the reliability-gated cause beyond the
    425-locus curated subset.
  - tier_assignment.json : {chr_pos: tier} for the display.

Tier (minN=10): A = reliable_cv & cluster_real & hp_consistent & adequate
                B = !reliable_cv & cluster_real & hp_consistent & adequate
                '' otherwise
"""
import csv, json

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")
MINN = 10


def num(r, k):
    try:
        v = float(r.get(k, "")); return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def cvmax(r):
    return max([num(r, k) or 0 for k in ("CramersV", "CramersV_HPFamily", "CramersV_HPFine")])


def minhpn(r):
    return min(num(r, "HP1FamilyN") or 0, num(r, "HP2FamilyN") or 0)


def cluster_real(r):
    return tb(r, "ClusterPermanovaValid") and (num(r, "ClusterPermanovaP") or 1) <= 0.05 and not tb(r, "ClusterDispersionWarn")


def hp_consistent(r):
    return tb(r, "LabelHPPermanovaValid") and (num(r, "LabelHPPermanovaP") or 1) <= 0.05 and not tb(r, "LabelHPDispersionWarn")


def tier(r):
    if (num(r, "NumReads") or 0) < 20 or not cluster_real(r) or not hp_consistent(r) or minhpn(r) < MINN:
        # over-strict (structure but inadequate/no cv) still flagged for the re-run set
        return ""
    return "A" if cvmax(r) >= 0.1 else "B"


def is_overstrict(r):
    return ((num(r, "NumReads") or 0) >= 20 and tb(r, "PassedGating")
            and cluster_real(r) and hp_consistent(r) and cvmax(r) < 0.1)


def main():
    assign = {}
    tierA = []
    os_tp, os_fp = [], []
    counts = {"A": 0, "B": 0, "": 0}
    for cls in ("tp", "fp"):
        rows = list(csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")))
        for r in rows:
            t = tier(r)
            counts[t] += 1
            assign[f"{r['Chr']}_{r['Pos']}"] = t
            if t == "A" and cls == "tp":
                tierA.append(r)
            if is_overstrict(r):
                (os_tp if cls == "tp" else os_fp).append((r["Chr"], r["Pos"]))

    # export Tier A high-confidence (TP)
    tierA.sort(key=lambda r: -(num(r, "LabelHPPermanovaF") or 0))
    with open(f"{DV}/tierA_high_confidence.tsv", "w") as f:
        f.write("Chr\tPos\tNumReads\tHP1N\tHP2N\tminHPN\tgatedCVmax\tDbeta\t"
                "ClusterPermF\tHPpermF\tLOH\tQuality_Tier\tLOH_Bed_Annotation\n")
        for r in tierA:
            f.write(f"{r['Chr']}\t{r['Pos']}\t{int(num(r,'NumReads') or 0)}\t"
                    f"{int(num(r,'HP1FamilyN') or 0)}\t{int(num(r,'HP2FamilyN') or 0)}\t{minhpn(r)}\t"
                    f"{cvmax(r):.3f}\t{num(r,'HPMergedDelta') or 0:+.3f}\t"
                    f"{num(r,'ClusterPermanovaF') or 0:.1f}\t{num(r,'LabelHPPermanovaF') or 0:.1f}\t"
                    f"{tb(r,'Potential_LOH')}\t{r.get('Quality_Tier','')}\t{r.get('LOH_Bed_Annotation','')}\n")

    # position files for the genome-wide raw re-run
    for name, lst in (("overstrict_tp_pos.txt", os_tp), ("overstrict_fp_pos.txt", os_fp)):
        seen, out = set(), []
        for c, p in lst:
            if (c, p) not in seen:
                seen.add((c, p)); out.append(f"{c}\t{p}")
        with open(f"{DV}/{name}", "w") as f:
            f.write("\n".join(out) + "\n")

    with open(f"{DV}/tier_assignment.json", "w") as f:
        json.dump(assign, f)

    print(f"[91] tiers: A={counts['A']} B={counts['B']} none={counts['']}")
    print(f"[91] Tier A TP exported = {len(tierA)} -> tierA_high_confidence.tsv")
    print(f"[91] over-strict pos: tp={len(os_tp)} fp={len(os_fp)} (for raw re-run)")
    print(f"[91] tier_assignment.json = {len(assign)} loci")


if __name__ == "__main__":
    main()
