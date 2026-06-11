#!/usr/bin/env python3
"""
95 - Combine the TP/FP-independent validation across the TWO real ASM signal modes,
since a locus can have HP-associated methylation as a mean-SHIFT (Δβ) OR as discrete
CLUSTERING (partition) — different signals (Q3: only 260 overlap). Validating a
clustering locus on the Δβ axis (or vice versa) under-counts it.

  Axis-1 SHIFT     (script 94): HP-permutation on |Δβ| + MWU + CpG split-half +
                   strand concordance + read bootstrap.  Tests a real HP MEAN shift.
  Axis-2 CLUSTER   (ISM native, summary): LabelHPPermanova — HP labels (SNV-derived)
                   permuted on the read×CpG distance matrix; valid & p<=0.05 & no
                   PERMDISP dispersion warning.  Tests reads cluster by HP. Non-circular
                   (HP from germline SNP phasing, not methylation).

  validated_ASM = SHIFT-pass OR CLUSTER-pass.

Output: display_v2/validation2.json  (+ updates per-locus mode label for the display)
"""
import csv, json

DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")
EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"


def num(r, k):
    try:
        v = float(r.get(k, "")); return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def cluster_pass(r):
    return (tb(r, "LabelHPPermanovaValid") and (num(r, "LabelHPPermanovaP") or 1) <= 0.05
            and not tb(r, "LabelHPDispersionWarn") and (num(r, "NumReads") or 0) >= 20)


def main():
    val = json.load(open(f"{DV}/validation.json"))
    shift = val["loci"]                       # key -> Δβ-axis result (has 'validated')
    tier = json.load(open(f"{DV}/tier_assignment.json"))
    man = {m["chr"] + "_" + str(m["pos"]): m for m in json.load(open(f"{DV}/manifest.json"))}

    # load cluster-axis (LabelHPPermanova) for the curated loci from existence summary
    cluster = {}
    for cls in ("tp", "fp"):
        for r in csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")):
            k = f"{r['Chr']}_{r['Pos']}"
            if k in man:
                cluster[k] = cluster_pass(r)

    out = {}
    for k, m in man.items():
        s = shift.get(k, {})
        sh = bool(s.get("validated"))
        cl = bool(cluster.get(k, False))
        out[k] = dict(shift_pass=sh, shift_testable=bool(s.get("ok")),
                      cluster_pass=cl, validated=sh or cl,
                      mode=("both" if sh and cl else "shift" if sh else "cluster" if cl else "none"),
                      perm_q=s.get("perm_q"), boot=s.get("boot_stability"),
                      strand_conc=s.get("strand_concordant"), strand_meas=s.get("strand_measurable"))

    # aggregate by tier
    def agg(keys):
        n = len(keys)
        return dict(n=n,
                    validated=sum(1 for k in keys if out[k]["validated"]),
                    shift=sum(1 for k in keys if out[k]["mode"] == "shift"),
                    cluster=sum(1 for k in keys if out[k]["mode"] == "cluster"),
                    both=sum(1 for k in keys if out[k]["mode"] == "both"),
                    none=sum(1 for k in keys if out[k]["mode"] == "none"))
    by_tier = {t or "none": agg([k for k in out if tier.get(k, "") == t]) for t in ("A", "B", "")}
    dbonly = [k for k in out if abs(man[k].get("db", 0)) >= 0.2 and tier.get(k, "") not in ("A", "B")]
    clusters_small_db = [k for k in out if tier.get(k, "") == "A" and abs(man[k].get("db", 0)) < 0.1]

    summary = dict(
        n=len(out),
        validated_any=sum(1 for k in out if out[k]["validated"]),
        mode_counts={m: sum(1 for k in out if out[k]["mode"] == m) for m in ("both", "shift", "cluster", "none")},
        by_tier=by_tier,
        dbeta_only=dict(n=len(dbonly), validated=sum(1 for k in dbonly if out[k]["validated"]),
                        via_shift=sum(1 for k in dbonly if out[k]["mode"] in ("shift", "both"))),
        tierA_smallDb=dict(n=len(clusters_small_db),
                           validated=sum(1 for k in clusters_small_db if out[k]["validated"]),
                           via_cluster=sum(1 for k in clusters_small_db if out[k]["mode"] in ("cluster", "both"))),
    )
    json.dump(dict(summary=summary, loci=out), open(f"{DV}/validation2.json", "w"))

    print("== two-axis independent validation (TP/FP-free) ==")
    print(f"  curated n={summary['n']}  validated(any axis)={summary['validated_any']} "
          f"({100*summary['validated_any']/summary['n']:.0f}%)")
    print(f"  modes: both={summary['mode_counts']['both']} shift-only={summary['mode_counts']['shift']} "
          f"cluster-only={summary['mode_counts']['cluster']} none={summary['mode_counts']['none']}")
    print("  by Tier:")
    for t in ("A", "B", "none"):
        a = by_tier[t]
        print(f"    {t:4s}: n={a['n']:4d} validated={a['validated']:4d} ({100*a['validated']/max(a['n'],1):.0f}%) "
              f"[both={a['both']} shift={a['shift']} cluster={a['cluster']} none={a['none']}]")
    sa = summary["tierA_smallDb"]
    print(f"  Tier A small-Δβ(<0.1) clustering loci: n={sa['n']} validated={sa['validated']} "
          f"(via cluster axis {sa['via_cluster']}) <- these pass via CLUSTER not Δβ")
    do = summary["dbeta_only"]
    print(f"  Δβ-only (no cluster structure): n={do['n']} validated={do['validated']} via shift {do['via_shift']} "
          f"<- real HP mean-shift, validated independently (not artifact; FP-leaning is a TP/FP property)")
    print(f"[95] wrote {DV}/validation2.json")


if __name__ == "__main__":
    main()
