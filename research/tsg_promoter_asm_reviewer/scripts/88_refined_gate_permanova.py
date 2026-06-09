#!/usr/bin/env python3
"""
88 - Refined "meaningful ASM" gate: recover latent loci the CramersV reliability
gate (min expected cell <5, Cochran) zeroed, using PERMANOVA (distance-based,
robust to sparse/imbalanced HP groups) instead.

Original gate (Significant) = PassedGating & global_p<=0.05 & CramersV>=0.1 & NumReads>=20
  -> CramersV is chi2-based; unreliable (=>0) when an HP group is tiny (e.g. 3 vs 51).

Refined recovery (latent) = NOT Significant
  & PassedGating & NumReads>=20
  & ClusterPermanovaValid & ClusterPermanovaP<=0.05         (real cluster structure)
  & LabelHPPermanovaValid & LabelHPPermanovaP<=0.05         (clusters ALIGN with HP)
  & min(HP1FamilyN, HP2FamilyN) >= MINN                     (avoid 2-3 read degeneracy)

Sweeps MINN. Reports recovery counts + TP/FP + enrichment vs the 47:1 base rate so we
can state honestly whether the recovered set is TP-discriminative or characterization-only.

Output: display_v2/refined_gate.json
"""
import csv, json

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
       "genome_survey_v2/cn_confound/cross_sample/display_v2/refined_gate.json")
MINNS = [0, 3, 5, 8, 10, 15]


def num(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def load(cls):
    return list(csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")))


def is_latent(r, minn):
    if tb(r, "Significant"):
        return False
    if not tb(r, "PassedGating"):
        return False
    if (num(r, "NumReads") or 0) < 20:
        return False
    if not (tb(r, "ClusterPermanovaValid") and (num(r, "ClusterPermanovaP") or 1) <= 0.05):
        return False
    if not (tb(r, "LabelHPPermanovaValid") and (num(r, "LabelHPPermanovaP") or 1) <= 0.05):
        return False
    mn = min(num(r, "HP1FamilyN") or 0, num(r, "HP2FamilyN") or 0)
    return mn >= minn


def main():
    tp, fp = load("tp"), load("fp")
    n_tp, n_fp = len(tp), len(fp)
    base = (sum(1 for r in tp if tb(r, "Significant")), sum(1 for r in fp if tb(r, "Significant")))
    base_ratio = n_tp / n_fp  # overall TP:FP among all analyzed loci
    print(f"analyzed: TP={n_tp} FP={n_fp} (base TP:FP={base_ratio:.1f}:1)")
    print(f"original Significant PASS: TP={base[0]} FP={base[1]}  (TP:FP={base[0]/max(base[1],1):.1f}:1)\n")
    print(f"{'minN':>5} {'recTP':>6} {'recFP':>6} {'rec TP:FP':>10} {'enrich×base':>12} "
          f"{'combTP':>7} {'combFP':>7} {'recLOH%':>8}")

    sweep = []
    for minn in MINNS:
        rtp = [r for r in tp if is_latent(r, minn)]
        rfp = [r for r in fp if is_latent(r, minn)]
        rec_ratio = len(rtp) / max(len(rfp), 1)
        enrich = rec_ratio / base_ratio
        comb_tp, comb_fp = base[0] + len(rtp), base[1] + len(rfp)
        loh_pct = 100 * sum(1 for r in rtp if tb(r, "Potential_LOH")) / max(len(rtp), 1)
        print(f"{minn:>5} {len(rtp):>6} {len(rfp):>6} {rec_ratio:>9.1f}:1 {enrich:>11.2f}× "
              f"{comb_tp:>7} {comb_fp:>7} {loh_pct:>7.0f}%")
        sweep.append(dict(minN=minn, rec_tp=len(rtp), rec_fp=len(rfp),
                          rec_ratio=round(rec_ratio, 2), enrich_vs_base=round(enrich, 2),
                          comb_tp=comb_tp, comb_fp=comb_fp, rec_loh_pct=round(loh_pct, 1)))

    out = dict(n_tp=n_tp, n_fp=n_fp, base_ratio=round(base_ratio, 1),
               orig_pass_tp=base[0], orig_pass_fp=base[1], sweep=sweep,
               note=("recovered set TP:FP vs base 47:1 -> enrich~1 means characterization-only "
                     "(NOT a TP/FP discriminator, consistent with the concluded methylation-filter NEGATIVE)"))
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n[88] wrote {OUT}")


if __name__ == "__main__":
    main()
