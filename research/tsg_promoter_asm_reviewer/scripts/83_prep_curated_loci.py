#!/usr/bin/env python3
"""
83 - Prepare the curated locus set for the redesigned display.

For HCC1395 TP and FP, select every locus the user needs to JUDGE:
  - PASS    : Significant==true  (kept by the gate)  -> confirm they are real
  - MISSED  : Significant==false AND |HPMergedDelta|>=0.2 (filtered out but
              big mean-shift)  -> let the user judge if any should be included
For each locus, record the gate inputs + the *binding reject reason* so the
display can show exactly WHY it was kept or dropped.

Outputs (under genome_survey_v2/cn_confound/cross_sample/display_v2/):
  curated_loci.json   : list of dicts (chr,pos,cls,category,cv*,db,p,reads,loh,reject)
  tp_pos.txt / fp_pos.txt : bcftools -T position files (chr<TAB>pos)
"""
import csv, json, os

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
       "genome_survey_v2/cn_confound/cross_sample/display_v2")
os.makedirs(OUT, exist_ok=True)

DB_HIGH = 0.20   # |Δβ| threshold for "high mean-shift" missed loci


def num(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def reject_reason(r):
    """Binding reason a locus failed the Significant gate (gate =
    passed_gating & global_p<=0.05 & cramers_v>=0.1 & num_reads>=20)."""
    reasons = []
    reads = num(r, "NumReads") or 0
    if reads < 20:
        reasons.append(f"覆蓋不足 reads={int(reads)}<20")
    if not tb(r, "PassedGating"):
        reasons.append("未過 gating（結構/離散）")
    gp = num(r, "GlobalP")
    if gp is not None and gp > 0.05:
        reasons.append(f"p 不顯著 global_p={gp:.3g}>0.05")
    # best CramersV across the 3 axes (ALT / HPFamily / HPFine)
    cvs = [num(r, "CramersV"), num(r, "CramersV_HPFamily"), num(r, "CramersV_HPFine")]
    cvmax = max([c for c in cvs if c is not None] or [0.0])
    if cvmax < 0.1:
        reasons.append(f"無分群關聯 CramersV_max={cvmax:.2f}<0.1")
    return reasons


def pack(r, cls, category):
    cvs = [num(r, "CramersV"), num(r, "CramersV_HPFamily"), num(r, "CramersV_HPFine")]
    cvmax = max([c for c in cvs if c is not None] or [0.0])
    return dict(
        chr=r.get("Chr"), pos=int(r.get("Pos")), cls=cls, category=category,
        minhpn=int(min(num(r, "HP1FamilyN") or 0, num(r, "HP2FamilyN") or 0)),
        cv_alt=round(num(r, "CramersV") or 0.0, 3),
        cv_hpfam=round(num(r, "CramersV_HPFamily") or 0.0, 3),
        cv_hpfine=round(num(r, "CramersV_HPFine") or 0.0, 3),
        cv_max=round(cvmax, 3),
        db=round(num(r, "HPMergedDelta") or 0.0, 3),
        gp=num(r, "GlobalP"),
        reads=int(num(r, "NumReads") or 0),
        ncpg=int(num(r, "NumCpGs") or 0),
        loh=tb(r, "Potential_LOH"),
        passed_gating=tb(r, "PassedGating"),
        dispersion_warn=tb(r, "ClusterDispersionWarn"),
        significant=tb(r, "Significant"),
        reject=("" if category == "PASS" else "；".join(reject_reason(r))),
    )


def main():
    curated, tp_pos, fp_pos = [], [], []
    counts = {}
    for cls in ("tp", "fp"):
        rows = list(csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")))
        for r in rows:
            sig = tb(r, "Significant")
            db = abs(num(r, "HPMergedDelta") or 0.0)
            if sig:
                cat = "PASS"
            elif db >= DB_HIGH:
                cat = "MISSED"
            else:
                continue
            d = pack(r, cls, cat)
            curated.append(d)
            (tp_pos if cls == "tp" else fp_pos).append((d["chr"], d["pos"]))
            counts[(cls, cat)] = counts.get((cls, cat), 0) + 1

    # de-dup positions (a chr:pos appears once per class)
    for name, lst in (("tp_pos.txt", tp_pos), ("fp_pos.txt", fp_pos)):
        seen, lines = set(), []
        for c, p in lst:
            if (c, p) not in seen:
                seen.add((c, p)); lines.append(f"{c}\t{p}")
        with open(f"{OUT}/{name}", "w") as f:
            f.write("\n".join(lines) + "\n")

    with open(f"{OUT}/curated_loci.json", "w") as f:
        json.dump(curated, f)

    print(f"[83] curated total = {len(curated)}")
    for k in sorted(counts):
        print(f"   {k[0].upper()} {k[1]:7s} = {counts[k]}")
    print(f"[83] tp_pos={len(tp_pos)}  fp_pos={len(fp_pos)}")
    print(f"[83] wrote {OUT}/curated_loci.json")


if __name__ == "__main__":
    main()
