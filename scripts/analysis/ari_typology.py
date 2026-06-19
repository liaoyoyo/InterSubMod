#!/usr/bin/env python3
"""
ari_typology.py — offline 3-case cluster<->label typology from the C2 ARI columns.

Reads a significance_summary.csv emitted by the C2 binary (has ARI_Cluster_HP,
ARI_Cluster_Allele, WithinHP_CleanMultigroup, VerificationClass) and derives:
  - ARI distribution per axis (sentinel -2.0 = not computed)
  - 3-case typology via ARI_max = max(ARI_Cluster_HP, ARI_Cluster_Allele):
      aligned (>=0.5) / partial [0.2,0.5) / no-align (<0.2, computed) / no-label-axis (both -2)
  - within-HP substructure (WithinHP_CleanMultigroup) as an orthogonal flag
  - crosstabs vs VerificationClass

Pure offline read-back — every number is a deterministic count.
Usage: python3 ari_typology.py significance_summary.csv [out.json]
"""
import csv
import json
import math
import sys

TAU_HIGH = 0.5   # ARI >= this -> aligned (substantial agreement)
TAU_LOW = 0.2    # ARI < this (but computed) -> no-align / structure-not-label
SENTINEL = -1.5  # values <= this are the "-2.0 not computed" sentinel


def ffloat(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return math.nan


def bucket(v):
    if math.isnan(v) or v <= SENTINEL:
        return "sentinel(-2)"
    if v < -0.2:
        return "[-1,-0.2)"
    if v < 0.2:
        return "[-0.2,0.2)"
    if v < 0.5:
        return "[0.2,0.5)"
    if v < 0.8:
        return "[0.5,0.8)"
    return "[0.8,1.0]"


def dist(vals):
    comp = [v for v in vals if not (math.isnan(v) or v <= SENTINEL)]
    b = {}
    for v in vals:
        k = bucket(v)
        b[k] = b.get(k, 0) + 1
    comp_sorted = sorted(comp)
    n = len(comp_sorted)
    return {
        "computed": n,
        "sentinel": sum(1 for v in vals if math.isnan(v) or v <= SENTINEL),
        "mean": round(sum(comp) / n, 4) if n else None,
        "median": round(comp_sorted[n // 2], 4) if n else None,
        "buckets": dict(sorted(b.items())),
    }


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: ari_typology.py significance_summary.csv [out.json]")
    path = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else None
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    n = len(rows)

    hp = [ffloat(r.get("ARI_Cluster_HP", "")) for r in rows]
    al = [ffloat(r.get("ARI_Cluster_Allele", "")) for r in rows]

    def is_sent(v):
        return math.isnan(v) or v <= SENTINEL

    typ = {"aligned_ge_0.5": 0, "partial_0.2_0.5": 0, "noalign_lt_0.2_computed": 0, "no_label_axis_both_sentinel": 0}
    typ_by_vc = {}
    within_true = 0
    within_by_typ = {}
    for i, r in enumerate(rows):
        h, a = hp[i], al[i]
        comp = [x for x in (h, a) if not is_sent(x)]
        if not comp:
            t = "no_label_axis_both_sentinel"
        else:
            amax = max(comp)
            if amax >= TAU_HIGH:
                t = "aligned_ge_0.5"
            elif amax >= TAU_LOW:
                t = "partial_0.2_0.5"
            else:
                t = "noalign_lt_0.2_computed"
        typ[t] += 1
        vc = r.get("VerificationClass", "").strip() or "(blank)"
        typ_by_vc.setdefault(t, {})
        typ_by_vc[t][vc] = typ_by_vc[t].get(vc, 0) + 1
        w = str(r.get("WithinHP_CleanMultigroup", "")).strip().lower() in ("true", "1")
        if w:
            within_true += 1
            within_by_typ[t] = within_by_typ.get(t, 0) + 1

    result = {
        "source_csv": path,
        "n": n,
        "ari_hp": dist(hp),
        "ari_allele": dist(al),
        "typology_by_ARImax": {k: {"n": v, "frac": round(v / n, 4) if n else None} for k, v in typ.items()},
        "typology_x_verificationclass": typ_by_vc,
        "within_hp_substructure": {
            "clean_multigroup_true": within_true,
            "frac": round(within_true / n, 4) if n else None,
            "by_typology": within_by_typ,
        },
        "thresholds": {"TAU_HIGH": TAU_HIGH, "TAU_LOW": TAU_LOW},
    }
    if out:
        with open(out, "w") as fh:
            json.dump(result, fh, indent=2)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    if out:
        print(f"\n[written] {out}")


if __name__ == "__main__":
    main()
