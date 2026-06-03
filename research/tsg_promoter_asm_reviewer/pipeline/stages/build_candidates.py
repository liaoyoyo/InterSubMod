"""
build_candidates.py — reproducible cis-candidate shortlist for the scan.

Selects HP-axis nonLOH loci with ASM signal (|mean_delta| >= dmin AND wilcoxon p < pmax)
from the dual-axis survey TSV, attaches REF/ALT from the caller VCF. Output feeds
scan_cis_candidates.py. characterization only.

    python3 stages/build_candidates.py --dualaxis <tp.tsv> --vcf <clairs_tp.vcf> --out <json>
"""
import csv
import gzip
import json
import argparse

BONF = 0.05 / 51091


def _fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def build(dualaxis_tsv, vcf, out, dmin=0.10, pmax=0.05):
    cand = {}
    with open(dualaxis_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["axis_type"] != "HP" or r["loh_status"] != "nonLOH":
                continue
            d, p = _fl(r["mean_delta"]), _fl(r["wilcoxon_p"])
            if d is None or p is None or abs(d) < dmin or p >= pmax:
                continue
            k = (r["chrom"], r["somatic_pos"])
            if k not in cand or abs(d) > abs(cand[k][0]):
                cand[k] = (d, p, "bonf" if p < BONF else "nominal")
    positions = {f"{c}\t{p}" for c, p in cand}
    refalt = {}
    op = gzip.open if vcf.endswith(".gz") else open
    with op(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 5:
                continue
            if f"{p[0]}\t{p[1]}" in positions and len(p[3]) == 1 and len(p[4]) == 1:
                refalt[(p[0], p[1])] = (p[3], p[4])
    out_loci = []
    for (c, pos), (d, pv, sig) in sorted(cand.items(), key=lambda x: -abs(x[1][0])):
        if (c, pos) in refalt:
            ref, alt = refalt[(c, pos)]
            out_loci.append({"chrom": c, "pos": int(pos), "ref": ref, "alt": alt,
                             "dbeta_hp_screen": round(d, 3), "sig": sig})
    json.dump(out_loci, open(out, "w"), indent=1)
    n_bonf = sum(1 for r in out_loci if r["sig"] == "bonf")
    print(f"{len(out_loci)} candidates ({n_bonf} Bonferroni-sig, {len(out_loci) - n_bonf} nominal) -> {out}")
    return out_loci


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dualaxis", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dmin", type=float, default=0.10)
    ap.add_argument("--pmax", type=float, default=0.05)
    a = ap.parse_args()
    build(a.dualaxis, a.vcf, a.out, a.dmin, a.pmax)


if __name__ == "__main__":
    main()
