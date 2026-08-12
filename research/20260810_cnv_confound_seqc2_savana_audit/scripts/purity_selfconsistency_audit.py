#!/usr/bin/env python3
"""Per-segment purity self-consistency audit across all SAVANA-called samples.

HCC1395 is the calibration case: SEQC2 truth plus the refit grid showed its
published fit (purity 0.76 / ploidy 1.83) is wrong and the truth is ~1.0 / ~2.95.
This script asks whether a purely internal test -- with no external truth --
would have caught that, and then applies the same test to the samples where no
truth exists.

The test
--------
For a segment with tumour total copy number n and tumour purity rho, the most
extreme allelic imbalance possible is complete LOH (minor copy number 0).  Even
then the observed B-allele fraction is diluted by the normal contamination,
which always contributes one copy per allele:

    BAF_max(n, rho) = (rho * n + (1 - rho)) / (rho * n + 2 * (1 - rho))

An observed mean BAF above BAF_max cannot be produced at that purity by any
allele-specific configuration.  A sample with many such segments has an
under-estimated purity, regardless of how self-consistent the fit looked to the
caller.

Because BAF_max grows with n, the test is evaluated per segment against that
segment's own fitted copy number, rather than against one global threshold.
"""

from __future__ import annotations

import json
import os
from collections import defaultdict

WORK = "/big7_disk/liaoyoyo2001/cnv_sv_work"
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT = os.path.join(OUT_DIR, "purity_selfconsistency_audit.json")

MIN_HET = 10  # minimum het SNPs for a segment's mean BAF to be usable
TOL = 0.02  # BAF tolerance before calling a violation

SAMPLES = {
    "HCC1395": f"{WORK}/HCC1395/savana_wgs/cna",
    "H1437": None,
    "H2009": None,
    "HCC1954": None,
    "HCC1937": None,
    "COLO829": None,
}


def find_cna_dir(sample):
    """Locate the SAVANA cna output directory for a sample."""
    if SAMPLES.get(sample):
        return SAMPLES[sample]
    base = os.path.join(WORK, sample)
    for root, dirs, files in os.walk(base):
        if any(f.endswith("_fitted_purity_ploidy.tsv") for f in files):
            return root
    return None


def baf_max(n, rho):
    """Highest mean BAF achievable at tumour CN n and purity rho (minor CN = 0)."""
    num = rho * n + (1.0 - rho)
    den = rho * n + 2.0 * (1.0 - rho)
    return num / den if den > 0 else None


def audit(sample, cna_dir, override_purity=None):
    fit_path = os.path.join(cna_dir, f"{sample}_fitted_purity_ploidy.tsv")
    seg_path = os.path.join(cna_dir, f"{sample}_segmented_absolute_copy_number.tsv")
    if not (os.path.exists(fit_path) and os.path.exists(seg_path)):
        return None

    with open(fit_path) as fh:
        hdr = fh.readline().split()
        vals = fh.readline().split()
    fit = dict(zip(hdr, vals))
    rho = float(override_purity if override_purity is not None else fit["purity"])

    segs = 0
    usable = 0
    violations = 0
    worst = 0.0
    excess_list = []
    bafs = []

    with open(seg_path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(h)}
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            segs += 1
            try:
                nhet = int(p[col["no_hetSNPs"]])
                baf = float(p[col["meanBAF"]])
                cn = float(p[col["copyNumber"]])
            except (ValueError, IndexError, KeyError):
                continue
            if nhet < MIN_HET:
                continue
            usable += 1
            bafs.append(baf)
            bmax = baf_max(max(cn, 0.01), rho)
            if bmax is None:
                continue
            excess = baf - bmax
            excess_list.append(excess)
            if excess > TOL:
                violations += 1
                worst = max(worst, excess)

    bafs.sort()
    n = len(bafs)
    excess_list.sort()
    m = len(excess_list)

    return {
        "sample": sample,
        "cna_dir": cna_dir,
        "fitted_purity": float(fit["purity"]),
        "fitted_ploidy": float(fit["ploidy"]),
        "purity_used_for_test": rho,
        "segments_total": segs,
        "segments_usable": usable,
        "baf_median": round(bafs[n // 2], 4) if n else None,
        "baf_p90": round(bafs[int(n * 0.9)], 4) if n else None,
        "baf_max": round(bafs[-1], 4) if n else None,
        "violations": violations,
        "violation_rate_percent": round(100.0 * violations / usable, 2) if usable else None,
        "worst_excess_baf": round(worst, 4),
        "median_excess_baf": round(excess_list[m // 2], 4) if m else None,
    }


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    results = {}
    for sample in SAMPLES:
        d = find_cna_dir(sample)
        if not d:
            results[sample] = {"status": "NO_SAVANA_CNA_OUTPUT_FOUND"}
            continue
        r = audit(sample, d)
        results[sample] = r if r else {"status": "MISSING_FILES", "cna_dir": d}

    # Calibration: HCC1395 re-tested at the grid-search optimum
    calib = None
    if results.get("HCC1395", {}).get("cna_dir"):
        calib = audit("HCC1395", results["HCC1395"]["cna_dir"], override_purity=1.0)

    out = {
        "generated_by": os.path.basename(__file__),
        "test": {
            "formula": "BAF_max(n, rho) = (rho*n + (1-rho)) / (rho*n + 2*(1-rho))",
            "interpretation": "observed mean BAF above BAF_max cannot be produced at that purity by any allele-specific configuration",
            "min_het_snps": MIN_HET,
            "tolerance": TOL,
        },
        "calibration_case": {
            "sample": "HCC1395",
            "why": "only sample with an external CN truth (SEQC2); refit grid puts truth at purity~1.0 / ploidy~2.95",
            "at_published_purity": results.get("HCC1395"),
            "at_refit_purity_1.0": calib,
        },
        "samples": results,
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"{'sample':10s} {'purity':>7s} {'ploidy':>7s} {'usable':>7s} {'viol%':>7s} {'BAFmed':>7s} {'worst':>7s}")
    for s, r in results.items():
        if not r or "fitted_purity" not in r:
            print(f"{s:10s} {str(r.get('status','?')):>40s}")
            continue
        print(
            f"{s:10s} {r['fitted_purity']:7.2f} {r['fitted_ploidy']:7.2f} "
            f"{r['segments_usable']:7d} {r['violation_rate_percent']:7.2f} "
            f"{r['baf_median']:7.4f} {r['worst_excess_baf']:7.4f}"
        )
    if calib:
        print(
            f"\ncalibration HCC1395 @ purity=1.0: violations={calib['violations']} "
            f"({calib['violation_rate_percent']}%), worst excess={calib['worst_excess_baf']}"
        )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
