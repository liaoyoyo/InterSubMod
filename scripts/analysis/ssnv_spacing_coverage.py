#!/usr/bin/env python3
"""
T-Q1-COVERAGE — HCC1395 sSNV spacing x ONT read length coverage.

Paper motivation quantification: how sparse are somatic sSNVs genome-wide, and
can a single ONT read directly link >=2 sSNVs? This is the quantitative basis for
"why sSNV alone is hard to reconstruct subclones" and the justification for using
methylation as a complementary per-read signal.

All numbers are written to a JSON file. NO numbers are printed-as-truth here;
the report is written in a separate step after reading the JSON back (project
section 13.0 anti-fabrication rule).

Inputs (canonical paired HCC1395 pileup pipeline; see scripts/run_vcf_all_snv.sh:46-48):
  - TP sSNV VCF : /big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz
  - FP sSNV VCF : /big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz
  - tumor BAM   : /big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam (symlink to
                  /big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam)

Definitions / scope choices (recorded in JSON for provenance):
  - sSNV set = TP (ClairS pileup calls confirmed against SEQC2 truth). This is the
    set ISM consumes as "somatic". FP set analysed separately for context only.
  - Autosomes (chr1-22) only; all records are PASS filter (verified empirically).
  - "read length" = REFERENCE-CONSUMED span (reference_end - reference_start), i.e.
    the genomic interval a read covers, because that is what determines how many
    genomic sSNV positions a read can overlap. Query (SEQ) length also recorded.
  - Genome size denominator = sum of autosome lengths from the BAM header (hg38).

Outputs JSON keys: see build at end of main().
"""

import argparse
import gzip
import json
import os
import random
import sys
from bisect import bisect_left, bisect_right
from datetime import datetime, timezone

import numpy as np
import pysam

TP_VCF = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz"
FP_VCF = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz"
TUMOR_BAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
WINDOW = 5000  # +/-5000 bp ISM methylation window (Config window; section C2 of cognition doc)


def parse_vcf(path):
    """Return dict chrom -> sorted list of POS (1-based), plus filter/contig stats."""
    by_chrom = {}
    n_total = 0
    n_pass = 0
    n_nonpass = 0
    n_nonautosome = 0
    filters_seen = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom = f[0]
            pos = int(f[1])
            filt = f[6]
            n_total += 1
            filters_seen[filt] = filters_seen.get(filt, 0) + 1
            if filt == "PASS" or filt == ".":
                n_pass += 1
            else:
                n_nonpass += 1
            if chrom not in AUTOSOMES:
                n_nonautosome += 1
                continue
            by_chrom.setdefault(chrom, []).append(pos)
    for c in by_chrom:
        by_chrom[c].sort()
    return {
        "by_chrom": by_chrom,
        "n_total_records": n_total,
        "n_pass": n_pass,
        "n_nonpass": n_nonpass,
        "n_nonautosome": n_nonautosome,
        "n_autosome_used": sum(len(v) for v in by_chrom.values()),
        "filters_seen": filters_seen,
    }


def gap_distribution(by_chrom):
    """Pairwise adjacent-sSNV gaps (within-chromosome only)."""
    gaps = []
    per_chrom_counts = {}
    for c, positions in by_chrom.items():
        per_chrom_counts[c] = len(positions)
        for i in range(1, len(positions)):
            gaps.append(positions[i] - positions[i - 1])
    gaps = np.asarray(gaps, dtype=np.int64)
    return gaps, per_chrom_counts


def pct(arr, qs):
    return {f"p{q}": float(np.percentile(arr, q)) for q in qs}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tp-vcf", default=TP_VCF)
    ap.add_argument("--fp-vcf", default=FP_VCF)
    ap.add_argument("--bam", default=TUMOR_BAM)
    ap.add_argument("--n-reads", type=int, default=200000,
                    help="number of mapped reads to sample for length distribution")
    ap.add_argument("--n-reads-empirical", type=int, default=50000,
                    help="number of reads to intersect with sSNV coords (empirical >=2 check)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    result = {
        "task": "T-Q1-COVERAGE",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "tp_vcf": os.path.realpath(args.tp_vcf),
            "fp_vcf": os.path.realpath(args.fp_vcf),
            "bam": os.path.realpath(args.bam),
            "bam_symlink": args.bam,
        },
        "scope": {
            "ssnv_set_primary": "TP (ClairS pileup confirmed vs SEQC2 truth)",
            "chromosomes": "autosomes chr1-22 only",
            "read_length_metric": "reference-consumed span (reference_end - reference_start)",
            "window_bp": WINDOW,
        },
    }

    # ---- (a) sSNV totals ----
    tp = parse_vcf(args.tp_vcf)
    fp = parse_vcf(args.fp_vcf)
    result["ssnv_totals"] = {
        "tp_n_total_records": tp["n_total_records"],
        "tp_n_pass": tp["n_pass"],
        "tp_n_nonpass": tp["n_nonpass"],
        "tp_n_autosome_used": tp["n_autosome_used"],
        "tp_n_nonautosome": tp["n_nonautosome"],
        "tp_filters_seen": tp["filters_seen"],
        "fp_n_total_records": fp["n_total_records"],
        "fp_n_autosome_used": fp["n_autosome_used"],
    }

    # ---- genome size from BAM header (autosomes) ----
    bam = pysam.AlignmentFile(args.bam, "rb")
    lengths = dict(zip(bam.references, bam.lengths))
    autosome_len = {c: lengths[c] for c in AUTOSOMES if c in lengths}
    genome_size = int(sum(autosome_len.values()))
    result["genome"] = {
        "autosome_total_bp": genome_size,
        "n_autosomes": len(autosome_len),
        "missing_autosomes": [c for c in AUTOSOMES if c not in lengths],
    }

    # ---- (b) sSNV spacing distribution (TP) ----
    gaps, per_chrom_counts = gap_distribution(tp["by_chrom"])
    n_ssnv = tp["n_autosome_used"]
    result["spacing"] = {
        "n_ssnv": n_ssnv,
        "n_adjacent_pairs": int(gaps.size),
        "gap_median_bp": float(np.median(gaps)),
        "gap_mean_bp": float(np.mean(gaps)),
        "gap_min_bp": int(np.min(gaps)),
        "gap_max_bp": int(np.max(gaps)),
        "gap_percentiles_bp": pct(gaps, [1, 5, 10, 25, 50, 75, 90, 95, 99]),
        "mean_density_ssnv_per_mb": float(n_ssnv / (genome_size / 1e6)),
        "expected_uniform_spacing_bp": float(genome_size / n_ssnv),
        "frac_gaps_le_50kb": float(np.mean(gaps <= 50000)),
        "frac_gaps_le_window_2x": float(np.mean(gaps <= 2 * WINDOW)),  # both within one +/-5kb window
        "per_chrom_ssnv_counts": per_chrom_counts,
    }

    # ---- (c) ONT read length distribution (sampled) ----
    ref_spans = []
    query_lens = []
    n_seen = 0
    n_target = args.n_reads
    # iterate mapped primary reads; stop after n_target
    for aln in bam.fetch():
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        rs = aln.reference_start
        re = aln.reference_end
        if rs is None or re is None:
            continue
        ref_spans.append(re - rs)
        ql = aln.query_length or aln.infer_read_length()
        if ql:
            query_lens.append(ql)
        n_seen += 1
        if n_seen >= n_target:
            break
    ref_spans = np.asarray(ref_spans, dtype=np.int64)
    query_lens = np.asarray(query_lens, dtype=np.int64)
    result["read_length"] = {
        "n_reads_sampled": int(ref_spans.size),
        "note": "first N mapped primary reads via fetch() (genome-ordered, chr1-start onward)",
        "ref_span_median_bp": float(np.median(ref_spans)),
        "ref_span_mean_bp": float(np.mean(ref_spans)),
        "ref_span_percentiles_bp": pct(ref_spans, [1, 5, 10, 25, 50, 75, 90, 95, 99]),
        "query_len_median_bp": float(np.median(query_lens)) if query_lens.size else None,
        "query_len_mean_bp": float(np.mean(query_lens)) if query_lens.size else None,
    }

    # ---- (d) P(read spans >=2 sSNV): analytic (Poisson) weighted by read length ----
    density_per_bp = n_ssnv / genome_size  # mean sSNV per bp (autosome avg)
    # For a read of reference span L placed uniformly, expected sSNV count lambda = density * L.
    # P(>=2) under Poisson = 1 - e^-lam - lam e^-lam ; P(>=1) = 1 - e^-lam.
    lam = density_per_bp * ref_spans.astype(np.float64)
    p_ge1 = 1.0 - np.exp(-lam)
    p_ge2 = 1.0 - np.exp(-lam) - lam * np.exp(-lam)
    result["span_analytic_poisson"] = {
        "ssnv_density_per_bp": float(density_per_bp),
        "ssnv_density_per_kb": float(density_per_bp * 1000),
        "lambda_at_median_read": float(density_per_bp * np.median(ref_spans)),
        "P_ge1_ssnv_read_weighted": float(np.mean(p_ge1)),
        "P_ge2_ssnv_read_weighted": float(np.mean(p_ge2)),
        "method": "Poisson with lambda=density*ref_span, averaged over sampled read-span distribution",
    }

    # ---- (d-empirical) actually intersect a batch of reads with sSNV coords ----
    # Build per-chrom sorted arrays for fast bisect.
    tp_arrays = {c: np.asarray(v, dtype=np.int64) for c, v in tp["by_chrom"].items()}
    counts_hist = {}  # n_ssnv_on_read -> n_reads
    n_emp = 0
    n_emp_target = args.n_reads_empirical
    emp_ge1 = 0
    emp_ge2 = 0
    for aln in bam.fetch():
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        c = aln.reference_name
        if c not in tp_arrays:
            continue
        rs = aln.reference_start
        re = aln.reference_end
        if rs is None or re is None:
            continue
        arr = tp_arrays[c]
        # sSNV POS are 1-based; reference_start 0-based, reference_end exclusive 0-based.
        # A read covers genomic positions [rs+1, re] in 1-based terms.
        lo = bisect_left(arr, rs + 1)
        hi = bisect_right(arr, re)
        k = hi - lo
        counts_hist[k] = counts_hist.get(k, 0) + 1
        if k >= 1:
            emp_ge1 += 1
        if k >= 2:
            emp_ge2 += 1
        n_emp += 1
        if n_emp >= n_emp_target:
            break
    result["span_empirical"] = {
        "n_reads_intersected": n_emp,
        "note": "reads overlapping sSNV coords by genomic interval; spans gap if read covers >=2 sSNV positions",
        "P_ge1_ssnv": float(emp_ge1 / n_emp) if n_emp else None,
        "P_ge2_ssnv": float(emp_ge2 / n_emp) if n_emp else None,
        "frac_reads_0_ssnv": float(counts_hist.get(0, 0) / n_emp) if n_emp else None,
        "ssnv_per_read_hist": {str(k): counts_hist[k] for k in sorted(counts_hist)},
        "max_ssnv_on_single_read": max(counts_hist) if counts_hist else 0,
    }

    # ---- (e) +/-5000bp window coverage vs read length distribution ----
    # The ISM window per sSNV is +/-5000 = 10000 bp wide. How does that compare
    # to read spans? Fraction of reads whose ref span >= window width (10kb), and
    # fraction >= half-window (5kb, can reach the SNV from one side fully).
    win_full = 2 * WINDOW
    result["window_coverage"] = {
        "window_full_width_bp": win_full,
        "frac_reads_span_ge_full_window": float(np.mean(ref_spans >= win_full)),
        "frac_reads_span_ge_half_window": float(np.mean(ref_spans >= WINDOW)),
        "median_windows_covered_per_read": float(np.median(ref_spans) / win_full),
        "note": ("a read of ref span L can in principle inform the +/-5kb methylation "
                 "context of a sSNV near its center; ratio L/window indicates how many "
                 "full windows a read spans"),
    }

    bam.close()

    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=2)
    sys.stderr.write(f"[ssnv_spacing_coverage] wrote {args.out}\n")


if __name__ == "__main__":
    main()
