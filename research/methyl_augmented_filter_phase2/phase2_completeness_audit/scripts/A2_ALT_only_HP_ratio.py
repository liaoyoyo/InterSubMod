#!/usr/bin/env python3
"""
A2: ALT-only HP=1:HP=2 ratio scan across HCC1395 baseline + HCC1395 V6 + V6 5-sample extension.

Background:
  PI 4-29 report (20260429_PI_Report_4_29_Errata_01.md) reported HCC1395 baseline
  HP=1:HP=2 = 17.3:1 (ALT-supporting reads, per-SNV-locus filter).
  20260519 ISM TPFP comparison report independently reported 1.696:1 (all primary reads,
  chr8+chr19 scope, per A1 methodology).
  These two numbers are NOT in conflict -- they measure two different things:
    - PI 4-29: "Of reads that carry the ALT allele at PASS SNV positions, what is HP1:HP2 ratio?"
      This filters to variant-supporting reads only, where Layer 1.5 priority bias matters most.
    - 20260519 / A1: "Of all primary alignments in chr8+chr19, what is HP1:HP2 ratio?"
      This includes both REF and ALT reads, dominated by REF/wild-type background.

  A2 reproduces the PI 4-29 cohort (ALT-only) and reports both cohorts side-by-side for
  5 samples to verify (a) baseline HCC1395 ALT-only ratio matches PI 4-29 ~17.3:1,
  (b) V6 attenuates the ratio, (c) 4 V6 extension samples behave consistently.

Method (ALT-only filter):
  For each ClairS-TO PASS SNV (chr, pos, ref, alt):
    fetch reads covering pos via pysam pileup or fetch+CIGAR walk
    determine the base each read carries at pos
    if base == ALT -> count HP tag toward ALT-only counter
    if base == REF -> count toward REF-only counter (audit, not reported)
  HP1:HP2 ALT-only ratio = sum(HP1 in ALT-supporting reads) / sum(HP2 in ALT-supporting reads).

Scope:
  chr8 + chr19 only (matches A1 cost envelope; chr8 LOH+HPSig 7.4x FP enrichment,
  chr19 priority-bug-sensitive germline-rich region).
  Per-locus fetch is O(coverage) not O(chr length), so chr8+chr19 ~3.5k SNV HCC1395 /
  up to 12.5k SNV H2009 finish in 5-15 min/BAM.

VCF source (frozen 2026-05-21):
  - HCC1395 baseline + V6 -> /big7_disk/.../longphase-to-mod/output/baseline/tumor_phased.vcf
                              (LongPhase-TO output preserving ClairS-TO PASS SNV positions/ALT)
    NB: HCC1395 V6 has no own tumor_phased.vcf in v6_germline_absent_revert/;
        the SNV cohort is shared (same ClairS-TO call set; only haplotag changes).
  - HCC1937 V6 -> v6_5sample_extension/HCC1937/tumor_phased.vcf
  - HCC1954 V6 -> v6_5sample_extension/HCC1954/tumor_phased.vcf
  - H1437 V6   -> v6_5sample_extension/H1437/tumor_phased.vcf
  - H2009 V6   -> v6_5sample_extension/H2009/tumor_phased.vcf

Outputs:
  - TSV : A2_ALT_only_HP_ratio_5sample.tsv
      columns: sample | bam_version | scan_scope | n_SNV_positions |
               ALT_reads_HP1 | ALT_reads_HP2 | ALT_only_ratio | all_reads_ratio
  - PNG : figures/A2_ALT_only_vs_all_reads_ratio.png (grouped bar 6 rows x 2 cohorts)
  - JSON: A2_ALT_only_HP_ratio_5sample.json  (per-BAM raw counts)
"""

from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path
from collections import Counter

import pysam
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Six-value HP categories (per Errata 01) -- mirror of A1
HP_CATEGORIES = [0, 1, 2, 11, 21, 33]

SCAN_CHROMS = ["chr8", "chr19"]

# BAM + matching VCF source (one VCF per sample; HCC1395 baseline/V6 share VCF)
SAMPLES = [
    # (sample, bam_version, bam_path, vcf_path)
    ("HCC1395", "baseline",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_phased.vcf"),
    ("HCC1395", "V6",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_phased.vcf"),  # shared SNV set
    ("HCC1937", "V6",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/HCC1937/tumor_tagged.bam",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/HCC1937/tumor_phased.vcf"),
    ("HCC1954", "V6",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/HCC1954/tumor_tagged.bam",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/HCC1954/tumor_phased.vcf"),
    ("H1437", "V6",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/H1437/tumor_tagged.bam",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/H1437/tumor_phased.vcf"),
    ("H2009", "V6",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/H2009/tumor_tagged.bam",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/H2009/tumor_phased.vcf"),
]


def parse_pass_snv(vcf_path: str, chroms: list[str]) -> list[tuple]:
    """Read PASS biallelic SNV in given chromosomes. Returns list of (chrom, pos1-based, ref, alt).

    Skips non-PASS, non-SNV (indel), multi-allelic (comma in ALT).
    """
    keep_chroms = set(chroms)
    out = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _, ref, alt, _, filt = fields[:7]
            if chrom not in keep_chroms:
                continue
            if filt != "PASS":
                continue
            if "," in alt:
                continue  # multi-allelic ignored
            if len(ref) != 1 or len(alt) != 1:
                continue  # SNV only
            out.append((chrom, int(pos), ref.upper(), alt.upper()))
    return out


def scan_bam_alt(bam_path: str, snv_loci: list[tuple]) -> dict:
    """For each SNV locus, fetch reads, classify by base at pos, accumulate HP tag.

    Returns dict with:
      alt_hp_counts: Counter over HP values for ALT-supporting reads
      ref_hp_counts: Counter over HP values for REF-supporting reads (audit)
      all_hp_counts: Counter over HP values for ALL primary reads covering ANY SNV locus
                     (this differs from A1 'all reads in chr8+chr19' -- it is the locus-restricted
                      'all reads at PASS SNV positions' cohort; per-read may double-count across
                      nearby SNVs but the metric is reads-overlapping-PASS-SNV not unique reads)
      n_loci_processed, n_loci_no_coverage, n_reads_no_base, n_reads_total
    """
    alt_hp: Counter = Counter()
    ref_hp: Counter = Counter()
    all_hp: Counter = Counter()
    n_loci_processed = 0
    n_loci_no_coverage = 0
    n_reads_no_base = 0
    n_reads_total = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, pos1, ref, alt in snv_loci:
            pos0 = pos1 - 1  # 0-based
            try:
                reads = bam.fetch(chrom, pos0, pos0 + 1, multiple_iterators=False)
            except ValueError:
                continue
            locus_had_read = False
            for read in reads:
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                # Find query base at pos0 via aligned_pairs (None if deletion / no overlap)
                base = None
                # get_aligned_pairs(matches_only=False, with_seq=False) is O(read_len);
                # for efficiency use get_reference_positions(full_length=True) then map
                try:
                    refpos = read.get_reference_positions(full_length=True)
                    # refpos is list aligned to query length; find index where refpos == pos0
                    # NB: deletions yield None -> read has no base at pos0
                    idx = None
                    for q_idx, rp in enumerate(refpos):
                        if rp == pos0:
                            idx = q_idx
                            break
                    if idx is None:
                        n_reads_no_base += 1
                        continue
                    seq = read.query_sequence
                    if seq is None or idx >= len(seq):
                        n_reads_no_base += 1
                        continue
                    base = seq[idx].upper()
                except Exception:
                    n_reads_no_base += 1
                    continue

                try:
                    hp = read.get_tag("HP")
                except KeyError:
                    hp = -1
                all_hp[hp] += 1
                n_reads_total += 1
                locus_had_read = True
                if base == alt:
                    alt_hp[hp] += 1
                elif base == ref:
                    ref_hp[hp] += 1
                # else: third-base / N -> ignore (rare)
            if locus_had_read:
                n_loci_processed += 1
            else:
                n_loci_no_coverage += 1

    return {
        "alt_hp_counts": dict(alt_hp),
        "ref_hp_counts": dict(ref_hp),
        "all_hp_counts": dict(all_hp),
        "n_loci_processed": n_loci_processed,
        "n_loci_no_coverage": n_loci_no_coverage,
        "n_reads_no_base": n_reads_no_base,
        "n_reads_total": n_reads_total,
    }


def summarize(scan: dict, sample: str, version: str, scope: str, n_snv: int,
              bam_path: str, vcf_path: str, elapsed_sec: float) -> dict:
    alt_hp = scan["alt_hp_counts"]
    all_hp = scan["all_hp_counts"]

    alt_hp1 = alt_hp.get(1, 0)
    alt_hp2 = alt_hp.get(2, 0)
    all_hp1 = all_hp.get(1, 0)
    all_hp2 = all_hp.get(2, 0)

    alt_ratio = round(alt_hp1 / alt_hp2, 4) if alt_hp2 > 0 else float("nan")
    all_ratio = round(all_hp1 / all_hp2, 4) if all_hp2 > 0 else float("nan")

    return {
        "sample": sample,
        "bam_version": version,
        "scan_scope": scope,
        "n_SNV_positions": n_snv,
        "n_SNV_with_coverage": scan["n_loci_processed"],
        "ALT_reads_HP1": alt_hp1,
        "ALT_reads_HP2": alt_hp2,
        "ALT_reads_HP0": alt_hp.get(0, 0),
        "ALT_reads_HP11": alt_hp.get(11, 0),
        "ALT_reads_HP21": alt_hp.get(21, 0),
        "ALT_reads_HP33": alt_hp.get(33, 0),
        "ALT_reads_untagged": alt_hp.get(-1, 0),
        "ALT_reads_total": sum(alt_hp.values()),
        "ALL_reads_HP1": all_hp1,
        "ALL_reads_HP2": all_hp2,
        "ALL_reads_total": sum(all_hp.values()),
        "ALT_only_ratio": alt_ratio,
        "all_reads_ratio": all_ratio,
        "elapsed_sec": round(elapsed_sec, 1),
        "bam_path": bam_path,
        "vcf_path": vcf_path,
    }


def write_tsv(rows: list[dict], out_path: Path) -> None:
    cols = [
        "sample", "bam_version", "scan_scope", "n_SNV_positions", "n_SNV_with_coverage",
        "ALT_reads_HP1", "ALT_reads_HP2", "ALT_reads_HP0",
        "ALT_reads_HP11", "ALT_reads_HP21", "ALT_reads_HP33",
        "ALT_reads_untagged", "ALT_reads_total",
        "ALL_reads_HP1", "ALL_reads_HP2", "ALL_reads_total",
        "ALT_only_ratio", "all_reads_ratio",
        "elapsed_sec", "bam_path", "vcf_path",
    ]
    with open(out_path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def plot_grouped_bar(rows: list[dict], out_png: Path) -> None:
    """Grouped bar: each (sample, bam_version) shows two bars (ALT-only ratio vs all-reads ratio)."""
    labels = [f"{r['sample']}\n{r['bam_version']}" for r in rows]
    alt_ratios = np.array([r["ALT_only_ratio"] for r in rows], dtype=float)
    all_ratios = np.array([r["all_reads_ratio"] for r in rows], dtype=float)

    x = np.arange(len(rows))
    width = 0.38

    fig, ax = plt.subplots(figsize=(11, 6))
    b1 = ax.bar(x - width / 2, alt_ratios, width, color="#d62728",
                label="ALT-only HP1:HP2 (PI 4-29 cohort)", edgecolor="white")
    b2 = ax.bar(x + width / 2, all_ratios, width, color="#1f77b4",
                label="All-reads HP1:HP2 (20260519 / A1 cohort)", edgecolor="white")

    # PI 4-29 baseline reference line
    ax.axhline(17.3, color="#888", linestyle=":", linewidth=1.2, alpha=0.9)
    ax.text(len(rows) - 0.5, 17.3 + 0.4, "PI 4-29 baseline ALT-only = 17.3",
            color="#555", fontsize=9, ha="right", va="bottom")
    ax.axhline(1.0, color="#bbb", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.text(len(rows) - 0.5, 1.05, "balanced 1:1", color="#666", fontsize=8, ha="right")

    # Annotate each bar with its value
    for rects, vals in [(b1, alt_ratios), (b2, all_ratios)]:
        for rect, v in zip(rects, vals):
            if np.isfinite(v):
                ax.annotate(f"{v:.2f}", xy=(rect.get_x() + rect.get_width() / 2, v),
                            xytext=(0, 3), textcoords="offset points",
                            ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("HP1 : HP2 read ratio (log-friendly linear scale)")
    ax.set_title("A2: ALT-only vs all-reads HP1:HP2 ratio per (sample x BAM version)\n"
                 "Scope: chr8+chr19 PASS SNV (HCC1395 baseline+V6 / 4 V6 extension)")
    ax.legend(loc="upper left", fontsize=9)
    ymax = max(np.nanmax(alt_ratios) if np.any(np.isfinite(alt_ratios)) else 1.0, 18.0)
    ax.set_ylim(0, ymax * 1.1)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()


def _dump_partial(rows, raw, scope, out_tsv, out_json, out_png, t_start):
    """Incremental write: dump TSV + JSON + PNG so we can resume safely."""
    write_tsv(rows, out_tsv)
    with open(out_json, "w") as f:
        json.dump({
            "scan_scope": scope, "hp_categories": HP_CATEGORIES, "rows": raw,
            "total_elapsed_sec": round(time.time() - t_start, 1),
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "partial": True,
        }, f, indent=2)
    try:
        plot_grouped_bar(rows, out_png)
    except Exception as e:
        print(f"  [partial-png-skip] {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-png", required=True)
    parser.add_argument("--out-json", required=True,
                        help="Per-BAM raw counts JSON for audit trail")
    parser.add_argument("--skip-samples", default="",
                        help="Comma list of 'sample/version' to skip "
                             "(e.g. 'HCC1395/baseline,HCC1395/V6'). Useful for resume.")
    parser.add_argument("--seed-json", default=None,
                        help="Pre-existing partial JSON to merge into output (resume mode)")
    parser.add_argument("--chroms", default=",".join(SCAN_CHROMS),
                        help="Comma list of chromosomes (default chr8,chr19). "
                             "Use 'chr19' to dodge H2009 chr8 9532-SNV timeout.")
    args = parser.parse_args()
    skip_set = {s.strip() for s in args.skip_samples.split(",") if s.strip()}
    # Override SCAN_CHROMS if user passed a subset
    new_chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    if new_chroms != SCAN_CHROMS:
        SCAN_CHROMS.clear()
        SCAN_CHROMS.extend(new_chroms)
        print(f"[chroms] scan scope overridden -> {SCAN_CHROMS}", flush=True)

    out_tsv = Path(args.out_tsv)
    out_png = Path(args.out_png)
    out_json = Path(args.out_json)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)

    scope = "+".join(SCAN_CHROMS)
    rows: list[dict] = []
    raw: list[dict] = []
    t_start = time.time()

    # Seed from prior partial JSON (resume mode)
    if args.seed_json and Path(args.seed_json).exists():
        print(f"[seed] loading prior rows from {args.seed_json}", flush=True)
        with open(args.seed_json) as f:
            seed = json.load(f)
        for r in seed.get("rows", []):
            # Reconstruct summary row from raw counts (so TSV format matches)
            scan = {
                "alt_hp_counts": {int(k): v for k, v in r["alt_hp_counts"].items()},
                "ref_hp_counts": {int(k): v for k, v in r["ref_hp_counts"].items()},
                "all_hp_counts": {int(k): v for k, v in r["all_hp_counts"].items()},
                "n_loci_processed": r["n_loci_processed"],
                "n_loci_no_coverage": r["n_loci_no_coverage"],
                "n_reads_no_base": r["n_reads_no_base"],
                "n_reads_total": r["n_reads_total"],
            }
            summary = summarize(scan, r["sample"], r["version"], r["scan_scope"],
                                r["n_SNV_positions"], r["bam_path"], r["vcf_path"],
                                r["elapsed_sec"])
            rows.append(summary)
            raw.append(r)
            skip_set.add(f"{r['sample']}/{r['version']}")
            print(f"  [seed] loaded {r['sample']}/{r['version']}", flush=True)

    # Cache parsed VCF (HCC1395 baseline + V6 share)
    vcf_cache: dict[str, list[tuple]] = {}

    for sample, version, bam_path, vcf_path in SAMPLES:
        key = f"{sample}/{version}"
        if key in skip_set:
            print(f"[skip] {key} (in skip-set or seeded)", flush=True)
            continue
        print(f"[scan] {sample}/{version} bam={bam_path}", flush=True)
        if not Path(bam_path).exists():
            print(f"  [SKIP] BAM not found", file=sys.stderr)
            continue
        if not Path(vcf_path).exists():
            print(f"  [SKIP] VCF not found: {vcf_path}", file=sys.stderr)
            continue

        if vcf_path not in vcf_cache:
            print(f"  [vcf] parsing {vcf_path}", flush=True)
            vcf_cache[vcf_path] = parse_pass_snv(vcf_path, SCAN_CHROMS)
            print(f"  [vcf] n_PASS_SNV = {len(vcf_cache[vcf_path])} ({scope})", flush=True)

        snv_loci = vcf_cache[vcf_path]
        t0 = time.time()
        scan = scan_bam_alt(bam_path, snv_loci)
        elapsed = time.time() - t0
        summary = summarize(scan, sample, version, scope, len(snv_loci),
                            bam_path, vcf_path, elapsed)
        rows.append(summary)
        raw.append({
            "sample": sample, "version": version, "bam_path": bam_path,
            "vcf_path": vcf_path, "scan_scope": scope,
            "n_SNV_positions": len(snv_loci),
            "alt_hp_counts": scan["alt_hp_counts"],
            "ref_hp_counts": scan["ref_hp_counts"],
            "all_hp_counts": scan["all_hp_counts"],
            "n_loci_processed": scan["n_loci_processed"],
            "n_loci_no_coverage": scan["n_loci_no_coverage"],
            "n_reads_no_base": scan["n_reads_no_base"],
            "n_reads_total": scan["n_reads_total"],
            "elapsed_sec": elapsed,
        })
        print(f"  -> n_SNV={len(snv_loci):,} ALT_HP1={summary['ALT_reads_HP1']:,} "
              f"ALT_HP2={summary['ALT_reads_HP2']:,} ALT_ratio={summary['ALT_only_ratio']} "
              f"ALL_ratio={summary['all_reads_ratio']} elapsed={elapsed:.1f}s",
              flush=True)
        # Incremental dump after each BAM to survive kills
        _dump_partial(rows, raw, scope, out_tsv, out_json, out_png, t_start)
        print(f"  [partial-dump] rows={len(rows)} -> {out_tsv.name}", flush=True)

    write_tsv(rows, out_tsv)
    with open(out_json, "w") as f:
        json.dump({
            "scan_scope": scope, "hp_categories": HP_CATEGORIES, "rows": raw,
            "total_elapsed_sec": round(time.time() - t_start, 1),
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        }, f, indent=2)
    plot_grouped_bar(rows, out_png)
    print(f"\n[done] TSV  -> {out_tsv}")
    print(f"[done] PNG  -> {out_png}")
    print(f"[done] JSON -> {out_json}")
    print(f"[done] total elapsed = {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
