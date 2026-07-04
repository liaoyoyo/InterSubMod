#!/usr/bin/env python3
"""Convert a SAVANA segmented absolute copy-number TSV into the 4-column
gain/loss/loh BED consumed by sm_region_integration.py (env SM_CNBED).

Downstream contract (sm_region_integration.py load_cn/cn_state):
  - whitespace-split, >=4 cols: chrom  start(int)  end(int)  label
  - cn_state(pos) returns the label of the first overlapping interval, else
    "neutral" (SM_CN_SPARSE=0) / "unknown" (SM_CN_SPARSE=1).
  - labels that matter downstream: {gain, loss, loh, neutral}; loh+neutral are
    counted as "clean" (topology-reliable) in shape_clean.

Label decision tree (one label per SAVANA segment -> non-overlapping BED, so
cn_state is unambiguous; cleaner than SEQC2's overlapping multi-label BED):

    copyNumber            >= --gain-cn (2.5)  -> gain   (amplified: multiplicity ambiguous)
    copyNumber            <= --loss-cn (1.5)  -> loss   (deletion)
    minorAlleleCopyNumber <  --loh-minor(0.5) -> loh    (copy-neutral allelic imbalance)
    otherwise                                 -> neutral (NOT emitted; cn_state default)

Thresholds default to the HCC1395-validated feasibility values (LOH minorCN<0.5
reproduced SEQC2 loh Jaccard 0.962). gain/loss are diploid-anchored direction
calls (weaker under aneuploid baseline, per feasibility doc) -> tunable via CLI.

Run downstream with SM_CN_SPARSE=0 (SAVANA is genome-wide, non-sparse).
"""
import argparse
import sys

DEFAULT_CHROMS = [f"chr{c}" for c in range(1, 23)] + ["chrX"]


def parse_float(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def convert(tsv_path, gain_cn, loss_cn, loh_minor, chroms, min_bins):
    rows = []
    counts = {"gain": 0, "loss": 0, "loh": 0, "neutral": 0}
    bp = {"gain": 0, "loss": 0, "loh": 0, "neutral": 0}
    skipped_na = 0
    chrom_set = set(chroms) if chroms else None

    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_chrom = header.index("chromosome")
            i_start = header.index("start")
            i_end = header.index("end")
            i_cn = header.index("copyNumber")
            i_minor = header.index("minorAlleleCopyNumber")
        except ValueError as e:
            sys.exit(f"ERROR: SAVANA header missing expected column: {e}\nheader={header}")
        i_bins = header.index("bin_count") if "bin_count" in header else None

        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) <= i_minor:
                continue
            chrom = p[i_chrom]
            if chrom_set is not None and chrom not in chrom_set:
                continue
            if not p[i_start].isdigit() or not p[i_end].isdigit():
                continue
            start, end = int(p[i_start]), int(p[i_end])
            if min_bins and i_bins is not None:
                bc = parse_float(p[i_bins])
                if bc is None or bc < min_bins:
                    continue
            cn = parse_float(p[i_cn])
            minor = parse_float(p[i_minor])
            if cn is None:
                skipped_na += 1
                continue

            if cn >= gain_cn:
                label = "gain"
            elif cn <= loss_cn:
                label = "loss"
            elif minor is not None and minor < loh_minor:
                label = "loh"
            else:
                label = "neutral"

            counts[label] += 1
            bp[label] += end - start
            if label != "neutral":                       # neutral -> cn_state default, not emitted
                rows.append((chrom, start, end, label))

    rows.sort(key=lambda r: (r[0], r[1]))
    return rows, counts, bp, skipped_na


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", help="SAVANA *_segmented_absolute_copy_number.tsv")
    ap.add_argument("-o", "--out", required=True, help="output SM_CNBED (.bed)")
    ap.add_argument("--gain-cn", type=float, default=2.5)
    ap.add_argument("--loss-cn", type=float, default=1.5)
    ap.add_argument("--loh-minor", type=float, default=0.5)
    ap.add_argument("--min-bins", type=float, default=0,
                    help="drop segments with bin_count < this (0=keep all)")
    ap.add_argument("--chroms", default=",".join(DEFAULT_CHROMS),
                    help="comma-separated chroms to keep (default autosomes+chrX); "
                         "'all' keeps every contig")
    args = ap.parse_args()

    chroms = None if args.chroms == "all" else args.chroms.split(",")
    rows, counts, bp, skipped_na = convert(
        args.tsv, args.gain_cn, args.loss_cn, args.loh_minor, chroms, args.min_bins)

    with open(args.out, "w") as fh:
        for chrom, start, end, label in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{label}\n")

    total = sum(counts.values())
    print(f"input : {args.tsv}")
    print(f"output: {args.out}  ({len(rows)} gain/loss/loh intervals)")
    print(f"thresholds: gain_cn>={args.gain_cn} loss_cn<={args.loss_cn} "
          f"loh_minor<{args.loh_minor} min_bins={args.min_bins}")
    print(f"segments classified (n={total}, skipped_NA={skipped_na}):")
    for k in ("gain", "loss", "loh", "neutral"):
        mb = bp[k] / 1e6
        print(f"  {k:8s} {counts[k]:>6}  ({mb:,.1f} Mb)"
              + ("  [not emitted; cn_state default]" if k == "neutral" else ""))


if __name__ == "__main__":
    main()
