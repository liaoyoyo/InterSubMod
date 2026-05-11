#!/usr/bin/env python3
"""
Analyze per-chr HP tag distributions across baseline/V5/V6/paired_T.

Aggregates 4 BAMs x 8 chrs.
Maps HP tags to canonical groups:
  HP1 group = {1, 11, "1", "1-1"}
  HP2 group = {2, 21, "2", "2-1"}
  HP3/somatic = {33, "3", 3}
  untagged   = {NA, "0", 0}

Computes per-chr per-version:
  - HP1_grp / HP2_grp ratio
  - L1 distance to paired_T's HP1/HP2 proportions
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/paired_priority_bug_audit/v6_quantification_evidence")
IN_FILE = OUT_DIR / "per_chr_hp_all.tsv"


def categorize(hp):
    s = str(hp).strip()
    if s in ("1", "11", "1-1"):
        return "HP1_grp"
    if s in ("2", "21", "2-1"):
        return "HP2_grp"
    if s in ("3", "33"):
        return "HP3_amb"
    if s in ("NA", "0", "."):
        return "untagged"
    return f"other:{s}"


def main():
    if not IN_FILE.exists():
        sys.exit(f"missing input: {IN_FILE}")

    # data[(version, chr)][category] = count
    data = defaultdict(lambda: defaultdict(int))
    with IN_FILE.open() as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            v = row["version"]
            c = row["chr"]
            try:
                n = int(row["count"])
            except ValueError:
                continue
            cat = categorize(row["hp_tag"])
            data[(v, c)][cat] += n

    versions = ["baseline", "V5", "V6", "paired_T"]
    chrs = sorted({c for v, c in data.keys()},
                  key=lambda x: (0 if x == "chrX" else int(x.replace("chr", "")), x))
    # standard order
    chrs = ["chr1", "chr3", "chr5", "chr8", "chr11", "chr17", "chr19", "chrX"]

    # Write per-chr per-version table
    out_tsv = OUT_DIR / "per_chr_hp_summary.tsv"
    with out_tsv.open("w") as f:
        f.write("chr\tversion\tHP1_grp\tHP2_grp\tHP3_amb\tuntagged\ttotal\tHP1_ratio\tHP1_to_HP2_ratio\n")
        for c in chrs:
            for v in versions:
                d = data.get((v, c), {})
                h1 = d.get("HP1_grp", 0)
                h2 = d.get("HP2_grp", 0)
                h3 = d.get("HP3_amb", 0)
                u = d.get("untagged", 0)
                total = h1 + h2 + h3 + u
                hp1_ratio = h1 / max(total, 1)
                ratio = h1 / max(h2, 1)
                f.write(f"{c}\t{v}\t{h1}\t{h2}\t{h3}\t{u}\t{total}\t{hp1_ratio:.4f}\t{ratio:.4f}\n")

    # Compute L1 distance to paired_T per chr per version
    # Use HP1/HP2 proportions among tagged HP1/HP2 reads (exclude HP3/untagged)
    out_dist = OUT_DIR / "per_chr_distance_to_paired.tsv"
    with out_dist.open("w") as f:
        f.write("chr\tversion\tHP1_prop\tpaired_HP1_prop\tL1_distance_to_paired\tHP1_to_HP2_ratio\tpaired_ratio\n")
        for c in chrs:
            d_p = data.get(("paired_T", c), {})
            ph1 = d_p.get("HP1_grp", 0)
            ph2 = d_p.get("HP2_grp", 0)
            ptot = ph1 + ph2
            paired_hp1_prop = ph1 / max(ptot, 1)
            paired_ratio = ph1 / max(ph2, 1)
            for v in versions:
                d = data.get((v, c), {})
                h1 = d.get("HP1_grp", 0)
                h2 = d.get("HP2_grp", 0)
                tot = h1 + h2
                hp1_prop = h1 / max(tot, 1)
                # L1 distance on HP1/HP2 distribution (2-vector)
                l1 = abs(hp1_prop - paired_hp1_prop) + abs((1 - hp1_prop) - (1 - paired_hp1_prop))
                ratio = h1 / max(h2, 1)
                f.write(f"{c}\t{v}\t{hp1_prop:.4f}\t{paired_hp1_prop:.4f}\t{l1:.4f}\t{ratio:.4f}\t{paired_ratio:.4f}\n")

    # Compute genome-wide aggregate (across all 8 chrs)
    out_genome = OUT_DIR / "genome_wide_hp_summary.tsv"
    with out_genome.open("w") as f:
        f.write("version\tHP1_grp\tHP2_grp\tHP3_amb\tuntagged\ttotal\tHP1_to_HP2_ratio\tHP1_prop_among_tagged\tL1_distance_to_paired\n")
        # paired ref
        ph1_tot = sum(data[(("paired_T", c))].get("HP1_grp", 0) for c in chrs)
        ph2_tot = sum(data[(("paired_T", c))].get("HP2_grp", 0) for c in chrs)
        paired_hp1_prop = ph1_tot / max(ph1_tot + ph2_tot, 1)
        for v in versions:
            h1 = sum(data[((v, c))].get("HP1_grp", 0) for c in chrs)
            h2 = sum(data[((v, c))].get("HP2_grp", 0) for c in chrs)
            h3 = sum(data[((v, c))].get("HP3_amb", 0) for c in chrs)
            u = sum(data[((v, c))].get("untagged", 0) for c in chrs)
            total = h1 + h2 + h3 + u
            ratio = h1 / max(h2, 1)
            hp1_prop = h1 / max(h1 + h2, 1)
            l1 = abs(hp1_prop - paired_hp1_prop) + abs((1 - hp1_prop) - (1 - paired_hp1_prop))
            f.write(f"{v}\t{h1}\t{h2}\t{h3}\t{u}\t{total}\t{ratio:.4f}\t{hp1_prop:.4f}\t{l1:.4f}\n")

    print(f"WRITE {out_tsv}")
    print(f"WRITE {out_dist}")
    print(f"WRITE {out_genome}")

    # Print summary to stdout
    print()
    print("=" * 90)
    print("GENOME-WIDE (8 representative chrs combined)")
    print("=" * 90)
    print(open(out_genome).read())

    print()
    print("=" * 90)
    print("PER-CHR L1 distance to paired_T (HP1/HP2 distribution)")
    print("=" * 90)
    print(open(out_dist).read())


if __name__ == "__main__":
    main()
