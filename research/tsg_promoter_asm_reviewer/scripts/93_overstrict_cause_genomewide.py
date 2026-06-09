#!/usr/bin/env python3
"""
93 - Genome-wide confirmation of the over-strict cause: read RAW CramersV from the
over-strict raw re-run (script 92, clustering/significance.json) for ALL 2788 TP +
131 FP over-strict loci, and split:
  reliability-gated : raw CramersV_max >= 0.3 but gated to <0.1 (sparse, Cochran)
  partition-mismatch: raw CramersV_max < 0.1 (even raw chi2 weak; PERMANOVA-only)
  mid               : 0.1-0.3
Confirms whether the 420/425 curated-subset finding holds genome-wide.

Output: display_v2/cause_genomewide.json
"""
import json, glob, os, csv

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
RAW = "/big7_disk/liaoyoyo2001/ism_overstrict_raw"
DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")


def raw_max(cls, chrom, pos):
    for d in glob.glob(f"{RAW}/HCC1395_{cls}/overstrict_{cls}/{chrom}/{chrom}_{pos}/*"):
        j = f"{d}/clustering/significance.json"
        if os.path.exists(j):
            s = json.load(open(j))
            return max(s.get("global_" + a, {}).get("cramers_v", 0)
                       for a in ("alt", "hp", "hp_family", "hp_fine"))
    return None


def main():
    res = {}
    for cls in ("tp", "fp"):
        pos = [tuple(l.split("\t")) for l in open(f"{DV}/overstrict_{cls}_pos.txt").read().split("\n") if l.strip()]
        rel = mis = mid = miss = 0
        rmlist = []
        for chrom, p in pos:
            rm = raw_max(cls, chrom, p)
            if rm is None:
                miss += 1; continue
            rmlist.append(rm)
            if rm >= 0.3:
                rel += 1
            elif rm < 0.1:
                mis += 1
            else:
                mid += 1
        n = len(pos)
        res[cls] = dict(n=n, found=n - miss, reliability_gated=rel, partition_mismatch=mis,
                        mid=mid, not_found=miss,
                        pct_reliability=round(100 * rel / max(n - miss, 1), 1),
                        median_rawCV=round(sorted(rmlist)[len(rmlist) // 2], 3) if rmlist else None)
    with open(f"{DV}/cause_genomewide.json", "w") as f:
        json.dump(res, f, indent=2)
    for cls in ("tp", "fp"):
        r = res[cls]
        print(f"{cls.upper()} over-strict n={r['n']} (found raw {r['found']}): "
              f"reliability-gated(raw≥0.3)={r['reliability_gated']} ({r['pct_reliability']}%)  "
              f"partition-mismatch(raw<0.1)={r['partition_mismatch']}  mid={r['mid']}  "
              f"median_rawCV={r['median_rawCV']}")
    print(f"\n[93] wrote {DV}/cause_genomewide.json")


if __name__ == "__main__":
    main()
