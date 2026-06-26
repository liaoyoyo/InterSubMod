"""[S5 SEQC2 truth overlay] 把 SEQC2 真值疊到我的 TP∪FP union census。
每 sSNV → {in_highconf, in_superset, in_hc_region} → truth_class:
  highconf            : 在 HighConf sSNV (金標準確認 somatic)
  superset_only       : 在 superSet 但非 HighConf (較寬真值)
  in_HC_not_called    : 在 HC callable region 內但 SEQC2 未收 (γ 類 — region 可評估、ONT 偵測、SEQC2 沒叫)
  outside_HC          : 不在 HC region (unevaluable)
交叉 my source(TP/FP) × truth_class × linked-somatic。輸出 sm_seqc2_overlay.json。
"""
import json
import gzip
from bisect import bisect_right
from collections import Counter, defaultdict

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
S = "/big8_disk/data/HCC1395/SEQC2"
IN = f"{A}/sm_linkage_genomewide.json"
LED = f"{A}/sm_completeness_ledger.json"
OUT = f"{A}/sm_seqc2_overlay.json"


def load_vcf_set(path):
    s = set()
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.split("\t")
            s.add((p[0], int(p[1])))
    return s


def load_hc_regions(path):
    starts = defaultdict(list)
    ends = defaultdict(list)
    with open(path) as fh:
        for ln in fh:
            p = ln.split("\t")
            if len(p) < 3:
                continue
            c = p[0]
            starts[c].append(int(p[1]))
            ends[c].append(int(p[2]))
    for c in starts:
        z = sorted(zip(starts[c], ends[c]))
        starts[c] = [a for a, _ in z]
        ends[c] = [b for _, b in z]
    return starts, ends


def in_region(starts, ends, chrom, pos):
    st = starts.get(chrom)
    if not st:
        return False
    i = bisect_right(st, pos) - 1
    return i >= 0 and ends[chrom][i] >= pos  # BED 0-based half-open; pos 1-based ≈ within


def main():
    d = json.load(open(IN))
    census = d["census"]
    hc = load_vcf_set(f"{S}/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz")
    ss = load_vcf_set(f"{S}/sSNV.MSDUKT.superSet.v1.2.vcf.gz")
    rs, re = load_hc_regions(f"{S}/High-Confidence_Regions_v1.2.bed")

    by_class = Counter()
    src_x_class = Counter()
    linked_somatic_x_class = Counter()
    gamma_like = []   # FP/in_HC_not_called + linked + somatic = γ-類 recovered
    for key, c in census.items():
        chrom, pos = c["chrom"], c["pos"]
        ih = (chrom, pos) in hc
        iss = (chrom, pos) in ss
        ireg = in_region(rs, re, chrom, pos)
        if ih:
            cls = "highconf"
        elif iss:
            cls = "superset_only"
        elif ireg:
            cls = "in_HC_not_called"
        else:
            cls = "outside_HC"
        by_class[cls] += 1
        src_x_class[(c["src"], cls)] += 1
        linked = (c["n_realized_links"] >= 1)
        if linked and c.get("somatic"):
            linked_somatic_x_class[cls] += 1
            if cls == "in_HC_not_called":
                gamma_like.append({"locus": key, "src": c["src"], "vaf": c.get("vaf"),
                                   "normal": [c.get("normal_ref"), c.get("normal_alt")],
                                   "max_coread": c.get("max_coread")})

    out = {
        "seqc2_counts": {"highconf": len(hc), "superset": len(ss)},
        "universe_by_truth_class": dict(by_class),
        "source_x_truth_class": {f"{k[0]}|{k[1]}": v for k, v in src_x_class.items()},
        "linked_somatic_x_truth_class": dict(linked_somatic_x_class),
        "gamma_like_in_HC_not_called_linked_somatic": {
            "count": len(gamma_like),
            "note": "在 HC callable region 內、SEQC2 未收、但 ONT normal-confirmed somatic 且單分子連鎖 = γ 類",
            "examples": sorted(gamma_like, key=lambda x: -(x["max_coread"] or 0))[:30],
        },
    }
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(f"SEQC2-OVERLAY by_class={dict(by_class)} linked_somatic_x_class={dict(linked_somatic_x_class)} "
          f"gamma_like={len(gamma_like)} -> {OUT}")


if __name__ == "__main__":
    main()
