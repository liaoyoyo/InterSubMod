#!/usr/bin/env python3
"""[P2 全基因組 sSNV 對連鎖地景] 枚舉所有 ≤MAXD 的 TP somatic SNV 對 → tumor BAM per-read 共現 2×2 →
分類 nested(克隆階層) / sibling(同HP互斥) / allelic(異HP互斥) / co_linked / independent + normal=REF 確認。
回答「nested 等克隆結構全基因組多普遍」。輸出 p2_genomewide_linkage.json。"""
import sys, json, gzip
from collections import Counter, defaultdict
from bisect import bisect_right
import pysam
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot/scripts")
import p2_snv_linkage as P2

A = P2.A
VD = P2.VD
TBAM = P2.TBAM
NBAM = P2.NBAM
MAXD = 8000
COREAD_MIN = 6


def load_snvs(chrom):
    out = []
    with gzip.open(f"{VD}/filtered_snv_tp_{chrom}.vcf.gz", "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.split("\t")
            out.append((int(p[1]), p[3], p[4]))
    return sorted(out)


def classify(rr, ra, ar, aa):
    if aa < 2:
        return "mutual_excl" if (ra + ar) >= 2 else "sparse"
    if ra == 0 and ar == 0:
        return "co_linked"
    if ar == 0:
        return "nested_a_in_b"
    if ra == 0:
        return "nested_b_in_a"
    return "independent"


def main():
    tb = pysam.AlignmentFile(TBAM, "rb")
    nb = pysam.AlignmentFile(NBAM, "rb")
    tally = Counter(); nested = []; sibling = []; n_pairs = 0; n_tested = 0
    for c in range(1, 23):
        chrom = f"chr{c}"
        snvs = load_snvs(chrom)
        pos = [s[0] for s in snvs]
        for i in range(len(snvs)):
            j = i + 1
            while j < len(snvs) and pos[j] - pos[i] <= MAXD:
                n_pairs += 1
                a, b = snvs[i], snvs[j]
                tcalls, thp = P2.per_read_allele(tb, chrom, [a, b])
                multi = {rn: c2 for rn, c2 in tcalls.items() if a[0] in c2 and b[0] in c2
                         and c2[a[0]] in ("REF", "ALT") and c2[b[0]] in ("REF", "ALT")}
                if len(multi) >= COREAD_MIN:
                    n_tested += 1
                    tab = Counter((c2[a[0]], c2[b[0]]) for c2 in multi.values())
                    rr, ra, ar, aa = tab[("REF", "REF")], tab[("REF", "ALT")], tab[("ALT", "REF")], tab[("ALT", "ALT")]
                    rel = classify(rr, ra, ar, aa)
                    # HP 一致性(只對 ALT reads)
                    hpa = Counter(thp[rn] for rn, c2 in tcalls.items() if c2.get(a[0]) == "ALT" and rn in thp)
                    hpb = Counter(thp[rn] for rn, c2 in tcalls.items() if c2.get(b[0]) == "ALT" and rn in thp)
                    da = hpa.most_common(1)[0][0] if hpa else None
                    db = hpb.most_common(1)[0][0] if hpb else None
                    same_hp = (da == db and da is not None)
                    rec = {"chrom": chrom, "a": a[0], "b": b[0], "dist": b[0] - a[0], "coread": len(multi),
                           "cells": {"RR": rr, "RA": ra, "AR": ar, "AA": aa}, "rel": rel, "hp_a": da, "hp_b": db, "same_hp": same_hp}
                    if rel in ("nested_a_in_b", "nested_b_in_a"):
                        tally["nested" + ("_sameHP" if same_hp else "_diffHP")] += 1
                        if same_hp:
                            nested.append(rec)
                    elif rel == "mutual_excl":
                        tally["sibling" if same_hp else "allelic"] += 1
                        if same_hp:
                            sibling.append(rec)
                    else:
                        tally[rel] += 1
                j += 1
        print(f"[{chrom}] pairs={n_pairs} tested={n_tested} tally={dict(tally)}", flush=True)
    tb.close(); nb.close()
    out = {"MAXD": MAXD, "n_pairs_within": n_pairs, "n_tested_coread>=6": n_tested, "tally": dict(tally),
           "nested_sameHP_examples": sorted(nested, key=lambda r: -r["coread"])[:40],
           "sibling_sameHP_examples": sorted(sibling, key=lambda r: -r["coread"])[:40]}
    json.dump(out, open(f"{A}/p2_genomewide_linkage.json", "w"), ensure_ascii=False, indent=1)
    print("DONE tally=" + json.dumps(out["tally"], ensure_ascii=False), flush=True)
    print(f"nested_sameHP={len(nested)} sibling_sameHP={len(sibling)}", flush=True)


if __name__ == "__main__":
    main()
