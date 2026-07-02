#!/usr/bin/env python3
"""[P2 multi-sSNV feasibility] 候選 anchor(=pos+5000) 周邊 TP somatic SNV 密度 → 哪些候選有 ≥2 sSNV 可連鎖
(同 read-length / window)。9 tumor-specific + 30 anchored + 全 1139 context。純讀 VCF + JSON。
輸出 p2_feasibility.json + 各候選鄰近 sSNV 清單(供 linkage)。"""
import json, gzip, os
import numpy as np
from bisect import bisect_left, bisect_right

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
VD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"

# ---- load all TP somatic SNV positions per chrom ----
SNV = {}
for c in range(1, 23):
    ch = f"chr{c}"; f = f"{VD}/filtered_snv_tp_{ch}.vcf.gz"
    if not os.path.exists(f):
        continue
    pos = []
    with gzip.open(f, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.split("\t")
            pos.append(int(p[1]))
    SNV[ch] = sorted(pos)
total = sum(len(v) for v in SNV.values())
print(f"TP somatic SNV loaded: {total}")

THR = [1000, 5000, 10000, 20000, 50000]


def neighbors(ch, anchor):
    arr = SNV.get(ch, [])
    out = {}
    for t in THR:
        lo = bisect_left(arr, anchor - t); hi = bisect_right(arr, anchor + t)
        # 排除 anchor 自身(±2bp 容差)
        n = sum(1 for p in arr[lo:hi] if abs(p - anchor) > 2)
        out[t] = n
    # 最近鄰距離 + 該距離內所有 sSNV pos (50kb)
    lo = bisect_left(arr, anchor - 50000); hi = bisect_right(arr, anchor + 50000)
    near = sorted([p for p in arr[lo:hi] if abs(p - anchor) > 2], key=lambda p: abs(p - anchor))
    nn = abs(near[0] - anchor) if near else None
    return out, nn, near[:6]


cands = json.load(open(f"{A}/cis_candidates_resolved.json"))
nine = {f"{r['chrom']}:{r['pos']}" for r in json.load(open(f"{A}/anchored_retest.json"))["still_candidate_on_geno_axis"]}
anch = {f"{r['chrom']}:{r['pos']}" for r in json.load(open(f"{A}/survivor_census.json"))["geno_anchored_loci"]}
# nact verdict 標記
nv = {f"{r['chrom']}:{r['pos']}": r["nact_verdict"] for r in json.load(open(f"{A}/nact_results.json"))}

rows = []
for c in cands:
    key = f"{c['chrom']}:{c['pos']}"
    anchor = int(c["pos"]) + 5000
    nb, nn, near = neighbors(c["chrom"], anchor)
    rows.append({"key": key, "set": c["set"], "anchor": anchor, "nact": nv.get(key),
                 "in9": key in nine, "in30": key in anch, "nbr": nb, "nearest": nn, "near_snvs": near})


def summarize(subset, name):
    s = {"n": len(subset)}
    for t in THR:
        s[f"has_nbr_{t}"] = sum(1 for r in subset if r["nbr"][t] >= 1)
    s["nearest_median"] = int(np.median([r["nearest"] for r in subset if r["nearest"] is not None])) if any(r["nearest"] for r in subset) else None
    return s


summ = {"total_tp_snv": total,
        "ALL_1139": summarize(rows, "all"),
        "nine_tumorspecific": summarize([r for r in rows if r["in9"]], "9"),
        "thirty_anchored": summarize([r for r in rows if r["in30"]], "30"),
        "candidate_102": summarize([r for r in rows if r["nact"] == "candidate_subclone"], "102"),
        "linkable_9": [{"key": r["key"], "set": r["set"], "nearest": r["nearest"], "n_10k": r["nbr"][10000], "near_snvs": r["near_snvs"]}
                       for r in rows if r["in9"] and r["nbr"][10000] >= 1],
        "linkable_30": [{"key": r["key"], "set": r["set"], "nearest": r["nearest"], "n_10k": r["nbr"][10000], "near_snvs": r["near_snvs"]}
                        for r in rows if r["in30"] and r["nbr"][10000] >= 1]}
json.dump({"summary": summ, "per_candidate": rows}, open(f"{A}/p2_feasibility.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(summ, ensure_ascii=False, indent=1))
print("[-> p2_feasibility.json]")
