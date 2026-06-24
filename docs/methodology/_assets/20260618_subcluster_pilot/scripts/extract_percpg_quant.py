#!/usr/bin/env python3
"""[差異甲基定位點量化] 全 34,736 走 region dirs, 每軸累積差異 CpG 分類(方向 hyper/hypo + |Δβ| bins),
每位點算 union 差異 CpG 數(跨軸去重) + per-axis n_sig → histogram 資料。逐 CpG 只數量+方向+強度(不存座標)。
輸出 percpg_cpg_classification.json(per軸全域分類) + percpg_per_locus.json(每位點, histogram 用)。"""
import os, glob, json, sys
from collections import Counter, defaultdict
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__)); import per_cpg_diff as P

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
KEEP = f"{WT}/output/_phylo_wg_full"
BINS = [(0.2, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 1.01)]; BINLAB = [".2-.3", ".3-.5", ".5-.7", ".7+"]
rec = {f"{r['chrom']}:{int(r['pos'])}": r["set"] for r in json.load(open(f"{A}/phylo_cpp_wg_full_records_v6.json"))}


def work(region):
    try:
        M, cpgs = P.load_meth(region)
        if not M:
            return None
        base = os.path.basename(region).split("_"); key = f"{base[0]}:{int(base[1])}"
        setl = rec.get(key, "?")
        gs = P.groupings(P.reads_meta(region))
        union = set(); axis_nsig = {}; per_axis_cls = {}
        for name, g in gs.items():
            tg = P.test_grouping(M, cpgs, g)
            union |= tg["sigset"]; axis_nsig[name] = tg["n_sig_mwu"]
            sd = tg["sig_dbetas"]
            per_axis_cls[name] = {"n_sig": len(sd), "n_hyper": sum(1 for x in sd if x > 0), "n_hypo": sum(1 for x in sd if x < 0),
                                  "bins": [sum(1 for x in sd if lo <= abs(x) < hi) for lo, hi in BINS]}
        # somatic-marker = HP1/HP1-1 或 HP2/HP2-1; subclone-marker = cluster_split
        som = max(axis_nsig.get("HP1_vs_HP1-1", 0), axis_nsig.get("HP2_vs_HP2-1", 0))
        sub = axis_nsig.get("cluster_split", 0)
        return (key, {"set": setl, "union": len(union), "axis_nsig": axis_nsig, "som_marker": som, "sub_marker": sub}, per_axis_cls)
    except Exception:
        return None


def main():
    regions = [os.path.dirname(os.path.dirname(p)) for p in glob.glob(f"{KEEP}/**/clustering/phylo_groups.tsv", recursive=True)]
    print(f"regions: {len(regions)}", flush=True)
    with Pool(24) as p:
        res = [r for r in p.map(work, regions, chunksize=20) if r]
    # per-locus
    per_locus = {k: v for k, v, _ in res}
    json.dump(per_locus, open(f"{A}/percpg_per_locus.json", "w"))
    # 全域 per-axis 分類累積
    glob_axis = defaultdict(lambda: {"n_sig": 0, "n_hyper": 0, "n_hypo": 0, "bins": [0, 0, 0, 0], "n_loci": 0})
    for _, _, cls in res:
        for ax, c in cls.items():
            ga = glob_axis[ax]; ga["n_sig"] += c["n_sig"]; ga["n_hyper"] += c["n_hyper"]; ga["n_hypo"] += c["n_hypo"]
            ga["n_loci"] += 1 if c["n_sig"] > 0 else 0
            for i in range(4): ga["bins"][i] += c["bins"][i]
    out = {"axes": {ax: {**ga, "bin_labels": BINLAB, "hyper_pct": round(100 * ga["n_hyper"] / ga["n_sig"], 1) if ga["n_sig"] else 0} for ax, ga in glob_axis.items()},
           "n_loci": len(per_locus),
           "total_sig_cpg": sum(ga["n_sig"] for ga in glob_axis.values())}
    json.dump(out, open(f"{A}/percpg_cpg_classification.json", "w"), ensure_ascii=False, indent=1)
    print("DONE total_sig_cpg=" + str(out["total_sig_cpg"]) + " axes=" + str(list(glob_axis.keys())), flush=True)


if __name__ == "__main__":
    main()
