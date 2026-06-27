#!/usr/bin/env python3
"""[差異 CpG 距 SNV 距離分佈] 全 34,736 走 region dirs, 每軸的顯著差異 CpG 與 SNV 錨點(window 中心 start+5000)
的距離分佈(bins 0-500/500-1k/1-2k/2-5k bp), 並對照 background(全 CpG 距離) → 差異是否局部富集近 SNV(cis)。
逐 CpG 座標來自 methylation.csv 欄頭(C++ 已輸出)。輸出 percpg_distance.json。"""
import os, glob, json, sys
from collections import defaultdict
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__)); import per_cpg_diff as P

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
KEEP = f"{WT}/output/_phylo_wg_full"
EDGES = [0, 500, 1000, 2000, 5001]; BL = ["0-500", "500-1k", "1-2k", "2-5k"]


def dbin(dist):
    for i in range(4):
        if EDGES[i] <= dist < EDGES[i + 1]:
            return i
    return 3


def work(region):
    try:
        M, cpgs = P.load_meth(region)
        if not M:
            return None
        base = os.path.basename(region).split("_"); start = int(base[1]); snv = start + 5000
        bg = [0, 0, 0, 0]
        for c in cpgs:
            bg[dbin(abs(c - snv))] += 1
        gs = P.groupings(P.reads_meta(region))
        per_axis = {}
        for name, g in gs.items():
            tg = P.test_grouping(M, cpgs, g)
            b = [0, 0, 0, 0]
            for j in tg["sigset"]:
                b[dbin(abs(cpgs[j] - snv))] += 1
            per_axis[name] = b
        return (bg, per_axis)
    except Exception:
        return None


def main():
    regions = [os.path.dirname(os.path.dirname(p)) for p in glob.glob(f"{KEEP}/**/clustering/phylo_groups.tsv", recursive=True)]
    print(f"regions: {len(regions)}", flush=True)
    with Pool(24) as p:
        res = [r for r in p.map(work, regions, chunksize=20) if r]
    bg = [0, 0, 0, 0]; axis = defaultdict(lambda: [0, 0, 0, 0])
    for b, pa in res:
        for i in range(4):
            bg[i] += b[i]
        for ax, ab in pa.items():
            for i in range(4):
                axis[ax][i] += ab[i]
    def pct(b):
        t = sum(b) or 1; return [round(100 * x / t, 1) for x in b]
    out = {"bins": BL, "background_all_cpg": {"counts": bg, "pct": pct(bg)},
           "axes": {ax: {"counts": ab, "pct": pct(ab),
                         "near_snv_enrich": round((ab[0] / max(1, sum(ab))) / (bg[0] / max(1, sum(bg))), 2)} for ax, ab in axis.items()},
           "n_loci": len(res)}
    json.dump(out, open(f"{A}/percpg_distance.json", "w"), ensure_ascii=False, indent=1)
    print("DONE " + json.dumps({"background_pct": out["background_all_cpg"]["pct"],
                                "axes_near_enrich": {a: v["near_snv_enrich"] for a, v in out["axes"].items()}}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
