#!/usr/bin/env python3
"""[每位點 tumor/normal + geometry 切群數 補充] 走保留的 region dirs(output/_phylo_wg_full),
每位點抽: n_tumor/n_normal(reads.tsv is_tumor) + geometry_ng(matrix 0.7×max 平切, size≥3)。
不重跑 binary(用保留 matrix/reads)。輸出 phylo_cpp_wg_full_records_full.json(records_cn + n_tumor/n_normal/geometry_ng)。"""
import os, csv, glob, json, sys
import numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
from scipy.cluster.hierarchy import fcluster
from multiprocessing import Pool

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
KEEP = f"{WT}/output/_phylo_wg_full"


def work(region):
    try:
        base = os.path.basename(region); parts = base.split("_")
        chrom, start = parts[0], int(parts[1])
        # T/N from reads.tsv
        nt = nn = 0
        rt = f"{region}/reads/reads.tsv"
        if os.path.exists(rt):
            rows = open(rt).read().splitlines()
            if rows:
                hdr = rows[0].split("\t")
                ti = hdr.index("is_tumor") if "is_tumor" in hdr else None
                if ti is not None:
                    for r in rows[1:]:
                        v = r.split("\t")[ti] if len(r.split("\t")) > ti else ""
                        if v == "1": nt += 1
                        elif v == "0": nn += 1
        # geometry_ng from matrix + phylo_groups.tsv (peeled subset)
        geom = None
        tsv = f"{region}/clustering/phylo_groups.tsv"; mat = f"{region}/distance/BERNOULLI/matrix.csv"
        if os.path.exists(tsv) and os.path.exists(mat):
            rid = [int(r["read_id"]) for r in csv.DictReader(open(tsv), delimiter="\t")]
            dids, Dm = CR.loadm(mat)
            idx = [dids.index(str(x)) for x in rid if str(x) in dids]
            if len(idx) >= 6:
                sub = Dm[np.ix_(idx, idx)]
                Z, _ = CR.linkZ(sub)
                lab = fcluster(Z, 0.7 * float(Z[:, 2].max()), "distance")
                geom = sum(1 for s in Counter(lab).values() if s >= 3)
        return (chrom, start, nt, nn, geom)
    except Exception:
        return None


def main():
    regions = [os.path.dirname(os.path.dirname(p)) for p in glob.glob(f"{KEEP}/**/clustering/phylo_groups.tsv", recursive=True)]
    print(f"regions: {len(regions)}", flush=True)
    with Pool(24) as p:
        res = [r for r in p.map(work, regions, chunksize=40) if r]
    amap = {(c, s): (nt, nn, g) for (c, s, nt, nn, g) in res}
    rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_cn.json"))
    hit = 0
    for r in rec:
        kk = (r["chrom"], int(r["pos"]))
        if kk in amap:
            nt, nn, g = amap[kk]
            r["n_tumor"] = nt; r["n_normal"] = nn; r["geometry_ng"] = g; hit += 1
        else:
            r["n_tumor"] = None; r["n_normal"] = None; r["geometry_ng"] = None
    json.dump(rec, open(f"{A}/phylo_cpp_wg_full_records_full.json", "w"))
    # stats
    g_div = sum(1 for r in rec if r.get("geometry_ng") and r["geometry_ng"] >= 2 and r["coarse_ng"] < 2)
    tn = [(r.get("n_tumor"), r.get("n_normal")) for r in rec if r.get("n_tumor") is not None]
    S = {"matched": hit, "n": len(rec),
         "geometry_divergence(geom≥2 but coarse<2)": g_div,
         "median_n_tumor": int(np.median([t for t, _ in tn])) if tn else None,
         "median_n_normal": int(np.median([n for _, n in tn])) if tn else None,
         "geom_ng_dist": dict(Counter(r.get("geometry_ng") for r in rec if r.get("geometry_ng")))}
    json.dump(S, open(f"{A}/tn_geometry_stats.json", "w"), indent=2)
    print("DONE " + json.dumps(S, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
