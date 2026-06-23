#!/usr/bin/env python3
"""[每位點甲基(β)摘要抽取] 走保留 region dirs(output/_phylo_wg_full), 每位點從 methylation/methylation.csv
(read×CpG β 矩陣) + phylo_groups.tsv(coarse_label) + reads.tsv(is_tumor/hp) 算甲基觀察指標。
不重跑 binary。輸出 phylo_cpp_wg_full_methylation.json(key=(chrom,pos))。模仿 extract_tn_geometry.py。

每位點指標:
  n_cpg / mean_beta / std_beta / frac_hypo(β<0.3) / frac_hyper(β>0.7)
  dbeta_group = coarse leaf 群間 mean β 之 max-min (驅動切群的甲基差 = cis-ASM/subclone 訊號; <2 群=None)
  dbeta_tn    = mean β(tumor) - mean β(normal)  (somatic 甲基偏移)
  dbeta_hp    = mean β(HP1) - mean β(HP2)        (allelic/cis-ASM; hp 開頭 1 vs 2)
"""
import os, csv, glob, json, sys
import numpy as np
from multiprocessing import Pool

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
KEEP = f"{WT}/output/_phylo_wg_full"


def _load_meth(path):
    """methylation.csv: row0=header(read_id,cpg1,...); 回 {read_id(str): np.array(β, NaN for missing)}, n_cpg."""
    rows = open(path).read().splitlines()
    if len(rows) < 2:
        return {}, 0
    ncpg = len(rows[0].split(",")) - 1
    M = {}
    for ln in rows[1:]:
        f = ln.split(",")
        if len(f) < 2:
            continue
        vals = np.array([float(x) if x not in ("", "NA", "nan", "NaN") else np.nan for x in f[1:]])
        M[f[0]] = vals
    return M, ncpg


def _mean(arrs):
    """flatten list of β arrays, nanmean (None if all-NaN/empty)."""
    if not arrs:
        return None
    v = np.concatenate(arrs)
    if v.size == 0 or np.all(np.isnan(v)):
        return None
    return float(np.nanmean(v))


def work(region):
    try:
        base = os.path.basename(region); parts = base.split("_")
        chrom, start = parts[0], int(parts[1])
        meth = f"{region}/methylation/methylation.csv"
        grp = f"{region}/clustering/phylo_groups.tsv"
        rt = f"{region}/reads/reads.tsv"
        if not os.path.exists(meth):
            return None
        M, ncpg = _load_meth(meth)
        if not M:
            return None
        allv = np.concatenate([M[r] for r in M])
        mean_beta = _mean([allv])
        if mean_beta is None:
            return None
        std_beta = float(np.nanstd(allv))
        frac_hypo = float(np.nanmean(allv < 0.3))
        frac_hyper = float(np.nanmean(allv > 0.7))
        # dbeta_group: split by coarse_label (peeled subset)
        dbeta_group = None
        if os.path.exists(grp):
            g2 = {}
            for row in csv.DictReader(open(grp), delimiter="\t"):
                if row.get("is_outlier") == "1":
                    continue
                g2.setdefault(row["coarse_label"], []).append(row["read_id"])
            gm = []
            for g, reads in g2.items():
                mb = _mean([M[r] for r in reads if r in M])
                if mb is not None:
                    gm.append(mb)
            if len(gm) >= 2:
                dbeta_group = float(max(gm) - min(gm))
        # dbeta_tn / dbeta_hp from reads.tsv
        dbeta_tn = dbeta_hp = None
        if os.path.exists(rt):
            rows = open(rt).read().splitlines()
            if rows:
                hdr = rows[0].split("\t")
                ti = hdr.index("is_tumor") if "is_tumor" in hdr else None
                hi = hdr.index("hp") if "hp" in hdr else None
                tum, nor, hp1, hp2 = [], [], [], []
                for r in rows[1:]:
                    c = r.split("\t")
                    rid = c[0]
                    if rid not in M:
                        continue
                    if ti is not None and len(c) > ti:
                        (tum if c[ti] == "1" else nor if c[ti] == "0" else []).append(M[rid])
                    if hi is not None and len(c) > hi:
                        h = c[hi]
                        if h.startswith("1"):
                            hp1.append(M[rid])
                        elif h.startswith("2"):
                            hp2.append(M[rid])
                mt, mn = _mean(tum), _mean(nor)
                if mt is not None and mn is not None:
                    dbeta_tn = float(mt - mn)
                m1, m2 = _mean(hp1), _mean(hp2)
                if m1 is not None and m2 is not None:
                    dbeta_hp = float(m1 - m2)
        return (chrom, start, {
            "n_cpg": ncpg, "mean_beta": round(mean_beta, 4), "std_beta": round(std_beta, 4),
            "frac_hypo": round(frac_hypo, 4), "frac_hyper": round(frac_hyper, 4),
            "dbeta_group": round(dbeta_group, 4) if dbeta_group is not None else None,
            "dbeta_tn": round(dbeta_tn, 4) if dbeta_tn is not None else None,
            "dbeta_hp": round(dbeta_hp, 4) if dbeta_hp is not None else None,
        })
    except Exception:
        return None


def main():
    regions = [os.path.dirname(os.path.dirname(p)) for p in
               glob.glob(f"{KEEP}/**/clustering/phylo_groups.tsv", recursive=True)]
    print(f"regions: {len(regions)}", flush=True)
    with Pool(24) as p:
        res = [r for r in p.map(work, regions, chunksize=40) if r]
    out = {f"{c}:{s}": d for (c, s, d) in res}
    json.dump(out, open(f"{A}/phylo_cpp_wg_full_methylation.json", "w"))
    # coverage / sanity stats
    nz_grp = sum(1 for d in out.values() if d["dbeta_group"] is not None)
    nz_tn = sum(1 for d in out.values() if d["dbeta_tn"] is not None)
    mbs = [d["mean_beta"] for d in out.values()]
    S = {"matched": len(out), "regions_walked": len(regions),
         "with_dbeta_group(>=2grp)": nz_grp, "with_dbeta_tn": nz_tn,
         "overall_mean_beta_median": round(float(np.median(mbs)), 4) if mbs else None,
         "mean_beta_min": round(min(mbs), 4) if mbs else None,
         "mean_beta_max": round(max(mbs), 4) if mbs else None}
    json.dump(S, open(f"{A}/methylation_stats.json", "w"), indent=2)
    print("DONE " + json.dumps(S, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
