#!/usr/bin/env python3
"""[chr17:48360161 綜合觀察 compute] BERNOULLI 距離(精確,DistanceMatrix.cpp 公式) + UPGMA + cut 分類
+ phylo coarse_label join(ISM 甲基分類) + 甲基cluster×基因型lineage 交叉 + per-CpG 基因軸歸因
(每 CpG 對齊哪個 sSNV REF/ALT / lineage L1-L2 / HP1-HP1-1)。輸出 chr17_tree_data.json。"""
import json, sys, csv
from collections import Counter, defaultdict
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import mannwhitneyu

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
D = json.load(open(f"{A}/chr17_subclone_data.json"))
reads = D["reads"]; cpgs = D["cpgs"]; K = len(cpgs); n = len(reads)
M = np.array([[np.nan if v is None else v for v in r["meth"]] for r in reads])  # n×K
LIN = {"L0_ancestral_root": "L0", "L1_alpha_only(ancestor)": "L1", "L2_alpha_beta(descendant)": "L2", "other": "other"}

# ---- region phylo coarse_label join ----
region = [c["region_dir"] for c in json.load(open(f"{A}/cis_candidates_resolved.json"))
          if c["chrom"] == "chr17" and c["pos"] == "48360161"][0]
coarse = {}
for row in csv.DictReader(open(f"{region}/clustering/phylo_groups.tsv"), delimiter="\t"):
    coarse[row["read_id"]] = row.get("coarse_label", "")
for r in reads:
    r["coarse"] = coarse.get(r["rid"], "NA")


# ---- BERNOULLI 距離 (精確公式) ----
def bern(pi, pj):
    mask = ~np.isnan(pi) & ~np.isnan(pj)
    if mask.sum() < 3:
        return None
    p = pi[mask]; q = pj[mask]
    w = (2 * np.abs(p - 0.5)) * (2 * np.abs(q - 0.5))
    delta = p * (1 - q) + (1 - p) * q
    sw = w.sum()
    return float((w * delta).sum() / sw) if sw > 1e-9 else 0.0


Dm = np.zeros((n, n)); miss = 0
for i in range(n):
    for j in range(i + 1, n):
        d = bern(M[i], M[j])
        if d is None:
            d = 0.5; miss += 1
        Dm[i, j] = Dm[j, i] = d
Z = linkage(squareform(Dm, checks=False), method="average")
order = [int(x) for x in leaves_list(Z)]
clu2 = [int(x) for x in fcluster(Z, 2, "maxclust")]
clu3 = [int(x) for x in fcluster(Z, 3, "maxclust")]

# ---- 交叉表: 甲基分類 × 基因型 lineage ----
def crosstab(labels):
    ct = defaultdict(Counter)
    for r, lab in zip(reads, labels):
        ct[str(lab)][LIN[r["lineage"]]] += 1
    return {k: dict(v) for k, v in ct.items()}


cross_coarse = crosstab([r["coarse"] for r in reads])
cross_clu2 = crosstab(clu2)
cross_clu3 = crosstab(clu3)


# ---- per-CpG 基因軸歸因 ----
def axis_label(r, ax):
    g = r["geno"]
    if ax == "ALT@48365089(α)":
        return {"REF": 0, "ALT": 1}.get(g.get("48365089"))
    if ax == "ALT@48365161(β2)":
        return {"REF": 0, "ALT": 1}.get(g.get("48365161"))
    if ax == "ALT@48362515(β1)":
        return {"REF": 0, "ALT": 1}.get(g.get("48362515"))
    if ax == "lineage_L1_vs_L2":
        return {"L1_alpha_only(ancestor)": 0, "L2_alpha_beta(descendant)": 1}.get(r["lineage"])
    if ax == "HP1_vs_HP1-1":
        return {"1": 0, "1-1": 1}.get(r["hp"])
    return None


AXES = ["ALT@48365089(α)", "ALT@48365161(β2)", "ALT@48362515(β1)", "lineage_L1_vs_L2", "HP1_vs_HP1-1"]


def bh(pv):
    pv = np.asarray(pv, float); m = len(pv)
    if m == 0:
        return pv
    o = np.argsort(pv); q = np.empty(m); c = 1.0
    for i in range(m - 1, -1, -1):
        c = min(c, pv[o[i]] * m / (i + 1)); q[o[i]] = c
    return q


percpg = []
axis_sig = {ax: 0 for ax in AXES}
for j in range(K):
    rowj = {"cpg": cpgs[j], "axes": {}}
    for ax in AXES:
        g0 = [M[i, j] for i in range(n) if axis_label(reads[i], ax) == 0 and not np.isnan(M[i, j])]
        g1 = [M[i, j] for i in range(n) if axis_label(reads[i], ax) == 1 and not np.isnan(M[i, j])]
        if len(g0) >= 3 and len(g1) >= 3:
            db = float(np.mean(g1) - np.mean(g0))
            try:
                p = float(mannwhitneyu(g0, g1, alternative="two-sided").pvalue)
            except ValueError:
                p = 1.0
            rowj["axes"][ax] = {"dbeta": round(db, 3), "p": round(p, 4), "n0": len(g0), "n1": len(g1)}
    percpg.append(rowj)
# BH per axis + best axis per CpG
for ax in AXES:
    idx = [k for k in range(K) if ax in percpg[k]["axes"]]
    qs = bh([percpg[k]["axes"][ax]["p"] for k in idx])
    for t, k in enumerate(idx):
        percpg[k]["axes"][ax]["q"] = round(float(qs[t]), 4)
        sig = qs[t] < 0.05 and abs(percpg[k]["axes"][ax]["dbeta"]) >= 0.2
        percpg[k]["axes"][ax]["sig"] = bool(sig)
        if sig:
            axis_sig[ax] += 1
for k in range(K):
    sigax = [(ax, percpg[k]["axes"][ax]) for ax in AXES if percpg[k]["axes"].get(ax, {}).get("sig")]
    percpg[k]["best_axis"] = max(sigax, key=lambda x: abs(x[1]["dbeta"]))[0] if sigax else None

best_counts = Counter(p["best_axis"] for p in percpg if p["best_axis"])
out = {"n_reads": n, "n_cpg": K, "dist_miss_pairs": miss,
       "read_order_upgma": order, "linkage_heights": [round(float(Z[i, 2]), 4) for i in range(len(Z))],
       "clusters_k2": clu2, "clusters_k3": clu3,
       "lineage_per_read": [LIN[r["lineage"]] for r in reads], "coarse_per_read": [r["coarse"] for r in reads],
       "cross_coarse_x_lineage": cross_coarse, "cross_clu2_x_lineage": cross_clu2, "cross_clu3_x_lineage": cross_clu3,
       "distance_matrix": [[round(float(Dm[i, j]), 3) for j in range(n)] for i in range(n)],
       "axis_sig_cpg_count": axis_sig, "best_axis_counts": dict(best_counts), "percpg_attrib": percpg}
json.dump(out, open(f"{A}/chr17_tree_data.json", "w"), ensure_ascii=False)
print(json.dumps({"甲基cluster(coarse)×lineage": cross_coarse, "UPGMA_cut2×lineage": cross_clu2,
                  "每軸顯著CpG數": axis_sig, "每CpG最佳歸因軸": dict(best_counts)}, ensure_ascii=False, indent=1))
print("[-> chr17_tree_data.json]")
