#!/usr/bin/env python3
"""
structure_vs_readcount_cn.py
條件分析：在「甲基判定有結構」的區 (structure_exists==True) 內，
分群狀況 (optimal_k / n_clusters / structure strength) 是否與 read 數 / SEQC2 CN / SAVANA CN 相關。
+ 結構『存在與否』是否被 coverage/CN 驅動。
輸入：region_table_tp.tsv (analyze_kism_vs_cn_perread.py 產) + SEQC2 gain/loss CN beds。
"""
import pandas as pd, numpy as np, json
import scipy.stats as st
from collections import defaultdict
from pathlib import Path

RES = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/20260621_kism_vs_cn_perread/_assets/region_table_tp.tsv"
D = "/big8_disk/data/HCC1395/SEQC2/CNV"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/20260621_kism_vs_cn_perread/_assets/structure_vs_readcount_cn.json"

res = pd.read_csv(RES, sep="\t")
seg = defaultdict(list)
for fn in ["ngs_benchmark_cnv_gain_cn.bed", "ngs_benchmark_cnv_loss_cn.bed"]:
    for ln in open(f"{D}/{fn}"):
        f = ln.split("\t")
        if len(f) >= 4:
            try: seg[f[0]].append((int(f[1]), int(f[2]), float(f[3])))
            except Exception: pass
for c in seg: seg[c].sort()
def seqc2_cn(chrom, pos):
    for s, e, cn in seg.get(chrom, []):
        if s <= pos <= e: return cn
    return 2.0
res["seqc2_cn"] = [seqc2_cn(c, int(p)) for c, p in zip(res["chr"], res["pos"])]

def sp(a, b):
    m = pd.DataFrame({"a": a, "b": b}).dropna()
    if len(m) < 10: return None
    rho, p = st.spearmanr(m["a"], m["b"])
    return {"rho": round(float(rho), 4), "p": float(p), "n": int(len(m))}

s = res[res["structure_exists"] == True].copy()
s["neglogp"] = -np.log10(s["permanova_p"].clip(lower=1e-6))
out = {"structured_n": int(len(s)), "within_structured": {
    "optimal_k_vs_nreads": sp(s["optimal_k"], s["n_reads"]),
    "optimal_k_vs_seqc2_cn": sp(s["optimal_k"], s["seqc2_cn"]),
    "optimal_k_vs_savana_cn": sp(s["optimal_k"], s["copyNumber"]),
    "n_clusters_vs_nreads": sp(s["n_clusters"], s["n_reads"]),
    "n_clusters_vs_seqc2_cn": sp(s["n_clusters"], s["seqc2_cn"]),
    "strength_vs_nreads": sp(s["neglogp"], s["n_reads"]),
    "strength_vs_seqc2_cn": sp(s["neglogp"], s["seqc2_cn"]),
}}
stru = res[res["structure_exists"]]; non = res[~res["structure_exists"]]
out["existence"] = {
    "nreads_median_structured": float(stru["n_reads"].median()),
    "nreads_median_unstructured": float(non["n_reads"].median()),
    "nreads_mwu_p": float(st.mannwhitneyu(stru["n_reads"], non["n_reads"])[1]),
    "seqc2_cn_median_structured": float(stru["seqc2_cn"].median()),
    "seqc2_cn_median_unstructured": float(non["seqc2_cn"].median()),
    "structure_rate_by_seqc2_cn": {str(cn): round(res[res["seqc2_cn"] == cn]["structure_exists"].mean()*100, 1)
                                   for cn in sorted(res["seqc2_cn"].unique()) if len(res[res["seqc2_cn"] == cn]) >= 50},
}
Path(OUT).write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
