#!/usr/bin/env python3
"""[逐 CpG 差異甲基核心 v2] site-first: 每 CpG 對一個二元分組測差異。
兩統計並列: Fisher 2×2(meth>0.5/unmeth × 組) + Mann-Whitney U(連續 β 秩和, rank-based 對 read 非獨立穩健)。
  (beta-binomial-DSS 需 replicate 估 dispersion → 單樣本 read-level 不適用, 改 MWU; pilot 已證壞掉的 quasi-φ 循環。)
BH-FDR 各自校正。sig = FDR<0.05 且 |Δβ|≥0.2(用 MWU)。
🔑 判別器(cis-ASM vs subclone): cluster-split sig CpG 集合 與 各 label sig CpG 集合的 containment overlap →
   cluster CpG ⊆ label = cis-ASM 複述; cluster CpG distinct 且 coherent run≥3 = subclone novel。
分組: HP1/HP2, HP1/HP1-1, HP2/HP2-1, ALT/REF, tumor/normal, cluster-split(主分叉 2 群)。純讀 region dir, 零 binary。"""
import os, csv
from collections import Counter
import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu

DTHR = 0.2; FDR = 0.05; RUNMIN = 3


def load_meth(region):
    mr = open(f"{region}/methylation/methylation.csv").read().strip().split("\n")
    if len(mr) < 2:
        return {}, []
    cpgs = [int(c) for c in mr[0].split(",")[1:]]
    M = {}
    for ln in mr[1:]:
        q = ln.split(",")
        M[q[0]] = np.array([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    return M, cpgs


def reads_meta(region):
    meta = {}
    rt = f"{region}/reads/reads.tsv"
    if os.path.exists(rt):
        rows = open(rt).read().splitlines(); hdr = rows[0].split("\t")
        ti = hdr.index("is_tumor") if "is_tumor" in hdr else None
        hi = hdr.index("hp") if "hp" in hdr else None
        for r in rows[1:]:
            c = r.split("\t"); rid = c[0]
            meta[rid] = {"hp": c[hi] if hi is not None and len(c) > hi else "",
                         "is_tumor": c[ti] if ti is not None and len(c) > ti else ""}
    grp = f"{region}/clustering/phylo_groups.tsv"
    if os.path.exists(grp):
        for row in csv.DictReader(open(grp), delimiter="\t"):
            rid = row["read_id"]; meta.setdefault(rid, {})
            meta[rid]["coarse"] = row.get("coarse_label", ""); meta[rid]["alt"] = row.get("alt_support", "")
    return meta


def groupings(meta):
    out = {}
    def binarize(name, fn):
        g = {rid: v for rid, m in meta.items() if (v := fn(m)) is not None}
        if g and sum(v == 0 for v in g.values()) >= 3 and sum(v == 1 for v in g.values()) >= 3:
            out[name] = g
    binarize("HP1_vs_HP2", lambda m: 0 if m.get("hp", "").startswith("1") else (1 if m.get("hp", "").startswith("2") else None))
    binarize("HP1_vs_HP1-1", lambda m: 0 if m.get("hp") == "1" else (1 if m.get("hp") == "1-1" else None))
    binarize("HP2_vs_HP2-1", lambda m: 0 if m.get("hp") == "2" else (1 if m.get("hp") == "2-1" else None))
    binarize("REF_vs_ALT", lambda m: 0 if m.get("alt") == "REF" else (1 if m.get("alt") == "ALT" else None))
    binarize("tumor_vs_normal", lambda m: 1 if m.get("is_tumor") == "1" else (0 if m.get("is_tumor") == "0" else None))
    # cluster-split = 主分叉(第一個有 ≥2 distinct prefix 的層級, 取 2 大 prefix)
    leaves = sorted(set(m["coarse"] for m in meta.values() if m.get("coarse") and m["coarse"] not in ("outlier", "other")))
    if len(leaves) >= 2:
        Ln = 2
        for L in range(2, 8):
            if len(set("-".join(x.split("-")[:L]) for x in leaves)) >= 2:
                Ln = L; break
        def pref(lab): return "-".join(lab.split("-")[:Ln]) if lab and lab not in ("outlier", "other") else None
        top2 = [k for k, _ in Counter(pref(m["coarse"]) for m in meta.values() if pref(m.get("coarse", ""))).most_common(2)]
        if len(top2) >= 2:
            binarize("cluster_split", lambda m: 0 if pref(m.get("coarse", "")) == top2[0] else (1 if pref(m.get("coarse", "")) == top2[1] else None))
    return out


def _bh(pv):
    pv = np.asarray(pv, float); n = len(pv); o = np.argsort(pv); q = np.empty(n); c = 1.0
    for i in range(n - 1, -1, -1):
        c = min(c, pv[o[i]] * n / (i + 1)); q[o[i]] = c
    return q


def test_grouping(M, cpgs, g):
    g0 = [r for r in M if g.get(r) == 0]; g1 = [r for r in M if g.get(r) == 1]
    K = len(cpgs); d = np.full(K, np.nan); pf = np.ones(K); pm = np.ones(K)
    A0 = np.array([M[r] for r in g0]); A1 = np.array([M[r] for r in g1])
    for j in range(K):
        v0 = A0[:, j][~np.isnan(A0[:, j])]; v1 = A1[:, j][~np.isnan(A1[:, j])]
        if len(v0) < 2 or len(v1) < 2:
            continue
        d[j] = np.mean(v1) - np.mean(v0)
        a, b = int((v0 > 0.5).sum()), int((v0 <= 0.5).sum()); c, e = int((v1 > 0.5).sum()), int((v1 <= 0.5).sum())
        try: pf[j] = fisher_exact([[a, b], [c, e]])[1]
        except Exception: pf[j] = 1.0
        try: pm[j] = mannwhitneyu(v0, v1, alternative="two-sided").pvalue
        except Exception: pm[j] = 1.0
    qf = _bh(pf); qm = _bh(pm)
    sigf = (qf < FDR) & (np.abs(d) >= DTHR); sigm = (qm < FDR) & (np.abs(d) >= DTHR)
    order = np.argsort(cpgs); best = cur = 0; sgn = 0
    for j in order:
        if sigm[j] and np.isfinite(d[j]) and (sgn == 0 or np.sign(d[j]) == sgn):
            cur += 1; sgn = np.sign(d[j]); best = max(best, cur)
        else:
            cur = 1 if (sigm[j] and np.isfinite(d[j])) else 0; sgn = np.sign(d[j]) if cur else 0
    return {"n_cpg": int(np.isfinite(d).sum()), "n_sig_fisher": int(sigf.sum()), "n_sig_mwu": int(sigm.sum()),
            "max_coherent_run": int(best), "sigset": set(np.where(sigm)[0].tolist()), "n0": len(g0), "n1": len(g1),
            "sig_dbetas": [round(float(d[j]), 3) for j in np.where(sigm)[0]],  # 顯著 CpG 的 Δβ(帶號) → 方向/強度
            "max_abs_dbeta": round(float(np.nanmax(np.abs(d))), 3) if np.isfinite(d).any() else None}


def analyze_region(region):
    M, cpgs = load_meth(region)
    if not M:
        return None
    res = {name: test_grouping(M, cpgs, g) for name, g in groupings(reads_meta(region)).items()}
    cs = res.get("cluster_split"); label_axes = {k: v for k, v in res.items() if k != "cluster_split"}
    overlap = {}; verdict = "no_cluster"
    if cs and cs["sigset"]:
        SC = cs["sigset"]
        for k, v in label_axes.items():
            overlap[k] = {"containment": round(len(SC & v["sigset"]) / len(SC), 3), "n_overlap": len(SC & v["sigset"]), "label_nsig": v["n_sig_mwu"]}
        best_cont = max((o["containment"] for o in overlap.values()), default=0.0)
        best_axis = max(overlap.items(), key=lambda kv: kv[1]["containment"], default=(None, None))[0]
        if best_cont >= 0.5:
            verdict = f"cis-ASM(cluster≈{best_axis})"
        elif cs["max_coherent_run"] >= RUNMIN:
            verdict = "subclone_novel"
        else:
            verdict = "weak_scattered"
    elif cs:
        verdict = "cluster_no_sigCpG"
    blab = max(label_axes.items(), key=lambda kv: kv[1]["n_sig_mwu"], default=(None, None))
    return {"n_cpg": len(cpgs), "verdict": verdict,
            "cluster_split_run": cs["max_coherent_run"] if cs else None, "cluster_split_nsig": cs["n_sig_mwu"] if cs else None,
            "best_label_axis": blab[0], "best_label_nsig": blab[1]["n_sig_mwu"] if blab[1] else 0, "overlap": overlap,
            "groupings": {k: {kk: vv for kk, vv in v.items() if kk != "sigset"} for k, v in res.items()}}
