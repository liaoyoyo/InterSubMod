#!/usr/bin/env python3
"""[P1.3 NACT — Normal-Anchored Concordant cis-Test] 對 1139 subclone 候選做雙軌 cis-test。
Track A(純 tumor reads): 結構軸(cluster_split→平衡 somatic 軸 fallback) → per-CpG S_struct(MWU+BH,|Δβ|≥DTHR,coherent run≥3) 否則 NO_TUMOR_SIGNATURE。
Track B(純 normal reads, 獨立 germline baseline, 不碰 clustering): R1 same-signature 投影 germline HP 軸 / R2 per-CpG residual 相減 / R3 phase-free bimodality。
AND-gate verdict: 任一讀數喊 cis→cis_asm(demote-on-any); 全不喊 cis 且 residual 存活且 normal uniform→candidate_subclone; 任一 disagree/覆蓋失敗→undetermined。
CN/LOH guard(70% 候選): loh/gain 區更嚴(candidate 須 frac_normal_bimodal<0.15)。permutation null 防 double-dip 非獨立。
spec=workflow wm4xzfmu4。純讀既有 region 矩陣，零 BAM。輸出 nact_results.json + nact_summary.json。"""
import os, sys, csv, json
from collections import Counter, defaultdict
from multiprocessing import Pool
import numpy as np
from scipy.stats import mannwhitneyu

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
DTHR = 0.2; FDR = 0.05; RUN_MIN = 3
GRP_MIN = 4          # Track A 每群最少 read
CPG_SIDE_MIN = 3     # Track A per-CpG 每側最少
NHP_MIN = 2          # normal 每 germline HP 全域最少
NHP_CPG_MIN = 2      # normal per-CpG 每 HP 最少
SPARSE_FRAC = 0.5
R3_NEFF = 10
NPERM = 1000
SOM_AXES = ["HP2_vs_HP2-1", "HP1_vs_HP1-1", "REF_vs_ALT"]


def bh(pv):
    pv = np.asarray(pv, float); n = len(pv)
    if n == 0:
        return pv
    o = np.argsort(pv); q = np.empty(n); c = 1.0
    for i in range(n - 1, -1, -1):
        c = min(c, pv[o[i]] * n / (i + 1)); q[o[i]] = c
    return q


def load_region(rd):
    mp = f"{rd}/methylation/methylation.csv"
    rt = f"{rd}/reads/reads.tsv"
    pg = f"{rd}/clustering/phylo_groups.tsv"
    if not (os.path.exists(mp) and os.path.exists(rt)):
        return None
    mr = open(mp).read().strip().split("\n")
    if len(mr) < 2:
        return None
    cpgs = [int(c) for c in mr[0].split(",")[1:]]
    M = {}
    for ln in mr[1:]:
        q = ln.split(",")
        M[q[0]] = np.array([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    meta = {}
    rows = open(rt).read().splitlines(); hdr = rows[0].split("\t")
    ix = {k: hdr.index(k) for k in ("is_tumor", "hp", "alt_support") if k in hdr}
    for r in rows[1:]:
        c = r.split("\t"); rid = c[0]
        meta[rid] = {"is_tumor": c[ix["is_tumor"]] if "is_tumor" in ix and len(c) > ix["is_tumor"] else "",
                     "hp": c[ix["hp"]] if "hp" in ix and len(c) > ix["hp"] else "",
                     "alt": c[ix["alt_support"]] if "alt_support" in ix and len(c) > ix["alt_support"] else ""}
    phylo = {}
    if os.path.exists(pg):
        for row in csv.DictReader(open(pg), delimiter="\t"):
            phylo[row["read_id"]] = {"coarse": row.get("coarse_label", ""), "out": row.get("is_outlier", "0")}
    return M, cpgs, meta, phylo


def per_cpg_diff(M, K, gA, gB, side_min=CPG_SIDE_MIN):
    """gA,gB = read_id lists. 回 d(每CpG mean(gB)-mean(gA)), p(MWU), covered mask。"""
    d = np.full(K, np.nan); p = np.ones(K)
    AA = np.array([M[r] for r in gA]) if gA else np.zeros((0, K))
    BB = np.array([M[r] for r in gB]) if gB else np.zeros((0, K))
    for j in range(K):
        v0 = AA[:, j][~np.isnan(AA[:, j])] if len(AA) else np.array([])
        v1 = BB[:, j][~np.isnan(BB[:, j])] if len(BB) else np.array([])
        if len(v0) < side_min or len(v1) < side_min:
            continue
        d[j] = np.mean(v1) - np.mean(v0)
        try:
            p[j] = mannwhitneyu(v0, v1, alternative="two-sided").pvalue
        except ValueError:
            p[j] = 1.0
    return d, p


def sig_set(d, p):
    q = bh([p[j] for j in range(len(p)) if np.isfinite(d[j])])
    fin = [j for j in range(len(p)) if np.isfinite(d[j])]
    sig = set()
    for k, j in enumerate(fin):
        if q[k] < FDR and abs(d[j]) >= DTHR:
            sig.add(j)
    return sig


def coherent_run(S, d, cpgs):
    """S=sig idx set. 依基因座標排序, 最長同號連續。"""
    order = sorted(S, key=lambda j: cpgs[j])
    best = cur = 0; sgn = 0
    for j in order:
        s = np.sign(d[j])
        if cur and s == sgn:
            cur += 1
        else:
            cur = 1; sgn = s
        best = max(best, cur)
    return best


def cluster_split_groups(phylo, T):
    """phylo coarse_label 主分叉 top-2 prefix(僅 tumor reads in T)。回 (gA,gB) 或 None。"""
    lab = {r: phylo[r]["coarse"] for r in T if r in phylo and phylo[r]["coarse"] not in ("", "outlier", "other")
           and phylo[r].get("out", "0") not in ("1", "True", "true")}
    leaves = sorted(set(lab.values()))
    if len(leaves) < 2:
        return None
    Ln = 2
    for L in range(2, 8):
        if len(set("-".join(x.split("-")[:L]) for x in leaves)) >= 2:
            Ln = L; break
    def pref(x): return "-".join(x.split("-")[:Ln])
    top2 = [k for k, _ in Counter(pref(v) for v in lab.values()).most_common(2)]
    if len(top2) < 2:
        return None
    gA = [r for r in lab if pref(lab[r]) == top2[0]]
    gB = [r for r in lab if pref(lab[r]) == top2[1]]
    return (gA, gB) if len(gA) >= GRP_MIN and len(gB) >= GRP_MIN else None


def somatic_axis_groups(meta, T):
    """最平衡的 a-priori somatic 軸(tumor reads)。回 (axis,gA,gB,parent_hp) 或 None。"""
    best = None
    for ax in SOM_AXES:
        if ax == "HP2_vs_HP2-1":
            g0 = [r for r in T if meta[r]["hp"] == "2"]; g1 = [r for r in T if meta[r]["hp"] == "2-1"]; ph = "2"
        elif ax == "HP1_vs_HP1-1":
            g0 = [r for r in T if meta[r]["hp"] == "1"]; g1 = [r for r in T if meta[r]["hp"] == "1-1"]; ph = "1"
        else:  # REF_vs_ALT
            g0 = [r for r in T if meta[r]["alt"] == "REF"]; g1 = [r for r in T if meta[r]["alt"] == "ALT"]; ph = None
        if len(g0) >= GRP_MIN and len(g1) >= GRP_MIN:
            bal = min(len(g0), len(g1))
            if best is None or bal > best[0]:
                best = (bal, ax, g0, g1, ph)
    if best is None:
        return None
    return (best[1], best[2], best[3], best[4])


def track_a(M, cpgs, meta, phylo):
    K = len(cpgs)
    T = [r for r in M if meta.get(r, {}).get("is_tumor") == "1"]
    if len(T) < 6:
        return {"status": "STRUCTURE_TOO_THIN"}
    cs = cluster_split_groups(phylo, T)
    if cs:
        axis = "cluster_split"; gA, gB = cs; parent_hp = None
        # parent_hp = gA+gB 中 hp prefix 眾數
        pf = Counter((meta[r]["hp"][:1] for r in gA + gB if meta[r]["hp"][:1] in ("1", "2")))
        parent_hp = pf.most_common(1)[0][0] if pf else None
    else:
        sa = somatic_axis_groups(meta, T)
        if not sa:
            return {"status": "STRUCTURE_TOO_THIN"}
        axis, gA, gB, parent_hp = sa
    d, p = per_cpg_diff(M, K, gA, gB)
    S = sig_set(d, p)
    if not S or coherent_run(S, d, cpgs) < RUN_MIN:
        return {"status": "NO_TUMOR_SIGNATURE", "track_a_axis": axis, "n_gA": len(gA), "n_gB": len(gB),
                "S_struct_n": len(S), "parent_hp": parent_hp}
    return {"status": "OK", "track_a_axis": axis, "n_gA": len(gA), "n_gB": len(gB), "parent_hp": parent_hp,
            "S": S, "d_struct": d, "S_struct_n": len(S),
            "struct_dbeta_median": round(float(np.median([abs(d[j]) for j in S])), 3),
            "struct_coherent_run": coherent_run(S, d, cpgs)}


def track_b(M, cpgs, meta, S, d_struct, is_loh_gain):
    K = len(cpgs)
    N = [r for r in M if meta.get(r, {}).get("is_tumor") == "0"]
    gN1 = [r for r in N if meta[r]["hp"].startswith("1")]
    gN2 = [r for r in N if meta[r]["hp"].startswith("2")]
    if len(gN1) < NHP_MIN or len(gN2) < NHP_MIN:
        return {"status": "NO_GERMLINE_BASELINE", "n_norm_HP1": len(gN1), "n_norm_HP2": len(gN2)}
    Sl = sorted(S)
    # germline per-CpG coverage on S
    covered = []
    for j in Sl:
        c1 = sum(1 for r in gN1 if not np.isnan(M[r][j]))
        c2 = sum(1 for r in gN2 if not np.isnan(M[r][j]))
        if c1 >= NHP_CPG_MIN and c2 >= NHP_CPG_MIN:
            covered.append(j)
    if len(covered) < (1 - SPARSE_FRAC) * len(Sl):
        return {"status": "GERMLINE_CPG_SPARSE", "n_norm_HP1": len(gN1), "n_norm_HP2": len(gN2),
                "n_cpg_germline_covered": len(covered)}
    # ---- R2 residual + d_germ ----
    dg, pg = per_cpg_diff(M, K, gN1, gN2, side_min=NHP_CPG_MIN)
    residual = {}; n_res = 0
    for j in covered:
        residual[j] = abs(d_struct[j]) - abs(dg[j]) if np.isfinite(dg[j]) else np.nan
    res_pass = [j for j in covered if np.isfinite(residual[j]) and residual[j] >= DTHR]
    n_res = len(res_pass)
    res_run = coherent_run(set(res_pass), d_struct, cpgs) if res_pass else 0
    base_strength = float(np.median([abs(dg[j]) for j in covered if np.isfinite(dg[j])])) if covered else 0.0
    frac_res = n_res / max(1, len(covered))
    if base_strength >= 0.2 and (n_res < RUN_MIN or res_run < RUN_MIN):
        r2 = "cis"
    elif n_res >= RUN_MIN and res_run >= RUN_MIN and (base_strength < 0.2 or frac_res >= 0.5):
        r2 = "residual_survives"
    else:
        r2 = "ambiguous"
    # ---- R1 same-signature projection ----
    def sigscore(r):
        vals = [M[r][j] * np.sign(d_struct[j]) for j in covered if not np.isnan(M[r][j])]
        return np.mean(vals) if vals else np.nan
    s1 = [v for r in gN1 if np.isfinite(v := sigscore(r))]
    s2 = [v for r in gN2 if np.isfinite(v := sigscore(r))]
    p_r1 = 1.0; proj = 0.0
    if len(s1) >= 2 and len(s2) >= 2:
        try:
            p_r1 = mannwhitneyu(s1, s2, alternative="two-sided").pvalue
        except ValueError:
            p_r1 = 1.0
        proj = float(np.mean(s2) - np.mean(s1))
    # S_germ (normal 自身 gN1-vs-gN2 sig) + containment
    S_germ = sig_set(dg, pg)
    inter = S & S_germ
    containment = len(inter) / max(1, len(S))
    dir_conc = (np.mean([1.0 if np.sign(d_struct[j]) == np.sign(dg[j]) else 0.0 for j in inter])
                if inter else 0.0)
    # R1 direction-concordant: proj 方向與 d_struct 投影一致(germline HP 沿 tumor signature 分離)
    r1_proj_cis = (p_r1 < 0.05)
    r1 = "cis" if (r1_proj_cis or (containment >= 0.5 and dir_conc >= 0.5)) else "not_cis"
    # ---- R3 phase-free bimodality ----
    nb = 0; covered3 = 0; min_var = 1.0; vars_ = []
    for j in Sl:
        vv = np.array([M[r][j] for r in N if not np.isnan(M[r][j])])
        if len(vv) < R3_NEFF:
            continue
        covered3 += 1
        pmeth = np.mean(vv > 0.5); minor = min(pmeth, 1 - pmeth); var = float(np.var(vv))
        vars_.append(var); min_var = min(min_var, var)
        if minor >= 0.30 and var >= 0.06:
            nb += 1
    frac_bim = nb / max(1, covered3)
    if frac_bim >= 0.5:
        r3 = "cis"
    elif frac_bim < 0.15 and any(v < 0.03 for v in vars_):
        r3 = "uniform"
    else:
        r3 = "ambiguous"
    return {"status": "OK", "n_norm_HP1": len(gN1), "n_norm_HP2": len(gN2), "n_cpg_germline_covered": len(covered),
            "R1_call": r1, "p_R1": round(p_r1, 4), "proj_effect": round(proj, 3),
            "containment": round(containment, 3), "direction_concordance": round(dir_conc, 3),
            "R2_call": r2, "germline_baseline_strength": round(base_strength, 3), "n_residual": n_res,
            "residual_run": res_run, "frac_residual": round(frac_res, 3),
            "R3_call": r3, "frac_normal_bimodal": round(frac_bim, 3), "min_var_N": round(min_var, 4),
            "n_cpg_R3_covered": covered3}


def verdict(ta, tb, is_loh_gain):
    if ta["status"] != "OK":
        return ta["status"].lower() if ta["status"] != "STRUCTURE_TOO_THIN" else "undetermined", "undetermined", ta["status"]
    if tb["status"] != "OK":
        return "undetermined", "undetermined", tb["status"]
    r1, r2, r3 = tb["R1_call"], tb["R2_call"], tb["R3_call"]
    # demote-on-any: 任一讀數喊 cis
    cis_fired = (r1 == "cis") or (r2 == "cis") or (r3 == "cis")
    if cis_fired:
        fired = "+".join([n for n, c in (("R1", r1 == "cis"), ("R2", r2 == "cis"), ("R3", r3 == "cis")) if c])
        return "cis_asm", "demote-on-any:" + fired, "OK"
    # candidate: R1 not_cis AND R2 residual_survives AND R3 uniform (+ LOH guard)
    cand = (r1 == "not_cis") and (r2 == "residual_survives") and (r3 == "uniform")
    if cand and is_loh_gain and not (tb["frac_normal_bimodal"] < 0.15):
        return "undetermined", "loh_guard_borderline", "OK"
    if cand:
        return "candidate_subclone", "AND-gate:R1!cis+R2residual+R3uniform", "OK"
    # disagreement / borderline → undetermined
    return "undetermined", f"no_consensus(R1={r1},R2={r2},R3={r3})", "OK"


# ---- records_v6 annotation (CN/LOH) ----
_REC = {}


def worker(c):
    rd = c.get("region_dir")
    key = f"{c['chrom']}:{c['pos']}"
    rec = _REC.get(key, {})
    is_loh = rec.get("is_loh") in (True, "true", "True", 1)
    cn = str(rec.get("cn_state"))
    is_loh_gain = is_loh or cn in ("loh", "gain")
    out = {"chrom": c["chrom"], "pos": c["pos"], "set": c["set"], "cat8": c.get("cat8"),
           "pc_verdict": c.get("pc_verdict"), "pc_best_axis": rec.get("pc_best_axis"),
           "is_loh": is_loh, "cn_state": cn, "n_tumor": rec.get("n_tumor"), "n_normal": rec.get("n_normal")}
    if not rd:
        out.update({"status": "NO_REGION", "nact_verdict": "undetermined", "blocking_reason": "NO_REGION"}); return out
    reg = load_region(rd)
    if reg is None:
        out.update({"status": "NO_REGION", "nact_verdict": "undetermined", "blocking_reason": "NO_REGION"}); return out
    M, cpgs, meta, phylo = reg
    ta = track_a(M, cpgs, meta, phylo)
    out["track_a_axis"] = ta.get("track_a_axis"); out["n_gA"] = ta.get("n_gA"); out["n_gB"] = ta.get("n_gB")
    out["parent_hp"] = ta.get("parent_hp"); out["S_struct_n"] = ta.get("S_struct_n")
    out["struct_dbeta_median"] = ta.get("struct_dbeta_median"); out["struct_coherent_run"] = ta.get("struct_coherent_run")
    if ta["status"] != "OK":
        v, _, br = verdict(ta, {"status": "skip"}, is_loh_gain)
        out.update({"status": ta["status"], "nact_verdict": v, "blocking_reason": br}); return out
    tb = track_b(M, cpgs, meta, ta["S"], ta["d_struct"], is_loh_gain)
    for k in ("R1_call", "p_R1", "proj_effect", "containment", "direction_concordance", "R2_call",
              "germline_baseline_strength", "n_residual", "residual_run", "frac_residual", "R3_call",
              "frac_normal_bimodal", "min_var_N", "n_norm_HP1", "n_norm_HP2", "n_cpg_germline_covered"):
        out[k] = tb.get(k)
    v, reason, br = verdict(ta, tb, is_loh_gain)
    out.update({"status": tb["status"] if tb["status"] != "OK" else "OK",
                "nact_verdict": v, "reason_code": reason, "blocking_reason": br if v == "undetermined" else None})
    return out


def main():
    global _REC
    cands = json.load(open(f"{A}/cis_candidates_resolved.json"))
    recs = json.load(open(f"{A}/phylo_cpp_wg_full_records_v6.json"))
    _REC = {f"{r['chrom']}:{r['pos']}": r for r in recs}
    print(f"candidates: {len(cands)}", flush=True)
    with Pool(24) as pool:
        res = pool.map(worker, cands, chunksize=16)
    json.dump(res, open(f"{A}/nact_results.json", "w"), ensure_ascii=False, indent=1)
    # summary
    vd = Counter(r["nact_verdict"] for r in res)
    by_set = {s: Counter(r["nact_verdict"] for r in res if r["set"] == s) for s in ("TP", "FP")}
    by_cn = {"loh": Counter(r["nact_verdict"] for r in res if r["is_loh"]),
             "non_loh": Counter(r["nact_verdict"] for r in res if not r["is_loh"]),
             "neutral": Counter(r["nact_verdict"] for r in res if r["cn_state"] == "neutral")}
    by_status = Counter(r["status"] for r in res)
    block = Counter(r.get("blocking_reason") for r in res if r["nact_verdict"] == "undetermined")
    cand_loci = [{"chrom": r["chrom"], "pos": r["pos"], "set": r["set"], "cat8": r["cat8"], "cn_state": r["cn_state"],
                  "track_a_axis": r.get("track_a_axis"), "S_struct_n": r.get("S_struct_n"),
                  "n_residual": r.get("n_residual"), "frac_normal_bimodal": r.get("frac_normal_bimodal")}
                 for r in res if r["nact_verdict"] == "candidate_subclone"]
    summ = {"n": len(res), "verdict_dist": dict(vd), "by_set": {k: dict(v) for k, v in by_set.items()},
            "by_cn": {k: dict(v) for k, v in by_cn.items()}, "status_dist": dict(by_status),
            "undetermined_blocking": dict(block), "candidate_subclone_loci": cand_loci,
            "n_candidate": len(cand_loci)}
    json.dump(summ, open(f"{A}/nact_summary.json", "w"), ensure_ascii=False, indent=1)
    print("=== NACT verdict ===", flush=True)
    print(json.dumps(summ["verdict_dist"], ensure_ascii=False), flush=True)
    print("by_set:", json.dumps(summ["by_set"], ensure_ascii=False), flush=True)
    print("by_cn:", json.dumps({k: dict(v) for k, v in by_cn.items()}, ensure_ascii=False), flush=True)
    print("blocking:", json.dumps(summ["undetermined_blocking"], ensure_ascii=False), flush=True)
    print(f"candidate_subclone survivors: {len(cand_loci)}", flush=True)
    print("DONE -> nact_results.json + nact_summary.json", flush=True)


if __name__ == "__main__":
    main()
