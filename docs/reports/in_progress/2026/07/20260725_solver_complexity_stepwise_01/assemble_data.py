#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_data.py — 彙整 solver 複雜度分析的全部已驗證數字 → data.json

來源(全部本輪實跑,可重算):
  1. complexity_bench.py   → 28 案例逐階段實測(28/28 鏡像驗證通過)
  2. complexity_theory.py  → 大整數精確理論上界
  3. layered_reconstruction_HCC1395.json → 真實問題規模分佈(engineering baseline)

§13-A:HTML 只從本檔注入數字,不手打。
"""
import collections
import json
import statistics as st
from math import comb
from pathlib import Path

_ARCHIVE = Path(__file__).parent / "sources"
_TMP = Path("/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/"
            "efb6e3d8-c4af-43d8-ac97-9dffdbec60ed/scratchpad")
# 優先用歸檔副本(可重現);scratchpad 僅作生成當下的 fallback
SCRATCH = _ARCHIVE if (_ARCHIVE / "complexity_bench.json").exists() else _TMP
REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
OUTDIR = REPO / "docs/reports/in_progress/2026/07/20260725_solver_complexity_stepwise_01"
REAL = REPO / "docs/methodology/_assets/20260618_subcluster_pilot/layered_reconstruction_HCC1395.json"

MAXK, EXTRA_CAP, BUDGET = 8, 4, 150_000

D = {"_meta": {
    "title": "布爾超立方體上「最小 group-Steiner 樹集合」求解 — 逐步驟計算量與理論上界",
    "generated_from": ["complexity_bench.py", "complexity_theory.py",
                       "layered_reconstruction_HCC1395.json"],
    "solver": "docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/"
              "tree_enumeration_solver.py (2026-07-06)",
    "params": {"MAX_SNV": MAXK, "extra_cap": EXTRA_CAP, "per_level_budget": BUDGET,
               "production_ANALYSIS_TREE_CAP": 0, "production_VERIFY_EVERY": 1},
}}

# ---------------- 1. 實測 bench ----------------
bench = json.load(open(SCRATCH / "complexity_bench.json"))
cases = [c for c in bench["cases"] if "error" not in c]
D["bench"] = {
    "n_cases": len(cases),
    "n_mirror_ok": sum(1 for c in cases if c.get("mirror_ok")),
    "n_capped": sum(1 for c in cases if c.get("capped")),
}

STAGES = [("t3_subcubes", "S3 建 subcube group"), ("t4_pool", "S4 建候選 pool"),
          ("t5_level_search", "S5 逐層最小解搜尋"), ("t6_analytic_count", "S6 分析式計數"),
          ("t7_materialize", "S7 實體化樹")]
tot = {k: sum(c[k] for c in cases) for k, _ in STAGES}
s_all = sum(tot.values())
D["stage_share"] = [{"key": k, "label": lab, "ms": round(tot[k] * 1000, 2),
                     "pct": round(tot[k] / s_all * 100, 2)} for k, lab in STAGES]
D["stage_share_total_ms"] = round(s_all * 1000, 2)

rates = [c["t5_level_search"] / c["candidates_tested"] * 1e6
         for c in cases if c["candidates_tested"] >= 1000]
D["throughput_us_per_candidate"] = {
    "min": round(min(rates), 2), "median": round(st.median(rates), 2),
    "max": round(max(rates), 2), "n_samples": len(rates)}

D["bench_cases"] = [{
    "case": c["case"], "k": c["k"], "f": c["n_full"], "p": c["n_partial"],
    "P": c["P_extra_pool"], "e_min": c["e_min"], "capped": bool(c["capped"]),
    "cand": c["candidates_tested"], "n_feasible_N": c["n_feasible_N"],
    "n_trees": c["n_trees"],
    "t5_ms": round(c["t5_level_search"] * 1000, 2),
    "t7_ms": round(c["t7_materialize"] * 1000, 2),
    "cap_kind": (None if not c["capped"] else
                 ("budget" if "budget" in (c["cap_reason"] or "") else "extra_cap_greedy")),
} for c in cases]

v45 = [{"case": c["case"], "t5_ms": round(c["t5_level_search"] * 1000, 2),
        "V4_ms": round(c["t_verify_V4"] * 1000, 2),
        "V5_ms": round(c["t_verify_V5"] * 1000, 2),
        "ratio": round((c["t_verify_V4"] + c["t_verify_V5"]) / c["t5_level_search"], 2)}
       for c in cases if c.get("t_verify_V5") is not None and c["t5_level_search"] > 0.001]
D["verify_cost"] = {"rows": v45,
                    "ratio_min": min(r["ratio"] for r in v45),
                    "ratio_max": max(r["ratio"] for r in v45)}
D["pool_build_max"] = {
    "max_subset_ops": max(c["pool_subset_ops"] for c in cases),
    "max_t4_ms": round(max(c["t4_pool"] for c in cases) * 1000, 2)}

# ---------------- 2. 理論上界 ----------------
def E_of_P(P):
    for e in range(0, min(EXTRA_CAP, P) + 1):
        if comb(P, e) > BUDGET:
            return e - 1
    return min(EXTRA_CAP, P)

curve = []
for P in range(0, 2 ** MAXK):
    E = E_of_P(P)
    tot_c = sum(comb(P, e) for e in range(0, E + 1)) if E >= 0 else 0
    curve.append({"P": P, "E": E, "cand": tot_c})
best = max(curve, key=lambda r: r["cand"])
D["step5_bound"] = {
    "formula": "Cand(P) = sum_{e=0}^{E(P)} C(P,e),  E(P)=min(4, min{e: C(P,e)>150000} - 1)",
    "global_max": best["cand"], "at_P": best["P"], "at_E": best["E"],
    "curve": curve,
    "max_P_allowing_e4": max(r["P"] for r in curve if r["E"] >= 4),
    "max_P_allowing_e3": max(r["P"] for r in curve if r["E"] >= 3),
    "sawtooth_pairs": [
        {"P": 45, "E": 4, "cand": next(r["cand"] for r in curve if r["P"] == 45)},
        {"P": 46, "E": 3, "cand": next(r["cand"] for r in curve if r["P"] == 46)},
        {"P": 97, "E": 3, "cand": next(r["cand"] for r in curve if r["P"] == 97)},
        {"P": 98, "E": 2, "cand": next(r["cand"] for r in curve if r["P"] == 98)},
        {"P": 255, "E": 2, "cand": next(r["cand"] for r in curve if r["P"] == 255)},
    ],
}
D["step5_bound"]["worst_walltime_s"] = {
    q: round(best["cand"] * D["throughput_us_per_candidate"][q] / 1e6, 3)
    for q in ("min", "median", "max")}

def A_max(k):
    v = 1
    for j in range(1, k + 1):
        v *= j ** comb(k, j)
    return v

D["A_max"] = [{"k": k, "vertices": 2 ** k, "value": str(A_max(k)),
               "digits": len(str(A_max(k))),
               "sci": f"{float(A_max(k)):.3e}" if len(str(A_max(k))) < 300 else "—"}
              for k in range(1, MAXK + 1)]
D["A_max_k8_factors"] = [{"j": j, "exp": comb(MAXK, j)} for j in range(1, MAXK + 1)]
D["reachable_bound"] = [{"f": f, "bound": str(MAXK ** (f + EXTRA_CAP)),
                         "digits": len(str(MAXK ** (f + EXTRA_CAP)))}
                        for f in [1, 2, 3, 5, 8, 12]]
D["state_space"] = [{"k": k, "vertices": 2 ** k, "max_P": 2 ** k - 1}
                    for k in range(2, MAXK + 1)]

# ---------------- 3. 真實分佈 ----------------
u = json.load(open(REAL))["detail"]
n = len(u)

def hist(key, cap=None):
    c = collections.Counter(x.get(key) for x in u)
    rows = []
    for v, ct in sorted(c.items(), key=lambda t: (t[0] is None, t[0])):
        if cap is not None and isinstance(v, int) and v > cap:
            continue
        rows.append({"v": v, "n": ct, "pct": round(ct / n * 100, 2)})
    return rows

def stats(key):
    vals = [x[key] for x in u if isinstance(x.get(key), (int, float))]
    return {"min": min(vals), "median": int(st.median(vals)),
            "mean": round(st.mean(vals), 2), "max": max(vals)}

nt = [x["n_trees"] for x in u if isinstance(x.get("n_trees"), int)]
cap_reasons = collections.Counter()
for x in u:
    if x.get("capped"):
        r = x.get("cap_reason") or ""
        cap_reasons["budget" if "budget" in r else "extra_cap_greedy"] += 1

D["real"] = {
    "source": "docs/methodology/_assets/20260618_subcluster_pilot/layered_reconstruction_HCC1395.json",
    "caveat": "2026-07-11 稽核判定為 upstream-mismatched engineering baseline"
              "(6/7 歷史 tagged BAM 受 --truth-bed 限制)。此處僅用於『計算規模量級』參考,"
              "不得作為正式 Results 的比例數字。",
    "n_units": n,
    "k": {"hist": hist("n_sSNV"), "stats": stats("n_sSNV")},
    "f": {"hist": hist("n_full_pops"), "stats": stats("n_full_pops")},
    "p": {"hist": hist("n_partial", cap=12), "stats": stats("n_partial")},
    "e_min": {"hist": hist("n_hidden", cap=11), "stats": stats("n_hidden")},
    "n_trees": {"stats": {"min": min(nt), "median": int(st.median(nt)),
                          "mean": round(st.mean(nt), 2), "max": max(nt)},
                "total": sum(nt),
                "tail": [{"thr": t, "n": sum(1 for v in nt if v >= t),
                          "pct": round(sum(1 for v in nt if v >= t) / len(nt) * 100, 2)}
                         for t in [1, 2, 4, 8, 32, 100, 1000]]},
    "capped": {"n": sum(1 for x in u if x.get("capped")),
               "pct": round(sum(1 for x in u if x.get("capped")) / n * 100, 2),
               "by_reason": dict(cap_reasons)},
    "L1_class": [{"cls": k2, "n": v, "pct": round(v / n * 100, 2)}
                 for k2, v in collections.Counter(x.get("L1_class") for x in u).most_common()],
}
D["theory_vs_real"] = {
    "A_max_k8_digits": len(str(A_max(8))),
    "real_max_trees": max(nt),
    "real_total_trees": sum(nt),
    "gap_orders_of_magnitude": len(str(A_max(8))) - len(str(max(nt))),
}

OUTDIR.mkdir(parents=True, exist_ok=True)
out = OUTDIR / "20260725_solver_complexity_stepwise_01.data.json"
with open(out, "w") as fh:
    json.dump(D, fh, indent=1, ensure_ascii=False)
print(f"wrote {out}")
print(f"  bench cases={D['bench']['n_cases']} mirror_ok={D['bench']['n_mirror_ok']}")
print(f"  step5 global_max={D['step5_bound']['global_max']:,} @P={D['step5_bound']['at_P']}")
print(f"  A_max(8) digits={D['theory_vs_real']['A_max_k8_digits']}")
print(f"  real units={n:,} max_trees={max(nt)} total_trees={sum(nt):,}")
