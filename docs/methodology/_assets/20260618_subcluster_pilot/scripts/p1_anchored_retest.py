#!/usr/bin/env python3
"""[P1.3 anchored 30 genotype-軸 re-test] 對 genotype-anchored 30 survivors,
強制 track_a_axis = best a-priori 軸(HP1/HP1-1 等)而非 cluster_split → 非循環 signature → NACT Track B + verdict。
測『anchored 在 genotype 軸上是 cis-ASM(characterization) 還是真 tumor-specific』。重用 p1_nact_cistest。純讀。"""
import sys, json
import numpy as np
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot/scripts")
import p1_nact_cistest as P

A = P.A
census = json.load(open(f"{A}/survivor_census.json"))
anchored = census["geno_anchored_loci"]
cands = json.load(open(f"{A}/cis_candidates_resolved.json"))
RD = {f"{c['chrom']}:{c['pos']}:{c['set']}": c["region_dir"] for c in cands}
recs = json.load(open(f"{A}/phylo_cpp_wg_full_records_v6.json"))
P._REC = {f"{r['chrom']}:{r['pos']}": r for r in recs}

AXMAP = {"HP1_vs_HP1-1": (lambda m: m["hp"] == "1", lambda m: m["hp"] == "1-1", "1"),
         "HP2_vs_HP2-1": (lambda m: m["hp"] == "2", lambda m: m["hp"] == "2-1", "2"),
         "REF_vs_ALT": (lambda m: m["alt"] == "REF", lambda m: m["alt"] == "ALT", None)}

out = []
for c in anchored:
    rd = RD.get(f"{c['chrom']}:{c['pos']}:{c['set']}")
    reg = P.load_region(rd) if rd else None
    if reg is None:
        out.append({**c, "retest": "NO_REGION"}); continue
    M, cpgs, meta, phylo = reg
    K = len(cpgs)
    T = [x for x in M if meta.get(x, {}).get("is_tumor") == "1"]
    ax = c["best_apriori_axis"]; f0, f1, parent_hp = AXMAP[ax]
    gA = [x for x in T if f0(meta[x])]; gB = [x for x in T if f1(meta[x])]
    if len(gA) < 4 or len(gB) < 4:
        out.append({**c, "retest": "AXIS_THIN"}); continue
    d, p = P.per_cpg_diff(M, K, gA, gB)
    S = P.sig_set(d, p)
    run = P.coherent_run(S, d, cpgs) if S else 0
    if not S or run < P.RUN_MIN:
        out.append({**c, "retest": "NO_GENO_SIGNATURE", "S_geno": len(S), "run": run}); continue
    rec = P._REC.get(f"{c['chrom']}:{c['pos']}", {})
    is_lg = rec.get("is_loh") is True or str(rec.get("cn_state")) in ("loh", "gain")
    tb = P.track_b(M, cpgs, meta, S, d, is_lg)
    ta = {"status": "OK", "S": S, "d_struct": d}
    v, reason, br = P.verdict(ta, tb, is_lg)
    out.append({**c, "retest": "OK", "S_geno": len(S), "run": run,
                "geno_dbeta_median": round(float(np.median([abs(d[j]) for j in S])), 3),
                "R1": tb.get("R1_call"), "R2": tb.get("R2_call"), "R3": tb.get("R3_call"),
                "containment": tb.get("containment"), "germline_baseline_strength": tb.get("germline_baseline_strength"),
                "geno_axis_verdict": v, "reason": reason})

from collections import Counter
status = Counter(o.get("retest") for o in out)
gv = Counter(o.get("geno_axis_verdict") for o in out if o.get("retest") == "OK")
summ = {"n_anchored": len(anchored), "retest_status": dict(status), "geno_axis_verdict": dict(gv),
        "still_candidate_on_geno_axis": [{"chrom": o["chrom"], "pos": o["pos"], "set": o["set"], "axis": o["best_apriori_axis"],
                                          "S_geno": o.get("S_geno"), "geno_dbeta_median": o.get("geno_dbeta_median")}
                                         for o in out if o.get("geno_axis_verdict") == "candidate_subclone"], "per_locus": out}
json.dump(summ, open(f"{A}/anchored_retest.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps({k: v for k, v in summ.items() if k != "per_locus"}, ensure_ascii=False, indent=1))
print("[-> anchored_retest.json]")
