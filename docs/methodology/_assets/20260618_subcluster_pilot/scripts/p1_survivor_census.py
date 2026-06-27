#!/usr/bin/env python3
"""[P1.3 survivor 決定性 census] 對 102 candidate_subclone survivors:
(1) genotype-homogeneity: cluster_split 用的 tumor reads 中 dominant (hp,alt) frac + 最佳 a-priori 軸平衡 → 分
    genotype_homogeneous(無 a-priori 軸,純 double-dip collider) vs genotype_anchored(有≥4/side 軸,可救少數);
(2) candidate vs cis_asm 對比(germline_baseline/struct_dbeta/containment/dir_concordance median) → 獨立驗證
    'tumor signature 閘不判別 + R1 結構性必然' claim;
(3) R3 frac_normal_bimodal 門檻敏感度(count@0.05-0.30) → 誠實 count fragility。
重用 p1_nact_cistest load_region/cluster_split_groups/somatic_axis_groups。純讀。輸出 survivor_census.json。"""
import sys, json
from collections import Counter
import numpy as np
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot/scripts")
import p1_nact_cistest as P

A = P.A
res = json.load(open(f"{A}/nact_results.json"))
cands = json.load(open(f"{A}/cis_candidates_resolved.json"))
RD = {f"{c['chrom']}:{c['pos']}:{c['set']}": c["region_dir"] for c in cands}
survivors = [r for r in res if r["nact_verdict"] == "candidate_subclone"]
cis = [r for r in res if r["nact_verdict"] == "cis_asm"]


def axis_balance(meta, T):
    """3 a-priori 軸的最佳 min(side); 回 (best_balance, axis)。"""
    best = (0, None)
    for ax, g0f, g1f in (("HP2_vs_HP2-1", lambda m: m["hp"] == "2", lambda m: m["hp"] == "2-1"),
                         ("HP1_vs_HP1-1", lambda m: m["hp"] == "1", lambda m: m["hp"] == "1-1"),
                         ("REF_vs_ALT", lambda m: m["alt"] == "REF", lambda m: m["alt"] == "ALT")):
        n0 = sum(1 for r in T if g0f(meta[r])); n1 = sum(1 for r in T if g1f(meta[r]))
        bal = min(n0, n1)
        if bal > best[0]:
            best = (bal, ax)
    return best


census = []
for r in survivors:
    rd = RD.get(f"{r['chrom']}:{r['pos']}:{r['set']}")
    reg = P.load_region(rd) if rd else None
    if reg is None:
        census.append({**{k: r[k] for k in ("chrom", "pos", "set", "cat8", "cn_state")}, "geno": "NO_REGION"}); continue
    M, cpgs, meta, phylo = reg
    T = [x for x in M if meta.get(x, {}).get("is_tumor") == "1"]
    cs = P.cluster_split_groups(phylo, T)
    grp = (cs[0] + cs[1]) if cs else T  # cluster_split 用的 reads, 否則全 tumor
    genos = [(meta[x]["hp"], meta[x]["alt"]) for x in grp]
    dom = Counter(genos).most_common(1)[0][1] if genos else 0
    dom_frac = dom / max(1, len(genos))
    bal, ax = axis_balance(meta, grp)
    geno = "anchored" if bal >= 4 else "homogeneous"
    census.append({"chrom": r["chrom"], "pos": r["pos"], "set": r["set"], "cat8": r["cat8"], "cn_state": r["cn_state"],
                   "track_a_axis": r.get("track_a_axis"), "n_grp": len(grp), "dominant_geno_frac": round(dom_frac, 3),
                   "best_apriori_balance": bal, "best_apriori_axis": ax, "geno": geno,
                   "containment": r.get("containment"), "direction_concordance": r.get("direction_concordance"),
                   "germline_baseline_strength": r.get("germline_baseline_strength"),
                   "struct_dbeta_median": r.get("struct_dbeta_median"), "frac_normal_bimodal": r.get("frac_normal_bimodal")})


def med(items, k):
    v = [x[k] for x in items if x.get(k) is not None]
    return round(float(np.median(v)), 3) if v else None


# candidate vs cis_asm 對比(獨立驗證)
cmp = {"candidate": {"n": len(survivors), "struct_dbeta_median": med(survivors, "struct_dbeta_median"),
                     "germline_baseline_strength": med(survivors, "germline_baseline_strength"),
                     "containment": med(survivors, "containment"), "direction_concordance": med(survivors, "direction_concordance"),
                     "frac_normal_bimodal": med(survivors, "frac_normal_bimodal")},
       "cis_asm": {"n": len(cis), "struct_dbeta_median": med(cis, "struct_dbeta_median"),
                   "germline_baseline_strength": med(cis, "germline_baseline_strength"),
                   "containment": med(cis, "containment"), "direction_concordance": med(cis, "direction_concordance"),
                   "frac_normal_bimodal": med(cis, "frac_normal_bimodal")}}
gb = [c["germline_baseline_strength"] for c in census if c.get("germline_baseline_strength") is not None]
# geno census
geno_counts = Counter(c["geno"] for c in census)
anchored = [c for c in census if c["geno"] == "anchored"]
# R3 門檻敏感度
fb = [r.get("frac_normal_bimodal") for r in survivors if r.get("frac_normal_bimodal") is not None]
thr_sens = {str(t): sum(1 for x in fb if x < t) for t in (0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30)}

out = {"n_survivors": len(survivors), "geno_census": dict(geno_counts),
       "geno_anchored_loci": [{k: c[k] for k in ("chrom", "pos", "set", "cat8", "cn_state", "best_apriori_axis", "best_apriori_balance", "dominant_geno_frac")} for c in anchored],
       "candidate_vs_cis_comparison": cmp,
       "survivor_germline_baseline": {"median": round(float(np.median(gb)), 3) if gb else None,
                                      "ge_0.2": sum(1 for x in gb if x >= 0.2), "lt_0.05": sum(1 for x in gb if x < 0.05), "n": len(gb)},
       "R3_threshold_sensitivity": thr_sens,
       "dominant_geno_frac_dist": {"median": round(float(np.median([c["dominant_geno_frac"] for c in census if "dominant_geno_frac" in c])), 3),
                                   "eq_1.0": sum(1 for c in census if c.get("dominant_geno_frac") == 1.0),
                                   "ge_0.9": sum(1 for c in census if (c.get("dominant_geno_frac") or 0) >= 0.9)},
       "per_survivor": census}
json.dump(out, open(f"{A}/survivor_census.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("per_survivor", "geno_anchored_loci")}, ensure_ascii=False, indent=1))
print(f"\ngeno_anchored loci ({len(anchored)}):")
for c in anchored:
    print(f"  {c['chrom']}:{c['pos']} {c['set']} {c['cat8']} axis={c['best_apriori_axis']} bal={c['best_apriori_balance']} dom_frac={c['dominant_geno_frac']}")
print("[-> survivor_census.json]")
