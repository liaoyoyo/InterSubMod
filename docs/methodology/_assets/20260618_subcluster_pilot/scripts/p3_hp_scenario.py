#!/usr/bin/env python3
"""[HP 情境驗證] 對 102 candidate + 9 tumor-specific + 4 confirmed, 統計 tumor reads HP 組成
→ 分類用戶三情境: H1(單 HP 多甲基群=double-dip) / H2(同 HP 內 2nd sSNV 互斥=subclone) / H3(HP1-1 與 HP2-1 並存=cross-parental allelic)。
+ 對 H3 loci 檢 HP1-1 vs HP2-1 甲基差是否=germline HP1 vs HP2 cis-ASM。純讀。"""
import sys, json
from collections import Counter
import numpy as np
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot/scripts")
import p1_nact_cistest as P

A = P.A
res = json.load(open(f"{A}/nact_results.json"))
cands = json.load(open(f"{A}/cis_candidates_resolved.json"))
RD = {f"{c['chrom']}:{c['pos']}:{c['set']}": c["region_dir"] for c in cands}
nine = {f"{r['chrom']}:{r['pos']}" for r in json.load(open(f"{A}/anchored_retest.json"))["still_candidate_on_geno_axis"]}
conf = set(json.load(open(f"{A}/p2_linkage.json")).keys())
survivors = [r for r in res if r["nact_verdict"] == "candidate_subclone"]


def hp_compose(rd):
    reg = P.load_region(rd)
    if reg is None:
        return None
    M, cpgs, meta, phylo = reg
    T = [x for x in M if meta.get(x, {}).get("is_tumor") == "1"]
    hpc = Counter(meta[x]["hp"] for x in T)
    return hpc, len(T)


def scenario(hpc):
    h1, h11, h2, h21, h0 = hpc.get("1", 0), hpc.get("1-1", 0), hpc.get("2", 0), hpc.get("2-1", 0), hpc.get("0", 0)
    # somatic-carrying reads
    som = {"1-1": h11, "2-1": h21}
    cross_parental = h11 >= 4 and h21 >= 4                 # H3: 兩 parental 都有 somatic
    within_par = (h1 >= 4 and h11 >= 4) or (h2 >= 4 and h21 >= 4)  # 同 parental REF/ALT 可分
    tot = sum(hpc.values())
    dom = max(hpc.values()) / tot if tot else 0
    homog = dom >= 0.85 and not within_par               # 單一群主導且無同-parental 軸
    return {"HP1": h1, "HP1-1": h11, "HP2": h2, "HP2-1": h21, "HP0": h0,
            "cross_parental_H3": cross_parental, "within_parental": within_par, "homogeneous_H1": homog, "dom_frac": round(dom, 2)}


rows = []
for r in survivors:
    key = f"{r['chrom']}:{r['pos']}"
    rd = RD.get(f"{key}:{r['set']}")
    hc = hp_compose(rd) if rd else None
    if not hc:
        continue
    sc = scenario(hc[0])
    rows.append({"key": key, "set": r["set"], "cat8": r["cat8"], "in9": key in nine, "in4conf": key in conf, **sc})

n = len(rows)
summ = {"n_survivors_scanned": n,
        "H3_cross_parental(HP1-1+HP2-1)": sum(1 for x in rows if x["cross_parental_H3"]),
        "within_parental_somatic_axis": sum(1 for x in rows if x["within_parental"]),
        "H1_homogeneous(double-dip risk)": sum(1 for x in rows if x["homogeneous_H1"]),
        "examples_H3": [{"key": x["key"], "HP1-1": x["HP1-1"], "HP2-1": x["HP2-1"]} for x in rows if x["cross_parental_H3"]][:8],
        "examples_within_parental": [{"key": x["key"], "in9": x["in9"], "HP1": x["HP1"], "HP1-1": x["HP1-1"], "HP2": x["HP2"], "HP2-1": x["HP2-1"]}
                                     for x in rows if x["within_parental"]][:8],
        "the_9_tumorspecific": [{"key": x["key"], "HP": {k: x[k] for k in ("HP1", "HP1-1", "HP2", "HP2-1", "HP0")},
                                 "within_par": x["within_parental"], "homog": x["homogeneous_H1"], "in4conf": x["in4conf"]}
                                for x in rows if x["in9"]]}
json.dump({"summary": summ, "per_survivor": rows}, open(f"{A}/hp_scenario.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(summ, ensure_ascii=False, indent=1))
print("[-> hp_scenario.json]")
