#!/usr/bin/env python3
"""[merge GW linkage parts] 合併 sm_linkage_gw_part_{1..5}.json → sm_linkage_genomewide.json。
chroms 互斥故 census 直接併；pairs concat；universe/tally/counts 累加。"""
import json
import glob
from collections import Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = f"{A}/sm_linkage_genomewide.json"

parts = sorted(glob.glob(f"{A}/sm_linkage_gw_part_*.json"))
assert parts, "no part files found"
pairs = []
census = {}
tally = Counter()
uni = {"n_tp": 0, "n_fp": 0, "n_union": 0, "per_chrom": {}}
n_groups_multi = n_singletons = 0
params = None
for p in parts:
    d = json.load(open(p))
    params = params or d["params"]
    pairs.extend(d["pairs"])
    census.update(d["census"])
    for k in ("n_tp", "n_fp", "n_union"):
        uni[k] += d["universe"][k]
    uni["per_chrom"].update(d["universe"]["per_chrom"])
    for k, v in d["tally_powered_somatic"].items():
        tally[k] += v
    n_groups_multi += d["n_groups_multi"]
    n_singletons += d["n_singletons"]

out = {"params": params, "universe": uni, "n_groups_multi": n_groups_multi,
       "n_singletons": n_singletons, "n_pairs_recorded": len(pairs),
       "tally_powered_somatic": dict(tally), "pairs": pairs, "census": census,
       "merged_from": [p.split("/")[-1] for p in parts]}
json.dump(out, open(OUT, "w"), ensure_ascii=False)
print(f"MERGED {len(parts)} parts: union={uni['n_union']} (TP={uni['n_tp']} FP={uni['n_fp']}) "
      f"census={len(census)} pairs={len(pairs)} tally={dict(tally)} -> {OUT}")
