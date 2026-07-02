#!/usr/bin/env python3
"""[v3 dashboard 標註] records_full → 每位點加最終判別 verdict 欄(讓 dashboard 可篩):
  wstate(S1–S6, 重用 classify) / structure_tier(1 confident/2 marginal/3 none) / geom_divergence / align_status。
純讀已落檔, 零 compute。輸出 phylo_cpp_wg_full_records_v3.json。"""
import json, os, sys
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__))
from phylo_weak_structure_classify import classify  # S1–S6 互斥決策流(同邏輯)

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
TIER = {"S1_CONFIRMED_MULTI_ALIGNED": 1, "S2_CONFIRMED_MULTI_UNALIGNED": 1,
        "S3_UNSTABLE_MULTI": 2, "S4_WEAK_FINE_ONLY": 2, "S5_WEAK_RESIDUAL": 2, "S6_NO_STRUCTURE_CLEAN": 3}
WSHORT = {"S1_CONFIRMED_MULTI_ALIGNED": "S1", "S2_CONFIRMED_MULTI_UNALIGNED": "S2", "S3_UNSTABLE_MULTI": "S3",
          "S4_WEAK_FINE_ONLY": "S4", "S5_WEAK_RESIDUAL": "S5", "S6_NO_STRUCTURE_CLEAN": "S6"}


def main():
    rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_full.json"))
    for r in rec:
        w = classify(r)
        r["wstate"] = WSHORT.get(w, w)
        r["structure_tier"] = TIER.get(w, 3)
        g = r.get("geometry_ng")
        r["geom_divergence"] = 1 if (g and g >= 2 and r["coarse_ng"] < 2) else 0
        r["align_status"] = ("aligned" if r["aligned"] else "unaligned") if r["coarse_ng"] >= 2 else "NA"
    json.dump(rec, open(f"{A}/phylo_cpp_wg_full_records_v3.json", "w"))
    out = {"n": len(rec),
           "wstate": dict(Counter(r["wstate"] for r in rec)),
           "structure_tier": dict(Counter(r["structure_tier"] for r in rec)),
           "align_status": dict(Counter(r["align_status"] for r in rec)),
           "geom_divergence": sum(r["geom_divergence"] for r in rec),
           "wstate_sumcheck": sum(Counter(r["wstate"] for r in rec).values()) == len(rec)}
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
