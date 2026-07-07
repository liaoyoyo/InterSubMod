#!/usr/bin/env python3
"""Merge the 8 chunk outputs of methyl_bimodality_cn_rate per sample -> rate by CN,
flag, Q5, AND determinacy x methyl-status cross-tab (complete, incl E_subcube)."""
import glob
import json
from collections import defaultdict, Counter

BASE = "/big7_disk/liaoyoyo2001/cnv_sv_work/methyl_lineage"
SAMPLES = ["H1437", "H2009", "HCC1954"]
ORDER = ["A_determined", "B_pairwise", "E_subcube", "A_ambiguous",
         "C_underdetermined", "recurrence_required", "incompatible", "other"]


def merge_sample(s):
    parts = sorted(glob.glob(f"{BASE}/cnrate_{s}_part*.json"))
    raw = defaultdict(lambda: [0, 0])
    flag = defaultdict(Counter)
    det = defaultdict(Counter)
    q5 = []
    nreg = 0
    for p in parts:
        d = json.load(open(p))
        nreg += d["summary"]["n_regions_multiSNV"]
        for cn, v in d.get("raw_rate", {}).items():
            raw[cn][0] += v[0]; raw[cn][1] += v[1]
        for cn, fl in d.get("flag_by_cn", {}).items():
            for k, n in fl.items():
                flag[cn][k] += n
        for dc, st in d.get("det_ct", {}).items():
            for k, n in st.items():
                det[dc][k] += n
        q5.extend(d.get("q5_candidates", []))
    rate = {cn: {"n_bimodal": v[0], "n_tested": v[1],
                 "bimodal_rate": round(v[0] / v[1], 4) if v[1] else None}
            for cn, v in sorted(raw.items())}
    merged = {"sample": s, "n_parts": len(parts), "n_regions_multiSNV": nreg,
              "bimodal_rate_by_cn": rate, "flag_by_cn": {k: dict(v) for k, v in flag.items()},
              "det_x_methyl": {k: dict(v) for k, v in det.items()},
              "n_q5": len(q5), "q5_candidates": q5}
    json.dump(merged, open(f"{BASE}/cnrate_{s}_merged.json", "w"), ensure_ascii=False, indent=1)
    return merged


def main():
    res = [merge_sample(s) for s in SAMPLES]
    pool = defaultdict(Counter)
    for r in res:
        for dc, st in r["det_x_methyl"].items():
            for k, n in st.items():
                pool[dc][k] += n
    json.dump({"pooled_det_x_methyl": {k: dict(v) for k, v in pool.items()},
               "per_sample": {r["sample"]: r["det_x_methyl"] for r in res}},
              open(f"{BASE}/determinacy_x_methyl_complete.json", "w"), ensure_ascii=False, indent=1)

    print("=== POOLED (3 CN 樣本) — determinacy(前判斷) × 甲基狀態 [完整含 E_subcube] ===")
    print(f"{'determinacy':22}{'total':>7}{'residual':>9}{'res%':>6}{'cn_gain':>8}{'hp':>4}{'none':>7}{'Q5':>4}")
    tot = Counter()
    for dc in ORDER:
        v = pool.get(dc)
        if not v:
            continue
        for k, n in v.items():
            tot[k] += n
        t = v.get("total", 0); rn = v.get("residual", 0)
        print(f"{dc:22}{t:>7}{rn:>9}{(100*rn/t if t else 0):>5.1f}%{v.get('cn_gain',0):>8}"
              f"{v.get('hp',0):>4}{v.get('none',0):>7}{v.get('q5',0):>4}")
    t = tot.get("total", 0)
    print(f"{'ALL':22}{t:>7}{tot.get('residual',0):>9}{(100*tot.get('residual',0)/t if t else 0):>5.1f}%"
          f"{tot.get('cn_gain',0):>8}{tot.get('hp',0):>4}{tot.get('none',0):>7}{tot.get('q5',0):>4}")
    print("\nresidual=甲基可助L3(非HP非CN-gain雙峰) | cn_gain=CN擴增混淆 | none=甲基無訊號 | Q5=難建樹輔助")
    hard = ["A_ambiguous", "C_underdetermined", "recurrence_required", "incompatible"]
    ht = sum(pool[d].get("total", 0) for d in hard)
    hr = sum(pool[d].get("residual", 0) for d in hard)
    hq = sum(pool[d].get("q5", 0) for d in hard)
    print(f"\n難建樹區合計(A_ambiguous+C_underdet+recurrence+incompatible): "
          f"total={ht} residual={hr}({100*hr/ht if ht else 0:.1f}%) Q5={hq}")


if __name__ == "__main__":
    main()
