#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
可驗證觀察(純數 read,零統計推論):單讀脆弱判定是否集中在 低頻VAF / 太遠 / 低coread。
弱 = 決定拓樸的那格只有 1 條 read。
  nested_a_in_b 決定格=RA; nested_b_in_a=AR; independent=min(RA,AR)。
輸出 weak vs strong 在 VAF/dist/coread 的分布 → 供定義 confidence 準則。
用法：python3 observe_weak_pairs.py
"""
import os, csv, json, statistics as st
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

def load(name):
    p = os.path.join(DATA, "lists", name)
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t")) if os.path.exists(p) else []

def deciding_off(d):
    """回 (config_kind, deciding_off_count)。"""
    ra, ar = int(d["RA"]), int(d["AR"]); cfg = d["config"]
    if cfg == "nested_a_in_b": return "nested", ra
    if cfg == "nested_b_in_a": return "nested", ar
    if cfg == "independent": return "independent", min(ra, ar)
    return cfg, None

pairs = []
for c in ("nested_a_in_b", "nested_b_in_a", "independent"):
    for hp in ("sameHP", "diffHP"):
        pairs += load(f"{c}_{hp}.tsv")

def minvaf(d):
    try: return min(float(d["vaf_a"]), float(d["vaf_b"]))
    except: return None

def bins(vals, edges, labels):
    c = Counter()
    for v in vals:
        if v is None: continue
        for e, l in zip(edges, labels):
            if v < e: c[l] += 1; break
        else: c[labels[-1]] += 1
    return dict(c)

def summarize(rows, tag):
    if not rows: return {}
    mv = [minvaf(d) for d in rows]; mv = [x for x in mv if x is not None]
    dist = [int(d["dist"]) for d in rows]; cor = [int(d["coread"]) for d in rows]
    return {
        "n": len(rows),
        "min_vaf_median": round(st.median(mv), 3) if mv else None,
        "min_vaf_bins(<0.05/<0.1/<0.25/<0.5/>=0.5)": bins(mv, [0.05, 0.1, 0.25, 0.5], ["<0.05", "0.05-0.1", "0.1-0.25", "0.25-0.5", ">=0.5"]),
        "dist_median": int(st.median(dist)),
        "dist_bins(<1k/<5k/<20k/<50k/>=50k)": bins(dist, [1000, 5000, 20000, 50000], ["<1k", "1-5k", "5-20k", "20-50k", ">=50k"]),
        "coread_median": int(st.median(cor)),
        "coread_bins(<6/6-10/11-20/21-50/>50)": bins(cor, [6, 11, 21, 51], ["2-5", "6-10", "11-20", "21-50", ">50"]),
    }

out = {}
for kind in ("nested", "independent"):
    sub = [d for d in pairs if deciding_off(d)[0] == kind]
    weak = [d for d in sub if deciding_off(d)[1] == 1]
    strong = [d for d in sub if (deciding_off(d)[1] or 0) >= 2]
    out[kind] = {
        "total": len(sub), "weak_single_read": len(weak), "weak_pct": round(100*len(weak)/len(sub), 1) if sub else None,
        "WEAK(決定格=1)": summarize(weak, "weak"),
        "STRONG(決定格>=2)": summarize(strong, "strong"),
    }

with open(os.path.join(DATA, "weak_pair_observation.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)

for kind in out:
    o = out[kind]
    print(f"===== {kind}: {o['total']:,} 對, 弱(單讀)={o['weak_single_read']:,} ({o['weak_pct']}%) =====")
    for grp in ("WEAK(決定格=1)", "STRONG(決定格>=2)"):
        s = o[grp]
        if not s: continue
        print(f"  [{grp}] n={s['n']:,}")
        print(f"    min_VAF 中位={s['min_vaf_median']}  bins={s['min_vaf_bins(<0.05/<0.1/<0.25/<0.5/>=0.5)']}")
        print(f"    dist 中位={s['dist_median']:,}bp  bins={s['dist_bins(<1k/<5k/<20k/<50k/>=50k)']}")
        print(f"    coread 中位={s['coread_median']}  bins={s['coread_bins(<6/6-10/11-20/21-50/>50)']}")
print("OK wrote weak_pair_observation.json")
