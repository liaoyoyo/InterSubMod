#!/usr/bin/env python3
"""Merge the 8 chunk outputs of methyl_bimodality_cn_rate per sample -> rate by CN + q5."""
import glob
import json
from collections import defaultdict, Counter

BASE = "/big7_disk/liaoyoyo2001/cnv_sv_work/methyl_lineage"
SAMPLES = ["H1437", "H2009", "HCC1954"]
out_all = []
for s in SAMPLES:
    parts = sorted(glob.glob(f"{BASE}/cnrate_{s}_part*.json"))
    raw = defaultdict(lambda: [0, 0])
    flag = defaultdict(Counter)
    q5 = []
    nreg = 0
    for p in parts:
        d = json.load(open(p))
        nreg += d["summary"]["n_regions_multiSNV"]
        for cn, v in d.get("raw_rate", {}).items():
            raw[cn][0] += v[0]
            raw[cn][1] += v[1]
        for cn, fl in d.get("flag_by_cn", {}).items():
            for k, n in fl.items():
                flag[cn][k] += n
        q5.extend(d.get("q5_candidates", []))
    rate = {cn: {"n_bimodal": v[0], "n_tested": v[1],
                 "bimodal_rate": round(v[0] / v[1], 4) if v[1] else None}
            for cn, v in sorted(raw.items())}
    merged = {"sample": s, "n_parts": len(parts), "n_regions_multiSNV": nreg,
              "bimodal_rate_by_cn": rate, "flag_by_cn": {k: dict(v) for k, v in flag.items()},
              "n_q5": len(q5), "q5_candidates": q5}
    json.dump(merged, open(f"{BASE}/cnrate_{s}_merged.json", "w"), ensure_ascii=False, indent=1)
    out_all.append(merged)

print(f"{'sample':9} {'gain_rate':>10} {'neutral_rate':>13} {'loh_rate':>9} {'loss_rate':>10} {'Q5':>5}  (bimodal/tested)")
for m in out_all:
    r = m["bimodal_rate_by_cn"]
    g = lambda cn: (f"{r[cn]['bimodal_rate']:.3f}({r[cn]['n_bimodal']}/{r[cn]['n_tested']})" if cn in r else "-")
    print(f"{m['sample']:9} {g('gain'):>10} {g('neutral'):>13} {g('loh'):>9} {g('loss'):>10} {m['n_q5']:>5}")
print("\n=== 因果判讀: gain_rate vs neutral_rate (gain>>neutral = CN-gain 因果提高雙峰率) ===")
for m in out_all:
    r = m["bimodal_rate_by_cn"]
    gr = r.get("gain", {}).get("bimodal_rate")
    nr = r.get("neutral", {}).get("bimodal_rate")
    if gr is not None and nr is not None:
        ratio = round(gr / nr, 2) if nr else None
        print(f"  {m['sample']}: gain {gr} / neutral {nr} = {ratio}x  -> {'CN-gain 顯著提高' if ratio and ratio>1.5 else '差異小' if ratio else '?'}")
