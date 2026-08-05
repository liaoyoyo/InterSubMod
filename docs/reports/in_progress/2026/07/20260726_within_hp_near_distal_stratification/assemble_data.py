#!/usr/bin/env python3
"""Assemble data.json for the within-HP near/distal stratification report.

Reads only the archived source JSONs under sources/ (produced by the analysis
runs on 2026-07-26). Performs no new statistics beyond re-deriving summary
counts from the stored per-unit records, so every number on the page is
traceable to a file on disk (CLAUDE.md 13-A).
"""
import json
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "sources")


def load(name):
    with open(os.path.join(SRC, name)) as fh:
        return json.load(fh)


step1_units = load("within_hp_test.json")
nd = load("within_hp_near_distal.json")
decay = load("within_hp_decay.json")
ctrl = load("tp_fp_control.json")

# ---- Step 1: marginal within-HP test (per-read mean over whole +-5kb window)
s1 = [u for u in step1_units if "p" in u]
sig = [u for u in s1 if u["p"] < 0.05]
big = [u for u in sig if abs(u["d"]) >= 0.1]
both_loci = {}
for u in s1:
    both_loci.setdefault(u["loc"], []).append(u)
both = [v for v in both_loci.values() if len(v) == 2]
step1 = {
    "n_unit": len(s1),
    "n_sig": len(sig),
    "pct_sig": len(sig) / len(s1) * 100,
    "n_big": len(big),
    "pct_big": len(big) / len(s1) * 100,
    "median_abs_d": st.median([abs(u["d"]) for u in s1]),
    "median_abs_d_sig": st.median([abs(u["d"]) for u in sig]),
    "n_both_hp": len(both),
    "n_both_sig": sum(1 for v in both if all(x["p"] < 0.05 for x in v)),
    "n_both_concordant": sum(1 for v in both if v[0]["d"] * v[1]["d"] > 0),
    "top": [
        {"loc": u["loc"], "hp": u["hp"], "d": u["d"], "p": u["p"],
         "n_alt": u["n_alt"], "n_ref": u["n_ref"]}
        for u in sorted(sig, key=lambda x: -abs(x["d"]))[:5]
    ],
}

# ---- Step 2: CpG-count matching kills the near>distal gap
step2 = dict(decay["matched"])

# ---- Step 3/4: distance-decay curves, observed vs permutation null, TP vs FP
step3 = {"tp": ctrl["tp"]["curve"], "fp": ctrl["fp"]["curve"]}
step4 = {
    "tp": dict(ctrl["tp"]["pooled"], scanned=ctrl["tp"]["scanned"],
               units=ctrl["tp"]["units"], pool=ctrl["tp"]["total_loci"]),
    "fp": dict(ctrl["fp"]["pooled"], scanned=ctrl["fp"]["scanned"],
               units=ctrl["fp"]["units"], pool=ctrl["fp"]["total_loci"]),
}
step4["delta_excess"] = step4["tp"]["exc"] - step4["fp"]["exc"]
step4["delta_ratio"] = step4["tp"]["ratio"] - step4["fp"]["ratio"]
step4["tp_yield"] = step4["tp"]["units"] / step4["tp"]["scanned"] * 100
step4["fp_yield"] = step4["fp"]["units"] / step4["fp"]["scanned"] * 100

# ---- near/distal significance rates (uncorrected, for the L3 audit trail)
u_nd = nd["unit"]
rates = {}
for tag in ("near", "distal"):
    s = [u for u in u_nd if tag in u]
    sg = [u for u in s if u[tag]["p"] < 0.05]
    rates[tag] = {
        "n": len(s), "n_sig": len(sg), "pct": len(sg) / len(s) * 100,
        "median_abs_d": st.median([abs(u[tag]["d"]) for u in s]),
    }

data = {
    "meta": {
        "title": "within-HP near vs distal 甲基分層",
        "date": "2026-07-26",
        "sample": "HCC1395",
        "archive": ("/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/"
                    "20260119_all-with-w5000_1"),
        "sources": sorted(os.listdir(SRC)),
    },
    "design": {
        "window": "±5,000 bp",
        "min_reads_per_group": 6,
        "min_cov_per_cpg": 4,
        "n_perm": 3,
        "hp_rule": "germline family = HP tag 破折號前第一碼（1-1 與 1 同族）",
        "stat": "Mann-Whitney U (two-sided)",
    },
    "step1": step1,
    "step2": step2,
    "step3": step3,
    "step4": step4,
    "rates": rates,
}

out = os.path.join(HERE, "20260726_within_hp_near_distal.data.json")
with open(out, "w") as fh:
    json.dump(data, fh, ensure_ascii=False, indent=1)
print("wrote", out)
for k in ("step1", "step2", "step4"):
    print(" ", k, json.dumps(data[k], ensure_ascii=False)[:220])
