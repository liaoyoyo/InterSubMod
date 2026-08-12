#!/usr/bin/env python3
"""Accumulate copy-number-clean deep edges across H1437 / H2009 / HCC1954.

HCC1395 alone yields only 72 clean deep edges, too few to carry a lineage
claim.  The three samples whose SAVANA fit passed the BAF self-consistency
audit (H1437 7.17%, H2009 9.81%, HCC1954 5.62% violations) are the only
candidates for adding more.

The neutral gate is the one validated on HCC1395 in
savana_neutral_callability.py: total CN ~ 2 AND segment mean BAF below the LOH
threshold.  On HCC1395 that gate reached 90.29% precision at 7.94% recall once
copy number was recalibrated -- deliberately conservative, which is the right
trade for an anchor set.

Two caveats are carried into the output rather than being resolved:
  * these three samples' purity/ploidy are SAVANA's published values.  They
    passed an internal consistency test but have no external truth, so the
    90.29% precision measured on HCC1395 is an expectation, not a measurement,
    for them.
  * at lower purity the BAF of an LOH segment is pulled toward 0.5, so the LOH
    gate is weaker on HCC1954 (purity 0.66) than on H1437/H2009 (0.95).  A
    stricter variant additionally requiring minor-allele CN >= 0.5 is reported
    alongside.
"""

from __future__ import annotations

import bisect
import json
import math
import os
from collections import Counter, defaultdict

WORK = "/big7_disk/liaoyoyo2001/cnv_sv_work"
ROUND = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/samples"
)
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT = os.path.join(OUT_DIR, "multisample_clean_deep_edges.json")

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
CN_TOL = 0.5
BAF_LOH = 0.75
MINOR_MIN = 0.5

SAMPLES = {
    "H1437": f"{WORK}/H1437/savana_wgs/cna_normalhet/H1437_segmented_absolute_copy_number.tsv",
    "H2009": f"{WORK}/H2009/savana_wgs/cna_normalhet/H2009_segmented_absolute_copy_number.tsv",
    "HCC1954": f"{WORK}/HCC1954/savana_wgs/cna_normalhet/HCC1954_segmented_absolute_copy_number.tsv",
}
PURITY = {"H1437": 0.95, "H2009": 0.95, "HCC1954": 0.66}


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {
        "count": k,
        "n": n,
        "percent": round(100 * p, 2),
        "lo": round(100 * max(0.0, c - h), 2),
        "hi": round(100 * min(1.0, c + h), 2),
    }


class Segs:
    def __init__(self, path):
        self.by_chrom = defaultdict(list)
        with open(path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            col = {c: i for i, c in enumerate(hdr)}
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                chrom = p[col["chromosome"]]
                if chrom not in AUTOSOMES:
                    continue
                def num(key):
                    try:
                        return float(p[col[key]])
                    except (ValueError, IndexError, KeyError):
                        return None
                try:
                    s, e = int(p[col["start"]]), int(p[col["end"]])
                except ValueError:
                    continue
                self.by_chrom[chrom].append(
                    (s, e, num("copyNumber"), num("minorAlleleCopyNumber"), num("meanBAF"))
                )
        self.starts = {}
        for c, v in self.by_chrom.items():
            v.sort()
            self.starts[c] = [t[0] for t in v]

    def at(self, chrom, pos):
        v = self.by_chrom.get(chrom)
        if not v:
            return None
        i = bisect.bisect_right(self.starts[chrom], pos) - 1
        if i < 0:
            return None
        s, e, cnv, minor, baf = v[i]
        if s <= pos < e:
            return {"cn": cnv, "minor": minor, "baf": baf}
        return None


def classify(seg, strict=False):
    """clean / altered / unknown under the validated gate."""
    if seg is None or seg["cn"] is None:
        return "unknown"
    if abs(seg["cn"] - 2.0) > CN_TOL:
        return "altered"
    baf = seg["baf"]
    if baf is None:
        return "unknown"
    if baf >= BAF_LOH:
        return "altered"
    if strict:
        if seg["minor"] is None or seg["minor"] < MINOR_MIN:
            return "altered"
    return "clean"


def analyse(sample, seg_path):
    segs = Segs(seg_path)
    topo = f"{ROUND}/{sample}/{sample}.topology.jsonl"
    rows = []
    with open(topo) as fh:
        for ln in fh:
            u = json.loads(ln)
            if u.get("unit_status") != "ranked":
                continue
            chrom = u.get("chrom")
            pos = u.get("active_positions") or []
            if chrom not in AUTOSOMES or not pos:
                continue
            k = int(u.get("active_bit_count") or 0)
            ttc = int(u.get("total_tree_count") or 0)
            edges = u.get("representative_best_edges") or []
            deep = [e for e in edges if e.get("parent_label") != "ROOT"]
            # every active site must sit on clean ground for the unit to count
            calls = [classify(segs.at(chrom, p)) for p in pos]
            calls_strict = [classify(segs.at(chrom, p), strict=True) for p in pos]
            state = (
                "clean"
                if all(c == "clean" for c in calls)
                else "unknown"
                if any(c == "unknown" for c in calls)
                else "altered"
            )
            state_strict = (
                "clean"
                if all(c == "clean" for c in calls_strict)
                else "unknown"
                if any(c == "unknown" for c in calls_strict)
                else "altered"
            )
            rows.append(
                {
                    "chrom": chrom,
                    "k": k,
                    "ttc": ttc,
                    "determined": ttc == 1 or bool(u.get("best_tree_unique")),
                    "structure_unique": ttc == 1,
                    "af_resolved": ttc > 1 and bool(u.get("best_tree_unique")),
                    "n_edges": len(edges),
                    "n_deep": len(deep),
                    "cn_state": state,
                    "cn_state_strict": state_strict,
                }
            )

    def summarise(sub, label):
        det = [r for r in sub if r["determined"]]
        return {
            "label": label,
            "units": len(sub),
            "determined": wilson(len(det), len(sub)) if sub else None,
            "structure_unique": sum(1 for r in det if r["structure_unique"]),
            "af_resolved": sum(1 for r in det if r["af_resolved"]),
            "edges": sum(r["n_edges"] for r in det),
            "deep_edges": sum(r["n_deep"] for r in det),
            "chromosomes": len({r["chrom"] for r in det}),
        }

    res = {"sample": sample, "purity_used_by_savana": PURITY[sample],
           "ranked_units_autosomal": len(rows)}
    for gate, key in [("cn_state", ""), ("cn_state_strict", "_strict")]:
        clean = [r for r in rows if r[gate] == "clean"]
        res[f"clean{key}"] = {
            "all": summarise(clean, f"{sample} clean all"),
            "k_ge_2": summarise([r for r in clean if r["k"] >= 2], f"{sample} clean k>=2"),
            "k_ge_3": summarise([r for r in clean if r["k"] >= 3], f"{sample} clean k>=3"),
            "chrom_distribution_k_ge_2": dict(
                Counter(r["chrom"] for r in clean if r["k"] >= 2).most_common()
            ),
        }
    res["state_counts"] = dict(Counter(r["cn_state"] for r in rows))
    res["state_counts_strict"] = dict(Counter(r["cn_state_strict"] for r in rows))
    altered = [r for r in rows if r["cn_state"] == "altered" and r["k"] >= 2]
    res["altered_k_ge_2_reference"] = summarise(altered, f"{sample} altered k>=2")
    return res


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    per_sample = {}
    for s, p in SAMPLES.items():
        if not os.path.exists(p):
            per_sample[s] = {"status": "MISSING_SAVANA_CN", "path": p}
            continue
        per_sample[s] = analyse(s, p)

    # pooled with HCC1395's SEQC2-based clean set
    hcc = json.load(open(os.path.join(OUT_DIR, "clean_ground_evolution_yield.json")))
    hcc_k2 = hcc["cn_neutral"]["k_ge_2"]
    hcc_k3 = hcc["cn_neutral"]["k_ge_3"]

    pooled = {
        "HCC1395_seqc2_truth": {
            "k_ge_2_units": hcc_k2["n_units"],
            "k_ge_2_determined": hcc_k2["determined_total"]["count"],
            "k_ge_2_deep_edges": hcc_k2["determined_deep_edges"],
            "k_ge_3_units": hcc_k3["n_units"],
            "k_ge_3_determined": hcc_k3["determined_total"]["count"],
            "k_ge_3_deep_edges": hcc_k3["determined_deep_edges"],
            "basis": "SEQC2 benchmark CN truth",
        }
    }
    tot_k2_units = hcc_k2["n_units"]
    tot_k2_deep = hcc_k2["determined_deep_edges"]
    tot_k3_units = hcc_k3["n_units"]
    tot_k3_deep = hcc_k3["determined_deep_edges"]
    for s, r in per_sample.items():
        if "clean" not in r:
            continue
        pooled[s] = {
            "k_ge_2_units": r["clean"]["k_ge_2"]["units"],
            "k_ge_2_determined": r["clean"]["k_ge_2"]["determined"]["count"]
            if r["clean"]["k_ge_2"]["determined"]
            else 0,
            "k_ge_2_deep_edges": r["clean"]["k_ge_2"]["deep_edges"],
            "k_ge_3_units": r["clean"]["k_ge_3"]["units"],
            "k_ge_3_determined": r["clean"]["k_ge_3"]["determined"]["count"]
            if r["clean"]["k_ge_3"]["determined"]
            else 0,
            "k_ge_3_deep_edges": r["clean"]["k_ge_3"]["deep_edges"],
            "basis": "SAVANA published fit + validated neutral gate",
        }
        tot_k2_units += r["clean"]["k_ge_2"]["units"]
        tot_k2_deep += r["clean"]["k_ge_2"]["deep_edges"]
        tot_k3_units += r["clean"]["k_ge_3"]["units"]
        tot_k3_deep += r["clean"]["k_ge_3"]["deep_edges"]

    out = {
        "generated_by": os.path.basename(__file__),
        "neutral_gate": {
            "total_cn": "abs(CN - 2) <= %.1f" % CN_TOL,
            "heterozygosity": "segment mean BAF < %.2f" % BAF_LOH,
            "strict_variant_adds": "minor-allele CN >= %.1f" % MINOR_MIN,
            "validated_on": "HCC1395 vs SEQC2: precision 90.29%, recall 7.94% after recalibration",
            "unit_rule": "every active site of the unit must fall on clean ground",
        },
        "caveats": [
            "H1437/H2009/HCC1954 use SAVANA's published purity/ploidy; they passed the internal BAF consistency audit but have no external CN truth, so the 90.29% precision is an expectation transferred from HCC1395, not a measurement on these samples",
            "at lower purity an LOH segment's BAF is pulled toward 0.5, weakening the LOH gate on HCC1954 (purity 0.66) relative to H1437/H2009 (0.95); the strict variant is reported alongside",
            "cross-sample pooling of deep edges counts independent local lineages, not one shared tree",
        ],
        "per_sample": per_sample,
        "pooled_clean_anchor_set": pooled,
        "pooled_totals": {
            "k_ge_2_units": tot_k2_units,
            "k_ge_2_deep_edges": tot_k2_deep,
            "k_ge_3_units": tot_k3_units,
            "k_ge_3_deep_edges": tot_k3_deep,
            "hcc1395_share_of_deep_edges_percent": round(
                100.0 * hcc_k2["determined_deep_edges"] / tot_k2_deep, 2
            )
            if tot_k2_deep
            else None,
        },
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    for s, r in per_sample.items():
        if "clean" not in r:
            print(f"{s}: {r.get('status')}")
            continue
        c2 = r["clean"]["k_ge_2"]
        c3 = r["clean"]["k_ge_3"]
        s2 = r["clean_strict"]["k_ge_2"]
        print(
            f"{s:9s} ranked={r['ranked_units_autosomal']:6d} states={r['state_counts']}\n"
            f"          clean k>=2: units={c2['units']:4d} det={c2['determined']['percent'] if c2['determined'] else 0:5.1f}% "
            f"deep={c2['deep_edges']:4d} chrom={c2['chromosomes']}\n"
            f"          clean k>=3: units={c3['units']:4d} deep={c3['deep_edges']:4d}\n"
            f"          strict k>=2: units={s2['units']:4d} deep={s2['deep_edges']:4d}"
        )
    print("\npooled:", json.dumps(out["pooled_totals"], indent=1))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
