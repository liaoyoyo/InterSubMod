#!/usr/bin/env python3
"""Join the three per-locus layers into one sSNV-lineage x methylation annotation.

Inputs (all produced/archived on 2026-07-26):
  methyl_locus_scan.tsv  9,414 locus x germline-HP units, within-HP statistics
  sig_scan.tsv           30,476 loci, the pipeline's own methylation clustering
  methyl_class_x_linkage_annotation.tsv  30,077 loci, sSNV linkage + CN

Every methylation flag is calibrated against the matched unimodal null carried
in the scan (sep_null), never against an absolute eyeballed cutoff.
"""
import csv, json, bisect, collections, statistics as st

S = "/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/efb6e3d8-c4af-43d8-ac97-9dffdbec60ed/scratchpad"
LNK = "docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/methyl_class_x_linkage_annotation.tsv"

units = list(csv.DictReader(open(S + "/methyl_locus_scan.tsv"), delimiter="\t"))
sig = {(r["chrom"], int(r["pos"])): r for r in csv.DictReader(open(S + "/sig_scan.tsv"), delimiter="\t")}
lnk = {}
bych = collections.defaultdict(list)
for r in csv.DictReader(open(LNK), delimiter="\t"):
    k = (r["chrom"], int(r["snv_pos"]))
    lnk[k] = r
    bych[r["chrom"]].append(int(r["snv_pos"]))
for c in bych:
    bych[c].sort()


def dens(c, p, w=50000):
    a = bych[c]
    return bisect.bisect_right(a, p + w) - bisect.bisect_left(a, p - w) - 1


# --- null calibration: 95th percentile of the matched unimodal null
nulls = sorted(float(u["sep_null"]) for u in units)
Q95 = nulls[int(0.95 * len(nulls))]

rows = []
for u in units:
    c, p = u["chrom"], int(u["pos"])
    sg, lk = sig.get((c, p)), lnk.get((c, p))
    if not lk:
        continue
    sep_a, sep_r = float(u["sep_alt"]), float(u["sep_ref"])
    dsp = abs(float(u["dsplit_alt"]))
    dar, par = float(u["d_altref"]), float(u["p_altref"])
    exc = float(u["exc_altref"])

    alt_multi = sep_a >= Q95 and dsp >= 0.20
    ref_multi = sep_r >= Q95
    shift = par < 0.05 and abs(dar) >= 0.10 and exc > 0
    strand_flag = float(u["strand_p"]) < 0.05

    if alt_multi and not ref_multi:
        mstate = "M2_ALT_ONLY_MULTI"
    elif alt_multi and ref_multi:
        mstate = "M3_REGION_MULTI"
    elif shift:
        mstate = "M1_ALTREF_SHIFT"
    else:
        mstate = "M0_NONE"

    npart = int(lk["n_linkage_partners"])
    maxco = int(lk["max_coread"] or 0)
    bs = lk["both_somatic"] == "True"
    if npart >= 1 and bs and maxco >= 6:
        lstate = "L1_BACKBONE"
    elif npart >= 1:
        lstate = "L2_WEAK_LINK"
    else:
        lstate = "L3_ISOLATED"

    rows.append({
        "chrom": c, "pos": p, "hp": u["hp"],
        "lineage_state": lstate, "n_partners": npart, "max_coread": maxco,
        "top_rel": lk["top_rel"], "both_somatic": lk["both_somatic"],
        "local_density_50kb": dens(c, p),
        "cn_state": lk["cn_state"], "is_loh": lk["is_loh"],
        "methyl_state": mstate,
        "flag_alt_multi": int(alt_multi), "flag_ref_multi": int(ref_multi),
        "flag_altref_shift": int(shift), "flag_strand_assoc": int(strand_flag),
        "n_alt": int(u["n_alt"]), "n_ref": int(u["n_ref"]), "n_common_cpg": int(u["n_common"]),
        "mean_alt": float(u["mean_alt"]), "mean_ref": float(u["mean_ref"]),
        "d_altref": dar, "p_altref": par, "excess_vs_null": exc,
        "sep_alt": sep_a, "dsplit_alt": float(u["dsplit_alt"]), "sep_ref": sep_r,
        "sep_null": float(u["sep_null"]),
        "pipeline_optimal_k": int(sg["optimal_k"]) if sg else 0,
        "pipeline_passed_gating": int(sg["passed_gating"]) if sg else 0,
        "v_alt": float(sg["v_alt"]) if sg else 0.0,
        "v_hp": float(sg["v_hp"]) if sg else 0.0,
        "cpg_span": "{0}-{1}".format(u["cpg_lo"], u["cpg_hi"]),
    })

cols = list(rows[0].keys())
OUT = "docs/reports/in_progress/2026/07/20260726_ssnv_lineage_x_methyl_annotation/ssnv_lineage_x_methyl_annotation.tsv"
import os
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")

# ------------------------------------------------------------------ summary
LS = ["L1_BACKBONE", "L2_WEAK_LINK", "L3_ISOLATED"]
MS = ["M2_ALT_ONLY_MULTI", "M3_REGION_MULTI", "M1_ALTREF_SHIFT", "M0_NONE"]
tab = {l: {m: 0 for m in MS} for l in LS}
for r in rows:
    tab[r["lineage_state"]][r["methyl_state"]] += 1

strat = []
DB = [(0, 0, "0"), (1, 1, "1"), (2, 3, "2-3"), (4, 8, "4-8"), (9, 10 ** 9, ">=9")]
for lo, hi, lab in DB:
    e = {"bin": lab}
    for l in LS:
        v = [r for r in rows if r["lineage_state"] == l and lo <= r["local_density_50kb"] <= hi]
        e[l] = {"n": len(v),
                "pct_alt_only": (sum(1 for r in v if r["methyl_state"] == "M2_ALT_ONLY_MULTI") / len(v) * 100) if v else None,
                "pct_shift": (sum(1 for r in v if r["flag_altref_shift"]) / len(v) * 100) if v else None}
    strat.append(e)

data = {
    "meta": {"date": "2026-07-26", "sample": "HCC1395",
             "archive": "20260119_all-with-w5000_1", "null_q95": round(Q95, 3),
             "n_units": len(rows), "n_loci": len(set((r["chrom"], r["pos"]) for r in rows)),
             "tsv": os.path.basename(OUT)},
    "inputs": {"scan_units": len(units), "sig_loci": len(sig), "linkage_loci": len(lnk)},
    "table": tab,
    "totals_lineage": {l: sum(tab[l].values()) for l in LS},
    "totals_methyl": {m: sum(tab[l][m] for l in LS) for m in MS},
    "strat": strat,
    "pipeline_k": collections.Counter(r["pipeline_optimal_k"] for r in rows),
    "gating": sum(r["pipeline_passed_gating"] for r in rows),
    "strand_flag": sum(r["flag_strand_assoc"] for r in rows),
    "ref_multi": sum(r["flag_ref_multi"] for r in rows),
    "alt_multi": sum(r["flag_alt_multi"] for r in rows),
    "cn": {cn: {"n": sum(1 for r in rows if r["cn_state"] == cn),
                "pct_alt_only": sum(1 for r in rows if r["cn_state"] == cn and r["methyl_state"] == "M2_ALT_ONLY_MULTI")
                / max(sum(1 for r in rows if r["cn_state"] == cn), 1) * 100}
           for cn in sorted(set(r["cn_state"] for r in rows))},
    "examples": sorted([r for r in rows if r["methyl_state"] == "M2_ALT_ONLY_MULTI"
                        and r["lineage_state"] == "L1_BACKBONE"],
                       key=lambda r: -abs(r["dsplit_alt"]))[:12],
}
data["pipeline_k"] = {str(k): v for k, v in sorted(data["pipeline_k"].items())}
json.dump(data, open(S + "/annotation_summary.json", "w"), ensure_ascii=False, indent=1)

print("null-95th (matched unimodal) =", round(Q95, 3))
print("wrote", OUT, len(rows), "units /", data["meta"]["n_loci"], "loci\n")
print("═══ sSNV-lineage 狀態 × 甲基標記狀態（單元數）═══")
hdr = "{0:<16}".format("") + "".join("{0:>20}".format(m) for m in MS) + "{0:>9}".format("合計")
print(hdr)
for l in LS:
    tot = sum(tab[l].values())
    line = "{0:<16}".format(l)
    for m in MS:
        n = tab[l][m]
        line += "{0:>13} {1:>5.1f}%".format(n, n / max(tot, 1) * 100)
    print(line + "{0:>9}".format(tot))
tot = len(rows)
line = "{0:<16}".format("合計")
for m in MS:
    n = data["totals_methyl"][m]
    line += "{0:>13} {1:>5.1f}%".format(n, n / tot * 100)
print(line + "{0:>9}".format(tot))
print("\n═══ 密度分層：M2（僅 ALT 多群）比例 ═══")
print("{0:<10}".format("密度") + "".join("{0:>22}".format(l) for l in LS))
for e in strat:
    line = "{0:<10}".format(e["bin"])
    for l in LS:
        v = e[l]
        line += "{0:>22}".format("n={0} {1}".format(v["n"], "-" if v["pct_alt_only"] is None else "{0:.1f}%".format(v["pct_alt_only"])))
    print(line)
print("\nREF 內多群（內建對照）= {0} ({1:.1f}%)   ALT 內多群 = {2} ({3:.1f}%)".format(
    data["ref_multi"], data["ref_multi"] / tot * 100, data["alt_multi"], data["alt_multi"] / tot * 100))
print("strand 關聯旗標 = {0} ({1:.1f}%)".format(data["strand_flag"], data["strand_flag"] / tot * 100))
print("管線 optimal_k 分布:", data["pipeline_k"], " passed_gating =", data["gating"])
