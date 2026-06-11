#!/usr/bin/env python3
"""
62_fast_cnv_validation_batch.py — FAST batch (existing files, zero BAM) for the
CNV-vs-mutation methylation validation. Three analyses:
  (1) CN-state class-contrast: HP-axis |Δβ| in CN=2-neutral-nonLOH vs gain vs cnLOH
      (does clustering shrink where copy-partition is structurally absent?)
  (2) falsifiers: signed dose-direction (artifact predicts +ρ) + cnLOH consistency
  (3) per-locus copy-vs-cis decision table over all 816 scanned loci
Outputs JSON; numbers read back by the report (no fabrication).
"""
import json, gzip
import numpy as np
from scipy import stats
from collections import Counter, defaultdict
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
OUT = {}

# ---------- load master forest (HP-axis records) ----------
rows = []
with open(f"{ROOT}/genome_survey_v2/cn_confound/master_o1_cn.tsv") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hdr)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[ix["axis_type"]] != "HP":
            continue
        def g(c):
            v = p[ix[c]]
            return v
        try:
            rows.append({"abs_delta": float(p[ix["abs_delta"]]), "wp": float(p[ix["wilcoxon_p"]]),
                         "median_cn": float(p[ix["median_cn"]]) if p[ix["median_cn"]] not in ("", "NA", "nan") else None,
                         "cn_class": p[ix["cn_class"]], "loh": p[ix["loh_status"]],
                         "cnloh": p[ix["cnloh_flag"]], "is_tp": p[ix["is_tp"]]})
        except (ValueError, KeyError):
            continue
print(f"HP-axis records: {len(rows)}")
print("cn_class:", dict(Counter(r["cn_class"] for r in rows)))
print("loh_status:", dict(Counter(r["loh"] for r in rows)))
print("cnloh_flag:", dict(Counter(r["cnloh"] for r in rows)))

# ---------- (1) CN-state class-contrast ----------
def stratum(rs, pred):
    v = [r["abs_delta"] for r in rs if pred(r)]
    return v
SIG = lambda r: r["wp"] < 0.05
is_cnloh = lambda r: r["cnloh"] in ("1", "True", "true", "TRUE")
neutral_nonloh = lambda r: r["cn_class"] == "neutral" and r["loh"] == "nonLOH" and not is_cnloh(r)
gain_nonloh = lambda r: r["cn_class"] == "gain" and r["loh"] == "nonLOH"
strata = {
    "neutral_nonLOH": neutral_nonloh,
    "gain_nonLOH": gain_nonloh,
    "cnLOH": is_cnloh,
    "gain_LOH": lambda r: r["cn_class"] == "gain" and r["loh"] == "LOH",
    "loss": lambda r: r["cn_class"] == "loss",
}
c1 = {"sig": {}, "all": {}}
for name, pred in strata.items():
    sig_v = stratum([r for r in rows if SIG(r)], pred)
    all_v = stratum(rows, pred)
    c1["sig"][name] = {"n": len(sig_v), "median_abs_delta": round(float(np.median(sig_v)), 4) if sig_v else None}
    c1["all"][name] = {"n": len(all_v), "median_abs_delta": round(float(np.median(all_v)), 4) if all_v else None}
# decisive MW: neutral_nonLOH vs gain_nonLOH (sig)
sn = stratum([r for r in rows if SIG(r)], neutral_nonloh)
sg = stratum([r for r in rows if SIG(r)], gain_nonloh)
mw = stats.mannwhitneyu(sn, sg, alternative="two-sided") if sn and sg else None
c1["decisive_MW_neutral_vs_gain_sig"] = {"neutral_n": len(sn), "neutral_median": round(float(np.median(sn)), 4) if sn else None,
                                         "gain_n": len(sg), "gain_median": round(float(np.median(sg)), 4) if sg else None,
                                         "MW_p": round(float(mw.pvalue), 4) if mw else None,
                                         "verdict": "INDISTINGUISHABLE -> dosage NOT driver" if mw and mw.pvalue > 0.05 else "DIFFERENT"}
OUT["analysis1_cn_class_contrast"] = c1

# ---------- (2) falsifiers ----------
sig_rows = [r for r in rows if SIG(r) and r["median_cn"] is not None]
ad = np.array([r["abs_delta"] for r in sig_rows]); cn = np.array([r["median_cn"] for r in sig_rows])
rho = stats.spearmanr(ad, cn)
# dosage artifact predicts POSITIVE rho (more |Δβ| at higher CN). Observed:
OUT["analysis2_falsifier"] = {
    "signed_dose_rho_absdelta_vs_CN": {"rho": round(float(rho.correlation), 4), "p": float(f"{rho.pvalue:.2e}"), "n": len(sig_rows),
        "verdict": "REVERSED/NULL -> dosage-artifact model REFUTED" if rho.correlation <= 0.1 else "positive -> compatible with dosage artifact"},
    "cnLOH_consistency": {"cnLOH_sig_median": c1["sig"]["cnLOH"]["median_abs_delta"], "cnLOH_sig_n": c1["sig"]["cnLOH"]["n"],
        "neutral_sig_median": c1["sig"]["neutral_nonLOH"]["median_abs_delta"],
        "note": "cnLOH = copy-neutral (CN=2) + single parental copy; if copy-driven its clustering should differ systematically from balanced neutral"},
}

# ---------- (3) per-locus copy-vs-cis decision table (816 loci) ----------
scan = json.load(open(f"{ROOT}/genome_survey_v2/cis_scan_full.json"))
cp = {r["locus"]: r for r in json.load(open(f"{ROOT}/genome_survey_v2/copy_partition_confirm.json"))["rows"]}
# CGI lookup
CGI = defaultdict(list)
with gzip.open(f"{ROOT}/data/cpgIslandExt_hg38.txt.gz", "rt") as f:
    for line in f:
        p = line.rstrip().split("\t")
        if len(p) >= 5:
            CGI[p[1]].append((int(p[2]) + 1, int(p[3])))
def cgi_dist(c, pos):
    best = None
    for s, e in CGI.get(c, []):
        d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        best = d if best is None else min(best, d)
    return best if best is not None else 10**9
BONF = 0.05 / 816
def verdict(r):
    loc = r["locus"]; chrom, pos = loc.split(":"); pos = int(pos)
    if r.get("mechanical_cis") in ("CREATES_CpG", "DESTROYS_CpG"):
        return "mechanical-artifact"
    if r.get("cis_tier") != "T3":
        return "not-cis (T0/T2)"
    pc = r.get("p_cis")
    if pc is None or pc >= BONF:
        return "T3-nominal (fails Bonferroni)"
    # T3 + Bonferroni-clean + mechanical NEUTRAL
    cgid = cgi_dist(chrom, pos)
    in_cp = loc in cp
    if in_cp:
        dw, dh = cp[loc].get("d_within"), cp[loc].get("d_HP")
        if dw is not None and dh and abs(dw) < 0.5 * abs(dh):
            return "copy-artifact (copy-test: copy-dominated)"
        return "cis-candidate (copy-test: allele-dominated)"
    # not copy-testable (pure-ALT tag) -> use CGI as weak prior
    return "cis-candidate-UNTESTABLE (pure-ALT tag; CGI-desert)" if cgid > 4000 else "cis-candidate-UNTESTABLE (pure-ALT tag; near-CGI)"
vt = Counter()
detail = []
for r in scan:
    v = verdict(r)
    vt[v] += 1
    if r.get("cis_tier") == "T3" and (r.get("p_cis") or 1) < BONF:
        chrom, pos = r["locus"].split(":")
        detail.append({"locus": r["locus"], "verdict": v, "total_cn": r.get("total_cn"), "loh": r.get("loh"),
                       "mechanical_cis": r.get("mechanical_cis"), "d_cis": r.get("d_cis"), "p_cis": r.get("p_cis"),
                       "cgi_dist": cgi_dist(chrom, int(pos))})
OUT["analysis3_decision_table"] = {"verdict_counts": dict(vt), "bonferroni_T3_detail": sorted(detail, key=lambda x: x["cgi_dist"])}

json.dump(OUT, open(f"{ROOT}/genome_survey_v2/fast_cnv_validation.json", "w"), indent=1)
print("\n=== written fast_cnv_validation.json ===")
print(json.dumps({"A1_decisive": c1["decisive_MW_neutral_vs_gain_sig"],
                  "A2_signed_rho": OUT["analysis2_falsifier"]["signed_dose_rho_absdelta_vs_CN"],
                  "A3_verdicts": dict(vt)}, ensure_ascii=False, indent=1))
