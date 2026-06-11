#!/usr/bin/env python3
"""
73 - Aggregate ISM significance_summary across 6 samples × {TP,FP,FN} (existence)
+ HCC1395 cis (normal-anchored). Deterministic; §13 layer-A.

EXISTENCE (ism_existence_scan/<sample>_<cls>/significance_summary.csv):
  per (sample,class): n, n_with_somatic_subhap, n_significant (ISM gate),
  significant_rate, n_hpmerged_sig, median CramersV_HPFine, median |HPMergedDelta|.
  -> TP vs FP vs FN ASM-existence comparison (powered, replaces the small-N script 69).

CIS (ism_cis_scan/HCC1395_<cls>_cis/significance_summary.csv):
  per class: HP_Residual_Delta distribution, n significant residual (cis candidates),
  n SampleASM_Sig. -> HCC1395 cis-ASM positions.

Output: research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/ism_aggregate.json
"""
import os, csv, json
import numpy as np

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
CIS = "/big7_disk/liaoyoyo2001/ism_cis_scan"
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
       "genome_survey_v2/cn_confound/cross_sample/ism_aggregate.json")

ALL = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
CANCER = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
          "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}
CLS = ["tp", "fp", "fn"]


def fnum(v):
    try:
        f = float(v)
        return f if not np.isnan(f) else None
    except (ValueError, TypeError):
        return None


def load(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return list(csv.DictReader(f))


def med(vals):
    a = np.array([v for v in vals if v is not None])
    return round(float(np.median(a)), 4) if a.size else None


def summarize_existence(rows):
    n = len(rows)
    sig = sum(1 for r in rows if str(r.get("Significant", "")).lower() == "true")
    hpm_sig = sum(1 for r in rows if str(r.get("HPMergedSig", "")).lower() == "true")
    # somatic subhap present = HPFineN_HP1S or HP2S >= 3
    subhap = 0
    for r in rows:
        n1s = fnum(r.get("HPFineN_HP1S")) or 0
        n2s = fnum(r.get("HPFineN_HP2S")) or 0
        if n1s >= 3 or n2s >= 3:
            subhap += 1
    # subhap-conditional significance (matched comparison: only loci WITH somatic subhap)
    sig_in_subhap = 0
    for r in rows:
        n1s = fnum(r.get("HPFineN_HP1S")) or 0
        n2s = fnum(r.get("HPFineN_HP2S")) or 0
        if (n1s >= 3 or n2s >= 3) and str(r.get("Significant", "")).lower() == "true":
            sig_in_subhap += 1
    cramers = [fnum(r.get("CramersV_HPFine")) for r in rows]
    hpmd = [fnum(r.get("HPMergedDelta")) for r in rows]
    hpmd_abs = [abs(x) for x in hpmd if x is not None]
    nreads = [fnum(r.get("NumReads")) for r in rows]
    return dict(
        n=n, n_with_subhap=subhap,
        n_significant=sig, significant_rate=round(sig / n, 4) if n else None,
        n_sig_in_subhap=sig_in_subhap,
        sig_rate_in_subhap=round(sig_in_subhap / subhap, 4) if subhap else None,
        n_hpmerged_sig=hpm_sig, hpmerged_sig_rate=round(hpm_sig / n, 4) if n else None,
        median_cramersV_HPFine=med(cramers),
        median_abs_HPMergedDelta=round(float(np.median(hpmd_abs)), 4) if hpmd_abs else None,
        median_NumReads=med(nreads),
    )


def summarize_cis(rows):
    n = len(rows)
    res_d = [fnum(r.get("HP_Residual_Delta")) for r in rows]
    res_p = [fnum(r.get("HP_Residual_P")) for r in rows]
    # cis candidate = residual significant (p<0.05) AND |residual delta| meaningful
    cis_cand = 0
    for r in rows:
        p = fnum(r.get("HP_Residual_P"))
        d = fnum(r.get("HP_Residual_Delta"))
        if p is not None and p < 0.05 and d is not None and abs(d) >= 0.10:
            cis_cand += 1
    sample_asm_sig = sum(1 for r in rows if str(r.get("SampleASM_Sig", "")).lower() == "true")
    rd_abs = [abs(x) for x in res_d if x is not None]
    return dict(
        n=n, n_residual_computed=sum(1 for x in res_d if x is not None),
        n_cis_candidate=cis_cand, cis_candidate_rate=round(cis_cand / n, 4) if n else None,
        median_abs_HP_Residual_Delta=round(float(np.median(rd_abs)), 4) if rd_abs else None,
        n_sampleASM_sig=sample_asm_sig,
    )


def main():
    existence = {}
    for s in ALL:
        existence[s] = {}
        for c in CLS:
            rows = load(f"{EX}/{s}_{c}/significance_summary.csv")
            existence[s][c] = summarize_existence(rows) if rows else None

    cis = {}
    for c in CLS:
        rows = load(f"{CIS}/HCC1395_{c}_cis/significance_summary.csv")
        cis[c] = summarize_cis(rows) if rows else None

    # TP vs FP discrimination (powered): significant_rate per class, pooled + per sample
    disc = {}
    for c in CLS:
        tot_sig = sum(existence[s][c]["n_significant"] for s in ALL if existence[s][c])
        tot_n = sum(existence[s][c]["n"] for s in ALL if existence[s][c])
        # subhap-conditional pooled (matched: only loci with somatic subhap)
        tot_subhap = sum(existence[s][c]["n_with_subhap"] for s in ALL if existence[s][c])
        tot_sig_subhap = sum(existence[s][c]["n_sig_in_subhap"] for s in ALL if existence[s][c])
        disc[c] = dict(pooled_significant_rate=round(tot_sig / tot_n, 4) if tot_n else None,
                       pooled_n=tot_n, pooled_n_sig=tot_sig,
                       pooled_n_with_subhap=tot_subhap,
                       pooled_sig_rate_in_subhap=round(tot_sig_subhap / tot_subhap, 4) if tot_subhap else None,
                       pooled_n_sig_in_subhap=tot_sig_subhap,
                       per_sample_rates={s: existence[s][c]["significant_rate"] for s in ALL
                                         if existence[s][c]})

    out = dict(
        meta=dict(script="73_aggregate_ism.py", samples=ALL, classes=CLS,
                  tool="unmodified ISM (build/bin/inter_sub_mod) significance_summary.csv",
                  note=("Existence = ISM native Significant gate (PassedGating AND global_p<=0.05 "
                        "AND CramersV>=0.1 AND NumReads>=20). Cis = HCC1395 normal-anchored "
                        "HP_Residual_Delta = tumor_HP_delta - normal_HP_delta (significant + "
                        "|delta|>=0.10 = cis candidate). TP/FP/FN = caller-vs-truth classes.")),
        existence=existence,
        discrimination=disc,
        cis_hcc1395=cis,
    )
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2,
                  default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
    print(f"[73] wrote {OUT}")
    print("[73] EXISTENCE significant-rate per (sample,class):")
    print(f"     {'sample':9s} {'TP':>16s} {'FP':>16s} {'FN':>16s}")
    for s in ALL:
        cells = []
        for c in CLS:
            e = existence[s][c]
            cells.append(f"{e['n_significant']}/{e['n']}({e['significant_rate']})" if e else "—")
        print(f"     {s:9s} {cells[0]:>16s} {cells[1]:>16s} {cells[2]:>16s}")
    print("[73] DISCRIMINATION (pooled significant-rate): " +
          " ".join(f"{c.upper()}={disc[c]['pooled_significant_rate']}({disc[c]['pooled_n_sig']}/{disc[c]['pooled_n']})" for c in CLS))
    print("[73] HCC1395 CIS:")
    for c in CLS:
        v = cis[c]
        if v:
            print(f"     {c.upper()}: n={v['n']} cis_candidate={v['n_cis_candidate']}({v['cis_candidate_rate']}) "
                  f"sampleASM_sig={v['n_sampleASM_sig']} med|residual|={v['median_abs_HP_Residual_Delta']}")


if __name__ == "__main__":
    main()
