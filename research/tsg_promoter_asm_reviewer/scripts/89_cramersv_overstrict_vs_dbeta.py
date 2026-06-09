#!/usr/bin/env python3
"""
89 - Three diagnostics on the CramersV gate vs Delta-beta, over the full HCC1395
TP/FP existence-scan summary (frozen baseline; gated values).

Q1  Over-strict?  How many TP have CLEAN HP-aligned cluster structure (PERMANOVA,
    distance-based, robust to sparsity; dispersion_warning excluded) yet gated
    CramersV<0.1 -> missed by the Significant gate. Root cause = min(HP1N,HP2N)
    (Cochran min-expected<5 makes chi2/CramersV unreliable -> zeroed).

Q2  Detection recipe: rank "real structure" by PERMANOVA + min-HP-N robustness,
    not by CramersV. Dump the missed-but-real TP list.

Q3  Why Delta-beta flags more/clearer structure: compare detection sets
    (Dbeta>=0.2 vs CramersV>=0.1 vs PERMANOVA), and the SPECIFICITY trade-off
    (TP:FP of each set vs the 47:1 base; Dbeta-only is FP-leaning).

Output: display_v2/diag89.json  +  missed_real_TP.tsv
"""
import csv, json

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")


def num(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def load(cls):
    return list(csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")))


def cvmax_gated(r):
    return max([num(r, k) or 0 for k in ("CramersV", "CramersV_HPFamily", "CramersV_HPFine")])


def minhpn(r):
    return min(num(r, "HP1FamilyN") or 0, num(r, "HP2FamilyN") or 0)


def clean_cluster(r):
    """Real cluster structure: PERMANOVA valid+sig, no dispersion artifact."""
    return (tb(r, "ClusterPermanovaValid") and (num(r, "ClusterPermanovaP") or 1) <= 0.05
            and not tb(r, "ClusterDispersionWarn"))


def hp_aligned(r):
    """Clusters align with haplotype labels (distance-based), no dispersion artifact."""
    return (tb(r, "LabelHPPermanovaValid") and (num(r, "LabelHPPermanovaP") or 1) <= 0.05
            and not tb(r, "LabelHPDispersionWarn"))


def main():
    tp, fp = load("tp"), load("fp")
    n_tp, n_fp = len(tp), len(fp)
    base = n_tp / n_fp
    out = {"n_tp": n_tp, "n_fp": n_fp, "base_ratio": round(base, 1)}

    # ---- Q1: over-strict (clean HP-aligned structure but CramersV gated to <0.1) ----
    def overstrict(rows):
        return [r for r in rows if (num(r, "NumReads") or 0) >= 20 and tb(r, "PassedGating")
                and clean_cluster(r) and hp_aligned(r) and cvmax_gated(r) < 0.1]
    os_tp, os_fp = overstrict(tp), overstrict(fp)
    sig_tp = [r for r in tp if tb(r, "Significant")]
    out["q1"] = dict(
        overstrict_tp=len(os_tp), overstrict_fp=len(os_fp),
        pct_of_tp=round(100 * len(os_tp) / n_tp, 2),
        sig_tp=len(sig_tp),
        ratio_overstrict_to_sig=round(len(os_tp) / max(len(sig_tp), 1), 2),
        # root cause: min HP group size distribution among over-strict TP
        cochran_fail_minHPlt5=sum(1 for r in os_tp if minhpn(r) < 5),
        minHP_5to9=sum(1 for r in os_tp if 5 <= minhpn(r) < 10),
        minHP_ge10=sum(1 for r in os_tp if minhpn(r) >= 10),
        loh_pct=round(100 * sum(1 for r in os_tp if tb(r, "Potential_LOH")) / max(len(os_tp), 1), 1),
        overstrict_tp_fp_ratio=round(len(os_tp) / max(len(os_fp), 1), 1),
    )

    # ---- Q2: dump missed-but-real TP (ranked by HP-PERMANOVA F then min-HP-N) ----
    os_tp_sorted = sorted(os_tp, key=lambda r: (-(num(r, "LabelHPPermanovaF") or 0), -minhpn(r)))
    with open(f"{DV}/missed_real_TP.tsv", "w") as f:
        f.write("Chr\tPos\tNumReads\tminHPN\tHP1N\tHP2N\tgatedCVmax\tDbeta\t"
                "ClusterPermF\tHPpermF\tLOH\toptimal_kproxy\n")
        for r in os_tp_sorted:
            f.write(f"{r['Chr']}\t{r['Pos']}\t{int(num(r,'NumReads') or 0)}\t{minhpn(r)}\t"
                    f"{int(num(r,'HP1FamilyN') or 0)}\t{int(num(r,'HP2FamilyN') or 0)}\t"
                    f"{cvmax_gated(r):.3f}\t{num(r,'HPMergedDelta') or 0:+.3f}\t"
                    f"{num(r,'ClusterPermanovaF') or 0:.1f}\t{num(r,'LabelHPPermanovaF') or 0:.1f}\t"
                    f"{tb(r,'Potential_LOH')}\t{r.get('HPFine_NGroups_CF','')}\n")

    # ---- Q3: detection-set comparison (Dbeta vs CramersV vs PERMANOVA) + specificity ----
    def sets(rows):
        db = [r for r in rows if abs(num(r, "HPMergedDelta") or 0) >= 0.2 and (num(r, "NumReads") or 0) >= 20]
        cv = [r for r in rows if cvmax_gated(r) >= 0.1 and (num(r, "NumReads") or 0) >= 20]
        pm = [r for r in rows if clean_cluster(r) and hp_aligned(r) and (num(r, "NumReads") or 0) >= 20]
        return db, cv, pm
    db_tp, cv_tp, pm_tp = sets(tp)
    db_fp, cv_fp, pm_fp = sets(fp)
    # overlaps within TP
    key = lambda r: (r["Chr"], r["Pos"])
    sdb, scv = set(map(key, db_tp)), set(map(key, cv_tp))
    out["q3"] = dict(
        tp_dbeta=len(db_tp), tp_cramersv=len(cv_tp), tp_permanova=len(pm_tp),
        dbeta_only_tp=len(sdb - scv), cramersv_only_tp=len(scv - sdb), both_tp=len(sdb & scv),
        # specificity: TP:FP of each detection set vs base 47:1
        dbeta_ratio=round(len(db_tp) / max(len(db_fp), 1), 1),
        cramersv_ratio=round(len(cv_tp) / max(len(cv_fp), 1), 1),
        permanova_ratio=round(len(pm_tp) / max(len(pm_fp), 1), 1),
        dbeta_enrich=round((len(db_tp) / max(len(db_fp), 1)) / base, 2),
        cramersv_enrich=round((len(cv_tp) / max(len(cv_fp), 1)) / base, 2),
        # among Dbeta-flagged TP: how many have NO HP-aligned clean structure (mean-shift only)?
        dbeta_no_hp_structure_tp=sum(1 for r in db_tp if not (clean_cluster(r) and hp_aligned(r))),
        dbeta_no_hp_structure_fp=sum(1 for r in db_fp if not (clean_cluster(r) and hp_aligned(r))),
    )

    with open(f"{DV}/diag89.json", "w") as f:
        json.dump(out, f, indent=2)

    # print summary
    q1, q3 = out["q1"], out["q3"]
    print(f"base TP:FP = {base:.1f}:1   (TP={n_tp} FP={n_fp})\n")
    print("== Q1 over-strict (clean HP-aligned PERMANOVA structure but gated CramersV<0.1) ==")
    print(f"   over-strict TP = {q1['overstrict_tp']} ({q1['pct_of_tp']}% of TP); FP = {q1['overstrict_fp']}"
          f"  -> set TP:FP = {q1['overstrict_tp_fp_ratio']}:1")
    print(f"   vs Significant-PASS TP = {q1['sig_tp']}  (over-strict adds {q1['ratio_overstrict_to_sig']}× more)")
    print(f"   root cause minHPN: <5(Cochran fail)={q1['cochran_fail_minHPlt5']}  "
          f"5-9={q1['minHP_5to9']}  >=10={q1['minHP_ge10']}   LOH={q1['loh_pct']}%")
    print("\n== Q3 detection sets (TP) ==")
    print(f"   Dbeta>=0.2: {q3['tp_dbeta']}   CramersV>=0.1: {q3['tp_cramersv']}   PERMANOVA: {q3['tp_permanova']}")
    print(f"   Dbeta-only={q3['dbeta_only_tp']}  CramersV-only={q3['cramersv_only_tp']}  both={q3['both_tp']}")
    print(f"   specificity TP:FP  Dbeta={q3['dbeta_ratio']}:1 ({q3['dbeta_enrich']}×base)  "
          f"CramersV={q3['cramersv_ratio']}:1 ({q3['cramersv_enrich']}×base)")
    print(f"   Dbeta-flagged TP with NO clean HP structure (mean-shift only) = {q3['dbeta_no_hp_structure_tp']}"
          f"  (FP counterpart {q3['dbeta_no_hp_structure_fp']})")
    print(f"\n[89] wrote {DV}/diag89.json + missed_real_TP.tsv ({len(os_tp_sorted)} rows)")


if __name__ == "__main__":
    main()
