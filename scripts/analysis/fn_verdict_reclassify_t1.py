#!/usr/bin/env python3
"""
fn_verdict_reclassify_t1.py — T1 of the ISM verdict false-negative audit.

Offline A/B: recompute VerificationClass with the label-side PERMANOVA wired into
label_sig (and, for an upgrade test, cluster PERMANOVA wired into cluster_sig),
count REAL flips out of Noise/Weak, confirm orig_Strong is never demoted, and
ISOLATE the subsets that may NOT be "the structure we want":
  - sample-axis-only (FN10 tumor-vs-normal confound) — EXCLUDED from the real rescue
  - dispersion-warned label-PERMANOVA (possible betadisper artifact)
  - germline-HP-driven (labelHP) vs allele-driven (labelAllele)

Derivation of the ORIGINAL booleans from VerificationClass (matches
SignificanceAnalyzer.cpp:326-339 exactly):
  Strong   -> label_sig=T, cluster_sig=T
  Subclone -> label_sig=F, cluster_sig=T
  Weak     -> label_sig=T, cluster_sig=F
  Noise    -> label_sig=F, cluster_sig=F

NO C++ rerun, NO fabrication — every number is a deterministic count over the CSV.

Usage: python3 fn_verdict_reclassify_t1.py OUT.json LABEL=CSV [LABEL=CSV ...]
"""
import csv
import json
import math
import sys


def fbool(s):
    return str(s).strip().lower() in ("true", "1", "yes", "t")


def ffloat(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return math.nan


def le(field, row, alpha=0.05):
    v = ffloat(row.get(field, ""))
    return (not math.isnan(v)) and v <= alpha


def sigperm(p_field, valid_field, row):
    return fbool(row.get(valid_field, "")) and le(p_field, row)


def cls(label_sig, cluster_sig):
    if label_sig and cluster_sig:
        return "Strong"
    if (not label_sig) and cluster_sig:
        return "Subclone"
    if label_sig and (not cluster_sig):
        return "Weak"
    return "Noise"


def audit(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    n = len(rows)

    # transition matrices
    trans_label = {}   # T1(i): label_sig' = label_sig OR label-PERMANOVA
    trans_full = {}    # T1(ii): also cluster_sig' = cluster_sig OR cluster-PERMANOVA
    orig_strong_demoted_label = 0
    orig_strong_demoted_full = 0

    # rescue breakdowns (out of Noise specifically)
    noise_total = 0
    noise_rescued_label = 0          # Noise -> not Noise under T1(i)
    noise_rescued_dispclean = 0      # of which, label-PERMANOVA not dispersion-warned
    noise_by_labelhp = 0
    noise_by_labelallele = 0
    noise_sample_only_excluded = 0   # Noise sig ONLY on sample axis (FN10) — NOT counted as real
    noise_rescued_full_upgrade = 0   # Noise -> Subclone/Strong via cluster PERMANOVA

    weak_total = 0
    weak_upgraded_full = 0           # Weak -> Strong via cluster PERMANOVA

    for r in rows:
        vc = r.get("VerificationClass", "").strip()
        orig_label = vc in ("Strong", "Weak")
        orig_cluster = vc in ("Strong", "Subclone")

        plhp = sigperm("LabelHPPermanovaP", "LabelHPPermanovaValid", r)
        plal = sigperm("LabelAllelePermanovaP", "LabelAllelePermanovaValid", r)
        pclu = sigperm("ClusterPermanovaP", "ClusterPermanovaValid", r)
        label_perm = plhp or plal
        disp_warn_label = fbool(r.get("LabelHPDispersionWarn", "")) or fbool(r.get("LabelAlleleDispersionWarn", ""))
        sample_sig = le("SampleASM_P", r)

        # T1(i): label rescue only (exclude sample confound axis)
        label_sig_new = orig_label or label_perm
        nc_label = cls(label_sig_new, orig_cluster)
        key = f"{vc}->{nc_label}"
        trans_label[key] = trans_label.get(key, 0) + 1
        if vc == "Strong" and nc_label != "Strong":
            orig_strong_demoted_label += 1

        # T1(ii): also upgrade cluster axis via cluster PERMANOVA
        cluster_sig_new = orig_cluster or pclu
        nc_full = cls(label_sig_new, cluster_sig_new)
        key2 = f"{vc}->{nc_full}"
        trans_full[key2] = trans_full.get(key2, 0) + 1
        if vc == "Strong" and nc_full != "Strong":
            orig_strong_demoted_full += 1

        if vc == "Noise":
            noise_total += 1
            if label_perm:
                noise_rescued_label += 1
                if not disp_warn_label:
                    noise_rescued_dispclean += 1
                if plhp:
                    noise_by_labelhp += 1
                if plal:
                    noise_by_labelallele += 1
            elif sample_sig:
                noise_sample_only_excluded += 1
            if nc_full in ("Subclone", "Strong"):
                noise_rescued_full_upgrade += 1
        if vc == "Weak":
            weak_total += 1
            if nc_full == "Strong":
                weak_upgraded_full += 1

    return {
        "n_rows": n,
        "T1i_label_rescue": {
            "transitions": dict(sorted(trans_label.items(), key=lambda kv: -kv[1])),
            "orig_Strong_demoted": orig_strong_demoted_label,
            "noise_total": noise_total,
            "noise_rescued_to_Weak_by_labelPERMANOVA": noise_rescued_label,
            "of_which_dispersion_clean": noise_rescued_dispclean,
            "noise_rescue_by_labelHP": noise_by_labelhp,
            "noise_rescue_by_labelAllele": noise_by_labelallele,
            "noise_sample_axis_ONLY_excluded_FN10": noise_sample_only_excluded,
            "noise_rescue_frac": round(noise_rescued_label / noise_total, 4) if noise_total else None,
            "noise_rescue_frac_dispclean": round(noise_rescued_dispclean / noise_total, 4) if noise_total else None,
        },
        "T1ii_full_upgrade_cluster_too": {
            "transitions": dict(sorted(trans_full.items(), key=lambda kv: -kv[1])),
            "orig_Strong_demoted": orig_strong_demoted_full,
            "noise_to_Subclone_or_Strong": noise_rescued_full_upgrade,
            "weak_total": weak_total,
            "weak_to_Strong": weak_upgraded_full,
        },
    }


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: fn_verdict_reclassify_t1.py OUT.json LABEL=CSV [...]")
    out = sys.argv[1]
    res = {}
    for spec in sys.argv[2:]:
        label, _, path = spec.partition("=")
        res[label] = {"source_csv": path, **audit(path)}
    with open(out, "w") as fh:
        json.dump(res, fh, indent=2)
    for label, d in res.items():
        a = d["T1i_label_rescue"]
        b = d["T1ii_full_upgrade_cluster_too"]
        print(f"\n===== {label} ({d['n_rows']} rows) =====")
        print(f"  [T1i label rescue] Noise {a['noise_total']}: rescued->Weak by label-PERMANOVA = "
              f"{a['noise_rescued_to_Weak_by_labelPERMANOVA']} ({a['noise_rescue_frac']}), "
              f"dispersion-clean = {a['of_which_dispersion_clean']} ({a['noise_rescue_frac_dispclean']})")
        print(f"     by arm: labelHP={a['noise_rescue_by_labelHP']} labelAllele={a['noise_rescue_by_labelAllele']} "
              f"| sample-ONLY excluded (FN10) = {a['noise_sample_axis_ONLY_excluded_FN10']}")
        print(f"     orig_Strong_demoted = {a['orig_Strong_demoted']}")
        print(f"     transitions: {a['transitions']}")
        print(f"  [T1ii +cluster upgrade] Noise->Subclone/Strong = {b['noise_to_Subclone_or_Strong']} ; "
              f"Weak->Strong = {b['weak_to_Strong']}/{b['weak_total']} ; orig_Strong_demoted = {b['orig_Strong_demoted']}")
    print(f"\n[written] {out}")


if __name__ == "__main__":
    main()
