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

The ORIGINAL booleans are read directly from schema-v2 LabelFirstSupport and
ClusterFirstSupport.  A historical four-state file is accepted only with the
explicit --allow-unversioned-v1 flag, then mapped as follows:
  Strong   -> label_sig=T, cluster_sig=T
  Subclone -> label_sig=F, cluster_sig=T
  Weak     -> label_sig=T, cluster_sig=F
  Noise    -> label_sig=F, cluster_sig=F

NO C++ rerun, NO fabrication — every number is a deterministic count over the CSV.

Usage: python3 fn_verdict_reclassify_t1.py [--allow-unversioned-v1] OUT.json LABEL=CSV [LABEL=CSV ...]
"""
import argparse
import json
import math
import sys
import warnings
from pathlib import Path

import pandas as pd

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
if str(WORKSPACE_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKSPACE_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    SchemaContractError,
    read_evidence,
    select_current_view,
    select_legacy_view,
)


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


def load_original_support(path, allow_unversioned_v1=False):
    df = pd.read_csv(path)
    if "VerificationSchemaVersion" in df.columns:
        current = select_current_view(df)
        if current.unknown_counts:
            raise SchemaContractError(f"FN verdict input: unknown current classes: {current.unknown_counts}")
        evidence = read_evidence(df)
        label_values = evidence["LabelFirstSupport"]
        cluster_values = evidence["ClusterFirstSupport"]
        metadata = {
            "schema_status": "V2_EVIDENCE",
            "selection_field": "LabelFirstSupport+ClusterFirstSupport",
            "warnings": [],
            "unknown_counts": {},
        }
    else:
        if not allow_unversioned_v1:
            raise SchemaContractError(
                "FN verdict input is unversioned; --allow-unversioned-v1 is required"
            )
        legacy = select_legacy_view(
            df,
            allow_unversioned_v1=True,
            unversioned_unknown_policy="fail",
        )
        label_values = legacy.values.isin(["Strong", "Weak"])
        cluster_values = legacy.values.isin(["Strong", "Subclone"])
        metadata = legacy.metadata()
        if legacy.field == "VerificationClass_Legacy":
            message = (
                "UNVERSIONED: using explicit VerificationClass_Legacy under "
                "--allow-unversioned-v1 authorization"
            )
            warnings.warn(message, UserWarning, stacklevel=2)
            metadata = {
                **metadata,
                "schema_status": "UNVERSIONED_V1_EXPLICIT_LEGACY",
                "warnings": [message],
            }
    result = df.copy()
    result["_OriginalLabelFirstSupport"] = label_values.astype(bool)
    result["_OriginalClusterFirstSupport"] = cluster_values.astype(bool)
    return result.to_dict(orient="records"), metadata


def audit(path, allow_unversioned_v1=False):
    rows, verification_metadata = load_original_support(path, allow_unversioned_v1)
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
        orig_label = bool(r["_OriginalLabelFirstSupport"])
        orig_cluster = bool(r["_OriginalClusterFirstSupport"])
        vc = cls(orig_label, orig_cluster)

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
        "verification_provenance": verification_metadata,
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


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-unversioned-v1", action="store_true")
    parser.add_argument("out_json")
    parser.add_argument("inputs", nargs="+", metavar="LABEL=CSV")
    return parser.parse_args()


def main():
    args = parse_args()
    out = args.out_json
    res = {}
    for spec in args.inputs:
        label, _, path = spec.partition("=")
        if not label or not path:
            raise SystemExit(f"invalid input specification: {spec!r}; expected LABEL=CSV")
        res[label] = {
            "source_csv": path,
            **audit(path, allow_unversioned_v1=args.allow_unversioned_v1),
        }
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
