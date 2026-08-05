#!/usr/bin/env python3
"""
fn_verdict_readback_audit.py

T1 + T2 of the ISM structure<->label false-negative audit (workflow wf_fe623383-64b).

Pure OFFLINE read-back over an existing significance_summary.csv — NO C++ rerun,
NO fabrication. Quantifies, per region:
  Stage-1 (per-gate attrition): how many regions each gate/floor drops, flag big ones.
  Stage-2 (verdict-ignored FN): among legacy L4 class in {Noise, Weak}, how many
           already carry a VALID, significant PERMANOVA (cluster / label-HP / label-allele)
           that the discrete 2x2 verdict ignored (FN1), broken down by axis, plus the
           subset significant ONLY on an axis the global gate cannot open (FN10:
           sample / hp_fine), plus the FN3 presentational-vs-decision split
           (CramersV==0 with a significant PERMANOVA: still PassedGating?).

Usage:
  python3 fn_verdict_readback_audit.py OUT.json LABEL1=CSV1 [LABEL2=CSV2 ...]
      [--allow-unversioned-v1]

Every number written is a deterministic count over the CSV — grep-able back to source.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    SchemaContractError,
    select_legacy_view,
)


def fbool(s):
    return str(s).strip().lower() in ("true", "1", "yes", "t")


def ffloat(s):
    try:
        v = float(s)
        return v
    except (TypeError, ValueError):
        return math.nan


def sig(p_field, valid_field, row, alpha=0.05):
    """valid && p<=alpha"""
    if valid_field is not None and not fbool(row.get(valid_field, "")):
        return False
    p = ffloat(row.get(p_field, ""))
    return (not math.isnan(p)) and p <= alpha


def audit_csv(path, allow_unversioned_v1=False):
    frame = pd.read_csv(path, low_memory=False)
    view = select_legacy_view(
        frame,
        allow_unversioned_v1=allow_unversioned_v1,
        unversioned_unknown_policy="exclude",
    )
    frame = frame.copy()
    frame["_VerificationClass_L4"] = view.values
    rows = frame.to_dict("records")
    n = len(rows)

    # ---- Stage 1: per-gate attrition ----
    vclass = {}
    cnt = {
        "passed_gating_false": 0,
        "cluster_permanova_invalid": 0,
        "label_hp_permanova_invalid": 0,
        "label_allele_permanova_invalid": 0,
        "numreads_lt10": 0,
        "numreads_lt20": 0,
        "numcpgs_lt10": 0,
        "cluster_dispersion_warn": 0,
        "cramersv_zero": 0,
        "significant_true": 0,
        "suggestfilter_true": 0,
        "potential_loh_true": 0,
    }
    for r in rows:
        vc = r.get("_VerificationClass_L4")
        if not pd.isna(vc):
            vc = str(vc)
            vclass[vc] = vclass.get(vc, 0) + 1
        if not fbool(r.get("PassedGating", "")):
            cnt["passed_gating_false"] += 1
        if not fbool(r.get("ClusterPermanovaValid", "")):
            cnt["cluster_permanova_invalid"] += 1
        if not fbool(r.get("LabelHPPermanovaValid", "")):
            cnt["label_hp_permanova_invalid"] += 1
        if not fbool(r.get("LabelAllelePermanovaValid", "")):
            cnt["label_allele_permanova_invalid"] += 1
        nr = ffloat(r.get("NumReads", ""))
        if not math.isnan(nr) and nr < 10:
            cnt["numreads_lt10"] += 1
        if not math.isnan(nr) and nr < 20:
            cnt["numreads_lt20"] += 1
        nc = ffloat(r.get("NumCpGs", ""))
        if not math.isnan(nc) and nc < 10:
            cnt["numcpgs_lt10"] += 1
        if fbool(r.get("ClusterDispersionWarn", "")):
            cnt["cluster_dispersion_warn"] += 1
        cv = ffloat(r.get("CramersV", ""))
        if not math.isnan(cv) and cv == 0.0:
            cnt["cramersv_zero"] += 1
        if fbool(r.get("Significant", "")):
            cnt["significant_true"] += 1
        if fbool(r.get("SuggestFilter", "")):
            cnt["suggestfilter_true"] += 1
        if fbool(r.get("Potential_LOH", "")):
            cnt["potential_loh_true"] += 1

    # ---- Stage 2: verdict-ignored significant PERMANOVA (FN1/FN10) ----
    # per-region significance flags on the tests the verdict does NOT consult
    def permsig_cluster(r):
        return sig("ClusterPermanovaP", "ClusterPermanovaValid", r)

    def permsig_labelhp(r):
        return sig("LabelHPPermanovaP", "LabelHPPermanovaValid", r)

    def permsig_labelallele(r):
        return sig("LabelAllelePermanovaP", "LabelAllelePermanovaValid", r)

    def sample_axis_sig(r):  # FN10 orphaned tumor-vs-normal axis
        p = ffloat(r.get("SampleASM_P", ""))
        return (not math.isnan(p)) and p <= 0.05

    def hpfine_sig(r):  # FN10 hp_fine dead-ended axis
        p = ffloat(r.get("GlobalP_HPFine", ""))
        return (not math.isnan(p)) and p <= 0.05

    stage2 = {}
    for target in ("Noise", "Weak"):
        sub = [
            r for r in rows
            if not pd.isna(r.get("_VerificationClass_L4"))
            and r.get("_VerificationClass_L4") == target
        ]
        m = len(sub)
        any_perm = sum(
            1 for r in sub
            if permsig_cluster(r) or permsig_labelhp(r) or permsig_labelallele(r)
        )
        only_cluster = sum(1 for r in sub if permsig_cluster(r))
        only_labelhp = sum(1 for r in sub if permsig_labelhp(r))
        only_labelallele = sum(1 for r in sub if permsig_labelallele(r))
        # FN10: significant on an axis the gate can't open, among this class
        sample_orphan = sum(1 for r in sub if sample_axis_sig(r))
        hpfine_orphan = sum(1 for r in sub if hpfine_sig(r))
        # label-side PERMANOVA sig (the un-gated rescue tests that the verdict ignores)
        labelside_perm = sum(1 for r in sub if permsig_labelhp(r) or permsig_labelallele(r))
        stage2[target] = {
            "class_total": m,
            "with_any_valid_sig_permanova": any_perm,
            "frac_any_perm": round(any_perm / m, 4) if m else None,
            "by_arm_cluster": only_cluster,
            "by_arm_labelHP": only_labelhp,
            "by_arm_labelAllele": only_labelallele,
            "labelside_perm_sig": labelside_perm,
            "sample_axis_sig_orphan_FN10": sample_orphan,
            "hpfine_axis_sig_orphan_FN10": hpfine_orphan,
        }

    # ---- T2: FN3 presentational-vs-decision (CramersV==0 with sig PERMANOVA) ----
    cv0 = [r for r in rows if (not math.isnan(ffloat(r.get("CramersV", "")))) and ffloat(r.get("CramersV", "")) == 0.0]
    cv0_with_sigperm = [
        r for r in cv0
        if permsig_cluster(r) or permsig_labelhp(r) or permsig_labelallele(r)
    ]
    cv0_sigperm_passed = sum(1 for r in cv0_with_sigperm if fbool(r.get("PassedGating", "")))
    cv0_sigperm_blocked = sum(1 for r in cv0_with_sigperm if not fbool(r.get("PassedGating", "")))
    fn3 = {
        "cramersv_zero_total": len(cv0),
        "cramersv_zero_with_sig_permanova": len(cv0_with_sigperm),
        "of_those_still_passed_gating_PRESENTATIONAL": cv0_sigperm_passed,
        "of_those_gate_blocked_DECISION_LEVEL": cv0_sigperm_blocked,
    }

    return {
        "n_rows": n,
        "verification_contract": view.metadata(),
        "verification_unknown_excluded_rows": int(sum(view.unknown_counts.values())),
        "stage1_per_gate_attrition": cnt,
        "verification_class_counts": dict(sorted(vclass.items(), key=lambda kv: -kv[1])),
        "stage2_verdict_ignored_permanova_FN1_FN10": stage2,
        "T2_cramersv_zero_FN3": fn3,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out_json")
    parser.add_argument("inputs", nargs="+", metavar="LABEL=CSV")
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help=(
            "Explicitly authorize historical unversioned VerificationClass input; "
            "unknown values are excluded from Noise/Weak cohorts and counted."
        ),
    )
    args = parser.parse_args()
    out_path = args.out_json
    result = {}
    for spec in args.inputs:
        label, _, path = spec.partition("=")
        if not label or not path:
            parser.error(f"input must be LABEL=CSV, observed {spec!r}")
        try:
            audit = audit_csv(path, allow_unversioned_v1=args.allow_unversioned_v1)
        except SchemaContractError as exc:
            raise SystemExit(f"[fn-verdict-audit][schema-contract] {path}: {exc}") from exc
        result[label] = {"source_csv": path, **audit}
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2)
    # human-readable echo to stdout
    for label, d in result.items():
        print(f"\n===== {label}  ({d['n_rows']} rows)  {d['source_csv']} =====")
        print("  Verification contract:", d["verification_contract"])
        print("  Unknown rows excluded:", d["verification_unknown_excluded_rows"])
        print("  VerificationClass:", d["verification_class_counts"])
        print("  Stage1 attrition :", d["stage1_per_gate_attrition"])
        for cls, s in d["stage2_verdict_ignored_permanova_FN1_FN10"].items():
            print(f"  Stage2 {cls:5s}: {s['with_any_valid_sig_permanova']}/{s['class_total']} "
                  f"(frac={s['frac_any_perm']}) carry a VALID sig PERMANOVA verdict ignored | "
                  f"labelside={s['labelside_perm_sig']} sampleOrphan={s['sample_axis_sig_orphan_FN10']} "
                  f"hpfineOrphan={s['hpfine_axis_sig_orphan_FN10']}")
        print("  T2 FN3 CramersV=0:", d["T2_cramersv_zero_FN3"])
    print(f"\n[written] {out_path}")


if __name__ == "__main__":
    main()
