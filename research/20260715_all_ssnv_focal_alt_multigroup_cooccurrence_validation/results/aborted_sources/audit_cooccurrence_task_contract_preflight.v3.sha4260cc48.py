#!/usr/bin/env python3
"""Preflight every reduced cooccurrence worker task before BAM analysis."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import scipy

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import analyze_methyl_ssnv_cooccurrence as ANALYZER


SCHEMA_NAME = "intersubmod.cooccurrence_task_contract_preflight"
SCHEMA_VERSION = "1.2.0"
EXPECTED_TASKS = 102_842
EXPECTED_ELIGIBLE = 919
EXPECTED_EVALUABLE = 1_867
EXPECTED_STATUS_COUNTS = {
    "ELIGIBLE_M2_RESIDUAL_UNEXPLAINED_AND_AXES_DETERMINATE": 919,
    "INELIGIBLE_NOT_M2_PHASE_ANCHORED_ROBUST": 808,
    "INELIGIBLE_SCREEN_HP_CONFOUND": 4,
    "INELIGIBLE_SCREEN_TECHNICAL_CONFOUND": 136,
    "NOT_EVALUABLE_M2_AXIS_INDETERMINATE": 100_974,
    "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM": 1,
}
EXPECTED_RAW_CONSTANT_COUNTS = {"hp_exact": 93_534, "hp_family": 93_534, "strand": 39}
EXPECTED_GATE_CONSTANT_COUNTS = {"hp_exact": 93_533, "hp_family": 93_533, "strand": 39}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def normalized_status(value: Any) -> str:
    return str(value).split(":", 1)[0]


def summarize_tasks(tasks: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    raw_constant_counts: Counter[str] = Counter()
    gate_constant_counts: Counter[str] = Counter()
    status_counts: Counter[str] = Counter()
    eligible = 0
    evaluable = 0
    group_limit_examples: list[dict[str, Any]] = []
    for task in tasks:
        assignment = task["assignment"]
        screen_row = task["screen_row"]
        levels = ANALYZER.m2_categorical_level_counts(assignment, screen_row)
        for axis, count in levels.items():
            if int(count) == 1:
                raw_constant_counts[axis] += 1
        gate = ANALYZER.m2_screen_eligibility(screen_row, levels)
        status = normalized_status(gate["status"])
        status_counts[status] += 1
        eligible += bool(gate["eligible"])
        evaluable += bool(gate["evaluable"])
        gate_constant_counts.update(str(axis) for axis in gate.get("constant_axes", []))
        if status == "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM":
            group_limit_examples.append(
                {
                    "sample": assignment["sample"],
                    "chrom": assignment["chrom"],
                    "pos": int(assignment["pos"]),
                    "ref": screen_row["ref"],
                    "alt": screen_row["alt"],
                    "raw_categorical_level_counts": levels,
                    "gate_constant_axes": list(gate.get("constant_axes", [])),
                }
            )
    return {
        "task_count": len(tasks),
        "eligible": eligible,
        "evaluable": evaluable,
        "evaluable_ineligible": evaluable - eligible,
        "status_counts": dict(sorted(status_counts.items())),
        "raw_assignment_constant_axis_site_counts": dict(
            sorted(raw_constant_counts.items())
        ),
        "gate_evaluated_constant_axis_site_counts": dict(
            sorted(gate_constant_counts.items())
        ),
        "group_limit_examples": group_limit_examples,
    }


def runtime_statistical_api_probe() -> dict[str, Any]:
    labels_3x2 = ["g1"] * 10 + ["g2"] * 10 + ["g3"] * 10
    categories_3x2 = (
        ["R"] * 8
        + ["A"] * 2
        + ["R"] * 5
        + ["A"] * 5
        + ["R"] * 2
        + ["A"] * 8
    )
    pearson = ANALYZER.categorical_association(labels_3x2, categories_3x2)
    raw_pearson = ANALYZER.chi2_contingency(pearson["table"], correction=False)
    marker_pearson = ANALYZER.LIB.methyl_group_marker_association(
        labels_3x2, categories_3x2
    )
    raw_marker_pearson = ANALYZER.LIB.chi2_contingency(
        marker_pearson["table"], correction=False
    )
    labels_2x2 = ["g1"] * 10 + ["g2"] * 10
    categories_2x2 = ["R"] * 7 + ["A"] * 3 + ["R"] * 3 + ["A"] * 7
    fisher = ANALYZER.categorical_association(labels_2x2, categories_2x2)
    return {
        "scipy_version": scipy.__version__,
        "chi2_result_type": type(raw_pearson).__name__,
        "chi2_tuple_index_1_p": float(raw_pearson[1]),
        "pearson_3x2": pearson,
        "marker_chi2_tuple_index_1_p": float(raw_marker_pearson[1]),
        "library_marker_pearson_3x2": marker_pearson,
        "fisher_2x2": fisher,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--intersubmod-root", type=Path, required=True)
    parser.add_argument("--independent-m2-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    started_at_utc = now_utc()
    if os.path.lexists(args.output):
        raise FileExistsError(f"Refusing to overwrite output: {args.output}")
    independent = json.loads(args.independent_m2_audit.read_text(encoding="utf-8"))
    if independent.get("pass") is not True:
        raise RuntimeError("Independent M2 audit did not pass")
    gate_artifact = artifact(Path(ANALYZER.M2_GATE.__file__))
    if independent["inputs"]["production_gate_source_reference_only"] != gate_artifact:
        raise RuntimeError("Independent M2 audit and current production gate identity differ")

    manifest, entries, run_validation = ANALYZER.load_manifest(
        args.manifest,
        args.intersubmod_root,
        None,
        False,
    )
    selected_samples = {entry["sample"] for entry in entries}
    assignments = ANALYZER.load_assignments(args.assignments, selected_samples)
    sites = ANALYZER.load_stable_site_rows(
        args.sites,
        assignments,
        selected_samples,
        require_candidate_screen_fields=True,
    )
    tasks = ANALYZER.build_tasks(
        entries,
        assignments,
        sites,
        ANALYZER.DEFAULT_TOP_MARKERS,
        exact_state_space_ceiling=ANALYZER.EXACT_STATE_SPACE_CEILING,
    )
    observed = summarize_tasks(tasks)
    runtime_probe = runtime_statistical_api_probe()
    checks = {
        "task_count": observed["task_count"] == EXPECTED_TASKS,
        "assignment_count": len(assignments) == EXPECTED_TASKS,
        "site_count": len(sites) == EXPECTED_TASKS,
        "manifest_dataset_count": len(entries) == 7,
        "validated_sample_count": len(run_validation) == 7,
        "eligible": observed["eligible"] == EXPECTED_ELIGIBLE,
        "evaluable": observed["evaluable"] == EXPECTED_EVALUABLE,
        "evaluable_ineligible": (
            observed["evaluable_ineligible"] == EXPECTED_EVALUABLE - EXPECTED_ELIGIBLE
        ),
        "status_counts": observed["status_counts"] == EXPECTED_STATUS_COUNTS,
        "raw_assignment_constant_counts": (
            observed["raw_assignment_constant_axis_site_counts"]
            == EXPECTED_RAW_CONSTANT_COUNTS
        ),
        "gate_evaluated_constant_counts": (
            observed["gate_evaluated_constant_axis_site_counts"]
            == EXPECTED_GATE_CONSTANT_COUNTS
        ),
        "single_group_limit_example": len(observed["group_limit_examples"]) == 1,
        "group_limit_is_hcc1954_chr5_751076": bool(
            observed["group_limit_examples"]
            and observed["group_limit_examples"][0]["sample"] == "HCC1954"
            and observed["group_limit_examples"][0]["chrom"] == "chr5"
            and observed["group_limit_examples"][0]["pos"] == 751_076
        ),
        "independent_counts_match": independent["counts"]
        == {
            "all_rows": 469_849,
            "eligible": EXPECTED_ELIGIBLE,
            "evaluable_ineligible": EXPECTED_EVALUABLE - EXPECTED_ELIGIBLE,
            "m1_stable_rows": EXPECTED_TASKS,
            "not_evaluable_axis_indeterminate": 100_974,
            "not_evaluable_group_count_gt10": 1,
        },
        "independent_gate_constant_counts_match": (
            independent["constant_axis_assignment_proof_counts"]
            == EXPECTED_GATE_CONSTANT_COUNTS
        ),
        "runtime_pearson_3x2_testable": (
            runtime_probe["pearson_3x2"]["testable"] is True
        ),
        "runtime_pearson_3x2_uses_chi_square": (
            runtime_probe["pearson_3x2"]["analytic_test"]
            == "pearson_chi_square"
        ),
        "runtime_pearson_p_matches_tuple_index": abs(
            float(runtime_probe["pearson_3x2"]["p_analytic"])
            - float(runtime_probe["chi2_tuple_index_1_p"])
        )
        < 1e-15,
        "runtime_library_marker_pearson_3x2_testable": (
            runtime_probe["library_marker_pearson_3x2"]["testable"] is True
        ),
        "runtime_library_marker_pearson_uses_chi_square": (
            runtime_probe["library_marker_pearson_3x2"]["analytic_test"]
            == "pearson_chi_square_descriptive_only"
        ),
        "runtime_library_marker_pearson_p_matches_tuple_index": abs(
            float(runtime_probe["library_marker_pearson_3x2"]["p_analytic"])
            - float(runtime_probe["marker_chi2_tuple_index_1_p"])
        )
        < 1e-15,
        "runtime_fisher_2x2_testable": (
            runtime_probe["fisher_2x2"]["testable"] is True
        ),
        "runtime_fisher_2x2_uses_fisher_exact": (
            runtime_probe["fisher_2x2"]["analytic_test"]
            == "fisher_exact_2x2"
        ),
    }
    finished_at_utc = now_utc()
    payload = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "started_at_utc": started_at_utc,
        "finished_at_utc": finished_at_utc,
        "created_at_utc": finished_at_utc,
        "task_type": "B_comprehensive_validation",
        "command": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        "inputs": {
            "manifest": artifact(args.manifest),
            "assignments": artifact(args.assignments),
            "sites": artifact(args.sites),
            "independent_m2_audit": artifact(args.independent_m2_audit),
        },
        "code": {
            "preflight": artifact(Path(__file__)),
            "analyzer": artifact(Path(ANALYZER.__file__)),
            "ssnv_cooccurrence_lib": artifact(Path(ANALYZER.LIB.__file__)),
            "m2_screen_gate": gate_artifact,
        },
        "observed": {
            **observed,
            "runtime_statistical_api_probe": runtime_probe,
        },
        "checks": checks,
        "interpretation": {
            "raw_vs_gate_constant_count_difference": (
                "The single K=11 site has constant HP categories in the assignment, but M2 "
                "returns the group-limit status before evaluating axes; raw counts therefore "
                "exceed gate-evaluated counts by exactly one for hp_exact and hp_family."
            ),
            "scientific_result": False,
            "role": (
                "worker-task contract, M2 gate integration, and runtime statistical API "
                "compatibility preflight only"
            ),
        },
        "pass": all(checks.values()),
        "pass_semantics": "execution_and_task_contract_integrity_only_not_scientific_result",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    args.output.chmod(0o444)
    print(json.dumps({"output": str(args.output.resolve()), "checks": checks, "pass": payload["pass"]}))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
