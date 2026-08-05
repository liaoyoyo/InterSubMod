#!/usr/bin/env python3
"""Preflight every reduced cooccurrence worker task before BAM analysis."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Iterator, Mapping, Sequence

import scipy

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import analyze_methyl_ssnv_cooccurrence as ANALYZER


SCHEMA_NAME = "intersubmod.cooccurrence_task_contract_preflight"
SCHEMA_VERSION = "2.0.0"
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


def source_paths() -> dict[str, Path]:
    return {
        "preflight": Path(__file__).resolve(),
        **ANALYZER.code_source_paths(),
    }


def capture_source_identity() -> dict[str, dict[str, Any]]:
    return {name: artifact(path) for name, path in source_paths().items()}


def capture_source_modes() -> dict[str, str]:
    return {
        name: oct(path.stat().st_mode & 0o777)
        for name, path in source_paths().items()
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


def audit_raw_identity_task(task: dict[str, Any]) -> dict[str, Any]:
    assignment = task["assignment"]
    screen_row = task["screen_row"]
    focal = ANALYZER.LIB.Variant(
        str(assignment["chrom"]),
        int(assignment["pos"]),
        str(screen_row["ref"]),
        str(screen_row["alt"]),
    )
    raw_audit: dict[str, Any] = {}
    recovered = ANALYZER.recover_site_reads(
        assignment,
        focal,
        [],
        task["entry"],
        raw_audit,
    )
    ANALYZER.validate_raw_identity_audit(
        raw_audit,
        expected_recovered_count=len(recovered),
    )
    observed_digest = ANALYZER.focal_alt_readset_digest(recovered)
    expected_digest = str(screen_row.get("alt_readset_sha256") or "")
    if observed_digest != expected_digest:
        raise RuntimeError(
            f"Raw preflight focal ALT readset mismatch at {assignment['sample']}:"
            f"{focal.chrom}:{focal.pos}: raw={observed_digest} expected={expected_digest}"
        )
    result = {
        "sample": str(assignment["sample"]),
        "chrom": focal.chrom,
        "pos": focal.pos,
        "n_all_focal_ref_alt_reads": len(recovered),
        "raw_identity_expected_projections": int(raw_audit["expected_projections"]),
        "raw_identity_matched_records": int(raw_audit["matched_records"]),
        "raw_identity_duplicate_projections_collapsed": int(
            raw_audit["duplicate_projections_collapsed"]
        ),
        "raw_identity_duplicate_extra_records_collapsed": int(
            raw_audit["duplicate_extra_records_collapsed"]
        ),
        "raw_identity_exact_duplicate_projections_collapsed": int(
            raw_audit["exact_duplicate_projections_collapsed"]
        ),
        "raw_identity_rg_only_duplicate_projections_collapsed": int(
            raw_audit["rg_only_duplicate_projections_collapsed"]
        ),
        "raw_identity_duplicate_projection_sha256": str(
            raw_audit["duplicate_projection_sha256"]
        ),
        "duplicate_projection_examples": raw_audit["duplicate_projection_examples"],
        "duplicate_projection_records": raw_audit["duplicate_projection_records"],
    }
    ANALYZER.validate_raw_identity_site_fields(
        result,
        expected_recovered_count=len(recovered),
    )
    return result


def audit_raw_identity_chunk(tasks: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    return [audit_raw_identity_task(task) for task in tasks]


def chunked(values: Sequence[Any], size: int) -> Iterator[list[Any]]:
    for offset in range(0, len(values), size):
        yield list(values[offset : offset + size])


def bounded_raw_identity_map(
    tasks: Sequence[dict[str, Any]],
    *,
    workers: int,
    chunk_size: int,
    max_pending_chunks: int,
) -> Iterator[dict[str, Any]]:
    task_chunks = iter(chunked(tasks, chunk_size))
    if workers == 1:
        for task_chunk in task_chunks:
            yield from audit_raw_identity_chunk(task_chunk)
        return
    with ProcessPoolExecutor(max_workers=workers) as executor:
        pending: dict[Any, int] = {}
        completed_by_index: dict[int, list[dict[str, Any]]] = {}
        next_index = 0
        emit_index = 0

        def submit_one() -> bool:
            nonlocal next_index
            try:
                task_chunk = next(task_chunks)
            except StopIteration:
                return False
            pending[executor.submit(audit_raw_identity_chunk, task_chunk)] = next_index
            next_index += 1
            return True

        while len(pending) < max_pending_chunks and submit_one():
            pass
        try:
            while pending:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    index = pending.pop(future)
                    completed_by_index[index] = future.result()
                    submit_one()
                while emit_index in completed_by_index:
                    yield from completed_by_index.pop(emit_index)
                    emit_index += 1
        except Exception:
            for future in pending:
                future.cancel()
            raise


def summarize_raw_identity(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    site_digest = hashlib.sha256()
    per_sample: dict[str, Counter[str]] = {}
    duplicate_examples: list[dict[str, Any]] = []
    for row in rows:
        site_digest.update(
            ANALYZER.json_text(
                ANALYZER.raw_identity_site_digest_payload(row)
            ).encode("utf-8")
        )
        site_digest.update(b"\n")
        counter = per_sample.setdefault(str(row["sample"]), Counter())
        counter["site_tasks"] += 1
        for source, target in (
            ("raw_identity_expected_projections", "expected_projection_occurrences"),
            ("raw_identity_matched_records", "matched_record_occurrences"),
            (
                "raw_identity_duplicate_projections_collapsed",
                "duplicate_projection_occurrences_collapsed",
            ),
            (
                "raw_identity_duplicate_extra_records_collapsed",
                "duplicate_extra_record_occurrences_collapsed",
            ),
            (
                "raw_identity_exact_duplicate_projections_collapsed",
                "exact_duplicate_projection_occurrences_collapsed",
            ),
            (
                "raw_identity_rg_only_duplicate_projections_collapsed",
                "rg_only_duplicate_projection_occurrences_collapsed",
            ),
        ):
            counter[target] += int(row[source])
        if int(row["raw_identity_duplicate_projections_collapsed"]) > 0:
            counter["sites_with_collapsed_duplicates"] += 1
            if len(duplicate_examples) < 100:
                duplicate_examples.append(
                    {
                        "sample": row["sample"],
                        "chrom": row["chrom"],
                        "pos": int(row["pos"]),
                        "duplicate_projections": int(
                            row["raw_identity_duplicate_projections_collapsed"]
                        ),
                        "projection_examples": row["duplicate_projection_examples"],
                    }
                )
    aggregate = Counter()
    for counter in per_sample.values():
        aggregate.update(counter)
    incident = next(
        (
            row
            for row in rows
            if row["sample"] == "HCC1395"
            and row["chrom"] == "chr1"
            and int(row["pos"]) == 13_200_585
        ),
        None,
    )
    return {
        "equivalence_contract": ANALYZER.RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
        "allowed_differing_auxiliary_tags": sorted(
            ANALYZER.RAW_DUPLICATE_ALLOWED_DIFFERING_TAGS
        ),
        "count_semantics": (
            "site-weighted projection occurrences; one long read can occur at multiple focal sites"
        ),
        "aggregate": dict(sorted(aggregate.items())),
        "per_sample": {
            sample: dict(sorted(counter.items()))
            for sample, counter in sorted(per_sample.items())
        },
        "site_weighted_audit_sha256": site_digest.hexdigest(),
        "duplicate_examples_first_100_sites": duplicate_examples,
        "v5_incident_site_replay": dict(incident) if incident is not None else None,
        "all_result_rows_passed_invariant_validation": True,
        "missing_projection_policy": ANALYZER.RAW_IDENTITY_MISSING_POLICY,
        "conflicting_analysis_payload_policy": ANALYZER.RAW_IDENTITY_CONFLICT_POLICY,
        "failure_counts_materialized": False,
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
    parser.add_argument("--primary-artifact-audit-pre", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--chunk-size", type=int, default=8)
    parser.add_argument("--max-pending-chunks", type=int, default=16)
    parser.add_argument("--progress-every", type=int, default=1_000)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    started_at_utc = now_utc()
    if os.path.lexists(args.output):
        raise FileExistsError(f"Refusing to overwrite output: {args.output}")
    if (
        args.workers < 1
        or args.chunk_size < 1
        or args.max_pending_chunks < args.workers
        or args.progress_every < 1
    ):
        raise ValueError("Invalid bounded raw-identity preflight settings")
    source_identity_before = capture_source_identity()
    source_modes_before = capture_source_modes()
    if set(source_modes_before.values()) != {"0o444"}:
        raise RuntimeError(f"Preflight sources are not mode 0444: {source_modes_before}")
    independent = json.loads(args.independent_m2_audit.read_text(encoding="utf-8"))
    if independent.get("pass") is not True:
        raise RuntimeError("Independent M2 audit did not pass")
    primary_audit = json.loads(args.primary_artifact_audit_pre.read_text(encoding="utf-8"))
    if primary_audit.get("pass") is not True:
        raise RuntimeError("Primary artifact audit did not pass")
    gate_artifact = source_identity_before["m2_screen_gate"]
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
    raw_identity_rows: list[dict[str, Any]] = []
    for row in bounded_raw_identity_map(
        tasks,
        workers=args.workers,
        chunk_size=args.chunk_size,
        max_pending_chunks=args.max_pending_chunks,
    ):
        raw_identity_rows.append(row)
        if len(raw_identity_rows) % args.progress_every == 0:
            print(
                json.dumps(
                    {
                        "stage": "raw_identity_preflight",
                        "completed": len(raw_identity_rows),
                        "total": len(tasks),
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
    raw_identity = summarize_raw_identity(raw_identity_rows)
    raw_aggregate = raw_identity["aggregate"]
    incident_replay = raw_identity["v5_incident_site_replay"] or {}
    source_identity_after = capture_source_identity()
    source_modes_after = capture_source_modes()
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
        "raw_identity_task_count": len(raw_identity_rows) == EXPECTED_TASKS,
        "raw_identity_all_samples": set(raw_identity["per_sample"])
        == {entry["sample"] for entry in entries},
        "raw_identity_nonempty_projection_population": int(
            raw_aggregate.get("expected_projection_occurrences", 0)
        )
        > 0,
        "raw_identity_all_result_rows_passed_invariant_validation": (
            raw_identity["all_result_rows_passed_invariant_validation"] is True
        ),
        "raw_identity_missing_policy_is_fail_closed": (
            raw_identity["missing_projection_policy"]
            == ANALYZER.RAW_IDENTITY_MISSING_POLICY
        ),
        "raw_identity_conflict_policy_is_fail_closed": (
            raw_identity["conflicting_analysis_payload_policy"]
            == ANALYZER.RAW_IDENTITY_CONFLICT_POLICY
        ),
        "raw_identity_failure_counts_not_materialized": (
            raw_identity["failure_counts_materialized"] is False
        ),
        "v5_incident_site_replayed": (
            incident_replay.get("sample") == "HCC1395"
            and incident_replay.get("chrom") == "chr1"
            and int(incident_replay.get("pos", -1)) == 13_200_585
            and int(
                incident_replay.get(
                    "raw_identity_rg_only_duplicate_projections_collapsed", -1
                )
            )
            == 2
            and int(
                incident_replay.get(
                    "raw_identity_exact_duplicate_projections_collapsed", -1
                )
            )
            == 0
        ),
        "source_identity_unchanged_during_preflight": (
            source_identity_before == source_identity_after
        ),
        "source_files_read_only_at_start": set(source_modes_before.values()) == {"0o444"},
        "source_files_read_only_at_finish": set(source_modes_after.values()) == {"0o444"},
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
            "primary_artifact_audit_pre": artifact(args.primary_artifact_audit_pre),
        },
        "code": {
            "source_identity_before": source_identity_before,
            "source_identity_after": source_identity_after,
            "source_modes_before": source_modes_before,
            "source_modes_after": source_modes_after,
        },
        "parameters": {
            "workers": args.workers,
            "chunk_size": args.chunk_size,
            "max_pending_chunks": args.max_pending_chunks,
            "progress_every": args.progress_every,
        },
        "observed": {
            **observed,
            "runtime_statistical_api_probe": runtime_probe,
            "raw_bam_identity_recovery": raw_identity,
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
                "worker-task contract, M2 gate integration, runtime statistical API, and "
                "all-task raw-BAM identity recovery preflight only"
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
