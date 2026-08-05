#!/usr/bin/env python3
"""Finalize a Task-B cooccurrence release after runner-level reconciliation."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
from typing import Any, Iterator, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
TOPIC_ROOT = SCRIPT_DIR.parent
ANALYZER_PATH = SCRIPT_DIR / "analyze_methyl_ssnv_cooccurrence.py"
ANALYZER_SPEC = importlib.util.spec_from_file_location(
    "cooccurrence_release_analyzer_contract", ANALYZER_PATH
)
if ANALYZER_SPEC is None or ANALYZER_SPEC.loader is None:
    raise RuntimeError(f"Cannot load analyzer contract: {ANALYZER_PATH}")
ANALYZER = importlib.util.module_from_spec(ANALYZER_SPEC)
sys.modules[ANALYZER_SPEC.name] = ANALYZER
ANALYZER_SPEC.loader.exec_module(ANALYZER)
SOURCE_AUTHORITY = ANALYZER.SOURCE_AUTHORITY

SCHEMA_NAME = "intersubmod.cooccurrence_release_receipt"
SCHEMA_VERSION = "1.0.0"
RAW_IDENTITY_RELEASE_CONTRACT = ANALYZER.RAW_IDENTITY_RELEASE_CONTRACT
ANALYSIS_SCOPE_IDENTITY_CONTRACT = ANALYZER.ANALYSIS_SCOPE_IDENTITY_CONTRACT
PREFLIGHT_CHECK_KEYS = ANALYZER.PREFLIGHT_CANONICAL_CHECK_KEYS
EXPECTED_STABLE_SITE_TASKS = ANALYZER.EXPECTED_STABLE_SITE_TASKS
CANONICAL_CODE_PATHS = ANALYZER.code_source_paths()
PREFLIGHT_SOURCE = SCRIPT_DIR / "audit_cooccurrence_task_contract_preflight.py"
RUNNER_SOURCE = SCRIPT_DIR / "run_cooccurrence_v6_source_locked.sh"
RESULT_ROOT = TOPIC_ROOT / "results"
WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
SCREEN_ROOT = (
    WORKSPACE_ROOT
    / "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
)
CANONICAL_MANIFEST = RESULT_ROOT / "all_ssnv_input_manifest.json"
CANONICAL_ASSIGNMENTS = SCREEN_ROOT / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
CANONICAL_SCREEN_SITES = SCREEN_ROOT / "all_ssnv_site_results.tsv.gz"
CANONICAL_PRIMARY_AUDIT = (
    RESULT_ROOT
    / "stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream.json"
)
CANONICAL_INDEPENDENT_M2_AUDIT = RESULT_ROOT / "independent_m2_gate_recount.v3.json"
CANONICAL_INTERSUBMOD_ROOT = WORKSPACE_ROOT / "intersubmod_all_ssnv_v2_verification_fix"
CANONICAL_PREFLIGHT = (
    RESULT_ROOT / "cooccurrence_task_contract_preflight.v9_command_parity_full_runtime.json"
)
CANONICAL_PRODUCER_DIR = (
    WORKSPACE_ROOT
    / "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity"
)
CANONICAL_SAMPLES = frozenset(
    {"HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"}
)
RELEASE_CHECK_KEYS = frozenset(
    {
        "preflight_schema_scope_pass",
        "preflight_command_and_parameters_exact",
        "producer_schema_scope_pass",
        "summary_schema_scope_pass",
        "release_contract_exact",
        "preflight_inputs_and_lock_reconciled",
        "upstream_audits_bound_to_current_screen_inputs",
        "upstream_audits_bound_to_source_authority",
        "producer_inputs_and_lock_reconciled",
        "producer_command_and_parameters_exact",
        "producer_output_artifacts_reconciled",
        "runtime_environment_and_probe_chain_reconciled",
        "python_cache_contract_reconciled",
        "producer_source_lock_reconciled",
        "preflight_source_chain_includes_finalizer_and_runner",
        "analysis_scope_contract_exact",
        "site_weighted_analysis_digest_independently_recomputed",
        "summary_counts_independently_recomputed",
        "all_full_scope_rows_and_counts_independently_recomputed",
        "chronology_reconciled",
    }
)


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


def load_json(path: Path, label: str) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"{label} is not a JSON object")
    return value


def require_artifact(record: Any, path: Path, label: str) -> None:
    if not isinstance(record, Mapping) or dict(record) != artifact(path):
        raise RuntimeError(f"{label} artifact identity drifted")


def require_mode_0444(path: Path, label: str) -> None:
    if oct(path.stat().st_mode & 0o777) != "0o444":
        raise RuntimeError(f"{label} is not mode 0444: {path}")


def require_exact_path(path: Path, expected: Path, label: str) -> None:
    if path.resolve() != expected.resolve():
        raise RuntimeError(
            f"{label} path is not canonical: observed={path.resolve()} expected={expected.resolve()}"
        )


def parse_int(row: Mapping[str, str], field: str) -> int:
    try:
        return int(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError(f"Invalid integer field {field}: {row.get(field)!r}") from error


def parse_json_list(row: Mapping[str, str], field: str) -> list[Any]:
    try:
        value = json.loads(row[field])
    except (KeyError, TypeError, json.JSONDecodeError) as error:
        raise RuntimeError(f"Invalid JSON list field {field}: {row.get(field)!r}") from error
    if not isinstance(value, list):
        raise RuntimeError(f"Expected JSON list in {field}")
    return value


def iter_tsv_rows(
    path: Path, expected_columns: Sequence[str]
) -> Iterator[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(expected_columns):
            raise RuntimeError(
                f"{path.name} header drifted: observed={reader.fieldnames} expected={list(expected_columns)}"
            )
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise RuntimeError(f"Malformed TSV row at {path}:{line_number}")
            yield row


def require_equal(observed: Any, expected: Any, label: str) -> None:
    if observed != expected:
        raise RuntimeError(f"{label} mismatch: observed={observed!r} expected={expected!r}")


def _site_key(row: Mapping[str, Any]) -> tuple[str, str, int]:
    return str(row["sample"]), str(row["chrom"]), int(row["pos"])


def recompute_output_contract(
    *,
    sites_path: Path,
    pairs_path: Path,
    duplicates_path: Path,
    oracle_path: Path,
    expected_stable_sites: int = EXPECTED_STABLE_SITE_TASKS,
) -> dict[str, Any]:
    site_rows: dict[tuple[str, str, int], dict[str, Any]] = {}
    site_digest = hashlib.sha256()
    raw_counts: Counter[str] = Counter()
    previous_site_key: tuple[str, str, int] | None = None
    candidate_count = 0
    m2_eligible_site_count = 0

    for row in iter_tsv_rows(sites_path, ANALYZER.SITE_COLUMNS):
        key = _site_key(row)
        if previous_site_key is not None and key <= previous_site_key:
            raise RuntimeError(f"Site rows are not unique and strictly sorted at {key}")
        previous_site_key = key
        candidate = ANALYZER.parse_bool(
            row["multi_marker_molecular_haplotype_base_candidate"]
        )
        m2_eligible = ANALYZER.parse_bool(row["m2_screen_eligible"])
        candidate_count += candidate
        m2_eligible_site_count += m2_eligible
        expected = parse_int(row, "raw_identity_expected_projections")
        matched = parse_int(row, "raw_identity_matched_records")
        duplicate_projections = parse_int(
            row, "raw_identity_duplicate_projections_collapsed"
        )
        duplicate_extras = parse_int(
            row, "raw_identity_duplicate_extra_records_collapsed"
        )
        exact_duplicates = parse_int(
            row, "raw_identity_exact_duplicate_projections_collapsed"
        )
        rg_duplicates = parse_int(
            row, "raw_identity_rg_only_duplicate_projections_collapsed"
        )
        typed_digest_row: dict[str, Any] = dict(row)
        for field, value in (
            ("n_all_focal_ref_alt_reads", parse_int(row, "n_all_focal_ref_alt_reads")),
            ("raw_identity_expected_projections", expected),
            ("raw_identity_matched_records", matched),
            ("raw_identity_duplicate_projections_collapsed", duplicate_projections),
            ("raw_identity_duplicate_extra_records_collapsed", duplicate_extras),
            ("raw_identity_exact_duplicate_projections_collapsed", exact_duplicates),
            ("raw_identity_rg_only_duplicate_projections_collapsed", rg_duplicates),
        ):
            typed_digest_row[field] = value
        payload = ANALYZER.raw_identity_site_digest_payload(typed_digest_row)
        site_digest.update(ANALYZER.json_text(payload).encode("utf-8"))
        site_digest.update(b"\n")
        raw_counts.update(
            {
                "expected_projection_occurrences": expected,
                "matched_record_occurrences": matched,
                "duplicate_projection_occurrences_collapsed": duplicate_projections,
                "duplicate_extra_record_occurrences_collapsed": duplicate_extras,
                "exact_duplicate_projection_occurrences_collapsed": exact_duplicates,
                "rg_only_duplicate_projection_occurrences_collapsed": rg_duplicates,
                "sites_with_collapsed_duplicates": int(duplicate_projections > 0),
            }
        )
        site_rows[key] = {
            "sample": key[0],
            "chrom": key[1],
            "pos": key[2],
            "ref": row["ref"],
            "alt": row["alt"],
            "n_all_focal_ref_alt_reads": parse_int(row, "n_all_focal_ref_alt_reads"),
            "n_partner_markers": parse_int(row, "n_partner_markers"),
            "n_pair_rows_reconciled": parse_int(row, "n_pair_rows_reconciled"),
            "pair_row_count_reconciliation_pass": ANALYZER.parse_bool(
                row["pair_row_count_reconciliation_pass"]
            ),
            "n_endpoint_a_testable_markers": parse_int(
                row, "n_endpoint_a_testable_markers"
            ),
            "n_endpoint_a_exact_identifiable_markers": parse_int(
                row, "n_endpoint_a_exact_identifiable_markers"
            ),
            "n_endpoint_a_exact_not_identifiable_markers": parse_int(
                row, "n_endpoint_a_exact_not_identifiable_markers"
            ),
            "n_endpoint_a_conditional_permutable_markers": parse_int(
                row, "n_endpoint_a_conditional_permutable_markers"
            ),
            "n_pair_bh_discoveries": parse_int(row, "n_pair_bh_discoveries"),
            "pair_bh_discovery_positions": parse_json_list(
                row, "pair_bh_discovery_positions"
            ),
            "n_pair_by_discoveries": parse_int(row, "n_pair_by_discoveries"),
            "pair_by_discovery_positions": parse_json_list(
                row, "pair_by_discovery_positions"
            ),
            "n_pair_by_confirmed": parse_int(row, "n_pair_by_confirmed"),
            "pair_by_confirmed_positions": parse_json_list(
                row, "pair_by_confirmed_positions"
            ),
            "n_spatially_separated_pair_by_20bp": parse_int(
                row, "n_spatially_separated_pair_by_20bp"
            ),
            "spatially_separated_pair_by_positions_20bp": parse_json_list(
                row, "spatially_separated_pair_by_positions_20bp"
            ),
            "top_marker_positions": parse_json_list(row, "top_marker_positions"),
            "n_top_marker_pair_by_confirmed": parse_int(
                row, "n_top_marker_pair_by_confirmed"
            ),
            "top_marker_pair_by_confirmed_positions": parse_json_list(
                row, "top_marker_pair_by_confirmed_positions"
            ),
            "joint_signature_complete_marker_effect_supported_positions": parse_json_list(
                row, "joint_signature_complete_marker_effect_supported_positions"
            ),
            "n_same_complete_read_effect_supported_top_pair_by": parse_int(
                row, "n_same_complete_read_effect_supported_top_pair_by"
            ),
            "same_complete_read_effect_supported_top_pair_by_positions": parse_json_list(
                row, "same_complete_read_effect_supported_top_pair_by_positions"
            ),
            "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp": parse_int(
                row,
                "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp",
            ),
            "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp": parse_json_list(
                row,
                "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp",
            ),
            "joint_signature_global_by_discovery": ANALYZER.parse_bool(
                row["joint_signature_global_by_discovery"]
            ),
            "m2_screen_eligible": m2_eligible,
            "candidate": candidate,
            "n_endpoint_a_pre_candidates": parse_int(
                row, "n_endpoint_a_pre_candidates"
            ),
            "n_confirmed_markers": parse_int(row, "n_confirmed_markers"),
            "confirmed_marker_positions": parse_json_list(
                row, "confirmed_marker_positions"
            ),
            "n_independent_confirmed_markers_20bp": parse_int(
                row, "n_independent_confirmed_markers_20bp"
            ),
            "genetically_anchored_multi_marker_candidate": ANALYZER.parse_bool(
                row["genetically_anchored_multi_marker_candidate"]
            ),
            "genetically_anchored_multi_marker_candidate_by_sensitivity": ANALYZER.parse_bool(
                row["genetically_anchored_multi_marker_candidate_by_sensitivity"]
            ),
            "raw_identity_duplicate_projections_collapsed": duplicate_projections,
            "raw_identity_duplicate_extra_records_collapsed": duplicate_extras,
            "raw_identity_exact_duplicate_projections_collapsed": exact_duplicates,
            "raw_identity_rg_only_duplicate_projections_collapsed": rg_duplicates,
            "raw_identity_duplicate_projection_sha256": row[
                "raw_identity_duplicate_projection_sha256"
            ],
            "raw_identity_analysis_site_payload_sha256": row[
                "raw_identity_analysis_site_payload_sha256"
            ],
            "alt_readset_sha256": row["alt_readset_sha256"],
            "alt_readset_representative": ANALYZER.parse_bool(
                row["alt_readset_representative"]
            ),
        }

    if len(site_rows) != expected_stable_sites:
        raise RuntimeError(
            f"Site row count is not expected scope: {len(site_rows)} != {expected_stable_sites}"
        )

    pair_counts: Counter[tuple[str, str, int]] = Counter()
    pair_metrics: dict[tuple[str, str, int], Counter[str]] = defaultdict(Counter)
    pair_positions: dict[tuple[str, str, int], dict[str, list[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    oracle_pair_positions: dict[tuple[str, str, int], list[int]] = defaultdict(list)
    previous_pair_key: tuple[Any, ...] | None = None
    focal_partner_pairs = 0
    pair_by_confirmed = 0
    callability_gate_pass = 0

    for row in iter_tsv_rows(pairs_path, ANALYZER.PAIR_COLUMNS):
        focal_key = (row["sample"], row["focal_chrom"], parse_int(row, "focal_pos"))
        partner_pos = parse_int(row, "partner_pos")
        sort_key = (
            row["sample"],
            row["focal_chrom"],
            focal_key[2],
            partner_pos,
            row["partner_ref"],
            row["partner_alt"],
        )
        if previous_pair_key is not None and sort_key <= previous_pair_key:
            raise RuntimeError(f"Pair rows are not unique and strictly sorted at {sort_key}")
        previous_pair_key = sort_key
        site = site_rows.get(focal_key)
        if site is None:
            raise RuntimeError(f"Pair focal is absent from site table: {focal_key}")
        require_equal(row["focal_ref"], site["ref"], "pair focal REF")
        require_equal(row["focal_alt"], site["alt"], "pair focal ALT")
        require_equal(row["partner_chrom"], row["focal_chrom"], "pair chromosome")
        distance = parse_int(row, "distance_bp")
        require_equal(distance, partner_pos - focal_key[2], "pair signed distance")
        if not 0 < abs(distance) <= ANALYZER.PAIR_WINDOW_BP:
            raise RuntimeError(f"Pair lies outside the canonical window: {sort_key}")
        require_equal(
            row["partner_universe_contract"],
            "all_frozen_latest_LongPhase-S_PASS_biallelic_sSNVs_within_focal_plus_or_minus_5000bp",
            "pair universe contract",
        )
        require_equal(
            parse_int(row, "n_all_focal_ref_alt_reads"),
            site["n_all_focal_ref_alt_reads"],
            "pair focal read denominator",
        )
        exact_by = ANALYZER.parse_bool(row["endpoint_a_exact_by_discovery"])
        effect = ANALYZER.parse_bool(row["endpoint_a_effect_gate_pass"])
        conditional = ANALYZER.parse_bool(
            row["endpoint_a_conditional_sensitivity_pass"]
        )
        callability = ANALYZER.parse_bool(row["callability_gate_pass"])
        formal = ANALYZER.parse_bool(row["endpoint_a_formal_pair_by_confirmed"])
        require_equal(
            formal,
            exact_by and effect and conditional and callability,
            "formal pair confirmation logic",
        )
        require_equal(
            ANALYZER.parse_bool(row["endpoint_a_confirmed_association"]),
            formal,
            "confirmed association alias",
        )
        require_equal(
            ANALYZER.parse_bool(row["endpoint_a_confirmed_by_sensitivity"]),
            formal,
            "confirmed sensitivity alias",
        )
        metrics = pair_metrics[focal_key]
        testable = ANALYZER.parse_bool(row["endpoint_a_testable"])
        permutable = ANALYZER.parse_bool(row["endpoint_a_permutable"])
        bh = ANALYZER.parse_bool(row["endpoint_a_exact_bh_discovery"])
        metrics["testable"] += testable
        metrics["exact_identifiable"] += row["endpoint_a_p_fixed_margins_exact"] != ""
        metrics["exact_not_identifiable"] += testable and row[
            "endpoint_a_p_fixed_margins_exact"
        ] == ""
        metrics["permutable"] += permutable
        metrics["bh"] += bh
        metrics["by"] += exact_by
        metrics["formal"] += formal
        metrics["exact_family"] += (
            row["endpoint_a_global_fdr_family_status"] == "ELIGIBLE_M2_EXACT_FAMILY"
        )
        pair_counts[focal_key] += 1
        focal_partner_pairs += 1
        pair_by_confirmed += formal
        callability_gate_pass += callability
        oracle_pair_positions[focal_key].append(partner_pos)
        if bh:
            pair_positions[focal_key]["bh"].append(partner_pos)
        if exact_by:
            pair_positions[focal_key]["by"].append(partner_pos)
        if formal:
            pair_positions[focal_key]["formal"].append(partner_pos)

    # Recompute every per-site pair field and the formal multi-marker candidate.
    candidate_count_recomputed = 0
    for key, site in site_rows.items():
        metrics = pair_metrics[key]
        positions = pair_positions[key]
        n_pairs = pair_counts[key]
        bh_positions = sorted(positions["bh"])
        by_positions = sorted(positions["by"])
        confirmed_positions = sorted(positions["formal"])
        spaced_positions = ANALYZER.LIB.spatially_separated_positions(
            confirmed_positions,
            minimum_separation=ANALYZER.SPACED_MARKER_MIN_SEPARATION_BP,
        )
        top_positions = {int(value) for value in site["top_marker_positions"]}
        top_confirmed = sorted(top_positions.intersection(confirmed_positions))
        top_spaced = ANALYZER.LIB.spatially_separated_positions(
            top_confirmed,
            minimum_separation=ANALYZER.SPACED_MARKER_MIN_SEPARATION_BP,
        )
        complete_effect = {
            int(value)
            for value in site[
                "joint_signature_complete_marker_effect_supported_positions"
            ]
        }
        same_complete = sorted(set(top_confirmed).intersection(complete_effect))
        same_complete_spaced = ANALYZER.LIB.spatially_separated_positions(
            same_complete,
            minimum_separation=ANALYZER.SPACED_MARKER_MIN_SEPARATION_BP,
        )
        candidate = bool(
            site["m2_screen_eligible"]
            and len(spaced_positions) >= 2
            and len(top_spaced) >= 2
            and len(same_complete_spaced) >= 2
            and site["joint_signature_global_by_discovery"]
        )
        expected_fields = {
            "n_partner_markers": n_pairs,
            "n_pair_rows_reconciled": n_pairs,
            "pair_row_count_reconciliation_pass": True,
            "n_endpoint_a_testable_markers": metrics["testable"],
            "n_endpoint_a_exact_identifiable_markers": metrics["exact_identifiable"],
            "n_endpoint_a_exact_not_identifiable_markers": metrics[
                "exact_not_identifiable"
            ],
            "n_endpoint_a_conditional_permutable_markers": metrics["permutable"],
            "n_pair_bh_discoveries": len(bh_positions),
            "pair_bh_discovery_positions": bh_positions,
            "n_pair_by_discoveries": len(by_positions),
            "pair_by_discovery_positions": by_positions,
            "n_pair_by_confirmed": len(confirmed_positions),
            "pair_by_confirmed_positions": confirmed_positions,
            "n_spatially_separated_pair_by_20bp": len(spaced_positions),
            "spatially_separated_pair_by_positions_20bp": spaced_positions,
            "n_top_marker_pair_by_confirmed": len(top_confirmed),
            "top_marker_pair_by_confirmed_positions": top_confirmed,
            "n_same_complete_read_effect_supported_top_pair_by": len(same_complete),
            "same_complete_read_effect_supported_top_pair_by_positions": same_complete,
            "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp": len(
                same_complete_spaced
            ),
            "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp": same_complete_spaced,
            "candidate": candidate,
            "n_endpoint_a_pre_candidates": len(bh_positions),
            "n_confirmed_markers": len(confirmed_positions),
            "confirmed_marker_positions": confirmed_positions,
            "n_independent_confirmed_markers_20bp": len(spaced_positions),
            "genetically_anchored_multi_marker_candidate": candidate,
            "genetically_anchored_multi_marker_candidate_by_sensitivity": candidate,
        }
        for field, expected in expected_fields.items():
            require_equal(site[field], expected, f"site-pair reconciliation {key} {field}")
        candidate_count_recomputed += candidate
    require_equal(candidate_count, candidate_count_recomputed, "candidate recomputation")

    duplicate_state: dict[tuple[str, str, int], dict[str, Any]] = defaultdict(
        lambda: {"projections": [], "extras": 0, "exact": 0, "rg": 0}
    )
    previous_duplicate_key: tuple[Any, ...] | None = None
    duplicate_rows = 0
    for row in iter_tsv_rows(duplicates_path, ANALYZER.RAW_IDENTITY_DUPLICATE_COLUMNS):
        key = _site_key(row)
        if key not in site_rows:
            raise RuntimeError(f"Duplicate audit references absent site: {key}")
        projection = ANALYZER.TAGS.projection_key(
            row["projection_read_name"],
            row["projection_chrom"],
            parse_int(row, "projection_start"),
            parse_int(row, "projection_end"),
            parse_int(row, "projection_mapq"),
            row["projection_strand"],
        )
        sort_key = (*key, *projection)
        if previous_duplicate_key is not None and sort_key <= previous_duplicate_key:
            raise RuntimeError(
                f"Duplicate audit rows are not unique and strictly sorted at {sort_key}"
            )
        previous_duplicate_key = sort_key
        record_count = parse_int(row, "record_count")
        if record_count < 2:
            raise RuntimeError("Duplicate audit record_count is below two")
        differing = parse_json_list(row, "differing_auxiliary_tags")
        classification = row["classification"]
        if classification == "EXACT_DUPLICATE":
            require_equal(differing, [], "exact duplicate differing tags")
        elif classification == "RG_ONLY_DUPLICATE":
            require_equal(differing, ["RG"], "RG-only duplicate differing tags")
        else:
            raise RuntimeError(f"Unexpected duplicate classification: {classification}")
        state = duplicate_state[key]
        state["projections"].append(projection)
        state["extras"] += record_count - 1
        state["exact"] += classification == "EXACT_DUPLICATE"
        state["rg"] += classification == "RG_ONLY_DUPLICATE"
        duplicate_rows += 1

    for key, site in site_rows.items():
        state = duplicate_state[key]
        require_equal(
            len(state["projections"]),
            site["raw_identity_duplicate_projections_collapsed"],
            f"duplicate projection count {key}",
        )
        require_equal(
            state["extras"],
            site["raw_identity_duplicate_extra_records_collapsed"],
            f"duplicate extra count {key}",
        )
        require_equal(
            state["exact"],
            site["raw_identity_exact_duplicate_projections_collapsed"],
            f"exact duplicate count {key}",
        )
        require_equal(
            state["rg"],
            site["raw_identity_rg_only_duplicate_projections_collapsed"],
            f"RG duplicate count {key}",
        )
        require_equal(
            ANALYZER.projection_digest(state["projections"]),
            site["raw_identity_duplicate_projection_sha256"],
            f"duplicate projection digest {key}",
        )

    oracle = load_json(oracle_path, "oracle")
    require_equal(
        oracle.get("schema_name"),
        f"{ANALYZER.SCHEMA_NAME}.oracle_cases",
        "oracle schema",
    )
    require_equal(oracle.get("schema_version"), ANALYZER.SCHEMA_VERSION, "oracle version")
    require_equal(oracle.get("partner_window_bp"), ANALYZER.PAIR_WINDOW_BP, "oracle window")
    focal_cases = oracle.get("focal_cases")
    if not isinstance(focal_cases, list) or len(focal_cases) != len(ANALYZER.ORACLE_FOCALS):
        raise RuntimeError("Oracle focal-case count drifted")
    canonical_oracle_scope = expected_stable_sites == EXPECTED_STABLE_SITE_TASKS
    if canonical_oracle_scope:
        canonical_focal_cases = focal_cases
    else:
        canonical_focal_cases = []
        expected_ids = {definition["case_id"] for definition in ANALYZER.ORACLE_FOCALS}
        observed_ids = set()
        for case in focal_cases:
            if not isinstance(case, Mapping):
                raise RuntimeError("Oracle focal case is not an object")
            case_id = case.get("case_id")
            if not isinstance(case_id, str) or case_id in observed_ids:
                raise RuntimeError("Synthetic oracle case IDs are missing or duplicated")
            observed_ids.add(case_id)
            require_equal(case.get("partner_window_oracle_pass"), True, "oracle pass flag")
        require_equal(observed_ids, expected_ids, "synthetic oracle case IDs")

    for definition, case in zip(ANALYZER.ORACLE_FOCALS, canonical_focal_cases):
        if not isinstance(case, Mapping):
            raise RuntimeError("Oracle focal case is not an object")
        key = (definition["sample"], definition["chrom"], int(definition["pos"]))
        actual_positions = sorted(oracle_pair_positions[key])
        expected_positions = list(definition["expected_partner_positions"])
        require_equal(case.get("case_id"), definition["case_id"], "oracle case id")
        require_equal(case.get("sample"), definition["sample"], "oracle sample")
        require_equal(case.get("focal"), {"chrom": key[1], "pos": key[2]}, "oracle focal")
        require_equal(case.get("present"), key in site_rows, "oracle presence")
        require_equal(case.get("expected_partner_positions"), expected_positions, "oracle expected positions")
        require_equal(case.get("observed_partner_positions"), actual_positions, "oracle observed positions")
        require_equal(actual_positions, expected_positions, "oracle canonical positions")
        require_equal(case.get("partner_window_oracle_pass"), True, "oracle pass flag")
        embedded_site = case.get("site")
        if not isinstance(embedded_site, Mapping):
            raise RuntimeError("Oracle embedded site is missing")
        require_equal(_site_key(embedded_site), key, "oracle embedded site key")
        require_equal(
            embedded_site.get("raw_identity_analysis_site_payload_sha256"),
            site_rows[key]["raw_identity_analysis_site_payload_sha256"],
            "oracle embedded site analysis digest",
        )
        embedded_pairs = case.get("pairs")
        if not isinstance(embedded_pairs, list):
            raise RuntimeError("Oracle embedded pairs are malformed")
        require_equal(
            sorted(int(pair["partner_pos"]) for pair in embedded_pairs),
            actual_positions,
            "oracle embedded pair positions",
        )

    shared = oracle.get("shared_readset_case")
    definition = ANALYZER.SHARED_READSET_ORACLE
    if not isinstance(shared, Mapping):
        raise RuntimeError("Shared-readset oracle is malformed")
    require_equal(shared.get("case_id"), definition["case_id"], "shared oracle id")
    if canonical_oracle_scope:
        shared_keys = [
            (definition["sample"], definition["chrom"], int(pos))
            for pos in definition["positions"]
        ]
        shared_sites = [site_rows[key] for key in shared_keys if key in site_rows]
        shared_digests = [site["alt_readset_sha256"] for site in shared_sites]
        same_readset = (
            len(shared_sites) == 2
            and len(set(shared_digests)) == 1
            and bool(shared_digests[0])
        )
        one_representative = (
            len(shared_sites) == 2
            and sum(site["alt_readset_representative"] for site in shared_sites) == 1
        )
        require_equal(shared.get("positions"), list(definition["positions"]), "shared positions")
        require_equal(
            shared.get("present_positions"),
            [key[2] for key in shared_keys if key in site_rows],
            "shared present positions",
        )
        require_equal(shared.get("same_alt_readset"), same_readset, "shared readset flag")
        require_equal(
            shared.get("one_alt_readset_representative"),
            one_representative,
            "shared representative flag",
        )
        if not same_readset or not one_representative:
            raise RuntimeError("Shared-readset oracle did not pass")
    else:
        require_equal(shared.get("same_alt_readset"), True, "synthetic shared readset flag")
        require_equal(
            shared.get("one_alt_readset_representative"),
            True,
            "synthetic shared representative flag",
        )

    return {
        "stable_sites": len(site_rows),
        "site_weighted_audit_sha256": site_digest.hexdigest(),
        "focal_partner_pairs": focal_partner_pairs,
        "pair_by_confirmed": pair_by_confirmed,
        "multi_marker_molecular_haplotype_base_candidates": candidate_count_recomputed,
        "m2_global_fdr_eligible_sites": m2_eligible_site_count,
        "m2_exact_global_fdr_family_pairs": sum(
            metrics["exact_family"] for metrics in pair_metrics.values()
        ),
        "pair_callability_gate_pass": callability_gate_pass,
        "raw_identity_expected_projection_occurrences": raw_counts[
            "expected_projection_occurrences"
        ],
        "raw_identity_matched_record_occurrences": raw_counts[
            "matched_record_occurrences"
        ],
        "raw_identity_sites_with_collapsed_duplicates": raw_counts[
            "sites_with_collapsed_duplicates"
        ],
        "raw_identity_duplicate_projection_occurrences_collapsed": raw_counts[
            "duplicate_projection_occurrences_collapsed"
        ],
        "raw_identity_duplicate_extra_record_occurrences_collapsed": raw_counts[
            "duplicate_extra_record_occurrences_collapsed"
        ],
        "raw_identity_exact_duplicate_projection_occurrences_collapsed": raw_counts[
            "exact_duplicate_projection_occurrences_collapsed"
        ],
        "raw_identity_rg_only_duplicate_projection_occurrences_collapsed": raw_counts[
            "rg_only_duplicate_projection_occurrences_collapsed"
        ],
        "raw_identity_sparse_duplicate_rows": duplicate_rows,
        "oracle_focal_cases": len(focal_cases),
        "oracle_shared_readset_pass": True,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--producer-receipt", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--duplicates", type=Path, required=True)
    parser.add_argument("--oracle", type=Path, required=True)
    parser.add_argument("--runner-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def canonical_command(args: argparse.Namespace) -> list[str]:
    return [
        *ANALYZER.canonical_python_prefix(args.output.resolve().parent),
        str(Path(__file__).resolve()),
        "--preflight",
        str(args.preflight.resolve()),
        "--producer-receipt",
        str(args.producer_receipt.resolve()),
        "--summary",
        str(args.summary.resolve()),
        "--sites",
        str(args.sites.resolve()),
        "--pairs",
        str(args.pairs.resolve()),
        "--duplicates",
        str(args.duplicates.resolve()),
        "--oracle",
        str(args.oracle.resolve()),
        "--runner-source",
        str(args.runner_source.resolve()),
        "--output",
        str(args.output.resolve()),
    ]


def observed_process_command() -> list[str]:
    raw = Path("/proc/self/cmdline").read_bytes()
    if not raw.endswith(b"\0"):
        raise RuntimeError("Process command line is unavailable or malformed")
    return [os.fsdecode(token) for token in raw[:-1].split(b"\0")]


def main(argv: Sequence[str] | None = None) -> int:
    if argv is not None:
        raise RuntimeError("Formal release finalizer is direct-CLI only")
    started_at_utc = now_utc()
    args = parse_args(argv)
    actual_command = observed_process_command()
    if actual_command != canonical_command(args):
        raise RuntimeError("Release finalizer command is not canonical")
    finalizer_python_cache = ANALYZER.capture_python_cache_contract(
        args.output.resolve().parent
    )
    if os.path.lexists(args.output):
        raise FileExistsError(f"Refusing to overwrite release receipt: {args.output}")
    source_authority = SOURCE_AUTHORITY.validate_release_source_authority(
        {
            "release_source_authority_validator",
            "cooccurrence_analyzer",
            "cooccurrence_preflight",
            "cooccurrence_release_finalizer",
            "cooccurrence_release_runner",
            "ssnv_cooccurrence_lib",
            "latest_tag_join",
            "m2_screen_gate",
            "independent_m2_auditor",
            "primary_artifact_auditor",
        }
    )
    canonical_output_paths = {
        "producer_receipt": CANONICAL_PRODUCER_DIR / ANALYZER.RECEIPT_OUTPUT_NAME,
        "summary": CANONICAL_PRODUCER_DIR / ANALYZER.SUMMARY_OUTPUT_NAME,
        "sites": CANONICAL_PRODUCER_DIR / ANALYZER.SITE_OUTPUT_NAME,
        "pairs": CANONICAL_PRODUCER_DIR / ANALYZER.PAIR_OUTPUT_NAME,
        "duplicates": CANONICAL_PRODUCER_DIR
        / ANALYZER.RAW_IDENTITY_DUPLICATE_OUTPUT_NAME,
        "oracle": CANONICAL_PRODUCER_DIR / ANALYZER.CASE_OUTPUT_NAME,
        "output": CANONICAL_PRODUCER_DIR / "release_receipt.json",
    }
    require_exact_path(args.preflight, CANONICAL_PREFLIGHT, "release preflight")
    require_exact_path(args.runner_source, RUNNER_SOURCE, "release runner")
    for label, expected in canonical_output_paths.items():
        require_exact_path(getattr(args, label), expected, f"release {label}")

    immutable_inputs = {
        "preflight": args.preflight,
        "producer_receipt": args.producer_receipt,
        "summary": args.summary,
        "sites": args.sites,
        "pairs": args.pairs,
        "duplicates": args.duplicates,
        "oracle": args.oracle,
    }
    for label, path in immutable_inputs.items():
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(f"Release {label} input is missing: {path}")
        require_mode_0444(path, f"release {label} input")
    require_mode_0444(Path(__file__).resolve(), "release finalizer source")
    require_mode_0444(args.runner_source.resolve(), "release runner source")
    for role, path in CANONICAL_CODE_PATHS.items():
        require_mode_0444(path.resolve(strict=True), f"cooccurrence {role} source")
    require_mode_0444(PREFLIGHT_SOURCE.resolve(strict=True), "preflight source")

    preflight = load_json(args.preflight, "preflight")
    receipt = load_json(args.producer_receipt, "producer receipt")
    summary = load_json(args.summary, "summary")
    checks: dict[str, bool] = {}
    preflight_checks = preflight.get("checks") or {}
    if not (
        preflight.get("schema_name") == ANALYZER.PREFLIGHT_SCHEMA_NAME
        and preflight.get("schema_version") == ANALYZER.PREFLIGHT_SCHEMA_VERSION
        and preflight.get("task_type") == "B_comprehensive_validation"
        and preflight.get("pass") is True
        and preflight.get("pass_semantics") == ANALYZER.PREFLIGHT_PASS_SEMANTICS
        and set(preflight_checks) == PREFLIGHT_CHECK_KEYS
        and all(value is True for value in preflight_checks.values())
    ):
        raise RuntimeError("Preflight schema/scope/check contract did not pass")
    require_equal(
        preflight.get("source_authority"),
        source_authority,
        "preflight source authority",
    )
    checks["preflight_schema_scope_pass"] = True
    expected_preflight_command = ANALYZER.canonical_preflight_command(
        manifest=CANONICAL_MANIFEST,
        assignments=CANONICAL_ASSIGNMENTS,
        sites=CANONICAL_SCREEN_SITES,
        intersubmod_root=CANONICAL_INTERSUBMOD_ROOT,
        independent_m2_audit=CANONICAL_INDEPENDENT_M2_AUDIT,
        primary_artifact_audit_pre=CANONICAL_PRIMARY_AUDIT,
        output=CANONICAL_PREFLIGHT,
    )
    require_equal(preflight.get("command"), expected_preflight_command, "preflight command")
    require_equal(
        preflight.get("parameters"),
        ANALYZER.PREFLIGHT_CANONICAL_PARAMETERS,
        "preflight parameters",
    )
    checks["preflight_command_and_parameters_exact"] = True

    if not (
        receipt.get("schema_name") == f"{ANALYZER.SCHEMA_NAME}.run_receipt"
        and receipt.get("schema_version") == ANALYZER.SCHEMA_VERSION
        and receipt.get("task_type") == "B_comprehensive_validation"
        and receipt.get("scope") == "all_manifest_samples"
        and receipt.get("full_scope") is True
        and receipt.get("pass") is True
        and receipt.get("release_status")
        == "PRODUCER_PASS_PENDING_RUNNER_RECONCILIATION"
        and receipt.get("pass_semantics")
        == "execution_integrity_only_pending_release_reconciliation"
    ):
        raise RuntimeError("Producer schema/scope/status contract did not pass")
    require_equal(
        receipt.get("source_authority"),
        source_authority,
        "producer source authority",
    )
    checks["producer_schema_scope_pass"] = True
    if not (
        summary.get("schema_name") == f"{ANALYZER.SCHEMA_NAME}.summary"
        and summary.get("schema_version") == ANALYZER.SCHEMA_VERSION
        and summary.get("task_type") == "B_comprehensive_validation"
        and summary.get("scope") == "all_manifest_samples"
        and summary.get("status") == "EXECUTION_PASS"
        and summary.get("pass") is True
    ):
        raise RuntimeError("Summary schema/scope/status contract did not pass")
    checks["summary_schema_scope_pass"] = True
    if not (
        preflight.get("raw_identity_release_contract")
        == receipt.get("raw_identity_release_contract")
        == summary.get("raw_identity_release_contract")
        == RAW_IDENTITY_RELEASE_CONTRACT
    ):
        raise RuntimeError("Raw-identity release contract drifted")
    checks["release_contract_exact"] = True

    preflight_inputs = preflight.get("inputs") or {}
    expected_preflight_inputs = {
        "manifest": CANONICAL_MANIFEST,
        "assignments": CANONICAL_ASSIGNMENTS,
        "sites": CANONICAL_SCREEN_SITES,
        "primary_artifact_audit_pre": CANONICAL_PRIMARY_AUDIT,
        "independent_m2_audit": CANONICAL_INDEPENDENT_M2_AUDIT,
    }
    require_equal(set(preflight_inputs), set(expected_preflight_inputs), "preflight input roles")
    for role, path in expected_preflight_inputs.items():
        require_artifact(preflight_inputs.get(role), path, f"preflight input {role}")
    preflight_input_lock = preflight.get("input_lock") or {}
    if not (
        preflight_input_lock.get("identity_before") == preflight_inputs
        and preflight_input_lock.get("identity_after") == preflight_inputs
        and preflight_input_lock.get("all_primary_inputs_unchanged") is True
    ):
        raise RuntimeError("Preflight input lock did not reconcile")
    checks["preflight_inputs_and_lock_reconciled"] = True

    independent = load_json(CANONICAL_INDEPENDENT_M2_AUDIT, "independent M2 audit")
    primary = load_json(CANONICAL_PRIMARY_AUDIT, "primary artifact audit")
    if not (
        independent.get("schema_name") == "intersubmod.independent_m2_gate_recount"
        and independent.get("schema_version") == "2.0.0"
        and independent.get("pass") is True
        and primary.get("schema_name") == "intersubmod.stable_primary_artifact_audit"
        and primary.get("schema_version") == "2.0.0"
        and primary.get("pass") is True
    ):
        raise RuntimeError("Independent M2 or primary audit schema/pass drifted")
    for audit_name, audit in (("independent M2", independent), ("primary", primary)):
        audit_inputs = audit.get("inputs") or {}
        for role, path in (
            ("stable_assignments", CANONICAL_ASSIGNMENTS),
            ("site_results", CANONICAL_SCREEN_SITES),
        ):
            record = audit_inputs.get(role)
            expected = artifact(path)
            if not isinstance(record, Mapping) or any(
                record.get(field) != value for field, value in expected.items()
            ):
                raise RuntimeError(f"{audit_name} {role} input identity drifted")
    checks["upstream_audits_bound_to_current_screen_inputs"] = True
    require_equal(
        (independent.get("inputs") or {}).get("claim_contract"),
        SOURCE_AUTHORITY.artifact(SOURCE_AUTHORITY.CLAIM_CONTRACT_PATH),
        "independent M2 canonical claim contract",
    )
    require_equal(
        independent.get("code"),
        {
            "independent_recount": artifact(
                SOURCE_AUTHORITY.EXPECTED_SOURCE_PATHS["independent_m2_auditor"]
            )
        },
        "independent M2 source authority binding",
    )
    require_equal(
        primary.get("source_authority"),
        source_authority,
        "primary audit source authority",
    )
    checks["upstream_audits_bound_to_source_authority"] = True

    producer_inputs = receipt.get("inputs") or {}
    expected_producer_artifact_inputs = {
        **expected_preflight_inputs,
        "preflight_receipt": CANONICAL_PREFLIGHT,
    }
    require_equal(
        set(producer_inputs),
        {*expected_producer_artifact_inputs, "intersubmod_runs"},
        "producer input roles",
    )
    for role, path in expected_producer_artifact_inputs.items():
        require_artifact(producer_inputs.get(role), path, f"producer input {role}")
    intersubmod_runs = producer_inputs.get("intersubmod_runs")
    if not isinstance(intersubmod_runs, Mapping) or set(intersubmod_runs) != set(
        CANONICAL_SAMPLES
    ):
        raise RuntimeError("Producer InterSubMod run-validation role set drifted")
    for sample, run in intersubmod_runs.items():
        if not isinstance(run, Mapping) or (run.get("validation") or {}).get("pass") is not True:
            raise RuntimeError(f"Producer InterSubMod validation did not pass: {sample}")
        receipt_record = run.get("receipt")
        if not isinstance(receipt_record, Mapping):
            raise RuntimeError(f"Producer InterSubMod receipt is malformed: {sample}")
        run_receipt_path = Path(str(receipt_record.get("path", "")))
        require_artifact(receipt_record, run_receipt_path, f"InterSubMod run {sample}")
    producer_locked_inputs = {
        role: producer_inputs[role] for role in expected_producer_artifact_inputs
    }
    producer_input_lock = receipt.get("input_lock") or {}
    if not (
        producer_input_lock.get("identity_before") == producer_locked_inputs
        and producer_input_lock.get("identity_after_compute") == producer_locked_inputs
        and producer_input_lock.get("identity_after_output") == producer_locked_inputs
        and producer_input_lock.get("all_primary_inputs_unchanged") is True
    ):
        raise RuntimeError("Producer input lock did not reconcile")
    checks["producer_inputs_and_lock_reconciled"] = True

    expected_producer_command = ANALYZER.canonical_producer_command(
        manifest=CANONICAL_MANIFEST,
        assignments=CANONICAL_ASSIGNMENTS,
        sites=CANONICAL_SCREEN_SITES,
        primary_artifact_audit_pre=CANONICAL_PRIMARY_AUDIT,
        independent_m2_audit=CANONICAL_INDEPENDENT_M2_AUDIT,
        preflight_receipt=CANONICAL_PREFLIGHT,
        intersubmod_root=CANONICAL_INTERSUBMOD_ROOT,
        output_dir=CANONICAL_PRODUCER_DIR,
    )
    require_equal(receipt.get("command"), expected_producer_command, "producer command")
    require_equal(
        receipt.get("parameters"),
        ANALYZER.PRODUCER_CANONICAL_PARAMETERS,
        "producer parameters",
    )
    checks["producer_command_and_parameters_exact"] = True

    producer_outputs = receipt.get("outputs") or {}
    require_equal(
        set(producer_outputs),
        {
            "summary_json",
            "site_tsv",
            "pair_tsv",
            "raw_identity_duplicate_audit_tsv",
            "case_json",
        },
        "producer output roles",
    )
    for role, path in (
        ("summary_json", args.summary),
        ("site_tsv", args.sites),
        ("pair_tsv", args.pairs),
        ("raw_identity_duplicate_audit_tsv", args.duplicates),
        ("case_json", args.oracle),
    ):
        require_artifact(producer_outputs.get(role), path, f"producer {role}")
    checks["producer_output_artifacts_reconciled"] = True

    preflight_runtime_lock = preflight.get("runtime_environment_lock") or {}
    producer_runtime_lock = receipt.get("runtime_environment_lock") or {}
    preflight_runtime = preflight_runtime_lock.get("identity_before")
    preflight_probe = (preflight.get("observed") or {}).get(
        "runtime_statistical_api_probe"
    )
    if not (
        preflight_runtime_lock.get("identity_after") == preflight_runtime
        and preflight_runtime_lock.get("direct_runtime_unchanged") is True
        and producer_runtime_lock.get("identity_before") == preflight_runtime
        and producer_runtime_lock.get("identity_after_compute") == preflight_runtime
        and producer_runtime_lock.get("identity_after_output") == preflight_runtime
        and producer_runtime_lock.get("statistical_api_probe_before") == preflight_probe
        and producer_runtime_lock.get("statistical_api_probe_after_compute")
        == preflight_probe
        and producer_runtime_lock.get("statistical_api_probe_after_output")
        == preflight_probe
        and producer_runtime_lock.get("direct_runtime_unchanged") is True
    ):
        raise RuntimeError("Runtime environment or statistical probe chain drifted")
    checks["runtime_environment_and_probe_chain_reconciled"] = True

    producer_python_cache = receipt.get("python_cache_lock") or {}
    if not (
        producer_python_cache.get("identity_before") == finalizer_python_cache
        and producer_python_cache.get("identity_after_compute")
        == finalizer_python_cache
        and producer_python_cache.get("identity_after_output")
        == finalizer_python_cache
        and producer_python_cache.get("canonical_cache_unchanged") is True
    ):
        raise RuntimeError("Producer/finalizer Python cache contract drifted")
    checks["python_cache_contract_reconciled"] = True

    producer_code = receipt.get("code") or {}
    require_equal(set(producer_code), set(CANONICAL_CODE_PATHS), "producer code roles")
    for role, path in CANONICAL_CODE_PATHS.items():
        require_artifact(producer_code.get(role), path, f"producer code {role}")
    source_lock = receipt.get("source_lock") or {}
    if not (
        source_lock.get("source_identity_before") == producer_code
        and source_lock.get("source_identity_after_compute") == producer_code
        and source_lock.get("source_identity_after_output") == producer_code
        and source_lock.get("all_sources_read_only_and_unchanged") is True
        and all(
            set((source_lock.get(field) or {})) == set(CANONICAL_CODE_PATHS)
            and set((source_lock.get(field) or {}).values()) == {"0o444"}
            for field in (
                "source_modes_before",
                "source_modes_after_compute",
                "source_modes_after_output",
            )
        )
    ):
        raise RuntimeError("Producer source lock did not reconcile")
    checks["producer_source_lock_reconciled"] = True
    preflight_code = preflight.get("code") or {}
    preflight_sources_before = preflight_code.get("source_identity_before") or {}
    preflight_sources = preflight_code.get("source_identity_after") or {}
    preflight_modes_before = preflight_code.get("source_modes_before") or {}
    preflight_modes_after = preflight_code.get("source_modes_after") or {}
    expected_preflight_source_roles = {"preflight", *CANONICAL_CODE_PATHS}
    if not (
        set(preflight_sources) == expected_preflight_source_roles
        and preflight_sources_before == preflight_sources
        and set(preflight_modes_before) == expected_preflight_source_roles
        and preflight_modes_before == preflight_modes_after
        and set(preflight_modes_before.values()) == {"0o444"}
        and preflight_sources.get("preflight") == artifact(PREFLIGHT_SOURCE)
        and all(
            preflight_sources.get(role) == producer_code.get(role)
            for role in CANONICAL_CODE_PATHS
        )
        and preflight_sources.get("release_finalizer") == artifact(Path(__file__).resolve())
        and preflight_sources.get("release_runner") == artifact(RUNNER_SOURCE)
    ):
        raise RuntimeError("Preflight source chain did not reconcile")
    checks["preflight_source_chain_includes_finalizer_and_runner"] = True

    recomputed = recompute_output_contract(
        sites_path=args.sites,
        pairs_path=args.pairs,
        duplicates_path=args.duplicates,
        oracle_path=args.oracle,
    )
    raw_preflight = (preflight.get("observed") or {}).get(
        "raw_bam_identity_recovery"
    ) or {}
    raw_summary = summary.get("raw_bam_identity_recovery_audit") or {}
    if not (
        raw_preflight.get("analysis_scope_identity_contract")
        == raw_summary.get("analysis_scope_identity_contract")
        == ANALYSIS_SCOPE_IDENTITY_CONTRACT
    ):
        raise RuntimeError("Analysis-scope identity contract drifted")
    checks["analysis_scope_contract_exact"] = True
    require_equal(
        recomputed["site_weighted_audit_sha256"],
        raw_preflight.get("site_weighted_audit_sha256"),
        "recomputed vs preflight site digest",
    )
    require_equal(
        recomputed["site_weighted_audit_sha256"],
        raw_summary.get("site_weighted_audit_sha256"),
        "recomputed vs summary site digest",
    )
    checks["site_weighted_analysis_digest_independently_recomputed"] = True

    summary_expectations = {
        "n_stable_sites": recomputed["stable_sites"],
        "n_focal_partner_pairs": recomputed["focal_partner_pairs"],
        "n_pair_by_confirmed": recomputed["pair_by_confirmed"],
        "n_pair_callability_gate_pass": recomputed["pair_callability_gate_pass"],
        "n_multi_marker_molecular_haplotype_base_candidates": recomputed[
            "multi_marker_molecular_haplotype_base_candidates"
        ],
    }
    for field, expected in summary_expectations.items():
        require_equal(summary.get(field), expected, f"summary {field}")
    raw_summary_expectations = {
        "n_site_tasks": recomputed["stable_sites"],
        "n_expected_projection_occurrences": recomputed[
            "raw_identity_expected_projection_occurrences"
        ],
        "n_matched_record_occurrences": recomputed[
            "raw_identity_matched_record_occurrences"
        ],
        "n_sites_with_collapsed_duplicates": recomputed[
            "raw_identity_sites_with_collapsed_duplicates"
        ],
        "n_duplicate_projection_occurrences_collapsed": recomputed[
            "raw_identity_duplicate_projection_occurrences_collapsed"
        ],
        "n_sparse_duplicate_rows": recomputed["raw_identity_sparse_duplicate_rows"],
        "n_duplicate_extra_record_occurrences_collapsed": recomputed[
            "raw_identity_duplicate_extra_record_occurrences_collapsed"
        ],
        "n_exact_duplicate_projection_occurrences_collapsed": recomputed[
            "raw_identity_exact_duplicate_projection_occurrences_collapsed"
        ],
        "n_rg_only_duplicate_projection_occurrences_collapsed": recomputed[
            "raw_identity_rg_only_duplicate_projection_occurrences_collapsed"
        ],
    }
    for field, expected in raw_summary_expectations.items():
        require_equal(raw_summary.get(field), expected, f"summary raw audit {field}")
    if not (
        raw_summary.get("all_site_results_passed_invariant_validation") is True
        and raw_summary.get("missing_projection_policy")
        == ANALYZER.RAW_IDENTITY_MISSING_POLICY
        and raw_summary.get("conflicting_analysis_payload_policy")
        == ANALYZER.RAW_IDENTITY_CONFLICT_POLICY
        and raw_summary.get("failure_counts_materialized") is False
    ):
        raise RuntimeError("Summary raw-identity fail-closed policy drifted")
    site_pair_summary = summary.get("site_pair_count_reconciliation") or {}
    require_equal(site_pair_summary.get("pass"), True, "summary site-pair pass")
    require_equal(
        site_pair_summary.get("n_sites_reconciled"),
        recomputed["stable_sites"],
        "summary reconciled sites",
    )
    require_equal(
        site_pair_summary.get("n_pair_rows_reconciled"),
        recomputed["focal_partner_pairs"],
        "summary reconciled pairs",
    )
    require_equal(
        site_pair_summary.get("n_m2_global_fdr_eligible_sites"),
        recomputed["m2_global_fdr_eligible_sites"],
        "summary M2 eligible sites",
    )
    require_equal(
        site_pair_summary.get("n_m2_exact_global_fdr_family_pairs"),
        recomputed["m2_exact_global_fdr_family_pairs"],
        "summary exact family pairs",
    )
    checks["summary_counts_independently_recomputed"] = True

    receipt_counts = receipt.get("counts") or {}
    receipt_expectations = {
        "stable_sites": recomputed["stable_sites"],
        "focal_partner_pairs": recomputed["focal_partner_pairs"],
        "m2_global_fdr_eligible_sites": recomputed["m2_global_fdr_eligible_sites"],
        "m2_exact_global_fdr_family_pairs": recomputed[
            "m2_exact_global_fdr_family_pairs"
        ],
        "pair_by_confirmed": recomputed["pair_by_confirmed"],
        "multi_marker_molecular_haplotype_base_candidates": recomputed[
            "multi_marker_molecular_haplotype_base_candidates"
        ],
        "raw_identity_expected_projection_occurrences": recomputed[
            "raw_identity_expected_projection_occurrences"
        ],
        "raw_identity_duplicate_projection_occurrences_collapsed": recomputed[
            "raw_identity_duplicate_projection_occurrences_collapsed"
        ],
        "raw_identity_duplicate_extra_record_occurrences_collapsed": recomputed[
            "raw_identity_duplicate_extra_record_occurrences_collapsed"
        ],
        "raw_identity_sparse_duplicate_rows": recomputed[
            "raw_identity_sparse_duplicate_rows"
        ],
    }
    for field, expected in receipt_expectations.items():
        require_equal(receipt_counts.get(field), expected, f"receipt count {field}")
    if not (
        receipt_counts.get("raw_identity_all_site_results_passed_invariant_validation")
        is True
        and receipt_counts.get("raw_identity_missing_projection_policy")
        == ANALYZER.RAW_IDENTITY_MISSING_POLICY
        and receipt_counts.get("raw_identity_conflicting_analysis_payload_policy")
        == ANALYZER.RAW_IDENTITY_CONFLICT_POLICY
        and receipt_counts.get("raw_identity_failure_counts_materialized") is False
    ):
        raise RuntimeError("Producer receipt fail-closed policy drifted")
    preflight_aggregate = raw_preflight.get("aggregate") or {}
    require_equal(
        preflight_aggregate.get("site_tasks"),
        recomputed["stable_sites"],
        "preflight site tasks",
    )
    require_equal(
        preflight_aggregate.get("expected_projection_occurrences"),
        recomputed["raw_identity_expected_projection_occurrences"],
        "preflight expected projections",
    )
    require_equal(
        preflight_aggregate.get("duplicate_projection_occurrences_collapsed"),
        recomputed["raw_identity_duplicate_projection_occurrences_collapsed"],
        "preflight duplicate projections",
    )
    require_equal(
        preflight_aggregate.get("duplicate_extra_record_occurrences_collapsed"),
        recomputed["raw_identity_duplicate_extra_record_occurrences_collapsed"],
        "preflight duplicate extras",
    )
    checks["all_full_scope_rows_and_counts_independently_recomputed"] = True

    try:
        preflight_finished = datetime.fromisoformat(str(preflight["finished_at_utc"]))
        producer_started = datetime.fromisoformat(str(receipt["started_at_utc"]))
        producer_finished = datetime.fromisoformat(str(receipt["finished_at_utc"]))
        release_started = datetime.fromisoformat(started_at_utc)
        if not preflight_finished <= producer_started <= producer_finished <= release_started:
            raise RuntimeError("Release chronology did not reconcile")
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError("Release chronology is malformed") from error
    checks["chronology_reconciled"] = True

    finished_at_utc = now_utc()
    if set(checks) != RELEASE_CHECK_KEYS or not all(checks.values()):
        raise RuntimeError("Release finalizer check-key contract drifted")
    code = {
        "release_finalizer": artifact(Path(__file__).resolve()),
        "release_runner": artifact(args.runner_source.resolve()),
    }
    payload = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "started_at_utc": started_at_utc,
        "finished_at_utc": finished_at_utc,
        "created_at_utc": finished_at_utc,
        "task_type": "B_comprehensive_validation",
        "scope": "all_manifest_samples",
        "raw_identity_release_contract": RAW_IDENTITY_RELEASE_CONTRACT,
        "producer_status": "EXECUTION_PASS",
        "release_status": "RELEASE_RECONCILIATION_PASS",
        "pass_semantics": "runner_reconciled_release_integrity_only_not_scientific_claim",
        "command": actual_command,
        "source_authority": source_authority,
        "python_cache_contract": finalizer_python_cache,
        "inputs": {name: artifact(path) for name, path in immutable_inputs.items()},
        "code": code,
        "source_modes": {name: "0o444" for name in code},
        "recomputed": recomputed,
        "checks": checks,
        "pass": True,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    args.output.chmod(0o444)
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "release_status": payload["release_status"],
                "pass": True,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
