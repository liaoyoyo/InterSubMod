#!/usr/bin/env python3
"""Fail-closed post-run gate for one M2 schema 1.2/2.0 pilot task.

The verifier never opens BAM/VCF/read-tag inputs and never launches a pilot.
It authenticates the already persisted extraction/ranking receipts and tables,
independently recomputes the primary unit/candidate gates, invokes the
independent profile-likelihood audit, and evaluates the resource thresholds
frozen in ``20260716_M2_schema1_2_2_0安全資源pilot規劃_01.md``.

Only an overall ``GO`` exits zero.  ``PROBE`` is deliberately non-passing for
full-run release, even though all evidence may be internally authentic.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import verify_full_m2_receipts as independent


DATASETS = independent.DATASETS
AUTOSOMES = independent.AUTOSOMES
EXTRACTION_DIRNAME = "extraction"
RANKING_DIRNAME = "ranking_bootstrap0"
UNIT_SOURCE_NAME = "m2_component_hp_rank_units.tsv.gz"
RESOURCE_GATE_PRODUCER = Path(__file__).resolve().parent / "run_full_m2_extraction.py"

EXPECTED_EXTRACTION_PARAMETERS = {
    "mapq_min": 20,
    "baseq_min": 20,
    "samtools_threads": 1,
    "excluded_flags_decimal": 0x4 | 0x100 | 0x200 | 0x400 | 0x800,
    "excluded_flags_hex": "0xf04",
    "bridge_thresholds": [1, 2, 3, 5],
    "component_linkage_bases": list(independent.LINKAGE_BASES),
    "primary_component_linkage_bases": list(independent.PRIMARY_BASES),
}
EXPECTED_RANKING_SCOPE = {
    "thresholds": [1, 2, 3, 5],
    "component_bases": list(independent.PRIMARY_BASES),
    "component_basis_mode": "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY",
    "hp_families": ["1", "2"],
    "phase_set_unit": "exact known PS for primary; missing PS excluded",
}
EXPECTED_RANKING_PARAMETERS = {
    "method_contract": independent.EXPECTED_METHOD_CONTRACT,
    "structural_exact_pattern_minread_grid": [1, 2, 3, 5],
    "primary_structural_exact_pattern_minread": 3,
    "scoring_minread": 1,
    "exact_k_max": 12,
    "max_vertex_sets": 256,
    "solver_time_limit_seconds_per_milp": 30.0,
    "minimum_bq_error_rate": 0.000001,
    "maximum_bq_error_rate": 0.25,
    "fixed_error_grid_conditional_binary_flip_probability": [0.005, 0.01, 0.02, 0.05],
    "conditional_candidate_ranking_bootstrap_replicates": 0,
    "conditional_candidate_ranking_bootstrap_base_seed": 20260716,
    "tie_tolerance_log_likelihood": 0.000001,
}


def expected_ranking_parameters(bootstrap_replicates: int) -> dict[str, Any]:
    """Return the frozen rank2 contract with the requested B0/B20 setting."""
    require(bootstrap_replicates >= 0, "expected bootstrap replicates must be nonnegative")
    return {
        **EXPECTED_RANKING_PARAMETERS,
        "conditional_candidate_ranking_bootstrap_replicates": bootstrap_replicates,
    }

REQUIRED_EXTRACTION_CHECKS = (
    "samtools_exit_zero",
    "raw_alignment_class_conserved",
    "raw_filter_funnel_conserved",
    "every_eligible_alignment_exact_sidecar_joined",
    "multi_region_alignment_unique",
    "canonical_molecule_unique",
    "eligible_alignment_rows_written",
    "site_call_reason_mass_conserved",
    "fixed_ra_mass_conserved",
    "alt_mass_conserved",
    "all_linkage_units_and_thresholds_conserve_their_active_sites",
    "primary_components_only_use_known_ps",
    "missing_ps_never_primary",
    "primary_membership_is_active_fixed_ra_only",
    "site_catalog_cardinality",
)
REQUIRED_RANKING_CHECKS = (
    "upstream_receipt_all_pass",
    "component_site_membership_conserved",
    "projection_funnel_conserved",
    "structural_scoring_separation_conserved",
    "all_units_have_semantic_hash",
    "all_patterns_have_semantic_hash",
    "k_gt_exact_limit_never_claimed_global_optimal",
    "all_informative_scoring_molecules_use_bq",
    "primary_units_are_exact_known_ps",
    "no_cross_ps_pattern_pooling",
    "known_ps_never_mixed",
    "missing_ps_separate_diagnostic",
    "partial_group_constraints_all_satisfied_when_evaluated",
    "effective_k_site_mass_conserved",
    "topology_partition_conserved",
    "k_route_partition_conserved",
    "all_candidates_have_semantic_hash",
    "all_converged_candidate_fits_have_global_kkt_certificate",
    "optimizer_abstain_candidates_are_not_claimed_winners",
    "all_responsibilities_have_semantic_hash",
    "candidate_relative_weights_conserve_per_unit",
    "posterior_responsibilities_conserve_per_pattern_candidate",
    "all_unit_evaluations_have_runtime_diagnostics",
    "runtime_diagnostics_are_finite_nonnegative_and_segmented",
    "runtime_scope_counts_match_unit_rows",
    "runtime_values_excluded_from_scientific_semantic_hashes",
    "runtime_invocation_flags_match_segment_contract",
)

PRIMARY_UNIT_REQUIRED_COLUMNS = (
    "dataset", "chrom", "threshold", "component_basis", "phase_set", "component_id",
    "family", "structural_exact_pattern_minread", "structural_minread_role",
    "k_observed_alt_active", "candidate_generation_complete",
    "candidate_generation_status", "distinct_vertex_sets_V", "selection_status",
    "top_vertex_set_ids",
    "conditional_candidate_ranking_bootstrap_replicates",
    "conditional_candidate_ranking_bootstrap_status",
    "conditional_candidate_ranking_bootstrap_seed",
    "conditional_candidate_ranking_bootstrap_top_vertex_set_frequency",
    "conditional_candidate_ranking_bootstrap_primary_top_set_frequency",
    "conditional_candidate_ranking_bootstrap_all_converged",
)
CANDIDATE_GATE_REQUIRED_COLUMNS = (
    "dataset", "chrom", "threshold", "component_basis", "phase_set", "component_id",
    "family", "structural_exact_pattern_minread", "fit_converged", "fit_monotone",
    "optimizer_status", "refinement_iterations", "global_log_likelihood_gap_bound",
    "simplex_residual", "is_winner", "is_tied_winner",
)

GLOBAL_LL_GAP_TOLERANCE = 1e-8
SIMPLEX_RESIDUAL_TOLERANCE = 1e-12
FOUR_HOURS_SECONDS = 4 * 60 * 60
EIGHT_HOURS_SECONDS = 8 * 60 * 60
FOUR_GIB_KIB = 4 * 1024 * 1024
EIGHT_GIB_KIB = 8 * 1024 * 1024


class PilotVerificationError(RuntimeError):
    """Raised when persisted pilot evidence violates its hard contract."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PilotVerificationError(message)


def _strict_bool(value: Any, label: str) -> bool:
    if value not in {"true", "false"}:
        raise PilotVerificationError(f"invalid boolean {label}: {value!r}")
    return value == "true"


def _require_columns(fieldnames: Sequence[str] | None, required: Iterable[str], label: str) -> None:
    fields = set(fieldnames or ())
    missing = sorted(set(required) - fields)
    require(not missing, f"{label} missing required columns: {missing}")


def _require_subset(actual: Mapping[str, Any], expected: Mapping[str, Any], label: str) -> None:
    for key, expected_value in expected.items():
        require(key in actual, f"{label} missing key: {key}")
        require(actual[key] == expected_value, f"{label} mismatch at {key}: {actual[key]!r}")


def _verify_named_checks(receipt: Mapping[str, Any], required: Sequence[str], label: str) -> None:
    checks = receipt.get("checks") or {}
    missing = sorted(set(required) - set(checks))
    failed = sorted(key for key, value in checks.items() if value is not True)
    require(not missing, f"{label} missing named checks: {missing}")
    require(not failed, f"{label} has failed/non-boolean checks: {failed}")


def _verify_source_identity(identity: Mapping[str, Any], label: str) -> dict[str, Any]:
    try:
        path = Path(str(identity["path"])).resolve(strict=True)
        declared_sha = str(identity["sha256"])
    except (KeyError, OSError, TypeError) as exc:
        raise PilotVerificationError(f"malformed/unavailable {label} identity") from exc
    actual_sha = independent.sha256_path(path)
    require(actual_sha == declared_sha, f"{label} source SHA mismatch: {path}")
    return {"path": str(path), "size_bytes": path.stat().st_size, "sha256": actual_sha}


def verify_pilot_resource_gate(
    path: Path,
    *,
    pilot_root: Path,
    output_root: Path,
    dataset: str,
    chrom: str,
    stage: str,
    gate_label: str,
) -> dict[str, Any]:
    """Independently authenticate one frozen-runner pilot process/disk gate."""
    expected_path = pilot_root / "resource_gates" / f"{gate_label}.json"
    try:
        require(
            path.resolve(strict=True) == expected_path.resolve(strict=True),
            f"pilot {gate_label} resource gate path swap",
        )
        gate, authenticated = independent._authenticate_json_with_sidecar(
            expected_path, f"pilot {gate_label} resource gate", immutable=True
        )
    except independent.VerificationError as exc:
        raise PilotVerificationError(str(exc)) from exc
    independent._exact_keys(
        gate, independent.RESOURCE_GATE_RECEIPT_KEYS,
        f"pilot {gate_label} resource gate",
    )
    require(
        gate["schema_name"] == independent.RESOURCE_GATE_SCHEMA_NAME
        and gate["schema_version"] == independent.RESOURCE_GATE_SCHEMA_VERSION
        and gate["stage"] == stage
        and gate["gate_scope"] == "pilot",
        f"pilot {gate_label} resource gate schema/stage/scope mismatch",
    )
    expected_gate_id = independent.semantic_json_sha256({
        key: value for key, value in gate.items()
        if key not in {"gate_id", "receipt_integrity"}
    })
    require(gate["gate_id"] == expected_gate_id, f"pilot {gate_label} gate_id mismatch")
    target = gate.get("target") or {}
    require(
        set(target) == {
            "task_type", "dataset", "chrom", "gate_label", "output_root",
            "manifest", "release_manifest",
        },
        f"pilot {gate_label} target exact-key mismatch",
    )
    require(
        target["task_type"] == "B_COMPREHENSIVE_VALIDATION_RELEASE_PILOT"
        and target["dataset"] == dataset
        and target["chrom"] == chrom
        and target["gate_label"] == gate_label,
        f"pilot {gate_label} target task/scope mismatch",
    )
    root = output_root.resolve(strict=True)
    observed_root = root.stat()
    require(
        target["output_root"] == {
            "path": str(root), "st_dev": int(observed_root.st_dev),
            "st_ino": int(observed_root.st_ino),
        },
        f"pilot {gate_label} output-root identity mismatch",
    )
    manifest = target.get("manifest") or {}
    try:
        manifest_path = Path(str(manifest["path"])).resolve(strict=True)
        manifest_sha = independent.sha256_path(manifest_path)
    except (KeyError, OSError, TypeError) as exc:
        raise PilotVerificationError(f"pilot {gate_label} manifest identity malformed") from exc
    require(
        manifest == {"path": str(manifest_path), "sha256": manifest_sha},
        f"pilot {gate_label} manifest identity mismatch",
    )
    canonical_manifest_identity = {
        "path": str(manifest_path),
        "size_bytes": int(manifest_path.stat().st_size),
        "sha256": manifest_sha,
    }
    release_identity = target.get("release_manifest")
    require(isinstance(release_identity, dict), f"pilot {gate_label} has no frozen release binding")
    try:
        release_path = Path(str(release_identity.get("path", ""))).resolve(strict=True)
        release_document, release_authenticated = independent._authenticate_json_with_sidecar(
            release_path, f"pilot {gate_label} release manifest", immutable=True
        )
    except independent.VerificationError as exc:
        raise PilotVerificationError(str(exc)) from exc
    release_sidecar = release_path.with_name(f"{release_path.name}.sha256")
    canonical_release_identity = {
        "path": release_authenticated["path"],
        "sha256": release_authenticated["sha256"],
        "semantic_sha256": release_authenticated["semantic_sha256"],
        "sidecar": {
            "path": str(release_sidecar.resolve(strict=True)),
            "sha256": release_authenticated["sidecar_sha256"],
        },
    }
    require(
        release_identity == canonical_release_identity,
        f"pilot {gate_label} complete release manifest identity mismatch",
    )
    require(
        release_document.get("schema_name") == independent.RELEASE_SCHEMA_NAME
        and release_document.get("schema_version") == independent.RELEASE_SCHEMA_VERSION
        and release_document.get("authority_mode") == "CANONICAL_V5_FROZEN"
        and release_document.get("validation_evidence_eligible") is True
        and release_document.get("all_pass") is True,
        f"pilot {gate_label} release is not a passing frozen canonical binding",
    )
    frozen_manifest_copy = (
        (release_document.get("canonical_manifest") or {}).get("immutable_copy") or {}
    )
    frozen_manifest_relative = Path(str(frozen_manifest_copy.get("path", "")))
    require(
        not frozen_manifest_relative.is_absolute()
        and ".." not in frozen_manifest_relative.parts
        and frozen_manifest_relative.parts,
        f"pilot {gate_label} frozen canonical manifest path is invalid",
    )
    try:
        frozen_manifest_path = (
            release_path.parent / frozen_manifest_relative
        ).resolve(strict=True)
    except OSError as exc:
        raise PilotVerificationError(
            f"pilot {gate_label} frozen canonical manifest is unavailable"
        ) from exc
    require(
        frozen_manifest_path == manifest_path
        and frozen_manifest_copy.get("sha256") == manifest_sha,
        f"pilot {gate_label} manifest is not the canonical copy in its frozen release binding",
    )
    producer = _verify_source_identity(
        gate.get("producer_source") or {}, f"pilot {gate_label} resource-gate producer"
    )
    expected_producer = RESOURCE_GATE_PRODUCER.resolve(strict=True)
    require(
        producer["path"] == str(expected_producer)
        and producer["sha256"] == independent.sha256_path(expected_producer),
        f"pilot {gate_label} gate was not produced by the adjacent frozen runner",
    )
    process = gate.get("process_snapshot") or {}
    process_core = {
        "process_count": process.get("process_count"),
        "root_count": process.get("root_count"),
        "representatives": process.get("representatives"),
    }
    require(
        process_core == {"process_count": 0, "root_count": 0, "representatives": []}
        and process.get("semantic_sha256") == independent.semantic_json_sha256(process_core),
        f"pilot {gate_label} is not an authenticated zero-conflict snapshot",
    )
    filesystem = gate.get("filesystem_snapshot") or {}
    filesystem_core = {
        key: filesystem.get(key)
        for key in (
            "probe_path", "target_output_root", "st_dev", "f_bavail", "f_frsize",
            "available_bytes", "required_reserve_bytes", "disk_pass",
        )
    }
    require(
        isinstance(filesystem_core["f_bavail"], int)
        and isinstance(filesystem_core["f_frsize"], int)
        and filesystem_core["available_bytes"]
        == filesystem_core["f_bavail"] * filesystem_core["f_frsize"]
        and filesystem_core["required_reserve_bytes"]
        == independent.RESOURCE_GATE_REQUIRED_RESERVE_BYTES
        and filesystem_core["available_bytes"]
        >= independent.RESOURCE_GATE_REQUIRED_RESERVE_BYTES
        and filesystem_core["disk_pass"] is True
        and filesystem_core["probe_path"] == str(root)
        and filesystem_core["target_output_root"] == str(root)
        and filesystem_core["st_dev"] == int(observed_root.st_dev)
        and filesystem.get("semantic_sha256")
        == independent.semantic_json_sha256(filesystem_core),
        f"pilot {gate_label} 300 GiB filesystem attestation mismatch",
    )
    require(
        gate.get("checks")
        == {"zero_conflict": True, "disk_reserve": True, "all_pass": True},
        f"pilot {gate_label} resource gate checks are not all PASS",
    )
    return {
        "path": authenticated["path"],
        "sha256": authenticated["sha256"],
        "semantic_sha256": independent.semantic_json_sha256(gate),
        "gate_id": gate["gate_id"],
        "stage": stage,
        "gate_label": gate_label,
        "producer": producer,
        "available_bytes": filesystem_core["available_bytes"],
        "required_reserve_bytes": filesystem_core["required_reserve_bytes"],
        "canonical_manifest_identity": canonical_manifest_identity,
        "release_manifest_identity": canonical_release_identity,
    }


def verify_extraction(
    extraction_dir: Path, dataset: str, chrom: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    receipt_path = extraction_dir / "receipt.json"
    try:
        receipt, receipt_sha = independent.authenticate_receipt(
            receipt_path,
            "intersubmod.lossless_read_linkage_chromosome_receipt",
            "1.2.0",
        )
    except independent.VerificationError as exc:
        raise PilotVerificationError(str(exc)) from exc
    scope = receipt.get("scope") or {}
    require(scope.get("dataset") == dataset and scope.get("chrom") == chrom, "extraction scope mismatch")
    require(
        isinstance(scope.get("n_sSNV"), int) and not isinstance(scope.get("n_sSNV"), bool)
        and int(scope["n_sSNV"]) > 0,
        "extraction n_sSNV must be a positive integer",
    )
    _require_subset(receipt.get("parameters") or {}, EXPECTED_EXTRACTION_PARAMETERS, "extraction parameters")
    _verify_named_checks(receipt, REQUIRED_EXTRACTION_CHECKS, "extraction")
    try:
        count_checks = independent.extraction_count_checks(
            receipt.get("counts") or {}, f"{dataset}/{chrom}"
        )
    except independent.VerificationError as exc:
        raise PilotVerificationError(str(exc)) from exc
    outputs = receipt.get("outputs") or {}
    required_matches = {
        suffix: [name for name in outputs if name.endswith(suffix)]
        for suffix in independent.REQUIRED_EXTRACTION_SUFFIXES
    }
    require(
        all(len(matches) == 1 for matches in required_matches.values()),
        f"required extraction outputs missing/duplicated: {required_matches}",
    )
    verified_outputs = {}
    for name, identity in outputs.items():
        try:
            verified_outputs[name] = independent.verify_identity(identity, extraction_dir)
        except independent.VerificationError as exc:
            raise PilotVerificationError(str(exc)) from exc
    provenance = receipt.get("provenance") or {}
    verified_provenance = {
        key: _verify_source_identity(provenance.get(key) or {}, f"extraction provenance/{key}")
        for key in ("manifest", "extractor")
    }
    return receipt, {
        "receipt_path": str(receipt_path.resolve()),
        "receipt_sha256": receipt_sha,
        "scope": scope,
        "parameters": {key: (receipt.get("parameters") or {})[key] for key in EXPECTED_EXTRACTION_PARAMETERS},
        "independent_count_checks": count_checks,
        "verified_provenance": verified_provenance,
        "n_outputs_verified": len(verified_outputs),
        "outputs": verified_outputs,
    }


def verify_ranking_contract(
    ranking_dir: Path,
    extraction_dir: Path,
    extraction_receipt: Mapping[str, Any],
    dataset: str,
    chrom: str,
    expected_bootstrap_replicates: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    receipt_path = ranking_dir / "receipt.json"
    try:
        receipt, receipt_sha = independent.authenticate_receipt(
            receipt_path,
            "intersubmod.m2_symbolic_patterns_vertex_rank_receipt",
            "2.0.0",
        )
    except independent.VerificationError as exc:
        raise PilotVerificationError(str(exc)) from exc
    scope = receipt.get("scope") or {}
    require(scope.get("dataset") == [dataset] and scope.get("chrom") == [chrom], "ranking scope mismatch")
    _require_subset(scope, EXPECTED_RANKING_SCOPE, "ranking scope")
    parameters = receipt.get("parameters") or {}
    expected_parameters = expected_ranking_parameters(expected_bootstrap_replicates)
    _require_subset(parameters, expected_parameters, "ranking parameters")
    optimizer_contract = parameters.get("optimizer_contract") or {}
    require(
        optimizer_contract.get("global_log_likelihood_gap_tolerance") == GLOBAL_LL_GAP_TOLERANCE,
        "ranking optimizer global-gap tolerance mismatch",
    )
    require(
        optimizer_contract.get("simplex_residual_tolerance") == SIMPLEX_RESIDUAL_TOLERANCE,
        "ranking optimizer simplex tolerance mismatch",
    )
    require(optimizer_contract.get("same_read_vaf_added") is False, "same-read VAF was added")
    _verify_named_checks(receipt, REQUIRED_RANKING_CHECKS, "ranking")

    provenance = receipt.get("provenance") or {}
    ranker_identity = provenance.get("ranker") or {}
    verified_ranker = _verify_source_identity(ranker_identity, "ranking provenance/ranker")
    try:
        independent.verify_child_method_contract_and_ranker_source(
            receipt,
            Path(verified_ranker["path"]),
            verified_ranker["sha256"],
            f"{dataset}/{chrom}",
        )
    except independent.VerificationError as exc:
        raise PilotVerificationError(str(exc)) from exc

    extraction_receipt_path = (extraction_dir / "receipt.json").resolve()
    extraction_receipt_sha = independent.sha256_path(extraction_receipt_path)
    upstream_provenance = provenance.get("upstream_extraction_receipt") or {}
    require(
        Path(str(upstream_provenance.get("path", ""))).resolve() == extraction_receipt_path
        and upstream_provenance.get("sha256") == extraction_receipt_sha,
        "ranking provenance is not bound to the extraction receipt",
    )
    upstream = receipt.get("upstream") or {}
    require(
        Path(str(upstream.get("path", ""))).resolve() == extraction_receipt_path
        and upstream.get("sha256") == extraction_receipt_sha,
        "ranking upstream block is not bound to the extraction receipt",
    )

    extraction_outputs = extraction_receipt.get("outputs") or {}
    expected_input_names = {
        name for name in extraction_outputs if name.endswith(independent.REQUIRED_EXTRACTION_SUFFIXES)
    }
    input_files = receipt.get("input_files") or {}
    require(set(input_files) == expected_input_names, "ranking input file set differs from extraction receipt")
    for name, identity in input_files.items():
        require(identity == extraction_outputs[name], f"ranking input identity differs: {name}")
        try:
            independent.verify_identity(identity, extraction_dir)
        except independent.VerificationError as exc:
            raise PilotVerificationError(str(exc)) from exc

    outputs = receipt.get("outputs") or {}
    for required_name in (
        independent.PATTERN_SOURCE_NAME,
        UNIT_SOURCE_NAME,
        independent.CANDIDATE_SOURCE_NAME,
        independent.RUNTIME_SOURCE_NAME,
    ):
        require(required_name in outputs, f"ranking output missing: {required_name}")
    verified_outputs = {}
    for name, identity in outputs.items():
        try:
            verified_outputs[name] = independent.verify_identity(identity, ranking_dir)
        except independent.VerificationError as exc:
            raise PilotVerificationError(str(exc)) from exc
    try:
        runtime_audit = independent.verify_child_runtime_diagnostics(receipt, dataset, chrom)
    except independent.VerificationError as exc:
        raise PilotVerificationError(str(exc)) from exc
    return receipt, {
        "receipt_path": str(receipt_path.resolve()),
        "receipt_sha256": receipt_sha,
        "scope": scope,
        "parameters": {key: parameters[key] for key in expected_parameters},
        "ranker": verified_ranker,
        "upstream_extraction_receipt_sha256": extraction_receipt_sha,
        "n_outputs_verified": len(verified_outputs),
        "outputs": verified_outputs,
        "runtime_diagnostics_primary_rows": len(
            runtime_audit["all_primary"]["unit_total_elapsed_seconds"]
        ),
    }


def _unit_key(row: Mapping[str, str]) -> tuple[str, ...]:
    return tuple(
        str(row[key])
        for key in (
            "dataset", "chrom", "component_basis", "phase_set", "threshold",
            "component_id", "family", "structural_exact_pattern_minread",
        )
    )


def _audit_bootstrap_row(
    row: Mapping[str, str],
    expected_replicates: int,
    expected_base_seed: int,
    generation_complete: bool,
) -> str:
    status = row["conditional_candidate_ranking_bootstrap_status"]
    replicates = int(row["conditional_candidate_ranking_bootstrap_replicates"])
    seed_raw = row["conditional_candidate_ranking_bootstrap_seed"]
    frequency_raw = row["conditional_candidate_ranking_bootstrap_top_vertex_set_frequency"]
    primary_raw = row["conditional_candidate_ranking_bootstrap_primary_top_set_frequency"]
    converged_raw = row["conditional_candidate_ranking_bootstrap_all_converged"]
    try:
        frequencies = json.loads(frequency_raw)
    except json.JSONDecodeError as exc:
        raise PilotVerificationError("invalid bootstrap top-vertex frequency JSON") from exc
    require(isinstance(frequencies, dict), "bootstrap top-vertex frequency is not an object")
    for key, value in frequencies.items():
        require(
            isinstance(key, str) and isinstance(value, (int, float))
            and not isinstance(value, bool) and math.isfinite(float(value))
            and 0.0 <= float(value) <= 1.0,
            "invalid bootstrap top-vertex frequency",
        )

    vertex_count_raw = row["distinct_vertex_sets_V"]
    vertex_count = int(vertex_count_raw) if vertex_count_raw else None
    if not generation_complete:
        require(
            status == "NOT_APPLICABLE_NO_COMPLETE_MULTI_VERTEX_CANDIDATE"
            and replicates == 0 and seed_raw == "" and frequencies == {}
            and primary_raw == "" and converged_raw == "",
            "non-complete unit violates bootstrap not-applicable contract",
        )
        return status
    require(vertex_count is not None and vertex_count >= 1, "complete unit lacks vertex-set count")
    if vertex_count == 1:
        require(
            status == "NOT_APPLICABLE_V1" and replicates == 0 and seed_raw == ""
            and frequencies == {} and primary_raw == "" and converged_raw == "",
            "V=1 unit violates bootstrap not-applicable contract",
        )
        return status
    if status == "NOT_RUN_PRIMARY_NONCONVERGENCE":
        require(
            replicates == 0 and seed_raw == "" and frequencies == {}
            and primary_raw == "" and converged_raw == "",
            "primary-nonconvergence bootstrap fields are inconsistent",
        )
        return status
    require(seed_raw != "" and 0 <= int(seed_raw) < 2**32, "bootstrap seed is missing/out of range")
    seed_payload = (
        f"{row['dataset']}\0{row['chrom']}\0{row['component_basis']}\0"
        f"{row['threshold']}\0{row['component_id']}\0{row['family']}"
    )
    independently_derived_seed = (
        int(hashlib.sha256(seed_payload.encode()).hexdigest()[:16], 16)
        ^ expected_base_seed
    ) % (2**32)
    require(int(seed_raw) == independently_derived_seed, "bootstrap derived seed mismatch")
    if expected_replicates == 0:
        require(
            status == "NOT_RUN_REQUESTED_0" and replicates == 0 and frequencies == {}
            and primary_raw == "" and converged_raw == "",
            "B0 unit violates bootstrap not-run contract",
        )
        return status
    require(
        status in {"COMPLETE", "ABSTAIN_NONCONVERGENCE"}
        and replicates == expected_replicates and frequencies,
        "requested bootstrap unit has wrong status/replicate/frequency fields",
    )
    try:
        primary_frequency = float(primary_raw)
    except ValueError as exc:
        raise PilotVerificationError("bootstrap primary frequency is invalid") from exc
    require(
        math.isfinite(primary_frequency) and 0.0 <= primary_frequency <= 1.0,
        "bootstrap primary frequency is outside [0,1]",
    )
    expected_converged = status == "COMPLETE"
    require(
        _strict_bool(converged_raw, "conditional bootstrap all_converged")
        == expected_converged,
        "bootstrap status and all_converged disagree",
    )
    return status


def audit_primary_units(
    unit_path: Path,
    dataset: str,
    chrom: str,
    primary_minread: int,
    exact_k_max: int,
    expected_bootstrap_replicates: int,
    expected_bootstrap_base_seed: int,
) -> dict[str, Any]:
    route_counts: Counter[str] = Counter()
    generation_status_counts: Counter[str] = Counter()
    incomplete_keys: set[tuple[str, ...]] = set()
    incomplete_not_abstain: list[tuple[str, ...]] = []
    bootstrap_status_counts: Counter[str] = Counter()
    n_primary = 0
    with gzip.open(unit_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_columns(reader.fieldnames, PRIMARY_UNIT_REQUIRED_COLUMNS, "unit table")
        for row in reader:
            if row["structural_minread_role"] != "PRIMARY":
                continue
            n_primary += 1
            require(row["dataset"] == dataset and row["chrom"] == chrom, "primary unit scope mismatch")
            require(
                row["component_basis"] in independent.PRIMARY_BASES
                and row["phase_set"] not in {"", ".", "None", "NA", "UNKNOWN"},
                "non-primary or missing-PS unit in primary rows",
            )
            require(int(row["structural_exact_pattern_minread"]) == primary_minread, "primary minread row mismatch")
            effective_k = int(row["k_observed_alt_active"])
            complete = _strict_bool(row["candidate_generation_complete"], "candidate_generation_complete")
            status = row["candidate_generation_status"]
            bootstrap_status_counts[_audit_bootstrap_row(
                row,
                expected_bootstrap_replicates,
                expected_bootstrap_base_seed,
                complete,
            )] += 1
            generation_status_counts[status] += 1
            if status == "EXACT_OPTIMAL_VERTEX_SETS_COMPLETE":
                require(complete and effective_k <= exact_k_max, "exact-complete unit violates exact-k contract")
                route_counts["EXACT_COMPLETE"] += 1
            elif status.startswith("EXACT_ENUMERATION_INCOMPLETE"):
                require(not complete and effective_k <= exact_k_max, "exact-incomplete unit violates exact-k contract")
                route_counts["EXACT_INCOMPLETE"] += 1
                key = _unit_key(row)
                incomplete_keys.add(key)
                if "ABSTAIN" not in row["selection_status"] or row["top_vertex_set_ids"] != "[]":
                    incomplete_not_abstain.append(key)
            elif status == "K_GT_EXACT_LIMIT":
                require(not complete and effective_k > exact_k_max, "k>limit unit violates route contract")
                require("ABSTAIN" in row["selection_status"], "k>limit unit is not abstained")
                route_counts["GT_EXACT_LIMIT"] += 1
            elif status in {"NO_INFORMATIVE_PATTERN", "NO_STRUCTURAL_ALT_PATTERN_AT_DECLARED_MINREAD"}:
                require(not complete, "no-structure unit marked complete")
                route_counts["NOT_RUN_NO_STRUCTURE"] += 1
            else:
                raise PilotVerificationError(f"unknown candidate generation status: {status}")
    require(n_primary > 0, "pilot has no primary unit rows")
    exact_eligible = route_counts["EXACT_COMPLETE"] + route_counts["EXACT_INCOMPLETE"]
    inference_units = exact_eligible + route_counts["GT_EXACT_LIMIT"]
    incomplete_rate = (
        route_counts["EXACT_INCOMPLETE"] / exact_eligible if exact_eligible else None
    )
    exact_limit_coverage = exact_eligible / inference_units if inference_units else None
    route_names = ("EXACT_COMPLETE", "EXACT_INCOMPLETE", "NOT_RUN_NO_STRUCTURE", "GT_EXACT_LIMIT")
    return {
        "n_primary_units": n_primary,
        "route_counts": {key: int(route_counts[key]) for key in route_names},
        "candidate_generation_status_counts": dict(generation_status_counts),
        "expected_bootstrap_replicates": expected_bootstrap_replicates,
        "bootstrap_status_counts": dict(sorted(bootstrap_status_counts.items())),
        "bootstrap_evaluated_units": (
            bootstrap_status_counts["COMPLETE"]
            + bootstrap_status_counts["ABSTAIN_NONCONVERGENCE"]
        ),
        "exact_eligible_units": exact_eligible,
        "inference_required_units": inference_units,
        "exact_enumeration_incomplete_rate": incomplete_rate,
        "exact_limit_coverage": exact_limit_coverage,
        "n_incomplete_not_abstain": len(incomplete_not_abstain),
        "incomplete_unit_keys": incomplete_keys,
        "definitions": {
            "exact_enumeration_incomplete_rate": "EXACT_INCOMPLETE / (EXACT_COMPLETE + EXACT_INCOMPLETE)",
            "exact_limit_coverage": "(EXACT_COMPLETE + EXACT_INCOMPLETE) / ((EXACT_COMPLETE + EXACT_INCOMPLETE) + GT_EXACT_LIMIT)",
            "no_structure_units": "excluded from coverage denominator because no structural inference was requested",
        },
    }


def _nearest_rank(values: Sequence[int], probability: float) -> int | None:
    if not values:
        return None
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, math.ceil(probability * len(ordered)) - 1)]


def audit_candidate_certificates(
    candidate_path: Path,
    dataset: str,
    chrom: str,
    primary_minread: int,
    incomplete_unit_keys: set[tuple[str, ...]],
) -> dict[str, Any]:
    n_candidates = n_certified = n_winners = n_uncertified_winners = 0
    n_primary_candidates = 0
    max_gap = max_simplex = 0.0
    refinements: list[int] = []
    candidates_in_incomplete_units = 0
    with gzip.open(candidate_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_columns(reader.fieldnames, CANDIDATE_GATE_REQUIRED_COLUMNS, "candidate table")
        for row in reader:
            n_candidates += 1
            require(row["dataset"] == dataset and row["chrom"] == chrom, "candidate scope mismatch")
            minread = int(row["structural_exact_pattern_minread"])
            require(minread in EXPECTED_RANKING_PARAMETERS["structural_exact_pattern_minread_grid"], "candidate minread outside grid")
            converged = _strict_bool(row["fit_converged"], "fit_converged")
            monotone = _strict_bool(row["fit_monotone"], "fit_monotone")
            winner = _strict_bool(row["is_winner"], "is_winner")
            _strict_bool(row["is_tied_winner"], "is_tied_winner")
            gap = float(row["global_log_likelihood_gap_bound"])
            simplex = float(row["simplex_residual"])
            iterations = int(row["refinement_iterations"])
            require(
                math.isfinite(gap) and gap >= 0.0 and math.isfinite(simplex)
                and simplex >= 0.0 and iterations >= 0,
                "candidate certificate contains invalid numeric value",
            )
            certified = (
                converged and monotone and row["optimizer_status"].startswith("CERTIFIED_")
                and gap <= GLOBAL_LL_GAP_TOLERANCE
                and simplex <= SIMPLEX_RESIDUAL_TOLERANCE
            )
            n_certified += int(certified)
            n_winners += int(winner)
            n_uncertified_winners += int(winner and not certified)
            max_gap = max(max_gap, gap)
            max_simplex = max(max_simplex, simplex)
            refinements.append(iterations)
            if minread == primary_minread:
                n_primary_candidates += 1
                candidates_in_incomplete_units += int(_unit_key(row) in incomplete_unit_keys)
    return {
        "n_candidates": n_candidates,
        "n_primary_candidates": n_primary_candidates,
        "n_globally_certified_candidates": n_certified,
        "globally_certified_fraction": n_certified / n_candidates if n_candidates else None,
        "n_winners": n_winners,
        "n_uncertified_winners": n_uncertified_winners,
        "n_candidates_in_primary_incomplete_units": candidates_in_incomplete_units,
        "max_global_log_likelihood_gap_bound": max_gap if n_candidates else None,
        "max_simplex_residual": max_simplex if n_candidates else None,
        "refinement_iterations": {
            "n": len(refinements),
            "p99_nearest_rank": _nearest_rank(refinements, 0.99),
            "max": max(refinements) if refinements else None,
        },
    }


def compare_table_audits_to_receipt(
    ranking_receipt: Mapping[str, Any], unit_audit: Mapping[str, Any], candidate_audit: Mapping[str, Any]
) -> dict[str, bool]:
    aggregate = ranking_receipt.get("aggregate") or {}
    diagnostics = ((ranking_receipt.get("candidate_evidence_counts") or {}).get("optimizer_diagnostics") or {})
    route_names = ("EXACT_COMPLETE", "EXACT_INCOMPLETE", "NOT_RUN_NO_STRUCTURE", "GT_EXACT_LIMIT")
    route_counts = {
        key: int((unit_audit.get("route_counts") or {}).get(key, 0)) for key in route_names
    }
    expected_routes = {
        key: int((aggregate.get("k_route_counts") or {}).get(key, 0))
        for key in route_names
    }
    checks = {
        "primary_unit_count_matches_receipt": int(unit_audit["n_primary_units"])
        == int(aggregate.get("n_component_hp_units", -1)),
        "k_route_counts_match_receipt": route_counts == expected_routes,
        "bootstrap_status_counts_match_receipt": (
            unit_audit.get("bootstrap_status_counts") or {}
        ) == (aggregate.get("conditional_candidate_ranking_bootstrap_status_counts") or {}),
        "bootstrap_partition_counts_match_receipt": (
            int((unit_audit.get("bootstrap_status_counts") or {}).get("COMPLETE", 0))
            == int(aggregate.get("conditional_candidate_ranking_bootstrap_complete_units", -1))
            and sum(
                int(value) for key, value in (unit_audit.get("bootstrap_status_counts") or {}).items()
                if str(key).startswith("NOT_")
            ) == int(aggregate.get("conditional_candidate_ranking_bootstrap_not_run_units", -1))
            and int((unit_audit.get("bootstrap_status_counts") or {}).get("ABSTAIN_NONCONVERGENCE", 0))
            == int(aggregate.get("conditional_candidate_ranking_bootstrap_nonconverged_units", -1))
        ),
        "candidate_fit_count_matches_receipt": int(candidate_audit["n_candidates"])
        == int(diagnostics.get("candidate_fits", -1)),
        "certified_candidate_count_matches_receipt": int(candidate_audit["n_globally_certified_candidates"])
        == int(diagnostics.get("globally_kkt_certified_candidate_fits", -1)),
        "max_gap_matches_receipt": math.isclose(
            float(candidate_audit["max_global_log_likelihood_gap_bound"] or 0.0),
            float(diagnostics.get("max_global_log_likelihood_gap_bound", math.nan)),
            rel_tol=1e-12,
            abs_tol=1e-15,
        ),
        "max_simplex_matches_receipt": math.isclose(
            float(candidate_audit["max_simplex_residual"] or 0.0),
            float(diagnostics.get("max_simplex_residual", math.nan)),
            rel_tol=1e-12,
            abs_tol=1e-15,
        ),
    }
    failed = sorted(key for key, value in checks.items() if not value)
    require(not failed, f"table audits differ from ranking receipt: {failed}")
    return checks


_TIME_PATTERNS = {
    "elapsed": re.compile(r"^\s*Elapsed \(wall clock\) time .*?:\s*(\S+)\s*$"),
    "rss": re.compile(r"^\s*Maximum resident set size \(kbytes\):\s*(\d+)\s*$"),
    "exit": re.compile(r"^\s*Exit status:\s*(-?\d+)\s*$"),
}


def _parse_elapsed(value: str) -> float:
    parts = value.split(":")
    require(1 <= len(parts) <= 3, f"invalid elapsed time: {value!r}")
    try:
        numeric = [float(part) for part in parts]
    except ValueError as exc:
        raise PilotVerificationError(f"invalid elapsed time: {value!r}") from exc
    require(all(math.isfinite(item) and item >= 0 for item in numeric), "elapsed time is not finite/nonnegative")
    if len(numeric) == 3:
        hours, minutes, seconds = numeric
    elif len(numeric) == 2:
        hours, minutes, seconds = 0.0, numeric[0], numeric[1]
    else:
        hours, minutes, seconds = 0.0, 0.0, numeric[0]
    require(minutes < 60 and seconds < 60, f"invalid elapsed clock fields: {value!r}")
    return 3600 * hours + 60 * minutes + seconds


def parse_gnu_time(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing GNU time output: {path}")
    found: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="strict").splitlines():
        for key, pattern in _TIME_PATTERNS.items():
            match = pattern.match(line)
            if match:
                found[key] = match.group(1)
    require("elapsed" in found and "rss" in found, f"GNU time output lacks wall/RSS: {path}")
    if "exit" in found:
        require(int(found["exit"]) == 0, f"GNU time reports non-zero exit: {path}")
    elapsed = _parse_elapsed(found["elapsed"])
    rss_kib = int(found["rss"])
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": independent.sha256_path(path),
        "wall_seconds": elapsed,
        "wall_hours": elapsed / 3600,
        "maximum_resident_set_size_kib": rss_kib,
        "maximum_resident_set_size_gib": rss_kib / (1024 * 1024),
        "exit_status": int(found.get("exit", "0")),
    }


def _wall_status(seconds: float) -> str:
    return "GO" if seconds <= FOUR_HOURS_SECONDS else "PROBE" if seconds <= EIGHT_HOURS_SECONDS else "NO_GO"


def _rss_status(kib: int) -> str:
    return "GO" if kib <= FOUR_GIB_KIB else "PROBE" if kib <= EIGHT_GIB_KIB else "NO_GO"


def evaluate_gates(
    extraction_time: Mapping[str, Any],
    ranking_time: Mapping[str, Any],
    unit_audit: Mapping[str, Any],
    candidate_audit: Mapping[str, Any],
) -> dict[str, Any]:
    incomplete_rate = unit_audit.get("exact_enumeration_incomplete_rate")
    coverage = unit_audit.get("exact_limit_coverage")
    incomplete_abstain = int(unit_audit.get("n_incomplete_not_abstain", -1)) == 0
    no_incomplete_candidates = int(candidate_audit.get("n_candidates_in_primary_incomplete_units", -1)) == 0
    if incomplete_rate is None or not incomplete_abstain or not no_incomplete_candidates:
        incomplete_status = "NO_GO"
    elif float(incomplete_rate) <= 0.01:
        incomplete_status = "GO"
    elif float(incomplete_rate) <= 0.05:
        incomplete_status = "PROBE"
    else:
        incomplete_status = "NO_GO"
    if coverage is None:
        coverage_status = "NO_GO"
    elif float(coverage) >= 0.90:
        coverage_status = "GO"
    elif float(coverage) >= 0.80:
        coverage_status = "PROBE"
    else:
        coverage_status = "NO_GO"
    n_candidates = int(candidate_audit.get("n_candidates", 0))
    certified_status = (
        "GO" if n_candidates > 0
        and int(candidate_audit.get("n_globally_certified_candidates", -1)) == n_candidates
        and int(candidate_audit.get("n_uncertified_winners", -1)) == 0
        else "NO_GO"
    )
    refinement = candidate_audit.get("refinement_iterations") or {}
    p99 = refinement.get("p99_nearest_rank")
    maximum = refinement.get("max")
    refinement_status = (
        "GO" if p99 is not None and maximum is not None and int(p99) <= 100 and int(maximum) <= 1000
        else "PROBE" if p99 is not None and maximum is not None
        else "NO_GO"
    )
    expected_bootstrap = int(unit_audit.get("expected_bootstrap_replicates", 0))
    bootstrap_counts = unit_audit.get("bootstrap_status_counts") or {}
    bootstrap_evaluated = int(unit_audit.get("bootstrap_evaluated_units", 0))
    bootstrap_nonconverged = int(bootstrap_counts.get("ABSTAIN_NONCONVERGENCE", 0))
    bootstrap_status = (
        "GO" if expected_bootstrap == 0
        else "GO" if bootstrap_evaluated > 0 and bootstrap_nonconverged == 0
        else "PROBE"
    )
    metrics = {
        "extraction_wall": {
            "value_hours": extraction_time["wall_hours"],
            "go_max_hours": 4,
            "probe_max_hours": 8,
            "status": _wall_status(float(extraction_time["wall_seconds"])),
        },
        "ranking_wall": {
            "value_hours": ranking_time["wall_hours"],
            "go_max_hours": 4,
            "probe_max_hours": 8,
            "status": _wall_status(float(ranking_time["wall_seconds"])),
        },
        "extraction_peak_rss": {
            "value_gib": extraction_time["maximum_resident_set_size_gib"],
            "go_max_gib": 4,
            "probe_max_gib": 8,
            "status": _rss_status(int(extraction_time["maximum_resident_set_size_kib"])),
        },
        "ranking_peak_rss": {
            "value_gib": ranking_time["maximum_resident_set_size_gib"],
            "go_max_gib": 4,
            "probe_max_gib": 8,
            "status": _rss_status(int(ranking_time["maximum_resident_set_size_kib"])),
        },
        "k_le_12_exact_enumeration_incomplete_rate": {
            "value": incomplete_rate,
            "all_incomplete_abstain": incomplete_abstain,
            "no_candidates_from_incomplete_units": no_incomplete_candidates,
            "go_max": 0.01,
            "probe_max": 0.05,
            "status": incomplete_status,
        },
        "exact_limit_coverage": {
            "value": coverage,
            "go_min": 0.90,
            "probe_min": 0.80,
            "status": coverage_status,
        },
        "globally_certified_candidate_fits": {
            "value": candidate_audit.get("globally_certified_fraction"),
            "required": 1.0,
            "status": certified_status,
        },
        "refinement_iterations": {
            "p99_nearest_rank": p99,
            "max": maximum,
            "go_p99_max": 100,
            "go_absolute_max": 1000,
            "status": refinement_status,
            "policy": "iteration excess alone is PROBE; the independent 8 h wall gate remains NO_GO",
        },
        "conditional_candidate_ranking_bootstrap": {
            "expected_replicates": expected_bootstrap,
            "evaluated_units": bootstrap_evaluated,
            "status_counts": bootstrap_counts,
            "status": bootstrap_status,
            "policy": (
                "B0 is not requested; B>0 requires at least one evaluated unit and zero "
                "ABSTAIN_NONCONVERGENCE units for GO, otherwise PROBE"
            ),
        },
    }
    statuses = [entry["status"] for entry in metrics.values()]
    verdict = "NO_GO" if "NO_GO" in statuses else "PROBE" if "PROBE" in statuses else "GO"
    return {"verdict": verdict, "metrics": metrics, "all_metrics_go": verdict == "GO"}


def verify_pilot(
    pilot_root: Path,
    dataset: str,
    chrom: str,
    *,
    ranking_dirname: str = RANKING_DIRNAME,
    expected_bootstrap_replicates: int = 0,
    extraction_resource_gate_receipt: Path | None = None,
    ranking_resource_gate_receipt: Path | None = None,
) -> dict[str, Any]:
    root = pilot_root.resolve()
    require(root.is_dir(), f"pilot root is not a directory: {root}")
    require(
        ranking_dirname not in {"", ".", ".."}
        and not Path(ranking_dirname).is_absolute()
        and Path(ranking_dirname).name == ranking_dirname,
        "ranking dirname must be one relative path component",
    )
    require(expected_bootstrap_replicates >= 0, "expected bootstrap replicates must be nonnegative")
    require(
        (ranking_dirname, expected_bootstrap_replicates)
        in {("ranking_bootstrap0", 0), ("ranking_bootstrap20", 20)},
        "release pilot ranking dirname/bootstrap contract must be B0 or B20",
    )
    extraction_dir = root / EXTRACTION_DIRNAME
    ranking_dir = root / ranking_dirname
    extraction_resource_gate_receipt = (
        root / "resource_gates" / "extraction.json"
        if extraction_resource_gate_receipt is None
        else extraction_resource_gate_receipt
    )
    ranking_resource_gate_receipt = (
        root / "resource_gates" / f"{ranking_dirname}.json"
        if ranking_resource_gate_receipt is None
        else ranking_resource_gate_receipt
    )
    extraction_gate_audit = verify_pilot_resource_gate(
        extraction_resource_gate_receipt,
        pilot_root=root,
        output_root=extraction_dir,
        dataset=dataset,
        chrom=chrom,
        stage="extraction",
        gate_label="extraction",
    )
    ranking_gate_audit = verify_pilot_resource_gate(
        ranking_resource_gate_receipt,
        pilot_root=root,
        output_root=ranking_dir,
        dataset=dataset,
        chrom=chrom,
        stage="ranking",
        gate_label=ranking_dirname,
    )
    extraction_receipt, extraction_audit = verify_extraction(extraction_dir, dataset, chrom)
    extraction_child_manifest_identity = (
        extraction_audit["verified_provenance"]["manifest"]
    )
    require(
        extraction_gate_audit["canonical_manifest_identity"]
        == ranking_gate_audit["canonical_manifest_identity"]
        == extraction_child_manifest_identity,
        "pilot extraction/ranking gates and extraction child provenance do not share the exact "
        "canonical manifest identity",
    )
    require(
        extraction_gate_audit["release_manifest_identity"]
        == ranking_gate_audit["release_manifest_identity"],
        "pilot extraction and ranking gates do not share the exact frozen release binding",
    )
    resource_gate_cross_binding = {
        "canonical_manifest_identity": extraction_child_manifest_identity,
        "release_manifest_identity": extraction_gate_audit["release_manifest_identity"],
        "extraction_gate_manifest_equals_ranking_gate_manifest": True,
        "both_gate_manifests_equal_extraction_child_provenance": True,
        "extraction_gate_release_equals_ranking_gate_release": True,
    }
    ranking_receipt, ranking_contract_audit = verify_ranking_contract(
        ranking_dir,
        extraction_dir,
        extraction_receipt,
        dataset,
        chrom,
        expected_bootstrap_replicates,
    )
    parameters = ranking_receipt["parameters"]
    unit_path = Path(str(ranking_receipt["outputs"][UNIT_SOURCE_NAME]["path"]))
    candidate_path = Path(str(ranking_receipt["outputs"][independent.CANDIDATE_SOURCE_NAME]["path"]))
    unit_audit = audit_primary_units(
        unit_path,
        dataset,
        chrom,
        int(parameters["primary_structural_exact_pattern_minread"]),
        int(parameters["exact_k_max"]),
        expected_bootstrap_replicates,
        int(parameters["conditional_candidate_ranking_bootstrap_base_seed"]),
    )
    candidate_audit = audit_candidate_certificates(
        candidate_path,
        dataset,
        chrom,
        int(parameters["primary_structural_exact_pattern_minread"]),
        unit_audit["incomplete_unit_keys"],
    )
    table_receipt_checks = compare_table_audits_to_receipt(
        ranking_receipt, unit_audit, candidate_audit
    )
    try:
        profile_audit = independent.verify_child_profile_likelihood(
            ranking_receipt,
            dataset,
            chrom,
            expected_ranking_parameters(expected_bootstrap_replicates),
        )
    except independent.VerificationError as exc:
        raise PilotVerificationError(f"independent profile likelihood audit failed: {exc}") from exc
    require(
        int(profile_audit["n_candidates"]) == int(candidate_audit["n_primary_candidates"]),
        "profile audit does not cover every primary candidate",
    )
    extraction_time = parse_gnu_time(root / "extraction.time.txt")
    ranking_time = parse_gnu_time(root / f"{ranking_dirname}.time.txt")
    gates = evaluate_gates(extraction_time, ranking_time, unit_audit, candidate_audit)
    unit_audit_public = {key: value for key, value in unit_audit.items() if key != "incomplete_unit_keys"}
    receipt = {
        "schema_name": "intersubmod.m2_single_task_pilot_independent_verification",
        "schema_version": "1.0.0",
        "task_type": "B_COMPREHENSIVE_VALIDATION_RELEASE_PILOT",
        "scope": {
            "pilot_root": str(root),
            "dataset": dataset,
            "chrom": chrom,
            "n_tasks": 1,
            "ranking_dirname": ranking_dirname,
            "expected_bootstrap_replicates": expected_bootstrap_replicates,
        },
        "verification_independence": {
            "reads_bam": False,
            "reads_vcf": False,
            "launches_pilot_or_samtools": False,
            "authenticates_receipts_sidecars_and_output_identities": True,
            "streams_persisted_unit_and_candidate_tables": True,
            "profile_likelihood_recomputed_from_patterns_states_pi": True,
            "claim_boundary": (
                "This receipt gates one completed task only; it does not establish 154-task completeness, "
                "full-run disk/wall extrapolation, biological clone counts, or uniquely identified parent "
                "edges. Resource receipts attest launch-time snapshots under the non-hostile same-UID "
                "assumption; they are not cryptographic OS process-ancestry proofs."
            ),
        },
        "extraction": extraction_audit,
        "ranking": ranking_contract_audit,
        "unit_audit": unit_audit_public,
        "candidate_certificate_audit": candidate_audit,
        "profile_likelihood_independent_recomputation": profile_audit,
        "table_to_receipt_checks": table_receipt_checks,
        "process_resources": {
            "extraction": extraction_time,
            ranking_dirname: ranking_time,
        },
        "resource_gate_receipts": {
            "extraction": extraction_gate_audit,
            ranking_dirname: ranking_gate_audit,
        },
        "resource_gate_cross_binding": resource_gate_cross_binding,
        "release_gate": gates,
        "checks": {
            "extraction_schema_1_2_sidecar_outputs_scope_parameters_verified": True,
            "ranking_schema_2_0_sidecar_outputs_scope_parameters_verified": True,
            "ranking_bound_to_exact_extraction_receipt_and_outputs": True,
            "all_child_named_checks_true": True,
            "unit_routes_and_receipt_aggregate_match": True,
            "all_claimed_candidate_certificates_streamed": True,
            "profile_likelihood_recomputed_independently": (
                profile_audit.get("all_profile_likelihoods_and_certificates_match") is True
                and profile_audit.get("all_relative_weights_match") is True
                and profile_audit.get("all_winner_tie_partitions_match") is True
            ),
            "gnu_time_wall_and_rss_parsed": True,
            "extraction_and_ranking_resource_gates_authenticated": True,
            "resource_gate_manifest_and_release_cross_binding_verified": True,
            "release_gate_is_go": gates["verdict"] == "GO",
        },
    }
    receipt["all_pass"] = bool(receipt["checks"]) and all(
        value is True for value in receipt["checks"].values()
    )
    return receipt


def _fsync_directory(path: Path) -> None:
    """Durably publish directory-entry and permission changes."""
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    directory_fd = os.open(path, flags)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def _write_exclusive_immutable(path: Path, payload: bytes) -> None:
    """Create one single-link 0444 file without replacing any prior evidence."""
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    file_fd = os.open(path, flags, 0o444)
    try:
        with os.fdopen(file_fd, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
            os.fchmod(handle.fileno(), 0o444)
    except BaseException:
        # A partial file deliberately remains as fail-closed evidence.  It must
        # never be deleted or reused inside the same formal evidence root.
        raise
    observed = os.lstat(path)
    if observed.st_nlink != 1 or (observed.st_mode & 0o7777) != 0o444:
        raise RuntimeError(f"immutable single-link publication failed: {path}")


def write_receipt(path: Path, receipt: dict[str, Any]) -> str:
    """Durably publish immutable JSON+sidecar without an overwrite path."""
    sidecar = path.with_name(f"{path.name}.sha256")
    if os.path.lexists(path) or os.path.lexists(sidecar):
        raise FileExistsError(
            f"refusing to overwrite or complete preseeded verification evidence: "
            f"{path} / {sidecar}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    receipt["receipt_integrity"] = {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": sidecar.name,
        "covers": path.name,
    }
    payload = (
        json.dumps(receipt, ensure_ascii=False, allow_nan=False, indent=2) + "\n"
    ).encode("utf-8")
    sha = hashlib.sha256(payload).hexdigest()
    _write_exclusive_immutable(path, payload)
    _fsync_directory(path.parent)
    _write_exclusive_immutable(
        sidecar, f"{sha}  {path.name}\n".encode("ascii")
    )
    _fsync_directory(path.parent)
    return sha


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-root", required=True, type=Path)
    parser.add_argument("--dataset", required=True, choices=DATASETS)
    parser.add_argument("--chrom", required=True, choices=AUTOSOMES)
    parser.add_argument("--ranking-dirname", default=RANKING_DIRNAME)
    parser.add_argument("--expected-bootstrap-replicates", type=int, default=0)
    parser.add_argument("--extraction-resource-gate-receipt", required=True, type=Path)
    parser.add_argument("--resource-gate-receipt", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    default_output_name = (
        "pilot_gate_verification_receipt.json"
        if args.ranking_dirname == RANKING_DIRNAME
        and args.expected_bootstrap_replicates == 0
        else f"pilot_gate_verification_receipt.{args.ranking_dirname}.json"
    )
    output = args.output or args.pilot_root / default_output_name
    failure: str | None = None
    try:
        receipt = verify_pilot(
            args.pilot_root,
            args.dataset,
            args.chrom,
            ranking_dirname=args.ranking_dirname,
            expected_bootstrap_replicates=args.expected_bootstrap_replicates,
            extraction_resource_gate_receipt=args.extraction_resource_gate_receipt,
            ranking_resource_gate_receipt=args.resource_gate_receipt,
        )
    except (
        PilotVerificationError,
        independent.VerificationError,
        OSError,
        ValueError,
        KeyError,
        TypeError,
        csv.Error,
        gzip.BadGzipFile,
        json.JSONDecodeError,
        UnicodeError,
    ) as exc:
        failure = f"{type(exc).__name__}: {exc}"
        receipt = {
            "schema_name": "intersubmod.m2_single_task_pilot_independent_verification",
            "schema_version": "1.0.0",
            "task_type": "B_COMPREHENSIVE_VALIDATION_RELEASE_PILOT",
            "scope": {
                "pilot_root": str(args.pilot_root.resolve()),
                "dataset": args.dataset,
                "chrom": args.chrom,
                "n_tasks": 1,
                "ranking_dirname": args.ranking_dirname,
                "expected_bootstrap_replicates": args.expected_bootstrap_replicates,
            },
            "verification_independence": {
                "reads_bam": False,
                "reads_vcf": False,
                "launches_pilot_or_samtools": False,
            },
            "failure": failure,
            "release_gate": {"verdict": "NO_GO", "metrics": {}, "all_metrics_go": False},
            "checks": {"verification_completed_without_contract_violation": False},
            "all_pass": False,
        }
    sha = write_receipt(output, receipt)
    print(json.dumps({
        "output": str(output.resolve()),
        "sha256": sha,
        "verdict": (receipt.get("release_gate") or {}).get("verdict", "NO_GO"),
        "all_pass": receipt["all_pass"],
        "failure": failure,
    }, ensure_ascii=False))
    raise SystemExit(0 if receipt["all_pass"] else 1)


if __name__ == "__main__":
    main()
