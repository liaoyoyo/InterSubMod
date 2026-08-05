#!/usr/bin/env python3
"""Build the all-sSNV report from the v29 schema-recovery final dataset."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import sqlite3
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence


RECOVERY_DIR = Path(__file__).resolve().parent
RECOVERY_FINALIZER = (
    RECOVERY_DIR / "finalize_task_b_result_release_schema_recovery_v29.py"
)
TITLE = "全 sSNV focal-ALT 甲基多群與分子單倍型證據驗證"
PORTABLE_TITLE = "全 sSNV 甲基分群與共現驗證"
WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
CANONICAL_PYTHON_CACHE_DIRNAME = (
    ".python_cache_m2v5_completion_v2_bound_bootstrap"
)
CANONICAL_PYTHON_CACHE_ROOT = WORKSPACE_ROOT / CANONICAL_PYTHON_CACHE_DIRNAME
EXPECTED_SITES = 469_849
EXPECTED_DATASETS = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
TRUTH_LABELS = ("TP", "FP", "UNASSESSED")
CLAIM_IDS = ("M1", "M2", "G1", "G2", "R1", "B1", "C1", "L1", "L2")
REPORT_STRATUM_CLAIMS = ("M1", "M2", "G1", "G2")
CLAIM_STATUSES = frozenset({"PASS", "FAIL", "NOT_EVALUABLE", "NOT_RUN"})
CLAIM_NAMES = {
    "M1": "focal-ALT operational stable methyl-multigroup screen flag",
    "M2": "residual robust epigenetic partition",
    "G1": "LongPhase-S callset-anchored local co-segregation",
    "G2": "multi-marker molecular-haplotype base candidate",
    "R1": "strict methyl-partition robustness",
    "B1": "background-controlled molecular-haplotype candidate",
    "C1": "CN/CCF-conditioned candidate",
    "L1": "cellular subclone supported",
    "L2": "lineage/order supported",
}
M1_DENOMINATOR_DEFINITION = (
    "all LongPhase-S recalibrated FILTER=PASS chr1-22 biallelic focal "
    "dataset-sites in scope"
)
M2_DENOMINATOR_DEFINITION = (
    "M1 PASS sites with 2-10 observed methyl groups where each measured axis is "
    "either an observed aligned confound, adequately powered for non-alignment, "
    "or assignment-proven constant"
)
M2_CATEGORICAL_PLANNING_LEVEL_CEILINGS = {
    "hp_exact": 7,
    "hp_family": 5,
    "strand": 2,
}
M2_ASSIGNMENT_OBSERVED_LEVELS_ROLE = (
    "constant-axis proof only; observed levels do not replace the source-locked "
    "planning level ceilings"
)
BACKGROUND_CONTROL_REPLICATION_GATE_CONTRACT = (
    "lenient_coarse_modal_K2_without_membership_ARI_requirement_v1"
)
BACKGROUND_CONTROL_RELATION_TO_PRIMARY_M1 = (
    "lenient_predicate_superset_of_ARI_qualified_predicate_on_same_background_payload"
)
FINAL_SCHEMA = "intersubmod.all_ssnv_final_report_dataset"
FINAL_SCHEMA_VERSION = "2.0.0"
FINAL_RECEIPT_SCHEMA = "intersubmod.all_ssnv_final_report_dataset_run_receipt"
FINAL_RECEIPT_SCHEMA_VERSION = "2.0.0"
REPORT_RECEIPT_SCHEMA = "intersubmod.all_ssnv_report_build_receipt"
REPORT_RECEIPT_SCHEMA_VERSION = "1.2.0"
MANIFEST_SCHEMA = "intersubmod.all_ssnv_focal_alt_input_manifest"
MANIFEST_SCHEMA_VERSION = "1.0.0"
SCREEN_SCHEMA = "intersubmod.all_ssnv_focal_alt_multigroup_screen"
SCREEN_SCHEMA_VERSION = "1.2.0"
COOCCURRENCE_SCHEMA = "intersubmod.methyl_ssnv_cooccurrence_validation.summary"
COOCCURRENCE_SCHEMA_VERSION = "2.0.0"
RECONCILIATION_SCHEMA = "intersubmod.all_ssnv_output_reconciliation"
RECONCILIATION_SCHEMA_VERSION = "1.1.0"
IMMUTABILITY_SCHEMA = "intersubmod.all_ssnv_frozen_input_immutability"
IMMUTABILITY_SCHEMA_VERSION = "1.0.0"
TREE_AUDIT_SCHEMA = "intersubmod.latest_tree_input_contract_audit"
TREE_AUDIT_SCHEMA_VERSION = "1.0.0"
REFERENCE_AUDIT_SCHEMA = "intersubmod.all_ssnv_extraction_reference_identity_audit"
REFERENCE_AUDIT_SCHEMA_VERSION = "1.0.0"
REFERENCE_AUDIT_PASS_SEMANTICS = (
    "reference_identity_and_receipt_binding_only_not_scientific_confirmation"
)
LAYERED_SCHEMA = "intersubmod.layered_input_manifest"
LAYERED_SCHEMA_VERSION = "3.0.0"
PASS_SEMANTICS = "integration_integrity_only_not_scientific_confirmation"
MAX_SNAPSHOT_ROWS = 2_000
MAX_CASE_PAIR_ROWS = 50
CANONICAL_ORACLE = ("HCC1395_DORADO", "chr5", 750_311)
TERMINAL_M2_GATE_CONTRACT = (
    "m2-measured-axis-v3_effect-permutation-and-power-evaluability"
)
M2_GATE_CONTRACT = (
    "m2-measured-axis-v4_asymmetric-positive-confound-and-observed-categorical-levels"
)
JOINT_GLOBAL_FDR_FAMILY = "ELIGIBLE_M2_JOINT_SIGNATURE_GLOBAL_FDR_FAMILY"
FOUR_STATE_FAMILYWISE_CONFIDENCE = 0.95
FOUR_STATE_RELATION_FAMILY_SIZE = 3
FOUR_STATE_PER_RELATION_CONFIDENCE = 1.0 - (
    (1.0 - FOUR_STATE_FAMILYWISE_CONFIDENCE) / FOUR_STATE_RELATION_FAMILY_SIZE
)
FOUR_STATE_MULTIPLICITY_METHOD = "bonferroni_three_relation_models"
FOUR_STATE_MINIMUM_ZERO_VIOLATION_DEPTH = 203
CLAIM_CONTRACT_SOURCE_ROLE = "terminal or release claim contract"
TERMINAL_CLAIM_CONTRACT_SHA256 = (
    "da94a50d0717174ff007b75f2edad2de79bf3aebf6b15df179eb736e8d8f526e"
)
RELEASE_CLAIM_CONTRACT_SHA256 = (
    "da8778af6bcc9b1e8b2887eb2bcc87eca83d909be32a819ec5eb8f5f12c9f2af"
)
CLAIM_CONTRACT_REQUIRED_FRAGMENTS = (
    "terminal-claim-contract-v3",
    TERMINAL_M2_GATE_CONTRACT,
    JOINT_GLOBAL_FDR_FAMILY,
    "joint_signature_q_global_by",
    "joint_signature_global_by_discovery",
    "bonferroni_three_relation_models",
    "minimum zero-violation depth = 203",
    (
        "among_G1_formal_pairs_max_endpoint_a_n_informative_then_abs_distance_"
        "then_partner_identity_without_four_state_result"
    ),
    "source-locked screen contract v2",
)
RELEASE_CLAIM_CONTRACT_REQUIRED_FRAGMENTS = (
    "terminal-claim-contract-v5",
    M2_GATE_CONTRACT,
    JOINT_GLOBAL_FDR_FAMILY,
    "joint_signature_q_global_by",
    "joint_signature_global_by_discovery",
    "bonferroni_three_relation_models",
    "minimum zero-violation depth = 203",
    "只支援2-10個methyl groups",
    "陽性混雜證據即使低於80% planning power仍可進入M2 FAIL",
    "observed category count=1證明constant",
    "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM",
    "coarse_ng>=2",
    "同一background payload下，此predicate是加入ARI條件之predicate的superset",
    (
        "among_G1_formal_pairs_max_endpoint_a_n_informative_then_abs_distance_"
        "then_partner_identity_without_four_state_result"
    ),
    "source-locked screen contract v2",
)
FINAL_CANDIDATE_TSV_REQUIRED_FIELDS = (
    "candidate_key",
    "sample",
    "chrom",
    "pos",
    "ref",
    "alt",
    "m2_screen_gate_contract",
    "m2_low_power_axes",
    "joint_signature_global_fdr_family_status",
    "joint_signature_q_global_bh",
    "joint_signature_q_global_by",
    "joint_signature_global_bh_discovery",
    "joint_signature_global_by_discovery",
    "b1_prespecified_pair_key",
    "b1_prespecified_pair_is_witness",
    "b1_uses_posthoc_compatible_pair_search",
)
FINAL_WITNESS_TSV_REQUIRED_FIELDS = (
    "witness_pair_key",
    "sample",
    "focal_chrom",
    "focal_pos",
    "focal_ref",
    "focal_alt",
    "partner_chrom",
    "partner_pos",
    "partner_ref",
    "partner_alt",
    "global_by_q",
    "four_state_familywise_confidence",
    "four_state_relation_family_size",
    "four_state_multiplicity_method",
    "four_state_minimum_zero_violation_depth",
    "same_pair_four_state_witness",
    "b1_prespecified_pair",
)

SiteKey = tuple[str, str, int, str, str]


class ReportContractError(RuntimeError):
    """Raised when a report input violates a frozen v2 contract."""


@dataclass(frozen=True)
class ReportInputs:
    final_dataset: Path
    final_receipt: Path
    candidate_catalog: Path
    candidate_witness_pairs: Path
    claim_contract: Path
    manifest: Path
    screen_summary: Path
    cooccurrence_sites: Path
    cooccurrence_pairs: Path
    cooccurrence_summary: Path
    output_reconciliation: Path
    post_immutability_audit: Path
    tree_input_audit: Path
    reference_identity_audit: Path
    earlier_fp_report: Path | None = None
    final_release_receipt: Path | None = None
    final_release_signature: Path | None = None
    final_release_public_key: Path | None = None


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def require_file(path: Path, label: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file() or path.stat().st_size <= 0:
        raise ReportContractError(f"Missing or empty {label}: {path}")
    return path


def load_json(path: Path, label: str) -> dict[str, Any]:
    path = require_file(path, label)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ReportContractError(f"Invalid JSON in {label}: {path}") from error
    if not isinstance(payload, dict):
        raise ReportContractError(f"{label} must contain a JSON object: {path}")
    return payload


def canonical_python_prefix() -> list[str]:
    return [
        sys.executable,
        "-I",
        "-X",
        f"pycache_prefix={CANONICAL_PYTHON_CACHE_ROOT}",
    ]


def validate_detached_final_release(inputs: ReportInputs) -> dict[str, Any]:
    provided = (
        inputs.final_release_receipt,
        inputs.final_release_signature,
        inputs.final_release_public_key,
    )
    if all(path is None for path in provided):
        return {"status": "NOT_INCLUDED_NON_RELEASE_TEST_FIXTURE", "pass": None}
    if any(path is None for path in provided):
        raise ReportContractError("Detached final release inputs are incomplete")
    script = RECOVERY_FINALIZER
    command = [
        *canonical_python_prefix(),
        str(script),
        "--verify",
        "--receipt",
        str(inputs.final_release_receipt.resolve()),
        "--signature",
        str(inputs.final_release_signature.resolve()),
    ]
    verification = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    if verification.returncode != 0:
        raise ReportContractError(
            "Detached final result signature verification failed: "
            f"{verification.stderr.strip()}"
        )
    try:
        payload = json.loads(verification.stdout)
    except json.JSONDecodeError as error:
        raise ReportContractError("Detached final result verifier output is invalid") from error
    if not isinstance(payload, dict) or payload.get("pass") is not True:
        raise ReportContractError("Detached final result verifier did not pass")
    declared_key = payload.get("public_key")
    verify_declared_artifact(
        declared_key,
        inputs.final_release_public_key,
        "detached final result public key",
    )
    return payload


def open_text(path: Path):
    path = require_file(path, "TSV input")
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def iter_tsv(path: Path, required: Iterable[str], label: str) -> Iterator[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        if not fields or len(fields) != len(set(fields)):
            raise ReportContractError(f"{label} has a missing or duplicate header")
        missing = sorted(set(required).difference(fields))
        if missing:
            raise ReportContractError(f"{label} missing columns: {missing}")
        for row in reader:
            yield dict(row)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def display_path(path: Path, repo_root: Path) -> str:
    resolved = path.expanduser().resolve()
    try:
        return f"InterSubMod/{resolved.relative_to(repo_root.resolve()).as_posix()}"
    except ValueError:
        return str(resolved)


def portable_path(path: Path, repo_root: Path) -> str:
    resolved = path.expanduser().resolve()
    try:
        return f"InterSubMod/{resolved.relative_to(repo_root.resolve()).as_posix()}"
    except ValueError:
        parts = resolved.parts
        if "big7_disk_output" in parts:
            return Path(*parts[parts.index("big7_disk_output") :]).as_posix()
        digest = hashlib.sha256(str(resolved).encode()).hexdigest()[:12]
        return f"external-data/{digest}_{resolved.name}"


def file_identity(role: str, path: Path, repo_root: Path) -> dict[str, Any]:
    path = require_file(path, role)
    return {
        "role": role,
        "path": portable_path(path, repo_root),
        "display_path": display_path(path, repo_root),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def require_schema(
    payload: Mapping[str, Any],
    schema_name: str,
    schema_version: str,
    label: str,
) -> None:
    if payload.get("schema_name") != schema_name:
        raise ReportContractError(
            f"Unexpected {label} schema: {payload.get('schema_name')!r} != {schema_name!r}"
        )
    if payload.get("schema_version") != schema_version:
        raise ReportContractError(
            f"Unexpected {label} schema version: {payload.get('schema_version')!r} "
            f"!= {schema_version!r}"
        )
    if payload.get("pass") is not True:
        raise ReportContractError(f"{label} did not PASS")


def to_int(value: Any, *, field: str) -> int:
    if isinstance(value, bool):
        raise ReportContractError(f"Boolean is not an integer for {field}")
    try:
        result = int(value)
    except (TypeError, ValueError) as error:
        raise ReportContractError(f"Invalid integer for {field}: {value!r}") from error
    return result


def optional_float(value: Any, *, field: str) -> float | None:
    if value is None or str(value).strip() in {"", "NA", "N/A", "null", "None"}:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ReportContractError(f"Invalid float for {field}: {value!r}") from error
    if not math.isfinite(result):
        raise ReportContractError(f"Non-finite float for {field}: {value!r}")
    return result


def to_bool(value: Any, *, field: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes"}:
        return True
    if normalized in {"0", "false", "no", ""}:
        return False
    raise ReportContractError(f"Invalid boolean for {field}: {value!r}")


def ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def format_count(value: int) -> str:
    return f"{value:,}"


def format_percent(value: float | None) -> str:
    return "NA" if value is None else f"{100 * value:.2f}%"


def site_key(row: Mapping[str, Any], *, pair: bool = False) -> SiteKey:
    prefix = "focal_" if pair else ""
    try:
        key = (
            str(row[f"{prefix}sample"] if pair else row["sample"]),
            str(row[f"{prefix}chrom"]),
            int(row[f"{prefix}pos"]),
            str(row[f"{prefix}ref"]).upper(),
            str(row[f"{prefix}alt"]).upper(),
        )
    except (KeyError, TypeError, ValueError) as error:
        if pair:
            try:
                key = (
                    str(row["sample"]),
                    str(row["focal_chrom"]),
                    int(row["focal_pos"]),
                    str(row["focal_ref"]).upper(),
                    str(row["focal_alt"]).upper(),
                )
            except (KeyError, TypeError, ValueError) as nested_error:
                raise ReportContractError(f"Malformed pair focal key: {row!r}") from nested_error
        else:
            raise ReportContractError(f"Malformed site key: {row!r}") from error
    if key[0] not in EXPECTED_DATASETS or not key[1] or key[2] <= 0 or not key[3] or not key[4]:
        raise ReportContractError(f"Invalid site key: {key!r}")
    return key


def verify_declared_artifact(reference: Any, path: Path, label: str) -> None:
    if not isinstance(reference, Mapping):
        raise ReportContractError(f"Missing artifact record for {label}")
    path = require_file(path, label)
    declared_path = reference.get("path")
    if not declared_path or Path(str(declared_path)).expanduser().resolve() != path:
        raise ReportContractError(f"{label} declared path does not match supplied path")
    if to_int(reference.get("size_bytes"), field=f"{label}.size_bytes") != path.stat().st_size:
        raise ReportContractError(f"{label} declared size does not match current file")
    declared_hash = reference.get("sha256")
    if not isinstance(declared_hash, str) or len(declared_hash) != 64:
        raise ReportContractError(f"{label} lacks declared SHA-256")
    if sha256(path) != declared_hash:
        raise ReportContractError(f"{label} SHA-256 mismatch")


def verify_embedded_artifact(reference: Any, label: str) -> Path:
    if not isinstance(reference, Mapping) or not reference.get("path"):
        raise ReportContractError(f"Missing artifact record for {label}")
    path = Path(str(reference["path"])).expanduser().resolve()
    verify_declared_artifact(reference, path, label)
    return path


def parse_utc_timestamp(value: Any, *, label: str) -> datetime:
    if not isinstance(value, str) or not value.strip():
        raise ReportContractError(f"{label} timestamp is missing")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise ReportContractError(f"{label} timestamp is invalid: {value!r}") from error
    if parsed.tzinfo is None:
        raise ReportContractError(f"{label} timestamp lacks timezone")
    return parsed.astimezone(timezone.utc)


def validate_final_receipt(
    payload: Mapping[str, Any], inputs: ReportInputs
) -> Path:
    require_schema(
        payload,
        FINAL_RECEIPT_SCHEMA,
        FINAL_RECEIPT_SCHEMA_VERSION,
        "final dataset receipt",
    )
    if payload.get("task_type") != "B_comprehensive_validation":
        raise ReportContractError("Final dataset receipt task type is not comprehensive validation")
    if payload.get("pass_semantics") != PASS_SEMANTICS:
        raise ReportContractError("Final dataset receipt pass semantics drift")
    command = payload.get("command")
    if (
        not isinstance(command, list)
        or not command
        or any(not isinstance(token, str) or not token.strip() for token in command)
    ):
        raise ReportContractError("Final dataset receipt command must be a non-empty string array")
    outputs = payload.get("outputs")
    if not isinstance(outputs, Mapping):
        raise ReportContractError("Final dataset receipt lacks outputs")
    verify_declared_artifact(
        outputs.get("final_report_dataset"),
        inputs.final_dataset,
        "final receipt output final_report_dataset",
    )
    verify_declared_artifact(
        outputs.get("candidate_catalog"),
        inputs.candidate_catalog,
        "final receipt output candidate_catalog",
    )
    verify_declared_artifact(
        outputs.get("candidate_witness_pairs"),
        inputs.candidate_witness_pairs,
        "final receipt output candidate_witness_pairs",
    )
    code = payload.get("code")
    if not isinstance(code, Mapping):
        raise ReportContractError("Final dataset receipt lacks code identities")
    return verify_embedded_artifact(
        code.get("final_report_dataset_builder"),
        "final_report_dataset_builder code",
    )


def serialized_tsv_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return str(value)


def validate_final_machine_tables(final: Mapping[str, Any], inputs: ReportInputs) -> None:
    candidate_rows = list(
        iter_tsv(
            inputs.candidate_catalog,
            FINAL_CANDIDATE_TSV_REQUIRED_FIELDS,
            "final candidate catalog TSV",
        )
    )
    witness_rows = list(
        iter_tsv(
            inputs.candidate_witness_pairs,
            FINAL_WITNESS_TSV_REQUIRED_FIELDS,
            "final candidate witness-pair TSV",
        )
    )
    candidate_json = final.get("candidate_catalog")
    witness_json = final.get("candidate_witness_pairs")
    counts = final.get("counts")
    if not isinstance(candidate_json, list) or not isinstance(witness_json, list):
        raise ReportContractError("Final dataset candidate machine tables are not arrays")
    if not isinstance(counts, Mapping):
        raise ReportContractError("Final dataset lacks candidate machine-table counts")
    expected_candidate_count = to_int(
        counts.get("candidate_catalog_rows"), field="candidate_catalog_rows"
    )
    expected_witness_count = to_int(
        counts.get("candidate_witness_pair_rows"), field="candidate_witness_pair_rows"
    )
    if len(candidate_rows) != len(candidate_json) or len(candidate_rows) != expected_candidate_count:
        raise ReportContractError("Candidate catalog TSV/JSON/count rows do not reconcile")
    if len(witness_rows) != len(witness_json) or len(witness_rows) != expected_witness_count:
        raise ReportContractError("Candidate witness-pair TSV/JSON/count rows do not reconcile")

    for label, rows, payload_rows, key_field, required_fields in (
        (
            "candidate catalog",
            candidate_rows,
            candidate_json,
            "candidate_key",
            FINAL_CANDIDATE_TSV_REQUIRED_FIELDS,
        ),
        (
            "candidate witness pairs",
            witness_rows,
            witness_json,
            "witness_pair_key",
            FINAL_WITNESS_TSV_REQUIRED_FIELDS,
        ),
    ):
        if any(not isinstance(row, Mapping) for row in payload_rows):
            raise ReportContractError(f"Final JSON {label} contains a non-object row")
        tsv_by_key = {row[key_field]: row for row in rows}
        json_by_key = {str(row.get(key_field, "")): row for row in payload_rows}
        if (
            len(tsv_by_key) != len(rows)
            or len(json_by_key) != len(payload_rows)
            or not all(tsv_by_key)
            or not all(json_by_key)
            or set(tsv_by_key) != set(json_by_key)
        ):
            raise ReportContractError(f"Final {label} key set is missing, duplicate, or inconsistent")
        for key, json_row in json_by_key.items():
            tsv_row = tsv_by_key[key]
            for field in required_fields:
                if tsv_row[field] != serialized_tsv_value(json_row.get(field)):
                    raise ReportContractError(
                        f"Final {label} TSV/JSON drift for key={key} field={field}"
                    )


def validate_all_final_input_artifacts(
    final: Mapping[str, Any], final_receipt: Mapping[str, Any]
) -> dict[str, Path]:
    records = final.get("input_artifacts")
    receipt_records = final_receipt.get("inputs")
    if not isinstance(records, Mapping) or not records:
        raise ReportContractError("Final dataset lacks complete input_artifacts")
    if receipt_records != records:
        raise ReportContractError("Final dataset and receipt input_artifacts differ")
    paths: dict[str, Path] = {}
    for role, reference in records.items():
        if not isinstance(role, str) or not role:
            raise ReportContractError("Final input artifact role is malformed")
        paths[role] = verify_embedded_artifact(
            reference, f"final input artifact {role}"
        )
    required_roles = {
        "manifest",
        "screen_sites",
        "screen_assignments",
        "screen_summary",
        "screen_receipt",
        "cooccurrence_sites",
        "cooccurrence_pairs",
        "cooccurrence_summary",
        "cooccurrence_receipt",
        "tumor_ref_sites",
        "tumor_ref_summary",
        "tumor_ref_receipt",
        "primary_artifact_audit_pre",
        "primary_artifact_audit_post",
    }
    missing = sorted(required_roles.difference(paths))
    if missing:
        raise ReportContractError(
            f"Final input artifact inventory is incomplete: {missing}"
        )
    return paths


def validate_recovery_provenance_artifacts(
    screen_receipt: Mapping[str, Any],
) -> dict[str, Path]:
    execution = screen_receipt.get("execution")
    if not isinstance(execution, Mapping):
        raise ReportContractError("Terminal screen receipt lacks execution metadata")
    is_recovery = execution.get("recovery_merge") is True
    recovery = screen_receipt.get("recovery")
    if not is_recovery:
        if recovery not in (None, {}):
            raise ReportContractError(
                "Non-recovery terminal screen contains stray recovery metadata"
            )
        return {}
    if not isinstance(recovery, Mapping):
        raise ReportContractError("Recovery terminal screen lacks recovery provenance")
    if recovery.get("serial_parallel_exact_equivalence_required") is not True:
        raise ReportContractError("Recovery screen did not require exact equivalence")
    checks = recovery.get("serial_parallel_exact_equivalence_checks")
    if not isinstance(checks, Mapping) or not checks or not all(
        value is True for value in checks.values()
    ):
        raise ReportContractError("Recovery exact-equivalence checks are incomplete")
    if (
        recovery.get("prefix_source_files_complete_closed_run") is not True
        or recovery.get("prefix_read_reached_eof") is not True
        or recovery.get("prefix_replacement_source_dependencies_exact") is not True
        or recovery.get("recovery_source_identity_unchanged_during_merge") is not True
    ):
        raise ReportContractError("Recovery source/provenance gates did not all PASS")
    pinned_sha = str(recovery.get("pinned_analyzer_sha256", ""))
    if re.fullmatch(r"[0-9a-f]{64}", pinned_sha) is None:
        raise ReportContractError("Recovery pinned analyzer SHA-256 is invalid")
    prefix = recovery.get("prefix_source_artifacts")
    if not isinstance(prefix, Mapping):
        raise ReportContractError("Recovery prefix artifact inventory is missing")
    references: dict[str, Any] = {
        "recovery prefix sites": prefix.get("sites"),
        "recovery prefix assignments": prefix.get("assignments"),
        "recovery prefix summary": prefix.get("summary"),
        "recovery prefix receipt": prefix.get("receipt"),
        "recovery prefix source-lock receipt": prefix.get("source_lock_receipt"),
        "recovery replacement summary": recovery.get("replacement_summary"),
        "recovery replacement receipt": recovery.get("replacement_receipt"),
        "recovery serial-parallel equivalence receipt": recovery.get(
            "serial_parallel_exact_equivalence_receipt"
        ),
    }
    return {
        role: verify_embedded_artifact(reference, role)
        for role, reference in references.items()
    }


def validate_claim_contract(path: Path, *, release_gate_pass: bool) -> str:
    path = require_file(path, "claim contract")
    observed_sha256 = sha256(path)
    expected_sha256 = (
        RELEASE_CLAIM_CONTRACT_SHA256
        if release_gate_pass
        else TERMINAL_CLAIM_CONTRACT_SHA256
    )
    if observed_sha256 != expected_sha256:
        raise ReportContractError(
            "Claim contract SHA-256 mismatch for release state: "
            f"release_gate_pass={release_gate_pass} expected={expected_sha256} "
            f"observed={observed_sha256}"
        )
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as error:
        raise ReportContractError(f"Cannot read claim contract: {path}") from error
    missing_ids = [
        claim_id
        for claim_id in CLAIM_IDS
        if re.search(rf"(?<![A-Za-z0-9_]){re.escape(claim_id)}(?![A-Za-z0-9_])", text)
        is None
    ]
    missing_names = [name for name in CLAIM_NAMES.values() if name not in text]
    required_fragments = (
        RELEASE_CLAIM_CONTRACT_REQUIRED_FRAGMENTS
        if release_gate_pass
        else CLAIM_CONTRACT_REQUIRED_FRAGMENTS
    )
    missing_fragments = [fragment for fragment in required_fragments if fragment not in text]
    if missing_ids or missing_names or missing_fragments:
        raise ReportContractError(
            "Claim contract does not contain the complete terminal claim vocabulary: "
            f"missing_ids={missing_ids}, missing_names={missing_names}, "
            f"missing_fragments={missing_fragments}"
        )
    return text


def validate_manifest(payload: Mapping[str, Any]) -> dict[str, int]:
    require_schema(payload, MANIFEST_SCHEMA, MANIFEST_SCHEMA_VERSION, "input manifest")
    samples = payload.get("samples")
    if not isinstance(samples, list) or tuple(row.get("sample") for row in samples) != EXPECTED_DATASETS:
        raise ReportContractError("Input manifest must contain the canonical 7-dataset order")
    totals = payload.get("totals")
    if not isinstance(totals, Mapping):
        raise ReportContractError("Input manifest lacks totals")
    if to_int(totals.get("all_ssnv"), field="manifest totals.all_ssnv") != EXPECTED_SITES:
        raise ReportContractError("Input manifest does not contain 469,849 sites")
    truth_total = sum(
        to_int(totals.get(f"truth_{label.lower()}"), field=f"manifest truth {label}")
        for label in TRUTH_LABELS
    )
    if truth_total != EXPECTED_SITES:
        raise ReportContractError("Manifest TP+FP+UNASSESSED does not reconcile to 469,849")
    sample_counts: dict[str, int] = {}
    sample_total = 0
    for row in samples:
        counts = row.get("counts")
        if not isinstance(counts, Mapping):
            raise ReportContractError(f"Manifest sample lacks counts: {row.get('sample')}")
        sample = str(row["sample"])
        population = to_int(counts.get("all_ssnv"), field=f"{sample}.all_ssnv")
        sample_truth = sum(
            to_int(counts.get(f"truth_{label.lower()}"), field=f"{sample}.{label}")
            for label in TRUTH_LABELS
        )
        if sample_truth != population:
            raise ReportContractError(f"Manifest truth strata do not reconcile for {sample}")
        sample_counts[sample] = population
        sample_total += population
    if sample_total != EXPECTED_SITES:
        raise ReportContractError("Manifest per-dataset totals do not reconcile")
    return sample_counts


def load_and_validate_layered_manifest(
    manifest: Mapping[str, Any],
) -> tuple[dict[str, Any], Path, dict[str, Any]]:
    reference = manifest.get("layered_root")
    if not isinstance(reference, Mapping) or not reference.get("path"):
        raise ReportContractError("Input manifest lacks the frozen layered_root artifact")
    path = Path(str(reference["path"])).expanduser().resolve()
    verify_declared_artifact(reference, path, "layered input manifest")
    payload = load_json(path, "layered input manifest")
    if payload.get("schema_name") != LAYERED_SCHEMA or payload.get("schema_version") != LAYERED_SCHEMA_VERSION:
        raise ReportContractError("Layered input manifest schema/version drift")
    if payload.get("dataset_count") != 7 or payload.get("biological_sample_count") != 6:
        raise ReportContractError("Layered input manifest dataset/sample scope drift")
    production = payload.get("production_summary")
    contract = payload.get("analysis_contract")
    if not isinstance(production, Mapping) or not isinstance(contract, Mapping):
        raise ReportContractError("Layered input manifest lacks production or analysis contract")
    gates = {
        "producer_7_of_7": production.get("all_pass") is True
        and production.get("passed_dataset_count") == 7,
        "clairs_raw_all_input": contract.get("longphase_input_contract")
        == "normalized_ClairS_raw_all",
        "longphase_pass_tree": contract.get("tree_input_contract")
        == "longphase_s_recalibrated_FILTER_PASS",
        "external_sidecar": contract.get("read_tag_mode") == "external_sidecar",
        "embedded_tags_ignored": contract.get("embedded_tag_policy") == "ignore",
        "exact_join_required": contract.get("require_exact_join") is True,
    }
    if not all(gates.values()):
        failed = sorted(name for name, passed in gates.items() if not passed)
        raise ReportContractError(f"Current ClairS/LongPhase-S flow contract failed: {failed}")
    samples = payload.get("samples")
    if not isinstance(samples, list) or tuple(row.get("sample") for row in samples) != EXPECTED_DATASETS:
        raise ReportContractError("Layered manifest lacks the canonical producer-receipt sample order")
    bam_rows: list[dict[str, Any]] = []
    for sample_row in samples:
        sample = str(sample_row["sample"])
        read_tags = sample_row.get("read_tags")
        if not isinstance(read_tags, Mapping):
            raise ReportContractError(f"Layered manifest lacks read-tag contract for {sample}")
        declared = read_tags.get("producer_capture_receipt")
        if not isinstance(declared, Mapping) or not declared.get("path"):
            raise ReportContractError(f"Layered manifest lacks producer capture receipt for {sample}")
        identity = declared.get("identity")
        if not isinstance(identity, Mapping):
            raise ReportContractError(f"Producer capture receipt lacks identity for {sample}")
        receipt_path = Path(str(declared["path"])).expanduser().resolve()
        verify_declared_artifact(
            {"path": str(receipt_path), **identity},
            receipt_path,
            f"{sample} producer capture receipt",
        )
        receipt = load_json(receipt_path, f"{sample} producer capture receipt")
        if (
            receipt.get("schema_name") != "intersubmod.longphase_raw_all_capture_receipt"
            or receipt.get("schema_version") != "2.0.0"
            or receipt.get("sample") != sample
        ):
            raise ReportContractError(f"Producer capture receipt identity drift for {sample}")
        policy = receipt.get("bam_output_policy")
        policy_gates = {
            "named_fifo": isinstance(policy, Mapping) and policy.get("transport") == "named_fifo",
            "not_persisted": isinstance(policy, Mapping) and policy.get("persisted_bam") is False,
            "no_regular_bam": isinstance(policy, Mapping) and policy.get("regular_bam_count") == 0,
            "fifo_at_closeout": isinstance(policy, Mapping)
            and policy.get("is_fifo_at_closeout") is True,
        }
        if not all(policy_gates.values()):
            failed = sorted(name for name, passed in policy_gates.items() if not passed)
            raise ReportContractError(f"Producer BAM-output policy failed for {sample}: {failed}")
        bam_rows.append(
            {
                "sample": sample,
                "transport": policy["transport"],
                "persisted_bam": policy["persisted_bam"],
                "regular_bam_count": policy["regular_bam_count"],
                "is_fifo_at_closeout": policy["is_fifo_at_closeout"],
                "receipt_path": str(receipt_path),
                "receipt_sha256": identity["sha256"],
                "status": "PASS",
            }
        )
    bam_audit = {
        "schema_name": "intersubmod.report_producer_bam_output_audit",
        "schema_version": "1.0.0",
        "dataset_count": len(bam_rows),
        "named_fifo_count": sum(row["transport"] == "named_fifo" for row in bam_rows),
        "persisted_bam_count": sum(bool(row["persisted_bam"]) for row in bam_rows),
        "regular_bam_count": sum(int(row["regular_bam_count"]) for row in bam_rows),
        "rows": bam_rows,
        "pass": len(bam_rows) == len(EXPECTED_DATASETS),
    }
    return payload, path, bam_audit


def validate_screen_summary(payload: Mapping[str, Any]) -> dict[str, Any]:
    require_schema(payload, SCREEN_SCHEMA, SCREEN_SCHEMA_VERSION, "current screen summary")
    scope = payload.get("scope")
    if not isinstance(scope, Mapping):
        raise ReportContractError("Screen summary lacks scope")
    if scope.get("full_469849") is not True:
        raise ReportContractError("Screen summary is not marked full_469849")
    selected = scope.get("selected_datasets", scope.get("selected_samples"))
    if tuple(selected or ()) != EXPECTED_DATASETS:
        raise ReportContractError("Screen summary does not cover the canonical 7 datasets")
    for field in ("expected_sites", "processed_sites"):
        if to_int(scope.get(field), field=f"screen scope.{field}") != EXPECTED_SITES:
            raise ReportContractError(f"Screen summary {field} is not 469,849")
    audit = payload.get("latest_hp_ps_terminal_join_audit")
    if not isinstance(audit, Mapping):
        raise ReportContractError("Screen summary lacks current terminal HP/PS audit")
    tag_gates = {
        "pass": audit.get("pass") is True,
        "same_run_sidecar": audit.get("authoritative_tag_source")
        == "same_run_LongPhase_S_external_HP_PS_sidecar",
        "embedded_tags_unused": audit.get("embedded_reads_tsv_hp_used_for_analysis") is False,
        "join_before_alt": audit.get("join_occurs_before_focal_ALT_selection") is True,
        "all_sites": to_int(audit.get("n_sites"), field="tag audit n_sites") == EXPECTED_SITES,
        "status_counts": audit.get("join_status_counts") == {"PASS": EXPECTED_SITES},
        "all_sites_pass": audit.get("all_sites_pass") is True,
        "every_row_joined": audit.get("every_reads_tsv_row_joined") is True,
        "unique_projection": audit.get("all_projection_identities_unique") is True,
        "zero_multimatch": to_int(
            audit.get("n_projection_multimatch_site_reads"), field="tag audit multimatch"
        )
        == 0,
    }
    joined = to_int(
        audit.get("n_exact_hp_ps_site_read_joins"), field="tag audit exact joins"
    )
    read_rows = to_int(audit.get("n_reads_tsv_site_rows"), field="tag audit reads rows")
    tag_gates["joined_count_reconciliation"] = joined == read_rows
    if not all(tag_gates.values()):
        failed = sorted(name for name, passed in tag_gates.items() if not passed)
        raise ReportContractError(f"Current terminal HP/PS audit failed closed: {failed}")
    pooled = payload.get("pooled_site_weighted")
    if not isinstance(pooled, Mapping) or to_int(pooled.get("n_sites"), field="screen pooled n") != EXPECTED_SITES:
        raise ReportContractError("Screen pooled site count does not reconcile")
    return dict(audit)


def validate_cooccurrence_summary(payload: Mapping[str, Any]) -> None:
    require_schema(
        payload,
        COOCCURRENCE_SCHEMA,
        COOCCURRENCE_SCHEMA_VERSION,
        "cooccurrence summary",
    )
    if tuple(payload.get("samples", ())) != EXPECTED_DATASETS:
        raise ReportContractError("Cooccurrence summary dataset order/content drift")
    reconciliation = payload.get("site_pair_count_reconciliation")
    if not isinstance(reconciliation, Mapping) or reconciliation.get("pass") is not True:
        raise ReportContractError("Cooccurrence site/pair reconciliation did not PASS")


def validate_output_reconciliation(payload: Mapping[str, Any]) -> None:
    require_schema(
        payload,
        RECONCILIATION_SCHEMA,
        RECONCILIATION_SCHEMA_VERSION,
        "output reconciliation",
    )
    if payload.get("status") != "EXECUTION_PASS":
        raise ReportContractError("Output reconciliation status is not EXECUTION_PASS")
    totals = payload.get("totals")
    if not isinstance(totals, Mapping):
        raise ReportContractError("Output reconciliation lacks totals")
    exact_fields = ("expected_sites", "reads_keys", "methylation_keys", "bernoulli_keys")
    if any(to_int(totals.get(field), field=f"reconciliation {field}") != EXPECTED_SITES for field in exact_fields):
        raise ReportContractError("Output reconciliation does not cover all 469,849 site artifacts")
    zero_fields = (
        "region_failures",
        "reads_missing",
        "reads_extra",
        "reads_duplicates",
        "reads_empty",
        "methylation_missing",
        "methylation_extra",
        "methylation_duplicates",
        "methylation_empty",
        "bernoulli_missing",
        "bernoulli_extra",
        "bernoulli_duplicates",
        "bernoulli_empty",
    )
    if any(to_int(totals.get(field, 0), field=f"reconciliation {field}") != 0 for field in zero_fields):
        raise ReportContractError("Output reconciliation contains missing/extra/duplicate/empty artifacts")
    datasets = payload.get("datasets", payload.get("samples"))
    if not isinstance(datasets, list) or tuple(row.get("sample", row.get("dataset")) for row in datasets) != EXPECTED_DATASETS:
        raise ReportContractError("Output reconciliation does not contain canonical 7-dataset rows")
    if not all(row.get("pass") is True for row in datasets):
        raise ReportContractError("At least one output reconciliation dataset failed")


def validate_post_immutability(
    payload: Mapping[str, Any], path: Path, screen_receipt: Mapping[str, Any]
) -> None:
    require_schema(
        payload,
        IMMUTABILITY_SCHEMA,
        IMMUTABILITY_SCHEMA_VERSION,
        "post-run immutability audit",
    )
    phase = str(payload.get("audit_phase", "")).lower()
    if phase not in {"post", "post_run", "post-run"} and "post" not in path.name.lower():
        raise ReportContractError("Immutability input is not identifiable as a post-run audit")
    totals = payload.get("totals")
    if not isinstance(totals, Mapping):
        raise ReportContractError("Post-run immutability audit lacks totals")
    n_samples = to_int(totals.get("n_samples"), field="immutability n_samples")
    n_sample_pass = to_int(totals.get("n_sample_pass"), field="immutability n_sample_pass")
    n_artifacts = to_int(totals.get("n_artifacts"), field="immutability n_artifacts")
    n_artifact_pass = to_int(
        totals.get("n_artifact_pass"), field="immutability n_artifact_pass"
    )
    if n_samples != 7 or n_sample_pass != 7 or n_artifacts <= 0 or n_artifacts != n_artifact_pass:
        raise ReportContractError("Post-run immutability counts did not reconcile")
    audit_created = parse_utc_timestamp(
        payload.get("created_at_utc"), label="post-run immutability creation"
    )
    screen_finished = parse_utc_timestamp(
        screen_receipt.get("finished_at_utc"), label="terminal screen finish"
    )
    if audit_created < screen_finished:
        raise ReportContractError(
            "Post-run immutability audit predates the selected terminal screen"
        )


def validate_tree_audit(payload: Mapping[str, Any]) -> None:
    require_schema(payload, TREE_AUDIT_SCHEMA, TREE_AUDIT_SCHEMA_VERSION, "tree-input audit")
    totals = payload.get("totals")
    if not isinstance(totals, Mapping) or to_int(totals.get("all_ssnv"), field="tree all_ssnv") != EXPECTED_SITES:
        raise ReportContractError("Tree-input audit does not cover 469,849 sites")
    checks = payload.get("top_level_checks")
    if not isinstance(checks, Mapping) or not checks or not all(value is True for value in checks.values()):
        raise ReportContractError("Tree-input audit top-level checks did not all PASS")
    samples = payload.get("samples")
    if not isinstance(samples, list) or tuple(row.get("sample") for row in samples) != EXPECTED_DATASETS:
        raise ReportContractError("Tree-input audit does not contain canonical 7-dataset rows")
    if not all(row.get("pass") is True for row in samples):
        raise ReportContractError("At least one tree-input dataset failed")


def validate_reference_identity_audit(payload: Mapping[str, Any]) -> None:
    require_schema(
        payload,
        REFERENCE_AUDIT_SCHEMA,
        REFERENCE_AUDIT_SCHEMA_VERSION,
        "extraction reference identity audit",
    )
    if payload.get("task_type") != "B_comprehensive_validation":
        raise ReportContractError(
            "Extraction reference identity audit task type is not comprehensive validation"
        )
    if payload.get("pass_semantics") != REFERENCE_AUDIT_PASS_SEMANTICS:
        raise ReportContractError("Extraction reference identity audit pass semantics drift")
    if "469,849" not in str(payload.get("scope", "")):
        raise ReportContractError(
            "Extraction reference identity audit scope does not identify 469,849 dataset-sites"
        )

    reference = payload.get("reference")
    if not isinstance(reference, Mapping):
        raise ReportContractError("Extraction reference identity audit lacks reference identity")
    if to_int(reference.get("size_bytes"), field="reference size_bytes") <= 0:
        raise ReportContractError("Extraction reference identity audit reference is empty")
    if re.fullmatch(r"[0-9a-f]{64}", str(reference.get("full_sha256", ""))) is None:
        raise ReportContractError("Extraction reference full SHA-256 is missing or malformed")
    fai = reference.get("fai")
    if not isinstance(fai, Mapping):
        raise ReportContractError("Extraction reference identity audit lacks FAI identity")
    if to_int(fai.get("size_bytes"), field="reference FAI size_bytes") <= 0:
        raise ReportContractError("Extraction reference FAI is empty")
    if re.fullmatch(r"[0-9a-f]{64}", str(fai.get("full_sha256", ""))) is None:
        raise ReportContractError("Extraction reference FAI full SHA-256 is missing or malformed")

    receipts = payload.get("sample_receipts")
    if not isinstance(receipts, list) or len(receipts) != len(EXPECTED_DATASETS):
        raise ReportContractError(
            "Extraction reference identity audit does not contain 7 sample receipts"
        )
    observed_samples = [str(row.get("sample", "")) for row in receipts if isinstance(row, Mapping)]
    if len(observed_samples) != len(receipts) or set(observed_samples) != set(EXPECTED_DATASETS):
        raise ReportContractError(
            "Extraction reference identity audit receipt dataset membership drift"
        )
    receipt_gates = (
        all(row.get("receipt_pass") is True for row in receipts),
        all(to_int(row.get("exit_code"), field="reference receipt exit_code") == 0 for row in receipts),
        all(row.get("reference_path_equal") is True for row in receipts),
        all(row.get("command_reference_path_equal") is True for row in receipts),
        all(re.fullmatch(r"[0-9a-f]{64}", str(row.get("sha256", ""))) is not None for row in receipts),
    )
    if not all(receipt_gates):
        raise ReportContractError(
            "Extraction reference identity audit sample receipt bindings did not all PASS"
        )

    checks = payload.get("checks")
    if not isinstance(checks, Mapping) or not checks or not all(value is True for value in checks.values()):
        raise ReportContractError(
            "Extraction reference identity audit checks did not all PASS"
        )
    limitations = payload.get("limitations")
    if not isinstance(limitations, list) or not limitations or not all(
        isinstance(value, str) and value.strip() for value in limitations
    ):
        raise ReportContractError(
            "Extraction reference identity audit must preserve its limitations"
        )


def validate_final_scope(payload: Mapping[str, Any]) -> None:
    require_schema(payload, FINAL_SCHEMA, FINAL_SCHEMA_VERSION, "final report dataset")
    if payload.get("task_type") != "B_comprehensive_validation":
        raise ReportContractError("Final dataset task type is not comprehensive validation")
    if payload.get("pass_semantics") != PASS_SEMANTICS:
        raise ReportContractError("Final dataset pass semantics drift")
    scope = payload.get("scope")
    if not isinstance(scope, Mapping):
        raise ReportContractError("Final dataset lacks scope")
    if tuple(scope.get("datasets", ())) != EXPECTED_DATASETS:
        raise ReportContractError("Final dataset does not contain canonical 7 datasets")
    if to_int(scope.get("dataset_count"), field="final dataset_count") != 7:
        raise ReportContractError("Final dataset_count is not 7")
    if to_int(scope.get("biological_sample_count"), field="final biological_sample_count") != 6:
        raise ReportContractError("Final biological_sample_count is not 6")
    for field in ("expected_screen_sites", "observed_screen_sites"):
        if to_int(scope.get(field), field=f"final scope.{field}") != EXPECTED_SITES:
            raise ReportContractError(f"Final {field} is not 469,849")
    technical = scope.get("technical_replicate")
    if not isinstance(technical, Mapping) or technical.get("counts_as_independent_biological_n") is not False:
        raise ReportContractError("HCC technical replicate is incorrectly counted as biological n")


def validate_final_input_links(final: Mapping[str, Any], inputs: ReportInputs) -> None:
    records = final.get("input_artifacts")
    if not isinstance(records, Mapping):
        raise ReportContractError("Final dataset lacks input_artifacts")
    expected = {
        "manifest": inputs.manifest,
        "screen_summary": inputs.screen_summary,
        "cooccurrence_sites": inputs.cooccurrence_sites,
        "cooccurrence_pairs": inputs.cooccurrence_pairs,
        "cooccurrence_summary": inputs.cooccurrence_summary,
    }
    for role, path in expected.items():
        verify_declared_artifact(records.get(role), path, role)


def normalized_metric(payload: Any, claim_id: str, label: str) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ReportContractError(f"Missing metric {label}/{claim_id}")
    values = {
        field: to_int(payload.get(field), field=f"{label}/{claim_id}/{field}")
        for field in ("numerator", "denominator", "not_evaluable", "not_run", "population")
    }
    if any(value < 0 for value in values.values()):
        raise ReportContractError(f"Negative metric count in {label}/{claim_id}")
    if values["numerator"] > values["denominator"]:
        raise ReportContractError(f"Numerator exceeds denominator in {label}/{claim_id}")
    if values["denominator"] + values["not_evaluable"] + values["not_run"] != values["population"]:
        raise ReportContractError(f"Status counts do not reconcile in {label}/{claim_id}")
    observed_ratio = optional_float(payload.get("ratio"), field=f"{label}/{claim_id}/ratio")
    expected_ratio = ratio(values["numerator"], values["denominator"])
    if observed_ratio is None and expected_ratio is not None:
        raise ReportContractError(f"Missing ratio in {label}/{claim_id}")
    if observed_ratio is not None and (
        expected_ratio is None
        or not math.isclose(observed_ratio, expected_ratio, rel_tol=0, abs_tol=1e-12)
    ):
        raise ReportContractError(f"Ratio does not reconcile in {label}/{claim_id}")
    definition = payload.get("denominator_definition")
    if not isinstance(definition, str) or not definition.strip():
        raise ReportContractError(f"Missing denominator definition in {label}/{claim_id}")
    return {**values, "ratio": observed_ratio, "denominator_definition": definition}


def normalize_metric_stratum(
    payload: Any,
    label: str,
    expected_population: int,
) -> dict[str, dict[str, Any]]:
    if not isinstance(payload, Mapping):
        raise ReportContractError(f"Missing metric stratum: {label}")
    if any(claim_id not in payload for claim_id in CLAIM_IDS):
        raise ReportContractError(f"Metric stratum lacks one or more claims: {label}")
    metrics = {
        claim_id: normalized_metric(payload[claim_id], claim_id, label)
        for claim_id in CLAIM_IDS
    }
    if any(metric["population"] != expected_population for metric in metrics.values()):
        raise ReportContractError(f"Metric population drift in {label}")
    return metrics


def _sum_metric_field(
    strata: Iterable[Mapping[str, Mapping[str, Any]]], claim_id: str, field: str
) -> int:
    return sum(to_int(stratum[claim_id][field], field=f"aggregate {claim_id}/{field}") for stratum in strata)


def validate_metrics(
    final: Mapping[str, Any], manifest: Mapping[str, Any]
) -> dict[str, Any]:
    funnel = final.get("funnel_metrics")
    if not isinstance(funnel, Mapping):
        raise ReportContractError("Final dataset lacks funnel_metrics")
    manifest_samples = {str(row["sample"]): row["counts"] for row in manifest["samples"]}
    manifest_totals = manifest["totals"]
    pooled = normalize_metric_stratum(funnel.get("pooled"), "pooled", EXPECTED_SITES)
    per_sample: dict[str, dict[str, dict[str, Any]]] = {}
    source_per_sample = funnel.get("per_sample")
    if not isinstance(source_per_sample, Mapping) or set(source_per_sample) != set(EXPECTED_DATASETS):
        raise ReportContractError("Final per_sample metric strata drift")
    for sample in EXPECTED_DATASETS:
        per_sample[sample] = normalize_metric_stratum(
            source_per_sample[sample],
            f"sample/{sample}",
            to_int(manifest_samples[sample]["all_ssnv"], field=f"manifest {sample}"),
        )
    truth: dict[str, dict[str, dict[str, Any]]] = {}
    source_truth = funnel.get("truth_strata")
    if not isinstance(source_truth, Mapping) or set(source_truth) != set(TRUTH_LABELS):
        raise ReportContractError("Final truth_strata key set drift")
    for truth_label in TRUTH_LABELS:
        truth[truth_label] = normalize_metric_stratum(
            source_truth[truth_label],
            f"truth/{truth_label}",
            to_int(
                manifest_totals[f"truth_{truth_label.lower()}"],
                field=f"manifest truth {truth_label}",
            ),
        )
    sample_truth: dict[str, dict[str, dict[str, Any]]] = {}
    source_sample_truth = funnel.get("sample_by_truth")
    expected_keys = tuple(f"{sample}|{truth_label}" for sample in EXPECTED_DATASETS for truth_label in TRUTH_LABELS)
    if not isinstance(source_sample_truth, Mapping) or set(source_sample_truth) != set(expected_keys):
        raise ReportContractError("Final sample_by_truth key set drift")
    for key in expected_keys:
        sample, truth_label = key.split("|", 1)
        sample_truth[key] = normalize_metric_stratum(
            source_sample_truth[key],
            f"sample_truth/{key}",
            to_int(
                manifest_samples[sample][f"truth_{truth_label.lower()}"],
                field=f"manifest {key}",
            ),
        )
    fields = ("numerator", "denominator", "not_evaluable", "not_run", "population")
    for claim_id in CLAIM_IDS:
        for field in fields:
            pooled_value = pooled[claim_id][field]
            if _sum_metric_field(per_sample.values(), claim_id, field) != pooled_value:
                raise ReportContractError(f"Per-sample {claim_id}/{field} does not reconcile to pooled")
            if _sum_metric_field(truth.values(), claim_id, field) != pooled_value:
                raise ReportContractError(f"Truth {claim_id}/{field} does not reconcile to pooled")
        for sample in EXPECTED_DATASETS:
            cells = [sample_truth[f"{sample}|{truth_label}"] for truth_label in TRUTH_LABELS]
            for field in fields:
                if _sum_metric_field(cells, claim_id, field) != per_sample[sample][claim_id][field]:
                    raise ReportContractError(
                        f"Sample-by-truth {sample}/{claim_id}/{field} does not reconcile"
                    )
        for truth_label in TRUTH_LABELS:
            cells = [sample_truth[f"{sample}|{truth_label}"] for sample in EXPECTED_DATASETS]
            for field in fields:
                if _sum_metric_field(cells, claim_id, field) != truth[truth_label][claim_id][field]:
                    raise ReportContractError(
                        f"Sample-by-truth {truth_label}/{claim_id}/{field} does not reconcile"
                    )
    status_counts = final.get("counts", {}).get("claim_status_counts")
    if not isinstance(status_counts, Mapping):
        raise ReportContractError("Final dataset lacks claim_status_counts")
    for claim_id in CLAIM_IDS:
        counts = status_counts.get(claim_id)
        if not isinstance(counts, Mapping) or any(status not in counts for status in CLAIM_STATUSES):
            raise ReportContractError(f"Missing claim status counts for {claim_id}")
        expected = {
            "PASS": pooled[claim_id]["numerator"],
            "FAIL": pooled[claim_id]["denominator"] - pooled[claim_id]["numerator"],
            "NOT_EVALUABLE": pooled[claim_id]["not_evaluable"],
            "NOT_RUN": pooled[claim_id]["not_run"],
        }
        observed = {status: to_int(counts[status], field=f"{claim_id}/{status}") for status in CLAIM_STATUSES}
        if observed != expected:
            raise ReportContractError(f"Claim status counts do not reconcile for {claim_id}")
    return {
        "pooled": pooled,
        "per_sample": per_sample,
        "truth_strata": truth,
        "sample_by_truth": sample_truth,
    }


def validate_m1_operational_screen(
    final: Mapping[str, Any], metrics: Mapping[str, Any]
) -> dict[str, Any]:
    payload = final.get("m1_operational_screen")
    if not isinstance(payload, Mapping):
        raise ReportContractError("Final dataset lacks m1_operational_screen")
    pooled = metrics["pooled"]["M1"]
    if (
        pooled["denominator"] != EXPECTED_SITES
        or pooled["population"] != EXPECTED_SITES
        or pooled["not_evaluable"] != 0
        or pooled["not_run"] != 0
        or pooled["denominator_definition"] != M1_DENOMINATOR_DEFINITION
    ):
        raise ReportContractError("M1 is not an all-469,849-site operational flag yield")
    n_all = to_int(payload.get("n_all_dataset_sites"), field="M1 n_all_dataset_sites")
    n_evaluable = to_int(payload.get("n_screen_evaluable"), field="M1 n_screen_evaluable")
    n_not_evaluable = to_int(
        payload.get("n_screen_not_evaluable"), field="M1 n_screen_not_evaluable"
    )
    n_flagged = to_int(
        payload.get("n_flagged_stable_multigroup"), field="M1 n_flagged"
    )
    n_not_flagged = to_int(payload.get("n_not_flagged_all"), field="M1 n_not_flagged")
    n_evaluable_not_flagged = to_int(
        payload.get("n_evaluable_not_flagged"), field="M1 n_evaluable_not_flagged"
    )
    if n_all != EXPECTED_SITES or n_evaluable + n_not_evaluable != n_all:
        raise ReportContractError("M1 screen evaluability counts do not reconcile")
    if n_flagged != pooled["numerator"] or n_flagged + n_not_flagged != n_all:
        raise ReportContractError("M1 flagged/not-flagged counts do not reconcile")
    if n_evaluable_not_flagged != n_evaluable - n_flagged or n_evaluable < n_flagged:
        raise ReportContractError("M1 evaluable/not-flagged counts do not reconcile")
    if payload.get("status_semantics") != "FLAGGED_VS_NOT_FLAGGED_OPERATIONAL_SCREEN_ONLY":
        raise ReportContractError("M1 operational status semantics drift")
    if payload.get("denominator_definition") != M1_DENOMINATOR_DEFINITION:
        raise ReportContractError("M1 operational denominator definition drift")
    if payload.get("global_null_validity_exported_for_nonstable_sites") is not False:
        raise ReportContractError("M1 global-null limitation is not preserved")
    if payload.get("biological_prevalence_estimate") is not None:
        raise ReportContractError("M1 incorrectly exposes a biological prevalence estimate")
    observed_all = optional_float(payload.get("flag_yield"), field="M1 flag_yield")
    observed_evaluable = optional_float(
        payload.get("flag_yield_among_screen_evaluable"),
        field="M1 flag_yield_among_screen_evaluable",
    )
    if observed_all != ratio(n_flagged, n_all):
        raise ReportContractError("M1 all-site flag yield does not reconcile")
    if observed_evaluable != ratio(n_flagged, n_evaluable):
        raise ReportContractError("M1 screen-evaluable secondary yield does not reconcile")
    return {
        "n_all_dataset_sites": n_all,
        "n_screen_evaluable": n_evaluable,
        "n_screen_not_evaluable": n_not_evaluable,
        "n_flagged_stable_multigroup": n_flagged,
        "n_not_flagged_all": n_not_flagged,
        "n_evaluable_not_flagged": n_evaluable_not_flagged,
        "flag_yield": observed_all,
        "flag_yield_among_screen_evaluable": observed_evaluable,
    }


def validate_m2_evaluability_contract(
    final: Mapping[str, Any], metrics: Mapping[str, Any]
) -> dict[str, Any]:
    payload = final.get("m2_evaluability_contract")
    if not isinstance(payload, Mapping):
        raise ReportContractError("Final dataset lacks m2_evaluability_contract")
    m1 = metrics["pooled"]["M1"]
    m2 = metrics["pooled"]["M2"]
    if payload.get("gate_contract") != M2_GATE_CONTRACT:
        raise ReportContractError("M2 evaluability gate contract drift")
    if (
        payload.get("denominator_definition") != M2_DENOMINATOR_DEFINITION
        or m2["denominator_definition"] != M2_DENOMINATOR_DEFINITION
    ):
        raise ReportContractError("M2 evaluability denominator definition drift")
    if (
        to_int(payload.get("minimum_supported_methyl_groups"), field="M2 minimum groups")
        != 2
        or to_int(payload.get("maximum_supported_methyl_groups"), field="M2 maximum groups")
        != 10
    ):
        raise ReportContractError("M2 supported methyl-group range must be 2-10")
    n_m1 = to_int(payload.get("n_m1_pass"), field="M2 n_m1_pass")
    n_evaluable = to_int(payload.get("n_m2_evaluable"), field="M2 n_evaluable")
    n_not_evaluable = to_int(
        payload.get("n_m2_not_evaluable"), field="M2 n_not_evaluable"
    )
    if (
        n_m1 != m1["numerator"]
        or n_evaluable != m2["denominator"]
        or n_not_evaluable != m2["not_evaluable"]
        or n_evaluable + n_not_evaluable != n_m1
    ):
        raise ReportContractError("M2 evaluability counts do not reconcile")
    reasons = payload.get("not_evaluable_reason_counts")
    if not isinstance(reasons, Mapping):
        raise ReportContractError("M2 not-evaluable reason census is missing")
    parsed_reasons = {
        str(reason): to_int(count, field=f"M2 reason {reason}")
        for reason, count in reasons.items()
    }
    if any(count <= 0 for count in parsed_reasons.values()):
        raise ReportContractError("M2 not-evaluable reason counts must be positive")
    if sum(parsed_reasons.values()) != n_not_evaluable:
        raise ReportContractError("M2 not-evaluable reason census does not reconcile")
    group_reason = "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM"
    n_group_limit = to_int(
        payload.get("n_group_count_exceeds_planning_model_maximum"),
        field="M2 group-limit count",
    )
    if n_group_limit != parsed_reasons.get(group_reason, 0):
        raise ReportContractError("M2 group-limit count does not reconcile")
    examples = payload.get("group_count_exceeds_planning_model_examples")
    if not isinstance(examples, list) or len(examples) > n_group_limit:
        raise ReportContractError("M2 group-limit examples are malformed")
    examples_complete = payload.get("group_count_exceeds_examples_complete")
    if examples_complete is not (n_group_limit <= 50):
        raise ReportContractError("M2 group-limit example completeness flag drift")
    expected_example_count = n_group_limit if examples_complete else 50
    if len(examples) != expected_example_count:
        raise ReportContractError("M2 group-limit example count drift")
    parsed_examples: list[dict[str, Any]] = []
    for example in examples:
        if not isinstance(example, Mapping):
            raise ReportContractError("M2 group-limit example is not an object")
        observed_groups = to_int(
            example.get("observed_methyl_groups"), field="M2 example observed groups"
        )
        if observed_groups <= 10 or example.get("maximum_supported_methyl_groups") != 10:
            raise ReportContractError("M2 group-limit example does not exceed 10 groups")
        expected_statuses = {
            "m1_status": "PASS",
            "m2_status": "NOT_EVALUABLE",
            "g1_status": "NOT_RUN",
            "g2_status": "NOT_RUN",
            "b1_status": "NOT_RUN",
        }
        if any(example.get(field) != value for field, value in expected_statuses.items()):
            raise ReportContractError("M2 group-limit downstream claim behavior drift")
        if not str(example.get("reason", "")).startswith(group_reason + ":"):
            raise ReportContractError("M2 group-limit example reason drift")
        parsed_examples.append(dict(example))
    expected_behavior = {
        "M1": "PASS retained",
        "M2": "NOT_EVALUABLE excluded from PASS/FAIL denominator",
        "G1": "NOT_RUN",
        "G2": "NOT_RUN",
        "B1": "NOT_RUN",
    }
    if payload.get("group_count_exceeds_claim_behavior") != expected_behavior:
        raise ReportContractError("M2 group-limit claim behavior contract drift")
    if (
        payload.get("categorical_planning_level_ceilings")
        != M2_CATEGORICAL_PLANNING_LEVEL_CEILINGS
        or payload.get("assignment_observed_levels_role")
        != M2_ASSIGNMENT_OBSERVED_LEVELS_ROLE
    ):
        raise ReportContractError("M2 categorical planning-level contract drift")
    independent = payload.get("independent_logic_audit")
    if independent is None:
        independent = {"status": "NOT_INCLUDED", "pass": None}
    if not isinstance(independent, Mapping):
        raise ReportContractError("M2 independent logic audit metadata is malformed")
    independent_status = str(independent.get("status", ""))
    if independent_status == "PASS_LOGIC_INDEPENDENT_RECOUNT":
        audit_counts = independent.get("counts")
        if not isinstance(audit_counts, Mapping):
            raise ReportContractError("M2 independent logic audit counts are missing")
        expected_audit_counts = {
            "all_rows": EXPECTED_SITES,
            "m1_stable_rows": n_m1,
            "eligible": m2["numerator"],
            "evaluable_ineligible": n_evaluable - m2["numerator"],
            "not_evaluable_axis_indeterminate": parsed_reasons.get(
                "NOT_EVALUABLE_M2_AXIS_INDETERMINATE", 0
            ),
            "not_evaluable_group_count_gt10": n_group_limit,
        }
        if {key: to_int(audit_counts.get(key), field=f"M2 audit {key}") for key in expected_audit_counts} != expected_audit_counts:
            raise ReportContractError("M2 independent logic audit counts do not reconcile")
        if (
            independent.get("production_gate_imported") is not False
            or independent.get("production_gate_functions_called") is not False
        ):
            raise ReportContractError("M2 independent logic audit imported or called production gate")
        aligned_below_power = to_int(
            independent.get(
                "n_evaluable_sites_with_aligned_axis_below_negative_evaluability_power"
            ),
            field="M2 aligned-below-power evaluable sites",
        )
        if aligned_below_power < 0:
            raise ReportContractError("M2 aligned-below-power count is negative")
    elif independent_status != "NOT_INCLUDED":
        raise ReportContractError("M2 independent logic audit status drift")
    return {
        "gate_contract": M2_GATE_CONTRACT,
        "denominator_definition": M2_DENOMINATOR_DEFINITION,
        "minimum_supported_methyl_groups": 2,
        "maximum_supported_methyl_groups": 10,
        "n_m1_pass": n_m1,
        "n_m2_evaluable": n_evaluable,
        "n_m2_not_evaluable": n_not_evaluable,
        "not_evaluable_reason_counts": parsed_reasons,
        "n_group_count_exceeds_planning_model_maximum": n_group_limit,
        "group_count_exceeds_planning_model_examples": parsed_examples,
        "group_count_exceeds_examples_complete": examples_complete,
        "independent_logic_audit": dict(independent),
    }


def aggregate_status(metric: Mapping[str, Any]) -> str:
    if metric["numerator"] > 0:
        return "PASS"
    if metric["denominator"] > 0:
        return "FAIL"
    if metric["not_evaluable"] > 0:
        return "NOT_EVALUABLE"
    return "NOT_RUN"


def validate_claim_ladder(
    final: Mapping[str, Any], metrics: Mapping[str, Any]
) -> list[dict[str, Any]]:
    source = final.get("claim_ladder")
    if not isinstance(source, list) or tuple(row.get("claim_id") for row in source) != CLAIM_IDS:
        raise ReportContractError("Claim ladder must contain ordered M1/M2/G1/G2/R1/B1/C1/L1/L2")
    rows: list[dict[str, Any]] = []
    for order, row in enumerate(source, 1):
        claim_id = str(row["claim_id"])
        if row.get("claim_name") != CLAIM_NAMES[claim_id]:
            raise ReportContractError(f"Claim name drift for {claim_id}")
        pooled = metrics["pooled"][claim_id]
        observed = {
            "dataset_numerator": to_int(row.get("dataset_numerator"), field=f"{claim_id} numerator"),
            "dataset_denominator": to_int(row.get("dataset_denominator"), field=f"{claim_id} denominator"),
            "dataset_not_evaluable": to_int(
                row.get("dataset_not_evaluable"), field=f"{claim_id} not_evaluable"
            ),
            "dataset_not_run": to_int(row.get("dataset_not_run"), field=f"{claim_id} not_run"),
        }
        expected = {
            "dataset_numerator": pooled["numerator"],
            "dataset_denominator": pooled["denominator"],
            "dataset_not_evaluable": pooled["not_evaluable"],
            "dataset_not_run": pooled["not_run"],
        }
        if observed != expected:
            raise ReportContractError(f"Claim ladder counts do not reconcile for {claim_id}")
        dataset_ratio = optional_float(row.get("dataset_ratio"), field=f"{claim_id} dataset_ratio")
        if dataset_ratio != pooled["ratio"]:
            raise ReportContractError(f"Claim ladder ratio does not reconcile for {claim_id}")
        status = str(row.get("status"))
        if status not in CLAIM_STATUSES or status != aggregate_status(pooled):
            raise ReportContractError(f"Claim ladder aggregate status drift for {claim_id}")
        biological_numerator = to_int(
            row.get("biological_numerator"), field=f"{claim_id} biological numerator"
        )
        biological_denominator = to_int(
            row.get("biological_denominator"), field=f"{claim_id} biological denominator"
        )
        biological_ratio = optional_float(
            row.get("biological_ratio"), field=f"{claim_id} biological ratio"
        )
        if biological_numerator > biological_denominator:
            raise ReportContractError(f"Biological numerator exceeds denominator for {claim_id}")
        if biological_ratio != ratio(biological_numerator, biological_denominator):
            raise ReportContractError(f"Biological ratio does not reconcile for {claim_id}")
        if row.get("automatic_upgrade_prohibited") is not True:
            raise ReportContractError(f"Automatic claim upgrade is not prohibited for {claim_id}")
        rows.append(
            {
                "claim_order": order,
                "claim_id": claim_id,
                "claim_name": CLAIM_NAMES[claim_id],
                **observed,
                "dataset_ratio": dataset_ratio,
                "biological_numerator": biological_numerator,
                "biological_denominator": biological_denominator,
                "biological_ratio": biological_ratio,
                "status": status,
                "denominator_definition": str(row.get("denominator_definition", "")),
                "guardrail": str(row.get("guardrail", "")),
            }
        )
    for claim_id in ("L1", "L2"):
        row = rows[CLAIM_IDS.index(claim_id)]
        if (
            row["dataset_numerator"] != 0
            or row["dataset_denominator"] != 0
            or row["status"] not in {"NOT_EVALUABLE", "NOT_RUN"}
        ):
            raise ReportContractError(
                f"{claim_id} must remain NOT_RUN/NOT_EVALUABLE without orthogonal evidence"
            )
    return rows


def validate_candidates(final: Mapping[str, Any], metrics: Mapping[str, Any]) -> list[dict[str, Any]]:
    candidates = final.get("candidate_catalog")
    witnesses = final.get("candidate_witness_pairs")
    if not isinstance(candidates, list) or not isinstance(witnesses, list):
        raise ReportContractError("Final dataset candidate arrays are malformed")
    g2_count = metrics["pooled"]["G2"]["numerator"]
    counts = final.get("counts")
    if not isinstance(counts, Mapping):
        raise ReportContractError("Final dataset lacks counts")
    for field in ("g2_candidates", "candidate_catalog_rows"):
        if to_int(counts.get(field), field=f"counts.{field}") != g2_count:
            raise ReportContractError(f"{field} does not reconcile to G2 PASS")
    if len(candidates) != g2_count:
        raise ReportContractError("Candidate catalog row count does not reconcile to G2 PASS")
    keys: set[SiteKey] = set()
    for row in candidates:
        if not isinstance(row, Mapping):
            raise ReportContractError("Candidate catalog contains a non-object row")
        key = site_key(row)
        if key in keys:
            raise ReportContractError(f"Duplicate candidate key: {key}")
        keys.add(key)
        if row.get("g2_status") != "PASS":
            raise ReportContractError(f"Candidate catalog row is not G2 PASS: {key}")
        if (
            row.get("m2_screen_gate_contract") != M2_GATE_CONTRACT
            or row.get("m2_screen_evaluable") is not True
            or not str(row.get("m2_screen_eligibility_status", "")).startswith(
                "ELIGIBLE_M2_RESIDUAL_UNEXPLAINED_AND_AXES_DETERMINATE"
            )
            or row.get("m2_indeterminate_axes") != []
            or row.get("m2_low_power_axes") != []
        ):
            raise ReportContractError(f"Candidate M2 measured-axis gate drift: {key}")
        joint_q_by = optional_float(
            row.get("joint_signature_q_global_by"),
            field=f"candidate {key} joint_signature_q_global_by",
        )
        if (
            row.get("joint_signature_global_fdr_family_status")
            != JOINT_GLOBAL_FDR_FAMILY
            or joint_q_by is None
            or joint_q_by > 0.05
            or row.get("joint_signature_global_by_discovery") is not True
        ):
            raise ReportContractError(f"Candidate G2 global-BY gate drift: {key}")
        witness_count = to_int(
            row.get("n_same_pair_four_state_witnesses"),
            field=f"candidate {key} prespecified B1 witness count",
        )
        if (
            witness_count not in {0, 1}
            or row.get("b1_uses_posthoc_compatible_pair_search") is not False
            or row.get("b1_prespecified_pair_is_witness") != (witness_count == 1)
            or (row.get("b1_status") == "PASS" and witness_count != 1)
        ):
            raise ReportContractError(f"Candidate B1 prespecified-pair gate drift: {key}")
        if row.get("b1_status") == "PASS" and (
            row.get("tumor_ref_status") != "PASS"
            or row.get("normal_control_status") != "PASS"
        ):
            raise ReportContractError(
                f"Candidate B1 PASS lacks passing tumor/normal background controls: {key}"
            )
    if g2_count == 0 and witnesses:
        raise ReportContractError("Zero G2 result contains candidate witness rows")
    return [dict(row) for row in candidates]


def validate_m2_axis_statistic_provenance(final: Mapping[str, Any]) -> dict[str, Any]:
    payload = final.get("m2_axis_statistic_provenance")
    if not isinstance(payload, Mapping):
        raise ReportContractError("Final dataset lacks M2 axis-statistic provenance")
    expected_checks = {
        "axis sample-size reconciliation",
        "499-permutation add-one p-value grid",
        "effect threshold classification",
        "80-percent planning-power evaluability",
        "assignment-derived categorical constant-axis proof",
        "asymmetric positive-confound versus negative-evaluability power decision",
    }
    if (
        payload.get("axis_effect_and_permutation_p_source")
        != "source_locked_focal_alt_multigroup_screen_producer"
        or payload.get("downstream_raw_read_axis_statistic_recomputed") is not False
        or payload.get("downstream_axis_classification_recomputed") is not True
        or set(payload.get("downstream_checks") or []) != expected_checks
        or payload.get("categorical_planning_level_ceilings")
        != M2_CATEGORICAL_PLANNING_LEVEL_CEILINGS
        or payload.get("assignment_observed_levels_role")
        != M2_ASSIGNMENT_OBSERVED_LEVELS_ROLE
        or "producer-derived" not in str(payload.get("claim_guardrail", ""))
    ):
        raise ReportContractError("M2 axis-statistic provenance contract drift")
    recovery = payload.get("screen_recovery_source_validation")
    if not isinstance(recovery, Mapping) or (
        recovery.get("mode") != "source_locked_prefix_plus_seed_parallel_replacement"
        or recovery.get("prefix_source_lock_pass") is not True
        or recovery.get("serial_parallel_exact_equivalence_pass") is not True
        or not re.fullmatch(r"[0-9a-f]{64}", str(recovery.get("pinned_analyzer_sha256", "")))
    ):
        raise ReportContractError("M2 screen-producer source validation drift")
    return dict(payload)


def validate_background_control_replication_gate(
    final: Mapping[str, Any],
) -> dict[str, Any]:
    payload = final.get("background_control_replication_gate")
    if not isinstance(payload, Mapping):
        raise ReportContractError(
            "Final dataset lacks background_control_replication_gate"
        )
    if payload.get("contract") != BACKGROUND_CONTROL_REPLICATION_GATE_CONTRACT:
        raise ReportContractError("Background-control replication gate contract drift")
    if payload.get("relation_to_primary_m1_replication_flags") != (
        BACKGROUND_CONTROL_RELATION_TO_PRIMARY_M1
    ):
        raise ReportContractError("Background-control relation to primary M1 drift")
    if payload.get("applies_to") != ["tumor_REF", "matched_normal_REF"]:
        raise ReportContractError("Background-control gate scope drift")
    if payload.get("required_conditions") != [
        "coarse_ng>=2",
        "modal_fraction>=0.7_via_unstable_false",
    ]:
        raise ReportContractError("Background-control required conditions drift")
    if payload.get("membership_ari_minimum_required") is not False:
        raise ReportContractError(
            "Background-control gate must disclose that membership ARI is not required"
        )
    if payload.get("b1_pass_direction") != (
        "requires_no_lenient_background_replication"
    ):
        raise ReportContractError("Background-control B1 pass direction drift")
    if payload.get("false_positive_direction") != (
        "cannot_increase_B1_passes_vs_ARI_qualified_predicate_on_same_background_payload"
    ):
        raise ReportContractError("Background-control false-positive direction drift")
    if payload.get("false_negative_direction") != (
        "may_conservatively_reduce_B1_passes_when_K_is_stable_but_membership_is_not"
    ):
        raise ReportContractError("Background-control false-negative direction drift")
    if payload.get("scientific_interpretation") != (
        "background nonreplication guardrail, not an exact primary-M1 replay"
    ):
        raise ReportContractError("Background-control scientific interpretation drift")
    return dict(payload)


def validate_tumor_ref_source_identity_attestation(
    final: Mapping[str, Any],
) -> dict[str, Any]:
    payload = final.get("tumor_ref_source_identity_attestation")
    if not isinstance(payload, Mapping):
        raise ReportContractError("Final dataset lacks tumor-REF source identity attestation")
    status = payload.get("status")
    if status == "NOT_INCLUDED_INTERMEDIATE_TERMINAL_BUILD":
        if (
            payload.get("release_gate_pass") is not False
            or payload.get("publishable_task_b_release") is not False
            or "intermediate" not in str(payload.get("interpretation", "")).lower()
        ):
            raise ReportContractError("Intermediate source-attestation status is malformed")
        return dict(payload)
    if status != "VERIFIED_BOUNDED_RETROSPECTIVE_SOURCE_IDENTITY":
        raise ReportContractError("Tumor-REF source identity attestation status drift")
    if (
        payload.get("release_gate_pass") is not True
        or payload.get("publishable_task_b_release") is not True
        or payload.get("audit_class")
        != "bounded_retrospective_source_file_identity"
        or payload.get("source_roles") != ["analyzer", "focal_alt_cluster_lib"]
    ):
        raise ReportContractError("Tumor-REF source identity release gate is not passing")
    source_sha = payload.get("source_sha256")
    if not isinstance(source_sha, Mapping) or set(source_sha) != {
        "analyzer",
        "focal_alt_cluster_lib",
    }:
        raise ReportContractError("Tumor-REF source identity SHA role set drift")
    if any(
        not re.fullmatch(r"[0-9a-f]{64}", str(source_sha[role]))
        for role in source_sha
    ):
        raise ReportContractError("Tumor-REF source identity SHA is malformed")
    limitation = str(payload.get("limitation", ""))
    if "not a prelaunch lock" not in limitation or "environment attestation" not in limitation:
        raise ReportContractError("Tumor-REF source identity limitation disclosure drift")
    receipt_path = verify_embedded_artifact(
        payload.get("receipt"), "tumor-REF bounded source identity receipt"
    )
    snapshot_path = verify_embedded_artifact(
        payload.get("snapshot"), "tumor-REF during-execution source snapshot"
    )
    result = dict(payload)
    result["_receipt_path"] = receipt_path
    result["_snapshot_path"] = snapshot_path
    return result


def validate_technical_replication(final: Mapping[str, Any]) -> dict[str, Any]:
    payload = final.get("technical_replication")
    if not isinstance(payload, Mapping):
        raise ReportContractError("Final dataset lacks HCC technical replication")
    numerator = to_int(payload.get("numerator"), field="technical numerator")
    denominator = to_int(payload.get("denominator"), field="technical denominator")
    observed_ratio = optional_float(payload.get("ratio"), field="technical ratio")
    if numerator > denominator or observed_ratio != ratio(numerator, denominator):
        raise ReportContractError("HCC technical replication ratio does not reconcile")
    status = str(payload.get("status"))
    expected_status = (
        "ANY_CONCORDANT_EXACT_PAIR_OBSERVED"
        if numerator
        else "NO_CONCORDANT_EXACT_PAIR_OBSERVED"
        if denominator
        else "NOT_EVALUABLE"
    )
    if status != expected_status:
        raise ReportContractError("HCC technical replication status does not reconcile")
    if (
        payload.get("biological_n") != 1
        or payload.get("independent_biological_replication_n") != 0
        or payload.get("replication_claim_status") != "NOT_EVALUABLE_BIOLOGICAL_N1"
        or payload.get("inferential_confidence_interval") is not None
        or payload.get("pair_independence_assumption_met") is not False
        or payload.get("required_for_b1") is not False
    ):
        raise ReportContractError("HCC technical replication biological interpretation drift")
    biological = final.get("biological_replication")
    if not isinstance(biological, Mapping) or biological.get("status") != "NOT_RUN":
        raise ReportContractError("Biological replication must remain NOT_RUN")
    return {
        "comparison": "HCC1395 vs HCC1395_DORADO",
        "status": status,
        "numerator": numerator,
        "denominator": denominator,
        "ratio": observed_ratio,
        "not_evaluable_one_platform_only": to_int(
            payload.get("not_evaluable_one_platform_only"),
            field="technical not_evaluable",
        ),
        "biological_n": 1,
        "independent_biological_replication_n": 0,
        "replication_claim_status": "NOT_EVALUABLE_BIOLOGICAL_N1",
        "denominator_definition": str(payload.get("denominator_definition", "")),
        "interpretation": (
            "descriptive technical concordance only; dependent pair opportunities and "
            "biological n=1 prohibit a replication PASS claim"
        ),
    }


def validate_focal_partner_truth_matrix(final: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = final.get("focal_partner_truth_matrix")
    if not isinstance(rows, list) or len(rows) != len(TRUTH_LABELS) ** 2:
        raise ReportContractError("Final focal x partner truth matrix must contain 9 rows")
    expected_order = [
        (focal_truth, partner_truth)
        for focal_truth in TRUTH_LABELS
        for partner_truth in TRUTH_LABELS
    ]
    normalized: list[dict[str, Any]] = []
    for row, expected_key in zip(rows, expected_order):
        if not isinstance(row, Mapping):
            raise ReportContractError("Focal x partner truth matrix contains a non-object row")
        observed_key = (
            str(row.get("focal_truth_label")),
            str(row.get("partner_truth_label")),
        )
        if observed_key != expected_key:
            raise ReportContractError("Focal x partner truth matrix order/content drift")
        n_all = to_int(row.get("n_all_pair_rows"), field=f"truth matrix {observed_key}/all")
        n_formal = to_int(
            row.get("n_g1_formal_pair_rows"), field=f"truth matrix {observed_key}/formal"
        )
        if n_all < 0 or n_formal < 0 or n_formal > n_all:
            raise ReportContractError("Focal x partner truth matrix count drift")
        normalized.append(
            {
                "focal_truth_label": observed_key[0],
                "partner_truth_label": observed_key[1],
                "n_all_pair_rows": n_all,
                "n_g1_formal_pair_rows": n_formal,
            }
        )
    counts = final.get("counts")
    if not isinstance(counts, Mapping) or sum(
        row["n_all_pair_rows"] for row in normalized
    ) != to_int(counts.get("cooccurrence_pair_rows"), field="cooccurrence_pair_rows"):
        raise ReportContractError("Focal x partner truth matrix does not reconcile to pair rows")
    return normalized


def json_cell(value: Any, *, field: str, default: Any) -> Any:
    if value is None or str(value).strip() in {"", "NA", "N/A", "null", "None"}:
        return default
    if isinstance(value, (dict, list)):
        return value
    try:
        return json.loads(str(value))
    except json.JSONDecodeError as error:
        raise ReportContractError(f"Invalid JSON in {field}: {value!r}") from error


def pair_rank(row: Mapping[str, Any]) -> tuple[Any, ...]:
    p_value = optional_float(
        row.get("endpoint_a_p_fixed_margins_exact"), field="endpoint_a_p_fixed_margins_exact"
    )
    effect = optional_float(row.get("endpoint_a_cramers_v"), field="endpoint_a_cramers_v")
    delta = optional_float(
        row.get("endpoint_a_delta_alt_fraction"), field="endpoint_a_delta_alt_fraction"
    )
    if p_value is None or effect is None or delta is None:
        raise ReportContractError("Pair ranking requires exact p, Cramer's V, and delta ALT fraction")
    return (
        p_value,
        -effect,
        -abs(delta),
        str(row["sample"]),
        str(row["focal_chrom"]),
        int(row["focal_pos"]),
        int(row["partner_pos"]),
        str(row["partner_ref"]),
        str(row["partner_alt"]),
    )


def well_explained_pair_rank(
    row: Mapping[str, Any], exact_effect_rank: tuple[Any, ...]
) -> tuple[Any, ...]:
    state_counts = json_cell(
        row.get("endpoint_b_state_counts"), field="endpoint_b_state_counts", default={}
    )
    if not isinstance(state_counts, Mapping):
        raise ReportContractError("endpoint_b_state_counts must be an object")
    called_depth = sum(
        to_int(state_counts.get(state, 0), field=f"well-explained four-state {state}")
        for state in ("RR", "AR", "RA", "AA")
    )
    relation = str(row.get("endpoint_b_relation_compatibility", "")).strip()
    relation_available = relation not in {"", "NA", "N/A", "NOT_EVALUABLE", "UNTESTABLE"}
    conditional_pass = to_bool(
        row.get("endpoint_a_conditional_sensitivity_pass"),
        field="endpoint_a_conditional_sensitivity_pass",
    )
    formal_pair = to_bool(
        row.get("endpoint_a_formal_pair_by_confirmed"),
        field="endpoint_a_formal_pair_by_confirmed",
    )
    return (
        0 if conditional_pass else 1,
        0 if formal_pair else 1,
        0 if relation_available else 1,
        -called_depth,
        *exact_effect_rank,
    )


def pair_display(row: Mapping[str, Any]) -> tuple[str, str]:
    focal = site_key(row, pair=True)
    focal_label = f"{focal[0]} {focal[1]}:{focal[2]} {focal[3]}>{focal[4]}"
    partner_label = (
        f"{row['partner_chrom']}:{int(row['partner_pos'])} "
        f"{str(row['partner_ref']).upper()}>{str(row['partner_alt']).upper()}"
    )
    return focal_label, partner_label


def case_registry_row(
    *,
    view_order: int,
    view_role: str,
    selection_definition: str,
    selection_status: str,
    pair: Mapping[str, Any] | None,
    preferred_target: str = "N/A",
    negative_result: bool = False,
) -> dict[str, Any]:
    focal_site, partner_site = (
        pair_display(pair) if pair is not None else ("N/A", "N/A")
    )
    return {
        "view_order": view_order,
        "view_role": view_role,
        "selection_definition": selection_definition,
        "selection_status": selection_status,
        "preferred_target": preferred_target,
        "selected_focal_site": focal_site,
        "selected_partner_site": partner_site,
        "negative_result": negative_result,
    }


def choose_case_pair(
    pair_path: Path,
    candidates: Sequence[Mapping[str, Any]],
) -> tuple[
    dict[str, str] | None,
    str,
    list[dict[str, Any]],
    dict[str, dict[str, str] | None],
]:
    required = {
        "sample",
        "focal_chrom",
        "focal_pos",
        "focal_ref",
        "focal_alt",
        "partner_chrom",
        "partner_pos",
        "partner_ref",
        "partner_alt",
        "endpoint_a_groups",
        "endpoint_a_table",
        "endpoint_a_p_fixed_margins_exact",
        "endpoint_a_q_global_by",
        "endpoint_a_cramers_v",
        "endpoint_a_delta_alt_fraction",
        "endpoint_a_conditional_status",
        "endpoint_a_p_conditional_perm",
        "endpoint_a_permutations",
        "endpoint_a_permutable",
        "endpoint_a_conditional_sensitivity_pass",
        "endpoint_a_formal_pair_by_confirmed",
        "callability_testable",
        "callability_q_global_by",
        "callability_cramers_v",
        "callability_noncallable_core_reads",
        "callability_gate_status",
        "callability_gate_pass",
        "endpoint_b_state_counts",
        "endpoint_b_n_called_depth",
        "endpoint_b_error_ceiling",
        "endpoint_b_error_model_confidence",
        "endpoint_b_familywise_confidence",
        "endpoint_b_relation_family_size",
        "endpoint_b_multiplicity_method",
        "endpoint_b_minimum_zero_violation_depth",
        "endpoint_b_focal_ancestor_violation_p_exact",
        "endpoint_b_focal_ancestor_violation_upper_bound",
        "endpoint_b_focal_ancestor_violation_threshold",
        "endpoint_b_focal_ancestor_violation_status",
        "endpoint_b_partner_ancestor_violation_p_exact",
        "endpoint_b_partner_ancestor_violation_upper_bound",
        "endpoint_b_partner_ancestor_violation_threshold",
        "endpoint_b_partner_ancestor_violation_status",
        "endpoint_b_branching_violation_p_exact",
        "endpoint_b_branching_violation_upper_bound",
        "endpoint_b_branching_violation_threshold",
        "endpoint_b_branching_violation_status",
        "endpoint_b_complete_four_state_testable",
        "endpoint_b_relation_compatibility",
        "endpoint_b_compatible_relation_models",
        "endpoint_b_n_compatible_relation_models",
        "focal_truth_label",
        "partner_truth_label",
        "focal_ssnv_branch",
        "partner_ssnv_branch",
        "focal_component_id",
        "partner_component_id",
        "topology_scope",
        "topology_region",
        "topology_order_status",
        "topology_claim_guardrail",
    }
    candidate_keys = {site_key(row) for row in candidates}
    candidate_best: dict[str, str] | None = None
    candidate_rank: tuple[Any, ...] | None = None
    canonical_best: dict[str, str] | None = None
    canonical_rank: tuple[Any, ...] | None = None
    extreme_best: dict[str, str] | None = None
    extreme_rank: tuple[Any, ...] | None = None
    explained_best: dict[str, str] | None = None
    explained_rank: tuple[Any, ...] | None = None
    negative_best: dict[str, str] | None = None
    negative_rank: tuple[Any, ...] | None = None
    for row in iter_tsv(pair_path, required, "cooccurrence pair table"):
        try:
            rank = pair_rank(row)
        except ReportContractError:
            continue
        focal = site_key(row, pair=True)
        if extreme_rank is None or rank < extreme_rank:
            extreme_best = row
            extreme_rank = rank
        explained = well_explained_pair_rank(row, rank)
        if explained_rank is None or explained < explained_rank:
            explained_best = row
            explained_rank = explained
        if focal[:3] == CANONICAL_ORACLE:
            if canonical_rank is None or rank < canonical_rank:
                canonical_best = row
                canonical_rank = rank
        if focal in candidate_keys and to_bool(
            row.get("endpoint_a_formal_pair_by_confirmed"),
            field="endpoint_a_formal_pair_by_confirmed",
        ):
            if candidate_rank is None or rank < candidate_rank:
                candidate_best = row
                candidate_rank = rank
        conditional_evaluable_fail = (
            to_bool(row.get("endpoint_a_permutable"), field="endpoint_a_permutable")
            and to_int(row.get("endpoint_a_permutations"), field="endpoint_a_permutations")
            == 999
            and not to_bool(
                row.get("endpoint_a_conditional_sensitivity_pass"),
                field="endpoint_a_conditional_sensitivity_pass",
            )
        )
        if conditional_evaluable_fail and not to_bool(
            row.get("endpoint_a_formal_pair_by_confirmed"),
            field="endpoint_a_formal_pair_by_confirmed",
        ):
            if negative_rank is None or rank < negative_rank:
                negative_best = row
                negative_rank = rank

    if candidate_keys and candidate_best is not None:
        primary = candidate_best
        mode = "G2_PASS_BASE_CANDIDATE"
    elif not candidate_keys and extreme_best is not None:
        primary = extreme_best
        mode = "NON_CONFIRMING_FALLBACK"
    else:
        primary = None
        mode = "NOT_APPLICABLE_NO_ELIGIBLE_PAIR"

    oracle_text = f"{CANONICAL_ORACLE[0]} {CANONICAL_ORACLE[1]}:{CANONICAL_ORACLE[2]}"
    if canonical_best is not None:
        canonical_status = "SELECTED_PRE_REGISTERED_ORACLE"
    else:
        canonical_status = "TARGET_UNAVAILABLE_NO_SUBSTITUTION"
    selected_pairs = {
        "primary_report_case": primary,
        "canonical_pre_registered_oracle": canonical_best,
        "extreme_global_exact_effect": extreme_best,
        "well_explained": explained_best,
        "evaluable_statistical_negative": negative_best,
    }
    registry = [
        case_registry_row(
            view_order=1,
            view_role="aggregate",
            selection_definition=(
                "All 469,849 dataset-site rows and claim-specific denominators; no focal pair selection."
            ),
            selection_status="AVAILABLE_ALL_DATASET_SITES",
            pair=None,
        ),
        case_registry_row(
            view_order=2,
            view_role="canonical_pre_registered_oracle",
            selection_definition=(
                "Prefer the pre-registered HCC1395_DORADO chr5:750311 focal site; "
                "otherwise use the deterministic global exact/effect fallback or N/A."
            ),
            selection_status=canonical_status,
            pair=canonical_best,
            preferred_target=oracle_text,
        ),
        case_registry_row(
            view_order=3,
            view_role="extreme_global_exact_effect",
            selection_definition=(
                "Minimum fixed-margins exact p, then maximum Cramer's V, then maximum absolute delta ALT fraction."
            ),
            selection_status=(
                "SELECTED_GLOBAL_EXTREME"
                if extreme_best is not None
                else "NOT_APPLICABLE_NO_ELIGIBLE_PAIR"
            ),
            pair=extreme_best,
        ),
        case_registry_row(
            view_order=4,
            view_role="well_explained",
            selection_definition=(
                "Prefer conditional/formal and interpretable four-state rows, then greater called depth, "
                "then the deterministic exact/effect rank."
            ),
            selection_status=(
                "SELECTED_WELL_EXPLAINED"
                if explained_best is not None
                else "NOT_APPLICABLE_NO_ELIGIBLE_PAIR"
            ),
            pair=explained_best,
        ),
        case_registry_row(
            view_order=5,
            view_role="evaluable_statistical_negative",
            selection_definition=(
                "Require an executed 999-permutation conditional test that failed the "
                "conditional sensitivity endpoint and did not form a formal pair."
            ),
            selection_status=(
                "SELECTED_EVALUATED_ENDPOINT_FAIL"
                if negative_best is not None
                else "NOT_APPLICABLE_NO_EVALUATED_FAIL_PAIR"
            ),
            pair=negative_best,
            negative_result=negative_best is not None,
        ),
    ]
    return primary, mode, registry, selected_pairs


def find_site_row(path: Path, key: SiteKey) -> dict[str, str]:
    required = {
        "sample",
        "chrom",
        "pos",
        "ref",
        "alt",
        "truth_label",
        "ssnv_branch",
        "component_id",
        "modal_assignment_ari_min",
        "hp_axis_confound",
        "technical_axis_confound",
        "residual_unexplained_multigroup",
        "joint_signature_n_complete_reads",
        "joint_signature_testable",
        "joint_signature_groups",
        "joint_signature_categories",
        "joint_signature_table",
        "joint_signature_conditional_status",
        "joint_signature_p_conditional_perm",
        "joint_signature_permutations",
        "joint_signature_permutable",
        "joint_signature_sensitivity_pass",
        "joint_signature_global_fdr_family_status",
        "joint_signature_q_global_by",
        "joint_signature_global_by_discovery",
        "multi_marker_molecular_haplotype_base_candidate",
    }
    found: dict[str, str] | None = None
    for row in iter_tsv(path, required, "cooccurrence site table"):
        if site_key(row) != key:
            continue
        if found is not None:
            raise ReportContractError(f"Duplicate selected cooccurrence site row: {key}")
        found = row
    if found is None:
        raise ReportContractError(f"Selected case is missing from cooccurrence site table: {key}")
    return found


def matrix_rows(
    groups: Any,
    categories: Any,
    table: Any,
    *,
    group_field: str,
    category_field: str,
    status: str,
) -> list[dict[str, Any]]:
    if not isinstance(groups, list) or not isinstance(categories, list) or not isinstance(table, list):
        return [
            {
                "row_order": 1,
                group_field: "NOT_EVALUABLE",
                category_field: "NOT_EVALUABLE",
                "count": None,
                "status": status,
            }
        ]
    if len(table) != len(groups) or any(not isinstance(row, list) or len(row) != len(categories) for row in table):
        raise ReportContractError("Case matrix dimensions do not reconcile")
    rows: list[dict[str, Any]] = []
    for group_index, group in enumerate(groups):
        for category_index, category in enumerate(categories):
            count = to_int(table[group_index][category_index], field="case matrix count")
            if count < 0:
                raise ReportContractError("Case matrix contains a negative count")
            rows.append(
                {
                    "row_order": len(rows) + 1,
                    group_field: str(group),
                    category_field: str(category),
                    "count": count,
                    "status": status,
                }
            )
    return rows


def four_state_model_rows(
    pair: Mapping[str, Any], *, view_role: str
) -> list[dict[str, Any]]:
    state_counts = json_cell(
        pair.get("endpoint_b_state_counts"), field="endpoint_b_state_counts", default={}
    )
    if not isinstance(state_counts, Mapping) or set(state_counts) != {
        "RR",
        "AR",
        "RA",
        "AA",
        "O",
        "X",
    }:
        raise ReportContractError("Four-state model requires exact RR/AR/RA/AA/O/X counts")
    counts = {
        state: to_int(state_counts[state], field=f"four-state model {state}")
        for state in state_counts
    }
    called_depth = sum(counts[state] for state in ("RR", "AR", "RA", "AA"))
    declared_depth = to_int(
        pair.get("endpoint_b_n_called_depth"), field="endpoint_b_n_called_depth"
    )
    if called_depth != declared_depth:
        raise ReportContractError("Four-state called-depth reconciliation failed")
    error_ceiling = optional_float(
        pair.get("endpoint_b_error_ceiling"), field="endpoint_b_error_ceiling"
    )
    confidence = optional_float(
        pair.get("endpoint_b_error_model_confidence"),
        field="endpoint_b_error_model_confidence",
    )
    familywise_confidence = optional_float(
        pair.get("endpoint_b_familywise_confidence"),
        field="endpoint_b_familywise_confidence",
    )
    family_size = to_int(
        pair.get("endpoint_b_relation_family_size"),
        field="endpoint_b_relation_family_size",
    )
    multiplicity_method = str(pair.get("endpoint_b_multiplicity_method"))
    minimum_zero_violation_depth = to_int(
        pair.get("endpoint_b_minimum_zero_violation_depth"),
        field="endpoint_b_minimum_zero_violation_depth",
    )
    if error_ceiling is None or confidence is None or familywise_confidence is None:
        raise ReportContractError("Four-state fixed-error model parameters are missing")
    if (
        not math.isclose(
            confidence, FOUR_STATE_PER_RELATION_CONFIDENCE, rel_tol=0, abs_tol=1e-12
        )
        or not math.isclose(
            familywise_confidence,
            FOUR_STATE_FAMILYWISE_CONFIDENCE,
            rel_tol=0,
            abs_tol=1e-12,
        )
        or family_size != FOUR_STATE_RELATION_FAMILY_SIZE
        or multiplicity_method != FOUR_STATE_MULTIPLICITY_METHOD
        or minimum_zero_violation_depth != FOUR_STATE_MINIMUM_ZERO_VIOLATION_DEPTH
    ):
        raise ReportContractError("Four-state simultaneous-confidence contract drift")
    compatible_models = json_cell(
        pair.get("endpoint_b_compatible_relation_models"),
        field="endpoint_b_compatible_relation_models",
        default=[],
    )
    if not isinstance(compatible_models, list):
        raise ReportContractError("Four-state compatible model list is malformed")
    if len(compatible_models) != to_int(
        pair.get("endpoint_b_n_compatible_relation_models"),
        field="endpoint_b_n_compatible_relation_models",
    ):
        raise ReportContractError("Four-state compatible model count drift")
    complete = to_bool(
        pair.get("endpoint_b_complete_four_state_testable"),
        field="endpoint_b_complete_four_state_testable",
    )
    specs = (
        (
            "FOCAL_ANCESTOR",
            counts["RA"],
            counts["RA"] + counts["AA"],
            counts["AR"] >= 3 and counts["AA"] >= 3,
            "endpoint_b_focal_ancestor_violation",
        ),
        (
            "PARTNER_ANCESTOR",
            counts["AR"],
            counts["AR"] + counts["AA"],
            counts["RA"] >= 3 and counts["AA"] >= 3,
            "endpoint_b_partner_ancestor_violation",
        ),
        (
            "BRANCHING",
            counts["AA"],
            counts["AR"] + counts["RA"] + counts["AA"],
            counts["AR"] >= 3 and counts["RA"] >= 3,
            "endpoint_b_branching_violation",
        ),
    )
    rows: list[dict[str, Any]] = []
    for order, (model, violations, denominator, split_support, prefix) in enumerate(
        specs, 1
    ):
        threshold = optional_float(
            pair.get(f"{prefix}_threshold"), field=f"{prefix}_threshold"
        )
        if threshold is None or not math.isclose(
            threshold, error_ceiling, rel_tol=0, abs_tol=1e-15
        ):
            raise ReportContractError("Four-state model threshold/error-ceiling drift")
        rows.append(
            {
                "view_role": view_role,
                "model_order": order,
                "model": model,
                "violations": violations,
                "denominator": denominator,
                "violation_rate": ratio(violations, denominator),
                "p_exact_greater": optional_float(
                    pair.get(f"{prefix}_p_exact"), field=f"{prefix}_p_exact"
                ),
                "upper_per_relation": optional_float(
                    pair.get(f"{prefix}_upper_bound"),
                    field=f"{prefix}_upper_bound",
                ),
                "error_ceiling": error_ceiling,
                "per_relation_confidence": confidence,
                "familywise_confidence": familywise_confidence,
                "relation_family_size": family_size,
                "multiplicity_method": multiplicity_method,
                "minimum_zero_violation_depth": minimum_zero_violation_depth,
                "status": str(pair.get(f"{prefix}_status")),
                "complete_four_state_depth": complete,
                "split_support": split_support,
                "compatible": model in compatible_models,
                "relation_status": str(
                    pair.get("endpoint_b_relation_compatibility")
                ),
                "claim_guardrail": (
                    "Fixed-error-model compatibility is pairwise post hoc context, "
                    "not proof of cellular ancestry or unique mutation order."
                ),
            }
        )
    return rows


def build_case_view_summary(
    *,
    view_order: int,
    view_role: str,
    selection_status: str,
    pair: Mapping[str, Any] | None,
    site: Mapping[str, Any] | None,
    candidate: Mapping[str, Any] | None,
    negative_result: bool,
) -> dict[str, Any]:
    if pair is None or site is None:
        return {
            "view_order": view_order,
            "view_role": view_role,
            "selection_status": selection_status,
            "evidence_status": "NOT_APPLICABLE",
            "negative_result": False,
            "negative_result_scope": "N/A is not a statistical or biological negative",
            "focal_site": "N/A",
            "partner_site": "N/A",
            "focal_truth_label": "N/A",
            "partner_truth_label": "N/A",
            "truth_pair": "N/A / N/A",
            "g2_status": "NOT_APPLICABLE",
            "exact_p": None,
            "global_by_q": None,
            "cramers_v": None,
            "delta_alt_fraction": None,
            "conditional_p": None,
            "conditional_permutations": None,
            "conditional_status": "NOT_APPLICABLE",
            "conditional_sensitivity_pass": None,
            "callability_gate_status": "NOT_APPLICABLE",
            "callability_gate_pass": None,
            "noncallable_core_reads": None,
            "modal_assignment_ari_min": None,
            "hp_axis_confound": None,
            "technical_axis_confound": None,
            "residual_unexplained_multigroup": None,
            "focal_branch": "N/A",
            "partner_branch": "N/A",
            "focal_component": "N/A",
            "partner_component": "N/A",
            "topology_scope": "N/A",
            "topology_region": "N/A",
            "topology_order_status": "N/A",
            "topology_claim_guardrail": "N/A",
            "normal_focal_status": "NOT_APPLICABLE",
            "normal_focal_called_reads": None,
            "normal_focal_alt_reads": None,
            "normal_partner_status": "NOT_EVALUATED_BY_DESIGN",
            "pair_reused_across_views": False,
            "same_pair_view_roles": "N/A",
        }
    focal_label, partner_label = pair_display(pair)
    g2 = to_bool(
        site.get("multi_marker_molecular_haplotype_base_candidate"),
        field="multi_marker_molecular_haplotype_base_candidate",
    )
    if negative_result and (
        not to_bool(pair.get("endpoint_a_permutable"), field="endpoint_a_permutable")
        or to_int(pair.get("endpoint_a_permutations"), field="endpoint_a_permutations")
        != 999
        or to_bool(
            pair.get("endpoint_a_conditional_sensitivity_pass"),
            field="endpoint_a_conditional_sensitivity_pass",
        )
    ):
        raise ReportContractError("Statistical negative case is not an evaluated endpoint FAIL")
    return {
        "view_order": view_order,
        "view_role": view_role,
        "selection_status": selection_status,
        "evidence_status": "AVAILABLE",
        "negative_result": negative_result,
        "negative_result_scope": (
            "evaluated conditional endpoint FAIL; not biological absence"
            if negative_result
            else "not designated as the evaluated statistical negative view"
        ),
        "focal_site": focal_label,
        "partner_site": partner_label,
        "focal_truth_label": str(pair.get("focal_truth_label")),
        "partner_truth_label": str(pair.get("partner_truth_label")),
        "truth_pair": (
            f"{pair.get('focal_truth_label')} / {pair.get('partner_truth_label')}"
        ),
        "g2_status": "PASS" if g2 else "FAIL",
        "exact_p": optional_float(
            pair.get("endpoint_a_p_fixed_margins_exact"), field="case exact p"
        ),
        "global_by_q": optional_float(
            pair.get("endpoint_a_q_global_by"), field="case BY q"
        ),
        "cramers_v": optional_float(
            pair.get("endpoint_a_cramers_v"), field="case Cramer's V"
        ),
        "delta_alt_fraction": optional_float(
            pair.get("endpoint_a_delta_alt_fraction"),
            field="case delta ALT fraction",
        ),
        "conditional_p": optional_float(
            pair.get("endpoint_a_p_conditional_perm"), field="case conditional p"
        ),
        "conditional_permutations": to_int(
            pair.get("endpoint_a_permutations"), field="case permutations"
        ),
        "conditional_status": str(pair.get("endpoint_a_conditional_status")),
        "conditional_sensitivity_pass": to_bool(
            pair.get("endpoint_a_conditional_sensitivity_pass"),
            field="case conditional sensitivity",
        ),
        "callability_gate_status": str(pair.get("callability_gate_status")),
        "callability_gate_pass": to_bool(
            pair.get("callability_gate_pass"), field="callability_gate_pass"
        ),
        "noncallable_core_reads": to_int(
            pair.get("callability_noncallable_core_reads"),
            field="callability_noncallable_core_reads",
        ),
        "modal_assignment_ari_min": optional_float(
            site.get("modal_assignment_ari_min"), field="modal_assignment_ari_min"
        ),
        "hp_axis_confound": to_bool(
            site.get("hp_axis_confound"), field="hp_axis_confound"
        ),
        "technical_axis_confound": to_bool(
            site.get("technical_axis_confound"), field="technical_axis_confound"
        ),
        "residual_unexplained_multigroup": to_bool(
            site.get("residual_unexplained_multigroup"),
            field="residual_unexplained_multigroup",
        ),
        "focal_branch": str(pair.get("focal_ssnv_branch")),
        "partner_branch": str(pair.get("partner_ssnv_branch")),
        "focal_component": str(pair.get("focal_component_id") or "N/A"),
        "partner_component": str(pair.get("partner_component_id") or "N/A"),
        "topology_scope": str(pair.get("topology_scope") or "N/A"),
        "topology_region": str(pair.get("topology_region") or "N/A"),
        "topology_order_status": str(pair.get("topology_order_status") or "N/A"),
        "topology_claim_guardrail": str(pair.get("topology_claim_guardrail")),
        "normal_focal_status": str(
            candidate.get("normal_control_status")
            if candidate is not None
            else "NOT_RUN_NON_G2_VIEW"
        ),
        "normal_focal_called_reads": (
            candidate.get("normal_called_reads") if candidate is not None else None
        ),
        "normal_focal_alt_reads": (
            candidate.get("normal_alt_reads") if candidate is not None else None
        ),
        "normal_partner_status": "NOT_EVALUATED_BY_DESIGN",
        "pair_reused_across_views": False,
        "same_pair_view_roles": "",
    }


def build_case_evidence(
    pair: Mapping[str, Any] | None,
    site: Mapping[str, Any] | None,
    mode: str,
    selection_registry: Sequence[Mapping[str, Any]],
    selected_pairs: Mapping[str, Mapping[str, Any] | None],
    selected_sites: Mapping[str, Mapping[str, Any] | None],
    candidates: Sequence[Mapping[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    if pair is None or site is None:
        if mode != "NOT_APPLICABLE_NO_ELIGIBLE_PAIR":
            raise ReportContractError("Missing pair/site for an applicable report case")
        not_applicable = "NOT_APPLICABLE_NO_ELIGIBLE_PAIR"
        return {
            "case_selection_registry": [dict(row) for row in selection_registry],
            "case_summary": [
                {
                    "case_rank": 1,
                    "case_role": "primary_report_case",
                    "selection_mode": not_applicable,
                    "evidence_status": "NOT_APPLICABLE",
                    "negative_result": False,
                    "focal_site": "N/A",
                    "partner_site": "N/A",
                    "truth_label": "N/A",
                    "g2_status": "NOT_APPLICABLE",
                    "exact_p": None,
                    "global_by_q": None,
                    "cramers_v": None,
                    "delta_alt_fraction": None,
                    "conditional_status": "NOT_APPLICABLE",
                    "conditional_sensitivity_pass": None,
                    "four_state_relation": "NOT_APPLICABLE",
                    "joint_signature_n_complete_reads": None,
                    "joint_signature_status": "NOT_APPLICABLE",
                    "joint_signature_sensitivity_pass": None,
                    "interpretation": (
                        "No sortable exact-p/effect pair is available. This presentation state is N/A, "
                        "not a negative biological result."
                    ),
                }
            ],
            "case_group_partner": [
                {
                    "row_order": 1,
                    "methyl_group": "N/A",
                    "partner_call": "N/A",
                    "count": None,
                    "status": not_applicable,
                }
            ],
            "case_four_state": [
                {
                    "state_order": 1,
                    "state": "N/A",
                    "count": None,
                    "relation_status": not_applicable,
                }
            ],
            "case_joint_signature": [
                {
                    "row_order": 1,
                    "methyl_group": "N/A",
                    "joint_signature": "N/A",
                    "count": None,
                    "status": not_applicable,
                }
            ],
            "case_four_state_models": [
                {
                    "view_role": "primary_report_case",
                    "model_order": 1,
                    "model": "N/A",
                    "violations": None,
                    "denominator": None,
                    "violation_rate": None,
                    "p_exact_greater": None,
                    "upper_per_relation": None,
                    "error_ceiling": None,
                    "per_relation_confidence": None,
                    "familywise_confidence": None,
                    "relation_family_size": None,
                    "multiplicity_method": "NOT_APPLICABLE",
                    "minimum_zero_violation_depth": None,
                    "status": not_applicable,
                    "complete_four_state_depth": None,
                    "split_support": None,
                    "compatible": None,
                    "relation_status": not_applicable,
                    "claim_guardrail": "N/A is not a biological negative result.",
                }
            ],
            "case_view_summary": [
                build_case_view_summary(
                    view_order=index,
                    view_role=role,
                    selection_status="NOT_APPLICABLE_NO_ELIGIBLE_PAIR",
                    pair=None,
                    site=None,
                    candidate=None,
                    negative_result=False,
                )
                for index, role in enumerate(selected_pairs, 1)
            ],
            "case_view_group_partner": [],
            "case_view_four_state": [],
            "case_view_four_state_models": [],
        }
    focal = site_key(pair, pair=True)
    partner = (
        str(pair["partner_chrom"]),
        int(pair["partner_pos"]),
        str(pair["partner_ref"]).upper(),
        str(pair["partner_alt"]).upper(),
    )
    g2_site = to_bool(
        site.get("multi_marker_molecular_haplotype_base_candidate"),
        field="multi_marker_molecular_haplotype_base_candidate",
    )
    joint_global_by = optional_float(
        site.get("joint_signature_q_global_by"),
        field="joint_signature_q_global_by",
    )
    joint_global_by_discovery = to_bool(
        site.get("joint_signature_global_by_discovery"),
        field="joint_signature_global_by_discovery",
    )
    if g2_site and (
        site.get("joint_signature_global_fdr_family_status")
        != JOINT_GLOBAL_FDR_FAMILY
        or joint_global_by is None
        or joint_global_by > 0.05
        or not joint_global_by_discovery
    ):
        raise ReportContractError("Selected G2 case does not satisfy the global-BY gate")
    if (mode == "G2_PASS_BASE_CANDIDATE") != g2_site:
        raise ReportContractError("Selected case mode does not reconcile to the site G2 status")
    candidate_by_key = {site_key(row): row for row in candidates}
    primary_view = build_case_view_summary(
        view_order=1,
        view_role="primary_report_case",
        selection_status=mode,
        pair=pair,
        site=site,
        candidate=candidate_by_key.get(focal),
        negative_result=False,
    )
    case_summary = {
        "case_rank": 1,
        "case_role": "primary_report_case",
        "selection_mode": mode,
        **{key: value for key, value in primary_view.items() if key != "view_order"},
        "truth_label": str(site.get("truth_label")),
        "four_state_relation": str(pair.get("endpoint_b_relation_compatibility")),
        "joint_signature_n_complete_reads": to_int(
            site.get("joint_signature_n_complete_reads"), field="joint complete reads"
        ),
        "joint_signature_status": str(site.get("joint_signature_conditional_status")),
        "joint_signature_p": optional_float(
            site.get("joint_signature_p_conditional_perm"),
            field="joint signature conditional p",
        ),
        "joint_signature_permutations": to_int(
            site.get("joint_signature_permutations"),
            field="joint signature permutations",
        ),
        "joint_signature_permutable": to_bool(
            site.get("joint_signature_permutable"),
            field="joint signature permutable",
        ),
        "joint_signature_sensitivity_pass": to_bool(
            site.get("joint_signature_sensitivity_pass"), field="joint signature sensitivity"
        ),
        "joint_signature_q_global_by": joint_global_by,
        "joint_signature_global_by_discovery": joint_global_by_discovery,
        "interpretation": (
            CLAIM_NAMES["G2"]
            if g2_site
            else "non-confirming exact-p/effect-ranked case; not a G2 result"
        ),
    }
    group_partner = matrix_rows(
        json_cell(pair.get("endpoint_a_groups"), field="endpoint_a_groups", default=[]),
        ["R", "A"],
        json_cell(pair.get("endpoint_a_table"), field="endpoint_a_table", default=[]),
        group_field="methyl_group",
        category_field="partner_call",
        status=str(pair.get("endpoint_a_conditional_status")),
    )
    state_counts = json_cell(
        pair.get("endpoint_b_state_counts"), field="endpoint_b_state_counts", default={}
    )
    if not isinstance(state_counts, Mapping):
        raise ReportContractError("endpoint_b_state_counts must be an object")
    four_state = [
        {
            "state_order": order,
            "state": state,
            "count": to_int(state_counts.get(state, 0), field=f"four-state {state}"),
            "relation_status": str(pair.get("endpoint_b_relation_compatibility")),
        }
        for order, state in enumerate(("RR", "AR", "RA", "AA", "O", "X"), 1)
    ]
    joint_signature = matrix_rows(
        json_cell(site.get("joint_signature_groups"), field="joint_signature_groups", default=[]),
        json_cell(
            site.get("joint_signature_categories"),
            field="joint_signature_categories",
            default=[],
        ),
        json_cell(site.get("joint_signature_table"), field="joint_signature_table", default=[]),
        group_field="methyl_group",
        category_field="joint_signature",
        status=str(site.get("joint_signature_conditional_status")),
    )
    registry_by_role = {
        str(row["view_role"]): row for row in selection_registry
    }
    view_order = {
        "primary_report_case": 1,
        "canonical_pre_registered_oracle": 2,
        "extreme_global_exact_effect": 3,
        "well_explained": 4,
        "evaluable_statistical_negative": 5,
    }
    case_view_summary: list[dict[str, Any]] = []
    case_view_group_partner: list[dict[str, Any]] = []
    case_view_four_state: list[dict[str, Any]] = []
    case_view_four_state_models: list[dict[str, Any]] = []
    pair_roles: dict[tuple[str, str], list[str]] = {}
    for role, selected_pair in selected_pairs.items():
        selected_site = selected_sites.get(role)
        registry_row = registry_by_role.get(role, {})
        selection_status = (
            mode if role == "primary_report_case" else str(registry_row.get("selection_status"))
        )
        negative = role == "evaluable_statistical_negative" and selected_pair is not None
        selected_key = site_key(selected_pair, pair=True) if selected_pair is not None else None
        summary_row = build_case_view_summary(
            view_order=view_order[role],
            view_role=role,
            selection_status=selection_status,
            pair=selected_pair,
            site=selected_site,
            candidate=(candidate_by_key.get(selected_key) if selected_key is not None else None),
            negative_result=negative,
        )
        case_view_summary.append(summary_row)
        if selected_pair is None:
            continue
        pair_identity = pair_display(selected_pair)
        pair_roles.setdefault(pair_identity, []).append(role)
        selected_group_rows = matrix_rows(
            json_cell(
                selected_pair.get("endpoint_a_groups"),
                field="endpoint_a_groups",
                default=[],
            ),
            ["R", "A"],
            json_cell(
                selected_pair.get("endpoint_a_table"),
                field="endpoint_a_table",
                default=[],
            ),
            group_field="methyl_group",
            category_field="partner_call",
            status=str(selected_pair.get("endpoint_a_conditional_status")),
        )
        case_view_group_partner.extend(
            {"view_role": role, **row} for row in selected_group_rows
        )
        selected_counts = json_cell(
            selected_pair.get("endpoint_b_state_counts"),
            field="endpoint_b_state_counts",
            default={},
        )
        if not isinstance(selected_counts, Mapping):
            raise ReportContractError("Selected view four-state counts are malformed")
        case_view_four_state.extend(
            {
                "view_role": role,
                "state_order": order,
                "state": state,
                "count": to_int(
                    selected_counts.get(state, 0), field=f"view four-state {state}"
                ),
                "relation_status": str(
                    selected_pair.get("endpoint_b_relation_compatibility")
                ),
            }
            for order, state in enumerate(("RR", "AR", "RA", "AA", "O", "X"), 1)
        )
        case_view_four_state_models.extend(
            four_state_model_rows(selected_pair, view_role=role)
        )
    for row in case_view_summary:
        if row["evidence_status"] != "AVAILABLE":
            continue
        roles = pair_roles.get((row["focal_site"], row["partner_site"]), [])
        row["pair_reused_across_views"] = len(roles) > 1
        row["same_pair_view_roles"] = ", ".join(roles)
    return {
        "case_selection_registry": [dict(row) for row in selection_registry],
        "case_summary": [case_summary],
        "case_group_partner": group_partner,
        "case_four_state": four_state,
        "case_joint_signature": joint_signature,
        "case_four_state_models": four_state_model_rows(
            pair, view_role="primary_report_case"
        ),
        "case_view_summary": case_view_summary,
        "case_view_group_partner": case_view_group_partner,
        "case_view_four_state": case_view_four_state,
        "case_view_four_state_models": case_view_four_state_models,
    }


def build_stratum_rows(metrics: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    strata = [
        *(('sample', sample, metrics["per_sample"][sample]) for sample in EXPECTED_DATASETS),
        *(("truth", label, metrics["truth_strata"][label]) for label in TRUTH_LABELS),
    ]
    for stratum_type, stratum, payload in strata:
        stratum_order = (
            EXPECTED_DATASETS.index(stratum) + 1
            if stratum_type == "sample"
            else TRUTH_LABELS.index(stratum) + 1
        )
        for claim_order, claim_id in enumerate(REPORT_STRATUM_CLAIMS, 1):
            metric = payload[claim_id]
            rows.append(
                {
                    "stratum_type": stratum_type,
                    "stratum_order": stratum_order,
                    "stratum": stratum,
                    "claim_order": claim_order,
                    "claim_id": claim_id,
                    "claim_name": CLAIM_NAMES[claim_id],
                    "numerator": metric["numerator"],
                    "denominator": metric["denominator"],
                    "ratio": metric["ratio"],
                    "not_evaluable": metric["not_evaluable"],
                    "not_run": metric["not_run"],
                    "population": metric["population"],
                    "denominator_definition": metric["denominator_definition"],
                }
            )
    return rows


def build_claim_chart_rows(claim_rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in claim_rows:
        for unit, numerator_field, denominator_field, ratio_field in (
            ("dataset sites", "dataset_numerator", "dataset_denominator", "dataset_ratio"),
            ("biological sites", "biological_numerator", "biological_denominator", "biological_ratio"),
        ):
            if row["claim_id"] == "M1":
                if unit == "dataset sites":
                    series_label = "M1 FLAGGED / all dataset-sites"
                    metric_definition = "FLAGGED / all dataset-sites"
                else:
                    series_label = "M1 descriptive biological-site aggregation"
                    metric_definition = (
                        "FLAGGED / descriptive biological-site aggregation"
                    )
                interpretation = "operational aggregation; not biological prevalence"
            else:
                series_label = f"{unit}: PASS / claim-specific evaluable denominator"
                metric_definition = "PASS / claim-specific evaluable denominator"
                interpretation = "claim-specific evidence ratio"
            rows.append(
                {
                    "claim_order": row["claim_order"],
                    "claim_id": row["claim_id"],
                    "statistical_unit": unit,
                    "series_label": series_label,
                    "numerator": row[numerator_field],
                    "denominator": row[denominator_field],
                    "ratio": row[ratio_field],
                    "status": row["status"],
                    "metric_definition": metric_definition,
                    "interpretation": interpretation,
                }
            )
    return rows


def build_overview_case_chart_rows(
    claim_chart_rows: Sequence[Mapping[str, Any]],
    case_group_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    rows = [
        {
            "panel": "claim ladder",
            "category": str(row["claim_id"]),
            "series": str(row["series_label"]),
            "proportion": row["ratio"],
            "numerator": row["numerator"],
            "denominator": row["denominator"],
            "status": row["status"],
            "metric_definition": row["metric_definition"],
            "interpretation": row["interpretation"],
        }
        for row in claim_chart_rows
        if row["claim_id"] in REPORT_STRATUM_CLAIMS
    ]
    by_group: dict[str, list[Mapping[str, Any]]] = {}
    for row in case_group_rows:
        count = row.get("count")
        if count is None:
            continue
        by_group.setdefault(str(row["methyl_group"]), []).append(row)
    for methyl_group, group_rows in by_group.items():
        denominator = sum(to_int(row["count"], field="case group count") for row in group_rows)
        if denominator <= 0:
            continue
        for row in group_rows:
            numerator = to_int(row["count"], field="case group count")
            rows.append(
                {
                    "panel": "actual case",
                    "category": f"case {methyl_group}",
                    "series": f"partner {row['partner_call']}",
                    "proportion": ratio(numerator, denominator),
                    "numerator": numerator,
                    "denominator": denominator,
                    "status": row["status"],
                    "metric_definition": (
                        "partner R or A reads / R+A reads within methyl group"
                    ),
                    "interpretation": (
                        "molecular read composition; not cellular clone identity"
                    ),
                }
            )
    return rows


def build_flow_rows(
    layered: Mapping[str, Any],
    tag_audit: Mapping[str, Any],
    producer_bam_audit: Mapping[str, Any],
) -> list[dict[str, Any]]:
    contract = layered["analysis_contract"]
    production = layered["production_summary"]
    return [
        {
            "step_order": 1,
            "flow_step": "normalized ClairS raw/all",
            "status": "PASS",
            "scope": f"{production['passed_dataset_count']}/7 datasets",
            "contract": contract["longphase_input_contract"],
        },
        {
            "step_order": 2,
            "flow_step": "LongPhase-S recalibration",
            "status": "PASS",
            "scope": "same production run",
            "contract": str(contract.get("tagging_semantics", "LongPhase-S recalibration")),
        },
        {
            "step_order": 3,
            "flow_step": "LongPhase-S BAM-named output transport",
            "status": "PASS",
            "scope": f"{producer_bam_audit['named_fifo_count']}/7 named FIFO; 0 persisted BAM",
            "contract": "no regular tagged BAM created or overwritten; tags persisted only in sidecars",
        },
        {
            "step_order": 4,
            "flow_step": "LongPhase-S recalibrated FILTER=PASS tree backbone",
            "status": "PASS",
            "scope": "7/7 datasets",
            "contract": contract["tree_input_contract"],
        },
        {
            "step_order": 5,
            "flow_step": "chr1-22 biallelic sSNV frozen screen",
            "status": "PASS",
            "scope": format_count(EXPECTED_SITES),
            "contract": "same PASS backbone subset",
        },
        {
            "step_order": 6,
            "flow_step": "raw BAM MM/ML plus current terminal HP/PS join",
            "status": "PASS",
            "scope": format_count(
                to_int(
                    tag_audit["n_exact_hp_ps_site_read_joins"],
                    field="terminal HP/PS joins",
                )
            ),
            "contract": "same-run external sidecar; embedded tags ignored",
        },
    ]


def build_audit_rows(
    screen_summary: Mapping[str, Any],
    reconciliation: Mapping[str, Any],
    immutability: Mapping[str, Any],
    tree_audit: Mapping[str, Any],
    reference_audit: Mapping[str, Any],
    producer_bam_audit: Mapping[str, Any],
    m2_provenance: Mapping[str, Any],
    source_attestation: Mapping[str, Any],
) -> list[dict[str, Any]]:
    terminal = screen_summary["latest_hp_ps_terminal_join_audit"]
    return [
        {
            "audit_order": 1,
            "audit": "current screen schema and terminal HP/PS",
            "schema_version": screen_summary["schema_version"],
            "status": "PASS",
            "scope": format_count(to_int(terminal["n_sites"], field="terminal sites")),
            "key_result": f"{format_count(to_int(terminal['n_exact_hp_ps_site_read_joins'], field='terminal joins'))} exact site-read joins",
        },
        {
            "audit_order": 2,
            "audit": "producer BAM-output persistence",
            "schema_version": producer_bam_audit["schema_version"],
            "status": "PASS",
            "scope": f"{producer_bam_audit['named_fifo_count']}/7 named FIFO",
            "key_result": "persisted BAM=0; regular tagged BAM=0; no overwrite",
        },
        {
            "audit_order": 3,
            "audit": "post-run output reconciliation",
            "schema_version": reconciliation["schema_version"],
            "status": "PASS",
            "scope": format_count(to_int(reconciliation["totals"]["expected_sites"], field="reconciled sites")),
            "key_result": "reads/methylation/BERNOULLI exact",
        },
        {
            "audit_order": 4,
            "audit": "post-run frozen input immutability",
            "schema_version": immutability["schema_version"],
            "status": "PASS",
            "scope": f"{immutability['totals']['n_artifact_pass']}/{immutability['totals']['n_artifacts']} artifacts",
            "key_result": "no observed input mutation",
        },
        {
            "audit_order": 5,
            "audit": "latest tree-input identity",
            "schema_version": tree_audit["schema_version"],
            "status": "PASS",
            "scope": "7/7 datasets",
            "key_result": "same-run LongPhase-S recalibrated FILTER=PASS",
        },
        {
            "audit_order": 6,
            "audit": "extraction reference identity",
            "schema_version": reference_audit["schema_version"],
            "status": "PASS",
            "scope": "7/7 extraction receipts",
            "key_result": "full reference SHA-256 plus frozen chunk binding",
        },
        {
            "audit_order": 7,
            "audit": "M2 screen-producer source and downstream reclassification",
            "schema_version": M2_GATE_CONTRACT,
            "status": "PASS_WITH_PROVENANCE_LIMIT",
            "scope": "source-locked producer plus terminal gate",
            "key_result": (
                f"pinned analyzer {m2_provenance['screen_recovery_source_validation']['pinned_analyzer_sha256'][:12]}...; "
                "raw-read axis statistic not independently recomputed"
            ),
        },
        {
            "audit_order": 8,
            "audit": "tumor-REF bounded retrospective source identity",
            "schema_version": "1.1.0",
            "status": (
                "PASS"
                if source_attestation.get("release_gate_pass") is True
                else "INTERMEDIATE_NOT_RELEASE_READY"
            ),
            "scope": (
                "analyzer + focal_alt_cluster_lib"
                if source_attestation.get("release_gate_pass") is True
                else "receipt not included"
            ),
            "key_result": (
                "source identities unchanged and producer manifest hash-bound"
                if source_attestation.get("release_gate_pass") is True
                else "final Task Type B release requires post-run source receipt"
            ),
        },
    ]


def sql_literal(value: Any) -> str:
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if isinstance(value, float) and not math.isfinite(value):
            raise ReportContractError("Cannot serialize non-finite SQL value")
        return repr(value)
    return "'" + str(value).replace("'", "''") + "'"


def snapshot_sql(rows: Sequence[Mapping[str, Any]], columns: Sequence[str], relation: str) -> str:
    quoted_columns = ", ".join(f'"{column}"' for column in columns)
    if rows:
        values = ",\n  ".join(
            "(" + ", ".join(sql_literal(row.get(column)) for column in columns) + ")"
            for row in rows
        )
        sql = (
            f'WITH "{relation}" ({quoted_columns}) AS (\n  VALUES {values}\n)\n'
            f'SELECT * FROM "{relation}";'
        )
    else:
        projections = ", ".join(f'NULL AS "{column}"' for column in columns)
        sql = f"SELECT {projections} WHERE 0;"
    try:
        observed = sqlite3.connect(":memory:").execute(sql).fetchall()
    except sqlite3.Error as error:
        raise ReportContractError(f"Generated snapshot SQL is invalid for {relation}") from error
    if len(observed) != len(rows):
        raise ReportContractError(f"Generated snapshot SQL row count drift for {relation}")
    return sql


def source_spec(
    source_id: str,
    label: str,
    identity: Mapping[str, Any],
    *,
    rows: Sequence[Mapping[str, Any]] | None = None,
    columns: Sequence[str] | None = None,
    relation: str | None = None,
    description: str,
    filters: Sequence[str] = (),
    metric_definitions: Sequence[str] = (),
    tables_used: Sequence[str] = (),
    generated_at: str,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "id": source_id,
        "label": label,
        "path": identity["path"],
        "size_bytes": identity["size_bytes"],
        "sha256": identity["sha256"],
    }
    query: dict[str, Any] = {
        "engine": "SQLite",
        "language": "sql",
        "description": description,
        "executed_at": generated_at,
        "id": f"{source_id}-{str(identity['sha256'])[:12]}",
        "tables_used": list(tables_used or (identity["path"],)),
        "filters": list(filters),
        "metric_definitions": [
            *metric_definitions,
            f"Authoritative source SHA-256: {identity['sha256']}",
        ],
    }
    if rows is not None and columns is not None and relation is not None:
        query["sql"] = snapshot_sql(rows, columns, relation)
    payload["query"] = query
    return payload


def markdown_table(rows: Sequence[Mapping[str, Any]], columns: Sequence[tuple[str, str]]) -> str:
    def cell(value: Any) -> str:
        if value is None:
            return "NA"
        if isinstance(value, bool):
            return "true" if value else "false"
        if isinstance(value, float):
            return f"{value:.6g}"
        return str(value).replace("|", "\\|").replace("\n", " ")

    header = "| " + " | ".join(label for _, label in columns) + " |"
    separator = "|" + "|".join("---" for _ in columns) + "|"
    body = [
        "| " + " | ".join(cell(row.get(field)) for field, _ in columns) + " |"
        for row in rows
    ]
    return "\n".join([header, separator, *body])


def biological_answer_for_claims(
    *, g1: Mapping[str, Any], g2: Mapping[str, Any]
) -> str:
    if g2["dataset_numerator"]:
        return (
            f"**目前最強的分子證據結論是 {format_count(g2['dataset_numerator'])} 個位點達到 "
            f"{CLAIM_NAMES['G2']}。** 這表示 focal ALT molecules 內有與多個 LongPhase-S PASS calls "
            "一致的局部分子結構；尚無正交 cellular identity，因此不能升級成 L1 或 L2。"
        )
    if g1["dataset_numerator"]:
        return (
            f"**目前沒有位點達到 {CLAIM_NAMES['G2']}。** 因此本輪最多支持 focal ALT reads 內的"
            "甲基異質性與局部共分離，不支持 cellular identity 或 lineage/order。"
        )
    return (
        f"**目前沒有位點達到 {CLAIM_NAMES['G2']}，G1 也為 0。** 因此本輪最多支持 focal ALT "
        "reads 內的甲基異質性或 residual epigenetic partition；沒有 read-level partner "
        "co-segregation 證據，更不支持 cellular identity 或 lineage/order。"
    )


def build_report_sections(
    claim_rows: Sequence[Mapping[str, Any]],
    stratum_rows: Sequence[Mapping[str, Any]],
    truth_matrix_rows: Sequence[Mapping[str, Any]],
    flow_rows: Sequence[Mapping[str, Any]],
    audit_rows: Sequence[Mapping[str, Any]],
    technical_rows: Sequence[Mapping[str, Any]],
    case_data: Mapping[str, Sequence[Mapping[str, Any]]],
    source_rows: Sequence[Mapping[str, Any]],
    history_rows: Sequence[Mapping[str, Any]],
    m1_operational: Mapping[str, Any],
    m2_evaluability: Mapping[str, Any],
    background_control_gate: Mapping[str, Any],
    source_attestation: Mapping[str, Any],
    *,
    include_inline_tables: bool = True,
) -> list[tuple[str, str]]:
    def inline_table(
        rows: Sequence[Mapping[str, Any]],
        columns: Sequence[tuple[str, str]],
        *,
        heading: str | None = None,
    ) -> str:
        if not include_inline_tables:
            return ""
        prefix = f"\n\n{heading}\n\n" if heading else "\n\n"
        return prefix + markdown_table(rows, columns)

    by_claim = {str(row["claim_id"]): row for row in claim_rows}
    m1 = by_claim["M1"]
    m2 = by_claim["M2"]
    g1 = by_claim["G1"]
    g2 = by_claim["G2"]
    r1 = by_claim["R1"]
    b1 = by_claim["B1"]
    c1 = by_claim["C1"]
    group_limit_examples = m2_evaluability[
        "group_count_exceeds_planning_model_examples"
    ]
    group_limit_count = m2_evaluability[
        "n_group_count_exceeds_planning_model_maximum"
    ]
    aligned_below_power_count = to_int(
        m2_evaluability["independent_logic_audit"].get(
            "n_evaluable_sites_with_aligned_axis_below_negative_evaluability_power",
            0,
        ),
        field="M2 aligned-below-power report count",
    )
    if group_limit_examples:
        group_limit_sites = ", ".join(
            f"{row['dataset']} {row['chrom']}:{row['pos']} {row['ref']}>{row['alt']} "
            f"({row['observed_methyl_groups']} groups)"
            for row in group_limit_examples
        )
    else:
        group_limit_sites = "none observed"
    if source_attestation.get("release_gate_pass") is True:
        source_attestation_statement = (
            "Tumor-REF 長計算另有 bounded retrospective source-file identity receipt PASS："
            "analyzer 與 focal library 在執行中觀察、完成後未漂移，並與 producer manifest hash 綁定。"
            "這不是 prelaunch lock，也不是完整執行環境 attestation。"
        )
    else:
        source_attestation_statement = (
            "此版本是 intermediate terminal build，尚未納入 tumor-REF post-run source identity receipt，"
            "不得作 Task Type B 最終發布版。"
        )
    biological_answer = biological_answer_for_claims(g1=g1, g2=g2)
    if c1["dataset_denominator"]:
        c1_summary = (
            f"C1={c1['dataset_numerator']}/{c1['dataset_denominator']} "
            f"({format_percent(c1['dataset_ratio'])})"
        )
    else:
        c1_summary = (
            f"C1={c1['dataset_numerator']}/0 evaluable（缺少 prespecified focal-partner joint CN/CCF "
            "model，屬結構性 `NOT_EVALUABLE`；不是 CN/CCF 否證 0 個候選）"
        )
    conclusion = (
        "## 證據結論先行\n\n"
        f"{biological_answer}\n\n"
        f"M1 為 **{format_count(m1['dataset_numerator'])}/{format_count(m1['dataset_denominator'])} "
        f"({format_percent(m1['dataset_ratio'])})**；M2 為 "
        f"**{format_count(m2['dataset_numerator'])}/{format_count(m2['dataset_denominator'])} "
        f"({format_percent(m2['dataset_ratio'])})**；G1 為 "
        f"**{format_count(g1['dataset_numerator'])}/{format_count(g1['dataset_denominator'])} "
        f"({format_percent(g1['dataset_ratio'])})**。後續 robustness/background/CN-CCF 通過數依序為 "
        f"R1={r1['dataset_numerator']}、B1={b1['dataset_numerator']}、{c1_summary}。\n\n"
        "M1 的比例是依本流程 frozen gate 計算的 **operational screen yield**；它不是腫瘤中的"
        "生物 prevalence，也不是已另行建立的 null-valid-only census。任何 prevalence 敘述都需要"
        "先補齊並明示 null-valid/invalid 的全母體分層。技術上可評估的位點為 "
        f"**{format_count(m1_operational['n_screen_evaluable'])}**，技術上不可評估為 "
        f"**{format_count(m1_operational['n_screen_not_evaluable'])}**；在技術可評估子集中的次要描述比例為 "
        f"**{format_count(m1_operational['n_flagged_stable_multigroup'])}/"
        f"{format_count(m1_operational['n_screen_evaluable'])} "
        f"({format_percent(m1_operational['flag_yield_among_screen_evaluable'])})**。\n\n"
        "本報告的 469,849 是 **dataset-site union**：同一座標若出現在不同 dataset（包含 "
        "HCC1395 與其 DORADO 技術重處理）會是不同列；它不是 469,849 個 unique biological loci。"
        f"\n\n**發布 provenance：** {source_attestation_statement}"
    )
    claims = (
        "## 九層 claim ladder 保留不可評估與未執行狀態\n\n"
        "M1 以全部 dataset-sites 表示 flagged/not-flagged 的 operational screen，技術可評估性另行保留，"
        "且 not-flagged 不代表生物陰性。M2 及其後各層使用自己的預先定義分母；"
        "`NOT_EVALUABLE` 與 `NOT_RUN` 不得併入 FAIL。"
        "L1/L2 沒有正交證據輸入，必須維持 `NOT_EVALUABLE`/`NOT_RUN`。"
        + inline_table(
            claim_rows,
            (
                ("claim_id", "ID"),
                ("claim_name", "正式名稱"),
                ("status", "狀態"),
                ("dataset_numerator", "分子"),
                ("dataset_denominator", "分母"),
                ("dataset_ratio", "比例"),
                ("dataset_not_evaluable", "NOT_EVALUABLE"),
                ("dataset_not_run", "NOT_RUN"),
            ),
        )
    )
    strata = (
        "## 每個資料集與 focal truth stratum 都保留 M1/M2/G1/G2 的實際分母\n\n"
        "HCC1395_DORADO 是 HCC1395 的技術重處理，不增加 biological n；site-level truth "
        "只指 focal site 且僅在分群後 post hoc 加入。Pair-level focal x partner truth 另表呈現，"
        "不能由 focal truth 推定 partner truth。"
        + inline_table(
            stratum_rows,
            (
                ("stratum_type", "分層"),
                ("stratum", "資料集/truth"),
                ("claim_id", "claim"),
                ("numerator", "分子"),
                ("denominator", "分母"),
                ("ratio", "比例"),
                ("not_evaluable", "NOT_EVALUABLE"),
                ("not_run", "NOT_RUN"),
            ),
        )
        + inline_table(
            truth_matrix_rows,
            (
                ("focal_truth_label", "focal truth"),
                ("partner_truth_label", "partner truth"),
                ("n_all_pair_rows", "all pair rows"),
                ("n_g1_formal_pair_rows", "G1 formal pair rows"),
            ),
            heading="### Pair-level focal x partner truth matrix",
        )
    )
    flow = (
        "## 最新 ClairS、LongPhase-S 與 terminal tag 契約全部通過\n\n"
        "Allele/MM/ML 來自 frozen raw BAM；HP/PS 只由同一次 LongPhase-S producer 的 external sidecar "
        "在 focal ALT 選取前回接。7/7 producer receipt 均證明 BAM-named output 是即時消費的 named FIFO，"
        "`persisted_bam=false`、regular tagged BAM=0；因此沒有產生或覆寫原本 tagged BAM。"
        "HP/PS 是 phase annotation，不是 cellular label。"
        + inline_table(
            flow_rows,
            (
                ("step_order", "順序"),
                ("flow_step", "流程"),
                ("status", "狀態"),
                ("scope", "範圍"),
                ("contract", "契約"),
            ),
        )
        + inline_table(
            audit_rows,
            (
                ("audit", "稽核"),
                ("schema_version", "schema"),
                ("status", "狀態"),
                ("scope", "範圍"),
                ("key_result", "關鍵結果"),
            ),
        )
    )
    technical = (
        "## HCC1395 雙流程只回答技術重現\n\n"
        "同一 biological specimen 的 exact focal+partner pair 可檢查流程穩定性，但不能當成獨立病人重現。"
        + inline_table(
            technical_rows,
            (
                ("comparison", "比較"),
                ("status", "狀態"),
                ("numerator", "分子"),
                ("denominator", "分母"),
                ("ratio", "比例"),
                ("biological_n", "biological n"),
                ("interpretation", "解釋"),
            ),
        )
    )
    views = (
        "## 四個視角分開回答整體、canonical、極端與可解釋個案\n\n"
        "Aggregate 不選 focal pair；其餘角色各自使用預先固定且 deterministic 的選取定義。"
        "Canonical target 缺失時明示 N/A，禁止用 extreme 代替。不同角色可以落在同一 pair，"
        "但會標示重複角色，不能因此把角色定義合併。`evaluable_statistical_negative` 僅表示"
        "已執行的 999-permutation endpoint 未通過，不是 negative biological result；其他 N/A "
        "也不能解讀為生物陰性。"
        + inline_table(
            case_data["case_selection_registry"],
            (
                ("view_role", "視角"),
                ("selection_definition", "獨立選取定義"),
                ("selection_status", "狀態"),
                ("preferred_target", "預定 target"),
                ("selected_focal_site", "focal"),
                ("selected_partner_site", "partner"),
                ("negative_result", "negative result"),
            ),
        )
    )
    case = case_data["case_summary"][0]
    if case["selection_mode"] == "NOT_APPLICABLE_NO_ELIGIBLE_PAIR":
        case_opening = (
            "Pair table 沒有具完整 exact p、Cramer's V 與 delta ALT fraction 的可排序列；"
            "以下保留結構化 **N/A rows**。這是 presentation not-applicable，"
            "**不是 negative biological result**，也不能推論沒有局部分子結構。"
        )
    elif case["selection_mode"] == "NON_CONFIRMING_FALLBACK":
        case_opening = (
            "G2=0；以下是依 fixed-margins exact p、Cramer's V、delta ALT fraction 排序最前的"
            "**未確認個案（non-confirming）**。它只用來顯示實際資料形狀，不得計為候選。"
        )
    else:
        case_opening = (
            f"以下是實際 {CLAIM_NAMES['G2']} 的最強 exact-p/effect witness；"
            "它仍只支持 bulk molecules 的局部結構。"
        )
    cases = (
        "## Actual case 同時呈現 group x partner、four-state 與 joint signature\n\n"
        f"{case_opening}"
        + inline_table(
            case_data["case_summary"],
            (
                ("selection_mode", "選取模式"),
                ("evidence_status", "證據狀態"),
                ("negative_result", "negative result"),
                ("focal_site", "focal"),
                ("partner_site", "partner"),
                ("g2_status", "G2"),
                ("exact_p", "exact p"),
                ("global_by_q", "global BY q"),
                ("cramers_v", "Cramer's V"),
                ("delta_alt_fraction", "delta ALT fraction"),
                ("conditional_status", "conditional"),
                ("callability_gate_status", "callability"),
                ("modal_assignment_ari_min", "ARI min"),
                ("hp_axis_confound", "HP confound"),
                ("technical_axis_confound", "technical confound"),
                ("residual_unexplained_multigroup", "residual"),
            ),
        )
        + inline_table(
            case_data["case_view_summary"],
            (
                ("view_role", "視角"),
                ("selection_status", "選取狀態"),
                ("evidence_status", "證據"),
                ("negative_result", "統計 negative"),
                ("negative_result_scope", "negative 範圍"),
                ("focal_site", "focal"),
                ("partner_site", "partner"),
                ("focal_truth_label", "focal truth"),
                ("partner_truth_label", "partner truth"),
                ("topology_order_status", "posthoc topology"),
                ("normal_focal_status", "normal focal"),
                ("normal_partner_status", "normal partner"),
                ("same_pair_view_roles", "共享 pair 角色"),
            ),
            heading="### 多視角實際個案與結構化 N/A",
        )
        + inline_table(
            case_data["case_group_partner"],
            (
                ("methyl_group", "methyl group"),
                ("partner_call", "partner call"),
                ("count", "reads"),
                ("status", "status"),
            ),
            heading="### Group x partner",
        )
        + inline_table(
            case_data["case_four_state"],
            (("state", "state"), ("count", "reads"), ("relation_status", "fixed-error-model status")),
            heading="### Four-state",
        )
        + inline_table(
            case_data["case_four_state_models"],
            (
                ("model", "模型"),
                ("violations", "violations"),
                ("denominator", "分母"),
                ("violation_rate", "比例"),
                ("p_exact_greater", "exact p"),
                ("upper_per_relation", "98.33% per-model upper"),
                ("error_ceiling", "固定 error ceiling"),
                ("familywise_confidence", "familywise confidence"),
                ("multiplicity_method", "multiplicity"),
                ("status", "狀態"),
                ("split_support", "split support"),
                ("compatible", "compatible"),
            ),
            heading="### Four-state fixed-error model（完整模型）",
        )
        + inline_table(
            case_data["case_joint_signature"],
            (
                ("methyl_group", "methyl group"),
                ("joint_signature", "joint signature"),
                ("count", "reads"),
                ("status", "status"),
            ),
            heading="### Joint signature",
        )
        + "\n\nMatched-normal 只在 focal candidate 評估；partner normal 在本流程"
        "固定標記為 `NOT_EVALUATED_BY_DESIGN`。Topology 僅為 pairwise posthoc context，"
        "不是 cellular lineage 或可識別 mutation order。"
    )
    ancestral_hypothesis = (
        "## 共同 ancestral ALT 模型的分子預測可被檢驗，但不能由單一位點定案\n\n"
        "使用者提出的情境在分子層級會產生可檢驗的模型預測：若 clone 1 與 clone 2 都承載 focal S1 ALT，"
        "而只有較晚狀態承載 partner S2 ALT，則限制在 S1-ALT molecules 後，S2 應同時出現 REF 與 ALT；"
        "若甲基群能輔助分離兩個 latent molecular states，則不同 methyl groups 對 S2 REF/ALT 或"
        "multi-marker signatures 會呈現可重算的富集。G1/G2正是這個模型預測的統計化版本；"
        "通過只表示資料與預測相容，不表示共同 ancestral ALT 已被直接觀察。\n\n"
        "但 S1-ALT 共同存在只說明讀段共享一個 allele state。即使加入 focal-REF four-state compatibility、"
        "至少多個 genetic markers、CN/purity/CCF 與正交 cellular identity時，同樣資料仍可由"
        "branching、cis-ASM、HP/PS、CN/LOH或其他cell-state機制產生；因此單一sSNV甲基多群不能"
        "推出linear ancestry，也不能把methyl group數直接當clone數。"
    )
    methods = (
        "## 判定方法與 denominator contract\n\n"
        "1. M1 母體是 469,849 個 LongPhase-S recalibrated `FILTER=PASS` chr1-22 biallelic "
        "dataset-site rows；這是 dataset-site union，不是去重後的 unique biological loci。\n"
        "2. M2 只在 M1 PASS、觀察 2-10 個 methyl groups，且 HP/technical measured axes 可判定時評估；"
        "超過 10 群保留 M1 PASS，但 M2 為 `NOT_EVALUABLE`且不進入 PASS/FAIL denominator。"
        f"本次超過上限者為 {format_count(group_limit_count)} 個：{group_limit_sites}。"
        "軸若 effect 達門檻且 permutation p<0.05，已構成 observed aligned confound；即使該軸低於"
        "用來支持陰性判定的 80% planning power，仍保守列為可評估的 M2 FAIL。"
        f"本次共有 {format_count(aligned_below_power_count)} 個可評估位點至少含一個這類軸。"
        "只有要判定 not-aligned 的軸才要求 planning power>=80% 與每群至少 5 reads；"
        "效應達門檻但 permutation p 未達 0.05，或未對齊但 power 不足，才標為 `NOT_EVALUABLE`。"
        "Categorical planning levels 固定為 HP-exact/family/strand=7/5/2；assignment-derived observed levels "
        "只用來證明缺失統計量的 axis 確為 constant，不取代 planning ceilings。axis effect 與 "
        "499-permutation add-one p 是 source-locked "
        "screen producer 的輸出；本 terminal gate 會獨立重算分類、樣本數、p-value 網格與 power，"
        "但沒有從原始 reads 再計一次 axis statistic。\n"
        "3. G1 只在至少一個 exact-testable partner 的 M2 sites 評估；先把所有 M2 exact-testable "
        "focal-partner pairs 納入同一全域 BY family，再套 effect、conditional 999-permutation "
        "sensitivity 與 callability。Permutation 可執行性不是 global FDR family membership 的"
        "事後篩選條件。\n"
        f"4. G2 僅稱 **{CLAIM_NAMES['G2']}**，要求 spatially separated markers、同批 complete-read "
        "joint signature，且 joint-signature 的完整可檢定 family 通過全域 BY；raw p<0.05 只保留為 sensitivity。\n"
        "5. R1/B1/C1 各自只在上游可評估集合內計算；B1 固定使用 G1 正式配對中資訊量最高的"
        "單一配對，選取時不查看 four-state 結果；三種 relation models 使用 Bonferroni，"
        "每模型 98.33% 單側 upper bound 維持 familywise 95% confidence。fixed-error ceiling "
        "無法識別時是 NOT_EVALUABLE。tumor-REF 與 matched-normal REF 的背景重現採較寬鬆的 "
        f"`{background_control_gate['contract']}`：要求 coarse K>=2 且 modal-stable，但不要求 "
        "membership ARI>=0.8。在同一 background payload 上，這個 lenient predicate 是加入 ARI "
        "條件後 predicate 的 superset；B1 只在較寬鬆背景 gate 仍不重現時通過。這不比較 ALT 與 REF "
        "的實際 flag 集合；因此不會增加"
        "B1 pass，只可能保守減少候選。它不是"
        "完全對稱的 primary-M1 replay。\n"
        "6. L1/L2 需要本流程未整合的 orthogonal cellular identity 與 identifiable order evidence。"
    )
    limits = (
        "## 限制、穩健性邊界與下一個可否證 gate\n\n"
        "LongPhase-S PASS call 不是正交 somatic truth；bulk read 共分離仍可能受 allele-specific methylation、"
        "CN/LOH、purity、mapping、read geometry 與未量測 cell state 影響。下一步應取得 exact-locus "
        "allele-specific CN/purity/CCF 與 single-cell、colony、spatial 或 multi-region 同一 cellular population 證據，"
        "再評估 L1/L2。M2 axis effect/p 的信任根是已鎖定 SHA-256 的 screen analyzer 與 focal library；"
        "本輪下游未做全位點 raw-read statistic 重算，因此結論必須保留這項 producer-derived provenance 限制。"
        "背景控制另使用 coarse/modal 的 lenient replication gate；沒有 membership ARI 的背景 partition "
        "仍會被算成重現並保守 veto B1，故背景 FAIL 不等同已證明存在與 focal ALT 相同的穩定分群。"
        f" {source_attestation_statement}"
    )
    recommendations = (
        "## 建議的下一步\n\n"
        "1. 以本次 G2/R1/B1 的 exact site/pair key 鎖定後續 matched-normal 與 CN/CCF input，"
        "不可換 pair 拼接 gate。\n"
        "2. 對仍可評估的候選取得 allele-specific CN、purity、multiplicity 與 CCF，並保留 model-unavailable。\n"
        "3. 只有在 single-cell、colony、spatial 或 multi-region 證據連到同一 cellular population 後，"
        "才評估 L1；L2 另需可識別的 order evidence。"
    )
    further_questions = (
        "## 尚待回答的問題\n\n"
        "- 哪些 G2/R1 分子結構在 matched-normal REF-only methyl background 下仍不重現？\n"
        "- allele-specific CN 與 purity/CCF 納入後，B1 候選有多少仍具一致的 cellular-fraction 解釋？\n"
        "- 技術重處理的一致性是否能在獨立 biological samples 重現，而非只在 biological n=1 成立？\n"
        "- 哪種正交 cellular assay 能對 exact focal/partner key 提供最直接且可否證的 L1 證據？"
    )
    literature = (
        "## 外部方法對照\n\n"
        "- [ClairS (Nature Methods, 2026)](https://www.nature.com/articles/s41592-026-03152-4) "
        "界定 ONT somatic calling；本報告仍把 truth label 與 caller PASS 分開。\n"
        "- [LongPhase-S preprint](https://doi.org/10.1101/2025.11.20.689492) 與 "
        "[source repository](https://github.com/CCU-Bioinformatics-Lab/longphase-s) 支持 recalibration/phase "
        "annotation 的工具定位；HP/PS 不被當成 clone truth。\n"
        "- [MethPhaser](https://pmc.ncbi.nlm.nih.gov/articles/PMC11193733/) 顯示長讀可連結 allele 與 methylation；"
        "這支持觀察 read-linked epigenetic structure，但不自動給出 cellular identity。\n"
        "- [Bulk subclonal reconstruction practical guide](https://pmc.ncbi.nlm.nih.gov/articles/PMC7867630/) "
        "與 [bulk clone inference pitfalls](https://pmc.ncbi.nlm.nih.gov/articles/PMC7044161/) 說明 purity、CN、"
        "取樣與模型不可識別性；因此本報告將 G2 與 L1/L2 分開。\n"
        "- [Single-nucleus methylome/RNA/exome integration](https://www.nature.com/articles/s41588-026-02642-7) "
        "代表可用於 L1 的正交 cellular evidence 類型，而非本次 bulk read 分群已完成的證據。"
    )
    sections: list[tuple[str, str]] = [
        ("title", f"# {TITLE}"),
        ("conclusion", conclusion),
        ("claims", claims),
        ("strata", strata),
        ("flow", flow),
        ("technical", technical),
        ("views", views),
        ("case", cases),
        ("ancestral_hypothesis", ancestral_hypothesis),
    ]
    if history_rows:
        sections.append(
            (
                "history",
                "## Earlier FP-only 結果只回答 specificity\n\n"
                "舊分析的研究問題是 FP specificity，不是全 sSNV prevalence；兩者的母體與分母不同，"
                "不得把 FP-only rate 當成本次 469,849-site prevalence，也不得直接作數值升降比較。"
                + inline_table(
                    history_rows,
                    (
                        ("scope", "歷史範圍"),
                        ("question", "研究問題"),
                        ("path", "來源"),
                        ("sha256", "SHA-256"),
                    ),
                ),
            )
        )
    sections.extend(
        [
            ("methods", methods),
            ("literature", literature),
            ("limits", limits),
            ("recommendations", recommendations),
            ("further_questions", further_questions),
            (
                "sources",
                "## 證據來源與 hash"
                + inline_table(
                    source_rows,
                    (
                        ("role", "角色"),
                        ("display_path", "路徑"),
                        ("size_bytes", "bytes"),
                        ("sha256", "SHA-256"),
                    ),
                ),
            ),
        ]
    )
    return sections


def markdown_block(block_id: str, body: str, source_id: str | None = None) -> dict[str, Any]:
    block: dict[str, Any] = {
        "id": block_id,
        "type": "markdown",
        "body": body,
        "layout": "full",
    }
    if source_id:
        block["sourceId"] = source_id
    return block


def chart_block(block_id: str, chart_id: str) -> dict[str, Any]:
    return {"id": block_id, "type": "chart", "chartId": chart_id, "layout": "full"}


def table_block(block_id: str, table_id: str) -> dict[str, Any]:
    return {"id": block_id, "type": "table", "tableId": table_id, "layout": "full"}


def build_portable_summary(
    datasets: Mapping[str, Sequence[Mapping[str, Any]]],
) -> str:
    by_claim = {
        str(row["claim_id"]): row for row in datasets["claim_ladder"]
    }
    claim_counts = "、".join(
        f"{claim_id} {format_count(by_claim[claim_id]['dataset_numerator'])}/"
        f"{format_count(by_claim[claim_id]['dataset_denominator'])}"
        for claim_id in ("M1", "M2", "G1", "G2")
    )
    case_summary = datasets["case_summary"][0]
    if case_summary["selection_mode"] == "NOT_APPLICABLE_NO_ELIGIBLE_PAIR":
        case_clause = "actual case=N/A（無 eligible pair）"
    else:
        case_label = (
            "G2 候選"
            if int(by_claim["G2"]["dataset_numerator"]) > 0
            else "non-confirming witness"
        )
        case_clause = f"右={case_label} 的 methyl-group partner R/A composition"
    return (
        f"**結論：**469,849 LongPhase-S PASS dataset-sites：{claim_counts}。"
        "G2 僅為 molecular-haplotype base candidate，非 cellular subclone/lineage；"
        "共同 ancestral ALT 假設可預測 focal-ALT reads 內 partner R/A 混合與methyl-group富集，"
        "但單一focal位點本身不能區分linear、branching或cellular clone；"
        "tree、terminal HP/PS、output/immutability audits 均 PASS。"
        f"圖左=claim rates；{case_clause}。"
        "完整分母、限制與 SHA-256 見 `report.md`/Sources。"
    )


def build_artifact(
    sections: Sequence[tuple[str, str]],
    datasets: Mapping[str, Sequence[Mapping[str, Any]]],
    sources: Sequence[Mapping[str, Any]],
    generated_at: str,
    include_history: bool,
) -> dict[str, Any]:
    charts = [
        {
            "id": "chart-overview-case",
            "title": "Claim funnel and actual case molecular composition",
            "subtitle": (
                "M1=FLAGGED/all dataset-sites（另列非 prevalence 的描述性 biological-site aggregation）；"
                "M2-G2=PASS/claim-specific evaluable denominator；case bars 為各 methyl group 內 partner R/A read 比例"
            ),
            "type": "bar",
            "dataset": "overview_case_chart",
            "sourceId": "src-overview-case-view",
            "intent": "comparison",
            "encodings": {
                "x": {
                    "field": "category",
                    "type": "ordinal",
                    "label": "claim or case methyl group",
                },
                "y": {
                    "field": "proportion",
                    "type": "quantitative",
                    "label": "proportion",
                    "format": "percent",
                },
                "color": {
                    "field": "series",
                    "type": "nominal",
                    "label": "evidence series",
                },
                "tooltip": [
                    {"field": "panel", "type": "nominal", "label": "panel"},
                    {"field": "numerator", "type": "quantitative", "label": "numerator"},
                    {"field": "denominator", "type": "quantitative", "label": "denominator"},
                    {"field": "status", "type": "nominal", "label": "status"},
                    {
                        "field": "metric_definition",
                        "type": "nominal",
                        "label": "metric definition",
                    },
                    {
                        "field": "interpretation",
                        "type": "nominal",
                        "label": "interpretation",
                    },
                ],
            },
            "settings": {
                "orientation": "vertical",
                "groupMode": "grouped",
                "sort": "none",
            },
            "legend": {"position": "bottom", "title": "evidence series"},
            "valueFormat": "percent",
            "layout": "full",
        },
        {
            "id": "chart-claim-ladder",
            "title": "九層 claim ladder 通過比例",
            "subtitle": (
                "M1 為全 dataset-sites operational FLAGGED yield（非 prevalence）；"
                "M2-L2 才是 PASS/claim-specific evaluable denominator"
            ),
            "type": "bar",
            "dataset": "claim_chart",
            "sourceId": "src-claims-view",
            "intent": "comparison",
            "encodings": {
                "x": {"field": "claim_id", "type": "ordinal", "label": "claim"},
                "y": {"field": "ratio", "type": "quantitative", "label": "通過比例", "format": "percent"},
                "color": {"field": "series_label", "type": "nominal", "label": "指標定義"},
                "tooltip": [
                    {"field": "numerator", "type": "quantitative", "label": "分子"},
                    {"field": "denominator", "type": "quantitative", "label": "分母"},
                    {"field": "status", "type": "nominal", "label": "狀態"},
                    {"field": "metric_definition", "type": "nominal", "label": "指標定義"},
                    {"field": "interpretation", "type": "nominal", "label": "解釋界線"},
                ],
            },
            "settings": {"orientation": "vertical", "groupMode": "grouped", "sort": "none"},
            "legend": {"position": "bottom", "title": "指標定義"},
            "valueFormat": "percent",
            "layout": "full",
        },
        {
            "id": "chart-strata",
            "title": "資料集與 focal truth strata 的 M1/M2/G1/G2 比例",
            "subtitle": (
                "M1=FLAGGED/all stratum dataset-sites（非 prevalence）；M2-G2=PASS/evaluable；"
                "site truth 僅指 focal 且為 post hoc"
            ),
            "type": "bar",
            "dataset": "stratum_claim_metrics",
            "sourceId": "src-strata-view",
            "intent": "comparison",
            "encodings": {
                "x": {"field": "stratum", "type": "nominal", "label": "資料集或 truth"},
                "y": {"field": "ratio", "type": "quantitative", "label": "通過比例", "format": "percent"},
                "color": {"field": "claim_id", "type": "nominal", "label": "claim"},
                "tooltip": [
                    {"field": "stratum_type", "type": "nominal", "label": "分層"},
                    {"field": "numerator", "type": "quantitative", "label": "分子"},
                    {"field": "denominator", "type": "quantitative", "label": "分母"},
                    {"field": "not_evaluable", "type": "quantitative", "label": "NOT_EVALUABLE"},
                    {"field": "not_run", "type": "quantitative", "label": "NOT_RUN"},
                ],
            },
            "settings": {"orientation": "vertical", "groupMode": "grouped", "sort": "none"},
            "legend": {"position": "bottom", "title": "claim"},
            "valueFormat": "percent",
            "layout": "full",
        },
        {
            "id": "chart-focal-partner-truth",
            "title": "Pair-level focal x partner truth composition",
            "subtitle": "Focal truth 不代表 partner truth；所有標籤僅作 post hoc annotation",
            "type": "bar",
            "dataset": "focal_partner_truth_matrix",
            "sourceId": "src-truth-matrix-view",
            "intent": "composition",
            "encodings": {
                "x": {
                    "field": "focal_truth_label",
                    "type": "ordinal",
                    "label": "focal truth",
                },
                "y": {
                    "field": "n_all_pair_rows",
                    "type": "quantitative",
                    "label": "pair rows",
                    "format": "number",
                },
                "color": {
                    "field": "partner_truth_label",
                    "type": "nominal",
                    "label": "partner truth",
                },
                "tooltip": [
                    {
                        "field": "n_g1_formal_pair_rows",
                        "type": "quantitative",
                        "label": "G1 formal pair rows",
                    }
                ],
            },
            "settings": {"orientation": "vertical", "groupMode": "stacked", "sort": "none"},
            "legend": {"position": "bottom", "title": "partner truth"},
            "layout": "full",
        },
        {
            "id": "chart-case-group-partner",
            "title": "Report case methyl group x partner read counts",
            "subtitle": "G2=0 時為預定 non-confirming witness；read counts 不代表 cellular clone",
            "type": "bar",
            "dataset": "case_group_partner",
            "sourceId": "src-case-group-view",
            "intent": "composition",
            "encodings": {
                "x": {"field": "methyl_group", "type": "nominal", "label": "methyl group"},
                "y": {"field": "count", "type": "quantitative", "label": "reads", "format": "number"},
                "color": {"field": "partner_call", "type": "nominal", "label": "partner call"},
                "tooltip": [{"field": "status", "type": "nominal", "label": "conditional status"}],
            },
            "settings": {"orientation": "vertical", "groupMode": "grouped", "sort": "none"},
            "legend": {"position": "bottom", "title": "partner call"},
            "layout": "full",
        },
    ]
    if not any(
        isinstance(row.get("count"), (int, float))
        and not isinstance(row.get("count"), bool)
        for row in datasets["case_group_partner"]
    ):
        charts = [chart for chart in charts if chart["id"] != "chart-case-group-partner"]
    tables = [
        {
            "id": "table-claim-ladder",
            "title": "九層 claim 狀態與精確分母",
            "subtitle": "M1=FLAGGED/NOT_FLAGGED operational；M2-L2 保留 claim-specific states",
            "dataset": "claim_ladder",
            "sourceId": "src-claims-view",
            "defaultSort": {"field": "claim_order", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "claim_order", "label": "order", "format": "number"},
                {"field": "claim_id", "label": "claim", "type": "text"},
                {"field": "claim_name", "label": "正式名稱", "type": "text"},
                {"field": "status", "label": "狀態", "type": "text"},
                {"field": "dataset_numerator", "label": "分子", "format": "number"},
                {"field": "dataset_denominator", "label": "分母", "format": "number"},
                {"field": "dataset_ratio", "label": "比例", "format": "percent"},
            ],
        },
        {
            "id": "table-strata",
            "title": "Per-sample 與 focal-truth-stratum 分子分母",
            "subtitle": "M1/M2/G1/G2；7 datasets 與 TP/FP/UNASSESSED",
            "dataset": "stratum_claim_metrics",
            "sourceId": "src-strata-view",
            "defaultSort": {"field": "stratum_type", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "stratum_type", "label": "分層", "type": "text"},
                {"field": "stratum", "label": "資料集/truth", "type": "text"},
                {"field": "claim_id", "label": "claim", "type": "text"},
                {"field": "numerator", "label": "分子", "format": "number"},
                {"field": "denominator", "label": "分母", "format": "number"},
                {"field": "ratio", "label": "比例", "format": "percent"},
                {"field": "not_evaluable", "label": "NOT_EVALUABLE", "format": "number"},
                {"field": "not_run", "label": "NOT_RUN", "format": "number"},
            ],
        },
        {
            "id": "table-m2-evaluability",
            "title": "M2 2-10群 planning-model 邊界",
            "subtitle": "超過10群保留M1 PASS；M2為NOT_EVALUABLE，G1/G2/B1不執行",
            "dataset": "m2_evaluability",
            "sourceId": "src-claims-view",
            "defaultSort": {"field": "record_type", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "record_type", "label": "record", "type": "text"},
                {"field": "site", "label": "site", "type": "text"},
                {"field": "count", "label": "count", "format": "number"},
                {
                    "field": "observed_methyl_groups",
                    "label": "observed groups",
                    "format": "number",
                },
                {
                    "field": "maximum_supported_methyl_groups",
                    "label": "model max",
                    "format": "number",
                },
                {
                    "field": "claim_behavior",
                    "label": "claim behavior",
                    "type": "text",
                },
                {"field": "reason", "label": "reason", "type": "text"},
            ],
        },
        {
            "id": "table-focal-partner-truth",
            "title": "Pair-level focal x partner truth matrix",
            "subtitle": "Site truth axis 與 partner truth axis 分離；truth 不參與候選選取",
            "dataset": "focal_partner_truth_matrix",
            "sourceId": "src-truth-matrix-view",
            "defaultSort": {"field": "focal_truth_label", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "focal_truth_label", "label": "focal truth", "type": "text"},
                {"field": "partner_truth_label", "label": "partner truth", "type": "text"},
                {"field": "n_all_pair_rows", "label": "all pair rows", "format": "number"},
                {
                    "field": "n_g1_formal_pair_rows",
                    "label": "G1 formal pair rows",
                    "format": "number",
                },
            ],
        },
        {
            "id": "table-flow",
            "title": "Current ClairS/LongPhase-S processing contract",
            "subtitle": "7 datasets；same-run PASS backbone 與 external terminal HP/PS join",
            "dataset": "flow",
            "sourceId": "src-flow-view",
            "defaultSort": {"field": "step_order", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "step_order", "label": "order", "format": "number"},
                {"field": "flow_step", "label": "流程", "type": "text"},
                {"field": "status", "label": "狀態", "type": "text"},
                {"field": "scope", "label": "範圍", "type": "text"},
                {"field": "contract", "label": "契約", "type": "text"},
            ],
        },
        {
            "id": "table-audits",
            "title": "Release-blocking audit results",
            "subtitle": (
                "screen/tag、output、post immutability、tree input、reference identity 均須 PASS"
            ),
            "dataset": "audit_status",
            "sourceId": "src-audit-view",
            "defaultSort": {"field": "audit_order", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "audit_order", "label": "order", "format": "number"},
                {"field": "audit", "label": "稽核", "type": "text"},
                {"field": "schema_version", "label": "schema", "type": "text"},
                {"field": "status", "label": "狀態", "type": "text"},
                {"field": "scope", "label": "範圍", "type": "text"},
                {"field": "key_result", "label": "關鍵結果", "type": "text"},
            ],
        },
        {
            "id": "table-technical",
            "title": "HCC1395 technical replication",
            "subtitle": "Exact shared focal+partner opportunities；biological n=1",
            "dataset": "technical_replication",
            "sourceId": "src-technical-view",
            "defaultSort": {"field": "comparison", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "comparison", "label": "比較", "type": "text"},
                {"field": "status", "label": "狀態", "type": "text"},
                {"field": "numerator", "label": "分子", "format": "number"},
                {"field": "denominator", "label": "分母", "format": "number"},
                {"field": "ratio", "label": "比例", "format": "percent"},
                {"field": "not_evaluable_one_platform_only", "label": "單平台 only", "format": "number"},
                {"field": "biological_n", "label": "biological n", "format": "number"},
            ],
        },
        {
            "id": "table-case-selection-registry",
            "title": "Four-view case selection registry",
            "subtitle": "Aggregate、pre-registered oracle、global extreme、well-explained 各自獨立選取",
            "dataset": "case_selection_registry",
            "sourceId": "src-case-selection-view",
            "defaultSort": {"field": "view_role", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "view_order", "label": "order", "format": "number"},
                {"field": "view_role", "label": "視角", "type": "text"},
                {"field": "selected_focal_site", "label": "focal", "type": "text"},
                {"field": "selected_partner_site", "label": "partner", "type": "text"},
                {"field": "negative_result", "label": "negative result", "type": "boolean"},
            ],
        },
        {
            "id": "table-case-summary",
            "title": "Actual case exact-p/effect evidence",
            "subtitle": "G2=0 時明示 NON_CONFIRMING_FALLBACK",
            "dataset": "case_summary",
            "sourceId": "src-case-summary-view",
            "defaultSort": {"field": "case_rank", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "case_rank", "label": "rank", "format": "number"},
                {"field": "focal_site", "label": "focal", "type": "text"},
                {"field": "partner_site", "label": "partner", "type": "text"},
                {"field": "g2_status", "label": "G2", "type": "text"},
                {"field": "exact_p", "label": "exact p", "format": "number"},
                {"field": "global_by_q", "label": "BY q", "format": "number"},
                {"field": "cramers_v", "label": "Cramer's V", "format": "number"},
                {"field": "delta_alt_fraction", "label": "delta ALT", "format": "number"},
            ],
        },
        {
            "id": "table-case-view-summary",
            "title": "多視角實際個案、統計 negative 與結構化 N/A",
            "subtitle": (
                "Canonical 不以 extreme 替代；negative 僅指已評估 endpoint FAIL，非生物陰性"
            ),
            "dataset": "case_view_summary",
            "sourceId": "src-case-views",
            "defaultSort": {"field": "view_role", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "view_role", "label": "視角", "type": "text"},
                {"field": "selection_status", "label": "選取狀態", "type": "text"},
                {"field": "negative_result", "label": "統計 negative", "type": "boolean"},
                {"field": "focal_site", "label": "focal", "type": "text"},
                {"field": "partner_site", "label": "partner", "type": "text"},
                {"field": "truth_pair", "label": "focal / partner truth", "type": "text"},
                {"field": "topology_order_status", "label": "posthoc topology", "type": "text"},
                {"field": "normal_focal_status", "label": "normal focal", "type": "text"},
            ],
        },
        {
            "id": "table-case-group-partner",
            "title": "Actual case group x partner counts",
            "subtitle": "Methyl group by partner REF/ALT read counts",
            "dataset": "case_group_partner",
            "sourceId": "src-case-group-view",
            "defaultSort": {"field": "row_order", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "row_order", "label": "order", "format": "number"},
                {"field": "methyl_group", "label": "methyl group", "type": "text"},
                {"field": "partner_call", "label": "partner call", "type": "text"},
                {"field": "count", "label": "reads", "format": "number"},
                {"field": "status", "label": "status", "type": "text"},
            ],
        },
        {
            "id": "table-case-four-state",
            "title": "Actual case four-state counts",
            "subtitle": "RR/AR/RA/AA/O/X under the fixed error model",
            "dataset": "case_four_state",
            "sourceId": "src-case-four-state-view",
            "defaultSort": {"field": "state_order", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "state_order", "label": "order", "format": "number"},
                {"field": "state", "label": "state", "type": "text"},
                {"field": "count", "label": "reads", "format": "number"},
                {"field": "relation_status", "label": "relation", "type": "text"},
            ],
        },
        {
            "id": "table-case-four-state-models",
            "title": "Primary case four-state fixed-error models",
            "subtitle": (
                "三個 relation models 使用 Bonferroni；每模型 98.33% upper 維持整體 95% confidence"
            ),
            "dataset": "case_four_state_models",
            "sourceId": "src-case-four-state-view",
            "defaultSort": {"field": "model", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "model", "label": "model", "type": "text"},
                {"field": "violations", "label": "violations", "format": "number"},
                {"field": "denominator", "label": "denominator", "format": "number"},
                {
                    "field": "upper_per_relation",
                    "label": "98.33% per-model upper",
                    "format": "percent",
                },
                {"field": "error_ceiling", "label": "ceiling", "format": "percent"},
                {
                    "field": "familywise_confidence",
                    "label": "familywise confidence",
                    "format": "percent",
                },
                {"field": "multiplicity_method", "label": "multiplicity", "type": "text"},
                {"field": "status", "label": "status", "type": "text"},
            ],
        },
        {
            "id": "table-case-view-four-state-models",
            "title": "Selected views four-state model comparison",
            "subtitle": "同一 pair 若被多個視角選中會保留角色標記；不視為獨立個案",
            "dataset": "case_view_four_state_models",
            "sourceId": "src-case-view-four-state",
            "defaultSort": {"field": "view_role", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "view_role", "label": "視角", "type": "text"},
                {"field": "model", "label": "model", "type": "text"},
                {"field": "violations", "label": "violations", "format": "number"},
                {"field": "denominator", "label": "denominator", "format": "number"},
                {
                    "field": "upper_per_relation",
                    "label": "98.33% per-model upper",
                    "format": "percent",
                },
                {
                    "field": "familywise_confidence",
                    "label": "familywise confidence",
                    "format": "percent",
                },
                {"field": "status", "label": "status", "type": "text"},
                {"field": "compatible", "label": "compatible", "type": "boolean"},
            ],
        },
        {
            "id": "table-case-joint-signature",
            "title": "Actual case joint signature",
            "subtitle": "Complete-read multi-marker signature by methyl group",
            "dataset": "case_joint_signature",
            "sourceId": "src-case-joint-view",
            "defaultSort": {"field": "row_order", "direction": "asc"},
            "density": "spacious",
            "layout": "full",
            "columns": [
                {"field": "row_order", "label": "order", "format": "number"},
                {"field": "methyl_group", "label": "methyl group", "type": "text"},
                {"field": "joint_signature", "label": "signature", "type": "text"},
                {"field": "count", "label": "reads", "format": "number"},
                {"field": "status", "label": "status", "type": "text"},
            ],
        },
        {
            "id": "table-sources",
            "title": "Input source inventory",
            "subtitle": "完整 portable path 與 SHA-256 保存在來源詳情及 Markdown 證據表",
            "dataset": "source_inventory",
            "sourceId": "src-source-view",
            "defaultSort": {"field": "role", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "role", "label": "角色", "type": "text"},
                {"field": "size_bytes", "label": "bytes", "format": "number"},
            ],
        },
    ]
    if include_history:
        tables.append(
            {
                "id": "table-history",
                "title": "Earlier FP-only scope boundary",
                "subtitle": "Specificity question；not an all-sSNV prevalence denominator",
                "dataset": "history",
                "sourceId": "src-history-view",
                "defaultSort": {"field": "scope", "direction": "asc"},
                "density": "spacious",
                "layout": "full",
                "columns": [
                    {"field": "scope", "label": "範圍", "type": "text"},
                    {"field": "question", "label": "研究問題", "type": "text"},
                ],
            }
        )
    section_map = dict(sections)
    required_sections = {
        "title",
        "conclusion",
        "claims",
        "strata",
        "flow",
        "technical",
        "views",
        "case",
        "ancestral_hypothesis",
        "methods",
        "literature",
        "limits",
        "recommendations",
        "further_questions",
        "sources",
    }
    missing_sections = sorted(required_sections.difference(section_map))
    if missing_sections:
        raise ReportContractError(
            f"Complete Markdown report sections are missing: {missing_sections}"
        )
    blocks: list[dict[str, Any]] = [
        markdown_block("block-title", f"# {PORTABLE_TITLE}"),
        markdown_block("block-portable-summary", build_portable_summary(datasets)),
        markdown_block("block-conclusion", section_map["conclusion"], "src-claims-view"),
        chart_block("block-overview-case-chart", "chart-overview-case"),
        markdown_block("block-claims", section_map["claims"], "src-claims-view"),
        chart_block("block-claim-chart", "chart-claim-ladder"),
        table_block("block-claim-table", "table-claim-ladder"),
        table_block("block-m2-evaluability", "table-m2-evaluability"),
        markdown_block("block-strata", section_map["strata"], "src-strata-view"),
        chart_block("block-strata-chart", "chart-strata"),
        table_block("block-strata-table", "table-strata"),
        chart_block("block-focal-partner-truth-chart", "chart-focal-partner-truth"),
        table_block("block-focal-partner-truth-table", "table-focal-partner-truth"),
        markdown_block("block-flow", section_map["flow"], "src-flow-view"),
        table_block("block-flow-table", "table-flow"),
        table_block("block-audit-table", "table-audits"),
        markdown_block("block-technical", section_map["technical"], "src-technical-view"),
        table_block("block-technical-table", "table-technical"),
        markdown_block("block-views", section_map["views"], "src-case-selection-view"),
        table_block("block-case-selection", "table-case-selection-registry"),
        markdown_block("block-case", section_map["case"], "src-case-summary-view"),
        table_block("block-case-summary", "table-case-summary"),
        table_block("block-case-view-summary", "table-case-view-summary"),
    ]
    if any(chart["id"] == "chart-case-group-partner" for chart in charts):
        blocks.append(chart_block("block-case-group-chart", "chart-case-group-partner"))
    blocks.extend(
        [
            table_block("block-case-group", "table-case-group-partner"),
            table_block("block-case-four-state", "table-case-four-state"),
            table_block("block-case-four-state-models", "table-case-four-state-models"),
            table_block(
                "block-case-view-four-state-models",
                "table-case-view-four-state-models",
            ),
            table_block("block-case-joint", "table-case-joint-signature"),
            markdown_block(
                "block-ancestral-hypothesis",
                section_map["ancestral_hypothesis"],
                "src-case-summary-view",
            ),
        ]
    )
    if include_history:
        blocks.extend(
            [
                markdown_block("block-history", section_map["history"], "src-history-view"),
                table_block("block-history-table", "table-history"),
            ]
        )
    blocks.extend(
        [
            markdown_block("block-methods", section_map["methods"], "src-claims-view"),
            markdown_block("block-literature", section_map["literature"]),
            markdown_block("block-limits", section_map["limits"], "src-claims-view"),
            markdown_block(
                "block-recommendations", section_map["recommendations"], "src-claims-view"
            ),
            markdown_block(
                "block-further-questions", section_map["further_questions"]
            ),
            markdown_block("block-sources", section_map["sources"], "src-source-view"),
            table_block("block-source-table", "table-sources"),
        ]
    )
    return {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": PORTABLE_TITLE,
            "description": (
                "Portable visual summary of the schema 2.0.0 all-sSNV claim ladder; "
                "complete evidence remains in the bounded snapshot, sources, and report.md."
            ),
            "generatedAt": generated_at,
            "sources": list(sources),
            "cards": [],
            "charts": charts,
            "tables": tables,
            "blocks": blocks,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "ready",
            "datasets": {name: list(rows) for name, rows in datasets.items()},
        },
        "sources": list(sources),
    }


def validate_artifact_shape(payload: Mapping[str, Any]) -> None:
    if payload.get("surface") != "report":
        raise ReportContractError("Artifact surface must be report")
    manifest = payload.get("manifest")
    snapshot = payload.get("snapshot")
    if not isinstance(manifest, Mapping) or not isinstance(snapshot, Mapping):
        raise ReportContractError("Artifact must contain manifest and snapshot objects")
    blocks = manifest.get("blocks")
    if not isinstance(blocks, list) or not blocks:
        raise ReportContractError("Artifact report has no blocks")
    first = blocks[0]
    if not isinstance(first, Mapping) or first.get("type") != "markdown":
        raise ReportContractError("First report block must be markdown")
    if first.get("body") != f"# {manifest.get('title')}":
        raise ReportContractError("First markdown # title must equal manifest title")
    charts = manifest.get("charts")
    chart_blocks = [block for block in blocks if block.get("type") == "chart"]
    if not isinstance(charts, list) or not charts or not chart_blocks:
        raise ReportContractError("Artifact must contain at least one native chart")
    datasets = snapshot.get("datasets")
    if not isinstance(datasets, Mapping) or any(not isinstance(rows, list) for rows in datasets.values()):
        raise ReportContractError("snapshot.datasets must be an object of row arrays")
    row_count = sum(len(rows) for rows in datasets.values())
    if row_count > MAX_SNAPSHOT_ROWS:
        raise ReportContractError(
            f"Artifact bounded snapshot exceeds {MAX_SNAPSHOT_ROWS} rows: {row_count}"
        )
    tables = manifest.get("tables")
    if not isinstance(tables, list) or not tables:
        raise ReportContractError("Artifact must contain native tables")
    sources = manifest.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ReportContractError("Artifact must contain canonical sources")
    source_ids = {source.get("id") for source in sources if isinstance(source, Mapping)}
    referenced_source_ids = {
        payload["sourceId"]
        for payload in (*blocks, *charts, *tables)
        if isinstance(payload, Mapping) and payload.get("sourceId")
    }
    unresolved = sorted(referenced_source_ids.difference(source_ids))
    if unresolved:
        raise ReportContractError(f"Artifact contains unresolved source IDs: {unresolved}")
    for table in tables:
        default_sort = table.get("defaultSort")
        columns = {column.get("field") for column in table.get("columns", [])}
        if (
            not isinstance(default_sort, Mapping)
            or default_sort.get("direction") not in {"asc", "desc"}
            or default_sort.get("field") not in columns
        ):
            raise ReportContractError(f"Table lacks a valid defaultSort: {table.get('id')}")
    registry = datasets.get("case_selection_registry")
    expected_roles = {
        "aggregate",
        "canonical_pre_registered_oracle",
        "extreme_global_exact_effect",
        "well_explained",
        "evaluable_statistical_negative",
    }
    if (
        not isinstance(registry, list)
        or {row.get("view_role") for row in registry if isinstance(row, Mapping)}
        != expected_roles
        or any(not isinstance(row, Mapping) for row in registry)
    ):
        raise ReportContractError("Artifact case selection registry contract failed")
    for row in registry:
        role = row["view_role"]
        expected_negative = (
            role == "evaluable_statistical_negative"
            and row["selection_status"] == "SELECTED_EVALUATED_ENDPOINT_FAIL"
        )
        if row.get("negative_result") is not expected_negative:
            raise ReportContractError(
                "Artifact statistical-negative role does not reconcile"
            )
    canonical = next(
        row for row in registry if row["view_role"] == "canonical_pre_registered_oracle"
    )
    if canonical["selection_status"] == "TARGET_UNAVAILABLE_NO_SUBSTITUTION" and (
        canonical["selected_focal_site"] != "N/A"
        or canonical["selected_partner_site"] != "N/A"
    ):
        raise ReportContractError("Unavailable canonical target was silently substituted")
    required_block_ids = {
        "block-conclusion",
        "block-claim-table",
        "block-strata-table",
        "block-audit-table",
        "block-case-selection",
        "block-case-view-summary",
        "block-case-four-state-models",
        "block-case-joint",
        "block-ancestral-hypothesis",
        "block-methods",
        "block-limits",
        "block-sources",
        "block-source-table",
    }
    observed_block_ids = {
        block.get("id") for block in blocks if isinstance(block, Mapping)
    }
    if not required_block_ids.issubset(observed_block_ids):
        raise ReportContractError("Portable artifact omits required visible report sections")


def build_outputs(
    inputs: ReportInputs,
    output_dir: Path,
    *,
    repo_root: Path,
    generated_at: str | None = None,
) -> dict[str, Path]:
    output_dir = output_dir.expanduser().resolve()
    if os.path.lexists(output_dir):
        raise ReportContractError(f"Output directory already exists; overwrite prohibited: {output_dir}")
    generated_at = generated_at or now_utc()
    final_release_validation = validate_detached_final_release(inputs)
    final = load_json(inputs.final_dataset, "final report dataset")
    final_receipt = load_json(inputs.final_receipt, "final dataset receipt")
    manifest = load_json(inputs.manifest, "input manifest")
    screen_summary = load_json(inputs.screen_summary, "screen summary")
    cooccurrence_summary = load_json(inputs.cooccurrence_summary, "cooccurrence summary")
    reconciliation = load_json(inputs.output_reconciliation, "output reconciliation")
    immutability = load_json(inputs.post_immutability_audit, "post-run immutability audit")
    tree_audit = load_json(inputs.tree_input_audit, "tree-input audit")
    reference_audit = load_json(
        inputs.reference_identity_audit, "extraction reference identity audit"
    )
    validate_final_scope(final)
    final_builder_path = validate_final_receipt(final_receipt, inputs)
    validate_final_machine_tables(final, inputs)
    final_input_paths = validate_all_final_input_artifacts(final, final_receipt)
    screen_receipt = load_json(
        final_input_paths["screen_receipt"], "terminal screen receipt"
    )
    if (
        screen_receipt.get("schema_name")
        != "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest"
        or screen_receipt.get("status") != "EXECUTION_PASS"
        or screen_receipt.get("pass") is not True
    ):
        raise ReportContractError("Terminal screen receipt is not an execution PASS")
    recovery_input_paths = validate_recovery_provenance_artifacts(screen_receipt)
    source_attestation = validate_tumor_ref_source_identity_attestation(final)
    validate_claim_contract(
        inputs.claim_contract,
        release_gate_pass=source_attestation.get("release_gate_pass") is True,
    )
    validate_manifest(manifest)
    layered, layered_path, producer_bam_audit = load_and_validate_layered_manifest(manifest)
    tag_audit = validate_screen_summary(screen_summary)
    validate_cooccurrence_summary(cooccurrence_summary)
    validate_output_reconciliation(reconciliation)
    validate_post_immutability(
        immutability, inputs.post_immutability_audit, screen_receipt
    )
    validate_tree_audit(tree_audit)
    validate_reference_identity_audit(reference_audit)
    validate_final_input_links(final, inputs)
    metrics = validate_metrics(final, manifest)
    m1_operational = validate_m1_operational_screen(final, metrics)
    m2_evaluability = validate_m2_evaluability_contract(final, metrics)
    if (
        source_attestation.get("release_gate_pass") is True
        and m2_evaluability["independent_logic_audit"].get("status")
        != "PASS_LOGIC_INDEPENDENT_RECOUNT"
    ):
        raise ReportContractError(
            "Source-attested release requires a passing logic-independent M2 recount"
        )
    background_control_gate = validate_background_control_replication_gate(final)
    claim_rows = validate_claim_ladder(final, metrics)
    m2_provenance = validate_m2_axis_statistic_provenance(final)
    candidates = validate_candidates(final, metrics)
    technical_rows = [validate_technical_replication(final)]
    truth_matrix_rows = validate_focal_partner_truth_matrix(final)
    pair, case_mode, selection_registry, selected_pairs = choose_case_pair(
        inputs.cooccurrence_pairs, candidates
    )
    selected_sites = {
        role: (
            find_site_row(inputs.cooccurrence_sites, site_key(selected, pair=True))
            if selected is not None
            else None
        )
        for role, selected in selected_pairs.items()
    }
    selected_site = selected_sites["primary_report_case"]
    case_data = build_case_evidence(
        pair,
        selected_site,
        case_mode,
        selection_registry,
        selected_pairs,
        selected_sites,
        candidates,
    )
    stratum_rows = build_stratum_rows(metrics)
    claim_chart_rows = build_claim_chart_rows(claim_rows)
    overview_case_chart_rows = build_overview_case_chart_rows(
        claim_chart_rows,
        case_data["case_group_partner"],
    )
    flow_rows = build_flow_rows(layered, tag_audit, producer_bam_audit)
    audit_rows = build_audit_rows(
        screen_summary,
        reconciliation,
        immutability,
        tree_audit,
        reference_audit,
        producer_bam_audit,
        m2_provenance,
        source_attestation,
    )

    repo_root = repo_root.expanduser().resolve()
    source_paths: list[tuple[str, Path]] = [
        ("final schema 2.0.0 dataset", inputs.final_dataset),
        ("final dataset producer receipt", inputs.final_receipt),
        ("final candidate catalog TSV", inputs.candidate_catalog),
        ("final candidate witness-pair TSV", inputs.candidate_witness_pairs),
        (CLAIM_CONTRACT_SOURCE_ROLE, inputs.claim_contract),
        ("final dataset builder code", final_builder_path),
        ("layered ClairS/LongPhase-S manifest", layered_path),
        ("post-run output reconciliation", inputs.output_reconciliation),
        ("post-run input immutability", inputs.post_immutability_audit),
        ("latest tree-input audit", inputs.tree_input_audit),
        ("extraction reference identity audit", inputs.reference_identity_audit),
    ]
    if final_release_validation.get("pass") is True:
        source_paths.extend(
            (
                ("signed Task-B final release receipt", inputs.final_release_receipt),
                ("Task-B final release Ed25519 signature", inputs.final_release_signature),
                ("Task-B final release Ed25519 public key", inputs.final_release_public_key),
            )
        )
    final_input_role_names = {
        "manifest": "frozen all-sSNV manifest",
        "screen_sites": "current terminal screen site table",
        "screen_assignments": "current terminal screen stable assignments",
        "screen_summary": "current screen summary",
        "screen_receipt": "current terminal screen receipt",
        "cooccurrence_sites": "cooccurrence site table",
        "cooccurrence_pairs": "cooccurrence pair table",
        "cooccurrence_summary": "cooccurrence summary",
        "cooccurrence_receipt": "cooccurrence producer receipt",
        "tumor_ref_sites": "tumor-REF control site table",
        "tumor_ref_summary": "tumor-REF control summary",
        "tumor_ref_receipt": "tumor-REF producer receipt",
        "primary_artifact_audit_pre": "stable-primary pre-consumer audit",
        "primary_artifact_audit_post": "stable-primary post-consumer audit",
        "tumor_ref_source_identity_receipt": (
            "tumor-REF bounded retrospective source identity receipt"
        ),
    }
    for role, path in final_input_paths.items():
        source_role = final_input_role_names.get(role)
        if source_role is None:
            schema = "non-JSON artifact"
            try:
                candidate_payload = json.loads(path.read_text(encoding="utf-8"))
                if isinstance(candidate_payload, Mapping):
                    schema = str(candidate_payload.get("schema_name", "JSON artifact"))
            except (OSError, UnicodeDecodeError, json.JSONDecodeError):
                pass
            source_role = f"final downstream input {role}: {schema}"
        source_paths.append((source_role, path))
    source_paths.extend(recovery_input_paths.items())
    if source_attestation.get("release_gate_pass") is True:
        source_paths.append(
            (
                "tumor-REF during-execution source identity snapshot",
                source_attestation["_snapshot_path"],
            )
        )
    source_paths.extend(
        (
            f"producer BAM-output receipt {row['sample']}",
            Path(str(row["receipt_path"])),
        )
        for row in producer_bam_audit["rows"]
    )
    if inputs.earlier_fp_report is not None:
        source_paths.append(("earlier FP-only report", inputs.earlier_fp_report))
    role_counts: dict[str, int] = {}
    unique_source_paths: list[tuple[str, Path]] = []
    for role, path in source_paths:
        role_counts[role] = role_counts.get(role, 0) + 1
        unique_role = role if role_counts[role] == 1 else f"{role} #{role_counts[role]}"
        unique_source_paths.append((unique_role, path))
    identities = {
        role: file_identity(role, path, repo_root)
        for role, path in unique_source_paths
    }
    source_rows = [
        {
            "role": identity["role"],
            "path": identity["path"],
            "display_path": identity["display_path"],
            "size_bytes": identity["size_bytes"],
            "sha256": identity["sha256"],
        }
        for identity in identities.values()
    ]
    history_rows: list[dict[str, Any]] = []
    if inputs.earlier_fp_report is not None:
        history = identities["earlier FP-only report"]
        history_rows.append(
            {
                "scope": "earlier FP-only",
                "question": "specificity; not all-sSNV prevalence",
                "path": history["path"],
                "sha256": history["sha256"],
            }
        )
    markdown_sections = build_report_sections(
        claim_rows,
        stratum_rows,
        truth_matrix_rows,
        flow_rows,
        audit_rows,
        technical_rows,
        case_data,
        source_rows,
        history_rows,
        m1_operational,
        m2_evaluability,
        background_control_gate,
        source_attestation,
    )
    artifact_source_rows = [
        {**row, "display_path": row["path"]}
        for row in source_rows
    ]
    artifact_sections = build_report_sections(
        claim_rows,
        stratum_rows,
        truth_matrix_rows,
        flow_rows,
        audit_rows,
        technical_rows,
        case_data,
        artifact_source_rows,
        history_rows,
        m1_operational,
        m2_evaluability,
        background_control_gate,
        source_attestation,
        include_inline_tables=False,
    )
    metadata = (
        "<!--\n"
        f"建立時間: {generated_at}\n"
        "目標: 依 schema 2.0.0 輸出 all-sSNV claim ladder 完整報告\n"
        "處理範圍: 7 datasets / 6 biological samples / chr1-22 / 469,849 sSNV\n"
        f"關聯檔案: {display_path(inputs.final_dataset, repo_root)}\n"
        "-->\n\n"
    )
    report_text = metadata + "\n\n".join(body for _, body in markdown_sections) + "\n"

    m2_evaluability_rows = [
        {
            "record_type": "summary",
            "site": "ALL_M1_PASS",
            "count": m2_evaluability[
                "n_group_count_exceeds_planning_model_maximum"
            ],
            "observed_methyl_groups": None,
            "maximum_supported_methyl_groups": 10,
            "m1_status": "PASS retained when >10 groups",
            "m2_status": "NOT_EVALUABLE when >10 groups",
            "g1_status": "NOT_RUN",
            "g2_status": "NOT_RUN",
            "b1_status": "NOT_RUN",
            "claim_behavior": "M1 PASS; M2 NE; G1/G2/B1 NOT_RUN",
            "reason": "2-10-group planning model boundary",
        },
        *[
            {
                "record_type": "observed_example",
                "site": (
                    f"{row['dataset']} {row['chrom']}:{row['pos']} "
                    f"{row['ref']}>{row['alt']}"
                ),
                "count": 1,
                "observed_methyl_groups": row["observed_methyl_groups"],
                "maximum_supported_methyl_groups": row[
                    "maximum_supported_methyl_groups"
                ],
                "m1_status": row["m1_status"],
                "m2_status": row["m2_status"],
                "g1_status": row["g1_status"],
                "g2_status": row["g2_status"],
                "b1_status": row["b1_status"],
                "claim_behavior": "M1 PASS; M2 NE; G1/G2/B1 NOT_RUN",
                "reason": row["reason"],
            }
            for row in m2_evaluability[
                "group_count_exceeds_planning_model_examples"
            ]
        ],
    ]
    datasets: dict[str, list[dict[str, Any]]] = {
        "claim_ladder": claim_rows,
        "claim_chart": claim_chart_rows,
        "overview_case_chart": overview_case_chart_rows,
        "stratum_claim_metrics": stratum_rows,
        "focal_partner_truth_matrix": truth_matrix_rows,
        "flow": flow_rows,
        "audit_status": audit_rows,
        "technical_replication": technical_rows,
        "m2_evaluability": m2_evaluability_rows,
        "producer_bam_output": [
            {
                **row,
                "receipt_path": identities[
                    f"producer BAM-output receipt {row['sample']}"
                ]["path"],
            }
            for row in producer_bam_audit["rows"]
        ],
        **{name: list(rows) for name, rows in case_data.items()},
        "source_inventory": [
            {field: row[field] for field in ("role", "path", "size_bytes", "sha256")}
            for row in source_rows
        ],
        "history": history_rows,
    }
    final_identity = identities["final schema 2.0.0 dataset"]
    layered_identity = identities["layered ClairS/LongPhase-S manifest"]
    site_identity = identities["cooccurrence site table"]
    pair_identity = identities["cooccurrence pair table"]
    tree_identity = identities["latest tree-input audit"]
    script_identity = file_identity("report builder", Path(__file__), repo_root)
    sources = [
        source_spec(
            "src-claims-view",
            "Validated schema 2.0.0 claim ladder",
            final_identity,
            rows=claim_rows,
            columns=tuple(claim_rows[0]),
            relation="claim_ladder_snapshot",
            description=(
                "Executed bounded VALUES projection after claim-status, denominator, "
                "and release-contract reconciliation."
            ),
            filters=("task_type=B comprehensive validation", "claims=M1 through L2"),
            metric_definitions=(
                "M1 dataset ratio = operational FLAGGED / all dataset-sites; NOT_FLAGGED includes technically non-evaluable sites and is not a biological negative",
                "M1 biological-site ratio = descriptive FLAGGED / biological-site aggregation; not biological prevalence",
                "M2-L2 ratio = PASS / claim-specific evaluable denominator",
            ),
            generated_at=generated_at,
        ),
        source_spec(
            "src-strata-view",
            "Validated per-sample and focal-truth-stratum metrics",
            final_identity,
            rows=stratum_rows,
            columns=tuple(stratum_rows[0]),
            relation="stratum_claim_snapshot",
            description="Executed bounded VALUES projection after pooled/sample/truth reconciliation.",
            filters=("claims=M1,M2,G1,G2", "truth labels are post hoc"),
            metric_definitions=(
                "M1 ratio = operational FLAGGED / all sites in the stratum; not biological prevalence",
                "M2-G2 ratio = PASS / claim-specific evaluable denominator",
            ),
            generated_at=generated_at,
        ),
        source_spec(
            "src-truth-matrix-view",
            "Validated pair-level focal x partner truth matrix",
            final_identity,
            rows=truth_matrix_rows,
            columns=tuple(truth_matrix_rows[0]),
            relation="focal_partner_truth_matrix_snapshot",
            description=(
                "Executed bounded VALUES projection separating focal-site truth from partner truth."
            ),
            filters=("truth labels are post hoc", "all pair rows and G1 formal pair rows"),
            metric_definitions=("site-level truth strata use focal truth only",),
            generated_at=generated_at,
        ),
        source_spec(
            "src-technical-view",
            "HCC1395 technical replication metrics",
            final_identity,
            rows=technical_rows,
            columns=tuple(technical_rows[0]),
            relation="technical_replication_snapshot",
            description="Executed bounded VALUES projection of exact shared-pair technical concordance.",
            filters=("HCC1395 and HCC1395_DORADO only",),
            metric_definitions=("biological n remains 1",),
            generated_at=generated_at,
        ),
        source_spec(
            "src-flow-view",
            "Current ClairS to LongPhase-S processing flow",
            layered_identity,
            rows=flow_rows,
            columns=tuple(flow_rows[0]),
            relation="processing_flow_snapshot",
            description="Executed bounded VALUES projection of the layered and terminal-tag flow contracts.",
            filters=("same-run LongPhase-S PASS", "external HP/PS sidecar"),
            tables_used=(
                layered_identity["path"],
                identities["current screen summary"]["path"],
                *(
                    identities[f"producer BAM-output receipt {sample}"]["path"]
                    for sample in EXPECTED_DATASETS
                ),
            ),
            generated_at=generated_at,
        ),
        source_spec(
            "src-case-selection-view",
            "Independent four-view case selection registry",
            pair_identity,
            rows=case_data["case_selection_registry"],
            columns=tuple(case_data["case_selection_registry"][0]),
            relation="case_selection_registry_snapshot",
            description="Executed bounded VALUES projection of four independently defined selection roles.",
            filters=("aggregate", "canonical oracle", "global extreme", "well-explained"),
            tables_used=(
                final_identity["path"],
                pair_identity["path"],
            ),
            generated_at=generated_at,
        ),
        source_spec(
            "src-case-summary-view",
            "Primary report-case summary",
            pair_identity,
            rows=case_data["case_summary"],
            columns=tuple(case_data["case_summary"][0]),
            relation="case_summary_snapshot",
            description="Executed bounded VALUES projection of the primary report case or structured N/A state.",
            filters=(f"selection_mode={case_mode}",),
            tables_used=(pair_identity["path"], site_identity["path"]),
            generated_at=generated_at,
        ),
        source_spec(
            "src-case-views",
            "Independent selected-case views and structured N/A states",
            pair_identity,
            rows=case_data["case_view_summary"],
            columns=tuple(case_data["case_view_summary"][0]),
            relation="case_view_summary_snapshot",
            description=(
                "Executed bounded VALUES projection of primary, canonical, extreme, "
                "well-explained, and evaluated-negative case roles."
            ),
            filters=(
                "canonical target is never substituted",
                "negative_result means evaluated endpoint FAIL only",
            ),
            metric_definitions=(
                "N/A is not a biological negative result",
                "topology is pairwise posthoc context only",
            ),
            tables_used=(pair_identity["path"], site_identity["path"], final_identity["path"]),
            generated_at=generated_at,
        ),
        source_spec(
            "src-case-group-view",
            "Primary case group by partner evidence",
            pair_identity,
            rows=case_data["case_group_partner"],
            columns=tuple(case_data["case_group_partner"][0]),
            relation="case_group_partner_snapshot",
            description="Executed bounded VALUES projection of methyl-group by partner R/A counts.",
            filters=(f"selection_mode={case_mode}",),
            metric_definitions=("Cells are read counts; N/A is not a negative result.",),
            generated_at=generated_at,
        ),
        source_spec(
            "src-overview-case-view",
            "Claim funnel and actual case molecular composition",
            final_identity,
            rows=overview_case_chart_rows,
            columns=tuple(overview_case_chart_rows[0]),
            relation="overview_case_chart_snapshot",
            description=(
                "Executed bounded projection of the M1 operational all-dataset-site flagged yield, "
                "the descriptive non-prevalence M1 biological-site aggregation, M2-G2 "
                "PASS/evaluable proportions, and case partner-read composition."
            ),
            filters=("claims=M1,M2,G1,G2", f"case_selection_mode={case_mode}"),
            metric_definitions=(
                "M1 dataset-site proportion = FLAGGED / all dataset-sites",
                "M1 biological-site proportion = descriptive FLAGGED / biological-site aggregation; not biological prevalence",
                "M2-G2 claim proportion = PASS / claim-specific evaluable denominator",
                "case proportion = partner R or A reads / R+A reads within methyl group",
            ),
            tables_used=(final_identity["path"], pair_identity["path"]),
            generated_at=generated_at,
        ),
        source_spec(
            "src-case-four-state-view",
            "Primary case four-state evidence",
            pair_identity,
            rows=case_data["case_four_state"],
            columns=tuple(case_data["case_four_state"][0]),
            relation="case_four_state_snapshot",
            description="Executed bounded VALUES projection of RR/AR/RA/AA/O/X fixed-error-model counts.",
            filters=(f"selection_mode={case_mode}",),
            metric_definitions=("Relation status is compatibility, not identified ancestry.",),
            generated_at=generated_at,
        ),
        source_spec(
            "src-case-view-four-state",
            "Selected-view four-state fixed-error models",
            pair_identity,
            rows=case_data["case_view_four_state_models"],
            columns=tuple(case_data["case_four_state_models"][0]),
            relation="case_view_four_state_models_snapshot",
            description=(
                "Executed bounded VALUES projection of focal-ancestor, partner-ancestor, "
                "and branching fixed-error-model checks for each selected case role."
            ),
            filters=("selected deterministic case roles only",),
            metric_definitions=(
                "2% fixed error ceiling and 95% one-sided upper bound",
                "model compatibility is not identified cellular ancestry",
            ),
            generated_at=generated_at,
        ),
        source_spec(
            "src-case-joint-view",
            "Primary case complete-read joint signature",
            site_identity,
            rows=case_data["case_joint_signature"],
            columns=tuple(case_data["case_joint_signature"][0]),
            relation="case_joint_signature_snapshot",
            description="Executed bounded VALUES projection of complete-read joint signatures by methyl group.",
            filters=(f"selection_mode={case_mode}",),
            metric_definitions=("Cells are complete-read counts; N/A is not a negative result.",),
            generated_at=generated_at,
        ),
        source_spec(
            "src-audit-view",
            "Current pipeline and release-blocking audits",
            tree_identity,
            rows=audit_rows,
            columns=tuple(audit_rows[0]),
            relation="audit_status_snapshot",
            description="Executed bounded VALUES projection of independently supplied PASS audits.",
            filters=("all release-blocking audits must PASS",),
            tables_used=tuple(
                identities[role]["path"]
                for role in (
                    "current screen summary",
                    "post-run output reconciliation",
                    "post-run input immutability",
                    "latest tree-input audit",
                    "extraction reference identity audit",
                )
            )
            + tuple(
                identities[f"producer BAM-output receipt {sample}"]["path"]
                for sample in EXPECTED_DATASETS
            ),
            generated_at=generated_at,
        ),
        source_spec(
            "src-source-view",
            "Report source inventory",
            script_identity,
            rows=datasets["source_inventory"],
            columns=tuple(datasets["source_inventory"][0]),
            relation="source_inventory_snapshot",
            description="Executed bounded VALUES projection of source paths, sizes, and hashes.",
            tables_used=tuple(row["path"] for row in datasets["source_inventory"]),
            generated_at=generated_at,
        ),
    ]
    if history_rows:
        history_identity = identities["earlier FP-only report"]
        sources.append(
            source_spec(
                "src-history-view",
                "Earlier FP-only scope reference",
                history_identity,
                rows=history_rows,
                columns=tuple(history_rows[0]),
                relation="history_scope_snapshot",
                description="Executed bounded VALUES projection of the historical scope boundary.",
                filters=("specificity question only", "not an all-sSNV prevalence estimate"),
                generated_at=generated_at,
            )
        )
    artifact_payload = build_artifact(
        artifact_sections,
        datasets,
        sources,
        generated_at,
        include_history=bool(history_rows),
    )
    validate_artifact_shape(artifact_payload)
    artifact_text = json.dumps(artifact_payload, ensure_ascii=False, indent=2) + "\n"
    snapshot_rows = sum(len(rows) for rows in datasets.values())

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=False, exist_ok=False)
    report_path = output_dir / "report.md"
    artifact_path = output_dir / "artifact.json"
    report_path.write_text(report_text, encoding="utf-8")
    artifact_path.write_text(artifact_text, encoding="utf-8")
    receipt_path = output_dir / "report_build_receipt.json"
    report_receipt = {
        "schema_name": REPORT_RECEIPT_SCHEMA,
        "schema_version": REPORT_RECEIPT_SCHEMA_VERSION,
        "created_at_utc": generated_at,
        "task_type": (
            "B_comprehensive_validation"
            if final_release_validation.get("pass") is True
            else "NON_RELEASE_TEST_FIXTURE_REPORT"
        ),
        "formal_task_b_release_eligible": final_release_validation.get("pass") is True,
        "pass_semantics": "report_build_integrity_only_not_scientific_confirmation",
        "inputs": {role: identity for role, identity in identities.items()},
        "claim_contract": identities[CLAIM_CONTRACT_SOURCE_ROLE],
        "code": {
            "report_builder": script_identity,
            "final_report_dataset_builder": identities["final dataset builder code"],
        },
        "outputs": {
            "report_md": file_identity("rendered Markdown report", report_path, repo_root),
            "artifact_json": file_identity(
                "canonical report artifact", artifact_path, repo_root
            ),
        },
        "snapshot_rows": snapshot_rows,
        "snapshot_dataset_rows": {
            name: len(rows) for name, rows in datasets.items()
        },
        "artifact_presentation_scope": (
            "complete_portable_technical_report_with_visible_narrative_native_charts_tables_"
            "four_views_evaluable_negative_structured_na_and_actual_case_evidence"
        ),
        "validations": {
            "detached_final_result_signature": (
                "PASS"
                if final_release_validation.get("pass") is True
                else "NOT_INCLUDED_NON_RELEASE_TEST_FIXTURE"
            ),
            "final_dataset_receipt_identity": "PASS",
            "claim_contract_complete_vocabulary": "PASS",
            "claim_contract_release_state_semantics": (
                "PASS_RELEASE_V5"
                if source_attestation.get("release_gate_pass") is True
                else "PASS_INTERMEDIATE_V3"
            ),
            "claim_contract_sha256_recorded": "PASS",
            "artifact_shape": "PASS",
            "extraction_reference_identity": "PASS",
            "terminal_screen_receipt": "PASS",
            "recovery_provenance_artifacts": (
                "PASS" if recovery_input_paths else "NOT_APPLICABLE_NON_RECOVERY_RUN"
            ),
            "post_immutability_is_fresh_after_terminal_screen": "PASS",
            "all_final_input_artifacts_sha256": "PASS",
            "tumor_ref_bounded_source_identity_release_gate": (
                "PASS"
                if source_attestation.get("release_gate_pass") is True
                else "NOT_INCLUDED_INTERMEDIATE_ONLY"
            ),
            "canonical_target_not_substituted": "PASS",
            "evaluable_negative_is_statistical_not_biological": "PASS",
            "output_overwrite_allowed": False,
        },
        "pass": True,
    }
    receipt_path.write_text(
        json.dumps(report_receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    outputs = {
        "report.md": report_path,
        "artifact.json": artifact_path,
        "report_build_receipt.json": receipt_path,
    }
    for path in outputs.values():
        path.chmod(0o444)
    print(
        json.dumps(
            {
                "output_dir": str(output_dir),
                "outputs": {
                    name: {
                        "path": str(path),
                        "size_bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                    for name, path in outputs.items()
                },
                "snapshot_rows": snapshot_rows,
                "report_sha256_from_memory": sha256_text(report_text),
                "html_generated": False,
                "pass": True,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-dataset", type=Path, required=True)
    parser.add_argument("--final-receipt", type=Path, required=True)
    parser.add_argument("--final-release-receipt", type=Path, required=True)
    parser.add_argument("--final-release-signature", type=Path, required=True)
    parser.add_argument("--final-release-public-key", type=Path, required=True)
    parser.add_argument("--candidate-catalog", type=Path, required=True)
    parser.add_argument("--candidate-witness-pairs", type=Path, required=True)
    parser.add_argument("--claim-contract", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--screen-summary", type=Path, required=True)
    parser.add_argument("--cooccurrence-sites", type=Path, required=True)
    parser.add_argument("--cooccurrence-pairs", type=Path, required=True)
    parser.add_argument("--cooccurrence-summary", type=Path, required=True)
    parser.add_argument("--output-reconciliation", type=Path, required=True)
    parser.add_argument(
        "--post-immutability-audit",
        "--immutability-audit",
        dest="post_immutability_audit",
        type=Path,
        required=True,
    )
    parser.add_argument("--tree-input-audit", type=Path, required=True)
    parser.add_argument("--reference-identity-audit", type=Path, required=True)
    parser.add_argument("--earlier-fp-report", type=Path)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    inputs = ReportInputs(
        final_dataset=args.final_dataset,
        final_receipt=args.final_receipt,
        final_release_receipt=args.final_release_receipt,
        final_release_signature=args.final_release_signature,
        final_release_public_key=args.final_release_public_key,
        candidate_catalog=args.candidate_catalog,
        candidate_witness_pairs=args.candidate_witness_pairs,
        claim_contract=args.claim_contract,
        manifest=args.manifest,
        screen_summary=args.screen_summary,
        cooccurrence_sites=args.cooccurrence_sites,
        cooccurrence_pairs=args.cooccurrence_pairs,
        cooccurrence_summary=args.cooccurrence_summary,
        output_reconciliation=args.output_reconciliation,
        post_immutability_audit=args.post_immutability_audit,
        tree_input_audit=args.tree_input_audit,
        reference_identity_audit=args.reference_identity_audit,
        earlier_fp_report=args.earlier_fp_report,
    )
    build_outputs(inputs, args.output_dir, repo_root=args.repo_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
