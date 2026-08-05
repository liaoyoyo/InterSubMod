#!/usr/bin/env python3
"""Build a full positional-singleton M1/M2 supplemental audit."""

from __future__ import annotations

import argparse
import ctypes
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import stat
import subprocess
import sys
from types import ModuleType
from typing import Any, Iterable, Mapping, Sequence
import uuid


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import release_source_authority as SOURCE_AUTHORITY  # noqa: E402


SCHEMA_NAME = "intersubmod.positional_singleton_methyl_multigroup_audit"
SCHEMA_VERSION = "2.0.0"
ASSIGNMENT_SCHEMA = "intersubmod.all_ssnv_stable_multigroup_read_assignment"
POSITIONAL_GAP_BP = 50_000
EXPECTED_DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
EXPECTED_TRUTH_LABELS = ("FP", "TP", "UNASSESSED")
EXPECTED_CHROMS = frozenset(f"chr{index}" for index in range(1, 23))
EXPECTED_ALLELES = frozenset("ACGT")
EXPECTED_BRANCHES = frozenset(
    {"max_snv_excluded", "positional_singleton", "retained"}
)
EXPECTED_TASK_TYPE = "B_comprehensive_validation"
EXPECTED_SCOPE = (
    "7 datasets; chr1-22; LongPhase-S recalibrated FILTER=PASS; "
    "biallelic sSNV"
)
EXPECTED_CLAIM_SHA256 = (
    "da8778af6bcc9b1e8b2887eb2bcc87eca83d909be32a819ec5eb8f5f12c9f2af"
)
EXPECTED_AUTHORITY_ID = SOURCE_AUTHORITY.AUTHORITY_ID
EXPECTED_AUTHORIZED_SOURCE_ROLES = frozenset(
    {
        "cn_ccf_annotator",
        "completion_runner",
        "cooccurrence_analyzer",
        "cooccurrence_preflight",
        "cooccurrence_release_finalizer",
        "cooccurrence_release_runner",
        "final_dataset_builder",
        "final_result_release_finalizer",
        "focal_alt_cluster_lib",
        "frozen_input_auditor",
        "independent_m2_auditor",
        "latest_tag_join",
        "m2_screen_gate",
        "matched_normal_analyzer",
        "matched_normal_runner",
        "portable_layout_qa",
        "portable_report_deliverer",
        "primary_artifact_auditor",
        "release_source_authority_validator",
        "report_builder",
        "ssnv_cooccurrence_lib",
        "strict_producer",
        "tumor_ref_source_verifier",
    }
)
EXPECTED_M2_CONTRACT = (
    "m2-measured-axis-v4_asymmetric-positive-confound-and-observed-"
    "categorical-levels"
)
EXPECTED_M2_CHECKS = frozenset(
    {
        "aligned_below_power_evaluable_matches_expected",
        "all_rows_match_expected",
        "axis_indeterminate_matches_expected",
        "eligible_matches_expected",
        "evaluable_ineligible_matches_expected",
        "group_count_gt10_matches_expected",
        "group_distribution_conserved",
        "group_limit_examples_reconcile",
        "m1_outcome_categories_conserved",
        "m1_stable_rows_match_expected",
        "stable_screen_assignment_keys_exact",
    }
)
EXPECTED_TREE_CHECKS = frozenset(
    {
        "input_manifest_pass",
        "layered_manifest_sha_matches_frozen_record",
        "same_seven_datasets",
        "producer_seven_of_seven_pass",
        "tree_input_contract_is_latest_longphase_s_pass",
        "read_tag_contract_requires_exact_external_join",
        "focal_scope_is_469849",
    }
)
EXPECTED_M2_AXIS_SPECS = (
    ("hp_exact", "categorical", "v", 0.30),
    ("hp_family", "categorical", "v", 0.30),
    ("strand", "categorical", "v", 0.30),
    ("start", "continuous", "eta2", 0.14),
    ("end", "continuous", "eta2", 0.14),
    ("length", "continuous", "eta2", 0.14),
    ("mapq", "continuous", "eta2", 0.14),
    ("cpg_called", "continuous", "eta2", 0.14),
)
EXPECTED_M2_CATEGORY_LEVEL_CEILINGS = {
    "hp_exact": 7,
    "hp_family": 5,
    "strand": 2,
}
EXPECTED_M2_GLOBAL_COUNTS = {
    "all_rows": 469_849,
    "eligible": 919,
    "evaluable_ineligible": 948,
    "m1_stable_rows": 102_842,
    "not_evaluable_axis_indeterminate": 100_974,
    "not_evaluable_group_count_gt10": 1,
}
EXPECTED_M2_LOGIC_INDEPENDENCE = {
    "assignment_categories_and_coarse_group_counts_recomputed": True,
    "production_gate_functions_called": False,
    "production_gate_imported": False,
    "screen_effect_and_p_values_reused_as_frozen_inputs": True,
}
EXPECTED_GLOBAL_TRUTH_COUNTS = {
    "truth_tp": 335_296,
    "truth_fp": 7_745,
    "truth_unassessed": 126_808,
}
REQUIRED_CLAIM_FRAGMENTS = (
    "`terminal-claim-contract-v5`",
    "Primary key：`(dataset, CHROM, POS, REF, ALT)`",
    EXPECTED_M2_CONTRACT,
    "L1前只能稱molecular-haplotype或read-level epigenetic evidence",
)
STRICT_NOT_EVALUATED_REASON = (
    "NOT_EVALUATED_AT_M1_SCREEN: R1 strict methyl-partition robustness is "
    "downstream-only and requires a G2 multi-marker molecular-haplotype base "
    "candidate; not evaluated is not failure"
)
EXPECTED_COUNTS = {
    "all_dataset_sites": 469_849,
    "m1_stable_all": 102_842,
    "positional_singleton": 50_432,
    "singleton_m1_evaluable": 48_347,
    "singleton_m1_not_evaluable": 2_085,
    "singleton_m1_flagged": 5_961,
    "singleton_m2_pass": 30,
    "singleton_m2_fail": 18,
    "singleton_m2_axis_indeterminate": 5_913,
    "singleton_m2_group_count_gt10": 0,
    "singleton_m2_not_run": 44_471,
    "singleton_min_finite_nearest_gap_bp": 50_003,
}
EXPECTED_BRANCH_COUNTS = {
    "max_snv_excluded": 225_268,
    "positional_singleton": 50_432,
    "retained": 194_149,
}
EXPECTED_M1_NOT_EVALUABLE_REASONS = {
    "INCOMPLETE_DISTANCE_BELOW_MINIMUM": 16,
    "INSUFFICIENT_MATRIX_JOINED_FOCAL_ALT_READS": 2_069,
}
EXPECTED_M2_FAIL_SUBTYPES = {
    "HP_AXIS_CONFOUND": 0,
    "NOT_PHASE_ANCHORED_ROBUST": 14,
    "TECHNICAL_AXIS_CONFOUND": 4,
}

SiteKey = tuple[str, str, int, str, str]


class SingletonAuditError(RuntimeError):
    """Raised when the frozen singleton contract does not reconcile."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(
    path: Path,
    *,
    include_mode: bool = True,
    public_path: Path | None = None,
) -> dict[str, Any]:
    resolved = path.expanduser().resolve(strict=True)
    info = resolved.stat()
    record: dict[str, Any] = {
        "path": str(
            public_path.expanduser().resolve()
            if public_path is not None
            else resolved
        ),
        "size_bytes": info.st_size,
        "mtime_ns": info.st_mtime_ns,
        "sha256": sha256(resolved),
    }
    if include_mode:
        record["mode"] = oct(stat.S_IMODE(info.st_mode))
    return record


def open_text(path: Path) -> Any:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def parse_bool(value: Any, *, field: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes"}:
        return True
    if normalized in {"0", "false", "no"}:
        return False
    raise SingletonAuditError(f"Invalid boolean for {field}: {value!r}")


def parse_optional_float(value: Any) -> float | None:
    if value in {None, ""}:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise SingletonAuditError(f"Invalid optional float: {value!r}") from error
    if not math.isfinite(parsed):
        raise SingletonAuditError(f"Non-finite optional float: {value!r}")
    return parsed


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SingletonAuditError(f"Unable to read {label}: {path}") from error
    if not isinstance(value, dict):
        raise SingletonAuditError(f"{label} must be a JSON object: {path}")
    return value


def require_schema(
    record: Mapping[str, Any],
    *,
    label: str,
    schema_name: str,
    schema_version: str,
) -> None:
    observed = (record.get("schema_name"), record.get("schema_version"))
    expected = (schema_name, schema_version)
    if observed != expected:
        raise SingletonAuditError(
            f"{label} schema drift: observed={observed!r}, expected={expected!r}"
        )


def require_true_checks(
    checks: Any, *, label: str, expected_keys: frozenset[str]
) -> None:
    if not isinstance(checks, dict) or frozenset(checks) != expected_keys:
        raise SingletonAuditError(f"{label} check key set drift")
    failed = sorted(key for key, value in checks.items() if value is not True)
    if failed:
        raise SingletonAuditError(f"{label} failed checks: {failed}")


def require_identity(
    path: Path,
    recorded: Any,
    *,
    label: str,
    require_size: bool = True,
) -> None:
    if not isinstance(recorded, dict):
        raise SingletonAuditError(f"{label} identity is not an object")
    observed = artifact(path)
    required = ("path", "sha256", "size_bytes") if require_size else ("path", "sha256")
    mismatches = {
        field: {"observed": observed[field], "recorded": recorded.get(field)}
        for field in required
        if observed[field] != recorded.get(field)
    }
    if mismatches:
        raise SingletonAuditError(f"{label} identity drift: {mismatches}")


def require_read_only(
    identities: Mapping[str, Mapping[str, Any]], names: Iterable[str]
) -> None:
    wrong = {
        name: identities[name].get("mode")
        for name in names
        if identities[name].get("mode") != "0o444"
    }
    if wrong:
        raise SingletonAuditError(f"Protected input mode drift: {wrong}")


def fsync_file(path: Path) -> None:
    with path.open("rb") as handle:
        os.fsync(handle.fileno())


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def rate_record(
    numerator: int, denominator: int, denominator_definition: str
) -> dict[str, Any]:
    if denominator <= 0:
        raise SingletonAuditError(
            f"Rate denominator must be positive: {denominator_definition}"
        )
    return {
        "numerator": numerator,
        "denominator": denominator,
        "denominator_definition": denominator_definition,
        "value": numerator / denominator,
    }


def rename_no_replace(source: Path, target: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise SingletonAuditError("renameat2 is required for no-replace publication")
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    result = renameat2(
        -100,
        os.fsencode(source),
        -100,
        os.fsencode(target),
        1,
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise FileExistsError(
            error_number,
            f"No-replace publication failed: {os.strerror(error_number)}",
            str(target),
        )


def verify_source_authority(
    *,
    paths: Mapping[str, Path],
    identities: Mapping[str, Mapping[str, Any]],
    primary_audit: Mapping[str, Any],
) -> dict[str, Any]:
    canonical_paths = {
        "source_authority": SOURCE_AUTHORITY.AUTHORITY_PATH,
        "source_authority_approval": SOURCE_AUTHORITY.APPROVAL_PATH,
        "source_authority_signature": SOURCE_AUTHORITY.APPROVAL_SIGNATURE_PATH,
        "source_authority_public_key": SOURCE_AUTHORITY.PUBLIC_KEY_PATH,
    }
    path_mismatches = {
        name: {
            "observed": str(paths[name].resolve()),
            "expected": str(expected.resolve()),
        }
        for name, expected in canonical_paths.items()
        if paths[name].resolve() != expected.resolve()
    }
    if path_mismatches:
        raise SingletonAuditError(
            f"Legacy or noncanonical source-authority path rejected: {path_mismatches}"
        )
    try:
        canonical = SOURCE_AUTHORITY.validate_release_source_authority()
    except (SOURCE_AUTHORITY.SourceAuthorityError, OSError, ValueError) as error:
        raise SingletonAuditError(
            f"Canonical v5 source-authority validation failed: {error}"
        ) from error
    if (
        canonical.get("pass") is not True
        or canonical.get("authority_id") != EXPECTED_AUTHORITY_ID
    ):
        raise SingletonAuditError("Canonical v5 source authority is not approved")
    expected_source_set = str(canonical["source_set_sha256"])
    expected_head = str(canonical["authorized_git_head"])
    protected = (
        "source_authority",
        "source_authority_approval",
        "source_authority_signature",
        "source_authority_public_key",
    )
    require_read_only(identities, protected)
    expected_records = {
        "source_authority": canonical["authority_manifest"],
        "source_authority_approval": canonical["detached_approval"],
        "source_authority_signature": canonical[
            "detached_approval_signature"
        ],
        "source_authority_public_key": canonical["external_public_key"],
    }
    identity_keys = ("path", "size_bytes", "sha256", "mode")
    identity_key_set = frozenset(identity_keys)

    def canonical_identity(record: Mapping[str, Any]) -> dict[str, Any]:
        return {key: record.get(key) for key in identity_keys}

    malformed_expected = {
        name: (
            sorted(expected)
            if isinstance(expected, Mapping)
            else type(expected).__name__
        )
        for name, expected in expected_records.items()
        if not isinstance(expected, Mapping) or set(expected) != identity_key_set
    }
    if malformed_expected:
        raise SingletonAuditError(
            "Canonical source-authority identity schema drift: "
            f"{malformed_expected}"
        )
    wrong_identities = {
        name: identities[name]
        for name, expected in expected_records.items()
        if canonical_identity(identities[name]) != dict(expected)
    }
    if wrong_identities:
        raise SingletonAuditError(
            f"Source-authority trust-anchor drift: {wrong_identities}"
        )
    try:
        with SOURCE_AUTHORITY.BoundArtifactReader() as reader:
            openssl_fd, _, openssl_artifact = reader.open(
                SOURCE_AUTHORITY.OPENSSL_PATH, include_mode=True
            )
            _, authority_bytes, bound_authority = reader.open(
                paths["source_authority"], include_mode=True
            )
            approval_fd, approval_bytes, bound_approval = reader.open(
                paths["source_authority_approval"], include_mode=True
            )
            signature_fd, _, bound_signature = reader.open(
                paths["source_authority_signature"], include_mode=True
            )
            public_key_fd, _, bound_public_key = reader.open(
                paths["source_authority_public_key"], include_mode=True
            )
            bound_identities = {
                "source_authority": bound_authority,
                "source_authority_approval": bound_approval,
                "source_authority_signature": bound_signature,
                "source_authority_public_key": bound_public_key,
            }
            if any(
                canonical_identity(bound_identities[name])
                != canonical_identity(identities[name])
                for name in protected
            ):
                raise SingletonAuditError(
                    "Source-authority identity changed before bound verification"
                )
            if openssl_artifact["sha256"] != SOURCE_AUTHORITY.OPENSSL_SHA256:
                raise SingletonAuditError(
                    "Source-authority OpenSSL verifier drift"
                )
            if not SOURCE_AUTHORITY.verify_ed25519_signature_fds(
                openssl_fd=openssl_fd,
                data_fd=approval_fd,
                public_key_fd=public_key_fd,
                signature_fd=signature_fd,
            ):
                raise SingletonAuditError(
                    "Source-authority Ed25519 signature verification failed"
                )
            try:
                authority = json.loads(authority_bytes.decode("utf-8"))
                approval = json.loads(approval_bytes.decode("utf-8"))
            except (UnicodeDecodeError, json.JSONDecodeError) as error:
                raise SingletonAuditError(
                    "Bound source-authority JSON is malformed"
                ) from error
            if not isinstance(authority, dict) or not isinstance(approval, dict):
                raise SingletonAuditError(
                    "Bound source-authority records must be JSON objects"
                )
            reader.retain_until_process_exit()
    except SOURCE_AUTHORITY.SourceAuthorityError as error:
        raise SingletonAuditError(
            f"Unable to bind source-authority trust artifacts: {error}"
        ) from error
    require_schema(
        authority,
        label="source authority",
        schema_name="intersubmod.release_source_authority",
        schema_version="1.0.0",
    )
    require_schema(
        approval,
        label="source authority approval",
        schema_name="intersubmod.release_source_authority.approval",
        schema_version="1.0.0",
    )
    for label, record in (("authority", authority), ("approval", approval)):
        if (
            record.get("authority_id") != EXPECTED_AUTHORITY_ID
            or record.get("approval_status") != "APPROVED_FOR_FULL_TASK_B_RUN"
        ):
            raise SingletonAuditError(f"Source {label} approval contract drift")
    require_identity(
        paths["source_authority"],
        approval.get("authority_manifest"),
        label="approval-bound authority manifest",
    )
    require_identity(
        paths["source_authority_public_key"],
        approval.get("public_key"),
        label="approval-bound public key",
    )
    if (
        authority.get("source_set_sha256") != expected_source_set
        or authority.get("repository", {}).get("git_head_at_authorization")
        != expected_head
    ):
        raise SingletonAuditError("Authorized source-set or Git HEAD drift")
    sources = authority.get("sources")
    if (
        not isinstance(sources, dict)
        or frozenset(sources) != EXPECTED_AUTHORIZED_SOURCE_ROLES
    ):
        raise SingletonAuditError("Authorized source role set drift")
    source_lines: list[str] = []
    for role, record in sorted(sources.items()):
        source_path = Path(str(record.get("path", ""))).resolve(strict=True)
        require_identity(source_path, record, label=f"authorized source {role}")
        observed = artifact(source_path)
        if observed["mode"] != "0o444" or record.get("mode") != "0o444":
            raise SingletonAuditError(f"Authorized source is not read-only: {role}")
        source_lines.append(
            "|".join(
                (
                    role,
                    observed["path"],
                    str(observed["size_bytes"]),
                    observed["sha256"],
                    observed["mode"],
                )
            )
        )
    observed_source_set = hashlib.sha256(
        "\n".join(source_lines).encode("utf-8")
    ).hexdigest()
    if observed_source_set != expected_source_set:
        raise SingletonAuditError(
            f"Live authorized source-set drift: {observed_source_set}"
        )
    require_identity(
        paths["claim_contract"],
        authority.get("claim_contract"),
        label="authority claim contract",
    )
    require_identity(
        paths["independent_m2_auditor"],
        sources["independent_m2_auditor"],
        label="authority independent M2 source",
    )
    reviews = approval.get("review_approvals")
    if (
        not isinstance(reviews, list)
        or len(reviews) < 2
        or any(
            review.get("verdict") != "APPROVE"
            or review.get("findings_closed") is not True
            or review.get("reviewed_git_head") != expected_head
            or review.get("reviewed_source_set_sha256")
            != expected_source_set
            for review in reviews
        )
    ):
        raise SingletonAuditError("External source-authority review drift")
    primary_authority = primary_audit.get("source_authority")
    if (
        not isinstance(primary_authority, dict)
        or primary_authority.get("authority_id") != EXPECTED_AUTHORITY_ID
        or primary_authority.get("source_set_sha256") != expected_source_set
        or primary_authority.get("pass") is not True
    ):
        raise SingletonAuditError("Primary audit lost v5 source-authority binding")
    for role, path_name, primary_name in (
        ("authority", "source_authority", "authority_manifest"),
        (
            "approval",
            "source_authority_approval",
            "detached_approval",
        ),
        (
            "signature",
            "source_authority_signature",
            "detached_approval_signature",
        ),
        (
            "public key",
            "source_authority_public_key",
            "external_public_key",
        ),
    ):
        require_identity(
            paths[path_name],
            primary_authority.get(primary_name),
            label=f"primary-audit {role}",
        )
    return {
        "authority_id": EXPECTED_AUTHORITY_ID,
        "authority_manifest": identities["source_authority"],
        "detached_approval": identities["source_authority_approval"],
        "detached_approval_signature": identities[
            "source_authority_signature"
        ],
        "public_key": identities["source_authority_public_key"],
        "source_set_sha256": observed_source_set,
        "authorized_git_head": expected_head,
        "authorized_source_roles": sorted(sources),
        "external_approvals": len(reviews),
        "pass": True,
    }


def screen_key(row: Mapping[str, Any]) -> SiteKey:
    return (
        str(row["dataset"]),
        str(row["chrom"]),
        int(row["pos"]),
        str(row["ref"]),
        str(row["alt"]),
    )


def key_text(key: SiteKey) -> str:
    return f"{key[0]}|{key[1]}|{key[2]}|{key[3]}|{key[4]}"


def biological_key_text(row: Mapping[str, Any]) -> str:
    return (
        f"{row['biological_id']}|{row['chrom']}|{row['pos']}|"
        f"{row['ref']}|{row['alt']}"
    )


def load_module(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise SingletonAuditError(f"Unable to load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_m2_module(module: ModuleType) -> dict[str, Any]:
    expected = {
        "SCHEMA_NAME": "intersubmod.independent_m2_gate_recount",
        "SCHEMA_VERSION": "2.0.0",
        "ASSIGNMENT_SCHEMA": ASSIGNMENT_SCHEMA,
        "GATE_CONTRACT": EXPECTED_M2_CONTRACT,
        "P_MAX": 0.05,
        "PERMUTATIONS": 499,
        "P_FLOOR": 0.002,
        "TARGET_POWER": 0.80,
        "POWER_ALPHA": 0.05,
        "MIN_GROUP_N": 5,
        "MAX_GROUPS": 10,
        "CATEGORY_LEVEL_CEILINGS": EXPECTED_M2_CATEGORY_LEVEL_CEILINGS,
        "AXIS_SPECS": EXPECTED_M2_AXIS_SPECS,
        "INDETERMINATE_STATUSES": frozenset(
            {"HIGH_EFFECT_P_NOT_PASS", "LOW_POWER", "MISSING"}
        ),
    }
    observed = {
        name: getattr(module, name, None)
        for name in expected
    }
    mismatches = {
        name: {"observed": observed[name], "expected": value}
        for name, value in expected.items()
        if observed[name] != value
    }
    if mismatches:
        raise SingletonAuditError(f"Independent M2 module contract drift: {mismatches}")
    required_functions = (
        "assignment_key",
        "classify_site",
        "hp_family",
        "parse_cluster_sizes",
    )
    missing_functions = [
        name for name in required_functions if not callable(getattr(module, name, None))
    ]
    if missing_functions:
        raise SingletonAuditError(
            f"Independent M2 module lacks functions: {missing_functions}"
        )
    return {
        "gate_contract": observed["GATE_CONTRACT"],
        "axis_specs": [list(value) for value in observed["AXIS_SPECS"]],
        "category_level_ceilings": observed["CATEGORY_LEVEL_CEILINGS"],
        "minimum_group_n": observed["MIN_GROUP_N"],
        "maximum_groups": observed["MAX_GROUPS"],
        "permutations": observed["PERMUTATIONS"],
    }


def expected_dataset_contracts(
    manifest: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    samples = manifest.get("samples")
    if not isinstance(samples, list) or len(samples) != len(EXPECTED_DATASETS):
        raise SingletonAuditError("Input manifest sample list is not exactly seven rows")
    contracts: dict[str, dict[str, Any]] = {}
    for sample in samples:
        if not isinstance(sample, dict):
            raise SingletonAuditError("Input manifest sample entry is not an object")
        dataset = str(sample.get("sample", ""))
        if dataset in contracts:
            raise SingletonAuditError(f"Duplicate input-manifest sample: {dataset}")
        checks = sample.get("checks")
        if not isinstance(checks, dict) or not checks or not all(
            value is True for value in checks.values()
        ):
            raise SingletonAuditError(f"Input-manifest sample checks failed: {dataset}")
        counts = sample.get("counts")
        if not isinstance(counts, dict):
            raise SingletonAuditError(f"Input-manifest counts missing: {dataset}")
        contracts[dataset] = {
            "biological_id": str(sample.get("biological_id", "")),
            "counts": counts,
        }
    if tuple(sorted(contracts)) != EXPECTED_DATASETS:
        raise SingletonAuditError(
            f"Input-manifest datasets drift: {sorted(contracts)}"
        )
    return contracts


def validate_source_chain(
    *,
    paths: Mapping[str, Path],
    identities: Mapping[str, Mapping[str, Any]],
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, dict[str, Any]],
]:
    manifest = load_json(paths["input_manifest"], label="input manifest")
    run_manifest = load_json(
        paths["screen_run_manifest"], label="screen run manifest"
    )
    tree_audit = load_json(
        paths["tree_contract_audit"], label="tree contract audit"
    )
    primary_audit = load_json(
        paths["primary_artifact_audit"], label="primary artifact audit"
    )
    m2_receipt = load_json(
        paths["independent_m2_receipt"], label="independent M2 receipt"
    )
    require_schema(
        manifest,
        label="input manifest",
        schema_name="intersubmod.all_ssnv_focal_alt_input_manifest",
        schema_version="1.0.0",
    )
    require_schema(
        run_manifest,
        label="screen run manifest",
        schema_name="intersubmod.all_ssnv_focal_alt_multigroup_run_manifest",
        schema_version="1.2.0",
    )
    require_schema(
        tree_audit,
        label="tree contract audit",
        schema_name="intersubmod.latest_tree_input_contract_audit",
        schema_version="1.0.0",
    )
    require_schema(
        primary_audit,
        label="primary artifact audit",
        schema_name="intersubmod.stable_primary_artifact_audit",
        schema_version="2.0.0",
    )
    require_schema(
        m2_receipt,
        label="independent M2 receipt",
        schema_name="intersubmod.independent_m2_gate_recount",
        schema_version="2.0.0",
    )
    if (
        manifest.get("pass") is not True
        or manifest.get("task_type") != EXPECTED_TASK_TYPE
        or manifest.get("scope") != EXPECTED_SCOPE
    ):
        raise SingletonAuditError("Input manifest release contract drift")
    if (
        run_manifest.get("pass") is not True
        or run_manifest.get("status") != "EXECUTION_PASS"
    ):
        raise SingletonAuditError("Screen run manifest is not EXECUTION_PASS")
    if tree_audit.get("pass") is not True:
        raise SingletonAuditError("Tree contract audit is not PASS")
    require_true_checks(
        tree_audit.get("top_level_checks"),
        label="tree contract audit",
        expected_keys=EXPECTED_TREE_CHECKS,
    )
    if (
        primary_audit.get("pass") is not True
        or primary_audit.get("task_type") != EXPECTED_TASK_TYPE
        or primary_audit.get("formal_task_b_release_eligible") is not True
    ):
        raise SingletonAuditError("Primary artifact audit is not release-eligible")
    source_authority = verify_source_authority(
        paths=paths,
        identities=identities,
        primary_audit=primary_audit,
    )
    if (
        m2_receipt.get("pass") is not True
        or m2_receipt.get("task_type") != EXPECTED_TASK_TYPE
        or m2_receipt.get("contract") != EXPECTED_M2_CONTRACT
    ):
        raise SingletonAuditError("Independent M2 receipt contract drift")
    require_true_checks(
        m2_receipt.get("checks"),
        label="independent M2 receipt",
        expected_keys=EXPECTED_M2_CHECKS,
    )
    if m2_receipt.get("counts") != EXPECTED_M2_GLOBAL_COUNTS:
        raise SingletonAuditError("Independent M2 receipt count drift")
    if m2_receipt.get("logic_independence") != EXPECTED_M2_LOGIC_INDEPENDENCE:
        raise SingletonAuditError("Independent M2 receipt lost logic independence")

    manifest_totals = manifest.get("totals")
    expected_totals = {
        "all_ssnv": EXPECTED_COUNTS["all_dataset_sites"],
        **EXPECTED_GLOBAL_TRUTH_COUNTS,
        "ledger_branches": EXPECTED_BRANCH_COUNTS,
    }
    if (
        manifest_totals != expected_totals
        or tree_audit.get("totals") != expected_totals
        or tree_audit.get("scope") != EXPECTED_SCOPE
    ):
        raise SingletonAuditError("Manifest/tree totals or scope drift")

    require_identity(
        paths["input_manifest"],
        run_manifest.get("input_manifest"),
        label="run-manifest input manifest",
    )
    require_identity(
        paths["input_manifest"],
        tree_audit.get("input_manifest"),
        label="tree-audit input manifest",
        require_size=False,
    )
    for role in ("site_results", "stable_assignments"):
        require_identity(
            paths[role],
            run_manifest.get("outputs", {}).get(role),
            label=f"screen-run {role}",
        )
        require_identity(
            paths[role],
            primary_audit.get("inputs", {}).get(role),
            label=f"primary-audit {role}",
        )
        require_identity(
            paths[role],
            m2_receipt.get("inputs", {}).get(role),
            label=f"independent-M2 {role}",
        )
    require_identity(
        paths["claim_contract"],
        m2_receipt.get("inputs", {}).get("claim_contract"),
        label="independent-M2 claim contract",
    )
    require_identity(
        paths["claim_contract"],
        primary_audit.get("source_authority", {}).get("claim_contract"),
        label="source-authority claim contract",
    )
    require_identity(
        paths["independent_m2_auditor"],
        m2_receipt.get("code", {}).get("independent_recount"),
        label="independent-M2 source",
    )
    if identities["claim_contract"]["sha256"] != EXPECTED_CLAIM_SHA256:
        raise SingletonAuditError("Claim-contract v5 SHA-256 drift")
    require_read_only(
        identities,
        (
            "claim_contract",
            "independent_m2_auditor",
            "independent_m2_receipt",
            "primary_artifact_audit",
        ),
    )

    run_counts = run_manifest.get("counts")
    if not isinstance(run_counts, dict) or {
        "expected_sites": run_counts.get("expected_sites"),
        "processed_sites": run_counts.get("processed_sites"),
        "stable_assignment_records": run_counts.get("stable_assignment_records"),
    } != {
        "expected_sites": EXPECTED_COUNTS["all_dataset_sites"],
        "processed_sites": EXPECTED_COUNTS["all_dataset_sites"],
        "stable_assignment_records": EXPECTED_COUNTS["m1_stable_all"],
    }:
        raise SingletonAuditError("Screen run count drift")
    if primary_audit.get("counts") != {
        "stable_sites": EXPECTED_COUNTS["m1_stable_all"],
        "assignment_records": EXPECTED_COUNTS["m1_stable_all"],
        "primary_artifacts_expected": 308_526,
        "primary_artifacts_verified": 308_526,
    }:
        raise SingletonAuditError("Primary artifact count drift")
    primary_verification = primary_audit.get("verification")
    if (
        not isinstance(primary_verification, dict)
        or primary_verification.get("stable_site_assignment_key_sets_exact") is not True
        or primary_verification.get("path_size_sha256_verified") is not True
        or primary_verification.get("errors") != 0
        or primary_verification.get("artifact_roles_exact")
        != ["distance_matrix", "methylation_matrix", "reads"]
    ):
        raise SingletonAuditError("Primary artifact verification drift")

    contracts = run_manifest.get("contracts")
    if (
        not isinstance(contracts, dict)
        or contracts.get("m1_stability_gate_contract")
        != "coarse_ng>=2 AND not unstable AND modal_assignment_ari_min>=0.8"
        or contracts.get("strict_methyl_partition_robustness_status")
        != "NOT_EVALUATED_AT_M1_SCREEN"
        or contracts.get("strict_confirm_status_legacy_alias") != "NOT_RUN"
        or contracts.get("strict_confirm_candidate_is_formal_r1_claim") is not False
        or contracts.get("latest_hp_ps_terminal_join", {}).get("pass") is not True
        or contracts.get("latest_hp_ps_terminal_join", {}).get("all_sites_pass")
        is not True
    ):
        raise SingletonAuditError("Screen run contract drift")

    claim_text = paths["claim_contract"].read_text(encoding="utf-8")
    missing_fragments = [
        fragment for fragment in REQUIRED_CLAIM_FRAGMENTS if fragment not in claim_text
    ]
    if missing_fragments:
        raise SingletonAuditError(
            f"Claim contract lacks required fragments: {missing_fragments}"
        )
    dataset_contracts = expected_dataset_contracts(manifest)
    completed = run_manifest.get("completed_dataset_runs")
    if not isinstance(completed, dict) or set(completed) != set(dataset_contracts):
        raise SingletonAuditError("Screen completed-dataset set drift")
    for dataset, contract in dataset_contracts.items():
        validation = completed[dataset].get("validation", {})
        expected_n = contract["counts"]["all_ssnv"]
        if any(
            validation.get(field) != expected_n
            for field in ("expected_vcf_sites", "reads_files", "bernoulli_matrix_files", "methylation_files")
        ):
            raise SingletonAuditError(
                f"Screen completed-dataset count drift: {dataset}"
            )
    return (
        manifest,
        run_manifest,
        tree_audit,
        primary_audit,
        m2_receipt,
        source_authority,
        dataset_contracts,
    )


def m1_recomputed(row: Mapping[str, Any]) -> bool:
    if not parse_bool(row["m1_evaluable"], field="m1_evaluable"):
        return False
    try:
        coarse_ng = int(row["coarse_ng"])
    except (TypeError, ValueError) as error:
        raise SingletonAuditError(
            f"Invalid coarse_ng for evaluable site {screen_key(row)}"
        ) from error
    unstable = parse_bool(row["unstable"], field="unstable")
    ari_min = parse_optional_float(row["modal_assignment_ari_min"])
    return coarse_ng >= 2 and not unstable and ari_min is not None and ari_min >= 0.8


def _neighbor_key(item: Mapping[str, Any] | None) -> str:
    return "" if item is None else key_text(item["key"])


def recompute_components(
    groups: Mapping[tuple[str, str], Sequence[Mapping[str, Any]]],
) -> tuple[dict[SiteKey, dict[str, Any]], dict[str, int]]:
    singleton_details: dict[SiteKey, dict[str, Any]] = {}
    counters: Counter[str] = Counter(
        {
            "component_rows": 0,
            "component_identity_mismatch": 0,
            "singleton_branch_mismatch": 0,
            "recomputed_singletons": 0,
            "singleton_positional_contract_failure": 0,
        }
    )
    finite_singleton_gaps: list[int] = []
    for (dataset, chrom), unsorted_items in sorted(groups.items()):
        items = sorted(
            unsorted_items,
            key=lambda item: (
                int(item["pos"]),
                str(item["key"][3]),
                str(item["key"][4]),
            ),
        )
        components: list[list[Mapping[str, Any]]] = []
        active: list[Mapping[str, Any]] = []
        previous_pos: int | None = None
        for item in items:
            pos = int(item["pos"])
            if active and previous_pos is not None and pos - previous_pos > POSITIONAL_GAP_BP:
                components.append(active)
                active = []
            active.append(item)
            previous_pos = pos
        if active:
            components.append(active)

        index_by_key = {item["key"]: index for index, item in enumerate(items)}
        for component in components:
            first_pos = int(component[0]["pos"])
            last_pos = int(component[-1]["pos"])
            expected_id = f"{chrom}:{first_pos}-{last_pos}"
            expected_size = len(component)
            for item in component:
                counters["component_rows"] += 1
                if (
                    str(item["component_id"]) != expected_id
                    or int(item["component_size"]) != expected_size
                ):
                    counters["component_identity_mismatch"] += 1
                branch = str(item["branch"])
                if (expected_size == 1) != (branch == "positional_singleton"):
                    counters["singleton_branch_mismatch"] += 1
                if expected_size != 1:
                    continue
                counters["recomputed_singletons"] += 1
                index = index_by_key[item["key"]]
                left = items[index - 1] if index > 0 else None
                right = items[index + 1] if index + 1 < len(items) else None
                left_gap = (
                    int(item["pos"]) - int(left["pos"]) if left is not None else None
                )
                right_gap = (
                    int(right["pos"]) - int(item["pos"]) if right is not None else None
                )
                finite = [gap for gap in (left_gap, right_gap) if gap is not None]
                nearest_gap = min(finite) if finite else None
                if nearest_gap is not None:
                    finite_singleton_gaps.append(nearest_gap)
                positional_pass = (
                    str(item["component_id"]) == expected_id
                    and int(item["component_size"]) == 1
                    and branch == "positional_singleton"
                    and not bool(item["selected_for_read_census"])
                    and (left_gap is None or left_gap > POSITIONAL_GAP_BP)
                    and (right_gap is None or right_gap > POSITIONAL_GAP_BP)
                )
                if not positional_pass:
                    counters["singleton_positional_contract_failure"] += 1
                singleton_details[item["key"]] = {
                    "dataset": dataset,
                    "chrom": chrom,
                    "left_neighbor_key": _neighbor_key(left),
                    "right_neighbor_key": _neighbor_key(right),
                    "left_gap_bp": left_gap,
                    "right_gap_bp": right_gap,
                    "nearest_gap_bp": nearest_gap,
                    "recomputed_component_id": expected_id,
                    "recomputed_component_size": 1,
                    "positional_contract_pass": positional_pass,
                }
    counters["minimum_finite_singleton_nearest_gap_bp"] = (
        min(finite_singleton_gaps) if finite_singleton_gaps else -1
    )
    return singleton_details, dict(counters)


def load_screen_index(
    path: Path,
    dataset_contracts: Mapping[str, Mapping[str, Any]],
) -> tuple[
    dict[tuple[str, str], list[dict[str, Any]]],
    set[SiteKey],
    dict[str, Any],
]:
    groups: defaultdict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    all_keys: set[SiteKey] = set()
    stable_keys: set[SiteKey] = set()
    branch_counts: Counter[str] = Counter()
    truth_counts: Counter[str] = Counter()
    dataset_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    dataset_biological_ids: dict[str, str] = {}
    all_rows = 0
    duplicate_keys = 0
    strict_contract_mismatches: Counter[str] = Counter()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "dataset",
            "sample",
            "biological_id",
            "truth_label",
            "chrom",
            "pos",
            "ref",
            "alt",
            "ssnv_branch",
            "component_id",
            "component_size",
            "ledger_selected_for_read_census",
            "stable_null_multigroup",
            "filter",
            "latest_tag_join_status",
            "m1_stability_gate_contract",
            "strict_methyl_partition_robustness_status",
            "strict_methyl_partition_robustness_not_evaluable_reason",
            "strict_confirm_status",
            "strict_confirm_candidate",
            "strict_confirm_candidate_is_formal_r1_claim",
            "strict_confirm_reason",
            "m1_evaluable",
            "m1_not_evaluable_reason",
            "analysis_status",
            "coarse_ng",
            "unstable",
            "modal_assignment_ari_min",
        }
        missing = sorted(required.difference(reader.fieldnames or ()))
        if missing:
            raise SingletonAuditError(f"Site TSV lacks required fields: {missing}")
        for row in reader:
            all_rows += 1
            key = screen_key(row)
            dataset, chrom, pos, ref, alt = key
            if dataset not in dataset_contracts:
                raise SingletonAuditError(f"Unknown dataset at row {all_rows}: {dataset}")
            if chrom not in EXPECTED_CHROMS:
                raise SingletonAuditError(f"Invalid chromosome at {key}")
            if pos <= 0:
                raise SingletonAuditError(f"Invalid position at {key}")
            if ref not in EXPECTED_ALLELES or alt not in EXPECTED_ALLELES or ref == alt:
                raise SingletonAuditError(f"Invalid biallelic SNV alleles at {key}")
            if key in all_keys:
                duplicate_keys += 1
            all_keys.add(key)
            biological_id = str(row["biological_id"])
            expected_biological_id = str(
                dataset_contracts[dataset]["biological_id"]
            )
            observed_biological = dataset_biological_ids.setdefault(
                dataset, biological_id
            )
            if (
                observed_biological != biological_id
                or biological_id != expected_biological_id
                or str(row["sample"]) != dataset
            ):
                raise SingletonAuditError(f"Dataset identity drift at {key}")
            branch = str(row["ssnv_branch"])
            truth = str(row["truth_label"])
            if branch not in EXPECTED_BRANCHES:
                raise SingletonAuditError(f"Unknown sSNV branch at {key}: {branch}")
            if truth not in EXPECTED_TRUTH_LABELS:
                raise SingletonAuditError(f"Unknown truth label at {key}: {truth}")
            if str(row["filter"]) != "PASS":
                raise SingletonAuditError(f"Non-PASS focal site at {key}")
            if str(row["latest_tag_join_status"]) != "PASS":
                raise SingletonAuditError(f"Latest HP/PS join is not PASS at {key}")
            if (
                str(row["m1_stability_gate_contract"])
                != "coarse_ng>=2 AND not unstable AND modal_assignment_ari_min>=0.8"
            ):
                raise SingletonAuditError(f"M1 contract drift at {key}")
            branch_counts[branch] += 1
            truth_counts[truth] += 1
            dataset_counts[dataset]["all_ssnv"] += 1
            dataset_counts[dataset][f"truth_{truth.lower()}"] += 1
            dataset_counts[dataset][f"branch_{branch}"] += 1
            stable = parse_bool(
                row["stable_null_multigroup"], field="stable_null_multigroup"
            )
            m1_evaluable = parse_bool(
                row["m1_evaluable"], field="m1_evaluable"
            )
            analysis_status = str(row["analysis_status"])
            if analysis_status not in {
                "evaluable",
                "incomplete_distance_below_minimum",
                "insufficient_alt_reads",
            }:
                raise SingletonAuditError(
                    f"Unknown M1 analysis status at {key}: {analysis_status}"
                )
            try:
                coarse_ng = int(row["coarse_ng"])
            except (TypeError, ValueError) as error:
                raise SingletonAuditError(
                    f"Invalid coarse_ng at {key}: {row['coarse_ng']!r}"
                ) from error
            if coarse_ng < 0:
                raise SingletonAuditError(f"Negative coarse_ng at {key}")
            ari_min = parse_optional_float(row["modal_assignment_ari_min"])
            unstable_value = (
                None
                if row["unstable"] in {None, ""}
                else parse_bool(row["unstable"], field="unstable")
            )
            m1_reason = str(row["m1_not_evaluable_reason"])
            if m1_evaluable:
                if (
                    analysis_status != "evaluable"
                    or unstable_value is None
                    or ari_min is None
                    or m1_reason
                ):
                    raise SingletonAuditError(
                        f"Evaluable M1 field contract drift at {key}"
                    )
            else:
                if (
                    analysis_status == "evaluable"
                    or stable
                    or m1_reason
                    not in {
                        "INCOMPLETE_DISTANCE_BELOW_MINIMUM",
                        "INSUFFICIENT_FOCAL_ALT_READS",
                        "INSUFFICIENT_MATRIX_JOINED_FOCAL_ALT_READS",
                    }
                ):
                    raise SingletonAuditError(
                        f"Non-evaluable M1 field contract drift at {key}"
                    )
            if stable != m1_recomputed(row):
                raise SingletonAuditError(f"Global M1 recomputation drift at {key}")
            strict_candidate = parse_bool(
                row["strict_confirm_candidate"], field="strict_confirm_candidate"
            )
            formal_r1 = parse_bool(
                row["strict_confirm_candidate_is_formal_r1_claim"],
                field="strict_confirm_candidate_is_formal_r1_claim",
            )
            strict_reason = str(row["strict_confirm_reason"])
            strict_robustness_reason = str(
                row["strict_methyl_partition_robustness_not_evaluable_reason"]
            )
            strict_expectations = {
                "strict_methyl_partition_robustness_status": (
                    str(row["strict_methyl_partition_robustness_status"])
                    == "NOT_EVALUATED_AT_M1_SCREEN"
                ),
                "strict_confirm_status": str(row["strict_confirm_status"])
                == "NOT_RUN",
                "strict_candidate_is_m1_alias": strict_candidate == stable,
                "strict_formal_r1_false": formal_r1 is False,
                "strict_reason_exact": strict_reason == STRICT_NOT_EVALUATED_REASON,
                "strict_robustness_reason_exact": (
                    strict_robustness_reason == STRICT_NOT_EVALUATED_REASON
                ),
                "strict_reasons_equal": strict_reason == strict_robustness_reason,
            }
            for name, passed in strict_expectations.items():
                if not passed:
                    strict_contract_mismatches[name] += 1
            if stable:
                stable_keys.add(key)
            groups[(dataset, key[1])].append(
                {
                    "key": key,
                    "pos": key[2],
                    "component_id": row["component_id"],
                    "component_size": int(row["component_size"]),
                    "branch": branch,
                    "selected_for_read_census": parse_bool(
                        row["ledger_selected_for_read_census"],
                        field="ledger_selected_for_read_census",
                    ),
                }
            )
    metadata = {
        "all_rows": all_rows,
        "unique_keys": len(all_keys),
        "duplicate_keys": duplicate_keys,
        "stable_keys": len(stable_keys),
        "branch_counts": {
            branch: int(branch_counts.get(branch, 0))
            for branch in sorted(EXPECTED_BRANCHES)
        },
        "truth_counts": {
            truth: int(truth_counts.get(truth, 0))
            for truth in EXPECTED_TRUTH_LABELS
        },
        "dataset_counts": {
            dataset: {
                field: int(dataset_counts[dataset].get(field, 0))
                for field in (
                    "all_ssnv",
                    "branch_max_snv_excluded",
                    "branch_positional_singleton",
                    "branch_retained",
                    "truth_fp",
                    "truth_tp",
                    "truth_unassessed",
                )
            }
            for dataset in sorted(dataset_contracts)
        },
        "dataset_biological_ids": dict(sorted(dataset_biological_ids.items())),
        "datasets": sorted(dataset_biological_ids),
        "strict_contract_mismatches": dict(
            sorted(strict_contract_mismatches.items())
        ),
    }
    return dict(groups), stable_keys, metadata


def load_assignment_proofs(
    path: Path, m2_module: ModuleType
) -> tuple[dict[SiteKey, dict[str, Any]], dict[str, int]]:
    proofs: dict[SiteKey, dict[str, Any]] = {}
    counters: Counter[str] = Counter()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            record = json.loads(line)
            if (
                record.get("schema_name") != ASSIGNMENT_SCHEMA
                or record.get("schema_version") != "1.0.0"
            ):
                raise SingletonAuditError(
                    f"Unexpected assignment schema at line {line_number}"
                )
            key = m2_module.assignment_key(record)
            if key in proofs:
                raise SingletonAuditError(f"Duplicate assignment key: {key}")
            core_reads = record.get("core_reads")
            if not isinstance(core_reads, list) or not core_reads:
                raise SingletonAuditError(f"Assignment lacks core reads: {key}")
            read_ids = [str(read.get("read_id")) for read in core_reads]
            read_names = [str(read.get("read_name")) for read in core_reads]
            if len(set(read_ids)) != len(read_ids) or len(set(read_names)) != len(
                read_names
            ):
                raise SingletonAuditError(f"Duplicate core read identity: {key}")
            if [str(value) for value in record.get("read_ids", [])] != read_ids:
                raise SingletonAuditError(f"Assignment read_ids drift: {key}")
            if [str(value) for value in record.get("read_names", [])] != read_names:
                raise SingletonAuditError(f"Assignment read_names drift: {key}")
            for read in core_reads:
                if any(
                    field not in read or read[field] is None
                    for field in ("label", "latest_hp", "strand")
                ):
                    raise SingletonAuditError(
                        f"Assignment core read lacks M2 proof: {key}"
                    )
            group_sizes = Counter(str(read["label"]) for read in core_reads)
            proofs[key] = {
                "group_sizes": dict(group_sizes),
                "category_levels": {
                    "hp_exact": len(
                        {str(read["latest_hp"]) for read in core_reads}
                    ),
                    "hp_family": len(
                        {m2_module.hp_family(read["latest_hp"]) for read in core_reads}
                    ),
                    "strand": len({str(read["strand"]) for read in core_reads}),
                },
                "core_read_n": len(core_reads),
                "min_group_n": min(group_sizes.values()),
            }
            counters["assignment_rows"] += 1
            counters["core_reads"] += len(core_reads)
    return proofs, dict(counters)


def dataset_role(dataset: str) -> tuple[str, int]:
    if dataset == "HCC1395":
        return "technical_pair_primary", 1
    if dataset == "HCC1395_DORADO":
        return "technical_replicate", 0
    return "biological_sample_primary", 1


def dense_breakdown(
    source: Mapping[str, Counter[str]], strata: Iterable[str]
) -> dict[str, dict[str, int]]:
    metrics = (
        "sites",
        "m1_evaluable",
        "m1_flagged",
        "m2_pass",
        "m2_fail",
        "m2_not_evaluable",
        "m2_not_run",
    )
    return {
        stratum: {metric: int(source.get(stratum, Counter()).get(metric, 0)) for metric in metrics}
        for stratum in strata
    }


def classify_fail_subtype(row: Mapping[str, Any]) -> str:
    if parse_bool(row["hp_axis_confound"], field="hp_axis_confound"):
        return "HP_AXIS_CONFOUND"
    if parse_bool(row["technical_axis_confound"], field="technical_axis_confound"):
        return "TECHNICAL_AXIS_CONFOUND"
    return "NOT_PHASE_ANCHORED_ROBUST"


def m2_reason_codes(
    row: Mapping[str, Any], bucket: str, axes: Mapping[str, Mapping[str, Any]]
) -> list[str]:
    if bucket == "eligible":
        return ["ALL_MEASURED_AXES_DETERMINATE_NO_ALIGNED_CONFOUND"]
    if bucket == "evaluable_ineligible":
        return [classify_fail_subtype(row)]
    if bucket == "not_evaluable_group_count_gt10":
        return ["GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM"]
    reasons = [
        f"{prefix}:{axis['status']}"
        for prefix, axis in sorted(axes.items())
        if axis["status"] in {"HIGH_EFFECT_P_NOT_PASS", "LOW_POWER", "MISSING"}
    ]
    return reasons or ["AXIS_INDETERMINATE"]


SITE_OUTPUT_FIELDS = (
    "dataset",
    "sample",
    "biological_id",
    "dataset_role",
    "biological_n_contribution",
    "chrom",
    "pos",
    "ref",
    "alt",
    "site_key",
    "biological_site_key",
    "truth_label",
    "ssnv_branch",
    "component_id",
    "component_size",
    "left_neighbor_key",
    "right_neighbor_key",
    "left_gap_bp",
    "right_gap_bp",
    "nearest_gap_bp",
    "recomputed_component_id",
    "recomputed_component_size",
    "positional_contract_pass",
    "m1_evaluable",
    "m1_status",
    "m1_reason",
    "analysis_status",
    "n_alt_matrix",
    "n_alt_after_peel",
    "n_alt_peeled",
    "coarse_ng",
    "unstable",
    "modal_assignment_ari_min",
    "stable_null_multigroup",
    "cluster_sizes",
    "m2_status",
    "m2_bucket",
    "m2_reason_codes",
    "methyl_group_count",
    "core_read_n",
    "min_group_n",
    "m2_aligned_axes",
    "m2_indeterminate_axes",
    "m2_constant_axes",
    "m2_low_power_axes",
    "m2_axis_details_json",
    "strict_confirm_candidate",
    "strict_confirm_candidate_deprecated_alias",
    "strict_confirm_candidate_is_formal_r1_claim",
    "strict_confirm_status",
    "strict_confirm_reason",
    "g1_status",
    "g1_reason",
    "g2_status",
    "g2_reason",
    "r1_status",
    "r1_reason",
    "final_link_status",
)


def build_audit(
    *,
    site_results: Path,
    stable_assignments: Path,
    input_manifest: Path,
    screen_run_manifest: Path,
    tree_contract_audit: Path,
    primary_artifact_audit: Path,
    source_authority: Path,
    source_authority_approval: Path,
    source_authority_signature: Path,
    source_authority_public_key: Path,
    independent_m2_auditor: Path,
    independent_m2_receipt: Path,
    claim_contract: Path,
    output_dir: Path,
) -> dict[str, Any]:
    output_dir = output_dir.expanduser().resolve()
    if os.path.lexists(output_dir):
        raise FileExistsError(f"Refusing to overwrite output directory: {output_dir}")
    input_paths = {
        "site_results": site_results.resolve(strict=True),
        "stable_assignments": stable_assignments.resolve(strict=True),
        "input_manifest": input_manifest.resolve(strict=True),
        "screen_run_manifest": screen_run_manifest.resolve(strict=True),
        "tree_contract_audit": tree_contract_audit.resolve(strict=True),
        "primary_artifact_audit": primary_artifact_audit.resolve(strict=True),
        "source_authority": source_authority.resolve(strict=True),
        "source_authority_approval": source_authority_approval.resolve(strict=True),
        "source_authority_signature": source_authority_signature.resolve(
            strict=True
        ),
        "source_authority_public_key": source_authority_public_key.resolve(
            strict=True
        ),
        "independent_m2_auditor": independent_m2_auditor.resolve(strict=True),
        "independent_m2_receipt": independent_m2_receipt.resolve(strict=True),
        "claim_contract": claim_contract.resolve(strict=True),
    }
    input_identity_before = {
        name: artifact(path) for name, path in input_paths.items()
    }
    supplemental_source_before = artifact(Path(__file__).resolve())
    (
        manifest,
        run_manifest,
        tree_audit,
        primary_audit,
        m2_receipt,
        source_authority_validation,
        dataset_contracts,
    ) = validate_source_chain(
        paths=input_paths,
        identities=input_identity_before,
    )

    m2_module = load_module(
        input_paths["independent_m2_auditor"],
        "positional_singleton_independent_m2_replay",
    )
    m2_contract_validation = validate_m2_module(m2_module)
    groups, stable_keys, screen_metadata = load_screen_index(
        input_paths["site_results"], dataset_contracts
    )
    singleton_details, positional_counts = recompute_components(groups)
    proofs, assignment_counts = load_assignment_proofs(
        input_paths["stable_assignments"], m2_module
    )
    assignment_keys = set(proofs)

    singleton_rows: list[dict[str, Any]] = []
    m2_pass_rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    m1_reason_counts: Counter[str] = Counter()
    m2_bucket_counts: Counter[str] = Counter()
    m2_fail_subtypes: Counter[str] = Counter()
    m2_axis_status_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    group_count_distribution: Counter[int] = Counter()
    by_dataset: defaultdict[str, Counter[str]] = defaultdict(Counter)
    by_truth: defaultdict[str, Counter[str]] = defaultdict(Counter)
    by_dataset_truth: defaultdict[str, Counter[str]] = defaultdict(Counter)
    alias_mismatch = 0
    formal_r1_true = 0
    strict_not_run_mismatch = 0
    m1_recompute_mismatch = 0
    singleton_stable_keys: set[SiteKey] = set()

    with open_text(input_paths["site_results"]) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["ssnv_branch"] != "positional_singleton":
                continue
            key = screen_key(row)
            positional = singleton_details.get(key)
            if positional is None:
                raise SingletonAuditError(f"Singleton lacks positional proof: {key}")
            counts["singleton_sites"] += 1
            truth = str(row["truth_label"])
            dataset = str(row["dataset"])
            strata = (by_dataset[dataset], by_truth[truth], by_dataset_truth[f"{dataset}|{truth}"])
            for target in strata:
                target["sites"] += 1

            m1_evaluable = parse_bool(row["m1_evaluable"], field="m1_evaluable")
            declared_m1 = parse_bool(
                row["stable_null_multigroup"], field="stable_null_multigroup"
            )
            if declared_m1 != m1_recomputed(row):
                m1_recompute_mismatch += 1
            if m1_evaluable:
                counts["m1_evaluable"] += 1
                for target in strata:
                    target["m1_evaluable"] += 1
            else:
                counts["m1_not_evaluable"] += 1
                m1_reason_counts[str(row["m1_not_evaluable_reason"])] += 1
            if declared_m1:
                counts["m1_flagged"] += 1
                singleton_stable_keys.add(key)
                for target in strata:
                    target["m1_flagged"] += 1

            strict_candidate = parse_bool(
                row["strict_confirm_candidate"], field="strict_confirm_candidate"
            )
            if strict_candidate != declared_m1:
                alias_mismatch += 1
            if parse_bool(
                row["strict_confirm_candidate_is_formal_r1_claim"],
                field="strict_confirm_candidate_is_formal_r1_claim",
            ):
                formal_r1_true += 1
            if row["strict_confirm_status"] != "NOT_RUN":
                strict_not_run_mismatch += 1

            m2_status = "NOT_RUN"
            m2_bucket = "not_run_m1_not_flagged"
            reasons = ["M1_NOT_FLAGGED"]
            axes: dict[str, Mapping[str, Any]] = {}
            group_count: int | str = ""
            core_read_n: int | str = ""
            min_group_n: int | str = ""
            aligned_axes: list[str] = []
            indeterminate_axes: list[str] = []
            constant_axes: list[str] = []
            low_power_axes: list[str] = []
            if declared_m1:
                proof = proofs.get(key)
                if proof is None:
                    raise SingletonAuditError(
                        f"Singleton M1 site lacks assignment proof: {key}"
                    )
                result = m2_module.classify_site(row, proof)
                m2_bucket = str(result["bucket"])
                axes = result["axes"]
                m2_bucket_counts[m2_bucket] += 1
                group_count = len(m2_module.parse_cluster_sizes(row["cluster_sizes"]))
                group_count_distribution[int(group_count)] += 1
                core_read_n = int(proof["core_read_n"])
                min_group_n = int(proof["min_group_n"])
                for prefix, axis in axes.items():
                    status = str(axis["status"])
                    m2_axis_status_counts[prefix][status] += 1
                aligned_axes = sorted(
                    prefix for prefix, axis in axes.items() if axis["aligned"]
                )
                indeterminate_axes = sorted(
                    prefix
                    for prefix, axis in axes.items()
                    if axis["status"]
                    in {"HIGH_EFFECT_P_NOT_PASS", "LOW_POWER", "MISSING"}
                )
                constant_axes = sorted(
                    prefix
                    for prefix, axis in axes.items()
                    if axis["status"] == "CONSTANT"
                )
                low_power_axes = sorted(
                    prefix
                    for prefix, axis in axes.items()
                    if axis["status"] == "LOW_POWER"
                )
                reasons = m2_reason_codes(row, m2_bucket, axes)
                if m2_bucket == "eligible":
                    m2_status = "PASS"
                    counts["m2_pass"] += 1
                    for target in strata:
                        target["m2_pass"] += 1
                elif m2_bucket == "evaluable_ineligible":
                    m2_status = "FAIL"
                    counts["m2_fail"] += 1
                    subtype = classify_fail_subtype(row)
                    m2_fail_subtypes[subtype] += 1
                    for target in strata:
                        target["m2_fail"] += 1
                else:
                    m2_status = "NOT_EVALUABLE"
                    counts["m2_not_evaluable"] += 1
                    for target in strata:
                        target["m2_not_evaluable"] += 1
            else:
                counts["m2_not_run"] += 1
                for target in strata:
                    target["m2_not_run"] += 1

            role, contribution = dataset_role(dataset)
            output_row = {
                "dataset": dataset,
                "sample": row["sample"],
                "biological_id": row["biological_id"],
                "dataset_role": role,
                "biological_n_contribution": contribution,
                "chrom": row["chrom"],
                "pos": row["pos"],
                "ref": row["ref"],
                "alt": row["alt"],
                "site_key": key_text(key),
                "biological_site_key": biological_key_text(row),
                "truth_label": truth,
                "ssnv_branch": row["ssnv_branch"],
                "component_id": row["component_id"],
                "component_size": row["component_size"],
                **positional,
                "m1_evaluable": str(m1_evaluable).lower(),
                "m1_status": "FLAGGED" if declared_m1 else "NOT_FLAGGED",
                "m1_reason": (
                    "OPERATIONAL_STABLE_MULTIGROUP"
                    if declared_m1
                    else row["m1_not_evaluable_reason"]
                    or "EVALUABLE_NOT_FLAGGED"
                ),
                "analysis_status": row["analysis_status"],
                "n_alt_matrix": row["n_alt_matrix"],
                "n_alt_after_peel": row["n_alt_after_peel"],
                "n_alt_peeled": row["n_alt_peeled"],
                "coarse_ng": row["coarse_ng"],
                "unstable": row["unstable"],
                "modal_assignment_ari_min": row["modal_assignment_ari_min"],
                "stable_null_multigroup": row["stable_null_multigroup"],
                "cluster_sizes": row["cluster_sizes"],
                "m2_status": m2_status,
                "m2_bucket": m2_bucket,
                "m2_reason_codes": json.dumps(reasons, separators=(",", ":")),
                "methyl_group_count": group_count,
                "core_read_n": core_read_n,
                "min_group_n": min_group_n,
                "m2_aligned_axes": ",".join(aligned_axes),
                "m2_indeterminate_axes": ",".join(indeterminate_axes),
                "m2_constant_axes": ",".join(constant_axes),
                "m2_low_power_axes": ",".join(low_power_axes),
                "m2_axis_details_json": json.dumps(
                    axes, ensure_ascii=True, sort_keys=True, separators=(",", ":")
                ),
                "strict_confirm_candidate": row["strict_confirm_candidate"],
                "strict_confirm_candidate_deprecated_alias": (
                    "M1_COMPATIBILITY_ALIAS_NOT_FORMAL_R1"
                ),
                "strict_confirm_candidate_is_formal_r1_claim": row[
                    "strict_confirm_candidate_is_formal_r1_claim"
                ],
                "strict_confirm_status": row["strict_confirm_status"],
                "strict_confirm_reason": row["strict_confirm_reason"],
                "g1_status": "NOT_RUN",
                "g1_reason": "NO_LOCAL_SSNV_PARTNER_POSITIONAL_SINGLETON",
                "g2_status": "NOT_RUN",
                "g2_reason": "G1_NOT_RUN_NO_LOCAL_SSNV_PARTNER",
                "r1_status": "NOT_RUN",
                "r1_reason": "G2_NOT_RUN_NO_MULTI_MARKER_BASE_CANDIDATE",
                "final_link_status": "PENDING_SIGNED_FINAL_DATASET_LINK",
            }
            singleton_rows.append(output_row)
            if m2_status == "PASS":
                m2_pass_rows.append(output_row)

    input_identity_after = {
        name: artifact(path) for name, path in input_paths.items()
    }
    supplemental_source_after = artifact(Path(__file__).resolve())
    expected_dataset_counts = {
        dataset: {
            "all_ssnv": int(contract["counts"]["all_ssnv"]),
            "branch_max_snv_excluded": int(
                contract["counts"]["ledger_branches"]["max_snv_excluded"]
            ),
            "branch_positional_singleton": int(
                contract["counts"]["ledger_branches"]["positional_singleton"]
            ),
            "branch_retained": int(
                contract["counts"]["ledger_branches"]["retained"]
            ),
            "truth_fp": int(contract["counts"]["truth_fp"]),
            "truth_tp": int(contract["counts"]["truth_tp"]),
            "truth_unassessed": int(contract["counts"]["truth_unassessed"]),
        }
        for dataset, contract in sorted(dataset_contracts.items())
    }
    checks = {
        "input_identity_unchanged": input_identity_before == input_identity_after,
        "supplemental_source_identity_unchanged": (
            supplemental_source_before == supplemental_source_after
        ),
        "source_chain_contracts_validated": True,
        "m2_replay_constants_exact": m2_contract_validation["gate_contract"]
        == EXPECTED_M2_CONTRACT,
        "all_rows_match_expected": screen_metadata["all_rows"]
        == EXPECTED_COUNTS["all_dataset_sites"],
        "all_primary_keys_unique": screen_metadata["duplicate_keys"] == 0
        and screen_metadata["unique_keys"] == screen_metadata["all_rows"],
        "datasets_exact": tuple(screen_metadata["datasets"]) == EXPECTED_DATASETS,
        "biological_sample_count_is_six": len(
            set(screen_metadata["dataset_biological_ids"].values())
        )
        == 6,
        "branch_counts_match_expected": screen_metadata["branch_counts"]
        == EXPECTED_BRANCH_COUNTS,
        "truth_counts_match_expected": screen_metadata["truth_counts"]
        == {
            "FP": manifest["totals"]["truth_fp"],
            "TP": manifest["totals"]["truth_tp"],
            "UNASSESSED": manifest["totals"]["truth_unassessed"],
        },
        "dataset_counts_match_bound_manifest": screen_metadata["dataset_counts"]
        == expected_dataset_counts,
        "all_screen_strict_r1_fields_exact": not screen_metadata[
            "strict_contract_mismatches"
        ],
        "component_identity_exact": positional_counts.get(
            "component_identity_mismatch", 0
        )
        == 0,
        "singleton_branch_exact": positional_counts.get(
            "singleton_branch_mismatch", 0
        )
        == 0,
        "singleton_positional_contract_exact": positional_counts.get(
            "singleton_positional_contract_failure", 0
        )
        == 0,
        "singleton_count_matches_expected": counts["singleton_sites"]
        == EXPECTED_COUNTS["positional_singleton"],
        "singleton_minimum_gap_matches_expected": positional_counts[
            "minimum_finite_singleton_nearest_gap_bp"
        ]
        == EXPECTED_COUNTS["singleton_min_finite_nearest_gap_bp"],
        "global_m1_stable_matches_expected": len(stable_keys)
        == EXPECTED_COUNTS["m1_stable_all"],
        "stable_assignment_keys_exact": assignment_keys == stable_keys,
        "singleton_assignment_intersection_exact": len(
            assignment_keys.intersection(singleton_details)
        )
        == EXPECTED_COUNTS["singleton_m1_flagged"],
        "m1_recompute_exact": m1_recompute_mismatch == 0,
        "singleton_m1_evaluable_matches_expected": counts["m1_evaluable"]
        == EXPECTED_COUNTS["singleton_m1_evaluable"],
        "singleton_m1_not_evaluable_matches_expected": counts["m1_not_evaluable"]
        == EXPECTED_COUNTS["singleton_m1_not_evaluable"],
        "singleton_m1_reason_census_matches_expected": dict(
            sorted(m1_reason_counts.items())
        )
        == EXPECTED_M1_NOT_EVALUABLE_REASONS,
        "singleton_m1_flagged_matches_expected": counts["m1_flagged"]
        == EXPECTED_COUNTS["singleton_m1_flagged"],
        "singleton_stable_key_set_exact": singleton_stable_keys
        == assignment_keys.intersection(singleton_details),
        "singleton_m2_pass_matches_expected": counts["m2_pass"]
        == EXPECTED_COUNTS["singleton_m2_pass"],
        "singleton_m2_fail_matches_expected": counts["m2_fail"]
        == EXPECTED_COUNTS["singleton_m2_fail"],
        "singleton_m2_axis_indeterminate_matches_expected": m2_bucket_counts[
            "not_evaluable_axis_indeterminate"
        ]
        == EXPECTED_COUNTS["singleton_m2_axis_indeterminate"],
        "singleton_m2_group_count_gt10_matches_expected": m2_bucket_counts[
            "not_evaluable_group_count_gt10"
        ]
        == EXPECTED_COUNTS["singleton_m2_group_count_gt10"],
        "singleton_m2_not_run_matches_expected": counts["m2_not_run"]
        == EXPECTED_COUNTS["singleton_m2_not_run"],
        "singleton_m2_conserved": counts["m1_flagged"]
        == counts["m2_pass"] + counts["m2_fail"] + counts["m2_not_evaluable"],
        "singleton_full_m2_status_conserved": counts["singleton_sites"]
        == counts["m2_pass"]
        + counts["m2_fail"]
        + counts["m2_not_evaluable"]
        + counts["m2_not_run"],
        "singleton_m2_fail_subtypes_match_expected": {
            key: m2_fail_subtypes.get(key, 0)
            for key in sorted(EXPECTED_M2_FAIL_SUBTYPES)
        }
        == EXPECTED_M2_FAIL_SUBTYPES,
        "legacy_strict_candidate_is_exact_m1_alias": alias_mismatch == 0,
        "legacy_alias_is_never_formal_r1": formal_r1_true == 0,
        "strict_r1_status_is_not_run": strict_not_run_mismatch == 0,
        "m2_pass_case_count_matches": len(m2_pass_rows) == counts["m2_pass"],
        "all_supplemental_downstream_claims_not_run": all(
            row["g1_status"] == row["g2_status"] == row["r1_status"] == "NOT_RUN"
            for row in singleton_rows
        ),
    }
    if not all(checks.values()):
        failed = sorted(key for key, passed in checks.items() if not passed)
        raise SingletonAuditError(f"Singleton audit checks failed: {failed}")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = output_dir.parent / (
        f".{output_dir.name}.staging.{os.getpid()}.{uuid.uuid4().hex}"
    )
    staging_dir.mkdir(mode=0o755)
    site_output = staging_dir / "positional_singleton_site_audit.tsv.gz"
    candidate_output = staging_dir / "positional_singleton_m2_pass_cases.tsv"
    final_site_output = output_dir / site_output.name
    final_candidate_output = output_dir / candidate_output.name
    with gzip.open(site_output, "xt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=SITE_OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(singleton_rows)
    candidate_fields = (
        "dataset",
        "biological_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "truth_label",
        "methyl_group_count",
        "core_read_n",
        "min_group_n",
        "m2_status",
        "m2_reason_codes",
        "m2_constant_axes",
        "m2_axis_details_json",
    )
    with candidate_output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=candidate_fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(m2_pass_rows)
    site_output.chmod(0o444)
    candidate_output.chmod(0o444)
    fsync_file(site_output)
    fsync_file(candidate_output)

    summary = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "B_comprehensive_validation",
        "scope": (
            f"all_{counts['singleton_sites']}_positional_singleton_dataset_sites_"
            f"across_{len(EXPECTED_DATASETS)}_datasets_chr1_22"
        ),
        "command": [str(Path(sys.executable).resolve()), *sys.argv],
        "inputs": input_identity_after,
        "code": {
            "supplemental_auditor": supplemental_source_after,
            "m2_replay_source": input_identity_after["independent_m2_auditor"],
        },
        "source_chain": {
            "main_task_b_source_authority": source_authority_validation,
            "screen_run_manifest_status": run_manifest["status"],
            "primary_artifact_audit_release_eligible": primary_audit[
                "formal_task_b_release_eligible"
            ],
            "primary_artifact_set_sha256": primary_audit["verification"][
                "artifact_set_sha256"
            ],
            "independent_m2_receipt_contract": m2_receipt["contract"],
            "tree_input_contract_pass": tree_audit["pass"],
            "claim_contract_sha256": input_identity_after["claim_contract"][
                "sha256"
            ],
        },
        "contracts": {
            "positional_singleton": (
                "same_dataset_chrom_adjacent_gap_gt_50000_component_size_1_v1"
            ),
            "m1": "coarse_ng>=2 AND not unstable AND modal_assignment_ari_min>=0.8",
            "m2": m2_contract_validation["gate_contract"],
            "m2_replay_parameters": m2_contract_validation,
            "statistical_unit": "(dataset,chrom,pos,ref,alt)",
            "claim_ceiling": "M2_read_level_residual_epigenetic_partition",
        },
        "screen_metadata": screen_metadata,
        "positional_recomputation": positional_counts,
        "assignment_recomputation": assignment_counts,
        "counts": dict(sorted(counts.items())),
        "m1_not_evaluable_reasons": dict(sorted(m1_reason_counts.items())),
        "m2_buckets": dict(sorted(m2_bucket_counts.items())),
        "m2_fail_subtypes": {
            key: m2_fail_subtypes.get(key, 0)
            for key in sorted(EXPECTED_M2_FAIL_SUBTYPES)
        },
        "m2_axis_status_census": {
            prefix: dict(sorted(values.items()))
            for prefix, values in sorted(m2_axis_status_counts.items())
        },
        "methyl_group_count_distribution": {
            str(group_count): count
            for group_count, count in sorted(group_count_distribution.items())
        },
        "status_census": {
            "m1": {
                "FLAGGED": counts["m1_flagged"],
                "NOT_FLAGGED": counts["singleton_sites"] - counts["m1_flagged"],
            },
            "m2": {
                "PASS": counts["m2_pass"],
                "FAIL": counts["m2_fail"],
                "NOT_EVALUABLE": counts["m2_not_evaluable"],
                "NOT_RUN": counts["m2_not_run"],
            },
            "g1": {"NOT_RUN": counts["singleton_sites"]},
            "g2": {"NOT_RUN": counts["singleton_sites"]},
            "r1": {"NOT_RUN": counts["singleton_sites"]},
        },
        "breakdowns": {
            "dataset": dense_breakdown(by_dataset, EXPECTED_DATASETS),
            "truth": dense_breakdown(by_truth, EXPECTED_TRUTH_LABELS),
            "dataset_truth": dense_breakdown(
                by_dataset_truth,
                (
                    f"{dataset}|{truth}"
                    for dataset in EXPECTED_DATASETS
                    for truth in EXPECTED_TRUTH_LABELS
                ),
            ),
        },
        "rates": {
            "singleton_fraction_of_all_dataset_sites": rate_record(
                counts["singleton_sites"],
                screen_metadata["all_rows"],
                f"all {screen_metadata['all_rows']:,} LongPhase-S recalibrated "
                "FILTER=PASS autosomal biallelic sSNV dataset-sites",
            ),
            "m1_evaluable_fraction_of_singletons": rate_record(
                counts["m1_evaluable"],
                counts["singleton_sites"],
                f"all {counts['singleton_sites']:,} positional-singleton "
                "dataset-sites",
            ),
            "m1_flag_yield_all_singletons": rate_record(
                counts["m1_flagged"],
                counts["singleton_sites"],
                f"all {counts['singleton_sites']:,} positional-singleton "
                "dataset-sites",
            ),
            "m1_flag_yield_evaluable_singletons": rate_record(
                counts["m1_flagged"],
                counts["m1_evaluable"],
                f"{counts['m1_evaluable']:,} M1-evaluable positional-singleton "
                "dataset-sites",
            ),
            "m2_evaluable_fraction_of_singleton_m1": rate_record(
                counts["m2_pass"] + counts["m2_fail"],
                counts["m1_flagged"],
                f"{counts['m1_flagged']:,} M1-flagged positional-singleton "
                "dataset-sites",
            ),
            "m2_pass_fraction_of_singleton_m1": rate_record(
                counts["m2_pass"],
                counts["m1_flagged"],
                f"{counts['m1_flagged']:,} M1-flagged positional-singleton "
                "dataset-sites",
            ),
            "m2_pass_fraction_of_m2_evaluable_singletons": rate_record(
                counts["m2_pass"],
                counts["m2_pass"] + counts["m2_fail"],
                f"{counts['m2_pass'] + counts['m2_fail']:,} selected "
                "M2-evaluable positional-singleton dataset-sites; not "
                "population-extrapolatable",
            ),
            "m2_pass_fraction_of_all_singletons": rate_record(
                counts["m2_pass"],
                counts["singleton_sites"],
                f"all {counts['singleton_sites']:,} positional-singleton "
                "dataset-sites",
            ),
        },
        "legacy_alias_audit": {
            "strict_confirm_candidate_equals_m1_mismatches": alias_mismatch,
            "formal_r1_true_rows": formal_r1_true,
            "strict_status_not_run_mismatches": strict_not_run_mismatch,
            "interpretation": (
                "strict_confirm_candidate is a deprecated M1 compatibility alias "
                "and is never a formal R1 result in this source screen"
            ),
        },
        "outputs": {
            "site_audit": artifact(site_output, public_path=final_site_output),
            "m2_pass_cases": artifact(
                candidate_output, public_path=final_candidate_output
            ),
        },
        "checks": checks,
        "limitations": [
            "Positional singleton means no same-dataset/chrom sSNV within the "
            "50 kb component contract; it is not a general read-sharing graph "
            "degree-zero definition.",
            "M1 is an operational screen yield, not biological prevalence.",
            "M2 evaluates eight measured axes and does not exclude unmeasured "
            "technical, biological, purity, or CN confounding.",
            "The M2-evaluable subset is highly selected; 30/48 must not be "
            "extrapolated to all M1 or singleton sites.",
            "HCC1395 and HCC1395_DORADO are one biological sample and cannot be "
            "counted as independent replication.",
            "Without a local partner, G1/G2 and formal R1 are NOT_RUN; M1/M2 do "
            "not establish cellular subclones, clone number, linear ancestry, "
            "or lineage order.",
        ],
        "pass_semantics": (
            "complete_positional_singleton_position_m1_m2_and_legacy_alias_"
            "reconciliation_not_cellular_subclone_confirmation"
        ),
        "pass": True,
    }
    summary_output = staging_dir / "positional_singleton_audit_summary.json"
    final_summary_output = output_dir / summary_output.name
    with summary_output.open("x", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    summary_output.chmod(0o444)
    fsync_file(summary_output)
    success_output = staging_dir / "_SUCCESS.json"
    final_success_output = output_dir / success_output.name
    success_record = {
        "schema_name": "intersubmod.atomic_release_marker",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "summary": artifact(summary_output, public_path=final_summary_output),
        "site_audit_sha256": summary["outputs"]["site_audit"]["sha256"],
        "m2_pass_cases_sha256": summary["outputs"]["m2_pass_cases"]["sha256"],
        "pass": True,
    }
    with success_output.open("x", encoding="utf-8") as handle:
        json.dump(success_record, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    success_output.chmod(0o444)
    fsync_file(success_output)
    fsync_directory(staging_dir)
    staging_dir.chmod(0o555)
    rename_no_replace(staging_dir, output_dir)
    fsync_directory(output_dir.parent)
    return {
        "output_dir": str(output_dir),
        "summary": str(final_summary_output),
        "success_marker": str(final_success_output),
        "site_rows": len(singleton_rows),
        "m2_pass_cases": len(m2_pass_rows),
        "rates": summary["rates"],
        "atomic_release": True,
        "pass": True,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-results", type=Path, required=True)
    parser.add_argument("--stable-assignments", type=Path, required=True)
    parser.add_argument("--input-manifest", type=Path, required=True)
    parser.add_argument("--screen-run-manifest", type=Path, required=True)
    parser.add_argument("--tree-contract-audit", type=Path, required=True)
    parser.add_argument("--primary-artifact-audit", type=Path, required=True)
    parser.add_argument("--source-authority", type=Path, required=True)
    parser.add_argument("--source-authority-approval", type=Path, required=True)
    parser.add_argument("--source-authority-signature", type=Path, required=True)
    parser.add_argument("--source-authority-public-key", type=Path, required=True)
    parser.add_argument("--independent-m2-auditor", type=Path, required=True)
    parser.add_argument("--independent-m2-receipt", type=Path, required=True)
    parser.add_argument("--claim-contract", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = build_audit(
        site_results=args.site_results.expanduser(),
        stable_assignments=args.stable_assignments.expanduser(),
        input_manifest=args.input_manifest.expanduser(),
        screen_run_manifest=args.screen_run_manifest.expanduser(),
        tree_contract_audit=args.tree_contract_audit.expanduser(),
        primary_artifact_audit=args.primary_artifact_audit.expanduser(),
        source_authority=args.source_authority.expanduser(),
        source_authority_approval=args.source_authority_approval.expanduser(),
        source_authority_signature=args.source_authority_signature.expanduser(),
        source_authority_public_key=args.source_authority_public_key.expanduser(),
        independent_m2_auditor=args.independent_m2_auditor.expanduser(),
        independent_m2_receipt=args.independent_m2_receipt.expanduser(),
        claim_contract=args.claim_contract.expanduser(),
        output_dir=args.output_dir.expanduser().resolve(),
    )
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
