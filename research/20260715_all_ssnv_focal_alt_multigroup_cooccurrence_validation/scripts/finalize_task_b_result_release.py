#!/usr/bin/env python3
"""Create or verify the detached-signature release receipt for the final Task-B dataset."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import build_all_ssnv_final_report_dataset as BUILDER  # noqa: E402
import release_source_authority as SOURCE_AUTHORITY  # noqa: E402


SCHEMA_NAME = "intersubmod.task_b_final_dataset_release_receipt"
SCHEMA_VERSION = "1.0.0"
FINAL_DATASET_DIR = BUILDER.CANONICAL_FINAL_DATASET_DIR
BUILDER_RECEIPT = FINAL_DATASET_DIR / "run_receipt.json"
RELEASE_RECEIPT = (
    SCRIPT_DIR.parent / "results" / "task_b_final_dataset_release_receipt.v1.json"
)
RELEASE_SIGNATURE = Path(str(RELEASE_RECEIPT) + ".ed25519.sig")
RESULT_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
    "20260719_all_ssnv_result_v5_post_reboot_bootstrap/ed25519_public.pem"
)
RESULT_PUBLIC_KEY_SHA256 = (
    "5b7d5d026835ec6ec677bcd886bc16ac71444117dabc44a084ce3ede3a4db5a9"
)
RESULT_PRIVATE_KEY = RESULT_PUBLIC_KEY.with_name(
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
REPORT_SCHEMA_NAME = "intersubmod.task_b_final_report_release_receipt"
REPORT_RELEASE_RECEIPT = (
    SCRIPT_DIR.parent / "results" / "task_b_final_report_release_receipt.v1.json"
)
REPORT_RELEASE_SIGNATURE = Path(str(REPORT_RELEASE_RECEIPT) + ".ed25519.sig")
REPORT_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
    "20260719_all_ssnv_report_v5_post_reboot_bootstrap/ed25519_public.pem"
)
REPORT_PUBLIC_KEY_SHA256 = (
    "572e27167e1eea4c39ca53546ba33868b2874ec7fa3ed1682db821b2e50fa439"
)
REPORT_PRIVATE_KEY = REPORT_PUBLIC_KEY.with_name(
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
CANONICAL_PYTHON_CACHE_DIRNAME = (
    ".python_cache_m2v5_completion_v2_bound_bootstrap"
)
CANONICAL_PYTHON_CACHE_ROOT = WORKSPACE_ROOT / CANONICAL_PYTHON_CACHE_DIRNAME
REPORT_DIR = WORKSPACE_ROOT / "all_ssnv_final_report_v5_m2v5_source_attested"
REPORT_OUTPUTS = {
    "report_markdown": REPORT_DIR / "report.md",
    "canonical_report_artifact": REPORT_DIR / "artifact.json",
    "report_build_receipt": REPORT_DIR / "report_build_receipt.json",
    "portable_html": REPORT_DIR / "all_ssnv_focal_alt_multigroup_cooccurrence_report.standalone.html",
    "portable_delivery_receipt": REPORT_DIR / "portable_report_delivery_receipt.json",
    "portable_official_verification_screenshot": REPORT_DIR / "portable_report_official_verification.png",
    "portable_desktop_screenshot": REPORT_DIR / "portable_report_desktop_1440x1000.png",
    "portable_mobile_screenshot": REPORT_DIR / "portable_report_mobile_390x844.png",
    "portable_desktop_qa": REPORT_DIR / "portable_report_desktop_qa.json",
    "portable_mobile_qa": REPORT_DIR / "portable_report_mobile_qa.json",
}
EXPECTED_OUTPUTS = {
    "final_report_dataset": FINAL_DATASET_DIR / "final_report_dataset.json",
    "per_sample_metrics": FINAL_DATASET_DIR / "per_sample_metrics.tsv",
    "candidate_catalog": FINAL_DATASET_DIR / "candidate_catalog.tsv",
    "candidate_witness_pairs": FINAL_DATASET_DIR / "candidate_witness_pairs.tsv",
    "claim_ladder": FINAL_DATASET_DIR / "claim_ladder.tsv",
}
CREATE_CHECKS = frozenset(
    {
        "direct_canonical_process",
        "source_authority_verified",
        "builder_receipt_schema_scope_pass",
        "builder_direct_command_exact",
        "builder_source_identity_exact",
        "formal_validation_gates_pass",
        "all_final_output_artifacts_exact_and_read_only",
        "final_dataset_scope_and_counts_exact",
        "detached_signature_target_fresh",
    }
)
REPORT_CREATE_CHECKS = frozenset(
    {
        "direct_canonical_process",
        "source_authority_verified",
        "signed_final_dataset_release_reverified",
        "report_builder_receipt_exact",
        "portable_delivery_receipt_exact",
        "desktop_layout_qa_pass",
        "mobile_layout_qa_pass",
        "all_report_artifacts_exact_and_read_only",
        "detached_report_signature_target_fresh",
    }
)


class FinalReleaseError(RuntimeError):
    """Raised when the final Task-B result cannot be released."""


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise FinalReleaseError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_nonfinite_constant(value: str) -> None:
    raise FinalReleaseError(f"Non-finite JSON constant: {value}")


def parse_json_bytes(data: bytes, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            data.decode("utf-8"),
            object_pairs_hook=reject_duplicate_keys,
            parse_constant=reject_nonfinite_constant,
        )
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        FinalReleaseError,
    ) as error:
        raise FinalReleaseError(f"{label} is invalid JSON") from error
    if not isinstance(value, dict):
        raise FinalReleaseError(f"{label} is not an object")
    return value


def observed_process_command() -> list[str]:
    raw = Path("/proc/self/cmdline").read_bytes()
    if not raw.endswith(b"\0"):
        raise FinalReleaseError("Process command line is unavailable or malformed")
    return [os.fsdecode(token) for token in raw[:-1].split(b"\0")]


def canonical_python_prefix() -> list[str]:
    return [
        sys.executable,
        "-I",
        "-X",
        f"pycache_prefix={CANONICAL_PYTHON_CACHE_ROOT}",
    ]


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        data = path.read_bytes()
    except OSError as error:
        raise FinalReleaseError(f"Invalid {label}: {path}") from error
    return parse_json_bytes(data, label=label)


def require_mode(path: Path, mode: str, label: str) -> None:
    observed = oct(path.resolve(strict=True).stat().st_mode & 0o777)
    if observed != mode:
        raise FinalReleaseError(f"{label} mode drift: {observed} != {mode}")


def require_artifact(record: Any, path: Path, label: str) -> dict[str, Any]:
    expected = SOURCE_AUTHORITY.artifact(path)
    if not isinstance(record, Mapping) or dict(record) != expected:
        raise FinalReleaseError(f"{label} artifact identity drift")
    return expected


def create_command() -> list[str]:
    return [
        *canonical_python_prefix(),
        str(Path(__file__).resolve()),
        "--create",
        "--final-dataset-dir",
        str(FINAL_DATASET_DIR.resolve()),
        "--builder-receipt",
        str(BUILDER_RECEIPT.resolve()),
        "--output",
        str(RELEASE_RECEIPT.resolve()),
    ]


def verify_command() -> list[str]:
    return [
        *canonical_python_prefix(),
        str(Path(__file__).resolve()),
        "--verify",
        "--receipt",
        str(RELEASE_RECEIPT.resolve()),
        "--signature",
        str(RELEASE_SIGNATURE.resolve()),
    ]


def create_report_command() -> list[str]:
    return [
        *canonical_python_prefix(),
        str(Path(__file__).resolve()),
        "--create-report",
        "--report-dir",
        str(REPORT_DIR.resolve()),
        "--dataset-release-receipt",
        str(RELEASE_RECEIPT.resolve()),
        "--dataset-release-signature",
        str(RELEASE_SIGNATURE.resolve()),
        "--output",
        str(REPORT_RELEASE_RECEIPT.resolve()),
    ]


def verify_report_command() -> list[str]:
    return [
        *canonical_python_prefix(),
        str(Path(__file__).resolve()),
        "--verify-report",
        "--receipt",
        str(REPORT_RELEASE_RECEIPT.resolve()),
        "--signature",
        str(REPORT_RELEASE_SIGNATURE.resolve()),
    ]


def validate_builder_release() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    source_authority = SOURCE_AUTHORITY.validate_release_source_authority(
        SOURCE_AUTHORITY.EXPECTED_SOURCE_PATHS
    )
    require_mode(BUILDER_RECEIPT, "0o444", "final builder receipt")
    receipt = load_json(BUILDER_RECEIPT, "final builder receipt")
    validations = receipt.get("validations")
    outputs = receipt.get("outputs")
    if (
        receipt.get("schema_name")
        != "intersubmod.all_ssnv_final_report_dataset_run_receipt"
        or receipt.get("schema_version") != BUILDER.SCHEMA_VERSION
        or receipt.get("task_type") != "B_comprehensive_validation"
        or receipt.get("formal_task_b_release_eligible") is not True
        or receipt.get("pass_semantics") != BUILDER.PASS_SEMANTICS
        or receipt.get("pass") is not True
        or receipt.get("source_authority") != source_authority
        or not isinstance(validations, Mapping)
        or not isinstance(outputs, Mapping)
    ):
        raise FinalReleaseError("Final builder receipt schema/scope did not pass")
    expected_builder_commands = BUILDER.canonical_task_b_final_builder_commands(
        FINAL_DATASET_DIR
    )
    if receipt.get("command") not in expected_builder_commands:
        raise FinalReleaseError("Final builder direct command drift")
    code = receipt.get("code")
    if not isinstance(code, Mapping) or set(code) != {"final_report_dataset_builder"}:
        raise FinalReleaseError("Final builder source role drift")
    require_artifact(
        code.get("final_report_dataset_builder"),
        Path(BUILDER.__file__).resolve(),
        "final builder source",
    )
    if receipt.get("source_modes") != {"final_report_dataset_builder": "0o444"}:
        raise FinalReleaseError("Final builder source mode drift")
    strict_replay = validations.get("strict_statistics_deterministic_replay")
    primary_recount = validations.get("independent_primary_artifact_recount")
    if (
        validations.get("task_type_b_full_scope") is not True
        or validations.get("formal_task_b_release_eligible") is not True
        or not isinstance(strict_replay, Mapping)
        or strict_replay.get("pass") is not True
        or strict_replay.get("all_strict_output_fields_replayed") is not True
        or not isinstance(primary_recount, Mapping)
        or primary_recount.get("pass") is not True
        or primary_recount.get("implementation_independence") is not True
        or int(primary_recount.get("stable_sites", -1)) != 102_842
        or int(primary_recount.get("primary_artifacts_verified", -1)) != 308_526
    ):
        raise FinalReleaseError("Final builder formal replay/recount gates did not pass")
    observed_outputs: dict[str, dict[str, Any]] = {}
    if set(outputs) != set(EXPECTED_OUTPUTS):
        raise FinalReleaseError("Final builder output role set drift")
    for role, path in EXPECTED_OUTPUTS.items():
        require_mode(path, "0o444", f"final output {role}")
        observed_outputs[role] = require_artifact(
            outputs.get(role), path, f"final output {role}"
        )
    dataset = load_json(EXPECTED_OUTPUTS["final_report_dataset"], "final dataset")
    counts = dataset.get("counts")
    if (
        dataset.get("schema_name") != "intersubmod.all_ssnv_final_report_dataset"
        or dataset.get("task_type") != "B_comprehensive_validation"
        or dataset.get("formal_task_b_release_eligible") is not True
        or dataset.get("pass") is not True
        or not isinstance(counts, Mapping)
        or int(counts.get("screen_sites", -1)) != 469_849
        or int(counts.get("stable_sites", -1)) != 102_842
    ):
        raise FinalReleaseError("Final dataset full-scope/count contract drift")
    return receipt, source_authority, observed_outputs


def create_release() -> dict[str, Any]:
    if observed_process_command() != create_command():
        raise FinalReleaseError("Result finalizer was not executed as the direct canonical process")
    if os.path.lexists(RELEASE_RECEIPT) or os.path.lexists(RELEASE_SIGNATURE):
        raise FileExistsError("Refusing stale final release receipt or signature")
    builder_receipt, source_authority, outputs = validate_builder_release()
    public_key = SOURCE_AUTHORITY.artifact(RESULT_PUBLIC_KEY, include_mode=True)
    if (
        public_key["sha256"] != RESULT_PUBLIC_KEY_SHA256
        or public_key["mode"] != "0o444"
    ):
        raise FinalReleaseError("Result-signing public key drift")
    checks = {name: True for name in CREATE_CHECKS}
    payload = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "task_type": "B_comprehensive_validation",
        "scope": "all_7_datasets_469849_sites_final_dataset",
        "release_status": "AWAITING_ONE_TIME_ED25519_RESULT_SIGNATURE",
        "pass_semantics": "terminal_result_identity_release_not_scientific_confirmation",
        "command": create_command(),
        "source_authority": source_authority,
        "inputs": {
            "builder_receipt": SOURCE_AUTHORITY.artifact(BUILDER_RECEIPT),
            "final_outputs": outputs,
        },
        "code": {
            "result_release_finalizer": SOURCE_AUTHORITY.artifact(
                Path(__file__).resolve()
            )
        },
        "source_modes": {"result_release_finalizer": "0o444"},
        "result_signature_contract": {
            "algorithm": "Ed25519",
            "public_key": public_key,
            "signed_artifact": str(RELEASE_RECEIPT.resolve()),
            "signature": str(RELEASE_SIGNATURE.resolve()),
            "private_key_lifecycle": "encrypted_one_time_key_chmod_000_after_signing",
            "private_key_path": str(RESULT_PRIVATE_KEY.resolve()),
        },
        "builder_counts": builder_receipt.get("counts"),
        "checks": checks,
        "pass": True,
    }
    RELEASE_RECEIPT.parent.mkdir(parents=True, exist_ok=True)
    with RELEASE_RECEIPT.open("x", encoding="utf-8") as handle:
        json.dump(
            payload,
            handle,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        handle.write("\n")
    RELEASE_RECEIPT.chmod(0o444)
    return payload


def validate_release_signature_artifacts() -> dict[str, Any]:
    with SOURCE_AUTHORITY.BoundArtifactReader() as reader:
        _, private_key = reader.open_metadata(
            RESULT_PRIVATE_KEY, include_mode=True
        )
        if private_key["mode"] != "0o0":
            raise FinalReleaseError(
                "retired one-time result private key mode drift"
            )
        openssl_fd, _, openssl_artifact = reader.open(
            SOURCE_AUTHORITY.OPENSSL_PATH, include_mode=True
        )
        if openssl_artifact["sha256"] != SOURCE_AUTHORITY.OPENSSL_SHA256:
            raise FinalReleaseError("Result signature OpenSSL verifier drift")
        receipt_fd, receipt_bytes, receipt_artifact = reader.open(
            RELEASE_RECEIPT, include_mode=True
        )
        key_fd, _, public_key = reader.open(RESULT_PUBLIC_KEY, include_mode=True)
        signature_fd, _, signature = reader.open(
            RELEASE_SIGNATURE, include_mode=True
        )
        if (
            public_key["sha256"] != RESULT_PUBLIC_KEY_SHA256
            or public_key["mode"] != "0o444"
            or receipt_artifact["mode"] != "0o444"
            or signature["mode"] != "0o444"
        ):
            raise FinalReleaseError("Result receipt/key/signature identity drift")
        if not SOURCE_AUTHORITY.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=receipt_fd,
            public_key_fd=key_fd,
            signature_fd=signature_fd,
        ):
            raise FinalReleaseError("Final result Ed25519 signature verification failed")
        payload = parse_json_bytes(
            receipt_bytes,
            label="Signed final release receipt",
        )
        reader.retain_until_process_exit()
    _, source_authority, outputs = validate_builder_release()
    expected_signature_contract = {
        "algorithm": "Ed25519",
        "public_key": public_key,
        "signed_artifact": str(RELEASE_RECEIPT.resolve()),
        "signature": str(RELEASE_SIGNATURE.resolve()),
        "private_key_lifecycle": "encrypted_one_time_key_chmod_000_after_signing",
        "private_key_path": str(private_key["path"]),
    }
    if (
        payload.get("schema_name") != SCHEMA_NAME
        or payload.get("schema_version") != SCHEMA_VERSION
        or payload.get("task_type") != "B_comprehensive_validation"
        or payload.get("release_status")
        != "AWAITING_ONE_TIME_ED25519_RESULT_SIGNATURE"
        or payload.get("source_authority") != source_authority
        or payload.get("checks") != {name: True for name in CREATE_CHECKS}
        or payload.get("inputs", {}).get("final_outputs") != outputs
        or payload.get("result_signature_contract") != expected_signature_contract
        or payload.get("pass") is not True
    ):
        raise FinalReleaseError("Signed final release receipt content drift")
    return {
        "schema_name": f"{SCHEMA_NAME}.verification",
        "schema_version": SCHEMA_VERSION,
        "receipt": receipt_artifact,
        "signature": signature,
        "public_key": public_key,
        "source_authority": source_authority,
        "signature_verified": True,
        "all_final_outputs_reverified": True,
        "pass": True,
    }


def verify_release_signature() -> dict[str, Any]:
    if observed_process_command() != verify_command():
        raise FinalReleaseError("Result signature verifier process is not canonical")
    return validate_release_signature_artifacts()


def require_artifact_fields(record: Any, path: Path, label: str) -> dict[str, Any]:
    expected = SOURCE_AUTHORITY.artifact(path)
    if not isinstance(record, Mapping) or any(
        record.get(field) != value for field, value in expected.items()
    ):
        raise FinalReleaseError(f"{label} artifact identity drift")
    return expected


def validate_report_release_outputs() -> tuple[
    dict[str, Any], dict[str, Any], dict[str, dict[str, Any]]
]:
    dataset_signature_verification = validate_release_signature_artifacts()
    source_authority = SOURCE_AUTHORITY.validate_release_source_authority(
        SOURCE_AUTHORITY.EXPECTED_SOURCE_PATHS
    )
    if dataset_signature_verification.get("source_authority") != source_authority:
        raise FinalReleaseError("Signed dataset/report source authority drift")

    observed_outputs: dict[str, dict[str, Any]] = {}
    for role, path in REPORT_OUTPUTS.items():
        require_mode(path, "0o444", f"report output {role}")
        observed_outputs[role] = SOURCE_AUTHORITY.artifact(path)

    report_receipt = load_json(REPORT_OUTPUTS["report_build_receipt"], "report build receipt")
    report_outputs = report_receipt.get("outputs")
    report_code = report_receipt.get("code")
    if (
        report_receipt.get("schema_name")
        != "intersubmod.all_ssnv_report_build_receipt"
        or report_receipt.get("schema_version") != "1.2.0"
        or report_receipt.get("task_type") != "B_comprehensive_validation"
        or report_receipt.get("formal_task_b_release_eligible") is not True
        or report_receipt.get("pass") is not True
        or not isinstance(report_outputs, Mapping)
        or not isinstance(report_code, Mapping)
    ):
        raise FinalReleaseError("Formal report builder receipt did not pass")
    require_artifact_fields(
        report_outputs.get("report_md"),
        REPORT_OUTPUTS["report_markdown"],
        "report Markdown",
    )
    require_artifact_fields(
        report_outputs.get("artifact_json"),
        REPORT_OUTPUTS["canonical_report_artifact"],
        "canonical report artifact",
    )
    require_artifact_fields(
        report_code.get("report_builder"),
        SCRIPT_DIR / "build_all_ssnv_report_artifact.py",
        "report builder source",
    )
    require_artifact_fields(
        report_code.get("final_report_dataset_builder"),
        Path(BUILDER.__file__).resolve(),
        "final dataset builder source in report receipt",
    )

    portable_receipt = load_json(
        REPORT_OUTPUTS["portable_delivery_receipt"], "portable delivery receipt"
    )
    if (
        portable_receipt.get("schema_name")
        != "intersubmod.portable_report_delivery_receipt"
        or portable_receipt.get("schema_version") != "1.0.0"
        or portable_receipt.get("status") != "PASS"
        or portable_receipt.get("pass") is not True
    ):
        raise FinalReleaseError("Portable report delivery receipt did not pass")
    require_artifact(
        portable_receipt.get("artifact"),
        REPORT_OUTPUTS["canonical_report_artifact"],
        "portable input artifact",
    )
    require_artifact(
        portable_receipt.get("output"),
        REPORT_OUTPUTS["portable_html"],
        "portable HTML",
    )

    for role in ("portable_desktop_qa", "portable_mobile_qa"):
        qa = load_json(REPORT_OUTPUTS[role], role)
        if (
            qa.get("pass") is not True
            or int(qa.get("overlapCount", -1)) != 0
            or qa.get("consoleErrors") != []
            or int(qa.get("documentScrollWidth", 1))
            > int(qa.get("documentClientWidth", 0)) + 1
            or int(qa.get("bodyScrollWidth", 1))
            > int(qa.get("bodyClientWidth", 0)) + 1
        ):
            raise FinalReleaseError(f"{role} layout QA did not pass")
    for role in (
        "portable_official_verification_screenshot",
        "portable_desktop_screenshot",
        "portable_mobile_screenshot",
    ):
        with REPORT_OUTPUTS[role].open("rb") as handle:
            if handle.read(8) != b"\x89PNG\r\n\x1a\n":
                raise FinalReleaseError(f"{role} is not a PNG artifact")
    return source_authority, dataset_signature_verification, observed_outputs


def create_report_release() -> dict[str, Any]:
    if observed_process_command() != create_report_command():
        raise FinalReleaseError("Report finalizer was not executed as the direct canonical process")
    if os.path.lexists(REPORT_RELEASE_RECEIPT) or os.path.lexists(
        REPORT_RELEASE_SIGNATURE
    ):
        raise FileExistsError("Refusing stale final report release receipt or signature")
    source_authority, dataset_verification, outputs = validate_report_release_outputs()
    public_key = SOURCE_AUTHORITY.artifact(REPORT_PUBLIC_KEY, include_mode=True)
    if (
        public_key["sha256"] != REPORT_PUBLIC_KEY_SHA256
        or public_key["mode"] != "0o444"
    ):
        raise FinalReleaseError("Report-signing public key drift")
    payload = {
        "schema_name": REPORT_SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "task_type": "B_comprehensive_validation",
        "scope": "all_7_datasets_469849_sites_final_report_and_portable_evidence",
        "release_status": "AWAITING_ONE_TIME_ED25519_REPORT_SIGNATURE",
        "pass_semantics": "terminal_report_identity_release_not_scientific_confirmation",
        "command": create_report_command(),
        "source_authority": source_authority,
        "inputs": {
            "signed_dataset_release": {
                "receipt": SOURCE_AUTHORITY.artifact(RELEASE_RECEIPT),
                "signature": SOURCE_AUTHORITY.artifact(RELEASE_SIGNATURE),
                "public_key": SOURCE_AUTHORITY.artifact(
                    RESULT_PUBLIC_KEY, include_mode=True
                ),
                "verification_pass": dataset_verification.get("pass") is True,
            },
            "report_outputs": outputs,
        },
        "code": {
            "report_release_finalizer": SOURCE_AUTHORITY.artifact(
                Path(__file__).resolve()
            )
        },
        "source_modes": {"report_release_finalizer": "0o444"},
        "report_signature_contract": {
            "algorithm": "Ed25519",
            "public_key": public_key,
            "signed_artifact": str(REPORT_RELEASE_RECEIPT.resolve()),
            "signature": str(REPORT_RELEASE_SIGNATURE.resolve()),
            "private_key_lifecycle": "encrypted_one_time_key_chmod_000_after_signing",
            "private_key_path": str(REPORT_PRIVATE_KEY.resolve()),
        },
        "checks": {name: True for name in REPORT_CREATE_CHECKS},
        "pass": True,
    }
    REPORT_RELEASE_RECEIPT.parent.mkdir(parents=True, exist_ok=True)
    with REPORT_RELEASE_RECEIPT.open("x", encoding="utf-8") as handle:
        json.dump(
            payload,
            handle,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        handle.write("\n")
    REPORT_RELEASE_RECEIPT.chmod(0o444)
    return payload


def validate_report_release_signature_artifacts() -> dict[str, Any]:
    with SOURCE_AUTHORITY.BoundArtifactReader() as reader:
        _, private_key = reader.open_metadata(
            REPORT_PRIVATE_KEY, include_mode=True
        )
        if private_key["mode"] != "0o0":
            raise FinalReleaseError(
                "Retired one-time report private key mode drift"
            )
        openssl_fd, _, openssl_artifact = reader.open(
            SOURCE_AUTHORITY.OPENSSL_PATH, include_mode=True
        )
        if openssl_artifact["sha256"] != SOURCE_AUTHORITY.OPENSSL_SHA256:
            raise FinalReleaseError("Report signature OpenSSL verifier drift")
        receipt_fd, receipt_bytes, receipt_artifact = reader.open(
            REPORT_RELEASE_RECEIPT, include_mode=True
        )
        key_fd, _, public_key = reader.open(REPORT_PUBLIC_KEY, include_mode=True)
        signature_fd, _, signature = reader.open(
            REPORT_RELEASE_SIGNATURE, include_mode=True
        )
        if (
            public_key["sha256"] != REPORT_PUBLIC_KEY_SHA256
            or public_key["mode"] != "0o444"
            or receipt_artifact["mode"] != "0o444"
            or signature["mode"] != "0o444"
        ):
            raise FinalReleaseError("Report receipt/key/signature identity drift")
        if not SOURCE_AUTHORITY.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=receipt_fd,
            public_key_fd=key_fd,
            signature_fd=signature_fd,
        ):
            raise FinalReleaseError("Final report Ed25519 signature verification failed")
        payload = parse_json_bytes(
            receipt_bytes,
            label="Signed report release receipt",
        )
        reader.retain_until_process_exit()
    source_authority, dataset_verification, outputs = validate_report_release_outputs()
    expected_signature_contract = {
        "algorithm": "Ed25519",
        "public_key": public_key,
        "signed_artifact": str(REPORT_RELEASE_RECEIPT.resolve()),
        "signature": str(REPORT_RELEASE_SIGNATURE.resolve()),
        "private_key_lifecycle": "encrypted_one_time_key_chmod_000_after_signing",
        "private_key_path": str(private_key["path"]),
    }
    signed_dataset_release = payload.get("inputs", {}).get("signed_dataset_release", {})
    if (
        payload.get("schema_name") != REPORT_SCHEMA_NAME
        or payload.get("schema_version") != SCHEMA_VERSION
        or payload.get("task_type") != "B_comprehensive_validation"
        or payload.get("release_status")
        != "AWAITING_ONE_TIME_ED25519_REPORT_SIGNATURE"
        or payload.get("source_authority") != source_authority
        or payload.get("checks") != {name: True for name in REPORT_CREATE_CHECKS}
        or payload.get("inputs", {}).get("report_outputs") != outputs
        or signed_dataset_release.get("receipt")
        != SOURCE_AUTHORITY.artifact(RELEASE_RECEIPT)
        or signed_dataset_release.get("signature")
        != SOURCE_AUTHORITY.artifact(RELEASE_SIGNATURE)
        or signed_dataset_release.get("public_key")
        != SOURCE_AUTHORITY.artifact(RESULT_PUBLIC_KEY, include_mode=True)
        or signed_dataset_release.get("verification_pass")
        != (dataset_verification.get("pass") is True)
        or payload.get("report_signature_contract") != expected_signature_contract
        or payload.get("pass") is not True
    ):
        raise FinalReleaseError("Signed final report release receipt content drift")
    return {
        "schema_name": f"{REPORT_SCHEMA_NAME}.verification",
        "schema_version": SCHEMA_VERSION,
        "receipt": receipt_artifact,
        "signature": signature,
        "public_key": public_key,
        "source_authority": source_authority,
        "signed_dataset_release_reverified": True,
        "signature_verified": True,
        "all_report_outputs_reverified": True,
        "pass": True,
    }


def verify_report_release_signature() -> dict[str, Any]:
    if observed_process_command() != verify_report_command():
        raise FinalReleaseError("Report signature verifier process is not canonical")
    return validate_report_release_signature_artifacts()


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--create", action="store_true")
    mode.add_argument("--verify", action="store_true")
    mode.add_argument("--create-report", action="store_true")
    mode.add_argument("--verify-report", action="store_true")
    parser.add_argument("--final-dataset-dir", type=Path)
    parser.add_argument("--builder-receipt", type=Path)
    parser.add_argument("--report-dir", type=Path)
    parser.add_argument("--dataset-release-receipt", type=Path)
    parser.add_argument("--dataset-release-signature", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--signature", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    if argv is not None:
        raise FinalReleaseError("Formal result release finalizer is direct-CLI only")
    args = parse_args()
    if args.create:
        if (
            args.final_dataset_dir.resolve() != FINAL_DATASET_DIR.resolve()
            or args.builder_receipt.resolve() != BUILDER_RECEIPT.resolve()
            or args.output.resolve() != RELEASE_RECEIPT.resolve()
            or args.receipt is not None
            or args.signature is not None
            or args.report_dir is not None
            or args.dataset_release_receipt is not None
            or args.dataset_release_signature is not None
        ):
            raise FinalReleaseError("Result finalizer create paths are not canonical")
        payload = create_release()
        print(json.dumps({"output": str(RELEASE_RECEIPT), **payload}, indent=2))
    elif args.verify:
        if (
            args.receipt.resolve() != RELEASE_RECEIPT.resolve()
            or args.signature.resolve() != RELEASE_SIGNATURE.resolve()
            or any(
                value is not None
                for value in (
                    args.final_dataset_dir,
                    args.builder_receipt,
                    args.output,
                    args.report_dir,
                    args.dataset_release_receipt,
                    args.dataset_release_signature,
                )
            )
        ):
            raise FinalReleaseError("Result signature verification paths are not canonical")
        print(json.dumps(verify_release_signature(), indent=2))
    elif args.create_report:
        if (
            args.report_dir.resolve() != REPORT_DIR.resolve()
            or args.dataset_release_receipt.resolve() != RELEASE_RECEIPT.resolve()
            or args.dataset_release_signature.resolve() != RELEASE_SIGNATURE.resolve()
            or args.output.resolve() != REPORT_RELEASE_RECEIPT.resolve()
            or any(
                value is not None
                for value in (
                    args.final_dataset_dir,
                    args.builder_receipt,
                    args.receipt,
                    args.signature,
                )
            )
        ):
            raise FinalReleaseError("Report finalizer create paths are not canonical")
        payload = create_report_release()
        print(json.dumps({"output": str(REPORT_RELEASE_RECEIPT), **payload}, indent=2))
    else:
        if (
            args.receipt.resolve() != REPORT_RELEASE_RECEIPT.resolve()
            or args.signature.resolve() != REPORT_RELEASE_SIGNATURE.resolve()
            or any(
                value is not None
                for value in (
                    args.final_dataset_dir,
                    args.builder_receipt,
                    args.report_dir,
                    args.dataset_release_receipt,
                    args.dataset_release_signature,
                    args.output,
                )
            )
        ):
            raise FinalReleaseError("Report signature verification paths are not canonical")
        print(json.dumps(verify_report_release_signature(), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
