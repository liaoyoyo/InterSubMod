#!/usr/bin/env python3
"""Create a signed-ready supplemental positional-singleton release receipt."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
from typing import Any
import uuid


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import release_source_authority as SOURCE_AUTHORITY  # noqa: E402


OPENSSL = Path("/usr/bin/openssl")
FINAL_DATASET_PUBLIC_KEY_SHA256 = (
    "5b7d5d026835ec6ec677bcd886bc16ac71444117dabc44a084ce3ede3a4db5a9"
)
FINAL_REPORT_PUBLIC_KEY_SHA256 = (
    "572e27167e1eea4c39ca53546ba33868b2874ec7fa3ed1682db821b2e50fa439"
)
SUPPLEMENTAL_PUBLIC_KEY_SHA256 = (
    "c8490467e7ad9c0477429d37c799dd4143773b781b6fc33ab08199a67dd99b91"
)
INPUT_NAMES = (
    "audit_summary",
    "audit_success",
    "report_receipt",
    "report_success",
    "final_dataset_receipt",
    "final_dataset_signature",
    "final_dataset_public_key",
    "final_report_receipt",
    "final_report_signature",
    "final_report_public_key",
    "supplemental_public_key",
)
JSON_INPUT_LABELS = {
    "audit_summary": "audit summary",
    "audit_success": "audit success marker",
    "report_receipt": "supplemental report receipt",
    "report_success": "supplemental report success marker",
    "final_dataset_receipt": "final dataset receipt",
    "final_report_receipt": "final report receipt",
}


class ReleaseError(RuntimeError):
    """Raised when a supplemental release contract is not satisfied."""


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ReleaseError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_nonfinite_constant(value: str) -> None:
    raise ReleaseError(f"Non-finite JSON constant: {value}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve(strict=True)
    info = resolved.stat()
    return {
        "path": str(resolved),
        "size_bytes": info.st_size,
        "sha256": sha256(resolved),
        "mode": oct(stat.S_IMODE(info.st_mode)),
    }


def parse_json_bytes(data: bytes, *, path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            data.decode("utf-8"),
            object_pairs_hook=reject_duplicate_keys,
            parse_constant=reject_nonfinite_constant,
        )
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        ReleaseError,
    ) as error:
        raise ReleaseError(f"Unable to read {label}: {path}") from error
    if not isinstance(value, dict):
        raise ReleaseError(f"{label} must be a JSON object")
    return value


def require_identity(
    observed: dict[str, Any], record: Any, *, label: str
) -> None:
    if not isinstance(record, dict):
        raise ReleaseError(f"{label} identity is missing")
    for field in ("path", "size_bytes", "sha256", "mode"):
        if observed[field] != record.get(field):
            raise ReleaseError(
                f"{label} identity drift for {field}: "
                f"{observed[field]!r} != {record.get(field)!r}"
            )


def validate_bound_inputs(
    reader: SOURCE_AUTHORITY.BoundArtifactReader,
    paths: dict[str, Path],
) -> tuple[
    dict[str, dict[str, Any]],
    dict[str, dict[str, Any]],
    dict[str, Any],
]:
    source_fd, _, source_artifact = reader.open(
        Path(__file__), include_mode=True
    )
    del source_fd
    openssl_fd, _, openssl_artifact = reader.open(
        SOURCE_AUTHORITY.OPENSSL_PATH, include_mode=True
    )
    if openssl_artifact["sha256"] != SOURCE_AUTHORITY.OPENSSL_SHA256:
        raise ReleaseError("Supplemental release OpenSSL verifier drift")

    descriptors: dict[str, int] = {}
    contents: dict[str, bytes] = {}
    identities: dict[str, dict[str, Any]] = {}
    for name in INPUT_NAMES:
        descriptor, data, record = reader.open(paths[name], include_mode=True)
        descriptors[name] = descriptor
        contents[name] = data
        identities[name] = record
    wrong_modes = {
        name: record["mode"]
        for name, record in identities.items()
        if record["mode"] != "0o444"
    }
    if wrong_modes:
        raise ReleaseError(f"Supplemental input protected mode drift: {wrong_modes}")

    payloads = {
        name: parse_json_bytes(
            contents[name],
            path=paths[name],
            label=label,
        )
        for name, label in JSON_INPUT_LABELS.items()
    }
    summary = payloads["audit_summary"]
    audit_success = payloads["audit_success"]
    report_receipt = payloads["report_receipt"]
    report_success = payloads["report_success"]
    dataset_receipt = payloads["final_dataset_receipt"]
    formal_report_receipt = payloads["final_report_receipt"]

    summary_checks = summary.get("checks")
    if (
        summary.get("schema_name")
        != "intersubmod.positional_singleton_methyl_multigroup_audit"
        or summary.get("schema_version") != "2.0.0"
        or summary.get("pass") is not True
        or not isinstance(summary_checks, dict)
        or not summary_checks
        or not all(value is True for value in summary_checks.values())
    ):
        raise ReleaseError("Positional-singleton audit contract drift")
    if (
        audit_success.get("schema_name")
        != "intersubmod.atomic_release_marker"
        or audit_success.get("pass") is not True
    ):
        raise ReleaseError("Audit success-marker contract drift")
    require_identity(
        identities["audit_summary"],
        audit_success.get("summary"),
        label="audit summary success binding",
    )
    if (
        report_receipt.get("schema_name")
        != "intersubmod.positional_singleton_methyl_multigroup_report"
        or report_receipt.get("schema_version") != "1.0.0"
        or report_receipt.get("pass") is not True
        or report_receipt.get("claim_ceiling")
        != "M2_read_level_residual_epigenetic_partition"
    ):
        raise ReleaseError("Supplemental report receipt contract drift")
    require_identity(
        identities["audit_summary"],
        report_receipt.get("inputs", {}).get("audit_summary"),
        label="supplemental report audit-summary binding",
    )
    if not all(
        report_receipt.get("formal_release_links", {}).get(field) is True
        for field in (
            "final_dataset_receipt_pass",
            "final_report_receipt_pass",
            "final_dataset_builder_receipt_pass",
        )
    ):
        raise ReleaseError("Supplemental report formal-release link drift")
    if (
        report_success.get("schema_name")
        != "intersubmod.atomic_release_marker"
        or report_success.get("pass") is not True
    ):
        raise ReleaseError("Supplemental report success-marker contract drift")
    require_identity(
        identities["report_receipt"],
        report_success.get("receipt"),
        label="supplemental report success binding",
    )
    if (
        dataset_receipt.get("schema_name")
        != "intersubmod.task_b_final_dataset_release_receipt"
        or dataset_receipt.get("pass") is not True
        or formal_report_receipt.get("schema_name")
        != "intersubmod.task_b_final_report_release_receipt"
        or formal_report_receipt.get("pass") is not True
    ):
        raise ReleaseError("Formal Task-B release contract drift")

    signature_specs = (
        (
            "final_dataset_receipt",
            "final_dataset_signature",
            "final_dataset_public_key",
            FINAL_DATASET_PUBLIC_KEY_SHA256,
            "final dataset release",
        ),
        (
            "final_report_receipt",
            "final_report_signature",
            "final_report_public_key",
            FINAL_REPORT_PUBLIC_KEY_SHA256,
            "final report release",
        ),
    )
    for receipt_name, signature_name, key_name, expected_key_sha, label in (
        signature_specs
    ):
        if identities[key_name]["sha256"] != expected_key_sha:
            raise ReleaseError(f"{label} public-key SHA-256 drift")
        if not SOURCE_AUTHORITY.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=descriptors[receipt_name],
            public_key_fd=descriptors[key_name],
            signature_fd=descriptors[signature_name],
        ):
            raise ReleaseError(f"{label} detached signature failed")
    if (
        identities["supplemental_public_key"]["sha256"]
        != SUPPLEMENTAL_PUBLIC_KEY_SHA256
    ):
        raise ReleaseError("Supplemental integrity public-key drift")

    return identities, payloads, source_artifact


def verify_signature(
    *,
    receipt: Path,
    signature: Path,
    public_key: Path,
    expected_public_key_sha256: str,
    label: str,
) -> None:
    try:
        with SOURCE_AUTHORITY.BoundArtifactReader() as reader:
            openssl_fd, _, openssl_artifact = reader.open(
                SOURCE_AUTHORITY.OPENSSL_PATH, include_mode=True
            )
            receipt_fd, _, receipt_artifact = reader.open(
                receipt, include_mode=True
            )
            signature_fd, _, signature_artifact = reader.open(
                signature, include_mode=True
            )
            public_key_fd, _, public_key_artifact = reader.open(
                public_key, include_mode=True
            )
            if openssl_artifact["sha256"] != SOURCE_AUTHORITY.OPENSSL_SHA256:
                raise ReleaseError(f"{label} OpenSSL verifier drift")
            if any(
                record["mode"] != "0o444"
                for record in (
                    receipt_artifact,
                    signature_artifact,
                    public_key_artifact,
                )
            ):
                raise ReleaseError(f"{label} protected mode drift")
            if public_key_artifact["sha256"] != expected_public_key_sha256:
                raise ReleaseError(f"{label} public-key SHA-256 drift")
            if not SOURCE_AUTHORITY.verify_ed25519_signature_fds(
                openssl_fd=openssl_fd,
                data_fd=receipt_fd,
                public_key_fd=public_key_fd,
                signature_fd=signature_fd,
            ):
                raise ReleaseError(f"{label} detached signature failed")
            reader.retain_until_process_exit()
    except SOURCE_AUTHORITY.SourceAuthorityError as error:
        raise ReleaseError(
            f"{label} bound signature verification failed: {error}"
        ) from error


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_create(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing path: {path}")
    pending = path.parent / f".{path.name}.pending.{os.getpid()}.{uuid.uuid4().hex}"
    encoded = (
        json.dumps(
            payload,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    descriptor = os.open(
        pending,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL,
        0o444,
    )
    try:
        offset = 0
        while offset < len(encoded):
            offset += os.write(descriptor, encoded[offset:])
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    try:
        os.link(pending, path)
        pending.unlink()
        fsync_directory(path.parent)
    except BaseException:
        pending.unlink(missing_ok=True)
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("create", "verify"),
        default="create",
    )
    parser.add_argument("--audit-summary", type=Path, required=True)
    parser.add_argument("--audit-success", type=Path, required=True)
    parser.add_argument("--report-receipt", type=Path, required=True)
    parser.add_argument("--report-success", type=Path, required=True)
    parser.add_argument("--final-dataset-receipt", type=Path, required=True)
    parser.add_argument("--final-dataset-signature", type=Path, required=True)
    parser.add_argument("--final-dataset-public-key", type=Path, required=True)
    parser.add_argument("--final-report-receipt", type=Path, required=True)
    parser.add_argument("--final-report-signature", type=Path, required=True)
    parser.add_argument("--final-report-public-key", type=Path, required=True)
    parser.add_argument("--supplemental-public-key", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--supplemental-signature", type=Path)
    parser.add_argument("--supplemental-private-key", type=Path)
    parser.add_argument("--verification-output", type=Path)
    return parser.parse_args()


def resolve_input_paths(args: argparse.Namespace) -> dict[str, Path]:
    return {
        name: getattr(args, name).expanduser().resolve(strict=True)
        for name in INPUT_NAMES
    }


def build_release_payload(
    *,
    output: Path,
    identities: dict[str, dict[str, Any]],
    payloads: dict[str, dict[str, Any]],
    source_artifact: dict[str, Any],
) -> dict[str, Any]:
    summary = payloads["audit_summary"]
    return {
        "schema_name": (
            "intersubmod.positional_singleton_supplemental_release_receipt"
        ),
        "schema_version": "1.0.0",
        "task_type": "B_comprehensive_validation",
        "scope": "all_50432_positional_singleton_dataset_sites",
        "command": [str(Path(sys.executable).resolve()), *sys.argv],
        "inputs": identities,
        "code": {"supplemental_release_finalizer": source_artifact},
        "counts": {
            "singleton_sites": summary["counts"]["singleton_sites"],
            "m1_evaluable": summary["counts"]["m1_evaluable"],
            "m1_flagged": summary["counts"]["m1_flagged"],
            "m2_pass": summary["counts"]["m2_pass"],
        },
        "claim_ceiling": "M2_read_level_residual_epigenetic_partition",
        "sealed_artifact_scope": (
            "After the required out-of-band detached signature is created, "
            "this receipt directly seals the audit summary, supplemental "
            "report receipt, both atomic success markers, both formal signed "
            "Task-B receipts, and all identities transitively recorded by "
            "the supplemental report receipt."
        ),
        "signature_contract": {
            "algorithm": "Ed25519",
            "state_at_receipt_creation": "UNSIGNED_PENDING_OUT_OF_BAND",
            "expected_signature_path": str(output) + ".ed25519.sig",
            "public_key": identities["supplemental_public_key"],
            "required_signature_mode": "0o444",
            "required_private_key_mode_after_signing": "0o0",
        },
        "authority_semantics": (
            "supplemental_integrity_seal_after_formal_task_b_release;"
            "not_v4_release_reauthorization;"
            "pass_requires_subsequent_detached_signature_and_verification"
        ),
        "pass": True,
    }


def create_release(args: argparse.Namespace) -> dict[str, Any]:
    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing path: {output}")
    paths = resolve_input_paths(args)
    with SOURCE_AUTHORITY.BoundArtifactReader() as reader:
        identities, payloads, source_artifact = validate_bound_inputs(
            reader, paths
        )
        payload = build_release_payload(
            output=output,
            identities=identities,
            payloads=payloads,
            source_artifact=source_artifact,
        )
        reader.retain_until_process_exit()
    atomic_create(output, payload)
    return {
        "output": str(output),
        "sha256": sha256(output),
        "mode": oct(stat.S_IMODE(output.stat().st_mode)),
        "release_state": "UNSIGNED_PENDING_OUT_OF_BAND",
        "pass": True,
    }


def verify_release(args: argparse.Namespace) -> dict[str, Any]:
    required = {
        "supplemental_signature": args.supplemental_signature,
        "supplemental_private_key": args.supplemental_private_key,
        "verification_output": args.verification_output,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise ReleaseError(
            f"Verify mode requires arguments: {', '.join(sorted(missing))}"
        )
    receipt_path = args.output.expanduser().resolve(strict=True)
    signature_path = args.supplemental_signature.expanduser().resolve(strict=True)
    private_key_path = (
        args.supplemental_private_key.expanduser().resolve(strict=True)
    )
    verification_output = args.verification_output.expanduser().resolve()
    if verification_output.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing path: {verification_output}"
        )
    expected_signature_path = Path(str(receipt_path) + ".ed25519.sig")
    if signature_path != expected_signature_path:
        raise ReleaseError(
            "Supplemental signature path differs from signed receipt contract"
        )

    paths = resolve_input_paths(args)
    with SOURCE_AUTHORITY.BoundArtifactReader() as reader:
        identities, payloads, source_artifact = validate_bound_inputs(
            reader, paths
        )
        openssl_fd, _, openssl_artifact = reader.open(
            SOURCE_AUTHORITY.OPENSSL_PATH, include_mode=True
        )
        receipt_fd, receipt_bytes, receipt_artifact = reader.open(
            receipt_path, include_mode=True
        )
        signature_fd, _, signature_artifact = reader.open(
            signature_path, include_mode=True
        )
        public_key_fd, _, public_key_artifact = reader.open(
            paths["supplemental_public_key"], include_mode=True
        )
        _, private_key_artifact = reader.open_metadata(
            private_key_path, include_mode=True
        )
        if openssl_artifact["sha256"] != SOURCE_AUTHORITY.OPENSSL_SHA256:
            raise ReleaseError("Supplemental post-sign OpenSSL verifier drift")
        if receipt_artifact["mode"] != "0o444":
            raise ReleaseError("Supplemental signed receipt mode drift")
        if signature_artifact["mode"] != "0o444":
            raise ReleaseError("Supplemental detached signature mode drift")
        if public_key_artifact != identities["supplemental_public_key"]:
            raise ReleaseError("Supplemental post-sign public-key identity drift")
        if private_key_artifact["mode"] != "0o0":
            raise ReleaseError("Supplemental one-time private key is not retired")
        if not SOURCE_AUTHORITY.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=receipt_fd,
            public_key_fd=public_key_fd,
            signature_fd=signature_fd,
        ):
            raise ReleaseError("Supplemental detached signature failed")
        receipt = parse_json_bytes(
            receipt_bytes,
            path=receipt_path,
            label="signed supplemental release receipt",
        )
        expected_signature_contract = {
            "algorithm": "Ed25519",
            "state_at_receipt_creation": "UNSIGNED_PENDING_OUT_OF_BAND",
            "expected_signature_path": str(expected_signature_path),
            "public_key": identities["supplemental_public_key"],
            "required_signature_mode": "0o444",
            "required_private_key_mode_after_signing": "0o0",
        }
        expected_counts = {
            name: payloads["audit_summary"]["counts"][name]
            for name in (
                "singleton_sites",
                "m1_evaluable",
                "m1_flagged",
                "m2_pass",
            )
        }
        command = receipt.get("command")
        if (
            receipt.get("schema_name")
            != "intersubmod.positional_singleton_supplemental_release_receipt"
            or receipt.get("schema_version") != "1.0.0"
            or receipt.get("task_type") != "B_comprehensive_validation"
            or receipt.get("scope")
            != "all_50432_positional_singleton_dataset_sites"
            or not isinstance(command, list)
            or str(receipt_path) not in command
            or receipt.get("inputs") != identities
            or receipt.get("code")
            != {"supplemental_release_finalizer": source_artifact}
            or receipt.get("counts") != expected_counts
            or receipt.get("claim_ceiling")
            != "M2_read_level_residual_epigenetic_partition"
            or receipt.get("signature_contract") != expected_signature_contract
            or receipt.get("pass") is not True
        ):
            raise ReleaseError("Signed supplemental release receipt contract drift")
        verification_payload = {
            "schema_name": (
                "intersubmod.positional_singleton_supplemental_"
                "post_signature_verification"
            ),
            "schema_version": "1.0.0",
            "task_type": "B_comprehensive_validation",
            "scope": "all_50432_positional_singleton_dataset_sites",
            "verified_receipt": receipt_artifact,
            "detached_signature": signature_artifact,
            "public_key": public_key_artifact,
            "retired_private_key": private_key_artifact,
            "live_inputs": identities,
            "code": {"supplemental_release_finalizer": source_artifact},
            "checks": {
                "receipt_contract_recomputed": True,
                "detached_ed25519_signature_verified": True,
                "one_time_private_key_retired": True,
                "all_bound_direct_input_identities_reverified": True,
                "all_bound_descriptors_retained_until_process_exit": True,
            },
            "reverification_scope": {
                "identity_count": len(identities),
                "scope": "direct_receipt_inputs_only",
                "transitive_artifacts_cryptographically_bound_but_not_reopened": True,
            },
            "verification_status": "SIGNED_SUPPLEMENTAL_RELEASE_REVERIFIED",
            "pass": True,
        }
        reader.retain_until_process_exit()
    atomic_create(verification_output, verification_payload)
    return {
        "output": str(verification_output),
        "sha256": sha256(verification_output),
        "mode": oct(stat.S_IMODE(verification_output.stat().st_mode)),
        "verification_status": "SIGNED_SUPPLEMENTAL_RELEASE_REVERIFIED",
        "pass": True,
    }


def main() -> None:
    args = parse_args()
    result = create_release(args) if args.mode == "create" else verify_release(args)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
