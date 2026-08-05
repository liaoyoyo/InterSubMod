#!/usr/bin/env python3
"""Archive the signed v14 generation after its fail-closed R-stage error."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
from typing import Any


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
RESULT_ROOT = TOPIC_ROOT / "results"
REVIEW_ROOT = TOPIC_ROOT / "reviews"
WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
FAILURE_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v14_signed_authority_r_metadata_only_relation_schema_mismatch"
)

PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")
REPLAYER = AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v14.py"
FORMAL_ENV = {
    "PATH": "/usr/bin:/bin",
    "HOME": "/tmp",
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "PYTHONHASHSEED": "0",
    "PYTHONNOUSERSITE": "1",
    "PYTHONDONTWRITEBYTECODE": "1",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "BLIS_NUM_THREADS": "1",
}
FORMAL_COMMAND = [str(PYTHON), "-I", "-B", str(REPLAYER)]

AUTHORITY_BUNDLE = RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v14.bundle"
VERIFICATION_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v14.json"
)
REVIEWS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v14.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v14.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v14.json"
    ),
}
SOURCES = {
    "authority_validator": AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v14.py",
    "continuation_verifier": AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v14.py",
    "runner_gate_replay": REPLAYER,
    "downstream_continuation": AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v14.py",
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v14.py",
    "regression_tests": AUDIT_ROOT / "schema_recovery_tests/test_tumor_ref_schema_recovery_v14.py",
    "ceremony_builder": AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v14.py",
    "archive_script": Path(__file__).resolve(strict=True),
}

AUTHORITY_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v14/ed25519_public.pem"
)
AUTHORITY_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v14/ed25519_private_one_time.pem"
)
CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v4/ed25519_public.pem"
)
CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v4/ed25519_private_one_time_resident.pem"
)

R_OUTPUTS = (
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v14.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v14.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v14.success_witness.json",
)
C_OUTPUTS = (
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v14.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v14.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v14.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v14.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v14.json",
)
UNUSED_COMPATIBILITY_OUTPUTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v14.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v14.json.ed25519.sig",
)
CANONICAL_DOWNSTREAM_OUTPUTS = (
    WORKSPACE_ROOT / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "all_ssnv_final_report_dataset_v5_m2v5_source_attested",
    WORKSPACE_ROOT / "all_ssnv_final_report_v5_m2v5_source_attested",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json",
    RESULT_ROOT / "positional_singleton_supplemental_release_receipt.v1.json",
)

EXPECTED_ERROR = (
    "Metadata-only artifact cannot be reopened for bytes: "
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260722_m2v5_terminal_v2/ed25519_private_one_time_resident.pem"
)


class ArchiveError(RuntimeError):
    pass


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def artifact_record(path: Path) -> dict[str, Any]:
    leaf = os.lstat(path)
    if not stat.S_ISREG(leaf.st_mode) or path.resolve(strict=True) != path:
        raise ArchiveError(f"Artifact is not one canonical regular file: {path}")
    data = path.read_bytes()
    after = os.lstat(path)
    if (leaf.st_dev, leaf.st_ino, leaf.st_size, leaf.st_mtime_ns, leaf.st_ctime_ns) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    ):
        raise ArchiveError(f"Artifact changed while hashing: {path}")
    return {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": sha256_bytes(data),
        "mode": oct(stat.S_IMODE(after.st_mode)),
    }


def metadata_record(path: Path) -> dict[str, Any]:
    leaf = os.lstat(path)
    if not stat.S_ISREG(leaf.st_mode) or path.resolve(strict=True) != path:
        raise ArchiveError(f"Metadata path is not one canonical regular file: {path}")
    return {
        "path": str(path),
        "size_bytes": leaf.st_size,
        "mode": oct(stat.S_IMODE(leaf.st_mode)),
    }


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ArchiveError(f"JSON root is not an object: {path}")
    return value


def write_exclusive(path: Path, data: bytes, mode: int = 0o444) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC, 0o600)
    try:
        offset = 0
        while offset < len(data):
            offset += os.write(descriptor, data[offset:])
        os.fsync(descriptor)
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def require_absent(paths: tuple[Path, ...], label: str) -> None:
    present = [str(path) for path in paths if path.exists() or path.is_symlink()]
    if present:
        raise ArchiveError(f"{label} unexpectedly present: {present}")


def move_no_replace(source: Path, destination: Path) -> None:
    if not source.exists() or source.is_symlink():
        raise ArchiveError(f"Archive source missing or invalid: {source}")
    if destination.exists() or destination.is_symlink():
        raise ArchiveError(f"Archive destination already exists: {destination}")
    os.rename(source, destination)


def record_at_final_path(staging_path: Path, final_path: Path) -> dict[str, Any]:
    record = artifact_record(staging_path)
    record["path"] = str(final_path)
    return record


def main() -> int:
    if FAILURE_ROOT.exists() or FAILURE_ROOT.is_symlink():
        raise ArchiveError(f"Failure archive already exists: {FAILURE_ROOT}")
    require_absent(R_OUTPUTS, "v14 R outputs before diagnostic replay")
    require_absent(C_OUTPUTS, "v14 C outputs before diagnostic replay")
    require_absent(CANONICAL_DOWNSTREAM_OUTPUTS, "canonical downstream outputs")
    if not AUTHORITY_BUNDLE.is_dir() or not VERIFICATION_RECEIPT.is_file():
        raise ArchiveError("v14 signed authority or V receipt is missing")
    if any(not path.is_file() for path in REVIEWS.values()):
        raise ArchiveError("One or more v14 formal reviews are missing")

    completed = subprocess.run(
        FORMAL_COMMAND,
        cwd=REPO_ROOT,
        env=FORMAL_ENV,
        capture_output=True,
        check=False,
        timeout=300,
    )
    stderr_text = completed.stderr.decode("utf-8", errors="strict")
    if completed.returncode != 1 or completed.stdout != b"" or EXPECTED_ERROR not in stderr_text:
        raise ArchiveError(
            "Formal v14 diagnostic replay did not reproduce the reviewed fail-closed result"
        )
    wait_match = re.search(
        r"Runner-gate worker did not produce waitpid normal exit zero: "
        r"pid=(\d+) waited_pid=(\d+) status=(\d+)",
        stderr_text,
    )
    if wait_match is None or wait_match.group(1) != wait_match.group(2):
        raise ArchiveError("Formal v14 waitpid failure evidence is incomplete")
    wait_status = int(wait_match.group(3))
    if (
        wait_status != 17920
        or not os.WIFEXITED(wait_status)
        or os.WEXITSTATUS(wait_status) != 70
    ):
        raise ArchiveError(f"Unexpected v14 worker wait status: {wait_status}")
    require_absent(R_OUTPUTS, "v14 R outputs after diagnostic replay")
    require_absent(C_OUTPUTS, "v14 C outputs after diagnostic replay")

    authority_payload = load_json(AUTHORITY_BUNDLE / "authority.json")
    review_payloads = {role: load_json(path) for role, path in REVIEWS.items()}
    if authority_payload.get("authority_id") != "20260723_tumor_ref_schema_recovery_v14":
        raise ArchiveError("v14 authority identity drift")
    if authority_payload.get("pass") is not True:
        raise ArchiveError("v14 authority does not declare pass=true")
    if any(payload.get("verdict") != "APPROVE" or payload.get("pass") is not True for payload in review_payloads.values()):
        raise ArchiveError("v14 formal review verdict drift")

    authority_private = metadata_record(AUTHORITY_PRIVATE_KEY)
    continuation_private = metadata_record(CONTINUATION_PRIVATE_KEY)
    if authority_private["mode"] != "0o0":
        raise ArchiveError("v14 authority private key was not retired")
    if continuation_private["mode"] != "0o400":
        raise ArchiveError("v14 continuation private key is not unused mode 0400")

    staging = FAILURE_ROOT.parent / f".{FAILURE_ROOT.name}.staging.{os.getpid()}"
    staging.mkdir(mode=0o755)
    try:
        move_no_replace(AUTHORITY_BUNDLE, staging / AUTHORITY_BUNDLE.name)
        move_no_replace(VERIFICATION_RECEIPT, staging / VERIFICATION_RECEIPT.name)
        for source in REVIEWS.values():
            move_no_replace(source, staging / source.name)

        stderr_path = staging / "formal_replay_stderr.log"
        write_exclusive(stderr_path, completed.stderr)

        archived = {
            "authority": record_at_final_path(
                staging / AUTHORITY_BUNDLE.name / "authority.json",
                FAILURE_ROOT / AUTHORITY_BUNDLE.name / "authority.json",
            ),
            "authority_signature": record_at_final_path(
                staging / AUTHORITY_BUNDLE.name / "authority.ed25519.sig",
                FAILURE_ROOT / AUTHORITY_BUNDLE.name / "authority.ed25519.sig",
            ),
            "authority_commit": record_at_final_path(
                staging / AUTHORITY_BUNDLE.name / "commit.json",
                FAILURE_ROOT / AUTHORITY_BUNDLE.name / "commit.json",
            ),
            "verification_receipt": record_at_final_path(
                staging / VERIFICATION_RECEIPT.name,
                FAILURE_ROOT / VERIFICATION_RECEIPT.name,
            ),
            "formal_replay_stderr": record_at_final_path(
                stderr_path,
                FAILURE_ROOT / stderr_path.name,
            ),
        }
        for role, source in REVIEWS.items():
            archived[f"review_{role}"] = record_at_final_path(
                staging / source.name,
                FAILURE_ROOT / source.name,
            )

        reviewed_hashes = next(iter(review_payloads.values()))
        evidence = {
            "schema_name": "intersubmod.tumor_ref_schema_recovery_signed_generation_failure",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "generation": "v14",
            "status": "SIGNED_AUTHORITY_AND_V_PASS_R_FAILED_BEFORE_FIRST_R_OUTPUT",
            "review_contract": {
                "reviewed_source_set_sha256": reviewed_hashes["reviewed_source_set_sha256"],
                "legacy_source_set_sha256": reviewed_hashes["legacy_source_set_sha256"],
                "prior_recovery_chain_sha256": reviewed_hashes["prior_recovery_chain_sha256"],
                "rejected_generations_sha256": reviewed_hashes["rejected_generations_sha256"],
                "scope_sha256": reviewed_hashes["scope_sha256"],
                "terminal_key_rotation_sha256": reviewed_hashes["terminal_key_rotation_sha256"],
            },
            "runtime": {
                "v_exit_code": 0,
                "r_supervisor_exit_code": completed.returncode,
                "r_worker_exit_code": 70,
                "r_waitpid_raw_status": wait_status,
                "r_worker_normal_exit": True,
                "r_worker_success": False,
                "c_started": False,
                "formal_command": FORMAL_COMMAND,
                "clean_environment": FORMAL_ENV,
                "stdout": {
                    "size_bytes": len(completed.stdout),
                    "sha256": sha256_bytes(completed.stdout),
                },
                "stderr": archived["formal_replay_stderr"],
                "error_type": "GateError",
                "error_message": EXPECTED_ERROR,
                "failure_location": (
                    "runner_gate_replay.recursively_validate_artifact_relations."
                    "schema_recovery_authority.terminal_key_state.legacy_private_key"
                ),
            },
            "diagnostic": {
                "record_shape": ["mode", "path", "size_bytes", "state"],
                "private_key_bytes_read": False,
                "root_cause": (
                    "The recursive relation classifier treated size_bytes without sha256 as a "
                    "content-readable artifact after the same path had been descriptor-bound "
                    "metadata-only. The reader correctly rejected that schema conflict."
                ),
                "required_fix": (
                    "Add an exact metadata-plus-size relation type that validates canonical path, "
                    "mode, size, inode binding, and registry compatibility without reading private "
                    "key bytes."
                ),
                "failure_before_first_r_output": True,
                "failure_before_c_start": True,
                "scientific_payload_changed": False,
            },
            "archived_artifacts": archived,
            "source_artifacts": {role: artifact_record(path) for role, path in SOURCES.items()},
            "reviews": {
                role: {
                    "reviewer": payload["reviewer"],
                    "reviewer_agent_id": payload["reviewer_agent_id"],
                    "verdict": payload["verdict"],
                    "findings_empty": not any(
                        payload.get(key)
                        for key in ("high_findings", "medium_findings", "low_findings")
                    ),
                }
                for role, payload in review_payloads.items()
            },
            "keys": {
                "authority_public_key": artifact_record(AUTHORITY_PUBLIC_KEY),
                "authority_private_key": {
                    **authority_private,
                    "retired_after_single_signature": True,
                },
                "continuation_public_key": artifact_record(CONTINUATION_PUBLIC_KEY),
                "continuation_private_key": {
                    **continuation_private,
                    "state": "UNUSED_NOT_RETIRED_NO_SIGNATURE_CREATED",
                },
            },
            "formal_output_state": {
                "original_v14_authority_review_v_slots_absent_after_archive": True,
                "r_receipt_created": False,
                "r_log_created": False,
                "r_success_witness_created": False,
                "continuation_receipt_created": False,
                "continuation_exit_attestation_created": False,
                "continuation_signature_created": False,
                "continuation_success_witness_created": False,
                "continuation_incident_created": False,
                "canonical_downstream_outputs_created": False,
                "final_dataset_release_receipt_created": False,
                "final_report_release_receipt_created": False,
                "singleton_supplemental_release_receipt_created": False,
            },
            "scientific_payload_changed": False,
            "claim_ceiling_changed": False,
            "required_resolution": (
                "Create append-only v15 with fresh authority and continuation keys, exact "
                "metadata-plus-size private-key relation handling, an actual-path regression, "
                "fresh independent reviews, and distinct V/R/C output slots."
            ),
            "pass": False,
        }
        evidence_data = json.dumps(
            evidence,
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        ).encode("ascii") + b"\n"
        write_exclusive(staging / "failure_evidence.json", evidence_data)
        staging_fd = os.open(staging, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        try:
            os.fsync(staging_fd)
        finally:
            os.close(staging_fd)
        os.rename(staging, FAILURE_ROOT)
        parent_fd = os.open(FAILURE_ROOT.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        try:
            os.fsync(parent_fd)
        finally:
            os.close(parent_fd)
    except BaseException:
        raise

    original_generation_slots = (
        AUTHORITY_BUNDLE,
        VERIFICATION_RECEIPT,
        *REVIEWS.values(),
        *UNUSED_COMPATIBILITY_OUTPUTS,
        *R_OUTPUTS,
        *C_OUTPUTS,
    )
    require_absent(original_generation_slots, "v14 original generation slots after archive")
    evidence_record = artifact_record(FAILURE_ROOT / "failure_evidence.json")
    print(
        json.dumps(
            {
                "archive": str(FAILURE_ROOT),
                "diagnostic_replay_exit_code": completed.returncode,
                "worker_wait_status": wait_status,
                "failure_evidence": evidence_record,
                "v14_original_slots_absent": True,
                "pass": True,
            },
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(
            json.dumps(
                {
                    "error": {"type": type(error).__name__, "message": str(error)},
                    "pass": False,
                },
                ensure_ascii=True,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)
