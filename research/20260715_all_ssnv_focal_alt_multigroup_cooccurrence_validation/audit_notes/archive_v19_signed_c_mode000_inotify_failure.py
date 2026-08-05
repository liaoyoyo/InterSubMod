#!/usr/bin/env python3
"""Archive signed v19 after C could not watch a retired mode-000 key."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
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
    / "20260723_v19_signed_authority_c_mode000_inotify_permission_denied"
)

PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")
STRACE = Path("/usr/bin/strace")
CONTINUATION = AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v19.py"
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
FORMAL_COMMAND = [str(PYTHON), "-I", "-B", str(CONTINUATION), "--supervise-and-sign"]
DIAGNOSTIC_COMMAND = [
    str(STRACE),
    "-f",
    "-qq",
    "-s",
    "20000",
    "-e",
    "trace=write",
    *FORMAL_COMMAND,
]

AUTHORITY_BUNDLE = RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v19.bundle"
VERIFICATION_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v19.json"
)
R_OUTPUTS = {
    "replay_receipt": RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v19.json",
    "replay_log": RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v19.log",
    "replay_success_witness": (
        RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v19.success_witness.json"
    ),
}
REVIEWS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v19.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v19.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v19.json"
    ),
}
SOURCES = {
    "authority_validator": AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v19.py",
    "continuation_verifier": AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v19.py",
    "runner_gate_replay": AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v19.py",
    "downstream_continuation": CONTINUATION,
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v19.py",
    "regression_tests": AUDIT_ROOT / "schema_recovery_tests/test_tumor_ref_schema_recovery_v19.py",
    "ceremony_builder": AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v19.py",
    "archive_script": Path(__file__).resolve(strict=True),
}
REVIEW_TRANSPORT = {
    "external_prompt": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v19_round2_prompt.md"
    ),
    "external_schema": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v19_round2.schema.json"
    ),
    "external_runner": (
        AUDIT_ROOT / "run_external_claude_schema_recovery_review_v19_round2.py"
    ),
    "external_envelope": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v19_round2.claude_cli.envelope.json"
    ),
    "external_stderr": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v19_round2.claude_cli.stderr.txt"
    ),
}

AUTHORITY_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v19/ed25519_public.pem"
)
AUTHORITY_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v19/ed25519_private_one_time.pem"
)
CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v9/ed25519_public.pem"
)
CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v9/ed25519_private_one_time_resident.pem"
)

C_OUTPUTS = (
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v19.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v19.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v19.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v19.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v19.json",
)
CANONICAL_DOWNSTREAM_OUTPUTS = (
    WORKSPACE_ROOT / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "all_ssnv_final_report_dataset_v5_m2v5_source_attested",
    WORKSPACE_ROOT / "all_ssnv_final_report_v5_m2v5_source_attested",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json.ed25519.sig",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json.ed25519.sig",
    RESULT_ROOT / "positional_singleton_supplemental_release_receipt.v1.json",
    RESULT_ROOT / "positional_singleton_supplemental_release_receipt.v1.json.ed25519.sig",
)
RETIRED_RELEASE_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
    "20260717_all_ssnv_v3/"
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
CHILD_ERROR_MESSAGE = (
    "inotify file watch failed for "
    + str(RETIRED_RELEASE_PRIVATE_KEY)
    + ": Permission denied"
)
CHILD_STDERR_SHA256 = (
    "47b76169652069e2cccdcae810a8b1cf9f4eb3bc34df804fde03f0d3030d0c33"
)


class ArchiveError(RuntimeError):
    pass


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def artifact_record(path: Path) -> dict[str, Any]:
    before = os.lstat(path)
    if (
        not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or path.resolve(strict=True) != path
    ):
        raise ArchiveError(f"Artifact is not one canonical regular file: {path}")
    data = path.read_bytes()
    after = os.lstat(path)
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    ) != (
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
    if (
        not stat.S_ISREG(leaf.st_mode)
        or leaf.st_nlink != 1
        or path.resolve(strict=True) != path
    ):
        raise ArchiveError(f"Metadata path is not one canonical regular file: {path}")
    return {
        "path": str(path),
        "size_bytes": leaf.st_size,
        "mode": oct(stat.S_IMODE(leaf.st_mode)),
        "link_count": leaf.st_nlink,
    }


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ArchiveError(f"JSON root is not an object: {path}")
    return value


def write_exclusive(path: Path, data: bytes, mode: int = 0o444) -> None:
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
        0o600,
    )
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
    present = [str(path) for path in paths if os.path.lexists(path)]
    if present:
        raise ArchiveError(f"{label} unexpectedly present: {present}")


def move_no_replace(source: Path, destination: Path) -> None:
    if not os.path.lexists(source) or source.is_symlink():
        raise ArchiveError(f"Archive source missing or invalid: {source}")
    if os.path.lexists(destination):
        raise ArchiveError(f"Archive destination already exists: {destination}")
    os.rename(source, destination)


def copy_no_replace(source: Path, destination: Path) -> None:
    before = os.lstat(source)
    if (
        not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or source.resolve(strict=True) != source
        or os.path.lexists(destination)
    ):
        raise ArchiveError(f"Archive copy source or destination is invalid: {source}")
    data = source.read_bytes()
    after = os.lstat(source)
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    ):
        raise ArchiveError(f"Archive copy source changed while reading: {source}")
    write_exclusive(destination, data)


def record_at_final_path(staging_path: Path, final_path: Path) -> dict[str, Any]:
    record = artifact_record(staging_path)
    record["path"] = str(final_path)
    return record


def symlink_record(alias: Path, target: Path) -> dict[str, Any]:
    before = os.lstat(alias)
    link_text = os.readlink(alias)
    after = os.lstat(alias)
    if (
        not stat.S_ISLNK(before.st_mode)
        or before.st_nlink != 1
        or before.st_dev != after.st_dev
        or before.st_ino != after.st_ino
        or before.st_mtime_ns != after.st_mtime_ns
        or link_text != target.name
        or alias.resolve(strict=True) != target
    ):
        raise ArchiveError(f"Tumor-REF summary alias contract drift: {alias}")
    return {
        "path": str(alias),
        "link_text": link_text,
        "resolved_target": str(target),
        "mode": oct(stat.S_IMODE(after.st_mode)),
        "link_count": after.st_nlink,
    }


def main() -> int:
    if os.path.lexists(FAILURE_ROOT):
        raise ArchiveError(f"Failure archive already exists: {FAILURE_ROOT}")
    require_absent(C_OUTPUTS, "v19 C outputs before diagnostic replay")
    require_absent(CANONICAL_DOWNSTREAM_OUTPUTS, "canonical downstream outputs")
    if not AUTHORITY_BUNDLE.is_dir() or not VERIFICATION_RECEIPT.is_file():
        raise ArchiveError("v19 signed authority or V receipt is missing")
    if any(not path.is_file() for path in R_OUTPUTS.values()):
        raise ArchiveError("One or more v19 R outputs are missing")
    if any(not path.is_file() for path in REVIEWS.values()):
        raise ArchiveError("One or more v19 formal reviews are missing")
    retired_release_private = metadata_record(RETIRED_RELEASE_PRIVATE_KEY)
    if (
        retired_release_private["mode"] != "0o0"
        or retired_release_private["size_bytes"] != 290
        or retired_release_private["link_count"] != 1
    ):
        raise ArchiveError("Retired release private-key metadata drift")

    completed = subprocess.run(
        FORMAL_COMMAND,
        cwd=REPO_ROOT,
        env=FORMAL_ENV,
        capture_output=True,
        check=False,
        timeout=300,
    )
    if completed.returncode != 1 or completed.stdout != b"":
        raise ArchiveError("Formal v19 C diagnostic did not fail closed as expected")
    failure = json.loads(completed.stderr.decode("ascii"))
    expected_parent_message = (
        "Supervised continuation child did not exit zero: "
        f"returncode=1 stderr_sha256={CHILD_STDERR_SHA256}"
    )
    if (
        failure.get("error")
        != {
            "message": expected_parent_message,
            "type": "ContinuationError",
        }
        or failure.get("pass") is not False
    ):
        raise ArchiveError("Formal v19 C diagnostic error contract drift")
    require_absent(C_OUTPUTS, "v19 C outputs after diagnostic replay")
    require_absent(CANONICAL_DOWNSTREAM_OUTPUTS, "canonical downstream outputs")

    diagnostic = subprocess.run(
        DIAGNOSTIC_COMMAND,
        cwd=REPO_ROOT,
        env=FORMAL_ENV,
        capture_output=True,
        check=False,
        timeout=300,
    )
    if (
        diagnostic.returncode != 1
        or diagnostic.stdout != b""
        or diagnostic.stderr.count(CHILD_ERROR_MESSAGE.encode("ascii")) != 1
        or CHILD_STDERR_SHA256.encode("ascii") not in diagnostic.stderr
    ):
        raise ArchiveError("v19 strace diagnostic did not capture the exact child error")
    require_absent(C_OUTPUTS, "v19 C outputs after strace diagnostic replay")
    require_absent(CANONICAL_DOWNSTREAM_OUTPUTS, "canonical downstream outputs")

    authority_payload = load_json(AUTHORITY_BUNDLE / "authority.json")
    verification_payload = load_json(VERIFICATION_RECEIPT)
    replay_payload = load_json(R_OUTPUTS["replay_receipt"])
    witness_payload = load_json(R_OUTPUTS["replay_success_witness"])
    review_payloads = {role: load_json(path) for role, path in REVIEWS.items()}
    if (
        authority_payload.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v19"
        or authority_payload.get("pass") is not True
        or verification_payload.get("pass") is not True
        or replay_payload.get("pass") is not True
        or witness_payload.get("pass") is not True
        or witness_payload.get("worker_wait", {}).get("exit_code") != 0
        or witness_payload.get("worker_wait", {}).get("normal_exit") is not True
    ):
        raise ArchiveError("v19 authority/V/R evidence drift")
    if any(
        payload.get("verdict") != "APPROVE" or payload.get("pass") is not True
        for payload in review_payloads.values()
    ):
        raise ArchiveError("v19 formal review verdict drift")

    authority_private = metadata_record(AUTHORITY_PRIVATE_KEY)
    continuation_private = metadata_record(CONTINUATION_PRIVATE_KEY)
    if authority_private["mode"] != "0o0":
        raise ArchiveError("v19 authority private key was not retired")
    if continuation_private["mode"] != "0o400":
        raise ArchiveError("v19 continuation private key is not unused mode 0400")

    staging = FAILURE_ROOT.parent / f".{FAILURE_ROOT.name}.staging.{os.getpid()}"
    staging.mkdir(mode=0o755)
    move_no_replace(AUTHORITY_BUNDLE, staging / AUTHORITY_BUNDLE.name)
    move_no_replace(VERIFICATION_RECEIPT, staging / VERIFICATION_RECEIPT.name)
    for source in R_OUTPUTS.values():
        move_no_replace(source, staging / source.name)
    for source in REVIEWS.values():
        move_no_replace(source, staging / source.name)

    stderr_path = staging / "formal_continuation_stderr.log"
    strace_path = staging / "diagnostic_strace_write_syscalls.log"
    write_exclusive(stderr_path, completed.stderr)
    write_exclusive(strace_path, diagnostic.stderr)
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
        "formal_continuation_stderr": record_at_final_path(
            stderr_path,
            FAILURE_ROOT / stderr_path.name,
        ),
        "diagnostic_strace_write_syscalls": record_at_final_path(
            strace_path,
            FAILURE_ROOT / strace_path.name,
        ),
    }
    for role, source in R_OUTPUTS.items():
        archived[role] = record_at_final_path(
            staging / source.name,
            FAILURE_ROOT / source.name,
        )
    for role, source in REVIEWS.items():
        archived[f"review_{role}"] = record_at_final_path(
            staging / source.name,
            FAILURE_ROOT / source.name,
        )

    reviewed = next(iter(review_payloads.values()))
    evidence = {
        "schema_name": "intersubmod.tumor_ref_schema_recovery_signed_generation_failure",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "generation": "v19",
        "status": (
            "SIGNED_AUTHORITY_V_AND_R_PASS_C_CHILD_FAILED_AFTER_FRESH_V_"
            "BEFORE_DOWNSTREAM"
        ),
        "review_contract": {
            key: reviewed[key]
            for key in (
                "reviewed_source_set_sha256",
                "legacy_source_set_sha256",
                "prior_recovery_chain_sha256",
                "rejected_generations_sha256",
                "scope_sha256",
                "terminal_key_rotation_sha256",
            )
        },
        "runtime": {
            "v_exit_code": 0,
            "r_supervisor_exit_code": 0,
            "r_worker_exit_code": 0,
            "c_supervisor_exit_code": completed.returncode,
            "c_child_started": True,
            "c_child_fresh_verifier_passed": True,
            "formal_command": FORMAL_COMMAND,
            "clean_environment": FORMAL_ENV,
            "stdout": {
                "size_bytes": len(completed.stdout),
                "sha256": sha256_bytes(completed.stdout),
            },
            "stderr": archived["formal_continuation_stderr"],
            "error_type": "ContinuationError",
            "error_message": expected_parent_message,
            "child_stderr_sha256": CHILD_STDERR_SHA256,
            "child_error_message": CHILD_ERROR_MESSAGE,
            "failure_location": (
                "downstream_continuation.run.PathMutationSentinel.__init__."
                "collect_protected_paths"
            ),
        },
        "diagnostic": {
            "retired_release_private_key": retired_release_private,
            "strace_command": DIAGNOSTIC_COMMAND,
            "strace_exit_code": diagnostic.returncode,
            "strace_stderr": archived["diagnostic_strace_write_syscalls"],
            "syscall_trace_captured_exact_child_error": True,
            "root_cause": (
                "v19 correctly preserved the already-retired release-authority "
                "private key at mode 000, but PathMutationSentinel attempted to add "
                "a direct inotify file watch to every collected protected path. "
                "Linux rejected the direct watch with EACCES after fresh V and before "
                "the downstream leaf started. The key bytes were never read and no "
                "scientific or final-release artifact was created."
            ),
            "required_fix": (
                "Keep the retired key mode 000 and metadata/O_PATH bound. For exactly "
                "mode-000 regular files where inotify_add_watch returns EACCES, use "
                "the already-required named parent-directory watch as the event "
                "source and retain terminal identity rechecks. Record the parent-only "
                "fallback count/hash in the signed mutation-monitor contract and add "
                "real regressions for replacement, chmod, and non-mode-000 failures."
            ),
            "failure_before_c_child_start": False,
            "failure_after_fresh_verifier_pass": True,
            "failure_before_downstream_leaf_start": True,
            "failure_before_first_c_output": True,
            "scientific_payload_changed": False,
        },
        "archived_artifacts": archived,
        "source_artifacts": {
            role: artifact_record(path) for role, path in SOURCES.items()
        },
        "review_transport_artifacts": {
            role: artifact_record(path) for role, path in REVIEW_TRANSPORT.items()
        },
        "reviews": {
            role: {
                "reviewer": payload["reviewer"],
                "reviewer_agent_id": payload["reviewer_agent_id"],
                "verdict": payload["verdict"],
                "high_findings": len(payload.get("high_findings", [])),
                "medium_findings": len(payload.get("medium_findings", [])),
                "unresolved_conditions": len(payload.get("unresolved_conditions", [])),
            }
            for role, payload in review_payloads.items()
        },
        "keys": {
            "authority_public_key": artifact_record(AUTHORITY_PUBLIC_KEY),
            "authority_private_key": {
                **authority_private,
                "state": "RETIRED_AFTER_SINGLE_AUTHORITY_SIGNATURE",
            },
            "continuation_public_key": artifact_record(CONTINUATION_PUBLIC_KEY),
            "continuation_private_key": {
                **continuation_private,
                "state": "UNUSED_NOT_RETIRED_NO_SIGNATURE_CREATED_MUST_NOT_BE_REUSED",
            },
        },
        "formal_output_state": {
            "authority_created": True,
            "verification_receipt_created": True,
            "replay_receipt_created": True,
            "replay_log_created": True,
            "replay_success_witness_created": True,
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
            "Create append-only v20 with fresh authority and terminal-v10 keys, "
            "strict mode-000 parent-watch fallback, live regressions, fresh "
            "independent reviews, and distinct V/R/C output slots."
        ),
        "pass": False,
    }
    write_exclusive(
        staging / "failure_evidence.json",
        json.dumps(evidence, ensure_ascii=True, indent=2, sort_keys=True).encode(
            "ascii"
        )
        + b"\n",
    )
    staging_fd = os.open(staging, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    try:
        os.fsync(staging_fd)
    finally:
        os.close(staging_fd)
    os.rename(staging, FAILURE_ROOT)
    parent_fd = os.open(
        FAILURE_ROOT.parent,
        os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC,
    )
    try:
        os.fsync(parent_fd)
    finally:
        os.close(parent_fd)

    original_slots = (
        AUTHORITY_BUNDLE,
        VERIFICATION_RECEIPT,
        *R_OUTPUTS.values(),
        *REVIEWS.values(),
        *C_OUTPUTS,
    )
    require_absent(original_slots, "v19 original generation slots after archive")
    require_absent(CANONICAL_DOWNSTREAM_OUTPUTS, "canonical downstream outputs")
    print(
        json.dumps(
            {
                "archive": str(FAILURE_ROOT),
                "diagnostic_continuation_exit_code": completed.returncode,
                "failure_evidence": artifact_record(
                    FAILURE_ROOT / "failure_evidence.json"
                ),
                "v19_original_slots_absent": True,
                "canonical_downstream_outputs_absent": True,
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
                    "error": {
                        "type": type(error).__name__,
                        "message": str(error),
                    },
                    "pass": False,
                },
                ensure_ascii=True,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)
