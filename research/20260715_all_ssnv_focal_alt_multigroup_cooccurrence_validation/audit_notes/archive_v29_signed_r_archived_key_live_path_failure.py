#!/usr/bin/env python3
"""Archive signed v29 after its fail-closed historical runner-key path error."""

from __future__ import annotations

from datetime import datetime, timezone
import ctypes
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
from typing import Any, Mapping


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
    / "20260724_v29_signed_authority_v_pass_r_archived_key_live_path_mismatch"
)

PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")
BASH = Path("/usr/bin/bash")
REPLAYER = AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v29.py"
RUNNER = TOPIC_ROOT / "scripts/run_m2v5_recovered_completion_chain.sh"
RUNNER_PREFIX_LINE_COUNT = 358
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

AUTHORITY_BUNDLE = RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v29.bundle"
VERIFICATION_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v29.json"
)
REVIEWS = {
    "mendel": REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v29.json",
    "nash": REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_nash.v29.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_external_claude_opus.v29.json"
    ),
}
REVIEW_TRANSPORT = {
    "mendel": AUDIT_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v29.multi_agent.transport.json",
    "nash": AUDIT_ROOT / "20260724_tumor_ref_schema_recovery_nash.v29.multi_agent.transport.json",
    "external_envelope": (
        AUDIT_ROOT / "20260724_external_claude_schema_recovery_review_v29.claude_cli.envelope.json"
    ),
    "external_stderr": (
        AUDIT_ROOT / "20260724_external_claude_schema_recovery_review_v29.claude_cli.stderr.txt"
    ),
    "external_schema": (
        AUDIT_ROOT / "20260724_external_claude_schema_recovery_review_v29.schema.json"
    ),
    "external_prompt": (
        AUDIT_ROOT / "20260724_external_claude_schema_recovery_review_v29_prompt.md"
    ),
}
SOURCES = {
    "authority_validator": AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v29.py",
    "continuation_verifier": AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v29.py",
    "runner_gate_replay": REPLAYER,
    "downstream_continuation": AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v29.py",
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v29.py",
    "regression_tests": AUDIT_ROOT / "schema_recovery_tests/test_tumor_ref_schema_recovery_v29.py",
    "ceremony_builder": AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v29.py",
    "dataset_builder": AUDIT_ROOT / "build_all_ssnv_final_report_dataset_schema_recovery_v29.py",
    "result_finalizer": AUDIT_ROOT / "finalize_task_b_result_release_schema_recovery_v29.py",
    "report_builder": AUDIT_ROOT / "build_all_ssnv_report_artifact_schema_recovery_v29.py",
    "review_publisher": AUDIT_ROOT / "publish_tumor_ref_schema_recovery_reviews_v29.py",
    "external_review_runner": AUDIT_ROOT / "run_external_claude_schema_recovery_review_v29.py",
    "archive_script": Path(__file__).resolve(strict=True),
}

R_OUTPUTS = (
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v29.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v29.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v29.success_witness.json",
)
C_OUTPUTS = (
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v29.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v29.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v29.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v29.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v29.json",
)
FINAL_OUTPUTS = (
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
)

EXPECTED_CLOSED_SIGNER_PIDS = (3446130, 3446180)
AT_FDCWD = -100
RENAME_NOREPLACE = 1

KEY_SPECS: dict[str, dict[str, Any]] = {
    "authority_v29": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260724_v29"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/archive/"
            "20260724_v29_signed_authority_v_pass_r_archived_key_live_path_failure_01"
        ),
        "expected_mode": "0o0",
        "state": "RETIRED_AFTER_SIGNED_V29_AUTHORITY_R_RUNTIME_FAILURE",
    },
    "terminal_v19": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260724_m2v5_terminal_v19"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
            "20260724_m2v5_terminal_v19_unused_after_signed_v29_r_failure_01"
        ),
        "expected_mode": "0o400",
        "state": "UNUSED_NO_TERMINAL_SIGNATURE_RETIRED_AFTER_V29_R_FAILURE",
    },
    "result_v6": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
            "20260724_all_ssnv_result_v6_v29_recovery"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/archive/"
            "20260724_all_ssnv_result_v6_unused_v29_r_failure_01"
        ),
        "expected_mode": "0o400",
        "state": "UNUSED_NO_SIGNATURE_KEY_UNHELD_RETIRED_AFTER_V29_R_FAILURE",
    },
    "report_v6": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
            "20260724_all_ssnv_report_v6_v29_recovery"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/archive/"
            "20260724_all_ssnv_report_v6_unused_v29_r_failure_01"
        ),
        "expected_mode": "0o400",
        "state": "UNUSED_NO_SIGNATURE_KEY_UNHELD_RETIRED_AFTER_V29_R_FAILURE",
    },
}

EXPECTED_PREFIX_STDERR = (
    "Required file is missing or empty: "
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
    "20260719_all_ssnv_result_v5_post_reboot_bootstrap/ed25519_public.pem\n"
).encode("ascii")
EXPECTED_PREFIX_STDERR_SHA256 = "18477fe5a82535b402f6f04e72ebcea147f5c5cee0475dba6e250ae75078c4a0"


class ArchiveError(RuntimeError):
    """Raised before any overwrite-prone or semantically invalid archive action."""


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def strict_json(path: Path) -> dict[str, Any]:
    def reject_duplicate(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ArchiveError(f"Duplicate JSON key in {path}: {key}")
            result[key] = value
        return result

    value = json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=reject_duplicate,
        parse_constant=lambda value: (_ for _ in ()).throw(
            ArchiveError(f"Non-finite JSON constant in {path}: {value}")
        ),
    )
    if not isinstance(value, dict):
        raise ArchiveError(f"JSON root is not an object: {path}")
    return value


def file_record(path: Path) -> dict[str, Any]:
    before = os.lstat(path)
    if not stat.S_ISREG(before.st_mode) or path.resolve(strict=True) != path:
        raise ArchiveError(f"Not one canonical regular file: {path}")
    data = path.read_bytes()
    after = os.lstat(path)
    identity = lambda item: (
        item.st_dev,
        item.st_ino,
        item.st_mode,
        item.st_size,
        item.st_mtime_ns,
        item.st_ctime_ns,
        item.st_nlink,
    )
    if identity(before) != identity(after) or after.st_nlink != 1:
        raise ArchiveError(f"Artifact identity changed or link count is not one: {path}")
    return {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": sha256_bytes(data),
        "mode": oct(stat.S_IMODE(after.st_mode)),
        "link_count": after.st_nlink,
    }


def metadata_record(path: Path) -> dict[str, Any]:
    observed = os.lstat(path)
    if not stat.S_ISREG(observed.st_mode) or path.resolve(strict=True) != path:
        raise ArchiveError(f"Not one canonical private-key file: {path}")
    return {
        "path": str(path),
        "size_bytes": observed.st_size,
        "mode": oct(stat.S_IMODE(observed.st_mode)),
        "link_count": observed.st_nlink,
    }


def write_exclusive(path: Path, data: bytes, mode: int = 0o444) -> None:
    descriptor = os.open(
        path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC, 0o600
    )
    try:
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise ArchiveError(f"Short write: {path}")
            offset += written
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
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise ArchiveError("renameat2 is required for no-replace archival")
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    result = renameat2(
        AT_FDCWD,
        os.fsencode(source),
        AT_FDCWD,
        os.fsencode(destination),
        RENAME_NOREPLACE,
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise ArchiveError(
            f"No-replace archive rename failed: {source} -> {destination}: "
            f"{os.strerror(error_number)}"
        )


def harden_tree(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        if path.is_dir():
            os.chmod(path, 0o755)
        elif path.is_file():
            os.chmod(path, 0o444)
        else:
            raise ArchiveError(f"Unsupported archive entry: {path}")
    os.chmod(root, 0o755)


def final_record(staging_path: Path, relative_path: Path) -> dict[str, Any]:
    record = file_record(staging_path)
    record["path"] = str(FAILURE_ROOT / relative_path)
    return record


def active_processes_for_key_roots() -> list[dict[str, Any]]:
    needles = [str(spec["source"]).encode("utf-8") for spec in KEY_SPECS.values()]
    active: list[dict[str, Any]] = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit() or int(proc.name) == os.getpid():
            continue
        try:
            command = (proc / "cmdline").read_bytes()
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        if any(needle in command for needle in needles):
            active.append(
                {
                    "pid": int(proc.name),
                    "cmdline_sha256": sha256_bytes(command),
                }
            )
    return active


def active_key_descriptor_holders() -> dict[str, list[dict[str, Any]]]:
    protected_identities: dict[tuple[int, int], str] = {}
    for spec in KEY_SPECS.values():
        source = Path(spec["source"])
        destination = Path(spec["destination"])
        root = source if source.is_dir() else destination
        if not root.is_dir():
            raise ArchiveError(
                f"Neither live nor archived key root exists during holder scan: {source}"
            )
        for path in (root, *sorted(root.iterdir())):
            observed = os.stat(path, follow_symlinks=False)
            protected_identities[(observed.st_dev, observed.st_ino)] = str(path)

    active: list[dict[str, Any]] = []
    uninspectable_same_uid: list[dict[str, Any]] = []
    current_uid = os.getuid()
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit() or int(proc.name) == os.getpid():
            continue
        try:
            process_stat = proc.stat()
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        if process_stat.st_uid != current_uid:
            continue
        try:
            descriptors = list((proc / "fd").iterdir())
        except (FileNotFoundError, ProcessLookupError):
            continue
        except PermissionError:
            try:
                command = (proc / "cmdline").read_bytes()
            except (FileNotFoundError, PermissionError, ProcessLookupError):
                command = b""
            uninspectable_same_uid.append(
                {
                    "pid": int(proc.name),
                    "cmdline_sha256": sha256_bytes(command),
                    "descriptor_scope": "UNINSPECTABLE_NON_DUMPABLE_PROCESS",
                }
            )
            continue
        for descriptor in descriptors:
            try:
                observed = os.stat(descriptor)
            except (FileNotFoundError, ProcessLookupError):
                continue
            except PermissionError as error:
                raise ArchiveError(
                    f"Cannot stat descriptor {descriptor}: {error}"
                ) from error
            identity = (observed.st_dev, observed.st_ino)
            if identity in protected_identities:
                active.append(
                    {
                        "pid": int(proc.name),
                        "fd": descriptor.name,
                        "protected_path": protected_identities[identity],
                    }
                )
    return {
        "active_protected_descriptors": active,
        "uninspectable_same_uid_processes": uninspectable_same_uid,
    }


def require_signers_closed_and_keys_unheld(stage: str) -> dict[str, Any]:
    reused = [pid for pid in EXPECTED_CLOSED_SIGNER_PIDS if Path(f"/proc/{pid}").exists()]
    if reused:
        raise ArchiveError(f"Expected closed v29 signer PID still exists at {stage}: {reused}")
    command_matches = active_processes_for_key_roots()
    descriptor_scan = active_key_descriptor_holders()
    descriptor_holders = descriptor_scan["active_protected_descriptors"]
    if command_matches or descriptor_holders:
        raise ArchiveError(
            f"v29 key process/descriptor remains active at {stage}: "
            f"cmdline={command_matches} descriptors={descriptor_holders}"
        )
    return {
        "stage": stage,
        "expected_closed_signer_pids": list(EXPECTED_CLOSED_SIGNER_PIDS),
        "expected_signer_pids_absent": True,
        "key_root_cmdline_matches": command_matches,
        **descriptor_scan,
        "trusted_boundary": (
            "linux_kernel_and_uncompromised_same_uid_research_account; non-dumpable "
            "same-UID descriptor tables may be uninspectable and are recorded, not "
            "claimed inspected; malicious same-UID processes remain out of scope"
        ),
        "pass": True,
    }


def complete_downstream_output_slots() -> tuple[Path, ...]:
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("v29_replayer_archive_contract", REPLAYER)
    if spec is None or spec.loader is None:
        raise ArchiveError("Unable to load frozen v29 replayer output-slot contract")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    slots = getattr(module, "DOWNSTREAM_OUTPUT_SLOTS", None)
    if not isinstance(slots, tuple) or not slots:
        raise ArchiveError("Frozen v29 replayer downstream output-slot contract drift")
    result = tuple(Path(path) for path in slots)
    if len(result) != len(set(result)):
        raise ArchiveError("Frozen v29 replayer downstream output slots contain duplicates")
    return result


def archive_key(role: str, spec: Mapping[str, Any]) -> dict[str, Any]:
    source = Path(spec["source"])
    destination = Path(spec["destination"])
    if destination.parent.resolve(strict=True) != destination.parent:
        raise ArchiveError(f"Non-canonical key archive parent: {destination.parent}")
    if not source.is_dir() or os.path.lexists(destination):
        raise ArchiveError(f"Invalid key archive source/destination: {role}")
    public_paths = sorted(source.glob("ed25519_public.pem"))
    private_paths = sorted(source.glob("ed25519_private*.pem"))
    if len(public_paths) != 1 or len(private_paths) != 1:
        raise ArchiveError(f"Unexpected key inventory: {role}")
    public_record = file_record(public_paths[0])
    private_record = metadata_record(private_paths[0])
    if (
        public_record["mode"] != "0o444"
        or private_record["mode"] != spec["expected_mode"]
        or private_record["link_count"] != 1
    ):
        raise ArchiveError(f"Key mode/link state drift: {role}")
    record = {
        "schema_name": "intersubmod.failed_formal_generation_key_archive",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "generation": "v29",
        "role": role,
        "state": spec["state"],
        "source_directory": str(source),
        "archive_directory": str(destination),
        "public_key": public_record,
        "private_key_metadata_only": private_record,
        "private_key_bytes_read": False,
        "key_reuse_forbidden": True,
        "pass": False,
    }
    write_exclusive(
        source / "FAILED_KEY_ARCHIVE_RECORD.v1.json",
        json.dumps(record, ensure_ascii=True, indent=2, sort_keys=True).encode("ascii")
        + b"\n",
    )
    move_no_replace(source, destination)
    return {
        "path": str(destination),
        "state": spec["state"],
        "public_key": {
            **public_record,
            "path": str(destination / "ed25519_public.pem"),
        },
        "private_key_pre_archive": private_record,
    }


def runner_prefix() -> bytes:
    data = RUNNER.read_bytes()
    lines = data.splitlines(keepends=True)
    if len(lines) < 633 or lines[358] not in (b"\n", b"\r\n"):
        raise ArchiveError("Runner physical-line boundary drift")
    return b"".join(lines[:RUNNER_PREFIX_LINE_COUNT])


def reproduce_failure() -> tuple[subprocess.CompletedProcess[bytes], subprocess.CompletedProcess[bytes]]:
    formal = subprocess.run(
        FORMAL_COMMAND,
        cwd=REPO_ROOT,
        env=FORMAL_ENV,
        capture_output=True,
        check=False,
        timeout=900,
    )
    formal_stderr = formal.stderr.decode("utf-8", errors="strict")
    if (
        formal.returncode != 1
        or formal.stdout != b""
        or EXPECTED_PREFIX_STDERR_SHA256 not in formal_stderr
        or "Runner-only preflight replay failed: exit=1" not in formal_stderr
    ):
        raise ArchiveError("Formal v29 R failure was not reproduced exactly")
    wait_match = re.search(
        r"Runner-gate worker did not produce waitpid normal exit zero: "
        r"pid=(\d+) waited_pid=(\d+) status=(\d+)",
        formal_stderr,
    )
    if (
        wait_match is None
        or wait_match.group(1) != wait_match.group(2)
        or int(wait_match.group(3)) != 17920
    ):
        raise ArchiveError("Formal v29 waitpid evidence drift")

    prefix = subprocess.run(
        [str(BASH), "-s"],
        input=runner_prefix(),
        cwd=REPO_ROOT,
        env=FORMAL_ENV,
        capture_output=True,
        check=False,
        timeout=900,
    )
    if (
        prefix.returncode != 1
        or prefix.stdout != b""
        or prefix.stderr != EXPECTED_PREFIX_STDERR
        or sha256_bytes(prefix.stderr) != EXPECTED_PREFIX_STDERR_SHA256
    ):
        raise ArchiveError("Exact runner-prefix missing-key failure drift")
    return formal, prefix


def main() -> int:
    if os.path.lexists(FAILURE_ROOT):
        raise ArchiveError(f"Failure archive already exists: {FAILURE_ROOT}")
    all_absent = tuple(
        dict.fromkeys(
            (
                *R_OUTPUTS,
                *C_OUTPUTS,
                *FINAL_OUTPUTS,
                *complete_downstream_output_slots(),
            )
        )
    )
    require_absent(all_absent, "v29 downstream output")
    required = (
        AUTHORITY_BUNDLE,
        VERIFICATION_RECEIPT,
        *REVIEWS.values(),
        *REVIEW_TRANSPORT.values(),
        *SOURCES.values(),
        *(Path(spec["source"]) for spec in KEY_SPECS.values()),
    )
    missing = [str(path) for path in required if not os.path.lexists(path)]
    if missing:
        raise ArchiveError(f"Required v29 failure evidence missing: {missing}")
    signer_holder_checks = [
        require_signers_closed_and_keys_unheld("before diagnostic replay")
    ]

    authority = strict_json(AUTHORITY_BUNDLE / "authority.json")
    verification = strict_json(VERIFICATION_RECEIPT)
    reviews = {role: strict_json(path) for role, path in REVIEWS.items()}
    if (
        authority.get("authority_id") != "20260724_tumor_ref_schema_recovery_v29"
        or authority.get("pass") is not True
        or verification.get("pass") is not True
        or any(
            payload.get("verdict") != "APPROVE" or payload.get("pass") is not True
            for payload in reviews.values()
        )
    ):
        raise ArchiveError("v29 authority/V/review contract drift")

    formal, prefix = reproduce_failure()
    require_absent(all_absent, "v29 downstream output after diagnostic replay")
    signer_holder_checks.append(
        require_signers_closed_and_keys_unheld("after diagnostic replay")
    )

    staging = FAILURE_ROOT.parent / f".{FAILURE_ROOT.name}.staging.{os.getpid()}"
    staging.mkdir(mode=0o755)
    (staging / "review_transport").mkdir(mode=0o755)
    move_no_replace(AUTHORITY_BUNDLE, staging / AUTHORITY_BUNDLE.name)
    move_no_replace(VERIFICATION_RECEIPT, staging / VERIFICATION_RECEIPT.name)
    for path in REVIEWS.values():
        move_no_replace(path, staging / path.name)
    for path in REVIEW_TRANSPORT.values():
        move_no_replace(path, staging / "review_transport" / path.name)
    write_exclusive(staging / "formal_replay_stderr.log", formal.stderr)
    write_exclusive(staging / "runner_prefix_stderr.log", prefix.stderr)

    signer_holder_checks.append(
        require_signers_closed_and_keys_unheld("immediately before key archival")
    )
    key_archives = {}
    for role, spec in KEY_SPECS.items():
        signer_holder_checks.append(
            require_signers_closed_and_keys_unheld(f"before archiving {role}")
        )
        key_archives[role] = archive_key(role, spec)
    harden_tree(staging)

    archived_artifacts = {
        "authority": final_record(
            staging / AUTHORITY_BUNDLE.name / "authority.json",
            Path(AUTHORITY_BUNDLE.name) / "authority.json",
        ),
        "authority_signature": final_record(
            staging / AUTHORITY_BUNDLE.name / "authority.ed25519.sig",
            Path(AUTHORITY_BUNDLE.name) / "authority.ed25519.sig",
        ),
        "authority_commit": final_record(
            staging / AUTHORITY_BUNDLE.name / "commit.json",
            Path(AUTHORITY_BUNDLE.name) / "commit.json",
        ),
        "verification_receipt": final_record(
            staging / VERIFICATION_RECEIPT.name,
            Path(VERIFICATION_RECEIPT.name),
        ),
        "formal_replay_stderr": final_record(
            staging / "formal_replay_stderr.log", Path("formal_replay_stderr.log")
        ),
        "runner_prefix_stderr": final_record(
            staging / "runner_prefix_stderr.log", Path("runner_prefix_stderr.log")
        ),
    }
    for role, path in REVIEWS.items():
        archived_artifacts[f"review_{role}"] = final_record(
            staging / path.name, Path(path.name)
        )
    for role, path in REVIEW_TRANSPORT.items():
        archived_artifacts[f"review_transport_{role}"] = final_record(
            staging / "review_transport" / path.name,
            Path("review_transport") / path.name,
        )

    review_contract = next(iter(reviews.values()))
    evidence = {
        "schema_name": "intersubmod.tumor_ref_schema_recovery_signed_generation_failure",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "generation": "v29",
        "status": "SIGNED_AUTHORITY_AND_V_PASS_R_FAILED_BEFORE_FIRST_R_OUTPUT",
        "task_type": "B_comprehensive_validation",
        "scope": "all_7_datasets_469849_sites",
        "review_contract": {
            key: review_contract[key]
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
            "r_supervisor_exit_code": formal.returncode,
            "r_worker_exit_code": 70,
            "r_waitpid_raw_status": 17920,
            "r_worker_normal_exit": True,
            "r_worker_success": False,
            "c_started": False,
            "formal_command": FORMAL_COMMAND,
            "clean_environment": FORMAL_ENV,
            "formal_stdout": {
                "size_bytes": len(formal.stdout),
                "sha256": sha256_bytes(formal.stdout),
            },
            "formal_stderr": archived_artifacts["formal_replay_stderr"],
            "runner_prefix_stderr": archived_artifacts["runner_prefix_stderr"],
            "runner_prefix_stderr_sha256": EXPECTED_PREFIX_STDERR_SHA256,
            "failure_location": "runner_gate_replay.worker_main.runner_prefix_subprocess",
            "error_type": "GateError",
        },
        "root_cause": {
            "primary_contract_failure": (
                "The v29 runner-only gate executed historical completion-runner physical lines "
                "1-358 unchanged. Those lines require the failed-v28 result/report keys at their "
                "former live paths, while v29 correctly quarantined those keys in append-only "
                "archives and bound fresh v6 result/report keys only in downstream continuation."
            ),
            "probe_gap": (
                "The v29 read-only probe checked replayer source self-binding and static helper "
                "contracts but did not execute the real runner prefix under the archived-key live "
                "state."
            ),
            "scientific_payload_changed": False,
            "claim_ceiling_changed": False,
        },
        "rejected_remediation": {
            "temporary_live_key_projection": False,
            "reason": (
                "The failed-v28 archive records set key_reuse_forbidden=true; copying the still-"
                "usable report-v5 private key back to its former live path would weaken quarantine "
                "semantics even if no signature were requested."
            ),
        },
        "required_resolution": (
            "Create append-only v30 with fresh authority, terminal, result, and report keys; "
            "perform an authority-bound deterministic archive-path rebase of only the three "
            "historical result/report key constants before runner-prefix execution; record both "
            "prefix digests and the exact mapping; and add a real archived-key-state regression."
        ),
        "archived_artifacts": archived_artifacts,
        "key_archives": key_archives,
        "signer_and_key_holder_checks": signer_holder_checks,
        "source_artifacts_left_immutable_at_original_paths": {
            role: file_record(path) for role, path in SOURCES.items()
        },
        "reviews": {
            role: {
                "reviewer": payload.get("reviewer"),
                "reviewer_agent_id": payload.get("reviewer_agent_id"),
                "verdict": payload.get("verdict"),
                "pass": payload.get("pass"),
            }
            for role, payload in reviews.items()
        },
        "terminal_output_state": {
            "r_log_created": False,
            "r_receipt_created": False,
            "r_success_witness_created": False,
            "c_started": False,
            "final_dataset_created": False,
            "final_report_created": False,
            "dataset_release_signature_created": False,
            "report_release_signature_created": False,
            "terminal_signature_created": False,
        },
        "scientific_payload_changed": False,
        "claim_ceiling_changed": False,
        "pass": False,
    }
    write_exclusive(
        staging / "failure_evidence.json",
        json.dumps(evidence, ensure_ascii=True, indent=2, sort_keys=True).encode("ascii")
        + b"\n",
    )
    summary = f"""<!--
建立時間: 2026-07-24
目標: 封存 v29 已簽 authority 與 V 通過後的 historical runner-key live-path 失敗
處理範圍: v29 authority/reviews/V/R failure/key quarantine；科學資料未重算、未改變
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/{FAILURE_ROOT.name}/failure_evidence.json
-->

# v29 Signed Authority R Failure Archive

- v29 authority 與 promotion verification（32/32）PASS；正式 R 在第一個 replay output 前 fail-closed。
- exact runner 第 1-358 行仍要求已封存 v28 result/report keys 的舊 live 路徑；第一個缺檔為 result-v5 public key。
- v29 probe/test 驗證了 replayer source 與 helper contract，但未在真實 archived-key live state 執行 runner prefix。
- 未建立 R receipt/log/witness、未啟動 C、未建立 dataset/report 或任何 v29 terminal/result/report signature。
- 不採用舊私鑰 live projection；authority、terminal-v19、result-v6、report-v6 keys 全數封存且禁止重用。
- v30 必須使用 authority-bound archive-path rebase、真實 archived-state regression 與四組全新 keys。
"""
    write_exclusive(staging / "SUMMARY.md", summary.encode("utf-8"))
    harden_tree(staging)
    move_no_replace(staging, FAILURE_ROOT)
    parent_fd = os.open(
        FAILURE_ROOT.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
    )
    try:
        os.fsync(parent_fd)
    finally:
        os.close(parent_fd)

    original_slots = (
        AUTHORITY_BUNDLE,
        VERIFICATION_RECEIPT,
        *REVIEWS.values(),
        *REVIEW_TRANSPORT.values(),
        *R_OUTPUTS,
        *C_OUTPUTS,
        *FINAL_OUTPUTS,
        *(Path(spec["source"]) for spec in KEY_SPECS.values()),
    )
    require_absent(original_slots, "v29 active slot after archive")
    result = {
        "archive": str(FAILURE_ROOT),
        "failure_evidence": file_record(FAILURE_ROOT / "failure_evidence.json"),
        "summary": file_record(FAILURE_ROOT / "SUMMARY.md"),
        "formal_replay_exit_code": formal.returncode,
        "runner_prefix_exit_code": prefix.returncode,
        "runner_prefix_stderr_sha256": sha256_bytes(prefix.stderr),
        "key_archives": key_archives,
        "v29_active_slots_absent": True,
        "pass": True,
    }
    print(json.dumps(result, ensure_ascii=True, indent=2, sort_keys=True))
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
