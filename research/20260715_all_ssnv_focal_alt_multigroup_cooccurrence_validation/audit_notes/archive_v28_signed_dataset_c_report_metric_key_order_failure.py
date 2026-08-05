#!/usr/bin/env python3
"""Archive v28 after the signed dataset passed but report schema validation failed."""

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
    / "20260724_v28_signed_dataset_c_report_metric_key_order_mismatch"
)

AUTHORITY_BUNDLE = RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v28.bundle"
RESULT_ARTIFACTS = {
    "verification_receipt": (
        RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v28.json"
    ),
    "runner_replay_receipt": RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.json",
    "runner_replay_log": RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.log",
    "runner_replay_success_witness": (
        RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.success_witness.json"
    ),
    "continuation_incident": (
        RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v28.json"
    ),
}
REVIEWS = {
    "mendel": REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v28.json",
    "nash": REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_nash.v28.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_external_claude_opus.v28.json"
    ),
}
REVIEW_TRANSPORT = {
    "mendel_transport": (
        AUDIT_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v28.multi_agent.transport.json"
    ),
    "nash_transport": (
        AUDIT_ROOT / "20260724_tumor_ref_schema_recovery_nash.v28.multi_agent.transport.json"
    ),
    "external_prompt": AUDIT_ROOT / "20260724_external_claude_schema_recovery_review_v28_prompt.md",
    "external_schema": AUDIT_ROOT / "20260724_external_claude_schema_recovery_review_v28.schema.json",
    "external_runner": AUDIT_ROOT / "run_external_claude_schema_recovery_review_v28.py",
    "external_envelope": (
        AUDIT_ROOT / "20260724_external_claude_schema_recovery_review_v28.claude_cli.envelope.json"
    ),
    "external_stderr": (
        AUDIT_ROOT / "20260724_external_claude_schema_recovery_review_v28.claude_cli.stderr.txt"
    ),
    "review_publisher": AUDIT_ROOT / "publish_tumor_ref_schema_recovery_reviews_v28.py",
}
SOURCES = {
    "authority_validator": AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v28.py",
    "ceremony_builder": AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v28.py",
    "continuation_verifier": AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v28.py",
    "downstream_continuation": AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v28.py",
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v28.py",
    "runner_gate_replay": AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v28.py",
    "regression_tests": AUDIT_ROOT / "schema_recovery_tests/test_tumor_ref_schema_recovery_v28.py",
    "final_dataset_builder": AUDIT_ROOT / "build_all_ssnv_final_report_dataset_schema_recovery_v28.py",
    "result_finalizer": AUDIT_ROOT / "finalize_task_b_result_release_schema_recovery_v28.py",
    "report_builder": AUDIT_ROOT / "build_all_ssnv_report_artifact_schema_recovery_v28.py",
    "archive_script": Path(__file__).resolve(strict=True),
}
DOWNSTREAM_OUTPUTS = {
    "strict_not_applicable": (
        WORKSPACE_ROOT / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5"
    ),
    "matched_normal_not_applicable": (
        WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5"
    ),
    "cn_ccf_not_applicable": (
        WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5"
    ),
    "primary_post_audit": (
        RESULT_ROOT / "stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json"
    ),
    "frozen_post_audit": (
        RESULT_ROOT / "frozen_input_immutability.post_m2v5_downstream_v3_source_authority_v5.json"
    ),
    "final_dataset": (
        WORKSPACE_ROOT / "all_ssnv_final_report_dataset_v5_m2v5_source_attested"
    ),
    "dataset_release_receipt": RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json",
    "dataset_release_signature": (
        RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json.ed25519.sig"
    ),
    "completion_log": (
        WORKSPACE_ROOT / "m2v5_recovered_completion_chain_v2_source_authority_v5.log"
    ),
    "completion_cache": (
        WORKSPACE_ROOT / ".python_cache_m2v5_completion_v2_bound_bootstrap"
    ),
}

FINAL_OUTPUTS_THAT_MUST_BE_ABSENT = (
    WORKSPACE_ROOT / "all_ssnv_final_report_v5_m2v5_source_attested",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v28.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v28.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v28.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v28.json",
)

KEY_SPECS = {
    "authority_v28": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260724_v28"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/archive/"
            "20260724_v28_signed_dataset_c_report_metric_key_order_failure_01"
        ),
        "private_glob": "ed25519_private*.pem",
        "expected_mode": "0o0",
        "state": "RETIRED_AFTER_SINGLE_SIGNED_AUTHORITY_THAT_FAILED_AT_C",
    },
    "terminal_v18": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260724_m2v5_terminal_v18"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
            "20260724_m2v5_terminal_v18_unused_after_signed_v28_c_failure_01"
        ),
        "private_glob": "ed25519_private*.pem",
        "expected_mode": "0o400",
        "state": "UNUSED_NO_TERMINAL_SIGNATURE_RETIRED_FROM_REUSE",
    },
    "result_v5": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
            "20260719_all_ssnv_result_v5_post_reboot_bootstrap"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/archive/"
            "20260724_all_ssnv_result_v5_consumed_v28_c_failure_01"
        ),
        "private_glob": "ed25519_private*.pem",
        "expected_mode": "0o0",
        "state": "CONSUMED_DATASET_SIGNATURE_PROVISIONAL_AFTER_C_FAILURE",
    },
    "report_v5": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
            "20260719_all_ssnv_report_v5_post_reboot_bootstrap"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/archive/"
            "20260724_all_ssnv_report_v5_unused_v28_c_failure_01"
        ),
        "private_glob": "ed25519_private*.pem",
        "expected_mode": "0o400",
        "state": "UNUSED_REPORT_KEY_RETIRED_FROM_FAILED_GENERATION",
    },
}

EXPECTED_ERROR = "ReportContractError: Final per_sample metric strata drift"
OPENSSL = Path("/usr/bin/openssl")


class ArchiveError(RuntimeError):
    """Raised before any unsupported or overwrite-prone archive operation."""


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ArchiveError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_nonfinite(value: str) -> None:
    raise ArchiveError(f"Non-finite JSON constant: {value}")


def strict_json(path: Path) -> dict[str, Any]:
    value = json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=reject_duplicate_keys,
        parse_constant=reject_nonfinite,
    )
    if not isinstance(value, dict):
        raise ArchiveError(f"JSON root is not an object: {path}")
    return value


def file_record(path: Path) -> dict[str, Any]:
    before = os.lstat(path)
    if (
        not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or path.resolve(strict=True) != path
    ):
        raise ArchiveError(f"Not one canonical regular file: {path}")
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
        raise ArchiveError(f"File changed while hashing: {path}")
    return {
        "mode": oct(stat.S_IMODE(after.st_mode)),
        "path": str(path),
        "sha256": sha256_bytes(data),
        "size_bytes": len(data),
    }


def metadata_record(path: Path) -> dict[str, Any]:
    observed = os.lstat(path)
    if not stat.S_ISREG(observed.st_mode) or observed.st_nlink != 1:
        raise ArchiveError(f"Invalid metadata-only file: {path}")
    return {
        "link_count": observed.st_nlink,
        "mode": oct(stat.S_IMODE(observed.st_mode)),
        "path": str(path),
        "size_bytes": observed.st_size,
    }


def write_exclusive(path: Path, data: bytes, mode: int = 0o444) -> None:
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
        0o600,
    )
    try:
        view = memoryview(data)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise ArchiveError(f"Short archive write: {path}")
            view = view[written:]
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
        raise ArchiveError(f"Archive destination exists: {destination}")
    os.rename(source, destination)


def harden_tree(root: Path) -> None:
    for current_root, directories, filenames in os.walk(root):
        current = Path(current_root)
        for directory in directories:
            os.chmod(current / directory, 0o755)
        for filename in filenames:
            os.chmod(current / filename, 0o444)


def tree_inventory(root: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ArchiveError(f"Archive tree contains symlink: {path}")
        if path.is_file():
            record = file_record(path)
            record["relative_path"] = str(path.relative_to(root))
            records.append(record)
    return records


def verify_dataset_signature() -> dict[str, Any]:
    receipt = DOWNSTREAM_OUTPUTS["dataset_release_receipt"]
    signature = DOWNSTREAM_OUTPUTS["dataset_release_signature"]
    public_key = KEY_SPECS["result_v5"]["source"] / "ed25519_public.pem"
    result = subprocess.run(
        [
            str(OPENSSL),
            "pkeyutl",
            "-verify",
            "-rawin",
            "-pubin",
            "-inkey",
            str(public_key),
            "-sigfile",
            str(signature),
            "-in",
            str(receipt),
        ],
        check=False,
        capture_output=True,
        env={"PATH": "/usr/bin:/bin", "HOME": "/tmp", "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8"},
    )
    if result.returncode != 0:
        raise ArchiveError("v28 provisional dataset signature did not verify")
    return {
        "public_key": file_record(public_key),
        "receipt": file_record(receipt),
        "signature": file_record(signature),
        "openssl_exit_code": result.returncode,
        "verified": True,
    }


def write_key_archive_record(role: str, spec: dict[str, Any]) -> dict[str, Any]:
    source = spec["source"]
    destination = spec["destination"]
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not source.is_dir() or os.path.lexists(destination):
        raise ArchiveError(f"Key archive source/destination invalid: {role}")
    private_paths = sorted(source.glob(spec["private_glob"]))
    public_paths = sorted(source.glob("ed25519_public.pem"))
    if len(private_paths) != 1 or len(public_paths) != 1:
        raise ArchiveError(f"Unexpected key inventory: {role}")
    private_metadata = metadata_record(private_paths[0])
    if private_metadata["mode"] != spec["expected_mode"]:
        raise ArchiveError(f"Private-key state drift: {role}")
    record = {
        "schema_name": "intersubmod.failed_formal_generation_key_archive",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "generation": "v28",
        "role": role,
        "state": spec["state"],
        "source_directory": str(source),
        "archive_directory": str(destination),
        "public_key": file_record(public_paths[0]),
        "private_key_metadata_only": private_metadata,
        "private_key_bytes_read": False,
        "key_reuse_forbidden": True,
        "pass": False,
    }
    write_exclusive(
        source / "FAILED_KEY_ARCHIVE_RECORD.v1.json",
        json.dumps(record, ensure_ascii=True, indent=2, sort_keys=True).encode("ascii") + b"\n",
    )
    move_no_replace(source, destination)
    return {"path": str(destination), "state": spec["state"], "private_key_pre_archive": private_metadata}


def validate_failure_state() -> dict[str, Any]:
    if os.path.lexists(FAILURE_ROOT):
        raise ArchiveError(f"Failure archive exists: {FAILURE_ROOT}")
    require_absent(FINAL_OUTPUTS_THAT_MUST_BE_ABSENT, "v28 terminal/report outputs")
    required = (
        AUTHORITY_BUNDLE,
        *RESULT_ARTIFACTS.values(),
        *REVIEWS.values(),
        *REVIEW_TRANSPORT.values(),
        *SOURCES.values(),
        *DOWNSTREAM_OUTPUTS.values(),
        *(spec["source"] for spec in KEY_SPECS.values()),
    )
    missing = [str(path) for path in required if not os.path.lexists(path)]
    if missing:
        raise ArchiveError(f"Required v28 failure evidence missing: {missing}")

    authority = strict_json(AUTHORITY_BUNDLE / "authority.json")
    verification = strict_json(RESULT_ARTIFACTS["verification_receipt"])
    replay = strict_json(RESULT_ARTIFACTS["runner_replay_receipt"])
    witness = strict_json(RESULT_ARTIFACTS["runner_replay_success_witness"])
    incident = strict_json(RESULT_ARTIFACTS["continuation_incident"])
    dataset_receipt = strict_json(DOWNSTREAM_OUTPUTS["dataset_release_receipt"])
    final_dataset = strict_json(DOWNSTREAM_OUTPUTS["final_dataset"] / "final_report_dataset.json")
    completion_log = DOWNSTREAM_OUTPUTS["completion_log"].read_text(encoding="utf-8", errors="strict")
    if (
        authority.get("authority_id") != "20260724_tumor_ref_schema_recovery_v28"
        or authority.get("pass") is not True
        or verification.get("pass") is not True
        or replay.get("pass") is not True
        or witness.get("pass") is not True
        or incident.get("pass") is not False
        or incident.get("release_authority_granted") is not False
        or incident.get("continuation_receipt_exists") is not False
        or incident.get("continuation_signature_exists") is not False
        or dataset_receipt.get("pass") is not True
        or dataset_receipt.get("scope") != "all_7_datasets_469849_sites_final_dataset"
        or final_dataset.get("pass") is not True
        or final_dataset.get("counts", {}).get("screen_sites") != 469849
        or EXPECTED_ERROR not in completion_log
    ):
        raise ArchiveError("v28 authority/V/R/C failure contract drift")
    reviews = {role: strict_json(path) for role, path in REVIEWS.items()}
    if any(
        payload.get("verdict") != "APPROVE" or payload.get("pass") is not True
        for payload in reviews.values()
    ):
        raise ArchiveError("v28 review contract drift")
    per_sample = final_dataset.get("funnel_metrics", {}).get("per_sample")
    if not isinstance(per_sample, dict) or tuple(per_sample) != tuple(sorted(per_sample)):
        raise ArchiveError("Observed sorted JSON per-sample key order drift")
    return {
        "authority": authority,
        "verification": verification,
        "replay": replay,
        "witness": witness,
        "incident": incident,
        "reviews": reviews,
        "dataset_signature": verify_dataset_signature(),
        "per_sample_keys": list(per_sample),
    }


def main() -> int:
    observed = validate_failure_state()
    source_records = {role: file_record(path) for role, path in SOURCES.items()}
    original_modes = {
        role: oct(stat.S_IMODE(os.lstat(path).st_mode))
        for role, path in DOWNSTREAM_OUTPUTS.items()
    }
    staging = FAILURE_ROOT.parent / f".{FAILURE_ROOT.name}.staging.{os.getpid()}"
    staging.mkdir(mode=0o755)
    (staging / "review_transport").mkdir(mode=0o755)
    (staging / "observed_downstream_outputs").mkdir(mode=0o755)

    move_no_replace(AUTHORITY_BUNDLE, staging / AUTHORITY_BUNDLE.name)
    for path in RESULT_ARTIFACTS.values():
        move_no_replace(path, staging / path.name)
    for path in REVIEWS.values():
        move_no_replace(path, staging / path.name)
    for path in REVIEW_TRANSPORT.values():
        move_no_replace(path, staging / "review_transport" / path.name)
    for path in DOWNSTREAM_OUTPUTS.values():
        move_no_replace(path, staging / "observed_downstream_outputs" / path.name)

    key_archives = {
        role: write_key_archive_record(role, spec)
        for role, spec in KEY_SPECS.items()
    }
    harden_tree(staging)
    evidence = {
        "schema_name": "intersubmod.tumor_ref_schema_recovery_signed_generation_failure",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "generation": "v28",
        "status": "SIGNED_AUTHORITY_V_R_PASS_C_FAILED_AFTER_SIGNED_DATASET_BEFORE_REPORT",
        "task_type": "B_comprehensive_validation",
        "scope": "all_7_datasets_469849_sites",
        "root_cause": {
            "error_type": "ReportContractError",
            "error_fragment": EXPECTED_ERROR,
            "failure_location": "build_all_ssnv_report_artifact.validate_metrics",
            "primary_contract_failure": (
                "The report builder treated JSON object insertion order as schema for "
                "per_sample, truth_strata, and sample_by_truth. The final dataset is "
                "canonically serialized with sorted keys, so the exact content passed "
                "but the order-only predicate failed."
            ),
            "scientific_payload_changed": False,
            "claim_ceiling_changed": False,
        },
        "verified_pre_failure_state": {
            "authority_pass": observed["authority"]["pass"],
            "verification_pass": observed["verification"]["pass"],
            "runner_replay_pass": observed["replay"]["pass"],
            "runner_success_witness_pass": observed["witness"]["pass"],
            "continuation_incident_pass": observed["incident"]["pass"],
            "dataset_signature": observed["dataset_signature"],
            "per_sample_keys_as_serialized": observed["per_sample_keys"],
            "reviews": {
                role: {"verdict": payload["verdict"], "pass": payload["pass"]}
                for role, payload in observed["reviews"].items()
            },
        },
        "original_downstream_modes": original_modes,
        "archived_inventory": tree_inventory(staging),
        "source_artifacts_left_immutable_at_original_paths": source_records,
        "key_archives": key_archives,
        "terminal_output_state": {
            "final_dataset_created": True,
            "dataset_release_receipt_created": True,
            "dataset_release_signature_created": True,
            "dataset_release_signature_independently_verified": True,
            "final_report_created": False,
            "report_release_receipt_created": False,
            "report_release_signature_created": False,
            "continuation_receipt_created": False,
            "continuation_signature_created": False,
            "continuation_success_witness_created": False,
        },
        "required_resolution": (
            "Use append-only v29 with exact-key-set validation for JSON mappings, add a "
            "canonical sorted-key round-trip regression, and use fresh authority, terminal, "
            "result-release, and report-release keys."
        ),
        "pass": False,
    }
    write_exclusive(
        staging / "failure_evidence.json",
        json.dumps(evidence, ensure_ascii=True, indent=2, sort_keys=True).encode("ascii") + b"\n",
    )
    summary = f"""<!--
建立時間: 2026-07-24
目標: 封存 v28 已簽 dataset 後的 report metric key-order schema 失敗
處理範圍: v28 authority/reviews/V/R/C incident/provisional dataset/key quarantine
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/{FAILURE_ROOT.name}/failure_evidence.json
-->

# v28 Signed Dataset C Failure Archive

- v28 authority、V、R 與 final dataset builder 均 PASS；dataset receipt 的 Ed25519 signature 亦獨立驗證成功。
- C 在 report builder `validate_metrics` fail-closed，尚未建立 final report 或 terminal continuation receipt。
- 根因是把 JSON object key insertion order 當作 schema；canonical JSON 會排序 keys，內容未改變。
- provisional dataset、receipt、signature 與所有本輪 downstream outputs 已完整封存，不提供 release authority。
- authority、terminal、result-release、report-release keys 全部封存並禁止重用。
- v29 必須改用 exact key-set 驗證並增加 canonical round-trip regression。
"""
    write_exclusive(staging / "SUMMARY.md", summary.encode("utf-8"))
    harden_tree(staging)
    os.rename(staging, FAILURE_ROOT)
    parent_fd = os.open(FAILURE_ROOT.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    try:
        os.fsync(parent_fd)
    finally:
        os.close(parent_fd)

    original_slots = (
        AUTHORITY_BUNDLE,
        *RESULT_ARTIFACTS.values(),
        *REVIEWS.values(),
        *REVIEW_TRANSPORT.values(),
        *DOWNSTREAM_OUTPUTS.values(),
        *(spec["source"] for spec in KEY_SPECS.values()),
    )
    require_absent(original_slots, "v28 active slots after archive")
    require_absent(FINAL_OUTPUTS_THAT_MUST_BE_ABSENT, "v28 report/terminal slots after archive")
    result = {
        "archive": str(FAILURE_ROOT),
        "failure_evidence": file_record(FAILURE_ROOT / "failure_evidence.json"),
        "summary": file_record(FAILURE_ROOT / "SUMMARY.md"),
        "key_archives": key_archives,
        "v28_active_slots_absent": True,
        "report_terminal_slots_absent": True,
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
                {"error": {"type": type(error).__name__, "message": str(error)}, "pass": False},
                ensure_ascii=True,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)
