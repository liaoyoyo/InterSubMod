#!/usr/bin/env python3
"""Archive signed v25 after the tumor-REF v1/v6 audit-path relation failure."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import stat
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
    / "20260724_v25_signed_authority_c_tumor_ref_pre_audit_path_relation_mismatch"
)

AUTHORITY_BUNDLE = RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v25.bundle"
RESULT_ARTIFACTS = {
    "verification_receipt": (
        RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v25.json"
    ),
    "runner_replay_receipt": (
        RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v25.json"
    ),
    "runner_replay_log": (
        RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v25.log"
    ),
    "runner_replay_success_witness": (
        RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v25.success_witness.json"
    ),
    "continuation_incident": (
        RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v25.json"
    ),
}
REVIEWS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v25.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v25.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v25.json"
    ),
}
REVIEW_TRANSPORT = {
    "mendel_transport": (
        AUDIT_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v25.multi_agent.transport.json"
    ),
    "nash_transport": (
        AUDIT_ROOT / "20260723_tumor_ref_schema_recovery_nash.v25.multi_agent.transport.json"
    ),
    "external_prompt": (
        AUDIT_ROOT / "20260723_external_claude_schema_recovery_review_v25_prompt.md"
    ),
    "external_schema": (
        AUDIT_ROOT / "20260723_external_claude_schema_recovery_review_v25.schema.json"
    ),
    "external_runner": (
        AUDIT_ROOT / "run_external_claude_schema_recovery_review_v25.py"
    ),
    "external_envelope": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v25.claude_cli.envelope.json"
    ),
    "external_stderr": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v25.claude_cli.stderr.txt"
    ),
    "review_publisher": (
        AUDIT_ROOT / "publish_tumor_ref_schema_recovery_reviews_v25.py"
    ),
}
SOURCES = {
    "authority_validator": AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v25.py",
    "ceremony_builder": AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v25.py",
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v25.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v25.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v25.py",
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v25.py"
    ),
    "regression_tests": (
        AUDIT_ROOT / "schema_recovery_tests/test_tumor_ref_schema_recovery_v25.py"
    ),
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
    "completion_log": (
        WORKSPACE_ROOT / "m2v5_recovered_completion_chain_v2_source_authority_v5.log"
    ),
    "completion_cache": (
        WORKSPACE_ROOT / ".python_cache_m2v5_completion_v2_bound_bootstrap"
    ),
}
GENERATED_CACHE = {
    "validator_python39_cache": (
        AUDIT_ROOT / "__pycache__/validate_tumor_ref_schema_recovery_authority_v25.cpython-39.pyc"
    ),
}

FINAL_OUTPUTS_THAT_MUST_BE_ABSENT = (
    WORKSPACE_ROOT / "all_ssnv_final_report_dataset_v5_m2v5_source_attested",
    WORKSPACE_ROOT / "all_ssnv_final_report_v5_m2v5_source_attested",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json.ed25519.sig",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json.ed25519.sig",
    RESULT_ROOT / "positional_singleton_supplemental_release_receipt.v1.json",
    RESULT_ROOT / "positional_singleton_supplemental_release_receipt.v1.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v25.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v25.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v25.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v25.json",
)

AUTHORITY_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260723_v25"
)
AUTHORITY_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/archive/"
    "20260723_v25_signed_c_tumor_ref_pre_audit_path_relation_failure_01"
)
TERMINAL_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v15"
)
TERMINAL_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
    "20260723_m2v5_terminal_v15_unused_after_signed_v25_c_failure_01"
)

LEGACY_AUDIT = RESULT_ROOT / "stable_primary_artifact_audit.v1_pre_downstream.json"
CURRENT_AUDIT = (
    RESULT_ROOT / "stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream.json"
)
TUMOR_REF_MANIFEST = (
    WORKSPACE_ROOT / "all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel/run_manifest.json"
)
EXPECTED_ERROR_FRAGMENT = (
    "tumor-REF pre-downstream primary artifact audit path mismatch: "
    f"declared='{LEGACY_AUDIT}' observed='{CURRENT_AUDIT}'"
)


class ArchiveError(RuntimeError):
    """Raised before any overwrite or unsupported archive state."""


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def strict_json(path: Path) -> dict[str, Any]:
    value = json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_keys,
        parse_constant=_reject_nonfinite,
    )
    if not isinstance(value, dict):
        raise ArchiveError(f"JSON root is not an object: {path}")
    return value


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ArchiveError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_nonfinite(value: str) -> None:
    raise ArchiveError(f"Non-finite JSON constant: {value}")


def file_record(path: Path, *, include_mode: bool = True) -> dict[str, Any]:
    before = os.lstat(path)
    if (
        not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or path.resolve(strict=True) != path
    ):
        raise ArchiveError(f"Not one canonical regular file: {path}")
    data = path.read_bytes()
    after = os.lstat(path)
    stable = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    ) == (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    )
    if not stable:
        raise ArchiveError(f"File changed while hashing: {path}")
    record: dict[str, Any] = {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": sha256_bytes(data),
    }
    if include_mode:
        record["mode"] = oct(stat.S_IMODE(after.st_mode))
    return record


def metadata_record(path: Path) -> dict[str, Any]:
    observed = os.lstat(path)
    if not stat.S_ISREG(observed.st_mode) or observed.st_nlink != 1:
        raise ArchiveError(f"Invalid metadata-only file: {path}")
    return {
        "path": str(path),
        "size_bytes": observed.st_size,
        "mode": oct(stat.S_IMODE(observed.st_mode)),
        "link_count": observed.st_nlink,
    }


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


def harden_tree(root: Path) -> None:
    for current_root, directories, files in os.walk(root):
        current = Path(current_root)
        for directory in directories:
            os.chmod(current / directory, 0o755)
        for filename in files:
            path = current / filename
            if path.name.startswith("ed25519_private"):
                continue
            os.chmod(path, 0o444)


def tree_inventory(root: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ArchiveError(f"Archive tree contains a symlink: {path}")
        if path.is_file():
            record = file_record(path)
            record["relative_path"] = str(path.relative_to(root))
            records.append(record)
    return records


def write_key_archive_record(
    source: Path,
    destination: Path,
    *,
    generation: str,
    state: str,
) -> None:
    if not source.is_dir() or os.path.lexists(destination):
        raise ArchiveError(f"Key archive source/destination invalid: {source}")
    private_paths = sorted(source.glob("ed25519_private*.pem"))
    public_paths = sorted(source.glob("ed25519_public.pem"))
    if len(private_paths) != 1 or len(public_paths) != 1:
        raise ArchiveError(f"Unexpected key inventory: {source}")
    record = {
        "schema_name": "intersubmod.failed_formal_generation_key_archive",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "generation": generation,
        "state": state,
        "source_directory": str(source),
        "archive_directory": str(destination),
        "public_key": file_record(public_paths[0]),
        "private_key_metadata_only": metadata_record(private_paths[0]),
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


def validate_failure_state() -> dict[str, Any]:
    if os.path.lexists(FAILURE_ROOT):
        raise ArchiveError(f"Failure archive already exists: {FAILURE_ROOT}")
    require_absent(FINAL_OUTPUTS_THAT_MUST_BE_ABSENT, "terminal/final outputs")
    required = (
        AUTHORITY_BUNDLE,
        *RESULT_ARTIFACTS.values(),
        *REVIEWS.values(),
        *REVIEW_TRANSPORT.values(),
        *SOURCES.values(),
        *DOWNSTREAM_OUTPUTS.values(),
        AUTHORITY_KEY_ROOT,
        TERMINAL_KEY_ROOT,
        LEGACY_AUDIT,
        CURRENT_AUDIT,
        TUMOR_REF_MANIFEST,
    )
    missing = [str(path) for path in required if not os.path.lexists(path)]
    if missing:
        raise ArchiveError(f"Required v25 failure evidence missing: {missing}")

    authority = strict_json(AUTHORITY_BUNDLE / "authority.json")
    verification = strict_json(RESULT_ARTIFACTS["verification_receipt"])
    replay = strict_json(RESULT_ARTIFACTS["runner_replay_receipt"])
    witness = strict_json(RESULT_ARTIFACTS["runner_replay_success_witness"])
    incident = strict_json(RESULT_ARTIFACTS["continuation_incident"])
    reviews = {role: strict_json(path) for role, path in REVIEWS.items()}
    completion_log = DOWNSTREAM_OUTPUTS["completion_log"].read_text(
        encoding="utf-8", errors="strict"
    )
    if (
        authority.get("authority_id") != "20260723_tumor_ref_schema_recovery_v25"
        or authority.get("pass") is not True
        or verification.get("pass") is not True
        or replay.get("pass") is not True
        or witness.get("pass") is not True
        or incident.get("pass") is not False
        or incident.get("release_authority_granted") is not False
        or incident.get("continuation_receipt_exists") is not False
        or incident.get("continuation_signature_exists") is not False
        or EXPECTED_ERROR_FRAGMENT not in completion_log
    ):
        raise ArchiveError("v25 authority/V/R/C failure contract drift")
    if any(
        payload.get("verdict") != "APPROVE" or payload.get("pass") is not True
        for payload in reviews.values()
    ):
        raise ArchiveError("v25 review contract drift")

    legacy = strict_json(LEGACY_AUDIT)
    current = strict_json(CURRENT_AUDIT)
    manifest = strict_json(TUMOR_REF_MANIFEST)
    declared = manifest.get("inputs", {}).get("primary_artifact_audit_pre", {})
    if (
        declared.get("path") != str(LEGACY_AUDIT)
        or declared.get("sha256") != file_record(LEGACY_AUDIT)["sha256"]
        or legacy.get("verification", {}).get("artifact_set_sha256")
        != current.get("verification", {}).get("artifact_set_sha256")
        or legacy.get("counts") != current.get("counts")
        or legacy.get("inputs", {}).get("site_results")
        != {
            key: value
            for key, value in current.get("inputs", {}).get("site_results", {}).items()
            if key in legacy.get("inputs", {}).get("site_results", {})
        }
        or legacy.get("inputs", {}).get("stable_assignments")
        != {
            key: value
            for key, value in current.get("inputs", {}).get("stable_assignments", {}).items()
            if key in legacy.get("inputs", {}).get("stable_assignments", {})
        }
    ):
        raise ArchiveError("Observed v1/v6 audit lineage evidence drift")
    return {
        "authority": authority,
        "verification": verification,
        "replay": replay,
        "witness": witness,
        "incident": incident,
        "reviews": reviews,
        "legacy_audit": legacy,
        "current_audit": current,
    }


def main() -> int:
    observed = validate_failure_state()
    original_modes = {
        role: (
            oct(stat.S_IMODE(os.lstat(path).st_mode))
            if path.is_file()
            else oct(stat.S_IMODE(os.lstat(path).st_mode))
        )
        for role, path in DOWNSTREAM_OUTPUTS.items()
    }
    source_records = {role: file_record(path) for role, path in SOURCES.items()}
    authority_private = metadata_record(
        AUTHORITY_KEY_ROOT / "ed25519_private_one_time.pem"
    )
    terminal_private = metadata_record(
        TERMINAL_KEY_ROOT / "ed25519_private_one_time_resident.pem"
    )
    if authority_private["mode"] != "0o0" or terminal_private["mode"] != "0o400":
        raise ArchiveError("v25 authority/terminal private-key state drift")

    staging = FAILURE_ROOT.parent / f".{FAILURE_ROOT.name}.staging.{os.getpid()}"
    staging.mkdir(mode=0o755)
    (staging / "review_transport").mkdir(mode=0o755)
    (staging / "observed_downstream_outputs").mkdir(mode=0o755)
    (staging / "generated_cache").mkdir(mode=0o755)

    move_no_replace(AUTHORITY_BUNDLE, staging / AUTHORITY_BUNDLE.name)
    for path in RESULT_ARTIFACTS.values():
        move_no_replace(path, staging / path.name)
    for path in REVIEWS.values():
        move_no_replace(path, staging / path.name)
    for path in REVIEW_TRANSPORT.values():
        move_no_replace(path, staging / "review_transport" / path.name)
    for path in DOWNSTREAM_OUTPUTS.values():
        move_no_replace(path, staging / "observed_downstream_outputs" / path.name)
    for path in GENERATED_CACHE.values():
        if os.path.lexists(path):
            move_no_replace(path, staging / "generated_cache" / path.name)

    write_key_archive_record(
        AUTHORITY_KEY_ROOT,
        AUTHORITY_KEY_ARCHIVE,
        generation="v25",
        state="RETIRED_AFTER_SINGLE_SIGNED_AUTHORITY_THAT_FAILED_AT_C",
    )
    write_key_archive_record(
        TERMINAL_KEY_ROOT,
        TERMINAL_KEY_ARCHIVE,
        generation="terminal_v15",
        state="UNUSED_NO_TERMINAL_SIGNATURE_CREATED_RETIRED_FROM_REUSE",
    )

    harden_tree(staging)
    archived_inventory = tree_inventory(staging)
    evidence = {
        "schema_name": "intersubmod.tumor_ref_schema_recovery_signed_generation_failure",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "generation": "v25",
        "status": "SIGNED_AUTHORITY_V_R_PASS_C_FAILED_AT_FINAL_DATASET_INTEGRATION",
        "task_type": "B_comprehensive_validation",
        "scope": "all_7_datasets_469849_sites",
        "root_cause": {
            "error_type": "ContractError",
            "error_fragment": EXPECTED_ERROR_FRAGMENT,
            "failure_location": "build_all_ssnv_final_report_dataset.load_tumor_ref",
            "primary_contract_failure": (
                "The immutable tumor-REF manifest declares the legacy v1 pre-audit, "
                "while the formal final builder requires the source-authorized v6 pre-audit."
            ),
            "secondary_contract_failure": (
                "The v25 Python continuation declared fresh source_authority_v6 output slots, "
                "but the composed completion-runner prefix still executed canonical v5 slots."
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
            "reviews": {
                role: {
                    "verdict": payload["verdict"],
                    "pass": payload["pass"],
                }
                for role, payload in observed["reviews"].items()
            },
        },
        "audit_lineage_observation": {
            "legacy_v1": file_record(LEGACY_AUDIT),
            "current_v6": file_record(CURRENT_AUDIT),
            "tumor_ref_manifest": file_record(TUMOR_REF_MANIFEST),
            "artifact_set_sha256": observed["legacy_audit"]["verification"][
                "artifact_set_sha256"
            ],
            "counts": observed["legacy_audit"]["counts"],
            "lineage_not_yet_authorized_by_v25": True,
        },
        "original_downstream_modes": original_modes,
        "archived_inventory": archived_inventory,
        "source_artifacts_left_immutable_at_original_paths": source_records,
        "key_archives": {
            "authority": {
                "path": str(AUTHORITY_KEY_ARCHIVE),
                "private_key_pre_archive": authority_private,
                "state": "RETIRED_NO_REUSE",
            },
            "terminal": {
                "path": str(TERMINAL_KEY_ARCHIVE),
                "private_key_pre_archive": terminal_private,
                "state": "UNUSED_BUT_RETIRED_FROM_REUSE",
            },
        },
        "terminal_output_state": {
            "final_dataset_created": False,
            "final_report_created": False,
            "dataset_release_receipt_created": False,
            "dataset_release_signature_created": False,
            "report_release_receipt_created": False,
            "report_release_signature_created": False,
            "supplemental_release_receipt_created": False,
            "continuation_receipt_created": False,
            "continuation_signature_created": False,
            "continuation_success_witness_created": False,
        },
        "required_resolution": (
            "Use append-only v26 with a separately reviewed recovery builder that validates "
            "the exact v1-to-v6 artifact lineage and split chronology; keep canonical v5 leaf "
            "producer paths aligned with the signed v7 source authority; use fresh authority "
            "and terminal keys."
        ),
        "pass": False,
    }
    write_exclusive(
        staging / "failure_evidence.json",
        json.dumps(evidence, ensure_ascii=True, indent=2, sort_keys=True).encode("ascii")
        + b"\n",
    )
    summary = f"""<!--
建立時間: 2026-07-24
目標: 封存 v25 已簽 authority 的 C-stage tumor-REF audit lineage 失敗
處理範圍: v25 authority/reviews/V/R/C incident/downstream observations/key quarantine
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/{FAILURE_ROOT.name}/failure_evidence.json
-->

# v25 Signed C Failure Archive

- v25 authority、V 與 R 均 PASS；C 在 final dataset integration fail-closed。
- 直接原因：tumor-REF manifest 綁定 v1 pre-audit，builder canonical input 為 v6 pre-audit。
- v1/v6 審查同一 `102,842` stable sites、`308,526` artifacts 與同一 artifact-set SHA，
  但 v25 沒有授權這個跨代 lineage relation，因此不可 bypass。
- 次要原因：v25 continuation 的 Python 預期槽與 rendered shell 實際槽不一致。
- 沒有 final dataset/report/supplemental receipt 或 signature；科學資料與 claim ceiling 未改變。
- authority key 已退役；未使用的 terminal-v15 key 也已封存並禁止重用。
"""
    write_exclusive(staging / "SUMMARY.md", summary.encode("utf-8"))
    harden_tree(staging)

    staging_fd = os.open(staging, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    try:
        os.fsync(staging_fd)
    finally:
        os.close(staging_fd)
    os.rename(staging, FAILURE_ROOT)
    parent_fd = os.open(
        FAILURE_ROOT.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
    )
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
        *GENERATED_CACHE.values(),
        AUTHORITY_KEY_ROOT,
        TERMINAL_KEY_ROOT,
    )
    require_absent(original_slots, "v25 active slots after archive")
    require_absent(FINAL_OUTPUTS_THAT_MUST_BE_ABSENT, "terminal/final outputs after archive")
    result = {
        "archive": str(FAILURE_ROOT),
        "failure_evidence": file_record(FAILURE_ROOT / "failure_evidence.json"),
        "summary": file_record(FAILURE_ROOT / "SUMMARY.md"),
        "authority_key_archive": str(AUTHORITY_KEY_ARCHIVE),
        "terminal_key_archive": str(TERMINAL_KEY_ARCHIVE),
        "v25_active_slots_absent": True,
        "terminal_final_slots_absent": True,
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
