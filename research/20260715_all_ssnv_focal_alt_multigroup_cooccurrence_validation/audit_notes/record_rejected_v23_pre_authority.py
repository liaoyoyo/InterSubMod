#!/usr/bin/env python3
"""Record the immutable pre-authority rejection evidence for recovery v23."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path


TOPIC_ROOT = Path(__file__).resolve().parent.parent
RESULT_ROOT = TOPIC_ROOT / "results"
REVIEW_ROOT = TOPIC_ROOT / "reviews"
ARCHIVE_ROOT = (
    TOPIC_ROOT
    / "audit_notes"
    / "rejected_pre_authority_reviews"
    / "20260723_v23_nonexistent_node_runtime_path"
)
AUTHORITY_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260723_v23"
)
TERMINAL_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v13"
)

SOURCE_PATHS = {
    "authority_validator": ARCHIVE_ROOT / "validate_tumor_ref_schema_recovery_authority_v23.py",
    "ceremony_builder": ARCHIVE_ROOT / "build_tumor_ref_schema_recovery_authority_v23.py",
    "continuation_verifier": ARCHIVE_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v23.py",
    "downstream_continuation": ARCHIVE_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v23.py",
    "readonly_probe": ARCHIVE_ROOT / "probe_tumor_ref_schema_recovery_sources_v23.py",
    "regression_tests": ARCHIVE_ROOT / "test_tumor_ref_schema_recovery_v23.py",
    "runner_gate_replay": ARCHIVE_ROOT / "replay_m2v5_runner_only_gates_recovery_v23.py",
}
REVIEW_PATHS = {
    "mendel": ARCHIVE_ROOT / "mendel_request_changes.json",
    "nash": ARCHIVE_ROOT / "nash_request_changes.json",
    "external_claude_opus": (
        ARCHIVE_ROOT
        / "20260723_external_claude_schema_recovery_review_v23.claude_cli.envelope.json"
    ),
}
FORMAL_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v23.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v23.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v23.json.ed25519.sig",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v23.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v23.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v23.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v23.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v23.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v23.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v23.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v23.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v23.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v23.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v23.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v23.json",
)


def artifact(path: Path) -> dict[str, object]:
    data = path.read_bytes()
    mode = path.stat().st_mode & 0o7777
    return {
        "mode": f"0o{mode:o}",
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(data).hexdigest(),
        "size_bytes": len(data),
    }


def create_immutable(path: Path, data: bytes) -> None:
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
    try:
        os.write(fd, data)
        os.fsync(fd)
    finally:
        os.close(fd)
    os.chmod(path, 0o444)


def main() -> int:
    if not ARCHIVE_ROOT.is_dir():
        raise RuntimeError(f"Archive root is absent: {ARCHIVE_ROOT}")
    occupied = [str(path) for path in FORMAL_OUTPUT_SLOTS if os.path.lexists(path)]
    if occupied:
        raise RuntimeError(f"Rejected v23 formal output slots are occupied: {occupied}")
    staging = list(
        RESULT_ROOT.glob(".tumor_ref_promotion_schema_recovery_authority.v23.bundle.staging.*")
    )
    if staging:
        raise RuntimeError(f"Rejected v23 staging paths are occupied: {staging}")

    source_set = {role: artifact(path) for role, path in SOURCE_PATHS.items()}
    reviews = {role: artifact(path) for role, path in REVIEW_PATHS.items()}
    for record in source_set.values():
        if record["mode"] != "0o444":
            raise RuntimeError(f"Rejected source is mutable: {record['path']}")
    for record in reviews.values():
        if record["mode"] != "0o444":
            raise RuntimeError(f"Rejected review is mutable: {record['path']}")

    mendel = json.loads(REVIEW_PATHS["mendel"].read_text(encoding="utf-8"))
    nash = json.loads(REVIEW_PATHS["nash"].read_text(encoding="utf-8"))
    external = json.loads(
        REVIEW_PATHS["external_claude_opus"].read_text(encoding="utf-8")
    )
    external_review = external["structured_output"]
    if (
        mendel.get("verdict") != "REQUEST_CHANGES"
        or nash.get("verdict") != "REQUEST_CHANGES"
        or external_review.get("verdict") != "APPROVE"
    ):
        raise RuntimeError("Strictest-review inputs do not match the captured verdicts")

    authority_private = artifact(AUTHORITY_KEY_ROOT / "ed25519_private_one_time.pem")
    authority_public = artifact(AUTHORITY_KEY_ROOT / "ed25519_public.pem")
    terminal_private = artifact(
        TERMINAL_KEY_ROOT / "ed25519_private_one_time_resident.pem"
    )
    terminal_public = artifact(TERMINAL_KEY_ROOT / "ed25519_public.pem")
    if authority_private["mode"] != "0o400" or terminal_private["mode"] != "0o400":
        raise RuntimeError("Rejected unused private keys are not quarantined at mode 0400")

    authority_private["state"] = "UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED"
    terminal_private["state"] = "UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED"
    evidence = {
        "schema_name": "intersubmod.tumor_ref_schema_recovery_pre_authority_candidate_rejection",
        "schema_version": "1.1.0",
        "candidate_generation": "v23",
        "status": "REJECTED_BEFORE_AUTHORITY_BY_STRICTEST_REVIEW_WINS",
        "successor_generation": "v24",
        "reviewed_source_set_sha256": "7304704c23a5cf73c6a4f83f28345a7999a57eddfe507b350f1e894affab25fd",
        "source_set": source_set,
        "reviews": {
            "mendel": {
                "artifact": reviews["mendel"],
                "reviewer_agent_id": mendel["reviewer_agent_id"],
                "verdict": mendel["verdict"],
                "attribution": "orchestrator-recorded multi-agent transport",
            },
            "nash": {
                "artifact": reviews["nash"],
                "reviewer_agent_id": nash["reviewer_agent_id"],
                "verdict": nash["verdict"],
                "attribution": "orchestrator-recorded multi-agent transport",
            },
            "external_claude_opus": {
                "artifact": reviews["external_claude_opus"],
                "reviewer_agent_id": external_review["reviewer_agent_id"],
                "verdict": external_review["verdict"],
                "attribution": "Claude CLI session envelope",
            },
        },
        "review_attribution_semantics": (
            "Internal reviewer IDs are orchestrator-recorded transport attribution only and "
            "are not cryptographic proof of reviewer authorship. The authority may authenticate "
            "frozen review bytes and the coordinator's decision, but must not claim "
            "reviewer-authenticated signatures."
        ),
        "strictest_review_summary": {
            "approve_count": 1,
            "request_changes_count": 2,
            "blocking_high_count": 2,
            "blocking_medium_count": 1,
            "nonblocking_low_count": len(external_review.get("low_findings", [])),
            "blocking_findings": [
                {
                    "id": "V23-H1",
                    "severity": "high",
                    "summary": (
                        "The mandatory Node runtime path points to absent v23.22.1 while its "
                        "pinned hash matches the existing v22.22.1 executable."
                    ),
                },
                {
                    "id": "V23-M1",
                    "severity": "medium",
                    "summary": (
                        "The probe and regression suite lacked a real existence and SHA-256 "
                        "check for every mandatory executable in GATE_INPUT_PATHS."
                    ),
                },
            ],
        },
        "pre_authority_state": {
            "authority_created": False,
            "authority_signature_created": False,
            "canonical_review_outputs_created": False,
            "verification_receipt_created": False,
            "replay_outputs_created": False,
            "continuation_outputs_created": False,
            "scientific_payload_changed": False,
            "canonical_receipt_bytes_changed": False,
            "claim_ceiling_changed": False,
        },
        "quarantined_unused_keys": {
            "authority": {
                "private_key": authority_private,
                "public_key": authority_public,
            },
            "terminal": {
                "private_key": terminal_private,
                "public_key": terminal_public,
            },
        },
        "formal_output_slots": {
            "expected_count": len(FORMAL_OUTPUT_SLOTS),
            "all_absent": True,
            "paths": [str(path.resolve()) for path in FORMAL_OUTPUT_SLOTS],
        },
        "scientific_claim_ceiling": {
            "confirmed_cellular_subclones": 0,
            "linear_ancestry_calls": 0,
            "permitted_claim": "latent molecular substructure candidates only",
        },
        "pass": True,
    }
    evidence_bytes = (
        json.dumps(evidence, ensure_ascii=False, indent=2, sort_keys=False) + "\n"
    ).encode("utf-8")
    create_immutable(ARCHIVE_ROOT / "rejection_evidence.json", evidence_bytes)

    summary = """<!--
建立時間: 2026-07-23
目標: 封存 recovery v23 在 authority 建立前被獨立審查拒絕的完整證據
處理範圍: 七個凍結來源、三方審查 transport、未使用金鑰與 15 個正式輸出槽位
關聯檔案: rejection_evidence.json
-->

# Recovery v23 pre-authority rejection

v23 依 strictest-review-wins 規則於 authority 建立前拒絕。Mendel 與 Nash
獨立指出 continuation 將必需的 Node runtime 指向不存在的 `v23.22.1`；
固定 SHA-256 實際對應仍存在的 `v22.22.1`。Mendel 另指出 558 項測試與
唯讀 probe 沒有實際驗證所有必需執行檔的存在性與 SHA-256。

外部 Claude Opus 回報 APPROVE，但不凌駕兩份 REQUEST_CHANGES。v23
authority、正式 review 與 V/R/C 輸出均未建立；兩組 v23 私鑰保持未使用並
標記永不重用。科學資料、469,849 位點、read-tag join 與結論上限均未改變。
"""
    create_immutable(ARCHIVE_ROOT / "SUMMARY.md", summary.encode("utf-8"))
    print(
        json.dumps(
            {
                "evidence": artifact(ARCHIVE_ROOT / "rejection_evidence.json"),
                "formal_output_slots_absent": True,
                "pass": True,
                "summary": artifact(ARCHIVE_ROOT / "SUMMARY.md"),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
