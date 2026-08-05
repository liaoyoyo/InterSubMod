#!/usr/bin/env python3
"""Record immutable pre-authority rejection evidence for recovery v24."""

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
    / "20260723_v24_runtime_execution_and_inventory_scope_findings"
)
AUTHORITY_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260723_v24"
)
TERMINAL_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v14"
)

SOURCE_PATHS = {
    "authority_validator": ARCHIVE_ROOT / "validate_tumor_ref_schema_recovery_authority_v24.py",
    "ceremony_builder": ARCHIVE_ROOT / "build_tumor_ref_schema_recovery_authority_v24.py",
    "continuation_verifier": ARCHIVE_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v24.py",
    "downstream_continuation": ARCHIVE_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v24.py",
    "readonly_probe": ARCHIVE_ROOT / "probe_tumor_ref_schema_recovery_sources_v24.py",
    "regression_tests": ARCHIVE_ROOT / "test_tumor_ref_schema_recovery_v24.py",
    "runner_gate_replay": ARCHIVE_ROOT / "replay_m2v5_runner_only_gates_recovery_v24.py",
}
REVIEW_PATHS = {
    "mendel": ARCHIVE_ROOT / "mendel_request_changes.json",
    "nash": ARCHIVE_ROOT / "nash_request_changes.json",
    "external_claude_opus": (
        ARCHIVE_ROOT
        / "20260723_external_claude_schema_recovery_review_v24.claude_cli.envelope.json"
    ),
}
TRANSPORT_PATHS = {
    "external_prompt": (
        ARCHIVE_ROOT / "20260723_external_claude_schema_recovery_review_v24_prompt.md"
    ),
    "external_schema": (
        ARCHIVE_ROOT / "20260723_external_claude_schema_recovery_review_v24.schema.json"
    ),
    "external_runner": ARCHIVE_ROOT / "run_external_claude_schema_recovery_review_v24.py",
    "external_stderr": (
        ARCHIVE_ROOT
        / "20260723_external_claude_schema_recovery_review_v24.claude_cli.stderr.txt"
    ),
}
FORMAL_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v24.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v24.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v24.json.ed25519.sig",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v24.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v24.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v24.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v24.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v24.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v24.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v24.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v24.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v24.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v24.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v24.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v24.json",
)


def artifact(path: Path) -> dict[str, object]:
    resolved = path.resolve(strict=True)
    data = resolved.read_bytes()
    mode = resolved.stat().st_mode & 0o7777
    return {
        "mode": f"0o{mode:o}",
        "path": str(resolved),
        "sha256": hashlib.sha256(data).hexdigest(),
        "size_bytes": len(data),
    }


def create_immutable(path: Path, data: bytes) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
    try:
        view = memoryview(data)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise RuntimeError(f"Short write: {path}")
            view = view[written:]
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
    finally:
        os.close(descriptor)


def main() -> int:
    if not ARCHIVE_ROOT.is_dir():
        raise RuntimeError(f"Archive root is absent: {ARCHIVE_ROOT}")
    occupied = [str(path) for path in FORMAL_OUTPUT_SLOTS if os.path.lexists(path)]
    if occupied:
        raise RuntimeError(f"Rejected v24 formal output slots are occupied: {occupied}")
    staging = list(
        RESULT_ROOT.glob(".tumor_ref_promotion_schema_recovery_authority.v24.bundle.staging.*")
    )
    if staging:
        raise RuntimeError(f"Rejected v24 staging paths are occupied: {staging}")

    source_set = {role: artifact(path) for role, path in SOURCE_PATHS.items()}
    reviews = {role: artifact(path) for role, path in REVIEW_PATHS.items()}
    transports = {role: artifact(path) for role, path in TRANSPORT_PATHS.items()}
    for record in (*source_set.values(), *reviews.values(), *transports.values()):
        if record["mode"] != "0o444":
            raise RuntimeError(f"Rejected v24 evidence is mutable: {record['path']}")

    mendel = json.loads(REVIEW_PATHS["mendel"].read_text(encoding="utf-8"))
    nash = json.loads(REVIEW_PATHS["nash"].read_text(encoding="utf-8"))
    external_envelope = json.loads(
        REVIEW_PATHS["external_claude_opus"].read_text(encoding="utf-8")
    )
    external = external_envelope["structured_output"]
    if (
        mendel.get("verdict") != "REQUEST_CHANGES"
        or nash.get("verdict") != "REQUEST_CHANGES"
        or external.get("verdict") != "APPROVE"
        or mendel.get("pass") is not False
        or nash.get("pass") is not False
        or external.get("pass") is not True
    ):
        raise RuntimeError("Strictest-review inputs do not match captured v24 verdicts")

    authority_private = artifact(AUTHORITY_KEY_ROOT / "ed25519_private_one_time.pem")
    authority_public = artifact(AUTHORITY_KEY_ROOT / "ed25519_public.pem")
    terminal_private = artifact(
        TERMINAL_KEY_ROOT / "ed25519_private_one_time_resident.pem"
    )
    terminal_public = artifact(TERMINAL_KEY_ROOT / "ed25519_public.pem")
    if authority_private["mode"] != "0o400" or terminal_private["mode"] != "0o400":
        raise RuntimeError("Rejected v24 unused private keys are not mode 0400")
    authority_private["state"] = "UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED"
    terminal_private["state"] = "UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED"

    evidence = {
        "schema_name": "intersubmod.tumor_ref_schema_recovery_pre_authority_candidate_rejection",
        "schema_version": "1.1.0",
        "candidate_generation": "v24",
        "status": "REJECTED_BEFORE_AUTHORITY_BY_STRICTEST_REVIEW_WINS",
        "successor_generation": "v25",
        "reviewed_source_set_sha256": (
            "dd4287bdd37a27ddf6bb5782806a814bebd39d2382d62c027f93c3cf81cf665b"
        ),
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
                "reviewer_agent_id": external["reviewer_agent_id"],
                "verdict": external["verdict"],
                "attribution": "Claude CLI session envelope",
            },
        },
        "review_transport_evidence": transports,
        "review_attribution_semantics": (
            "Reviewer IDs are orchestrator-recorded transport attribution only and are "
            "not cryptographic proof of reviewer authorship."
        ),
        "strictest_review_summary": {
            "approve_count": 1,
            "request_changes_count": 2,
            "blocking_high_count": 0,
            "blocking_medium_count": 4,
            "unique_blocking_finding_count": 3,
            "blocking_findings": [
                {
                    "id": "V24-M1",
                    "severity": "medium",
                    "summary": (
                        "Node and Chromium were identity-checked but not actually executed "
                        "by regression tests; missing/substituted runtime rejection was absent."
                    ),
                },
                {
                    "id": "V24-M2",
                    "severity": "medium",
                    "summary": (
                        "The occupied-state parameterization covered only 334 of 349 "
                        "forbidden slots while the authority claimed all 349."
                    ),
                },
                {
                    "id": "V24-M3",
                    "severity": "medium",
                    "summary": (
                        "RECOVERY_SCOPE omitted rejected v23 and AUTHORITY_CHECKS claimed "
                        "22 staging-pattern regressions although the active contract had 23."
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
            "authority": {"private_key": authority_private, "public_key": authority_public},
            "terminal": {"private_key": terminal_private, "public_key": terminal_public},
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
    create_immutable(
        ARCHIVE_ROOT / "rejection_evidence.json",
        (json.dumps(evidence, ensure_ascii=False, indent=2) + "\n").encode("utf-8"),
    )
    summary = """<!--
建立時間: 2026-07-23
目標: 封存 recovery v24 在 authority 建立前被獨立審查拒絕的完整證據
處理範圍: 七個凍結來源、三方審查 transport、未使用金鑰與 15 個正式輸出槽位
關聯檔案: rejection_evidence.json
-->

# Recovery v24 pre-authority rejection

v24 依 strictest-review-wins 規則於 authority 建立前拒絕。Mendel 與 Nash
一致指出 regression 未實際執行 Node 與 Chromium，亦缺少 missing/substituted
runtime 的 fail-closed 測試。Nash 另指出 occupied-state regression 只覆蓋
334/349 slots，scope 漏列 v23，且 authority 宣稱 22 patterns 而實際為 23。

外部 Claude Opus 回報 APPROVE，但不凌駕兩份 REQUEST_CHANGES。v24
authority、正式 review 與 V/R/C 輸出均未建立；authority-v24 與 terminal-v14
私鑰保持未使用並標記永不重用。科學資料、469,849 位點、read-tag join 與
結論上限均未改變。
"""
    create_immutable(ARCHIVE_ROOT / "SUMMARY.md", summary.encode("utf-8"))
    print(
        json.dumps(
            {
                "evidence": artifact(ARCHIVE_ROOT / "rejection_evidence.json"),
                "formal_output_slots_absent": True,
                "summary": artifact(ARCHIVE_ROOT / "SUMMARY.md"),
                "pass": True,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
