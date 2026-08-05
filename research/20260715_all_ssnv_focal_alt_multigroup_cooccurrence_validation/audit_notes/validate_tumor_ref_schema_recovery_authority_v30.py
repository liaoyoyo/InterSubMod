#!/usr/bin/env python3
"""Validate the append-only tumor-REF transition-relation recovery authority."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
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

PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")
PYTHON_TARGET = Path(
    "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python3.11"
)
EXPECTED_PYTHON_RUNTIME = {
    "path": str(PYTHON_TARGET),
    "size_bytes": 25_409_784,
    "sha256": "777797a57eb75b28f530191628d26b14afada9ced2cb51c0ecae1eb62796062e",
    "mode": "0o775",
}
OPENSSL = Path("/usr/bin/openssl")
AUTHORITY_ID = "20260724_tumor_ref_schema_recovery_v30"
SCHEMA_NAME = "intersubmod.tumor_ref_schema_recovery_authority"
SCHEMA_VERSION = "1.0.0"
REVIEW_SCHEMA_NAME = "intersubmod.tumor_ref_schema_recovery_review"
REVIEW_ATTRIBUTION_SEMANTICS = (
    "Review artifacts are immutable authority-bound coordinator records attributed "
    "to distinct multi-agent or Claude CLI transport identifiers. They are not "
    "cryptographic signatures by the named reviewers and do not prove reviewer "
    "authorship independently of the coordinating session."
)
AUTHORITY_STATUS = (
    "APPROVED_FOR_TUMOR_REF_V1_V6_AUDIT_LINEAGE_SPLIT_CHRONOLOGY_"
    "AND_HISTORICAL_CURRENT_RUNTIME_PROJECTION_SEPARATION_"
    "AND_BOUND_FINAL_DATASET_REPORT_RECOVERY_V30"
)
HISTORICAL_V19_V21_AUTHORITY_STATUS = (
    "APPROVED_FOR_TRANSITION_ALIAS_METADATA_PLUS_SIZE_"
    "AND_TERMINAL_KEY_RECOVERY_ONLY"
)
AUTHORITY_PASS_SEMANTICS = (
    "authorizes only transition-context-aware historical identity handling, exact "
    "six-key executable-alias binding, exact metadata-plus-size relation validation "
    "without reading private-key bytes, exact legacy eight-key stat relation "
    "validation without inferring link_count or registering historical records as "
    "live identities, strict metadata-to-metadata-plus-size enrichment, exact "
    "tumor-REF summary alias-to-canonical-target binding with descriptor protection, "
    "a signed two-directory inotify inventory, transient-restore regression, and "
    "parent-watch-before-snapshot mode-000 fallback with setup identity recheck, "
    "plus exact immutable tumor-REF v1 to source-authorized v6 primary-audit "
    "lineage and split producer chronology, canonical-path and descriptor-bound "
    "builder/probe/test execution, "
    "descriptor-bound Python execution with canonical argv0, and a distinct "
    "fresh four-role bootstrap with separately bound pre-authority and failed-v8 "
    "archives, an exact reviewer digest binding over the prior failed signed-"
    "recovery aggregate, and a fresh v21 key for recovery-v30 terminal outputs while "
    "preserving the legacy v2 signed key contract and quarantining the unused "
    "failed-v16 v6, failed-v17 v7, failed-v18 v8, failed-v19 v9, rejected-v20 "
    "v10, failed-v21 v11, rejected-v22 v12, rejected-v23 v13, rejected-v24 "
    "v14, failed-v25 v15, failed-v26 v16, rejected-v27 v17, failed-v28 "
    "v18, and failed-v29 v19 keys; reviewer "
    "identifiers are "
    "transport attribution rather than cryptographic reviewer-authorship proof; does not "
    "alter scientific "
    "payload, canonical receipt bytes, or biological claim ceiling"
)

KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260724_v30_four_role_rebootstrap"
)
PUBLIC_KEY = KEY_ROOT / "ed25519_public.pem"
PRIVATE_KEY = KEY_ROOT / "ed25519_private_one_time.pem"
EXPECTED_PUBLIC_KEY_SHA256 = (
    "a5b0e0b2c2a9f220d988597b47c8eb1d5446de401a932102948d829ffd0611ed"
)
AUTHORIZED_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260722_m2v5_terminal_v2"
)
AUTHORIZED_CONTINUATION_PUBLIC_KEY = (
    AUTHORIZED_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
AUTHORIZED_CONTINUATION_PRIVATE_KEY = (
    AUTHORIZED_CONTINUATION_KEY_ROOT / "ed25519_private_one_time_resident.pem"
)
EXPECTED_AUTHORIZED_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "a98370415594ebc9e8eb166bb70cab3550067b8a83db7d73ad23fac0891cb8b5"
)
RECOVERY_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260724_m2v5_terminal_v21"
)
RECOVERY_CONTINUATION_PUBLIC_KEY = (
    RECOVERY_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
RECOVERY_CONTINUATION_PRIVATE_KEY = (
    RECOVERY_CONTINUATION_KEY_ROOT / "ed25519_private_one_time_resident.pem"
)
EXPECTED_RECOVERY_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "3ea7624ed42caba9bd51ade25a4c9a037f0b84689b4e2a3563c8205bbb136fcd"
)
AUTHORITY_BUNDLE = RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v30.bundle"
AUTHORITY = AUTHORITY_BUNDLE / "authority.json"
AUTHORITY_SIGNATURE = AUTHORITY_BUNDLE / "authority.ed25519.sig"
AUTHORITY_COMMIT = AUTHORITY_BUNDLE / "commit.json"

VALIDATOR = AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v30.py"
RECOVERY_VERIFIER = AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v30.py"
RECOVERY_REPLAYER = AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v30.py"
RECOVERY_CONTINUATION = (
    AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v30.py"
)
READONLY_PROBE = AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v30.py"
REGRESSION_TEST = (
    AUDIT_ROOT
    / "schema_recovery_tests"
    / "test_tumor_ref_schema_recovery_v30.py"
)
CEREMONY_BUILDER = AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v30.py"
FOUR_ROLE_BOOTSTRAP_BUILDER = AUDIT_ROOT / "bootstrap_v30_four_role_keys.py"
RECOVERY_FINAL_DATASET_BUILDER = (
    AUDIT_ROOT / "build_all_ssnv_final_report_dataset_schema_recovery_v29.py"
)
RECOVERY_RESULT_REPORT_FINALIZER = (
    AUDIT_ROOT / "finalize_task_b_result_release_schema_recovery_v30.py"
)
RECOVERY_REPORT_BUILDER = (
    AUDIT_ROOT / "build_all_ssnv_report_artifact_schema_recovery_v29.py"
)
READONLY_PROBE_COMMAND = [
    "/usr/bin/env",
    "-i",
    "PATH=/usr/bin:/bin",
    "HOME=/tmp",
    "LANG=C.UTF-8",
    "LC_ALL=C.UTF-8",
    "PYTHONHASHSEED=0",
    "PYTHONNOUSERSITE=1",
    "PYTHONDONTWRITEBYTECODE=1",
    "OMP_NUM_THREADS=1",
    "OPENBLAS_NUM_THREADS=1",
    "MKL_NUM_THREADS=1",
    "NUMEXPR_NUM_THREADS=1",
    "BLIS_NUM_THREADS=1",
    str(PYTHON),
    "-I",
    "-B",
    str(READONLY_PROBE),
]

LEGACY_VERIFIER = AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_v2.py"
LEGACY_REPLAYER = AUDIT_ROOT / "replay_m2v5_runner_only_gates_v1.py"
LEGACY_CONTINUATION = AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_v1.py"

SOURCE_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json"
)
CANONICAL_RECEIPT = (
    WORKSPACE_ROOT
    / "tumor_ref_recovery_source_identity_v1"
    / "post_run_source_identity.receipt.json"
)
PROMOTION_AUTHORIZATION = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_authorization.v3.json"
)
PREPARE_ATTESTATION = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_prepare_exit_attestation.v1.json"
)
PREPARE_SIGNATURE = Path(str(PREPARE_ATTESTATION) + ".ed25519.sig")
PROMOTION_COMPLETION = RESULT_ROOT / "tumor_ref_source_receipt_promotion.v3.json"
PROMOTE_ATTESTATION = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_promote_exit_attestation.v1.json"
)
PROMOTE_SIGNATURE = Path(str(PROMOTE_ATTESTATION) + ".ed25519.sig")

RECOVERY_VERIFICATION_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v30.json"
)
RECOVERY_REPLAY_RECEIPT = RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.json"
RECOVERY_REPLAY_LOG = RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.log"
RECOVERY_REPLAY_SUCCESS_WITNESS = (
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.success_witness.json"
)
RECOVERY_CONTINUATION_RECEIPT = (
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v30.json"
)
RECOVERY_EXIT_ATTESTATION = (
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v30.json"
)
RECOVERY_EXIT_SIGNATURE = Path(str(RECOVERY_EXIT_ATTESTATION) + ".ed25519.sig")
RECOVERY_SUCCESS_WITNESS = (
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v30.json"
)
RECOVERY_INCIDENT = RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v30.json"

REVIEW_PATHS = {
    "mendel": REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v30.json",
    "nash": REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_nash.v30.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_external_claude_opus.v30.json"
    ),
}
EXPECTED_REVIEWERS = {
    "mendel": "Mendel",
    "nash": "Nash",
    "external_claude_opus": "External Claude Opus",
}
EXPECTED_REVIEWER_AGENT_IDS = {
    "mendel": "019f929e-b3db-7c81-b483-ca9c4bcdf155",
    "nash": "019f929f-3d98-7390-bec5-294f49ae2c56",
    "external_claude_opus": None,
}

RECOVERY_SOURCE_PATHS = {
    "authority_validator": VALIDATOR,
    "continuation_verifier": RECOVERY_VERIFIER,
    "runner_gate_replay": RECOVERY_REPLAYER,
    "downstream_continuation": RECOVERY_CONTINUATION,
    "readonly_probe": READONLY_PROBE,
    "regression_tests": REGRESSION_TEST,
    "ceremony_builder": CEREMONY_BUILDER,
    "four_role_bootstrap_builder": FOUR_ROLE_BOOTSTRAP_BUILDER,
    "recovery_final_dataset_builder": RECOVERY_FINAL_DATASET_BUILDER,
    "recovery_result_report_finalizer": RECOVERY_RESULT_REPORT_FINALIZER,
    "recovery_report_builder": RECOVERY_REPORT_BUILDER,
}
LEGACY_SOURCE_PATHS = {
    "continuation_verifier": LEGACY_VERIFIER,
    "runner_gate_replay": LEGACY_REPLAYER,
    "downstream_continuation": LEGACY_CONTINUATION,
}
ORIGINAL_CHAIN_PATHS = {
    "promotion_authorization": PROMOTION_AUTHORIZATION,
    "prepare_exit_attestation": PREPARE_ATTESTATION,
    "prepare_exit_attestation_signature": PREPARE_SIGNATURE,
    "promotion_completion": PROMOTION_COMPLETION,
    "promote_exit_attestation": PROMOTE_ATTESTATION,
    "promote_exit_attestation_signature": PROMOTE_SIGNATURE,
    "source_receipt": SOURCE_RECEIPT,
    "canonical_receipt": CANONICAL_RECEIPT,
}

PRIOR_RECOVERY_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260722_v1/ed25519_public.pem"
)
PRIOR_RECOVERY_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260722_v1/ed25519_private_one_time.pem"
)
PRIOR_RECOVERY_CHAIN_PATHS = {
    "authority": RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v1.json",
    "authority_signature": RESULT_ROOT
    / "tumor_ref_promotion_schema_recovery_authority.v1.json.ed25519.sig",
    "public_key": PRIOR_RECOVERY_PUBLIC_KEY,
    "review_mendel": REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_mendel.v1.json",
    "review_nash": REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_nash.v1.json",
    "review_external_claude_opus": REVIEW_ROOT
    / "20260722_tumor_ref_schema_recovery_external_claude_opus.v1.json",
    "authority_validator": AUDIT_ROOT
    / "validate_tumor_ref_schema_recovery_authority_v1.py",
    "runner_gate_replay": AUDIT_ROOT
    / "replay_m2v5_runner_only_gates_recovery_v1.py",
    "downstream_continuation": AUDIT_ROOT
    / "continue_m2v5_after_tumor_ref_promotion_recovery_v1.py",
    "verification_receipt": RESULT_ROOT
    / "tumor_ref_source_receipt_promotion_verification.recovery.v1.json",
    "runner_failure_evidence": AUDIT_ROOT
    / "20260722_rrec_v1_historical_transition_live_binding_failure.v1.json",
}

REJECTED_V2_EVIDENCE = (
    AUDIT_ROOT / "20260722_recovery_v2_formal_rejection_and_key_retirement.v1.json"
)
REJECTED_V2_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260722_v2/ed25519_private_one_time.pem"
)
REJECTED_V2_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v2.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v2.json.ed25519.sig",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_mendel.v2.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_nash.v2.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_external_claude_opus.v2.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v2.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v2.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v2.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v2.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v2.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v2.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v2.json",
)
EXPECTED_REJECTED_V2_EVIDENCE_SHA256 = (
    "f0c64e46ab34872149ca0cae830e79f48d818b591f8100160ac0cce44ed5ba84"
)
REJECTED_V3_EVIDENCE = (
    AUDIT_ROOT / "20260722_recovery_v3_formal_rejection_and_key_retirement.v1.json"
)
REJECTED_V3_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260722_v3/ed25519_private_one_time.pem"
)
REJECTED_V3_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v3.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v3.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v3.json.ed25519.sig",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_mendel.v3.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_nash.v3.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_external_claude_opus.v3.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v3.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v3.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v3.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v3.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v3.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v3.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v3.json",
)
EXPECTED_REJECTED_V3_EVIDENCE_SHA256 = (
    "c31009328e2130449422e01a5fd766446f6673e6f2625be0f4f380e3f41e4ef5"
)
REJECTED_V4_EVIDENCE = (
    AUDIT_ROOT / "20260722_recovery_v4_formal_rejection_and_key_retirement.v1.json"
)
REJECTED_V4_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260722_v4/ed25519_private_one_time.pem"
)
REJECTED_V4_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v4.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v4.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v4.json.ed25519.sig",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_mendel.v4.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_nash.v4.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_external_claude_opus.v4.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v4.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v4.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v4.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v4.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v4.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v4.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v4.json",
)
EXPECTED_REJECTED_V4_EVIDENCE_SHA256 = (
    "a1044ae1a0580b9a6587e30ddfeea22afeac84128cae0d8a28ee3a64619b9fb3"
)
REJECTED_V5_EVIDENCE = (
    AUDIT_ROOT / "20260722_recovery_v5_formal_rejection_and_key_retirement.v1.json"
)
REJECTED_V5_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260722_v5/ed25519_private_one_time.pem"
)
REJECTED_V5_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v5.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v5.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v5.json.ed25519.sig",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_mendel.v5.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_nash.v5.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_external_claude_opus.v5.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v5.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v5.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v5.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v5.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v5.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v5.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v5.json",
)
EXPECTED_REJECTED_V5_EVIDENCE_SHA256 = (
    "0495af623ba463b822f4f823ce28d49745c11c056ae391a36deeda06e0e78047"
)
REJECTED_V6_EVIDENCE = (
    AUDIT_ROOT / "20260722_recovery_v6_formal_rejection_and_key_retirement.v1.json"
)
REJECTED_V6_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260722_v6/ed25519_private_one_time.pem"
)
REJECTED_V6_ARCHIVED_REVIEWS = {
    "mendel": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260722_recovery_v6_post_review_probe_rejected"
        / "reviews"
        / "20260722_tumor_ref_schema_recovery_mendel.v6.json"
    ),
    "nash": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260722_recovery_v6_post_review_probe_rejected"
        / "reviews"
        / "20260722_tumor_ref_schema_recovery_nash.v6.json"
    ),
    "external_claude_opus": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260722_recovery_v6_post_review_probe_rejected"
        / "reviews"
        / "20260722_tumor_ref_schema_recovery_external_claude_opus.v6.json"
    ),
}
REJECTED_V6_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v6.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v6.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v6.json.ed25519.sig",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_mendel.v6.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_nash.v6.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_external_claude_opus.v6.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v6.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v6.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v6.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v6.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v6.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v6.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v6.json",
)
EXPECTED_REJECTED_V6_EVIDENCE_SHA256 = (
    "2a8a4a779d1e6df8a31e90adb64ae37e89d075ce78a3f0cba8942aea4359c9ab"
)
REJECTED_V7_EVIDENCE = (
    AUDIT_ROOT / "20260722_recovery_v7_formal_rejection_and_key_retirement.v1.json"
)
REJECTED_V7_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260722_v7/ed25519_private_one_time.pem"
)
REJECTED_V7_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v7.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v7.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v7.json.ed25519.sig",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_mendel.v7.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_nash.v7.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_external_claude_opus.v7.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v7.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v7.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v7.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v7.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v7.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v7.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v7.json",
)
EXPECTED_REJECTED_V7_EVIDENCE_SHA256 = (
    "07858f6f990cb9d72b9fca421d3e064b68d8641ec056c96a25c8adacb98c5956"
)
REJECTED_V8_EVIDENCE = (
    AUDIT_ROOT / "20260722_recovery_v8_formal_rejection_and_key_retirement.v1.json"
)
REJECTED_V8_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260722_v8/ed25519_private_one_time.pem"
)
REJECTED_V8_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v8.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v8.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v8.json.ed25519.sig",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_mendel.v8.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_nash.v8.json",
    REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_external_claude_opus.v8.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v8.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v8.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v8.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v8.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v8.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v8.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v8.json",
)
EXPECTED_REJECTED_V8_EVIDENCE_SHA256 = (
    "fe1724c88f2e168322fa16858db7450afd7145d8a0cbda097f019bfa052f6cd9"
)

FAILED_V9_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v9_signed_authority_runtime_17_23_slot_mismatch"
)
FAILED_V9_BUNDLE = (
    FAILED_V9_ROOT / "tumor_ref_promotion_schema_recovery_authority.v9.bundle"
)
FAILED_V9_EVIDENCE = FAILED_V9_ROOT / "failure_evidence.json"
FAILED_V9_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260722_v9/ed25519_public.pem"
)
FAILED_V9_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260722_v9/ed25519_private_one_time.pem"
)
FAILED_V9_ARCHIVED_ARTIFACT_PATHS = {
    "authority": FAILED_V9_BUNDLE / "authority.json",
    "commit": FAILED_V9_BUNDLE / "commit.json",
    "external_review": (
        FAILED_V9_ROOT
        / "20260722_tumor_ref_schema_recovery_external_claude_opus.v9.json"
    ),
    "mendel_review": (
        FAILED_V9_ROOT / "20260722_tumor_ref_schema_recovery_mendel.v9.json"
    ),
    "nash_review": (
        FAILED_V9_ROOT / "20260722_tumor_ref_schema_recovery_nash.v9.json"
    ),
    "signature": FAILED_V9_BUNDLE / "authority.ed25519.sig",
}
FAILED_V9_ORIGINAL_REVIEW_PATHS = {
    "mendel": REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_mendel.v9.json",
    "nash": REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_nash.v9.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260722_tumor_ref_schema_recovery_external_claude_opus.v9.json"
    ),
}
FAILED_V9_ARCHIVED_REVIEW_KEYS = {
    "mendel": "mendel_review",
    "nash": "nash_review",
    "external_claude_opus": "external_review",
}
FAILED_V9_EXPECTED_REVIEWERS = {
    "mendel": ("Mendel", "019f8bc7-c78d-79e2-8bef-75a8bbea8ff8"),
    "nash": ("Nash", "019f8b75-d075-7a22-848c-f9fed3fae06e"),
    "external_claude_opus": (
        "External Claude Opus",
        "7724f023-2491-47ed-9a1e-7856ae194add",
    ),
}
FAILED_V9_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v9.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v9.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v9.json.ed25519.sig",
    *FAILED_V9_ORIGINAL_REVIEW_PATHS.values(),
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v9.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v9.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v9.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v9.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v9.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v9.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v9.json",
)
EXPECTED_FAILED_V9_EVIDENCE_SHA256 = (
    "024d58b94c61c7b94ff3e842d47039fb228cb69e09e29806a000caef7f041bcc"
)

FAILED_V10_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
)
FAILED_V10_BUNDLE = (
    FAILED_V10_ROOT / "tumor_ref_promotion_schema_recovery_authority.v10.bundle"
)
FAILED_V10_EVIDENCE = FAILED_V10_ROOT / "failure_evidence.json"
FAILED_V10_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v10/ed25519_public.pem"
)
FAILED_V10_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v10/ed25519_private_one_time.pem"
)
FAILED_V10_SOURCE_PATHS = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v10.py"
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v10.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v10.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v10.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v10.py",
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v10.py"
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v10.py"
    ),
}
FAILED_V10_ARCHIVED_ARTIFACT_PATHS = {
    "authority": FAILED_V10_BUNDLE / "authority.json",
    "authority_signature": FAILED_V10_BUNDLE / "authority.ed25519.sig",
    "authority_commit": FAILED_V10_BUNDLE / "commit.json",
    "review_mendel": (
        FAILED_V10_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v10.json"
    ),
    "review_nash": (
        FAILED_V10_ROOT / "20260723_tumor_ref_schema_recovery_nash.v10.json"
    ),
    "review_external": (
        FAILED_V10_ROOT
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v10.json"
    ),
    "verification_receipt": (
        FAILED_V10_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v10.json"
    ),
    "replay_receipt": (
        FAILED_V10_ROOT / "m2v5_runner_only_gate_replay.recovery.v10.json"
    ),
    "replay_log": FAILED_V10_ROOT / "m2v5_runner_only_gate_replay.recovery.v10.log",
}
FAILED_V10_ORIGINAL_REVIEW_PATHS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v10.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v10.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v10.json"
    ),
}
FAILED_V10_ARCHIVED_REVIEW_KEYS = {
    "mendel": "review_mendel",
    "nash": "review_nash",
    "external_claude_opus": "review_external",
}
FAILED_V10_EXPECTED_REVIEWERS = {
    "mendel": ("Mendel", "019f8bc7-c78d-79e2-8bef-75a8bbea8ff8"),
    "nash": ("Nash", "019f8b75-d075-7a22-848c-f9fed3fae06e"),
    "external_claude_opus": (
        "External Claude Opus",
        "7724f023-2491-47ed-9a1e-7856ae194add",
    ),
}
FAILED_V10_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v10.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v10.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v10.json.ed25519.sig",
    *FAILED_V10_ORIGINAL_REVIEW_PATHS.values(),
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v10.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v10.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v10.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v10.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v10.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v10.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v10.json",
)
EXPECTED_FAILED_V10_EVIDENCE_SHA256 = (
    "25768002648008ad0ae6c459731f33ece68d0c74f484b4084de0fa89a6227701"
)

FAILED_V11_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
)
FAILED_V11_BUNDLE = (
    FAILED_V11_ROOT / "tumor_ref_promotion_schema_recovery_authority.v11.bundle"
)
FAILED_V11_EVIDENCE = FAILED_V11_ROOT / "failure_evidence.json"
FAILED_V11_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v11/ed25519_public.pem"
)
FAILED_V11_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v11/ed25519_private_one_time.pem"
)
FAILED_V11_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260722_m2v5_terminal_v2/ed25519_public.pem"
)
FAILED_V11_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260722_m2v5_terminal_v2/ed25519_private_one_time_resident.pem"
)
FAILED_V11_SOURCE_PATHS = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v11.py"
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v11.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v11.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v11.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v11.py",
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v11.py"
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v11.py"
    ),
}
FAILED_V11_ARCHIVED_ARTIFACT_PATHS = {
    "authority": FAILED_V11_BUNDLE / "authority.json",
    "authority_signature": FAILED_V11_BUNDLE / "authority.ed25519.sig",
    "authority_commit": FAILED_V11_BUNDLE / "commit.json",
    "review_mendel": (
        FAILED_V11_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v11.json"
    ),
    "review_nash": (
        FAILED_V11_ROOT / "20260723_tumor_ref_schema_recovery_nash.v11.json"
    ),
    "review_external": (
        FAILED_V11_ROOT
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v11.json"
    ),
    "verification_receipt": (
        FAILED_V11_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v11.json"
    ),
    "replay_receipt": (
        FAILED_V11_ROOT / "m2v5_runner_only_gate_replay.recovery.v11.json"
    ),
    "replay_log": FAILED_V11_ROOT / "m2v5_runner_only_gate_replay.recovery.v11.log",
}
FAILED_V11_ORIGINAL_REVIEW_PATHS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v11.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v11.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v11.json"
    ),
}
FAILED_V11_ARCHIVED_REVIEW_KEYS = {
    "mendel": "review_mendel",
    "nash": "review_nash",
    "external_claude_opus": "review_external",
}
FAILED_V11_EXPECTED_REVIEWERS = {
    "mendel": ("Mendel", "019f8bc7-c78d-79e2-8bef-75a8bbea8ff8"),
    "nash": ("Nash", "019f8b75-d075-7a22-848c-f9fed3fae06e"),
    "external_claude_opus": (
        "External Claude Opus",
        "7724f023-2491-47ed-9a1e-7856ae194add",
    ),
}
FAILED_V11_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v11.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v11.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v11.json.ed25519.sig",
    *FAILED_V11_ORIGINAL_REVIEW_PATHS.values(),
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v11.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v11.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v11.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v11.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v11.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v11.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v11.json",
)
EXPECTED_FAILED_V11_EVIDENCE_SHA256 = (
    "bc76d6026a20857145feaeb188a4cc5f933ed8b262e2f3f7d107d4bde13f7619"
)

FAILED_V12_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v12_signed_authority_c_python_executable_alias_relation_mismatch"
)
FAILED_V12_BUNDLE = (
    FAILED_V12_ROOT / "tumor_ref_promotion_schema_recovery_authority.v12.bundle"
)
FAILED_V12_EVIDENCE = FAILED_V12_ROOT / "failure_evidence.json"
FAILED_V12_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v12/ed25519_public.pem"
)
FAILED_V12_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v12/ed25519_private_one_time.pem"
)
FAILED_V12_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260722_m2v5_terminal_v2/ed25519_public.pem"
)
FAILED_V12_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260722_m2v5_terminal_v2/ed25519_private_one_time_resident.pem"
)
FAILED_V12_SOURCE_PATHS = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v12.py"
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v12.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v12.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v12.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v12.py",
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v12.py"
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v12.py"
    ),
}
FAILED_V12_ARCHIVED_ARTIFACT_PATHS = {
    "authority": FAILED_V12_BUNDLE / "authority.json",
    "authority_signature": FAILED_V12_BUNDLE / "authority.ed25519.sig",
    "authority_commit": FAILED_V12_BUNDLE / "commit.json",
    "review_mendel": (
        FAILED_V12_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v12.json"
    ),
    "review_nash": (
        FAILED_V12_ROOT / "20260723_tumor_ref_schema_recovery_nash.v12.json"
    ),
    "review_external": (
        FAILED_V12_ROOT
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v12.json"
    ),
    "verification_receipt": (
        FAILED_V12_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v12.json"
    ),
    "replay_receipt": (
        FAILED_V12_ROOT / "m2v5_runner_only_gate_replay.recovery.v12.json"
    ),
    "replay_log": FAILED_V12_ROOT / "m2v5_runner_only_gate_replay.recovery.v12.log",
}
FAILED_V12_ORIGINAL_REVIEW_PATHS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v12.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v12.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v12.json"
    ),
}
FAILED_V12_ARCHIVED_REVIEW_KEYS = {
    "mendel": "review_mendel",
    "nash": "review_nash",
    "external_claude_opus": "review_external",
}
FAILED_V12_EXPECTED_REVIEWERS = {
    "mendel": ("Mendel", "019f8cf9-af38-7d71-8f02-b1eb09742504"),
    "nash": ("Nash", "019f8cf9-b625-7762-a018-b4cb0000bbd7"),
    "external_claude_opus": (
        "External Claude Opus",
        "7724f023-2491-47ed-9a1e-7856ae194add",
    ),
}
FAILED_V12_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v12.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v12.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v12.json.ed25519.sig",
    *FAILED_V12_ORIGINAL_REVIEW_PATHS.values(),
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v12.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v12.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v12.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v12.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v12.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v12.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v12.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v12.json",
)
EXPECTED_FAILED_V12_EVIDENCE_SHA256 = (
    "a5655e8c19fa0ef52e008bf045aefa09a6be1172a5adcfe8e153e62a0cbac62b"
)

REJECTED_V13_ROOT = (
    AUDIT_ROOT
    / "rejected_pre_authority_reviews"
    / "20260723_v13_candidate_nash_medium_findings"
)
REJECTED_V13_SUPERSEDED_EVIDENCE = REJECTED_V13_ROOT / "rejection_evidence.json"
EXPECTED_REJECTED_V13_SUPERSEDED_EVIDENCE_SHA256 = (
    "a6dad5b9c6ec7d26a6390751b8296ce0d7eea1f9dfc419633458a88a36dc05fa"
)
REJECTED_V13_EVIDENCE = REJECTED_V13_ROOT / "rejection_evidence.v2.json"
EXPECTED_REJECTED_V13_EVIDENCE_SHA256 = (
    "b42de4167900069fbe563244b9058a8291fe1baf097fb0a5f5e77dee74ee02eb"
)
REJECTED_V13_KEY_ARCHIVE_LEDGER = REJECTED_V13_ROOT / "key_archive_ledger.v1.jsonl"
EXPECTED_REJECTED_V13_KEY_ARCHIVE_LEDGER_SHA256 = (
    "0f4b9e0d6d48187cd7b63bafa46472efb83d8f1f3d790d4503e8423102799321"
)
REJECTED_V13_SOURCE_PATHS = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v13.py"
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v13.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v13.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v13.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v13.py",
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v13.py"
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v13.py"
    ),
}
REJECTED_V13_AUTHORITY_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "archive/20260723_v13_unused_pre_authority_review_rejection_01"
)
REJECTED_V13_TERMINAL_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "archive/20260723_m2v5_terminal_v3_unused_pre_authority_review_rejection_01"
)
REJECTED_V13_ORIGINAL_AUTHORITY_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v13"
)
REJECTED_V13_ORIGINAL_TERMINAL_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v3"
)
REJECTED_V13_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v13.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v13.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v13.json.ed25519.sig",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v13.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v13.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v13.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v13.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v13.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v13.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v13.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v13.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v13.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v13.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v13.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v13.json",
)

FAILED_V14_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v14_signed_authority_r_metadata_only_relation_schema_mismatch"
)
FAILED_V14_BUNDLE = (
    FAILED_V14_ROOT / "tumor_ref_promotion_schema_recovery_authority.v14.bundle"
)
FAILED_V14_EVIDENCE = FAILED_V14_ROOT / "failure_evidence.json"
EXPECTED_FAILED_V14_EVIDENCE_SHA256 = (
    "8e3f5c13be270279a11f4e3e11d708a5815d787e5b0eccaf524885616166af9c"
)
FAILED_V14_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v14/ed25519_public.pem"
)
FAILED_V14_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v14/ed25519_private_one_time.pem"
)
FAILED_V14_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v4/ed25519_public.pem"
)
FAILED_V14_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v4/ed25519_private_one_time_resident.pem"
)
FAILED_V14_SOURCE_PATHS = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v14.py"
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v14.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v14.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v14.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v14.py",
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v14.py"
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v14.py"
    ),
    "archive_script": AUDIT_ROOT / "archive_v14_signed_r_metadata_relation_failure.py",
}
FAILED_V14_ARCHIVED_ARTIFACT_PATHS = {
    "authority": FAILED_V14_BUNDLE / "authority.json",
    "authority_signature": FAILED_V14_BUNDLE / "authority.ed25519.sig",
    "authority_commit": FAILED_V14_BUNDLE / "commit.json",
    "verification_receipt": (
        FAILED_V14_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v14.json"
    ),
    "formal_replay_stderr": FAILED_V14_ROOT / "formal_replay_stderr.log",
    "review_mendel": (
        FAILED_V14_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v14.json"
    ),
    "review_nash": (
        FAILED_V14_ROOT / "20260723_tumor_ref_schema_recovery_nash.v14.json"
    ),
    "review_external_claude_opus": (
        FAILED_V14_ROOT
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v14.json"
    ),
}
FAILED_V14_ORIGINAL_REVIEW_PATHS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v14.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v14.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v14.json"
    ),
}
FAILED_V14_ARCHIVED_REVIEW_KEYS = {
    "mendel": "review_mendel",
    "nash": "review_nash",
    "external_claude_opus": "review_external_claude_opus",
}
FAILED_V14_EXPECTED_REVIEWERS = {
    "mendel": ("Mendel", "019f8d4b-8e42-7830-9efe-d92e5d197cf8"),
    "nash": ("Nash", "019f8d4b-9388-7e53-9f3d-dd134beed980"),
    "external_claude_opus": (
        "External Claude Opus",
        "828d4dcb-1aed-4ade-bbf8-433ee06a895e",
    ),
}
FAILED_V14_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v14.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v14.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v14.json.ed25519.sig",
    *FAILED_V14_ORIGINAL_REVIEW_PATHS.values(),
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v14.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v14.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v14.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v14.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v14.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v14.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v14.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v14.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v14.json",
)

REJECTED_V15_ROOT = (
    AUDIT_ROOT
    / "rejected_pre_authority_reviews"
    / "20260723_v15_candidate_runtime_schema_and_provenance_findings"
)
REJECTED_V15_EVIDENCE = REJECTED_V15_ROOT / "rejection_evidence.json"
EXPECTED_REJECTED_V15_EVIDENCE_SHA256 = (
    "40b329951a199fdf33e2efa7c73b2d7779d12d5c707d237ed7cc9e1734930297"
)
REJECTED_V15_KEY_ARCHIVE_LEDGER = REJECTED_V15_ROOT / "key_archive_ledger.v1.jsonl"
EXPECTED_REJECTED_V15_KEY_ARCHIVE_LEDGER_SHA256 = (
    "6409605672db27706ca1666225a8dae20ae31a61f4893414f2cf4fb284b9f9bc"
)
REJECTED_V15_REVIEW_PATHS = {
    "mendel": REJECTED_V15_ROOT / "mendel_request_changes.json",
    "nash": REJECTED_V15_ROOT / "nash_request_changes.json",
    "external_claude_opus": REJECTED_V15_ROOT / "external_claude_opus_approve.json",
}
REJECTED_V15_EXTERNAL_ENVELOPE = (
    AUDIT_ROOT
    / "20260723_external_claude_schema_recovery_review_v15.attempt2.claude_cli.envelope.json"
)
REJECTED_V15_EXTERNAL_STDERR = (
    AUDIT_ROOT
    / "20260723_external_claude_schema_recovery_review_v15.attempt2.claude_cli.stderr.txt"
)
REJECTED_V15_SOURCE_PATHS = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v15.py"
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v15.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v15.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v15.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v15.py",
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v15.py"
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v15.py"
    ),
}
REJECTED_V15_AUTHORITY_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "archive/20260723_v15_unused_pre_authority_review_rejection_01"
)
REJECTED_V15_TERMINAL_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "archive/20260723_m2v5_terminal_v5_unused_pre_authority_review_rejection_01"
)
REJECTED_V15_ORIGINAL_AUTHORITY_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v15"
)
REJECTED_V15_ORIGINAL_TERMINAL_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v5"
)
REJECTED_V15_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v15.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v15.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v15.json.ed25519.sig",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v15.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v15.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v15.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v15.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v15.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v15.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v15.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v15.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v15.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v15.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v15.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v15.json",
)

REJECTED_V19_ROUND1_ROOT = (
    AUDIT_ROOT
    / "rejected_pre_authority_reviews"
    / "20260723_v19_candidate_mutation_monitor_provenance_finding"
)
REJECTED_V19_ROUND1_EVIDENCE = REJECTED_V19_ROUND1_ROOT / "rejection_evidence.json"
EXPECTED_REJECTED_V19_ROUND1_EVIDENCE_SHA256 = (
    "d99750a48971e1e872661b6f6435a8b1adfb0a7f0dbdeeb6e47fe938e0c02b2f"
)
REJECTED_V19_ROUND1_REVIEW_PATHS = {
    "mendel": REJECTED_V19_ROUND1_ROOT / "mendel_request_changes.json",
    "nash": REJECTED_V19_ROUND1_ROOT / "nash_request_changes.json",
}
REJECTED_V19_ROUND1_EXTERNAL_ENVELOPE = (
    REJECTED_V19_ROUND1_ROOT
    / "20260723_external_claude_schema_recovery_review_v19.claude_cli.envelope.json"
)
REJECTED_V19_ROUND1_EXTERNAL_STDERR = (
    REJECTED_V19_ROUND1_ROOT
    / "20260723_external_claude_schema_recovery_review_v19.claude_cli.stderr.txt"
)
REJECTED_V19_ROUND1_TRANSPORT_PATHS = {
    "external_prompt": (
        REJECTED_V19_ROUND1_ROOT
        / "20260723_external_claude_schema_recovery_review_v19_prompt.md"
    ),
    "external_schema": (
        REJECTED_V19_ROUND1_ROOT
        / "20260723_external_claude_schema_recovery_review_v19.schema.json"
    ),
    "external_runner": (
        REJECTED_V19_ROUND1_ROOT
        / "run_external_claude_schema_recovery_review_v19.py"
    ),
}
REJECTED_V19_ROUND1_SOURCE_PATHS = {
    "authority_validator": (
        REJECTED_V19_ROUND1_ROOT
        / "source_snapshot/validate_tumor_ref_schema_recovery_authority_v19.py"
    ),
    "continuation_verifier": (
        REJECTED_V19_ROUND1_ROOT
        / "source_snapshot/verify_tumor_ref_receipt_promotion_recovery_v19.py"
    ),
    "runner_gate_replay": (
        REJECTED_V19_ROUND1_ROOT
        / "source_snapshot/replay_m2v5_runner_only_gates_recovery_v19.py"
    ),
    "downstream_continuation": (
        REJECTED_V19_ROUND1_ROOT
        / "source_snapshot/continue_m2v5_after_tumor_ref_promotion_recovery_v19.py"
    ),
    "readonly_probe": (
        REJECTED_V19_ROUND1_ROOT
        / "source_snapshot/probe_tumor_ref_schema_recovery_sources_v19.py"
    ),
    "regression_tests": (
        REJECTED_V19_ROUND1_ROOT
        / "source_snapshot/schema_recovery_tests/test_tumor_ref_schema_recovery_v19.py"
    ),
    "ceremony_builder": (
        REJECTED_V19_ROUND1_ROOT
        / "source_snapshot/build_tumor_ref_schema_recovery_authority_v19.py"
    ),
}
EXPECTED_REJECTED_V19_ROUND1_SOURCE_IDENTITIES = {
    "authority_validator": (303_580, "618382f2dd3b6da74a6a965a8d5236acb86cfaff477e35d7e36fd0936936c92d"),
    "continuation_verifier": (127_906, "a24abd4b90b997fc8d5ff16115a195922969260c8cb39aa7c682cde7fea06c4d"),
    "runner_gate_replay": (152_178, "b3486543089a72d6d2a987f054b5261fd307e26ef6301d5959883b487eea9cd9"),
    "downstream_continuation": (361_428, "ade2f9b9331b709487e06ca9a72bf46a248f7dc196edc7f01ec655d2472ccc2f"),
    "readonly_probe": (82_617, "26ccf2c1f43c83a9fdc1b30234116eafd90939d8b0bd7d486c184b3ff1961829"),
    "regression_tests": (144_833, "e05ab8bf210a231658bc2240d12427709b8b100fdd1bfc8359e6e05518a6983c"),
    "ceremony_builder": (53_317, "36ac487de3a4e687086eafdcb2170f8c5ae45209ee71f6ea25e250b926cc1129"),
}

REJECTED_V29_ROUND1_ROOT = (
    AUDIT_ROOT
    / "rejected_pre_authority_reviews"
    / "20260724_v29_round1_failed_v28_terminal_binding_gap"
)
REJECTED_V29_ROUND1_EVIDENCE = REJECTED_V29_ROUND1_ROOT / "rejection_evidence.json"
EXPECTED_REJECTED_V29_ROUND1_EVIDENCE_SHA256 = (
    "6d44a75b32737c0e16a7c05c3801ade734977336982a7f08260b70610442cafa"
)
REJECTED_V29_ROUND1_SUMMARY = REJECTED_V29_ROUND1_ROOT / "SUMMARY.md"
EXPECTED_REJECTED_V29_ROUND1_SUMMARY_SHA256 = (
    "56122d97ee422d9d80d5508ed0655cd7c7fb9be75887b2b4bc38ac0be8c4fa92"
)
REJECTED_V29_ROUND1_PAYLOAD_RELATIVE_PATHS = (
    "source_snapshot/build_all_ssnv_final_report_dataset_schema_recovery_v29.py",
    "source_snapshot/build_all_ssnv_report_artifact_schema_recovery_v29.py",
    "source_snapshot/build_tumor_ref_schema_recovery_authority_v29.py",
    "source_snapshot/continue_m2v5_after_tumor_ref_promotion_recovery_v29.py",
    "source_snapshot/finalize_task_b_result_release_schema_recovery_v29.py",
    "source_snapshot/probe_tumor_ref_schema_recovery_sources_v29.py",
    "source_snapshot/replay_m2v5_runner_only_gates_recovery_v29.py",
    "source_snapshot/schema_recovery_tests/test_tumor_ref_schema_recovery_v29.py",
    "source_snapshot/validate_tumor_ref_schema_recovery_authority_v29.py",
    "source_snapshot/verify_tumor_ref_receipt_promotion_recovery_v29.py",
    (
        "review_transport/20260724_external_claude_schema_recovery_review_v29."
        "claude_cli.envelope.json"
    ),
    (
        "review_transport/20260724_external_claude_schema_recovery_review_v29."
        "claude_cli.stderr.txt"
    ),
    "review_transport/20260724_external_claude_schema_recovery_review_v29.schema.json",
    "review_transport/20260724_external_claude_schema_recovery_review_v29_prompt.md",
    "review_transport/publish_tumor_ref_schema_recovery_reviews_v29.py",
    "review_transport/run_external_claude_schema_recovery_review_v29.py",
    "mendel_initial_approve.json",
    "mendel_corrected_request_changes.json",
    "nash_request_changes.json",
    "orchestrator_ast_reproduction.json",
)
REJECTED_V29_ROUND1_PAYLOAD_PATHS = {
    relative_path: REJECTED_V29_ROUND1_ROOT / relative_path
    for relative_path in REJECTED_V29_ROUND1_PAYLOAD_RELATIVE_PATHS
}

FAILED_V16_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v16_signed_authority_c_legacy_eight_key_stat_schema_mismatch"
)
FAILED_V16_BUNDLE = (
    FAILED_V16_ROOT / "tumor_ref_promotion_schema_recovery_authority.v16.bundle"
)
FAILED_V16_EVIDENCE = FAILED_V16_ROOT / "failure_evidence.json"
EXPECTED_FAILED_V16_EVIDENCE_SHA256 = (
    "c6ba5ced3fa16b06bb2e4dbd974c8db3e8f4c9f1c9304f1021da061835bf0c25"
)
FAILED_V16_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v16/ed25519_public.pem"
)
FAILED_V16_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v16/ed25519_private_one_time.pem"
)
FAILED_V16_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v6/ed25519_public.pem"
)
FAILED_V16_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v6/ed25519_private_one_time_resident.pem"
)
FAILED_V16_SOURCE_PATHS = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v16.py"
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v16.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v16.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v16.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v16.py",
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v16.py"
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v16.py"
    ),
    "archive_script": (
        AUDIT_ROOT / "archive_v16_signed_c_legacy_stat_schema_failure.py"
    ),
}
FAILED_V16_REVIEW_TRANSPORT_PATHS = {
    "external_prompt": (
        AUDIT_ROOT / "20260723_external_claude_schema_recovery_review_v16_prompt.md"
    ),
    "external_schema": (
        AUDIT_ROOT / "20260723_external_claude_schema_recovery_review_v16.schema.json"
    ),
    "external_attempt1_envelope": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v16.claude_cli.envelope.json"
    ),
    "external_attempt1_stderr": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v16.claude_cli.stderr.txt"
    ),
    "external_attempt2_envelope": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v16.attempt2.claude_cli.envelope.json"
    ),
    "external_attempt2_stderr": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v16.attempt2.claude_cli.stderr.txt"
    ),
}
FAILED_V16_ARCHIVED_ARTIFACT_PATHS = {
    "authority": FAILED_V16_BUNDLE / "authority.json",
    "authority_signature": FAILED_V16_BUNDLE / "authority.ed25519.sig",
    "authority_commit": FAILED_V16_BUNDLE / "commit.json",
    "verification_receipt": (
        FAILED_V16_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v16.json"
    ),
    "replay_receipt": (
        FAILED_V16_ROOT / "m2v5_runner_only_gate_replay.recovery.v16.json"
    ),
    "replay_log": (
        FAILED_V16_ROOT / "m2v5_runner_only_gate_replay.recovery.v16.log"
    ),
    "replay_success_witness": (
        FAILED_V16_ROOT
        / "m2v5_runner_only_gate_replay.recovery.v16.success_witness.json"
    ),
    "formal_continuation_stderr": (
        FAILED_V16_ROOT / "formal_continuation_stderr.log"
    ),
    "review_mendel": (
        FAILED_V16_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v16.json"
    ),
    "review_nash": (
        FAILED_V16_ROOT / "20260723_tumor_ref_schema_recovery_nash.v16.json"
    ),
    "review_external_claude_opus": (
        FAILED_V16_ROOT
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v16.json"
    ),
}
FAILED_V16_ORIGINAL_REVIEW_PATHS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v16.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v16.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v16.json"
    ),
}
FAILED_V16_ARCHIVED_REVIEW_KEYS = {
    "mendel": "review_mendel",
    "nash": "review_nash",
    "external_claude_opus": "review_external_claude_opus",
}
FAILED_V16_EXPECTED_REVIEWERS = {
    "mendel": ("Mendel", "019f8df0-28ad-78d0-b62c-81f512f90c8d"),
    "nash": ("Nash", "019f8df0-9010-71e0-9f39-e30df506f59f"),
    "external_claude_opus": (
        "External Claude Opus",
        "981fda5f-bebe-47cd-83b1-7a152421d474",
    ),
}
FAILED_V16_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v16.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v16.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v16.json.ed25519.sig",
    *FAILED_V16_ORIGINAL_REVIEW_PATHS.values(),
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v16.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v16.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v16.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v16.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v16.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v16.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v16.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v16.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v16.json",
)

FAILED_V17_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v17_signed_authority_c_metadata_enrichment_relation_conflict"
)
FAILED_V17_BUNDLE = (
    FAILED_V17_ROOT / "tumor_ref_promotion_schema_recovery_authority.v17.bundle"
)
FAILED_V17_EVIDENCE = FAILED_V17_ROOT / "failure_evidence.json"
EXPECTED_FAILED_V17_EVIDENCE_SHA256 = (
    "b237a8d03b4a483abd47b3239d5d8d13442794835c090dfe4122b60834b025a9"
)
FAILED_V17_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v17/ed25519_public.pem"
)
FAILED_V17_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v17/ed25519_private_one_time.pem"
)
FAILED_V17_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v7/ed25519_public.pem"
)
FAILED_V17_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v7/ed25519_private_one_time_resident.pem"
)
FAILED_V17_SOURCE_PATHS = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v17.py"
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v17.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v17.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v17.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v17.py",
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v17.py"
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v17.py"
    ),
    "archive_script": (
        AUDIT_ROOT / "archive_v17_signed_c_metadata_enrichment_failure.py"
    ),
}
FAILED_V17_REVIEW_TRANSPORT_PATHS = {
    "external_prompt": (
        AUDIT_ROOT / "20260723_external_claude_schema_recovery_review_v17_prompt.md"
    ),
    "external_schema": (
        AUDIT_ROOT / "20260723_external_claude_schema_recovery_review_v17.schema.json"
    ),
    "external_runner": (
        AUDIT_ROOT / "run_external_claude_schema_recovery_review_v17.py"
    ),
    "external_envelope": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v17.claude_cli.envelope.json"
    ),
    "external_stderr": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v17.claude_cli.stderr.txt"
    ),
}
FAILED_V17_ARCHIVED_ARTIFACT_PATHS = {
    "authority": FAILED_V17_BUNDLE / "authority.json",
    "authority_signature": FAILED_V17_BUNDLE / "authority.ed25519.sig",
    "authority_commit": FAILED_V17_BUNDLE / "commit.json",
    "verification_receipt": (
        FAILED_V17_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v17.json"
    ),
    "replay_receipt": (
        FAILED_V17_ROOT / "m2v5_runner_only_gate_replay.recovery.v17.json"
    ),
    "replay_log": (
        FAILED_V17_ROOT / "m2v5_runner_only_gate_replay.recovery.v17.log"
    ),
    "replay_success_witness": (
        FAILED_V17_ROOT
        / "m2v5_runner_only_gate_replay.recovery.v17.success_witness.json"
    ),
    "formal_continuation_stderr": (
        FAILED_V17_ROOT / "formal_continuation_stderr.log"
    ),
    "review_mendel": (
        FAILED_V17_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v17.json"
    ),
    "review_nash": (
        FAILED_V17_ROOT / "20260723_tumor_ref_schema_recovery_nash.v17.json"
    ),
    "review_external_claude_opus": (
        FAILED_V17_ROOT
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v17.json"
    ),
}
FAILED_V17_TRACE_ROOT = FAILED_V17_ROOT / "diagnostic_write_traces"
FAILED_V17_ORIGINAL_REVIEW_PATHS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v17.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v17.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v17.json"
    ),
}
FAILED_V17_ARCHIVED_REVIEW_KEYS = {
    "mendel": "review_mendel",
    "nash": "review_nash",
    "external_claude_opus": "review_external_claude_opus",
}
FAILED_V17_EXPECTED_REVIEWERS = {
    "mendel": ("Mendel", "019f8e32-7009-7dd3-8049-58994abdb27b"),
    "nash": ("Nash", "019f8e32-3b58-7391-8da9-271165963b06"),
    "external_claude_opus": (
        "External Claude Opus",
        "c0b814e1-1046-49d0-864d-ea37f933c35b",
    ),
}
FAILED_V17_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v17.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v17.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v17.json.ed25519.sig",
    *FAILED_V17_ORIGINAL_REVIEW_PATHS.values(),
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v17.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v17.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v17.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v17.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v17.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v17.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v17.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v17.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v17.json",
)

FAILED_V18_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v18_signed_authority_c_tumor_ref_summary_alias_noncanonical"
)
FAILED_V18_BUNDLE = (
    FAILED_V18_ROOT / "tumor_ref_promotion_schema_recovery_authority.v18.bundle"
)
FAILED_V18_EVIDENCE = FAILED_V18_ROOT / "failure_evidence.json"
EXPECTED_FAILED_V18_EVIDENCE_SHA256 = (
    "d4cca03adc3382f48177cf8cec705189c66ffbc063c5533bcf794629ccaead2b"
)
FAILED_V18_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v18/ed25519_public.pem"
)
FAILED_V18_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v18/ed25519_private_one_time.pem"
)
FAILED_V18_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v8/ed25519_public.pem"
)
FAILED_V18_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v8/ed25519_private_one_time_resident.pem"
)
FAILED_V18_SOURCE_PATHS = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v18.py"
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v18.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v18.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v18.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v18.py",
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v18.py"
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v18.py"
    ),
    "archive_script": (
        AUDIT_ROOT / "archive_v18_signed_c_tumor_ref_summary_alias_failure.py"
    ),
}
FAILED_V18_REVIEW_TRANSPORT_PATHS = {
    "external_prompt": (
        AUDIT_ROOT / "20260723_external_claude_schema_recovery_review_v18_prompt.md"
    ),
    "external_schema": (
        AUDIT_ROOT / "20260723_external_claude_schema_recovery_review_v18.schema.json"
    ),
    "external_runner": (
        AUDIT_ROOT / "run_external_claude_schema_recovery_review_v18.py"
    ),
    "external_envelope": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v18.claude_cli.envelope.json"
    ),
    "external_stderr": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v18.claude_cli.stderr.txt"
    ),
}
FAILED_V18_ARCHIVED_ARTIFACT_PATHS = {
    "authority": FAILED_V18_BUNDLE / "authority.json",
    "authority_signature": FAILED_V18_BUNDLE / "authority.ed25519.sig",
    "authority_commit": FAILED_V18_BUNDLE / "commit.json",
    "verification_receipt": (
        FAILED_V18_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v18.json"
    ),
    "replay_receipt": (
        FAILED_V18_ROOT / "m2v5_runner_only_gate_replay.recovery.v18.json"
    ),
    "replay_log": (
        FAILED_V18_ROOT / "m2v5_runner_only_gate_replay.recovery.v18.log"
    ),
    "replay_success_witness": (
        FAILED_V18_ROOT
        / "m2v5_runner_only_gate_replay.recovery.v18.success_witness.json"
    ),
    "formal_continuation_stderr": (
        FAILED_V18_ROOT / "formal_continuation_stderr.log"
    ),
    "diagnostic_readme": (
        FAILED_V18_ROOT / "diagnostic_write_traces" / "README.md"
    ),
    "review_mendel": (
        FAILED_V18_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v18.json"
    ),
    "review_nash": (
        FAILED_V18_ROOT / "20260723_tumor_ref_schema_recovery_nash.v18.json"
    ),
    "review_external_claude_opus": (
        FAILED_V18_ROOT
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v18.json"
    ),
}
FAILED_V18_TRACE_ROOT = FAILED_V18_ROOT / "diagnostic_write_traces"
FAILED_V18_TRACE_GLOB = "v18_c_write_trace.*"
FAILED_V18_ORIGINAL_REVIEW_PATHS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v18.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v18.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v18.json"
    ),
}
FAILED_V18_ARCHIVED_REVIEW_KEYS = {
    "mendel": "review_mendel",
    "nash": "review_nash",
    "external_claude_opus": "review_external_claude_opus",
}
FAILED_V18_EXPECTED_REVIEWERS = {
    "mendel": ("Mendel", "019f8e71-cab2-7ce1-aa6f-986efd4efac0"),
    "nash": ("Nash", "019f8e71-c521-7b00-b7df-82b1c2b3a6a1"),
    "external_claude_opus": (
        "External Claude Opus",
        "bf87dba8-975f-485b-9954-1f28aad9dd42",
    ),
}
FAILED_V18_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v18.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v18.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v18.json.ed25519.sig",
    *FAILED_V18_ORIGINAL_REVIEW_PATHS.values(),
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v18.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v18.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v18.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v18.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v18.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v18.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v18.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v18.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v18.json",
)

FAILED_V19_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v19_signed_authority_c_mode000_inotify_permission_denied"
)
FAILED_V19_BUNDLE = (
    FAILED_V19_ROOT / "tumor_ref_promotion_schema_recovery_authority.v19.bundle"
)
FAILED_V19_EVIDENCE = FAILED_V19_ROOT / "failure_evidence.json"
EXPECTED_FAILED_V19_EVIDENCE_SHA256 = (
    "45cb20e1f2346a278fef33ec83baeca527ad726e8c96526a7977215a1437268b"
)
FAILED_V19_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v19/ed25519_public.pem"
)
FAILED_V19_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v19/ed25519_private_one_time.pem"
)
FAILED_V19_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v9/ed25519_public.pem"
)
FAILED_V19_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v9/ed25519_private_one_time_resident.pem"
)
FAILED_V19_SOURCE_PATHS = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v19.py"
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v19.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v19.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v19.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v19.py",
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v19.py"
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v19.py"
    ),
    "archive_script": (
        AUDIT_ROOT / "archive_v19_signed_c_mode000_inotify_failure.py"
    ),
}
FAILED_V19_REVIEW_TRANSPORT_PATHS = {
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
FAILED_V19_ARCHIVED_ARTIFACT_PATHS = {
    "authority": FAILED_V19_BUNDLE / "authority.json",
    "authority_signature": FAILED_V19_BUNDLE / "authority.ed25519.sig",
    "authority_commit": FAILED_V19_BUNDLE / "commit.json",
    "verification_receipt": (
        FAILED_V19_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v19.json"
    ),
    "replay_receipt": (
        FAILED_V19_ROOT / "m2v5_runner_only_gate_replay.recovery.v19.json"
    ),
    "replay_log": (
        FAILED_V19_ROOT / "m2v5_runner_only_gate_replay.recovery.v19.log"
    ),
    "replay_success_witness": (
        FAILED_V19_ROOT
        / "m2v5_runner_only_gate_replay.recovery.v19.success_witness.json"
    ),
    "formal_continuation_stderr": (
        FAILED_V19_ROOT / "formal_continuation_stderr.log"
    ),
    "diagnostic_strace_write_syscalls": (
        FAILED_V19_ROOT / "diagnostic_strace_write_syscalls.log"
    ),
    "review_mendel": (
        FAILED_V19_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v19.json"
    ),
    "review_nash": (
        FAILED_V19_ROOT / "20260723_tumor_ref_schema_recovery_nash.v19.json"
    ),
    "review_external_claude_opus": (
        FAILED_V19_ROOT
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v19.json"
    ),
}
FAILED_V19_ORIGINAL_REVIEW_PATHS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v19.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v19.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v19.json"
    ),
}
FAILED_V19_ARCHIVED_REVIEW_KEYS = {
    "mendel": "review_mendel",
    "nash": "review_nash",
    "external_claude_opus": "review_external_claude_opus",
}
FAILED_V19_EXPECTED_REVIEWERS = {
    "mendel": ("Mendel", "019f8ea2-2b9a-7271-90b2-0725bcc354a6"),
    "nash": ("Nash", "019f8ea2-261f-7520-8a18-0d692b5002e0"),
    "external_claude_opus": (
        "External Claude Opus",
        "db7885d1-6d8c-4bdb-8fcc-7bb3f6f03124",
    ),
}
FAILED_V19_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v19.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v19.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v19.json.ed25519.sig",
    *FAILED_V19_ORIGINAL_REVIEW_PATHS.values(),
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v19.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v19.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v19.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v19.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v19.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v19.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v19.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v19.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v19.json",
)
FAILED_V19_CHILD_ERROR_MESSAGE = (
    "inotify file watch failed for "
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
    "20260717_all_ssnv_v3/"
    "ed25519_private_encrypted_unrecoverable_after_signing.pem: Permission denied"
)
FAILED_V19_CHILD_STDERR_SHA256 = (
    "47b76169652069e2cccdcae810a8b1cf9f4eb3bc34df804fde03f0d3030d0c33"
)
REJECTED_V20_ROOT = (
    AUDIT_ROOT
    / "rejected_pre_authority_reviews"
    / "20260723_v20_candidate_scope_gate_and_setup_toctou_findings"
)
REJECTED_V20_EVIDENCE = REJECTED_V20_ROOT / "rejection_evidence.json"
EXPECTED_REJECTED_V20_EVIDENCE_SHA256 = (
    "8b851d38304df5aa6f350bd934d5d82c7eacc287faaa3aee9c28cd73a75b0a34"
)
REJECTED_V20_SOURCE_PATHS = {
    "authority_validator": (
        REJECTED_V20_ROOT / "validate_tumor_ref_schema_recovery_authority_v20.py"
    ),
    "ceremony_builder": (
        REJECTED_V20_ROOT / "build_tumor_ref_schema_recovery_authority_v20.py"
    ),
    "continuation_verifier": (
        REJECTED_V20_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v20.py"
    ),
    "downstream_continuation": (
        REJECTED_V20_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v20.py"
    ),
    "readonly_probe": (
        REJECTED_V20_ROOT / "probe_tumor_ref_schema_recovery_sources_v20.py"
    ),
    "regression_tests": REJECTED_V20_ROOT / "test_tumor_ref_schema_recovery_v20.py",
    "runner_gate_replay": (
        REJECTED_V20_ROOT / "replay_m2v5_runner_only_gates_recovery_v20.py"
    ),
}
EXPECTED_REJECTED_V20_SOURCE_IDENTITIES = {
    "authority_validator": (
        341_929,
        "4e34624e5054b683da9946cf7b843ced01c9fb4f9aab7800df34c2bddd080731",
    ),
    "ceremony_builder": (
        53_672,
        "b8f7b47fa0937e39dd071162c2a0b34c62ba33940488364c600772dbdf7ed080",
    ),
    "continuation_verifier": (
        127_906,
        "7805b8e14982031d391ddc98d51736b39bb47a9a77f95150afb66bc7a6dd4ce7",
    ),
    "downstream_continuation": (
        366_859,
        "687e92e6d0ed1c4ee8777435b558da03de914b280d7dcbf4483ad66e8d8b1fb4",
    ),
    "readonly_probe": (
        86_719,
        "432fd8f676f89456a753d5bed86388c1dbfec4aa8d9e9e9c39c94519a0b46a2f",
    ),
    "regression_tests": (
        151_124,
        "f970de837c22c51e068ebd35996d4dbf09ec59eedc5f9676c0c7e5e66a85f3f5",
    ),
    "runner_gate_replay": (
        152_178,
        "f1de1c5e9a4893a47d401ff1b4c41d469e5c686b3ec80286ad0a7afb1fbc1e51",
    ),
}
REJECTED_V20_INTERNAL_REVIEW_PATHS = {
    "mendel": REJECTED_V20_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v20.json",
    "nash": REJECTED_V20_ROOT / "20260723_tumor_ref_schema_recovery_nash.v20.json",
}
REJECTED_V20_EXTERNAL_ENVELOPE = (
    REJECTED_V20_ROOT
    / "20260723_external_claude_schema_recovery_review_v20.claude_cli.envelope.json"
)
REJECTED_V20_EXTERNAL_STDERR = (
    REJECTED_V20_ROOT
    / "20260723_external_claude_schema_recovery_review_v20.claude_cli.stderr.txt"
)
REJECTED_V20_TRANSPORT_PATHS = {
    "external_prompt": (
        REJECTED_V20_ROOT / "20260723_external_claude_schema_recovery_review_v20_prompt.md"
    ),
    "external_schema": (
        REJECTED_V20_ROOT / "20260723_external_claude_schema_recovery_review_v20.schema.json"
    ),
    "external_runner": (
        REJECTED_V20_ROOT / "run_external_claude_schema_recovery_review_v20.py"
    ),
}
REJECTED_V20_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v20/ed25519_public.pem"
)
REJECTED_V20_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v20/ed25519_private_one_time.pem"
)
REJECTED_V20_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v10/ed25519_public.pem"
)
REJECTED_V20_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v10/ed25519_private_one_time_resident.pem"
)
REJECTED_V20_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v20.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v20.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v20.json.ed25519.sig",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v20.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v20.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v20.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v20.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v20.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v20.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v20.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v20.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v20.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v20.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v20.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v20.json",
)
FAILED_V21_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v21_signed_authority_c_pre_audit_canonical_launcher_relation_mismatch"
)
FAILED_V21_BUNDLE = (
    FAILED_V21_ROOT / "tumor_ref_promotion_schema_recovery_authority.v21.bundle"
)
FAILED_V21_EVIDENCE = FAILED_V21_ROOT / "failure_evidence.json"
EXPECTED_FAILED_V21_EVIDENCE_SHA256 = (
    "1180a6db140929298e97ebfd84aa3e1a1efece07b1ceab3d123056c99dcd7109"
)
FAILED_V21_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v21/ed25519_public.pem"
)
FAILED_V21_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v21/ed25519_private_one_time.pem"
)
FAILED_V21_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v11/ed25519_public.pem"
)
FAILED_V21_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v11/ed25519_private_one_time_resident.pem"
)
FAILED_V21_SOURCE_PATHS = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v21.py"
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v21.py"
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v21.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v21.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v21.py",
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v21.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v21.py"
    ),
}
EXPECTED_FAILED_V21_SOURCE_SHA256 = {
    "authority_validator": "b0761a96d1a96a32805e0a14412b5afb9bb403fb06bc300a9f507b73c3afeb75",
    "ceremony_builder": "f6cbbd41f191114f1b8c7baa79156a05f33cf82fbdd83053bfde06ff39003074",
    "continuation_verifier": "e95d646da009ad4457c0c0f71d74f57fd207208d4cf541dec90e3476e729cde6",
    "downstream_continuation": "82b2b49e1fd3675fb962735ceb3c041ef7f16ad29e8eac3d4027674789bc50eb",
    "readonly_probe": "8158a737279078fe6192dd6cf599dbfb16b3e8f068f8018689b924cd01eb3d9a",
    "regression_tests": "43026ee29f2741d7dc0edc0fc54bbdd759713e10a5710422fbd28c34adc08dc6",
    "runner_gate_replay": "4f181e07f9ae2175d3cbb92e80e8f17478d820941fa9f3518c793bbd75007a08",
}
FAILED_V21_ARCHIVED_ARTIFACT_PATHS = {
    "authority": FAILED_V21_BUNDLE / "authority.json",
    "authority_signature": FAILED_V21_BUNDLE / "authority.ed25519.sig",
    "authority_commit": FAILED_V21_BUNDLE / "commit.json",
    "verification_receipt": (
        FAILED_V21_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v21.json"
    ),
    "replay_receipt": (
        FAILED_V21_ROOT / "m2v5_runner_only_gate_replay.recovery.v21.json"
    ),
    "replay_log": FAILED_V21_ROOT / "m2v5_runner_only_gate_replay.recovery.v21.log",
    "replay_success_witness": (
        FAILED_V21_ROOT
        / "m2v5_runner_only_gate_replay.recovery.v21.success_witness.json"
    ),
    "continuation_incident": (
        FAILED_V21_ROOT / "m2v5_downstream_continuation_incident.recovery.v21.json"
    ),
    "review_mendel": (
        FAILED_V21_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v21.json"
    ),
    "review_nash": (
        FAILED_V21_ROOT / "20260723_tumor_ref_schema_recovery_nash.v21.json"
    ),
    "review_external_claude_opus": (
        FAILED_V21_ROOT
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v21.json"
    ),
    "strict_not_applicable": (
        FAILED_V21_ROOT
        / "observed_downstream_outputs"
        / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5"
        / "not_applicable_receipt.json"
    ),
    "matched_normal_not_applicable": (
        FAILED_V21_ROOT
        / "observed_downstream_outputs"
        / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5"
        / "not_applicable_receipt.json"
    ),
    "cn_ccf_receipt": (
        FAILED_V21_ROOT
        / "observed_downstream_outputs"
        / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5"
        / "receipt.json"
    ),
    "cn_ccf_tsv": (
        FAILED_V21_ROOT
        / "observed_downstream_outputs"
        / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5"
        / "candidate_cn_ccf_annotations.tsv.gz"
    ),
    "primary_post_audit": (
        FAILED_V21_ROOT
        / "observed_downstream_outputs"
        / "stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json"
    ),
    "frozen_post_audit": (
        FAILED_V21_ROOT
        / "observed_downstream_outputs"
        / "frozen_input_immutability.post_m2v5_downstream_v3_source_authority_v5.json"
    ),
    "continuation_log": (
        FAILED_V21_ROOT
        / "observed_downstream_outputs"
        / "m2v5_recovered_completion_chain_v2_source_authority_v5.log"
    ),
}
EXPECTED_FAILED_V21_ARCHIVED_SHA256 = {
    "authority": "895aaa71b625dd3b7ce193daa97061ff5b3acb7ef1eccdab8f053247b7ce2acc",
    "authority_signature": "729cfc146c44253e0eb621f5ca4cc1995402e238d32d341b91a3f2efda4d8e42",
    "authority_commit": "b8d4295c198122a5fec2d040bf67ac891a0a45e83550bb87dfcd50d5c419e714",
    "verification_receipt": "c61e931b02dbf3fa3a1aaf39b2b7de4f0d95352af9896843eb1611dc4c818869",
    "replay_receipt": "ae4594eee049247a8a51d1d5a8d5000e34640f048ad23d5d33ed80cba13dfb1c",
    "replay_log": "90262907e8e0f548fa2feb4c2e0c8bbfc9ec04272e0b2e907bcc7111e7f6a62f",
    "replay_success_witness": "203452a7e715cb7b61a5fea5928c947b650f905d2995218905644e0d2f6904ee",
    "continuation_incident": "1712d805d5f618aed9400084f17305dbe9d188be56c1109dd424e171f47f3c44",
    "review_mendel": "263c5a6645a714159787e6d92259ad5448bda7cbd829e8cb83fadacc21a20557",
    "review_nash": "04e401cd775d5397645e2eb091bf43ae4ae710c248005892ef1113fdab15f1f7",
    "review_external_claude_opus": "497cb4fbcb30e7bd30a33dd293fd61c6d9bd91ca6e48b7e2f195a393ecb927b9",
    "strict_not_applicable": "6d630e1dfb98a4c5258ac2ac62d4c12d70f864aa7a87834ef65c80fdfb230ab4",
    "matched_normal_not_applicable": "7e8d7a158e4d303c248af248c18b6af9741202a98d94b93f902c08c478dda268",
    "cn_ccf_receipt": "13e53a7e914e7486f9ccd70ac52fbd2f0b888da4740207b2805d4e325093f3c0",
    "cn_ccf_tsv": "3af165f7d000ad89e37e04e021f632f42aeba80909660c72e1911f496835d5db",
    "primary_post_audit": "0a778de4578cc2d25e072341b421e850cdc97d4eabbc3bb466cb17dec605231a",
    "frozen_post_audit": "8d41d3935745f81f4d9ac32b36f6b34df2e4eefc3f248a34055abbe993109226",
    "continuation_log": "4054b06eb7d3428ba36c61ee95199261b777dfa6f8f62ae909849a1dd5905dbe",
}
FAILED_V21_REVIEW_TRANSPORT_PATHS = {
    "external_envelope": (
        FAILED_V21_ROOT
        / "review_transport"
        / "20260723_external_claude_schema_recovery_review_v21.claude_cli.envelope.json"
    ),
    "external_stderr": (
        FAILED_V21_ROOT
        / "review_transport"
        / "20260723_external_claude_schema_recovery_review_v21.claude_cli.stderr.txt"
    ),
    "external_schema": (
        FAILED_V21_ROOT
        / "review_transport"
        / "20260723_external_claude_schema_recovery_review_v21.schema.json"
    ),
    "external_prompt": (
        FAILED_V21_ROOT
        / "review_transport"
        / "20260723_external_claude_schema_recovery_review_v21_prompt.md"
    ),
    "external_runner": (
        FAILED_V21_ROOT
        / "review_transport"
        / "run_external_claude_schema_recovery_review_v21.py"
    ),
}
EXPECTED_FAILED_V21_TRANSPORT_SHA256 = {
    "external_envelope": "132dee162a9764376bd3fed5181e946df8b49079d22af1e5ef5ab4b63791df88",
    "external_stderr": "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a",
    "external_schema": "c968d34429743d4ab30254573560b74e14b9ab8db9876a3bd8b5eed556645c06",
    "external_prompt": "28f1e84b5d70dfff17b019ec70cb49ef2235106ef37c2495c2fe8f48ede5a1cf",
    "external_runner": "26252eaf4fdd2802f62d6dbb86d7c255e3a7db4220c5c6b9ce1a9ae68dbd4ddf",
}
FAILED_V21_ORIGINAL_REVIEW_PATHS = {
    "mendel": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v21.json",
    "nash": REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v21.json",
    "external_claude_opus": (
        REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v21.json"
    ),
}
FAILED_V21_EXPECTED_REVIEWERS = {
    "mendel": ("Mendel", "019f8f0e-b59e-7811-ada3-b540856095af"),
    "nash": ("Nash", "019f8f0e-7eb7-7f43-b0df-80d451fe37cf"),
    "external_claude_opus": (
        "External Claude Opus",
        "ab3ab34d-76f2-4da6-ad2a-af5fb88aee2d",
    ),
}
FAILED_V21_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v21.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v21.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v21.json.ed25519.sig",
    *FAILED_V21_ORIGINAL_REVIEW_PATHS.values(),
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v21.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v21.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v21.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v21.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v21.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v21.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v21.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v21.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v21.json",
)
FAILED_V21_ORIGINAL_OBSERVED_OUTPUT_SLOTS = (
    WORKSPACE_ROOT / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT
    / "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5",
    RESULT_ROOT / "stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json",
    RESULT_ROOT / "frozen_input_immutability.post_m2v5_downstream_v3_source_authority_v5.json",
    WORKSPACE_ROOT / "m2v5_recovered_completion_chain_v2_source_authority_v5.log",
    WORKSPACE_ROOT / ".python_cache_m2v5_completion_v2_bound_bootstrap",
)
REJECTED_V22_ROOT = (
    AUDIT_ROOT
    / "rejected_pre_authority_reviews"
    / "20260723_v22_candidate_builder_self_binding_review_provenance_v5_watch_findings"
)
REJECTED_V22_EVIDENCE = REJECTED_V22_ROOT / "rejection_evidence.json"
EXPECTED_REJECTED_V22_EVIDENCE_SHA256 = (
    "62d1690831d8279aecb717d4bd6f7058708ea6cbcae23dfed556064ae7cf5cf4"
)
REJECTED_V22_SOURCE_PATHS = {
    "authority_validator": (
        REJECTED_V22_ROOT / "validate_tumor_ref_schema_recovery_authority_v22.py"
    ),
    "ceremony_builder": (
        REJECTED_V22_ROOT / "build_tumor_ref_schema_recovery_authority_v22.py"
    ),
    "continuation_verifier": (
        REJECTED_V22_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v22.py"
    ),
    "downstream_continuation": (
        REJECTED_V22_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v22.py"
    ),
    "readonly_probe": (
        REJECTED_V22_ROOT / "probe_tumor_ref_schema_recovery_sources_v22.py"
    ),
    "regression_tests": REJECTED_V22_ROOT / "test_tumor_ref_schema_recovery_v22.py",
    "runner_gate_replay": (
        REJECTED_V22_ROOT / "replay_m2v5_runner_only_gates_recovery_v22.py"
    ),
}
EXPECTED_REJECTED_V22_SOURCE_IDENTITIES = {
    "authority_validator": (396_507, "803ce84531c4007a4cf1d69d09e79a789ce557ad1fd5fc0c92dcee0dee3b5119"),
    "ceremony_builder": (53_672, "56b8ee1c8e529aec2245a83cad2734bc8e4d05e5c4c02f6c35ec3fdfc4cb7d31"),
    "continuation_verifier": (129_383, "af6557b49214aab8177887e01424f6c54f2c517c289719d2b9f8a32928d12ed1"),
    "downstream_continuation": (380_779, "d1f4f375677826ba65117d3ae811dad65f736196aa6c23f8a137045e8d1fe9f1"),
    "readonly_probe": (92_975, "1ac085cdf1afbdb177050747019405794449933076194c6a624a0b11e532259c"),
    "regression_tests": (160_411, "2f58d243cc77343501efd64e3fa60f681c2090e22456778b59fdd01031011f7f"),
    "runner_gate_replay": (153_508, "2f891a82553dc34bd8e241f33a13ca6f36c29213ab065678ca1efa05a8599916"),
}
REJECTED_V22_REVIEW_PATHS = {
    "mendel": REJECTED_V22_ROOT / "mendel_request_changes.json",
    "nash": REJECTED_V22_ROOT / "nash_request_changes.json",
    "external_claude_opus": (
        REJECTED_V22_ROOT
        / "20260723_external_claude_schema_recovery_review_v22.claude_cli.envelope.json"
    ),
}
REJECTED_V22_TRANSPORT_PATHS = {
    "external_prompt": (
        REJECTED_V22_ROOT / "20260723_external_claude_schema_recovery_review_v22_prompt.md"
    ),
    "external_schema": (
        REJECTED_V22_ROOT / "20260723_external_claude_schema_recovery_review_v22.schema.json"
    ),
    "external_runner": (
        REJECTED_V22_ROOT / "run_external_claude_schema_recovery_review_v22.py"
    ),
    "external_stderr": (
        REJECTED_V22_ROOT
        / "20260723_external_claude_schema_recovery_review_v22.claude_cli.stderr.txt"
    ),
}
REJECTED_V22_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v22/ed25519_public.pem"
)
REJECTED_V22_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v22/ed25519_private_one_time.pem"
)
REJECTED_V22_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v12/ed25519_public.pem"
)
REJECTED_V22_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v12/ed25519_private_one_time_resident.pem"
)
REJECTED_V22_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v22.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v22.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v22.json.ed25519.sig",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v22.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v22.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v22.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v22.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v22.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v22.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v22.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v22.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v22.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v22.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v22.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v22.json",
)
REJECTED_V23_ROOT = (
    AUDIT_ROOT
    / "rejected_pre_authority_reviews"
    / "20260723_v23_nonexistent_node_runtime_path"
)
REJECTED_V23_EVIDENCE = REJECTED_V23_ROOT / "rejection_evidence.json"
EXPECTED_REJECTED_V23_EVIDENCE_SHA256 = (
    "2c6625c0bc3c711e6fba690f91bbffd50c30eb5a885c6b6b6f0a7f1b44a01858"
)
REJECTED_V23_SOURCE_PATHS = {
    "authority_validator": (
        REJECTED_V23_ROOT / "validate_tumor_ref_schema_recovery_authority_v23.py"
    ),
    "ceremony_builder": (
        REJECTED_V23_ROOT / "build_tumor_ref_schema_recovery_authority_v23.py"
    ),
    "continuation_verifier": (
        REJECTED_V23_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v23.py"
    ),
    "downstream_continuation": (
        REJECTED_V23_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v23.py"
    ),
    "readonly_probe": (
        REJECTED_V23_ROOT / "probe_tumor_ref_schema_recovery_sources_v23.py"
    ),
    "regression_tests": REJECTED_V23_ROOT / "test_tumor_ref_schema_recovery_v23.py",
    "runner_gate_replay": (
        REJECTED_V23_ROOT / "replay_m2v5_runner_only_gates_recovery_v23.py"
    ),
}
EXPECTED_REJECTED_V23_SOURCE_IDENTITIES = {
    "authority_validator": (424_582, "9232148c6dadf015658a14312c44f2c1893aa4ce13c54aacab273c9844ccabe1"),
    "ceremony_builder": (60_231, "f2f97256f1e6419fe98dc5908181765885c679e5fddd06e830d3876fa528e115"),
    "continuation_verifier": (129_383, "7c9138d503eb58e4bcb8eb264e6acd28bb8cd02dc2790267c62c306ce8edfdfd"),
    "downstream_continuation": (383_965, "327934cbd6827566a32ffa6e54ea44b584140677c40ec9dd596f41737c2f5b9c"),
    "readonly_probe": (105_946, "8ea403f2bc0ce60b02614c0af02e80cd57e78329a4a83410bdedcc069d494426"),
    "regression_tests": (167_122, "de26903ace2627cecdc5b87938be1fd130e8139025c4bd0746d0aa049cf7a17f"),
    "runner_gate_replay": (153_508, "8e03bb6fdd0b5cc2cd15e9aaaa20daeaa2ecb6baa81fca89c18efbd2baac3a9b"),
}
REJECTED_V23_REVIEW_PATHS = {
    "mendel": REJECTED_V23_ROOT / "mendel_request_changes.json",
    "nash": REJECTED_V23_ROOT / "nash_request_changes.json",
    "external_claude_opus": (
        REJECTED_V23_ROOT
        / "20260723_external_claude_schema_recovery_review_v23.claude_cli.envelope.json"
    ),
}
REJECTED_V23_TRANSPORT_PATHS = {
    "external_prompt": (
        REJECTED_V23_ROOT / "20260723_external_claude_schema_recovery_review_v23_prompt.md"
    ),
    "external_schema": (
        REJECTED_V23_ROOT / "20260723_external_claude_schema_recovery_review_v23.schema.json"
    ),
    "external_runner": (
        REJECTED_V23_ROOT / "run_external_claude_schema_recovery_review_v23.py"
    ),
    "external_stderr": (
        REJECTED_V23_ROOT
        / "20260723_external_claude_schema_recovery_review_v23.claude_cli.stderr.txt"
    ),
}
REJECTED_V23_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v23/ed25519_public.pem"
)
REJECTED_V23_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v23/ed25519_private_one_time.pem"
)
REJECTED_V23_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v13/ed25519_public.pem"
)
REJECTED_V23_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v13/ed25519_private_one_time_resident.pem"
)
REJECTED_V23_OUTPUT_SLOTS = (
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
REJECTED_V24_ROOT = (
    AUDIT_ROOT
    / "rejected_pre_authority_reviews"
    / "20260723_v24_runtime_execution_and_inventory_scope_findings"
)
REJECTED_V24_EVIDENCE = REJECTED_V24_ROOT / "rejection_evidence.json"
EXPECTED_REJECTED_V24_EVIDENCE_SHA256 = (
    "d8ec62cd282190737fd951eafb1eec89d54e9dec36b0af70584fba243e52857b"
)
REJECTED_V24_SOURCE_PATHS = {
    "authority_validator": (
        REJECTED_V24_ROOT / "validate_tumor_ref_schema_recovery_authority_v24.py"
    ),
    "ceremony_builder": (
        REJECTED_V24_ROOT / "build_tumor_ref_schema_recovery_authority_v24.py"
    ),
    "continuation_verifier": (
        REJECTED_V24_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v24.py"
    ),
    "downstream_continuation": (
        REJECTED_V24_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v24.py"
    ),
    "readonly_probe": (
        REJECTED_V24_ROOT / "probe_tumor_ref_schema_recovery_sources_v24.py"
    ),
    "regression_tests": REJECTED_V24_ROOT / "test_tumor_ref_schema_recovery_v24.py",
    "runner_gate_replay": (
        REJECTED_V24_ROOT / "replay_m2v5_runner_only_gates_recovery_v24.py"
    ),
}
EXPECTED_REJECTED_V24_SOURCE_IDENTITIES = {
    "authority_validator": (442_595, "8e2c5b5d4363583082a95cb742dff838e4ee9ac3fb9ff800bfe4dcb4b5387ba2"),
    "ceremony_builder": (60_648, "703623d52452bde8f8ae80cc86dbbb9cce6133e6b4cdf870efdcafedb5f02e9a"),
    "continuation_verifier": (129_383, "eea1927f92d09217e9cb648a791c36f0cb0cb1c9ac4d4cf94c0be4b506fad02d"),
    "downstream_continuation": (388_261, "6968b819e3a2a43f7d3d7aee3d75268279ff31af2ebd1d19a3ea8b3854ab30d2"),
    "readonly_probe": (109_466, "127bf610fe30a4a34fd76c079a0085a263bfed304f446c77f3543425e9a7174c"),
    "regression_tests": (169_789, "6451daf6bf0dfb74834fffe0bc7cea2ce45e53fb2ab834ee4fb0eccf033df9b3"),
    "runner_gate_replay": (153_508, "4640cccade1e8166d13132cb75347f2bb77db857cedaeafe62d917c62b07a3d8"),
}
REJECTED_V24_REVIEW_PATHS = {
    "mendel": REJECTED_V24_ROOT / "mendel_request_changes.json",
    "nash": REJECTED_V24_ROOT / "nash_request_changes.json",
    "external_claude_opus": (
        REJECTED_V24_ROOT
        / "20260723_external_claude_schema_recovery_review_v24.claude_cli.envelope.json"
    ),
}
REJECTED_V24_TRANSPORT_PATHS = {
    "external_prompt": (
        REJECTED_V24_ROOT / "20260723_external_claude_schema_recovery_review_v24_prompt.md"
    ),
    "external_schema": (
        REJECTED_V24_ROOT / "20260723_external_claude_schema_recovery_review_v24.schema.json"
    ),
    "external_runner": (
        REJECTED_V24_ROOT / "run_external_claude_schema_recovery_review_v24.py"
    ),
    "external_stderr": (
        REJECTED_V24_ROOT
        / "20260723_external_claude_schema_recovery_review_v24.claude_cli.stderr.txt"
    ),
}
REJECTED_V24_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v24/ed25519_public.pem"
)
REJECTED_V24_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260723_v24/ed25519_private_one_time.pem"
)
REJECTED_V24_CONTINUATION_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v14/ed25519_public.pem"
)
REJECTED_V24_CONTINUATION_PRIVATE_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v14/ed25519_private_one_time_resident.pem"
)
REJECTED_V24_OUTPUT_SLOTS = (
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
FAILED_V25_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260724_v25_signed_authority_c_tumor_ref_pre_audit_path_relation_mismatch"
)
FAILED_V25_EVIDENCE = FAILED_V25_ROOT / "failure_evidence.json"
EXPECTED_FAILED_V25_EVIDENCE_SHA256 = (
    "b49d2dac5b610ec34f214a54bfd0c159f795b08ab693761a874e09b284e393be"
)
FAILED_V25_BUNDLE = (
    FAILED_V25_ROOT / "tumor_ref_promotion_schema_recovery_authority.v25.bundle"
)
FAILED_V25_SOURCE_PATHS = {
    "authority_validator": AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v25.py",
    "ceremony_builder": AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v25.py",
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v25.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v25.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v25.py",
    "regression_tests": (
        AUDIT_ROOT / "schema_recovery_tests/test_tumor_ref_schema_recovery_v25.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v25.py"
    ),
}
EXPECTED_FAILED_V25_SOURCE_SHA256 = {
    "authority_validator": "58d19165741dfec488fa0603dcbb01118ed0f9addddf982abbdeb16b609ed2cc",
    "ceremony_builder": "9da7c2a2b221450908685270924a1ab65c14a0dda6df096664c9550dd40fde7c",
    "continuation_verifier": "f5e6f00b936b7bd315ce8652c441d121f6ca0d1e5c1d679f8e946efc86f062f0",
    "downstream_continuation": "151d6810a4f3dc27a50ba5e99aa23ad1f297ce9fb9de2053b816f6b55f2e8e66",
    "readonly_probe": "c73f67e2dfae33770165b691702332e2fed73fb497fbfdd381aae056d798e1fb",
    "regression_tests": "d5f753e85269f55477a8680a3b15581219560b41e8ba1e9a0aacb86feb2284e4",
    "runner_gate_replay": "871c675a6c5967b0593212ffb6b5db5749879a7ef7b821f94eaa0f09b2e1ecfb",
}
FAILED_V25_AUTHORITY_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/archive/"
    "20260723_v25_signed_c_tumor_ref_pre_audit_path_relation_failure_01"
)
FAILED_V25_PUBLIC_KEY = FAILED_V25_AUTHORITY_KEY_ARCHIVE / "ed25519_public.pem"
FAILED_V25_PRIVATE_KEY = (
    FAILED_V25_AUTHORITY_KEY_ARCHIVE / "ed25519_private_one_time.pem"
)
FAILED_V25_TERMINAL_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
    "20260723_m2v5_terminal_v15_unused_after_signed_v25_c_failure_01"
)
FAILED_V25_TERMINAL_PUBLIC_KEY = (
    FAILED_V25_TERMINAL_KEY_ARCHIVE / "ed25519_public.pem"
)
FAILED_V25_TERMINAL_PRIVATE_KEY = (
    FAILED_V25_TERMINAL_KEY_ARCHIVE / "ed25519_private_one_time_resident.pem"
)
FAILED_V25_ORIGINAL_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v25.bundle",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v25.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v25.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v25.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v25.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v25.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v25.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v25.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v25.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v25.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v25.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v25.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v25.json",
    *FAILED_V21_ORIGINAL_OBSERVED_OUTPUT_SLOTS,
)
FAILED_V26_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260724_v26_signed_authority_c_signed_promotion_runtime_role_projection_mismatch"
)
FAILED_V26_EVIDENCE = FAILED_V26_ROOT / "failure_evidence.json"
EXPECTED_FAILED_V26_EVIDENCE_SHA256 = (
    "0ee3e2923e6cad144b3decff0c8005fb845af3c1b934ba3cc2256b0751961732"
)
FAILED_V26_BUNDLE = (
    FAILED_V26_ROOT / "tumor_ref_promotion_schema_recovery_authority.v26.bundle"
)
FAILED_V26_SOURCE_PATHS = {
    "authority_validator": AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v26.py",
    "ceremony_builder": AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v26.py",
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v26.py"
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v26.py"
    ),
    "readonly_probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v26.py",
    "regression_tests": (
        AUDIT_ROOT / "schema_recovery_tests/test_tumor_ref_schema_recovery_v26.py"
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v26.py"
    ),
    "recovery_final_dataset_builder": (
        AUDIT_ROOT / "build_all_ssnv_final_report_dataset_schema_recovery_v26.py"
    ),
    "recovery_result_report_finalizer": (
        AUDIT_ROOT / "finalize_task_b_result_release_schema_recovery_v26.py"
    ),
    "recovery_report_builder": (
        AUDIT_ROOT / "build_all_ssnv_report_artifact_schema_recovery_v26.py"
    ),
}
EXPECTED_FAILED_V26_SOURCE_SHA256 = {
    "authority_validator": "ead51f4b3134268f7ef163ede38827d38d1f07383e2d7c44f949d971b2c315ab",
    "ceremony_builder": "4ab964b97b63e2903c9c3e431545942743768b5485fda0fbeaabd686f29cf840",
    "continuation_verifier": "e32a0d57d49af05a123ae0225f8b6a39fd221e9efccfadf96445a62d08e9b885",
    "downstream_continuation": "a79d3630570f523068b1220b20aa78369b4ba7b665672a6b01fcbcdb47eedd54",
    "readonly_probe": "4fcf5df25224d2e18a663bc846957cc4c43913162e70d792cef85d7fbe25b187",
    "regression_tests": "1670c508bb0ba2be3b5e8439b0465b54c8415a240f1d4a5b928a443a76a87f99",
    "runner_gate_replay": "2f19d2cdcc96c65ed96c2d0aacbe7fcc24cc303adc01f70d77b9e3150b72355b",
    "recovery_final_dataset_builder": "cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f",
    "recovery_result_report_finalizer": "ab758306b80c26c17c169285273f99e235ba2c16496cd1b9a2a8fdce974f73d1",
    "recovery_report_builder": "2bf11448db06c7729e94bd28bfee2f2d83269728eedd09dac185a87b10488527",
}
FAILED_V26_REVIEW_PATHS = {
    "mendel": FAILED_V26_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v26.json",
    "nash": FAILED_V26_ROOT / "20260724_tumor_ref_schema_recovery_nash.v26.json",
    "external_claude_opus": (
        FAILED_V26_ROOT
        / "20260724_tumor_ref_schema_recovery_external_claude_opus.v26.json"
    ),
}
EXPECTED_FAILED_V26_REVIEW_SHA256 = {
    "mendel": "bf5a7690fccbc7d82c46c792f3447ffead507ab3ab4269ddd7b105e3ddd9da78",
    "nash": "9b3caa622f146df192261ab6613d6fd35edc13707d72c8c5eb7836c7177c56d3",
    "external_claude_opus": "8881218b0595ed370c9e018ff05263e1e6e4456d49cdaf09a87639ec16fa7603",
}
FAILED_V26_EXECUTION_PATHS = {
    "verification_receipt": (
        FAILED_V26_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v26.json"
    ),
    "replay_receipt": (
        FAILED_V26_ROOT / "m2v5_runner_only_gate_replay.recovery.v26.json"
    ),
    "replay_log": FAILED_V26_ROOT / "m2v5_runner_only_gate_replay.recovery.v26.log",
    "replay_success_witness": (
        FAILED_V26_ROOT
        / "m2v5_runner_only_gate_replay.recovery.v26.success_witness.json"
    ),
}
EXPECTED_FAILED_V26_EXECUTION_SHA256 = {
    "verification_receipt": "b0658144a819c6d19cc1f446f15909155d533886b3a51a7ee06be8b01c59f28e",
    "replay_receipt": "dccabc1cd95a40bea1a7ca52fbfa7ef0496d6fb19341e514eb8e5de50442939e",
    "replay_log": "f0ae55837a7d1f402de7678c598694131d7b1c1f1c160758d2df7dddb08f8c41",
    "replay_success_witness": "6a2764f739c27669152b946a7725c67481b4e8cc609b56a1f33375b4966421b0",
}
FAILED_V26_AUTHORITY_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/archive/"
    "20260724_v26_signed_c_runtime_role_projection_failure_01"
)
FAILED_V26_PUBLIC_KEY = FAILED_V26_AUTHORITY_KEY_ARCHIVE / "ed25519_public.pem"
FAILED_V26_PRIVATE_KEY = (
    FAILED_V26_AUTHORITY_KEY_ARCHIVE / "ed25519_private_one_time.pem"
)
FAILED_V26_AUTHORITY_KEY_ARCHIVE_RECORD = (
    FAILED_V26_AUTHORITY_KEY_ARCHIVE / "FAILED_KEY_ARCHIVE_RECORD.v1.json"
)
FAILED_V26_TERMINAL_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
    "20260724_m2v5_terminal_v16_unused_after_signed_v26_c_failure_01"
)
FAILED_V26_TERMINAL_PUBLIC_KEY = (
    FAILED_V26_TERMINAL_KEY_ARCHIVE / "ed25519_public.pem"
)
FAILED_V26_TERMINAL_PRIVATE_KEY = (
    FAILED_V26_TERMINAL_KEY_ARCHIVE / "ed25519_private_one_time_resident.pem"
)
FAILED_V26_TERMINAL_KEY_ARCHIVE_RECORD = (
    FAILED_V26_TERMINAL_KEY_ARCHIVE / "FAILED_KEY_ARCHIVE_RECORD.v1.json"
)
FAILED_V26_ORIGINAL_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v26.bundle",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v26.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v26.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v26.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v26.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v26.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v26.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v26.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v26.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v26.json",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v26.json",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_nash.v26.json",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_external_claude_opus.v26.json",
)
REJECTED_V27_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260724_v27_pre_authority_runtime_contract_and_review_transport_rejection"
)
REJECTED_V27_EVIDENCE = REJECTED_V27_ROOT / "rejection_evidence.json"
EXPECTED_REJECTED_V27_EVIDENCE_SHA256 = (
    "3c73542487bfa9a5bc5a107c8b6bc50d8e1af19062ab1925f7821b00a0cce604"
)
REJECTED_V27_SOURCE_PATHS = {
    "authority_validator": (
        REJECTED_V27_ROOT / "sources/validate_tumor_ref_schema_recovery_authority_v27.py"
    ),
    "ceremony_builder": (
        REJECTED_V27_ROOT / "sources/build_tumor_ref_schema_recovery_authority_v27.py"
    ),
    "continuation_verifier": (
        REJECTED_V27_ROOT / "sources/verify_tumor_ref_receipt_promotion_recovery_v27.py"
    ),
    "downstream_continuation": (
        REJECTED_V27_ROOT / "sources/continue_m2v5_after_tumor_ref_promotion_recovery_v27.py"
    ),
    "readonly_probe": (
        REJECTED_V27_ROOT / "sources/probe_tumor_ref_schema_recovery_sources_v27.py"
    ),
    "regression_tests": (
        REJECTED_V27_ROOT / "sources/test_tumor_ref_schema_recovery_v27.py"
    ),
    "runner_gate_replay": (
        REJECTED_V27_ROOT / "sources/replay_m2v5_runner_only_gates_recovery_v27.py"
    ),
    "recovery_final_dataset_builder": (
        REJECTED_V27_ROOT
        / "sources/build_all_ssnv_final_report_dataset_schema_recovery_v27.py"
    ),
    "recovery_result_report_finalizer": (
        REJECTED_V27_ROOT / "sources/finalize_task_b_result_release_schema_recovery_v27.py"
    ),
    "recovery_report_builder": (
        REJECTED_V27_ROOT / "sources/build_all_ssnv_report_artifact_schema_recovery_v27.py"
    ),
}
EXPECTED_REJECTED_V27_SOURCE_IDENTITIES = {
    "authority_validator": (497_952, "12d194b97f6dac302d3e5234c858ceff3ff00170e5d99344bc9af8f29f601801"),
    "ceremony_builder": (61_065, "23137d0f174e435c7cdfcd3c9b4b35efd91063104ff06c17ff0df75c6aa8e5c0"),
    "continuation_verifier": (129_377, "5e83ed0add71ba832bdc6c5b3729caba198d1c1bc9d831426526d1059016af25"),
    "downstream_continuation": (403_619, "1005bde75fcbc35fbe53457195bda6e55e01a8df52a3974ab61776c69a170b95"),
    "readonly_probe": (134_418, "8be86c609904d4d891db872a1339b34d33f3dfab3a61002b60d47f8360ea2a75"),
    "regression_tests": (194_893, "a9bb4e654643e82758a6687a6c37f0fb323cc480e0d44c1a9f9c8ee1fbd05891"),
    "runner_gate_replay": (153_502, "5302d3307048fd24df3bd82015b4729248c3226a057c1cf873e6487569ff4aa9"),
    "recovery_final_dataset_builder": (351_388, "cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f"),
    "recovery_result_report_finalizer": (33_500, "35508a5049e5bba2b1e72e8b08a2e9d0f35665cfb9208a082549923d6a2b3fc7"),
    "recovery_report_builder": (238_722, "f16cc51bbec3ee23e7598f8ded3eb50f85e13e1634d5c795c8335024cbfeac32"),
}
REJECTED_V27_REVIEW_PATHS = {
    "mendel": REJECTED_V27_ROOT / "reviews/mendel_request_changes.json",
    "nash": REJECTED_V27_ROOT / "reviews/nash_request_changes.json",
    "external_claude_opus": REJECTED_V27_ROOT / "reviews/external_claude_approve.json",
}
EXPECTED_REJECTED_V27_REVIEW_IDENTITIES = {
    "mendel": (3_949, "50bbe7cfb2b2ca32e402240c698666d4a45462527083e081cc65c8e2525f95fb"),
    "nash": (4_623, "cede19b5565cf64f5569eda2a7319236d9da3e0cdaa109cc87ec6feb853f74b7"),
    "external_claude_opus": (6_413, "6a97e4a693a6b2b3b79b15281b09276d2ae58cf1fb4f87fedadcde32e7d4a54e"),
}
REJECTED_V27_TRANSPORT_PATHS = {
    "external_prompt": REJECTED_V27_ROOT / "review_transport/20260724_external_claude_schema_recovery_review_v27_prompt.md",
    "external_schema": REJECTED_V27_ROOT / "review_transport/20260724_external_claude_schema_recovery_review_v27.schema.json",
    "external_runner": REJECTED_V27_ROOT / "review_transport/run_external_claude_schema_recovery_review_v27.py",
    "review_publisher": REJECTED_V27_ROOT / "review_transport/publish_tumor_ref_schema_recovery_reviews_v27.py",
    "external_envelope": REJECTED_V27_ROOT / "review_transport/20260724_external_claude_schema_recovery_review_v27.claude_cli.envelope.json",
    "external_stderr": REJECTED_V27_ROOT / "review_transport/20260724_external_claude_schema_recovery_review_v27.claude_cli.stderr.txt",
}
EXPECTED_REJECTED_V27_TRANSPORT_IDENTITIES = {
    "external_prompt": (8_652, "4df1602a33fcd8f5842c5b50361a632027514db2332711332b6cf9774a1d4983"),
    "external_schema": (4_332, "a03238bd84fed5008a5429865e154753066d8fb1f303bea500a20fc53442cc54"),
    "external_runner": (4_761, "4bbc029342a767bfe1a61d53d2c6314627c32a9e996c0160f818b54f389fe1ca"),
    "review_publisher": (8_542, "882d3ac0f061381920236b7ab17db09319fda4550c1f1626cf2ff274a9d36a1f"),
    "external_envelope": (15_658, "3a9ca228970554f506420e5ca4fccc0724a19c8c6e2ae0a42a3454b0db5a2f5d"),
    "external_stderr": (302, "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a"),
}
REJECTED_V27_AUTHORITY_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/archive/"
    "20260724_v27_unused_pre_authority_review_rejection_01"
)
REJECTED_V27_ORIGINAL_AUTHORITY_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260724_v27"
)
REJECTED_V27_TERMINAL_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
    "20260724_m2v5_terminal_v17_unused_pre_authority_review_rejection_01"
)
REJECTED_V27_ORIGINAL_TERMINAL_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260724_m2v5_terminal_v17"
)
REJECTED_V27_PUBLIC_KEY = REJECTED_V27_AUTHORITY_KEY_ARCHIVE / "ed25519_public.pem"
REJECTED_V27_PRIVATE_KEY = (
    REJECTED_V27_AUTHORITY_KEY_ARCHIVE / "ed25519_private_one_time.pem"
)
REJECTED_V27_AUTHORITY_KEY_ARCHIVE_RECORD = (
    REJECTED_V27_AUTHORITY_KEY_ARCHIVE / "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
)
REJECTED_V27_CONTINUATION_PUBLIC_KEY = (
    REJECTED_V27_TERMINAL_KEY_ARCHIVE / "ed25519_public.pem"
)
REJECTED_V27_CONTINUATION_PRIVATE_KEY = (
    REJECTED_V27_TERMINAL_KEY_ARCHIVE / "ed25519_private_one_time_resident.pem"
)
REJECTED_V27_TERMINAL_KEY_ARCHIVE_RECORD = (
    REJECTED_V27_TERMINAL_KEY_ARCHIVE / "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
)
REJECTED_V27_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v27.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v27.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v27.json.ed25519.sig",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v27.json",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_nash.v27.json",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_external_claude_opus.v27.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v27.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v27.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v27.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v27.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v27.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v27.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v27.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v27.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v27.json",
)
FAILED_V28_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260724_v28_signed_dataset_c_report_metric_key_order_mismatch"
)
FAILED_V28_EVIDENCE = FAILED_V28_ROOT / "failure_evidence.json"
EXPECTED_FAILED_V28_EVIDENCE_SHA256 = (
    "a49fd9339448f75834e18b71f6b43a257ab99aaf26d833f8e9a2644d4e65fbcb"
)
FAILED_V28_SUMMARY = FAILED_V28_ROOT / "SUMMARY.md"
FAILED_V28_BUNDLE = (
    FAILED_V28_ROOT / "tumor_ref_promotion_schema_recovery_authority.v28.bundle"
)
FAILED_V28_SOURCE_PATHS = {
    "archive_script": AUDIT_ROOT
    / "archive_v28_signed_dataset_c_report_metric_key_order_failure.py",
    "authority_validator": VALIDATOR.with_name(
        "validate_tumor_ref_schema_recovery_authority_v28.py"
    ),
    "ceremony_builder": CEREMONY_BUILDER.with_name(
        "build_tumor_ref_schema_recovery_authority_v28.py"
    ),
    "continuation_verifier": RECOVERY_VERIFIER.with_name(
        "verify_tumor_ref_receipt_promotion_recovery_v28.py"
    ),
    "downstream_continuation": RECOVERY_CONTINUATION.with_name(
        "continue_m2v5_after_tumor_ref_promotion_recovery_v28.py"
    ),
    "readonly_probe": READONLY_PROBE.with_name(
        "probe_tumor_ref_schema_recovery_sources_v28.py"
    ),
    "regression_tests": REGRESSION_TEST.with_name(
        "test_tumor_ref_schema_recovery_v28.py"
    ),
    "runner_gate_replay": RECOVERY_REPLAYER.with_name(
        "replay_m2v5_runner_only_gates_recovery_v28.py"
    ),
    "final_dataset_builder": RECOVERY_FINAL_DATASET_BUILDER.with_name(
        "build_all_ssnv_final_report_dataset_schema_recovery_v28.py"
    ),
    "result_finalizer": RECOVERY_RESULT_REPORT_FINALIZER.with_name(
        "finalize_task_b_result_release_schema_recovery_v28.py"
    ),
    "report_builder": RECOVERY_REPORT_BUILDER.with_name(
        "build_all_ssnv_report_artifact_schema_recovery_v28.py"
    ),
}
FAILED_V28_REVIEW_PATHS = {
    "mendel": FAILED_V28_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v28.json",
    "nash": FAILED_V28_ROOT / "20260724_tumor_ref_schema_recovery_nash.v28.json",
    "external_claude_opus": (
        FAILED_V28_ROOT
        / "20260724_tumor_ref_schema_recovery_external_claude_opus.v28.json"
    ),
}
FAILED_V28_EXECUTION_PATHS = {
    "verification_receipt": (
        FAILED_V28_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v28.json"
    ),
    "replay_receipt": FAILED_V28_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.json",
    "replay_log": FAILED_V28_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.log",
    "replay_success_witness": (
        FAILED_V28_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.success_witness.json"
    ),
    "continuation_incident": (
        FAILED_V28_ROOT / "m2v5_downstream_continuation_incident.recovery.v28.json"
    ),
}
FAILED_V28_PROVISIONAL_DATASET_PATHS = {
    "final_dataset": (
        FAILED_V28_ROOT
        / "observed_downstream_outputs"
        / "all_ssnv_final_report_dataset_v5_m2v5_source_attested"
        / "final_report_dataset.json"
    ),
    "builder_receipt": (
        FAILED_V28_ROOT
        / "observed_downstream_outputs"
        / "all_ssnv_final_report_dataset_v5_m2v5_source_attested"
        / "run_receipt.json"
    ),
    "dataset_release_receipt": (
        FAILED_V28_ROOT
        / "observed_downstream_outputs"
        / "task_b_final_dataset_release_receipt.v1.json"
    ),
    "dataset_release_signature": (
        FAILED_V28_ROOT
        / "observed_downstream_outputs"
        / "task_b_final_dataset_release_receipt.v1.json.ed25519.sig"
    ),
}
FAILED_V28_AUTHORITY_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/archive/"
    "20260724_v28_signed_dataset_c_report_metric_key_order_failure_01"
)
FAILED_V28_TERMINAL_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
    "20260724_m2v5_terminal_v18_unused_after_signed_v28_c_failure_01"
)
FAILED_V28_RESULT_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/archive/"
    "20260724_all_ssnv_result_v5_consumed_v28_c_failure_01"
)
FAILED_V28_REPORT_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/archive/"
    "20260724_all_ssnv_report_v5_unused_v28_c_failure_01"
)
FAILED_V28_KEY_ARCHIVES = {
    "authority_v28": {
        "root": FAILED_V28_AUTHORITY_KEY_ARCHIVE,
        "public": FAILED_V28_AUTHORITY_KEY_ARCHIVE / "ed25519_public.pem",
        "private": FAILED_V28_AUTHORITY_KEY_ARCHIVE / "ed25519_private_one_time.pem",
        "private_mode": 0o000,
        "public_sha256": "0a4d3a99befa388255c506e5a2a77c2acfd831a8b586d7a5b140585015584e3e",
        "record_sha256": "cdfb48e11b5adaf1a76340e6246d1edb73cf96ad9c7591e1401ece368a4a8c49",
        "state": "RETIRED_AFTER_SINGLE_SIGNED_AUTHORITY_THAT_FAILED_AT_C",
    },
    "terminal_v18": {
        "root": FAILED_V28_TERMINAL_KEY_ARCHIVE,
        "public": FAILED_V28_TERMINAL_KEY_ARCHIVE / "ed25519_public.pem",
        "private": (
            FAILED_V28_TERMINAL_KEY_ARCHIVE / "ed25519_private_one_time_resident.pem"
        ),
        "private_mode": 0o400,
        "public_sha256": "ded65f60581476b08c530c04d795c55bcb393558f393acb0f171629dce47add7",
        "record_sha256": "456076174d9322335eba48b3f070601b2e80ef8d77566fc0d84c1810a55c3036",
        "state": "UNUSED_NO_TERMINAL_SIGNATURE_RETIRED_FROM_REUSE",
    },
    "result_v5": {
        "root": FAILED_V28_RESULT_KEY_ARCHIVE,
        "public": FAILED_V28_RESULT_KEY_ARCHIVE / "ed25519_public.pem",
        "private": (
            FAILED_V28_RESULT_KEY_ARCHIVE
            / "ed25519_private_encrypted_unrecoverable_after_signing.pem"
        ),
        "private_mode": 0o000,
        "public_sha256": "5b7d5d026835ec6ec677bcd886bc16ac71444117dabc44a084ce3ede3a4db5a9",
        "record_sha256": "a25b458c0118248dff43fb7fbc77d14aea20b8df76b28305646aee7b917d9d26",
        "state": "CONSUMED_DATASET_SIGNATURE_PROVISIONAL_AFTER_C_FAILURE",
    },
    "report_v5": {
        "root": FAILED_V28_REPORT_KEY_ARCHIVE,
        "public": FAILED_V28_REPORT_KEY_ARCHIVE / "ed25519_public.pem",
        "private": (
            FAILED_V28_REPORT_KEY_ARCHIVE
            / "ed25519_private_encrypted_unrecoverable_after_signing.pem"
        ),
        "private_mode": 0o400,
        "public_sha256": "572e27167e1eea4c39ca53546ba33868b2874ec7fa3ed1682db821b2e50fa439",
        "record_sha256": "5e9d9114233b742626a1d300aaf721bfe3c4a4527717d39db1a930a4fe8f3576",
        "state": "UNUSED_REPORT_KEY_RETIRED_FROM_FAILED_GENERATION",
    },
}
FAILED_V28_ORIGINAL_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v28.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v28.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v28.json.ed25519.sig",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v28.json",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_nash.v28.json",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_external_claude_opus.v28.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v28.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v28.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v28.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v28.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v28.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v28.json",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json.ed25519.sig",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json.ed25519.sig",
)
FAILED_V29_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260724_v29_signed_authority_v_pass_r_archived_key_live_path_mismatch"
)
FAILED_V29_EVIDENCE = FAILED_V29_ROOT / "failure_evidence.json"
EXPECTED_FAILED_V29_EVIDENCE_SHA256 = (
    "4765ee8326b7d209a5b09b9f6676bd5945bd842e558f36eb1464d62e9474bf4c"
)
FAILED_V29_SUMMARY = FAILED_V29_ROOT / "SUMMARY.md"
EXPECTED_FAILED_V29_SUMMARY_SHA256 = (
    "a6013d073738c6b891fcc8bda531ebc76df6c0af745b00f0ce1135d23c76cc39"
)
FAILED_V29_KEY_ARCHIVES = {
    "authority_v29": {
        "root": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "archive/20260724_v29_signed_authority_v_pass_r_archived_key_live_path_failure_01"
        ),
        "private_name": "ed25519_private_one_time.pem",
        "private_mode": 0o000,
        "public_sha256": "819c0ee9d6c6729ad41f3ea1e8071b90ddf506981c68b6927b8abc98e78e9cda",
        "record_sha256": "c7e2c9cce9cd12643fb104c37524beccc99556155ea15b824dccb90ff983ca72",
        "state": "RETIRED_AFTER_SIGNED_V29_AUTHORITY_R_RUNTIME_FAILURE",
    },
    "terminal_v19": {
        "root": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260724_m2v5_terminal_v19_unused_after_signed_v29_r_failure_01"
        ),
        "private_name": "ed25519_private_one_time_resident.pem",
        "private_mode": 0o400,
        "public_sha256": "04d6acab01d56b0bfe25726a242904afd38bbc3ee47d4e3db29e9eb154e23e8b",
        "record_sha256": "66dd4099a57ce9aa94bf0d465fd4dd49dc6e648e5f9287d7737bd2fc1a5311c9",
        "state": "UNUSED_NO_TERMINAL_SIGNATURE_RETIRED_AFTER_V29_R_FAILURE",
    },
    "result_v6": {
        "root": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
            "archive/20260724_all_ssnv_result_v6_unused_v29_r_failure_01"
        ),
        "private_name": "ed25519_private_one_time.pem",
        "private_mode": 0o400,
        "public_sha256": "84438d0a91200108ee06ad7600a3c5428804f37567c351ba33036843ae837864",
        "record_sha256": "22fe8a7bdc18839c425c7b8abe881262ee450a7f6b2284e7e44a164484854e07",
        "state": "UNUSED_NO_SIGNATURE_KEY_UNHELD_RETIRED_AFTER_V29_R_FAILURE",
    },
    "report_v6": {
        "root": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
            "archive/20260724_all_ssnv_report_v6_unused_v29_r_failure_01"
        ),
        "private_name": "ed25519_private_one_time.pem",
        "private_mode": 0o400,
        "public_sha256": "79a684d855ee2d0010691c2a42439389d0e9148d0b84157e04b322c188df6c59",
        "record_sha256": "6d02b12e986ce5b41d2dd18d0051f3425540e4b642ee427e442582fdf147ecfe",
        "state": "UNUSED_NO_SIGNATURE_KEY_UNHELD_RETIRED_AFTER_V29_R_FAILURE",
    },
}
FAILED_V29_ORIGINAL_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v29.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v29.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v29.json.ed25519.sig",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v29.json",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_nash.v29.json",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_external_claude_opus.v29.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v29.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v29.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v29.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v29.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v29.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v29.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v29.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v29.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v29.json",
)
V30_BOOTSTRAP_PREPARE = AUDIT_ROOT / "20260724_v30_four_role_key_bootstrap_prepare.json"
V30_BOOTSTRAP_PROGRESS = AUDIT_ROOT / "20260724_v30_four_role_key_bootstrap_progress.jsonl"
V30_BOOTSTRAP_RECEIPT = AUDIT_ROOT / "20260724_v30_four_role_key_bootstrap_receipt.json"
V30_BOOTSTRAP_SUCCESS = AUDIT_ROOT / "20260724_v30_four_role_key_bootstrap_success.json"
EXPECTED_V30_BOOTSTRAP_SHA256 = {
    "prepare": "3fb0232160a517d1b202462d64a73643e5932ed541845d58cca778240460b91b",
    "progress": "328776bf5124c6e11d15236891e8a3e393bc05b4d90a3ac7f6b7f5d0d1c46343",
    "receipt": "6c86539fd8c15bb3e98299474fcd27ebe426bf06d99366d4ea36f0d9212fe4be",
    "success": "3d4d09eb64f561fa6e2c1d5bd544f8e371aa44eb44a5daf47ac53b69843dd898",
}
V30_RESULT_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
    "20260724_all_ssnv_result_v9_v30_recovery/ed25519_public.pem"
)
V30_REPORT_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
    "20260724_all_ssnv_report_v9_v30_recovery/ed25519_public.pem"
)
EXPECTED_V30_RESULT_PUBLIC_KEY_SHA256 = (
    "0d985d9afce029c06275b6932d51db050f807db50ed27c3bf66cd6f9e201f267"
)
EXPECTED_V30_REPORT_PUBLIC_KEY_SHA256 = (
    "8eaf44b95e216320b45ca4109cf34b2c86b8fd4dcf06bb58025c1e27126bf5b5"
)
TUMOR_REF_SUMMARY_ALIAS = (
    WORKSPACE_ROOT
    / "all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel"
    / "summary.json"
)
TUMOR_REF_SUMMARY_TARGET = (
    WORKSPACE_ROOT
    / "all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel"
    / "all_ssnv_tumor_ref_control_summary.json"
)

CEREMONY_NONGENERATIONAL_FORBIDDEN_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.v3.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.v1.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.v1.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v1.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v1.log",
    RESULT_ROOT / "m2v5_downstream_continuation.v1.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.v1.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.v1.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v1.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v1.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v1.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v1.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_incident.v3.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_authorization.v3.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion.v3.json.ed25519.sig",
    WORKSPACE_ROOT / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v6",
    WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v6",
    WORKSPACE_ROOT
    / "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v6",
    WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v6",
    WORKSPACE_ROOT / "all_ssnv_final_report_dataset_v5_m2v5_source_attested",
    WORKSPACE_ROOT / "all_ssnv_final_report_v5_m2v5_source_attested",
)
ACTIVE_RECOVERY_OUTPUT_SLOTS = (
    AUTHORITY_BUNDLE,
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v30.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v30.json.ed25519.sig",
    RECOVERY_VERIFICATION_RECEIPT,
    RECOVERY_REPLAY_RECEIPT,
    RECOVERY_REPLAY_LOG,
    RECOVERY_REPLAY_SUCCESS_WITNESS,
    RECOVERY_CONTINUATION_RECEIPT,
    RECOVERY_EXIT_ATTESTATION,
    RECOVERY_EXIT_SIGNATURE,
    RECOVERY_SUCCESS_WITNESS,
    RECOVERY_INCIDENT,
)
CEREMONY_FORBIDDEN_OUTPUT_SLOTS = tuple(
    dict.fromkeys(
        (
            *CEREMONY_NONGENERATIONAL_FORBIDDEN_OUTPUT_SLOTS,
            *REJECTED_V2_OUTPUT_SLOTS,
            *REJECTED_V3_OUTPUT_SLOTS,
            *REJECTED_V4_OUTPUT_SLOTS,
            *REJECTED_V5_OUTPUT_SLOTS,
            *REJECTED_V6_OUTPUT_SLOTS,
            *REJECTED_V7_OUTPUT_SLOTS,
            *REJECTED_V8_OUTPUT_SLOTS,
            *FAILED_V9_OUTPUT_SLOTS,
            *FAILED_V10_OUTPUT_SLOTS,
            *FAILED_V11_OUTPUT_SLOTS,
            *FAILED_V12_OUTPUT_SLOTS,
            *REJECTED_V13_OUTPUT_SLOTS,
            *FAILED_V14_OUTPUT_SLOTS,
            *REJECTED_V15_OUTPUT_SLOTS,
            *FAILED_V16_OUTPUT_SLOTS,
            *FAILED_V17_OUTPUT_SLOTS,
            *FAILED_V18_OUTPUT_SLOTS,
            *FAILED_V19_OUTPUT_SLOTS,
            *REJECTED_V20_OUTPUT_SLOTS,
            *FAILED_V21_OUTPUT_SLOTS,
            *FAILED_V21_ORIGINAL_OBSERVED_OUTPUT_SLOTS,
            *REJECTED_V22_OUTPUT_SLOTS,
            *REJECTED_V23_OUTPUT_SLOTS,
            *REJECTED_V24_OUTPUT_SLOTS,
            *FAILED_V25_ORIGINAL_OUTPUT_SLOTS,
            *FAILED_V26_ORIGINAL_OUTPUT_SLOTS,
            *REJECTED_V27_OUTPUT_SLOTS,
            *FAILED_V28_ORIGINAL_OUTPUT_SLOTS,
            *FAILED_V29_ORIGINAL_OUTPUT_SLOTS,
            *ACTIVE_RECOVERY_OUTPUT_SLOTS,
        )
    )
)
CEREMONY_STAGING_PATTERNS = tuple(
    f".tumor_ref_promotion_schema_recovery_authority.v{generation}.bundle.staging.*"
    for generation in range(2, 31)
)
EXPECTED_CEREMONY_FORBIDDEN_OUTPUT_SLOT_COUNT = 439
EXPECTED_REGRESSION_SUMMARY = "733 passed"
CEREMONY_FORBIDDEN_OUTPUT_SLOTS_SHA256 = hashlib.sha256(
    json.dumps(
        sorted(str(path) for path in CEREMONY_FORBIDDEN_OUTPUT_SLOTS),
        ensure_ascii=True,
        separators=(",", ":"),
    ).encode("ascii")
).hexdigest()
CEREMONY_STAGING_PATTERNS_SHA256 = hashlib.sha256(
    json.dumps(
        list(CEREMONY_STAGING_PATTERNS),
        ensure_ascii=True,
        separators=(",", ":"),
    ).encode("ascii")
).hexdigest()

EXPECTED_LEGACY_SHA256 = {
    "continuation_verifier": "03ff3b32368efffafa35491e04621508a46134d36407f060c3da12f90f2432a8",
    "runner_gate_replay": "10f1607aca3ef1a99da7fd77dcd6a207e0ba7003a6e3547a35b28926771fd694",
    "downstream_continuation": "f7b77bd16bd86ae1cbd97e85eebb38a882998b09bc9228fa5b045abfc0ffcfbd",
}
EXPECTED_ORIGINAL_CHAIN_SHA256 = {
    "promotion_authorization": "88cc8b726a5f33a0dbb9776751f4c3e0e6cc2d9a7483731e1fb81d6cbaf9d58e",
    "prepare_exit_attestation": "195c71f6bb5d15e07939511c2df9dc3d0c9d3e3e2793076dce6cfa25771dda75",
    "prepare_exit_attestation_signature": "df0a9209b51ac7ab48aee49417069afb60f373d05a7408060c8e826692343037",
    "promotion_completion": "62945e218753da216edc6d67f67669b8ad0bc986f767379894f4fbfb72459ae8",
    "promote_exit_attestation": "c8018b0fee52999787ddf286f0d3e86e789075c3e67aa7cd936d325c91d2fe75",
    "promote_exit_attestation_signature": "29372ff37faaa60d11777b57e572c591a9fd23a4b6fbf224e6c84a0ef80237ea",
    "source_receipt": "d9fd6ecfc92b1e76ef6b68932faaab1e3d8e8540f521e5a7065f760bdd228e19",
    "canonical_receipt": "d9fd6ecfc92b1e76ef6b68932faaab1e3d8e8540f521e5a7065f760bdd228e19",
}
EXPECTED_PRIOR_RECOVERY_CHAIN_SHA256 = {
    "authority": "af7c2cf68b7edbfe3883cd65b273213985bb9058e1799984860c3a09e1725fda",
    "authority_signature": "ab2026b3ae7fce9c2e602e06dd95760f1d386a5ec90525471600f50ba5b60243",
    "public_key": "1614f15414b3719187f2fcf281794f206e6f12b7c32fd9ab7d532ef167ba1f34",
    "review_mendel": "7341d17279611b2c1a6a97815b68bd2f73fa4fc9d9a56e6f6b802cb1efba9502",
    "review_nash": "3631d4a9c087a64cddf80e9184949fb5ac0b7cf5fd6bb780248c08dc38c1c1e2",
    "review_external_claude_opus": "f51fe4042a40244ae866b78a64438fe77bda19642d8a93a7b4e9f15caad0c556",
    "authority_validator": "1dfb3f3efd62068535788a24bf711cabb1ff425a7b09458dcbd182193cb8ac89",
    "runner_gate_replay": "dee756a0ca55884ca88260b3b5f1af71b397b26952a80ad1813855d03733ed3c",
    "downstream_continuation": "81175713f1a124d4d60347581723bd8ae8f79f3e4cf9539e14edd1354dbab94d",
    "verification_receipt": "b1f9b22d570ca68f91fd3039c2d2f9c1956c1bd769f4b30c983ddf6a122a13a5",
    "runner_failure_evidence": "7af8f793c1a8ca534445e598608ca3edf0d01f237d4130f43537b6f8c3fb7cb0",
}

HISTORICAL_EXECUTION_KEYS = [
    "ctime_ns",
    "device",
    "inode",
    "mode",
    "mtime_ns",
    "path",
    "sha256",
    "size_bytes",
]
CURRENT_IDENTITY_KEYS = [
    "ctime_ns",
    "device",
    "inode",
    "link_count",
    "mode",
    "mtime_ns",
    "path",
    "sha256",
    "size_bytes",
]
RECOVERY_SCOPE = {
    "defect_class": (
        "tumor_ref_primary_audit_lineage_plus_transition_alias_terminal_ceremony_"
        "metadata_and_legacy_stat_relation_schema_recovery"
    ),
    "defect_pointers": [
        (
            "promotion_authorization.evidence.focal_source_identity_transition."
            "during_execution"
        ),
        "replay_runtime_artifacts.python_executable_alias_relation",
        "runner_gate_replay.terminal_receipt_and_success_witness",
        "authority_builder.final_absence_mutation_watch",
        (
            "promotion_verification.schema_recovery_authority.terminal_key_state."
            "legacy_private_key"
        ),
        (
            "promotion_authorization_and_completion.historical_focal_source_relation."
            "exact_legacy_eight_key_stat_schema"
        ),
        (
            "downstream_continuation.PathMutationSentinel."
            "retired_mode000_named_parent_watch_fallback"
        ),
        (
            "downstream_continuation.compose_downstream_script."
            "bound_fd_canonical_argv0_wrapper"
        ),
        (
            "final_dataset.tumor_ref.primary_artifact_audit_pre."
            "immutable_v1_to_source_authorized_v6_lineage"
        ),
        (
            "final_dataset.primary_artifact_audit_window."
            "tumor_ref_legacy_epoch_plus_current_downstream_epoch"
        ),
    ],
    "historical_execution_identity_keys": HISTORICAL_EXECUTION_KEYS,
    "current_identity_keys": CURRENT_IDENTITY_KEYS,
    "historical_link_count_observation": "NOT_RECORDED_NOT_INFERRED",
    "current_link_count_requirement": 1,
    "prior_authority_id": "20260722_tumor_ref_schema_recovery_v1",
    "prior_runner_output_created": False,
    "allowed_changes": [
        "transition_context_aware_historical_relation_validation",
        "exact_six_key_executable_alias_relation_binding",
        "exact_typed_relation_schema_and_per_reader_path_registry",
        "replay_process_lifetime_alias_and_target_inotify_monitor",
        "runner_supervisor_waitpid_normal_exit_zero_success_witness",
        "persistent_final_absence_inotify_watch_through_pre_rename_check",
        "exact_metadata_plus_size_relation_without_private_key_byte_read",
        "exact_legacy_eight_key_stat_relation_without_inferred_link_count",
        "historical_relation_excluded_from_live_identity_registry",
        "exact_tumor_ref_summary_alias_to_canonical_target_binding",
        "signed_portable_asset_and_tumor_ref_all_event_directory_inventory",
        "transient_alias_replace_and_restore_inotify_regression",
        "mode000_eacces_named_parent_watch_with_identity_recheck",
        "bound_fd_python_runtime_with_canonical_argv0",
        "canonical_builder_path_plus_bound_source_fd_reused_through_publish",
        "readonly_probe_and_regression_test_execute_from_bound_source_fds",
        "reviewer_transport_attribution_without_cryptographic_authorship_claim",
        (
            "continuous_watch_for_v21_original_v5_and_rejected_v22_v23_v24_"
            "output_slots"
        ),
        "actual_bound_fd_execution_of_all_four_mandatory_runtimes",
        "missing_and_substituted_mandatory_runtime_fail_closed_regressions",
        "fresh_v30_terminal_v21_key_rotation_without_mutating_legacy_signed_contract",
        "review_digest_binding_of_prior_failed_signed_recovery_aggregate",
        "formal_c_current_release_key_rebase_and_v30_finalizer_binding",
        "exact_tumor_ref_v1_to_v6_primary_audit_data_plane_lineage_validation",
        "split_tumor_ref_legacy_epoch_and_current_downstream_audit_chronology",
        "descriptor_bound_recovery_final_dataset_and_report_runtime_sources",
        "failed_v28_terminal_key_projection_exact_cross_module_binding",
        "rejected_v29_round1_review_correction_and_reproduction_chain_binding",
    ],
    "historical_record_is_live_identity": False,
    "current_record_remains_descriptor_bound": True,
    "historical_ctime_strictly_precedes_current": True,
    "rejected_unsigned_generations": [
        "v2",
        "v3",
        "v4",
        "v5",
        "v6",
        "v7",
        "v8",
        "v13",
        "v15",
        "v19_round1",
        "v20",
        "v22",
        "v23",
        "v24",
        "v27",
        "v29_round1",
    ],
    "rejected_v2_authority_created": False,
    "rejected_v2_private_key_retired": True,
    "rejected_v3_authority_created": False,
    "rejected_v3_private_key_retired": True,
    "rejected_v4_authority_created": False,
    "rejected_v4_private_key_retired": True,
    "rejected_v5_authority_created": False,
    "rejected_v5_private_key_retired": True,
    "rejected_v6_authority_created": False,
    "rejected_v6_private_key_retired": True,
    "rejected_v6_reviews_archived": True,
    "rejected_v7_authority_created": False,
    "rejected_v7_private_key_retired": True,
    "rejected_v8_authority_created": False,
    "rejected_v8_private_key_retired": True,
    "rejected_v13_authority_created": False,
    "rejected_v13_formal_reviews_published": False,
    "rejected_v13_keys_archived_unused_never_signed": True,
    "rejected_v13_outputs_remain_absent": True,
    "rejected_v15_authority_created": False,
    "rejected_v15_formal_reviews_published": False,
    "rejected_v15_keys_archived_unused_never_signed": True,
    "rejected_v15_outputs_remain_absent": True,
    "rejected_v19_round1_evidence_sources_and_reviews_hash_bound": True,
    "rejected_v19_round1_strictest_review_wins": True,
    "rejected_v19_round1_findings_corrected_before_authority": True,
    "tumor_ref_all_event_directory_set_signed_and_revalidated": True,
    "tumor_ref_transient_alias_restore_regression_passes": True,
    "rejected_v19_round1_authority_created": False,
    "rejected_v19_round1_signature_created": False,
    "rejected_v19_round1_strictest_review_wins": True,
    "rejected_v19_round1_findings_corrected_in_active_candidate": True,
    "failed_signed_generations": [
        "v9",
        "v10",
        "v11",
        "v12",
        "v14",
        "v16",
        "v17",
        "v18",
        "v19",
        "v21",
        "v25",
        "v26",
        "v28",
        "v29",
    ],
    "failed_v9_authority_created": True,
    "failed_v9_authority_archived": True,
    "failed_v9_private_key_retired": True,
    "failed_v9_runtime_outputs_created": False,
    "failed_v10_authority_created": True,
    "failed_v10_authority_archived": True,
    "failed_v10_private_key_retired": True,
    "failed_v10_verification_receipt_created": True,
    "failed_v10_replay_receipt_and_log_created": True,
    "failed_v10_continuation_outputs_created": False,
    "failed_v10_canonical_downstream_outputs_created": False,
    "failed_v11_authority_created": True,
    "failed_v11_authority_archived": True,
    "failed_v11_private_key_retired": True,
    "failed_v11_verification_receipt_created": True,
    "failed_v11_replay_receipt_and_log_created": True,
    "failed_v11_continuation_outputs_created": False,
    "failed_v11_canonical_downstream_outputs_created": False,
    "failed_v11_continuation_key_unused_not_retired": True,
    "failed_v12_authority_created": True,
    "failed_v12_authority_archived": True,
    "failed_v12_private_key_retired": True,
    "failed_v12_verification_receipt_created": True,
    "failed_v12_replay_receipt_and_log_created": True,
    "failed_v12_continuation_outputs_created": False,
    "failed_v12_canonical_downstream_outputs_created": False,
    "failed_v12_continuation_key_unused_not_retired": True,
    "failed_v14_authority_created": True,
    "failed_v14_authority_archived": True,
    "failed_v14_private_key_retired": True,
    "failed_v14_verification_receipt_created": True,
    "failed_v14_replay_outputs_created": False,
    "failed_v14_continuation_outputs_created": False,
    "failed_v14_canonical_downstream_outputs_created": False,
    "failed_v14_continuation_key_unused_not_retired": True,
    "failed_v16_authority_created": True,
    "failed_v16_authority_archived": True,
    "failed_v16_private_key_retired": True,
    "failed_v16_verification_receipt_created": True,
    "failed_v16_replay_receipt_and_log_created": True,
    "failed_v16_replay_success_witness_created": True,
    "failed_v16_continuation_child_started": False,
    "failed_v16_continuation_outputs_created": False,
    "failed_v16_canonical_downstream_outputs_created": False,
    "failed_v16_continuation_key_unused_not_retired_must_not_reuse": True,
    "failed_v16_legacy_eight_key_stat_schema_mismatch_reproduced": True,
    "failed_v17_authority_created": True,
    "failed_v17_authority_archived": True,
    "failed_v17_private_key_retired": True,
    "failed_v17_verification_receipt_created": True,
    "failed_v17_replay_receipt_and_log_created": True,
    "failed_v17_replay_success_witness_created": True,
    "failed_v17_continuation_child_started": True,
    "failed_v17_fresh_verifier_passed": True,
    "failed_v17_continuation_outputs_created": False,
    "failed_v17_canonical_downstream_outputs_created": False,
    "failed_v17_continuation_key_unused_not_retired_must_not_reuse": True,
    "failed_v17_metadata_enrichment_conflict_reproduced": True,
    "failed_v18_authority_created": True,
    "failed_v18_authority_archived": True,
    "failed_v18_private_key_retired": True,
    "failed_v18_verification_receipt_created": True,
    "failed_v18_replay_receipt_and_log_created": True,
    "failed_v18_replay_success_witness_created": True,
    "failed_v18_continuation_child_started": True,
    "failed_v18_fresh_verifier_passed": True,
    "failed_v18_continuation_outputs_created": False,
    "failed_v18_canonical_downstream_outputs_created": False,
    "failed_v18_continuation_key_unused_not_retired_must_not_reuse": True,
    "failed_v18_tumor_ref_summary_alias_noncanonical_reproduced": True,
    "failed_v19_authority_created": True,
    "failed_v19_authority_archived": True,
    "failed_v19_private_key_retired": True,
    "failed_v19_verification_receipt_created": True,
    "failed_v19_replay_receipt_and_log_created": True,
    "failed_v19_replay_success_witness_created": True,
    "failed_v19_continuation_child_started": True,
    "failed_v19_fresh_verifier_passed": True,
    "failed_v19_continuation_outputs_created": False,
    "failed_v19_canonical_downstream_outputs_created": False,
    "failed_v19_continuation_key_unused_not_retired_must_not_reuse": True,
    "failed_v19_mode000_inotify_eacces_reproduced": True,
    "failed_v21_authority_created": True,
    "failed_v21_authority_archived": True,
    "failed_v21_private_key_retired": True,
    "failed_v21_verification_receipt_created": True,
    "failed_v21_replay_receipt_and_log_created": True,
    "failed_v21_replay_success_witness_created": True,
    "failed_v21_continuation_child_started": True,
    "failed_v21_fresh_verifier_passed": True,
    "failed_v21_intermediate_downstream_outputs_archived": True,
    "failed_v21_final_dataset_created": False,
    "failed_v21_continuation_outputs_created": False,
    "failed_v21_continuation_key_unused_not_retired_must_not_reuse": True,
    "failed_v21_canonical_launcher_relation_mismatch_reproduced": True,
    "failed_v25_authority_created": True,
    "failed_v25_authority_archived": True,
    "failed_v25_private_key_retired": True,
    "failed_v25_verification_receipt_created": True,
    "failed_v25_replay_receipt_and_log_created": True,
    "failed_v25_replay_success_witness_created": True,
    "failed_v25_continuation_child_started": True,
    "failed_v25_final_dataset_created": False,
    "failed_v25_continuation_outputs_created": False,
    "failed_v25_terminal_key_archived_unused_must_not_reuse": True,
    "failed_v25_tumor_ref_v1_v6_audit_path_mismatch_reproduced": True,
    "failed_v26_authority_created": True,
    "failed_v26_authority_archived": True,
    "failed_v26_private_key_retired": True,
    "failed_v26_verification_receipt_created": True,
    "failed_v26_replay_receipt_and_log_created": True,
    "failed_v26_replay_success_witness_created": True,
    "failed_v26_continuation_child_started": True,
    "failed_v26_downstream_producers_started": False,
    "failed_v26_final_dataset_created": False,
    "failed_v26_continuation_outputs_created": False,
    "failed_v26_terminal_key_archived_unused_must_not_reuse": True,
    "failed_v26_signed_runtime_projection_mismatch_reproduced": True,
    "rejected_v27_authority_created": False,
    "rejected_v27_formal_reviews_published": False,
    "rejected_v27_strictest_reproducible_finding_wins": True,
    "rejected_v27_authority_and_terminal_keys_archived_unused_must_not_reuse": True,
    "rejected_v27_outputs_remain_absent": True,
    "rejected_v27_runtime_hash_transport_and_staging_findings_corrected": True,
    "failed_v28_authority_created": True,
    "failed_v28_authority_archived": True,
    "failed_v28_authority_private_key_retired": True,
    "failed_v28_verification_and_replay_passed": True,
    "failed_v28_provisional_dataset_created_and_signed": True,
    "failed_v28_final_report_and_terminal_outputs_created": False,
    "failed_v28_all_four_keys_archived_must_not_reuse": True,
    "failed_v28_json_mapping_order_bug_reproduced": True,
    "rejected_v29_round1_authority_created": False,
    "rejected_v29_round1_formal_reviews_published": False,
    "rejected_v29_round1_any_key_consumed": False,
    "rejected_v29_round1_strictest_reproducible_review_wins": True,
    "rejected_v29_round1_terminal_projection_gap_reproduced": True,
    "rejected_v29_round1_findings_corrected_in_active_candidate": True,
    "failed_v29_authority_created_and_signed": True,
    "failed_v29_verification_passed": True,
    "failed_v29_replay_failed_before_first_output": True,
    "failed_v29_downstream_continuation_started": False,
    "failed_v29_all_four_keys_archived_must_not_reuse": True,
    "failed_v29_archived_key_live_path_gap_reproduced": True,
    "v30_archive_rebase_restricted_to_three_historical_key_paths": True,
    "v30_four_fresh_key_bootstrap_hash_bound": True,
    "rejected_v20_authority_created": False,
    "rejected_v20_formal_reviews_published": False,
    "rejected_v20_strictest_review_wins": True,
    "rejected_v20_terminal_key_unused_not_retired_must_not_reuse": True,
    "rejected_v20_outputs_remain_absent": True,
    "rejected_v20_scope_test_and_setup_toctou_findings_corrected_in_active_candidate": True,
    "rejected_v22_authority_created": False,
    "rejected_v22_formal_reviews_published": False,
    "rejected_v22_strictest_review_wins": True,
    "rejected_v22_authority_and_terminal_keys_unused_must_not_reuse": True,
    "rejected_v22_outputs_remain_absent": True,
    "rejected_v22_findings_corrected_in_active_candidate": True,
    "rejected_v23_authority_created": False,
    "rejected_v23_formal_reviews_published": False,
    "rejected_v23_strictest_review_wins": True,
    "rejected_v23_authority_and_terminal_keys_unused_must_not_reuse": True,
    "rejected_v23_outputs_remain_absent": True,
    "rejected_v23_findings_corrected_in_active_candidate": True,
    "rejected_v24_authority_created": False,
    "rejected_v24_formal_reviews_published": False,
    "rejected_v24_strictest_review_wins": True,
    "rejected_v24_authority_and_terminal_keys_unused_must_not_reuse": True,
    "rejected_v24_outputs_remain_absent": True,
    "rejected_v24_findings_corrected_in_active_candidate": True,
    "legacy_v2_terminal_key_remains_unused_not_retired": True,
    "fresh_v21_terminal_key_required": True,
    "failed_v16_v6_terminal_key_quarantined": True,
    "failed_v17_v7_terminal_key_quarantined": True,
    "failed_v18_v8_terminal_key_quarantined": True,
    "failed_v19_v9_terminal_key_quarantined": True,
    "rejected_v20_v10_terminal_key_quarantined": True,
    "failed_v21_v11_terminal_key_quarantined": True,
    "rejected_v22_v12_terminal_key_quarantined": True,
    "rejected_v23_v13_terminal_key_quarantined": True,
    "rejected_v24_v14_terminal_key_quarantined": True,
    "failed_v25_v15_terminal_key_quarantined": True,
    "failed_v26_v16_terminal_key_quarantined": True,
    "rejected_v27_v17_terminal_key_quarantined": True,
    "failed_v28_v18_terminal_key_quarantined": True,
    "failed_v29_v19_terminal_key_quarantined": True,
    "all_terminal_keys_must_be_pairwise_distinct": True,
    "legacy_signed_terminal_contract_changed": False,
    "historical_authorized_output_slot_count": 17,
    "current_recovery_output_slot_count": 25,
    "signed_crosslink_inventory_contract": (
        "historical_17_authorized_slots_distinct_from_current_25_absence_slots"
    ),
    "ceremony_forbidden_output_slot_count": 439,
    "ceremony_forbidden_output_slots_sha256": CEREMONY_FORBIDDEN_OUTPUT_SLOTS_SHA256,
    "ceremony_staging_pattern_count": 29,
    "ceremony_staging_patterns_sha256": CEREMONY_STAGING_PATTERNS_SHA256,
    "terminal_absence_directories_fd_leased": True,
    "cooperating_writer_lease": "v30_public_key_flock_exclusive",
    "terminal_absence_scan": "two_pass_directory_generation_stable",
    "terminal_absence_event_watch": (
        "inotify_create_move_delete_attrib_overflow_fail_closed"
    ),
    "terminal_absence_parent_tokens": (
        "result_review_workspace_rechecked_in_blocked_signal_critical_section"
    ),
    "terminal_publish_critical_section": (
        "catchable_signals_blocked_stage_member_tokens_then_prebound_renameat2"
    ),
    "readonly_probe_no_output_writes_scope": (
        "protected_recovery_and_downstream_namespaces_only"
    ),
    "continuation_signature_key_lifecycle": (
        "stage_sealed_unlinked_memfd_then_retire_and_verify_disk_private_key_"
        "before_signing_close_memfd_before_signature_publication"
    ),
    "authority_signature_key_lifecycle": (
        "stage_sealed_unlinked_memfd_then_retire_and_verify_disk_private_key_"
        "before_signing_close_memfd_before_durable_staging"
    ),
    "continuation_signature_publication_semantics": (
        "no_replace_signature_link_is_provisional_then_precommit_verifier_passes_"
        "before_no_replace_durable_success_witness_final_verifier_requires_witness"
    ),
    "continuation_success_witness_publication_semantics": (
        "pre_fsynced_unlinked_inode_then_catchable_signals_blocked_final_no_replace_"
        "link_and_immediate_no_return_exit"
    ),
    "continuation_final_verifier_lease": (
        "shared_nonblocking_flock_proves_exclusive_writer_released_and_remains_excluded"
    ),
    "authority_publication": "atomic_directory_rename_noreplace",
    "builder_validator_binding": (
        "canonical_builder_path_and_proc_cmdline_bound_source_fd_plus_single_fd_"
        "validator_hash_compile_exec"
    ),
    "review_attribution_semantics": REVIEW_ATTRIBUTION_SEMANTICS,
    "authority_validator_source_record": "bootstrap_fd_record_reused_without_reopen",
    "post_sign_failure_policy": (
        "staged_signing_key_arms_BaseException_key_retirement_before_disk_retirement"
    ),
    "post_publish_validation": "retained_stage_and_member_fds_without_validator_reopen",
    "pre_rename_staging_recheck": (
        "signature_commit_then_terminal_directory_inode_member_bytes_recheck"
    ),
    "terminal_input_rechecks": ["before_sign", "before_atomic_publish"],
    "scientific_payload_changed": False,
    "canonical_receipt_bytes_changed": False,
    "legacy_signed_artifacts_changed": False,
    "claim_ceiling_changed": False,
}

LEGACY_COMMANDS = {
    "verify": [str(PYTHON), "-I", "-B", str(LEGACY_VERIFIER), "--verify"],
    "verify_and_record": [
        str(PYTHON), "-I", "-B", str(LEGACY_VERIFIER), "--verify-and-record"
    ],
    "runner_gate_replay": [str(PYTHON), "-I", "-B", str(LEGACY_REPLAYER)],
    "downstream_continuation": [
        str(PYTHON), "-I", "-B", str(LEGACY_CONTINUATION), "--supervise-and-sign"
    ],
}
RECOVERY_COMMANDS = {
    "verify": [str(PYTHON), "-I", "-B", str(RECOVERY_VERIFIER), "--verify"],
    "verify_and_record": [
        str(PYTHON), "-I", "-B", str(RECOVERY_VERIFIER), "--verify-and-record"
    ],
    "runner_gate_replay": [str(PYTHON), "-I", "-B", str(RECOVERY_REPLAYER)],
    "downstream_continuation": [
        str(PYTHON),
        "-I",
        "-B",
        str(RECOVERY_CONTINUATION),
        "--supervise-and-sign",
    ],
    "supervised_child": [
        str(PYTHON), "-I", "-B", str(RECOVERY_CONTINUATION), "--supervised-child"
    ],
    "verify_signed_terminal": [
        str(PYTHON),
        "-I",
        "-B",
        str(RECOVERY_CONTINUATION),
        "--verify-signed-terminal",
    ],
    "verify_signed_terminal_prewitness": [
        str(PYTHON),
        "-I",
        "-B",
        str(RECOVERY_CONTINUATION),
        "--verify-signed-terminal-prewitness",
    ],
}
SCHEMAS = {
    "legacy_verification": "intersubmod.tumor_ref_source_receipt_promotion_verification@3.0.0",
    "recovery_verification": (
        "intersubmod.tumor_ref_source_receipt_promotion_verification.recovery@1.0.0"
    ),
    "legacy_replay": "intersubmod.m2v5_runner_only_gate_replay@1.0.0",
    "recovery_replay": "intersubmod.m2v5_runner_only_gate_replay.recovery@1.0.0",
    "legacy_continuation": "intersubmod.m2v5_downstream_continuation@1.0.0",
    "recovery_continuation": (
        "intersubmod.m2v5_downstream_continuation.recovery@1.0.0"
    ),
}
RECEIPT_PATHS = {
    "legacy_verification": str(
        RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.v3.json"
    ),
    "recovery_verification": str(RECOVERY_VERIFICATION_RECEIPT),
    "legacy_replay": str(RESULT_ROOT / "m2v5_runner_only_gate_replay.v1.json"),
    "recovery_replay": str(RECOVERY_REPLAY_RECEIPT),
    "recovery_replay_log": str(RECOVERY_REPLAY_LOG),
    "recovery_replay_success_witness": str(RECOVERY_REPLAY_SUCCESS_WITNESS),
    "legacy_continuation": str(RESULT_ROOT / "m2v5_downstream_continuation.v1.json"),
    "recovery_continuation": str(RECOVERY_CONTINUATION_RECEIPT),
    "recovery_exit_attestation": str(RECOVERY_EXIT_ATTESTATION),
    "recovery_exit_signature": str(RECOVERY_EXIT_SIGNATURE),
    "recovery_success_witness": str(RECOVERY_SUCCESS_WITNESS),
    "recovery_incident": str(RECOVERY_INCIDENT),
}
AUTHORITY_CHECKS = {
    "authority_signature_verified": True,
    "prior_recovery_authority_signature_verified": True,
    "prior_recovery_private_key_remains_retired": True,
    "prior_runner_failure_evidence_hash_bound": True,
    "prior_runner_created_zero_replay_outputs": True,
    "historical_identity_schema_exactly_eight_keys": True,
    "historical_link_count_not_inferred": True,
    "current_identity_still_requires_nine_keys_and_link_count_one": True,
    "historical_transition_record_is_not_rebound_as_live_identity": True,
    "current_transition_record_remains_descriptor_bound": True,
    "historical_ctime_strictly_precedes_current": True,
    "executable_alias_relation_exact_six_key_schema": True,
    "executable_alias_and_target_descriptor_bound": True,
    "executable_alias_target_bytes_mode_verified_before_execution": True,
    "legacy_v2_terminal_key_contract_preserved": True,
    "legacy_v2_terminal_private_key_remains_mode_0400": True,
    "fresh_v21_terminal_key_is_distinct_and_mode_0400_before_signing": True,
    "fresh_v21_terminal_key_only_authorizes_recovery_v30_outputs": True,
    "failed_v29_signed_authority_v_pass_r_failure_hash_bound": True,
    "failed_v29_four_key_archives_hash_bound_and_not_reused": True,
    "v30_four_fresh_key_bootstrap_witnesses_hash_bound": True,
    "v30_four_role_public_keys_pairwise_distinct": True,
    "bound_python_fd_bytes_retain_canonical_argv0": True,
    "metadata_enrichment_upgrades_registry_and_preserves_conflict_detection": True,
    "tumor_ref_summary_alias_and_canonical_target_descriptor_bound": True,
    "tumor_ref_summary_alias_link_text_and_target_identity_verified": True,
    "tumor_ref_directory_transient_mutation_watch_required": True,
    "mode000_retired_key_named_parent_watch_fallback_required": True,
    "rejected_v2_evidence_hash_bound": True,
    "rejected_v2_private_key_retired": True,
    "rejected_v2_outputs_remain_absent": True,
    "rejected_v3_evidence_hash_bound": True,
    "rejected_v3_private_key_retired": True,
    "rejected_v3_outputs_remain_absent": True,
    "rejected_v4_evidence_hash_bound": True,
    "rejected_v4_private_key_retired": True,
    "rejected_v4_outputs_remain_absent": True,
    "rejected_v5_evidence_hash_bound": True,
    "rejected_v5_private_key_retired": True,
    "rejected_v5_outputs_remain_absent": True,
    "rejected_v6_evidence_hash_bound": True,
    "rejected_v6_private_key_retired": True,
    "rejected_v6_outputs_remain_absent": True,
    "rejected_v6_reviews_archived_and_hash_bound": True,
    "rejected_v7_evidence_hash_bound": True,
    "rejected_v7_private_key_retired": True,
    "rejected_v7_outputs_remain_absent": True,
    "rejected_v8_evidence_hash_bound": True,
    "rejected_v8_private_key_retired": True,
    "rejected_v8_outputs_remain_absent": True,
    "rejected_v13_evidence_and_key_archive_ledger_hash_bound": True,
    "rejected_v13_sources_and_review_adjudication_bound": True,
    "rejected_v13_keys_archived_unused_never_signed": True,
    "rejected_v13_outputs_remain_absent": True,
    "rejected_v15_evidence_reviews_and_key_archive_ledger_hash_bound": True,
    "rejected_v15_sources_and_strict_adjudication_bound": True,
    "rejected_v15_keys_archived_unused_never_signed": True,
    "rejected_v15_outputs_remain_absent": True,
    "failed_v9_evidence_hash_bound": True,
    "failed_v9_authority_signature_verified": True,
    "failed_v9_atomic_commit_verified": True,
    "failed_v9_reviews_archived_and_hash_bound": True,
    "failed_v9_private_key_retired": True,
    "failed_v9_original_outputs_remain_absent": True,
    "failed_v10_evidence_hash_bound": True,
    "failed_v10_authority_signature_verified": True,
    "failed_v10_atomic_commit_verified": True,
    "failed_v10_reviews_v_and_r_archived_and_hash_bound": True,
    "failed_v10_private_key_retired": True,
    "failed_v10_original_outputs_remain_absent": True,
    "failed_v11_evidence_hash_bound": True,
    "failed_v11_authority_signature_verified": True,
    "failed_v11_atomic_commit_verified": True,
    "failed_v11_reviews_v_and_r_archived_and_hash_bound": True,
    "failed_v11_private_key_retired": True,
    "failed_v11_continuation_key_unused_not_retired": True,
    "failed_v11_original_outputs_remain_absent": True,
    "failed_v11_basic_to_full_stat_mismatch_reproduced": True,
    "failed_v12_authority_signature_and_atomic_commit_verified": True,
    "failed_v12_source_and_review_artifacts_bound": True,
    "failed_v12_original_outputs_remain_absent": True,
    "failed_v12_executable_alias_relation_mismatch_reproduced": True,
    "failed_v16_evidence_hash_bound": True,
    "failed_v16_authority_signature_and_atomic_commit_verified": True,
    "failed_v16_source_review_v_and_r_artifacts_bound": True,
    "failed_v16_authority_private_key_retired": True,
    "failed_v16_continuation_key_unused_not_retired_must_not_reuse": True,
    "failed_v16_original_outputs_remain_absent": True,
    "failed_v16_legacy_eight_key_stat_mismatch_reproduced": True,
    "failed_v17_evidence_hash_bound": True,
    "failed_v17_authority_signature_and_atomic_commit_verified": True,
    "failed_v17_source_review_v_and_r_artifacts_bound": True,
    "failed_v17_authority_private_key_retired": True,
    "failed_v17_continuation_key_unused_not_retired_must_not_reuse": True,
    "failed_v17_original_outputs_remain_absent": True,
    "failed_v17_metadata_enrichment_conflict_reproduced": True,
    "failed_v17_syscall_trace_evidence_bound": True,
    "failed_v18_evidence_hash_bound": True,
    "failed_v18_authority_signature_and_atomic_commit_verified": True,
    "failed_v18_source_review_v_and_r_artifacts_bound": True,
    "failed_v18_authority_private_key_retired": True,
    "failed_v18_continuation_key_unused_not_retired_must_not_reuse": True,
    "failed_v18_original_outputs_remain_absent": True,
    "failed_v18_tumor_ref_summary_alias_noncanonical_reproduced": True,
    "failed_v18_syscall_trace_evidence_bound": True,
    "failed_v19_evidence_hash_bound": True,
    "failed_v19_authority_signature_and_atomic_commit_verified": True,
    "failed_v19_source_review_v_and_r_artifacts_bound": True,
    "failed_v19_authority_private_key_retired": True,
    "failed_v19_continuation_key_unused_not_retired_must_not_reuse": True,
    "failed_v19_original_outputs_remain_absent": True,
    "failed_v19_mode000_inotify_eacces_reproduced": True,
    "failed_v19_syscall_trace_evidence_bound": True,
    "failed_v21_evidence_hash_bound": True,
    "failed_v21_authority_signature_and_atomic_commit_verified": True,
    "failed_v21_source_review_v_r_incident_and_downstream_artifacts_bound": True,
    "failed_v21_authority_private_key_retired": True,
    "failed_v21_continuation_key_unused_not_retired_must_not_reuse": True,
    "failed_v21_original_outputs_remain_absent": True,
    "failed_v21_canonical_launcher_relation_mismatch_reproduced": True,
    "failed_v25_evidence_hash_bound": True,
    "failed_v25_authority_signature_verified": True,
    "failed_v25_atomic_commit_verified": True,
    "failed_v25_sources_hash_bound": True,
    "failed_v25_original_outputs_remain_absent": True,
    "failed_v25_tumor_ref_v1_v6_audit_path_mismatch_reproduced": True,
    "failed_v26_evidence_hash_bound": True,
    "failed_v26_authority_signature_verified": True,
    "failed_v26_atomic_commit_verified": True,
    "failed_v26_sources_reviews_v_and_r_hash_bound": True,
    "failed_v26_authority_private_key_retired": True,
    "failed_v26_terminal_key_archived_unused_must_not_reuse": True,
    "failed_v26_original_outputs_remain_absent": True,
    "failed_v26_signed_runtime_projection_mismatch_reproduced": True,
    "rejected_v27_evidence_sources_reviews_transport_and_keys_hash_bound": True,
    "rejected_v27_strictest_reproducible_finding_wins": True,
    "rejected_v27_authority_and_terminal_keys_archived_unused_and_quarantined": True,
    "rejected_v27_outputs_remain_absent": True,
    "failed_v28_evidence_and_full_archive_inventory_hash_bound": True,
    "failed_v28_authority_and_dataset_signatures_verified": True,
    "failed_v28_sources_reviews_v_r_incident_and_provisional_dataset_bound": True,
    "failed_v28_all_four_keys_archived_and_quarantined": True,
    "failed_v28_original_outputs_remain_absent": True,
    "failed_v28_json_mapping_order_failure_reproduced": True,
    "rejected_v29_round1_evidence_summary_and_exact_inventory_hash_bound": True,
    "rejected_v29_round1_review_correction_chain_hash_bound": True,
    "rejected_v29_round1_terminal_projection_gap_reproduced": True,
    "rejected_v29_round1_same_version_keys_unconsumed_before_correction": True,
    "rejected_v29_round1_findings_corrected_in_active_candidate": True,
    "historical_signed_runtime_projection_exact_11_roles": True,
    "current_reviewed_runtime_contract_exact_14_roles": True,
    "recovery_only_runtime_roles_excluded_from_historical_signed_projection": True,
    "tumor_ref_v1_v6_primary_audit_data_plane_lineage_exact": True,
    "tumor_ref_split_audit_chronology_fail_closed": True,
    "recovery_final_dataset_and_report_sources_frozen_and_bound": True,
    "rejected_v20_evidence_sources_reviews_and_transport_hash_bound": True,
    "rejected_v20_strictest_review_wins": True,
    "rejected_v20_authority_and_terminal_keys_unused_and_quarantined": True,
    "rejected_v20_outputs_remain_absent": True,
    "rejected_v22_evidence_sources_reviews_transport_and_keys_hash_bound": True,
    "rejected_v22_strictest_review_wins": True,
    "rejected_v22_authority_and_terminal_keys_unused_and_quarantined": True,
    "rejected_v22_outputs_remain_absent": True,
    "rejected_v23_evidence_sources_reviews_transport_and_keys_hash_bound": True,
    "rejected_v23_strictest_review_wins": True,
    "rejected_v23_authority_and_terminal_keys_unused_and_quarantined": True,
    "rejected_v23_outputs_remain_absent": True,
    "rejected_v24_evidence_sources_reviews_transport_and_keys_hash_bound": True,
    "rejected_v24_strictest_review_wins": True,
    "rejected_v24_authority_and_terminal_keys_unused_and_quarantined": True,
    "rejected_v24_outputs_remain_absent": True,
    "review_artifacts_are_transport_attributed_not_reviewer_signed": True,
    "reviewer_cryptographic_authorship_not_claimed": True,
    "v21_original_v5_intermediate_slots_continuously_watched": True,
    "v21_archived_artifacts_exact_keys_and_observed_outputs_crosslinked": True,
    "parent_watch_precedes_leaf_snapshot_and_setup_identity_recheck": True,
    "all_439_forbidden_slots_have_occupied_state_regressions": True,
    "all_29_staging_patterns_have_occupied_state_regressions": True,
    "setup_window_transient_chmod_and_replacement_regressions": True,
    "historical_17_slot_crosslinks_are_distinct_from_current_25_slot_absence_gate": True,
    "cooperating_writer_lock_held_through_post_publish_check": True,
    "forbidden_namespace_two_pass_generation_stable": True,
    "forbidden_namespace_inotify_mutation_watch_quiet": True,
    "terminal_stage_and_member_generations_rechecked": True,
    "all_absence_parent_generations_rechecked_before_renameat2": True,
    "no_python_callback_between_terminal_tokens_and_renameat2": True,
    "continuation_incident_write_requires_writer_lease": True,
    "full_forbidden_namespace_rechecked_before_sign": True,
    "full_forbidden_namespace_rechecked_after_sign": True,
    "full_forbidden_namespace_rechecked_before_key_retirement": True,
    "full_forbidden_namespace_rechecked_after_external_verification": True,
    "full_forbidden_namespace_rechecked_immediately_before_atomic_publish": True,
    "builder_validator_hash_compile_exec_share_one_fd": True,
    "builder_canonical_path_proc_cmdline_and_bound_source_fd_rechecked": True,
    "readonly_probe_executes_bound_probe_bytes_with_canonical_python_argv0": True,
    "pytest_executes_bound_regression_source_bytes": True,
    "all_four_mandatory_gate_runtimes_are_actually_executed_from_bound_fds": True,
    "missing_or_substituted_mandatory_gate_runtime_fails_closed": True,
    "authority_validator_record_reuses_bootstrap_fd_without_reopen": True,
    "signing_attempt_BaseException_path_retires_key": True,
    "authority_private_key_retired_before_signature_bytes_exist": True,
    "authority_signing_uses_sealed_unlinked_memfd": True,
    "authority_signing_memfd_closed_before_durable_staging": True,
    "continuation_private_key_retired_before_terminal_signing": True,
    "continuation_signing_uses_sealed_unlinked_memfd": True,
    "continuation_signature_provisional_until_independent_verification": True,
    "continuation_post_link_failure_incident_blocks_release_authority": True,
    "continuation_durable_success_witness_required_for_release_authority": True,
    "continuation_signature_without_success_witness_fails_closed": True,
    "continuation_success_witness_link_is_final_no_return_commit": True,
    "runner_replay_success_witness_requires_worker_waitpid_normal_exit_zero": True,
    "runner_replay_partial_log_receipt_witness_tuple_fails_closed": True,
    "continuation_final_verifier_requires_released_writer_lease": True,
    "builder_post_publish_does_not_reopen_validator_path": True,
    "all_ceremony_inputs_terminally_rechecked_before_sign": True,
    "all_ceremony_inputs_terminally_rechecked_before_atomic_publish": True,
    "authority_signature_commit_published_as_atomic_bundle": True,
    "staging_identity_and_bytes_rechecked_after_external_verification": True,
    "legacy_signed_chain_unchanged": True,
    "legacy_vrc_sources_unchanged": True,
    "recovery_vrc_sources_frozen_and_bound": True,
    "recovery_outputs_use_v30_paths_without_overwriting_v1_through_v29": True,
    "probe_tests_and_ceremony_builder_frozen_and_bound": True,
    "recovery_receipt_paths_are_distinct": True,
    "three_distinct_transport_attributed_review_artifacts_approve": True,
    "review_json_duplicate_keys_and_nonfinite_values_rejected": True,
    "review_json_float_literals_rejected": True,
    "review_json_nested_scalar_types_are_exact": True,
    "retired_private_key_fd_leased_and_terminally_rechecked": True,
    "scientific_payload_and_claim_ceiling_unchanged": True,
}

_FD_LEASES: list[tuple[Path, int, os.stat_result, str | None]] = []
_DIRECTORY_LEASES: list[tuple[Path, int, os.stat_result]] = []


class RecoveryAuthorityError(RuntimeError):
    pass


def _strict_equal(left: Any, right: Any) -> bool:
    if type(left) is not type(right):
        return False
    if isinstance(left, Mapping):
        return set(left) == set(right) and all(
            _strict_equal(left[key], right[key]) for key in left
        )
    if isinstance(left, (list, tuple)):
        return len(left) == len(right) and all(
            _strict_equal(left_value, right_value)
            for left_value, right_value in zip(left, right)
        )
    return left == right


def _json_sha256(value: Any) -> str:
    data = json.dumps(
        value,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")
    return hashlib.sha256(data).hexdigest()


def _validate_ceremony_absence_contract() -> dict[str, Any]:
    paths = sorted(str(path) for path in CEREMONY_FORBIDDEN_OUTPUT_SLOTS)
    expected_parents = {RESULT_ROOT, REVIEW_ROOT, WORKSPACE_ROOT}
    if (
        len(paths) != EXPECTED_CEREMONY_FORBIDDEN_OUTPUT_SLOT_COUNT
        or len(paths) != len(set(paths))
        or {path.parent for path in CEREMONY_FORBIDDEN_OUTPUT_SLOTS}
        != expected_parents
        or any(path in CEREMONY_FORBIDDEN_OUTPUT_SLOTS for path in REVIEW_PATHS.values())
        or hashlib.sha256(
            json.dumps(
                paths,
                ensure_ascii=True,
                separators=(",", ":"),
            ).encode("ascii")
        ).hexdigest()
        != CEREMONY_FORBIDDEN_OUTPUT_SLOTS_SHA256
        or hashlib.sha256(
            json.dumps(
                list(CEREMONY_STAGING_PATTERNS),
                ensure_ascii=True,
                separators=(",", ":"),
            ).encode("ascii")
        ).hexdigest()
        != CEREMONY_STAGING_PATTERNS_SHA256
    ):
        raise RecoveryAuthorityError("Ceremony forbidden-output inventory drift")
    return {
        "forbidden_output_slot_count": len(paths),
        "forbidden_output_slots_sha256": CEREMONY_FORBIDDEN_OUTPUT_SLOTS_SHA256,
        "staging_pattern_count": len(CEREMONY_STAGING_PATTERNS),
        "staging_patterns_sha256": CEREMONY_STAGING_PATTERNS_SHA256,
        "leased_parent_directories": sorted(str(path) for path in expected_parents),
        "pass": True,
    }


def _exact_keys(value: Any, expected: set[str], label: str) -> None:
    if not isinstance(value, Mapping) or set(value) != expected:
        observed = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise RecoveryAuthorityError(
            f"{label} keys are not exact: observed={observed}, expected={sorted(expected)}"
        )


def _open_bound(path: Path, *, expected_mode: str = "0o444") -> tuple[bytes, dict[str, Any], int]:
    if not path.is_absolute() or path.resolve(strict=True) != path:
        raise RecoveryAuthorityError(f"Non-canonical authority path: {path}")
    before = os.stat(path, follow_symlinks=False)
    if not stat.S_ISREG(before.st_mode):
        raise RecoveryAuthorityError(f"Authority input is not a regular file: {path}")
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    opened = os.fstat(descriptor)
    data = b"".join(
        os.pread(descriptor, min(8 * 1024 * 1024, opened.st_size - offset), offset)
        for offset in range(0, opened.st_size, 8 * 1024 * 1024)
    )
    after = os.stat(path, follow_symlinks=False)
    if (
        (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns)
        != (opened.st_dev, opened.st_ino, opened.st_size, opened.st_mtime_ns, opened.st_ctime_ns)
        or (opened.st_dev, opened.st_ino, opened.st_size, opened.st_mtime_ns, opened.st_ctime_ns)
        != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns)
        or len(data) != opened.st_size
        or opened.st_nlink != 1
        or oct(stat.S_IMODE(opened.st_mode)) != expected_mode
    ):
        os.close(descriptor)
        raise RecoveryAuthorityError(f"Authority input identity/mode drift: {path}")
    record = {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": expected_mode,
    }
    _FD_LEASES.append((path, descriptor, opened, record["sha256"]))
    return data, record, descriptor


def _records(paths: Mapping[str, Path]) -> dict[str, dict[str, Any]]:
    return {role: _open_bound(path)[1] for role, path in paths.items()}


def _open_retired_private_key_bound(path: Path) -> tuple[dict[str, str], int]:
    if not path.is_absolute() or path.resolve(strict=True) != path:
        raise RecoveryAuthorityError(f"Non-canonical retired private-key path: {path}")
    before = os.stat(path, follow_symlinks=False)
    if (
        not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or stat.S_IMODE(before.st_mode) != 0
        or not hasattr(os, "O_PATH")
    ):
        raise RecoveryAuthorityError("Recovery private key is not retired mode 0000")
    flags = os.O_PATH | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    opened = os.fstat(descriptor)
    after = os.stat(path, follow_symlinks=False)
    identity = lambda value: (
        value.st_dev,
        value.st_ino,
        value.st_size,
        value.st_mode,
        value.st_nlink,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )
    if identity(before) != identity(opened) or identity(opened) != identity(after):
        os.close(descriptor)
        raise RecoveryAuthorityError("Recovery private-key metadata changed while binding")
    _FD_LEASES.append((path, descriptor, opened, None))
    return {"path": str(path), "mode": "0o0"}, descriptor


def _open_private_key_metadata_bound(
    path: Path, *, expected_mode: int, state: str
) -> tuple[dict[str, Any], int]:
    if not path.is_absolute() or path.resolve(strict=True) != path:
        raise RecoveryAuthorityError(f"Non-canonical private-key path: {path}")
    before = os.stat(path, follow_symlinks=False)
    if (
        not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or stat.S_IMODE(before.st_mode) != expected_mode
        or not hasattr(os, "O_PATH")
    ):
        raise RecoveryAuthorityError(f"Private-key metadata state drift: {path}")
    flags = os.O_PATH | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    opened = os.fstat(descriptor)
    after = os.stat(path, follow_symlinks=False)
    identity = lambda value: (
        value.st_dev,
        value.st_ino,
        value.st_size,
        value.st_mode,
        value.st_nlink,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )
    if identity(before) != identity(opened) or identity(opened) != identity(after):
        os.close(descriptor)
        raise RecoveryAuthorityError("Private-key metadata changed while binding")
    _FD_LEASES.append((path, descriptor, opened, None))
    return {
        "path": str(path),
        "size_bytes": opened.st_size,
        "mode": oct(expected_mode),
        "state": state,
    }, descriptor


def _validate_terminal_key_rotation(
    *, expected_recovery_private_mode: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    if expected_recovery_private_mode not in {"0o400", "0o0"}:
        raise RecoveryAuthorityError("Unknown recovery continuation private-key mode")
    authorized_public_data, authorized_public_record, _ = _open_bound(
        AUTHORIZED_CONTINUATION_PUBLIC_KEY
    )
    authorized_private_record, _ = _open_private_key_metadata_bound(
        AUTHORIZED_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="LEGACY_SIGNED_CONTRACT_UNUSED_NOT_RETIRED",
    )
    failed_v16_public_data, failed_v16_public_record, _ = _open_bound(
        FAILED_V16_CONTINUATION_PUBLIC_KEY
    )
    failed_v16_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V16_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="FAILED_V16_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    failed_v17_public_data, failed_v17_public_record, _ = _open_bound(
        FAILED_V17_CONTINUATION_PUBLIC_KEY
    )
    failed_v17_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V17_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="FAILED_V17_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    failed_v18_public_data, failed_v18_public_record, _ = _open_bound(
        FAILED_V18_CONTINUATION_PUBLIC_KEY
    )
    failed_v18_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V18_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="FAILED_V18_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    failed_v19_public_data, failed_v19_public_record, _ = _open_bound(
        FAILED_V19_CONTINUATION_PUBLIC_KEY
    )
    failed_v19_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V19_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="FAILED_V19_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    rejected_v20_public_data, rejected_v20_public_record, _ = _open_bound(
        REJECTED_V20_CONTINUATION_PUBLIC_KEY
    )
    rejected_v20_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V20_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V20_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    failed_v21_public_data, failed_v21_public_record, _ = _open_bound(
        FAILED_V21_CONTINUATION_PUBLIC_KEY
    )
    failed_v21_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V21_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="FAILED_V21_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    rejected_v22_public_data, rejected_v22_public_record, _ = _open_bound(
        REJECTED_V22_CONTINUATION_PUBLIC_KEY
    )
    rejected_v22_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V22_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V22_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    rejected_v23_public_data, rejected_v23_public_record, _ = _open_bound(
        REJECTED_V23_CONTINUATION_PUBLIC_KEY
    )
    rejected_v23_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V23_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V23_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    rejected_v24_public_data, rejected_v24_public_record, _ = _open_bound(
        REJECTED_V24_CONTINUATION_PUBLIC_KEY
    )
    rejected_v24_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V24_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V24_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    failed_v25_public_data, failed_v25_public_record, _ = _open_bound(
        FAILED_V25_TERMINAL_PUBLIC_KEY
    )
    failed_v25_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V25_TERMINAL_PRIVATE_KEY,
        expected_mode=0o400,
        state="FAILED_V25_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
    )
    failed_v26_public_data, failed_v26_public_record, _ = _open_bound(
        FAILED_V26_TERMINAL_PUBLIC_KEY
    )
    failed_v26_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V26_TERMINAL_PRIVATE_KEY,
        expected_mode=0o400,
        state="FAILED_V26_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
    )
    rejected_v27_public_data, rejected_v27_public_record, _ = _open_bound(
        REJECTED_V27_CONTINUATION_PUBLIC_KEY
    )
    rejected_v27_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V27_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V27_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
    )
    failed_v28_terminal_spec = FAILED_V28_KEY_ARCHIVES["terminal_v18"]
    failed_v28_public_data, failed_v28_public_record, _ = _open_bound(
        failed_v28_terminal_spec["public"]
    )
    failed_v28_private_record, _ = _open_private_key_metadata_bound(
        failed_v28_terminal_spec["private"],
        expected_mode=0o400,
        state="FAILED_V28_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
    )
    failed_v29_terminal_spec = FAILED_V29_KEY_ARCHIVES["terminal_v19"]
    failed_v29_terminal_public_key = (
        failed_v29_terminal_spec["root"] / "ed25519_public.pem"
    )
    failed_v29_terminal_private_key = (
        failed_v29_terminal_spec["root"] / failed_v29_terminal_spec["private_name"]
    )
    failed_v29_public_data, failed_v29_public_record, _ = _open_bound(
        failed_v29_terminal_public_key
    )
    failed_v29_private_record, _ = _open_private_key_metadata_bound(
        failed_v29_terminal_private_key,
        expected_mode=0o400,
        state="FAILED_V29_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
    )
    recovery_public_data, recovery_public_record, _ = _open_bound(
        RECOVERY_CONTINUATION_PUBLIC_KEY
    )
    recovery_private_record, _ = _open_private_key_metadata_bound(
        RECOVERY_CONTINUATION_PRIVATE_KEY,
        expected_mode=int(expected_recovery_private_mode, 8),
        state=(
            "RECOVERY_V30_READY_FOR_SINGLE_SIGNATURE"
            if expected_recovery_private_mode == "0o400"
            else "RECOVERY_V30_RETIRED_AFTER_SINGLE_SIGNATURE"
        ),
    )
    if (
        len(authorized_public_data) != 113
        or authorized_public_record["sha256"]
        != EXPECTED_AUTHORIZED_CONTINUATION_PUBLIC_KEY_SHA256
        or len(failed_v16_public_data) != 113
        or failed_v16_public_record["sha256"]
        != "066949c0c36be413cd2fb60670e5a2fbc583ab9a3a4264ebf4d3766aba39867f"
        or len(failed_v17_public_data) != 113
        or failed_v17_public_record["sha256"]
        != "4b83e655f1a7a778691e27dc3df2257a230001c702afdf6e703e578c706e0b03"
        or len(failed_v18_public_data) != 113
        or failed_v18_public_record["sha256"]
        != "9fa7667e8076cee90a93fee44cfe08bc45a68e9ce28d01507b5c057463102a93"
        or len(failed_v19_public_data) != 113
        or failed_v19_public_record["sha256"]
        != "34d11ff6a699b96aaa8624b30f45fbc8845e6958dc04122b19c413a1a2c12c2d"
        or len(rejected_v20_public_data) != 113
        or rejected_v20_public_record["sha256"]
        != "09794c2ce162af3bf2f3117f6d11dea0c4bd626cbe50946267609058c6c0c291"
        or len(failed_v21_public_data) != 113
        or failed_v21_public_record["sha256"]
        != "825111aa7dcd25e60c6357243e228422d2cf07855d682aae30d90ebdcd5559d2"
        or len(rejected_v22_public_data) != 113
        or rejected_v22_public_record["sha256"]
        != "94ebaec5e5fc994dc75ffcf189f124726154fccc89ff1ea32add8398f06bf5b2"
        or len(rejected_v23_public_data) != 113
        or rejected_v23_public_record["sha256"]
        != "d050c0dfea29b469e271967ea7759d4cbb36c3cbec493010163be5938e81e54a"
        or len(rejected_v24_public_data) != 113
        or rejected_v24_public_record["sha256"]
        != "dcff9fbb753ebdda8525597354be18be7b3fa3d22a4f40a020cd9f677736cc8a"
        or len(failed_v25_public_data) != 113
        or failed_v25_public_record["sha256"]
        != "b0056c3f60d7a7204d782ac1cea31e1e9411200b371295e38b9afc3d5f67a1d1"
        or len(failed_v26_public_data) != 113
        or failed_v26_public_record["sha256"]
        != "b61e7f75ba3e418098c3c7f9fb19da0b261380d7985833eb6a59b99d6e1aeaee"
        or len(rejected_v27_public_data) != 113
        or rejected_v27_public_record["sha256"]
        != "355979601138f8dca29534db58e3862f30f30b49526c01d10a97b01ba91c26f9"
        or len(failed_v28_public_data) != 113
        or failed_v28_public_record["sha256"]
        != "ded65f60581476b08c530c04d795c55bcb393558f393acb0f171629dce47add7"
        or len(failed_v29_public_data) != 113
        or failed_v29_public_record["sha256"]
        != failed_v29_terminal_spec["public_sha256"]
        or len(recovery_public_data) != 113
        or recovery_public_record["sha256"]
        != EXPECTED_RECOVERY_CONTINUATION_PUBLIC_KEY_SHA256
        or len(
            {
                authorized_public_record["sha256"],
                failed_v16_public_record["sha256"],
                failed_v17_public_record["sha256"],
                failed_v18_public_record["sha256"],
                failed_v19_public_record["sha256"],
                rejected_v20_public_record["sha256"],
                failed_v21_public_record["sha256"],
                rejected_v22_public_record["sha256"],
                rejected_v23_public_record["sha256"],
                rejected_v24_public_record["sha256"],
                failed_v25_public_record["sha256"],
                failed_v26_public_record["sha256"],
                rejected_v27_public_record["sha256"],
                failed_v28_public_record["sha256"],
                failed_v29_public_record["sha256"],
                recovery_public_record["sha256"],
            }
        )
        != 16
        or len(
            {
                AUTHORIZED_CONTINUATION_PRIVATE_KEY,
                FAILED_V16_CONTINUATION_PRIVATE_KEY,
                FAILED_V17_CONTINUATION_PRIVATE_KEY,
                FAILED_V18_CONTINUATION_PRIVATE_KEY,
                FAILED_V19_CONTINUATION_PRIVATE_KEY,
                REJECTED_V20_CONTINUATION_PRIVATE_KEY,
                FAILED_V21_CONTINUATION_PRIVATE_KEY,
                REJECTED_V22_CONTINUATION_PRIVATE_KEY,
                REJECTED_V23_CONTINUATION_PRIVATE_KEY,
                REJECTED_V24_CONTINUATION_PRIVATE_KEY,
                FAILED_V25_TERMINAL_PRIVATE_KEY,
                FAILED_V26_TERMINAL_PRIVATE_KEY,
                REJECTED_V27_CONTINUATION_PRIVATE_KEY,
                failed_v28_terminal_spec["private"],
                failed_v29_terminal_private_key,
                RECOVERY_CONTINUATION_PRIVATE_KEY,
            }
        )
        != 16
    ):
        raise RecoveryAuthorityError("Continuation terminal key rotation identity drift")

    contract = {
        "legacy_signed_contract": {
            "public_key": authorized_public_record,
            "private_key": {
                "path": str(AUTHORIZED_CONTINUATION_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_remain_unretired": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.v1.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.v1.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.v1.json.ed25519.sig"
                ),
            },
            "status": "PRESERVED_NOT_EXECUTED",
        },
        "failed_v16_contract": {
            "public_key": failed_v16_public_record,
            "private_key": {
                "path": str(FAILED_V16_CONTINUATION_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v16.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v16.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v16.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v16.json"
                ),
            },
            "status": "FAILED_BEFORE_C_CHILD_UNUSED_KEY_QUARANTINED",
        },
        "failed_v17_contract": {
            "public_key": failed_v17_public_record,
            "private_key": {
                "path": str(FAILED_V17_CONTINUATION_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v17.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v17.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v17.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v17.json"
                ),
            },
            "status": "FAILED_AFTER_FRESH_V_BEFORE_DOWNSTREAM_UNUSED_KEY_QUARANTINED",
        },
        "failed_v18_contract": {
            "public_key": failed_v18_public_record,
            "private_key": {
                "path": str(FAILED_V18_CONTINUATION_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v18.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v18.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v18.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v18.json"
                ),
            },
            "status": "FAILED_AFTER_FRESH_V_AND_R_BEFORE_DOWNSTREAM_UNUSED_KEY_QUARANTINED",
        },
        "failed_v19_contract": {
            "public_key": failed_v19_public_record,
            "private_key": {
                "path": str(FAILED_V19_CONTINUATION_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v19.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v19.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v19.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v19.json"
                ),
            },
            "status": (
                "FAILED_AFTER_FRESH_V_AND_R_BEFORE_DOWNSTREAM_MODE000_INOTIFY_"
                "EACCES_UNUSED_KEY_QUARANTINED"
            ),
        },
        "rejected_v20_contract": {
            "public_key": rejected_v20_public_record,
            "private_key": {
                "path": str(REJECTED_V20_CONTINUATION_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v20.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v20.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v20.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v20.json"
                ),
            },
            "status": "REJECTED_PRE_AUTHORITY_UNUSED_KEY_QUARANTINED",
        },
        "failed_v21_contract": {
            "public_key": failed_v21_public_record,
            "private_key": {
                "path": str(FAILED_V21_CONTINUATION_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v21.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v21.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v21.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v21.json"
                ),
                "incident": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_incident.recovery.v21.json"
                ),
            },
            "status": (
                "FAILED_AFTER_INTERMEDIATE_DOWNSTREAM_OUTPUTS_BEFORE_FINAL_DATASET_"
                "UNUSED_KEY_QUARANTINED"
            ),
        },
        "rejected_v22_contract": {
            "public_key": rejected_v22_public_record,
            "private_key": {
                "path": str(REJECTED_V22_CONTINUATION_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v22.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v22.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v22.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v22.json"
                ),
            },
            "status": "REJECTED_PRE_AUTHORITY_UNUSED_KEY_QUARANTINED",
        },
        "rejected_v23_contract": {
            "public_key": rejected_v23_public_record,
            "private_key": {
                "path": str(REJECTED_V23_CONTINUATION_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v23.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v23.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v23.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v23.json"
                ),
            },
            "status": "REJECTED_PRE_AUTHORITY_UNUSED_KEY_QUARANTINED",
        },
        "rejected_v24_contract": {
            "public_key": rejected_v24_public_record,
            "private_key": {
                "path": str(REJECTED_V24_CONTINUATION_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v24.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v24.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v24.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v24.json"
                ),
            },
            "status": "REJECTED_PRE_AUTHORITY_UNUSED_KEY_QUARANTINED",
        },
        "failed_v25_contract": {
            "public_key": failed_v25_public_record,
            "private_key": {
                "path": str(FAILED_V25_TERMINAL_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v25.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v25.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v25.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v25.json"
                ),
            },
            "status": "SIGNED_AUTHORITY_C_FAILED_UNUSED_TERMINAL_KEY_ARCHIVED",
        },
        "failed_v26_contract": {
            "public_key": failed_v26_public_record,
            "private_key": {
                "path": str(FAILED_V26_TERMINAL_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v26.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v26.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v26.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v26.json"
                ),
            },
            "status": "SIGNED_AUTHORITY_C_FAILED_UNUSED_TERMINAL_KEY_ARCHIVED",
        },
        "rejected_v27_contract": {
            "public_key": rejected_v27_public_record,
            "private_key": {
                "path": str(REJECTED_V27_CONTINUATION_PRIVATE_KEY),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v27.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v27.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v27.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v27.json"
                ),
            },
            "status": "REJECTED_PRE_AUTHORITY_UNUSED_TERMINAL_KEY_ARCHIVED",
        },
        "failed_v28_contract": {
            "public_key": failed_v28_public_record,
            "private_key": {
                "path": str(failed_v28_terminal_spec["private"]),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v28.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v28.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v28.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v28.json"
                ),
            },
            "status": "SIGNED_AUTHORITY_C_REPORT_FAILED_UNUSED_TERMINAL_KEY_ARCHIVED",
        },
        "failed_v29_contract": {
            "public_key": failed_v29_public_record,
            "private_key": {
                "path": str(failed_v29_terminal_private_key),
                "required_mode": "0o400",
                "must_not_be_reused": True,
            },
            "terminal_outputs": {
                "execution_receipt": str(
                    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v29.json"
                ),
                "signed_artifact": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v29.json"
                ),
                "signature": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_exit_attestation.recovery.v29.json.ed25519.sig"
                ),
                "success_witness": str(
                    RESULT_ROOT
                    / "m2v5_downstream_continuation_supervisor_success.recovery.v29.json"
                ),
            },
            "status": "SIGNED_AUTHORITY_R_FAILED_UNUSED_TERMINAL_KEY_ARCHIVED",
        },
        "recovery_v30_contract": {
            "public_key": recovery_public_record,
            "private_key": {
                "path": str(RECOVERY_CONTINUATION_PRIVATE_KEY),
                "required_pre_signature_mode": "0o400",
                "required_post_signature_mode": "0o0",
            },
            "terminal_outputs": {
                "execution_receipt": str(RECOVERY_CONTINUATION_RECEIPT),
                "signed_artifact": str(RECOVERY_EXIT_ATTESTATION),
                "signature": str(RECOVERY_EXIT_SIGNATURE),
                "success_witness": str(RECOVERY_SUCCESS_WITNESS),
                "incident": str(RECOVERY_INCIDENT),
            },
            "status": "AUTHORIZED_BY_SIGNED_V30_RECOVERY_AUTHORITY_ONLY",
        },
        "keys_are_pairwise_distinct": True,
        "legacy_signed_payload_changed": False,
    }
    current_state = {
        "legacy_public_key": authorized_public_record,
        "legacy_private_key": authorized_private_record,
        "failed_v16_public_key": failed_v16_public_record,
        "failed_v16_private_key": failed_v16_private_record,
        "failed_v17_public_key": failed_v17_public_record,
        "failed_v17_private_key": failed_v17_private_record,
        "failed_v18_public_key": failed_v18_public_record,
        "failed_v18_private_key": failed_v18_private_record,
        "failed_v19_public_key": failed_v19_public_record,
        "failed_v19_private_key": failed_v19_private_record,
        "rejected_v20_public_key": rejected_v20_public_record,
        "rejected_v20_private_key": rejected_v20_private_record,
        "failed_v21_public_key": failed_v21_public_record,
        "failed_v21_private_key": failed_v21_private_record,
        "rejected_v22_public_key": rejected_v22_public_record,
        "rejected_v22_private_key": rejected_v22_private_record,
        "rejected_v23_public_key": rejected_v23_public_record,
        "rejected_v23_private_key": rejected_v23_private_record,
        "rejected_v24_public_key": rejected_v24_public_record,
        "rejected_v24_private_key": rejected_v24_private_record,
        "failed_v25_public_key": failed_v25_public_record,
        "failed_v25_private_key": failed_v25_private_record,
        "failed_v26_public_key": failed_v26_public_record,
        "failed_v26_private_key": failed_v26_private_record,
        "rejected_v27_public_key": rejected_v27_public_record,
        "rejected_v27_private_key": rejected_v27_private_record,
        "failed_v28_public_key": failed_v28_public_record,
        "failed_v28_private_key": failed_v28_private_record,
        "failed_v29_public_key": failed_v29_public_record,
        "failed_v29_private_key": failed_v29_private_record,
        "recovery_public_key": recovery_public_record,
        "recovery_private_key": recovery_private_record,
    }
    return contract, current_state


def _require_leases(*, required_paths: set[Path] | None = None) -> None:
    leased_paths: set[Path] = set()
    for path, descriptor, expected, expected_sha256 in _FD_LEASES:
        leased_paths.add(path)
        opened = os.fstat(descriptor)
        live = os.stat(path, follow_symlinks=False)
        identity = lambda value: (
            value.st_dev,
            value.st_ino,
            value.st_size,
            value.st_mode,
            value.st_nlink,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
        if identity(opened) != identity(expected) or identity(live) != identity(expected):
            raise RecoveryAuthorityError(f"Authority descriptor lease drift: {path}")
        if expected_sha256 is not None:
            data = b"".join(
                os.pread(descriptor, min(8 * 1024 * 1024, opened.st_size - offset), offset)
                for offset in range(0, opened.st_size, 8 * 1024 * 1024)
            )
            if len(data) != opened.st_size or hashlib.sha256(data).hexdigest() != expected_sha256:
                raise RecoveryAuthorityError(f"Authority descriptor content drift: {path}")
    for path, descriptor, expected in _DIRECTORY_LEASES:
        leased_paths.add(path)
        opened = os.fstat(descriptor)
        live = os.stat(path, follow_symlinks=False)
        identity = lambda value: (
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_nlink,
        )
        if identity(opened) != identity(expected) or identity(live) != identity(expected):
            raise RecoveryAuthorityError(f"Authority directory lease drift: {path}")
    if required_paths is not None and not required_paths.issubset(leased_paths):
        missing = sorted(str(path) for path in required_paths - leased_paths)
        raise RecoveryAuthorityError(f"Authority descriptor lease coverage gap: {missing}")


def _open_bundle_directory() -> int:
    path = AUTHORITY_BUNDLE
    if not path.is_absolute() or path.resolve(strict=True) != path:
        raise RecoveryAuthorityError("Recovery authority bundle path is not canonical")
    before = os.stat(path, follow_symlinks=False)
    if (
        not stat.S_ISDIR(before.st_mode)
        or stat.S_IMODE(before.st_mode) != 0o700
        or before.st_nlink != 2
    ):
        raise RecoveryAuthorityError("Recovery authority bundle metadata drift")
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    opened = os.fstat(descriptor)
    after = os.stat(path, follow_symlinks=False)
    identity = lambda value: (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
    )
    if identity(before) != identity(opened) or identity(opened) != identity(after):
        os.close(descriptor)
        raise RecoveryAuthorityError("Recovery authority bundle changed while binding")
    expected_members = {"authority.json", "authority.ed25519.sig", "commit.json"}
    if set(os.listdir(descriptor)) != expected_members:
        os.close(descriptor)
        raise RecoveryAuthorityError("Recovery authority bundle member set drift")
    _DIRECTORY_LEASES.append((path, descriptor, opened))
    return descriptor


def _load_json(data: bytes, label: str) -> Mapping[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        output: dict[str, Any] = {}
        for key, value in pairs:
            if key in output:
                raise RecoveryAuthorityError(f"Duplicate JSON key in {label}: {key}")
            output[key] = value
        return output

    def reject_nonfinite(value: str) -> None:
        raise RecoveryAuthorityError(f"Non-finite JSON value in {label}: {value}")

    def reject_float(value: str) -> None:
        raise RecoveryAuthorityError(f"Floating JSON value in {label}: {value}")

    try:
        value = json.loads(
            data,
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_nonfinite,
            parse_float=reject_float,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RecoveryAuthorityError(f"Invalid JSON for {label}") from error
    if not isinstance(value, Mapping):
        raise RecoveryAuthorityError(f"{label} is not a JSON object")
    return value


def _load_archived_transport_json(data: bytes, label: str) -> Mapping[str, Any]:
    """Load hash-pinned CLI transport metadata without treating costs as evidence."""

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        output: dict[str, Any] = {}
        for key, value in pairs:
            if key in output:
                raise RecoveryAuthorityError(f"Duplicate JSON key in {label}: {key}")
            output[key] = value
        return output

    def reject_nonfinite(value: str) -> None:
        raise RecoveryAuthorityError(f"Non-finite JSON value in {label}: {value}")

    try:
        value = json.loads(
            data,
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_nonfinite,
            parse_float=lambda value: value,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RecoveryAuthorityError(f"Invalid JSON for {label}") from error
    if not isinstance(value, Mapping):
        raise RecoveryAuthorityError(f"{label} is not a JSON object")
    return value


def _bundle_member(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "name": Path(record["path"]).name,
        "size_bytes": record["size_bytes"],
        "sha256": record["sha256"],
        "mode": record["mode"],
    }


def _validate_bundle_commit(
    data: bytes,
    authority_record: Mapping[str, Any],
    signature_record: Mapping[str, Any],
    retired_private_key: Mapping[str, Any],
) -> Mapping[str, Any]:
    payload = _load_json(data, "recovery authority atomic commit")
    _exact_keys(
        payload,
        {
            "schema_name",
            "schema_version",
            "transaction_id",
            "members",
            "authority",
            "signature",
            "retired_private_key",
            "pass",
        },
        "recovery authority atomic commit",
    )
    transaction_id = payload.get("transaction_id")
    if (
        payload.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_atomic_commit"
        or payload.get("schema_version") != SCHEMA_VERSION
        or not isinstance(transaction_id, str)
        or len(transaction_id) != 32
        or any(character not in "0123456789abcdef" for character in transaction_id)
        or not _strict_equal(
            payload.get("members"),
            ["authority.ed25519.sig", "authority.json", "commit.json"],
        )
        or not _strict_equal(payload.get("authority"), _bundle_member(authority_record))
        or not _strict_equal(payload.get("signature"), _bundle_member(signature_record))
        or not _strict_equal(payload.get("retired_private_key"), retired_private_key)
        or payload.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Recovery authority atomic commit contract drift")
    return payload


def _verify_signature(data_fd: int, public_fd: int, signature_fd: int) -> None:
    command = [
        str(OPENSSL),
        "pkeyutl",
        "-verify",
        "-pubin",
        "-inkey",
        f"/proc/self/fd/{public_fd}",
        "-rawin",
        "-in",
        f"/proc/self/fd/{data_fd}",
        "-sigfile",
        f"/proc/self/fd/{signature_fd}",
    ]
    completed = subprocess.run(
        command,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        close_fds=True,
        pass_fds=(data_fd, public_fd, signature_fd),
        env={"PATH": "/usr/bin:/bin", "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8"},
    )
    if completed.returncode != 0:
        raise RecoveryAuthorityError("Recovery authority Ed25519 signature verification failed")


def _validate_pinned_records(
    legacy_sources: Mapping[str, Mapping[str, Any]],
    original_chain: Mapping[str, Mapping[str, Any]],
) -> None:
    for role, expected_sha256 in EXPECTED_LEGACY_SHA256.items():
        if legacy_sources[role]["sha256"] != expected_sha256:
            raise RecoveryAuthorityError(f"Legacy source SHA drift: {role}")
    for role, expected_sha256 in EXPECTED_ORIGINAL_CHAIN_SHA256.items():
        if original_chain[role]["sha256"] != expected_sha256:
            raise RecoveryAuthorityError(f"Original promotion-chain SHA drift: {role}")
    source = original_chain["source_receipt"]
    canonical = original_chain["canonical_receipt"]
    if (
        source["size_bytes"] != 6_733
        or canonical["size_bytes"] != 6_733
        or source["sha256"] != canonical["sha256"]
    ):
        raise RecoveryAuthorityError("Source/canonical tumor-REF receipt bytes drift")


def _validate_prior_recovery_chain(
    records: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    if set(records) != set(PRIOR_RECOVERY_CHAIN_PATHS):
        raise RecoveryAuthorityError("Prior recovery-chain artifact set drift")
    for role, expected_sha256 in EXPECTED_PRIOR_RECOVERY_CHAIN_SHA256.items():
        if records[role]["sha256"] != expected_sha256:
            raise RecoveryAuthorityError(f"Prior recovery-chain SHA drift: {role}")

    authority_data, authority_record, authority_fd = _open_bound(
        PRIOR_RECOVERY_CHAIN_PATHS["authority"]
    )
    signature_data, signature_record, signature_fd = _open_bound(
        PRIOR_RECOVERY_CHAIN_PATHS["authority_signature"]
    )
    public_data, public_record, public_fd = _open_bound(PRIOR_RECOVERY_PUBLIC_KEY)
    if (
        not _strict_equal(authority_record, records["authority"])
        or not _strict_equal(signature_record, records["authority_signature"])
        or not _strict_equal(public_record, records["public_key"])
        or len(signature_data) != 64
        or len(public_data) != 113
    ):
        raise RecoveryAuthorityError("Prior recovery authority artifact identity drift")
    _verify_signature(authority_fd, public_fd, signature_fd)

    authority_payload = _load_json(authority_data, "prior recovery authority")
    if (
        authority_payload.get("authority_id")
        != "20260722_tumor_ref_schema_recovery_v1"
        or authority_payload.get("pass") is not True
        or authority_payload.get("recovery_sources", {})
        .get("runner_gate_replay", {})
        .get("sha256")
        != EXPECTED_PRIOR_RECOVERY_CHAIN_SHA256["runner_gate_replay"]
    ):
        raise RecoveryAuthorityError("Prior recovery authority payload drift")

    failure_data = _open_bound(PRIOR_RECOVERY_CHAIN_PATHS["runner_failure_evidence"])[0]
    failure_payload = _load_json(failure_data, "prior runner failure evidence")
    output_state = failure_payload.get("output_state_after_failure")
    if (
        failure_payload.get("status") != "FAIL_CLOSED_ZERO_REPLAY_OUTPUT"
        or failure_payload.get("execution", {}).get("exit_code") != 1
        or type(failure_payload.get("execution", {}).get("exit_code")) is not int
        or not _strict_equal(
            output_state,
            {
                "downstream_outputs_created": False,
                "replay_log_created": False,
                "replay_receipt_created": False,
            },
        )
        or failure_payload.get("root_cause", {}).get("defect_class")
        != "historical_transition_record_rebound_as_live_identity"
    ):
        raise RecoveryAuthorityError("Prior runner failure evidence contract drift")

    retired_private, _ = _open_retired_private_key_bound(PRIOR_RECOVERY_PRIVATE_KEY)
    return {
        "artifacts": {role: dict(record) for role, record in records.items()},
        "authority_signature_verified": True,
        "retired_private_key": retired_private,
        "runner_output_created": False,
    }


def _validate_rejected_generation(
    *,
    generation: str,
    evidence_path: Path,
    expected_size: int,
    expected_sha256: str,
    expected_source_set_sha256: str,
    private_key: Path,
    output_slots: tuple[Path, ...],
    decisive_reviewer: str,
    expected_high: list[str],
    expected_medium: list[str],
    expected_review_artifacts_created: bool = False,
    expected_archived_reviews: Mapping[str, Path] | None = None,
    expected_post_review_probe_failure: bool = False,
) -> dict[str, Any]:
    data, record, _ = _open_bound(evidence_path)
    if record["size_bytes"] != expected_size or record["sha256"] != expected_sha256:
        raise RecoveryAuthorityError(f"Rejected {generation} evidence identity drift")
    payload = _load_json(data, f"rejected {generation} recovery generation")
    review = payload.get("formal_reviews", {}).get(decisive_reviewer, {})
    if (
        payload.get("status") != "REJECTED_BEFORE_SIGNING"
        or payload.get("authority_output_created") is not False
        or payload.get("review_artifact_files_created")
        is not expected_review_artifacts_created
        or payload.get(f"{generation}_runtime_outputs_created") is not False
        or payload.get("source_set_sha256") != expected_source_set_sha256
        or review.get("verdict") != "REJECT"
        or review.get("high_findings") != expected_high
        or review.get("medium_findings") != expected_medium
    ):
        raise RecoveryAuthorityError(f"Rejected {generation} evidence contract drift")
    archived_review_records: dict[str, dict[str, Any]] = {}
    if expected_archived_reviews is not None:
        archived_review_records = _records(expected_archived_reviews)
        expected_post_review = {
            "direct_regression_exit_code": 1,
            "failed_test": "test_v6_full_builder_reuses_bootstrap_validator_record",
            "failure_class": "state_dependent_missing_review_expectation",
            "forbidden_output_slots_created": 0,
            "probe_exit_code": 1,
            "regression_summary": "45 passed, 1 failed",
            "review_evidence_state": "all_present",
        }
        if (
            payload.get("review_artifacts_archived") is not True
            or payload.get("active_review_slots_absent") is not True
            or payload.get("signing_attempted") is not False
            or not _strict_equal(
                payload.get("archived_review_artifacts"), archived_review_records
            )
            or not _strict_equal(
                payload.get("probe_evidence", {}).get("post_review"),
                expected_post_review,
            )
            or not expected_post_review_probe_failure
        ):
            raise RecoveryAuthorityError(
                f"Rejected {generation} archived-review contract drift"
            )
    elif expected_post_review_probe_failure:
        raise RecoveryAuthorityError(
            f"Rejected {generation} post-review failure lacks archived reviews"
        )
    occupied = [str(path) for path in output_slots if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            f".tumor_ref_promotion_schema_recovery_authority.{generation}.bundle.staging.*"
        )
    )
    if occupied or staging:
        raise RecoveryAuthorityError(
            f"Rejected {generation} output slot became occupied: {occupied + staging}"
        )
    retired_private, _ = _open_retired_private_key_bound(private_key)
    return {
        "evidence": record,
        "private_key": retired_private,
        "authority_created": False,
        "review_artifacts_created": expected_review_artifacts_created,
        "archived_review_artifacts": archived_review_records,
        "runtime_outputs_created": False,
        "pass": True,
    }


def _validate_rejected_v13_pre_authority() -> dict[str, Any]:
    _, superseded_record, _ = _open_bound(REJECTED_V13_SUPERSEDED_EVIDENCE)
    evidence_data, evidence_record, _ = _open_bound(REJECTED_V13_EVIDENCE)
    ledger_data, ledger_record, _ = _open_bound(REJECTED_V13_KEY_ARCHIVE_LEDGER)
    if (
        superseded_record["size_bytes"] != 5_830
        or superseded_record["sha256"]
        != EXPECTED_REJECTED_V13_SUPERSEDED_EVIDENCE_SHA256
        or evidence_record["size_bytes"] != 6_694
        or evidence_record["sha256"] != EXPECTED_REJECTED_V13_EVIDENCE_SHA256
        or ledger_record["size_bytes"] != 2_783
        or ledger_record["sha256"]
        != EXPECTED_REJECTED_V13_KEY_ARCHIVE_LEDGER_SHA256
    ):
        raise RecoveryAuthorityError("Rejected v13 evidence/ledger identity drift")
    evidence = _load_json(evidence_data, "rejected v13 pre-authority evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "task_type",
            "superseded_evidence",
            "correction",
            "candidate",
            "review_chronology",
            "adjudication",
            "formal_state",
            "key_archive",
            "replacement_generation",
            "status",
            "pass",
        },
        "rejected v13 pre-authority evidence",
    )
    source_records = _records(REJECTED_V13_SOURCE_PATHS)
    candidate = evidence.get("candidate")
    formal_state = evidence.get("formal_state")
    adjudication = evidence.get("adjudication")
    replacement = evidence.get("replacement_generation")
    key_archive = evidence.get("key_archive")
    reviews = evidence.get("review_chronology")
    correction = evidence.get("correction")
    if (
        evidence.get("schema_name")
        != "intersubmod.schema_recovery_pre_authority_rejection"
        or evidence.get("task_type") != "B_comprehensive_validation"
        or not isinstance(evidence.get("created_at_utc"), str)
        or not evidence["created_at_utc"].endswith("+00:00")
        or not _strict_equal(evidence.get("superseded_evidence"), superseded_record)
        or not _strict_equal(
            correction,
            {
                "field": "candidate.reviewed_source_set_sha256",
                "old_value": (
                    "7bf4358958d52b286b268c4b26158eff6d7c4ab514974157e677a196ca71882e"
                ),
                "new_value": _json_sha256(source_records),
                "reason": (
                    "replace_transcription_value_with_canonical_sorted_json_sha256_"
                    "of_the_seven_exact_source_records"
                ),
                "source_records_changed": False,
                "scientific_payload_changed": False,
                "claim_ceiling_changed": False,
            },
        )
        or evidence.get("schema_version") != "1.1.0"
        or not isinstance(candidate, Mapping)
        or set(candidate) != {"authority_id", "reviewed_source_set_sha256", "sources"}
        or candidate.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v13"
        or candidate.get("reviewed_source_set_sha256") != _json_sha256(source_records)
        or not _strict_equal(candidate.get("sources"), source_records)
        or not isinstance(reviews, list)
        or len(reviews) != 3
        or [review.get("reviewer") for review in reviews]
        != ["Mendel", "Nash", "External Claude Opus"]
        or reviews[0].get("adjudicated_verdict") != "REQUEST_CHANGES"
        or reviews[1].get("verdict") != "REQUEST_CHANGES"
        or reviews[2].get("verdict") != "APPROVE_WITH_LOW_FINDINGS"
        or not isinstance(adjudication, Mapping)
        or adjudication.get("policy") != "strictest_independent_review_wins"
        or type(adjudication.get("blocking_finding_count")) is not int
        or adjudication["blocking_finding_count"] != 2
        or adjudication.get("verdict") != "REJECTED_PRE_AUTHORITY"
        or adjudication.get("pass") is not False
        or not isinstance(formal_state, Mapping)
        or formal_state.get("authority_bundle_created") is not False
        or formal_state.get("formal_reviews_published") is not False
        or formal_state.get("authority_signature_created") is not False
        or formal_state.get("downstream_outputs_created") is not False
        or type(formal_state.get("formal_v13_slot_count")) is not int
        or formal_state["formal_v13_slot_count"] != len(REJECTED_V13_OUTPUT_SLOTS)
        or formal_state.get("formal_v13_slots_all_absent") is not True
        or not isinstance(key_archive, Mapping)
        or not _strict_equal(key_archive.get("ledger"), ledger_record)
        or not isinstance(replacement, Mapping)
        or replacement.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v14"
        or replacement.get("authority_public_key_sha256")
        != "91f3b81a0dfab1911b492269dc40ef150eed76bf3b16c4143ba541d16ffdc8a3"
        or replacement.get("terminal_public_key_sha256")
        != "091785646a7cadf2295f97668838b8279a2bb3b8a798e317dcf1ec6aba33d427"
        or evidence.get("status") != "REJECTED_PRE_AUTHORITY_ARCHIVED"
        or evidence.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v13 evidence contract drift")

    archive_contracts = (
        (
            "schema_recovery_authority_v13",
            REJECTED_V13_AUTHORITY_KEY_ARCHIVE,
            "ed25519_public.pem",
            "ed25519_private_one_time.pem",
            "99e8836c17adfa0559ad92b8605f2c8b1fd7e0a2064de207e34862b4f48ee56c",
        ),
        (
            "downstream_terminal_v3",
            REJECTED_V13_TERMINAL_KEY_ARCHIVE,
            "ed25519_public.pem",
            "ed25519_private_one_time_resident.pem",
            "a1bbafabe577ae3a05dffc0e61566d51ad7436baf980590f9fdac51179ca0d94",
        ),
    )
    archived_keys: dict[str, Any] = {}
    for role, root, public_name, private_name, public_sha256 in archive_contracts:
        if (
            root.resolve(strict=True) != root
            or stat.S_IMODE(os.lstat(root).st_mode) != 0o700
            or sorted(path.name for path in root.iterdir())
            != sorted(
                [public_name, private_name, "UNUSED_KEY_ARCHIVE_RECORD.v1.json"]
            )
        ):
            raise RecoveryAuthorityError(f"Rejected v13 key archive drift: {role}")
        _, public_record, public_fd = _open_bound(root / public_name)
        _, private_record, private_fd = _open_bound(
            root / private_name, expected_mode="0o400"
        )
        status_data, status_record, _ = _open_bound(
            root / "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
        )
        status_payload = _load_json(status_data, f"rejected v13 key archive {role}")
        if (
            public_record["sha256"] != public_sha256
            or public_record["mode"] != "0o444"
            or private_record["mode"] != "0o400"
            or os.fstat(public_fd).st_nlink != 1
            or os.fstat(private_fd).st_nlink != 1
            or status_record["mode"] != "0o444"
            or status_payload.get("role") != role
            or status_payload.get("archive") != str(root)
            or status_payload.get("signature_created") is not False
            or status_payload.get("must_never_be_used") is not True
        ):
            raise RecoveryAuthorityError(
                f"Rejected v13 archived key identity/lifecycle drift: {role}"
            )
        archived_keys[role] = {
            "archive": str(root),
            "public_key": public_record,
            "unused_private_key": private_record,
            "archive_record": status_record,
            "signature_created": False,
            "pass": True,
        }

    occupied = [
        str(path) for path in REJECTED_V13_OUTPUT_SLOTS if os.path.lexists(path)
    ]
    if (
        os.path.lexists(REJECTED_V13_ORIGINAL_AUTHORITY_KEY_ROOT)
        or os.path.lexists(REJECTED_V13_ORIGINAL_TERMINAL_KEY_ROOT)
        or occupied
    ):
        raise RecoveryAuthorityError(
            f"Rejected v13 original key/output state drift: {occupied}"
        )
    ledger_events = [
        _load_json(line, f"rejected v13 archive ledger event {index}")
        for index, line in enumerate(ledger_data.splitlines(), start=1)
    ]
    if (
        len(ledger_events) != 6
        or [event.get("event") for event in ledger_events]
        != [
            "PREFLIGHT_COMPLETE",
            "ARCHIVE_MOVE_START",
            "ARCHIVE_MOVE_COMPLETE",
            "ARCHIVE_MOVE_START",
            "ARCHIVE_MOVE_COMPLETE",
            "ARCHIVE_BATCH_COMPLETE",
        ]
        or ledger_events[1].get("role") != "schema_recovery_authority_v13"
        or ledger_events[2].get("role") != "schema_recovery_authority_v13"
        or ledger_events[3].get("role") != "downstream_terminal_v3"
        or ledger_events[4].get("role") != "downstream_terminal_v3"
        or not _strict_equal(
            ledger_events[-1],
            {
                "entry_count": 2,
                "event": "ARCHIVE_BATCH_COMPLETE",
                "status": "PASS",
            },
        )
    ):
        raise RecoveryAuthorityError("Rejected v13 archive ledger event drift")
    return {
        "superseded_evidence": superseded_record,
        "evidence": evidence_record,
        "key_archive_ledger": ledger_record,
        "sources": source_records,
        "archived_unused_keys": archived_keys,
        "authority_created": False,
        "formal_reviews_published": False,
        "signatures_created": False,
        "outputs_remain_absent": True,
        "blocking_finding_count": 2,
        "pass": True,
    }


def _validate_rejected_v15_pre_authority() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(REJECTED_V15_EVIDENCE)
    ledger_data, ledger_record, _ = _open_bound(REJECTED_V15_KEY_ARCHIVE_LEDGER)
    review_records = _records(REJECTED_V15_REVIEW_PATHS)
    envelope_data, envelope_record, _ = _open_bound(REJECTED_V15_EXTERNAL_ENVELOPE)
    _, stderr_record, _ = _open_bound(REJECTED_V15_EXTERNAL_STDERR)
    if (
        evidence_record["size_bytes"] != 6_290
        or evidence_record["sha256"] != EXPECTED_REJECTED_V15_EVIDENCE_SHA256
        or ledger_record["size_bytes"] != 2_801
        or ledger_record["sha256"]
        != EXPECTED_REJECTED_V15_KEY_ARCHIVE_LEDGER_SHA256
        or envelope_record["size_bytes"] != 18_191
        or envelope_record["sha256"]
        != "86671481a742a61b899eb2fe59a278d27e754de9c261d9084349fd2a320bd332"
        or stderr_record["size_bytes"] != 302
        or stderr_record["sha256"]
        != "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a"
    ):
        raise RecoveryAuthorityError("Rejected v15 evidence identity drift")

    evidence = _load_json(evidence_data, "rejected v15 pre-authority evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "task_type",
            "candidate",
            "review_chronology",
            "adjudication",
            "formal_state",
            "key_archive",
            "replacement_generation",
            "scientific_payload_changed",
            "claim_ceiling_changed",
            "status",
            "pass",
        },
        "rejected v15 pre-authority evidence",
    )
    source_records = _records(REJECTED_V15_SOURCE_PATHS)
    candidate = evidence.get("candidate")
    chronology = evidence.get("review_chronology")
    adjudication = evidence.get("adjudication")
    formal_state = evidence.get("formal_state")
    key_archive = evidence.get("key_archive")
    replacement = evidence.get("replacement_generation")
    if (
        evidence.get("schema_name")
        != "intersubmod.schema_recovery_pre_authority_rejection"
        or evidence.get("schema_version") != "1.0.0"
        or evidence.get("task_type") != "B_comprehensive_validation"
        or not isinstance(evidence.get("created_at_utc"), str)
        or not evidence["created_at_utc"].endswith("+00:00")
        or not isinstance(candidate, Mapping)
        or set(candidate) != {"authority_id", "reviewed_source_set_sha256", "sources"}
        or candidate.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v15"
        or candidate.get("reviewed_source_set_sha256") != _json_sha256(source_records)
        or not _strict_equal(candidate.get("sources"), source_records)
        or not isinstance(chronology, list)
        or len(chronology) != 3
        or [item.get("reviewer") for item in chronology]
        != ["Mendel", "Nash", "External Claude Opus"]
        or [item.get("verdict") for item in chronology]
        != ["REQUEST_CHANGES", "REQUEST_CHANGES", "APPROVE"]
        or [item.get("blocking_finding_count") for item in chronology] != [4, 3, 0]
        or chronology[2].get("review_completed_at_utc_not_used_for_authority")
        is not True
        or not isinstance(adjudication, Mapping)
        or adjudication.get("policy") != "strictest_independent_review_wins"
        or adjudication.get("verdict") != "REJECTED_PRE_AUTHORITY"
        or adjudication.get("pass") is not False
        or not isinstance(formal_state, Mapping)
        or formal_state.get("authority_bundle_created") is not False
        or formal_state.get("formal_reviews_published") is not False
        or formal_state.get("authority_signature_created") is not False
        or formal_state.get("downstream_outputs_created") is not False
        or formal_state.get("formal_v15_slot_count") != len(REJECTED_V15_OUTPUT_SLOTS)
        or type(formal_state.get("formal_v15_slot_count")) is not int
        or formal_state.get("formal_v15_slots_all_absent") is not True
        or not isinstance(key_archive, Mapping)
        or not _strict_equal(key_archive.get("ledger"), ledger_record)
        or not isinstance(replacement, Mapping)
        or replacement.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v16"
        or replacement.get("authority_public_key_sha256")
        != "540b64ed3615618efed89637069f772787fdd025acadbbd27b6e334423d2345e"
        or replacement.get("terminal_public_key_sha256")
        != "066949c0c36be413cd2fb60670e5a2fbc583ab9a3a4264ebf4d3766aba39867f"
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("status") != "REJECTED_PRE_AUTHORITY_ARCHIVED"
        or evidence.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v15 evidence contract drift")

    review_payloads = {
        role: _load_json(
            _open_bound(path)[0], f"rejected v15 candidate review {role}"
        )
        for role, path in REJECTED_V15_REVIEW_PATHS.items()
    }
    expected_reviews = {
        "mendel": (
            "Mendel",
            "019f8dd9-c4f7-7831-aa31-0ba6ce92d6b4",
            "REQUEST_CHANGES",
            False,
            2,
            2,
        ),
        "nash": (
            "Nash",
            "019f8dd9-eb02-7e31-905c-14a7fb0349f9",
            "REQUEST_CHANGES",
            False,
            1,
            2,
        ),
        "external_claude_opus": (
            "External Claude Opus",
            "8540c6b8-509c-404c-9575-07662444c433",
            "APPROVE",
            True,
            0,
            0,
        ),
    }
    for index, role in enumerate(("mendel", "nash", "external_claude_opus")):
        payload = review_payloads[role]
        reviewer, reviewer_id, verdict, passed, high_count, medium_count = (
            expected_reviews[role]
        )
        probe = payload.get("readonly_probe")
        if (
            payload.get("schema_name")
            != "intersubmod.schema_recovery_candidate_review"
            or payload.get("candidate_authority_id")
            != "20260723_tumor_ref_schema_recovery_v15"
            or payload.get("reviewer") != reviewer
            or payload.get("reviewer_agent_id") != reviewer_id
            or payload.get("verdict") != verdict
            or payload.get("pass") is not passed
            or not isinstance(payload.get("high_findings"), list)
            or len(payload["high_findings"]) != high_count
            or not isinstance(payload.get("medium_findings"), list)
            or len(payload["medium_findings"]) != medium_count
            or not isinstance(probe, Mapping)
            or probe.get("exit_code") != 0
            or probe.get("regression_summary") != "350 passed"
            or probe.get("forbidden_output_slots_checked") != 206
            or probe.get("staging_patterns_checked") != 14
            or probe.get("no_output_writes") is not True
            or chronology[index].get("reviewer_id") != reviewer_id
            or chronology[index].get("evidence_sha256")
            != review_records[role]["sha256"]
        ):
            raise RecoveryAuthorityError(
                f"Rejected v15 candidate review contract drift: {role}"
            )

    external_review = review_payloads["external_claude_opus"]
    # The hash-pinned CLI transport envelope contains duration floats; it is
    # historical transport evidence, not an authority or formal-review payload.
    envelope = json.loads(envelope_data.decode("utf-8"))
    structured = envelope.get("structured_output")
    if (
        not _strict_equal(external_review.get("raw_envelope"), envelope_record)
        or external_review.get("reviewed_source_set_sha256")
        != _json_sha256(source_records)
        or not isinstance(structured, Mapping)
        or structured.get("reviewer") != "External Claude Opus"
        or structured.get("reviewer_agent_id")
        != "8540c6b8-509c-404c-9575-07662444c433"
        or structured.get("verdict") != "APPROVE"
        or structured.get("reviewed_source_set_sha256")
        != _json_sha256(source_records)
        or structured.get("pass") is not True
        or envelope.get("is_error") is not False
    ):
        raise RecoveryAuthorityError("Rejected v15 external review envelope drift")

    archive_contracts = (
        (
            "schema_recovery_authority_v15",
            REJECTED_V15_AUTHORITY_KEY_ARCHIVE,
            REJECTED_V15_ORIGINAL_AUTHORITY_KEY_ROOT,
            "ed25519_public.pem",
            "ed25519_private_one_time.pem",
            "797c9c174b72ca588d6a63b16c4ab1f0b1bf465763a4ee030e367e1e807aaf4d",
        ),
        (
            "downstream_terminal_v5",
            REJECTED_V15_TERMINAL_KEY_ARCHIVE,
            REJECTED_V15_ORIGINAL_TERMINAL_KEY_ROOT,
            "ed25519_public.pem",
            "ed25519_private_one_time_resident.pem",
            "17a0b67f22b706da120ca46ce37e9104cc53ba2fd524fd37c92b261668e00f84",
        ),
    )
    archived_keys: dict[str, Any] = {}
    for role, root, original_root, public_name, private_name, public_sha256 in (
        archive_contracts
    ):
        if (
            root.resolve(strict=True) != root
            or stat.S_IMODE(os.lstat(root).st_mode) != 0o700
            or os.path.lexists(original_root)
            or sorted(path.name for path in root.iterdir())
            != sorted([public_name, private_name, "UNUSED_KEY_ARCHIVE_RECORD.v1.json"])
        ):
            raise RecoveryAuthorityError(f"Rejected v15 key archive drift: {role}")
        _, public_record, public_fd = _open_bound(root / public_name)
        private_record, private_fd = _open_private_key_metadata_bound(
            root / private_name,
            expected_mode=0o400,
            state="archived_unused_never_signed",
        )
        status_data, status_record, _ = _open_bound(
            root / "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
        )
        status_payload = _load_json(status_data, f"rejected v15 key archive {role}")
        evidence_key = key_archive.get(role)
        if (
            public_record["sha256"] != public_sha256
            or public_record["size_bytes"] != 113
            or private_record["size_bytes"] != 119
            or os.fstat(public_fd).st_nlink != 1
            or os.fstat(private_fd).st_nlink != 1
            or not isinstance(evidence_key, Mapping)
            or evidence_key.get("archive") != str(root)
            or evidence_key.get("archive_record_sha256") != status_record["sha256"]
            or evidence_key.get("public_key_sha256") != public_sha256
            or evidence_key.get("signature_created") is not False
            or evidence_key.get("status") != "ARCHIVED_UNUSED_NEVER_SIGNED"
            or status_payload.get("role") != role
            or status_payload.get("archive") != str(root)
            or status_payload.get("source") != str(original_root)
            or status_payload.get("signature_created") is not False
            or status_payload.get("must_never_be_used") is not True
            or status_payload.get("reason")
            != "REJECTED_PRE_AUTHORITY_BY_INDEPENDENT_HIGH_AND_MEDIUM_FINDINGS"
        ):
            raise RecoveryAuthorityError(
                f"Rejected v15 archived key identity/lifecycle drift: {role}"
            )
        archived_keys[role] = {
            "archive": str(root),
            "public_key": public_record,
            "unused_private_key_metadata": private_record,
            "archive_record": status_record,
            "signature_created": False,
            "private_key_bytes_read": False,
            "pass": True,
        }

    occupied = [
        str(path) for path in REJECTED_V15_OUTPUT_SLOTS if os.path.lexists(path)
    ]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v15.bundle.staging.*"
        )
    )
    ledger_events = [
        _load_json(line, f"rejected v15 archive ledger event {index}")
        for index, line in enumerate(ledger_data.splitlines(), start=1)
    ]
    if (
        occupied
        or staging
        or len(ledger_events) != 6
        or [event.get("event") for event in ledger_events]
        != [
            "PREFLIGHT_COMPLETE",
            "ARCHIVE_MOVE_START",
            "ARCHIVE_MOVE_COMPLETE",
            "ARCHIVE_MOVE_START",
            "ARCHIVE_MOVE_COMPLETE",
            "ARCHIVE_BATCH_COMPLETE",
        ]
        or ledger_events[0].get("formal_v15_slot_count")
        != len(REJECTED_V15_OUTPUT_SLOTS)
        or ledger_events[0].get("formal_v15_slots_all_absent") is not True
        or ledger_events[1].get("role") != "schema_recovery_authority_v15"
        or ledger_events[2].get("role") != "schema_recovery_authority_v15"
        or ledger_events[3].get("role") != "downstream_terminal_v5"
        or ledger_events[4].get("role") != "downstream_terminal_v5"
        or not _strict_equal(
            ledger_events[-1],
            {
                "entry_count": 2,
                "event": "ARCHIVE_BATCH_COMPLETE",
                "status": "PASS",
            },
        )
    ):
        raise RecoveryAuthorityError(
            f"Rejected v15 archive/output contract drift: {occupied + staging}"
        )
    return {
        "evidence": evidence_record,
        "key_archive_ledger": ledger_record,
        "candidate_reviews": review_records,
        "external_review_envelope": envelope_record,
        "external_review_stderr": stderr_record,
        "sources": source_records,
        "archived_unused_keys": archived_keys,
        "authority_created": False,
        "formal_reviews_published": False,
        "signatures_created": False,
        "outputs_remain_absent": True,
        "strictest_review_wins": True,
        "pass": True,
    }


def _validate_rejected_v19_round1_pre_authority() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(REJECTED_V19_ROUND1_EVIDENCE)
    evidence = _load_json(evidence_data, "rejected v19 round1 evidence")
    if (
        evidence_record["size_bytes"] != 4_015
        or evidence_record["sha256"]
        != EXPECTED_REJECTED_V19_ROUND1_EVIDENCE_SHA256
        or evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_pre_authority_candidate_rejection"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("candidate") != "v19_round1"
        or evidence.get("status")
        != "REJECTED_BEFORE_AUTHORITY_BY_STRICTEST_REVIEW_WINS"
        or evidence.get("reviewed_source_set_sha256")
        != "30c6b0c47ecaf59d00dff31ba99bc536a5c75cde32dda15a8f20a7b52e1b99e9"
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v19 round1 evidence drift")

    source_records = _records(REJECTED_V19_ROUND1_SOURCE_PATHS)
    for role, record in source_records.items():
        expected_size, expected_sha256 = (
            EXPECTED_REJECTED_V19_ROUND1_SOURCE_IDENTITIES[role]
        )
        if (
            record["size_bytes"] != expected_size
            or record["sha256"] != expected_sha256
            or record["mode"] != "0o444"
        ):
            raise RecoveryAuthorityError(
                f"Rejected v19 round1 source identity drift: {role}"
            )

    expected_review_ids = {
        "mendel": "019f8ea2-2b9a-7271-90b2-0725bcc354a6",
        "nash": "019f8ea2-261f-7520-8a18-0d692b5002e0",
    }
    expected_medium_counts = {"mendel": 1, "nash": 2}
    review_records: dict[str, Any] = {}
    for role, path in REJECTED_V19_ROUND1_REVIEW_PATHS.items():
        review_data, review_record, _ = _open_bound(path)
        review = _load_json(review_data, f"rejected v19 round1 review {role}")
        evidence_review = evidence.get("reviews", {}).get(role, {})
        if (
            review.get("reviewer") != ("Mendel" if role == "mendel" else "Nash")
            or review.get("reviewer_agent_id") != expected_review_ids[role]
            or review.get("verdict") != "REQUEST_CHANGES"
            or review.get("reviewed_source_set_sha256")
            != evidence.get("reviewed_source_set_sha256")
            or len(review.get("medium_findings", [])) != expected_medium_counts[role]
            or review.get("pass") is not False
            or evidence_review.get("path") != str(path)
            or evidence_review.get("mode") != review_record["mode"]
            or evidence_review.get("size_bytes") != review_record["size_bytes"]
            or evidence_review.get("sha256") != review_record["sha256"]
            or evidence_review.get("verdict") != "REQUEST_CHANGES"
            or evidence_review.get("medium_finding_count")
            != expected_medium_counts[role]
        ):
            raise RecoveryAuthorityError(
                f"Rejected v19 round1 review drift: {role}"
            )
        review_records[role] = review_record

    envelope_data, envelope_record, _ = _open_bound(
        REJECTED_V19_ROUND1_EXTERNAL_ENVELOPE
    )
    envelope = json.loads(envelope_data)
    if not isinstance(envelope, Mapping):
        raise RecoveryAuthorityError("Rejected v19 round1 envelope is not an object")
    _, stderr_record, _ = _open_bound(REJECTED_V19_ROUND1_EXTERNAL_STDERR)
    transport_records = _records(REJECTED_V19_ROUND1_TRANSPORT_PATHS)
    expected_transport = {
        "external_prompt": (
            6_008,
            "311d3f580e6278cd6494f571a2e39e3076a8e4913b1d1f081ada49941a172a2a",
        ),
        "external_schema": (
            3_731,
            "12e0a38274072847bf5e3ddceab4d2c229bf48a29bc3077df4e511dcdff74d4d",
        ),
        "external_runner": (
            4_765,
            "4a5688ffb92beca8443a3c7ee8f4318aa4d09f9f4ae81c3ad4f72801405273b1",
        ),
    }
    for role, record in transport_records.items():
        if (record["size_bytes"], record["sha256"]) != expected_transport[role]:
            raise RecoveryAuthorityError(
                f"Rejected v19 round1 external transport drift: {role}"
            )
    aborted = evidence.get("aborted_external_transport", {})
    if (
        envelope_record["size_bytes"] != 1_246
        or envelope_record["sha256"]
        != "edd852da8473926d47c9bf535bc995e7608cb9c7cfd5c0ba04b33060b7e88ab6"
        or stderr_record["size_bytes"] != 302
        or stderr_record["sha256"]
        != "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a"
        or envelope.get("session_id") != "006a1050-2d27-4f63-9ee0-80a42fb8f26d"
        or envelope.get("subtype") != "error_during_execution"
        or envelope.get("terminal_reason") != "aborted_streaming"
        or aborted.get("formal_verdict_produced") is not False
    ):
        raise RecoveryAuthorityError("Rejected v19 round1 external abort drift")

    authority_state = evidence.get("authority_state", {})
    if (
        authority_state.get("authority_created") is not False
        or authority_state.get("signature_created") is not False
        or authority_state.get("canonical_reviews_published") is not False
        or authority_state.get("authority_public_key_sha256")
        != "d494f66e8aea206e37e0e803ccb4e0ceb9cf2b244e71eccba6b370f96ddee2e0"
        or authority_state.get("authority_private_key_mode") != "0o400"
        or authority_state.get("terminal_v9_private_key_mode") != "0o400"
        or len(evidence.get("blocking_findings", [])) != 2
        or len(evidence.get("required_corrections", [])) != 3
    ):
        raise RecoveryAuthorityError("Rejected v19 round1 lifecycle drift")
    return {
        "evidence": evidence_record,
        "sources": source_records,
        "candidate_reviews": review_records,
        "aborted_external_envelope": envelope_record,
        "aborted_external_stderr": stderr_record,
        "external_transport": transport_records,
        "authority_created": False,
        "signature_created": False,
        "strictest_review_wins": True,
        "findings_corrected_in_active_candidate": True,
        "pass": True,
    }


def _validate_rejected_v20_pre_authority() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(REJECTED_V20_EVIDENCE)
    evidence = _load_json(evidence_data, "rejected v20 pre-authority evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "candidate_generation",
            "status",
            "successor_generation",
            "reviewed_source_set_sha256",
            "source_set",
            "reviews",
            "strictest_review_summary",
            "pre_authority_state",
            "quarantined_unused_keys",
            "formal_output_slots",
            "scientific_claim_ceiling",
            "pass",
        },
        "rejected v20 pre-authority evidence",
    )
    if (
        evidence_record["size_bytes"] != 10_996
        or evidence_record["sha256"] != EXPECTED_REJECTED_V20_EVIDENCE_SHA256
        or evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_pre_authority_candidate_rejection"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("candidate_generation") != "v20"
        or evidence.get("status")
        != "REJECTED_BEFORE_AUTHORITY_BY_STRICTEST_REVIEW_WINS"
        or evidence.get("successor_generation") != "v21"
        or evidence.get("reviewed_source_set_sha256")
        != "595e4169194b87cace651f07d532e0a395be4d5fb2f434cdba519c07bec01467"
        or evidence.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v20 evidence identity/contract drift")

    source_records = _records(REJECTED_V20_SOURCE_PATHS)
    for role, record in source_records.items():
        expected_size, expected_sha256 = EXPECTED_REJECTED_V20_SOURCE_IDENTITIES[role]
        if (
            record["mode"] != "0o444"
            or record["size_bytes"] != expected_size
            or record["sha256"] != expected_sha256
        ):
            raise RecoveryAuthorityError(f"Rejected v20 source identity drift: {role}")
    if not _strict_equal(evidence.get("source_set"), source_records):
        raise RecoveryAuthorityError("Rejected v20 evidence/source inventory drift")

    expected_internal_reviews = {
        "mendel": {
            "reviewer": "Mendel",
            "reviewer_agent_id": "019f8ef0-be1f-7803-886c-cf7175fb17cc",
            "verdict": "REQUEST_CHANGES",
            "pass": False,
            "high_count": 0,
            "medium_count": 2,
            "sha256": "934991ba2b99681f38d5f7270e0c0fc8bc9c17864096f768aa56a4de2b52a5c3",
            "size_bytes": 2_906,
        },
        "nash": {
            "reviewer": "Nash",
            "reviewer_agent_id": "019f8ef0-b87c-78c2-92cc-cbd0b406f55d",
            "verdict": "REQUEST_CHANGES",
            "pass": False,
            "high_count": 1,
            "medium_count": 0,
            "sha256": "ae96aa680b0edcd48ed8a6483dfb6c1b668c39c743436ada8dfcc2cdf031cd31",
            "size_bytes": 2_913,
        },
    }
    review_records: dict[str, Any] = {}
    for role, path in REJECTED_V20_INTERNAL_REVIEW_PATHS.items():
        review_data, review_record, _ = _open_bound(path)
        review = _load_json(review_data, f"rejected v20 review {role}")
        expected = expected_internal_reviews[role]
        evidence_review = evidence.get("reviews", {}).get(role, {})
        if (
            review_record["mode"] != "0o444"
            or review_record["size_bytes"] != expected["size_bytes"]
            or review_record["sha256"] != expected["sha256"]
            or review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != SCHEMA_VERSION
            or review.get("reviewer") != expected["reviewer"]
            or review.get("reviewer_agent_id") != expected["reviewer_agent_id"]
            or review.get("verdict") != expected["verdict"]
            or review.get("pass") is not expected["pass"]
            or review.get("reviewed_source_set_sha256")
            != evidence["reviewed_source_set_sha256"]
            or len(review.get("high_findings", [])) != expected["high_count"]
            or len(review.get("medium_findings", [])) != expected["medium_count"]
            or evidence_review.get("reviewer_agent_id")
            != expected["reviewer_agent_id"]
            or evidence_review.get("verdict") != expected["verdict"]
            or not _strict_equal(evidence_review.get("artifact"), review_record)
        ):
            raise RecoveryAuthorityError(f"Rejected v20 review drift: {role}")
        review_records[role] = review_record

    envelope_data, envelope_record, _ = _open_bound(REJECTED_V20_EXTERNAL_ENVELOPE)
    envelope = json.loads(envelope_data)
    if not isinstance(envelope, Mapping):
        raise RecoveryAuthorityError("Rejected v20 external envelope is not an object")
    structured = envelope.get("structured_output")
    evidence_external = evidence.get("reviews", {}).get("external_claude_opus", {})
    _, stderr_record, _ = _open_bound(REJECTED_V20_EXTERNAL_STDERR)
    transport_records = _records(REJECTED_V20_TRANSPORT_PATHS)
    expected_transport = {
        "external_prompt": (
            6_627,
            "580b8d4e64a98478fb1136897764c7180b8f605bc031fb59f04eee6d4a2e1bf7",
        ),
        "external_schema": (
            3_731,
            "30ecc4da8bdc1237980ff401d474c36b6f7e1b57ed269d8f5a523db59b716635",
        ),
        "external_runner": (
            4_749,
            "4b370eec50f1e6417e9e074749ee59b4fcbd2570b56f4e425ad06279c1ff5eca",
        ),
    }
    for role, record in transport_records.items():
        if (record["size_bytes"], record["sha256"]) != expected_transport[role]:
            raise RecoveryAuthorityError(f"Rejected v20 transport drift: {role}")
    if (
        envelope_record["mode"] != "0o444"
        or envelope_record["size_bytes"] != 14_215
        or envelope_record["sha256"]
        != "00e192975175e4220a35de309e88a779de51818ba7d692450778b73f2a536a26"
        or stderr_record["size_bytes"] != 302
        or stderr_record["sha256"]
        != "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a"
        or envelope.get("is_error") is not False
        or envelope.get("session_id") != "84c32c61-b227-4378-b2fb-1169c0821db9"
        or not isinstance(structured, Mapping)
        or structured.get("reviewer") != "External Claude Opus"
        or structured.get("reviewer_agent_id")
        != "84c32c61-b227-4378-b2fb-1169c0821db9"
        or structured.get("verdict") != "APPROVE"
        or structured.get("pass") is not True
        or structured.get("reviewed_source_set_sha256")
        != evidence["reviewed_source_set_sha256"]
        or evidence_external.get("reviewer_agent_id")
        != "84c32c61-b227-4378-b2fb-1169c0821db9"
        or evidence_external.get("verdict") != "APPROVE"
        or not _strict_equal(evidence_external.get("artifact"), envelope_record)
    ):
        raise RecoveryAuthorityError("Rejected v20 external review drift")

    authority_public_data, authority_public_record, _ = _open_bound(
        REJECTED_V20_PUBLIC_KEY
    )
    authority_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V20_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V20_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    terminal_public_data, terminal_public_record, _ = _open_bound(
        REJECTED_V20_CONTINUATION_PUBLIC_KEY
    )
    terminal_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V20_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V20_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    evidence_keys = evidence.get("quarantined_unused_keys", {})
    evidence_authority = evidence_keys.get("authority", {})
    evidence_terminal = evidence_keys.get("terminal", {})
    if (
        len(authority_public_data) != 113
        or authority_public_record["sha256"]
        != "122fb38c395c37d1a4a2786c385110397db30f4b2db0ae3e4944f55355656fa9"
        or len(terminal_public_data) != 113
        or terminal_public_record["sha256"]
        != "09794c2ce162af3bf2f3117f6d11dea0c4bd626cbe50946267609058c6c0c291"
        or authority_private_record["size_bytes"] != 119
        or terminal_private_record["size_bytes"] != 119
        or evidence_authority.get("public_key") != authority_public_record
        or evidence_terminal.get("public_key") != terminal_public_record
        or evidence_authority.get("private_key", {}).get("path")
        != str(REJECTED_V20_PRIVATE_KEY)
        or evidence_authority.get("private_key", {}).get("mode") != "0o400"
        or evidence_authority.get("private_key", {}).get("size_bytes") != 119
        or evidence_authority.get("private_key", {}).get("sha256")
        != "cc6317a6e2fef86b3151315774ae5bfc44a3ff0282d23506ea91b4e5a3f1e259"
        or evidence_authority.get("private_key", {}).get("state")
        != "UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED"
        or evidence_terminal.get("private_key", {}).get("path")
        != str(REJECTED_V20_CONTINUATION_PRIVATE_KEY)
        or evidence_terminal.get("private_key", {}).get("mode") != "0o400"
        or evidence_terminal.get("private_key", {}).get("size_bytes") != 119
        or evidence_terminal.get("private_key", {}).get("sha256")
        != "3c708f24690697a3b714987cde2e05cdf0c9708d004f5022f3b84d151251c330"
        or evidence_terminal.get("private_key", {}).get("state")
        != "UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED"
    ):
        raise RecoveryAuthorityError("Rejected v20 key quarantine drift")

    output_contract = evidence.get("formal_output_slots", {})
    occupied = [str(path) for path in REJECTED_V20_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v20.bundle.staging.*"
        )
    )
    pre_authority = evidence.get("pre_authority_state", {})
    review_summary = evidence.get("strictest_review_summary", {})
    claim_ceiling = evidence.get("scientific_claim_ceiling", {})
    if (
        occupied
        or staging
        or output_contract.get("expected_count") != len(REJECTED_V20_OUTPUT_SLOTS)
        or output_contract.get("all_absent") is not True
        or output_contract.get("paths")
        != [str(path) for path in REJECTED_V20_OUTPUT_SLOTS]
        or pre_authority.get("authority_created") is not False
        or pre_authority.get("authority_signature_created") is not False
        or pre_authority.get("canonical_review_outputs_created") is not False
        or pre_authority.get("verification_receipt_created") is not False
        or pre_authority.get("replay_outputs_created") is not False
        or pre_authority.get("continuation_outputs_created") is not False
        or pre_authority.get("scientific_payload_changed") is not False
        or pre_authority.get("canonical_receipt_bytes_changed") is not False
        or pre_authority.get("claim_ceiling_changed") is not False
        or review_summary.get("approve_count") != 1
        or review_summary.get("request_changes_count") != 2
        or review_summary.get("blocking_high_count") != 1
        or review_summary.get("blocking_medium_count") != 2
        or len(review_summary.get("blocking_findings", [])) != 3
        or claim_ceiling.get("confirmed_cellular_subclones") != 0
        or claim_ceiling.get("linear_ancestry_calls") != 0
        or claim_ceiling.get("permitted_claim")
        != "latent molecular substructure candidates only"
    ):
        raise RecoveryAuthorityError(
            f"Rejected v20 lifecycle/output contract drift: {occupied + staging}"
        )
    return {
        "evidence": evidence_record,
        "sources": source_records,
        "candidate_reviews": review_records,
        "external_review_envelope": envelope_record,
        "external_review_stderr": stderr_record,
        "external_transport": transport_records,
        "authority_public_key": authority_public_record,
        "authority_private_key_metadata": authority_private_record,
        "terminal_public_key": terminal_public_record,
        "terminal_private_key_metadata": terminal_private_record,
        "authority_created": False,
        "signature_created": False,
        "formal_reviews_published": False,
        "outputs_remain_absent": True,
        "strictest_review_wins": True,
        "blocking_high_count": 1,
        "blocking_medium_count": 2,
        "findings_corrected_in_active_candidate": True,
        "pass": True,
    }


def _validate_rejected_v22_pre_authority() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(REJECTED_V22_EVIDENCE)
    evidence = _load_json(evidence_data, "rejected v22 pre-authority evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "candidate_generation",
            "status",
            "successor_generation",
            "reviewed_source_set_sha256",
            "source_set",
            "reviews",
            "review_attribution_semantics",
            "strictest_review_summary",
            "pre_authority_state",
            "quarantined_unused_keys",
            "formal_output_slots",
            "scientific_claim_ceiling",
            "pass",
        },
        "rejected v22 pre-authority evidence",
    )
    if (
        evidence_record["size_bytes"] != 11_084
        or evidence_record["sha256"] != EXPECTED_REJECTED_V22_EVIDENCE_SHA256
        or evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_pre_authority_candidate_rejection"
        or evidence.get("schema_version") != "1.1.0"
        or evidence.get("candidate_generation") != "v22"
        or evidence.get("status")
        != "REJECTED_BEFORE_AUTHORITY_BY_STRICTEST_REVIEW_WINS"
        or evidence.get("successor_generation") != "v23"
        or evidence.get("reviewed_source_set_sha256")
        != "c404f20adcfcfb61622a58cbf299099c6091e17f618c4df43556f064f9b4e7cd"
        or evidence.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v22 evidence identity/contract drift")

    source_records = _records(REJECTED_V22_SOURCE_PATHS)
    for role, record in source_records.items():
        expected_size, expected_sha256 = EXPECTED_REJECTED_V22_SOURCE_IDENTITIES[role]
        if (
            record["mode"] != "0o444"
            or record["size_bytes"] != expected_size
            or record["sha256"] != expected_sha256
        ):
            raise RecoveryAuthorityError(f"Rejected v22 source identity drift: {role}")
    if not _strict_equal(evidence.get("source_set"), source_records):
        raise RecoveryAuthorityError("Rejected v22 evidence/source inventory drift")

    review_records = _records(REJECTED_V22_REVIEW_PATHS)
    expected_review_identities = {
        "mendel": (1_962, "b6f768e3f312aed83a67d27db90070a966d02b6d0769ed33f3e4bc027a4ab720"),
        "nash": (1_701, "fd0a15c889d7af9b05d97cc569e15f8cb46102c75d83ce4bf90bf5c4cc565592"),
        "external_claude_opus": (12_118, "91b35d7a81a1e4aee56f85a06e71ccd40c8fcfda1a3b12446152d8f9dc7ae288"),
    }
    expected_reviewer_ids = {
        "mendel": "019f8f63-bee6-79f1-98dc-51a541e50fcf",
        "nash": "019f8f63-c4c3-7151-9fcb-ffe28df5a8cb",
        "external_claude_opus": "6f44c5f2-b8fb-49bb-a5cc-4987b3b30009",
    }
    for role, record in review_records.items():
        if (
            record["mode"] != "0o444"
            or (record["size_bytes"], record["sha256"])
            != expected_review_identities[role]
            or not _strict_equal(
                evidence.get("reviews", {}).get(role, {}).get("artifact"), record
            )
            or evidence.get("reviews", {}).get(role, {}).get("reviewer_agent_id")
            != expected_reviewer_ids[role]
        ):
            raise RecoveryAuthorityError(f"Rejected v22 review identity drift: {role}")

    for role in ("mendel", "nash"):
        review = _load_json(
            _open_bound(REJECTED_V22_REVIEW_PATHS[role])[0],
            f"rejected v22 captured review {role}",
        )
        if (
            review.get("schema_name")
            != "intersubmod.tumor_ref_schema_recovery_review_capture"
            or review.get("reviewer_agent_id") != expected_reviewer_ids[role]
            or review.get("verdict") != "REQUEST_CHANGES"
            or review.get("captured_by") != "coordinating_codex_session"
            or review.get("attribution_semantics")
            != "Orchestrator-recorded multi-agent transport attribution; not cryptographic proof of reviewer authorship."
            or len(review.get("high_findings", [])) != 1
            or len(review.get("medium_findings", [])) != (1 if role == "mendel" else 2)
            or review.get("pass") is not False
        ):
            raise RecoveryAuthorityError(f"Rejected v22 captured review drift: {role}")

    envelope = _load_archived_transport_json(
        _open_bound(REJECTED_V22_REVIEW_PATHS["external_claude_opus"])[0],
        "rejected v22 external envelope",
    )
    external_result = _load_json(
        envelope.get("result", ""), "rejected v22 external structured result"
    )
    if (
        envelope.get("session_id") != expected_reviewer_ids["external_claude_opus"]
        or envelope.get("is_error") is not False
        or external_result.get("reviewer") != "External Claude Opus"
        or external_result.get("reviewer_agent_id")
        != expected_reviewer_ids["external_claude_opus"]
        or external_result.get("verdict") != "APPROVE"
        or external_result.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v22 external review drift")

    transport_records = _records(REJECTED_V22_TRANSPORT_PATHS)
    expected_transport = {
        "external_prompt": (6_528, "8188d95b9af3f62e5813d843faa2d84ca33e81e0ec135130f99e564d2320d751"),
        "external_schema": (3_731, "094c902b2f544b821190577d1cc26c34b7b4e19b007fd833623267fcc7401592"),
        "external_runner": (4_749, "9202d2a91cbba0b397044e995fcc583c17358bfb748f7b17986d4ed8a50ed49e"),
        "external_stderr": (302, "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a"),
    }
    for role, record in transport_records.items():
        if (record["size_bytes"], record["sha256"]) != expected_transport[role]:
            raise RecoveryAuthorityError(f"Rejected v22 transport drift: {role}")

    authority_public_data, authority_public_record, _ = _open_bound(
        REJECTED_V22_PUBLIC_KEY
    )
    authority_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V22_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V22_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    terminal_public_data, terminal_public_record, _ = _open_bound(
        REJECTED_V22_CONTINUATION_PUBLIC_KEY
    )
    terminal_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V22_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V22_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    evidence_keys = evidence.get("quarantined_unused_keys", {})
    if (
        len(authority_public_data) != 113
        or authority_public_record["sha256"]
        != "a09adf97ee4b4afb2f90a0d12faf3db4ef1a3a0218676db333a6af6591c76051"
        or len(terminal_public_data) != 113
        or terminal_public_record["sha256"]
        != "94ebaec5e5fc994dc75ffcf189f124726154fccc89ff1ea32add8398f06bf5b2"
        or authority_private_record["size_bytes"] != 119
        or terminal_private_record["size_bytes"] != 119
        or evidence_keys.get("authority", {}).get("public_key")
        != authority_public_record
        or evidence_keys.get("terminal", {}).get("public_key")
        != terminal_public_record
        or evidence_keys.get("authority", {}).get("private_key", {}).get("sha256")
        != "bdcc806d6da6879e0625c995a1da4a780f1a79d36b16f75c287781d38e925fc1"
        or evidence_keys.get("terminal", {}).get("private_key", {}).get("sha256")
        != "435a3d58320cdf6cdb2d4159b6d1a71178fbcb211b9c38cebbeaa69a2c7837b5"
    ):
        raise RecoveryAuthorityError("Rejected v22 key quarantine drift")

    output_contract = evidence.get("formal_output_slots", {})
    occupied = [str(path) for path in REJECTED_V22_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v22.bundle.staging.*"
        )
    )
    pre_authority = evidence.get("pre_authority_state", {})
    review_summary = evidence.get("strictest_review_summary", {})
    claim_ceiling = evidence.get("scientific_claim_ceiling", {})
    if (
        occupied
        or staging
        or output_contract.get("expected_count") != len(REJECTED_V22_OUTPUT_SLOTS)
        or output_contract.get("all_absent") is not True
        or output_contract.get("paths") != [str(path) for path in REJECTED_V22_OUTPUT_SLOTS]
        or any(pre_authority.get(field) is not False for field in pre_authority)
        or review_summary.get("approve_count") != 1
        or review_summary.get("request_changes_count") != 2
        or review_summary.get("blocking_high_count") != 2
        or review_summary.get("blocking_medium_count") != 3
        or claim_ceiling.get("confirmed_cellular_subclones") != 0
        or claim_ceiling.get("linear_ancestry_calls") != 0
        or claim_ceiling.get("permitted_claim")
        != "latent molecular substructure candidates only"
        or "not cryptographic proof of reviewer authorship"
        not in evidence.get("review_attribution_semantics", "")
    ):
        raise RecoveryAuthorityError(
            f"Rejected v22 lifecycle/output contract drift: {occupied + staging}"
        )
    return {
        "evidence": evidence_record,
        "sources": source_records,
        "candidate_reviews": review_records,
        "external_transport": transport_records,
        "authority_public_key": authority_public_record,
        "authority_private_key_metadata": authority_private_record,
        "terminal_public_key": terminal_public_record,
        "terminal_private_key_metadata": terminal_private_record,
        "authority_created": False,
        "signature_created": False,
        "formal_reviews_published": False,
        "outputs_remain_absent": True,
        "strictest_review_wins": True,
        "review_authorship_cryptographically_proven": False,
        "findings_corrected_in_active_candidate": True,
        "pass": True,
    }


def _validate_rejected_v23_pre_authority() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(REJECTED_V23_EVIDENCE)
    evidence = _load_json(evidence_data, "rejected v23 pre-authority evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "candidate_generation",
            "status",
            "successor_generation",
            "reviewed_source_set_sha256",
            "source_set",
            "reviews",
            "review_attribution_semantics",
            "strictest_review_summary",
            "pre_authority_state",
            "quarantined_unused_keys",
            "formal_output_slots",
            "scientific_claim_ceiling",
            "pass",
        },
        "rejected v23 pre-authority evidence",
    )
    if (
        evidence_record["size_bytes"] != 11_042
        or evidence_record["sha256"] != EXPECTED_REJECTED_V23_EVIDENCE_SHA256
        or evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_pre_authority_candidate_rejection"
        or evidence.get("schema_version") != "1.1.0"
        or evidence.get("candidate_generation") != "v23"
        or evidence.get("status")
        != "REJECTED_BEFORE_AUTHORITY_BY_STRICTEST_REVIEW_WINS"
        or evidence.get("successor_generation") != "v24"
        or evidence.get("reviewed_source_set_sha256")
        != "7304704c23a5cf73c6a4f83f28345a7999a57eddfe507b350f1e894affab25fd"
        or evidence.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v23 evidence identity/contract drift")

    source_records = _records(REJECTED_V23_SOURCE_PATHS)
    for role, record in source_records.items():
        expected_size, expected_sha256 = EXPECTED_REJECTED_V23_SOURCE_IDENTITIES[role]
        if (
            record["mode"] != "0o444"
            or record["size_bytes"] != expected_size
            or record["sha256"] != expected_sha256
        ):
            raise RecoveryAuthorityError(f"Rejected v23 source identity drift: {role}")
    if not _strict_equal(evidence.get("source_set"), source_records):
        raise RecoveryAuthorityError("Rejected v23 evidence/source inventory drift")

    review_records = _records(REJECTED_V23_REVIEW_PATHS)
    expected_review_identities = {
        "mendel": (1_806, "b588cea29c4532dff26bd6956c588a87cd1e0f9ed4dec11b9355263f945beb0a"),
        "nash": (1_408, "d479bbc2c5a99d4a6559fe6c944a362ef04dd3b6c912c8e610d94da56b0bc3ca"),
        "external_claude_opus": (16_635, "54fa7c1e3cfcbfee52d6f9c6146fa8a3e4f65e63e9a2a4d10cbc4b9ef8a877b3"),
    }
    expected_reviewer_ids = {
        "mendel": "019f8f95-df75-70c3-a43f-58a7b2b9661d",
        "nash": "019f8f95-e5dd-7f33-96c9-ba74b5f96743",
        "external_claude_opus": "8ddaf10e-14e4-4e03-ae73-92c0c7051e69",
    }
    for role, record in review_records.items():
        if (
            record["mode"] != "0o444"
            or (record["size_bytes"], record["sha256"])
            != expected_review_identities[role]
            or not _strict_equal(
                evidence.get("reviews", {}).get(role, {}).get("artifact"), record
            )
            or evidence.get("reviews", {}).get(role, {}).get("reviewer_agent_id")
            != expected_reviewer_ids[role]
        ):
            raise RecoveryAuthorityError(f"Rejected v23 review identity drift: {role}")

    for role, medium_count in (("mendel", 1), ("nash", 0)):
        review = _load_json(
            _open_bound(REJECTED_V23_REVIEW_PATHS[role])[0],
            f"rejected v23 captured review {role}",
        )
        if (
            review.get("schema_name")
            != "intersubmod.tumor_ref_schema_recovery_review_capture"
            or review.get("reviewer_agent_id") != expected_reviewer_ids[role]
            or review.get("verdict") != "REQUEST_CHANGES"
            or review.get("captured_by") != "coordinating_codex_session"
            or review.get("attribution_semantics")
            != "Orchestrator-recorded multi-agent transport attribution; not cryptographic proof of reviewer authorship."
            or len(review.get("high_findings", [])) != 1
            or len(review.get("medium_findings", [])) != medium_count
            or review.get("pass") is not False
        ):
            raise RecoveryAuthorityError(f"Rejected v23 captured review drift: {role}")

    envelope = _load_archived_transport_json(
        _open_bound(REJECTED_V23_REVIEW_PATHS["external_claude_opus"])[0],
        "rejected v23 external envelope",
    )
    external_result = envelope.get("structured_output", {})
    if (
        envelope.get("session_id") != expected_reviewer_ids["external_claude_opus"]
        or envelope.get("is_error") is not False
        or not isinstance(external_result, Mapping)
        or external_result.get("reviewer") != "External Claude Opus"
        or external_result.get("reviewer_agent_id")
        != expected_reviewer_ids["external_claude_opus"]
        or external_result.get("verdict") != "APPROVE"
        or external_result.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v23 external review drift")

    transport_records = _records(REJECTED_V23_TRANSPORT_PATHS)
    expected_transport = {
        "external_prompt": (7_053, "71d11af166543e9a22dafc2cb5b0c4cf62edd9c96658c2e3d71493788ae1edf9"),
        "external_schema": (4_332, "2f087ef7a97157ab48a7969c00fef70936fd8639ed656ed7844f01fb4ad63ce7"),
        "external_runner": (4_749, "1cd5f457405de359f6a2c4f350b270e79596dfd4a66d369a61c26bf4b402d8a7"),
        "external_stderr": (302, "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a"),
    }
    for role, record in transport_records.items():
        if (record["size_bytes"], record["sha256"]) != expected_transport[role]:
            raise RecoveryAuthorityError(f"Rejected v23 transport drift: {role}")

    authority_public_data, authority_public_record, _ = _open_bound(
        REJECTED_V23_PUBLIC_KEY
    )
    authority_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V23_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V23_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    terminal_public_data, terminal_public_record, _ = _open_bound(
        REJECTED_V23_CONTINUATION_PUBLIC_KEY
    )
    terminal_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V23_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V23_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    evidence_keys = evidence.get("quarantined_unused_keys", {})
    if (
        len(authority_public_data) != 113
        or authority_public_record["sha256"]
        != "468fb1526e08b4e2cee44ee27d3cd20438beec2298769fca2a413e164759fbe1"
        or len(terminal_public_data) != 113
        or terminal_public_record["sha256"]
        != "d050c0dfea29b469e271967ea7759d4cbb36c3cbec493010163be5938e81e54a"
        or authority_private_record["size_bytes"] != 119
        or terminal_private_record["size_bytes"] != 119
        or evidence_keys.get("authority", {}).get("public_key")
        != authority_public_record
        or evidence_keys.get("terminal", {}).get("public_key")
        != terminal_public_record
        or evidence_keys.get("authority", {}).get("private_key", {}).get("sha256")
        != "a83d4369c3e3a06f421559a5a9deec38b1311ba2565e5d4296aef5e7feeea6c2"
        or evidence_keys.get("terminal", {}).get("private_key", {}).get("sha256")
        != "84c75838d2563f1a2668929e829f5b30a288cb37fd217d70c1c0651a70f80b03"
    ):
        raise RecoveryAuthorityError("Rejected v23 key quarantine drift")

    output_contract = evidence.get("formal_output_slots", {})
    occupied = [str(path) for path in REJECTED_V23_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v23.bundle.staging.*"
        )
    )
    pre_authority = evidence.get("pre_authority_state", {})
    review_summary = evidence.get("strictest_review_summary", {})
    claim_ceiling = evidence.get("scientific_claim_ceiling", {})
    if (
        occupied
        or staging
        or output_contract.get("expected_count") != len(REJECTED_V23_OUTPUT_SLOTS)
        or output_contract.get("all_absent") is not True
        or output_contract.get("paths") != [str(path) for path in REJECTED_V23_OUTPUT_SLOTS]
        or any(pre_authority.get(field) is not False for field in pre_authority)
        or review_summary.get("approve_count") != 1
        or review_summary.get("request_changes_count") != 2
        or review_summary.get("blocking_high_count") != 2
        or review_summary.get("blocking_medium_count") != 1
        or claim_ceiling.get("confirmed_cellular_subclones") != 0
        or claim_ceiling.get("linear_ancestry_calls") != 0
        or claim_ceiling.get("permitted_claim")
        != "latent molecular substructure candidates only"
        or "not cryptographic proof of reviewer authorship"
        not in evidence.get("review_attribution_semantics", "")
    ):
        raise RecoveryAuthorityError(
            f"Rejected v23 lifecycle/output contract drift: {occupied + staging}"
        )
    return {
        "evidence": evidence_record,
        "sources": source_records,
        "candidate_reviews": review_records,
        "external_transport": transport_records,
        "authority_public_key": authority_public_record,
        "authority_private_key_metadata": authority_private_record,
        "terminal_public_key": terminal_public_record,
        "terminal_private_key_metadata": terminal_private_record,
        "authority_created": False,
        "signature_created": False,
        "formal_reviews_published": False,
        "outputs_remain_absent": True,
        "strictest_review_wins": True,
        "review_authorship_cryptographically_proven": False,
        "findings_corrected_in_active_candidate": True,
        "pass": True,
    }


def _validate_rejected_v24_pre_authority() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(REJECTED_V24_EVIDENCE)
    evidence = _load_json(evidence_data, "rejected v24 pre-authority evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "candidate_generation",
            "status",
            "successor_generation",
            "reviewed_source_set_sha256",
            "source_set",
            "reviews",
            "review_transport_evidence",
            "review_attribution_semantics",
            "strictest_review_summary",
            "pre_authority_state",
            "quarantined_unused_keys",
            "formal_output_slots",
            "scientific_claim_ceiling",
            "pass",
        },
        "rejected v24 pre-authority evidence",
    )
    if (
        evidence_record["size_bytes"] != 13_158
        or evidence_record["sha256"] != EXPECTED_REJECTED_V24_EVIDENCE_SHA256
        or evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_pre_authority_candidate_rejection"
        or evidence.get("schema_version") != "1.1.0"
        or evidence.get("candidate_generation") != "v24"
        or evidence.get("status")
        != "REJECTED_BEFORE_AUTHORITY_BY_STRICTEST_REVIEW_WINS"
        or evidence.get("successor_generation") != "v25"
        or evidence.get("reviewed_source_set_sha256")
        != "dd4287bdd37a27ddf6bb5782806a814bebd39d2382d62c027f93c3cf81cf665b"
        or evidence.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v24 evidence identity/contract drift")

    source_records = _records(REJECTED_V24_SOURCE_PATHS)
    for role, record in source_records.items():
        expected_size, expected_sha256 = EXPECTED_REJECTED_V24_SOURCE_IDENTITIES[role]
        if (
            record["mode"] != "0o444"
            or record["size_bytes"] != expected_size
            or record["sha256"] != expected_sha256
        ):
            raise RecoveryAuthorityError(f"Rejected v24 source identity drift: {role}")
    if not _strict_equal(evidence.get("source_set"), source_records):
        raise RecoveryAuthorityError("Rejected v24 evidence/source inventory drift")

    review_records = _records(REJECTED_V24_REVIEW_PATHS)
    expected_review_identities = {
        "mendel": (3_979, "c3e11dc27f36bb1c42ea1c27010edaa5aa8361f6ebc3a36339485771a2cb8e3a"),
        "nash": (4_826, "404f40e7111b6b46f3c1292ccd8cf78c057312d250872ddf38c8b3eef80d672c"),
        "external_claude_opus": (19_758, "e4b67862fd33d2b37c00a6281e804f563e307a66cbcc27a48f0078b22619f340"),
    }
    expected_reviewer_ids = {
        "mendel": "019f8f95-df75-70c3-a43f-58a7b2b9661d",
        "nash": "019f8f95-e5dd-7f33-96c9-ba74b5f96743",
        "external_claude_opus": "8ddaf10e-14e4-4e03-ae73-92c0c7051e69",
    }
    for role, record in review_records.items():
        if (
            record["mode"] != "0o444"
            or (record["size_bytes"], record["sha256"])
            != expected_review_identities[role]
            or not _strict_equal(
                evidence.get("reviews", {}).get(role, {}).get("artifact"), record
            )
            or evidence.get("reviews", {}).get(role, {}).get("reviewer_agent_id")
            != expected_reviewer_ids[role]
        ):
            raise RecoveryAuthorityError(f"Rejected v24 review identity drift: {role}")

    for role, medium_count in (("mendel", 1), ("nash", 3)):
        review = _load_json(
            _open_bound(REJECTED_V24_REVIEW_PATHS[role])[0],
            f"rejected v24 captured review {role}",
        )
        if (
            review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != "1.0.0"
            or review.get("reviewer_agent_id") != expected_reviewer_ids[role]
            or review.get("reviewed_source_set_sha256")
            != evidence.get("reviewed_source_set_sha256")
            or review.get("verdict") != "REQUEST_CHANGES"
            or review.get("high_findings") != []
            or len(review.get("medium_findings", [])) != medium_count
            or review.get("pass") is not False
        ):
            raise RecoveryAuthorityError(f"Rejected v24 captured review drift: {role}")

    envelope = _load_archived_transport_json(
        _open_bound(REJECTED_V24_REVIEW_PATHS["external_claude_opus"])[0],
        "rejected v24 external envelope",
    )
    external_result = envelope.get("structured_output", {})
    if (
        envelope.get("session_id") != expected_reviewer_ids["external_claude_opus"]
        or envelope.get("is_error") is not False
        or not isinstance(external_result, Mapping)
        or external_result.get("reviewer") != "External Claude Opus"
        or external_result.get("reviewer_agent_id")
        != expected_reviewer_ids["external_claude_opus"]
        or external_result.get("reviewed_source_set_sha256")
        != evidence.get("reviewed_source_set_sha256")
        or external_result.get("verdict") != "APPROVE"
        or external_result.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v24 external review drift")

    transport_records = _records(REJECTED_V24_TRANSPORT_PATHS)
    expected_transport = {
        "external_prompt": (7_217, "89683e8a831780d24a0e28c85307701aeb6b5bf369a2de768b1f0622307cdcfa"),
        "external_schema": (4_332, "e47ddd8f41086038a9b64bafb3f0c0fee11db2dbc03e1738ecfa6cd4a93ad053"),
        "external_runner": (4_748, "4a0b1c2960141798b44486c5b6b3eb9bd89ea024fae7792f64a2bfe942d265e3"),
        "external_stderr": (302, "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a"),
    }
    for role, record in transport_records.items():
        if (
            record["mode"] != "0o444"
            or (record["size_bytes"], record["sha256"])
            != expected_transport[role]
            or not _strict_equal(
                evidence.get("review_transport_evidence", {}).get(role), record
            )
        ):
            raise RecoveryAuthorityError(f"Rejected v24 transport drift: {role}")

    authority_public_data, authority_public_record, _ = _open_bound(
        REJECTED_V24_PUBLIC_KEY
    )
    authority_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V24_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V24_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    terminal_public_data, terminal_public_record, _ = _open_bound(
        REJECTED_V24_CONTINUATION_PUBLIC_KEY
    )
    terminal_private_record, _ = _open_private_key_metadata_bound(
        REJECTED_V24_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="REJECTED_V24_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
    )
    evidence_keys = evidence.get("quarantined_unused_keys", {})
    if (
        len(authority_public_data) != 113
        or authority_public_record["sha256"]
        != "c15734463d661be0ead9a8292a3c0c23bd3391233882606df0e930eabb933879"
        or len(terminal_public_data) != 113
        or terminal_public_record["sha256"]
        != "dcff9fbb753ebdda8525597354be18be7b3fa3d22a4f40a020cd9f677736cc8a"
        or authority_private_record["size_bytes"] != 119
        or terminal_private_record["size_bytes"] != 119
        or evidence_keys.get("authority", {}).get("public_key")
        != authority_public_record
        or evidence_keys.get("terminal", {}).get("public_key")
        != terminal_public_record
        or evidence_keys.get("authority", {}).get("private_key", {}).get("sha256")
        != "5e36389dc4f73474f15bfbb737a4af4d0d2366e09d33315a48461d26c3011e1f"
        or evidence_keys.get("terminal", {}).get("private_key", {}).get("sha256")
        != "a12456bace4618c8e93c2c252885e88e9a80f41a4a766db9aebb4e16d00e3d31"
    ):
        raise RecoveryAuthorityError("Rejected v24 key quarantine drift")

    output_contract = evidence.get("formal_output_slots", {})
    occupied = [str(path) for path in REJECTED_V24_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v24.bundle.staging.*"
        )
    )
    pre_authority = evidence.get("pre_authority_state", {})
    review_summary = evidence.get("strictest_review_summary", {})
    claim_ceiling = evidence.get("scientific_claim_ceiling", {})
    if (
        occupied
        or staging
        or output_contract.get("expected_count") != len(REJECTED_V24_OUTPUT_SLOTS)
        or output_contract.get("all_absent") is not True
        or output_contract.get("paths") != [str(path) for path in REJECTED_V24_OUTPUT_SLOTS]
        or any(pre_authority.get(field) is not False for field in pre_authority)
        or review_summary.get("approve_count") != 1
        or review_summary.get("request_changes_count") != 2
        or review_summary.get("blocking_high_count") != 0
        or review_summary.get("blocking_medium_count") != 4
        or review_summary.get("unique_blocking_finding_count") != 3
        or claim_ceiling.get("confirmed_cellular_subclones") != 0
        or claim_ceiling.get("linear_ancestry_calls") != 0
        or claim_ceiling.get("permitted_claim")
        != "latent molecular substructure candidates only"
        or "not cryptographic proof of reviewer authorship"
        not in evidence.get("review_attribution_semantics", "")
    ):
        raise RecoveryAuthorityError(
            f"Rejected v24 lifecycle/output contract drift: {occupied + staging}"
        )
    return {
        "evidence": evidence_record,
        "sources": source_records,
        "candidate_reviews": review_records,
        "external_transport": transport_records,
        "authority_public_key": authority_public_record,
        "authority_private_key_metadata": authority_private_record,
        "terminal_public_key": terminal_public_record,
        "terminal_private_key_metadata": terminal_private_record,
        "authority_created": False,
        "signature_created": False,
        "formal_reviews_published": False,
        "outputs_remain_absent": True,
        "strictest_review_wins": True,
        "review_authorship_cryptographically_proven": False,
        "findings_corrected_in_active_candidate": True,
        "pass": True,
    }


def _validate_rejected_v27_pre_authority() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(REJECTED_V27_EVIDENCE)
    evidence = _load_json(evidence_data, "rejected v27 pre-authority evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "task_type",
            "generation",
            "status",
            "adjudication",
            "root_causes",
            "source_artifacts",
            "review_artifacts",
            "readonly_probe",
            "key_archives",
            "formal_output_contract",
            "scientific_scope",
            "pass",
        },
        "rejected v27 pre-authority evidence",
    )
    if (
        evidence_record["size_bytes"] != 4_778
        or evidence_record["sha256"] != EXPECTED_REJECTED_V27_EVIDENCE_SHA256
        or evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_pre_authority_rejection"
        or evidence.get("schema_version") != "1.0.0"
        or evidence.get("task_type") != "B_comprehensive_validation"
        or evidence.get("generation") != "v27"
        or evidence.get("status") != "REJECTED_PRE_AUTHORITY_ARCHIVED"
        or evidence.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v27 evidence identity/contract drift")

    adjudication = evidence.get("adjudication")
    root_causes = evidence.get("root_causes")
    if (
        not isinstance(adjudication, Mapping)
        or adjudication.get("verdict") != "REJECTED_PRE_AUTHORITY"
        or adjudication.get("policy")
        != "any_reproducible_high_finding_blocks_signature_even_if_another_reviewer_approves"
        or adjudication.get("mendel") != "REQUEST_CHANGES"
        or adjudication.get("nash") != "REQUEST_CHANGES"
        or adjudication.get("external_claude_opus") != "APPROVE"
        or adjudication.get("signature_created") is not False
        or root_causes
        != [
            "recovery_result_report_finalizer_and_recovery_report_builder_runtime_contracts_retained_v26_sha256_values",
            "formal_internal_review_transport_required_multi_agent_v1_but_requested_reviews_declared_multi_agent",
            "claimed_all_26_staging_pattern_regressions_but_parameterized_test_covered_only_indices_0_through_23",
        ]
    ):
        raise RecoveryAuthorityError("Rejected v27 adjudication drift")

    source_records = _records(REJECTED_V27_SOURCE_PATHS)
    evidence_sources = evidence.get("source_artifacts")
    if not isinstance(evidence_sources, Mapping) or set(evidence_sources) != set(
        REJECTED_V27_SOURCE_PATHS
    ):
        raise RecoveryAuthorityError("Rejected v27 source role set drift")
    for role, record in source_records.items():
        expected_size, expected_sha256 = EXPECTED_REJECTED_V27_SOURCE_IDENTITIES[role]
        evidence_source = evidence_sources.get(role)
        expected_relative_path = str(REJECTED_V27_SOURCE_PATHS[role].relative_to(REJECTED_V27_ROOT))
        if (
            record["mode"] != "0o444"
            or (record["size_bytes"], record["sha256"])
            != (expected_size, expected_sha256)
            or not isinstance(evidence_source, Mapping)
            or evidence_source.get("mode") != "0o444"
            or evidence_source.get("path") != expected_relative_path
            or evidence_source.get("size_bytes") != expected_size
            or evidence_source.get("sha256") != expected_sha256
        ):
            raise RecoveryAuthorityError(f"Rejected v27 source identity drift: {role}")

    review_records = _records(REJECTED_V27_REVIEW_PATHS)
    expected_review_contracts = {
        "mendel": ("Mendel", "019f909e-07a6-7af1-93fb-8bf5b2a4b69b", "REQUEST_CHANGES", False, 2, 0),
        "nash": ("Nash", "019f909e-0e85-71b3-8cdb-736098200420", "REQUEST_CHANGES", False, 1, 2),
        "external_claude_opus": ("External Claude Opus", "949fbc0f-1af7-4186-9ac6-eefbeb7e76be", "APPROVE", True, 0, 0),
    }
    evidence_reviews = evidence.get("review_artifacts")
    for role, record in review_records.items():
        expected_size, expected_sha256 = EXPECTED_REJECTED_V27_REVIEW_IDENTITIES[role]
        reviewer, reviewer_id, verdict, passed, high_count, medium_count = (
            expected_review_contracts[role]
        )
        payload = _load_json(
            _open_bound(REJECTED_V27_REVIEW_PATHS[role])[0],
            f"rejected v27 review {role}",
        )
        evidence_review = evidence_reviews.get(role) if isinstance(evidence_reviews, Mapping) else None
        attribution = payload.get("attribution")
        if (
            record["mode"] != "0o444"
            or (record["size_bytes"], record["sha256"])
            != (expected_size, expected_sha256)
            or not isinstance(evidence_review, Mapping)
            or evidence_review.get("sha256") != expected_sha256
            or evidence_review.get("verdict") != verdict
            or payload.get("schema_name") != REVIEW_SCHEMA_NAME
            or payload.get("schema_version") != SCHEMA_VERSION
            or payload.get("reviewer") != reviewer
            or payload.get("reviewer_agent_id") != reviewer_id
            or payload.get("verdict") != verdict
            or payload.get("pass") is not passed
            or len(payload.get("high_findings", [])) != high_count
            or len(payload.get("medium_findings", [])) != medium_count
            or not isinstance(attribution, Mapping)
            or attribution.get("cryptographic_reviewer_authorship_proven") is not False
            or attribution.get("transport_id") != reviewer_id
            or attribution.get("transport")
            != ("claude_cli_session_envelope" if role == "external_claude_opus" else "multi_agent")
        ):
            raise RecoveryAuthorityError(f"Rejected v27 review drift: {role}")

    transport_records = _records(REJECTED_V27_TRANSPORT_PATHS)
    for role, record in transport_records.items():
        if (
            record["mode"] != "0o444"
            or (record["size_bytes"], record["sha256"])
            != EXPECTED_REJECTED_V27_TRANSPORT_IDENTITIES[role]
        ):
            raise RecoveryAuthorityError(f"Rejected v27 transport drift: {role}")
    envelope = _load_archived_transport_json(
        _open_bound(REJECTED_V27_TRANSPORT_PATHS["external_envelope"])[0],
        "rejected v27 external envelope",
    )
    if (
        envelope.get("session_id") != expected_review_contracts["external_claude_opus"][1]
        or not _strict_equal(
            envelope.get("structured_output"),
            _load_json(
                _open_bound(REJECTED_V27_REVIEW_PATHS["external_claude_opus"])[0],
                "rejected v27 external review",
            ),
        )
    ):
        raise RecoveryAuthorityError("Rejected v27 external envelope drift")

    archive_contracts = (
        (
            "schema_recovery_authority_v27",
            REJECTED_V27_AUTHORITY_KEY_ARCHIVE,
            REJECTED_V27_ORIGINAL_AUTHORITY_KEY_ROOT,
            "ed25519_public.pem",
            "ed25519_private_one_time.pem",
            "db8c8b5c5ea1a7e4235efb59966c6d41981ed9b6085af6226c07226b56ab3cb3",
            "1e653a445197726e8ffa2d66361708814a3a9b6e27bb7c129d9c2c314f7d139a",
        ),
        (
            "downstream_terminal_v17",
            REJECTED_V27_TERMINAL_KEY_ARCHIVE,
            REJECTED_V27_ORIGINAL_TERMINAL_KEY_ROOT,
            "ed25519_public.pem",
            "ed25519_private_one_time_resident.pem",
            "355979601138f8dca29534db58e3862f30f30b49526c01d10a97b01ba91c26f9",
            "59ee02af13e34d0c2394cae45f915972d86e0078223c9c35d5af454a28a0c3cc",
        ),
    )
    archived_keys: dict[str, Any] = {}
    for role, root, original_root, public_name, private_name, public_sha256, record_sha256 in archive_contracts:
        if (
            root.resolve(strict=True) != root
            or stat.S_IMODE(os.lstat(root).st_mode) != 0o700
            or os.path.lexists(original_root)
            or sorted(path.name for path in root.iterdir())
            != sorted([public_name, private_name, "UNUSED_KEY_ARCHIVE_RECORD.v1.json"])
        ):
            raise RecoveryAuthorityError(f"Rejected v27 key archive drift: {role}")
        _, public_record, public_fd = _open_bound(root / public_name)
        private_record, private_fd = _open_private_key_metadata_bound(
            root / private_name,
            expected_mode=0o400,
            state="REJECTED_V27_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
        )
        status_data, status_record, _ = _open_bound(
            root / "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
        )
        status_payload = _load_json(status_data, f"rejected v27 key archive {role}")
        if (
            public_record["size_bytes"] != 113
            or public_record["sha256"] != public_sha256
            or private_record["size_bytes"] != 119
            or os.fstat(public_fd).st_nlink != 1
            or os.fstat(private_fd).st_nlink != 1
            or status_record["sha256"] != record_sha256
            or status_payload.get("role") != role
            or status_payload.get("archive") != str(root)
            or status_payload.get("source") != str(original_root)
            or status_payload.get("signature_created") is not False
            or status_payload.get("must_never_be_used") is not True
            or status_payload.get("reason")
            != "REJECTED_PRE_AUTHORITY_BY_INDEPENDENT_HIGH_AND_MEDIUM_FINDINGS"
        ):
            raise RecoveryAuthorityError(f"Rejected v27 key lifecycle drift: {role}")
        archived_keys[role] = {
            "archive": str(root),
            "public_key": public_record,
            "unused_private_key_metadata": private_record,
            "archive_record": status_record,
            "private_key_bytes_read": False,
            "pass": True,
        }

    occupied = [str(path) for path in REJECTED_V27_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v27.bundle.staging.*"
        )
    )
    probe = evidence.get("readonly_probe")
    formal = evidence.get("formal_output_contract")
    science = evidence.get("scientific_scope")
    if (
        occupied
        or staging
        or not isinstance(probe, Mapping)
        or probe.get("exit_code") != 0
        or probe.get("forbidden_output_slots_checked") != 390
        or probe.get("staging_patterns_absent") != 26
        or probe.get("regression_summary") != "656 passed"
        or probe.get("no_output_writes") is not True
        or not isinstance(formal, Mapping)
        or any(value is not False for value in formal.values())
        or not isinstance(science, Mapping)
        or science.get("datasets") != 7
        or science.get("same_run_longphase_s_recalibrated_filter_pass_ssnv") != 469_849
        or science.get("scientific_payload_changed") is not False
        or science.get("claim_ceiling_changed") is not False
        or science.get("confirmed_cellular_subclones") != 0
        or science.get("linear_ancestry_calls") != 0
    ):
        raise RecoveryAuthorityError(
            f"Rejected v27 lifecycle/output contract drift: {occupied + staging}"
        )
    return {
        "evidence": evidence_record,
        "sources": source_records,
        "candidate_reviews": review_records,
        "review_transport": transport_records,
        "archived_unused_keys": archived_keys,
        "authority_created": False,
        "signature_created": False,
        "formal_reviews_published": False,
        "outputs_remain_absent": True,
        "strictest_reproducible_finding_wins": True,
        "review_authorship_cryptographically_proven": False,
        "findings_corrected_in_active_candidate": True,
        "pass": True,
    }


def _validate_rejected_v29_round1_pre_authority() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(REJECTED_V29_ROUND1_EVIDENCE)
    summary_data, summary_record, _ = _open_bound(REJECTED_V29_ROUND1_SUMMARY)
    evidence = _load_json(evidence_data, "rejected v29 round1 evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "candidate",
            "task_type",
            "status",
            "reviewed_source_set_sha256",
            "inventory_file_count",
            "inventory",
            "review_chain",
            "blocking_reproduction",
            "required_corrections",
            "authority_state",
            "active_public_key_sha256",
            "scientific_payload_changed",
            "claim_ceiling_changed",
            "claim_ceiling",
            "pass",
        },
        "rejected v29 round1 evidence",
    )
    if (
        evidence_record["size_bytes"] != 6_363
        or evidence_record["sha256"]
        != EXPECTED_REJECTED_V29_ROUND1_EVIDENCE_SHA256
        or summary_record["size_bytes"] != 3_045
        or summary_record["sha256"] != EXPECTED_REJECTED_V29_ROUND1_SUMMARY_SHA256
        or EXPECTED_REJECTED_V29_ROUND1_EVIDENCE_SHA256.encode("ascii")
        not in summary_data
        or evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_pre_authority_candidate_rejection"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("candidate") != "v29_round1"
        or evidence.get("task_type") != "B_comprehensive_validation"
        or evidence.get("status")
        != "REJECTED_BEFORE_AUTHORITY_BY_STRICTEST_REPRODUCIBLE_REVIEW_WINS"
        or evidence.get("reviewed_source_set_sha256")
        != "0087b195ce2b7bfb0495f4f3e2c879b11484a645218a318402d14b47879a019e"
        or evidence.get("inventory_file_count") != 20
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v29 round1 evidence identity drift")

    inventory = evidence.get("inventory")
    if (
        not isinstance(inventory, Mapping)
        or set(inventory) != set(REJECTED_V29_ROUND1_PAYLOAD_RELATIVE_PATHS)
    ):
        raise RecoveryAuthorityError("Rejected v29 round1 inventory scope drift")
    expected_archive_paths = set(REJECTED_V29_ROUND1_PAYLOAD_RELATIVE_PATHS) | {
        "rejection_evidence.json",
        "SUMMARY.md",
    }
    observed_archive_paths = {
        path.relative_to(REJECTED_V29_ROUND1_ROOT).as_posix()
        for path in REJECTED_V29_ROUND1_ROOT.rglob("*")
        if path.is_file()
    }
    if observed_archive_paths != expected_archive_paths:
        raise RecoveryAuthorityError("Rejected v29 round1 archive file-set drift")

    payload_data: dict[str, bytes] = {}
    payload_records: dict[str, dict[str, Any]] = {}
    for relative_path, path in REJECTED_V29_ROUND1_PAYLOAD_PATHS.items():
        data, record, _ = _open_bound(path)
        expected_identity = inventory.get(relative_path)
        if (
            not isinstance(expected_identity, list)
            or len(expected_identity) != 2
            or type(expected_identity[0]) is not int
            or not isinstance(expected_identity[1], str)
            or record["size_bytes"] != expected_identity[0]
            or record["sha256"] != expected_identity[1]
        ):
            raise RecoveryAuthorityError(
                f"Rejected v29 round1 payload identity drift: {relative_path}"
            )
        payload_data[relative_path] = data
        payload_records[relative_path] = record

    mendel_initial = _load_json(
        payload_data["mendel_initial_approve.json"],
        "rejected v29 round1 Mendel initial review",
    )
    mendel_corrected = _load_json(
        payload_data["mendel_corrected_request_changes.json"],
        "rejected v29 round1 Mendel corrected review",
    )
    nash = _load_json(
        payload_data["nash_request_changes.json"],
        "rejected v29 round1 Nash review",
    )
    reproduction = _load_json(
        payload_data["orchestrator_ast_reproduction.json"],
        "rejected v29 round1 projection reproduction",
    )
    reviewed_digest = evidence["reviewed_source_set_sha256"]
    if (
        mendel_initial.get("reviewer_agent_id")
        != "019f9140-b2d4-7cc3-9f31-5b06645e66a5"
        or mendel_initial.get("verdict") != "APPROVE"
        or mendel_initial.get("reviewed_source_set_sha256") != reviewed_digest
        or mendel_initial.get("superseded") is not True
        or mendel_initial.get("authority_permission_revoked") is not True
        or mendel_initial.get("pass") is not False
        or mendel_corrected.get("reviewer_agent_id")
        != "019f9140-b2d4-7cc3-9f31-5b06645e66a5"
        or mendel_corrected.get("verdict") != "REQUEST_CHANGES"
        or mendel_corrected.get("reviewed_source_set_sha256") != reviewed_digest
        or len(mendel_corrected.get("high_findings", [])) != 1
        or len(mendel_corrected.get("medium_findings", [])) != 2
        or len(mendel_corrected.get("low_findings", [])) != 1
        or mendel_corrected.get("authority_allowed") is not False
        or mendel_corrected.get("pass") is not False
        or nash.get("reviewer_agent_id")
        != "019f9140-e2e8-76b2-9d63-cadbb2b5e964"
        or nash.get("verdict") != "REQUEST_CHANGES"
        or nash.get("reviewed_source_set_sha256") != reviewed_digest
        or len(nash.get("high_findings", [])) != 1
        or len(nash.get("medium_findings", [])) != 2
        or len(nash.get("low_findings", [])) != 1
        or nash.get("authority_allowed") is not False
        or nash.get("pass") is not False
    ):
        raise RecoveryAuthorityError("Rejected v29 round1 internal review drift")

    observed = reproduction.get("observed", {})
    lifecycle = reproduction.get("lifecycle_state_at_reproduction", {})
    if (
        reproduction.get("reviewed_source_set_sha256") != reviewed_digest
        or observed.get("validator_state_field_count") != 30
        or observed.get("continuation_expected_state_field_count") != 28
        or observed.get("validator_public_key_count") != 15
        or observed.get("continuation_expected_public_key_count") != 14
        or observed.get("validator_minus_continuation")
        != ["failed_v28_private_key", "failed_v28_public_key"]
        or observed.get("continuation_minus_validator") != []
        or observed.get("strict_json_equal_after_adding_one_validator_only_key")
        is not False
        or lifecycle.get("v29_authority_exists") is not False
        or lifecycle.get("v29_formal_reviews_exist") is not False
        or lifecycle.get("any_v29_key_consumed") is not False
        or reproduction.get("scientific_payload_changed") is not False
        or reproduction.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Rejected v29 round1 reproduction drift")

    envelope = json.loads(
        payload_data[
            "review_transport/20260724_external_claude_schema_recovery_review_v29."
            "claude_cli.envelope.json"
        ]
    )
    if not isinstance(envelope, Mapping):
        raise RecoveryAuthorityError("Rejected v29 round1 envelope is not an object")
    structured_output = envelope.get("structured_output")
    review_chain = evidence.get("review_chain", {})
    if (
        envelope.get("session_id") != "2b0e25af-34e4-4ae4-b171-b41ea89a37b8"
        or envelope.get("terminal_reason") != "completed"
        or not isinstance(structured_output, Mapping)
        or structured_output.get("reviewer_agent_id")
        != "2b0e25af-34e4-4ae4-b171-b41ea89a37b8"
        or structured_output.get("verdict") != "APPROVE"
        or structured_output.get("reviewed_source_set_sha256") != reviewed_digest
        or structured_output.get("pass") is not True
        or review_chain.get("strictest_review_wins") is not True
        or review_chain.get("external_claude_opus", {}).get(
            "formal_review_published"
        )
        is not False
        or review_chain.get("external_claude_opus", {}).get(
            "superseded_by_strictest_reproducible_finding"
        )
        is not True
    ):
        raise RecoveryAuthorityError("Rejected v29 round1 external review drift")

    blocking = evidence.get("blocking_reproduction", {})
    authority_state = evidence.get("authority_state", {})
    claim_ceiling = evidence.get("claim_ceiling", {})
    active_public_keys = evidence.get("active_public_key_sha256", {})
    if (
        blocking.get("validator_state_field_count") != 30
        or blocking.get("continuation_expected_state_field_count") != 28
        or blocking.get("validator_public_key_count") != 15
        or blocking.get("continuation_expected_public_key_count") != 14
        or blocking.get("validator_only_keys")
        != ["failed_v28_private_key", "failed_v28_public_key"]
        or blocking.get("strict_json_equal") is not False
        or len(evidence.get("required_corrections", [])) != 4
        or authority_state.get("authority_created") is not False
        or authority_state.get("authority_signature_created") is not False
        or authority_state.get("formal_reviews_published") is not False
        or authority_state.get("any_active_v29_key_consumed") is not False
        or authority_state.get("all_active_v29_private_keys_mode") != "0o400"
        or authority_state.get("all_active_v29_private_keys_link_count") != 1
        or authority_state.get("same_version_round2_key_reuse_allowed") is not True
        or active_public_keys
        != {
            "authority": "819c0ee9d6c6729ad41f3ea1e8071b90ddf506981c68b6927b8abc98e78e9cda",
            "terminal": "04d6acab01d56b0bfe25726a242904afd38bbc3ee47d4e3db29e9eb154e23e8b",
            "result": "84438d0a91200108ee06ad7600a3c5428804f37567c351ba33036843ae837864",
            "report": "79a684d855ee2d0010691c2a42439389d0e9148d0b84157e04b322c188df6c59",
        }
        or claim_ceiling.get("latent_molecular_substructure_candidates") is not True
        or claim_ceiling.get("confirmed_cellular_subclone_count") != 0
        or claim_ceiling.get("linear_ancestry_count") != 0
    ):
        raise RecoveryAuthorityError("Rejected v29 round1 lifecycle/claim drift")

    return {
        "evidence": evidence_record,
        "summary": summary_record,
        "payload_inventory": payload_records,
        "payload_inventory_count": len(payload_records),
        "reviewed_source_set_sha256": reviewed_digest,
        "mendel_initial_verdict": "APPROVE_SUPERSEDED",
        "mendel_corrected_verdict": "REQUEST_CHANGES",
        "nash_verdict": "REQUEST_CHANGES",
        "external_claude_verdict": "APPROVE_SUPERSEDED",
        "strictest_reproducible_review_wins": True,
        "authority_created_before_correction": False,
        "any_key_consumed_before_correction": False,
        "findings_corrected_in_active_candidate": True,
        "scientific_payload_changed": False,
        "claim_ceiling_changed": False,
        "pass": True,
    }


def _validate_rejected_generations() -> dict[str, Any]:
    return {
        "v2": _validate_rejected_generation(
            generation="v2",
            evidence_path=REJECTED_V2_EVIDENCE,
            expected_size=4_482,
            expected_sha256=EXPECTED_REJECTED_V2_EVIDENCE_SHA256,
            expected_source_set_sha256=(
                "90619369caa25359a6b5eda3ebd80e4705d50c2dbd4570ea68a3a4356c62fc74"
            ),
            private_key=REJECTED_V2_PRIVATE_KEY,
            output_slots=REJECTED_V2_OUTPUT_SLOTS,
            decisive_reviewer="nash",
            expected_high=["H1_builder_executes_second_unhashed_validator_read"],
            expected_medium=[
                "M1_preflight_inputs_not_rechecked_before_sign_and_publish"
            ],
        ),
        "v3": _validate_rejected_generation(
            generation="v3",
            evidence_path=REJECTED_V3_EVIDENCE,
            expected_size=4_473,
            expected_sha256=EXPECTED_REJECTED_V3_EVIDENCE_SHA256,
            expected_source_set_sha256=(
                "0533756fde39575a72504f0d6a93b01a52043a524895f0d019a1db36d0e6210a"
            ),
            private_key=REJECTED_V3_PRIVATE_KEY,
            output_slots=REJECTED_V3_OUTPUT_SLOTS,
            decisive_reviewer="nash",
            expected_high=["H1_post_sign_key_lifecycle_gap"],
            expected_medium=[
                "M1_single_open_contract_false",
                "M2_atomic_publish_TOCTOU",
            ],
        ),
        "v4": _validate_rejected_generation(
            generation="v4",
            evidence_path=REJECTED_V4_EVIDENCE,
            expected_size=4_676,
            expected_sha256=EXPECTED_REJECTED_V4_EVIDENCE_SHA256,
            expected_source_set_sha256=(
                "cd6a9f55f9086bdc94a6f3fad9b094680ea56427ae188a0bb4ac318804ae4e4d"
            ),
            private_key=REJECTED_V4_PRIVATE_KEY,
            output_slots=REJECTED_V4_OUTPUT_SLOTS,
            decisive_reviewer="mendel",
            expected_high=["H1_post_sign_BaseException_bypasses_key_retirement"],
            expected_medium=[
                "M1_bootstrap_validator_is_reopened_after_publish",
                "M2_staging_recheck_is_not_terminal",
            ],
        ),
        "v5": _validate_rejected_generation(
            generation="v5",
            evidence_path=REJECTED_V5_EVIDENCE,
            expected_size=4_322,
            expected_sha256=EXPECTED_REJECTED_V5_EVIDENCE_SHA256,
            expected_source_set_sha256=(
                "6cdbe386816fdfc22a6605d2a97110fb3cc4bcb2990c35c04188c2b840e02950"
            ),
            private_key=REJECTED_V5_PRIVATE_KEY,
            output_slots=REJECTED_V5_OUTPUT_SLOTS,
            decisive_reviewer="mendel",
            expected_high=[],
            expected_medium=["M1_missing_behavioral_post_publish_reopen_test"],
        ),
        "v6": _validate_rejected_generation(
            generation="v6",
            evidence_path=REJECTED_V6_EVIDENCE,
            expected_size=6_506,
            expected_sha256=EXPECTED_REJECTED_V6_EVIDENCE_SHA256,
            expected_source_set_sha256=(
                "87d68c542d9461d59d96a6d32b006467e4401f9a80e2c65ad878dbeffd632622"
            ),
            private_key=REJECTED_V6_PRIVATE_KEY,
            output_slots=REJECTED_V6_OUTPUT_SLOTS,
            decisive_reviewer="post_review_probe",
            expected_high=["H1_state_dependent_regression_breaks_post_review_probe"],
            expected_medium=[],
            expected_review_artifacts_created=True,
            expected_archived_reviews=REJECTED_V6_ARCHIVED_REVIEWS,
            expected_post_review_probe_failure=True,
        ),
        "v7": _validate_rejected_generation(
            generation="v7",
            evidence_path=REJECTED_V7_EVIDENCE,
            expected_size=5_732,
            expected_sha256=EXPECTED_REJECTED_V7_EVIDENCE_SHA256,
            expected_source_set_sha256=(
                "9eb4f27bd7ca496bbc489bbf188477fad79677aab95fa14618156d1d62ce0d24"
            ),
            private_key=REJECTED_V7_PRIVATE_KEY,
            output_slots=REJECTED_V7_OUTPUT_SLOTS,
            decisive_reviewer="mendel",
            expected_high=[],
            expected_medium=[
                "M1_forbidden_output_absence_not_terminally_rechecked"
            ],
        ),
        "v8": _validate_rejected_generation(
            generation="v8",
            evidence_path=REJECTED_V8_EVIDENCE,
            expected_size=5_695,
            expected_sha256=EXPECTED_REJECTED_V8_EVIDENCE_SHA256,
            expected_source_set_sha256=(
                "f391fd7687ea1cb5ae97931f03f64e93b5e53e490c6ebd911615e0b06aa01f90"
            ),
            private_key=REJECTED_V8_PRIVATE_KEY,
            output_slots=REJECTED_V8_OUTPUT_SLOTS,
            decisive_reviewer="nash",
            expected_high=[],
            expected_medium=[
                "M1_stage_identity_not_terminal_before_rename",
                "M2_forbidden_inventory_scan_has_intra_scan_TOCTOU",
            ],
        ),
        "v13": _validate_rejected_v13_pre_authority(),
        "v15": _validate_rejected_v15_pre_authority(),
        "v19_round1": _validate_rejected_v19_round1_pre_authority(),
        "v20": _validate_rejected_v20_pre_authority(),
        "v22": _validate_rejected_v22_pre_authority(),
        "v23": _validate_rejected_v23_pre_authority(),
        "v24": _validate_rejected_v24_pre_authority(),
        "v27": _validate_rejected_v27_pre_authority(),
        "v29_round1": _validate_rejected_v29_round1_pre_authority(),
    }


def _validate_failed_v9_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V9_EVIDENCE)
    if (
        evidence_record["size_bytes"] != 4_741
        or evidence_record["sha256"] != EXPECTED_FAILED_V9_EVIDENCE_SHA256
    ):
        raise RecoveryAuthorityError("Failed v9 evidence identity drift")
    evidence = _load_json(evidence_data, "failed signed v9 recovery evidence")
    _exact_keys(
        evidence,
        {
            "archived_artifacts",
            "authority_created",
            "authority_original_path",
            "authority_original_path_absent_after_archive",
            "created_at_utc",
            "diagnostic",
            "downstream_outputs_created",
            "generation",
            "private_key",
            "public_key",
            "recovery_verification_receipt_created",
            "resolution_contract",
            "runtime_error",
            "runtime_exit_code",
            "schema_name",
            "schema_version",
            "scientific_payload_changed",
            "source_set_sha256",
            "status",
            "pass",
            "pass_semantics",
        },
        "failed signed v9 recovery evidence",
    )

    archived_records = _records(FAILED_V9_ARCHIVED_ARTIFACT_PATHS)
    for role, record in archived_records.items():
        _exact_keys(record, {"path", "size_bytes", "sha256", "mode"}, f"v9 archive {role}")
    if not _strict_equal(evidence.get("archived_artifacts"), archived_records):
        raise RecoveryAuthorityError("Failed v9 archived-artifact records drift")

    public_data, public_record, public_fd = _open_bound(FAILED_V9_PUBLIC_KEY)
    retired_private, _ = _open_retired_private_key_bound(FAILED_V9_PRIVATE_KEY)
    expected_private_evidence = {
        **retired_private,
        "retired_after_single_signature": True,
    }
    if (
        len(public_data) != 113
        or public_record["sha256"]
        != "1aaa9f2fc71fa854f042fe8a8c6e99cc0a1d72d4224f4e85953e5d747c348b7b"
        or not _strict_equal(evidence.get("public_key"), public_record)
        or not _strict_equal(evidence.get("private_key"), expected_private_evidence)
    ):
        raise RecoveryAuthorityError("Failed v9 key identity/retirement drift")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V9_ARCHIVED_ARTIFACT_PATHS["authority"]
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V9_ARCHIVED_ARTIFACT_PATHS["signature"]
    )
    commit_data, commit_record, _ = _open_bound(
        FAILED_V9_ARCHIVED_ARTIFACT_PATHS["commit"]
    )
    if (
        len(signature_data) != 64
        or not _strict_equal(authority_record, archived_records["authority"])
        or not _strict_equal(signature_record, archived_records["signature"])
        or not _strict_equal(commit_record, archived_records["commit"])
    ):
        raise RecoveryAuthorityError("Failed v9 signed bundle identity drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    _validate_bundle_commit(
        commit_data,
        authority_record,
        signature_record,
        retired_private,
    )

    authority = _load_json(authority_data, "failed signed v9 authority")
    source_set_sha256 = _json_sha256(authority.get("recovery_sources"))
    if (
        authority.get("authority_id") != "20260722_tumor_ref_schema_recovery_v9"
        or authority.get("status")
        != "APPROVED_FOR_TRANSITION_RELATION_RECOVERY_ONLY"
        or authority.get("pass") is not True
        or not _strict_equal(authority.get("public_key"), public_record)
        or not _strict_equal(authority.get("retired_private_key"), retired_private)
        or source_set_sha256
        != "5d2c99f9cb84a055821e56650728822409d962a9097ac08c1eceea98599e9d8f"
        or evidence.get("source_set_sha256") != source_set_sha256
    ):
        raise RecoveryAuthorityError("Failed v9 authority payload drift")

    reviewer_ids: set[str] = set()
    authority_review_records = authority.get("review_artifacts")
    if not isinstance(authority_review_records, Mapping):
        raise RecoveryAuthorityError("Failed v9 authority review records missing")
    review_expected_keys = {
        "schema_name",
        "schema_version",
        "reviewer",
        "reviewer_agent_id",
        "review_completed_at_utc",
        "verdict",
        "reviewed_source_set_sha256",
        "legacy_source_set_sha256",
        "prior_recovery_chain_sha256",
        "rejected_generations_sha256",
        "trusted_recovery_public_key",
        "scope_sha256",
        "readonly_probe",
        "high_findings",
        "medium_findings",
        "low_findings",
        "unresolved_conditions",
        "summary",
        "pass",
    }
    for role, (reviewer, reviewer_id) in FAILED_V9_EXPECTED_REVIEWERS.items():
        archive_key = FAILED_V9_ARCHIVED_REVIEW_KEYS[role]
        archived_record = archived_records[archive_key]
        original_record = authority_review_records.get(role)
        if (
            not isinstance(original_record, Mapping)
            or original_record.get("path") != str(FAILED_V9_ORIGINAL_REVIEW_PATHS[role])
            or any(
                original_record.get(field) != archived_record[field]
                for field in ("size_bytes", "sha256", "mode")
            )
        ):
            raise RecoveryAuthorityError(f"Failed v9 review archive transition drift: {role}")
        review_data = _open_bound(FAILED_V9_ARCHIVED_ARTIFACT_PATHS[archive_key])[0]
        review = _load_json(review_data, f"failed v9 archived review {role}")
        _exact_keys(review, review_expected_keys, f"failed v9 archived review {role}")
        if (
            review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != SCHEMA_VERSION
            or review.get("reviewer") != reviewer
            or review.get("reviewer_agent_id") != reviewer_id
            or review.get("verdict") != "APPROVE"
            or review.get("reviewed_source_set_sha256") != source_set_sha256
            or review.get("legacy_source_set_sha256")
            != _json_sha256(authority.get("legacy_sources"))
            or review.get("prior_recovery_chain_sha256")
            != _json_sha256(authority.get("prior_recovery_chain"))
            or review.get("rejected_generations_sha256")
            != _json_sha256(authority.get("rejected_unsigned_generations"))
            or not _strict_equal(review.get("trusted_recovery_public_key"), public_record)
            or review.get("scope_sha256") != _json_sha256(authority.get("scope"))
            or review.get("readonly_probe", {}).get("forbidden_output_slots_checked")
            != 121
            or review.get("readonly_probe", {}).get("regression_summary")
            != "191 passed"
            or review.get("readonly_probe", {}).get("no_output_writes") is not True
            or not _strict_equal(review.get("high_findings"), [])
            or not _strict_equal(review.get("medium_findings"), [])
            or not _strict_equal(review.get("unresolved_conditions"), [])
            or review.get("pass") is not True
        ):
            raise RecoveryAuthorityError(f"Failed v9 archived review contract drift: {role}")
        reviewer_ids.add(reviewer_id)
    if len(reviewer_ids) != len(FAILED_V9_EXPECTED_REVIEWERS):
        raise RecoveryAuthorityError("Failed v9 reviewer identities are not distinct")

    expected_diagnostic = {
        "comparison": "signed_authorization.downstream_output_absence",
        "expected_current_inventory_length": 23,
        "first_changed_shared_index": None,
        "observed_signed_historical_length": 17,
        "root_cause": (
            "historical cross-link incorrectly compared against the expanded current "
            "recovery output inventory"
        ),
    }
    expected_runtime_error = {
        "message": "Signed authorization exact schema/cross-link drift",
        "type": "VerificationError",
    }
    occupied = [str(path) for path in FAILED_V9_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v9.bundle.staging.*"
        )
    )
    if (
        evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_failed_signed_generation"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("generation") != "v9"
        or evidence.get("status") != "SIGNED_AUTHORITY_RUNTIME_REJECTED_ARCHIVED"
        or evidence.get("authority_created") is not True
        or evidence.get("authority_original_path") != str(FAILED_V9_OUTPUT_SLOTS[0])
        or evidence.get("authority_original_path_absent_after_archive") is not True
        or evidence.get("downstream_outputs_created") is not False
        or evidence.get("recovery_verification_receipt_created") is not False
        or evidence.get("runtime_exit_code") != 1
        or type(evidence.get("runtime_exit_code")) is not int
        or not _strict_equal(evidence.get("diagnostic"), expected_diagnostic)
        or not _strict_equal(evidence.get("runtime_error"), expected_runtime_error)
        or evidence.get("resolution_contract")
        != (
            "v10 must compare signed historical cross-links against the exact 17-slot "
            "authorized inventory while independently enforcing absence of the full current "
            "23-slot recovery inventory"
        )
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("pass") is not True
        or evidence.get("pass_semantics")
        != (
            "records and cryptographically preserves the failed v9 authority generation; "
            "does not authorize reuse, continuation, or any scientific claim"
        )
        or occupied
        or staging
    ):
        raise RecoveryAuthorityError(
            f"Failed v9 runtime/archive contract drift: {occupied + staging}"
        )

    return {
        "evidence": evidence_record,
        "archived_artifacts": archived_records,
        "public_key": public_record,
        "retired_private_key": retired_private,
        "authority_signature_verified": True,
        "atomic_commit_verified": True,
        "authority_created": True,
        "runtime_outputs_created": False,
        "historical_authorized_output_slot_count": 17,
        "current_recovery_output_slot_count": 23,
        "source_set_sha256": source_set_sha256,
        "pass": True,
    }


def _validate_failed_v10_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V10_EVIDENCE)
    if (
        evidence_record["size_bytes"] != 6_190
        or evidence_record["sha256"] != EXPECTED_FAILED_V10_EVIDENCE_SHA256
    ):
        raise RecoveryAuthorityError("Failed v10 evidence identity drift")
    evidence = _load_json(evidence_data, "failed signed v10 recovery evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "generation",
            "status",
            "source_set_sha256",
            "runtime",
            "diagnostic",
            "archived_artifacts",
            "keys",
            "original_slots_absent_after_archive",
            "continuation_receipt_created",
            "continuation_exit_attestation_created",
            "continuation_signature_created",
            "continuation_incident_created",
            "canonical_downstream_outputs_created",
            "scientific_payload_changed",
            "claim_ceiling_changed",
            "required_resolution",
            "pass",
        },
        "failed signed v10 recovery evidence",
    )

    archived_records = _records(FAILED_V10_ARCHIVED_ARTIFACT_PATHS)
    for role, record in archived_records.items():
        _exact_keys(
            record,
            {"path", "size_bytes", "sha256", "mode"},
            f"v10 archive {role}",
        )
    if not _strict_equal(evidence.get("archived_artifacts"), archived_records):
        raise RecoveryAuthorityError("Failed v10 archived-artifact records drift")

    public_data, public_record, public_fd = _open_bound(FAILED_V10_PUBLIC_KEY)
    retired_private, _ = _open_retired_private_key_bound(FAILED_V10_PRIVATE_KEY)
    expected_key_evidence = {
        "public_key": public_record,
        "private_key": {
            "path": str(FAILED_V10_PRIVATE_KEY),
            "size_bytes": 119,
            "mode": "0o0",
            "retired_after_single_signature": True,
        },
    }
    if (
        len(public_data) != 113
        or public_record["sha256"]
        != "ac17fe53932a426216b12ed875b1ef84c946e0362f8a34359bff2d23bb51bc4e"
        or not _strict_equal(evidence.get("keys"), expected_key_evidence)
    ):
        raise RecoveryAuthorityError("Failed v10 key identity/retirement drift")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V10_ARCHIVED_ARTIFACT_PATHS["authority"]
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V10_ARCHIVED_ARTIFACT_PATHS["authority_signature"]
    )
    commit_data, commit_record, _ = _open_bound(
        FAILED_V10_ARCHIVED_ARTIFACT_PATHS["authority_commit"]
    )
    if len(signature_data) != 64:
        raise RecoveryAuthorityError("Failed v10 signature size drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    _validate_bundle_commit(
        commit_data,
        authority_record,
        signature_record,
        retired_private,
    )

    authority = _load_json(authority_data, "failed signed v10 authority")
    source_set_sha256 = _json_sha256(authority.get("recovery_sources"))
    live_v10_sources = _records(FAILED_V10_SOURCE_PATHS)
    if (
        authority.get("authority_id") != "20260723_tumor_ref_schema_recovery_v10"
        or authority.get("status")
        != "APPROVED_FOR_TRANSITION_RELATION_RECOVERY_ONLY"
        or authority.get("pass") is not True
        or source_set_sha256
        != "21cd1f719147a536abe91b3a54f419dbd50e9c0f794d4fd37583519e472bcd3a"
        or evidence.get("source_set_sha256") != source_set_sha256
        or not _strict_equal(authority.get("recovery_sources"), live_v10_sources)
        or not _strict_equal(authority.get("public_key"), public_record)
        or not _strict_equal(authority.get("retired_private_key"), retired_private)
    ):
        raise RecoveryAuthorityError("Failed v10 authority/source payload drift")

    review_expected_keys = {
        "schema_name",
        "schema_version",
        "reviewer",
        "reviewer_agent_id",
        "review_completed_at_utc",
        "verdict",
        "reviewed_source_set_sha256",
        "legacy_source_set_sha256",
        "prior_recovery_chain_sha256",
        "rejected_generations_sha256",
        "trusted_recovery_public_key",
        "scope_sha256",
        "readonly_probe",
        "high_findings",
        "medium_findings",
        "low_findings",
        "unresolved_conditions",
        "summary",
        "pass",
    }
    authority_reviews = authority.get("review_artifacts")
    if not isinstance(authority_reviews, Mapping):
        raise RecoveryAuthorityError("Failed v10 authority review records missing")
    reviewer_ids: set[str] = set()
    for role, (reviewer, reviewer_id) in FAILED_V10_EXPECTED_REVIEWERS.items():
        archive_key = FAILED_V10_ARCHIVED_REVIEW_KEYS[role]
        archived_record = archived_records[archive_key]
        original_record = authority_reviews.get(role)
        if (
            not isinstance(original_record, Mapping)
            or original_record.get("path")
            != str(FAILED_V10_ORIGINAL_REVIEW_PATHS[role])
            or any(
                original_record.get(field) != archived_record[field]
                for field in ("size_bytes", "sha256", "mode")
            )
        ):
            raise RecoveryAuthorityError(
                f"Failed v10 review archive transition drift: {role}"
            )
        review = _load_json(
            _open_bound(FAILED_V10_ARCHIVED_ARTIFACT_PATHS[archive_key])[0],
            f"failed v10 archived review {role}",
        )
        _exact_keys(review, review_expected_keys, f"failed v10 archived review {role}")
        if (
            review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != SCHEMA_VERSION
            or review.get("reviewer") != reviewer
            or review.get("reviewer_agent_id") != reviewer_id
            or review.get("verdict") != "APPROVE"
            or review.get("reviewed_source_set_sha256") != source_set_sha256
            or review.get("legacy_source_set_sha256")
            != _json_sha256(authority.get("legacy_sources"))
            or review.get("prior_recovery_chain_sha256")
            != _json_sha256(authority.get("prior_recovery_chain"))
            or review.get("rejected_generations_sha256")
            != _json_sha256(authority.get("rejected_unsigned_generations"))
            or not _strict_equal(review.get("trusted_recovery_public_key"), public_record)
            or review.get("scope_sha256") != _json_sha256(authority.get("scope"))
            or review.get("readonly_probe", {}).get("forbidden_output_slots_checked")
            != 134
            or review.get("readonly_probe", {}).get("regression_summary")
            != "213 passed"
            or review.get("readonly_probe", {}).get("no_output_writes") is not True
            or not _strict_equal(review.get("high_findings"), [])
            or not _strict_equal(review.get("medium_findings"), [])
            or not _strict_equal(review.get("unresolved_conditions"), [])
            or review.get("pass") is not True
        ):
            raise RecoveryAuthorityError(
                f"Failed v10 archived review contract drift: {role}"
            )
        reviewer_ids.add(reviewer_id)
    if len(reviewer_ids) != len(FAILED_V10_EXPECTED_REVIEWERS):
        raise RecoveryAuthorityError("Failed v10 reviewer identities are not distinct")

    verification = _load_json(
        _open_bound(FAILED_V10_ARCHIVED_ARTIFACT_PATHS["verification_receipt"])[0],
        "failed v10 verification receipt",
    )
    replay = _load_json(
        _open_bound(FAILED_V10_ARCHIVED_ARTIFACT_PATHS["replay_receipt"])[0],
        "failed v10 replay receipt",
    )
    verification_authority = verification.get("schema_recovery_authority")
    original_authority_path = str(
        RESULT_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v10.bundle"
        / "authority.json"
    )
    if (
        verification.get("schema_name")
        != "intersubmod.tumor_ref_source_receipt_promotion_verification.recovery"
        or verification.get("mode") != "verify-and-record"
        or verification.get("pass") is not True
        or not isinstance(verification_authority, Mapping)
        or verification_authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v10"
        or verification_authority.get("authority", {}).get("path")
        != original_authority_path
        or any(
            verification_authority.get("authority", {}).get(field)
            != archived_records["authority"][field]
            for field in ("size_bytes", "sha256", "mode")
        )
        or replay.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery"
        or replay.get("status")
        != "RUNNER_PREFIX_REPLAY_EVIDENCE_PASS_NON_AUTHORITATIVE"
        or replay.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Failed v10 V/R evidence contract drift")

    expected_runtime = {
        "v_exit_code": 0,
        "r_exit_code": 0,
        "c_exit_code": 1,
        "c_child_stderr_sha256": (
            "980809848fe8e8e0af55a7139b4710f8910bd8af90726055a29b52693b7c46e3"
        ),
        "c_error_type": "ContinuationError",
        "c_error_message": "Promotion verification payload contract drift",
    }
    expected_diagnostic = {
        "field": "historical_incident_disclosure.archive_receipt",
        "producer_value": (
            "authorization.evidence.failed_attempt_archive.receipt artifact record"
        ),
        "consumer_expected_value": (
            "authorization.evidence.failed_attempt_archive full object"
        ),
        "all_other_strict_verification_comparisons_passed": True,
        "failure_before_first_downstream_write": True,
    }
    occupied = [str(path) for path in FAILED_V10_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v10.bundle.staging.*"
        )
    )
    if (
        evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_signed_generation_failure"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("generation") != "v10"
        or evidence.get("status")
        != "SIGNED_AUTHORITY_V_AND_R_PASS_C_FAILED_BEFORE_DOWNSTREAM"
        or not _strict_equal(evidence.get("runtime"), expected_runtime)
        or not _strict_equal(evidence.get("diagnostic"), expected_diagnostic)
        or evidence.get("original_slots_absent_after_archive") is not True
        or evidence.get("continuation_receipt_created") is not False
        or evidence.get("continuation_exit_attestation_created") is not False
        or evidence.get("continuation_signature_created") is not False
        or evidence.get("continuation_incident_created") is not False
        or evidence.get("canonical_downstream_outputs_created") is not False
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("required_resolution")
        != (
            "Full append-only v11 with fresh source set, authority key, reviews, V/R/C "
            "slots, and producer-consumer parity regression."
        )
        or evidence.get("pass") is not False
        or occupied
        or staging
    ):
        raise RecoveryAuthorityError(
            f"Failed v10 runtime/archive contract drift: {occupied + staging}"
        )

    return {
        "evidence": evidence_record,
        "archived_artifacts": archived_records,
        "public_key": public_record,
        "retired_private_key": retired_private,
        "authority_signature_verified": True,
        "atomic_commit_verified": True,
        "authority_created": True,
        "verification_receipt_created": True,
        "replay_receipt_and_log_created": True,
        "continuation_outputs_created": False,
        "canonical_downstream_outputs_created": False,
        "source_set_sha256": source_set_sha256,
        "pass": True,
    }


def _validate_failed_v11_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V11_EVIDENCE)
    if (
        evidence_record["size_bytes"] != 9_812
        or evidence_record["sha256"] != EXPECTED_FAILED_V11_EVIDENCE_SHA256
    ):
        raise RecoveryAuthorityError("Failed v11 evidence identity drift")
    evidence = _load_json(evidence_data, "failed signed v11 recovery evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "generation",
            "status",
            "source_set_sha256",
            "runtime",
            "diagnostic",
            "archived_artifacts",
            "keys",
            "forbidden_output_slots",
            "original_slots_absent_after_archive",
            "continuation_receipt_created",
            "continuation_exit_attestation_created",
            "continuation_signature_created",
            "continuation_incident_created",
            "canonical_downstream_outputs_created",
            "scientific_payload_changed",
            "claim_ceiling_changed",
            "required_resolution",
            "pass",
        },
        "failed signed v11 recovery evidence",
    )

    archived_records = _records(FAILED_V11_ARCHIVED_ARTIFACT_PATHS)
    for role, record in archived_records.items():
        _exact_keys(
            record,
            {"path", "size_bytes", "sha256", "mode"},
            f"v11 archive {role}",
        )
    if not _strict_equal(evidence.get("archived_artifacts"), archived_records):
        raise RecoveryAuthorityError("Failed v11 archived-artifact records drift")

    public_data, public_record, public_fd = _open_bound(FAILED_V11_PUBLIC_KEY)
    retired_private, _ = _open_retired_private_key_bound(FAILED_V11_PRIVATE_KEY)
    continuation_public_data, continuation_public_record, _ = _open_bound(
        FAILED_V11_CONTINUATION_PUBLIC_KEY
    )
    continuation_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V11_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="UNUSED_NOT_RETIRED",
    )
    expected_key_evidence = {
        "public_key": public_record,
        "private_key": {
            "path": str(FAILED_V11_PRIVATE_KEY),
            "size_bytes": 119,
            "mode": "0o0",
            "retired_after_single_signature": True,
        },
        "continuation_key_state": {
            "public_key": continuation_public_record,
            "private_key": continuation_private_record,
        },
    }
    if (
        len(public_data) != 113
        or public_record["sha256"]
        != "2c1e3b609cb0f58dd23233d257ef66f3950ee5953dcda18474ce209525180b7b"
        or len(continuation_public_data) != 113
        or continuation_public_record["sha256"]
        != "a98370415594ebc9e8eb166bb70cab3550067b8a83db7d73ad23fac0891cb8b5"
        or not _strict_equal(evidence.get("keys"), expected_key_evidence)
    ):
        raise RecoveryAuthorityError("Failed v11 key identity/lifecycle drift")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V11_ARCHIVED_ARTIFACT_PATHS["authority"]
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V11_ARCHIVED_ARTIFACT_PATHS["authority_signature"]
    )
    commit_data, commit_record, _ = _open_bound(
        FAILED_V11_ARCHIVED_ARTIFACT_PATHS["authority_commit"]
    )
    if len(signature_data) != 64:
        raise RecoveryAuthorityError("Failed v11 signature size drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    _validate_bundle_commit(
        commit_data,
        authority_record,
        signature_record,
        retired_private,
    )

    authority = _load_json(authority_data, "failed signed v11 authority")
    source_set_sha256 = _json_sha256(authority.get("recovery_sources"))
    live_v11_sources = _records(FAILED_V11_SOURCE_PATHS)
    if (
        authority.get("authority_id") != "20260723_tumor_ref_schema_recovery_v11"
        or authority.get("status")
        != "APPROVED_FOR_TRANSITION_RELATION_RECOVERY_ONLY"
        or authority.get("pass") is not True
        or source_set_sha256
        != "b7e4963db70e0b1be3624818f152597e13febd3f0e592d85080fc7470b61487a"
        or evidence.get("source_set_sha256") != source_set_sha256
        or not _strict_equal(authority.get("recovery_sources"), live_v11_sources)
        or not _strict_equal(authority.get("public_key"), public_record)
        or not _strict_equal(authority.get("retired_private_key"), retired_private)
    ):
        raise RecoveryAuthorityError("Failed v11 authority/source payload drift")

    review_expected_keys = {
        "schema_name",
        "schema_version",
        "reviewer",
        "reviewer_agent_id",
        "review_completed_at_utc",
        "verdict",
        "reviewed_source_set_sha256",
        "legacy_source_set_sha256",
        "prior_recovery_chain_sha256",
        "rejected_generations_sha256",
        "trusted_recovery_public_key",
        "scope_sha256",
        "readonly_probe",
        "high_findings",
        "medium_findings",
        "low_findings",
        "unresolved_conditions",
        "summary",
        "pass",
    }
    authority_reviews = authority.get("review_artifacts")
    if not isinstance(authority_reviews, Mapping):
        raise RecoveryAuthorityError("Failed v11 authority review records missing")
    reviewer_ids: set[str] = set()
    for role, (reviewer, reviewer_id) in FAILED_V11_EXPECTED_REVIEWERS.items():
        archive_key = FAILED_V11_ARCHIVED_REVIEW_KEYS[role]
        archived_record = archived_records[archive_key]
        original_record = authority_reviews.get(role)
        if (
            not isinstance(original_record, Mapping)
            or original_record.get("path")
            != str(FAILED_V11_ORIGINAL_REVIEW_PATHS[role])
            or any(
                original_record.get(field) != archived_record[field]
                for field in ("size_bytes", "sha256", "mode")
            )
        ):
            raise RecoveryAuthorityError(
                f"Failed v11 review archive transition drift: {role}"
            )
        review = _load_json(
            _open_bound(FAILED_V11_ARCHIVED_ARTIFACT_PATHS[archive_key])[0],
            f"failed v11 archived review {role}",
        )
        _exact_keys(review, review_expected_keys, f"failed v11 archived review {role}")
        if (
            review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != SCHEMA_VERSION
            or review.get("reviewer") != reviewer
            or review.get("reviewer_agent_id") != reviewer_id
            or review.get("verdict") != "APPROVE"
            or review.get("reviewed_source_set_sha256") != source_set_sha256
            or review.get("legacy_source_set_sha256")
            != _json_sha256(authority.get("legacy_sources"))
            or review.get("prior_recovery_chain_sha256")
            != _json_sha256(authority.get("prior_recovery_chain"))
            or review.get("rejected_generations_sha256")
            != _json_sha256(authority.get("rejected_unsigned_generations"))
            or not _strict_equal(review.get("trusted_recovery_public_key"), public_record)
            or review.get("scope_sha256") != _json_sha256(authority.get("scope"))
            or review.get("readonly_probe", {}).get("forbidden_output_slots_checked")
            != 147
            or review.get("readonly_probe", {}).get("regression_summary")
            != "232 passed"
            or review.get("readonly_probe", {}).get("no_output_writes") is not True
            or not _strict_equal(review.get("high_findings"), [])
            or not _strict_equal(review.get("medium_findings"), [])
            or not _strict_equal(review.get("unresolved_conditions"), [])
            or review.get("pass") is not True
        ):
            raise RecoveryAuthorityError(
                f"Failed v11 archived review contract drift: {role}"
            )
        reviewer_ids.add(reviewer_id)
    if len(reviewer_ids) != len(FAILED_V11_EXPECTED_REVIEWERS):
        raise RecoveryAuthorityError("Failed v11 reviewer identities are not distinct")

    verification = _load_json(
        _open_bound(FAILED_V11_ARCHIVED_ARTIFACT_PATHS["verification_receipt"])[0],
        "failed v11 verification receipt",
    )
    replay = _load_json(
        _open_bound(FAILED_V11_ARCHIVED_ARTIFACT_PATHS["replay_receipt"])[0],
        "failed v11 replay receipt",
    )
    original_bundle = RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v11.bundle"

    def require_embedded_authority(
        payload: Mapping[str, Any], *, runtime_role: str
    ) -> None:
        embedded = payload.get("schema_recovery_authority")
        if not isinstance(embedded, Mapping):
            raise RecoveryAuthorityError(f"Failed v11 {runtime_role} authority missing")
        expected_originals = {
            "authority": original_bundle / "authority.json",
            "signature": original_bundle / "authority.ed25519.sig",
            "commit": original_bundle / "commit.json",
        }
        expected_archived = {
            "authority": archived_records["authority"],
            "signature": archived_records["authority_signature"],
            "commit": archived_records["authority_commit"],
        }
        for field, original_path in expected_originals.items():
            observed = embedded.get(field)
            archived = expected_archived[field]
            if (
                not isinstance(observed, Mapping)
                or observed.get("path") != str(original_path)
                or any(
                    observed.get(key) != archived[key]
                    for key in ("size_bytes", "sha256", "mode")
                )
            ):
                raise RecoveryAuthorityError(
                    f"Failed v11 {runtime_role} embedded {field} drift"
                )
        if (
            embedded.get("authority_id")
            != "20260723_tumor_ref_schema_recovery_v11"
            or embedded.get("authority_bundle") != str(original_bundle)
            or embedded.get("runtime_role") != runtime_role
            or not _strict_equal(
                embedded.get("runtime_source"), live_v11_sources[runtime_role]
            )
            or not _strict_equal(embedded.get("public_key"), public_record)
            or not _strict_equal(embedded.get("retired_private_key"), retired_private)
            or embedded.get("pass") is not True
        ):
            raise RecoveryAuthorityError(
                f"Failed v11 {runtime_role} embedded authority drift"
            )

    require_embedded_authority(verification, runtime_role="continuation_verifier")
    require_embedded_authority(replay, runtime_role="runner_gate_replay")
    v_source = verification.get("source_receipt")
    v_canonical = verification.get("canonical_receipt")
    trust_chain = replay.get("promotion_trust_chain")
    if not isinstance(trust_chain, Mapping):
        raise RecoveryAuthorityError("Failed v11 replay trust chain missing")
    r_source = trust_chain.get("historical_source_receipt")
    r_canonical = trust_chain.get("canonical_source_receipt")
    basic_keys = {"path", "size_bytes", "sha256", "mode"}
    full_keys = set(CURRENT_IDENTITY_KEYS)
    missing_stat_keys = {"device", "inode", "link_count", "mtime_ns", "ctime_ns"}
    if (
        verification.get("schema_name")
        != "intersubmod.tumor_ref_source_receipt_promotion_verification.recovery"
        or verification.get("mode") != "verify-and-record"
        or verification.get("pass") is not True
        or replay.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery"
        or replay.get("status")
        != "RUNNER_PREFIX_REPLAY_EVIDENCE_PASS_NON_AUTHORITATIVE"
        or replay.get("pass") is not True
        or not isinstance(v_source, Mapping)
        or not isinstance(v_canonical, Mapping)
        or not isinstance(r_source, Mapping)
        or not isinstance(r_canonical, Mapping)
        or set(v_source) != full_keys
        or set(v_canonical) != full_keys
        or set(r_source) != basic_keys
        or set(r_canonical) != basic_keys
        or full_keys - set(r_source) != missing_stat_keys
        or full_keys - set(r_canonical) != missing_stat_keys
        or not _strict_equal(r_source, {key: v_source[key] for key in basic_keys})
        or not _strict_equal(r_canonical, {key: v_canonical[key] for key in basic_keys})
    ):
        raise RecoveryAuthorityError("Failed v11 V/R projection mismatch evidence drift")

    expected_runtime = {
        "v_exit_code": 0,
        "r_exit_code": 0,
        "c_exit_code": 1,
        "c_supervised_child_returncode": 1,
        "c_child_stderr": {
            "encoding": "utf-8",
            "size_bytes": 230,
            "sha256": "f4d25c4f1857e3081e7ba85f9e82e1428b1d3599c678771bb1834f196a7b12a6",
        },
        "c_supervisor_stderr": {
            "encoding": "utf-8",
            "size_bytes": 312,
            "sha256": "bf4ba5073ad488ad9cfffb896e60fad369c4a769d169ab99861450b01b0df79c",
        },
        "c_error_type": "ContinuationError",
        "c_error_message": "Runner-only replay is not exact non-authoritative evidence",
    }
    expected_diagnostic = {
        "fields": [
            "promotion_trust_chain.historical_source_receipt",
            "promotion_trust_chain.canonical_source_receipt",
        ],
        "producer_projection_keys": sorted(basic_keys),
        "consumer_expected_projection_keys": sorted(full_keys),
        "missing_in_producer": sorted(missing_stat_keys),
        "strict_mismatch_count": 2,
        "all_previous_short_circuit_comparisons_passed": True,
        "failure_before_first_downstream_write": True,
    }
    occupied = [str(path) for path in FAILED_V11_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v11.bundle.staging.*"
        )
    )
    if (
        evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_signed_generation_failure"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("generation") != "v11"
        or evidence.get("status")
        != "SIGNED_AUTHORITY_V_AND_R_PASS_C_FAILED_BEFORE_DOWNSTREAM"
        or not _strict_equal(evidence.get("runtime"), expected_runtime)
        or not _strict_equal(evidence.get("diagnostic"), expected_diagnostic)
        or not _strict_equal(
            evidence.get("forbidden_output_slots"),
            [str(path) for path in FAILED_V11_OUTPUT_SLOTS],
        )
        or evidence.get("original_slots_absent_after_archive") is not True
        or evidence.get("continuation_receipt_created") is not False
        or evidence.get("continuation_exit_attestation_created") is not False
        or evidence.get("continuation_signature_created") is not False
        or evidence.get("continuation_incident_created") is not False
        or evidence.get("canonical_downstream_outputs_created") is not False
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("required_resolution")
        != (
            "Full append-only v12 with fresh source set, authority key, reviews, V/R/C "
            "slots, full-stat replay records, and producer-consumer parity regression."
        )
        or evidence.get("pass") is not False
        or occupied
        or staging
    ):
        raise RecoveryAuthorityError(
            f"Failed v11 runtime/archive contract drift: {occupied + staging}"
        )

    return {
        "evidence": evidence_record,
        "archived_artifacts": archived_records,
        "public_key": public_record,
        "retired_private_key": retired_private,
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": continuation_private_record,
        "authority_signature_verified": True,
        "atomic_commit_verified": True,
        "authority_created": True,
        "verification_receipt_created": True,
        "replay_receipt_and_log_created": True,
        "continuation_outputs_created": False,
        "canonical_downstream_outputs_created": False,
        "basic_to_full_stat_mismatch_reproduced": True,
        "source_set_sha256": source_set_sha256,
        "pass": True,
    }


def _validate_failed_v12_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V12_EVIDENCE)
    if (
        evidence_record["size_bytes"] != 11_521
        or evidence_record["sha256"] != EXPECTED_FAILED_V12_EVIDENCE_SHA256
    ):
        raise RecoveryAuthorityError("Failed v12 evidence identity drift")
    evidence = _load_json(evidence_data, "failed signed v12 recovery evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "generation",
            "status",
            "source_set_sha256",
            "scope_sha256",
            "runtime",
            "diagnostic",
            "archived_artifacts",
            "source_artifacts",
            "keys",
            "formal_output_state",
            "scientific_payload_changed",
            "claim_ceiling_changed",
            "required_resolution",
            "pass",
        },
        "failed signed v12 recovery evidence",
    )

    archived_records = _records(FAILED_V12_ARCHIVED_ARTIFACT_PATHS)
    source_records = _records(FAILED_V12_SOURCE_PATHS)
    if not _strict_equal(evidence.get("archived_artifacts"), archived_records):
        raise RecoveryAuthorityError("Failed v12 archived-artifact records drift")
    if not _strict_equal(evidence.get("source_artifacts"), source_records):
        raise RecoveryAuthorityError("Failed v12 source-artifact records drift")

    public_data, public_record, public_fd = _open_bound(FAILED_V12_PUBLIC_KEY)
    retired_private, _ = _open_retired_private_key_bound(FAILED_V12_PRIVATE_KEY)
    continuation_public_data, continuation_public_record, _ = _open_bound(
        FAILED_V12_CONTINUATION_PUBLIC_KEY
    )
    continuation_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V12_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="UNUSED_NOT_RETIRED",
    )
    expected_key_evidence = {
        "authority_public_key": public_record,
        "authority_private_key": {
            "path": str(FAILED_V12_PRIVATE_KEY),
            "size_bytes": 119,
            "mode": "0o0",
            "retired_after_single_signature": True,
        },
        "continuation_key_state": {
            "public_key": continuation_public_record,
            "private_key": continuation_private_record,
        },
    }
    if (
        len(public_data) != 113
        or public_record["sha256"]
        != "b6c0734543c9608f8af830c89bdde071a36b3cc38ab9d5c53399fd63a2d46d3d"
        or len(continuation_public_data) != 113
        or continuation_public_record["sha256"]
        != "a98370415594ebc9e8eb166bb70cab3550067b8a83db7d73ad23fac0891cb8b5"
        or not _strict_equal(evidence.get("keys"), expected_key_evidence)
    ):
        raise RecoveryAuthorityError("Failed v12 key identity/lifecycle drift")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V12_ARCHIVED_ARTIFACT_PATHS["authority"]
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V12_ARCHIVED_ARTIFACT_PATHS["authority_signature"]
    )
    commit_data, commit_record, _ = _open_bound(
        FAILED_V12_ARCHIVED_ARTIFACT_PATHS["authority_commit"]
    )
    if len(signature_data) != 64:
        raise RecoveryAuthorityError("Failed v12 signature size drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    _validate_bundle_commit(
        commit_data,
        authority_record,
        signature_record,
        retired_private,
    )

    authority = _load_json(authority_data, "failed signed v12 authority")
    source_set_sha256 = _json_sha256(authority.get("recovery_sources"))
    scope_sha256 = _json_sha256(authority.get("scope"))
    if (
        authority.get("authority_id") != "20260723_tumor_ref_schema_recovery_v12"
        or authority.get("status")
        != "APPROVED_FOR_TRANSITION_RELATION_RECOVERY_ONLY"
        or authority.get("pass") is not True
        or source_set_sha256
        != "e738d06ad0730854a0d72ac29ad9f7dd1d45abba1d4c5ae7aeda3465d791ce93"
        or scope_sha256
        != "56d03544995cd01919e94761816e76b141ff6714e8d51fe2b86e80af71770f5a"
        or evidence.get("source_set_sha256") != source_set_sha256
        or evidence.get("scope_sha256") != scope_sha256
        or not _strict_equal(authority.get("recovery_sources"), source_records)
        or not _strict_equal(authority.get("public_key"), public_record)
        or not _strict_equal(authority.get("retired_private_key"), retired_private)
    ):
        raise RecoveryAuthorityError("Failed v12 authority/source payload drift")

    review_expected_keys = {
        "schema_name",
        "schema_version",
        "reviewer",
        "reviewer_agent_id",
        "review_completed_at_utc",
        "verdict",
        "reviewed_source_set_sha256",
        "legacy_source_set_sha256",
        "prior_recovery_chain_sha256",
        "rejected_generations_sha256",
        "trusted_recovery_public_key",
        "scope_sha256",
        "readonly_probe",
        "high_findings",
        "medium_findings",
        "low_findings",
        "unresolved_conditions",
        "summary",
        "pass",
    }
    authority_reviews = authority.get("review_artifacts")
    if not isinstance(authority_reviews, Mapping):
        raise RecoveryAuthorityError("Failed v12 authority review records missing")
    reviewer_ids: set[str] = set()
    for role, (reviewer, reviewer_id) in FAILED_V12_EXPECTED_REVIEWERS.items():
        archive_key = FAILED_V12_ARCHIVED_REVIEW_KEYS[role]
        archived_record = archived_records[archive_key]
        original_record = authority_reviews.get(role)
        if (
            not isinstance(original_record, Mapping)
            or original_record.get("path")
            != str(FAILED_V12_ORIGINAL_REVIEW_PATHS[role])
            or any(
                original_record.get(field) != archived_record[field]
                for field in ("size_bytes", "sha256", "mode")
            )
        ):
            raise RecoveryAuthorityError(
                f"Failed v12 review archive transition drift: {role}"
            )
        review = _load_json(
            _open_bound(FAILED_V12_ARCHIVED_ARTIFACT_PATHS[archive_key])[0],
            f"failed v12 archived review {role}",
        )
        _exact_keys(review, review_expected_keys, f"failed v12 archived review {role}")
        probe = review.get("readonly_probe")
        if (
            review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != SCHEMA_VERSION
            or review.get("reviewer") != reviewer
            or review.get("reviewer_agent_id") != reviewer_id
            or review.get("verdict") != "APPROVE"
            or review.get("reviewed_source_set_sha256") != source_set_sha256
            or review.get("legacy_source_set_sha256")
            != _json_sha256(authority.get("legacy_sources"))
            or review.get("prior_recovery_chain_sha256")
            != _json_sha256(authority.get("prior_recovery_chain"))
            or review.get("rejected_generations_sha256")
            != _json_sha256(authority.get("rejected_unsigned_generations"))
            or not _strict_equal(review.get("trusted_recovery_public_key"), public_record)
            or review.get("scope_sha256") != scope_sha256
            or not isinstance(probe, Mapping)
            or probe.get("forbidden_output_slots_checked") != 161
            or probe.get("regression_summary") != "255 passed"
            or probe.get("no_output_writes") is not True
            or not _strict_equal(review.get("high_findings"), [])
            or not _strict_equal(review.get("medium_findings"), [])
            or not _strict_equal(review.get("unresolved_conditions"), [])
            or review.get("pass") is not True
        ):
            raise RecoveryAuthorityError(
                f"Failed v12 archived review contract drift: {role}"
            )
        reviewer_ids.add(reviewer_id)
    if len(reviewer_ids) != len(FAILED_V12_EXPECTED_REVIEWERS):
        raise RecoveryAuthorityError("Failed v12 reviewer identities are not distinct")

    verification = _load_json(
        _open_bound(FAILED_V12_ARCHIVED_ARTIFACT_PATHS["verification_receipt"])[0],
        "failed v12 verification receipt",
    )
    replay = _load_json(
        _open_bound(FAILED_V12_ARCHIVED_ARTIFACT_PATHS["replay_receipt"])[0],
        "failed v12 replay receipt",
    )
    original_bundle = RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v12.bundle"

    def require_embedded_authority(
        payload: Mapping[str, Any], *, runtime_role: str
    ) -> None:
        embedded = payload.get("schema_recovery_authority")
        if not isinstance(embedded, Mapping):
            raise RecoveryAuthorityError(f"Failed v12 {runtime_role} authority missing")
        for field, archived_key in (
            ("authority", "authority"),
            ("signature", "authority_signature"),
            ("commit", "authority_commit"),
        ):
            observed = embedded.get(field)
            archived = archived_records[archived_key]
            original_name = {
                "authority": "authority.json",
                "signature": "authority.ed25519.sig",
                "commit": "commit.json",
            }[field]
            if (
                not isinstance(observed, Mapping)
                or observed.get("path") != str(original_bundle / original_name)
                or any(
                    observed.get(key) != archived[key]
                    for key in ("size_bytes", "sha256", "mode")
                )
            ):
                raise RecoveryAuthorityError(
                    f"Failed v12 {runtime_role} embedded {field} drift"
                )
        if (
            embedded.get("authority_id")
            != "20260723_tumor_ref_schema_recovery_v12"
            or embedded.get("authority_bundle") != str(original_bundle)
            or embedded.get("runtime_role") != runtime_role
            or not _strict_equal(
                embedded.get("runtime_source"), source_records[runtime_role]
            )
            or not _strict_equal(embedded.get("public_key"), public_record)
            or not _strict_equal(embedded.get("retired_private_key"), retired_private)
            or embedded.get("pass") is not True
        ):
            raise RecoveryAuthorityError(
                f"Failed v12 {runtime_role} embedded authority drift"
            )

    require_embedded_authority(verification, runtime_role="continuation_verifier")
    require_embedded_authority(replay, runtime_role="runner_gate_replay")
    declared_alias = evidence["diagnostic"]["declared_alias_record"]
    alias_paths = (
        ("completion_runner", "authorization_bound_runtime_contract", "runtime_artifacts", "python"),
        ("promotion_verifier", "fd_artifact_binding", "python", "artifact"),
        ("runtime", "python"),
    )
    for json_path in alias_paths:
        node: Any = replay
        for field in json_path:
            if not isinstance(node, Mapping):
                raise RecoveryAuthorityError("Failed v12 replay alias path is not a mapping")
            node = node.get(field)
        if not _strict_equal(node, declared_alias):
            raise RecoveryAuthorityError("Failed v12 replay alias record drift")

    alias_path = Path(declared_alias["path"])
    resolved_path = Path(declared_alias["resolved_path"])
    alias_stat = os.lstat(alias_path)
    target_data, target_record, _ = _open_bound(
        resolved_path, expected_mode=declared_alias["mode"]
    )
    expected_live_alias = {
        "path": str(alias_path),
        "device": alias_stat.st_dev,
        "inode": alias_stat.st_ino,
        "link_count": alias_stat.st_nlink,
        "size_bytes": alias_stat.st_size,
        "mode": oct(stat.S_IMODE(alias_stat.st_mode)),
        "file_type": "symbolic_link",
        "symlink_target": os.readlink(alias_path),
    }
    expected_live_target = {
        "path": str(resolved_path),
        "device": os.stat(resolved_path, follow_symlinks=False).st_dev,
        "inode": os.stat(resolved_path, follow_symlinks=False).st_ino,
        "link_count": os.stat(resolved_path, follow_symlinks=False).st_nlink,
        "size_bytes": len(target_data),
        "mode": target_record["mode"],
        "file_type": "regular_file",
        "sha256": target_record["sha256"],
    }
    if (
        not stat.S_ISLNK(alias_stat.st_mode)
        or alias_path.resolve(strict=True) != resolved_path
        or not _strict_equal(evidence["diagnostic"]["live_alias"], expected_live_alias)
        or not _strict_equal(
            evidence["diagnostic"]["live_resolved_target"], expected_live_target
        )
    ):
        raise RecoveryAuthorityError("Failed v12 live executable alias evidence drift")

    expected_runtime = {
        "v_exit_code": 0,
        "r_exit_code": 0,
        "c_exit_code": 1,
        "c_supervised_child_returncode": 1,
        "c_child_stderr": {
            "encoding": "utf-8",
            "size_bytes": 265,
            "sha256": "1075f008453f30fb03e64b730fe16e4a9421ab669de7ffe18ef0c86e44fc3d34",
        },
        "c_supervisor_stderr": {
            "encoding": "utf-8",
            "size_bytes": 312,
            "sha256": "37532fa119f3e183006f8b6e68df2fc6638d9618e451d97d1e1426b52c606006",
        },
        "c_error_type": "ContinuationError",
        "c_error_message": (
            "Declared artifact relation drift: "
            "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python"
        ),
        "diagnostic_provenance": (
            "The exact child error JSON was recovered from the formal child write syscall; "
            "its reconstructed bytes match the stderr SHA-256 reported by the formal supervisor."
        ),
    }
    expected_output_state = {
        "original_v12_authority_review_v_r_slots_absent_after_archive": True,
        "continuation_receipt_created": False,
        "continuation_exit_attestation_created": False,
        "continuation_signature_created": False,
        "continuation_success_witness_created": False,
        "continuation_incident_created": False,
        "canonical_downstream_outputs_created": False,
        "final_dataset_release_receipt_created": False,
        "final_report_release_receipt_created": False,
    }
    occupied = [str(path) for path in FAILED_V12_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v12.bundle.staging.*"
        )
    )
    if (
        evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_signed_generation_failure"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("generation") != "v12"
        or evidence.get("status")
        != "SIGNED_AUTHORITY_V_AND_R_PASS_C_FAILED_BEFORE_DOWNSTREAM"
        or not _strict_equal(evidence.get("runtime"), expected_runtime)
        or evidence["diagnostic"].get("failure_before_first_downstream_write")
        is not True
        or evidence["diagnostic"].get("scientific_payload_changed") is not False
        or not _strict_equal(evidence.get("formal_output_state"), expected_output_state)
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("required_resolution")
        != (
            "Create append-only v13 with a fresh authority key and fresh continuation key, "
            "explicit executable-alias relation binding, process-lifetime alias and target "
            "leases, adversarial symlink regression tests, fresh independent reviews, and "
            "new V/R/C output slots."
        )
        or evidence.get("pass") is not False
        or occupied
        or staging
    ):
        raise RecoveryAuthorityError(
            f"Failed v12 runtime/archive contract drift: {occupied + staging}"
        )

    return {
        "evidence": evidence_record,
        "archived_artifacts": archived_records,
        "public_key": public_record,
        "retired_private_key": retired_private,
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": continuation_private_record,
        "authority_signature_verified": True,
        "atomic_commit_verified": True,
        "authority_created": True,
        "verification_receipt_created": True,
        "replay_receipt_and_log_created": True,
        "continuation_outputs_created": False,
        "canonical_downstream_outputs_created": False,
        "executable_alias_relation_mismatch_reproduced": True,
        "source_set_sha256": source_set_sha256,
        "pass": True,
    }


def _validate_failed_v14_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V14_EVIDENCE)
    if (
        evidence_record["size_bytes"] != 12_657
        or evidence_record["sha256"] != EXPECTED_FAILED_V14_EVIDENCE_SHA256
    ):
        raise RecoveryAuthorityError("Failed v14 evidence identity drift")
    evidence = _load_json(evidence_data, "failed signed v14 recovery evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "generation",
            "status",
            "review_contract",
            "runtime",
            "diagnostic",
            "archived_artifacts",
            "source_artifacts",
            "reviews",
            "keys",
            "formal_output_state",
            "scientific_payload_changed",
            "claim_ceiling_changed",
            "required_resolution",
            "pass",
        },
        "failed signed v14 recovery evidence",
    )

    archived_records = _records(FAILED_V14_ARCHIVED_ARTIFACT_PATHS)
    source_records = _records(FAILED_V14_SOURCE_PATHS)
    if not _strict_equal(evidence.get("archived_artifacts"), archived_records):
        raise RecoveryAuthorityError("Failed v14 archived-artifact records drift")
    if not _strict_equal(evidence.get("source_artifacts"), source_records):
        raise RecoveryAuthorityError("Failed v14 source-artifact records drift")
    if any(record.get("mode") != "0o444" for record in source_records.values()):
        raise RecoveryAuthorityError("Failed v14 source set is not immutable mode 0444")

    public_data, public_record, public_fd = _open_bound(FAILED_V14_PUBLIC_KEY)
    retired_private, _ = _open_retired_private_key_bound(FAILED_V14_PRIVATE_KEY)
    continuation_public_data, continuation_public_record, _ = _open_bound(
        FAILED_V14_CONTINUATION_PUBLIC_KEY
    )
    continuation_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V14_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="UNUSED_NOT_RETIRED_NO_SIGNATURE_CREATED",
    )
    expected_keys = {
        "authority_public_key": public_record,
        "authority_private_key": {
            "path": str(FAILED_V14_PRIVATE_KEY),
            "size_bytes": 119,
            "mode": "0o0",
            "retired_after_single_signature": True,
        },
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": continuation_private_record,
    }
    if (
        len(public_data) != 113
        or public_record["sha256"]
        != "91f3b81a0dfab1911b492269dc40ef150eed76bf3b16c4143ba541d16ffdc8a3"
        or len(continuation_public_data) != 113
        or continuation_public_record["sha256"]
        != "091785646a7cadf2295f97668838b8279a2bb3b8a798e317dcf1ec6aba33d427"
        or not _strict_equal(evidence.get("keys"), expected_keys)
    ):
        raise RecoveryAuthorityError("Failed v14 key identity/lifecycle drift")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V14_ARCHIVED_ARTIFACT_PATHS["authority"]
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V14_ARCHIVED_ARTIFACT_PATHS["authority_signature"]
    )
    commit_data, commit_record, _ = _open_bound(
        FAILED_V14_ARCHIVED_ARTIFACT_PATHS["authority_commit"]
    )
    if len(signature_data) != 64:
        raise RecoveryAuthorityError("Failed v14 signature size drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    _validate_bundle_commit(commit_data, authority_record, signature_record, retired_private)

    authority = _load_json(authority_data, "failed signed v14 authority")
    authority_sources = authority.get("recovery_sources")
    expected_authority_sources = {
        role: record for role, record in source_records.items() if role != "archive_script"
    }
    review_contract = evidence.get("review_contract")
    if not isinstance(review_contract, Mapping):
        raise RecoveryAuthorityError("Failed v14 review contract missing")
    source_set_sha256 = _json_sha256(authority_sources)
    scope_sha256 = _json_sha256(authority.get("scope"))
    if (
        authority.get("authority_id") != "20260723_tumor_ref_schema_recovery_v14"
        or authority.get("status")
        != "APPROVED_FOR_TRANSITION_ALIAS_AND_TERMINAL_KEY_RECOVERY_ONLY"
        or authority.get("pass") is not True
        or not _strict_equal(authority_sources, expected_authority_sources)
        or not _strict_equal(authority.get("public_key"), public_record)
        or not _strict_equal(authority.get("retired_private_key"), retired_private)
        or review_contract.get("reviewed_source_set_sha256") != source_set_sha256
        or review_contract.get("scope_sha256") != scope_sha256
        or review_contract.get("legacy_source_set_sha256")
        != _json_sha256(authority.get("legacy_sources"))
        or review_contract.get("prior_recovery_chain_sha256")
        != _json_sha256(authority.get("prior_recovery_chain"))
        or review_contract.get("rejected_generations_sha256")
        != _json_sha256(authority.get("rejected_unsigned_generations"))
        or review_contract.get("terminal_key_rotation_sha256")
        != _json_sha256(authority.get("terminal_key_rotation"))
    ):
        raise RecoveryAuthorityError("Failed v14 authority/source contract drift")

    authority_reviews = authority.get("review_artifacts")
    evidence_reviews = evidence.get("reviews")
    if not isinstance(authority_reviews, Mapping) or not isinstance(evidence_reviews, Mapping):
        raise RecoveryAuthorityError("Failed v14 review records missing")
    reviewer_ids: set[str] = set()
    for role, (reviewer, reviewer_id) in FAILED_V14_EXPECTED_REVIEWERS.items():
        archive_key = FAILED_V14_ARCHIVED_REVIEW_KEYS[role]
        archived_record = archived_records[archive_key]
        original_record = authority_reviews.get(role)
        if (
            not isinstance(original_record, Mapping)
            or original_record.get("path")
            != str(FAILED_V14_ORIGINAL_REVIEW_PATHS[role])
            or any(
                original_record.get(field) != archived_record[field]
                for field in ("size_bytes", "sha256", "mode")
            )
        ):
            raise RecoveryAuthorityError(
                f"Failed v14 review archive transition drift: {role}"
            )
        review = _load_json(
            _open_bound(FAILED_V14_ARCHIVED_ARTIFACT_PATHS[archive_key])[0],
            f"failed v14 archived review {role}",
        )
        probe = review.get("readonly_probe")
        expected_review_evidence = {
            "reviewer": reviewer,
            "reviewer_agent_id": reviewer_id,
            "verdict": "APPROVE",
            "findings_empty": True,
        }
        if (
            review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != SCHEMA_VERSION
            or review.get("reviewer") != reviewer
            or review.get("reviewer_agent_id") != reviewer_id
            or review.get("verdict") != "APPROVE"
            or review.get("reviewed_source_set_sha256") != source_set_sha256
            or review.get("scope_sha256") != scope_sha256
            or review.get("terminal_key_rotation_sha256")
            != review_contract.get("terminal_key_rotation_sha256")
            or not isinstance(probe, Mapping)
            or probe.get("forbidden_output_slots_checked") != 191
            or probe.get("regression_summary") != "329 passed"
            or probe.get("no_output_writes") is not True
            or any(review.get(key) for key in ("high_findings", "medium_findings", "low_findings"))
            or review.get("unresolved_conditions") != []
            or review.get("pass") is not True
            or not _strict_equal(evidence_reviews.get(role), expected_review_evidence)
        ):
            raise RecoveryAuthorityError(f"Failed v14 review contract drift: {role}")
        reviewer_ids.add(reviewer_id)
    if len(reviewer_ids) != len(FAILED_V14_EXPECTED_REVIEWERS):
        raise RecoveryAuthorityError("Failed v14 reviewer identities are not distinct")

    verification = _load_json(
        _open_bound(FAILED_V14_ARCHIVED_ARTIFACT_PATHS["verification_receipt"])[0],
        "failed v14 verification receipt",
    )
    embedded_authority = verification.get("schema_recovery_authority")
    if (
        verification.get("schema_name")
        != "intersubmod.tumor_ref_source_receipt_promotion_verification.recovery"
        or verification.get("mode") != "verify-and-record"
        or verification.get("pass") is not True
        or not isinstance(embedded_authority, Mapping)
        or embedded_authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v14"
        or embedded_authority.get("runtime_role") != "continuation_verifier"
        or embedded_authority.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Failed v14 V receipt contract drift")

    stderr_data = _open_bound(
        FAILED_V14_ARCHIVED_ARTIFACT_PATHS["formal_replay_stderr"]
    )[0]
    stderr_text = stderr_data.decode("utf-8")
    expected_error = (
        "Metadata-only artifact cannot be reopened for bytes: "
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
        "20260722_m2v5_terminal_v2/ed25519_private_one_time_resident.pem"
    )
    runtime = evidence.get("runtime")
    diagnostic = evidence.get("diagnostic")
    output_state = evidence.get("formal_output_state")
    expected_output_state = {
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
    }
    occupied = [str(path) for path in FAILED_V14_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v14.bundle.staging.*"
        )
    )
    if (
        evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_signed_generation_failure"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("generation") != "v14"
        or evidence.get("status")
        != "SIGNED_AUTHORITY_AND_V_PASS_R_FAILED_BEFORE_FIRST_R_OUTPUT"
        or not isinstance(runtime, Mapping)
        or runtime.get("v_exit_code") != 0
        or runtime.get("r_supervisor_exit_code") != 1
        or runtime.get("r_worker_exit_code") != 70
        or runtime.get("r_waitpid_raw_status") != 17920
        or runtime.get("r_worker_normal_exit") is not True
        or runtime.get("r_worker_success") is not False
        or runtime.get("c_started") is not False
        or not _strict_equal(runtime.get("stderr"), archived_records["formal_replay_stderr"])
        or runtime.get("error_type") != "GateError"
        or runtime.get("error_message") != expected_error
        or expected_error not in stderr_text
        or "status=17920" not in stderr_text
        or not isinstance(diagnostic, Mapping)
        or diagnostic.get("record_shape") != ["mode", "path", "size_bytes", "state"]
        or diagnostic.get("private_key_bytes_read") is not False
        or diagnostic.get("failure_before_first_r_output") is not True
        or diagnostic.get("failure_before_c_start") is not True
        or diagnostic.get("scientific_payload_changed") is not False
        or not _strict_equal(output_state, expected_output_state)
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("pass") is not False
        or occupied
        or staging
    ):
        raise RecoveryAuthorityError(
            f"Failed v14 runtime/archive contract drift: {occupied + staging}"
        )

    return {
        "evidence": evidence_record,
        "archived_artifacts": archived_records,
        "public_key": public_record,
        "retired_private_key": retired_private,
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": continuation_private_record,
        "authority_signature_verified": True,
        "atomic_commit_verified": True,
        "authority_created": True,
        "verification_receipt_created": True,
        "replay_outputs_created": False,
        "continuation_outputs_created": False,
        "canonical_downstream_outputs_created": False,
        "metadata_only_relation_schema_mismatch_reproduced": True,
        "source_set_sha256": source_set_sha256,
        "pass": True,
    }


def _validate_failed_v16_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V16_EVIDENCE)
    if (
        evidence_record["size_bytes"] != 17_151
        or evidence_record["sha256"] != EXPECTED_FAILED_V16_EVIDENCE_SHA256
    ):
        raise RecoveryAuthorityError("Failed v16 evidence identity drift")
    evidence = _load_json(evidence_data, "failed signed v16 recovery evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "generation",
            "status",
            "review_contract",
            "runtime",
            "diagnostic",
            "archived_artifacts",
            "source_artifacts",
            "review_transport_artifacts",
            "reviews",
            "keys",
            "formal_output_state",
            "scientific_payload_changed",
            "claim_ceiling_changed",
            "required_resolution",
            "pass",
        },
        "failed signed v16 recovery evidence",
    )

    archived_records = _records(FAILED_V16_ARCHIVED_ARTIFACT_PATHS)
    source_records = _records(FAILED_V16_SOURCE_PATHS)
    transport_records = _records(FAILED_V16_REVIEW_TRANSPORT_PATHS)
    if not _strict_equal(evidence.get("archived_artifacts"), archived_records):
        raise RecoveryAuthorityError("Failed v16 archived-artifact records drift")
    if not _strict_equal(evidence.get("source_artifacts"), source_records):
        raise RecoveryAuthorityError("Failed v16 source-artifact records drift")
    if not _strict_equal(
        evidence.get("review_transport_artifacts"), transport_records
    ):
        raise RecoveryAuthorityError("Failed v16 review-transport records drift")
    if any(record.get("mode") != "0o444" for record in source_records.values()):
        raise RecoveryAuthorityError("Failed v16 source set is not immutable mode 0444")
    if any(record.get("mode") != "0o444" for record in transport_records.values()):
        raise RecoveryAuthorityError(
            "Failed v16 review transport is not immutable mode 0444"
        )

    public_data, public_record, public_fd = _open_bound(FAILED_V16_PUBLIC_KEY)
    retired_private, retired_private_fd = _open_retired_private_key_bound(
        FAILED_V16_PRIVATE_KEY
    )
    retired_private_stat = os.fstat(retired_private_fd)
    continuation_public_data, continuation_public_record, _ = _open_bound(
        FAILED_V16_CONTINUATION_PUBLIC_KEY
    )
    continuation_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V16_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="UNUSED_NOT_RETIRED_NO_SIGNATURE_CREATED_MUST_NOT_BE_REUSED",
    )
    expected_keys = {
        "authority_public_key": public_record,
        "authority_private_key": {
            "path": str(FAILED_V16_PRIVATE_KEY),
            "size_bytes": 119,
            "mode": "0o0",
            "link_count": 1,
            "state": "RETIRED_AFTER_SINGLE_AUTHORITY_SIGNATURE",
        },
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": {
            **continuation_private_record,
            "link_count": 1,
        },
    }
    if (
        len(public_data) != 113
        or public_record["sha256"]
        != "540b64ed3615618efed89637069f772787fdd025acadbbd27b6e334423d2345e"
        or retired_private_stat.st_size != 119
        or retired_private_stat.st_nlink != 1
        or len(continuation_public_data) != 113
        or continuation_public_record["sha256"]
        != "066949c0c36be413cd2fb60670e5a2fbc583ab9a3a4264ebf4d3766aba39867f"
        or not _strict_equal(evidence.get("keys"), expected_keys)
    ):
        raise RecoveryAuthorityError("Failed v16 key identity/lifecycle drift")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V16_ARCHIVED_ARTIFACT_PATHS["authority"]
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V16_ARCHIVED_ARTIFACT_PATHS["authority_signature"]
    )
    commit_data, commit_record, _ = _open_bound(
        FAILED_V16_ARCHIVED_ARTIFACT_PATHS["authority_commit"]
    )
    if len(signature_data) != 64:
        raise RecoveryAuthorityError("Failed v16 signature size drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    _validate_bundle_commit(
        commit_data,
        authority_record,
        signature_record,
        retired_private,
    )

    authority = _load_json(authority_data, "failed signed v16 authority")
    authority_sources = authority.get("recovery_sources")
    expected_authority_sources = {
        role: record
        for role, record in source_records.items()
        if role != "archive_script"
    }
    review_contract = evidence.get("review_contract")
    if not isinstance(review_contract, Mapping):
        raise RecoveryAuthorityError("Failed v16 review contract missing")
    source_set_sha256 = _json_sha256(authority_sources)
    scope_sha256 = _json_sha256(authority.get("scope"))
    if (
        authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v16"
        or authority.get("status")
        != "APPROVED_FOR_TRANSITION_ALIAS_METADATA_PLUS_SIZE_AND_TERMINAL_KEY_RECOVERY_ONLY"
        or authority.get("pass") is not True
        or not _strict_equal(authority_sources, expected_authority_sources)
        or not _strict_equal(authority.get("public_key"), public_record)
        or not _strict_equal(authority.get("retired_private_key"), retired_private)
        or review_contract.get("reviewed_source_set_sha256") != source_set_sha256
        or review_contract.get("scope_sha256") != scope_sha256
        or review_contract.get("legacy_source_set_sha256")
        != _json_sha256(authority.get("legacy_sources"))
        or review_contract.get("prior_recovery_chain_sha256")
        != _json_sha256(authority.get("prior_recovery_chain"))
        or review_contract.get("rejected_generations_sha256")
        != _json_sha256(authority.get("rejected_unsigned_generations"))
        or review_contract.get("terminal_key_rotation_sha256")
        != _json_sha256(authority.get("terminal_key_rotation"))
    ):
        raise RecoveryAuthorityError("Failed v16 authority/source contract drift")

    authority_reviews = authority.get("review_artifacts")
    evidence_reviews = evidence.get("reviews")
    if not isinstance(authority_reviews, Mapping) or not isinstance(
        evidence_reviews, Mapping
    ):
        raise RecoveryAuthorityError("Failed v16 review records missing")
    reviewer_ids: set[str] = set()
    for role, (reviewer, reviewer_id) in FAILED_V16_EXPECTED_REVIEWERS.items():
        archive_key = FAILED_V16_ARCHIVED_REVIEW_KEYS[role]
        archived_record = archived_records[archive_key]
        original_record = authority_reviews.get(role)
        if (
            not isinstance(original_record, Mapping)
            or original_record.get("path")
            != str(FAILED_V16_ORIGINAL_REVIEW_PATHS[role])
            or any(
                original_record.get(field) != archived_record[field]
                for field in ("size_bytes", "sha256", "mode")
            )
        ):
            raise RecoveryAuthorityError(
                f"Failed v16 review archive transition drift: {role}"
            )
        review = _load_json(
            _open_bound(FAILED_V16_ARCHIVED_ARTIFACT_PATHS[archive_key])[0],
            f"failed v16 archived review {role}",
        )
        probe = review.get("readonly_probe")
        expected_review_evidence = {
            "reviewer": reviewer,
            "reviewer_agent_id": reviewer_id,
            "verdict": "APPROVE",
            "high_findings": 0,
            "medium_findings": 0,
            "unresolved_conditions": 0,
        }
        if (
            review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != SCHEMA_VERSION
            or review.get("reviewer") != reviewer
            or review.get("reviewer_agent_id") != reviewer_id
            or review.get("verdict") != "APPROVE"
            or review.get("reviewed_source_set_sha256") != source_set_sha256
            or review.get("scope_sha256") != scope_sha256
            or review.get("terminal_key_rotation_sha256")
            != review_contract.get("terminal_key_rotation_sha256")
            or not isinstance(probe, Mapping)
            or probe.get("forbidden_output_slots_checked") != 221
            or probe.get("regression_summary") != "380 passed"
            or probe.get("no_output_writes") is not True
            or review.get("high_findings") != []
            or review.get("medium_findings") != []
            or review.get("unresolved_conditions") != []
            or review.get("pass") is not True
            or not _strict_equal(
                evidence_reviews.get(role), expected_review_evidence
            )
        ):
            raise RecoveryAuthorityError(f"Failed v16 review contract drift: {role}")
        reviewer_ids.add(reviewer_id)
    if len(reviewer_ids) != len(FAILED_V16_EXPECTED_REVIEWERS):
        raise RecoveryAuthorityError("Failed v16 reviewer identities are not distinct")

    verification = _load_json(
        _open_bound(FAILED_V16_ARCHIVED_ARTIFACT_PATHS["verification_receipt"])[0],
        "failed v16 verification receipt",
    )
    replay = _load_json(
        _open_bound(FAILED_V16_ARCHIVED_ARTIFACT_PATHS["replay_receipt"])[0],
        "failed v16 replay receipt",
    )
    witness = _load_json(
        _open_bound(
            FAILED_V16_ARCHIVED_ARTIFACT_PATHS["replay_success_witness"]
        )[0],
        "failed v16 replay success witness",
    )
    embedded_authority = verification.get("schema_recovery_authority")
    worker_wait = witness.get("worker_wait")
    if (
        verification.get("schema_name")
        != "intersubmod.tumor_ref_source_receipt_promotion_verification.recovery"
        or verification.get("mode") != "verify-and-record"
        or verification.get("pass") is not True
        or not isinstance(embedded_authority, Mapping)
        or embedded_authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v16"
        or embedded_authority.get("runtime_role") != "continuation_verifier"
        or embedded_authority.get("pass") is not True
        or replay.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery"
        or replay.get("pass") is not True
        or witness.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery.success_witness"
        or witness.get("pass") is not True
        or not isinstance(worker_wait, Mapping)
        or worker_wait.get("exit_code") != 0
        or worker_wait.get("normal_exit") is not True
        or worker_wait.get("raw_wait_status") != 0
    ):
        raise RecoveryAuthorityError("Failed v16 V/R receipt contract drift")

    stderr_data = _open_bound(
        FAILED_V16_ARCHIVED_ARTIFACT_PATHS["formal_continuation_stderr"]
    )[0]
    stderr_payload = _load_json(stderr_data, "failed v16 continuation stderr")
    runtime = evidence.get("runtime")
    diagnostic = evidence.get("diagnostic")
    output_state = evidence.get("formal_output_state")
    expected_output_state = {
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
    }
    expected_legacy_path = (
        REPO_ROOT
        / "research/20260715_single_fp_alt_multicluster_subclone_validation/"
        "scripts/focal_alt_cluster_lib.py"
    )
    expected_record_shape = [
        "ctime_ns",
        "device",
        "inode",
        "mode",
        "mtime_ns",
        "path",
        "sha256",
        "size_bytes",
    ]

    def legacy_records(value: Any) -> list[Mapping[str, Any]]:
        records: list[Mapping[str, Any]] = []
        if isinstance(value, Mapping):
            if (
                value.get("path") == str(expected_legacy_path)
                and sorted(value) == expected_record_shape
            ):
                records.append(value)
            for child in value.values():
                records.extend(legacy_records(child))
        elif isinstance(value, list):
            for child in value:
                records.extend(legacy_records(child))
        return records

    authorization = _load_json(
        _open_bound(ORIGINAL_CHAIN_PATHS["promotion_authorization"])[0],
        "promotion authorization for failed v16 diagnostic",
    )
    completion = _load_json(
        _open_bound(ORIGINAL_CHAIN_PATHS["promotion_completion"])[0],
        "promotion completion for failed v16 diagnostic",
    )
    actual_legacy_records = legacy_records(authorization) + legacy_records(completion)
    occupied = [str(path) for path in FAILED_V16_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v16.bundle.staging.*"
        )
    )
    if (
        evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_signed_generation_failure"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("generation") != "v16"
        or evidence.get("status")
        != "SIGNED_AUTHORITY_V_AND_R_PASS_C_FAILED_BEFORE_CHILD_START"
        or not isinstance(runtime, Mapping)
        or runtime.get("v_exit_code") != 0
        or runtime.get("r_supervisor_exit_code") != 0
        or runtime.get("r_worker_exit_code") != 0
        or runtime.get("c_supervisor_exit_code") != 1
        or runtime.get("c_child_started") is not False
        or not _strict_equal(
            runtime.get("stderr"), archived_records["formal_continuation_stderr"]
        )
        or runtime.get("error_type") != "ContinuationError"
        or runtime.get("error_message")
        != "Ambiguous identity-like relation schema"
        or stderr_payload.get("error")
        != {
            "message": "Ambiguous identity-like relation schema",
            "type": "ContinuationError",
        }
        or stderr_payload.get("pass") is not False
        or not isinstance(diagnostic, Mapping)
        or diagnostic.get("legacy_relation_path") != str(expected_legacy_path)
        or diagnostic.get("record_shape") != expected_record_shape
        or diagnostic.get("link_count_present_in_historical_schema") is not False
        or diagnostic.get("failure_before_c_child_start") is not True
        or diagnostic.get("failure_before_first_c_output") is not True
        or diagnostic.get("scientific_payload_changed") is not False
        or len(actual_legacy_records) != 2
        or not _strict_equal(actual_legacy_records[0], actual_legacy_records[1])
        or not _strict_equal(diagnostic.get("record"), actual_legacy_records[0])
        or not _strict_equal(output_state, expected_output_state)
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("pass") is not False
        or occupied
        or staging
    ):
        raise RecoveryAuthorityError(
            f"Failed v16 runtime/archive contract drift: {occupied + staging}"
        )

    return {
        "evidence": evidence_record,
        "archived_artifacts": archived_records,
        "review_transport_artifacts": transport_records,
        "public_key": public_record,
        "retired_private_key": retired_private,
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": continuation_private_record,
        "authority_signature_verified": True,
        "atomic_commit_verified": True,
        "authority_created": True,
        "verification_receipt_created": True,
        "replay_receipt_and_log_created": True,
        "replay_success_witness_created": True,
        "continuation_outputs_created": False,
        "canonical_downstream_outputs_created": False,
        "legacy_eight_key_stat_schema_mismatch_reproduced": True,
        "source_set_sha256": source_set_sha256,
        "pass": True,
    }


def _validate_failed_v17_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V17_EVIDENCE)
    if (
        evidence_record["size_bytes"] != 40_696
        or evidence_record["sha256"] != EXPECTED_FAILED_V17_EVIDENCE_SHA256
    ):
        raise RecoveryAuthorityError("Failed v17 evidence identity drift")
    evidence = _load_json(evidence_data, "failed signed v17 recovery evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "generation",
            "status",
            "review_contract",
            "runtime",
            "diagnostic",
            "archived_artifacts",
            "source_artifacts",
            "review_transport_artifacts",
            "reviews",
            "keys",
            "formal_output_state",
            "scientific_payload_changed",
            "claim_ceiling_changed",
            "required_resolution",
            "pass",
        },
        "failed signed v17 recovery evidence",
    )

    archived_records = _records(FAILED_V17_ARCHIVED_ARTIFACT_PATHS)
    trace_paths = sorted(FAILED_V17_TRACE_ROOT.iterdir())
    trace_data_records = [_open_bound(path)[:2] for path in trace_paths]
    trace_records = [record for _, record in trace_data_records]
    complete_archived_records = {
        **archived_records,
        "diagnostic_write_traces": trace_records,
    }
    source_records = _records(FAILED_V17_SOURCE_PATHS)
    transport_records = _records(FAILED_V17_REVIEW_TRANSPORT_PATHS)
    if not _strict_equal(
        evidence.get("archived_artifacts"), complete_archived_records
    ):
        raise RecoveryAuthorityError("Failed v17 archived-artifact records drift")
    if not _strict_equal(evidence.get("source_artifacts"), source_records):
        raise RecoveryAuthorityError("Failed v17 source-artifact records drift")
    if not _strict_equal(
        evidence.get("review_transport_artifacts"), transport_records
    ):
        raise RecoveryAuthorityError("Failed v17 review-transport records drift")
    if any(record.get("mode") != "0o444" for record in source_records.values()):
        raise RecoveryAuthorityError("Failed v17 source set is not immutable mode 0444")
    if any(record.get("mode") != "0o444" for record in transport_records.values()):
        raise RecoveryAuthorityError(
            "Failed v17 review transport is not immutable mode 0444"
        )
    if (
        len(trace_records) != 49
        or any(record.get("mode") != "0o444" for record in trace_records)
    ):
        raise RecoveryAuthorityError("Failed v17 syscall trace inventory drift")

    public_data, public_record, public_fd = _open_bound(FAILED_V17_PUBLIC_KEY)
    retired_private, retired_private_fd = _open_retired_private_key_bound(
        FAILED_V17_PRIVATE_KEY
    )
    retired_private_stat = os.fstat(retired_private_fd)
    continuation_public_data, continuation_public_record, _ = _open_bound(
        FAILED_V17_CONTINUATION_PUBLIC_KEY
    )
    continuation_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V17_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="UNUSED_NOT_RETIRED_NO_SIGNATURE_CREATED_MUST_NOT_BE_REUSED",
    )
    expected_keys = {
        "authority_public_key": public_record,
        "authority_private_key": {
            "path": str(FAILED_V17_PRIVATE_KEY),
            "size_bytes": 119,
            "mode": "0o0",
            "link_count": 1,
            "state": "RETIRED_AFTER_SINGLE_AUTHORITY_SIGNATURE",
        },
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": {
            **continuation_private_record,
            "link_count": 1,
        },
    }
    if (
        len(public_data) != 113
        or public_record["sha256"]
        != "ef3413162e1cb9850fb966a1505d71cf185fe79812415b27a93f4d5887154eac"
        or retired_private_stat.st_size != 119
        or retired_private_stat.st_nlink != 1
        or len(continuation_public_data) != 113
        or continuation_public_record["sha256"]
        != "4b83e655f1a7a778691e27dc3df2257a230001c702afdf6e703e578c706e0b03"
        or not _strict_equal(evidence.get("keys"), expected_keys)
    ):
        raise RecoveryAuthorityError("Failed v17 key identity/lifecycle drift")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V17_ARCHIVED_ARTIFACT_PATHS["authority"]
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V17_ARCHIVED_ARTIFACT_PATHS["authority_signature"]
    )
    commit_data, _, _ = _open_bound(
        FAILED_V17_ARCHIVED_ARTIFACT_PATHS["authority_commit"]
    )
    if len(signature_data) != 64:
        raise RecoveryAuthorityError("Failed v17 signature size drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    _validate_bundle_commit(
        commit_data,
        authority_record,
        signature_record,
        retired_private,
    )

    authority = _load_json(authority_data, "failed signed v17 authority")
    authority_sources = authority.get("recovery_sources")
    expected_authority_sources = {
        role: record
        for role, record in source_records.items()
        if role != "archive_script"
    }
    review_contract = evidence.get("review_contract")
    if not isinstance(review_contract, Mapping):
        raise RecoveryAuthorityError("Failed v17 review contract missing")
    source_set_sha256 = _json_sha256(authority_sources)
    scope_sha256 = _json_sha256(authority.get("scope"))
    if (
        authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v17"
        or authority.get("status")
        != "APPROVED_FOR_TRANSITION_ALIAS_METADATA_PLUS_SIZE_AND_TERMINAL_KEY_RECOVERY_ONLY"
        or authority.get("pass") is not True
        or not _strict_equal(authority_sources, expected_authority_sources)
        or not _strict_equal(authority.get("public_key"), public_record)
        or not _strict_equal(authority.get("retired_private_key"), retired_private)
        or review_contract.get("reviewed_source_set_sha256") != source_set_sha256
        or review_contract.get("scope_sha256") != scope_sha256
        or review_contract.get("legacy_source_set_sha256")
        != _json_sha256(authority.get("legacy_sources"))
        or review_contract.get("prior_recovery_chain_sha256")
        != _json_sha256(authority.get("prior_recovery_chain"))
        or review_contract.get("rejected_generations_sha256")
        != _json_sha256(authority.get("rejected_unsigned_generations"))
        or review_contract.get("terminal_key_rotation_sha256")
        != _json_sha256(authority.get("terminal_key_rotation"))
    ):
        raise RecoveryAuthorityError("Failed v17 authority/source contract drift")

    authority_reviews = authority.get("review_artifacts")
    evidence_reviews = evidence.get("reviews")
    if not isinstance(authority_reviews, Mapping) or not isinstance(
        evidence_reviews, Mapping
    ):
        raise RecoveryAuthorityError("Failed v17 review records missing")
    reviewer_ids: set[str] = set()
    for role, (reviewer, reviewer_id) in FAILED_V17_EXPECTED_REVIEWERS.items():
        archive_key = FAILED_V17_ARCHIVED_REVIEW_KEYS[role]
        archived_record = archived_records[archive_key]
        original_record = authority_reviews.get(role)
        if (
            not isinstance(original_record, Mapping)
            or original_record.get("path")
            != str(FAILED_V17_ORIGINAL_REVIEW_PATHS[role])
            or any(
                original_record.get(field) != archived_record[field]
                for field in ("size_bytes", "sha256", "mode")
            )
        ):
            raise RecoveryAuthorityError(
                f"Failed v17 review archive transition drift: {role}"
            )
        review = _load_json(
            _open_bound(FAILED_V17_ARCHIVED_ARTIFACT_PATHS[archive_key])[0],
            f"failed v17 archived review {role}",
        )
        probe = review.get("readonly_probe")
        expected_review_evidence = {
            "reviewer": reviewer,
            "reviewer_agent_id": reviewer_id,
            "verdict": "APPROVE",
            "high_findings": 0,
            "medium_findings": 0,
            "unresolved_conditions": 0,
        }
        if (
            review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != SCHEMA_VERSION
            or review.get("reviewer") != reviewer
            or review.get("reviewer_agent_id") != reviewer_id
            or review.get("verdict") != "APPROVE"
            or review.get("reviewed_source_set_sha256") != source_set_sha256
            or review.get("scope_sha256") != scope_sha256
            or review.get("terminal_key_rotation_sha256")
            != review_contract.get("terminal_key_rotation_sha256")
            or not isinstance(probe, Mapping)
            or probe.get("forbidden_output_slots_checked") != 236
            or probe.get("regression_summary") != "403 passed"
            or probe.get("no_output_writes") is not True
            or review.get("high_findings") != []
            or review.get("medium_findings") != []
            or review.get("unresolved_conditions") != []
            or review.get("pass") is not True
            or not _strict_equal(
                evidence_reviews.get(role), expected_review_evidence
            )
        ):
            raise RecoveryAuthorityError(f"Failed v17 review contract drift: {role}")
        reviewer_ids.add(reviewer_id)
    if len(reviewer_ids) != len(FAILED_V17_EXPECTED_REVIEWERS):
        raise RecoveryAuthorityError("Failed v17 reviewer identities are not distinct")

    verification = _load_json(
        _open_bound(FAILED_V17_ARCHIVED_ARTIFACT_PATHS["verification_receipt"])[0],
        "failed v17 verification receipt",
    )
    replay = _load_json(
        _open_bound(FAILED_V17_ARCHIVED_ARTIFACT_PATHS["replay_receipt"])[0],
        "failed v17 replay receipt",
    )
    witness = _load_json(
        _open_bound(
            FAILED_V17_ARCHIVED_ARTIFACT_PATHS["replay_success_witness"]
        )[0],
        "failed v17 replay success witness",
    )
    embedded_authority = verification.get("schema_recovery_authority")
    worker_wait = witness.get("worker_wait")
    if (
        verification.get("schema_name")
        != "intersubmod.tumor_ref_source_receipt_promotion_verification.recovery"
        or verification.get("mode") != "verify-and-record"
        or verification.get("pass") is not True
        or not isinstance(embedded_authority, Mapping)
        or embedded_authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v17"
        or embedded_authority.get("runtime_role") != "continuation_verifier"
        or embedded_authority.get("pass") is not True
        or replay.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery"
        or replay.get("pass") is not True
        or witness.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery.success_witness"
        or witness.get("pass") is not True
        or not isinstance(worker_wait, Mapping)
        or worker_wait.get("exit_code") != 0
        or worker_wait.get("normal_exit") is not True
        or worker_wait.get("raw_wait_status") != 0
    ):
        raise RecoveryAuthorityError("Failed v17 V/R receipt contract drift")

    stderr_data = _open_bound(
        FAILED_V17_ARCHIVED_ARTIFACT_PATHS["formal_continuation_stderr"]
    )[0]
    stderr_payload = _load_json(stderr_data, "failed v17 continuation stderr")
    runtime = evidence.get("runtime")
    diagnostic = evidence.get("diagnostic")
    output_state = evidence.get("formal_output_state")
    expected_output_state = {
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
    }
    legacy_private_key = Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
        "20260722_m2v5_terminal_v2/ed25519_private_one_time_resident.pem"
    )

    def matching_relations(value: Any) -> list[Mapping[str, Any]]:
        records: list[Mapping[str, Any]] = []
        if isinstance(value, Mapping):
            if value.get("path") == str(legacy_private_key):
                records.append(value)
            for child in value.values():
                records.extend(matching_relations(child))
        elif isinstance(value, list):
            for child in value:
                records.extend(matching_relations(child))
        return records

    relation_records = matching_relations(replay)
    relation_shapes = sorted(sorted(record) for record in relation_records)
    expected_shapes = sorted(
        [
            ["mode", "path"],
            ["mode", "path", "size_bytes", "state"],
            ["must_remain_unretired", "path", "required_mode"],
        ]
    )
    child_error_message = (
        "Conflicting declared relation schema or value: "
        + str(legacy_private_key)
    )
    child_error_occurrences = sum(
        child_error_message.encode("ascii") in data
        for data, _ in trace_data_records
    )
    expected_parent_message = (
        "Supervised continuation child did not exit zero: returncode=1 "
        "stderr_sha256=6350b2e1a2c8f0a680571b25c63dc052d55580f01fac4b77a7edbcad7dcbede7"
    )
    occupied = [str(path) for path in FAILED_V17_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v17.bundle.staging.*"
        )
    )
    if (
        evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_signed_generation_failure"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("generation") != "v17"
        or evidence.get("status")
        != "SIGNED_AUTHORITY_V_AND_R_PASS_C_CHILD_FAILED_AFTER_FRESH_V_BEFORE_DOWNSTREAM"
        or not isinstance(runtime, Mapping)
        or runtime.get("v_exit_code") != 0
        or runtime.get("r_supervisor_exit_code") != 0
        or runtime.get("r_worker_exit_code") != 0
        or runtime.get("c_supervisor_exit_code") != 1
        or runtime.get("c_child_started") is not True
        or runtime.get("c_child_fresh_verifier_passed") is not True
        or runtime.get("error_type") != "ContinuationError"
        or runtime.get("error_message") != expected_parent_message
        or runtime.get("child_error_message") != child_error_message
        or runtime.get("child_stderr_sha256")
        != "6350b2e1a2c8f0a680571b25c63dc052d55580f01fac4b77a7edbcad7dcbede7"
        or not _strict_equal(
            runtime.get("stderr"), archived_records["formal_continuation_stderr"]
        )
        or stderr_payload.get("error")
        != {"message": expected_parent_message, "type": "ContinuationError"}
        or stderr_payload.get("pass") is not False
        or not isinstance(diagnostic, Mapping)
        or diagnostic.get("legacy_private_key_path") != str(legacy_private_key)
        or not _strict_equal(diagnostic.get("relation_records"), relation_records)
        or not _strict_equal(diagnostic.get("relation_shapes"), relation_shapes)
        or relation_shapes != expected_shapes
        or diagnostic.get("syscall_trace_count") != 49
        or diagnostic.get("syscall_trace_captured_exact_child_error") is not True
        or child_error_occurrences != 1
        or diagnostic.get("failure_before_c_child_start") is not False
        or diagnostic.get("failure_after_fresh_verifier_pass") is not True
        or diagnostic.get("failure_before_downstream_leaf_start") is not True
        or diagnostic.get("failure_before_first_c_output") is not True
        or diagnostic.get("scientific_payload_changed") is not False
        or not _strict_equal(output_state, expected_output_state)
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("pass") is not False
        or occupied
        or staging
    ):
        raise RecoveryAuthorityError(
            f"Failed v17 runtime/archive contract drift: {occupied + staging}"
        )

    return {
        "evidence": evidence_record,
        "archived_artifacts": complete_archived_records,
        "review_transport_artifacts": transport_records,
        "public_key": public_record,
        "retired_private_key": retired_private,
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": continuation_private_record,
        "authority_signature_verified": True,
        "atomic_commit_verified": True,
        "authority_created": True,
        "verification_receipt_created": True,
        "replay_receipt_and_log_created": True,
        "replay_success_witness_created": True,
        "continuation_child_started": True,
        "fresh_verifier_passed": True,
        "continuation_outputs_created": False,
        "canonical_downstream_outputs_created": False,
        "metadata_enrichment_conflict_reproduced": True,
        "syscall_trace_evidence_bound": True,
        "source_set_sha256": source_set_sha256,
        "pass": True,
    }


def _validate_failed_v18_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V18_EVIDENCE)
    if (
        evidence_record["size_bytes"] != 40_663
        or evidence_record["sha256"] != EXPECTED_FAILED_V18_EVIDENCE_SHA256
    ):
        raise RecoveryAuthorityError("Failed v18 evidence identity drift")
    evidence = _load_json(evidence_data, "failed signed v18 recovery evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "generation",
            "status",
            "review_contract",
            "runtime",
            "diagnostic",
            "archived_artifacts",
            "source_artifacts",
            "review_transport_artifacts",
            "reviews",
            "keys",
            "formal_output_state",
            "scientific_payload_changed",
            "claim_ceiling_changed",
            "required_resolution",
            "pass",
        },
        "failed signed v18 recovery evidence",
    )

    archived_records = _records(FAILED_V18_ARCHIVED_ARTIFACT_PATHS)
    trace_paths = sorted(FAILED_V18_TRACE_ROOT.glob(FAILED_V18_TRACE_GLOB))
    trace_data_records = [_open_bound(path)[:2] for path in trace_paths]
    trace_records = [record for _, record in trace_data_records]
    complete_archived_records = {
        **archived_records,
        "diagnostic_write_traces": trace_records,
    }
    source_records = _records(FAILED_V18_SOURCE_PATHS)
    declared_transport = evidence.get("review_transport_artifacts")
    if not isinstance(declared_transport, Mapping):
        raise RecoveryAuthorityError("Failed v18 review-transport records missing")
    transport_records: dict[str, dict[str, Any]] = {}
    for role, path in FAILED_V18_REVIEW_TRANSPORT_PATHS.items():
        declared = declared_transport.get(role)
        if (
            not isinstance(declared, Mapping)
            or declared.get("mode") not in {"0o444", "0o664"}
        ):
            raise RecoveryAuthorityError(
                f"Failed v18 review-transport mode contract drift: {role}"
            )
        transport_records[role] = _open_bound(
            path, expected_mode=str(declared["mode"])
        )[1]
    if not _strict_equal(
        evidence.get("archived_artifacts"), complete_archived_records
    ):
        raise RecoveryAuthorityError("Failed v18 archived-artifact records drift")
    if not _strict_equal(evidence.get("source_artifacts"), source_records):
        raise RecoveryAuthorityError("Failed v18 source-artifact records drift")
    if not _strict_equal(
        evidence.get("review_transport_artifacts"), transport_records
    ):
        raise RecoveryAuthorityError("Failed v18 review-transport records drift")
    if any(record.get("mode") != "0o444" for record in source_records.values()):
        raise RecoveryAuthorityError("Failed v18 source set is not immutable mode 0444")
    if (
        len(trace_records) != 52
        or any(record.get("mode") != "0o444" for record in trace_records)
        or archived_records["diagnostic_readme"].get("mode") != "0o444"
    ):
        raise RecoveryAuthorityError("Failed v18 diagnostic inventory drift")

    public_data, public_record, public_fd = _open_bound(FAILED_V18_PUBLIC_KEY)
    retired_private, retired_private_fd = _open_retired_private_key_bound(
        FAILED_V18_PRIVATE_KEY
    )
    retired_private_stat = os.fstat(retired_private_fd)
    continuation_public_data, continuation_public_record, _ = _open_bound(
        FAILED_V18_CONTINUATION_PUBLIC_KEY
    )
    continuation_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V18_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="UNUSED_NOT_RETIRED_NO_SIGNATURE_CREATED_MUST_NOT_BE_REUSED",
    )
    expected_keys = {
        "authority_public_key": public_record,
        "authority_private_key": {
            "path": str(FAILED_V18_PRIVATE_KEY),
            "size_bytes": 119,
            "mode": "0o0",
            "link_count": 1,
            "state": "RETIRED_AFTER_SINGLE_AUTHORITY_SIGNATURE",
        },
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": {
            **continuation_private_record,
            "link_count": 1,
        },
    }
    if (
        len(public_data) != 113
        or public_record["sha256"]
        != "57062811b5b67f514f1454cac3e96fe55f7ab94be2fbbbace43431f68fbb65e4"
        or retired_private_stat.st_size != 119
        or retired_private_stat.st_nlink != 1
        or len(continuation_public_data) != 113
        or continuation_public_record["sha256"]
        != "9fa7667e8076cee90a93fee44cfe08bc45a68e9ce28d01507b5c057463102a93"
        or not _strict_equal(evidence.get("keys"), expected_keys)
    ):
        raise RecoveryAuthorityError("Failed v18 key identity/lifecycle drift")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V18_ARCHIVED_ARTIFACT_PATHS["authority"]
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V18_ARCHIVED_ARTIFACT_PATHS["authority_signature"]
    )
    commit_data, _, _ = _open_bound(
        FAILED_V18_ARCHIVED_ARTIFACT_PATHS["authority_commit"]
    )
    if len(signature_data) != 64:
        raise RecoveryAuthorityError("Failed v18 signature size drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    _validate_bundle_commit(
        commit_data,
        authority_record,
        signature_record,
        retired_private,
    )

    authority = _load_json(authority_data, "failed signed v18 authority")
    authority_sources = authority.get("recovery_sources")
    expected_authority_sources = {
        role: record
        for role, record in source_records.items()
        if role != "archive_script"
    }
    review_contract = evidence.get("review_contract")
    if not isinstance(review_contract, Mapping):
        raise RecoveryAuthorityError("Failed v18 review contract missing")
    source_set_sha256 = _json_sha256(authority_sources)
    scope_sha256 = _json_sha256(authority.get("scope"))
    if (
        authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v18"
        or authority.get("status")
        != "APPROVED_FOR_TRANSITION_ALIAS_METADATA_PLUS_SIZE_AND_TERMINAL_KEY_RECOVERY_ONLY"
        or authority.get("pass") is not True
        or not _strict_equal(authority_sources, expected_authority_sources)
        or not _strict_equal(authority.get("public_key"), public_record)
        or not _strict_equal(authority.get("retired_private_key"), retired_private)
        or review_contract.get("reviewed_source_set_sha256") != source_set_sha256
        or review_contract.get("scope_sha256") != scope_sha256
        or review_contract.get("legacy_source_set_sha256")
        != _json_sha256(authority.get("legacy_sources"))
        or review_contract.get("prior_recovery_chain_sha256")
        != _json_sha256(authority.get("prior_recovery_chain"))
        or review_contract.get("rejected_generations_sha256")
        != _json_sha256(authority.get("rejected_unsigned_generations"))
        or review_contract.get("terminal_key_rotation_sha256")
        != _json_sha256(authority.get("terminal_key_rotation"))
    ):
        raise RecoveryAuthorityError("Failed v18 authority/source contract drift")

    authority_reviews = authority.get("review_artifacts")
    evidence_reviews = evidence.get("reviews")
    if not isinstance(authority_reviews, Mapping) or not isinstance(
        evidence_reviews, Mapping
    ):
        raise RecoveryAuthorityError("Failed v18 review records missing")
    reviewer_ids: set[str] = set()
    for role, (reviewer, reviewer_id) in FAILED_V18_EXPECTED_REVIEWERS.items():
        archive_key = FAILED_V18_ARCHIVED_REVIEW_KEYS[role]
        archived_record = archived_records[archive_key]
        original_record = authority_reviews.get(role)
        if (
            not isinstance(original_record, Mapping)
            or original_record.get("path")
            != str(FAILED_V18_ORIGINAL_REVIEW_PATHS[role])
            or any(
                original_record.get(field) != archived_record[field]
                for field in ("size_bytes", "sha256", "mode")
            )
        ):
            raise RecoveryAuthorityError(
                f"Failed v18 review archive transition drift: {role}"
            )
        review = _load_json(
            _open_bound(FAILED_V18_ARCHIVED_ARTIFACT_PATHS[archive_key])[0],
            f"failed v18 archived review {role}",
        )
        probe = review.get("readonly_probe")
        expected_review_evidence = {
            "reviewer": reviewer,
            "reviewer_agent_id": reviewer_id,
            "verdict": "APPROVE",
            "high_findings": 0,
            "medium_findings": 0,
            "unresolved_conditions": 0,
        }
        if (
            review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != SCHEMA_VERSION
            or review.get("reviewer") != reviewer
            or review.get("reviewer_agent_id") != reviewer_id
            or review.get("verdict") != "APPROVE"
            or review.get("reviewed_source_set_sha256") != source_set_sha256
            or review.get("scope_sha256") != scope_sha256
            or review.get("terminal_key_rotation_sha256")
            != review_contract.get("terminal_key_rotation_sha256")
            or not isinstance(probe, Mapping)
            or probe.get("forbidden_output_slots_checked") != 251
            or probe.get("regression_summary") != "427 passed"
            or probe.get("no_output_writes") is not True
            or review.get("high_findings") != []
            or review.get("medium_findings") != []
            or review.get("unresolved_conditions") != []
            or review.get("pass") is not True
            or not _strict_equal(
                evidence_reviews.get(role), expected_review_evidence
            )
        ):
            raise RecoveryAuthorityError(f"Failed v18 review contract drift: {role}")
        reviewer_ids.add(reviewer_id)
    if len(reviewer_ids) != len(FAILED_V18_EXPECTED_REVIEWERS):
        raise RecoveryAuthorityError("Failed v18 reviewer identities are not distinct")

    verification = _load_json(
        _open_bound(FAILED_V18_ARCHIVED_ARTIFACT_PATHS["verification_receipt"])[0],
        "failed v18 verification receipt",
    )
    replay = _load_json(
        _open_bound(FAILED_V18_ARCHIVED_ARTIFACT_PATHS["replay_receipt"])[0],
        "failed v18 replay receipt",
    )
    witness = _load_json(
        _open_bound(
            FAILED_V18_ARCHIVED_ARTIFACT_PATHS["replay_success_witness"]
        )[0],
        "failed v18 replay success witness",
    )
    embedded_authority = verification.get("schema_recovery_authority")
    worker_wait = witness.get("worker_wait")
    if (
        verification.get("schema_name")
        != "intersubmod.tumor_ref_source_receipt_promotion_verification.recovery"
        or verification.get("mode") != "verify-and-record"
        or verification.get("pass") is not True
        or not isinstance(embedded_authority, Mapping)
        or embedded_authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v18"
        or embedded_authority.get("runtime_role") != "continuation_verifier"
        or embedded_authority.get("pass") is not True
        or replay.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery"
        or replay.get("pass") is not True
        or witness.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery.success_witness"
        or witness.get("pass") is not True
        or not isinstance(worker_wait, Mapping)
        or worker_wait.get("exit_code") != 0
        or worker_wait.get("normal_exit") is not True
        or worker_wait.get("raw_wait_status") != 0
    ):
        raise RecoveryAuthorityError("Failed v18 V/R receipt contract drift")

    alias_before = os.lstat(TUMOR_REF_SUMMARY_ALIAS)
    alias_flags = os.O_PATH | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0)
    alias_fd = os.open(TUMOR_REF_SUMMARY_ALIAS, alias_flags)
    alias_opened = os.fstat(alias_fd)
    alias_text = os.readlink(TUMOR_REF_SUMMARY_ALIAS)
    alias_after = os.lstat(TUMOR_REF_SUMMARY_ALIAS)
    target_data, target_record, _ = _open_bound(
        TUMOR_REF_SUMMARY_TARGET, expected_mode="0o664"
    )
    alias_identity = lambda value: (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )
    expected_alias = {
        "path": str(TUMOR_REF_SUMMARY_ALIAS),
        "link_text": TUMOR_REF_SUMMARY_TARGET.name,
        "resolved_target": str(TUMOR_REF_SUMMARY_TARGET),
        "mode": "0o777",
        "link_count": 1,
    }
    if (
        not stat.S_ISLNK(alias_opened.st_mode)
        or alias_identity(alias_before) != alias_identity(alias_opened)
        or alias_identity(alias_opened) != alias_identity(alias_after)
        or alias_text != TUMOR_REF_SUMMARY_TARGET.name
        or TUMOR_REF_SUMMARY_ALIAS.resolve(strict=True) != TUMOR_REF_SUMMARY_TARGET
        or len(target_data) != 89_279
        or target_record["sha256"]
        != "9e777c9e011d6f235dc0bbc7182325469a29cc063b341a9a9f6ee669509aeebd"
        or target_record["mode"] != "0o664"
    ):
        os.close(alias_fd)
        raise RecoveryAuthorityError("Failed v18 tumor-REF alias/target drift")
    _FD_LEASES.append((TUMOR_REF_SUMMARY_ALIAS, alias_fd, alias_opened, None))

    stderr_data = _open_bound(
        FAILED_V18_ARCHIVED_ARTIFACT_PATHS["formal_continuation_stderr"]
    )[0]
    stderr_payload = _load_json(stderr_data, "failed v18 continuation stderr")
    runtime = evidence.get("runtime")
    diagnostic = evidence.get("diagnostic")
    output_state = evidence.get("formal_output_state")
    expected_output_state = {
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
    }
    child_error_message = (
        "Gate input path is not canonical: " + str(TUMOR_REF_SUMMARY_ALIAS)
    )
    child_error_occurrences = sum(
        child_error_message.encode("ascii") in data
        for data, _ in trace_data_records
    )
    expected_parent_message = (
        "Supervised continuation child did not exit zero: returncode=1 "
        "stderr_sha256=035cb2d90b4db9bdf0dd490a40f75d605bd9a404d1d838e5753b837a5c60b2da"
    )
    occupied = [str(path) for path in FAILED_V18_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v18.bundle.staging.*"
        )
    )
    if (
        evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_signed_generation_failure"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("generation") != "v18"
        or evidence.get("status")
        != "SIGNED_AUTHORITY_V_AND_R_PASS_C_CHILD_FAILED_AFTER_FRESH_V_BEFORE_DOWNSTREAM"
        or evidence.get("required_resolution")
        != "Create append-only v19 with fresh authority and terminal-v9 keys, strict tumor-REF alias/target binding, live regressions, fresh independent reviews, and distinct V/R/C output slots."
        or not isinstance(runtime, Mapping)
        or runtime.get("v_exit_code") != 0
        or runtime.get("r_supervisor_exit_code") != 0
        or runtime.get("r_worker_exit_code") != 0
        or runtime.get("c_supervisor_exit_code") != 1
        or runtime.get("c_child_started") is not True
        or runtime.get("c_child_fresh_verifier_passed") is not True
        or runtime.get("error_type") != "ContinuationError"
        or runtime.get("error_message") != expected_parent_message
        or runtime.get("child_error_message") != child_error_message
        or runtime.get("child_stderr_sha256")
        != "035cb2d90b4db9bdf0dd490a40f75d605bd9a404d1d838e5753b837a5c60b2da"
        or not _strict_equal(
            runtime.get("stderr"), archived_records["formal_continuation_stderr"]
        )
        or stderr_payload.get("error")
        != {"message": expected_parent_message, "type": "ContinuationError"}
        or stderr_payload.get("pass") is not False
        or not isinstance(diagnostic, Mapping)
        or not _strict_equal(
            diagnostic.get("tumor_ref_summary_alias"), expected_alias
        )
        or not _strict_equal(
            diagnostic.get("tumor_ref_summary_target"), target_record
        )
        or diagnostic.get("syscall_trace_count") != 52
        or diagnostic.get("syscall_trace_captured_exact_child_error") is not True
        or child_error_occurrences != 1
        or diagnostic.get("failure_before_c_child_start") is not False
        or diagnostic.get("failure_after_fresh_verifier_pass") is not True
        or diagnostic.get("failure_before_downstream_leaf_start") is not True
        or diagnostic.get("failure_before_first_c_output") is not True
        or diagnostic.get("scientific_payload_changed") is not False
        or not _strict_equal(output_state, expected_output_state)
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("pass") is not False
        or occupied
        or staging
    ):
        raise RecoveryAuthorityError(
            f"Failed v18 runtime/archive contract drift: {occupied + staging}"
        )

    return {
        "evidence": evidence_record,
        "archived_artifacts": complete_archived_records,
        "review_transport_artifacts": transport_records,
        "public_key": public_record,
        "retired_private_key": retired_private,
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": continuation_private_record,
        "tumor_ref_summary_alias": expected_alias,
        "tumor_ref_summary_target": target_record,
        "authority_signature_verified": True,
        "atomic_commit_verified": True,
        "authority_created": True,
        "verification_receipt_created": True,
        "replay_receipt_and_log_created": True,
        "replay_success_witness_created": True,
        "continuation_child_started": True,
        "fresh_verifier_passed": True,
        "continuation_outputs_created": False,
        "canonical_downstream_outputs_created": False,
        "tumor_ref_summary_alias_noncanonical_reproduced": True,
        "syscall_trace_evidence_bound": True,
        "source_set_sha256": source_set_sha256,
        "pass": True,
    }


def _validate_failed_v19_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V19_EVIDENCE)
    if (
        evidence_record["size_bytes"] != 18_407
        or evidence_record["sha256"] != EXPECTED_FAILED_V19_EVIDENCE_SHA256
    ):
        raise RecoveryAuthorityError("Failed v19 evidence identity drift")
    evidence = _load_json(evidence_data, "failed signed v19 recovery evidence")
    _exact_keys(
        evidence,
        {
            "schema_name",
            "schema_version",
            "created_at_utc",
            "generation",
            "status",
            "review_contract",
            "runtime",
            "diagnostic",
            "archived_artifacts",
            "source_artifacts",
            "review_transport_artifacts",
            "reviews",
            "keys",
            "formal_output_state",
            "scientific_payload_changed",
            "claim_ceiling_changed",
            "required_resolution",
            "pass",
        },
        "failed signed v19 recovery evidence",
    )

    archived_records = _records(FAILED_V19_ARCHIVED_ARTIFACT_PATHS)
    source_records = _records(FAILED_V19_SOURCE_PATHS)
    declared_transport = evidence.get("review_transport_artifacts")
    if not isinstance(declared_transport, Mapping):
        raise RecoveryAuthorityError("Failed v19 review-transport records missing")
    transport_records: dict[str, dict[str, Any]] = {}
    for role, path in FAILED_V19_REVIEW_TRANSPORT_PATHS.items():
        declared = declared_transport.get(role)
        if not isinstance(declared, Mapping) or declared.get("mode") != "0o444":
            raise RecoveryAuthorityError(
                f"Failed v19 review-transport mode contract drift: {role}"
            )
        transport_records[role] = _open_bound(path)[1]
    if not _strict_equal(evidence.get("archived_artifacts"), archived_records):
        raise RecoveryAuthorityError("Failed v19 archived-artifact records drift")
    if not _strict_equal(evidence.get("source_artifacts"), source_records):
        raise RecoveryAuthorityError("Failed v19 source-artifact records drift")
    if not _strict_equal(
        evidence.get("review_transport_artifacts"), transport_records
    ):
        raise RecoveryAuthorityError("Failed v19 review-transport records drift")
    if any(record.get("mode") != "0o444" for record in source_records.values()):
        raise RecoveryAuthorityError("Failed v19 source set is not immutable mode 0444")

    public_data, public_record, public_fd = _open_bound(FAILED_V19_PUBLIC_KEY)
    retired_private, retired_private_fd = _open_retired_private_key_bound(
        FAILED_V19_PRIVATE_KEY
    )
    retired_private_stat = os.fstat(retired_private_fd)
    continuation_public_data, continuation_public_record, _ = _open_bound(
        FAILED_V19_CONTINUATION_PUBLIC_KEY
    )
    continuation_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V19_CONTINUATION_PRIVATE_KEY,
        expected_mode=0o400,
        state="UNUSED_NOT_RETIRED_NO_SIGNATURE_CREATED_MUST_NOT_BE_REUSED",
    )
    expected_keys = {
        "authority_public_key": public_record,
        "authority_private_key": {
            "path": str(FAILED_V19_PRIVATE_KEY),
            "size_bytes": 119,
            "mode": "0o0",
            "link_count": 1,
            "state": "RETIRED_AFTER_SINGLE_AUTHORITY_SIGNATURE",
        },
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": {
            **continuation_private_record,
            "link_count": 1,
        },
    }
    if (
        len(public_data) != 113
        or public_record["sha256"]
        != "d494f66e8aea206e37e0e803ccb4e0ceb9cf2b244e71eccba6b370f96ddee2e0"
        or retired_private_stat.st_size != 119
        or retired_private_stat.st_nlink != 1
        or len(continuation_public_data) != 113
        or continuation_public_record["sha256"]
        != "34d11ff6a699b96aaa8624b30f45fbc8845e6958dc04122b19c413a1a2c12c2d"
        or not _strict_equal(evidence.get("keys"), expected_keys)
    ):
        raise RecoveryAuthorityError("Failed v19 key identity/lifecycle drift")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V19_ARCHIVED_ARTIFACT_PATHS["authority"]
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V19_ARCHIVED_ARTIFACT_PATHS["authority_signature"]
    )
    commit_data, _, _ = _open_bound(
        FAILED_V19_ARCHIVED_ARTIFACT_PATHS["authority_commit"]
    )
    if len(signature_data) != 64:
        raise RecoveryAuthorityError("Failed v19 signature size drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    _validate_bundle_commit(
        commit_data,
        authority_record,
        signature_record,
        retired_private,
    )

    authority = _load_json(authority_data, "failed signed v19 authority")
    authority_sources = authority.get("recovery_sources")
    expected_authority_sources = {
        role: record
        for role, record in source_records.items()
        if role != "archive_script"
    }
    review_contract = evidence.get("review_contract")
    if not isinstance(review_contract, Mapping):
        raise RecoveryAuthorityError("Failed v19 review contract missing")
    source_set_sha256 = _json_sha256(authority_sources)
    scope_sha256 = _json_sha256(authority.get("scope"))
    if (
        authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v19"
        or authority.get("status") != HISTORICAL_V19_V21_AUTHORITY_STATUS
        or authority.get("pass") is not True
        or not _strict_equal(authority_sources, expected_authority_sources)
        or not _strict_equal(authority.get("public_key"), public_record)
        or not _strict_equal(authority.get("retired_private_key"), retired_private)
        or review_contract.get("reviewed_source_set_sha256") != source_set_sha256
        or review_contract.get("scope_sha256") != scope_sha256
        or review_contract.get("legacy_source_set_sha256")
        != _json_sha256(authority.get("legacy_sources"))
        or review_contract.get("prior_recovery_chain_sha256")
        != _json_sha256(authority.get("prior_recovery_chain"))
        or review_contract.get("rejected_generations_sha256")
        != _json_sha256(authority.get("rejected_unsigned_generations"))
        or review_contract.get("terminal_key_rotation_sha256")
        != _json_sha256(authority.get("terminal_key_rotation"))
    ):
        raise RecoveryAuthorityError("Failed v19 authority/source contract drift")

    authority_reviews = authority.get("review_artifacts")
    evidence_reviews = evidence.get("reviews")
    if not isinstance(authority_reviews, Mapping) or not isinstance(
        evidence_reviews, Mapping
    ):
        raise RecoveryAuthorityError("Failed v19 review records missing")
    reviewer_ids: set[str] = set()
    for role, (reviewer, reviewer_id) in FAILED_V19_EXPECTED_REVIEWERS.items():
        archive_key = FAILED_V19_ARCHIVED_REVIEW_KEYS[role]
        archived_record = archived_records[archive_key]
        original_record = authority_reviews.get(role)
        if (
            not isinstance(original_record, Mapping)
            or original_record.get("path")
            != str(FAILED_V19_ORIGINAL_REVIEW_PATHS[role])
            or any(
                original_record.get(field) != archived_record[field]
                for field in ("size_bytes", "sha256", "mode")
            )
        ):
            raise RecoveryAuthorityError(
                f"Failed v19 review archive transition drift: {role}"
            )
        review = _load_json(
            _open_bound(FAILED_V19_ARCHIVED_ARTIFACT_PATHS[archive_key])[0],
            f"failed v19 archived review {role}",
        )
        probe = review.get("readonly_probe")
        expected_review_evidence = {
            "reviewer": reviewer,
            "reviewer_agent_id": reviewer_id,
            "verdict": "APPROVE",
            "high_findings": 0,
            "medium_findings": 0,
            "unresolved_conditions": 0,
        }
        if (
            review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != SCHEMA_VERSION
            or review.get("reviewer") != reviewer
            or review.get("reviewer_agent_id") != reviewer_id
            or review.get("verdict") != "APPROVE"
            or review.get("reviewed_source_set_sha256") != source_set_sha256
            or review.get("scope_sha256") != scope_sha256
            or review.get("terminal_key_rotation_sha256")
            != review_contract.get("terminal_key_rotation_sha256")
            or not isinstance(probe, Mapping)
            or probe.get("forbidden_output_slots_checked") != 266
            or probe.get("regression_summary") != "436 passed"
            or probe.get("no_output_writes") is not True
            or review.get("high_findings") != []
            or review.get("medium_findings") != []
            or review.get("unresolved_conditions") != []
            or review.get("pass") is not True
            or not _strict_equal(
                evidence_reviews.get(role), expected_review_evidence
            )
        ):
            raise RecoveryAuthorityError(f"Failed v19 review contract drift: {role}")
        reviewer_ids.add(reviewer_id)
    if len(reviewer_ids) != len(FAILED_V19_EXPECTED_REVIEWERS):
        raise RecoveryAuthorityError("Failed v19 reviewer identities are not distinct")

    verification = _load_json(
        _open_bound(FAILED_V19_ARCHIVED_ARTIFACT_PATHS["verification_receipt"])[0],
        "failed v19 verification receipt",
    )
    replay = _load_json(
        _open_bound(FAILED_V19_ARCHIVED_ARTIFACT_PATHS["replay_receipt"])[0],
        "failed v19 replay receipt",
    )
    witness = _load_json(
        _open_bound(
            FAILED_V19_ARCHIVED_ARTIFACT_PATHS["replay_success_witness"]
        )[0],
        "failed v19 replay success witness",
    )
    embedded_authority = verification.get("schema_recovery_authority")
    worker_wait = witness.get("worker_wait")
    if (
        verification.get("schema_name")
        != "intersubmod.tumor_ref_source_receipt_promotion_verification.recovery"
        or verification.get("mode") != "verify-and-record"
        or verification.get("pass") is not True
        or not isinstance(embedded_authority, Mapping)
        or embedded_authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v19"
        or embedded_authority.get("runtime_role") != "continuation_verifier"
        or embedded_authority.get("pass") is not True
        or replay.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery"
        or replay.get("pass") is not True
        or witness.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery.success_witness"
        or witness.get("pass") is not True
        or not isinstance(worker_wait, Mapping)
        or worker_wait.get("exit_code") != 0
        or worker_wait.get("normal_exit") is not True
        or worker_wait.get("raw_wait_status") != 0
    ):
        raise RecoveryAuthorityError("Failed v19 V/R receipt contract drift")

    stderr_data = _open_bound(
        FAILED_V19_ARCHIVED_ARTIFACT_PATHS["formal_continuation_stderr"]
    )[0]
    stderr_payload = _load_json(stderr_data, "failed v19 continuation stderr")
    strace_data = _open_bound(
        FAILED_V19_ARCHIVED_ARTIFACT_PATHS["diagnostic_strace_write_syscalls"]
    )[0]
    runtime = evidence.get("runtime")
    diagnostic = evidence.get("diagnostic")
    output_state = evidence.get("formal_output_state")
    expected_output_state = {
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
    }
    expected_parent_message = (
        "Supervised continuation child did not exit zero: returncode=1 "
        f"stderr_sha256={FAILED_V19_CHILD_STDERR_SHA256}"
    )
    occupied = [str(path) for path in FAILED_V19_OUTPUT_SLOTS if os.path.lexists(path)]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v19.bundle.staging.*"
        )
    )
    retired_input = diagnostic.get("retired_release_private_key") if isinstance(
        diagnostic, Mapping
    ) else None
    if (
        evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_signed_generation_failure"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("generation") != "v19"
        or evidence.get("status")
        != "SIGNED_AUTHORITY_V_AND_R_PASS_C_CHILD_FAILED_AFTER_FRESH_V_BEFORE_DOWNSTREAM"
        or evidence.get("required_resolution")
        != "Create append-only v20 with fresh authority and terminal-v10 keys, strict mode-000 parent-watch fallback, live regressions, fresh independent reviews, and distinct V/R/C output slots."
        or not isinstance(runtime, Mapping)
        or runtime.get("v_exit_code") != 0
        or runtime.get("r_supervisor_exit_code") != 0
        or runtime.get("r_worker_exit_code") != 0
        or runtime.get("c_supervisor_exit_code") != 1
        or runtime.get("c_child_started") is not True
        or runtime.get("c_child_fresh_verifier_passed") is not True
        or runtime.get("error_type") != "ContinuationError"
        or runtime.get("error_message") != expected_parent_message
        or runtime.get("child_error_message") != FAILED_V19_CHILD_ERROR_MESSAGE
        or runtime.get("child_stderr_sha256") != FAILED_V19_CHILD_STDERR_SHA256
        or not _strict_equal(
            runtime.get("stderr"), archived_records["formal_continuation_stderr"]
        )
        or stderr_payload.get("error")
        != {"message": expected_parent_message, "type": "ContinuationError"}
        or stderr_payload.get("pass") is not False
        or not isinstance(diagnostic, Mapping)
        or diagnostic.get("strace_exit_code") != 1
        or diagnostic.get("syscall_trace_captured_exact_child_error") is not True
        or strace_data.count(FAILED_V19_CHILD_ERROR_MESSAGE.encode("ascii")) != 1
        or FAILED_V19_CHILD_STDERR_SHA256.encode("ascii") not in strace_data
        or not isinstance(retired_input, Mapping)
        or retired_input.get("mode") != "0o0"
        or retired_input.get("size_bytes") != 290
        or retired_input.get("link_count") != 1
        or diagnostic.get("failure_before_c_child_start") is not False
        or diagnostic.get("failure_after_fresh_verifier_pass") is not True
        or diagnostic.get("failure_before_downstream_leaf_start") is not True
        or diagnostic.get("failure_before_first_c_output") is not True
        or diagnostic.get("scientific_payload_changed") is not False
        or not _strict_equal(output_state, expected_output_state)
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("pass") is not False
        or occupied
        or staging
    ):
        raise RecoveryAuthorityError(
            f"Failed v19 runtime/archive contract drift: {occupied + staging}"
        )

    return {
        "evidence": evidence_record,
        "archived_artifacts": archived_records,
        "review_transport_artifacts": transport_records,
        "public_key": public_record,
        "retired_private_key": retired_private,
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": continuation_private_record,
        "authority_signature_verified": True,
        "atomic_commit_verified": True,
        "authority_created": True,
        "verification_receipt_created": True,
        "replay_receipt_and_log_created": True,
        "replay_success_witness_created": True,
        "continuation_child_started": True,
        "fresh_verifier_passed": True,
        "continuation_outputs_created": False,
        "canonical_downstream_outputs_created": False,
        "mode000_inotify_eacces_reproduced": True,
        "syscall_trace_evidence_bound": True,
        "source_set_sha256": source_set_sha256,
        "pass": True,
    }


def _validate_failed_v21_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V21_EVIDENCE)
    if evidence_record["sha256"] != EXPECTED_FAILED_V21_EVIDENCE_SHA256:
        raise RecoveryAuthorityError("Failed v21 evidence identity drift")
    evidence = _load_json(evidence_data, "failed signed v21 recovery evidence")
    _exact_keys(
        evidence,
        {
            "archived_artifacts",
            "claim_ceiling_changed",
            "created_at_utc",
            "diagnostic",
            "formal_output_state",
            "generation",
            "keys",
            "observed_downstream_outputs",
            "pass",
            "required_resolution",
            "review_contract",
            "reviews",
            "schema_name",
            "schema_version",
            "scientific_payload_changed",
            "source_artifacts",
            "status",
        },
        "failed signed v21 recovery evidence",
    )
    if (
        evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery.failed_formal_run"
        or evidence.get("schema_version") != SCHEMA_VERSION
        or evidence.get("generation") != "v21"
        or evidence.get("status")
        != "ARCHIVED_FAILED_AFTER_INTERMEDIATE_DOWNSTREAM_OUTPUTS_BEFORE_FINAL_DATASET"
        or evidence.get("required_resolution")
        != "Create append-only v22 with fresh authority and terminal-v12 keys, FD-bound canonical-argv0 Python execution, focused real-exec regressions, fresh independent reviews, and distinct V/R/C output slots."
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("pass") is not False
    ):
        raise RecoveryAuthorityError("Failed v21 top-level evidence contract drift")

    source_records = _records(FAILED_V21_SOURCE_PATHS)
    if (
        {role: record["sha256"] for role, record in source_records.items()}
        != EXPECTED_FAILED_V21_SOURCE_SHA256
        or evidence.get("source_artifacts") != EXPECTED_FAILED_V21_SOURCE_SHA256
    ):
        raise RecoveryAuthorityError("Failed v21 source identity drift")
    archived_records = _records(FAILED_V21_ARCHIVED_ARTIFACT_PATHS)
    if (
        {role: record["sha256"] for role, record in archived_records.items()}
        != EXPECTED_FAILED_V21_ARCHIVED_SHA256
    ):
        raise RecoveryAuthorityError("Failed v21 archived artifact identity drift")
    transport_records = _records(FAILED_V21_REVIEW_TRANSPORT_PATHS)
    if (
        {role: record["sha256"] for role, record in transport_records.items()}
        != EXPECTED_FAILED_V21_TRANSPORT_SHA256
    ):
        raise RecoveryAuthorityError("Failed v21 review transport identity drift")

    declared_archived = evidence.get("archived_artifacts")
    expected_declared_archived_roles = {
        "authority",
        "authority_commit",
        "authority_signature",
        "continuation_incident",
        "replay_log",
        "replay_receipt",
        "replay_success_witness",
        "verification_receipt",
    }
    if (
        not isinstance(declared_archived, Mapping)
        or set(declared_archived) != expected_declared_archived_roles
    ):
        raise RecoveryAuthorityError("Failed v21 archived evidence is missing")
    for role, declared in declared_archived.items():
        if role not in archived_records or not isinstance(declared, Mapping):
            raise RecoveryAuthorityError(f"Failed v21 archived role drift: {role}")
        observed = archived_records[role]
        if (
            declared.get("mode") != observed["mode"]
            or declared.get("size_bytes") != observed["size_bytes"]
            or declared.get("sha256") != observed["sha256"]
            or declared.get("path") != FAILED_V21_ARCHIVED_ARTIFACT_PATHS[role].relative_to(
                FAILED_V21_ROOT
            ).as_posix()
        ):
            raise RecoveryAuthorityError(
                f"Failed v21 archived declaration drift: {role}"
            )

    declared_observed = evidence.get("observed_downstream_outputs")
    observed_role_map = {
        "strict_confirmation_receipt": "strict_not_applicable",
        "matched_normal_receipt": "matched_normal_not_applicable",
        "candidate_cn_ccf_receipt": "cn_ccf_receipt",
        "primary_artifact_audit": "primary_post_audit",
        "frozen_input_audit": "frozen_post_audit",
    }
    if (
        not isinstance(declared_observed, Mapping)
        or set(declared_observed) != set(observed_role_map)
    ):
        raise RecoveryAuthorityError("Failed v21 observed-output inventory drift")
    for declared_role, archive_role in observed_role_map.items():
        declared = declared_observed[declared_role]
        observed = archived_records[archive_role]
        if (
            not isinstance(declared, Mapping)
            or declared.get("path")
            != FAILED_V21_ARCHIVED_ARTIFACT_PATHS[archive_role].relative_to(
                FAILED_V21_ROOT
            ).as_posix()
            or declared.get("sha256") != observed["sha256"]
        ):
            raise RecoveryAuthorityError(
                f"Failed v21 observed-output cross-link drift: {declared_role}"
            )

    public_data, public_record, public_fd = _open_bound(FAILED_V21_PUBLIC_KEY)
    retired_private, retired_private_fd = _open_retired_private_key_bound(
        FAILED_V21_PRIVATE_KEY
    )
    retired_private_stat = os.fstat(retired_private_fd)
    continuation_public_data, continuation_public_record, _ = _open_bound(
        FAILED_V21_CONTINUATION_PUBLIC_KEY
    )
    continuation_private_record, continuation_private_fd = (
        _open_private_key_metadata_bound(
            FAILED_V21_CONTINUATION_PRIVATE_KEY,
            expected_mode=0o400,
            state="UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        )
    )
    continuation_private_stat = os.fstat(continuation_private_fd)
    expected_keys = {
        "authority_private_key": {
            "link_count": 1,
            "mode": "0o0",
            "path": str(FAILED_V21_PRIVATE_KEY),
            "size_bytes": 119,
            "state": "RETIRED_AFTER_SINGLE_AUTHORITY_SIGNATURE",
        },
        "authority_public_key": public_record,
        "continuation_private_key": {
            "link_count": 1,
            "mode": "0o400",
            "path": str(FAILED_V21_CONTINUATION_PRIVATE_KEY),
            "size_bytes": 119,
            "state": "UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        },
        "continuation_public_key": continuation_public_record,
    }
    if (
        len(public_data) != 113
        or public_record["sha256"]
        != "b66aa240810122b7d9be0ff89840339a7a58563163a55b5d2a837178dffebafb"
        or retired_private_stat.st_size != 119
        or retired_private_stat.st_nlink != 1
        or len(continuation_public_data) != 113
        or continuation_public_record["sha256"]
        != "825111aa7dcd25e60c6357243e228422d2cf07855d682aae30d90ebdcd5559d2"
        or continuation_private_record["size_bytes"] != 119
        or continuation_private_stat.st_nlink != 1
        or not _strict_equal(evidence.get("keys"), expected_keys)
    ):
        raise RecoveryAuthorityError("Failed v21 key identity/lifecycle drift")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V21_ARCHIVED_ARTIFACT_PATHS["authority"]
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V21_ARCHIVED_ARTIFACT_PATHS["authority_signature"]
    )
    commit_data, _, _ = _open_bound(
        FAILED_V21_ARCHIVED_ARTIFACT_PATHS["authority_commit"]
    )
    if len(signature_data) != 64:
        raise RecoveryAuthorityError("Failed v21 signature size drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    _validate_bundle_commit(
        commit_data, authority_record, signature_record, retired_private
    )
    authority = _load_json(authority_data, "failed signed v21 authority")
    review_contract = evidence.get("review_contract")
    authority_sources = authority.get("recovery_sources")
    if (
        not isinstance(review_contract, Mapping)
        or authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v21"
        or authority.get("status") != HISTORICAL_V19_V21_AUTHORITY_STATUS
        or authority.get("pass") is not True
        or not _strict_equal(authority_sources, source_records)
        or not _strict_equal(authority.get("public_key"), public_record)
        or not _strict_equal(authority.get("retired_private_key"), retired_private)
        or review_contract.get("reviewed_source_set_sha256")
        != _json_sha256(authority_sources)
        or review_contract.get("scope_sha256") != _json_sha256(authority.get("scope"))
        or review_contract.get("legacy_source_set_sha256")
        != _json_sha256(authority.get("legacy_sources"))
        or review_contract.get("prior_recovery_chain_sha256")
        != _json_sha256(authority.get("prior_recovery_chain"))
        or review_contract.get("rejected_generations_sha256")
        != _json_sha256(authority.get("rejected_unsigned_generations"))
        or review_contract.get("terminal_key_rotation_sha256")
        != _json_sha256(authority.get("terminal_key_rotation"))
    ):
        raise RecoveryAuthorityError("Failed v21 authority/source contract drift")

    authority_reviews = authority.get("review_artifacts")
    evidence_reviews = evidence.get("reviews")
    archived_review_keys = {
        "mendel": "review_mendel",
        "nash": "review_nash",
        "external_claude_opus": "review_external_claude_opus",
    }
    if not isinstance(authority_reviews, Mapping) or not isinstance(
        evidence_reviews, Mapping
    ):
        raise RecoveryAuthorityError("Failed v21 review records are missing")
    reviewer_ids: set[str] = set()
    for role, (reviewer, reviewer_id) in FAILED_V21_EXPECTED_REVIEWERS.items():
        archive_key = archived_review_keys[role]
        archive_record = archived_records[archive_key]
        original_record = authority_reviews.get(role)
        if (
            not isinstance(original_record, Mapping)
            or original_record.get("path") != str(FAILED_V21_ORIGINAL_REVIEW_PATHS[role])
            or any(
                original_record.get(field) != archive_record[field]
                for field in ("size_bytes", "sha256", "mode")
            )
        ):
            raise RecoveryAuthorityError(
                f"Failed v21 review archive transition drift: {role}"
            )
        review = _load_json(
            _open_bound(FAILED_V21_ARCHIVED_ARTIFACT_PATHS[archive_key])[0],
            f"failed v21 archived review {role}",
        )
        declared_review = evidence_reviews.get(role)
        if (
            review.get("schema_name") != REVIEW_SCHEMA_NAME
            or review.get("schema_version") != SCHEMA_VERSION
            or review.get("reviewer") != reviewer
            or review.get("reviewer_agent_id") != reviewer_id
            or review.get("verdict") != "APPROVE"
            or review.get("reviewed_source_set_sha256")
            != review_contract.get("reviewed_source_set_sha256")
            or review.get("scope_sha256") != review_contract.get("scope_sha256")
            or review.get("terminal_key_rotation_sha256")
            != review_contract.get("terminal_key_rotation_sha256")
            or review.get("high_findings") != []
            or review.get("medium_findings") != []
            or review.get("unresolved_conditions") != []
            or review.get("pass") is not True
            or not isinstance(declared_review, Mapping)
            or declared_review.get("reviewer_agent_id") != reviewer_id
            or declared_review.get("sha256") != archive_record["sha256"]
            or declared_review.get("verdict") != "APPROVE"
            or declared_review.get("high_findings") != 0
            or declared_review.get("medium_findings") != 0
            or declared_review.get("unresolved_conditions") != 0
        ):
            raise RecoveryAuthorityError(f"Failed v21 review contract drift: {role}")
        reviewer_ids.add(reviewer_id)
    if len(reviewer_ids) != 3:
        raise RecoveryAuthorityError("Failed v21 reviewer identities are not distinct")

    verification = _load_json(
        _open_bound(FAILED_V21_ARCHIVED_ARTIFACT_PATHS["verification_receipt"])[0],
        "failed v21 verification receipt",
    )
    replay = _load_json(
        _open_bound(FAILED_V21_ARCHIVED_ARTIFACT_PATHS["replay_receipt"])[0],
        "failed v21 replay receipt",
    )
    witness = _load_json(
        _open_bound(FAILED_V21_ARCHIVED_ARTIFACT_PATHS["replay_success_witness"])[0],
        "failed v21 replay success witness",
    )
    incident = _load_json(
        _open_bound(FAILED_V21_ARCHIVED_ARTIFACT_PATHS["continuation_incident"])[0],
        "failed v21 continuation incident",
    )
    embedded_authority = verification.get("schema_recovery_authority")
    worker_wait = witness.get("worker_wait")
    expected_error = (
        "Supervised continuation child did not exit zero: returncode=1 "
        "stderr_sha256=58a2fa58167d886c83d0465a56da2fc5513d985e47b70cf1c3c04d998657748b"
    )
    if (
        verification.get("schema_name")
        != "intersubmod.tumor_ref_source_receipt_promotion_verification.recovery"
        or verification.get("mode") != "verify-and-record"
        or verification.get("pass") is not True
        or not isinstance(embedded_authority, Mapping)
        or embedded_authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v21"
        or embedded_authority.get("runtime_role") != "continuation_verifier"
        or embedded_authority.get("pass") is not True
        or replay.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery"
        or replay.get("pass") is not True
        or witness.get("schema_name")
        != "intersubmod.m2v5_runner_only_gate_replay.recovery.success_witness"
        or witness.get("pass") is not True
        or not isinstance(worker_wait, Mapping)
        or worker_wait.get("normal_exit") is not True
        or worker_wait.get("exit_code") != 0
        or worker_wait.get("raw_wait_status") != 0
        or incident.get("schema_name")
        != "intersubmod.m2v5_downstream_continuation.recovery.incident"
        or incident.get("status")
        != "DOWNSTREAM_ARTIFACTS_EXIST_WITHOUT_CLEAN_TERMINAL_RETURN"
        or incident.get("error")
        != {"message": expected_error, "type": "ContinuationError"}
        or incident.get("release_authority_granted") is not False
        or incident.get("continuation_receipt_exists") is not False
        or incident.get("continuation_signature_exists") is not False
        or incident.get("supervisor_success_witness_exists") is not False
        or incident.get("pass") is not False
    ):
        raise RecoveryAuthorityError("Failed v21 V/R/incident contract drift")

    strict_receipt = _load_json(
        _open_bound(FAILED_V21_ARCHIVED_ARTIFACT_PATHS["strict_not_applicable"])[0],
        "failed v21 strict receipt",
    )
    normal_receipt = _load_json(
        _open_bound(
            FAILED_V21_ARCHIVED_ARTIFACT_PATHS["matched_normal_not_applicable"]
        )[0],
        "failed v21 matched-normal receipt",
    )
    cn_receipt = _load_json(
        _open_bound(FAILED_V21_ARCHIVED_ARTIFACT_PATHS["cn_ccf_receipt"])[0],
        "failed v21 CN/CCF receipt",
    )
    primary_audit = _load_json(
        _open_bound(FAILED_V21_ARCHIVED_ARTIFACT_PATHS["primary_post_audit"])[0],
        "failed v21 primary audit",
    )
    frozen_audit = _load_json(
        _open_bound(FAILED_V21_ARCHIVED_ARTIFACT_PATHS["frozen_post_audit"])[0],
        "failed v21 frozen-input audit",
    )
    continuation_log = _open_bound(
        FAILED_V21_ARCHIVED_ARTIFACT_PATHS["continuation_log"]
    )[0]
    if (
        strict_receipt.get("status") != "NOT_APPLICABLE"
        or strict_receipt.get("reason") != "ZERO_SELECTED_CANDIDATES"
        or strict_receipt.get("counts", {}).get("n_selected_candidates") != 0
        or strict_receipt.get("pass") is not True
        or normal_receipt.get("status") != "NOT_APPLICABLE"
        or normal_receipt.get("reason") != "ZERO_SELECTED_CANDIDATES"
        or normal_receipt.get("counts", {}).get("n_selected_candidates") != 0
        or normal_receipt.get("pass") is not True
        or cn_receipt.get("status") != "NOT_APPLICABLE"
        or cn_receipt.get("pass") is not True
        or primary_audit.get("counts", {}).get("stable_sites") != 102_842
        or primary_audit.get("counts", {}).get("primary_artifacts_verified")
        != 308_526
        or primary_audit.get("verification", {}).get("artifact_set_sha256")
        != "195e77d571576f37debf1627edb57f9dc7edb935849bd0ae6e29b323b380b2ca"
        or primary_audit.get("pass") is not True
        or frozen_audit.get("totals")
        != {
            "n_samples": 7,
            "n_sample_pass": 7,
            "n_artifacts": 77,
            "n_artifact_pass": 77,
        }
        or frozen_audit.get("pass") is not True
        or continuation_log.count(
            b"ContractError: pre-downstream primary artifact audit command is not canonical"
        )
        != 1
    ):
        raise RecoveryAuthorityError("Failed v21 downstream evidence contract drift")

    diagnostic = evidence.get("diagnostic")
    output_state = evidence.get("formal_output_state")
    if (
        not isinstance(diagnostic, Mapping)
        or diagnostic.get("actual_downstream_python_launcher") != "/proc/self/fd/705"
        or diagnostic.get("canonical_python_launcher") != str(PYTHON)
        or diagnostic.get("failure_after_downstream_leaf_start") is not True
        or diagnostic.get("failure_before_final_dataset_creation") is not True
        or diagnostic.get("scientific_payload_changed") is not False
        or not isinstance(output_state, Mapping)
        or output_state.get("authority_created_and_archived") is not True
        or output_state.get("canonical_downstream_intermediate_outputs_created_and_archived")
        is not True
        or output_state.get("final_dataset_directory_created") is not False
        or output_state.get("continuation_receipt_created") is not False
        or output_state.get("continuation_signature_created") is not False
    ):
        raise RecoveryAuthorityError("Failed v21 diagnostic/output-state drift")

    occupied = [
        str(path)
        for path in (
            *FAILED_V21_OUTPUT_SLOTS,
            *FAILED_V21_ORIGINAL_OBSERVED_OUTPUT_SLOTS,
        )
        if os.path.lexists(path)
    ]
    staging = sorted(
        str(path)
        for path in RESULT_ROOT.glob(
            ".tumor_ref_promotion_schema_recovery_authority.v21.bundle.staging.*"
        )
    )
    if occupied or staging:
        raise RecoveryAuthorityError(
            f"Failed v21 original output slots are occupied: {occupied + staging}"
        )
    cache_archive = (
        FAILED_V21_ROOT
        / "observed_downstream_outputs"
        / ".python_cache_m2v5_completion_v2_bound_bootstrap"
    )
    if not cache_archive.is_dir():
        raise RecoveryAuthorityError("Failed v21 Python cache archive is absent")

    return {
        "evidence": evidence_record,
        "archived_artifacts": archived_records,
        "review_transport_artifacts": transport_records,
        "public_key": public_record,
        "retired_private_key": retired_private,
        "continuation_public_key": continuation_public_record,
        "continuation_private_key": continuation_private_record,
        "authority_signature_verified": True,
        "atomic_commit_verified": True,
        "verification_receipt_created": True,
        "replay_receipt_log_and_success_witness_created": True,
        "continuation_incident_created": True,
        "intermediate_downstream_outputs_archived": True,
        "original_observed_downstream_slots_absent": True,
        "final_dataset_created": False,
        "continuation_terminal_outputs_created": False,
        "canonical_launcher_relation_mismatch_reproduced": True,
        "source_set_sha256": _json_sha256(source_records),
        "scientific_payload_changed": False,
        "claim_ceiling_changed": False,
        "pass": True,
    }


def _validate_failed_v25_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V25_EVIDENCE)
    if evidence_record["sha256"] != EXPECTED_FAILED_V25_EVIDENCE_SHA256:
        raise RecoveryAuthorityError("Failed v25 evidence identity drift")
    evidence = _load_json(evidence_data, "failed v25 evidence")
    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V25_BUNDLE / "authority.json"
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V25_BUNDLE / "authority.ed25519.sig"
    )
    commit_data, commit_record, _ = _open_bound(FAILED_V25_BUNDLE / "commit.json")
    public_data, public_record, public_fd = _open_bound(FAILED_V25_PUBLIC_KEY)
    if (
        len(signature_data) != 64
        or len(public_data) != 113
        or public_record["sha256"]
        != "8b5dc9f1715a3b8b0dafee138a66135f1a29a8bcc324ae15f28a0fd7f3ea54a6"
    ):
        raise RecoveryAuthorityError("Failed v25 authority key/signature identity drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    authority = _load_json(authority_data, "failed v25 signed authority")
    commit = _load_json(commit_data, "failed v25 authority commit")
    source_records = _records(FAILED_V25_SOURCE_PATHS)
    if {
        role: record["sha256"] for role, record in source_records.items()
    } != EXPECTED_FAILED_V25_SOURCE_SHA256:
        raise RecoveryAuthorityError("Failed v25 source identity drift")
    retired_private, _ = _open_retired_private_key_bound(FAILED_V25_PRIVATE_KEY)
    terminal_public_data, terminal_public_record, _ = _open_bound(
        FAILED_V25_TERMINAL_PUBLIC_KEY
    )
    terminal_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V25_TERMINAL_PRIVATE_KEY,
        expected_mode=0o400,
        state="FAILED_V25_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
    )
    expected_transaction_id = "8bda9ac7f4aae520426aa217c10fc86b"
    expected_members = [
        "authority.ed25519.sig",
        "authority.json",
        "commit.json",
    ]
    expected_retired_private_key = {
        "mode": "0o0",
        "path": evidence.get("key_archives", {})
        .get("authority", {})
        .get("private_key_pre_archive", {})
        .get("path"),
    }
    expected_commit = {
        "authority": {
            "mode": authority_record["mode"],
            "name": "authority.json",
            "sha256": authority_record["sha256"],
            "size_bytes": authority_record["size_bytes"],
        },
        "members": expected_members,
        "pass": True,
        "retired_private_key": expected_retired_private_key,
        "schema_name": "intersubmod.tumor_ref_schema_recovery_atomic_commit",
        "schema_version": "1.0.0",
        "signature": {
            "mode": signature_record["mode"],
            "name": "authority.ed25519.sig",
            "sha256": signature_record["sha256"],
            "size_bytes": signature_record["size_bytes"],
        },
        "transaction_id": expected_transaction_id,
    }
    expected_authority_bundle = {
        "members": expected_members,
        "path": str(FAILED_V25_ORIGINAL_OUTPUT_SLOTS[0]),
        "publication": "atomic_directory_rename_noreplace",
        "transaction_id": expected_transaction_id,
    }
    if (
        len(terminal_public_data) != 113
        or terminal_public_record["sha256"]
        != "b0056c3f60d7a7204d782ac1cea31e1e9411200b371295e38b9afc3d5f67a1d1"
        or authority.get("authority_id")
        != "20260723_tumor_ref_schema_recovery_v25"
        or authority.get("pass") is not True
        or not _strict_equal(commit, expected_commit)
        or commit_record["sha256"]
        != "51d1ec8e0555bf5ff62cc068db98fab84ccf7e6f13143bd671fac579c67c44ff"
        or not _strict_equal(authority.get("authority_bundle"), expected_authority_bundle)
        or not _strict_equal(
            authority.get("retired_private_key"), expected_retired_private_key
        )
        or evidence.get("generation") != "v25"
        or evidence.get("status")
        != "SIGNED_AUTHORITY_V_R_PASS_C_FAILED_AT_FINAL_DATASET_INTEGRATION"
        or evidence.get("root_cause", {}).get("scientific_payload_changed") is not False
        or evidence.get("root_cause", {}).get("claim_ceiling_changed") is not False
        or evidence.get("audit_lineage_observation", {}).get("artifact_set_sha256")
        != "195e77d571576f37debf1627edb57f9dc7edb935849bd0ae6e29b323b380b2ca"
        or evidence.get("audit_lineage_observation", {}).get("counts")
        != {
            "assignment_records": 102_842,
            "primary_artifacts_expected": 308_526,
            "primary_artifacts_verified": 308_526,
            "stable_sites": 102_842,
        }
        or evidence.get("terminal_output_state")
        != {
            "continuation_receipt_created": False,
            "continuation_signature_created": False,
            "continuation_success_witness_created": False,
            "dataset_release_receipt_created": False,
            "dataset_release_signature_created": False,
            "final_dataset_created": False,
            "final_report_created": False,
            "report_release_receipt_created": False,
            "report_release_signature_created": False,
            "supplemental_release_receipt_created": False,
        }
        or evidence.get("pass") is not False
        or retired_private.get("mode") != "0o0"
        or terminal_private_record.get("mode") != "0o400"
        or any(os.path.lexists(path) for path in FAILED_V25_ORIGINAL_OUTPUT_SLOTS)
    ):
        raise RecoveryAuthorityError("Failed v25 signed-generation contract drift")
    return {
        "failure_evidence": evidence_record,
        "authority": authority_record,
        "authority_signature": signature_record,
        "authority_commit": commit_record,
        "authority_public_key": public_record,
        "authority_retired_private_key": retired_private,
        "terminal_public_key": terminal_public_record,
        "terminal_unused_private_key": terminal_private_record,
        "sources": source_records,
        "failure": {
            "stage": "C_FINAL_DATASET_INTEGRATION",
            "contract": "tumor_ref_v1_declared_vs_v6_canonical_audit_path",
            "scientific_payload_changed": False,
            "terminal_signature_created": False,
        },
        "pass": True,
    }


def _validate_failed_v26_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V26_EVIDENCE)
    if evidence_record["sha256"] != EXPECTED_FAILED_V26_EVIDENCE_SHA256:
        raise RecoveryAuthorityError("Failed v26 evidence identity drift")
    evidence = _load_json(evidence_data, "failed v26 evidence")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V26_BUNDLE / "authority.json"
    )
    signature_data, signature_record, signature_fd = _open_bound(
        FAILED_V26_BUNDLE / "authority.ed25519.sig"
    )
    commit_data, commit_record, _ = _open_bound(FAILED_V26_BUNDLE / "commit.json")
    public_data, public_record, public_fd = _open_bound(FAILED_V26_PUBLIC_KEY)
    if (
        len(signature_data) != 64
        or len(public_data) != 113
        or public_record["sha256"]
        != "36fde24fb209c62465d74b4f42a70ca2361bab3c0e326cac25539e0f42049a41"
    ):
        raise RecoveryAuthorityError("Failed v26 authority key/signature identity drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    authority = _load_json(authority_data, "failed v26 signed authority")
    commit = _load_json(commit_data, "failed v26 authority commit")

    source_records = _records(FAILED_V26_SOURCE_PATHS)
    source_hashes = {
        role: record["sha256"] for role, record in source_records.items()
    }
    review_records = _records(FAILED_V26_REVIEW_PATHS)
    review_hashes = {
        role: record["sha256"] for role, record in review_records.items()
    }
    execution_records = _records(FAILED_V26_EXECUTION_PATHS)
    execution_hashes = {
        role: record["sha256"] for role, record in execution_records.items()
    }
    if (
        source_hashes != EXPECTED_FAILED_V26_SOURCE_SHA256
        or review_hashes != EXPECTED_FAILED_V26_REVIEW_SHA256
        or execution_hashes != EXPECTED_FAILED_V26_EXECUTION_SHA256
        or {
            role: record.get("sha256")
            for role, record in authority.get("recovery_sources", {}).items()
        }
        != EXPECTED_FAILED_V26_SOURCE_SHA256
        or {
            role: record.get("sha256")
            for role, record in authority.get("review_artifacts", {}).items()
        }
        != EXPECTED_FAILED_V26_REVIEW_SHA256
    ):
        raise RecoveryAuthorityError("Failed v26 source/review/execution identity drift")

    retired_private, _ = _open_retired_private_key_bound(FAILED_V26_PRIVATE_KEY)
    terminal_public_data, terminal_public_record, _ = _open_bound(
        FAILED_V26_TERMINAL_PUBLIC_KEY
    )
    terminal_private_record, _ = _open_private_key_metadata_bound(
        FAILED_V26_TERMINAL_PRIVATE_KEY,
        expected_mode=0o400,
        state="FAILED_V26_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
    )
    authority_archive_data, authority_archive_record, _ = _open_bound(
        FAILED_V26_AUTHORITY_KEY_ARCHIVE_RECORD
    )
    terminal_archive_data, terminal_archive_record, _ = _open_bound(
        FAILED_V26_TERMINAL_KEY_ARCHIVE_RECORD
    )
    authority_archive = _load_json(
        authority_archive_data, "failed v26 authority key archive record"
    )
    terminal_archive = _load_json(
        terminal_archive_data, "failed v26 terminal key archive record"
    )

    expected_transaction_id = "90766a0bc37bafd1ed19be014c4115d0"
    expected_members = [
        "authority.ed25519.sig",
        "authority.json",
        "commit.json",
    ]
    expected_retired_private_key = {
        "mode": "0o0",
        "path": (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260724_v26/ed25519_private_one_time.pem"
        ),
    }
    expected_commit = {
        "authority": {
            "mode": authority_record["mode"],
            "name": "authority.json",
            "sha256": authority_record["sha256"],
            "size_bytes": authority_record["size_bytes"],
        },
        "members": expected_members,
        "pass": True,
        "retired_private_key": expected_retired_private_key,
        "schema_name": "intersubmod.tumor_ref_schema_recovery_atomic_commit",
        "schema_version": "1.0.0",
        "signature": {
            "mode": signature_record["mode"],
            "name": "authority.ed25519.sig",
            "sha256": signature_record["sha256"],
            "size_bytes": signature_record["size_bytes"],
        },
        "transaction_id": expected_transaction_id,
    }
    expected_authority_bundle = {
        "members": expected_members,
        "path": str(FAILED_V26_ORIGINAL_OUTPUT_SLOTS[0]),
        "publication": "atomic_directory_rename_noreplace",
        "transaction_id": expected_transaction_id,
    }

    verification = _load_json(
        _open_bound(FAILED_V26_EXECUTION_PATHS["verification_receipt"])[0],
        "failed v26 verification receipt",
    )
    replay = _load_json(
        _open_bound(FAILED_V26_EXECUTION_PATHS["replay_receipt"])[0],
        "failed v26 replay receipt",
    )
    replay_witness = _load_json(
        _open_bound(FAILED_V26_EXECUTION_PATHS["replay_success_witness"])[0],
        "failed v26 replay success witness",
    )
    review_payloads = {
        role: _load_json(
            _open_bound(path)[0], f"failed v26 review {role}"
        )
        for role, path in FAILED_V26_REVIEW_PATHS.items()
    }
    root_cause = evidence.get("root_cause")
    impact = evidence.get("impact")
    if (
        authority_record["sha256"]
        != "26e632c9334a1438d754e03397760aa3c59624a12649762d9752eeb13c6f033e"
        or signature_record["sha256"]
        != "ce9078a9280b81b4ae89637394a6e25b23c2a15ab973ca90559817cfafc2b51a"
        or commit_record["sha256"]
        != "f581ceedbd5ed02221e3171f826c048c2b3f8811a2de2a76a87f01935c906de4"
        or not _strict_equal(commit, expected_commit)
        or authority.get("authority_id")
        != "20260724_tumor_ref_schema_recovery_v26"
        or authority.get("pass") is not True
        or not _strict_equal(authority.get("authority_bundle"), expected_authority_bundle)
        or not _strict_equal(
            authority.get("retired_private_key"), expected_retired_private_key
        )
        or evidence.get("generation") != "v26"
        or evidence.get("failed_stage")
        != "C_pre_downstream_validate_promotion_chain"
        or evidence.get("error")
        != {
            "message": "Signed promotion authorization exact contract drift",
            "type": "ContinuationError",
        }
        or root_cause
        != {
            "current_only_runtime_roles": [
                "recovery_final_dataset_builder",
                "recovery_report_builder",
                "recovery_result_report_finalizer",
            ],
            "current_runtime_role_count": 14,
            "signed_authorization_runtime_role_count": 11,
            "signed_only_runtime_roles": [],
            "type": (
                "historical_signed_projection_contaminated_by_"
                "recovery_only_runtime_roles"
            ),
        }
        or not isinstance(impact, Mapping)
        or impact.get("downstream_producers_started") is not False
        or impact.get("scientific_payload_changed") is not False
        or impact.get("confirmed_cellular_subclone_count_changed") is not False
        or impact.get("linear_ancestry_count_changed") is not False
        or evidence.get("pass") is not False
        or verification.get("pass") is not True
        or replay.get("pass") is not True
        or replay_witness.get("pass") is not True
        or any(
            payload.get("verdict") != "APPROVE"
            or payload.get("pass") is not True
            or payload.get("high_findings") != []
            or payload.get("medium_findings") != []
            or payload.get("unresolved_conditions") != []
            for payload in review_payloads.values()
        )
        or authority_archive_record["sha256"]
        != "a3957acd555181a4a4b575db7dc6f0190db6cb882f020f2f30a445aa6eb33372"
        or terminal_archive_record["sha256"]
        != "ba9b0d191fb49aa1c5fe897009fd91ce1592c119b96133e70880354f9f73477e"
        or authority_archive.get("generation") != "v26"
        or authority_archive.get("key_reuse_forbidden") is not True
        or authority_archive.get("private_key_bytes_read") is not False
        or authority_archive.get("state")
        != "RETIRED_AFTER_SINGLE_SIGNED_AUTHORITY_THAT_FAILED_AT_C"
        or terminal_archive.get("generation") != "terminal_v16"
        or terminal_archive.get("key_reuse_forbidden") is not True
        or terminal_archive.get("private_key_bytes_read") is not False
        or terminal_archive.get("state")
        != "UNUSED_NO_TERMINAL_SIGNATURE_CREATED_RETIRED_FROM_REUSE"
        or retired_private.get("mode") != "0o0"
        or len(terminal_public_data) != 113
        or terminal_public_record["sha256"]
        != "b61e7f75ba3e418098c3c7f9fb19da0b261380d7985833eb6a59b99d6e1aeaee"
        or terminal_private_record.get("mode") != "0o400"
        or any(os.path.lexists(path) for path in FAILED_V26_ORIGINAL_OUTPUT_SLOTS)
    ):
        raise RecoveryAuthorityError("Failed v26 signed-generation contract drift")
    return {
        "failure_evidence": evidence_record,
        "authority": authority_record,
        "authority_signature": signature_record,
        "authority_commit": commit_record,
        "authority_public_key": public_record,
        "authority_retired_private_key": retired_private,
        "authority_key_archive_record": authority_archive_record,
        "terminal_public_key": terminal_public_record,
        "terminal_unused_private_key": terminal_private_record,
        "terminal_key_archive_record": terminal_archive_record,
        "sources": source_records,
        "reviews": review_records,
        "verification_and_replay": execution_records,
        "failure": {
            "stage": "C_PRE_DOWNSTREAM_VALIDATE_PROMOTION_CHAIN",
            "contract": "historical_signed_runtime_projection_exact_11_roles",
            "current_reviewed_runtime_role_count": 14,
            "scientific_payload_changed": False,
            "terminal_signature_created": False,
        },
        "pass": True,
    }


def _validate_failed_v28_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V28_EVIDENCE)
    if (
        evidence_record["size_bytes"] != 588_302
        or evidence_record["sha256"] != EXPECTED_FAILED_V28_EVIDENCE_SHA256
    ):
        raise RecoveryAuthorityError("Failed v28 evidence identity drift")
    evidence = _load_json(evidence_data, "failed v28 evidence")
    summary_record = _open_bound(FAILED_V28_SUMMARY)[1]
    if summary_record["sha256"] != (
        "a4b98ed12d2f7566e8e9eed293d0090a3fb9d269bd986f02e87670053c7453c8"
    ):
        raise RecoveryAuthorityError("Failed v28 summary identity drift")

    inventory = evidence.get("archived_inventory")
    if not isinstance(inventory, list) or len(inventory) != 742:
        raise RecoveryAuthorityError("Failed v28 archived inventory cardinality drift")
    inventory_relative_paths: set[str] = set()
    archive_root = FAILED_V28_ROOT.resolve(strict=True)
    for expected in inventory:
        _exact_keys(
            expected,
            {"mode", "path", "relative_path", "sha256", "size_bytes"},
            "failed v28 archived inventory row",
        )
        relative_path = str(expected["relative_path"])
        if relative_path in inventory_relative_paths:
            raise RecoveryAuthorityError("Failed v28 archived inventory duplicate path")
        inventory_relative_paths.add(relative_path)
        archived_path = (FAILED_V28_ROOT / relative_path).resolve(strict=True)
        if not archived_path.is_relative_to(archive_root):
            raise RecoveryAuthorityError("Failed v28 archived inventory path escape")
        observed = _open_bound(archived_path)[1]
        if {
            "mode": observed["mode"],
            "sha256": observed["sha256"],
            "size_bytes": observed["size_bytes"],
        } != {
            "mode": expected["mode"],
            "sha256": expected["sha256"],
            "size_bytes": expected["size_bytes"],
        }:
            raise RecoveryAuthorityError(
                f"Failed v28 archived artifact identity drift: {relative_path}"
            )

    source_records = _records(FAILED_V28_SOURCE_PATHS)
    if not _strict_equal(
        evidence.get("source_artifacts_left_immutable_at_original_paths"),
        source_records,
    ):
        raise RecoveryAuthorityError("Failed v28 frozen source identities drift")

    authority_data, authority_record, authority_fd = _open_bound(
        FAILED_V28_BUNDLE / "authority.json"
    )
    authority_signature_data, authority_signature_record, authority_signature_fd = (
        _open_bound(FAILED_V28_BUNDLE / "authority.ed25519.sig")
    )
    authority_commit_data, authority_commit_record, _ = _open_bound(
        FAILED_V28_BUNDLE / "commit.json"
    )
    authority_public_data, authority_public_record, authority_public_fd = _open_bound(
        FAILED_V28_KEY_ARCHIVES["authority_v28"]["public"]
    )
    if len(authority_signature_data) != 64 or len(authority_public_data) != 113:
        raise RecoveryAuthorityError("Failed v28 authority signature/key size drift")
    _verify_signature(authority_fd, authority_public_fd, authority_signature_fd)
    authority = _load_json(authority_data, "failed v28 signed authority")
    authority_commit = _load_json(authority_commit_data, "failed v28 authority commit")

    review_records = _records(FAILED_V28_REVIEW_PATHS)
    review_payloads = {
        role: _load_json(_open_bound(path)[0], f"failed v28 review {role}")
        for role, path in FAILED_V28_REVIEW_PATHS.items()
    }
    execution_records = _records(FAILED_V28_EXECUTION_PATHS)
    execution_payloads = {
        role: _load_json(_open_bound(path)[0], f"failed v28 execution {role}")
        for role, path in FAILED_V28_EXECUTION_PATHS.items()
        if path.suffix == ".json"
    }
    provisional_records = _records(FAILED_V28_PROVISIONAL_DATASET_PATHS)
    final_dataset = _load_archived_transport_json(
        _open_bound(FAILED_V28_PROVISIONAL_DATASET_PATHS["final_dataset"])[0],
        "failed v28 provisional final dataset",
    )
    dataset_release_data, dataset_release_record, dataset_release_fd = _open_bound(
        FAILED_V28_PROVISIONAL_DATASET_PATHS["dataset_release_receipt"]
    )
    dataset_signature_data, dataset_signature_record, dataset_signature_fd = _open_bound(
        FAILED_V28_PROVISIONAL_DATASET_PATHS["dataset_release_signature"]
    )
    result_public_data, result_public_record, result_public_fd = _open_bound(
        FAILED_V28_KEY_ARCHIVES["result_v5"]["public"]
    )
    if len(dataset_signature_data) != 64 or len(result_public_data) != 113:
        raise RecoveryAuthorityError("Failed v28 dataset signature/key size drift")
    _verify_signature(dataset_release_fd, result_public_fd, dataset_signature_fd)
    dataset_release = _load_json(
        dataset_release_data, "failed v28 provisional dataset release"
    )

    key_archive_records: dict[str, Any] = {}
    for role, spec in FAILED_V28_KEY_ARCHIVES.items():
        public_data, public_record, _ = _open_bound(spec["public"])
        archive_record_data, archive_record, _ = _open_bound(
            spec["root"] / "FAILED_KEY_ARCHIVE_RECORD.v1.json"
        )
        archive_payload = _load_json(
            archive_record_data, f"failed v28 {role} key archive record"
        )
        private_mode = int(spec["private_mode"])
        if private_mode == 0:
            private_record, _ = _open_retired_private_key_bound(spec["private"])
        else:
            private_record, _ = _open_private_key_metadata_bound(
                spec["private"],
                expected_mode=private_mode,
                state=f"FAILED_V28_{role.upper()}_ARCHIVED_MUST_NOT_BE_REUSED",
            )
        if (
            len(public_data) != 113
            or public_record["sha256"] != spec["public_sha256"]
            or archive_record["sha256"] != spec["record_sha256"]
            or archive_payload.get("generation") != "v28"
            or archive_payload.get("role") != role
            or archive_payload.get("state") != spec["state"]
            or archive_payload.get("key_reuse_forbidden") is not True
            or archive_payload.get("private_key_bytes_read") is not False
            or archive_payload.get("archive_directory") != str(spec["root"])
            or os.path.lexists(archive_payload.get("source_directory", ""))
        ):
            raise RecoveryAuthorityError(f"Failed v28 {role} key archive drift")
        key_archive_records[role] = {
            "public_key": public_record,
            "private_key": private_record,
            "archive_record": archive_record,
            "state": spec["state"],
        }

    claim_counts = final_dataset.get("counts", {}).get("claim_status_counts", {})
    expected_pass_counts = {"M1": 102_842, "M2": 919, "G1": 7, "G2": 0}
    if (
        authority_record["sha256"]
        != "01e75a6eb54b0510001abbd72907da1b269db4e80d5842e9f5b35163d18454b6"
        or authority_signature_record["sha256"]
        != "069dd80d3c97bf7e940ef9d9d85490e993cfdad1c3cde632f7e76b9302856ce5"
        or authority_commit_record["sha256"]
        != "a307a4dabac00457aca3c607eecacf7a2091eb484a6db5d562987b0985b5392e"
        or authority.get("authority_id")
        != "20260724_tumor_ref_schema_recovery_v28"
        or authority.get("pass") is not True
        or authority_commit.get("pass") is not True
        or any(
            payload.get("verdict") != "APPROVE"
            or payload.get("pass") is not True
            or payload.get("high_findings") != []
            or payload.get("medium_findings") != []
            or payload.get("unresolved_conditions") != []
            for payload in review_payloads.values()
        )
        or execution_payloads["verification_receipt"].get("pass") is not True
        or execution_payloads["replay_receipt"].get("pass") is not True
        or execution_payloads["replay_success_witness"].get("pass") is not True
        or execution_payloads["continuation_incident"].get("pass") is not False
        or execution_payloads["continuation_incident"].get("error", {}).get("type")
        != "ContinuationError"
        or final_dataset.get("schema_name")
        != "intersubmod.all_ssnv_final_report_dataset"
        or final_dataset.get("task_type") != "B_comprehensive_validation"
        or final_dataset.get("scope", {}).get("observed_screen_sites") != 469_849
        or final_dataset.get("pass") is not True
        or {
            claim_id: claim_counts.get(claim_id, {}).get("PASS")
            for claim_id in expected_pass_counts
        }
        != expected_pass_counts
        or dataset_release_record["sha256"]
        != "e2e6896aee1b6e82756a9dabd8e743dc86df7281720990d9a98a7a7496466543"
        or dataset_signature_record["sha256"]
        != "4f0f2c7c910aa32fc7fc78d52f576a5c4dbea36b88584f87ca7a7cc0730b8425"
        or dataset_release.get("schema_name")
        != "intersubmod.task_b_final_dataset_release_receipt"
        or dataset_release.get("pass") is not True
        or evidence.get("generation") != "v28"
        or evidence.get("status")
        != "SIGNED_AUTHORITY_V_R_PASS_C_FAILED_AFTER_SIGNED_DATASET_BEFORE_REPORT"
        or evidence.get("root_cause", {}).get("error_fragment")
        != "ReportContractError: Final per_sample metric strata drift"
        or evidence.get("root_cause", {}).get("scientific_payload_changed") is not False
        or evidence.get("root_cause", {}).get("claim_ceiling_changed") is not False
        or evidence.get("verified_pre_failure_state", {})
        .get("dataset_signature", {})
        .get("verified")
        is not True
        or evidence.get("terminal_output_state", {}).get("final_report_created") is not False
        or evidence.get("terminal_output_state", {}).get(
            "continuation_receipt_created"
        )
        is not False
        or evidence.get("pass") is not False
        or any(os.path.lexists(path) for path in FAILED_V28_ORIGINAL_OUTPUT_SLOTS)
    ):
        raise RecoveryAuthorityError("Failed v28 signed-generation contract drift")

    return {
        "failure_evidence": evidence_record,
        "summary": summary_record,
        "archived_inventory_count": len(inventory),
        "archived_inventory_sha256": _json_sha256(inventory),
        "authority": authority_record,
        "authority_signature": authority_signature_record,
        "authority_commit": authority_commit_record,
        "authority_public_key": authority_public_record,
        "authority_signature_verified": True,
        "sources": source_records,
        "reviews": review_records,
        "verification_replay_and_incident": execution_records,
        "provisional_dataset": provisional_records,
        "dataset_release_signature_verified": True,
        "key_archives": key_archive_records,
        "failure": {
            "stage": "C_REPORT_BUILD_AFTER_SIGNED_DATASET",
            "contract": "json_mapping_exact_key_set_not_insertion_order",
            "scientific_payload_changed": False,
            "terminal_signature_created": False,
            "report_signature_created": False,
        },
        "pass": True,
    }


def _validate_failed_v29_signed_recovery() -> dict[str, Any]:
    evidence_data, evidence_record, _ = _open_bound(FAILED_V29_EVIDENCE)
    summary_data, summary_record, _ = _open_bound(FAILED_V29_SUMMARY)
    if (
        evidence_record["size_bytes"] != 91_001
        or evidence_record["sha256"] != EXPECTED_FAILED_V29_EVIDENCE_SHA256
        or summary_record["size_bytes"] != 1_233
        or summary_record["sha256"] != EXPECTED_FAILED_V29_SUMMARY_SHA256
        or not summary_data.startswith(b"<!--")
    ):
        raise RecoveryAuthorityError("Failed v29 failure evidence identity drift")
    evidence = _load_json(evidence_data, "failed v29 failure evidence")

    archived = evidence.get("archived_artifacts")
    sources = evidence.get("source_artifacts_left_immutable_at_original_paths")
    if not isinstance(archived, Mapping) or not isinstance(sources, Mapping):
        raise RecoveryAuthorityError("Failed v29 archived/source artifact inventory missing")
    archived_records: dict[str, Any] = {}
    for role, expected in archived.items():
        if not isinstance(expected, Mapping) or expected.get("link_count") != 1:
            raise RecoveryAuthorityError(f"Failed v29 archived record drift: {role}")
        observed = _open_bound(Path(str(expected.get("path"))))[1]
        if any(
            observed.get(key) != expected.get(key)
            for key in ("path", "size_bytes", "sha256", "mode")
        ):
            raise RecoveryAuthorityError(f"Failed v29 archived identity drift: {role}")
        archived_records[str(role)] = observed
    source_records: dict[str, Any] = {}
    for role, expected in sources.items():
        if not isinstance(expected, Mapping) or expected.get("link_count") != 1:
            raise RecoveryAuthorityError(f"Failed v29 source record drift: {role}")
        observed = _open_bound(Path(str(expected.get("path"))))[1]
        if any(
            observed.get(key) != expected.get(key)
            for key in ("path", "size_bytes", "sha256", "mode")
        ):
            raise RecoveryAuthorityError(f"Failed v29 source identity drift: {role}")
        source_records[str(role)] = observed

    authority_fd = _open_bound(
        Path(archived["authority"]["path"])
    )[2]
    signature_data, _, signature_fd = _open_bound(
        Path(archived["authority_signature"]["path"])
    )
    authority_public_data, authority_public_record, authority_public_fd = _open_bound(
        FAILED_V29_KEY_ARCHIVES["authority_v29"]["root"] / "ed25519_public.pem"
    )
    _verify_signature(authority_fd, authority_public_fd, signature_fd)
    authority = _load_json(
        _open_bound(Path(archived["authority"]["path"]))[0],
        "failed v29 signed authority",
    )
    authority_commit = _load_json(
        _open_bound(Path(archived["authority_commit"]["path"]))[0],
        "failed v29 authority commit",
    )
    if (
        len(signature_data) != 64
        or len(authority_public_data) != 113
        or authority_public_record["sha256"]
        != FAILED_V29_KEY_ARCHIVES["authority_v29"]["public_sha256"]
        or authority.get("authority_id") != "20260724_tumor_ref_schema_recovery_v29"
        or authority.get("pass") is not True
        or authority_commit.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Failed v29 signed authority contract drift")

    key_archive_records: dict[str, Any] = {}
    for role, spec in FAILED_V29_KEY_ARCHIVES.items():
        public_data, public_record, _ = _open_bound(spec["root"] / "ed25519_public.pem")
        archive_data, archive_record, _ = _open_bound(
            spec["root"] / "FAILED_KEY_ARCHIVE_RECORD.v1.json"
        )
        private_path = spec["root"] / spec["private_name"]
        if spec["private_mode"] == 0:
            private_record, _ = _open_retired_private_key_bound(private_path)
        else:
            private_record, _ = _open_private_key_metadata_bound(
                private_path,
                expected_mode=spec["private_mode"],
                state=f"FAILED_V29_{role.upper()}_ARCHIVED_MUST_NOT_BE_REUSED",
            )
        archive_payload = _load_json(archive_data, f"failed v29 {role} key archive")
        if (
            len(public_data) != 113
            or public_record["sha256"] != spec["public_sha256"]
            or archive_record["sha256"] != spec["record_sha256"]
            or archive_payload.get("generation") != "v29"
            or archive_payload.get("role") != role
            or archive_payload.get("state") != spec["state"]
            or archive_payload.get("key_reuse_forbidden") is not True
            or archive_payload.get("private_key_bytes_read") is not False
            or archive_payload.get("archive_directory") != str(spec["root"])
            or os.path.lexists(archive_payload.get("source_directory", ""))
        ):
            raise RecoveryAuthorityError(f"Failed v29 {role} key archive drift")
        key_archive_records[role] = {
            "public_key": public_record,
            "private_key": private_record,
            "archive_record": archive_record,
            "state": spec["state"],
        }

    review_summary = evidence.get("reviews")
    review_contract = evidence.get("review_contract")
    for role, artifact_role in (
        ("mendel", "review_mendel"),
        ("nash", "review_nash"),
        ("external_claude_opus", "review_external_claude_opus"),
    ):
        review = _load_json(
            _open_bound(Path(archived[artifact_role]["path"]))[0],
            f"failed v29 review {role}",
        )
        if (
            not isinstance(review_summary, Mapping)
            or review.get("verdict") != "APPROVE"
            or review.get("pass") is not True
            or review.get("reviewed_source_set_sha256")
            != review_contract.get("reviewed_source_set_sha256")
            or review_summary.get(role, {}).get("reviewer_agent_id")
            != review.get("reviewer_agent_id")
        ):
            raise RecoveryAuthorityError(f"Failed v29 review contract drift: {role}")

    verification = _load_json(
        _open_bound(Path(archived["verification_receipt"]["path"]))[0],
        "failed v29 verification receipt",
    )
    runtime = evidence.get("runtime")
    terminal_state = evidence.get("terminal_output_state")
    root_cause = evidence.get("root_cause")
    if (
        evidence.get("schema_name")
        != "intersubmod.tumor_ref_schema_recovery_signed_generation_failure"
        or evidence.get("generation") != "v29"
        or evidence.get("status")
        != "SIGNED_AUTHORITY_AND_V_PASS_R_FAILED_BEFORE_FIRST_R_OUTPUT"
        or evidence.get("scientific_payload_changed") is not False
        or evidence.get("claim_ceiling_changed") is not False
        or evidence.get("pass") is not False
        or verification.get("pass") is not True
        or runtime.get("v_exit_code") != 0
        or runtime.get("r_worker_exit_code") != 70
        or runtime.get("r_supervisor_exit_code") != 1
        or runtime.get("runner_prefix_stderr_sha256")
        != "18477fe5a82535b402f6f04e72ebcea147f5c5cee0475dba6e250ae75078c4a0"
        or terminal_state
        != {
            "c_started": False,
            "dataset_release_signature_created": False,
            "final_dataset_created": False,
            "final_report_created": False,
            "r_log_created": False,
            "r_receipt_created": False,
            "r_success_witness_created": False,
            "report_release_signature_created": False,
            "terminal_signature_created": False,
        }
        or not isinstance(root_cause, Mapping)
        or root_cause.get("scientific_payload_changed") is not False
        or root_cause.get("claim_ceiling_changed") is not False
        or evidence.get("rejected_remediation", {}).get("temporary_live_key_projection")
        is not False
        or any(os.path.lexists(path) for path in FAILED_V29_ORIGINAL_OUTPUT_SLOTS)
    ):
        raise RecoveryAuthorityError("Failed v29 signed-generation contract drift")
    return {
        "failure_evidence": evidence_record,
        "summary": summary_record,
        "archived_artifacts": archived_records,
        "sources": source_records,
        "authority_public_key": authority_public_record,
        "authority_signature_verified": True,
        "key_archives": key_archive_records,
        "failure": {
            "stage": "R_RUNNER_PREFIX_BEFORE_FIRST_R_OUTPUT",
            "root_cause": "historical_live_key_paths_absent_after_required_quarantine",
            "scientific_payload_changed": False,
            "claim_ceiling_changed": False,
        },
        "pass": True,
    }


def _validate_v30_key_bootstrap() -> dict[str, Any]:
    def bind_expected_record(expected: Any, label: str) -> dict[str, Any]:
        if not isinstance(expected, Mapping):
            raise RecoveryAuthorityError(f"V30 bootstrap {label} record is missing")
        expected_keys = set(expected)
        allowed_keys = {
            "path",
            "size_bytes",
            "sha256",
            "mode",
            "link_count",
            "device",
            "inode",
        }
        if not expected_keys <= allowed_keys or not {
            "path",
            "size_bytes",
            "sha256",
            "mode",
        } <= expected_keys:
            raise RecoveryAuthorityError(
                f"V30 bootstrap {label} record schema drift"
            )
        data, basic, descriptor = _open_bound(
            Path(expected["path"]), expected_mode=expected["mode"]
        )
        opened = os.fstat(descriptor)
        observed_all = {
            **basic,
            "link_count": opened.st_nlink,
            "device": opened.st_dev,
            "inode": opened.st_ino,
        }
        observed = {key: observed_all[key] for key in expected}
        if len(data) != expected["size_bytes"] or not _strict_equal(
            observed, expected
        ):
            raise RecoveryAuthorityError(
                f"V30 bootstrap {label} bound identity drift"
            )
        return observed

    def bind_expected_private_metadata(expected: Any, label: str) -> dict[str, Any]:
        if not isinstance(expected, Mapping):
            raise RecoveryAuthorityError(
                f"V30 bootstrap {label} private metadata is missing"
            )
        expected_keys = set(expected)
        allowed_keys = {
            "path",
            "size_bytes",
            "mode",
            "link_count",
            "device",
            "inode",
        }
        if not expected_keys <= allowed_keys or not {
            "path",
            "size_bytes",
            "mode",
            "link_count",
        } <= expected_keys:
            raise RecoveryAuthorityError(
                f"V30 bootstrap {label} private metadata schema drift"
            )
        path = Path(expected["path"])
        if not path.is_absolute() or path.resolve(strict=True) != path:
            raise RecoveryAuthorityError(
                f"V30 bootstrap {label} private path is not canonical"
            )
        before = os.stat(path, follow_symlinks=False)
        expected_mode = int(expected["mode"], 8)
        if (
            not stat.S_ISREG(before.st_mode)
            or stat.S_IMODE(before.st_mode) != expected_mode
            or before.st_nlink != 1
            or not hasattr(os, "O_PATH")
        ):
            raise RecoveryAuthorityError(
                f"V30 bootstrap {label} private metadata state drift"
            )
        flags = os.O_PATH | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(path, flags)
        opened = os.fstat(descriptor)
        after = os.stat(path, follow_symlinks=False)
        identity = lambda value: (
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_nlink,
            value.st_size,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
        if identity(before) != identity(opened) or identity(opened) != identity(after):
            os.close(descriptor)
            raise RecoveryAuthorityError(
                f"V30 bootstrap {label} private metadata changed while binding"
            )
        observed_all = {
            "path": str(path),
            "size_bytes": opened.st_size,
            "mode": oct(stat.S_IMODE(opened.st_mode)),
            "link_count": opened.st_nlink,
            "device": opened.st_dev,
            "inode": opened.st_ino,
        }
        observed = {key: observed_all[key] for key in expected}
        if not _strict_equal(observed, expected):
            os.close(descriptor)
            raise RecoveryAuthorityError(
                f"V30 bootstrap {label} private metadata drift"
            )
        _FD_LEASES.append((path, descriptor, opened, None))
        return observed

    records: dict[str, Any] = {}
    payloads: dict[str, Any] = {}
    for role, path in (
        ("prepare", V30_BOOTSTRAP_PREPARE),
        ("receipt", V30_BOOTSTRAP_RECEIPT),
        ("success", V30_BOOTSTRAP_SUCCESS),
    ):
        data, record, _ = _open_bound(path)
        if record["sha256"] != EXPECTED_V30_BOOTSTRAP_SHA256[role]:
            raise RecoveryAuthorityError(f"V30 bootstrap {role} identity drift")
        records[role] = record
        payloads[role] = _load_json(data, f"v30 bootstrap {role}")

    progress_data, progress_record, _ = _open_bound(V30_BOOTSTRAP_PROGRESS)
    if progress_record["sha256"] != EXPECTED_V30_BOOTSTRAP_SHA256["progress"]:
        raise RecoveryAuthorityError("V30 bootstrap progress identity drift")
    progress_events = [
        _load_json(line, f"v30 bootstrap progress line {index}")
        for index, line in enumerate(progress_data.splitlines(), start=1)
    ]

    public_records: dict[str, Any] = {}
    for role, path, expected_sha256 in (
        ("authority", PUBLIC_KEY, EXPECTED_PUBLIC_KEY_SHA256),
        (
            "terminal",
            RECOVERY_CONTINUATION_PUBLIC_KEY,
            EXPECTED_RECOVERY_CONTINUATION_PUBLIC_KEY_SHA256,
        ),
        ("result", V30_RESULT_PUBLIC_KEY, EXPECTED_V30_RESULT_PUBLIC_KEY_SHA256),
        ("report", V30_REPORT_PUBLIC_KEY, EXPECTED_V30_REPORT_PUBLIC_KEY_SHA256),
    ):
        data, record, descriptor = _open_bound(path)
        opened = os.fstat(descriptor)
        if len(data) != 113 or record["sha256"] != expected_sha256:
            raise RecoveryAuthorityError(f"V30 bootstrap {role} public key drift")
        public_records[role] = {**record, "link_count": opened.st_nlink}
    prepare = payloads["prepare"]
    receipt = payloads["receipt"]
    success = payloads["success"]
    transaction_id = prepare.get("transaction_id")
    expected_progress_events = [
        "TRANSACTION_STARTED",
        "SIGNAL_FENCE_ENTERED",
        "AUTHORITY_STAGING_KEY_CREATED",
        "TERMINAL_STAGING_KEY_CREATED",
        "AUTHORITY_KEY_ROOT_PUBLISHED",
        "TERMINAL_KEY_ROOT_PUBLISHED",
        "BOOTSTRAP_RECEIPT_PUBLISHED_AWAITING_SUCCESS_WITNESS",
    ]
    prepare_record = {**records["prepare"], "link_count": 1}
    progress_record_with_link = {**progress_record, "link_count": 1}
    receipt_record = {**records["receipt"], "link_count": 1}
    pre_authority_archive = prepare.get("pre_authority_archive")
    failed_v8_archive = prepare.get("failed_v8_pre_ready_archive")
    if (
        prepare.get("schema_name")
        != "intersubmod.v30_four_role_key_bootstrap_prepare"
        or prepare.get("generation") != "v30"
        or prepare.get("pass") is not False
        or prepare.get("release_authority_granted") is not False
        or receipt.get("schema_name") != "intersubmod.v30_four_role_key_bootstrap"
        or receipt.get("generation") != "v30"
        or receipt.get("pass") is not True
        or receipt.get("release_authority_granted") is not False
        or success.get("schema_name")
        != "intersubmod.v30_four_role_key_bootstrap_success"
        or success.get("generation") != "v30"
        or success.get("pass") is not True
        or success.get("release_authority_granted") is not False
        or not isinstance(transaction_id, str)
        or len(transaction_id) != 32
        or transaction_id != receipt.get("transaction_id")
        or transaction_id != success.get("transaction_id")
        or len(progress_events) != len(expected_progress_events)
        or [event.get("event") for event in progress_events]
        != expected_progress_events
        or any(event.get("transaction_id") != transaction_id for event in progress_events)
        or receipt.get("four_role_public_keys_pairwise_distinct") is not True
        or len({record["sha256"] for record in public_records.values()}) != 4
        or receipt.get("authority", {}).get("public_key") != public_records["authority"]
        or receipt.get("terminal", {}).get("public_key") != public_records["terminal"]
        or receipt.get("result", {}).get("public_key") != public_records["result"]
        or receipt.get("report", {}).get("public_key") != public_records["report"]
        or receipt.get("prepare_receipt") != prepare_record
        or success.get("authority_public_key") != public_records["authority"]
        or success.get("terminal_public_key") != public_records["terminal"]
        or success.get("result", {}).get("public_key") != public_records["result"]
        or success.get("report", {}).get("public_key") != public_records["report"]
        or success.get("prepare_receipt") != prepare_record
        or success.get("progress") != progress_record_with_link
        or success.get("receipt") != receipt_record
        or not _strict_equal(receipt.get("result"), success.get("result"))
        or not _strict_equal(receipt.get("report"), success.get("report"))
        or not _strict_equal(pre_authority_archive, receipt.get("pre_authority_archive"))
        or not _strict_equal(pre_authority_archive, success.get("pre_authority_archive"))
        or not _strict_equal(failed_v8_archive, receipt.get("failed_v8_pre_ready_archive"))
        or not _strict_equal(failed_v8_archive, success.get("failed_v8_pre_ready_archive"))
        or success.get("checks", {}).get("four_role_public_keys_pairwise_distinct")
        is not True
        or success.get("checks", {}).get("pre_authority_four_role_archive_verified")
        is not True
        or success.get("checks", {}).get("failed_v8_pre_ready_archive_verified")
        is not True
        or success.get("checks", {}).get("progress_ledger_fsynced_before_witness")
        is not True
        or not isinstance(pre_authority_archive, Mapping)
        or pre_authority_archive.get("pass") is not True
        or pre_authority_archive.get("all_original_roots_absent") is not True
        or pre_authority_archive.get("key_reuse_forbidden") is not True
        or pre_authority_archive.get("partial_witness_absent") is not True
        or not isinstance(failed_v8_archive, Mapping)
        or failed_v8_archive.get("pass") is not True
        or failed_v8_archive.get("source_root_absent") is not True
        or failed_v8_archive.get("key_reuse_forbidden") is not True
        or failed_v8_archive.get("partial_witness_absent") is not True
        or receipt.get("private_key_bytes_recorded_in_receipt") is not False
    ):
        raise RecoveryAuthorityError("V30 four-key bootstrap contract drift")

    for role in ("authority", "terminal"):
        role_payload = receipt.get(role)
        if (
            not isinstance(role_payload, Mapping)
            or role_payload.get("private_key_bytes_recorded_in_receipt") is not False
        ):
            raise RecoveryAuthorityError(
                f"V30 bootstrap {role} private-key receipt drift"
            )
        bind_expected_private_metadata(
            role_payload.get("private_key_metadata_only"), f"{role} private key"
        )
    for role in ("result", "report"):
        role_payload = receipt.get(role)
        if (
            not isinstance(role_payload, Mapping)
            or role_payload.get("pass") is not True
            or role_payload.get("private_key_bytes_read") is not False
        ):
            raise RecoveryAuthorityError(
                f"V30 bootstrap {role} signer receipt drift"
            )
        bind_expected_private_metadata(
            role_payload.get("private_key_metadata_only"), f"{role} private key"
        )
        for nested_role in (
            "prepare_receipt",
            "ready_receipt",
            "signer_source",
            "wrapper_source",
        ):
            bind_expected_record(
                role_payload.get(nested_role), f"{role} {nested_role}"
            )

    for nested_role in ("prepare", "progress", "receipt"):
        bind_expected_record(
            pre_authority_archive.get(nested_role),
            f"pre-authority archive {nested_role}",
        )
    archives = pre_authority_archive.get("archives")
    if not isinstance(archives, Mapping) or set(archives) != {
        "authority_v30_initial",
        "terminal_v20_initial",
        "result_v7_initial",
        "report_v7_initial",
    }:
        raise RecoveryAuthorityError("V30 pre-authority archive role drift")
    for role, archive in archives.items():
        if (
            not isinstance(archive, Mapping)
            or archive.get("pass") is not True
            or archive.get("key_reuse_forbidden") is not True
        ):
            raise RecoveryAuthorityError(
                f"V30 pre-authority archive {role} contract drift"
            )
        bind_expected_record(
            archive.get("archive_record"), f"pre-authority {role} archive record"
        )
        bind_expected_record(
            archive.get("public_key"), f"pre-authority {role} public key"
        )
        bind_expected_private_metadata(
            archive.get("private_key_metadata_only"),
            f"pre-authority {role} private key",
        )

    for nested_role in (
        "prepare",
        "progress",
        "receipt",
        "failure_record",
        "previous_prepare",
        "public_key",
    ):
        bind_expected_record(
            failed_v8_archive.get(nested_role), f"failed-v8 archive {nested_role}"
        )
    bind_expected_private_metadata(
        failed_v8_archive.get("private_key_metadata_only"),
        "failed-v8 archive private key",
    )
    return {
        "prepare": records["prepare"],
        "progress": progress_record,
        "receipt": records["receipt"],
        "success": records["success"],
        "public_keys": public_records,
        "pre_authority_archive": pre_authority_archive,
        "failed_v8_pre_ready_archive": failed_v8_archive,
        "four_role_public_keys_pairwise_distinct": True,
        "release_authority_granted": False,
        "pass": True,
    }


def _validate_prior_failed_signed_recoveries() -> dict[str, Any]:
    return {
        "v9": _validate_failed_v9_signed_recovery(),
        "v10": _validate_failed_v10_signed_recovery(),
        "v11": _validate_failed_v11_signed_recovery(),
        "v12": _validate_failed_v12_signed_recovery(),
        "v14": _validate_failed_v14_signed_recovery(),
        "v16": _validate_failed_v16_signed_recovery(),
        "v17": _validate_failed_v17_signed_recovery(),
        "v18": _validate_failed_v18_signed_recovery(),
        "v19": _validate_failed_v19_signed_recovery(),
        "v21": _validate_failed_v21_signed_recovery(),
        "v25": _validate_failed_v25_signed_recovery(),
        "v26": _validate_failed_v26_signed_recovery(),
        "v28": _validate_failed_v28_signed_recovery(),
        "v29": _validate_failed_v29_signed_recovery(),
    }


def _is_canonical_uuid(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 36
        and value[8] == value[13] == value[18] == value[23] == "-"
        and all(
            character in "0123456789abcdef"
            for index, character in enumerate(value)
            if index not in {8, 13, 18, 23}
        )
    )


def _validate_review(
    role: str,
    payload: Mapping[str, Any],
    recovery_sources: Mapping[str, Any],
    legacy_sources: Mapping[str, Any],
    prior_recovery_chain: Mapping[str, Any],
    rejected_generations: Mapping[str, Any],
    prior_failed_signed_recovery: Mapping[str, Any],
    fresh_key_bootstrap: Mapping[str, Any],
    public_key: Mapping[str, Any],
    terminal_key_rotation: Mapping[str, Any],
) -> str:
    expected_keys = {
        "schema_name",
        "schema_version",
        "reviewer",
        "reviewer_agent_id",
        "attribution",
        "review_completed_at_utc",
        "verdict",
        "reviewed_source_set_sha256",
        "legacy_source_set_sha256",
        "prior_recovery_chain_sha256",
        "rejected_generations_sha256",
        "prior_failed_signed_recovery_sha256",
        "fresh_key_bootstrap_sha256",
        "trusted_recovery_public_key",
        "scope_sha256",
        "terminal_key_rotation_sha256",
        "readonly_probe",
        "high_findings",
        "medium_findings",
        "low_findings",
        "unresolved_conditions",
        "summary",
        "pass",
    }
    _exact_keys(payload, expected_keys, f"recovery review {role}")
    reviewer_id = payload.get("reviewer_agent_id")
    completed_at = payload.get("review_completed_at_utc")
    probe = payload.get("readonly_probe")
    attribution = payload.get("attribution")
    low_findings = payload.get("low_findings")
    summary = payload.get("summary")
    if (
        payload.get("schema_name") != REVIEW_SCHEMA_NAME
        or payload.get("schema_version") != SCHEMA_VERSION
        or payload.get("reviewer") != EXPECTED_REVIEWERS[role]
        or (
            role != "external_claude_opus"
            and reviewer_id != EXPECTED_REVIEWER_AGENT_IDS[role]
        )
        or (role == "external_claude_opus" and not _is_canonical_uuid(reviewer_id))
        or not isinstance(attribution, Mapping)
        or set(attribution)
        != {
            "type",
            "transport",
            "transport_id",
            "cryptographic_reviewer_authorship_proven",
            "semantics",
        }
        or attribution.get("type") != "orchestrator_recorded_transport"
        or attribution.get("transport")
        != ("claude_cli_session_envelope" if role == "external_claude_opus" else "multi_agent_v1")
        or attribution.get("transport_id") != reviewer_id
        or attribution.get("cryptographic_reviewer_authorship_proven") is not False
        or attribution.get("semantics") != REVIEW_ATTRIBUTION_SEMANTICS
        or not isinstance(completed_at, str)
        or not (completed_at.endswith("Z") or completed_at.endswith("+00:00"))
        or payload.get("verdict") != "APPROVE"
        or payload.get("reviewed_source_set_sha256")
        != _json_sha256(recovery_sources)
        or payload.get("legacy_source_set_sha256") != _json_sha256(legacy_sources)
        or payload.get("prior_recovery_chain_sha256")
        != _json_sha256(prior_recovery_chain)
        or payload.get("rejected_generations_sha256")
        != _json_sha256(rejected_generations)
        or payload.get("prior_failed_signed_recovery_sha256")
        != _json_sha256(prior_failed_signed_recovery)
        or payload.get("fresh_key_bootstrap_sha256")
        != _json_sha256(fresh_key_bootstrap)
        or not _strict_equal(payload.get("trusted_recovery_public_key"), public_key)
        or payload.get("scope_sha256") != _json_sha256(RECOVERY_SCOPE)
        or payload.get("terminal_key_rotation_sha256")
        != _json_sha256(terminal_key_rotation)
        or not isinstance(probe, Mapping)
        or set(probe)
        != {
            "command",
            "exit_code",
            "no_output_writes",
            "status",
            "forbidden_output_slots_checked",
            "regression_summary",
        }
        or probe.get("exit_code") != 0
        or probe.get("no_output_writes") is not True
        or probe.get("status") != "PASS"
        or probe.get("forbidden_output_slots_checked")
        != EXPECTED_CEREMONY_FORBIDDEN_OUTPUT_SLOT_COUNT
        or probe.get("regression_summary") != EXPECTED_REGRESSION_SUMMARY
        or type(probe.get("exit_code")) is not int
        or not _strict_equal(probe.get("command"), READONLY_PROBE_COMMAND)
        or not _strict_equal(payload.get("high_findings"), [])
        or not _strict_equal(payload.get("medium_findings"), [])
        or not isinstance(low_findings, list)
        or any(not isinstance(finding, str) or not finding for finding in low_findings)
        or not _strict_equal(payload.get("unresolved_conditions"), [])
        or not isinstance(summary, str)
        or not summary.strip()
        or payload.get("pass") is not True
    ):
        raise RecoveryAuthorityError(f"Recovery review did not exactly approve: {role}")
    return reviewer_id


def _build_validation_evidence(
    *,
    expected_runtime_role: str,
    recovery_sources: Mapping[str, Any],
    authority_record: Mapping[str, Any],
    signature_record: Mapping[str, Any],
    commit_record: Mapping[str, Any],
    public_record: Mapping[str, Any],
    expected_private: Mapping[str, Any],
    terminal_key_rotation: Mapping[str, Any],
    terminal_key_state: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "authority_id": AUTHORITY_ID,
        "authority": dict(authority_record),
        "signature": dict(signature_record),
        "commit": dict(commit_record),
        "authority_bundle": str(AUTHORITY_BUNDLE),
        "public_key": dict(public_record),
        "retired_private_key": dict(expected_private),
        "terminal_key_rotation": dict(terminal_key_rotation),
        "terminal_key_state": dict(terminal_key_state),
        "authority_validator": dict(recovery_sources["authority_validator"]),
        "runtime_role": expected_runtime_role,
        "runtime_source": dict(recovery_sources[expected_runtime_role]),
        "historical_link_count_observation": "NOT_RECORDED_NOT_INFERRED",
        "descriptor_lease_count": len(_FD_LEASES),
        "pass": True,
    }


def _validate_ceremony_execution_bindings(
    value: Mapping[str, Any], recovery_sources: Mapping[str, Any]
) -> None:
    _exact_keys(
        value,
        {
            "builder",
            "readonly_probe_and_regression_tests",
            "review_attribution_semantics",
            "pass",
        },
        "ceremony execution bindings",
    )
    builder = value.get("builder")
    probe = value.get("readonly_probe_and_regression_tests")
    if not isinstance(builder, Mapping) or not isinstance(probe, Mapping):
        raise RecoveryAuthorityError("Ceremony execution bindings are not mappings")
    _exact_keys(
        builder,
        {
            "source",
            "canonical_file",
            "sys_argv0",
            "proc_cmdline_contains_exact_canonical_script_once",
            "bound_before_validator_load",
            "same_path_inode_as_bound_source",
            "threat_model",
            "pass",
        },
        "ceremony builder execution binding",
    )
    _exact_keys(
        probe,
        {
            "probe_source",
            "python_runtime",
            "launch",
            "regression_test_source",
            "pass",
        },
        "probe/test execution binding",
    )
    regression = probe.get("regression_test_source")
    if isinstance(regression, Mapping):
        _exact_keys(
            regression,
            {
                "source",
                "python_runtime",
                "execution",
                "canonical_python_argv0",
                "pass",
            },
            "regression test execution binding",
        )
    if (
        not isinstance(regression, Mapping)
        or not _strict_equal(builder.get("source"), recovery_sources["ceremony_builder"])
        or builder.get("canonical_file") != str(CEREMONY_BUILDER)
        or builder.get("sys_argv0") != str(CEREMONY_BUILDER)
        or builder.get("proc_cmdline_contains_exact_canonical_script_once") is not True
        or builder.get("bound_before_validator_load") is not True
        or builder.get("same_path_inode_as_bound_source") is not True
        or builder.get("threat_model")
        != "trusted_same_uid_account_no_malicious_runtime_code_injection"
        or builder.get("pass") is not True
        or not _strict_equal(probe.get("probe_source"), recovery_sources["readonly_probe"])
        or not _strict_equal(probe.get("python_runtime"), EXPECTED_PYTHON_RUNTIME)
        or probe.get("launch")
        != "bound_python_fd_exec_a_canonical_alias_bound_probe_source_fd"
        or not _strict_equal(regression.get("source"), recovery_sources["regression_tests"])
        or not _strict_equal(regression.get("python_runtime"), EXPECTED_PYTHON_RUNTIME)
        or regression.get("execution")
        != (
            "pytest_preloaded_module_compiled_from_bound_source_fd_"
            "via_bound_python_fd"
        )
        or regression.get("canonical_python_argv0") != str(PYTHON)
        or regression.get("pass") is not True
        or probe.get("pass") is not True
        or value.get("review_attribution_semantics") != REVIEW_ATTRIBUTION_SEMANTICS
        or value.get("pass") is not True
    ):
        raise RecoveryAuthorityError("Ceremony execution binding contract drift")


def _validate_authority_payload(
    payload: Mapping[str, Any],
    *,
    transaction_id: str,
    public_record: Mapping[str, Any],
    expected_private: Mapping[str, Any],
    original_chain: Mapping[str, Any],
    prior_recovery_chain: Mapping[str, Any],
    prior_failed_signed_recovery: Mapping[str, Any],
    fresh_key_bootstrap: Mapping[str, Any],
    rejected_generations: Mapping[str, Any],
    legacy_sources: Mapping[str, Any],
    recovery_sources: Mapping[str, Any],
    review_records: Mapping[str, Any],
    terminal_key_rotation: Mapping[str, Any],
) -> None:
    expected_keys = {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "authority_id",
        "task_type",
        "status",
        "authority_bundle",
        "scope",
        "terminal_key_rotation",
        "public_key",
        "retired_private_key",
        "original_signed_chain",
        "prior_recovery_chain",
        "prior_failed_signed_recovery",
        "fresh_key_bootstrap",
        "rejected_unsigned_generations",
        "legacy_sources",
        "recovery_sources",
        "legacy_commands",
        "recovery_commands",
        "schemas",
        "receipt_paths",
        "review_artifacts",
        "review_attribution_semantics",
        "ceremony_execution_bindings",
        "checks",
        "pass",
        "pass_semantics",
    }
    _exact_keys(payload, expected_keys, "recovery authority")
    ceremony_execution_bindings = payload.get("ceremony_execution_bindings")
    if not isinstance(ceremony_execution_bindings, Mapping):
        raise RecoveryAuthorityError("Ceremony execution bindings are missing")
    _validate_ceremony_execution_bindings(
        ceremony_execution_bindings, recovery_sources
    )
    if (
        payload.get("schema_name") != SCHEMA_NAME
        or payload.get("schema_version") != SCHEMA_VERSION
        or payload.get("authority_id") != AUTHORITY_ID
        or payload.get("task_type") != "B_comprehensive_validation"
        or payload.get("status") != AUTHORITY_STATUS
        or not isinstance(payload.get("created_at_utc"), str)
        or not payload["created_at_utc"].endswith("+00:00")
        or not _strict_equal(
            payload.get("authority_bundle"),
            {
                "path": str(AUTHORITY_BUNDLE),
                "publication": "atomic_directory_rename_noreplace",
                "transaction_id": transaction_id,
                "members": [
                    "authority.ed25519.sig",
                    "authority.json",
                    "commit.json",
                ],
            },
        )
        or not _strict_equal(payload.get("scope"), RECOVERY_SCOPE)
        or not _strict_equal(
            payload.get("terminal_key_rotation"), terminal_key_rotation
        )
        or not _strict_equal(payload.get("public_key"), public_record)
        or not _strict_equal(payload.get("retired_private_key"), expected_private)
        or not _strict_equal(payload.get("original_signed_chain"), original_chain)
        or not _strict_equal(payload.get("prior_recovery_chain"), prior_recovery_chain)
        or not _strict_equal(
            payload.get("prior_failed_signed_recovery"), prior_failed_signed_recovery
        )
        or not _strict_equal(payload.get("fresh_key_bootstrap"), fresh_key_bootstrap)
        or not _strict_equal(
            payload.get("rejected_unsigned_generations"), rejected_generations
        )
        or not _strict_equal(payload.get("legacy_sources"), legacy_sources)
        or not _strict_equal(payload.get("recovery_sources"), recovery_sources)
        or not _strict_equal(payload.get("legacy_commands"), LEGACY_COMMANDS)
        or not _strict_equal(payload.get("recovery_commands"), RECOVERY_COMMANDS)
        or not _strict_equal(payload.get("schemas"), SCHEMAS)
        or not _strict_equal(payload.get("receipt_paths"), RECEIPT_PATHS)
        or not _strict_equal(payload.get("review_artifacts"), review_records)
        or payload.get("review_attribution_semantics") != REVIEW_ATTRIBUTION_SEMANTICS
        or not _strict_equal(payload.get("checks"), AUTHORITY_CHECKS)
        or payload.get("pass") is not True
        or payload.get("pass_semantics") != AUTHORITY_PASS_SEMANTICS
    ):
        raise RecoveryAuthorityError("Recovery authority exact contract drift")


def validate_recovery_authority(
    *,
    expected_runtime_role: str,
    expected_recovery_continuation_private_mode: str = "0o400",
) -> dict[str, Any]:
    if expected_runtime_role not in RECOVERY_SOURCE_PATHS:
        raise RecoveryAuthorityError(f"Unknown recovery runtime role: {expected_runtime_role}")

    _validate_ceremony_absence_contract()
    recovery_sources = _records(RECOVERY_SOURCE_PATHS)
    legacy_sources = _records(LEGACY_SOURCE_PATHS)
    original_chain = _records(ORIGINAL_CHAIN_PATHS)
    _validate_pinned_records(legacy_sources, original_chain)
    prior_recovery_records = _records(PRIOR_RECOVERY_CHAIN_PATHS)
    prior_recovery_chain = _validate_prior_recovery_chain(prior_recovery_records)
    rejected_generations = _validate_rejected_generations()
    prior_failed_signed_recovery = _validate_prior_failed_signed_recoveries()
    fresh_key_bootstrap = _validate_v30_key_bootstrap()
    terminal_key_rotation, terminal_key_state = _validate_terminal_key_rotation(
        expected_recovery_private_mode=expected_recovery_continuation_private_mode
    )

    public_data, public_record, public_fd = _open_bound(PUBLIC_KEY)
    if (
        len(public_data) != 113
        or public_record["sha256"] != EXPECTED_PUBLIC_KEY_SHA256
    ):
        raise RecoveryAuthorityError("Pinned recovery public key identity drift")
    expected_private, private_fd = _open_retired_private_key_bound(PRIVATE_KEY)
    del private_fd

    _open_bundle_directory()
    authority_data, authority_record, authority_fd = _open_bound(AUTHORITY)
    signature_data, signature_record, signature_fd = _open_bound(AUTHORITY_SIGNATURE)
    commit_data, commit_record, _ = _open_bound(AUTHORITY_COMMIT)
    if len(signature_data) != 64:
        raise RecoveryAuthorityError("Recovery authority signature size drift")
    _verify_signature(authority_fd, public_fd, signature_fd)
    commit_payload = _validate_bundle_commit(
        commit_data, authority_record, signature_record, expected_private
    )

    review_records: dict[str, dict[str, Any]] = {}
    review_payloads: dict[str, Mapping[str, Any]] = {}
    for role, path in REVIEW_PATHS.items():
        review_data, review_record, _ = _open_bound(path)
        review_records[role] = review_record
        review_payloads[role] = _load_json(review_data, f"recovery review {role}")
    reviewer_ids = {
        _validate_review(
            role,
            review_payloads[role],
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_record,
            terminal_key_rotation,
        )
        for role in REVIEW_PATHS
    }
    if len(reviewer_ids) != len(REVIEW_PATHS):
        raise RecoveryAuthorityError("Recovery reviewer identities are not distinct")

    payload = _load_json(authority_data, "recovery authority")
    _validate_authority_payload(
        payload,
        transaction_id=commit_payload["transaction_id"],
        public_record=public_record,
        expected_private=expected_private,
        original_chain=original_chain,
        prior_recovery_chain=prior_recovery_chain,
        prior_failed_signed_recovery=prior_failed_signed_recovery,
        fresh_key_bootstrap=fresh_key_bootstrap,
        rejected_generations=rejected_generations,
        legacy_sources=legacy_sources,
        recovery_sources=recovery_sources,
        review_records=review_records,
        terminal_key_rotation=terminal_key_rotation,
    )

    if len(HISTORICAL_EXECUTION_KEYS) != 8 or "link_count" in HISTORICAL_EXECUTION_KEYS:
        raise RecoveryAuthorityError("Historical identity schema is not exact eight-key")
    if len(CURRENT_IDENTITY_KEYS) != 9 or "link_count" not in CURRENT_IDENTITY_KEYS:
        raise RecoveryAuthorityError("Current identity schema is not exact nine-key")
    if set(RECEIPT_PATHS.values()) & {
        str(RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.v3.json"),
        str(RESULT_ROOT / "m2v5_runner_only_gate_replay.v1.json"),
        str(RESULT_ROOT / "m2v5_downstream_continuation.v1.json"),
    } != {
        RECEIPT_PATHS["legacy_verification"],
        RECEIPT_PATHS["legacy_replay"],
        RECEIPT_PATHS["legacy_continuation"],
    }:
        raise RecoveryAuthorityError("Legacy/recovery receipt path separation drift")

    required_paths = (
        set(RECOVERY_SOURCE_PATHS.values())
        | set(LEGACY_SOURCE_PATHS.values())
        | set(ORIGINAL_CHAIN_PATHS.values())
        | set(PRIOR_RECOVERY_CHAIN_PATHS.values())
        | set(REVIEW_PATHS.values())
        | {
            PUBLIC_KEY,
            PRIVATE_KEY,
            AUTHORIZED_CONTINUATION_PUBLIC_KEY,
            AUTHORIZED_CONTINUATION_PRIVATE_KEY,
            RECOVERY_CONTINUATION_PUBLIC_KEY,
            RECOVERY_CONTINUATION_PRIVATE_KEY,
            PRIOR_RECOVERY_PRIVATE_KEY,
            REJECTED_V2_EVIDENCE,
            REJECTED_V2_PRIVATE_KEY,
            REJECTED_V3_EVIDENCE,
            REJECTED_V3_PRIVATE_KEY,
            REJECTED_V4_EVIDENCE,
            REJECTED_V4_PRIVATE_KEY,
            REJECTED_V5_EVIDENCE,
            REJECTED_V5_PRIVATE_KEY,
            REJECTED_V6_EVIDENCE,
            REJECTED_V6_PRIVATE_KEY,
            *REJECTED_V6_ARCHIVED_REVIEWS.values(),
            REJECTED_V7_EVIDENCE,
            REJECTED_V7_PRIVATE_KEY,
            REJECTED_V8_EVIDENCE,
            REJECTED_V8_PRIVATE_KEY,
            FAILED_V9_EVIDENCE,
            FAILED_V9_PUBLIC_KEY,
            FAILED_V9_PRIVATE_KEY,
            *FAILED_V9_ARCHIVED_ARTIFACT_PATHS.values(),
            FAILED_V10_EVIDENCE,
            FAILED_V10_PUBLIC_KEY,
            FAILED_V10_PRIVATE_KEY,
            *FAILED_V10_SOURCE_PATHS.values(),
            *FAILED_V10_ARCHIVED_ARTIFACT_PATHS.values(),
            FAILED_V11_EVIDENCE,
            FAILED_V11_PUBLIC_KEY,
            FAILED_V11_PRIVATE_KEY,
            FAILED_V11_CONTINUATION_PUBLIC_KEY,
            FAILED_V11_CONTINUATION_PRIVATE_KEY,
            *FAILED_V11_SOURCE_PATHS.values(),
            *FAILED_V11_ARCHIVED_ARTIFACT_PATHS.values(),
            FAILED_V12_EVIDENCE,
            FAILED_V12_PUBLIC_KEY,
            FAILED_V12_PRIVATE_KEY,
            FAILED_V12_CONTINUATION_PUBLIC_KEY,
            FAILED_V12_CONTINUATION_PRIVATE_KEY,
            *FAILED_V12_SOURCE_PATHS.values(),
            *FAILED_V12_ARCHIVED_ARTIFACT_PATHS.values(),
            FAILED_V14_EVIDENCE,
            FAILED_V14_PUBLIC_KEY,
            FAILED_V14_PRIVATE_KEY,
            FAILED_V14_CONTINUATION_PUBLIC_KEY,
            FAILED_V14_CONTINUATION_PRIVATE_KEY,
            *FAILED_V14_SOURCE_PATHS.values(),
            *FAILED_V14_ARCHIVED_ARTIFACT_PATHS.values(),
            REJECTED_V15_EVIDENCE,
            REJECTED_V15_KEY_ARCHIVE_LEDGER,
            REJECTED_V19_ROUND1_EVIDENCE,
            *REJECTED_V19_ROUND1_REVIEW_PATHS.values(),
            REJECTED_V19_ROUND1_EXTERNAL_ENVELOPE,
            REJECTED_V19_ROUND1_EXTERNAL_STDERR,
            *REJECTED_V19_ROUND1_TRANSPORT_PATHS.values(),
            *REJECTED_V19_ROUND1_SOURCE_PATHS.values(),
            REJECTED_V29_ROUND1_EVIDENCE,
            REJECTED_V29_ROUND1_SUMMARY,
            *REJECTED_V29_ROUND1_PAYLOAD_PATHS.values(),
            *REJECTED_V15_REVIEW_PATHS.values(),
            REJECTED_V15_EXTERNAL_ENVELOPE,
            REJECTED_V15_EXTERNAL_STDERR,
            *REJECTED_V15_SOURCE_PATHS.values(),
            REJECTED_V15_AUTHORITY_KEY_ARCHIVE / "ed25519_public.pem",
            REJECTED_V15_AUTHORITY_KEY_ARCHIVE / "ed25519_private_one_time.pem",
            REJECTED_V15_AUTHORITY_KEY_ARCHIVE / "UNUSED_KEY_ARCHIVE_RECORD.v1.json",
            REJECTED_V15_TERMINAL_KEY_ARCHIVE / "ed25519_public.pem",
            REJECTED_V15_TERMINAL_KEY_ARCHIVE
            / "ed25519_private_one_time_resident.pem",
            REJECTED_V15_TERMINAL_KEY_ARCHIVE / "UNUSED_KEY_ARCHIVE_RECORD.v1.json",
            FAILED_V16_EVIDENCE,
            FAILED_V16_PUBLIC_KEY,
            FAILED_V16_PRIVATE_KEY,
            FAILED_V16_CONTINUATION_PUBLIC_KEY,
            FAILED_V16_CONTINUATION_PRIVATE_KEY,
            *FAILED_V16_SOURCE_PATHS.values(),
            *FAILED_V16_REVIEW_TRANSPORT_PATHS.values(),
            *FAILED_V16_ARCHIVED_ARTIFACT_PATHS.values(),
            FAILED_V17_EVIDENCE,
            FAILED_V17_PUBLIC_KEY,
            FAILED_V17_PRIVATE_KEY,
            FAILED_V17_CONTINUATION_PUBLIC_KEY,
            FAILED_V17_CONTINUATION_PRIVATE_KEY,
            *FAILED_V17_SOURCE_PATHS.values(),
            *FAILED_V17_REVIEW_TRANSPORT_PATHS.values(),
            *FAILED_V17_ARCHIVED_ARTIFACT_PATHS.values(),
            *FAILED_V17_TRACE_ROOT.iterdir(),
            FAILED_V18_EVIDENCE,
            FAILED_V18_PUBLIC_KEY,
            FAILED_V18_PRIVATE_KEY,
            FAILED_V18_CONTINUATION_PUBLIC_KEY,
            FAILED_V18_CONTINUATION_PRIVATE_KEY,
            *FAILED_V18_SOURCE_PATHS.values(),
            *FAILED_V18_REVIEW_TRANSPORT_PATHS.values(),
            *FAILED_V18_ARCHIVED_ARTIFACT_PATHS.values(),
            *FAILED_V18_TRACE_ROOT.glob(FAILED_V18_TRACE_GLOB),
            FAILED_V19_EVIDENCE,
            FAILED_V19_PUBLIC_KEY,
            FAILED_V19_PRIVATE_KEY,
            FAILED_V19_CONTINUATION_PUBLIC_KEY,
            FAILED_V19_CONTINUATION_PRIVATE_KEY,
            *FAILED_V19_SOURCE_PATHS.values(),
            *FAILED_V19_REVIEW_TRANSPORT_PATHS.values(),
            *FAILED_V19_ARCHIVED_ARTIFACT_PATHS.values(),
            REJECTED_V20_EVIDENCE,
            *REJECTED_V20_SOURCE_PATHS.values(),
            *REJECTED_V20_INTERNAL_REVIEW_PATHS.values(),
            REJECTED_V20_EXTERNAL_ENVELOPE,
            REJECTED_V20_EXTERNAL_STDERR,
            *REJECTED_V20_TRANSPORT_PATHS.values(),
            REJECTED_V20_PUBLIC_KEY,
            REJECTED_V20_PRIVATE_KEY,
            REJECTED_V20_CONTINUATION_PUBLIC_KEY,
            REJECTED_V20_CONTINUATION_PRIVATE_KEY,
            FAILED_V21_EVIDENCE,
            FAILED_V21_PUBLIC_KEY,
            FAILED_V21_PRIVATE_KEY,
            FAILED_V21_CONTINUATION_PUBLIC_KEY,
            FAILED_V21_CONTINUATION_PRIVATE_KEY,
            *FAILED_V21_SOURCE_PATHS.values(),
            *FAILED_V21_ARCHIVED_ARTIFACT_PATHS.values(),
            *FAILED_V21_REVIEW_TRANSPORT_PATHS.values(),
            REJECTED_V22_EVIDENCE,
            *REJECTED_V22_SOURCE_PATHS.values(),
            *REJECTED_V22_REVIEW_PATHS.values(),
            *REJECTED_V22_TRANSPORT_PATHS.values(),
            REJECTED_V22_PUBLIC_KEY,
            REJECTED_V22_PRIVATE_KEY,
            REJECTED_V22_CONTINUATION_PUBLIC_KEY,
            REJECTED_V22_CONTINUATION_PRIVATE_KEY,
            REJECTED_V23_EVIDENCE,
            *REJECTED_V23_SOURCE_PATHS.values(),
            *REJECTED_V23_REVIEW_PATHS.values(),
            *REJECTED_V23_TRANSPORT_PATHS.values(),
            REJECTED_V23_PUBLIC_KEY,
            REJECTED_V23_PRIVATE_KEY,
            REJECTED_V23_CONTINUATION_PUBLIC_KEY,
            REJECTED_V23_CONTINUATION_PRIVATE_KEY,
            REJECTED_V24_EVIDENCE,
            *REJECTED_V24_SOURCE_PATHS.values(),
            *REJECTED_V24_REVIEW_PATHS.values(),
            *REJECTED_V24_TRANSPORT_PATHS.values(),
            REJECTED_V24_PUBLIC_KEY,
            REJECTED_V24_PRIVATE_KEY,
            REJECTED_V24_CONTINUATION_PUBLIC_KEY,
            REJECTED_V24_CONTINUATION_PRIVATE_KEY,
            FAILED_V25_EVIDENCE,
            FAILED_V25_BUNDLE / "authority.json",
            FAILED_V25_BUNDLE / "authority.ed25519.sig",
            FAILED_V25_BUNDLE / "commit.json",
            FAILED_V25_PUBLIC_KEY,
            FAILED_V25_PRIVATE_KEY,
            FAILED_V25_TERMINAL_PUBLIC_KEY,
            FAILED_V25_TERMINAL_PRIVATE_KEY,
            *FAILED_V25_SOURCE_PATHS.values(),
            FAILED_V26_EVIDENCE,
            FAILED_V26_BUNDLE / "authority.json",
            FAILED_V26_BUNDLE / "authority.ed25519.sig",
            FAILED_V26_BUNDLE / "commit.json",
            FAILED_V26_PUBLIC_KEY,
            FAILED_V26_PRIVATE_KEY,
            FAILED_V26_AUTHORITY_KEY_ARCHIVE_RECORD,
            FAILED_V26_TERMINAL_PUBLIC_KEY,
            FAILED_V26_TERMINAL_PRIVATE_KEY,
            FAILED_V26_TERMINAL_KEY_ARCHIVE_RECORD,
            *FAILED_V26_SOURCE_PATHS.values(),
            *FAILED_V26_REVIEW_PATHS.values(),
            *FAILED_V26_EXECUTION_PATHS.values(),
            REJECTED_V27_EVIDENCE,
            *REJECTED_V27_SOURCE_PATHS.values(),
            *REJECTED_V27_REVIEW_PATHS.values(),
            *REJECTED_V27_TRANSPORT_PATHS.values(),
            REJECTED_V27_PUBLIC_KEY,
            REJECTED_V27_PRIVATE_KEY,
            REJECTED_V27_AUTHORITY_KEY_ARCHIVE_RECORD,
            REJECTED_V27_CONTINUATION_PUBLIC_KEY,
            REJECTED_V27_CONTINUATION_PRIVATE_KEY,
            REJECTED_V27_TERMINAL_KEY_ARCHIVE_RECORD,
            *tuple(path for path in FAILED_V28_ROOT.rglob("*") if path.is_file()),
            *tuple(
                path
                for spec in FAILED_V28_KEY_ARCHIVES.values()
                for path in (
                    spec["public"],
                    spec["private"],
                    spec["root"] / "FAILED_KEY_ARCHIVE_RECORD.v1.json",
                )
            ),
            *tuple(path for path in FAILED_V29_ROOT.rglob("*") if path.is_file()),
            *tuple(
                path
                for spec in FAILED_V29_KEY_ARCHIVES.values()
                for path in (
                    spec["root"] / "ed25519_public.pem",
                    spec["root"] / spec["private_name"],
                    spec["root"] / "FAILED_KEY_ARCHIVE_RECORD.v1.json",
                )
            ),
            V30_BOOTSTRAP_PREPARE,
            V30_BOOTSTRAP_PROGRESS,
            V30_BOOTSTRAP_RECEIPT,
            V30_BOOTSTRAP_SUCCESS,
            V30_RESULT_PUBLIC_KEY,
            V30_REPORT_PUBLIC_KEY,
            TUMOR_REF_SUMMARY_ALIAS,
            TUMOR_REF_SUMMARY_TARGET,
            AUTHORITY_BUNDLE,
            AUTHORITY,
            AUTHORITY_SIGNATURE,
            AUTHORITY_COMMIT,
        }
    )
    _require_leases(required_paths=required_paths)
    return _build_validation_evidence(
        expected_runtime_role=expected_runtime_role,
        recovery_sources=recovery_sources,
        authority_record=authority_record,
        signature_record=signature_record,
        commit_record=commit_record,
        public_record=public_record,
        expected_private=expected_private,
        terminal_key_rotation=terminal_key_rotation,
        terminal_key_state=terminal_key_state,
    )


def static_contract_probe() -> dict[str, Any]:
    ceremony_absence_contract = _validate_ceremony_absence_contract()
    recovery_sources = _records(RECOVERY_SOURCE_PATHS)
    legacy_sources = _records(LEGACY_SOURCE_PATHS)
    original_chain = _records(ORIGINAL_CHAIN_PATHS)
    _validate_pinned_records(legacy_sources, original_chain)
    prior_recovery_records = _records(PRIOR_RECOVERY_CHAIN_PATHS)
    prior_recovery_chain = _validate_prior_recovery_chain(prior_recovery_records)
    rejected_generations = _validate_rejected_generations()
    prior_failed_signed_recovery = _validate_prior_failed_signed_recoveries()
    fresh_key_bootstrap = _validate_v30_key_bootstrap()
    terminal_key_rotation, terminal_key_state = _validate_terminal_key_rotation(
        expected_recovery_private_mode="0o400"
    )
    public_data, public_record, _ = _open_bound(PUBLIC_KEY)
    if len(public_data) != 113 or public_record["sha256"] != EXPECTED_PUBLIC_KEY_SHA256:
        raise RecoveryAuthorityError("Pinned recovery public key identity drift")
    _require_leases(
        required_paths=(
            set(RECOVERY_SOURCE_PATHS.values())
            | set(LEGACY_SOURCE_PATHS.values())
            | set(ORIGINAL_CHAIN_PATHS.values())
            | set(PRIOR_RECOVERY_CHAIN_PATHS.values())
            | {
                PUBLIC_KEY,
                AUTHORIZED_CONTINUATION_PUBLIC_KEY,
                AUTHORIZED_CONTINUATION_PRIVATE_KEY,
                RECOVERY_CONTINUATION_PUBLIC_KEY,
                RECOVERY_CONTINUATION_PRIVATE_KEY,
                PRIOR_RECOVERY_PRIVATE_KEY,
                REJECTED_V2_EVIDENCE,
                REJECTED_V2_PRIVATE_KEY,
                REJECTED_V3_EVIDENCE,
                REJECTED_V3_PRIVATE_KEY,
                REJECTED_V4_EVIDENCE,
                REJECTED_V4_PRIVATE_KEY,
                REJECTED_V5_EVIDENCE,
                REJECTED_V5_PRIVATE_KEY,
                REJECTED_V6_EVIDENCE,
                REJECTED_V6_PRIVATE_KEY,
                *REJECTED_V6_ARCHIVED_REVIEWS.values(),
                REJECTED_V7_EVIDENCE,
                REJECTED_V7_PRIVATE_KEY,
                REJECTED_V8_EVIDENCE,
                REJECTED_V8_PRIVATE_KEY,
                FAILED_V9_EVIDENCE,
                FAILED_V9_PUBLIC_KEY,
                FAILED_V9_PRIVATE_KEY,
                *FAILED_V9_ARCHIVED_ARTIFACT_PATHS.values(),
                FAILED_V10_EVIDENCE,
                FAILED_V10_PUBLIC_KEY,
                FAILED_V10_PRIVATE_KEY,
                *FAILED_V10_SOURCE_PATHS.values(),
                *FAILED_V10_ARCHIVED_ARTIFACT_PATHS.values(),
                FAILED_V11_EVIDENCE,
                FAILED_V11_PUBLIC_KEY,
                FAILED_V11_PRIVATE_KEY,
                FAILED_V11_CONTINUATION_PUBLIC_KEY,
                FAILED_V11_CONTINUATION_PRIVATE_KEY,
                *FAILED_V11_SOURCE_PATHS.values(),
                *FAILED_V11_ARCHIVED_ARTIFACT_PATHS.values(),
                FAILED_V12_EVIDENCE,
                FAILED_V12_PUBLIC_KEY,
                FAILED_V12_PRIVATE_KEY,
                FAILED_V12_CONTINUATION_PUBLIC_KEY,
                FAILED_V12_CONTINUATION_PRIVATE_KEY,
                *FAILED_V12_SOURCE_PATHS.values(),
                *FAILED_V12_ARCHIVED_ARTIFACT_PATHS.values(),
                FAILED_V14_EVIDENCE,
                FAILED_V14_PUBLIC_KEY,
                FAILED_V14_PRIVATE_KEY,
                FAILED_V14_CONTINUATION_PUBLIC_KEY,
                FAILED_V14_CONTINUATION_PRIVATE_KEY,
                *FAILED_V14_SOURCE_PATHS.values(),
                *FAILED_V14_ARCHIVED_ARTIFACT_PATHS.values(),
                REJECTED_V15_EVIDENCE,
                REJECTED_V15_KEY_ARCHIVE_LEDGER,
                *REJECTED_V15_REVIEW_PATHS.values(),
                REJECTED_V15_EXTERNAL_ENVELOPE,
                REJECTED_V15_EXTERNAL_STDERR,
                *REJECTED_V15_SOURCE_PATHS.values(),
                REJECTED_V15_AUTHORITY_KEY_ARCHIVE / "ed25519_public.pem",
                REJECTED_V15_AUTHORITY_KEY_ARCHIVE
                / "ed25519_private_one_time.pem",
                REJECTED_V15_AUTHORITY_KEY_ARCHIVE
                / "UNUSED_KEY_ARCHIVE_RECORD.v1.json",
                REJECTED_V15_TERMINAL_KEY_ARCHIVE / "ed25519_public.pem",
                REJECTED_V15_TERMINAL_KEY_ARCHIVE
                / "ed25519_private_one_time_resident.pem",
                REJECTED_V15_TERMINAL_KEY_ARCHIVE
                / "UNUSED_KEY_ARCHIVE_RECORD.v1.json",
                REJECTED_V19_ROUND1_EVIDENCE,
                *REJECTED_V19_ROUND1_REVIEW_PATHS.values(),
                REJECTED_V19_ROUND1_EXTERNAL_ENVELOPE,
                REJECTED_V19_ROUND1_EXTERNAL_STDERR,
                *REJECTED_V19_ROUND1_TRANSPORT_PATHS.values(),
                *REJECTED_V19_ROUND1_SOURCE_PATHS.values(),
                REJECTED_V29_ROUND1_EVIDENCE,
                REJECTED_V29_ROUND1_SUMMARY,
                *REJECTED_V29_ROUND1_PAYLOAD_PATHS.values(),
                FAILED_V16_EVIDENCE,
                FAILED_V16_PUBLIC_KEY,
                FAILED_V16_PRIVATE_KEY,
                FAILED_V16_CONTINUATION_PUBLIC_KEY,
                FAILED_V16_CONTINUATION_PRIVATE_KEY,
                *FAILED_V16_SOURCE_PATHS.values(),
                *FAILED_V16_REVIEW_TRANSPORT_PATHS.values(),
                *FAILED_V16_ARCHIVED_ARTIFACT_PATHS.values(),
                FAILED_V17_EVIDENCE,
                FAILED_V17_PUBLIC_KEY,
                FAILED_V17_PRIVATE_KEY,
                FAILED_V17_CONTINUATION_PUBLIC_KEY,
                FAILED_V17_CONTINUATION_PRIVATE_KEY,
                *FAILED_V17_SOURCE_PATHS.values(),
                *FAILED_V17_REVIEW_TRANSPORT_PATHS.values(),
                *FAILED_V17_ARCHIVED_ARTIFACT_PATHS.values(),
                *FAILED_V17_TRACE_ROOT.iterdir(),
                FAILED_V18_EVIDENCE,
                FAILED_V18_PUBLIC_KEY,
                FAILED_V18_PRIVATE_KEY,
                FAILED_V18_CONTINUATION_PUBLIC_KEY,
                FAILED_V18_CONTINUATION_PRIVATE_KEY,
                *FAILED_V18_SOURCE_PATHS.values(),
                *FAILED_V18_REVIEW_TRANSPORT_PATHS.values(),
                *FAILED_V18_ARCHIVED_ARTIFACT_PATHS.values(),
                *FAILED_V18_TRACE_ROOT.glob(FAILED_V18_TRACE_GLOB),
                FAILED_V19_EVIDENCE,
                FAILED_V19_PUBLIC_KEY,
                FAILED_V19_PRIVATE_KEY,
                FAILED_V19_CONTINUATION_PUBLIC_KEY,
                FAILED_V19_CONTINUATION_PRIVATE_KEY,
                *FAILED_V19_SOURCE_PATHS.values(),
                *FAILED_V19_REVIEW_TRANSPORT_PATHS.values(),
                *FAILED_V19_ARCHIVED_ARTIFACT_PATHS.values(),
                REJECTED_V20_EVIDENCE,
                *REJECTED_V20_SOURCE_PATHS.values(),
                *REJECTED_V20_INTERNAL_REVIEW_PATHS.values(),
                REJECTED_V20_EXTERNAL_ENVELOPE,
                REJECTED_V20_EXTERNAL_STDERR,
                *REJECTED_V20_TRANSPORT_PATHS.values(),
                REJECTED_V20_PUBLIC_KEY,
                REJECTED_V20_PRIVATE_KEY,
                REJECTED_V20_CONTINUATION_PUBLIC_KEY,
                REJECTED_V20_CONTINUATION_PRIVATE_KEY,
                FAILED_V21_EVIDENCE,
                FAILED_V21_PUBLIC_KEY,
                FAILED_V21_PRIVATE_KEY,
                FAILED_V21_CONTINUATION_PUBLIC_KEY,
                FAILED_V21_CONTINUATION_PRIVATE_KEY,
                *FAILED_V21_SOURCE_PATHS.values(),
                *FAILED_V21_ARCHIVED_ARTIFACT_PATHS.values(),
                *FAILED_V21_REVIEW_TRANSPORT_PATHS.values(),
                REJECTED_V22_EVIDENCE,
                *REJECTED_V22_SOURCE_PATHS.values(),
                *REJECTED_V22_REVIEW_PATHS.values(),
                *REJECTED_V22_TRANSPORT_PATHS.values(),
                REJECTED_V22_PUBLIC_KEY,
                REJECTED_V22_PRIVATE_KEY,
                REJECTED_V22_CONTINUATION_PUBLIC_KEY,
                REJECTED_V22_CONTINUATION_PRIVATE_KEY,
                REJECTED_V23_EVIDENCE,
                *REJECTED_V23_SOURCE_PATHS.values(),
                *REJECTED_V23_REVIEW_PATHS.values(),
                *REJECTED_V23_TRANSPORT_PATHS.values(),
                REJECTED_V23_PUBLIC_KEY,
                REJECTED_V23_PRIVATE_KEY,
                REJECTED_V23_CONTINUATION_PUBLIC_KEY,
                REJECTED_V23_CONTINUATION_PRIVATE_KEY,
                REJECTED_V24_EVIDENCE,
                *REJECTED_V24_SOURCE_PATHS.values(),
                *REJECTED_V24_REVIEW_PATHS.values(),
                *REJECTED_V24_TRANSPORT_PATHS.values(),
                REJECTED_V24_PUBLIC_KEY,
                REJECTED_V24_PRIVATE_KEY,
                REJECTED_V24_CONTINUATION_PUBLIC_KEY,
                REJECTED_V24_CONTINUATION_PRIVATE_KEY,
                FAILED_V25_EVIDENCE,
                FAILED_V25_BUNDLE / "authority.json",
                FAILED_V25_BUNDLE / "authority.ed25519.sig",
                FAILED_V25_BUNDLE / "commit.json",
                FAILED_V25_PUBLIC_KEY,
                FAILED_V25_PRIVATE_KEY,
                FAILED_V25_TERMINAL_PUBLIC_KEY,
                FAILED_V25_TERMINAL_PRIVATE_KEY,
                *FAILED_V25_SOURCE_PATHS.values(),
                FAILED_V26_EVIDENCE,
                FAILED_V26_BUNDLE / "authority.json",
                FAILED_V26_BUNDLE / "authority.ed25519.sig",
                FAILED_V26_BUNDLE / "commit.json",
                FAILED_V26_PUBLIC_KEY,
                FAILED_V26_PRIVATE_KEY,
                FAILED_V26_AUTHORITY_KEY_ARCHIVE_RECORD,
                FAILED_V26_TERMINAL_PUBLIC_KEY,
                FAILED_V26_TERMINAL_PRIVATE_KEY,
                FAILED_V26_TERMINAL_KEY_ARCHIVE_RECORD,
                *FAILED_V26_SOURCE_PATHS.values(),
                *FAILED_V26_REVIEW_PATHS.values(),
                *FAILED_V26_EXECUTION_PATHS.values(),
                REJECTED_V27_EVIDENCE,
                *REJECTED_V27_SOURCE_PATHS.values(),
                *REJECTED_V27_REVIEW_PATHS.values(),
                *REJECTED_V27_TRANSPORT_PATHS.values(),
                REJECTED_V27_PUBLIC_KEY,
                REJECTED_V27_PRIVATE_KEY,
                REJECTED_V27_AUTHORITY_KEY_ARCHIVE_RECORD,
                REJECTED_V27_CONTINUATION_PUBLIC_KEY,
                REJECTED_V27_CONTINUATION_PRIVATE_KEY,
                REJECTED_V27_TERMINAL_KEY_ARCHIVE_RECORD,
                *tuple(path for path in FAILED_V28_ROOT.rglob("*") if path.is_file()),
                *tuple(
                    path
                    for spec in FAILED_V28_KEY_ARCHIVES.values()
                    for path in (
                        spec["public"],
                        spec["private"],
                        spec["root"] / "FAILED_KEY_ARCHIVE_RECORD.v1.json",
                    )
                ),
                *tuple(path for path in FAILED_V29_ROOT.rglob("*") if path.is_file()),
                *tuple(
                    path
                    for spec in FAILED_V29_KEY_ARCHIVES.values()
                    for path in (
                        spec["root"] / "ed25519_public.pem",
                        spec["root"] / spec["private_name"],
                        spec["root"] / "FAILED_KEY_ARCHIVE_RECORD.v1.json",
                    )
                ),
                V30_BOOTSTRAP_PREPARE,
                V30_BOOTSTRAP_PROGRESS,
                V30_BOOTSTRAP_RECEIPT,
                V30_BOOTSTRAP_SUCCESS,
                V30_RESULT_PUBLIC_KEY,
                V30_REPORT_PUBLIC_KEY,
                TUMOR_REF_SUMMARY_ALIAS,
                TUMOR_REF_SUMMARY_TARGET,
            }
        )
    )
    return {
        "authority_id": AUTHORITY_ID,
        "historical_execution_identity_keys": HISTORICAL_EXECUTION_KEYS,
        "current_identity_keys": CURRENT_IDENTITY_KEYS,
        "recovery_sources": recovery_sources,
        "legacy_sources": legacy_sources,
        "original_signed_chain": original_chain,
        "prior_recovery_chain": prior_recovery_chain,
        "prior_failed_signed_recovery": prior_failed_signed_recovery,
        "fresh_key_bootstrap": fresh_key_bootstrap,
        "rejected_unsigned_generations": rejected_generations,
        "ceremony_absence_contract": ceremony_absence_contract,
        "public_key": public_record,
        "terminal_key_rotation": terminal_key_rotation,
        "terminal_key_state": terminal_key_state,
        "pass": True,
    }


def main() -> int:
    try:
        invoked_source = Path(sys.argv[0]).resolve(strict=True)
        if invoked_source != VALIDATOR:
            raise RecoveryAuthorityError("Validator must execute canonical bound source")
        if sys.argv[1:] == ["--static-contract"]:
            result = static_contract_probe()
        elif len(sys.argv) == 3 and sys.argv[1] == "--validate-role":
            result = validate_recovery_authority(expected_runtime_role=sys.argv[2])
        else:
            raise RecoveryAuthorityError(
                "Usage: validator --static-contract | --validate-role ROLE"
            )
        print(json.dumps(result, ensure_ascii=True, indent=2, sort_keys=True))
        return 0
    except Exception as error:
        print(
            json.dumps(
                {
                    "schema_name": f"{SCHEMA_NAME}.error",
                    "schema_version": SCHEMA_VERSION,
                    "error": {"type": type(error).__name__, "message": str(error)},
                    "pass": False,
                },
                ensure_ascii=True,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
