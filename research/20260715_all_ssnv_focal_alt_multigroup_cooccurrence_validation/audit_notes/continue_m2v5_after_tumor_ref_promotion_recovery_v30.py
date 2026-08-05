#!/usr/bin/env python3
"""Continue the reviewed M2v5 Task-B runner after signed tumor-REF promotion.

The runner is intentionally fail-closed.  It treats the historical runner-only
replay as non-authoritative evidence, rebuilds the gate decision from bound file
descriptors, executes an exact static bridge around the already-materialized
tumor-REF receipt, and publishes one terminal receipt only after every downstream
release signature and relation has been independently revalidated.
"""

from __future__ import annotations

import ctypes
from datetime import datetime, timezone
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import stat
import subprocess
import sys
import types
from typing import Any, Callable, Mapping, NoReturn
import weakref


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
SCRIPT_ROOT = TOPIC_ROOT / "scripts"
RESULT_ROOT = TOPIC_ROOT / "results"
WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)

PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")
QA_PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/bin/python")
PRIMARY_PYTHON_RUNTIME = Path(
    "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python3.11"
)
QA_PYTHON_RUNTIME = Path("/bip7_disk/liaoyoyo2001/miniconda3/bin/python3.9")
BASH = Path("/usr/bin/bash")
OPENSSL = Path("/usr/bin/openssl")
NODE = Path("/bip7_disk/liaoyoyo2001/.nvm/versions/node/v22.22.1/bin/node")
QA_CHROMIUM = Path(
    "/bip7_disk/liaoyoyo2001/.cache/ms-playwright/"
    "chromium_headless_shell-1217/chrome-headless-shell-linux64/chrome-headless-shell"
)

CONTINUATION_RUNNER = (
    AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v30.py"
)
AUTHORIZED_CONTINUATION_RUNNER = (
    AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_v1.py"
)
PROMOTION_TOOL = AUDIT_ROOT / "promote_tumor_ref_source_receipt_v2.py"
CONTINUATION_VERIFIER = (
    AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v30.py"
)
AUTHORIZED_CONTINUATION_VERIFIER = (
    AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_v2.py"
)
RUNNER_GATE_REPLAY = (
    AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v30.py"
)
AUTHORIZED_RUNNER_GATE_REPLAY = (
    AUDIT_ROOT / "replay_m2v5_runner_only_gates_v1.py"
)
SCHEMA_RECOVERY_VALIDATOR = (
    AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v30.py"
)
EXPECTED_SCHEMA_RECOVERY_VALIDATOR_SHA256 = (
    "74c3e445fec567cf81184611b983d260b1b60344b8def63fd2feac099004daf1"
)
RECOVERY_OUTPUT_LEASE_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260724_v30_four_role_rebootstrap/ed25519_public.pem"
)
EXPECTED_RECOVERY_OUTPUT_LEASE_PUBLIC_KEY_SHA256 = (
    "a5b0e0b2c2a9f220d988597b47c8eb1d5446de401a932102948d829ffd0611ed"
)
COMPLETION_RUNNER = SCRIPT_ROOT / "run_m2v5_recovered_completion_chain.sh"
RELEASE_SOURCE_AUTHORITY = SCRIPT_ROOT / "release_source_authority.py"
FINAL_RELEASE_FINALIZER = (
    AUDIT_ROOT / "finalize_task_b_result_release_schema_recovery_v30.py"
)
FINAL_DATASET_BUILDER = (
    AUDIT_ROOT / "build_all_ssnv_final_report_dataset_schema_recovery_v29.py"
)
REPORT_BUILDER = (
    AUDIT_ROOT / "build_all_ssnv_report_artifact_schema_recovery_v29.py"
)
PORTABLE_PLUGIN_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/"
    "data-analytics/0.2.8-13ceeea1f599"
)
PORTABLE_PLUGIN_SCRIPT_ROOT = (
    PORTABLE_PLUGIN_ROOT / "skills" / "build-report" / "scripts"
)
PORTABLE_BUILDER_MODULE = PORTABLE_PLUGIN_SCRIPT_ROOT / "build_portable_artifact.mjs"
PORTABLE_CHART_EXTRACTOR_MODULE = (
    PORTABLE_PLUGIN_SCRIPT_ROOT / "extract_portable_chart_svgs.mjs"
)
PORTABLE_VERIFIER_MODULE = PORTABLE_PLUGIN_SCRIPT_ROOT / "verify_portable_artifact.mjs"
PORTABLE_BROWSER_HELPERS_MODULE = (
    PORTABLE_PLUGIN_SCRIPT_ROOT / "portable_browser_helpers.mjs"
)
PORTABLE_BROWSER_CLI_MODULE = PORTABLE_PLUGIN_SCRIPT_ROOT / "portable_browser_cli.mjs"
PORTABLE_SERVER_BUNDLE = PORTABLE_PLUGIN_ROOT / "mcp" / "server.cjs"
PORTABLE_ASSET_ROOT = PORTABLE_PLUGIN_ROOT / "assets"
PORTABLE_READER_ASSET_PARTS = tuple(
    PORTABLE_ASSET_ROOT / f"portable-artifact-reader.html.gz.b64.part{index:03d}"
    for index in range(1, 4)
)

AUTHORIZATION = RESULT_ROOT / "tumor_ref_source_receipt_promotion_authorization.v3.json"
PREPARE_EXIT_ATTESTATION = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_prepare_exit_attestation.v1.json"
)
AUTHORIZATION_SIGNATURE = Path(str(PREPARE_EXIT_ATTESTATION) + ".ed25519.sig")
LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE = Path(str(AUTHORIZATION) + ".ed25519.sig")
COMPLETION = RESULT_ROOT / "tumor_ref_source_receipt_promotion.v3.json"
PROMOTE_EXIT_ATTESTATION = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_promote_exit_attestation.v1.json"
)
COMPLETION_SIGNATURE = Path(str(PROMOTE_EXIT_ATTESTATION) + ".ed25519.sig")
LEGACY_COMPLETION_PAYLOAD_SIGNATURE = Path(str(COMPLETION) + ".ed25519.sig")
PROMOTION_VERIFICATION_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v30.json"
)
AUTHORIZED_PROMOTION_VERIFICATION_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.v3.json"
)
RUNNER_GATE_RECEIPT = RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.json"
RUNNER_GATE_LOG = RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.log"
RUNNER_GATE_SUCCESS_WITNESS = (
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.success_witness.json"
)
RUNNER_GATE_SUCCESS_WITNESS_SCHEMA = (
    "intersubmod.m2v5_runner_only_gate_replay.recovery.success_witness"
)
AUTHORIZED_RUNNER_GATE_RECEIPT = RESULT_ROOT / "m2v5_runner_only_gate_replay.v1.json"
PROMOTION_INCIDENT = RESULT_ROOT / "tumor_ref_source_receipt_promotion_incident.v3.json"
CONTINUATION_RECEIPT = RESULT_ROOT / "m2v5_downstream_continuation.recovery.v30.json"
CONTINUATION_EXIT_ATTESTATION = (
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v30.json"
)
CONTINUATION_SIGNATURE = Path(
    str(CONTINUATION_EXIT_ATTESTATION) + ".ed25519.sig"
)
CONTINUATION_SUCCESS_WITNESS = (
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v30.json"
)
CONTINUATION_INCIDENT = (
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v30.json"
)
AUTHORIZED_CONTINUATION_RECEIPT = (
    RESULT_ROOT / "m2v5_downstream_continuation.v1.json"
)
AUTHORIZED_CONTINUATION_EXIT_ATTESTATION = (
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.v1.json"
)
AUTHORIZED_CONTINUATION_SIGNATURE = Path(
    str(AUTHORIZED_CONTINUATION_EXIT_ATTESTATION) + ".ed25519.sig"
)

HISTORICAL_SOURCE_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json"
)
CANONICAL_SOURCE_RECEIPT = (
    WORKSPACE_ROOT
    / "tumor_ref_recovery_source_identity_v1"
    / "post_run_source_identity.receipt.json"
)

AUTHORIZATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_promotion_authorization/"
    "20260722_tumor_ref_receipt_promotion_v4"
)
AUTHORIZATION_PUBLIC_KEY = AUTHORIZATION_KEY_ROOT / "ed25519_public.pem"
AUTHORIZATION_PRIVATE_KEY = AUTHORIZATION_KEY_ROOT / (
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
COMPLETION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_promotion_completion/"
    "20260722_tumor_ref_receipt_promotion_v4"
)
COMPLETION_PUBLIC_KEY = COMPLETION_KEY_ROOT / "ed25519_public.pem"
COMPLETION_PRIVATE_KEY = COMPLETION_KEY_ROOT / (
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260722_m2v5_terminal_v2"
)
CONTINUATION_PUBLIC_KEY = CONTINUATION_KEY_ROOT / "ed25519_public.pem"
CONTINUATION_PRIVATE_KEY = CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
FAILED_V16_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v6"
)
FAILED_V16_CONTINUATION_PUBLIC_KEY = (
    FAILED_V16_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
FAILED_V16_CONTINUATION_PRIVATE_KEY = FAILED_V16_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
FAILED_V17_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v7"
)
FAILED_V17_CONTINUATION_PUBLIC_KEY = (
    FAILED_V17_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
FAILED_V17_CONTINUATION_PRIVATE_KEY = FAILED_V17_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
FAILED_V18_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v8"
)
FAILED_V18_CONTINUATION_PUBLIC_KEY = (
    FAILED_V18_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
FAILED_V18_CONTINUATION_PRIVATE_KEY = FAILED_V18_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
FAILED_V19_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v9"
)
FAILED_V19_CONTINUATION_PUBLIC_KEY = (
    FAILED_V19_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
FAILED_V19_CONTINUATION_PRIVATE_KEY = FAILED_V19_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
REJECTED_V20_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v10"
)
REJECTED_V20_CONTINUATION_PUBLIC_KEY = (
    REJECTED_V20_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
REJECTED_V20_CONTINUATION_PRIVATE_KEY = REJECTED_V20_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
FAILED_V21_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v11"
)
FAILED_V21_CONTINUATION_PUBLIC_KEY = (
    FAILED_V21_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
FAILED_V21_CONTINUATION_PRIVATE_KEY = FAILED_V21_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
REJECTED_V22_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v12"
)
REJECTED_V22_CONTINUATION_PUBLIC_KEY = (
    REJECTED_V22_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
REJECTED_V22_CONTINUATION_PRIVATE_KEY = REJECTED_V22_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
REJECTED_V23_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v13"
)
REJECTED_V23_CONTINUATION_PUBLIC_KEY = (
    REJECTED_V23_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
REJECTED_V23_CONTINUATION_PRIVATE_KEY = REJECTED_V23_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
REJECTED_V24_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260723_m2v5_terminal_v14"
)
REJECTED_V24_CONTINUATION_PUBLIC_KEY = (
    REJECTED_V24_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
REJECTED_V24_CONTINUATION_PRIVATE_KEY = REJECTED_V24_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
FAILED_V25_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
    "20260723_m2v5_terminal_v15_unused_after_signed_v25_c_failure_01"
)
FAILED_V25_CONTINUATION_PUBLIC_KEY = (
    FAILED_V25_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
FAILED_V25_CONTINUATION_PRIVATE_KEY = FAILED_V25_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
FAILED_V26_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
    "20260724_m2v5_terminal_v16_unused_after_signed_v26_c_failure_01"
)
FAILED_V26_CONTINUATION_PUBLIC_KEY = (
    FAILED_V26_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
FAILED_V26_CONTINUATION_PRIVATE_KEY = FAILED_V26_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
REJECTED_V27_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
    "20260724_m2v5_terminal_v17_unused_pre_authority_review_rejection_01"
)
REJECTED_V27_CONTINUATION_PUBLIC_KEY = (
    REJECTED_V27_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
REJECTED_V27_CONTINUATION_PRIVATE_KEY = REJECTED_V27_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
FAILED_V28_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
    "20260724_m2v5_terminal_v18_unused_after_signed_v28_c_failure_01"
)
FAILED_V28_CONTINUATION_PUBLIC_KEY = (
    FAILED_V28_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
FAILED_V28_CONTINUATION_PRIVATE_KEY = FAILED_V28_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
FAILED_V29_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/archive/"
    "20260724_m2v5_terminal_v19_unused_after_signed_v29_r_failure_01"
)
FAILED_V29_CONTINUATION_PUBLIC_KEY = (
    FAILED_V29_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
FAILED_V29_CONTINUATION_PRIVATE_KEY = FAILED_V29_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)
FAILED_V28_RESULT_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/archive/"
    "20260724_all_ssnv_result_v5_consumed_v28_c_failure_01"
)
FAILED_V28_REPORT_KEY_ARCHIVE = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/archive/"
    "20260724_all_ssnv_report_v5_unused_v28_c_failure_01"
)
FAILED_V28_RESULT_PUBLIC_KEY = FAILED_V28_RESULT_KEY_ARCHIVE / "ed25519_public.pem"
FAILED_V28_RESULT_KEY_ARCHIVE_RECORD = (
    FAILED_V28_RESULT_KEY_ARCHIVE / "FAILED_KEY_ARCHIVE_RECORD.v1.json"
)
FAILED_V28_REPORT_PUBLIC_KEY = FAILED_V28_REPORT_KEY_ARCHIVE / "ed25519_public.pem"
FAILED_V28_REPORT_PRIVATE_KEY = (
    FAILED_V28_REPORT_KEY_ARCHIVE
    / "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
FAILED_V28_REPORT_KEY_ARCHIVE_RECORD = (
    FAILED_V28_REPORT_KEY_ARCHIVE / "FAILED_KEY_ARCHIVE_RECORD.v1.json"
)
RECOVERY_CONTINUATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260724_m2v5_terminal_v21"
)
RECOVERY_CONTINUATION_PUBLIC_KEY = (
    RECOVERY_CONTINUATION_KEY_ROOT / "ed25519_public.pem"
)
RECOVERY_CONTINUATION_PRIVATE_KEY = RECOVERY_CONTINUATION_KEY_ROOT / (
    "ed25519_private_one_time_resident.pem"
)

SUPERVISOR_CAPABILITY_FD = 198
LINUX_MFD_CLOEXEC = 0x0001
LINUX_MFD_ALLOW_SEALING = 0x0002
LINUX_F_ADD_SEALS = 1033
LINUX_F_GET_SEALS = 1034
LINUX_F_SEAL_SEAL = 0x0001
LINUX_F_SEAL_SHRINK = 0x0002
LINUX_F_SEAL_GROW = 0x0004
LINUX_F_SEAL_WRITE = 0x0008
SUPERVISOR_CAPABILITY_SCHEMA = (
    "intersubmod.m2v5_downstream_continuation.supervisor_capability"
)
SUPERVISOR_CAPABILITY_SEAL_NAMES = (
    "F_SEAL_SEAL",
    "F_SEAL_SHRINK",
    "F_SEAL_GROW",
    "F_SEAL_WRITE",
)
SUPERVISOR_CAPABILITY_PROTOCOL = {
    "mechanism": "sealed_memfd_with_256_bit_nonce_and_live_parent_binding",
    "inherited_descriptor": SUPERVISOR_CAPABILITY_FD,
    "required_memfd_mode": "0o400",
    "required_memfd_link_count": 0,
    "required_seals": list(SUPERVISOR_CAPABILITY_SEAL_NAMES),
    "direct_supervised_child_invocation_authorized": False,
}
CONTINUATION_THREAT_MODEL = {
    "trusted_boundary": "linux_kernel_and_uncompromised_same_uid_research_account",
    "protects_against": [
        "accidental_direct_child_invocation",
        "partial_or_nonzero_child_execution",
        "pathname_replacement_and_unsigned_artifact_tampering",
    ],
    "malicious_same_uid_process_in_scope": False,
    "same_uid_limitation": (
        "A malicious process running as the key-owning UID can read or chmod a mode-0400 "
        "file key and therefore is outside this process-provenance claim."
    ),
    "hostile_same_uid_requires": (
        "an isolated service account plus privileged signer, TPM/HSM, or equivalent "
        "hardware-backed policy enforcement"
    ),
}

V7_AUTHORITY = (
    REPO_ROOT
    / "docs"
    / "provenance"
    / "source_authorities"
    / "20260722_all_ssnv_focal_alt_release_source_authority.v7.json"
)
V7_APPROVAL = V7_AUTHORITY.with_name(
    V7_AUTHORITY.name.removesuffix(".json") + ".approval.json"
)
V7_SIGNATURE = Path(str(V7_APPROVAL) + ".ed25519.sig")
V7_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
    "20260722_all_ssnv_v10_strict_command_parity_bootstrap/ed25519_public.pem"
)
V7_PRIVATE_KEY = V7_PUBLIC_KEY.with_name(
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)

SCREEN = WORKSPACE_ROOT / "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
SCREEN_SITES = SCREEN / "all_ssnv_site_results.tsv.gz"
SCREEN_ASSIGNMENTS = SCREEN / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
SCREEN_SUMMARY = SCREEN / "all_ssnv_summary.json"
MANIFEST = RESULT_ROOT / "all_ssnv_input_manifest.json"
PRIMARY_PRE = RESULT_ROOT / "stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream.json"
LEGACY_PRIMARY_PRE = RESULT_ROOT / "stable_primary_artifact_audit.v1_pre_downstream.json"
COOCCURRENCE = (
    WORKSPACE_ROOT
    / "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity"
)
TUMOR_REF = WORKSPACE_ROOT / "all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel"
TUMOR_REF_SUMMARY = TUMOR_REF / "all_ssnv_tumor_ref_control_summary.json"
TUMOR_REF_SUMMARY_ALIAS = TUMOR_REF / "summary.json"
SOURCE_SNAPSHOT = (
    WORKSPACE_ROOT
    / "tumor_ref_recovery_source_identity_v1"
    / "observed_during_execution.snapshot.json"
)
INDEPENDENT_M2_AUDIT = RESULT_ROOT / "independent_m2_gate_recount.v3.json"
COOCCURRENCE_PREFLIGHT = (
    RESULT_ROOT / "cooccurrence_task_contract_preflight.v9_command_parity_full_runtime.json"
)
CLAIM_CONTRACT = TOPIC_ROOT / "claim-contract-v5.md"
FINAL_RELEASE_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
    "20260724_all_ssnv_result_v9_v30_recovery/ed25519_public.pem"
)
FINAL_RELEASE_PRIVATE_KEY = FINAL_RELEASE_PUBLIC_KEY.with_name(
    "ed25519_private_one_time.pem"
)
REPORT_RELEASE_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
    "20260724_all_ssnv_report_v9_v30_recovery/ed25519_public.pem"
)
REPORT_RELEASE_PRIVATE_KEY = REPORT_RELEASE_PUBLIC_KEY.with_name(
    "ed25519_private_one_time.pem"
)
PRIVATE_KEY_PATHS = frozenset(
    {
        AUTHORIZATION_PRIVATE_KEY,
        COMPLETION_PRIVATE_KEY,
        CONTINUATION_PRIVATE_KEY,
        FAILED_V16_CONTINUATION_PRIVATE_KEY,
        FAILED_V17_CONTINUATION_PRIVATE_KEY,
        FAILED_V18_CONTINUATION_PRIVATE_KEY,
        FAILED_V19_CONTINUATION_PRIVATE_KEY,
        REJECTED_V20_CONTINUATION_PRIVATE_KEY,
        FAILED_V21_CONTINUATION_PRIVATE_KEY,
        REJECTED_V22_CONTINUATION_PRIVATE_KEY,
        REJECTED_V23_CONTINUATION_PRIVATE_KEY,
        REJECTED_V24_CONTINUATION_PRIVATE_KEY,
        FAILED_V25_CONTINUATION_PRIVATE_KEY,
        FAILED_V26_CONTINUATION_PRIVATE_KEY,
        REJECTED_V27_CONTINUATION_PRIVATE_KEY,
        FAILED_V28_CONTINUATION_PRIVATE_KEY,
        FAILED_V29_CONTINUATION_PRIVATE_KEY,
        RECOVERY_CONTINUATION_PRIVATE_KEY,
        FINAL_RELEASE_PRIVATE_KEY,
        REPORT_RELEASE_PRIVATE_KEY,
        FAILED_V28_REPORT_PRIVATE_KEY,
    }
)
TERMINAL_KEY_ROTATION_KEYS = frozenset(
    {
        "legacy_signed_contract",
        "failed_v16_contract",
        "failed_v17_contract",
        "failed_v18_contract",
        "failed_v19_contract",
        "rejected_v20_contract",
        "failed_v21_contract",
        "rejected_v22_contract",
        "rejected_v23_contract",
        "rejected_v24_contract",
        "failed_v25_contract",
        "failed_v26_contract",
        "rejected_v27_contract",
        "failed_v28_contract",
        "failed_v29_contract",
        "recovery_v30_contract",
        "keys_are_pairwise_distinct",
        "legacy_signed_payload_changed",
    }
)
TERMINAL_KEY_STATE_KEYS = frozenset(
    {
        "legacy_public_key",
        "legacy_private_key",
        "failed_v16_public_key",
        "failed_v16_private_key",
        "failed_v17_public_key",
        "failed_v17_private_key",
        "failed_v18_public_key",
        "failed_v18_private_key",
        "failed_v19_public_key",
        "failed_v19_private_key",
        "rejected_v20_public_key",
        "rejected_v20_private_key",
        "failed_v21_public_key",
        "failed_v21_private_key",
        "rejected_v22_public_key",
        "rejected_v22_private_key",
        "rejected_v23_public_key",
        "rejected_v23_private_key",
        "rejected_v24_public_key",
        "rejected_v24_private_key",
        "failed_v25_public_key",
        "failed_v25_private_key",
        "failed_v26_public_key",
        "failed_v26_private_key",
        "rejected_v27_public_key",
        "rejected_v27_private_key",
        "failed_v28_public_key",
        "failed_v28_private_key",
        "failed_v29_public_key",
        "failed_v29_private_key",
        "recovery_public_key",
        "recovery_private_key",
    }
)
EARLIER_FP_REPORT = (
    REPO_ROOT
    / "research"
    / "20260715_single_fp_alt_multicluster_subclone_validation"
    / "20260715_單一FP_sSNV_ALT_read甲基多群與subclone假說全量驗證_01.md"
)
INTERSUBMOD_BIN = REPO_ROOT / "build" / "bin" / "inter_sub_mod"
REFERENCE = Path("/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta")

STRICT = WORKSPACE_ROOT / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5"
MATCHED_RUN = WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5"
MATCHED_ANALYSIS = (
    WORKSPACE_ROOT / "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5"
)
CN_CCF = WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5"
PRIMARY_POST = RESULT_ROOT / "stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json"
FROZEN_POST = (
    RESULT_ROOT / "frozen_input_immutability.post_m2v5_downstream_v3_source_authority_v5.json"
)
FINAL_DATASET = WORKSPACE_ROOT / "all_ssnv_final_report_dataset_v5_m2v5_source_attested"
FINAL_REPORT = WORKSPACE_ROOT / "all_ssnv_final_report_v5_m2v5_source_attested"
ORIGINAL_RUNNER_LOG = WORKSPACE_ROOT / "m2v5_recovered_completion_chain_v2_source_authority_v5.log"
ORIGINAL_RUNNER_CACHE = WORKSPACE_ROOT / ".python_cache_m2v5_completion_v2_bound_bootstrap"
FINAL_DATASET_RELEASE_RECEIPT = RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json"
FINAL_DATASET_RELEASE_SIGNATURE = Path(str(FINAL_DATASET_RELEASE_RECEIPT) + ".ed25519.sig")
FINAL_REPORT_RELEASE_RECEIPT = RESULT_ROOT / "task_b_final_report_release_receipt.v1.json"
FINAL_REPORT_RELEASE_SIGNATURE = Path(str(FINAL_REPORT_RELEASE_RECEIPT) + ".ed25519.sig")

PORTABLE_HTML = FINAL_REPORT / "all_ssnv_focal_alt_multigroup_cooccurrence_report.standalone.html"
REPORT_OUTPUTS = {
    "report_markdown": FINAL_REPORT / "report.md",
    "canonical_report_artifact": FINAL_REPORT / "artifact.json",
    "report_build_receipt": FINAL_REPORT / "report_build_receipt.json",
    "portable_html": PORTABLE_HTML,
    "portable_delivery_receipt": FINAL_REPORT / "portable_report_delivery_receipt.json",
    "portable_official_verification_screenshot": FINAL_REPORT / "portable_report_official_verification.png",
    "portable_desktop_screenshot": FINAL_REPORT / "portable_report_desktop_1440x1000.png",
    "portable_mobile_screenshot": FINAL_REPORT / "portable_report_mobile_390x844.png",
    "portable_desktop_qa": FINAL_REPORT / "portable_report_desktop_qa.json",
    "portable_mobile_qa": FINAL_REPORT / "portable_report_mobile_qa.json",
}
FINAL_DATASET_OUTPUTS = {
    "final_report_dataset": FINAL_DATASET / "final_report_dataset.json",
    "per_sample_metrics": FINAL_DATASET / "per_sample_metrics.tsv",
    "candidate_catalog": FINAL_DATASET / "candidate_catalog.tsv",
    "candidate_witness_pairs": FINAL_DATASET / "candidate_witness_pairs.tsv",
    "claim_ladder": FINAL_DATASET / "claim_ladder.tsv",
}

BASE_DOWNSTREAM_OUTPUT_SLOTS = (
    STRICT,
    MATCHED_RUN,
    MATCHED_ANALYSIS,
    CN_CCF,
    PRIMARY_POST,
    FROZEN_POST,
    FINAL_DATASET,
    FINAL_REPORT,
    FINAL_DATASET_RELEASE_RECEIPT,
    FINAL_DATASET_RELEASE_SIGNATURE,
    FINAL_REPORT_RELEASE_RECEIPT,
    FINAL_REPORT_RELEASE_SIGNATURE,
    ORIGINAL_RUNNER_LOG,
    ORIGINAL_RUNNER_CACHE,
)
AUTHORIZED_DOWNSTREAM_OUTPUT_SLOTS = (
    WORKSPACE_ROOT / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT
    / "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5",
    RESULT_ROOT / "stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json",
    RESULT_ROOT / "frozen_input_immutability.post_m2v5_downstream_v3_source_authority_v5.json",
    FINAL_DATASET,
    FINAL_REPORT,
    FINAL_DATASET_RELEASE_RECEIPT,
    FINAL_DATASET_RELEASE_SIGNATURE,
    FINAL_REPORT_RELEASE_RECEIPT,
    FINAL_REPORT_RELEASE_SIGNATURE,
    WORKSPACE_ROOT / "m2v5_recovered_completion_chain_v2_source_authority_v5.log",
    WORKSPACE_ROOT / ".python_cache_m2v5_completion_v2_bound_bootstrap",
    AUTHORIZED_CONTINUATION_RECEIPT,
    AUTHORIZED_CONTINUATION_EXIT_ATTESTATION,
    AUTHORIZED_CONTINUATION_SIGNATURE,
)
CURRENT_DOWNSTREAM_OUTPUT_SLOTS = BASE_DOWNSTREAM_OUTPUT_SLOTS + (
    AUTHORIZED_CONTINUATION_RECEIPT,
    AUTHORIZED_CONTINUATION_EXIT_ATTESTATION,
    AUTHORIZED_CONTINUATION_SIGNATURE,
)
DOWNSTREAM_OUTPUT_SLOTS = CURRENT_DOWNSTREAM_OUTPUT_SLOTS + (
    CONTINUATION_RECEIPT,
    CONTINUATION_EXIT_ATTESTATION,
    CONTINUATION_SIGNATURE,
    CONTINUATION_SUCCESS_WITNESS,
    CONTINUATION_INCIDENT,
)
REQUIRED_CHILD_OUTPUT_ROOTS = {
    "strict": STRICT,
    "matched_run": MATCHED_RUN,
    "cn_ccf": CN_CCF,
    "final_dataset": FINAL_DATASET,
    "final_report": FINAL_REPORT,
}
OPTIONAL_CHILD_OUTPUT_ROOTS = {"matched_analysis": MATCHED_ANALYSIS}
REQUIRED_CHILD_OUTPUT_FILES = {
    "primary_post_audit": PRIMARY_POST,
    "frozen_post_audit": FROZEN_POST,
}

EXPECTED_ENVIRONMENT = {
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
PROMOTION_CACHE_ROOTS = {
    "--prepare-authorization": WORKSPACE_ROOT / ".python_cache_tumor_ref_promotion_v2_prepare",
    "--promote": WORKSPACE_ROOT / ".python_cache_tumor_ref_promotion_v2_promote",
}

EXPECTED_PYTHON_SHA256 = "777797a57eb75b28f530191628d26b14afada9ced2cb51c0ecae1eb62796062e"
EXPECTED_QA_PYTHON_SHA256 = (
    "cbd71e11d993aa3b7cc3ad553b97733f8054c1465ed7764dbbd386b9897bd5bc"
)
EXPECTED_BASH_SHA256 = "59474588a312b6b6e73e5a42a59bf71e62b55416b6c9d5e4a6e1c630c2a9ecd4"
EXPECTED_OPENSSL_SHA256 = "38c064f53a6619364c7947fc11cad09e1fc4218fc2c3e9016a7786b1341ecd08"
EXPECTED_NODE_SHA256 = "243fd8938011479f41b3de101842150fa990f33fbbb3f7aabd330857f2d79e1d"
EXPECTED_QA_CHROMIUM_SHA256 = "22d72b589f908111fa29ff4e61ead8a0601e925177f59c70e2e851c22e9b1886"
EXPECTED_PORTABLE_BUILDER_SHA256 = (
    "0b86883810cf23c81f563a48a7708283a6a4cadf11228ad261efb64426ef728e"
)
EXPECTED_PORTABLE_CHART_EXTRACTOR_SHA256 = (
    "77e78593d8971563b6394334128dd3ee243f2182d3d087591d3a5499ece9283b"
)
EXPECTED_PORTABLE_VERIFIER_SHA256 = (
    "b495c4cc34113fb2918118eac302f4e7e2152c1b1b9b63e3646cea97ecbf9b3f"
)
EXPECTED_PORTABLE_BROWSER_HELPERS_SHA256 = (
    "84aa4f8a2a11376ebee6942d2d7e10a083d16c65dfdb73114d7b34e51c27f69d"
)
EXPECTED_PORTABLE_BROWSER_CLI_SHA256 = (
    "aac3b12fc12c7ad2e044533f881791dd3b23bd0ee31ddf6682dd7f6de99e6596"
)
EXPECTED_PORTABLE_SERVER_BUNDLE_SHA256 = (
    "eff59c6085d2ab6b6153c80a03749e764e160f8c6711da8433f7bd6762e1db66"
)
EXPECTED_PORTABLE_READER_ASSET_SHA256 = (
    "154f1d561bab28174f88d71ae710599709e5a5eda64dee08b025a699449dbbfc",
    "9459ed4b76bd825daba9637564dd3122ca617a88dd16531d3e9bca5c7ced080e",
    "cdfe2e6787faa61d37f043df6996d5f988c07360bb2b3d8af2aa3c18e8db3ac0",
)
REVIEWED_RUNTIME_SOURCE_CONTRACTS = {
    "recovery_final_dataset_builder": {
        "path": FINAL_DATASET_BUILDER,
        "sha256": "cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f",
        "mode": "0o444",
    },
    "recovery_result_report_finalizer": {
        "path": FINAL_RELEASE_FINALIZER,
        "sha256": "b9645cf4f57653f078357421592cef19a168db8389ee68ae90898b9c8c63d318",
        "mode": "0o444",
    },
    "recovery_report_builder": {
        "path": REPORT_BUILDER,
        "sha256": "fe19be151bfd72978f87f4a003e8fb8732f54c21e49f702ea5c575b92f4ee9ae",
        "mode": "0o444",
    },
    "primary_python_runtime": {
        "path": PRIMARY_PYTHON_RUNTIME,
        "sha256": EXPECTED_PYTHON_SHA256,
        "mode": "0o775",
    },
    "qa_python_runtime": {
        "path": QA_PYTHON_RUNTIME,
        "sha256": EXPECTED_QA_PYTHON_SHA256,
        "mode": "0o755",
    },
    "portable_builder_module": {
        "path": PORTABLE_BUILDER_MODULE,
        "sha256": EXPECTED_PORTABLE_BUILDER_SHA256,
        "mode": "0o644",
    },
    "portable_chart_extractor_module": {
        "path": PORTABLE_CHART_EXTRACTOR_MODULE,
        "sha256": EXPECTED_PORTABLE_CHART_EXTRACTOR_SHA256,
        "mode": "0o644",
    },
    "portable_verifier_module": {
        "path": PORTABLE_VERIFIER_MODULE,
        "sha256": EXPECTED_PORTABLE_VERIFIER_SHA256,
        "mode": "0o644",
    },
    "portable_browser_helpers_module": {
        "path": PORTABLE_BROWSER_HELPERS_MODULE,
        "sha256": EXPECTED_PORTABLE_BROWSER_HELPERS_SHA256,
        "mode": "0o644",
    },
    "portable_browser_cli_module": {
        "path": PORTABLE_BROWSER_CLI_MODULE,
        "sha256": EXPECTED_PORTABLE_BROWSER_CLI_SHA256,
        "mode": "0o644",
    },
    "portable_server_bundle": {
        "path": PORTABLE_SERVER_BUNDLE,
        "sha256": EXPECTED_PORTABLE_SERVER_BUNDLE_SHA256,
        "mode": "0o644",
    },
    **{
        f"portable_reader_asset_part{index:03d}": {
            "path": path,
            "sha256": EXPECTED_PORTABLE_READER_ASSET_SHA256[index - 1],
            "mode": "0o644",
        }
        for index, path in enumerate(PORTABLE_READER_ASSET_PARTS, start=1)
    },
}
HISTORICAL_SIGNED_RUNTIME_SOURCE_ROLES = frozenset(
    {
        "primary_python_runtime",
        "qa_python_runtime",
        "portable_builder_module",
        "portable_chart_extractor_module",
        "portable_verifier_module",
        "portable_browser_helpers_module",
        "portable_browser_cli_module",
        "portable_server_bundle",
        "portable_reader_asset_part001",
        "portable_reader_asset_part002",
        "portable_reader_asset_part003",
    }
)
RECOVERY_ONLY_RUNTIME_SOURCE_ROLES = frozenset(
    {
        "recovery_final_dataset_builder",
        "recovery_result_report_finalizer",
        "recovery_report_builder",
    }
)
MANDATORY_GATE_EXECUTABLE_CONTRACTS = {
    "primary_python_runtime": REVIEWED_RUNTIME_SOURCE_CONTRACTS[
        "primary_python_runtime"
    ],
    "qa_python_runtime": REVIEWED_RUNTIME_SOURCE_CONTRACTS["qa_python_runtime"],
    "node": {
        "path": NODE,
        "sha256": EXPECTED_NODE_SHA256,
        "mode": "0o755",
    },
    "qa_chromium": {
        "path": QA_CHROMIUM,
        "sha256": EXPECTED_QA_CHROMIUM_SHA256,
        "mode": "0o755",
    },
}
EXPECTED_AUTHORIZATION_PUBLIC_KEY_SHA256 = (
    "e638b29a1d9c207dfa9849ff20a3867dacf702f7c1a899f1968ea1401a839677"
)
EXPECTED_COMPLETION_PUBLIC_KEY_SHA256 = (
    "6c55230f368db6b1fb0eedbc6c2362c0df7ee41a2699764fadb1f424e7036032"
)
EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "a98370415594ebc9e8eb166bb70cab3550067b8a83db7d73ad23fac0891cb8b5"
)
EXPECTED_FAILED_V16_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "066949c0c36be413cd2fb60670e5a2fbc583ab9a3a4264ebf4d3766aba39867f"
)
EXPECTED_FAILED_V17_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "4b83e655f1a7a778691e27dc3df2257a230001c702afdf6e703e578c706e0b03"
)
EXPECTED_FAILED_V18_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "9fa7667e8076cee90a93fee44cfe08bc45a68e9ce28d01507b5c057463102a93"
)
EXPECTED_FAILED_V19_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "34d11ff6a699b96aaa8624b30f45fbc8845e6958dc04122b19c413a1a2c12c2d"
)
EXPECTED_REJECTED_V20_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "09794c2ce162af3bf2f3117f6d11dea0c4bd626cbe50946267609058c6c0c291"
)
EXPECTED_FAILED_V21_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "825111aa7dcd25e60c6357243e228422d2cf07855d682aae30d90ebdcd5559d2"
)
EXPECTED_REJECTED_V22_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "94ebaec5e5fc994dc75ffcf189f124726154fccc89ff1ea32add8398f06bf5b2"
)
EXPECTED_REJECTED_V23_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "d050c0dfea29b469e271967ea7759d4cbb36c3cbec493010163be5938e81e54a"
)
EXPECTED_REJECTED_V24_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "dcff9fbb753ebdda8525597354be18be7b3fa3d22a4f40a020cd9f677736cc8a"
)
EXPECTED_FAILED_V25_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "b0056c3f60d7a7204d782ac1cea31e1e9411200b371295e38b9afc3d5f67a1d1"
)
EXPECTED_FAILED_V26_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "b61e7f75ba3e418098c3c7f9fb19da0b261380d7985833eb6a59b99d6e1aeaee"
)
EXPECTED_REJECTED_V27_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "355979601138f8dca29534db58e3862f30f30b49526c01d10a97b01ba91c26f9"
)
EXPECTED_FAILED_V28_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "ded65f60581476b08c530c04d795c55bcb393558f393acb0f171629dce47add7"
)
EXPECTED_FAILED_V29_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "04d6acab01d56b0bfe25726a242904afd38bbc3ee47d4e3db29e9eb154e23e8b"
)
EXPECTED_RECOVERY_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "3ea7624ed42caba9bd51ade25a4c9a037f0b84689b4e2a3563c8205bbb136fcd"
)
EXPECTED_V7_PUBLIC_KEY_SHA256 = (
    "cecb2287cf87f8bd948c7390bd5f7059578966e0fc917de420572ddd82d8f21b"
)
EXPECTED_FINAL_RELEASE_PUBLIC_KEY_SHA256 = (
    "0d985d9afce029c06275b6932d51db050f807db50ed27c3bf66cd6f9e201f267"
)
EXPECTED_REPORT_RELEASE_PUBLIC_KEY_SHA256 = (
    "8eaf44b95e216320b45ca4109cf34b2c86b8fd4dcf06bb58025c1e27126bf5b5"
)
EXPECTED_RELEASE_SOURCE_AUTHORITY_SHA256 = (
    "36621dc735d17cf3322e34aeb5d72247869935e7b87306cdb600ff6d522c52af"
)
EXPECTED_SOURCE_RECEIPT_SHA256 = (
    "d9fd6ecfc92b1e76ef6b68932faaab1e3d8e8540f521e5a7065f760bdd228e19"
)
EXPECTED_SOURCE_RECEIPT_SIZE = 6_733
EXPECTED_RUNNER_SHA256 = "2f5ff29333fccb091e1f812f67cd0cbd8b6e835c804c95fbbe6d32a83a96c6ba"
EXPECTED_RUNNER_LINE_COUNT = 638
EXPECTED_RUNNER_PREFIX_SHA256 = (
    "a5409c1a25f11e7a85b1bc9305689e4a52cdd4e81281cef7fb6658780bb97814"
)
EXPECTED_RUNNER_PREFIX_SIZE = 22_591
EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SHA256 = (
    "92c6667a96ebb1c126242e8a8bd4bc7983993fbb5e3fac0fad3b2f6588d3c24d"
)
EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SIZE = 22_623
EXPECTED_FAILED_V28_RESULT_PUBLIC_KEY_SHA256 = (
    "5b7d5d026835ec6ec677bcd886bc16ac71444117dabc44a084ce3ede3a4db5a9"
)
EXPECTED_FAILED_V28_REPORT_PUBLIC_KEY_SHA256 = (
    "572e27167e1eea4c39ca53546ba33868b2874ec7fa3ed1682db821b2e50fa439"
)
EXPECTED_FAILED_V28_RESULT_ARCHIVE_RECORD_SHA256 = (
    "a25b458c0118248dff43fb7fbc77d14aea20b8df76b28305646aee7b917d9d26"
)
EXPECTED_FAILED_V28_REPORT_ARCHIVE_RECORD_SHA256 = (
    "5e9d9114233b742626a1d300aaf721bfe3c4a4527717d39db1a930a4fe8f3576"
)
RUNNER_ARCHIVE_REBASE = (
    (
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
        "20260719_all_ssnv_result_v5_post_reboot_bootstrap/ed25519_public.pem",
        str(FAILED_V28_RESULT_PUBLIC_KEY),
        "failed_v28_result_public_key",
    ),
    (
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
        "20260719_all_ssnv_report_v5_post_reboot_bootstrap/ed25519_public.pem",
        str(FAILED_V28_REPORT_PUBLIC_KEY),
        "failed_v28_report_public_key",
    ),
    (
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
        "20260719_all_ssnv_report_v5_post_reboot_bootstrap/"
        "ed25519_private_encrypted_unrecoverable_after_signing.pem",
        str(FAILED_V28_REPORT_PRIVATE_KEY),
        "failed_v28_report_private_key_metadata_only",
    ),
)
RUNNER_CURRENT_RELEASE_KEY_REBASE = (
    (
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
        "20260719_all_ssnv_result_v5_post_reboot_bootstrap/ed25519_public.pem",
        str(FINAL_RELEASE_PUBLIC_KEY),
        "current_result_public_key",
    ),
    (
        "5b7d5d026835ec6ec677bcd886bc16ac71444117dabc44a084ce3ede3a4db5a9",
        EXPECTED_FINAL_RELEASE_PUBLIC_KEY_SHA256,
        "current_result_public_key_sha256",
    ),
    (
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
        "20260719_all_ssnv_report_v5_post_reboot_bootstrap/ed25519_public.pem",
        str(REPORT_RELEASE_PUBLIC_KEY),
        "current_report_public_key",
    ),
    (
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
        "20260719_all_ssnv_report_v5_post_reboot_bootstrap/"
        "ed25519_private_encrypted_unrecoverable_after_signing.pem",
        str(REPORT_RELEASE_PRIVATE_KEY),
        "current_report_private_key",
    ),
    (
        "572e27167e1eea4c39ca53546ba33868b2874ec7fa3ed1682db821b2e50fa439",
        EXPECTED_REPORT_RELEASE_PUBLIC_KEY_SHA256,
        "current_report_public_key_sha256",
    ),
)
FAILED_V28_LIVE_KEY_PATHS = tuple(
    Path(original) for original, _, _ in RUNNER_ARCHIVE_REBASE
)
EXPECTED_RUNNER_LEAF_SHA256 = (
    "757ddf3fac80be9ca35d1de1eeaaa1ff70f5108a91916abdc2a17a6f94818077"
)
EXPECTED_V7_AUTHORITY_ID = (
    "20260722_all_ssnv_focal_alt_task_b_release_v7_strict_command_parity"
)
EXPECTED_V7_SOURCE_SET_SHA256 = (
    "f5ba8a9e786971c4261b51e283a0f9df6e807a8aca695bf5041271fee5420f58"
)
AUTHORIZATION_ID = "20260722_tumor_ref_receipt_promotion_v3"

AUTHORIZATION_SCHEMA = "intersubmod.tumor_ref_source_receipt_promotion_authorization"
COMPLETION_SCHEMA = "intersubmod.tumor_ref_source_receipt_promotion"
VERIFICATION_SCHEMA = (
    "intersubmod.tumor_ref_source_receipt_promotion_verification.recovery"
)
PREPARE_EXIT_ATTESTATION_SCHEMA = {
    "name": "intersubmod.tumor_ref_source_receipt_promotion.prepare_exit_attestation",
    "version": "1.0.0",
}
PROMOTE_EXIT_ATTESTATION_SCHEMA = {
    "name": "intersubmod.tumor_ref_source_receipt_promotion.promote_exit_attestation",
    "version": "1.0.0",
}
REPLAY_SCHEMA = "intersubmod.m2v5_runner_only_gate_replay.recovery"
CONTINUATION_SCHEMA = "intersubmod.m2v5_downstream_continuation.recovery"
SCHEMA_VERSION = "3.0.0"
RECOVERY_SCHEMA_VERSION = "1.0.0"

BASIC_ARTIFACT_KEYS = frozenset({"path", "size_bytes", "sha256", "mode"})
EXECUTABLE_ALIAS_RECORD_KEYS = frozenset(
    {"path", "resolved_path", "symlink_target", "size_bytes", "sha256", "mode"}
)
STAT_ARTIFACT_KEYS = frozenset(
    {
        "path",
        "size_bytes",
        "sha256",
        "mode",
        "device",
        "inode",
        "link_count",
        "mtime_ns",
        "ctime_ns",
    }
)
LEGACY_STAT_ARTIFACT_KEYS = frozenset(
    {
        "path",
        "size_bytes",
        "sha256",
        "mode",
        "device",
        "inode",
        "mtime_ns",
        "ctime_ns",
    }
)
METADATA_ARTIFACT_KEYS = frozenset({"path", "mode"})
METADATA_PLUS_SIZE_STATE_KEYS = frozenset(
    {"path", "size_bytes", "mode", "state"}
)
IDENTITY_MARKER_KEYS = frozenset(
    (
        set(BASIC_ARTIFACT_KEYS)
        | set(EXECUTABLE_ALIAS_RECORD_KEYS)
        | set(STAT_ARTIFACT_KEYS)
        | set(METADATA_ARTIFACT_KEYS)
    )
    - {"path"}
)
CONTINUATION_RECEIPT_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "task_type",
        "status",
        "scope",
        "canonical_command",
        "clean_environment",
        "supervisor_capability",
        "code",
        "legacy_signed_code",
        "schema_recovery_authority",
        "historical_identity_contract",
        "self_binding",
        "signed_promotion",
        "fresh_promotion_verification",
        "runner_only_replay",
        "authority_semantics",
        "v7_validation",
        "parent_gate_validation",
        "mutation_monitor",
        "composed_runner_stdin",
        "child_execution",
        "terminal_releases",
        "terminal_signature_contract",
        "governance",
        "checks",
        "pass",
        "pass_semantics",
    }
)
CONTINUATION_EXIT_ATTESTATION_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "task_type",
        "status",
        "scope",
        "supervisor_command",
        "supervised_child_command",
        "clean_environment",
        "supervisor_source",
        "supervisor_source_binding",
        "python_runtime",
        "supervisor_capability",
        "child_wait",
        "execution_receipt",
        "execution_receipt_pass_semantics",
        "terminal_signature_contract",
        "checks",
        "pass",
        "pass_semantics",
    }
)
CONTINUATION_SUCCESS_WITNESS_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "task_type",
        "status",
        "scope",
        "supervisor_command",
        "prewitness_verifier_command",
        "continuation_source",
        "continuation_source_binding",
        "execution_receipt",
        "exit_attestation",
        "continuation_signature",
        "continuation_public_key",
        "retired_continuation_private_key",
        "supervisor_capability",
        "supervised_child_wait",
        "complete_child_output_inventory",
        "prewitness_verification",
        "checks",
        "release_authority_granted",
        "pass",
        "pass_semantics",
    }
)
SUPERVISOR_CAPABILITY_RECORD_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "protocol",
        "token_size_bytes",
        "token_sha256",
        "nonce_sha256",
        "supervisor_pid",
        "supervisor_start_ticks",
        "supervisor_command",
        "supervised_child_command",
        "supervisor_source",
        "python_runtime",
        "pass",
    }
)
SUPERVISOR_CAPABILITY_TOKEN_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "nonce",
        "supervisor_pid",
        "supervisor_start_ticks",
        "supervisor_command",
        "supervised_child_command",
        "supervisor_source",
        "python_runtime",
    }
)
PRIVATE_KEY_PRE_SIGNATURE_KEYS = frozenset(
    {
        "path",
        "mode",
        "device",
        "inode",
        "link_count",
        "mtime_ns",
        "ctime_ns",
    }
)
AUTHORIZATION_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "authorization_id",
        "task_type",
        "status",
        "scope",
        "commands",
        "reviewed_sources",
        "execution_binding",
        "producer_session",
        "producer_exit_attestation",
        "review_artifacts",
        "trusted_signing_keys",
        "source_receipt",
        "canonical_target_path",
        "evidence",
        "authorization_signature_contract",
        "completion_signature_contract",
        "continuation_gate",
        "historical_incident_disclosure",
        "downstream_output_absence",
        "checks",
        "pass",
        "pass_semantics",
    }
)
COMPLETION_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "authorization_id",
        "task_type",
        "status",
        "scope",
        "command",
        "authorization",
        "code",
        "execution_binding",
        "producer_session",
        "producer_exit_attestation",
        "source_receipt",
        "canonical_receipt",
        "evidence",
        "review_artifacts",
        "builder_gate_replay",
        "historical_incident_disclosure",
        "governance",
        "signature_contract",
        "checks",
        "pass",
        "pass_semantics",
    }
)
VERIFICATION_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "task_type",
        "mode",
        "authorization",
        "completion",
        "source_receipt",
        "canonical_receipt",
        "trusted_signing_keys",
        "private_key_retirement",
        "producer_exit_attestations",
        "historical_incident_disclosure",
        "release_source_authority_validator",
        "schema_recovery_authority",
        "historical_identity_contract",
        "checks",
        "retained_critical_input_descriptor_count",
        "governance",
        "pass",
        "pass_semantics",
    }
)
PRODUCER_EXIT_ATTESTATION_KEYS = frozenset(
    {
        "authorization_id",
        "phase",
        "supervisor_command",
        "producer_session",
        "child_wait",
        "produced_artifacts",
        "producer_source",
        "python_runtime",
        "signature_contract",
        "checks",
        "pass",
        "pass_semantics",
        "created_at_utc",
    }
)
PRODUCER_SESSION_KEYS = frozenset(
    {
        "protocol",
        "schema_version",
        "mode",
        "nonce_sha256",
        "parent_pid",
        "parent_start_time_ticks",
        "child_pid",
        "child_start_time_ticks",
        "fork_no_exec",
    }
)
PRODUCER_CHILD_WAIT_KEYS = frozenset(
    {
        "raw_wait_status",
        "waited_pid",
        "wifexited",
        "exit_status",
        "wifsignaled",
        "terminating_signal",
        "core_dumped",
    }
)
PREPARE_EXIT_CHECKS = {
    "parent_forked_child_without_exec": True,
    "parent_waitpid_returned_exact_child": True,
    "child_wait_was_normal_exit_zero": True,
    "producer_session_exactly_reopened": True,
    "authorization_reopened_with_full_stat_identity": True,
    "authorization_payload_strictly_revalidated": True,
    "legacy_authorization_payload_signature_absent": True,
    "legacy_completion_payload_signature_absent": True,
    "attestation_publication_is_atomic_no_replace": True,
}
PROMOTE_EXIT_CHECKS = {
    "parent_forked_child_without_exec": True,
    "parent_waitpid_returned_exact_child": True,
    "child_wait_was_normal_exit_zero": True,
    "producer_session_exactly_reopened": True,
    "prepare_attestation_signature_independently_verified": True,
    "prepare_attestation_authorization_binding_independently_verified": True,
    "completion_receipt_reopened_with_full_stat_identity": True,
    "canonical_receipt_reopened_with_full_stat_identity": True,
    "completion_payload_strictly_revalidated": True,
    "legacy_authorization_payload_signature_absent": True,
    "legacy_completion_payload_signature_absent": True,
    "attestation_publication_is_atomic_no_replace": True,
}
PREPARE_EXIT_PASS_SEMANTICS = (
    "parent_observed_fork_no_exec_child_normal_exit_zero_and_strictly_revalidated_"
    "the_authorization_output; this_attestation_requires_its_detached_authorization_"
    "key_signature_before_promotion"
)
PROMOTE_EXIT_PASS_SEMANTICS = (
    "parent_observed_fork_no_exec_child_normal_exit_zero_and_strictly_revalidated_"
    "the_completion_and_canonical_outputs_after_independently_verifying_the_signed_"
    "prepare_attestation; downstream_authority_still_requires_the_detached_completion_"
    "key_signature_and_all_continuation_gates"
)
REPLAY_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "task_type",
        "status",
        "scope",
        "canonical_command",
        "clean_environment",
        "code",
        "legacy_signed_code",
        "schema_recovery_authority",
        "historical_identity_contract",
        "failed_v28_key_archive_rebase_inputs",
        "runner_archive_rebase_contract",
        "runtime",
        "promotion_verifier",
        "promotion_trust_chain",
        "v7_validation",
        "completion_runner",
        "runner_only_replay",
        "output_slot_checks",
        "failed_v28_live_key_path_checks",
        "log",
        "claims",
        "downstream_not_executed",
        "downstream_authorized_after_this_gate",
        "pass",
        "pass_semantics",
    }
)
DATASET_RELEASE_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "task_type",
        "scope",
        "release_status",
        "pass_semantics",
        "command",
        "source_authority",
        "inputs",
        "code",
        "source_modes",
        "result_signature_contract",
        "builder_counts",
        "checks",
        "pass",
    }
)
REPORT_RELEASE_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "task_type",
        "scope",
        "release_status",
        "pass_semantics",
        "command",
        "source_authority",
        "inputs",
        "code",
        "source_modes",
        "report_signature_contract",
        "checks",
        "pass",
    }
)

AUTHORIZATION_CHECKS = {
    "three_distinct_review_artifacts_exact_source_bound": True,
    "producer_session_fork_no_exec_contract_recorded": True,
    "authorization_public_key_precommitted": True,
    "completion_public_key_precommitted": True,
    "authorization_private_key_pre_signature_mode_0400": True,
    "authorization_payload_signature_forbidden": True,
    "prepare_exit_attestation_signature_required": True,
    "completion_private_key_pre_signature_mode_0400": True,
    "continuation_public_key_precommitted": True,
    "continuation_private_key_pre_signature_mode_0400": True,
    "source_receipt_exact_and_builder_accepted": True,
    "historical_incident_and_prior_policy_hash_bound": True,
    "failed_attempt_archive_independently_replayed": True,
    "v3_transitive_signature_binding_verified": True,
    "v7_full_source_authority_verified": True,
    "all_downstream_output_slots_absent": True,
    "scientific_payload_changed": False,
}
COMPLETION_CHECKS = {
    "signed_prepare_exit_attestation_verified_over_bound_bytes": True,
    "prepare_attestation_authorization_binding_verified": True,
    "authorization_private_key_retired_mode_000": True,
    "completion_private_key_pre_signature_mode_0400": True,
    "completion_payload_signature_forbidden": True,
    "promote_exit_attestation_signature_required": True,
    "three_distinct_review_artifacts_exact_source_bound": True,
    "source_receipt_bound_no_follow": True,
    "canonical_staging_regular_mode_0444_link_count_zero": True,
    "completion_receipt_published_before_canonical_commit": True,
    "write_ahead_receipt_is_durable_incomplete_transaction_marker": True,
    "canonical_target_created_exclusively_after_completion_receipt": True,
    "canonical_target_regular_file": True,
    "canonical_target_link_count_one": True,
    "canonical_target_parent_fsynced_after_link": True,
    "canonical_target_post_link_path_inode_rebound": True,
    "source_and_target_distinct_inodes": True,
    "source_and_target_bytes_equal": True,
    "source_and_target_size_equal": True,
    "source_and_target_sha256_equal": True,
    "source_and_target_paths_bound_through_canonical_commit": True,
    "v3_transitive_signature_binding_verified": True,
    "v7_full_source_authority_verified": True,
    "final_builder_gate_passed_on_source_and_staged_target": True,
    "failed_attempt_evidence_replayed_and_preserved": True,
    "all_downstream_output_slots_absent": True,
    "completion_runner_end_to_end_pass_claimed": False,
    "original_verifier_success_claimed": False,
    "scientific_payload_changed": False,
}
DATASET_RELEASE_CHECKS = frozenset(
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
REPORT_RELEASE_CHECKS = frozenset(
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

CONTINUATION_POLICY = (
    "No strict, matched-normal, CN/CCF, final-dataset, or final-report command "
    "may run until both detached promotion signatures verify, the promotion verification "
    "receipt is recorded, the exact original runner lines 1-358 are bound and a "
    "deterministic three-path archive-rebased derivative executes successfully as "
    "non-authoritative evidence, and the signed downstream continuation runner "
    "independently revalidates every gate input from bound descriptors and proves an "
    "inherited sealed-memfd capability from its live supervisor parent. Final "
    "release authority additionally requires a machine-generated supervisor "
    "attestation that wait status was a normal exit zero, a detached signature "
    "over that attestation, independent signed-terminal verification, and the explicit "
    "trusted-account threat model recorded by the signed authorization."
)
AUTHORIZED_CONTINUATION_POLICY = (
    "No strict, matched-normal, CN/CCF, final-dataset, or final-report command "
    "may run until both detached promotion signatures verify, the promotion verification "
    "receipt is recorded, the original runner lines 1-358 replay successfully as "
    "non-authoritative evidence, and the signed downstream continuation runner "
    "independently revalidates every gate input from bound descriptors and proves an "
    "inherited sealed-memfd capability from its live supervisor parent. Final "
    "release authority additionally requires a machine-generated supervisor "
    "attestation that wait status was a normal exit zero, a detached signature "
    "over that attestation, independent signed-terminal verification, and the explicit "
    "trusted-account threat model recorded by the signed authorization."
)
TERMINAL_AUTHORITY_SEMANTICS = {
    "runner_only_replay_is_non_authoritative_evidence": True,
    "runner_only_replay_downstream_authorized_after_this_gate": False,
    "actual_execution_authority": (
        "signed_prepare_and_promote_exit_attestations_plus_payload_bindings_"
        "plus_fresh_verifier_and_"
        "bound_parent_gate_validation_and_sealed_parent_capability"
    ),
    "threat_model": CONTINUATION_THREAT_MODEL,
}
TERMINAL_GOVERNANCE = {
    "historical_runner_lines_359_402_executed": False,
    "canonical_source_receipt_recreated": False,
    "canonical_source_receipt_revalidated_from_signed_record": True,
    "shell_pathname_gate_used_as_authority": False,
    "all_gate_input_descriptors_retained_until_process_exit": True,
    "transient_input_path_mutation_monitor_retained_until_process_exit": True,
    "promotion_incident_absent_at_every_terminal_gate": True,
    "terminal_receipt_is_provisional_no_replace_validation_marker": True,
    "supervisor_wait_status_zero_before_exit_attestation_required": True,
    "sealed_memfd_parent_capability_required_before_downstream_execution": True,
    "detached_terminal_signature_covers_machine_exit_attestation": True,
    "published_terminal_signature_is_provisional_until_independent_verifier": True,
    "post_link_failure_may_leave_provisional_signature": True,
    "post_link_operations_do_not_mutate_signed_payload": True,
    "lease_holding_supervisor_failure_requires_continuation_incident": True,
    "durable_supervisor_success_witness_required_for_release_authority": True,
    "success_witness_published_only_after_precommit_verifier_passes": True,
    "success_witness_link_is_final_no_return_supervisor_commit": True,
    "final_verifier_requires_writer_lease_release": True,
    "signature_without_success_witness_never_grants_release_authority": True,
    "unsigned_terminal_receipt_grants_release_authority": False,
    "signed_dataset_and_report_receipts_are_provisional_without_this_terminal_receipt": True,
    "continuation_incident_absence_is_required_for_release_authority": True,
    "fallible_operations_after_success_witness_link": 0,
    "partial_outputs_claimed_complete": False,
    "threat_model": CONTINUATION_THREAT_MODEL,
}
TERMINAL_CHECKS = {
    "executing_source_matches_bound_recursive_module_code": True,
    "function_class_method_and_critical_global_binding_exact": True,
    "signed_promotion_authorization_exact": True,
    "signed_promotion_completion_exact": True,
    "schema_recovery_authority_verified": True,
    "legacy_signed_bindings_preserved": True,
    "prepare_producer_exit_attestation_signed_and_exact": True,
    "promote_producer_exit_attestation_signed_and_exact": True,
    "promotion_private_keys_retired_mode_000": True,
    "fresh_and_recorded_promotion_verifier_agree": True,
    "runner_replay_consumed_as_non_authoritative_evidence": True,
    "runner_replay_downstream_authority_false": True,
    "runner_gate_inputs_independently_bound_and_validated": True,
    "transient_path_replacement_detected_fail_closed": True,
    "cooccurrence_core_predicates_independently_recounted": True,
    "static_bridge_skipped_only_lines_359_402": True,
    "supervisor_capability_memfd_bound_to_live_parent": True,
    "downstream_child_exit_zero": True,
    "dataset_release_signature_and_relations_reverified": True,
    "report_release_signature_and_relations_reverified": True,
    "all_declared_release_artifacts_mode_0444_link_count_one": True,
    "complete_child_output_tree_inventory_bound_and_monitored": True,
    "final_release_private_keys_retired_and_rebound_at_terminal_precommit": True,
    "all_source_and_input_paths_stable_before_terminal_commit": True,
}
TERMINAL_PASS_SEMANTICS = (
    "provisional_terminal_execution_and_release_chain_validation_only; machine "
    "supervisor exit-zero attestation, detached continuation signature, retired "
    "terminal private key, durable supervisor-success witness, exact live artifacts, "
    "and absence of both incidents remain mandatory; not scientific confirmation"
)
EXIT_ATTESTATION_CHECKS = {
    "supervisor_source_and_interpreter_bound_before_child_launch": True,
    "sealed_memfd_capability_bound_to_supervisor_and_child_receipt": True,
    "supervised_child_started_by_this_process": True,
    "wait_status_observed_by_parent_process": True,
    "child_exited_normally_with_status_zero": True,
    "execution_receipt_absent_before_child_launch": True,
    "execution_receipt_published_by_child_no_replace": True,
    "execution_receipt_full_nested_contract_reverified_after_wait": True,
    "all_live_release_artifacts_reverified_after_wait": True,
    "continuation_incident_absent_before_attestation_commit": True,
    "attestation_signature_target_absent_before_signing": True,
    "persistent_private_key_retired_before_signature_bytes_exist": True,
    "sealed_unlinked_memfd_signing_key_closed_before_publication": True,
    "published_signature_provisional_until_independent_verification": True,
    "post_link_supervisor_failure_requires_continuation_incident": True,
    "durable_supervisor_success_witness_required_for_release_authority": True,
}
EXIT_ATTESTATION_PASS_SEMANTICS = (
    "machine_supervisor_wait_exit_zero_attestation_integrity_only; the detached "
    "signature remains provisional even after no-replace publication; retired one-time "
    "private key, exact live artifacts, independent signed-terminal verification, and "
    "incident absence remain mandatory; not scientific confirmation"
)
SUCCESS_WITNESS_CHECKS = {
    "signed_terminal_precommit_verifier_exit_zero": True,
    "prewitness_verifier_granted_no_release_authority": True,
    "execution_receipt_exit_attestation_and_signature_exact": True,
    "continuation_private_key_retired_mode_000": True,
    "complete_child_output_inventory_exact": True,
    "all_bound_inputs_and_release_artifacts_stable_before_commit": True,
    "continuation_incident_absent_before_commit": True,
    "witness_target_absent_before_no_replace_commit": True,
    "catchable_signals_blocked_before_final_witness_link": True,
    "witness_link_is_final_no_return_supervisor_commit_point": True,
    "final_verifier_requires_released_writer_lease": True,
}
SUCCESS_WITNESS_PASS_SEMANTICS = (
    "durable_clean_supervisor_completion_witness_only; independent final "
    "signed-terminal verification and incident absence remain mandatory; "
    "not scientific confirmation"
)
PREWITNESS_VERIFICATION_PASS_SEMANTICS = (
    "signed_terminal_precommit_integrity_only_no_release_authority; durable "
    "supervisor success witness is still absent and mandatory"
)
REPLAY_PASS_SEMANTICS = (
    "non_authoritative_archive_rebased_pathname_replay_evidence_only; signed downstream "
    "continuation must independently revalidate all gate inputs from bound "
    "descriptors before any downstream command"
)
VERIFICATION_CHECKS = {
    "clean_direct_runtime": True,
    "executing_verifier_definitions_match_bound_reviewed_source": True,
    "authority_module_executed_from_bound_source_bytes": True,
    "v7_release_source_authority_verified_with_process_lifetime_fd_lease": True,
    "authorization_schema_and_exact_keys_verified": True,
    "authorization_signature_verified_from_bound_fds": True,
    "prepare_exit_attestation_exact_schema_verified": True,
    "prepare_child_waitpid_normal_exit_zero_verified": True,
    "authorization_signature_verified_over_prepare_exit_attestation": True,
    "legacy_authorization_payload_signature_absent": True,
    "authorization_private_key_retired_mode_000": True,
    "completion_schema_and_exact_keys_verified": True,
    "completion_signature_verified_from_bound_fds": True,
    "promote_exit_attestation_exact_schema_verified": True,
    "promote_child_waitpid_normal_exit_zero_verified": True,
    "completion_signature_verified_over_promote_exit_attestation": True,
    "legacy_completion_payload_signature_absent": True,
    "promote_attestation_binds_signed_prepare_chain": True,
    "completion_private_key_retired_mode_000": True,
    "three_distinct_review_artifacts_verified": True,
    "four_core_plus_eleven_runtime_sources_verified": True,
    "portable_reader_asset_membership_exact": True,
    "source_and_canonical_bytes_size_sha256_equal": True,
    "source_and_canonical_independent_inodes": True,
    "canonical_regular_mode_0444_link_count_one": True,
    "historical_incident_disclosure_hash_bound": True,
    "exact_artifact_cross_links_verified": True,
    "critical_artifact_paths_rebound_before_success": True,
    "critical_descriptor_leases_retained_until_process_exit": True,
    "signed_schema_recovery_authority_verified": True,
    "historical_eight_key_identity_validated_without_inferred_link_count": True,
    "current_nine_key_identity_and_link_count_one_verified": True,
}
VERIFICATION_GOVERNANCE = {
    "original_completion_runner_end_to_end_pass_claimed": False,
    "original_verifier_reexecuted_successfully_claimed": False,
    "runner_only_gate_replay_evidence_required_before_downstream": True,
    "signed_downstream_continuation_revalidation_required": True,
    "verification_scope": "promotion_integrity_and_continuation_gate_only",
    "schema_recovery_scope": "historical_identity_schema_compatibility_only",
    "biological_claim_ceiling_changed": False,
}
VERIFICATION_PASS_SEMANTICS = (
    "signed byte-identical promotion chain verified through the scoped historical "
    "eight-key schema recovery; link_count was not recorded and was not inferred; "
    "this does not claim the original completion runner or verifier succeeded"
)

PROMOTION_VERIFIER_KEYS = frozenset(
    {
        "legacy_signed_command",
        "recovery_authorized_command",
        "actual_argv",
        "normalized_argv",
        "pass_fds",
        "fd_artifact_binding",
        "source_descriptor_bound",
        "exit_code",
        "stdout",
        "stderr",
        "verification_receipt",
        "fresh_verification_mode",
        "recorded_verification_mode",
        "fresh_and_recorded_evidence_agree",
        "pass",
    }
)
PROMOTION_TRUST_CHAIN_KEYS = frozenset(
    {
        "authorization",
        "prepare_exit_attestation",
        "authorization_signature",
        "authorization_signature_target",
        "authorization_public_key",
        "retired_authorization_private_key",
        "promotion_receipt",
        "promote_exit_attestation",
        "promotion_signature",
        "promotion_signature_target",
        "completion_public_key",
        "retired_completion_private_key",
        "continuation_public_key",
        "pre_signature_continuation_private_key",
        "historical_source_receipt",
        "canonical_source_receipt",
        "reviewed_sources",
        "portable_reader_asset_membership",
        "authorization_signature_verified",
        "completion_signature_verified",
        "prepare_child_waitpid_normal_exit_zero_verified",
        "promote_child_waitpid_normal_exit_zero_verified",
        "legacy_payload_signatures_absent",
        "promote_attestation_binds_signed_prepare_chain",
        "canonical_bytes_equal_historical_source",
        "canonical_inode_distinct_from_historical_source",
    }
)
REPLAY_RUNTIME_KEYS = frozenset({"python", "bash", "openssl"})
REPLAY_SUCCESS_WITNESS_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "task_type",
        "status",
        "scope",
        "canonical_command",
        "clean_environment",
        "supervisor_source",
        "worker_wait",
        "worker_outputs",
        "checks",
        "pass",
        "pass_semantics",
    }
)
REPLAY_SUCCESS_WITNESS_CHECKS = {
    "supervisor_source_descriptor_bound": True,
    "worker_waitpid_normal_exit_zero": True,
    "worker_receipt_source_matches_supervisor_source": True,
    "worker_log_and_receipt_immutable_and_bound": True,
    "worker_receipt_pass": True,
}
REPLAY_SUCCESS_WITNESS_PASS_SEMANTICS = (
    "supervisor_observed_worker_waitpid_normal_exit_zero_and_bound_the_"
    "immutable_replay_log_receipt_and_reviewed_worker_source"
)
REPLAY_COMPLETION_RUNNER_KEYS = frozenset(
    {
        "artifact",
        "observed_total_physical_lines",
        "included_line_start",
        "included_line_end",
        "first_excluded_line",
        "first_excluded_line_is_blank",
        "first_downstream_executable_line",
        "first_downstream_executable_line_prefix",
        "stdin_size_bytes",
        "stdin_sha256",
        "executed_stdin_size_bytes",
        "executed_stdin_sha256",
        "authorization_bound_runtime_contract",
    }
)
REPLAY_RUNNER_ONLY_KEYS = frozenset(
    {
        "command",
        "stdin",
        "exit_code",
        "stdout",
        "stderr",
        "shell_observation_checks",
        "authoritative_gate_validation",
        "pass",
    }
)
REPLAY_OUTPUT_SLOT_KEYS = frozenset(
    {"before_verifier", "after_verifier", "after_runner_replay", "all_absent"}
)
REPLAY_CLAIM_KEYS = frozenset(
    {
        "original_completion_runner_end_to_end_pass_claimed",
        "original_tumor_ref_verifier_successfully_reexecuted_claimed",
        "promotion_integrity_verified",
        "archive_rebased_runner_lines_1_358_executed_successfully",
        "runner_prefix_archive_rebased_three_paths",
        "failed_v28_private_key_bytes_read",
        "failed_v28_live_key_paths_created",
        "shell_pathname_gate_result_is_authoritative",
        "runner_line_359_or_later_executed",
        "strict_executed",
        "matched_normal_executed",
        "cn_ccf_executed",
        "final_dataset_or_report_executed",
        "downstream_authorized_by_this_receipt",
    }
)

# Linux inotify masks used to retain evidence of transient path replacement.
IN_MODIFY = 0x00000002
IN_ATTRIB = 0x00000004
IN_CLOSE_WRITE = 0x00000008
IN_MOVED_FROM = 0x00000040
IN_MOVED_TO = 0x00000080
IN_CREATE = 0x00000100
IN_DELETE = 0x00000200
IN_DELETE_SELF = 0x00000400
IN_MOVE_SELF = 0x00000800
IN_UNMOUNT = 0x00002000
IN_Q_OVERFLOW = 0x00004000
IN_IGNORED = 0x00008000
INOTIFY_MUTATION_MASK = (
    IN_MODIFY
    | IN_ATTRIB
    | IN_CLOSE_WRITE
    | IN_MOVED_FROM
    | IN_MOVED_TO
    | IN_CREATE
    | IN_DELETE
    | IN_DELETE_SELF
    | IN_MOVE_SELF
    | IN_UNMOUNT
    | IN_Q_OVERFLOW
    | IN_IGNORED
)

PRIMARY_PYTHON_WRAPPER = "run_bound_primary_python"
QA_PYTHON_WRAPPER = "run_bound_qa_python"
RUNTIME_FD_BRIDGE_TEMPLATE = b'''\

# Execute immutable descriptor bytes while retaining the signed canonical argv[0].
readonly INTERSUBMOD_PRIMARY_PYTHON_FD="__PRIMARY_PYTHON_FD__"
readonly INTERSUBMOD_QA_PYTHON_FD="__QA_PYTHON_FD__"

run_bound_primary_python() (
    exec -a "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python" \
        "/proc/self/fd/${INTERSUBMOD_PRIMARY_PYTHON_FD}" "$@"
)

run_bound_qa_python() (
    exec -a "/bip7_disk/liaoyoyo2001/miniconda3/bin/python" \
        "/proc/self/fd/${INTERSUBMOD_QA_PYTHON_FD}" "$@"
)

'''

BRIDGE_BYTES = b'''\

# Reviewed continuation bridge: canonical tumor-REF identity already exists.
readonly PROMOTED_SOURCE_RECEIPT_PATH="/big7_disk/liaoyoyo2001/big7_disk_output/"\
"synthesis/observation_workspaces/"\
"20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"\
"tumor_ref_recovery_source_identity_v1/post_run_source_identity.receipt.json"
readonly PROMOTED_SOURCE_RECEIPT_SHA256="d9fd6ecfc92b1e76ef6b68932faaab1e3d8e8540f521e5a7065f760bdd228e19"
require_file "${SOURCE_RECEIPT}"
[[ "$(realpath "${SOURCE_RECEIPT}")" == "${PROMOTED_SOURCE_RECEIPT_PATH}" ]]
[[ "$(stat -c '%s' "${SOURCE_RECEIPT}")" == "6733" ]]
[[ "$(sha256sum "${SOURCE_RECEIPT}" | awk '{print $1}')" == "${PROMOTED_SOURCE_RECEIPT_SHA256}" ]]
[[ "$(stat -c '%a' "${SOURCE_RECEIPT}")" == "444" ]]
require_json_pass "${SOURCE_RECEIPT}"

for path in \
    "${STRICT}" \
    "${MATCHED_RUN}" \
    "${MATCHED_ANALYSIS}" \
    "${CN_CCF}" \
    "${PRIMARY_POST}" \
    "${FROZEN_POST}" \
    "${FINAL_DATASET}" \
    "${REPORT}" \
    "${FINAL_RELEASE_RECEIPT}" \
    "${FINAL_RELEASE_SIGNATURE}" \
    "${REPORT_RELEASE_RECEIPT}" \
    "${REPORT_RELEASE_SIGNATURE}" \
    "${LOG_PATH}" \
    "${PYTHON_CACHE_ROOT}"; do
    require_absent "${path}"
done

umask 0077
/usr/bin/mkdir --mode=0700 -- "${PYTHON_CACHE_ROOT}"
[[ "$(stat -c '%a' "${PYTHON_CACHE_ROOT}")" == "700" ]]
umask 0022

exec > >(tee "${LOG_PATH}") 2>&1

printf 'Input manifest: %s\n' "${MANIFEST}"
printf 'Input screen: %s\n' "${SCREEN}"
printf 'Input cooccurrence: %s\n' "${COOCCURRENCE}"
printf 'Input recovered tumor-REF: %s\n' "${TUMOR_REF}"
printf 'Input independent M2 audit: %s\n' "${INDEPENDENT_M2_AUDIT}"
printf 'Input cooccurrence runtime preflight: %s\n' "${COOCCURRENCE_PREFLIGHT}"
printf 'Promoted source identity receipt: %s\n' "${SOURCE_RECEIPT}"
printf 'Output final dataset: %s\n' "${FINAL_DATASET}"
printf 'Output final report: %s\n' "${REPORT}"

'''

EPILOGUE_BYTES = b'''\

/usr/bin/chmod 0444 -- "${LOG_PATH}"
[[ "$(stat -c '%a' "${LOG_PATH}")" == "444" ]]
'''

GATE_INPUT_PATHS = {
    "manifest": MANIFEST,
    "screen_sites": SCREEN_SITES,
    "screen_assignments": SCREEN_ASSIGNMENTS,
    "screen_summary": SCREEN_SUMMARY,
    "primary_pre": PRIMARY_PRE,
    "legacy_primary_pre": LEGACY_PRIMARY_PRE,
    "cooccurrence_site": COOCCURRENCE / "methyl_ssnv_site_results.tsv.gz",
    "cooccurrence_pair": COOCCURRENCE / "methyl_ssnv_pair_results.tsv.gz",
    "cooccurrence_duplicate_audit": COOCCURRENCE / "raw_identity_duplicate_audit.tsv.gz",
    "cooccurrence_oracle": COOCCURRENCE / "oracle_cases.json",
    "cooccurrence_summary": COOCCURRENCE / "summary.json",
    "cooccurrence_run_receipt": COOCCURRENCE / "run_receipt.json",
    "cooccurrence_release_receipt": COOCCURRENCE / "release_receipt.json",
    "tumor_ref_sites": TUMOR_REF / "all_ssnv_tumor_ref_control_site_results.tsv.gz",
    "tumor_ref_summary": TUMOR_REF_SUMMARY,
    "tumor_ref_manifest": TUMOR_REF / "run_manifest.json",
    "source_snapshot": SOURCE_SNAPSHOT,
    "independent_m2_audit": INDEPENDENT_M2_AUDIT,
    "cooccurrence_preflight": COOCCURRENCE_PREFLIGHT,
    "claim_contract": CLAIM_CONTRACT,
    "v7_authority": V7_AUTHORITY,
    "v7_approval": V7_APPROVAL,
    "v7_signature": V7_SIGNATURE,
    "v7_public_key": V7_PUBLIC_KEY,
    "v7_private_key": V7_PRIVATE_KEY,
    "recovery_final_dataset_builder": FINAL_DATASET_BUILDER,
    "recovery_result_report_finalizer": FINAL_RELEASE_FINALIZER,
    "recovery_report_builder": REPORT_BUILDER,
    "final_release_public_key": FINAL_RELEASE_PUBLIC_KEY,
    "report_release_public_key": REPORT_RELEASE_PUBLIC_KEY,
    "earlier_fp_report": EARLIER_FP_REPORT,
    "intersubmod_binary": INTERSUBMOD_BIN,
    "reference": REFERENCE,
    "reference_fai": Path(str(REFERENCE) + ".fai"),
    "primary_python_runtime": PRIMARY_PYTHON_RUNTIME,
    "qa_python_runtime": QA_PYTHON_RUNTIME,
    "node": NODE,
    "qa_chromium": QA_CHROMIUM,
    "portable_builder_module": PORTABLE_BUILDER_MODULE,
    "portable_chart_extractor_module": PORTABLE_CHART_EXTRACTOR_MODULE,
    "portable_verifier_module": PORTABLE_VERIFIER_MODULE,
    "portable_browser_helpers_module": PORTABLE_BROWSER_HELPERS_MODULE,
    "portable_browser_cli_module": PORTABLE_BROWSER_CLI_MODULE,
    "portable_server_bundle": PORTABLE_SERVER_BUNDLE,
    **{
        f"portable_reader_asset_part{index:03d}": path
        for index, path in enumerate(PORTABLE_READER_ASSET_PARTS, start=1)
    },
}

AUTHORIZATION_EVIDENCE_KEYS = frozenset(
    {
        "source_snapshot",
        "tumor_ref_manifest",
        "focal_source_identity_transition",
        "prior_policy_review",
        "failed_attempt_archive",
        "historical_v3_authority",
        "current_v7_authority",
        "current_v7_runtime_validation",
        "final_builder",
        "final_builder_gate_execution",
        "python_runtime",
        "python_cache_contract",
        "final_builder_source_path_gate",
    }
)

CRITICAL_GLOBAL_NAMES = (
    "REPO_ROOT",
    "TOPIC_ROOT",
    "AUDIT_ROOT",
    "SCRIPT_ROOT",
    "RESULT_ROOT",
    "WORKSPACE_ROOT",
    "PYTHON",
    "QA_PYTHON",
    "PRIMARY_PYTHON_RUNTIME",
    "QA_PYTHON_RUNTIME",
    "BASH",
    "OPENSSL",
    "NODE",
    "QA_CHROMIUM",
    "CONTINUATION_RUNNER",
    "AUTHORIZED_CONTINUATION_RUNNER",
    "PROMOTION_TOOL",
    "CONTINUATION_VERIFIER",
    "AUTHORIZED_CONTINUATION_VERIFIER",
    "RUNNER_GATE_REPLAY",
    "AUTHORIZED_RUNNER_GATE_REPLAY",
    "SCHEMA_RECOVERY_VALIDATOR",
    "EXPECTED_SCHEMA_RECOVERY_VALIDATOR_SHA256",
    "RECOVERY_OUTPUT_LEASE_PUBLIC_KEY",
    "EXPECTED_RECOVERY_OUTPUT_LEASE_PUBLIC_KEY_SHA256",
    "COMPLETION_RUNNER",
    "RELEASE_SOURCE_AUTHORITY",
    "PORTABLE_PLUGIN_ROOT",
    "PORTABLE_PLUGIN_SCRIPT_ROOT",
    "PORTABLE_BUILDER_MODULE",
    "PORTABLE_CHART_EXTRACTOR_MODULE",
    "PORTABLE_VERIFIER_MODULE",
    "PORTABLE_BROWSER_HELPERS_MODULE",
    "PORTABLE_BROWSER_CLI_MODULE",
    "PORTABLE_SERVER_BUNDLE",
    "PORTABLE_ASSET_ROOT",
    "PORTABLE_READER_ASSET_PARTS",
    "AUTHORIZATION",
    "PREPARE_EXIT_ATTESTATION",
    "AUTHORIZATION_SIGNATURE",
    "LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE",
    "COMPLETION",
    "PROMOTE_EXIT_ATTESTATION",
    "COMPLETION_SIGNATURE",
    "LEGACY_COMPLETION_PAYLOAD_SIGNATURE",
    "PROMOTION_VERIFICATION_RECEIPT",
    "AUTHORIZED_PROMOTION_VERIFICATION_RECEIPT",
    "RUNNER_GATE_RECEIPT",
    "AUTHORIZED_RUNNER_GATE_RECEIPT",
    "RUNNER_GATE_LOG",
    "RUNNER_GATE_SUCCESS_WITNESS",
    "RUNNER_GATE_SUCCESS_WITNESS_SCHEMA",
    "PROMOTION_INCIDENT",
    "CONTINUATION_RECEIPT",
    "CONTINUATION_EXIT_ATTESTATION",
    "CONTINUATION_SIGNATURE",
    "CONTINUATION_SUCCESS_WITNESS",
    "CONTINUATION_INCIDENT",
    "AUTHORIZED_CONTINUATION_RECEIPT",
    "AUTHORIZED_CONTINUATION_EXIT_ATTESTATION",
    "AUTHORIZED_CONTINUATION_SIGNATURE",
    "HISTORICAL_SOURCE_RECEIPT",
    "CANONICAL_SOURCE_RECEIPT",
    "AUTHORIZATION_PUBLIC_KEY",
    "AUTHORIZATION_PRIVATE_KEY",
    "COMPLETION_PUBLIC_KEY",
    "COMPLETION_PRIVATE_KEY",
    "CONTINUATION_PUBLIC_KEY",
    "CONTINUATION_PRIVATE_KEY",
    "FAILED_V16_CONTINUATION_KEY_ROOT",
    "FAILED_V16_CONTINUATION_PUBLIC_KEY",
    "FAILED_V16_CONTINUATION_PRIVATE_KEY",
    "FAILED_V17_CONTINUATION_KEY_ROOT",
    "FAILED_V17_CONTINUATION_PUBLIC_KEY",
    "FAILED_V17_CONTINUATION_PRIVATE_KEY",
    "FAILED_V18_CONTINUATION_KEY_ROOT",
    "FAILED_V18_CONTINUATION_PUBLIC_KEY",
    "FAILED_V18_CONTINUATION_PRIVATE_KEY",
    "FAILED_V19_CONTINUATION_KEY_ROOT",
    "FAILED_V19_CONTINUATION_PUBLIC_KEY",
    "FAILED_V19_CONTINUATION_PRIVATE_KEY",
    "REJECTED_V20_CONTINUATION_KEY_ROOT",
    "REJECTED_V20_CONTINUATION_PUBLIC_KEY",
    "REJECTED_V20_CONTINUATION_PRIVATE_KEY",
    "FAILED_V21_CONTINUATION_KEY_ROOT",
    "FAILED_V21_CONTINUATION_PUBLIC_KEY",
    "FAILED_V21_CONTINUATION_PRIVATE_KEY",
    "REJECTED_V22_CONTINUATION_KEY_ROOT",
    "REJECTED_V22_CONTINUATION_PUBLIC_KEY",
    "REJECTED_V22_CONTINUATION_PRIVATE_KEY",
    "REJECTED_V23_CONTINUATION_KEY_ROOT",
    "REJECTED_V23_CONTINUATION_PUBLIC_KEY",
    "REJECTED_V23_CONTINUATION_PRIVATE_KEY",
    "REJECTED_V24_CONTINUATION_KEY_ROOT",
    "REJECTED_V24_CONTINUATION_PUBLIC_KEY",
    "REJECTED_V24_CONTINUATION_PRIVATE_KEY",
    "FAILED_V25_CONTINUATION_KEY_ROOT",
    "FAILED_V25_CONTINUATION_PUBLIC_KEY",
    "FAILED_V25_CONTINUATION_PRIVATE_KEY",
    "FAILED_V26_CONTINUATION_KEY_ROOT",
    "FAILED_V26_CONTINUATION_PUBLIC_KEY",
    "FAILED_V26_CONTINUATION_PRIVATE_KEY",
    "REJECTED_V27_CONTINUATION_KEY_ROOT",
    "REJECTED_V27_CONTINUATION_PUBLIC_KEY",
    "REJECTED_V27_CONTINUATION_PRIVATE_KEY",
    "FAILED_V28_CONTINUATION_KEY_ROOT",
    "FAILED_V28_CONTINUATION_PUBLIC_KEY",
    "FAILED_V28_CONTINUATION_PRIVATE_KEY",
    "FAILED_V29_CONTINUATION_KEY_ROOT",
    "FAILED_V29_CONTINUATION_PUBLIC_KEY",
    "FAILED_V29_CONTINUATION_PRIVATE_KEY",
    "FAILED_V28_RESULT_KEY_ARCHIVE",
    "FAILED_V28_REPORT_KEY_ARCHIVE",
    "FAILED_V28_RESULT_PUBLIC_KEY",
    "FAILED_V28_RESULT_KEY_ARCHIVE_RECORD",
    "FAILED_V28_REPORT_PUBLIC_KEY",
    "FAILED_V28_REPORT_PRIVATE_KEY",
    "FAILED_V28_REPORT_KEY_ARCHIVE_RECORD",
    "RECOVERY_CONTINUATION_KEY_ROOT",
    "RECOVERY_CONTINUATION_PUBLIC_KEY",
    "RECOVERY_CONTINUATION_PRIVATE_KEY",
    "PRIVATE_KEY_PATHS",
    "TERMINAL_KEY_ROTATION_KEYS",
    "TERMINAL_KEY_STATE_KEYS",
    "SUPERVISOR_CAPABILITY_FD",
    "LINUX_MFD_CLOEXEC",
    "LINUX_MFD_ALLOW_SEALING",
    "LINUX_F_ADD_SEALS",
    "LINUX_F_GET_SEALS",
    "LINUX_F_SEAL_SEAL",
    "LINUX_F_SEAL_SHRINK",
    "LINUX_F_SEAL_GROW",
    "LINUX_F_SEAL_WRITE",
    "SUPERVISOR_CAPABILITY_SCHEMA",
    "SUPERVISOR_CAPABILITY_SEAL_NAMES",
    "SUPERVISOR_CAPABILITY_PROTOCOL",
    "CONTINUATION_THREAT_MODEL",
    "V7_AUTHORITY",
    "V7_APPROVAL",
    "V7_SIGNATURE",
    "V7_PUBLIC_KEY",
    "V7_PRIVATE_KEY",
    "BASE_DOWNSTREAM_OUTPUT_SLOTS",
    "CURRENT_DOWNSTREAM_OUTPUT_SLOTS",
    "DOWNSTREAM_OUTPUT_SLOTS",
    "AUTHORIZED_DOWNSTREAM_OUTPUT_SLOTS",
    "REQUIRED_CHILD_OUTPUT_ROOTS",
    "OPTIONAL_CHILD_OUTPUT_ROOTS",
    "REQUIRED_CHILD_OUTPUT_FILES",
    "EXPECTED_ENVIRONMENT",
    "EXPECTED_PYTHON_SHA256",
    "EXPECTED_QA_PYTHON_SHA256",
    "EXPECTED_BASH_SHA256",
    "EXPECTED_OPENSSL_SHA256",
    "EXPECTED_NODE_SHA256",
    "EXPECTED_QA_CHROMIUM_SHA256",
    "EXPECTED_PORTABLE_BUILDER_SHA256",
    "EXPECTED_PORTABLE_CHART_EXTRACTOR_SHA256",
    "EXPECTED_PORTABLE_VERIFIER_SHA256",
    "EXPECTED_PORTABLE_BROWSER_HELPERS_SHA256",
    "EXPECTED_PORTABLE_BROWSER_CLI_SHA256",
    "EXPECTED_PORTABLE_SERVER_BUNDLE_SHA256",
    "EXPECTED_PORTABLE_READER_ASSET_SHA256",
    "REVIEWED_RUNTIME_SOURCE_CONTRACTS",
    "HISTORICAL_SIGNED_RUNTIME_SOURCE_ROLES",
    "RECOVERY_ONLY_RUNTIME_SOURCE_ROLES",
    "MANDATORY_GATE_EXECUTABLE_CONTRACTS",
    "EXPECTED_AUTHORIZATION_PUBLIC_KEY_SHA256",
    "EXPECTED_COMPLETION_PUBLIC_KEY_SHA256",
    "EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_FAILED_V16_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_FAILED_V17_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_FAILED_V18_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_FAILED_V19_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_REJECTED_V20_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_FAILED_V21_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_REJECTED_V22_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_REJECTED_V23_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_REJECTED_V24_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_FAILED_V25_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_FAILED_V26_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_REJECTED_V27_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_FAILED_V28_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_FAILED_V29_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_RECOVERY_CONTINUATION_PUBLIC_KEY_SHA256",
    "EXPECTED_RELEASE_SOURCE_AUTHORITY_SHA256",
    "EXPECTED_SOURCE_RECEIPT_SHA256",
    "EXPECTED_SOURCE_RECEIPT_SIZE",
    "EXPECTED_RUNNER_SHA256",
    "EXPECTED_RUNNER_LINE_COUNT",
    "EXPECTED_RUNNER_PREFIX_SHA256",
    "EXPECTED_RUNNER_PREFIX_SIZE",
    "EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SHA256",
    "EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SIZE",
    "EXPECTED_FAILED_V28_RESULT_PUBLIC_KEY_SHA256",
    "EXPECTED_FAILED_V28_REPORT_PUBLIC_KEY_SHA256",
    "EXPECTED_FAILED_V28_RESULT_ARCHIVE_RECORD_SHA256",
    "EXPECTED_FAILED_V28_REPORT_ARCHIVE_RECORD_SHA256",
    "RUNNER_ARCHIVE_REBASE",
    "RUNNER_CURRENT_RELEASE_KEY_REBASE",
    "FAILED_V28_LIVE_KEY_PATHS",
    "EXPECTED_RUNNER_LEAF_SHA256",
    "EXPECTED_V7_AUTHORITY_ID",
    "EXPECTED_V7_SOURCE_SET_SHA256",
    "AUTHORIZATION_ID",
    "AUTHORIZATION_SCHEMA",
    "COMPLETION_SCHEMA",
    "VERIFICATION_SCHEMA",
    "REPLAY_SCHEMA",
    "CONTINUATION_SCHEMA",
    "PREPARE_EXIT_ATTESTATION_SCHEMA",
    "PROMOTE_EXIT_ATTESTATION_SCHEMA",
    "SCHEMA_VERSION",
    "RECOVERY_SCHEMA_VERSION",
    "EXECUTABLE_ALIAS_RECORD_KEYS",
    "AUTHORIZATION_KEYS",
    "COMPLETION_KEYS",
    "PRODUCER_EXIT_ATTESTATION_KEYS",
    "PRODUCER_SESSION_KEYS",
    "PRODUCER_CHILD_WAIT_KEYS",
    "PREPARE_EXIT_CHECKS",
    "PROMOTE_EXIT_CHECKS",
    "PREPARE_EXIT_PASS_SEMANTICS",
    "PROMOTE_EXIT_PASS_SEMANTICS",
    "CONTINUATION_RECEIPT_KEYS",
    "CONTINUATION_EXIT_ATTESTATION_KEYS",
    "CONTINUATION_SUCCESS_WITNESS_KEYS",
    "SUPERVISOR_CAPABILITY_RECORD_KEYS",
    "SUPERVISOR_CAPABILITY_TOKEN_KEYS",
    "PRIVATE_KEY_PRE_SIGNATURE_KEYS",
    "VERIFICATION_KEYS",
    "REPLAY_KEYS",
    "AUTHORIZATION_CHECKS",
    "COMPLETION_CHECKS",
    "CONTINUATION_POLICY",
    "AUTHORIZED_CONTINUATION_POLICY",
    "TERMINAL_AUTHORITY_SEMANTICS",
    "TERMINAL_GOVERNANCE",
    "TERMINAL_CHECKS",
    "TERMINAL_PASS_SEMANTICS",
    "EXIT_ATTESTATION_CHECKS",
    "EXIT_ATTESTATION_PASS_SEMANTICS",
    "SUCCESS_WITNESS_CHECKS",
    "SUCCESS_WITNESS_PASS_SEMANTICS",
    "PREWITNESS_VERIFICATION_PASS_SEMANTICS",
    "REPLAY_PASS_SEMANTICS",
    "VERIFICATION_CHECKS",
    "VERIFICATION_GOVERNANCE",
    "VERIFICATION_PASS_SEMANTICS",
    "PROMOTION_VERIFIER_KEYS",
    "PROMOTION_TRUST_CHAIN_KEYS",
    "REPLAY_RUNTIME_KEYS",
    "REPLAY_SUCCESS_WITNESS_KEYS",
    "REPLAY_SUCCESS_WITNESS_CHECKS",
    "REPLAY_SUCCESS_WITNESS_PASS_SEMANTICS",
    "REPLAY_COMPLETION_RUNNER_KEYS",
    "REPLAY_RUNNER_ONLY_KEYS",
    "REPLAY_OUTPUT_SLOT_KEYS",
    "REPLAY_CLAIM_KEYS",
    "INOTIFY_MUTATION_MASK",
    "PRIMARY_PYTHON_WRAPPER",
    "QA_PYTHON_WRAPPER",
    "RUNTIME_FD_BRIDGE_TEMPLATE",
    "BRIDGE_BYTES",
    "EPILOGUE_BYTES",
    "GATE_INPUT_PATHS",
)

_BOOTSTRAP_LEASES: list[tuple[Path, int, os.stat_result]] = []
_EXECUTABLE_ALIAS_LEASES: dict[tuple[int, Path], "BoundExecutableAlias"] = {}
_DECLARED_RELATION_REGISTRY: weakref.WeakKeyDictionary[
    Any, dict[str, tuple[str, dict[str, Any]]]
] = weakref.WeakKeyDictionary()
_STREAM_LEASES: list["StreamBoundArtifactSet"] = []
_MUTATION_SENTINELS: list["PathMutationSentinel"] = []
_COMMIT_FDS: list[int] = []
_VERIFIER_LEASE_FDS: list[int] = []


class ContinuationError(RuntimeError):
    """Raised when the signed continuation contract cannot be proven exactly."""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def is_canonical_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def strict_json_equal(observed: Any, expected: Any) -> bool:
    """Compare decoded JSON without Python's bool/int equivalence."""

    if type(observed) is not type(expected):
        return False
    if isinstance(expected, dict):
        return set(observed) == set(expected) and all(
            strict_json_equal(observed[key], expected[key]) for key in expected
        )
    if isinstance(expected, list):
        return len(observed) == len(expected) and all(
            strict_json_equal(left, right)
            for left, right in zip(observed, expected, strict=True)
        )
    return observed == expected


def strict_json(data: bytes, label: str) -> dict[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        output: dict[str, Any] = {}
        for key, value in pairs:
            if key in output:
                raise ContinuationError(f"Duplicate JSON key in {label}: {key}")
            output[key] = value
        return output

    def reject_nonfinite(value: str) -> None:
        raise ContinuationError(f"Non-finite JSON value in {label}: {value}")

    try:
        payload = json.loads(
            data,
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_nonfinite,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContinuationError(f"Invalid JSON in {label}") from error
    if not isinstance(payload, dict):
        raise ContinuationError(f"{label} is not a JSON object")
    return payload


def encode_json(payload: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True, allow_nan=False)
        + "\n"
    ).encode("utf-8")


def require_exact_keys(value: Any, expected: frozenset[str], label: str) -> None:
    if not isinstance(value, Mapping) or set(value) != set(expected):
        observed = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise ContinuationError(
            f"{label} exact schema drift: observed={observed} expected={sorted(expected)}"
        )


def require_utc_timestamp(value: Any, label: str) -> datetime:
    if not isinstance(value, str):
        raise ContinuationError(f"{label} is not a timestamp string")
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise ContinuationError(f"{label} is not ISO-8601") from error
    if (
        parsed.tzinfo is None
        or parsed.utcoffset() != timezone.utc.utcoffset(parsed)
        or parsed.microsecond != 0
        or value != parsed.isoformat(timespec="seconds")
    ):
        raise ContinuationError(f"{label} is not canonical UTC")
    return parsed


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def same_stat(left: os.stat_result, right: os.stat_result) -> bool:
    return all(
        getattr(left, field) == getattr(right, field)
        for field in (
            "st_dev",
            "st_ino",
            "st_mode",
            "st_nlink",
            "st_size",
            "st_mtime_ns",
            "st_ctime_ns",
        )
    )


def stable_code_payload(code: types.CodeType) -> dict[str, Any]:
    constants: list[Any] = []
    for value in code.co_consts:
        if isinstance(value, types.CodeType):
            constants.append({"code": stable_code_payload(value)})
        elif isinstance(value, frozenset):
            constants.append({"frozenset": sorted(repr(item) for item in value)})
        elif isinstance(value, bytes):
            constants.append({"bytes_hex": value.hex()})
        else:
            constants.append(repr(value))
    return {
        "argcount": code.co_argcount,
        "posonlyargcount": code.co_posonlyargcount,
        "kwonlyargcount": code.co_kwonlyargcount,
        "nlocals": code.co_nlocals,
        "flags": code.co_flags,
        "firstlineno": code.co_firstlineno,
        "code_hex": code.co_code.hex(),
        "linetable_hex": code.co_linetable.hex(),
        "exceptiontable_hex": code.co_exceptiontable.hex(),
        "names": list(code.co_names),
        "varnames": list(code.co_varnames),
        "freevars": list(code.co_freevars),
        "cellvars": list(code.co_cellvars),
        "constants": constants,
    }


def stable_code_digest(code: types.CodeType) -> str:
    return sha256_bytes(
        json.dumps(stable_code_payload(code), sort_keys=True).encode("utf-8")
    )


def code_inventory(namespace: Mapping[str, Any], filename: str) -> dict[str, str]:
    inventory: dict[str, str] = {}
    for name, value in namespace.items():
        if isinstance(value, types.FunctionType) and value.__code__.co_filename == filename:
            inventory[name] = stable_code_digest(value.__code__)
        elif isinstance(value, type):
            for method_name, method in vars(value).items():
                function = (
                    method.__func__
                    if isinstance(method, (staticmethod, classmethod))
                    else method
                )
                if (
                    isinstance(function, types.FunctionType)
                    and function.__code__.co_filename == filename
                ):
                    inventory[f"{name}.{method_name}"] = stable_code_digest(
                        function.__code__
                    )
    return dict(sorted(inventory.items()))


def executing_module_code() -> types.CodeType:
    frame = sys._getframe()
    while frame is not None:
        if frame.f_code.co_name == "<module>" and frame.f_globals is globals():
            return frame.f_code
        frame = frame.f_back
    raise ContinuationError("Executing continuation module frame is unavailable")


def frozen_value(value: Any) -> Any:
    if isinstance(value, Path):
        return {"path": str(value)}
    if isinstance(value, Mapping):
        return {
            str(key): frozen_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [frozen_value(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted(
            (frozen_value(item) for item in value),
            key=lambda item: json.dumps(item, sort_keys=True),
        )
    if isinstance(value, bytes):
        return {"bytes_hex": value.hex()}
    if value is None or type(value) in (bool, int, float, str):
        return value
    raise ContinuationError(f"Unsupported critical-global type: {type(value).__name__}")


def verify_executing_source(self_data: bytes) -> dict[str, Any]:
    module_code = executing_module_code()
    filename = module_code.co_filename
    shadow: dict[str, Any] = {
        "__name__": "m2v5_downstream_continuation_bound_shadow",
        "__file__": filename,
        "__package__": "",
    }
    reviewed_code = compile(self_data, filename, "exec")
    exec(reviewed_code, shadow)
    live_inventory = code_inventory(globals(), filename)
    shadow_inventory = code_inventory(shadow, filename)
    live_critical = {
        name: frozen_value(globals()[name]) for name in CRITICAL_GLOBAL_NAMES
    }
    shadow_critical = {
        name: frozen_value(shadow[name]) for name in CRITICAL_GLOBAL_NAMES
    }
    checks = {
        "recursive_module_code_equal": (
            stable_code_digest(module_code) == stable_code_digest(reviewed_code)
        ),
        "function_and_class_method_inventory_nonempty": bool(live_inventory),
        "function_and_class_method_inventory_equal": live_inventory == shadow_inventory,
        "critical_globals_equal": live_critical == shadow_critical,
    }
    if not all(checks.values()):
        raise ContinuationError(
            "Executing continuation source differs from bound bytes: "
            + json.dumps(checks, sort_keys=True)
        )
    return {
        "source_sha256": sha256_bytes(self_data),
        "module_code_sha256": stable_code_digest(module_code),
        "function_and_class_method_count": len(live_inventory),
        "function_and_class_method_inventory_sha256": sha256_bytes(
            json.dumps(live_inventory, sort_keys=True).encode("utf-8")
        ),
        "critical_globals_sha256": sha256_bytes(
            json.dumps(live_critical, sort_keys=True).encode("utf-8")
        ),
        "pass": True,
    }


def recompute_execution_binding(data: bytes, filename: str) -> dict[str, Any]:
    namespace: dict[str, Any] = {
        "__name__": "promotion_execution_binding_independent_shadow",
        "__file__": filename,
        "__package__": "",
    }
    module_code = compile(data, filename, "exec")
    exec(module_code, namespace)
    inventory = code_inventory(namespace, filename)
    return {
        "source_sha256": sha256_bytes(data),
        "function_and_method_count": len(inventory),
        "inventory_sha256": sha256_bytes(
            json.dumps(inventory, sort_keys=True).encode("utf-8")
        ),
        "module_code_sha256": stable_code_digest(module_code),
        "critical_constants_equal": True,
        "canonical_commands_equal": True,
        "pass": True,
    }


def observed_command() -> list[str]:
    data = Path("/proc/self/cmdline").read_bytes()
    if not data.endswith(b"\0"):
        raise ContinuationError("Process command line is unavailable")
    return [os.fsdecode(token) for token in data[:-1].split(b"\0")]


def canonical_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(CONTINUATION_RUNNER),
        "--supervise-and-sign",
    ]


def authorized_canonical_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(AUTHORIZED_CONTINUATION_RUNNER),
        "--supervise-and-sign",
    ]


def supervised_child_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(CONTINUATION_RUNNER),
        "--supervised-child",
    ]


def authorized_supervised_child_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(AUTHORIZED_CONTINUATION_RUNNER),
        "--supervised-child",
    ]


def signed_terminal_verification_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(CONTINUATION_RUNNER),
        "--verify-signed-terminal",
    ]


def signed_terminal_prewitness_verification_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(CONTINUATION_RUNNER),
        "--verify-signed-terminal-prewitness",
    ]


def authorized_signed_terminal_verification_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(AUTHORIZED_CONTINUATION_RUNNER),
        "--verify-signed-terminal",
    ]


def require_clean_runtime(*, mode: str) -> None:
    commands = {
        "supervisor": canonical_command(),
        "child": supervised_child_command(),
        "signed_verifier": signed_terminal_verification_command(),
        "signed_prewitness_verifier": (
            signed_terminal_prewitness_verification_command()
        ),
    }
    argv = {
        "supervisor": [str(CONTINUATION_RUNNER), "--supervise-and-sign"],
        "child": [str(CONTINUATION_RUNNER), "--supervised-child"],
        "signed_verifier": [str(CONTINUATION_RUNNER), "--verify-signed-terminal"],
        "signed_prewitness_verifier": [
            str(CONTINUATION_RUNNER),
            "--verify-signed-terminal-prewitness",
        ],
    }
    if mode not in commands:
        raise ContinuationError(f"Unknown continuation runtime mode: {mode}")
    expected_command = commands[mode]
    expected_argv = argv[mode]
    observed = observed_command()
    actual_argv = list(sys.argv)
    fd_bound_source = False
    if mode in {"signed_verifier", "signed_prewitness_verifier"} and len(
        actual_argv
    ) == 2:
        source_token = actual_argv[0]
        prefix = "/proc/self/fd/"
        if source_token.startswith(prefix):
            try:
                source_fd = int(source_token.removeprefix(prefix))
            except ValueError as error:
                raise ContinuationError(
                    "Signed verifier source FD token is not canonical"
                ) from error
            source_stat = os.fstat(source_fd)
            canonical_stat = os.stat(CONTINUATION_RUNNER, follow_symlinks=False)
            fd_bound_source = (
                source_fd > 2
                and source_token == f"{prefix}{source_fd}"
                and same_stat(source_stat, canonical_stat)
                and actual_argv[1] == expected_argv[1]
                and observed
                == [
                    str(PYTHON),
                    "-I",
                    "-B",
                    source_token,
                    expected_argv[1],
                ]
            )
    if not fd_bound_source and (
        observed != expected_command or actual_argv != expected_argv
    ):
        raise ContinuationError("Continuation was not invoked by its canonical direct command")
    if dict(os.environ) != EXPECTED_ENVIRONMENT:
        raise ContinuationError("Continuation environment is not the exact clean environment")
    if Path(sys.executable) != PYTHON:
        raise ContinuationError("Python launcher token drift")
    if not (
        sys.flags.isolated == 1
        and sys.flags.ignore_environment == 1
        and sys.flags.no_user_site == 1
        and sys.flags.dont_write_bytecode == 1
    ):
        raise ContinuationError("Python isolation flags are incomplete")


def process_start_ticks(pid: int) -> int:
    try:
        data = Path(f"/proc/{pid}/stat").read_text(encoding="ascii")
        close = data.rfind(")")
        fields = data[close + 2 :].split()
        value = int(fields[19])
    except (OSError, UnicodeError, ValueError, IndexError) as error:
        raise ContinuationError(f"Cannot bind process start ticks for pid {pid}") from error
    if close <= 0 or value <= 0:
        raise ContinuationError(f"Invalid process start ticks for pid {pid}")
    return value


def process_command(pid: int) -> list[str]:
    try:
        data = Path(f"/proc/{pid}/cmdline").read_bytes()
    except OSError as error:
        raise ContinuationError(f"Cannot bind process command line for pid {pid}") from error
    if not data.endswith(b"\0"):
        raise ContinuationError(f"Process command line is incomplete for pid {pid}")
    return [os.fsdecode(token) for token in data[:-1].split(b"\0")]


def supervisor_capability_seal_mask() -> int:
    return (
        LINUX_F_SEAL_SEAL
        | LINUX_F_SEAL_SHRINK
        | LINUX_F_SEAL_GROW
        | LINUX_F_SEAL_WRITE
    )


def linux_memfd_create(name: str) -> int:
    libc = ctypes.CDLL(None, use_errno=True)
    function = getattr(libc, "memfd_create", None)
    if function is None:
        raise ContinuationError("libc memfd_create is unavailable")
    function.argtypes = [ctypes.c_char_p, ctypes.c_uint]
    function.restype = ctypes.c_int
    descriptor = int(
        function(
            name.encode("ascii"),
            LINUX_MFD_CLOEXEC | LINUX_MFD_ALLOW_SEALING,
        )
    )
    if descriptor < 0:
        error_number = ctypes.get_errno()
        raise ContinuationError(
            f"libc memfd_create failed: errno={error_number}"
        )
    return descriptor


def supervisor_capability_record(
    token: Mapping[str, Any], token_data: bytes
) -> dict[str, Any]:
    nonce = token.get("nonce")
    try:
        nonce_bytes = bytes.fromhex(nonce) if isinstance(nonce, str) else b""
    except ValueError as error:
        raise ContinuationError("Supervisor capability nonce is not hexadecimal") from error
    if len(nonce_bytes) != 32 or nonce != nonce_bytes.hex():
        raise ContinuationError("Supervisor capability nonce is not canonical 256-bit data")
    record = {
        "schema_name": SUPERVISOR_CAPABILITY_SCHEMA,
        "schema_version": "1.0.0",
        "protocol": SUPERVISOR_CAPABILITY_PROTOCOL,
        "token_size_bytes": len(token_data),
        "token_sha256": sha256_bytes(token_data),
        "nonce_sha256": sha256_bytes(nonce_bytes),
        "supervisor_pid": token.get("supervisor_pid"),
        "supervisor_start_ticks": token.get("supervisor_start_ticks"),
        "supervisor_command": token.get("supervisor_command"),
        "supervised_child_command": token.get("supervised_child_command"),
        "supervisor_source": token.get("supervisor_source"),
        "python_runtime": token.get("python_runtime"),
        "pass": True,
    }
    validate_supervisor_capability_record(record)
    return record


def validate_supervisor_capability_record(
    value: Any, chain: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    require_exact_keys(
        value,
        SUPERVISOR_CAPABILITY_RECORD_KEYS,
        "supervisor capability record",
    )
    if (
        value.get("schema_name") != SUPERVISOR_CAPABILITY_SCHEMA
        or value.get("schema_version") != "1.0.0"
        or value.get("protocol") != SUPERVISOR_CAPABILITY_PROTOCOL
        or type(value.get("token_size_bytes")) is not int
        or value["token_size_bytes"] <= 0
        or not is_canonical_sha256(value.get("token_sha256"))
        or not is_canonical_sha256(value.get("nonce_sha256"))
        or type(value.get("supervisor_pid")) is not int
        or value["supervisor_pid"] <= 1
        or type(value.get("supervisor_start_ticks")) is not int
        or value["supervisor_start_ticks"] <= 0
        or not strict_json_equal(value.get("supervisor_command"), canonical_command())
        or not strict_json_equal(
            value.get("supervised_child_command"), supervised_child_command()
        )
        or not isinstance(value.get("supervisor_source"), Mapping)
        or not isinstance(value.get("python_runtime"), Mapping)
        or value.get("pass") is not True
    ):
        raise ContinuationError("Supervisor capability record contract drift")
    if chain is not None and (
        not strict_json_equal(value.get("supervisor_source"), chain["self_record"])
        or not strict_json_equal(value.get("python_runtime"), chain["python_record"])
    ):
        raise ContinuationError("Supervisor capability source/runtime authority drift")
    return dict(value)


def reserve_supervisor_capability_descriptor() -> dict[str, Any]:
    try:
        os.fstat(SUPERVISOR_CAPABILITY_FD)
    except OSError as error:
        if error.errno != errno.EBADF:
            raise
    else:
        raise ContinuationError(
            "Reserved supervisor capability descriptor is occupied before reservation"
        )
    flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0)
    source = os.open("/dev/null", flags)
    try:
        if source == SUPERVISOR_CAPABILITY_FD:
            os.set_inheritable(source, False)
            descriptor = source
            source = -1
        else:
            descriptor = os.dup2(
                source, SUPERVISOR_CAPABILITY_FD, inheritable=False
            )
    finally:
        if source >= 0:
            os.close(source)
    opened = os.fstat(descriptor)
    if descriptor != SUPERVISOR_CAPABILITY_FD or not stat.S_ISCHR(opened.st_mode):
        raise ContinuationError("Supervisor capability descriptor reservation drift")
    return {"descriptor": descriptor, "stat": opened}


def create_supervisor_capability(
    chain: Mapping[str, Any], reservation: Mapping[str, Any]
) -> tuple[int, dict[str, Any]]:
    if (
        not isinstance(reservation, Mapping)
        or set(reservation) != {"descriptor", "stat"}
        or type(reservation.get("descriptor")) is not int
        or reservation["descriptor"] != SUPERVISOR_CAPABILITY_FD
        or not isinstance(reservation.get("stat"), os.stat_result)
        or not same_stat(
            reservation["stat"], os.fstat(SUPERVISOR_CAPABILITY_FD)
        )
        or not stat.S_ISCHR(reservation["stat"].st_mode)
    ):
        raise ContinuationError("Supervisor capability reservation was not retained")

    pid = os.getpid()
    token = {
        "schema_name": SUPERVISOR_CAPABILITY_SCHEMA,
        "schema_version": "1.0.0",
        "nonce": os.getrandom(32).hex(),
        "supervisor_pid": pid,
        "supervisor_start_ticks": process_start_ticks(pid),
        "supervisor_command": canonical_command(),
        "supervised_child_command": supervised_child_command(),
        "supervisor_source": chain["self_record"],
        "python_runtime": chain["python_record"],
    }
    token_data = encode_json(token)
    descriptor = linux_memfd_create("intersubmod-supervisor-capability")
    retained = -1
    try:
        offset = 0
        while offset < len(token_data):
            written = os.write(descriptor, token_data[offset:])
            if written <= 0:
                raise ContinuationError("Short write to supervisor capability memfd")
            offset += written
        os.fchmod(descriptor, 0o400)
        os.fsync(descriptor)
        seal_mask = supervisor_capability_seal_mask()
        fcntl.fcntl(descriptor, LINUX_F_ADD_SEALS, seal_mask)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or stat.S_IMODE(opened.st_mode) != 0o400
            or opened.st_nlink != 0
            or fcntl.fcntl(descriptor, LINUX_F_GET_SEALS) != seal_mask
            or os.pread(descriptor, opened.st_size, 0) != token_data
        ):
            raise ContinuationError("Supervisor capability memfd sealing drift")
        if descriptor == SUPERVISOR_CAPABILITY_FD:
            os.set_inheritable(descriptor, True)
            retained = descriptor
            descriptor = -1
        else:
            retained = os.dup2(
                descriptor, SUPERVISOR_CAPABILITY_FD, inheritable=True
            )
        if retained != SUPERVISOR_CAPABILITY_FD:
            raise ContinuationError("Supervisor capability descriptor assignment drift")
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    record = supervisor_capability_record(token, token_data)
    validate_supervisor_capability_record(record, chain)
    return retained, record


def bind_supervisor_capability() -> dict[str, Any]:
    descriptor = SUPERVISOR_CAPABILITY_FD
    try:
        opened = os.fstat(descriptor)
    except OSError as error:
        raise ContinuationError(
            "Supervised child lacks the inherited supervisor capability descriptor"
        ) from error
    seal_mask = supervisor_capability_seal_mask()
    data = os.pread(descriptor, opened.st_size, 0)
    if (
        not stat.S_ISREG(opened.st_mode)
        or stat.S_IMODE(opened.st_mode) != 0o400
        or opened.st_nlink != 0
        or fcntl.fcntl(descriptor, LINUX_F_GET_SEALS) != seal_mask
        or len(data) != opened.st_size
    ):
        raise ContinuationError("Inherited supervisor capability memfd drift")
    token = strict_json(data, "supervisor capability token")
    require_exact_keys(
        token, SUPERVISOR_CAPABILITY_TOKEN_KEYS, "supervisor capability token"
    )
    parent = os.getppid()
    parent_start = process_start_ticks(parent)
    self_executable = os.stat("/proc/self/exe")
    parent_executable = os.stat(f"/proc/{parent}/exe")
    if (
        token.get("schema_name") != SUPERVISOR_CAPABILITY_SCHEMA
        or token.get("schema_version") != "1.0.0"
        or type(token.get("supervisor_pid")) is not int
        or token.get("supervisor_pid") != parent
        or type(token.get("supervisor_start_ticks")) is not int
        or token.get("supervisor_start_ticks") != parent_start
        or not strict_json_equal(process_command(parent), canonical_command())
        or (self_executable.st_dev, self_executable.st_ino)
        != (parent_executable.st_dev, parent_executable.st_ino)
        or not strict_json_equal(token.get("supervisor_command"), canonical_command())
        or not strict_json_equal(
            token.get("supervised_child_command"), supervised_child_command()
        )
    ):
        raise ContinuationError("Supervisor capability live-parent binding drift")
    record = supervisor_capability_record(token, data)
    return validate_supervisor_capability_record(record)


def bootstrap_release_authority() -> tuple[Any, dict[str, Any]]:
    path = RELEASE_SOURCE_AUTHORITY
    if path.resolve(strict=True) != path or stat.S_ISLNK(os.lstat(path).st_mode):
        raise ContinuationError("Release authority source path is not canonical")
    flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    try:
        before = os.fstat(descriptor)
        data = b"".join(
            os.pread(descriptor, min(8 * 1024 * 1024, before.st_size - offset), offset)
            for offset in range(0, before.st_size, 8 * 1024 * 1024)
        )
        after = os.fstat(descriptor)
        live = os.lstat(path)
        if (
            not stat.S_ISREG(before.st_mode)
            or stat.S_IMODE(before.st_mode) != 0o444
            or len(data) != before.st_size
            or not same_stat(before, after)
            or not same_stat(after, live)
            or sha256_bytes(data) != EXPECTED_RELEASE_SOURCE_AUTHORITY_SHA256
        ):
            raise ContinuationError("Release authority bootstrap identity drift")
    except Exception:
        os.close(descriptor)
        raise
    module = types.ModuleType("bound_release_source_authority_for_m2v5_continuation")
    module.__file__ = str(path)
    module.__package__ = ""
    exec(compile(data, str(path), "exec"), module.__dict__)
    if not hasattr(module, "BoundArtifactReader") or not hasattr(
        module, "verify_ed25519_signature_fds"
    ):
        os.close(descriptor)
        raise ContinuationError("Bound release authority API is incomplete")
    _BOOTSTRAP_LEASES.append((path, descriptor, after))
    return module, {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": sha256_bytes(data),
        "mode": "0o444",
    }


def require_bootstrap_stable() -> None:
    for path, descriptor, opened in _BOOTSTRAP_LEASES:
        if not same_stat(opened, os.fstat(descriptor)) or not same_stat(
            opened, os.lstat(path)
        ):
            raise ContinuationError(f"Bootstrap path binding drift: {path}")
    for binding in _EXECUTABLE_ALIAS_LEASES.values():
        binding.require_stable()


def open_bound(
    reader: Any,
    path: Path,
    *,
    expected_mode: str | None = "0o444",
) -> tuple[int, bytes, dict[str, Any], os.stat_result]:
    descriptor, data, record = reader.open(path, include_mode=True)
    opened = os.fstat(descriptor)
    if expected_mode is not None and record.get("mode") != expected_mode:
        raise ContinuationError(
            f"Unexpected mode for {path}: {record.get('mode')} != {expected_mode}"
        )
    if opened.st_nlink != 1:
        raise ContinuationError(f"Trust artifact link count is not one: {path}")
    return descriptor, data, record, opened


class ExecutableAliasMutationWatch:
    """Fail on transient mutation of an executable alias or its target."""

    _EVENT_HEADER_SIZE = 16

    def __init__(self, alias_path: Path, target_path: Path) -> None:
        libc = ctypes.CDLL(None, use_errno=True)
        init = getattr(libc, "inotify_init1", None)
        add_watch = getattr(libc, "inotify_add_watch", None)
        if init is None or add_watch is None:
            raise ContinuationError("Linux inotify is required for alias mutation monitoring")
        init.argtypes = (ctypes.c_int,)
        init.restype = ctypes.c_int
        add_watch.argtypes = (ctypes.c_int, ctypes.c_char_p, ctypes.c_uint32)
        add_watch.restype = ctypes.c_int
        self._fd = init(os.O_NONBLOCK | os.O_CLOEXEC)
        if self._fd < 0:
            raise ContinuationError("inotify_init1 failed for executable alias")
        self._file_watches: set[int] = set()
        self._parent_watches: dict[int, set[str]] = {}
        try:
            for parent, name in (
                (alias_path.parent, alias_path.name),
                (target_path.parent, target_path.name),
            ):
                watch = add_watch(
                    self._fd, os.fsencode(parent), INOTIFY_MUTATION_MASK
                )
                if watch < 0:
                    raise ContinuationError(
                        f"inotify parent watch failed for executable alias: {parent}"
                    )
                self._parent_watches.setdefault(watch, set()).add(name)
            target_watch = add_watch(
                self._fd, os.fsencode(target_path), INOTIFY_MUTATION_MASK
            )
            if target_watch < 0:
                raise ContinuationError(
                    f"inotify target watch failed for executable alias: {target_path}"
                )
            self._file_watches.add(target_watch)
            self.assert_clean()
        except BaseException:
            self.close()
            raise

    def assert_clean(self) -> None:
        while True:
            try:
                data = os.read(self._fd, 1 << 16)
            except BlockingIOError:
                return
            if not data:
                return
            offset = 0
            while offset < len(data):
                if len(data) - offset < self._EVENT_HEADER_SIZE:
                    raise ContinuationError("Truncated executable alias inotify event")
                watch = int.from_bytes(data[offset : offset + 4], "little", signed=True)
                mask = int.from_bytes(data[offset + 4 : offset + 8], "little")
                name_size = int.from_bytes(data[offset + 12 : offset + 16], "little")
                end = offset + self._EVENT_HEADER_SIZE + name_size
                if end > len(data):
                    raise ContinuationError("Malformed executable alias inotify event")
                raw_name = data[offset + self._EVENT_HEADER_SIZE : end]
                name = os.fsdecode(raw_name.split(b"\0", 1)[0])
                terminal_mask = 0x00004000 | 0x00008000 | 0x00000400 | 0x00000800
                if (
                    mask & terminal_mask
                    or watch in self._file_watches
                    or name in self._parent_watches.get(watch, set())
                ):
                    raise ContinuationError("Executable alias mutation event observed")
                offset = end

    def close(self) -> None:
        if self._fd >= 0:
            os.close(self._fd)
            self._fd = -1


class BoundExecutableAlias:
    """Retain a symlink inode and its reviewed regular-file target."""

    def __init__(
        self,
        *,
        alias_path: Path,
        alias_fd: int,
        alias_stat: os.stat_result,
        link_text: str,
        resolved_path: Path,
        target_fd: int,
        target_stat: os.stat_result,
        mutation_watch: ExecutableAliasMutationWatch,
        record: Mapping[str, Any],
    ) -> None:
        self.alias_path = alias_path
        self.alias_fd = alias_fd
        self.alias_stat = alias_stat
        self.link_text = link_text
        self.resolved_path = resolved_path
        self.target_fd = target_fd
        self.target_stat = target_stat
        self.mutation_watch = mutation_watch
        self.record = dict(record)

    @classmethod
    def open(
        cls, reader: Any, record: Mapping[str, Any]
    ) -> "BoundExecutableAlias":
        require_exact_keys(record, EXECUTABLE_ALIAS_RECORD_KEYS, "executable alias record")
        if (
            not isinstance(record.get("path"), str)
            or not isinstance(record.get("resolved_path"), str)
            or not isinstance(record.get("symlink_target"), str)
            or type(record.get("size_bytes")) is not int
            or record["size_bytes"] < 0
            or not is_canonical_sha256(record.get("sha256"))
            or not isinstance(record.get("mode"), str)
        ):
            raise ContinuationError("Executable alias record scalar-type drift")
        alias_path = Path(record["path"])
        resolved_path = Path(record["resolved_path"])
        link_text = record["symlink_target"]
        if (
            not alias_path.is_absolute()
            or not resolved_path.is_absolute()
            or alias_path.parent.resolve(strict=True) != alias_path.parent
            or resolved_path.resolve(strict=True) != resolved_path
            or not link_text
        ):
            raise ContinuationError("Executable alias path contract drift")

        mutation_watch = ExecutableAliasMutationWatch(alias_path, resolved_path)
        alias_fd = -1
        try:
            before = os.lstat(alias_path)
            if not stat.S_ISLNK(before.st_mode) or before.st_nlink != 1:
                raise ContinuationError(f"Executable alias identity drift: {alias_path}")
            if not hasattr(os, "O_PATH"):
                raise ContinuationError("O_PATH is required for executable alias binding")
            alias_fd = os.open(
                alias_path,
                os.O_PATH | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0),
            )
            opened_alias = os.fstat(alias_fd)
            descriptor_link_text = os.readlink("", dir_fd=alias_fd)
            descriptor_resolved = (
                Path(descriptor_link_text)
                if Path(descriptor_link_text).is_absolute()
                else alias_path.parent / descriptor_link_text
            ).resolve(strict=True)
            if not stat.S_ISLNK(opened_alias.st_mode) or not same_stat(before, opened_alias):
                raise ContinuationError(f"Executable alias changed while opening: {alias_path}")
            if descriptor_link_text != link_text or descriptor_resolved != resolved_path:
                raise ContinuationError(f"Executable alias descriptor drift: {alias_path}")
            target_fd, _, target_record, target_stat = open_bound(
                reader,
                resolved_path,
                expected_mode=record["mode"],
            )
            expected_target = {
                "path": str(resolved_path),
                "size_bytes": record["size_bytes"],
                "sha256": record["sha256"],
                "mode": record["mode"],
            }
            if not strict_json_equal(target_record, expected_target):
                raise ContinuationError(f"Executable alias target bytes drift: {alias_path}")
            after_alias = os.lstat(alias_path)
            through_alias = os.stat(alias_path)
            if (
                not same_stat(before, after_alias)
                or os.readlink(alias_path) != link_text
                or os.readlink("", dir_fd=alias_fd) != link_text
                or descriptor_resolved != resolved_path
                or not same_stat(through_alias, target_stat)
            ):
                raise ContinuationError(f"Executable alias target binding drift: {alias_path}")
            mutation_watch.assert_clean()
            return cls(
                alias_path=alias_path,
                alias_fd=alias_fd,
                alias_stat=opened_alias,
                link_text=link_text,
                resolved_path=resolved_path,
                target_fd=target_fd,
                target_stat=target_stat,
                mutation_watch=mutation_watch,
                record=record,
            )
        except BaseException:
            if alias_fd >= 0:
                os.close(alias_fd)
            mutation_watch.close()
            raise

    def require_stable(self) -> None:
        self.mutation_watch.assert_clean()
        try:
            descriptor_alias = os.fstat(self.alias_fd)
            live_alias = os.lstat(self.alias_path)
            through_alias = os.stat(self.alias_path)
            target_descriptor = os.fstat(self.target_fd)
            live_target = os.lstat(self.resolved_path)
        except OSError as error:
            raise ContinuationError(
                f"Executable alias disappeared after binding: {self.alias_path}"
            ) from error
        if (
            not same_stat(self.alias_stat, descriptor_alias)
            or not same_stat(self.alias_stat, live_alias)
            or not stat.S_ISLNK(live_alias.st_mode)
            or os.readlink(self.alias_path) != self.link_text
            or os.readlink("", dir_fd=self.alias_fd) != self.link_text
            or not same_stat(self.target_stat, target_descriptor)
            or not same_stat(self.target_stat, live_target)
            or not same_stat(self.target_stat, through_alias)
        ):
            raise ContinuationError(f"Executable alias lease drift: {self.alias_path}")
        self.mutation_watch.assert_clean()


def bind_executable_alias_relation(
    node: Mapping[str, Any], reader: Any
) -> dict[str, Any]:
    alias_path = Path(str(node.get("path", "")))
    cache_key = (id(reader), alias_path)
    existing = _EXECUTABLE_ALIAS_LEASES.get(cache_key)
    if existing is not None:
        if not strict_json_equal(existing.record, dict(node)):
            raise ContinuationError(f"Conflicting executable alias relation: {alias_path}")
        existing.require_stable()
        return dict(existing.record)
    binding = BoundExecutableAlias.open(reader, node)
    _EXECUTABLE_ALIAS_LEASES[cache_key] = binding
    return dict(binding.record)


def open_private_metadata(
    reader: Any, path: Path, expected_mode: str
) -> tuple[int, dict[str, Any], os.stat_result]:
    descriptor, record = reader.open_metadata(path, include_mode=True)
    opened = os.fstat(descriptor)
    if not strict_json_equal(
        record, {"path": str(path.resolve(strict=True)), "mode": expected_mode}
    ):
        raise ContinuationError(f"Private-key mode contract drift: {path}")
    return descriptor, record, opened


def bind_running_python(reader: Any) -> tuple[int, dict[str, Any], os.stat_result]:
    if PYTHON.resolve(strict=True) != PRIMARY_PYTHON_RUNTIME:
        raise ContinuationError("Primary Python launcher target drift")
    python_fd, _, python_record, python_stat = open_bound(
        reader, PRIMARY_PYTHON_RUNTIME, expected_mode="0o775"
    )
    process_exe = os.stat("/proc/self/exe")
    if (
        python_record["sha256"] != EXPECTED_PYTHON_SHA256
        or (python_stat.st_dev, python_stat.st_ino)
        != (process_exe.st_dev, process_exe.st_ino)
    ):
        raise ContinuationError("Running Python interpreter identity drift")
    return python_fd, python_record, python_stat


def bind_reviewed_runtime_sources(
    reader: Any,
    *,
    primary_python_fd: int,
    primary_python_record: Mapping[str, Any],
) -> tuple[dict[str, dict[str, Any]], dict[str, int]]:
    if QA_PYTHON.resolve(strict=True) != QA_PYTHON_RUNTIME:
        raise ContinuationError("QA Python launcher target drift")
    records = {"primary_python_runtime": dict(primary_python_record)}
    descriptors = {"primary_python_runtime": primary_python_fd}
    primary_contract = REVIEWED_RUNTIME_SOURCE_CONTRACTS["primary_python_runtime"]
    if (
        primary_python_record.get("sha256") != primary_contract["sha256"]
        or primary_python_record.get("mode") != primary_contract["mode"]
        or primary_python_record.get("path") != str(primary_contract["path"])
    ):
        raise ContinuationError("Primary Python reviewed runtime contract drift")
    for role, contract in REVIEWED_RUNTIME_SOURCE_CONTRACTS.items():
        if role == "primary_python_runtime":
            continue
        descriptor, _, record, _ = open_bound(
            reader,
            contract["path"],
            expected_mode=str(contract["mode"]),
        )
        if (
            record.get("sha256") != contract["sha256"]
            or record.get("path") != str(contract["path"])
        ):
            raise ContinuationError(f"Reviewed runtime source identity drift: {role}")
        records[role] = record
        descriptors[role] = descriptor
    return records, descriptors


def historical_signed_runtime_projection(
    reviewed_runtime_sources: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, Any]]:
    current_roles = set(REVIEWED_RUNTIME_SOURCE_CONTRACTS)
    historical_roles = set(HISTORICAL_SIGNED_RUNTIME_SOURCE_ROLES)
    recovery_only_roles = set(RECOVERY_ONLY_RUNTIME_SOURCE_ROLES)
    if (
        historical_roles & recovery_only_roles
        or historical_roles | recovery_only_roles != current_roles
        or set(reviewed_runtime_sources) != current_roles
    ):
        raise ContinuationError("Historical/current runtime role partition drift")
    projection = {
        role: dict(reviewed_runtime_sources[role])
        for role in sorted(HISTORICAL_SIGNED_RUNTIME_SOURCE_ROLES)
    }
    if set(projection) != historical_roles:
        raise ContinuationError("Historical signed runtime projection drift")
    return projection


def with_stat(record: Mapping[str, Any], opened: os.stat_result) -> dict[str, Any]:
    return {
        **dict(record),
        "device": opened.st_dev,
        "inode": opened.st_ino,
        "link_count": opened.st_nlink,
        "mtime_ns": opened.st_mtime_ns,
        "ctime_ns": opened.st_ctime_ns,
    }


def basic_projection(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: record[key]
        for key in ("path", "size_bytes", "sha256")
        if key in record
    }


def validate_replay_source_identity_contract(
    trust_chain: Mapping[str, Any],
    expected_source_record: Mapping[str, Any],
    expected_canonical_record: Mapping[str, Any],
) -> None:
    for field, expected in (
        ("historical_source_receipt", expected_source_record),
        ("canonical_source_receipt", expected_canonical_record),
    ):
        observed = trust_chain.get(field)
        require_exact_keys(
            observed,
            STAT_ARTIFACT_KEYS,
            f"runner replay {field} full-stat record",
        )
        require_exact_keys(
            expected,
            STAT_ARTIFACT_KEYS,
            f"live {field} full-stat record",
        )
        if not strict_json_equal(observed, expected):
            raise ContinuationError(f"Runner replay {field} live identity drift")


def require_artifact_equal(
    observed: Any, expected: Mapping[str, Any], label: str
) -> None:
    if not isinstance(observed, Mapping) or not strict_json_equal(
        dict(observed), dict(expected)
    ):
        raise ContinuationError(f"{label} artifact identity drift")


class StreamBoundArtifactSet:
    """Stream-hash large gate inputs while retaining every descriptor."""

    def __init__(self) -> None:
        self._opened: dict[Path, tuple[int, os.stat_result, dict[str, Any], bytes | None]] = {}
        self._aliases: dict[
            Path, tuple[int, os.stat_result, str, Path, os.stat_result]
        ] = {}
        self._retained = False

    @property
    def descriptor_count(self) -> int:
        return len(self._opened) + len(self._aliases)

    def open(
        self,
        path: Path,
        *,
        capture_limit: int = 16 * 1024 * 1024,
        require_nonempty: bool = True,
        expected_mode: str | None = None,
        expected_link_count: int | None = None,
    ) -> tuple[int, dict[str, Any], bytes | None]:
        resolved = path.resolve(strict=True)
        if resolved != path or stat.S_ISLNK(os.lstat(path).st_mode):
            raise ContinuationError(f"Gate input path is not canonical: {path}")
        existing = self._opened.get(resolved)
        if existing is not None:
            return existing[0], dict(existing[2]), existing[3]
        descriptor = os.open(
            resolved,
            os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0),
        )
        try:
            before = os.fstat(descriptor)
            if not stat.S_ISREG(before.st_mode) or (require_nonempty and before.st_size <= 0):
                raise ContinuationError(f"Gate input is not a nonempty regular file: {resolved}")
            if expected_mode is not None and oct(stat.S_IMODE(before.st_mode)) != expected_mode:
                raise ContinuationError(f"Gate input mode drift: {resolved}")
            if expected_link_count is not None and before.st_nlink != expected_link_count:
                raise ContinuationError(f"Gate input link-count drift: {resolved}")
            digest = hashlib.sha256()
            captured: list[bytes] | None = [] if before.st_size <= capture_limit else None
            offset = 0
            while offset < before.st_size:
                block = os.pread(
                    descriptor,
                    min(8 * 1024 * 1024, before.st_size - offset),
                    offset,
                )
                if not block:
                    raise ContinuationError(f"Short read from gate input: {resolved}")
                digest.update(block)
                if captured is not None:
                    captured.append(block)
                offset += len(block)
            after = os.fstat(descriptor)
            live = os.lstat(resolved)
            if offset != before.st_size or not same_stat(before, after) or not same_stat(
                after, live
            ):
                raise ContinuationError(f"Gate input changed while hashing: {resolved}")
            record = {
                "path": str(resolved),
                "size_bytes": after.st_size,
                "sha256": digest.hexdigest(),
                "mode": oct(stat.S_IMODE(after.st_mode)),
            }
            data = b"".join(captured) if captured is not None else None
            self._opened[resolved] = (descriptor, after, record, data)
            return descriptor, dict(record), data
        except Exception:
            os.close(descriptor)
            raise

    def open_alias(
        self,
        alias: Path,
        target: Path,
        *,
        expected_link_text: str,
    ) -> dict[str, Any]:
        if not hasattr(os, "O_PATH"):
            raise ContinuationError("O_PATH is required for symlink gate binding")
        if target not in self._opened:
            raise ContinuationError(f"Alias target was not bound first: {target}")
        if (
            alias.parent.resolve(strict=True) != alias.parent
            or target.resolve(strict=True) != target
            or alias.resolve(strict=True) != target
            or Path(expected_link_text).is_absolute()
            or Path(expected_link_text).name != expected_link_text
            or expected_link_text != target.name
        ):
            raise ContinuationError(f"Gate alias target contract drift: {alias}")
        before = os.lstat(alias)
        descriptor = os.open(
            alias,
            os.O_PATH | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0),
        )
        try:
            opened = os.fstat(descriptor)
            link_text = os.readlink(alias)
            after = os.lstat(alias)
            target_opened = self._opened[target][1]
            target_live = os.lstat(target)
            if (
                not stat.S_ISLNK(before.st_mode)
                or before.st_nlink != 1
                or not same_stat(before, opened)
                or not same_stat(opened, after)
                or link_text != expected_link_text
                or alias.resolve(strict=True) != target
                or not same_stat(target_opened, target_live)
            ):
                raise ContinuationError(f"Gate alias identity drift: {alias}")
            self._aliases[alias] = (
                descriptor,
                opened,
                link_text,
                target,
                target_opened,
            )
            return {
                "path": str(alias),
                "link_text": link_text,
                "resolved_target": str(target),
                "mode": oct(stat.S_IMODE(opened.st_mode)),
                "link_count": opened.st_nlink,
            }
        except Exception:
            os.close(descriptor)
            raise

    def open_metadata(self, path: Path, expected_mode: str) -> tuple[int, dict[str, Any]]:
        if not hasattr(os, "O_PATH"):
            raise ContinuationError("O_PATH is required for metadata-only gate binding")
        resolved = path.resolve(strict=True)
        if resolved != path or stat.S_ISLNK(os.lstat(path).st_mode):
            raise ContinuationError(f"Metadata gate path is not canonical: {path}")
        descriptor = os.open(
            resolved, os.O_PATH | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0)
        )
        opened = os.fstat(descriptor)
        live = os.lstat(resolved)
        if (
            not stat.S_ISREG(opened.st_mode)
            or not same_stat(opened, live)
            or oct(stat.S_IMODE(opened.st_mode)) != expected_mode
        ):
            os.close(descriptor)
            raise ContinuationError(f"Metadata gate identity drift: {resolved}")
        record = {"path": str(resolved), "mode": expected_mode}
        self._opened[resolved] = (descriptor, opened, record, None)
        return descriptor, dict(record)

    def data_for(self, path: Path) -> bytes:
        resolved = path.resolve(strict=True)
        if resolved not in self._opened or self._opened[resolved][3] is None:
            raise ContinuationError(f"Gate input bytes were not captured: {resolved}")
        return self._opened[resolved][3]  # type: ignore[return-value]

    def record_for(self, path: Path) -> dict[str, Any]:
        resolved = path.resolve(strict=True)
        if resolved not in self._opened:
            raise ContinuationError(f"Gate input was not bound: {resolved}")
        return dict(self._opened[resolved][2])

    def alias_record_for(self, alias: Path) -> dict[str, Any]:
        if alias not in self._aliases:
            raise ContinuationError(f"Gate alias was not bound: {alias}")
        _, opened, link_text, target, _ = self._aliases[alias]
        return {
            "path": str(alias),
            "link_text": link_text,
            "resolved_target": str(target),
            "mode": oct(stat.S_IMODE(opened.st_mode)),
            "link_count": opened.st_nlink,
        }

    def fd_for(self, path: Path) -> int:
        resolved = path.resolve(strict=True)
        if resolved not in self._opened:
            raise ContinuationError(f"Gate input was not bound: {resolved}")
        return self._opened[resolved][0]

    def require_paths_still_bound(self) -> None:
        for path, (descriptor, opened, _, _) in self._opened.items():
            if not same_stat(opened, os.fstat(descriptor)) or not same_stat(
                opened, os.lstat(path)
            ):
                raise ContinuationError(f"Bound gate input drifted: {path}")
        for alias, (descriptor, opened, link_text, target, target_opened) in (
            self._aliases.items()
        ):
            if (
                not same_stat(opened, os.fstat(descriptor))
                or not same_stat(opened, os.lstat(alias))
                or os.readlink(alias) != link_text
                or alias.resolve(strict=True) != target
                or not same_stat(target_opened, os.lstat(target))
            ):
                raise ContinuationError(f"Bound gate alias drifted: {alias}")

    def retain_until_process_exit(self) -> None:
        if self._retained:
            raise ContinuationError("Stream artifact set was already retained")
        self.require_paths_still_bound()
        _STREAM_LEASES.append(self)
        self._retained = True

    def close(self) -> None:
        if self._retained:
            return
        while self._aliases:
            _, (descriptor, _, _, _, _) = self._aliases.popitem()
            os.close(descriptor)
        while self._opened:
            _, (descriptor, _, _, _) = self._opened.popitem()
            os.close(descriptor)


class PathMutationSentinel:
    """Detect even transient replacement or mutation of protected input paths."""

    _EVENT_HEADER_SIZE = 16

    @staticmethod
    def _allows_parent_only_mode000_fallback(
        *, error_number: int, file_mode: int
    ) -> bool:
        return error_number == errno.EACCES and file_mode == 0o000

    def __init__(
        self,
        paths: set[Path],
        *,
        all_event_directories: set[Path] | None = None,
    ) -> None:
        libc = ctypes.CDLL(None, use_errno=True)
        init = getattr(libc, "inotify_init1", None)
        add_watch = getattr(libc, "inotify_add_watch", None)
        if init is None or add_watch is None:
            raise ContinuationError("Linux inotify is required for mutation monitoring")
        init.argtypes = (ctypes.c_int,)
        init.restype = ctypes.c_int
        add_watch.argtypes = (ctypes.c_int, ctypes.c_char_p, ctypes.c_uint32)
        add_watch.restype = ctypes.c_int
        descriptor = init(os.O_NONBLOCK | os.O_CLOEXEC)
        if descriptor < 0:
            error_number = ctypes.get_errno()
            raise ContinuationError("inotify_init1 failed: " + os.strerror(error_number))
        self._fd = descriptor
        self._file_watches: dict[int, Path] = {}
        self._parent_watches: dict[int, tuple[Path, set[str]]] = {}
        self._all_event_parent_watches: set[int] = set()
        self._parent_only_mode000_paths: set[Path] = set()
        self._retained = False
        self._paths = tuple(sorted(paths))
        if any(
            not path.is_absolute()
            or Path(os.path.normpath(str(path))) != path
            for path in self._paths
        ):
            os.close(self._fd)
            raise ContinuationError("Mutation sentinel paths must be normalized absolute paths")
        try:
            parent_watch_descriptors: dict[Path, int] = {}
            for directory in sorted({path.parent for path in self._paths}):
                parent_wd = add_watch(
                    self._fd, os.fsencode(directory), INOTIFY_MUTATION_MASK
                )
                if parent_wd < 0:
                    error_number = ctypes.get_errno()
                    raise ContinuationError(
                        f"inotify parent watch failed for {directory}: "
                        + os.strerror(error_number)
                    )
                directory_stat = os.lstat(directory)
                if (
                    directory.resolve(strict=True) != directory
                    or stat.S_ISLNK(directory_stat.st_mode)
                    or not stat.S_ISDIR(directory_stat.st_mode)
                ):
                    raise ContinuationError(
                        f"Mutation sentinel parent is not canonical: {directory}"
                    )
                parent_watch_descriptors[directory] = parent_wd
                existing = self._parent_watches.get(parent_wd)
                if existing is not None and existing[0] != directory:
                    raise ContinuationError("Mutation sentinel parent watch alias drift")
                self._parent_watches.setdefault(parent_wd, (directory, set()))
            for path in self._paths:
                self._parent_watches[parent_watch_descriptors[path.parent]][1].add(
                    path.name
                )

            for directory in sorted(all_event_directories or set()):
                if (
                    not directory.is_absolute()
                    or Path(os.path.normpath(str(directory))) != directory
                ):
                    raise ContinuationError(
                        f"Mutation sentinel directory is not normalized: {directory}"
                    )
                directory_wd = add_watch(
                    self._fd,
                    os.fsencode(directory),
                    INOTIFY_MUTATION_MASK,
                )
                if directory_wd < 0:
                    error_number = ctypes.get_errno()
                    raise ContinuationError(
                        f"inotify directory watch failed for {directory}: "
                        + os.strerror(error_number)
                    )
                directory_stat = os.lstat(directory)
                if (
                    directory.resolve(strict=True) != directory
                    or stat.S_ISLNK(directory_stat.st_mode)
                    or not stat.S_ISDIR(directory_stat.st_mode)
                ):
                    raise ContinuationError(
                        f"Mutation sentinel directory is not canonical: {directory}"
                    )
                existing = self._parent_watches.get(directory_wd)
                if existing is not None and existing[0] != directory:
                    raise ContinuationError("Mutation sentinel directory watch alias drift")
                self._parent_watches.setdefault(directory_wd, (directory, set()))
                self._all_event_parent_watches.add(directory_wd)

            initial_stats: dict[Path, os.stat_result] = {}
            for path in self._paths:
                if path.resolve(strict=True) != path:
                    raise ContinuationError(
                        f"Mutation sentinel path is not canonical: {path}"
                    )
                leaf = os.lstat(path)
                if stat.S_ISLNK(leaf.st_mode) or not stat.S_ISREG(leaf.st_mode):
                    raise ContinuationError(
                        f"Mutation sentinel path is not a regular file: {path}"
                    )
                initial_stats[path] = leaf

            for path, leaf in initial_stats.items():
                file_wd = add_watch(
                    self._fd, os.fsencode(path), INOTIFY_MUTATION_MASK
                )
                if file_wd < 0:
                    error_number = ctypes.get_errno()
                    if self._allows_parent_only_mode000_fallback(
                        error_number=error_number,
                        file_mode=stat.S_IMODE(leaf.st_mode),
                    ):
                        self._parent_only_mode000_paths.add(path)
                    else:
                        raise ContinuationError(
                            f"inotify file watch failed for {path}: "
                            + os.strerror(error_number)
                        )
                else:
                    self._file_watches[file_wd] = path

            for path, initial in initial_stats.items():
                if path.resolve(strict=True) != path or not same_stat(
                    initial, os.lstat(path)
                ):
                    raise ContinuationError(
                        f"Mutation sentinel setup identity drifted: {path}"
                    )
            self._setup_identity_recheck_count = len(initial_stats)
            self.assert_clean()
        except Exception:
            os.close(self._fd)
            raise

    @property
    def protected_path_count(self) -> int:
        return len(self._paths)

    @property
    def protected_path_set_sha256(self) -> str:
        encoded = ("\n".join(str(path) for path in self._paths) + "\n").encode(
            "utf-8"
        )
        return sha256_bytes(encoded)

    @property
    def direct_file_watch_count(self) -> int:
        return len(self._paths) - len(self._parent_only_mode000_paths)

    @property
    def parent_only_mode000_path_count(self) -> int:
        return len(self._parent_only_mode000_paths)

    @property
    def parent_only_mode000_path_set_sha256(self) -> str:
        encoded = (
            "\n".join(str(path) for path in sorted(self._parent_only_mode000_paths))
            + "\n"
        ).encode("utf-8")
        return sha256_bytes(encoded)

    @property
    def setup_identity_recheck_count(self) -> int:
        return self._setup_identity_recheck_count

    @property
    def setup_identity_recheck_path_set_sha256(self) -> str:
        return self.protected_path_set_sha256

    def assert_clean(self) -> None:
        violations: list[str] = []
        while True:
            try:
                block = os.read(self._fd, 1024 * 1024)
            except BlockingIOError:
                break
            if not block:
                break
            offset = 0
            while offset < len(block):
                if len(block) - offset < self._EVENT_HEADER_SIZE:
                    raise ContinuationError("Truncated inotify event header")
                wd = int.from_bytes(
                    block[offset : offset + 4], sys.byteorder, signed=True
                )
                mask = int.from_bytes(
                    block[offset + 4 : offset + 8], sys.byteorder
                )
                name_length = int.from_bytes(
                    block[offset + 12 : offset + 16], sys.byteorder
                )
                event_end = offset + self._EVENT_HEADER_SIZE + name_length
                if event_end > len(block):
                    raise ContinuationError("Truncated inotify event payload")
                raw_name = block[offset + self._EVENT_HEADER_SIZE : event_end]
                name = raw_name.split(b"\0", 1)[0].decode(
                    "utf-8", errors="surrogateescape"
                )
                if mask & IN_Q_OVERFLOW:
                    violations.append("inotify_queue_overflow")
                elif wd in self._file_watches:
                    violations.append(
                        f"file:{self._file_watches[wd]}:mask={mask:#x}"
                    )
                elif wd in self._parent_watches:
                    parent, protected_names = self._parent_watches[wd]
                    parent_self_event = not name and mask & (
                        IN_DELETE_SELF | IN_MOVE_SELF | IN_UNMOUNT | IN_IGNORED
                    )
                    if (
                        wd in self._all_event_parent_watches
                        or name in protected_names
                        or parent_self_event
                    ):
                        violations.append(f"parent:{parent}/{name}:mask={mask:#x}")
                offset = event_end
        if violations:
            raise ContinuationError(
                "Protected path mutation was observed: " + "; ".join(violations[:8])
            )

    def retain_until_process_exit(self) -> None:
        if self._retained:
            raise ContinuationError("Mutation sentinel was already retained")
        self.assert_clean()
        _MUTATION_SENTINELS.append(self)
        self._retained = True


def collect_protected_paths(*values: Any) -> set[Path]:
    """Collect existing artifact paths from signed records and gate projections."""

    protected: set[Path] = set()
    excluded = set(PRIVATE_KEY_PATHS)

    def visit(value: Any) -> None:
        if isinstance(value, Mapping):
            path_value = value.get("path")
            if isinstance(path_value, str) and os.path.lexists(path_value):
                path = Path(path_value).resolve(strict=True)
                if path not in excluded and path.is_file():
                    protected.add(path)
            for child in value.values():
                visit(child)
        elif isinstance(value, (list, tuple)):
            for child in value:
                visit(child)

    for value in values:
        visit(value)
    protected.update(
        path.resolve(strict=True)
        for path in (
            CONTINUATION_RUNNER,
            AUTHORIZED_CONTINUATION_RUNNER,
            PROMOTION_TOOL,
            CONTINUATION_VERIFIER,
            AUTHORIZED_CONTINUATION_VERIFIER,
            RUNNER_GATE_REPLAY,
            AUTHORIZED_RUNNER_GATE_REPLAY,
            SCHEMA_RECOVERY_VALIDATOR,
            COMPLETION_RUNNER,
            RELEASE_SOURCE_AUTHORITY,
            AUTHORIZATION,
            PREPARE_EXIT_ATTESTATION,
            AUTHORIZATION_SIGNATURE,
            COMPLETION,
            PROMOTE_EXIT_ATTESTATION,
            COMPLETION_SIGNATURE,
            PROMOTION_VERIFICATION_RECEIPT,
            RUNNER_GATE_RECEIPT,
            RUNNER_GATE_LOG,
            HISTORICAL_SOURCE_RECEIPT,
            CANONICAL_SOURCE_RECEIPT,
        )
    )
    return protected


class AuthorizedModeTransition:
    """Bind a one-time key before and after its only permitted chmod 0400->0000."""

    def __init__(self, authority_module: Any, path: Path) -> None:
        self.path = path
        self.reader = authority_module.BoundArtifactReader()
        self.descriptor, self.pre_record = self.reader.open_metadata(
            path, include_mode=True
        )
        self.pre_stat = os.fstat(self.descriptor)
        self.post_descriptor: int | None = None
        self.post_stat: os.stat_result | None = None
        if not strict_json_equal(
            self.pre_record,
            {"path": str(path.resolve(strict=True)), "mode": "0o400"},
        ):
            self.reader.close()
            raise ContinuationError(f"One-time key is not in pre-signing mode 0400: {path}")
        self.reader.retain_until_process_exit()

    def verify_retired(self, post_reader: Any) -> dict[str, Any]:
        descriptor, record = post_reader.open_metadata(self.path, include_mode=True)
        post = os.fstat(descriptor)
        live_pre_fd = os.fstat(self.descriptor)
        stable_fields = ("st_dev", "st_ino", "st_nlink", "st_size", "st_mtime_ns")
        if (
            not strict_json_equal(
                record,
                {"path": str(self.path.resolve(strict=True)), "mode": "0o0"},
            )
            or any(getattr(post, field) != getattr(self.pre_stat, field) for field in stable_fields)
            or any(getattr(live_pre_fd, field) != getattr(post, field) for field in stable_fields)
            or stat.S_IMODE(post.st_mode) != 0
            or stat.S_IMODE(live_pre_fd.st_mode) != 0
            or post.st_ctime_ns < self.pre_stat.st_ctime_ns
        ):
            raise ContinuationError(f"One-time key retirement transition drift: {self.path}")
        self.post_descriptor = descriptor
        self.post_stat = post
        return dict(record)

    def require_retired_stable(self) -> None:
        if self.post_descriptor is None or self.post_stat is None:
            raise ContinuationError(
                f"One-time key retirement was not established: {self.path}"
            )
        pre_live = os.fstat(self.descriptor)
        post_live = os.fstat(self.post_descriptor)
        path_live = os.lstat(self.path)
        if (
            not same_stat(self.post_stat, pre_live)
            or not same_stat(self.post_stat, post_live)
            or not same_stat(self.post_stat, path_live)
            or stat.S_IMODE(path_live.st_mode) != 0
            or path_live.st_nlink != 1
            or self.path.resolve(strict=True) != self.path
        ):
            raise ContinuationError(
                f"One-time key retirement changed after verification: {self.path}"
            )


def bind_declared_relations(
    value: Any,
    reader: Any,
    *,
    historical_mode_records: set[tuple[tuple[str, ...], str, str]],
    relations: dict[str, dict[str, Any]] | None = None,
) -> dict[str, dict[str, Any]]:
    output = {} if relations is None else relations

    def validate_scalar_types(kind: str, node: Mapping[str, Any]) -> None:
        if not isinstance(node.get("path"), str) or not node["path"]:
            raise ContinuationError(f"{kind} relation path type drift")
        if "mode" in node:
            mode = node.get("mode")
            if (
                not isinstance(mode, str)
                or not mode.startswith("0o")
                or len(mode) <= 2
                or any(character not in "01234567" for character in mode[2:])
                or oct(int(mode, 8)) != mode
            ):
                raise ContinuationError(f"{kind} relation mode type drift")
        if "size_bytes" in node and (
            type(node.get("size_bytes")) is not int or node["size_bytes"] < 0
        ):
            raise ContinuationError(f"{kind} relation size type drift")
        if "sha256" in node and not is_canonical_sha256(node.get("sha256")):
            raise ContinuationError(f"{kind} relation SHA-256 type drift")
        if kind == "executable_alias" and (
            not isinstance(node.get("resolved_path"), str)
            or not node["resolved_path"]
            or not isinstance(node.get("symlink_target"), str)
            or not node["symlink_target"]
        ):
            raise ContinuationError("Executable alias relation scalar-type drift")
        if kind == "metadata_plus_size" and (
            node.get("mode") not in {"0o0", "0o400"}
            or not isinstance(node.get("state"), str)
            or not node["state"]
        ):
            raise ContinuationError(
                "Metadata-plus-size relation scalar-type drift"
            )
        if kind in {"stat", "legacy_stat"}:
            fields = ["device", "inode", "mtime_ns", "ctime_ns"]
            if kind == "stat":
                fields.append("link_count")
            for field in fields:
                if type(node.get(field)) is not int or node[field] < 0:
                    raise ContinuationError(
                        f"Declared stat relation integer type drift: {field}"
                    )

    def classify(node: Mapping[str, Any]) -> str | None:
        keys = frozenset(node)
        if keys == EXECUTABLE_ALIAS_RECORD_KEYS:
            kind = "executable_alias"
        elif keys == BASIC_ARTIFACT_KEYS:
            kind = "basic"
        elif keys == STAT_ARTIFACT_KEYS:
            kind = "stat"
        elif keys == LEGACY_STAT_ARTIFACT_KEYS:
            kind = "legacy_stat"
        elif keys == METADATA_ARTIFACT_KEYS:
            kind = "metadata"
        elif keys == METADATA_PLUS_SIZE_STATE_KEYS:
            kind = "metadata_plus_size"
        elif "path" in keys and keys & IDENTITY_MARKER_KEYS:
            raise ContinuationError("Ambiguous identity-like relation schema")
        else:
            return None
        validate_scalar_types(kind, node)
        return kind

    def register(kind: str, node: Mapping[str, Any]) -> None:
        path_text = str(node["path"])
        declared = dict(node)
        registry = _DECLARED_RELATION_REGISTRY.setdefault(reader, {})
        existing = registry.get(path_text)
        if existing is not None:
            existing_kind, existing_declared = existing
            metadata_kinds = {"metadata", "metadata_plus_size"}
            compatible_metadata_enrichment = (
                {kind, existing_kind} == metadata_kinds
                and all(
                    strict_json_equal(existing_declared.get(key), declared.get(key))
                    for key in METADATA_ARTIFACT_KEYS
                )
            )
            if compatible_metadata_enrichment:
                if kind == "metadata_plus_size":
                    registry[path_text] = (kind, declared)
                return
            content_kinds = {"basic", "stat", "legacy_stat"}
            compatible_content_enrichment = (
                kind in content_kinds
                and existing_kind in content_kinds
                and kind != existing_kind
                and all(
                    strict_json_equal(existing_declared.get(key), declared.get(key))
                    for key in BASIC_ARTIFACT_KEYS
                )
            )
            if compatible_content_enrichment:
                return
            if existing_kind != kind or not strict_json_equal(
                existing_declared, declared
            ):
                raise ContinuationError(
                    f"Conflicting declared relation schema or value: {path_text}"
                )
            return
        registry[path_text] = (kind, declared)

    def visit(node: Any, location: tuple[str, ...]) -> None:
        if isinstance(node, Mapping):
            kind = classify(node)
            historical_relation = (
                kind is not None
                and "mode" in node
                and (location, str(node["path"]), str(node["mode"]))
                in historical_mode_records
            )
            if kind is not None and not historical_relation:
                register(kind, node)
            if kind == "executable_alias":
                live = bind_executable_alias_relation(node, reader)
                output[str(node["path"])] = live
            elif kind in {"basic", "stat", "legacy_stat"}:
                path = Path(node["path"])
                declared_mode = node["mode"]
                if not historical_relation:
                    descriptor, _, live, opened = open_bound(
                        reader,
                        path,
                        expected_mode=declared_mode if declared_mode else None,
                    )
                    del descriptor
                    expected_basic = {
                        key: node[key]
                        for key in ("path", "size_bytes", "sha256", "mode")
                    }
                    observed_basic = {key: live[key] for key in expected_basic}
                    if not strict_json_equal(observed_basic, expected_basic):
                        raise ContinuationError(f"Declared artifact relation drift: {path}")
                    if kind in {"stat", "legacy_stat"}:
                        observed_stat = with_stat(live, opened)
                        expected_keys = (
                            STAT_ARTIFACT_KEYS
                            if kind == "stat"
                            else LEGACY_STAT_ARTIFACT_KEYS
                        )
                        expected_stat = {key: node[key] for key in expected_keys}
                        observed_stat = {
                            key: observed_stat[key] for key in expected_keys
                        }
                        if not strict_json_equal(observed_stat, expected_stat):
                            raise ContinuationError(f"Declared stat relation drift: {path}")
                    output[str(path)] = dict(live)
            elif kind == "metadata":
                path_text = node["path"]
                mode = node["mode"]
                if not historical_relation:
                    _, record, opened = open_private_metadata(
                        reader, Path(path_text), mode
                    )
                    registered_kind, _ = _DECLARED_RELATION_REGISTRY[reader][
                        path_text
                    ]
                    output[path_text] = (
                        {**record, "size_bytes": opened.st_size}
                        if registered_kind == "metadata_plus_size"
                        else record
                    )
            elif kind == "metadata_plus_size":
                path_text = node["path"]
                mode = node["mode"]
                if not historical_relation:
                    _, record, opened = open_private_metadata(
                        reader, Path(path_text), mode
                    )
                    if opened.st_size != node["size_bytes"]:
                        raise ContinuationError(
                            f"Declared metadata-plus-size relation drift: {path_text}"
                        )
                    output[path_text] = {
                        **record,
                        "size_bytes": opened.st_size,
                    }
            for key, child in node.items():
                visit(child, (*location, str(key)))
        elif isinstance(node, list):
            for index, child in enumerate(node):
                visit(child, (*location, str(index)))

    visit(value, ())
    return output


def promotion_command(mode: str) -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        "-X",
        f"pycache_prefix={PROMOTION_CACHE_ROOTS[mode]}",
        str(PROMOTION_TOOL),
        mode,
    ]


def verifier_command() -> list[str]:
    return [str(PYTHON), "-I", "-B", str(CONTINUATION_VERIFIER), "--verify"]


def authorized_verifier_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(AUTHORIZED_CONTINUATION_VERIFIER),
        "--verify",
    ]


def replay_command() -> list[str]:
    return [str(PYTHON), "-I", "-B", str(RUNNER_GATE_REPLAY)]


def authorized_replay_command() -> list[str]:
    return [str(PYTHON), "-I", "-B", str(AUTHORIZED_RUNNER_GATE_REPLAY)]


def validate_v7(authority_module: Any) -> dict[str, Any]:
    validation = authority_module.validate_release_source_authority(
        authority_module.EXPECTED_SOURCE_PATHS
    )
    lease = validation.get("validation_binding_lease")
    if (
        validation.get("pass") is not True
        or validation.get("authority_id") != EXPECTED_V7_AUTHORITY_ID
        or validation.get("source_set_sha256") != EXPECTED_V7_SOURCE_SET_SHA256
        or not isinstance(validation.get("all_source_roles_verified"), list)
        or len(validation["all_source_roles_verified"]) != 23
        or not isinstance(lease, Mapping)
        or lease.get("policy") != "BOUND_FILE_DESCRIPTORS_RETAINED_UNTIL_PROCESS_EXIT"
        or type(lease.get("descriptor_count")) is not int
        or lease["descriptor_count"] <= 0
    ):
        raise ContinuationError("v7 release authority validation drift")
    return validation


def require_absent(path: Path, label: str) -> None:
    if os.path.lexists(path):
        raise ContinuationError(f"{label} must be absent: {path}")


def require_incident_absent() -> None:
    require_absent(PROMOTION_INCIDENT, "Promotion incident")
    require_absent(CONTINUATION_INCIDENT, "Downstream continuation incident")
    require_absent(
        LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
        "Legacy authorization-payload signature",
    )
    require_absent(
        LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
        "Legacy completion-payload signature",
    )


def validate_producer_session(
    value: Any, label: str, *, expected_mode: str
) -> dict[str, Any]:
    require_exact_keys(value, PRODUCER_SESSION_KEYS, f"{label} producer_session")
    if (
        value.get("protocol") != "intersubmod.parent_fork_no_exec_child"
        or value.get("schema_version") != "1.0.0"
        or value.get("mode") != expected_mode
        or not is_canonical_sha256(value.get("nonce_sha256"))
        or type(value.get("parent_pid")) is not int
        or value["parent_pid"] <= 0
        or type(value.get("parent_start_time_ticks")) is not int
        or value["parent_start_time_ticks"] <= 0
        or type(value.get("child_pid")) is not int
        or value["child_pid"] <= 0
        or value["child_pid"] == value["parent_pid"]
        or type(value.get("child_start_time_ticks")) is not int
        or value["child_start_time_ticks"] <= 0
        or value.get("fork_no_exec") is not True
    ):
        raise ContinuationError(f"{label} producer_session contract drift")
    return dict(value)


def validate_producer_child_wait(
    value: Any, session: Mapping[str, Any], label: str
) -> dict[str, Any]:
    require_exact_keys(value, PRODUCER_CHILD_WAIT_KEYS, f"{label} child_wait")
    raw_status = value.get("raw_wait_status")
    if type(raw_status) is not int or raw_status < 0:
        raise ContinuationError(f"{label} child_wait raw status type drift")
    exited = os.WIFEXITED(raw_status)
    signaled = os.WIFSIGNALED(raw_status)
    decoded = {
        "waited_pid": value.get("waited_pid"),
        "raw_wait_status": raw_status,
        "wifexited": bool(exited),
        "exit_status": os.WEXITSTATUS(raw_status) if exited else None,
        "wifsignaled": bool(signaled),
        "terminating_signal": os.WTERMSIG(raw_status) if signaled else None,
        "core_dumped": bool(os.WCOREDUMP(raw_status)) if signaled else False,
    }
    if (
        not strict_json_equal(value, decoded)
        or value["raw_wait_status"] != 0
        or type(value.get("waited_pid")) is not int
        or value["waited_pid"] != session["child_pid"]
        or value.get("wifexited") is not True
        or type(value.get("exit_status")) is not int
        or value["exit_status"] != 0
        or value.get("wifsignaled") is not False
        or value.get("terminating_signal") is not None
        or value.get("core_dumped") is not False
    ):
        raise ContinuationError(f"{label} child_wait does not prove normal exit zero")
    return dict(value)


def producer_exit_attestation_announcement(
    path: Path,
    schema: Mapping[str, str],
    public_key: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "path": str(path),
        "schema": dict(schema),
        "signing_key": dict(public_key),
    }


def producer_exit_signature_contract(
    *,
    path: Path,
    schema: Mapping[str, str],
    signature: Path,
    public_key: Mapping[str, Any],
    private_key: Path,
    legacy_payload_signature: Path,
) -> dict[str, Any]:
    return {
        "attestation_schema": dict(schema),
        "algorithm": "Ed25519",
        "signed_artifact": str(path),
        "signature": str(signature),
        "public_key": dict(public_key),
        "private_key_lifecycle": {
            "pre_signature": {"path": str(private_key), "mode": "0o400"},
            "required_post_signature_mode": "0o0",
            "signing_must_follow_attestation_publication": True,
        },
        "legacy_payload_signature": {
            "path": str(legacy_payload_signature),
            "required_absent": True,
        },
    }


def validate_producer_exit_attestation(
    value: Any,
    *,
    label: str,
    expected_schema: Mapping[str, str],
    expected_phase: str,
    expected_command: list[str],
    expected_session: Mapping[str, Any],
    expected_artifacts: Mapping[str, Any],
    expected_source: Mapping[str, Any],
    expected_runtime: Mapping[str, Any],
    expected_signature_contract: Mapping[str, Any],
    expected_checks: Mapping[str, Any],
    expected_pass_semantics: str,
) -> tuple[dict[str, Any], datetime]:
    require_exact_keys(value, PRODUCER_EXIT_ATTESTATION_KEYS, label)
    created = require_utc_timestamp(value.get("created_at_utc"), f"{label} timestamp")
    session = validate_producer_session(
        value.get("producer_session"), label, expected_mode=expected_command[-1]
    )
    validate_producer_child_wait(value.get("child_wait"), session, label)
    if (
        value.get("authorization_id") != AUTHORIZATION_ID
        or value.get("phase") != expected_phase
        or not strict_json_equal(value.get("supervisor_command"), expected_command)
        or not strict_json_equal(session, dict(expected_session))
        or not strict_json_equal(value.get("produced_artifacts"), expected_artifacts)
        or not strict_json_equal(value.get("producer_source"), expected_source)
        or not strict_json_equal(value.get("python_runtime"), expected_runtime)
        or not strict_json_equal(
            value.get("signature_contract"), expected_signature_contract
        )
        or not strict_json_equal(
            value.get("signature_contract", {}).get("attestation_schema"),
            dict(expected_schema),
        )
        or not strict_json_equal(value.get("checks"), expected_checks)
        or value.get("pass") is not True
        or value.get("pass_semantics") != expected_pass_semantics
    ):
        raise ContinuationError(f"{label} exact contract drift")
    return dict(value), created


def validate_promotion_chain(
    authority_module: Any,
    reader: Any,
    v7_validation: Mapping[str, Any],
    bootstrap_record: Mapping[str, Any],
    *,
    continuation_private_live_mode: str = "0o400",
) -> dict[str, Any]:
    require_incident_absent()
    self_fd, self_data, self_record, self_stat = open_bound(reader, CONTINUATION_RUNNER)
    if self_record["mode"] != "0o444" or self_stat.st_nlink != 1:
        raise ContinuationError("Continuation source is not immutable mode 0444")
    self_binding = verify_executing_source(self_data)
    _, _, authorized_self_record, _ = open_bound(
        reader, AUTHORIZED_CONTINUATION_RUNNER
    )

    promotion_fd, promotion_data, promotion_record, promotion_stat = open_bound(
        reader, PROMOTION_TOOL
    )
    del promotion_fd
    promotion_full_record = with_stat(promotion_record, promotion_stat)
    verifier_fd, verifier_data, verifier_record, _ = open_bound(
        reader, CONTINUATION_VERIFIER
    )
    _, _, authorized_verifier_record, _ = open_bound(
        reader, AUTHORIZED_CONTINUATION_VERIFIER
    )
    replay_fd, replay_source_data, replay_record, _ = open_bound(reader, RUNNER_GATE_REPLAY)
    _, _, authorized_replay_record, _ = open_bound(
        reader, AUTHORIZED_RUNNER_GATE_REPLAY
    )
    del replay_source_data
    runner_fd, runner_data, runner_record, _ = open_bound(reader, COMPLETION_RUNNER)
    authority_fd, authority_data, authority_record, _ = open_bound(
        reader, RELEASE_SOURCE_AUTHORITY
    )
    del authority_fd
    if not strict_json_equal(authority_record, bootstrap_record) or sha256_bytes(authority_data) != (
        EXPECTED_RELEASE_SOURCE_AUTHORITY_SHA256
    ):
        raise ContinuationError("Executed release authority differs from bound source")
    _, recovery_validator_data, recovery_validator_record, _ = open_bound(
        reader, SCHEMA_RECOVERY_VALIDATOR
    )
    if recovery_validator_record.get("sha256") != (
        EXPECTED_SCHEMA_RECOVERY_VALIDATOR_SHA256
    ):
        raise ContinuationError("Schema-recovery authority validator identity drift")
    recovery_module = types.ModuleType(
        "bound_tumor_ref_schema_recovery_authority_for_continuation_v2"
    )
    recovery_module.__file__ = str(SCHEMA_RECOVERY_VALIDATOR)
    recovery_module.__package__ = ""
    exec(
        compile(
            recovery_validator_data,
            str(SCHEMA_RECOVERY_VALIDATOR),
            "exec",
            dont_inherit=True,
        ),
        recovery_module.__dict__,
    )
    if not hasattr(recovery_module, "validate_recovery_authority"):
        raise ContinuationError("Schema-recovery authority validator API is incomplete")
    recovery_authority = recovery_module.validate_recovery_authority(
        expected_runtime_role="downstream_continuation",
        expected_recovery_continuation_private_mode=continuation_private_live_mode,
    )
    if (
        recovery_authority.get("pass") is not True
        or recovery_authority.get("authority_id")
        != "20260724_tumor_ref_schema_recovery_v30"
        or not strict_json_equal(
            recovery_authority.get("runtime_source"), self_record
        )
    ):
        raise ContinuationError("Schema-recovery authority does not bind this continuation")

    openssl_fd, _, openssl_record, _ = open_bound(reader, OPENSSL, expected_mode="0o755")
    if openssl_record["sha256"] != EXPECTED_OPENSSL_SHA256:
        raise ContinuationError("OpenSSL verifier identity drift")
    python_fd, python_record, python_stat = bind_running_python(reader)
    python_full_record = with_stat(python_record, python_stat)
    reviewed_runtime_sources, reviewed_runtime_fds = bind_reviewed_runtime_sources(
        reader,
        primary_python_fd=python_fd,
        primary_python_record=python_record,
    )
    historical_runtime_sources = historical_signed_runtime_projection(
        reviewed_runtime_sources
    )

    auth_fd, auth_data, auth_basic, auth_stat = open_bound(reader, AUTHORIZATION)
    auth_record = with_stat(auth_basic, auth_stat)
    prepare_attestation_fd, prepare_attestation_data, prepare_attestation_basic, prepare_attestation_stat = (
        open_bound(reader, PREPARE_EXIT_ATTESTATION)
    )
    prepare_attestation_record = with_stat(
        prepare_attestation_basic, prepare_attestation_stat
    )
    auth_sig_fd, _, auth_sig_basic, auth_sig_stat = open_bound(
        reader, AUTHORIZATION_SIGNATURE
    )
    auth_sig_record = with_stat(auth_sig_basic, auth_sig_stat)
    auth_public_fd, _, auth_public_record, _ = open_bound(
        reader, AUTHORIZATION_PUBLIC_KEY
    )
    auth_private_fd, auth_private_record, _ = open_private_metadata(
        reader, AUTHORIZATION_PRIVATE_KEY, "0o0"
    )
    del auth_private_fd
    completion_fd, completion_data, completion_basic, completion_stat = open_bound(
        reader, COMPLETION
    )
    completion_record = with_stat(completion_basic, completion_stat)
    promote_attestation_fd, promote_attestation_data, promote_attestation_basic, promote_attestation_stat = (
        open_bound(reader, PROMOTE_EXIT_ATTESTATION)
    )
    promote_attestation_record = with_stat(
        promote_attestation_basic, promote_attestation_stat
    )
    completion_sig_fd, _, completion_sig_basic, completion_sig_stat = open_bound(
        reader, COMPLETION_SIGNATURE
    )
    completion_sig_record = with_stat(completion_sig_basic, completion_sig_stat)
    completion_public_fd, _, completion_public_record, _ = open_bound(
        reader, COMPLETION_PUBLIC_KEY
    )
    completion_private_fd, completion_private_record, _ = open_private_metadata(
        reader, COMPLETION_PRIVATE_KEY, "0o0"
    )
    del completion_private_fd
    continuation_public_fd, _, continuation_public_record, continuation_public_stat = (
        open_bound(reader, CONTINUATION_PUBLIC_KEY)
    )
    del continuation_public_fd
    (
        continuation_private_fd,
        continuation_private_live_record,
        continuation_private_stat,
    ) = open_private_metadata(
        reader, CONTINUATION_PRIVATE_KEY, "0o400"
    )
    del continuation_private_fd
    continuation_private_record = {
        "path": str(CONTINUATION_PRIVATE_KEY),
        "mode": "0o400",
    }
    (
        recovery_continuation_public_fd,
        _,
        recovery_continuation_public_record,
        recovery_continuation_public_stat,
    ) = open_bound(reader, RECOVERY_CONTINUATION_PUBLIC_KEY)
    del recovery_continuation_public_fd
    (
        recovery_continuation_private_fd,
        recovery_continuation_private_live_record,
        recovery_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        RECOVERY_CONTINUATION_PRIVATE_KEY,
        continuation_private_live_mode,
    )
    del recovery_continuation_private_fd
    recovery_continuation_private_record = {
        "path": str(RECOVERY_CONTINUATION_PRIVATE_KEY),
        "mode": "0o400",
    }
    (
        failed_v16_continuation_public_fd,
        _,
        failed_v16_continuation_public_record,
        failed_v16_continuation_public_stat,
    ) = open_bound(reader, FAILED_V16_CONTINUATION_PUBLIC_KEY)
    del failed_v16_continuation_public_fd
    (
        failed_v16_continuation_private_fd,
        failed_v16_continuation_private_live_record,
        failed_v16_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        FAILED_V16_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del failed_v16_continuation_private_fd
    (
        failed_v17_continuation_public_fd,
        _,
        failed_v17_continuation_public_record,
        failed_v17_continuation_public_stat,
    ) = open_bound(reader, FAILED_V17_CONTINUATION_PUBLIC_KEY)
    del failed_v17_continuation_public_fd
    (
        failed_v17_continuation_private_fd,
        failed_v17_continuation_private_live_record,
        failed_v17_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        FAILED_V17_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del failed_v17_continuation_private_fd
    (
        failed_v18_continuation_public_fd,
        _,
        failed_v18_continuation_public_record,
        failed_v18_continuation_public_stat,
    ) = open_bound(reader, FAILED_V18_CONTINUATION_PUBLIC_KEY)
    del failed_v18_continuation_public_fd
    (
        failed_v18_continuation_private_fd,
        failed_v18_continuation_private_live_record,
        failed_v18_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        FAILED_V18_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del failed_v18_continuation_private_fd
    (
        failed_v19_continuation_public_fd,
        _,
        failed_v19_continuation_public_record,
        failed_v19_continuation_public_stat,
    ) = open_bound(reader, FAILED_V19_CONTINUATION_PUBLIC_KEY)
    del failed_v19_continuation_public_fd
    (
        failed_v19_continuation_private_fd,
        failed_v19_continuation_private_live_record,
        failed_v19_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        FAILED_V19_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del failed_v19_continuation_private_fd
    (
        rejected_v20_continuation_public_fd,
        _,
        rejected_v20_continuation_public_record,
        rejected_v20_continuation_public_stat,
    ) = open_bound(reader, REJECTED_V20_CONTINUATION_PUBLIC_KEY)
    del rejected_v20_continuation_public_fd
    (
        rejected_v20_continuation_private_fd,
        rejected_v20_continuation_private_live_record,
        rejected_v20_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        REJECTED_V20_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del rejected_v20_continuation_private_fd
    (
        failed_v21_continuation_public_fd,
        _,
        failed_v21_continuation_public_record,
        failed_v21_continuation_public_stat,
    ) = open_bound(reader, FAILED_V21_CONTINUATION_PUBLIC_KEY)
    del failed_v21_continuation_public_fd
    (
        failed_v21_continuation_private_fd,
        failed_v21_continuation_private_live_record,
        failed_v21_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        FAILED_V21_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del failed_v21_continuation_private_fd
    (
        rejected_v22_continuation_public_fd,
        _,
        rejected_v22_continuation_public_record,
        rejected_v22_continuation_public_stat,
    ) = open_bound(reader, REJECTED_V22_CONTINUATION_PUBLIC_KEY)
    del rejected_v22_continuation_public_fd
    (
        rejected_v22_continuation_private_fd,
        rejected_v22_continuation_private_live_record,
        rejected_v22_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        REJECTED_V22_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del rejected_v22_continuation_private_fd
    (
        rejected_v23_continuation_public_fd,
        _,
        rejected_v23_continuation_public_record,
        rejected_v23_continuation_public_stat,
    ) = open_bound(reader, REJECTED_V23_CONTINUATION_PUBLIC_KEY)
    del rejected_v23_continuation_public_fd
    (
        rejected_v23_continuation_private_fd,
        rejected_v23_continuation_private_live_record,
        rejected_v23_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        REJECTED_V23_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del rejected_v23_continuation_private_fd
    (
        rejected_v24_continuation_public_fd,
        _,
        rejected_v24_continuation_public_record,
        rejected_v24_continuation_public_stat,
    ) = open_bound(reader, REJECTED_V24_CONTINUATION_PUBLIC_KEY)
    del rejected_v24_continuation_public_fd
    (
        rejected_v24_continuation_private_fd,
        rejected_v24_continuation_private_live_record,
        rejected_v24_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        REJECTED_V24_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del rejected_v24_continuation_private_fd
    (
        failed_v25_continuation_public_fd,
        _,
        failed_v25_continuation_public_record,
        failed_v25_continuation_public_stat,
    ) = open_bound(reader, FAILED_V25_CONTINUATION_PUBLIC_KEY)
    del failed_v25_continuation_public_fd
    (
        failed_v25_continuation_private_fd,
        failed_v25_continuation_private_live_record,
        failed_v25_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        FAILED_V25_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del failed_v25_continuation_private_fd
    (
        failed_v26_continuation_public_fd,
        _,
        failed_v26_continuation_public_record,
        failed_v26_continuation_public_stat,
    ) = open_bound(reader, FAILED_V26_CONTINUATION_PUBLIC_KEY)
    del failed_v26_continuation_public_fd
    (
        failed_v26_continuation_private_fd,
        failed_v26_continuation_private_live_record,
        failed_v26_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        FAILED_V26_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del failed_v26_continuation_private_fd
    (
        rejected_v27_continuation_public_fd,
        _,
        rejected_v27_continuation_public_record,
        rejected_v27_continuation_public_stat,
    ) = open_bound(reader, REJECTED_V27_CONTINUATION_PUBLIC_KEY)
    del rejected_v27_continuation_public_fd
    (
        rejected_v27_continuation_private_fd,
        rejected_v27_continuation_private_live_record,
        rejected_v27_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        REJECTED_V27_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del rejected_v27_continuation_private_fd
    (
        failed_v28_continuation_public_fd,
        _,
        failed_v28_continuation_public_record,
        failed_v28_continuation_public_stat,
    ) = open_bound(reader, FAILED_V28_CONTINUATION_PUBLIC_KEY)
    del failed_v28_continuation_public_fd
    (
        failed_v28_continuation_private_fd,
        failed_v28_continuation_private_live_record,
        failed_v28_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        FAILED_V28_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del failed_v28_continuation_private_fd
    (
        failed_v29_continuation_public_fd,
        _,
        failed_v29_continuation_public_record,
        failed_v29_continuation_public_stat,
    ) = open_bound(reader, FAILED_V29_CONTINUATION_PUBLIC_KEY)
    del failed_v29_continuation_public_fd
    (
        failed_v29_continuation_private_fd,
        failed_v29_continuation_private_live_record,
        failed_v29_continuation_private_stat,
    ) = open_private_metadata(
        reader,
        FAILED_V29_CONTINUATION_PRIVATE_KEY,
        "0o400",
    )
    del failed_v29_continuation_private_fd
    expected_terminal_key_state = {
        "legacy_public_key": continuation_public_record,
        "legacy_private_key": {
            "path": str(CONTINUATION_PRIVATE_KEY),
            "size_bytes": continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "LEGACY_SIGNED_CONTRACT_UNUSED_NOT_RETIRED",
        },
        "failed_v16_public_key": failed_v16_continuation_public_record,
        "failed_v16_private_key": {
            "path": str(FAILED_V16_CONTINUATION_PRIVATE_KEY),
            "size_bytes": failed_v16_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "FAILED_V16_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        },
        "failed_v17_public_key": failed_v17_continuation_public_record,
        "failed_v17_private_key": {
            "path": str(FAILED_V17_CONTINUATION_PRIVATE_KEY),
            "size_bytes": failed_v17_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "FAILED_V17_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        },
        "failed_v18_public_key": failed_v18_continuation_public_record,
        "failed_v18_private_key": {
            "path": str(FAILED_V18_CONTINUATION_PRIVATE_KEY),
            "size_bytes": failed_v18_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "FAILED_V18_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        },
        "failed_v19_public_key": failed_v19_continuation_public_record,
        "failed_v19_private_key": {
            "path": str(FAILED_V19_CONTINUATION_PRIVATE_KEY),
            "size_bytes": failed_v19_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "FAILED_V19_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        },
        "rejected_v20_public_key": rejected_v20_continuation_public_record,
        "rejected_v20_private_key": {
            "path": str(REJECTED_V20_CONTINUATION_PRIVATE_KEY),
            "size_bytes": rejected_v20_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "REJECTED_V20_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        },
        "failed_v21_public_key": failed_v21_continuation_public_record,
        "failed_v21_private_key": {
            "path": str(FAILED_V21_CONTINUATION_PRIVATE_KEY),
            "size_bytes": failed_v21_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "FAILED_V21_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        },
        "rejected_v22_public_key": rejected_v22_continuation_public_record,
        "rejected_v22_private_key": {
            "path": str(REJECTED_V22_CONTINUATION_PRIVATE_KEY),
            "size_bytes": rejected_v22_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "REJECTED_V22_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        },
        "rejected_v23_public_key": rejected_v23_continuation_public_record,
        "rejected_v23_private_key": {
            "path": str(REJECTED_V23_CONTINUATION_PRIVATE_KEY),
            "size_bytes": rejected_v23_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "REJECTED_V23_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        },
        "rejected_v24_public_key": rejected_v24_continuation_public_record,
        "rejected_v24_private_key": {
            "path": str(REJECTED_V24_CONTINUATION_PRIVATE_KEY),
            "size_bytes": rejected_v24_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "REJECTED_V24_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        },
        "failed_v25_public_key": failed_v25_continuation_public_record,
        "failed_v25_private_key": {
            "path": str(FAILED_V25_CONTINUATION_PRIVATE_KEY),
            "size_bytes": failed_v25_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "FAILED_V25_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
        },
        "failed_v26_public_key": failed_v26_continuation_public_record,
        "failed_v26_private_key": {
            "path": str(FAILED_V26_CONTINUATION_PRIVATE_KEY),
            "size_bytes": failed_v26_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "FAILED_V26_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
        },
        "rejected_v27_public_key": rejected_v27_continuation_public_record,
        "rejected_v27_private_key": {
            "path": str(REJECTED_V27_CONTINUATION_PRIVATE_KEY),
            "size_bytes": rejected_v27_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "REJECTED_V27_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
        },
        "failed_v28_public_key": failed_v28_continuation_public_record,
        "failed_v28_private_key": {
            "path": str(FAILED_V28_CONTINUATION_PRIVATE_KEY),
            "size_bytes": failed_v28_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "FAILED_V28_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
        },
        "failed_v29_public_key": failed_v29_continuation_public_record,
        "failed_v29_private_key": {
            "path": str(FAILED_V29_CONTINUATION_PRIVATE_KEY),
            "size_bytes": failed_v29_continuation_private_stat.st_size,
            "mode": "0o400",
            "state": "FAILED_V29_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
        },
        "recovery_public_key": recovery_continuation_public_record,
        "recovery_private_key": {
            "path": str(RECOVERY_CONTINUATION_PRIVATE_KEY),
            "size_bytes": recovery_continuation_private_stat.st_size,
            "mode": continuation_private_live_mode,
            "state": (
                "RECOVERY_V30_READY_FOR_SINGLE_SIGNATURE"
                if continuation_private_live_mode == "0o400"
                else "RECOVERY_V30_RETIRED_AFTER_SINGLE_SIGNATURE"
            ),
        },
    }
    if (
        set(expected_terminal_key_state) != TERMINAL_KEY_STATE_KEYS
        or not isinstance(recovery_authority.get("terminal_key_rotation"), Mapping)
        or set(recovery_authority["terminal_key_rotation"])
        != TERMINAL_KEY_ROTATION_KEYS
        or not isinstance(recovery_authority.get("terminal_key_state"), Mapping)
        or set(recovery_authority["terminal_key_state"]) != TERMINAL_KEY_STATE_KEYS
        or auth_sig_stat.st_size != 64
        or completion_sig_stat.st_size != 64
        or auth_public_record["sha256"] != EXPECTED_AUTHORIZATION_PUBLIC_KEY_SHA256
        or completion_public_record["sha256"] != EXPECTED_COMPLETION_PUBLIC_KEY_SHA256
        or continuation_public_record["sha256"]
        != EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256
        or failed_v16_continuation_public_record["sha256"]
        != EXPECTED_FAILED_V16_CONTINUATION_PUBLIC_KEY_SHA256
        or failed_v17_continuation_public_record["sha256"]
        != EXPECTED_FAILED_V17_CONTINUATION_PUBLIC_KEY_SHA256
        or failed_v18_continuation_public_record["sha256"]
        != EXPECTED_FAILED_V18_CONTINUATION_PUBLIC_KEY_SHA256
        or failed_v19_continuation_public_record["sha256"]
        != EXPECTED_FAILED_V19_CONTINUATION_PUBLIC_KEY_SHA256
        or rejected_v20_continuation_public_record["sha256"]
        != EXPECTED_REJECTED_V20_CONTINUATION_PUBLIC_KEY_SHA256
        or failed_v21_continuation_public_record["sha256"]
        != EXPECTED_FAILED_V21_CONTINUATION_PUBLIC_KEY_SHA256
        or rejected_v22_continuation_public_record["sha256"]
        != EXPECTED_REJECTED_V22_CONTINUATION_PUBLIC_KEY_SHA256
        or rejected_v23_continuation_public_record["sha256"]
        != EXPECTED_REJECTED_V23_CONTINUATION_PUBLIC_KEY_SHA256
        or rejected_v24_continuation_public_record["sha256"]
        != EXPECTED_REJECTED_V24_CONTINUATION_PUBLIC_KEY_SHA256
        or failed_v25_continuation_public_record["sha256"]
        != EXPECTED_FAILED_V25_CONTINUATION_PUBLIC_KEY_SHA256
        or failed_v26_continuation_public_record["sha256"]
        != EXPECTED_FAILED_V26_CONTINUATION_PUBLIC_KEY_SHA256
        or rejected_v27_continuation_public_record["sha256"]
        != EXPECTED_REJECTED_V27_CONTINUATION_PUBLIC_KEY_SHA256
        or failed_v28_continuation_public_record["sha256"]
        != EXPECTED_FAILED_V28_CONTINUATION_PUBLIC_KEY_SHA256
        or failed_v29_continuation_public_record["sha256"]
        != EXPECTED_FAILED_V29_CONTINUATION_PUBLIC_KEY_SHA256
        or recovery_continuation_public_record["sha256"]
        != EXPECTED_RECOVERY_CONTINUATION_PUBLIC_KEY_SHA256
        or len(
            {
                continuation_public_record["sha256"],
                failed_v16_continuation_public_record["sha256"],
                failed_v17_continuation_public_record["sha256"],
                failed_v18_continuation_public_record["sha256"],
                failed_v19_continuation_public_record["sha256"],
                rejected_v20_continuation_public_record["sha256"],
                failed_v21_continuation_public_record["sha256"],
                rejected_v22_continuation_public_record["sha256"],
                rejected_v23_continuation_public_record["sha256"],
                rejected_v24_continuation_public_record["sha256"],
                failed_v25_continuation_public_record["sha256"],
                failed_v26_continuation_public_record["sha256"],
                rejected_v27_continuation_public_record["sha256"],
                failed_v28_continuation_public_record["sha256"],
                failed_v29_continuation_public_record["sha256"],
                recovery_continuation_public_record["sha256"],
            }
        )
        != 16
        or len(
            {
                auth_public_record["sha256"],
                completion_public_record["sha256"],
                continuation_public_record["sha256"],
            }
        )
        != 3
        or continuation_public_stat.st_nlink != 1
        or continuation_private_stat.st_nlink != 1
        or failed_v16_continuation_public_stat.st_nlink != 1
        or failed_v16_continuation_private_stat.st_nlink != 1
        or failed_v17_continuation_public_stat.st_nlink != 1
        or failed_v17_continuation_private_stat.st_nlink != 1
        or failed_v18_continuation_public_stat.st_nlink != 1
        or failed_v18_continuation_private_stat.st_nlink != 1
        or failed_v19_continuation_public_stat.st_nlink != 1
        or failed_v19_continuation_private_stat.st_nlink != 1
        or rejected_v20_continuation_public_stat.st_nlink != 1
        or rejected_v20_continuation_private_stat.st_nlink != 1
        or failed_v21_continuation_public_stat.st_nlink != 1
        or failed_v21_continuation_private_stat.st_nlink != 1
        or rejected_v22_continuation_public_stat.st_nlink != 1
        or rejected_v22_continuation_private_stat.st_nlink != 1
        or rejected_v23_continuation_public_stat.st_nlink != 1
        or rejected_v23_continuation_private_stat.st_nlink != 1
        or rejected_v24_continuation_public_stat.st_nlink != 1
        or rejected_v24_continuation_private_stat.st_nlink != 1
        or failed_v25_continuation_public_stat.st_nlink != 1
        or failed_v25_continuation_private_stat.st_nlink != 1
        or failed_v26_continuation_public_stat.st_nlink != 1
        or failed_v26_continuation_private_stat.st_nlink != 1
        or rejected_v27_continuation_public_stat.st_nlink != 1
        or rejected_v27_continuation_private_stat.st_nlink != 1
        or failed_v28_continuation_public_stat.st_nlink != 1
        or failed_v28_continuation_private_stat.st_nlink != 1
        or failed_v29_continuation_public_stat.st_nlink != 1
        or failed_v29_continuation_private_stat.st_nlink != 1
        or recovery_continuation_public_stat.st_nlink != 1
        or recovery_continuation_private_stat.st_nlink != 1
        or not strict_json_equal(
            recovery_authority.get("terminal_key_state"),
            expected_terminal_key_state,
        )
    ):
        raise ContinuationError("Promotion signature/key identity drift")
    if not authority_module.verify_ed25519_signature_fds(
        openssl_fd=openssl_fd,
        data_fd=prepare_attestation_fd,
        public_key_fd=auth_public_fd,
        signature_fd=auth_sig_fd,
    ):
        raise ContinuationError("Promotion authorization signature verification failed")
    if not authority_module.verify_ed25519_signature_fds(
        openssl_fd=openssl_fd,
        data_fd=promote_attestation_fd,
        public_key_fd=completion_public_fd,
        signature_fd=completion_sig_fd,
    ):
        raise ContinuationError("Promotion completion signature verification failed")

    source_fd, source_data, source_basic, source_stat = open_bound(
        reader, HISTORICAL_SOURCE_RECEIPT
    )
    del source_fd
    canonical_fd, canonical_data, canonical_basic, canonical_stat = open_bound(
        reader, CANONICAL_SOURCE_RECEIPT
    )
    del canonical_fd
    source_record = with_stat(source_basic, source_stat)
    canonical_record = with_stat(canonical_basic, canonical_stat)
    if (
        len(source_data) != EXPECTED_SOURCE_RECEIPT_SIZE
        or source_basic["sha256"] != EXPECTED_SOURCE_RECEIPT_SHA256
        or canonical_basic["sha256"] != EXPECTED_SOURCE_RECEIPT_SHA256
        or source_data != canonical_data
        or source_stat.st_nlink != 1
        or canonical_stat.st_nlink != 1
        or (source_stat.st_dev, source_stat.st_ino)
        == (canonical_stat.st_dev, canonical_stat.st_ino)
    ):
        raise ContinuationError("Historical/canonical tumor-REF receipt identity drift")
    source_payload = strict_json(source_data, "historical source receipt")
    if (
        source_payload.get("schema_name")
        != "intersubmod.retrospective_running_source_identity.receipt"
        or source_payload.get("schema_version") != "1.2.0"
        or source_payload.get("task_type") != "B_comprehensive_validation"
        or source_payload.get("pass") is not True
    ):
        raise ContinuationError("Historical source receipt semantic contract drift")

    reviewed_sources = {
        "promotion_tool": promotion_record,
        "continuation_verifier": authorized_verifier_record,
        "runner_gate_replay": authorized_replay_record,
        "downstream_continuation": authorized_self_record,
        **historical_runtime_sources,
    }
    trusted_keys = {
        "authorization_public_key": auth_public_record,
        "completion_public_key": completion_public_record,
        "continuation_public_key": continuation_public_record,
    }
    authorization = strict_json(auth_data, "signed promotion authorization")
    completion = strict_json(completion_data, "signed promotion completion")
    require_exact_keys(authorization, AUTHORIZATION_KEYS, "promotion authorization")
    require_exact_keys(completion, COMPLETION_KEYS, "promotion completion")
    prepare_session = validate_producer_session(
        authorization.get("producer_session"),
        "promotion authorization",
        expected_mode="--prepare-authorization",
    )
    promote_session = validate_producer_session(
        completion.get("producer_session"),
        "promotion completion",
        expected_mode="--promote",
    )
    if prepare_session["nonce_sha256"] == promote_session["nonce_sha256"]:
        raise ContinuationError("Prepare/promote producer sessions reused one nonce")
    auth_time = require_utc_timestamp(
        authorization.get("created_at_utc"), "authorization created_at_utc"
    )
    completion_time = require_utc_timestamp(
        completion.get("created_at_utc"), "completion created_at_utc"
    )
    if completion_time < auth_time:
        raise ContinuationError("Promotion completion predates authorization")

    expected_commands = {
        "prepare_authorization": promotion_command("--prepare-authorization"),
        "promote": promotion_command("--promote"),
        "continuation_verify": authorized_verifier_command(),
        "runner_gate_replay": authorized_replay_command(),
        "downstream_continuation": authorized_canonical_command(),
    }
    expected_continuation_gate = {
        "verifier": authorized_verifier_record,
        "verification_receipt": str(AUTHORIZED_PROMOTION_VERIFICATION_RECEIPT),
        "runner_gate_replay": authorized_replay_record,
        "runner_gate_receipt": str(AUTHORIZED_RUNNER_GATE_RECEIPT),
        "downstream_continuation": authorized_self_record,
        "downstream_continuation_command": authorized_canonical_command(),
        "supervised_child_command": authorized_supervised_child_command(),
        "supervisor_capability_protocol": SUPERVISOR_CAPABILITY_PROTOCOL,
        "threat_model": CONTINUATION_THREAT_MODEL,
        "terminal_signature_contract": {
            "algorithm": "Ed25519",
            "public_key": continuation_public_record,
            "private_key": continuation_private_record,
            "execution_receipt": str(AUTHORIZED_CONTINUATION_RECEIPT),
            "signed_artifact": str(AUTHORIZED_CONTINUATION_EXIT_ATTESTATION),
            "signature": str(AUTHORIZED_CONTINUATION_SIGNATURE),
            "required_private_key_pre_signature_mode": "0o400",
            "required_private_key_post_signature_mode": "0o0",
            "supervisor_command": authorized_canonical_command(),
            "supervised_child_command": authorized_supervised_child_command(),
            "signed_terminal_verification_command": (
                authorized_signed_terminal_verification_command()
            ),
        },
        "policy": CONTINUATION_POLICY,
    }
    expected_auth_contract = producer_exit_signature_contract(
        path=PREPARE_EXIT_ATTESTATION,
        schema=PREPARE_EXIT_ATTESTATION_SCHEMA,
        signature=AUTHORIZATION_SIGNATURE,
        public_key=auth_public_record,
        private_key=AUTHORIZATION_PRIVATE_KEY,
        legacy_payload_signature=LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
    )
    expected_completion_contract = producer_exit_signature_contract(
        path=PROMOTE_EXIT_ATTESTATION,
        schema=PROMOTE_EXIT_ATTESTATION_SCHEMA,
        signature=COMPLETION_SIGNATURE,
        public_key=completion_public_record,
        private_key=COMPLETION_PRIVATE_KEY,
        legacy_payload_signature=LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
    )
    expected_disclosure = {
        "original_completion_runner_end_to_end_pass_claimed": False,
        "original_post_run_verifier_success_claimed": False,
        "canonical_receipt_will_be_created_by": (
            "signed_v3_exit_attested_byte_identical_promotion"
        ),
        "failed_attempt_preserved": True,
        "narrow_policy_supersession": (
            "Only the prior requirement that this canonical path be newly emitted by the "
            "historical completion-run verifier is superseded. Receipt content checks and "
            "all downstream release gates remain mandatory."
        ),
    }
    if (
        authorization.get("schema_name") != AUTHORIZATION_SCHEMA
        or authorization.get("schema_version") != SCHEMA_VERSION
        or authorization.get("authorization_id") != AUTHORIZATION_ID
        or authorization.get("task_type") != "B_comprehensive_validation"
        or authorization.get("status") != "AUTHORIZED_FOR_SINGLE_BYTE_IDENTICAL_PROMOTION"
        or authorization.get("scope")
        != "canonical_path_publication_only_no_scientific_payload_change"
        or not strict_json_equal(authorization.get("commands"), expected_commands)
        or not strict_json_equal(
            authorization.get("reviewed_sources"), reviewed_sources
        )
        or not strict_json_equal(
            authorization.get("trusted_signing_keys"), trusted_keys
        )
        or not strict_json_equal(authorization.get("source_receipt"), source_record)
        or authorization.get("canonical_target_path") != str(CANONICAL_SOURCE_RECEIPT)
        or not strict_json_equal(
            authorization.get("producer_exit_attestation"),
            producer_exit_attestation_announcement(
                PREPARE_EXIT_ATTESTATION,
                PREPARE_EXIT_ATTESTATION_SCHEMA,
                auth_public_record,
            ),
        )
        or set(authorization.get("evidence", {})) != set(AUTHORIZATION_EVIDENCE_KEYS)
        or not strict_json_equal(
            authorization.get("authorization_signature_contract"), expected_auth_contract
        )
        or not strict_json_equal(
            authorization.get("completion_signature_contract"),
            expected_completion_contract,
        )
        or not strict_json_equal(
            authorization.get("continuation_gate"), expected_continuation_gate
        )
        or not strict_json_equal(
            authorization.get("historical_incident_disclosure"), expected_disclosure
        )
        or not strict_json_equal(
            authorization.get("downstream_output_absence"),
            [str(path) for path in AUTHORIZED_DOWNSTREAM_OUTPUT_SLOTS],
        )
        or not strict_json_equal(authorization.get("checks"), AUTHORIZATION_CHECKS)
        or authorization.get("pass") is not True
        or authorization.get("pass_semantics")
        != (
            "authorization_payload_only; parent exit-attestation publication and detached "
            "authorization-key signature remain mandatory; no canonical publication has "
            "occurred and this payload does not permit downstream execution"
        )
    ):
        raise ContinuationError("Signed promotion authorization exact contract drift")
    recomputed_binding = recompute_execution_binding(promotion_data, str(PROMOTION_TOOL))
    if not strict_json_equal(authorization.get("execution_binding"), recomputed_binding):
        raise ContinuationError("Promotion execution binding was not independently reproduced")
    prepare_attestation_payload = strict_json(
        prepare_attestation_data, "signed prepare producer exit attestation"
    )
    prepare_attestation_payload, prepare_attestation_time = (
        validate_producer_exit_attestation(
            prepare_attestation_payload,
            label="prepare producer exit attestation",
            expected_schema=PREPARE_EXIT_ATTESTATION_SCHEMA,
            expected_phase="prepare_authorization",
            expected_command=promotion_command("--prepare-authorization"),
            expected_session=prepare_session,
            expected_artifacts={"authorization": auth_record},
            expected_source=promotion_full_record,
            expected_runtime=python_full_record,
            expected_signature_contract=expected_auth_contract,
            expected_checks=PREPARE_EXIT_CHECKS,
            expected_pass_semantics=PREPARE_EXIT_PASS_SEMANTICS,
        )
    )
    if prepare_attestation_time < auth_time:
        raise ContinuationError("Prepare producer exit attestation predates authorization")

    expected_authorization_link = {
        "artifact": auth_record,
        "producer_exit_attestation": prepare_attestation_record,
        "signature": auth_sig_record,
        "public_key": auth_public_record,
        "retired_private_key": auth_private_record,
        "producer_session": authorization["producer_session"],
        "signature_verified": True,
    }
    expected_completion_signature_contract = expected_completion_contract
    expected_builder_gate = {
        "source_path_status": "VERIFIED_BOUNDED_RETROSPECTIVE_SOURCE_IDENTITY",
        "canonical_path_status": "VERIFIED_BOUNDED_RETROSPECTIVE_SOURCE_IDENTITY",
        "release_gate_pass": True,
        "publishable_task_b_release": True,
        "builder_loaded_from_v7_bound_bytes": True,
        "bytecode_disabled": True,
    }
    expected_completion_disclosure = {
        "original_completion_runner_end_to_end_pass_claimed": False,
        "original_post_run_verifier_success_claimed": False,
        "canonical_receipt_created_by": (
            "signed_v3_exit_attested_byte_identical_promotion"
        ),
        "failed_attempt_archive": authorization["evidence"]["failed_attempt_archive"],
        "observed_failure": (
            "The historical runner attempt stopped because full-stat equality rejected a "
            "later 0664-to-0444 source permission tightening. No strict/matched/CN/final "
            "downstream output was produced by that attempt."
        ),
        "claim_limit": (
            "This write-ahead receipt is valid only with the exact live canonical copy, "
            "detached completion signature, fresh verifier, replay evidence, and signed "
            "continuation terminal receipt; it does not claim that the historical runner "
            "or its post-run verifier completed successfully."
        ),
    }
    governance = completion.get("governance")
    if not isinstance(governance, Mapping) or set(governance) != {
        "prior_review",
        "narrow_supersession",
        "continuation_verifier",
        "runner_gate_replay",
        "downstream_continuation",
        "downstream_continuation_command",
        "supervised_child_command",
        "continuation_verification_receipt",
        "runner_gate_receipt",
        "downstream_output_slots_absent_at_promotion",
        "canonical_publication_order",
        "canonical_alone_grants_downstream_authority",
        "write_ahead_receipt_marks_transaction_incomplete_until",
    }:
        raise ContinuationError("Promotion completion governance schema drift")
    if (
        completion.get("schema_name") != COMPLETION_SCHEMA
        or completion.get("schema_version") != SCHEMA_VERSION
        or completion.get("authorization_id") != AUTHORIZATION_ID
        or completion.get("task_type") != "B_comprehensive_validation"
        or completion.get("status")
        != "PROVISIONAL_WRITE_AHEAD_CANONICAL_PUBLICATION_INTENT"
        or completion.get("scope")
        != "canonical_path_publication_only_no_scientific_payload_change"
        or not strict_json_equal(
            completion.get("command"), promotion_command("--promote")
        )
        or not strict_json_equal(
            completion.get("producer_exit_attestation"),
            producer_exit_attestation_announcement(
                PROMOTE_EXIT_ATTESTATION,
                PROMOTE_EXIT_ATTESTATION_SCHEMA,
                completion_public_record,
            ),
        )
        or not strict_json_equal(
            completion.get("authorization"), expected_authorization_link
        )
        or not strict_json_equal(completion.get("code"), reviewed_sources)
        or not strict_json_equal(
            completion.get("execution_binding"), recomputed_binding
        )
        or not strict_json_equal(completion.get("source_receipt"), source_record)
        or not strict_json_equal(
            completion.get("canonical_receipt"),
            {
                key: canonical_record[key]
                for key in ("path", "size_bytes", "sha256", "mode")
            },
        )
        or not strict_json_equal(
            completion.get("evidence"), authorization.get("evidence")
        )
        or not strict_json_equal(
            completion.get("review_artifacts"), authorization.get("review_artifacts")
        )
        or not strict_json_equal(
            completion.get("builder_gate_replay"), expected_builder_gate
        )
        or not strict_json_equal(
            completion.get("historical_incident_disclosure"),
            expected_completion_disclosure,
        )
        or not strict_json_equal(
            completion.get("signature_contract"),
            expected_completion_signature_contract,
        )
        or not strict_json_equal(completion.get("checks"), COMPLETION_CHECKS)
        or completion.get("pass") is not True
        or completion.get("pass_semantics")
        != (
            "completion_payload_only; parent promote exit-attestation publication and "
            "detached completion-key signature remain mandatory, followed by continuation "
            "verification, runner-gate evidence, and the signed downstream terminal receipt"
        )
        or not strict_json_equal(
            governance.get("continuation_verifier"), authorized_verifier_record
        )
        or not strict_json_equal(
            governance.get("runner_gate_replay"), authorized_replay_record
        )
        or not strict_json_equal(
            governance.get("downstream_continuation"), authorized_self_record
        )
        or not strict_json_equal(
            governance.get("downstream_continuation_command"),
            authorized_canonical_command(),
        )
        or not strict_json_equal(
            governance.get("supervised_child_command"),
            authorized_supervised_child_command(),
        )
        or governance.get("continuation_verification_receipt")
        != str(AUTHORIZED_PROMOTION_VERIFICATION_RECEIPT)
        or governance.get("runner_gate_receipt")
        != str(AUTHORIZED_RUNNER_GATE_RECEIPT)
        or not strict_json_equal(
            governance.get("downstream_output_slots_absent_at_promotion"),
            [str(path) for path in AUTHORIZED_DOWNSTREAM_OUTPUT_SLOTS],
        )
        or governance.get("canonical_publication_order")
        != "provisional_completion_receipt_then_canonical_no_replace_commit"
        or governance.get("canonical_alone_grants_downstream_authority") is not False
        or governance.get("write_ahead_receipt_marks_transaction_incomplete_until")
        != "detached_completion_signature_plus_live_canonical_verification"
        or governance.get("narrow_supersession")
        != (
            "Only the canonical-emitter requirement is superseded; all content checks, "
            "source authorities, downstream gate replay, and release signatures remain."
        )
    ):
        raise ContinuationError("Signed promotion completion exact contract drift")
    promote_attestation_payload = strict_json(
        promote_attestation_data, "signed promote producer exit attestation"
    )
    promote_attestation_payload, promote_attestation_time = (
        validate_producer_exit_attestation(
            promote_attestation_payload,
            label="promote producer exit attestation",
            expected_schema=PROMOTE_EXIT_ATTESTATION_SCHEMA,
            expected_phase="promote",
            expected_command=promotion_command("--promote"),
            expected_session=promote_session,
            expected_artifacts={
                "completion_receipt": completion_record,
                "canonical_receipt": canonical_record,
                "prepare_exit_attestation": prepare_attestation_record,
                "prepare_exit_attestation_signature": auth_sig_record,
            },
            expected_source=promotion_full_record,
            expected_runtime=python_full_record,
            expected_signature_contract=expected_completion_contract,
            expected_checks=PROMOTE_EXIT_CHECKS,
            expected_pass_semantics=PROMOTE_EXIT_PASS_SEMANTICS,
        )
    )
    if (
        promote_attestation_time < completion_time
        or promote_attestation_time < prepare_attestation_time
    ):
        raise ContinuationError("Promote producer exit attestation timestamp order drift")

    historical_modes = {
        (
            (
                "authorization_signature_contract",
                "private_key_lifecycle",
                "pre_signature",
            ),
            str(AUTHORIZATION_PRIVATE_KEY),
            "0o400",
        ),
        (
            (
                "completion_signature_contract",
                "private_key_lifecycle",
                "pre_signature",
            ),
            str(COMPLETION_PRIVATE_KEY),
            "0o400",
        ),
        (
            ("continuation_gate", "terminal_signature_contract", "private_key"),
            str(CONTINUATION_PRIVATE_KEY),
            "0o400",
        ),
        (
            ("signature_contract", "private_key_lifecycle", "pre_signature"),
            str(COMPLETION_PRIVATE_KEY),
            "0o400",
        ),
    }
    focal_during = (
        authorization.get("evidence", {})
        .get("focal_source_identity_transition", {})
        .get("during_execution")
    )
    if isinstance(focal_during, Mapping):
        historical_modes.add(
            (
                (
                    "evidence",
                    "focal_source_identity_transition",
                    "during_execution",
                ),
                str(focal_during.get("path")),
                str(focal_during.get("mode")),
            )
        )
    bound_relations = bind_declared_relations(
        authorization,
        reader,
        historical_mode_records=historical_modes,
    )
    bind_declared_relations(
        completion,
        reader,
        historical_mode_records=historical_modes,
        relations=bound_relations,
    )
    current_v7 = authorization["evidence"]["current_v7_runtime_validation"]
    if (
        current_v7.get("authority_id") != EXPECTED_V7_AUTHORITY_ID
        or current_v7.get("source_set_sha256") != EXPECTED_V7_SOURCE_SET_SHA256
        or current_v7.get("pass") is not True
    ):
        raise ContinuationError("Signed promotion v7 evidence drift")
    return {
        "self_fd": self_fd,
        "self_record": self_record,
        "authorized_self_record": authorized_self_record,
        "self_binding": self_binding,
        "promotion_record": promotion_record,
        "promotion_full_record": promotion_full_record,
        "verifier_fd": verifier_fd,
        "verifier_data": verifier_data,
        "verifier_record": verifier_record,
        "authorized_verifier_record": authorized_verifier_record,
        "replay_fd": replay_fd,
        "replay_record": replay_record,
        "authorized_replay_record": authorized_replay_record,
        "runner_fd": runner_fd,
        "runner_data": runner_data,
        "runner_record": runner_record,
        "authority_record": authority_record,
        "schema_recovery_validator_record": recovery_validator_record,
        "schema_recovery_authority": recovery_authority,
        "schema_recovery_module": recovery_module,
        "openssl_fd": openssl_fd,
        "openssl_record": openssl_record,
        "python_fd": python_fd,
        "python_record": python_record,
        "python_full_record": python_full_record,
        "qa_python_fd": reviewed_runtime_fds["qa_python_runtime"],
        "qa_python_record": reviewed_runtime_sources["qa_python_runtime"],
        "reviewed_runtime_sources": reviewed_runtime_sources,
        "reviewed_runtime_fds": reviewed_runtime_fds,
        "reviewed_sources": reviewed_sources,
        "authorization": authorization,
        "authorization_record": auth_record,
        "prepare_exit_attestation_record": prepare_attestation_record,
        "prepare_exit_attestation": prepare_attestation_payload,
        "authorization_signature_record": auth_sig_record,
        "authorization_public_key_record": auth_public_record,
        "authorization_private_key_record": auth_private_record,
        "completion": completion,
        "completion_record": completion_record,
        "promote_exit_attestation_record": promote_attestation_record,
        "promote_exit_attestation": promote_attestation_payload,
        "completion_signature_record": completion_sig_record,
        "completion_public_key_record": completion_public_record,
        "completion_private_key_record": completion_private_record,
        "continuation_public_key_record": continuation_public_record,
        "continuation_public_key_stat": continuation_public_stat,
        "continuation_private_key_record": continuation_private_record,
        "continuation_private_key_live_record": continuation_private_live_record,
        "continuation_private_key_stat": continuation_private_stat,
        "failed_v16_continuation_public_key_record": (
            failed_v16_continuation_public_record
        ),
        "failed_v16_continuation_public_key_stat": (
            failed_v16_continuation_public_stat
        ),
        "failed_v16_continuation_private_key_live_record": (
            failed_v16_continuation_private_live_record
        ),
        "failed_v16_continuation_private_key_stat": (
            failed_v16_continuation_private_stat
        ),
        "failed_v17_continuation_public_key_record": (
            failed_v17_continuation_public_record
        ),
        "failed_v17_continuation_public_key_stat": (
            failed_v17_continuation_public_stat
        ),
        "failed_v17_continuation_private_key_live_record": (
            failed_v17_continuation_private_live_record
        ),
        "failed_v17_continuation_private_key_stat": (
            failed_v17_continuation_private_stat
        ),
        "failed_v18_continuation_public_key_record": (
            failed_v18_continuation_public_record
        ),
        "failed_v18_continuation_public_key_stat": (
            failed_v18_continuation_public_stat
        ),
        "failed_v18_continuation_private_key_live_record": (
            failed_v18_continuation_private_live_record
        ),
        "failed_v18_continuation_private_key_stat": (
            failed_v18_continuation_private_stat
        ),
        "failed_v19_continuation_public_key_record": (
            failed_v19_continuation_public_record
        ),
        "failed_v19_continuation_public_key_stat": (
            failed_v19_continuation_public_stat
        ),
        "failed_v19_continuation_private_key_live_record": (
            failed_v19_continuation_private_live_record
        ),
        "failed_v19_continuation_private_key_stat": (
            failed_v19_continuation_private_stat
        ),
        "rejected_v20_continuation_public_key_record": (
            rejected_v20_continuation_public_record
        ),
        "rejected_v20_continuation_public_key_stat": (
            rejected_v20_continuation_public_stat
        ),
        "rejected_v20_continuation_private_key_live_record": (
            rejected_v20_continuation_private_live_record
        ),
        "rejected_v20_continuation_private_key_stat": (
            rejected_v20_continuation_private_stat
        ),
        "failed_v21_continuation_public_key_record": (
            failed_v21_continuation_public_record
        ),
        "failed_v21_continuation_public_key_stat": (
            failed_v21_continuation_public_stat
        ),
        "failed_v21_continuation_private_key_live_record": (
            failed_v21_continuation_private_live_record
        ),
        "failed_v21_continuation_private_key_stat": (
            failed_v21_continuation_private_stat
        ),
        "rejected_v22_continuation_public_key_record": (
            rejected_v22_continuation_public_record
        ),
        "rejected_v22_continuation_public_key_stat": (
            rejected_v22_continuation_public_stat
        ),
        "rejected_v22_continuation_private_key_live_record": (
            rejected_v22_continuation_private_live_record
        ),
        "rejected_v22_continuation_private_key_stat": (
            rejected_v22_continuation_private_stat
        ),
        "rejected_v23_continuation_public_key_record": (
            rejected_v23_continuation_public_record
        ),
        "rejected_v23_continuation_public_key_stat": (
            rejected_v23_continuation_public_stat
        ),
        "rejected_v23_continuation_private_key_live_record": (
            rejected_v23_continuation_private_live_record
        ),
        "rejected_v23_continuation_private_key_stat": (
            rejected_v23_continuation_private_stat
        ),
        "rejected_v24_continuation_public_key_record": (
            rejected_v24_continuation_public_record
        ),
        "rejected_v24_continuation_public_key_stat": (
            rejected_v24_continuation_public_stat
        ),
        "rejected_v24_continuation_private_key_live_record": (
            rejected_v24_continuation_private_live_record
        ),
        "rejected_v24_continuation_private_key_stat": (
            rejected_v24_continuation_private_stat
        ),
        "failed_v25_continuation_public_key_record": (
            failed_v25_continuation_public_record
        ),
        "failed_v25_continuation_public_key_stat": (
            failed_v25_continuation_public_stat
        ),
        "failed_v25_continuation_private_key_live_record": (
            failed_v25_continuation_private_live_record
        ),
        "failed_v25_continuation_private_key_stat": (
            failed_v25_continuation_private_stat
        ),
        "failed_v26_continuation_public_key_record": (
            failed_v26_continuation_public_record
        ),
        "failed_v26_continuation_public_key_stat": (
            failed_v26_continuation_public_stat
        ),
        "failed_v26_continuation_private_key_live_record": (
            failed_v26_continuation_private_live_record
        ),
        "failed_v26_continuation_private_key_stat": (
            failed_v26_continuation_private_stat
        ),
        "rejected_v27_continuation_public_key_record": (
            rejected_v27_continuation_public_record
        ),
        "rejected_v27_continuation_public_key_stat": (
            rejected_v27_continuation_public_stat
        ),
        "rejected_v27_continuation_private_key_live_record": (
            rejected_v27_continuation_private_live_record
        ),
        "rejected_v27_continuation_private_key_stat": (
            rejected_v27_continuation_private_stat
        ),
        "failed_v28_continuation_public_key_record": (
            failed_v28_continuation_public_record
        ),
        "failed_v28_continuation_public_key_stat": (
            failed_v28_continuation_public_stat
        ),
        "failed_v28_continuation_private_key_live_record": (
            failed_v28_continuation_private_live_record
        ),
        "failed_v28_continuation_private_key_stat": (
            failed_v28_continuation_private_stat
        ),
        "failed_v29_continuation_public_key_record": (
            failed_v29_continuation_public_record
        ),
        "failed_v29_continuation_public_key_stat": (
            failed_v29_continuation_public_stat
        ),
        "failed_v29_continuation_private_key_live_record": (
            failed_v29_continuation_private_live_record
        ),
        "failed_v29_continuation_private_key_stat": (
            failed_v29_continuation_private_stat
        ),
        "recovery_continuation_public_key_record": (
            recovery_continuation_public_record
        ),
        "recovery_continuation_public_key_stat": recovery_continuation_public_stat,
        "recovery_continuation_private_key_record": (
            recovery_continuation_private_record
        ),
        "recovery_continuation_private_key_live_record": (
            recovery_continuation_private_live_record
        ),
        "recovery_continuation_private_key_stat": recovery_continuation_private_stat,
        "source_record": source_record,
        "canonical_record": canonical_record,
        "trusted_keys": trusted_keys,
        "v7_validation": dict(v7_validation),
        "bound_relation_count": len(bound_relations),
    }


def validate_verification_payload(
    payload: Mapping[str, Any], *, expected_mode: str, chain: Mapping[str, Any]
) -> datetime:
    require_exact_keys(payload, VERIFICATION_KEYS, "promotion verification payload")
    created = require_utc_timestamp(
        payload.get("created_at_utc"), f"promotion verification {expected_mode} timestamp"
    )
    checks = payload.get("checks")
    governance = payload.get("governance")
    incident_disclosure = {
        "incident_id": "20260722_completion_attempt2_mode_tightening_conflict",
        "archive_receipt": chain["authorization"]["evidence"][
            "failed_attempt_archive"
        ]["receipt"],
        "failure_class": (
            "post_execution_source_mode_tightening_conflicts_with_full_stat_identity_equality"
        ),
        "failed_before_downstream_output": True,
        "failed_evidence_preserved": True,
        "original_completion_runner_end_to_end_pass_claimed": False,
        "original_verifier_reexecuted_successfully_claimed": False,
        "scientific_payload_changed": False,
    }
    def verification_artifact(record: Mapping[str, Any]) -> dict[str, Any]:
        return {key: record[key] for key in BASIC_ARTIFACT_KEYS}

    expected_authorization = {
        "artifact": verification_artifact(chain["authorization_record"]),
        "signature": verification_artifact(chain["authorization_signature_record"]),
        "signature_verified": True,
    }
    expected_completion = {
        "artifact": verification_artifact(chain["completion_record"]),
        "signature": verification_artifact(chain["completion_signature_record"]),
        "signature_verified": True,
    }
    expected_exit_attestations = {
        "prepare": {
            "artifact": chain["prepare_exit_attestation_record"],
            "signature": verification_artifact(
                chain["authorization_signature_record"]
            ),
            "public_key": chain["authorization_public_key_record"],
            "signature_verified": True,
        },
        "promote": {
            "artifact": chain["promote_exit_attestation_record"],
            "signature": verification_artifact(chain["completion_signature_record"]),
            "public_key": chain["completion_public_key_record"],
            "signature_verified": True,
        },
    }
    expected_retirement = {
        "authorization": chain["authorization_private_key_record"],
        "completion": chain["completion_private_key_record"],
    }
    recovery_authority = payload.get("schema_recovery_authority")
    historical_identity = payload.get("historical_identity_contract")
    chain_recovery_authority = chain["schema_recovery_authority"]
    retained_count = payload.get("retained_critical_input_descriptor_count")
    if (
        payload.get("schema_name") != VERIFICATION_SCHEMA
        or payload.get("schema_version") != RECOVERY_SCHEMA_VERSION
        or payload.get("task_type") != "B_comprehensive_validation"
        or payload.get("mode") != expected_mode
        or not strict_json_equal(
            payload.get("authorization"), expected_authorization
        )
        or not strict_json_equal(payload.get("completion"), expected_completion)
        or not strict_json_equal(
            payload.get("producer_exit_attestations"), expected_exit_attestations
        )
        or not strict_json_equal(payload.get("source_receipt"), chain["source_record"])
        or not strict_json_equal(
            payload.get("canonical_receipt"), chain["canonical_record"]
        )
        or not strict_json_equal(
            payload.get("trusted_signing_keys"), chain["trusted_keys"]
        )
        or not strict_json_equal(
            payload.get("private_key_retirement"), expected_retirement
        )
        or not strict_json_equal(
            payload.get("historical_incident_disclosure"), incident_disclosure
        )
        or not strict_json_equal(
            payload.get("release_source_authority_validator"),
            chain["authority_record"],
        )
        or not isinstance(recovery_authority, Mapping)
        or recovery_authority.get("authority_id")
        != "20260724_tumor_ref_schema_recovery_v30"
        or recovery_authority.get("runtime_role") != "continuation_verifier"
        or not strict_json_equal(
            recovery_authority.get("runtime_source"), chain["verifier_record"]
        )
        or recovery_authority.get("authority")
        != chain_recovery_authority.get("authority")
        or recovery_authority.get("signature")
        != chain_recovery_authority.get("signature")
        or recovery_authority.get("public_key")
        != chain_recovery_authority.get("public_key")
        or recovery_authority.get("historical_link_count_observation")
        != "NOT_RECORDED_NOT_INFERRED"
        or recovery_authority.get("pass") is not True
        or not isinstance(historical_identity, Mapping)
        or set(historical_identity)
        != {
            "during_execution_keys",
            "current_identity_keys",
            "historical_link_count_observation",
            "current_link_count_requirement",
            "scientific_payload_changed",
        }
        or historical_identity.get("during_execution_keys")
        != ["ctime_ns", "device", "inode", "mode", "mtime_ns", "path", "sha256", "size_bytes"]
        or historical_identity.get("current_identity_keys")
        != [
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
        or historical_identity.get("historical_link_count_observation")
        != "NOT_RECORDED_NOT_INFERRED"
        or historical_identity.get("current_link_count_requirement") != 1
        or historical_identity.get("scientific_payload_changed") is not False
        or not strict_json_equal(checks, VERIFICATION_CHECKS)
        or not strict_json_equal(governance, VERIFICATION_GOVERNANCE)
        or type(retained_count) is not int
        or retained_count <= 0
        or payload.get("pass") is not True
        or payload.get("pass_semantics") != VERIFICATION_PASS_SEMANTICS
    ):
        raise ContinuationError("Promotion verification payload contract drift")
    return created


def run_fresh_verifier(
    reader: Any,
    chain: Mapping[str, Any],
) -> dict[str, Any]:
    _, recorded_data, recorded_record, _ = open_bound(
        reader, PROMOTION_VERIFICATION_RECEIPT
    )
    recorded = strict_json(recorded_data, "recorded promotion verification")
    recorded_time = validate_verification_payload(
        recorded, expected_mode="verify-and-record", chain=chain
    )
    completion_time = require_utc_timestamp(
        chain["completion"].get("created_at_utc"), "signed completion timestamp"
    )
    if recorded_time < completion_time:
        raise ContinuationError("Recorded promotion verification predates completion")
    reader.require_paths_still_bound()
    actual_argv = [
        f"/proc/self/fd/{chain['python_fd']}",
        "-I",
        "-B",
        f"/proc/self/fd/{chain['verifier_fd']}",
        "--verify",
    ]
    result = subprocess.run(
        actual_argv,
        cwd=REPO_ROOT,
        env=EXPECTED_ENVIRONMENT,
        pass_fds=(int(chain["python_fd"]), int(chain["verifier_fd"])),
        capture_output=True,
        check=False,
        timeout=600,
    )
    reader.require_paths_still_bound()
    require_incident_absent()
    if result.returncode != 0 or result.stderr:
        raise ContinuationError(
            "Fresh promotion verifier failed: "
            f"exit={result.returncode} stderr_sha256={sha256_bytes(result.stderr)}"
        )
    fresh = strict_json(result.stdout, "fresh promotion verifier stdout")
    fresh_time = validate_verification_payload(
        fresh, expected_mode="verify", chain=chain
    )
    comparable = set(VERIFICATION_KEYS) - {
        "created_at_utc",
        "mode",
        "retained_critical_input_descriptor_count",
    }
    if any(
        not strict_json_equal(fresh.get(key), recorded.get(key)) for key in comparable
    ):
        raise ContinuationError("Fresh and recorded promotion verification evidence differ")
    if fresh_time < recorded_time:
        raise ContinuationError("Fresh promotion verification predates recorded verification")
    historical_modes: set[tuple[tuple[str, ...], str, str]] = set()
    relations = bind_declared_relations(
        recorded, reader, historical_mode_records=historical_modes
    )
    bind_declared_relations(
        fresh,
        reader,
        historical_mode_records=historical_modes,
        relations=relations,
    )
    return {
        "recorded_payload": recorded,
        "recorded_record": recorded_record,
        "recorded_time": recorded_time,
        "fresh_payload": fresh,
        "fresh_time": fresh_time,
        "fresh_stdout": {
            "size_bytes": len(result.stdout),
            "sha256": sha256_bytes(result.stdout),
        },
        "fresh_stderr": {
            "size_bytes": len(result.stderr),
            "sha256": sha256_bytes(result.stderr),
        },
        "fresh_actual_argv": actual_argv,
        "fresh_normalized_argv": [
            "/proc/self/fd/<bound-python-fd>",
            "-I",
            "-B",
            "/proc/self/fd/<bound-verifier-fd>",
            "--verify",
        ],
        "fresh_pass_fds": {
            str(chain["python_fd"]): chain["python_record"],
            str(chain["verifier_fd"]): chain["verifier_record"],
        },
        "bound_relation_count": len(relations),
    }


def load_recorded_verification(
    reader: Any,
    chain: Mapping[str, Any],
) -> dict[str, Any]:
    _, recorded_data, recorded_record, _ = open_bound(
        reader, PROMOTION_VERIFICATION_RECEIPT
    )
    recorded = strict_json(recorded_data, "recorded promotion verification")
    recorded_time = validate_verification_payload(
        recorded, expected_mode="verify-and-record", chain=chain
    )
    completion_time = require_utc_timestamp(
        chain["completion"].get("created_at_utc"), "signed completion timestamp"
    )
    if recorded_time < completion_time:
        raise ContinuationError("Recorded promotion verification predates completion")
    historical_modes: set[tuple[tuple[str, ...], str, str]] = set()
    relations = bind_declared_relations(
        recorded, reader, historical_mode_records=historical_modes
    )
    return {
        "recorded_payload": recorded,
        "recorded_record": recorded_record,
        "recorded_time": recorded_time,
        "bound_relation_count": len(relations),
    }


def require_replay_python_alias_consistency(
    authoritative_target: Mapping[str, Any],
    runtime_alias: Any,
    fd_binding_alias: Any,
    completion_contract_alias: Any,
) -> None:
    require_exact_keys(
        authoritative_target,
        BASIC_ARTIFACT_KEYS,
        "authoritative primary Python target",
    )
    expected_alias = {
        "path": str(PYTHON),
        "resolved_path": str(PRIMARY_PYTHON_RUNTIME),
        "symlink_target": PRIMARY_PYTHON_RUNTIME.name,
        "size_bytes": authoritative_target["size_bytes"],
        "sha256": authoritative_target["sha256"],
        "mode": authoritative_target["mode"],
    }
    if authoritative_target["path"] != str(PRIMARY_PYTHON_RUNTIME):
        raise ContinuationError("Authoritative primary Python target path drift")
    for label, observed in (
        ("runtime", runtime_alias),
        ("promotion FD binding", fd_binding_alias),
        ("completion runner contract", completion_contract_alias),
    ):
        if (
            not isinstance(observed, Mapping)
            or set(observed) != set(EXECUTABLE_ALIAS_RECORD_KEYS)
            or not strict_json_equal(observed, expected_alias)
        ):
            raise ContinuationError(f"Runner replay {label} Python alias drift")


def validate_replay_success_witness(
    witness: Any,
    *,
    chain: Mapping[str, Any],
    receipt_record: Mapping[str, Any],
    log_record: Mapping[str, Any],
    replay_time: datetime,
) -> dict[str, Any]:
    require_exact_keys(
        witness,
        REPLAY_SUCCESS_WITNESS_KEYS,
        "runner replay success witness",
    )
    witness_time = require_utc_timestamp(
        witness.get("created_at_utc"), "runner replay success witness timestamp"
    )
    worker_wait = witness.get("worker_wait")
    worker_outputs = witness.get("worker_outputs")
    if (
        witness.get("schema_name") != RUNNER_GATE_SUCCESS_WITNESS_SCHEMA
        or witness.get("schema_version") != "1.0.0"
        or witness.get("task_type") != "B_comprehensive_validation"
        or witness.get("status") != "WORKER_WAITPID_NORMAL_EXIT_ZERO_ATTESTED"
        or witness.get("scope") != "completion_runner_physical_lines_1_358_only"
        or not strict_json_equal(witness.get("canonical_command"), replay_command())
        or not strict_json_equal(witness.get("clean_environment"), EXPECTED_ENVIRONMENT)
        or not strict_json_equal(
            witness.get("supervisor_source"), chain["replay_record"]
        )
        or not isinstance(worker_wait, Mapping)
        or set(worker_wait)
        != {"pid", "raw_wait_status", "normal_exit", "exit_code", "term_signal"}
        or type(worker_wait.get("pid")) is not int
        or worker_wait["pid"] <= 0
        or type(worker_wait.get("raw_wait_status")) is not int
        or worker_wait["raw_wait_status"] != 0
        or worker_wait.get("normal_exit") is not True
        or type(worker_wait.get("exit_code")) is not int
        or worker_wait["exit_code"] != 0
        or worker_wait.get("term_signal") is not None
        or not isinstance(worker_outputs, Mapping)
        or set(worker_outputs) != {"log", "receipt"}
        or not strict_json_equal(worker_outputs.get("log"), log_record)
        or not strict_json_equal(worker_outputs.get("receipt"), receipt_record)
        or not strict_json_equal(
            witness.get("checks"), REPLAY_SUCCESS_WITNESS_CHECKS
        )
        or witness.get("pass") is not True
        or witness.get("pass_semantics")
        != REPLAY_SUCCESS_WITNESS_PASS_SEMANTICS
        or witness_time < replay_time
    ):
        raise ContinuationError("Runner replay success witness contract drift")
    return {
        "timestamp": witness["created_at_utc"],
        "worker_wait": dict(worker_wait),
    }


def validate_runner_archive_rebase_evidence(
    reader: Any,
    replay: Mapping[str, Any],
    completion_runner: Mapping[str, Any],
    runner_only: Mapping[str, Any],
) -> dict[str, Any]:
    _, result_archive_data, result_archive_record, result_archive_stat = open_bound(
        reader, FAILED_V28_RESULT_KEY_ARCHIVE_RECORD
    )
    _, _, result_public_record, result_public_stat = open_bound(
        reader, FAILED_V28_RESULT_PUBLIC_KEY
    )
    _, report_archive_data, report_archive_record, report_archive_stat = open_bound(
        reader, FAILED_V28_REPORT_KEY_ARCHIVE_RECORD
    )
    _, _, report_public_record, report_public_stat = open_bound(
        reader, FAILED_V28_REPORT_PUBLIC_KEY
    )
    report_private_fd, report_private_record, report_private_stat = open_private_metadata(
        reader, FAILED_V28_REPORT_PRIVATE_KEY, "0o400"
    )
    del report_private_fd

    if (
        result_archive_record.get("sha256")
        != EXPECTED_FAILED_V28_RESULT_ARCHIVE_RECORD_SHA256
        or report_archive_record.get("sha256")
        != EXPECTED_FAILED_V28_REPORT_ARCHIVE_RECORD_SHA256
        or result_public_record.get("sha256")
        != EXPECTED_FAILED_V28_RESULT_PUBLIC_KEY_SHA256
        or report_public_record.get("sha256")
        != EXPECTED_FAILED_V28_REPORT_PUBLIC_KEY_SHA256
        or result_archive_stat.st_nlink != 1
        or report_archive_stat.st_nlink != 1
        or result_public_stat.st_nlink != 1
        or report_public_stat.st_nlink != 1
        or report_private_stat.st_nlink != 1
        or report_private_stat.st_size != 290
    ):
        raise ContinuationError("Failed-v28 archive-rebase artifact identity drift")

    result_archive = strict_json(result_archive_data, "failed-v28 result key archive")
    report_archive = strict_json(report_archive_data, "failed-v28 report key archive")
    expected_common = {
        "schema_name": "intersubmod.failed_formal_generation_key_archive",
        "schema_version": "1.0.0",
        "generation": "v28",
        "key_reuse_forbidden": True,
        "private_key_bytes_read": False,
        "pass": False,
    }
    for role, payload, archive_root, public_record, state in (
        (
            "result_v5",
            result_archive,
            FAILED_V28_RESULT_KEY_ARCHIVE,
            result_public_record,
            "CONSUMED_DATASET_SIGNATURE_PROVISIONAL_AFTER_C_FAILURE",
        ),
        (
            "report_v5",
            report_archive,
            FAILED_V28_REPORT_KEY_ARCHIVE,
            report_public_record,
            "UNUSED_REPORT_KEY_RETIRED_FROM_FAILED_GENERATION",
        ),
    ):
        if (
            any(payload.get(key) != value for key, value in expected_common.items())
            or payload.get("role") != role
            or payload.get("archive_directory") != str(archive_root)
            or payload.get("state") != state
            or not isinstance(payload.get("public_key"), Mapping)
            or payload["public_key"].get("sha256") != public_record.get("sha256")
            or payload["public_key"].get("size_bytes") != 113
            or payload["public_key"].get("mode") != "0o444"
        ):
            raise ContinuationError(f"Failed-v28 {role} archive contract drift")
    report_private_declared = report_archive.get("private_key_metadata_only")
    if (
        not isinstance(report_private_declared, Mapping)
        or report_private_declared.get("size_bytes") != 290
        or report_private_declared.get("mode") != "0o400"
        or report_private_declared.get("link_count") != 1
    ):
        raise ContinuationError("Failed-v28 report private-key metadata drift")

    expected_private_record = {
        **report_private_record,
        "size_bytes": report_private_stat.st_size,
    }
    expected_archive_inputs = {
        "result_key_archive_record": result_archive_record,
        "result_public_key": result_public_record,
        "report_key_archive_record": report_archive_record,
        "report_public_key": report_public_record,
        "report_private_key_metadata_only": expected_private_record,
        "report_private_key_link_count": report_private_stat.st_nlink,
        "private_key_bytes_read": False,
        "key_reuse_forbidden_preserved": True,
        "archive_paths_remain_canonical": True,
        "pass": True,
    }
    archive_inputs = replay.get("failed_v28_key_archive_rebase_inputs")
    if not strict_json_equal(archive_inputs, expected_archive_inputs):
        raise ContinuationError("Runner replay archive input contract drift")

    expected_mappings = [
        {
            "role": role,
            "original_path": original,
            "archive_path": archive,
            "replacement_count": 1,
        }
        for original, archive, role in RUNNER_ARCHIVE_REBASE
    ]
    expected_rebase_contract = {
        "contract": "deterministic_three_path_archive_rebase_v1",
        "original_prefix": {
            "size_bytes": EXPECTED_RUNNER_PREFIX_SIZE,
            "sha256": EXPECTED_RUNNER_PREFIX_SHA256,
        },
        "archive_rebased_prefix": {
            "size_bytes": EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SIZE,
            "sha256": EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SHA256,
        },
        "mappings": expected_mappings,
        "mapping_count": 3,
        "private_key_bytes_read": False,
        "host_live_key_paths_created": False,
        "pass": True,
    }
    if not strict_json_equal(
        replay.get("runner_archive_rebase_contract"), expected_rebase_contract
    ):
        raise ContinuationError("Runner archive-rebase transform contract drift")

    expected_live_absence = {str(path): True for path in FAILED_V28_LIVE_KEY_PATHS}
    expected_live_checks = {
        "before_verifier": expected_live_absence,
        "after_verifier": expected_live_absence,
        "after_runner_replay": expected_live_absence,
        "all_absent": True,
    }
    if (
        any(os.path.lexists(path) for path in FAILED_V28_LIVE_KEY_PATHS)
        or not strict_json_equal(
            replay.get("failed_v28_live_key_path_checks"), expected_live_checks
        )
    ):
        raise ContinuationError("Failed-v28 live key path absence contract drift")

    expected_stdin = {
        "source": str(COMPLETION_RUNNER),
        "physical_lines": "1-358",
        "original_size_bytes": EXPECTED_RUNNER_PREFIX_SIZE,
        "original_sha256": EXPECTED_RUNNER_PREFIX_SHA256,
        "archive_rebased_size_bytes": EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SIZE,
        "archive_rebased_sha256": EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SHA256,
        "transform_contract": "deterministic_three_path_archive_rebase_v1",
    }
    if (
        completion_runner.get("stdin_size_bytes") != EXPECTED_RUNNER_PREFIX_SIZE
        or completion_runner.get("stdin_sha256") != EXPECTED_RUNNER_PREFIX_SHA256
        or completion_runner.get("executed_stdin_size_bytes")
        != EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SIZE
        or completion_runner.get("executed_stdin_sha256")
        != EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SHA256
        or not strict_json_equal(runner_only.get("stdin"), expected_stdin)
    ):
        raise ContinuationError("Runner archive-rebase dual stdin identity drift")
    return {
        "archive_inputs": expected_archive_inputs,
        "transform": expected_rebase_contract,
        "live_path_checks": expected_live_checks,
        "private_key_bytes_read": False,
    }


def validate_replay_evidence(
    reader: Any,
    chain: Mapping[str, Any],
    verification: Mapping[str, Any],
) -> dict[str, Any]:
    _, receipt_data, receipt_record, _ = open_bound(reader, RUNNER_GATE_RECEIPT)
    _, log_data, log_record, _ = open_bound(reader, RUNNER_GATE_LOG)
    _, witness_data, witness_record, _ = open_bound(
        reader, RUNNER_GATE_SUCCESS_WITNESS
    )
    replay = strict_json(receipt_data, "runner-only replay receipt")
    witness = strict_json(witness_data, "runner replay success witness")
    require_exact_keys(replay, REPLAY_KEYS, "runner-only replay receipt")
    replay_time = require_utc_timestamp(
        replay.get("created_at_utc"), "runner-only replay timestamp"
    )
    witness_validation = validate_replay_success_witness(
        witness,
        chain=chain,
        receipt_record=receipt_record,
        log_record=log_record,
        replay_time=replay_time,
    )
    completion_runner = replay.get("completion_runner")
    runner_only = replay.get("runner_only_replay")
    claims = replay.get("claims")
    code = replay.get("code")
    legacy_signed_code = replay.get("legacy_signed_code")
    recovery_authority = replay.get("schema_recovery_authority")
    historical_identity = replay.get("historical_identity_contract")
    promotion_verifier = replay.get("promotion_verifier")
    trust_chain = replay.get("promotion_trust_chain")
    replay_v7 = replay.get("v7_validation")
    output_checks = replay.get("output_slot_checks")
    runtime = replay.get("runtime")
    require_exact_keys(runtime, REPLAY_RUNTIME_KEYS, "runner-only replay runtime")
    require_exact_keys(
        promotion_verifier,
        PROMOTION_VERIFIER_KEYS,
        "runner-only replay promotion verifier",
    )
    require_exact_keys(
        trust_chain,
        PROMOTION_TRUST_CHAIN_KEYS,
        "runner-only replay promotion trust chain",
    )
    validate_replay_source_identity_contract(
        trust_chain,
        chain["source_record"],
        chain["canonical_record"],
    )
    require_exact_keys(
        completion_runner,
        REPLAY_COMPLETION_RUNNER_KEYS,
        "runner-only replay completion runner",
    )
    require_exact_keys(
        runner_only, REPLAY_RUNNER_ONLY_KEYS, "runner-only replay execution"
    )
    require_exact_keys(
        output_checks, REPLAY_OUTPUT_SLOT_KEYS, "runner-only replay output checks"
    )
    require_exact_keys(claims, REPLAY_CLAIM_KEYS, "runner-only replay claims")
    archive_rebase = validate_runner_archive_rebase_evidence(
        reader, replay, completion_runner, runner_only
    )
    expected_shell_observation_checks = {
        "executing_replayer_definitions_match_bound_reviewed_source": True,
        "exact_original_runner_prefix_bytes_lines_1_358_bound": True,
        "deterministic_three_path_archive_rebase_exact": True,
        "failed_v28_key_archive_records_and_modes_bound": True,
        "private_key_bytes_not_read": True,
        "host_live_key_paths_not_created": True,
        "shell_pathname_checks_exited_zero_non_authoritative": True,
        "parent_did_not_treat_shell_result_as_downstream_authority": True,
        "signed_downstream_continuation_must_revalidate_bound_inputs": True,
        "runner_prefix_exit_zero_under_set_Eeuo_pipefail": True,
    }
    fd_bindings = promotion_verifier.get("fd_artifact_binding", {})
    if not isinstance(fd_bindings, Mapping) or set(fd_bindings) != {
        "python",
        "verifier",
    }:
        raise ContinuationError("Runner replay FD binding schema drift")
    python_fd_binding = fd_bindings.get("python", {})
    verifier_fd_binding = fd_bindings.get("verifier", {})
    runner_runtime_contract = completion_runner.get(
        "authorization_bound_runtime_contract"
    )
    runner_runtime_artifacts = (
        runner_runtime_contract.get("runtime_artifacts")
        if isinstance(runner_runtime_contract, Mapping)
        else None
    )
    runtime_python = runtime.get("python")
    if (
        not isinstance(python_fd_binding, Mapping)
        or set(python_fd_binding) != {"fd", "artifact"}
        or not isinstance(verifier_fd_binding, Mapping)
        or set(verifier_fd_binding) != {"fd", "artifact"}
    ):
        raise ContinuationError("Runner replay nested FD binding schema drift")
    if not isinstance(runner_runtime_artifacts, Mapping):
        raise ContinuationError("Runner replay runtime artifact contract drift")
    require_replay_python_alias_consistency(
        chain["python_record"],
        runtime_python,
        python_fd_binding.get("artifact"),
        runner_runtime_artifacts.get("python"),
    )
    python_fd_value = python_fd_binding.get("fd")
    verifier_fd_value = verifier_fd_binding.get("fd")
    expected_absence = {str(path): True for path in DOWNSTREAM_OUTPUT_SLOTS}
    expected_trust_chain = {
        "authorization": chain["authorization_record"],
        "prepare_exit_attestation": chain["prepare_exit_attestation_record"],
        "authorization_signature": chain["authorization_signature_record"],
        "authorization_signature_target": str(PREPARE_EXIT_ATTESTATION),
        "authorization_public_key": chain["authorization_public_key_record"],
        "retired_authorization_private_key": chain[
            "authorization_private_key_record"
        ],
        "promotion_receipt": chain["completion_record"],
        "promote_exit_attestation": chain["promote_exit_attestation_record"],
        "promotion_signature": chain["completion_signature_record"],
        "promotion_signature_target": str(PROMOTE_EXIT_ATTESTATION),
        "completion_public_key": chain["completion_public_key_record"],
        "retired_completion_private_key": chain["completion_private_key_record"],
        "continuation_public_key": chain["continuation_public_key_record"],
        "pre_signature_continuation_private_key": chain[
            "continuation_private_key_record"
        ],
        "historical_source_receipt": chain["source_record"],
        "canonical_source_receipt": chain["canonical_record"],
        "reviewed_sources": chain["reviewed_sources"],
        "portable_reader_asset_membership": [
            str(path) for path in PORTABLE_READER_ASSET_PARTS
        ],
        "authorization_signature_verified": True,
        "completion_signature_verified": True,
        "prepare_child_waitpid_normal_exit_zero_verified": True,
        "promote_child_waitpid_normal_exit_zero_verified": True,
        "legacy_payload_signatures_absent": True,
        "promote_attestation_binds_signed_prepare_chain": True,
        "canonical_bytes_equal_historical_source": True,
        "canonical_inode_distinct_from_historical_source": True,
    }
    if (
        replay.get("schema_name") != REPLAY_SCHEMA
        or replay.get("schema_version") != "1.0.0"
        or replay.get("task_type") != "B_comprehensive_validation"
        or replay.get("status")
        != "ARCHIVE_REBASED_RUNNER_PREFIX_REPLAY_EVIDENCE_PASS_NON_AUTHORITATIVE"
        or replay.get("scope") != "completion_runner_physical_lines_1_358_only"
        or not strict_json_equal(replay.get("canonical_command"), replay_command())
        or not strict_json_equal(replay.get("clean_environment"), EXPECTED_ENVIRONMENT)
        or not strict_json_equal(replay.get("log"), log_record)
        or replay.get("downstream_not_executed") is not True
        or replay.get("downstream_authorized_after_this_gate") is not False
        or replay.get("pass") is not True
        or replay.get("pass_semantics") != REPLAY_PASS_SEMANTICS
        or not isinstance(code, Mapping)
        or set(code)
        != {
            "runner_gate_replay",
            "continuation_verifier",
            "downstream_continuation",
            "release_source_authority_validator",
            "schema_recovery_authority_validator",
        }
        or not strict_json_equal(
            legacy_signed_code,
            {
                "runner_gate_replay": chain["authorized_replay_record"],
                "continuation_verifier": chain["authorized_verifier_record"],
                "downstream_continuation": chain["authorized_self_record"],
            },
        )
        or not isinstance(recovery_authority, Mapping)
        or recovery_authority.get("authority_id")
        != "20260724_tumor_ref_schema_recovery_v30"
        or recovery_authority.get("runtime_role") != "runner_gate_replay"
        or not strict_json_equal(
            recovery_authority.get("runtime_source"), chain["replay_record"]
        )
        or recovery_authority.get("authority")
        != chain["schema_recovery_authority"].get("authority")
        or recovery_authority.get("signature")
        != chain["schema_recovery_authority"].get("signature")
        or recovery_authority.get("pass") is not True
        or not strict_json_equal(
            historical_identity,
            {
                "historical_link_count_observation": "NOT_RECORDED_NOT_INFERRED",
                "scientific_payload_changed": False,
                "claim_ceiling_changed": False,
            },
        )
        or not strict_json_equal(
            promotion_verifier.get("verification_receipt"),
            verification["recorded_record"],
        )
        or promotion_verifier.get("fresh_and_recorded_evidence_agree") is not True
        or not strict_json_equal(
            promotion_verifier.get("legacy_signed_command"),
            authorized_verifier_command(),
        )
        or not strict_json_equal(
            promotion_verifier.get("recovery_authorized_command"),
            verifier_command(),
        )
        or type(python_fd_value) is not int
        or python_fd_value <= 2
        or type(verifier_fd_value) is not int
        or verifier_fd_value <= 2
        or python_fd_value == verifier_fd_value
        or not strict_json_equal(
            promotion_verifier.get("pass_fds"),
            [python_fd_value, verifier_fd_value],
        )
        or not strict_json_equal(
            promotion_verifier.get("actual_argv"),
            [
                f"/proc/self/fd/{python_fd_value}",
                "-I",
                "-B",
                f"/proc/self/fd/{verifier_fd_value}",
                "--verify",
            ],
        )
        or not strict_json_equal(
            verifier_fd_binding.get("artifact"), chain["verifier_record"]
        )
        or not strict_json_equal(
            promotion_verifier.get("normalized_argv"),
            [
                "/proc/self/fd/<bound-python-fd>",
                "-I",
                "-B",
                "/proc/self/fd/<bound-verifier-fd>",
                "--verify",
            ],
        )
        or type(promotion_verifier.get("exit_code")) is not int
        or promotion_verifier.get("exit_code") != 0
        or promotion_verifier.get("source_descriptor_bound") is not True
        or promotion_verifier.get("pass") is not True
        or not strict_json_equal(trust_chain, expected_trust_chain)
        or not isinstance(replay_v7, Mapping)
        or replay_v7.get("authority_id") != EXPECTED_V7_AUTHORITY_ID
        or replay_v7.get("source_set_sha256") != EXPECTED_V7_SOURCE_SET_SHA256
        or replay_v7.get("pass") is not True
        or output_checks.get("all_absent") is not True
        or not strict_json_equal(output_checks.get("before_verifier"), expected_absence)
        or not strict_json_equal(output_checks.get("after_verifier"), expected_absence)
        or not strict_json_equal(
            output_checks.get("after_runner_replay"), expected_absence
        )
        or not strict_json_equal(
            completion_runner.get("artifact"), chain["runner_record"]
        )
        or type(completion_runner.get("observed_total_physical_lines")) is not int
        or completion_runner["observed_total_physical_lines"]
        != EXPECTED_RUNNER_LINE_COUNT
        or type(completion_runner.get("included_line_start")) is not int
        or completion_runner["included_line_start"] != 1
        or type(completion_runner.get("included_line_end")) is not int
        or completion_runner["included_line_end"] != 358
        or not strict_json_equal(
            runner_only.get("command"), ["/proc/self/fd/<bound-bash-fd>", "-s"]
        )
        or type(runner_only.get("exit_code")) is not int
        or runner_only.get("exit_code") != 0
        or runner_only.get("pass") is not True
        or runner_only.get("authoritative_gate_validation") is not False
        or not strict_json_equal(
            runner_only.get("shell_observation_checks"),
            expected_shell_observation_checks,
        )
        or claims.get("original_completion_runner_end_to_end_pass_claimed") is not False
        or claims.get("original_tumor_ref_verifier_successfully_reexecuted_claimed")
        is not False
        or claims.get("promotion_integrity_verified") is not True
        or claims.get("archive_rebased_runner_lines_1_358_executed_successfully")
        is not True
        or claims.get("runner_prefix_archive_rebased_three_paths") is not True
        or claims.get("failed_v28_private_key_bytes_read") is not False
        or claims.get("failed_v28_live_key_paths_created") is not False
        or claims.get("shell_pathname_gate_result_is_authoritative") is not False
        or claims.get("downstream_authorized_by_this_receipt") is not False
        or claims.get("runner_line_359_or_later_executed") is not False
        or claims.get("strict_executed") is not False
        or claims.get("matched_normal_executed") is not False
        or claims.get("cn_ccf_executed") is not False
        or claims.get("final_dataset_or_report_executed") is not False
    ):
        raise ContinuationError("Runner-only replay is not exact non-authoritative evidence")
    if (
        not strict_json_equal(code.get("runner_gate_replay"), chain["replay_record"])
        or not strict_json_equal(
            code.get("continuation_verifier"), chain["verifier_record"]
        )
        or not strict_json_equal(
            code.get("downstream_continuation"), chain["self_record"]
        )
        or not strict_json_equal(
            code.get("release_source_authority_validator"), chain["authority_record"]
        )
        or not strict_json_equal(
            code.get("schema_recovery_authority_validator"),
            chain["schema_recovery_validator_record"],
        )
    ):
        raise ContinuationError("Runner-only replay code relations drift")
    if replay_time < verification["recorded_time"]:
        raise ContinuationError("Runner-only replay predates recorded promotion verification")
    historical_modes: set[tuple[tuple[str, ...], str, str]] = set()
    relations = bind_declared_relations(
        replay, reader, historical_mode_records=historical_modes
    )
    return {
        "receipt": receipt_record,
        "log": log_record,
        "success_witness": witness_record,
        "worker_wait": witness_validation["worker_wait"],
        "success_witness_timestamp": witness_validation["timestamp"],
        "timestamp": replay.get("created_at_utc"),
        "bound_relation_count": len(relations),
        "log_size_bytes": len(log_data),
        "archive_rebase": archive_rebase,
        "non_authoritative": True,
        "downstream_authorized_after_this_gate": False,
    }


def require_basic_declared_artifact(
    observed: Any, live: Mapping[str, Any], label: str
) -> None:
    expected = basic_projection(live)
    if not isinstance(observed, Mapping) or not strict_json_equal(
        dict(observed), expected
    ):
        raise ContinuationError(f"{label} declared identity drift")


def validate_mandatory_gate_executable_records(
    gates: StreamBoundArtifactSet,
) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for role, contract in MANDATORY_GATE_EXECUTABLE_CONTRACTS.items():
        path = contract["path"]
        if GATE_INPUT_PATHS.get(role) != path:
            raise ContinuationError(
                f"Mandatory gate executable path registry drift: {role}"
            )
        try:
            record = gates.record_for(path)
        except (FileNotFoundError, KeyError, OSError, ContinuationError) as error:
            raise ContinuationError(
                f"Mandatory gate executable is missing or unbound: {role}"
            ) from error
        if (
            set(record) != {"path", "size_bytes", "sha256", "mode"}
            or record.get("path") != str(path)
            or type(record.get("size_bytes")) is not int
            or record["size_bytes"] <= 0
            or record.get("mode") != contract["mode"]
            or record.get("sha256") != contract["sha256"]
        ):
            raise ContinuationError(
                f"Mandatory gate executable identity drift: {role}"
            )
        records[role] = record
    return records


def validate_reviewed_runtime_gate_registry() -> dict[str, Path]:
    expected = {
        role: contract["path"]
        for role, contract in REVIEWED_RUNTIME_SOURCE_CONTRACTS.items()
    }
    runtime_paths = set(expected.values())
    observed = {
        role: path for role, path in GATE_INPUT_PATHS.items() if path in runtime_paths
    }
    if observed != expected:
        raise ContinuationError("Reviewed runtime gate path registry drift")
    return expected


def validate_gate_inputs(
    gates: StreamBoundArtifactSet,
    reviewed_runtime_sources: Mapping[str, Any],
) -> dict[str, Any]:
    validate_reviewed_runtime_gate_registry()
    expected_modes = {
        role: str(contract["mode"])
        for role, contract in {
            **REVIEWED_RUNTIME_SOURCE_CONTRACTS,
            **MANDATORY_GATE_EXECUTABLE_CONTRACTS,
        }.items()
    }
    for role, path in GATE_INPUT_PATHS.items():
        if role == "v7_private_key":
            gates.open_metadata(path, "0o0")
        else:
            gates.open(
                path,
                expected_mode=expected_modes.get(role),
                expected_link_count=1,
            )
    tumor_ref_summary_alias = gates.open_alias(
        TUMOR_REF_SUMMARY_ALIAS,
        TUMOR_REF_SUMMARY,
        expected_link_text=TUMOR_REF_SUMMARY.name,
    )
    runtime_hashes = {
        contract["path"]: contract["sha256"]
        for contract in REVIEWED_RUNTIME_SOURCE_CONTRACTS.values()
    }
    for path, expected_sha256 in runtime_hashes.items():
        if gates.record_for(path)["sha256"] != expected_sha256:
            raise ContinuationError(f"Reviewed downstream runtime identity drift: {path}")
    gate_runtime_sources = {
        role: gates.record_for(contract["path"])
        for role, contract in REVIEWED_RUNTIME_SOURCE_CONTRACTS.items()
    }
    if not strict_json_equal(gate_runtime_sources, dict(reviewed_runtime_sources)):
        raise ContinuationError("Signed review/runtime gate source relation drift")
    mandatory_gate_executables = validate_mandatory_gate_executable_records(gates)
    observed_asset_parts = tuple(
        sorted(
            path
            for path in PORTABLE_ASSET_ROOT.iterdir()
            if path.name.startswith("portable-artifact-reader.html.gz.b64.part")
        )
    )
    if observed_asset_parts != PORTABLE_READER_ASSET_PARTS:
        raise ContinuationError("Portable reader asset part membership drift")
    if mandatory_gate_executables["node"]["sha256"] != EXPECTED_NODE_SHA256:
        raise ContinuationError("Node identity drift")
    if (
        mandatory_gate_executables["qa_chromium"]["sha256"]
        != EXPECTED_QA_CHROMIUM_SHA256
    ):
        raise ContinuationError("Chromium identity drift")
    if (
        gates.record_for(FINAL_RELEASE_PUBLIC_KEY)["sha256"]
        != EXPECTED_FINAL_RELEASE_PUBLIC_KEY_SHA256
        or gates.record_for(REPORT_RELEASE_PUBLIC_KEY)["sha256"]
        != EXPECTED_REPORT_RELEASE_PUBLIC_KEY_SHA256
    ):
        raise ContinuationError("Terminal release public-key identity drift")

    payloads = {
        "manifest": strict_json(gates.data_for(MANIFEST), "all-sSNV manifest"),
        "primary": strict_json(gates.data_for(PRIMARY_PRE), "primary artifact audit"),
        "legacy_primary": strict_json(
            gates.data_for(LEGACY_PRIMARY_PRE),
            "legacy tumor-REF primary artifact audit",
        ),
        "run": strict_json(
            gates.data_for(COOCCURRENCE / "run_receipt.json"),
            "cooccurrence run receipt",
        ),
        "release": strict_json(
            gates.data_for(COOCCURRENCE / "release_receipt.json"),
            "cooccurrence release receipt",
        ),
        "summary": strict_json(
            gates.data_for(COOCCURRENCE / "summary.json"),
            "cooccurrence summary",
        ),
        "preflight": strict_json(
            gates.data_for(COOCCURRENCE_PREFLIGHT), "cooccurrence preflight"
        ),
        "tumor_ref": strict_json(
            gates.data_for(TUMOR_REF / "run_manifest.json"), "tumor-REF run manifest"
        ),
        "snapshot": strict_json(gates.data_for(SOURCE_SNAPSHOT), "source snapshot"),
        "m2": strict_json(gates.data_for(INDEPENDENT_M2_AUDIT), "independent M2 audit"),
    }
    for label in (
        "manifest",
        "primary",
        "legacy_primary",
        "tumor_ref",
        "snapshot",
        "m2",
    ):
        if payloads[label].get("pass") is not True:
            raise ContinuationError(f"Required upstream receipt is not pass=true: {label}")
    if (
        payloads["legacy_primary"].get("schema_name")
        != "intersubmod.stable_primary_artifact_audit"
        or payloads["legacy_primary"].get("schema_version") != "1.0.0"
    ):
        raise ContinuationError("Legacy tumor-REF primary-audit schema drift")
    require_basic_declared_artifact(
        payloads["tumor_ref"].get("inputs", {}).get(
            "primary_artifact_audit_pre"
        ),
        gates.record_for(LEGACY_PRIMARY_PRE),
        "tumor-REF legacy primary-audit input",
    )
    run = payloads["run"]
    release = payloads["release"]
    summary = payloads["summary"]
    preflight = payloads["preflight"]
    if (
        run.get("schema_name")
        != "intersubmod.methyl_ssnv_cooccurrence_validation.run_receipt"
        or run.get("schema_version") != "2.0.0"
        or run.get("task_type") != "B_comprehensive_validation"
        or run.get("scope") != "all_manifest_samples"
        or run.get("full_scope") is not True
        or run.get("pass") is not True
        or release.get("schema_name") != "intersubmod.cooccurrence_release_receipt"
        or release.get("schema_version") != "1.0.0"
        or release.get("release_status") != "RELEASE_RECONCILIATION_PASS"
        or release.get("pass") is not True
        or summary.get("schema_name")
        != "intersubmod.methyl_ssnv_cooccurrence_validation.summary"
        or summary.get("schema_version") != "2.0.0"
        or summary.get("pass") is not True
        or preflight.get("schema_name")
        != "intersubmod.cooccurrence_task_contract_preflight"
        or preflight.get("schema_version") != "2.0.0"
        or preflight.get("task_type") != "B_comprehensive_validation"
        or preflight.get("pass") is not True
    ):
        raise ContinuationError("Cooccurrence schema/scope/release contract drift")

    output_paths = {
        "site_tsv": COOCCURRENCE / "methyl_ssnv_site_results.tsv.gz",
        "pair_tsv": COOCCURRENCE / "methyl_ssnv_pair_results.tsv.gz",
        "summary_json": COOCCURRENCE / "summary.json",
        "raw_identity_duplicate_audit_tsv": (
            COOCCURRENCE / "raw_identity_duplicate_audit.tsv.gz"
        ),
        "case_json": COOCCURRENCE / "oracle_cases.json",
    }
    for role, path in output_paths.items():
        require_basic_declared_artifact(
            run.get("outputs", {}).get(role),
            gates.record_for(path),
            f"cooccurrence run output {role}",
        )
    require_basic_declared_artifact(
        run.get("inputs", {}).get("preflight_receipt"),
        gates.record_for(COOCCURRENCE_PREFLIGHT),
        "cooccurrence preflight input",
    )
    producer_input_roles = {
        "manifest": MANIFEST,
        "assignments": SCREEN_ASSIGNMENTS,
        "sites": SCREEN_SITES,
        "primary_artifact_audit_pre": PRIMARY_PRE,
        "independent_m2_audit": INDEPENDENT_M2_AUDIT,
    }
    for role, path in producer_input_roles.items():
        require_basic_declared_artifact(
            run.get("inputs", {}).get(role),
            gates.record_for(path),
            f"cooccurrence producer input {role}",
        )
        require_basic_declared_artifact(
            preflight.get("inputs", {}).get(role),
            gates.record_for(path),
            f"cooccurrence preflight input {role}",
        )
    release_input_roles = {
        "preflight": COOCCURRENCE_PREFLIGHT,
        "producer_receipt": COOCCURRENCE / "run_receipt.json",
        "summary": COOCCURRENCE / "summary.json",
        "sites": COOCCURRENCE / "methyl_ssnv_site_results.tsv.gz",
        "pairs": COOCCURRENCE / "methyl_ssnv_pair_results.tsv.gz",
        "duplicates": COOCCURRENCE / "raw_identity_duplicate_audit.tsv.gz",
        "oracle": COOCCURRENCE / "oracle_cases.json",
    }
    for role, path in release_input_roles.items():
        require_basic_declared_artifact(
            release.get("inputs", {}).get(role),
            gates.record_for(path),
            f"cooccurrence release input {role}",
        )

    expected_counts = {
        "stable_sites": 102_842,
        "focal_partner_pairs": 407_738,
        "m2_exact_global_fdr_family_pairs": 147,
        "m2_global_fdr_eligible_sites": 919,
        "multi_marker_molecular_haplotype_base_candidates": 0,
        "pair_by_confirmed": 7,
        "raw_identity_sparse_duplicate_rows": 196_706,
        "raw_identity_expected_projection_occurrences": 9_356_980,
    }
    counts = run.get("counts")
    recomputed = release.get("recomputed")
    reconciliation = summary.get("site_pair_count_reconciliation")
    raw = summary.get("raw_bam_identity_recovery_audit")
    observed_raw = preflight.get("observed", {}).get("raw_bam_identity_recovery")
    release_checks = release.get("checks")
    if (
        not isinstance(counts, Mapping)
        or any(
            not strict_json_equal(counts.get(key), value)
            for key, value in expected_counts.items()
        )
        or counts.get("raw_identity_all_site_results_passed_invariant_validation")
        is not True
        or counts.get("raw_identity_missing_projection_policy")
        != "hard_fail_before_site_result"
        or counts.get("raw_identity_conflicting_analysis_payload_policy")
        != "hard_fail_before_site_result"
        or counts.get("raw_identity_failure_counts_materialized") is not False
        or not isinstance(recomputed, Mapping)
        or not strict_json_equal(recomputed.get("stable_sites"), 102_842)
        or not strict_json_equal(recomputed.get("focal_partner_pairs"), 407_738)
        or not strict_json_equal(
            recomputed.get("multi_marker_molecular_haplotype_base_candidates"), 0
        )
        or not isinstance(reconciliation, Mapping)
        or not strict_json_equal(reconciliation.get("n_sites_reconciled"), 102_842)
        or not strict_json_equal(
            reconciliation.get("n_pair_rows_reconciled"), 407_738
        )
        or reconciliation.get("pass") is not True
        or not isinstance(raw, Mapping)
        or raw.get("equivalence_contract") != "sam_core_and_all_aux_tags_except_RG_exact_v1"
        or not strict_json_equal(raw.get("allowed_differing_auxiliary_tags"), ["RG"])
        or raw.get("all_site_results_passed_invariant_validation") is not True
        or raw.get("missing_projection_policy") != "hard_fail_before_site_result"
        or raw.get("conflicting_analysis_payload_policy")
        != "hard_fail_before_site_result"
        or raw.get("failure_counts_materialized") is not False
        or not isinstance(observed_raw, Mapping)
        or observed_raw.get("site_weighted_audit_sha256")
        != raw.get("site_weighted_audit_sha256")
        or not strict_json_equal(
            observed_raw.get("aggregate", {}).get("site_tasks"), 102_842
        )
        or not strict_json_equal(
            observed_raw.get("aggregate", {}).get("expected_projection_occurrences"),
            9_356_980,
        )
        or not isinstance(release_checks, Mapping)
        or not release_checks
        or any(value is not True for value in release_checks.values())
    ):
        raise ContinuationError("Cooccurrence independent core predicate recount drift")
    preflight_source_after = preflight.get("code", {}).get("source_identity_after")
    expected_preflight_sources = dict(run.get("code", {}))
    expected_preflight_sources["preflight"] = (
        preflight_source_after.get("preflight")
        if isinstance(preflight_source_after, Mapping)
        else None
    )
    if not strict_json_equal(preflight_source_after, expected_preflight_sources):
        raise ContinuationError("Cooccurrence preflight/run source identity reconciliation drift")
    if (
        not strict_json_equal(
            run.get("source_lock", {}).get("source_identity_before"), run.get("code")
        )
        or not strict_json_equal(
            run.get("source_lock", {}).get("source_identity_after_compute"),
            run.get("code"),
        )
        or not strict_json_equal(
            run.get("source_lock", {}).get("source_identity_after_output"),
            run.get("code"),
        )
        or run.get("source_lock", {}).get("all_sources_read_only_and_unchanged") is not True
    ):
        raise ContinuationError("Cooccurrence producer source lock drift")
    gates.require_paths_still_bound()
    return {
        "descriptor_count": gates.descriptor_count,
        "records": {role: gates.record_for(path) for role, path in GATE_INPUT_PATHS.items()},
        "alias_records": {"tumor_ref_summary": tumor_ref_summary_alias},
        "reviewed_runtime_sources": gate_runtime_sources,
        "cooccurrence_counts": {key: counts[key] for key in expected_counts},
        "release_status": release["release_status"],
        "raw_identity_site_weighted_audit_sha256": raw["site_weighted_audit_sha256"],
        "pass": True,
    }


def render_runtime_fd_bridge(primary_python_fd: int, qa_python_fd: int) -> bytes:
    if (
        type(primary_python_fd) is not int
        or type(qa_python_fd) is not int
        or primary_python_fd <= 2
        or qa_python_fd <= 2
        or primary_python_fd == qa_python_fd
    ):
        raise ContinuationError("Runtime Python descriptor bridge FD drift")
    rendered = RUNTIME_FD_BRIDGE_TEMPLATE
    fd_tokens = {
        b"__PRIMARY_PYTHON_FD__": os.fsencode(str(primary_python_fd)),
        b"__QA_PYTHON_FD__": os.fsencode(str(qa_python_fd)),
    }
    for token, value in fd_tokens.items():
        if rendered.count(token) != 1:
            raise ContinuationError("Runtime Python descriptor bridge token drift")
        rendered = rendered.replace(token, value)
    if any(token in rendered for token in fd_tokens):
        raise ContinuationError("Runtime Python descriptor bridge rendering incomplete")
    return rendered


def compose_downstream_script(
    runner_data: bytes,
    *,
    primary_python_fd: int,
    qa_python_fd: int,
    primary_python_record: Mapping[str, Any],
    qa_python_record: Mapping[str, Any],
) -> tuple[bytes, dict[str, Any]]:
    if sha256_bytes(runner_data) != EXPECTED_RUNNER_SHA256:
        raise ContinuationError("Completion runner full-source SHA drift")
    lines = runner_data.splitlines(keepends=True)
    if len(lines) != EXPECTED_RUNNER_LINE_COUNT or not runner_data.endswith(b"\n"):
        raise ContinuationError("Completion runner physical-line contract drift")
    source_prefix = b"".join(lines[:358])
    excluded = b"".join(lines[358:402])
    leaf = b"".join(lines[402:])
    if (
        sha256_bytes(source_prefix) != EXPECTED_RUNNER_PREFIX_SHA256
        or sha256_bytes(leaf) != EXPECTED_RUNNER_LEAF_SHA256
        or lines[358] not in (b"\n", b"\r\n")
        or not lines[359].startswith(b"for path in \\")
        or not lines[394].startswith(b"run_step \"Verify bounded tumor-REF")
        or not lines[402].startswith(b"run_step \"Run strict multi-seed")
    ):
        raise ContinuationError("Completion runner bridge boundaries drift")
    if (
        type(primary_python_fd) is not int
        or type(qa_python_fd) is not int
        or primary_python_fd <= 2
        or qa_python_fd <= 2
        or primary_python_fd == qa_python_fd
        or primary_python_record.get("path") != str(PRIMARY_PYTHON_RUNTIME)
        or qa_python_record.get("path") != str(QA_PYTHON_RUNTIME)
    ):
        raise ContinuationError("Completion runner Python descriptor contract drift")
    rewrite_contracts = {
        "primary_python": {
            "launcher": str(PYTHON),
            "wrapper_function": PRIMARY_PYTHON_WRAPPER,
            "canonical_argv0": str(PYTHON),
            "bound_executable_fd_path": f"/proc/self/fd/{primary_python_fd}",
            "fd": primary_python_fd,
            "runtime": dict(primary_python_record),
        },
        "qa_python": {
            "launcher": str(QA_PYTHON),
            "wrapper_function": QA_PYTHON_WRAPPER,
            "canonical_argv0": str(QA_PYTHON),
            "bound_executable_fd_path": f"/proc/self/fd/{qa_python_fd}",
            "fd": qa_python_fd,
            "runtime": dict(qa_python_record),
        },
    }
    executed_prefix = source_prefix
    for label, contract in rewrite_contracts.items():
        launcher = os.fsencode(contract["launcher"])
        replacement = os.fsencode(contract["wrapper_function"])
        if executed_prefix.count(launcher) != 1:
            raise ContinuationError(f"Completion runner {label} launcher count drift")
        executed_prefix = executed_prefix.replace(launcher, replacement)
    recovery_script_rewrites = {
        b"${SCRIPT_ROOT}/finalize_task_b_result_release.py": os.fsencode(
            str(FINAL_RELEASE_FINALIZER)
        ),
        b"${SCRIPT_ROOT}/build_all_ssnv_final_report_dataset.py": os.fsencode(
            str(FINAL_DATASET_BUILDER)
        ),
        b"${SCRIPT_ROOT}/build_all_ssnv_report_artifact.py": os.fsencode(
            str(REPORT_BUILDER)
        ),
    }
    executed_leaf = leaf
    for source_token, recovery_path in recovery_script_rewrites.items():
        source_count = executed_prefix.count(source_token) + executed_leaf.count(
            source_token
        )
        if source_count != 1:
            raise ContinuationError(
                "Completion runner recovery-script token count drift: "
                f"{os.fsdecode(source_token)} count={source_count}"
            )
        executed_prefix = executed_prefix.replace(source_token, recovery_path)
        executed_leaf = executed_leaf.replace(source_token, recovery_path)
    current_release_key_rewrites = {}
    for source, replacement, role in RUNNER_CURRENT_RELEASE_KEY_REBASE:
        source_token = os.fsencode(source)
        replacement_token = os.fsencode(replacement)
        source_count = executed_prefix.count(source_token) + executed_leaf.count(
            source_token
        )
        if source_count != 1:
            raise ContinuationError(
                "Completion runner current-release key token count drift: "
                f"{role} count={source_count}"
            )
        executed_prefix = executed_prefix.replace(source_token, replacement_token)
        executed_leaf = executed_leaf.replace(source_token, replacement_token)
        current_release_key_rewrites[role] = {
            "source": source,
            "executed_value": replacement,
            "source_occurrence_count": source_count,
        }
    for source, replacement, role in RUNNER_CURRENT_RELEASE_KEY_REBASE:
        source_token = os.fsencode(source)
        replacement_token = os.fsencode(replacement)
        if source_token in executed_prefix or source_token in executed_leaf:
            raise ContinuationError(
                f"Completion runner stale current-release key token survived: {role}"
            )
        replacement_count = executed_prefix.count(
            replacement_token
        ) + executed_leaf.count(replacement_token)
        if replacement_count != 1:
            raise ContinuationError(
                "Completion runner current-release replacement count drift: "
                f"{role} count={replacement_count}"
            )
        current_release_key_rewrites[role][
            "executed_occurrence_count"
        ] = replacement_count
    rendered_runtime_bridge = render_runtime_fd_bridge(
        primary_python_fd, qa_python_fd
    )
    composed = (
        executed_prefix
        + rendered_runtime_bridge
        + BRIDGE_BYTES
        + executed_leaf
        + EPILOGUE_BYTES
    )
    return composed, {
        "source": {
            "path": str(COMPLETION_RUNNER),
            "size_bytes": len(runner_data),
            "sha256": sha256_bytes(runner_data),
            "physical_lines": len(lines),
        },
        "included_prefix": {
            "physical_lines": "1-358",
            "size_bytes": len(source_prefix),
            "sha256": sha256_bytes(source_prefix),
            "source_bytes_executed_after_wrapper_rewrite": False,
        },
        "runtime_fd_canonical_argv0_wrappers": rewrite_contracts,
        "executed_runtime_bound_prefix": {
            "physical_lines": "1-358",
            "size_bytes": len(executed_prefix),
            "sha256": sha256_bytes(executed_prefix),
            "all_python_launchers_use_bound_descriptors_with_canonical_argv0": True,
        },
        "recovery_script_rewrites": {
            os.fsdecode(source_token): {
                "executed_path": os.fsdecode(recovery_path),
                "source_occurrence_count": 1,
            }
            for source_token, recovery_path in sorted(
                recovery_script_rewrites.items(), key=lambda item: item[0]
            )
        },
        "current_release_key_rewrites": current_release_key_rewrites,
        "excluded_historical_segment": {
            "physical_lines": "359-402",
            "size_bytes": len(excluded),
            "sha256": sha256_bytes(excluded),
            "executed": False,
            "reason": "canonical_source_receipt_already_exists_under_signed_promotion",
        },
        "reviewed_static_bridge": {
            "size_bytes": len(BRIDGE_BYTES),
            "sha256": sha256_bytes(BRIDGE_BYTES),
        },
        "reviewed_runtime_fd_bridge_template": {
            "size_bytes": len(RUNTIME_FD_BRIDGE_TEMPLATE),
            "sha256": sha256_bytes(RUNTIME_FD_BRIDGE_TEMPLATE),
        },
        "rendered_runtime_fd_bridge": {
            "size_bytes": len(rendered_runtime_bridge),
            "sha256": sha256_bytes(rendered_runtime_bridge),
            "placeholder_count": 0,
            "execution_semantics": "bound_fd_bytes_with_canonical_argv0",
        },
        "included_downstream_leaf": {
            "physical_lines": "403-638",
            "size_bytes": len(leaf),
            "sha256": sha256_bytes(leaf),
            "executed_size_bytes": len(executed_leaf),
            "executed_sha256": sha256_bytes(executed_leaf),
            "source_bytes_executed_after_recovery_script_rewrite": False,
        },
        "reviewed_static_epilogue": {
            "purpose": "tighten_completion_log_to_mode_0444_before_child_exit",
            "size_bytes": len(EPILOGUE_BYTES),
            "sha256": sha256_bytes(EPILOGUE_BYTES),
        },
        "stdin": {
            "size_bytes": len(composed),
            "sha256": sha256_bytes(composed),
        },
    }


def run_downstream_leaf(
    reader: Any,
    chain: Mapping[str, Any],
    script: bytes,
) -> dict[str, Any]:
    bash_fd, _, bash_record, _ = open_bound(reader, BASH, expected_mode="0o755")
    if bash_record["sha256"] != EXPECTED_BASH_SHA256:
        raise ContinuationError("Bash runtime identity drift")
    primary_python_fd = int(chain["python_fd"])
    qa_python_fd = int(chain["qa_python_fd"])
    if len({bash_fd, primary_python_fd, qa_python_fd}) != 3:
        raise ContinuationError("Downstream runtime descriptors are not distinct")
    reader.require_paths_still_bound()
    require_bootstrap_stable()
    require_incident_absent()
    result = subprocess.run(
        [f"/proc/self/fd/{bash_fd}", "-s"],
        input=script,
        cwd=REPO_ROOT,
        env=EXPECTED_ENVIRONMENT,
        pass_fds=(bash_fd, primary_python_fd, qa_python_fd),
        check=False,
    )
    if result.returncode != 0:
        raise ContinuationError(f"Downstream completion leaf failed with exit {result.returncode}")
    return {
        "actual_argv": [f"/proc/self/fd/{bash_fd}", "-s"],
        "normalized_argv": ["/proc/self/fd/<bound-bash-fd>", "-s"],
        "pass_fds": {
            str(bash_fd): bash_record,
            str(primary_python_fd): chain["python_record"],
            str(qa_python_fd): chain["qa_python_record"],
        },
        "exit_code": result.returncode,
    }


def expected_finalizer_command(mode: str, *, python_launcher: str) -> list[str]:
    prefix = [
        python_launcher,
        "-I",
        "-X",
        f"pycache_prefix={ORIGINAL_RUNNER_CACHE}",
        str(FINAL_RELEASE_FINALIZER),
    ]
    if mode == "create-dataset":
        return [
            *prefix,
            "--create",
            "--final-dataset-dir",
            str(FINAL_DATASET),
            "--builder-receipt",
            str(FINAL_DATASET / "run_receipt.json"),
            "--output",
            str(FINAL_DATASET_RELEASE_RECEIPT),
        ]
    if mode == "create-report":
        return [
            *prefix,
            "--create-report",
            "--report-dir",
            str(FINAL_REPORT),
            "--dataset-release-receipt",
            str(FINAL_DATASET_RELEASE_RECEIPT),
            "--dataset-release-signature",
            str(FINAL_DATASET_RELEASE_SIGNATURE),
            "--output",
            str(FINAL_REPORT_RELEASE_RECEIPT),
        ]
    raise ContinuationError(f"Unknown finalizer command mode: {mode}")


def enumerate_child_output_trees() -> dict[str, dict[str, Any]]:
    inventory: dict[str, dict[str, Any]] = {}
    root_specs = (
        *((role, path, True) for role, path in REQUIRED_CHILD_OUTPUT_ROOTS.items()),
        *((role, path, False) for role, path in OPTIONAL_CHILD_OUTPUT_ROOTS.items()),
    )
    for role, root, required in root_specs:
        exists = os.path.lexists(root)
        if required and not exists:
            raise ContinuationError(f"Required child output root is absent: {root}")
        if not exists:
            inventory[role] = {
                "path": str(root),
                "required": False,
                "exists": False,
                "directories": [],
                "files": [],
            }
            continue
        root_stat = os.lstat(root)
        if (
            stat.S_ISLNK(root_stat.st_mode)
            or not stat.S_ISDIR(root_stat.st_mode)
            or root.resolve(strict=True) != root
        ):
            raise ContinuationError(f"Child output root is not canonical: {root}")
        directories: list[Path] = []
        files: list[Path] = []
        for current, dir_names, file_names in os.walk(root, topdown=True, followlinks=False):
            current_path = Path(current)
            current_stat = os.lstat(current_path)
            if (
                stat.S_ISLNK(current_stat.st_mode)
                or not stat.S_ISDIR(current_stat.st_mode)
                or current_path.resolve(strict=True) != current_path
            ):
                raise ContinuationError(
                    f"Child output directory is not canonical: {current_path}"
                )
            directories.append(current_path)
            dir_names.sort()
            file_names.sort()
            for name in dir_names:
                candidate = current_path / name
                candidate_stat = os.lstat(candidate)
                if stat.S_ISLNK(candidate_stat.st_mode) or not stat.S_ISDIR(
                    candidate_stat.st_mode
                ):
                    raise ContinuationError(
                        f"Child output tree contains a non-directory entry: {candidate}"
                    )
            for name in file_names:
                candidate = current_path / name
                candidate_stat = os.lstat(candidate)
                if (
                    stat.S_ISLNK(candidate_stat.st_mode)
                    or not stat.S_ISREG(candidate_stat.st_mode)
                    or candidate.resolve(strict=True) != candidate
                ):
                    raise ContinuationError(
                        f"Child output tree contains a non-regular file: {candidate}"
                    )
                files.append(candidate)
        if not files:
            raise ContinuationError(f"Child output root contains no files: {root}")
        inventory[role] = {
            "path": str(root),
            "required": required,
            "exists": True,
            "directories": sorted(directories),
            "files": sorted(files),
        }
    return inventory


def bind_complete_child_output_inventory(
    outputs: StreamBoundArtifactSet,
) -> tuple[dict[str, Any], set[Path], set[Path], Callable[[], None]]:
    path_inventory = enumerate_child_output_trees()
    bound_paths: set[Path] = set()
    watched_directories: set[Path] = set()
    public_roots: dict[str, Any] = {}
    for role, item in path_inventory.items():
        if not item["exists"]:
            public_roots[role] = {
                "path": item["path"],
                "required": False,
                "exists": False,
                "file_count": 0,
                "files": [],
            }
            continue
        watched_directories.update(item["directories"])
        records: list[dict[str, Any]] = []
        root = Path(item["path"])
        for path in item["files"]:
            outputs.open(
                path,
                capture_limit=64 * 1024 * 1024,
                expected_mode="0o444",
                expected_link_count=1,
            )
            bound_paths.add(path)
            records.append(
                {
                    "relative_path": str(path.relative_to(root)),
                    "artifact": basic_projection(outputs.record_for(path)),
                }
            )
        public_roots[role] = {
            "path": item["path"],
            "required": item["required"],
            "exists": True,
            "file_count": len(records),
            "files": records,
        }
    standalone: dict[str, Any] = {}
    for role, path in REQUIRED_CHILD_OUTPUT_FILES.items():
        outputs.open(
            path,
            capture_limit=64 * 1024 * 1024,
            expected_mode="0o444",
            expected_link_count=1,
        )
        bound_paths.add(path)
        standalone[role] = basic_projection(outputs.record_for(path))

    expected_tree_paths = {
        role: {
            "directories": tuple(item["directories"]),
            "files": tuple(item["files"]),
            "exists": item["exists"],
        }
        for role, item in path_inventory.items()
    }

    def require_child_inventory_stable() -> None:
        live = enumerate_child_output_trees()
        live_tree_paths = {
            role: {
                "directories": tuple(item["directories"]),
                "files": tuple(item["files"]),
                "exists": item["exists"],
            }
            for role, item in live.items()
        }
        if live_tree_paths != expected_tree_paths:
            raise ContinuationError("Complete child output tree inventory changed")
        outputs.require_paths_still_bound()

    inventory_payload = {
        "roots": public_roots,
        "standalone_files": standalone,
        "bound_file_count": len(bound_paths),
        "watched_directory_count": len(watched_directories),
        "inventory_sha256": sha256_bytes(
            json.dumps(
                {"roots": public_roots, "standalone_files": standalone},
                ensure_ascii=True,
                sort_keys=True,
                separators=(",", ":"),
            ).encode("ascii")
        ),
        "all_bound_files_regular_mode_0444_link_count_one": True,
        "pass": True,
    }
    return (
        inventory_payload,
        bound_paths,
        watched_directories,
        require_child_inventory_stable,
    )


def validate_terminal_releases(
    authority_module: Any,
    trust_reader: Any,
    v7_validation: Mapping[str, Any],
    reviewed_runtime_sources: Mapping[str, Any],
    primary_python_launcher: str,
    result_key_transition: AuthorizedModeTransition | None,
    report_key_transition: AuthorizedModeTransition | None,
) -> tuple[dict[str, Any], Callable[[], None]]:
    outputs = StreamBoundArtifactSet()
    (
        child_output_inventory,
        child_output_paths,
        child_output_directories,
        require_child_inventory_stable,
    ) = bind_complete_child_output_inventory(outputs)
    output_reader = authority_module.BoundArtifactReader()
    for path in (
        FINAL_DATASET / "run_receipt.json",
        *FINAL_DATASET_OUTPUTS.values(),
        *REPORT_OUTPUTS.values(),
        ORIGINAL_RUNNER_LOG,
    ):
        outputs.open(
            path,
            capture_limit=64 * 1024 * 1024,
            expected_mode="0o444",
            expected_link_count=1,
        )

    openssl_fd, _, openssl_record, _ = open_bound(
        output_reader, OPENSSL, expected_mode="0o755"
    )
    if openssl_record["sha256"] != EXPECTED_OPENSSL_SHA256:
        raise ContinuationError("Terminal OpenSSL identity drift")
    dataset_fd, dataset_receipt_data, dataset_receipt_record, dataset_receipt_stat = open_bound(
        output_reader, FINAL_DATASET_RELEASE_RECEIPT
    )
    dataset_sig_fd, _, dataset_sig_record, dataset_sig_stat = open_bound(
        output_reader, FINAL_DATASET_RELEASE_SIGNATURE
    )
    result_public_fd, _, result_public_record, _ = open_bound(
        output_reader, FINAL_RELEASE_PUBLIC_KEY
    )
    report_fd, report_receipt_data, report_receipt_record, report_receipt_stat = open_bound(
        output_reader, FINAL_REPORT_RELEASE_RECEIPT
    )
    report_sig_fd, _, report_sig_record, report_sig_stat = open_bound(
        output_reader, FINAL_REPORT_RELEASE_SIGNATURE
    )
    report_public_fd, _, report_public_record, _ = open_bound(
        output_reader, REPORT_RELEASE_PUBLIC_KEY
    )
    if (
        dataset_sig_stat.st_size != 64
        or report_sig_stat.st_size != 64
        or result_public_record["sha256"] != EXPECTED_FINAL_RELEASE_PUBLIC_KEY_SHA256
        or report_public_record["sha256"] != EXPECTED_REPORT_RELEASE_PUBLIC_KEY_SHA256
    ):
        raise ContinuationError("Terminal release signature/key identity drift")
    if not authority_module.verify_ed25519_signature_fds(
        openssl_fd=openssl_fd,
        data_fd=dataset_fd,
        public_key_fd=result_public_fd,
        signature_fd=dataset_sig_fd,
    ):
        raise ContinuationError("Final dataset release signature failed")
    if not authority_module.verify_ed25519_signature_fds(
        openssl_fd=openssl_fd,
        data_fd=report_fd,
        public_key_fd=report_public_fd,
        signature_fd=report_sig_fd,
    ):
        raise ContinuationError("Final report release signature failed")
    if result_key_transition is None:
        _, retired_result_private, retired_result_stat = open_private_metadata(
            output_reader, FINAL_RELEASE_PRIVATE_KEY, "0o0"
        )
        if retired_result_stat.st_nlink != 1:
            raise ContinuationError("Final result private-key link count drift")
    else:
        retired_result_private = result_key_transition.verify_retired(output_reader)
    if report_key_transition is None:
        _, retired_report_private, retired_report_stat = open_private_metadata(
            output_reader, REPORT_RELEASE_PRIVATE_KEY, "0o0"
        )
        if retired_report_stat.st_nlink != 1:
            raise ContinuationError("Final report private-key link count drift")
    else:
        retired_report_private = report_key_transition.verify_retired(output_reader)

    dataset_release = strict_json(dataset_receipt_data, "final dataset release receipt")
    report_release = strict_json(report_receipt_data, "final report release receipt")
    require_exact_keys(dataset_release, DATASET_RELEASE_KEYS, "final dataset release")
    require_exact_keys(report_release, REPORT_RELEASE_KEYS, "final report release")
    builder_receipt = strict_json(
        outputs.data_for(FINAL_DATASET / "run_receipt.json"), "final builder receipt"
    )
    final_dataset = strict_json(
        outputs.data_for(FINAL_DATASET_OUTPUTS["final_report_dataset"]), "final dataset"
    )
    report_build = strict_json(
        outputs.data_for(REPORT_OUTPUTS["report_build_receipt"]), "report build receipt"
    )
    portable = strict_json(
        outputs.data_for(REPORT_OUTPUTS["portable_delivery_receipt"]),
        "portable delivery receipt",
    )
    desktop_qa = strict_json(
        outputs.data_for(REPORT_OUTPUTS["portable_desktop_qa"]), "desktop QA"
    )
    mobile_qa = strict_json(
        outputs.data_for(REPORT_OUTPUTS["portable_mobile_qa"]), "mobile QA"
    )

    finalizer_record = trust_reader.open(
        FINAL_RELEASE_FINALIZER, include_mode=True
    )[2]
    builder_record = trust_reader.open(FINAL_DATASET_BUILDER, include_mode=True)[2]
    report_builder_record = trust_reader.open(REPORT_BUILDER, include_mode=True)[2]
    final_outputs = {
        role: basic_projection(outputs.record_for(path))
        for role, path in FINAL_DATASET_OUTPUTS.items()
    }
    report_outputs = {
        role: basic_projection(outputs.record_for(path))
        for role, path in REPORT_OUTPUTS.items()
    }
    if set(reviewed_runtime_sources) != set(REVIEWED_RUNTIME_SOURCE_CONTRACTS):
        raise ContinuationError("Reviewed runtime source role set drift")
    expected_portable_runtime = {
        "plugin_root": str(PORTABLE_PLUGIN_ROOT),
        "builder": basic_projection(
            reviewed_runtime_sources["portable_builder_module"]
        ),
        "chart_extractor": basic_projection(
            reviewed_runtime_sources["portable_chart_extractor_module"]
        ),
        "verifier": basic_projection(
            reviewed_runtime_sources["portable_verifier_module"]
        ),
    }
    expected_dataset_signature_contract = {
        "algorithm": "Ed25519",
        "public_key": result_public_record,
        "signed_artifact": str(FINAL_DATASET_RELEASE_RECEIPT),
        "signature": str(FINAL_DATASET_RELEASE_SIGNATURE),
        "private_key_lifecycle": "encrypted_one_time_key_chmod_000_after_signing",
        "private_key_path": str(FINAL_RELEASE_PRIVATE_KEY),
    }
    expected_report_signature_contract = {
        "algorithm": "Ed25519",
        "public_key": report_public_record,
        "signed_artifact": str(FINAL_REPORT_RELEASE_RECEIPT),
        "signature": str(FINAL_REPORT_RELEASE_SIGNATURE),
        "private_key_lifecycle": "encrypted_one_time_key_chmod_000_after_signing",
        "private_key_path": str(REPORT_RELEASE_PRIVATE_KEY),
    }
    if (
        dataset_release.get("schema_name")
        != "intersubmod.task_b_final_dataset_release_receipt"
        or dataset_release.get("schema_version") != "1.0.0"
        or dataset_release.get("task_type") != "B_comprehensive_validation"
        or dataset_release.get("scope") != "all_7_datasets_469849_sites_final_dataset"
        or dataset_release.get("release_status")
        != "AWAITING_ONE_TIME_ED25519_RESULT_SIGNATURE"
        or dataset_release.get("pass_semantics")
        != "terminal_result_identity_release_not_scientific_confirmation"
        or not strict_json_equal(
            dataset_release.get("command"),
            expected_finalizer_command(
                "create-dataset", python_launcher=primary_python_launcher
            ),
        )
        or not strict_json_equal(dataset_release.get("source_authority"), v7_validation)
        or not strict_json_equal(
            dataset_release.get("inputs", {}).get("builder_receipt"),
            basic_projection(outputs.record_for(FINAL_DATASET / "run_receipt.json")),
        )
        or not strict_json_equal(
            dataset_release.get("inputs", {}).get("final_outputs"), final_outputs
        )
        or not strict_json_equal(
            dataset_release.get("code"),
            {"result_release_finalizer": basic_projection(finalizer_record)},
        )
        or not strict_json_equal(
            dataset_release.get("source_modes"),
            {"result_release_finalizer": "0o444"},
        )
        or not strict_json_equal(
            dataset_release.get("result_signature_contract"),
            expected_dataset_signature_contract,
        )
        or not strict_json_equal(
            dataset_release.get("checks"),
            {name: True for name in DATASET_RELEASE_CHECKS},
        )
        or dataset_release.get("pass") is not True
    ):
        raise ContinuationError("Signed final dataset release relation drift")
    if (
        builder_receipt.get("schema_name")
        != "intersubmod.all_ssnv_final_report_dataset_run_receipt"
        or builder_receipt.get("schema_version") != "2.0.0"
        or builder_receipt.get("task_type") != "B_comprehensive_validation"
        or builder_receipt.get("formal_task_b_release_eligible") is not True
        or builder_receipt.get("pass_semantics")
        != "integration_integrity_only_not_scientific_confirmation"
        or not strict_json_equal(builder_receipt.get("source_authority"), v7_validation)
        or not strict_json_equal(
            builder_receipt.get("code"),
            {"final_report_dataset_builder": basic_projection(builder_record)},
        )
        or not strict_json_equal(
            builder_receipt.get("source_modes"),
            {"final_report_dataset_builder": "0o444"},
        )
        or not strict_json_equal(builder_receipt.get("outputs"), final_outputs)
        or builder_receipt.get("pass") is not True
        or final_dataset.get("schema_name")
        != "intersubmod.all_ssnv_final_report_dataset"
        or final_dataset.get("task_type") != "B_comprehensive_validation"
        or final_dataset.get("formal_task_b_release_eligible") is not True
        or not strict_json_equal(
            final_dataset.get("counts", {}).get("screen_sites"), 469_849
        )
        or not strict_json_equal(
            final_dataset.get("counts", {}).get("stable_sites"), 102_842
        )
        or final_dataset.get("pass") is not True
    ):
        raise ContinuationError("Final builder/dataset formal relation drift")
    signed_dataset = report_release.get("inputs", {}).get("signed_dataset_release")
    if (
        report_release.get("schema_name")
        != "intersubmod.task_b_final_report_release_receipt"
        or report_release.get("schema_version") != "1.0.0"
        or report_release.get("task_type") != "B_comprehensive_validation"
        or report_release.get("scope")
        != "all_7_datasets_469849_sites_final_report_and_portable_evidence"
        or report_release.get("release_status")
        != "AWAITING_ONE_TIME_ED25519_REPORT_SIGNATURE"
        or report_release.get("pass_semantics")
        != "terminal_report_identity_release_not_scientific_confirmation"
        or not strict_json_equal(
            report_release.get("command"),
            expected_finalizer_command(
                "create-report", python_launcher=primary_python_launcher
            ),
        )
        or not strict_json_equal(report_release.get("source_authority"), v7_validation)
        or not strict_json_equal(
            report_release.get("inputs", {}).get("report_outputs"), report_outputs
        )
        or not isinstance(signed_dataset, Mapping)
        or not strict_json_equal(
            signed_dataset.get("receipt"), basic_projection(dataset_receipt_record)
        )
        or not strict_json_equal(
            signed_dataset.get("signature"), basic_projection(dataset_sig_record)
        )
        or not strict_json_equal(signed_dataset.get("public_key"), result_public_record)
        or signed_dataset.get("verification_pass") is not True
        or not strict_json_equal(
            report_release.get("code"),
            {"report_release_finalizer": basic_projection(finalizer_record)},
        )
        or not strict_json_equal(
            report_release.get("source_modes"),
            {"report_release_finalizer": "0o444"},
        )
        or not strict_json_equal(
            report_release.get("report_signature_contract"),
            expected_report_signature_contract,
        )
        or not strict_json_equal(
            report_release.get("checks"),
            {name: True for name in REPORT_RELEASE_CHECKS},
        )
        or report_release.get("pass") is not True
    ):
        raise ContinuationError("Signed final report release relation drift")
    if (
        report_build.get("schema_name") != "intersubmod.all_ssnv_report_build_receipt"
        or report_build.get("schema_version") != "1.2.0"
        or report_build.get("task_type") != "B_comprehensive_validation"
        or report_build.get("formal_task_b_release_eligible") is not True
        or not strict_json_equal(
            report_build.get("code", {}).get("report_builder"),
            basic_projection(report_builder_record),
        )
        or not strict_json_equal(
            report_build.get("code", {}).get("final_report_dataset_builder"),
            basic_projection(builder_record),
        )
        or not strict_json_equal(
            report_build.get("outputs", {}).get("report_md"),
            basic_projection(outputs.record_for(REPORT_OUTPUTS["report_markdown"])),
        )
        or not strict_json_equal(
            report_build.get("outputs", {}).get("artifact_json"),
            basic_projection(outputs.record_for(REPORT_OUTPUTS["canonical_report_artifact"])),
        )
        or report_build.get("pass") is not True
        or portable.get("schema_name") != "intersubmod.portable_report_delivery_receipt"
        or portable.get("schema_version") != "1.0.0"
        or portable.get("status") != "PASS"
        or not strict_json_equal(
            portable.get("artifact"),
            basic_projection(outputs.record_for(REPORT_OUTPUTS["canonical_report_artifact"])),
        )
        or not strict_json_equal(
            portable.get("output"),
            basic_projection(outputs.record_for(REPORT_OUTPUTS["portable_html"])),
        )
        or not strict_json_equal(
            portable.get("official_runtime"), expected_portable_runtime
        )
        or portable.get("pass") is not True
    ):
        raise ContinuationError("Report build/portable relation drift")
    for label, qa in (("desktop", desktop_qa), ("mobile", mobile_qa)):
        numeric_fields = (
            "overlapCount",
            "documentScrollWidth",
            "documentClientWidth",
            "bodyScrollWidth",
            "bodyClientWidth",
        )
        if (
            not isinstance(qa, Mapping)
            or any(type(qa.get(field)) is not int for field in numeric_fields)
            or qa.get("pass") is not True
            or qa["overlapCount"] != 0
            or not strict_json_equal(qa.get("consoleErrors"), [])
            or qa["documentScrollWidth"] > qa["documentClientWidth"] + 1
            or qa["bodyScrollWidth"] > qa["bodyClientWidth"] + 1
        ):
            raise ContinuationError(f"Portable {label} QA drift")
    for role in (
        "portable_official_verification_screenshot",
        "portable_desktop_screenshot",
        "portable_mobile_screenshot",
    ):
        if os.pread(outputs.fd_for(REPORT_OUTPUTS[role]), 8, 0) != b"\x89PNG\r\n\x1a\n":
            raise ContinuationError(f"Report PNG signature drift: {role}")

    terminal_paths = {
        FINAL_DATASET / "run_receipt.json",
        *FINAL_DATASET_OUTPUTS.values(),
        *REPORT_OUTPUTS.values(),
        ORIGINAL_RUNNER_LOG,
        FINAL_DATASET_RELEASE_RECEIPT,
        FINAL_DATASET_RELEASE_SIGNATURE,
        FINAL_RELEASE_PUBLIC_KEY,
        FINAL_REPORT_RELEASE_RECEIPT,
        FINAL_REPORT_RELEASE_SIGNATURE,
        REPORT_RELEASE_PUBLIC_KEY,
        FAILED_V16_CONTINUATION_PUBLIC_KEY,
        FAILED_V17_CONTINUATION_PUBLIC_KEY,
        FAILED_V18_CONTINUATION_PUBLIC_KEY,
        FAILED_V19_CONTINUATION_PUBLIC_KEY,
        REJECTED_V20_CONTINUATION_PUBLIC_KEY,
        FAILED_V21_CONTINUATION_PUBLIC_KEY,
        REJECTED_V22_CONTINUATION_PUBLIC_KEY,
        REJECTED_V23_CONTINUATION_PUBLIC_KEY,
        REJECTED_V24_CONTINUATION_PUBLIC_KEY,
        FAILED_V25_CONTINUATION_PUBLIC_KEY,
        FAILED_V26_CONTINUATION_PUBLIC_KEY,
        REJECTED_V27_CONTINUATION_PUBLIC_KEY,
        FAILED_V28_CONTINUATION_PUBLIC_KEY,
        FAILED_V29_CONTINUATION_PUBLIC_KEY,
        RECOVERY_CONTINUATION_PUBLIC_KEY,
        AUTHORIZATION,
        PREPARE_EXIT_ATTESTATION,
        AUTHORIZATION_SIGNATURE,
        AUTHORIZATION_PUBLIC_KEY,
        COMPLETION,
        PROMOTE_EXIT_ATTESTATION,
        COMPLETION_SIGNATURE,
        COMPLETION_PUBLIC_KEY,
        HISTORICAL_SOURCE_RECEIPT,
        CANONICAL_SOURCE_RECEIPT,
        *child_output_paths,
    }
    require_child_inventory_stable()
    outputs.require_paths_still_bound()
    output_reader.require_paths_still_bound()
    trust_reader.require_paths_still_bound()
    if result_key_transition is not None:
        result_key_transition.require_retired_stable()
    if report_key_transition is not None:
        report_key_transition.require_retired_stable()
    require_bootstrap_stable()
    require_incident_absent()
    terminal_sentinel = PathMutationSentinel(
        terminal_paths,
        all_event_directories=child_output_directories,
    )
    outputs.retain_until_process_exit()
    output_reader.retain_until_process_exit()
    terminal_sentinel.retain_until_process_exit()

    def require_terminal_release_stable() -> None:
        require_child_inventory_stable()
        outputs.require_paths_still_bound()
        output_reader.require_paths_still_bound()
        trust_reader.require_paths_still_bound()
        if result_key_transition is not None:
            result_key_transition.require_retired_stable()
        if report_key_transition is not None:
            report_key_transition.require_retired_stable()
        terminal_sentinel.assert_clean()
        require_bootstrap_stable()
        require_incident_absent()

    public_result = {
        "dataset_release": {
            "receipt": with_stat(dataset_receipt_record, dataset_receipt_stat),
            "signature": dataset_sig_record,
            "public_key": result_public_record,
            "retired_private_key": retired_result_private,
            "signature_verified": True,
            "relations_reverified": True,
        },
        "report_release": {
            "receipt": with_stat(report_receipt_record, report_receipt_stat),
            "signature": report_sig_record,
            "public_key": report_public_record,
            "retired_private_key": retired_report_private,
            "signature_verified": True,
            "relations_reverified": True,
        },
        "final_dataset_counts": {
            "screen_sites": 469_849,
            "stable_sites": 102_842,
        },
        "report_output_count": len(report_outputs),
        "bound_output_descriptor_count": outputs.descriptor_count,
        "complete_child_output_inventory": child_output_inventory,
        "completion_log": basic_projection(outputs.record_for(ORIGINAL_RUNNER_LOG)),
        "all_declared_release_artifacts_mode_0444_link_count_one": True,
        "retired_private_keys_mode_000_link_count_one": True,
        "terminal_mutation_sentinel": {
            "protected_path_count": terminal_sentinel.protected_path_count,
            "protected_path_set_sha256": terminal_sentinel.protected_path_set_sha256,
            "retained_until_process_exit": True,
            "pass": True,
        },
        "pass": True,
    }
    return public_result, require_terminal_release_stable


def link_fd_no_replace(source_fd: int, parent_fd: int, target_name: str) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    linkat = getattr(libc, "linkat", None)
    if linkat is None:
        raise ContinuationError("linkat is required for terminal no-replace commit")
    linkat.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
    )
    linkat.restype = ctypes.c_int
    result = linkat(
        -100,
        os.fsencode(f"/proc/self/fd/{source_fd}"),
        parent_fd,
        os.fsencode(target_name),
        0x400,
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise ContinuationError(
            "Terminal atomic no-replace commit failed: " + os.strerror(error_number)
        )


def commit_terminal_receipt(
    payload: Mapping[str, Any],
    precommit: Callable[[], None],
    *,
    path: Path = CONTINUATION_RECEIPT,
    postcommit: Callable[[], None] | None = None,
) -> dict[str, Any]:
    data = encode_json(payload)
    require_absent(path, "Terminal receipt")
    if path.parent.resolve(strict=True) != path.parent or stat.S_ISLNK(
        os.lstat(path.parent).st_mode
    ):
        raise ContinuationError("Terminal receipt parent is not canonical")
    parent_flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0)
    if hasattr(os, "O_NOFOLLOW"):
        parent_flags |= os.O_NOFOLLOW
    parent_fd = os.open(path.parent, parent_flags)
    parent_opened = os.fstat(parent_fd)
    if not hasattr(os, "O_TMPFILE"):
        raise ContinuationError("O_TMPFILE is required for terminal staged commit")
    descriptor = os.open(
        ".",
        os.O_RDWR | os.O_TMPFILE | os.O_CLOEXEC,
        0o400,
        dir_fd=parent_fd,
    )
    offset = 0
    while offset < len(data):
        written = os.write(descriptor, data[offset:])
        if written <= 0:
            raise ContinuationError("Short write to terminal receipt staging inode")
        offset += written
    os.fchmod(descriptor, 0o444)
    os.fsync(descriptor)
    opened = os.fstat(descriptor)
    observed = b"".join(
        os.pread(descriptor, min(1024 * 1024, opened.st_size - position), position)
        for position in range(0, opened.st_size, 1024 * 1024)
    )
    if (
        not stat.S_ISREG(opened.st_mode)
        or opened.st_nlink != 0
        or stat.S_IMODE(opened.st_mode) != 0o444
        or observed != data
    ):
        raise ContinuationError("Terminal receipt staging validation failed")
    precommit()
    parent_live = os.stat(path.parent, follow_symlinks=False)
    if (
        parent_opened.st_dev,
        parent_opened.st_ino,
        parent_opened.st_mode,
    ) != (parent_live.st_dev, parent_live.st_ino, parent_live.st_mode):
        raise ContinuationError("Continuation receipt parent path rebound")
    require_absent(path, "Terminal receipt")
    _COMMIT_FDS.extend((parent_fd, descriptor))
    link_fd_no_replace(descriptor, parent_fd, path.name)
    committed = os.fstat(descriptor)
    linked = os.stat(path.name, dir_fd=parent_fd, follow_symlinks=False)
    reopened_flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        reopened_flags |= os.O_NOFOLLOW
    reopened_fd = os.open(path.name, reopened_flags, dir_fd=parent_fd)
    try:
        reopened = os.fstat(reopened_fd)
        reopened_data = b"".join(
            os.pread(
                reopened_fd,
                min(1024 * 1024, reopened.st_size - position),
                position,
            )
            for position in range(0, reopened.st_size, 1024 * 1024)
        )
        if (
            committed.st_nlink != 1
            or not same_stat(committed, linked)
            or not same_stat(committed, reopened)
            or reopened_data != data
        ):
            raise ContinuationError("Terminal receipt post-link path/inode drift")
    finally:
        os.close(reopened_fd)
    os.fsync(parent_fd)
    parent_after = os.stat(path.parent, follow_symlinks=False)
    if (
        parent_opened.st_dev,
        parent_opened.st_ino,
        parent_opened.st_mode,
    ) != (parent_after.st_dev, parent_after.st_ino, parent_after.st_mode):
        raise ContinuationError("Terminal receipt parent changed after durable commit")
    if path.resolve(strict=True) != path:
        raise ContinuationError("Terminal receipt path is not canonical after commit")
    durable_fd = os.open(path.name, reopened_flags, dir_fd=parent_fd)
    try:
        durable_source = os.fstat(descriptor)
        durable_opened = os.fstat(durable_fd)
        durable_linked = os.stat(path.name, dir_fd=parent_fd, follow_symlinks=False)
        durable_data = b"".join(
            os.pread(
                durable_fd,
                min(1024 * 1024, durable_opened.st_size - position),
                position,
            )
            for position in range(0, durable_opened.st_size, 1024 * 1024)
        )
        if (
            durable_source.st_nlink != 1
            or not same_stat(durable_source, durable_opened)
            or not same_stat(durable_source, durable_linked)
            or durable_data != data
        ):
            raise ContinuationError("Terminal receipt post-fsync path/inode drift")
    finally:
        os.close(durable_fd)
    if postcommit is not None:
        postcommit()
    committed_stat = os.fstat(descriptor)
    return {
        "descriptor": descriptor,
        "data": data,
        "artifact": with_stat(
            {
                "path": str(path),
                "size_bytes": len(data),
                "sha256": sha256_bytes(data),
                "mode": "0o444",
            },
            committed_stat,
        ),
    }


def catchable_signal_set() -> set[signal.Signals]:
    return {
        member
        for member in signal.valid_signals()
        if member not in {signal.SIGKILL, signal.SIGSTOP}
    }


def commit_terminal_success_witness(
    payload: Mapping[str, Any], precommit: Callable[[], None]
) -> NoReturn:
    """Publish the success witness as the supervisor's final no-return action."""
    data = encode_json(payload)
    path = CONTINUATION_SUCCESS_WITNESS
    require_absent(path, "Continuation supervisor-success witness")
    if path.parent.resolve(strict=True) != path.parent or stat.S_ISLNK(
        os.lstat(path.parent).st_mode
    ):
        raise ContinuationError("Success-witness parent is not canonical")
    parent_flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0)
    if hasattr(os, "O_NOFOLLOW"):
        parent_flags |= os.O_NOFOLLOW
    parent_fd = os.open(path.parent, parent_flags)
    parent_opened = os.fstat(parent_fd)
    if not hasattr(os, "O_TMPFILE"):
        raise ContinuationError("O_TMPFILE is required for success-witness commit")
    descriptor = os.open(
        ".",
        os.O_RDWR | os.O_TMPFILE | os.O_CLOEXEC,
        0o400,
        dir_fd=parent_fd,
    )
    _COMMIT_FDS.extend((parent_fd, descriptor))
    offset = 0
    while offset < len(data):
        written = os.write(descriptor, data[offset:])
        if written <= 0:
            raise ContinuationError("Short write to success-witness staging inode")
        offset += written
    os.fchmod(descriptor, 0o444)
    os.fsync(descriptor)
    staged = os.fstat(descriptor)
    staged_data = b"".join(
        os.pread(descriptor, min(1024 * 1024, staged.st_size - position), position)
        for position in range(0, staged.st_size, 1024 * 1024)
    )
    if (
        not stat.S_ISREG(staged.st_mode)
        or staged.st_nlink != 0
        or stat.S_IMODE(staged.st_mode) != 0o444
        or staged_data != data
    ):
        raise ContinuationError("Success-witness staging validation failed")

    precommit()
    parent_live = os.stat(path.parent, follow_symlinks=False)
    if not same_stat(parent_opened, parent_live):
        raise ContinuationError("Success-witness parent path rebound")
    require_absent(path, "Continuation supervisor-success witness")
    os.fsync(parent_fd)

    previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, catchable_signal_set())
    try:
        precommit()
        parent_terminal = os.stat(path.parent, follow_symlinks=False)
        staged_terminal = os.fstat(descriptor)
        if (
            not same_stat(parent_opened, parent_terminal)
            or not same_stat(staged, staged_terminal)
            or os.pread(descriptor, staged.st_size, 0) != data
        ):
            raise ContinuationError("Success-witness terminal token drift")
        require_absent(path, "Continuation supervisor-success witness")
        link_fd_no_replace(descriptor, parent_fd, path.name)
    except BaseException:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        raise

    # linkat success is the authority commit. A crash can only remove the
    # directory entry after reboot, which is fail-closed; no fallible work follows.
    os._exit(0)


def run() -> None:
    require_clean_runtime(mode="child")
    supervisor_capability = bind_supervisor_capability()
    require_absent(CONTINUATION_RECEIPT, "Continuation terminal receipt")
    require_absent(CONTINUATION_SIGNATURE, "Continuation terminal signature")
    require_incident_absent()
    for path in DOWNSTREAM_OUTPUT_SLOTS:
        require_absent(path, "Downstream output before continuation")

    authority_module, bootstrap_record = bootstrap_release_authority()
    v7_validation = validate_v7(authority_module)
    reader = authority_module.BoundArtifactReader()
    chain = validate_promotion_chain(
        authority_module, reader, v7_validation, bootstrap_record
    )
    validate_supervisor_capability_record(supervisor_capability, chain)
    continuation_public_record = chain["recovery_continuation_public_key_record"]
    continuation_public_stat = chain["recovery_continuation_public_key_stat"]
    continuation_private_record = chain["recovery_continuation_private_key_record"]
    continuation_private_stat = chain["recovery_continuation_private_key_stat"]
    if (
        continuation_public_record["sha256"]
        != EXPECTED_RECOVERY_CONTINUATION_PUBLIC_KEY_SHA256
        or continuation_public_stat.st_nlink != 1
        or continuation_private_stat.st_nlink != 1
    ):
        raise ContinuationError("Continuation terminal signing key contract drift")
    verification = run_fresh_verifier(reader, chain)
    replay = validate_replay_evidence(reader, chain, verification)

    gate_inputs = StreamBoundArtifactSet()
    gate_validation = validate_gate_inputs(
        gate_inputs, chain["reviewed_runtime_sources"]
    )
    report_key_transition = AuthorizedModeTransition(
        authority_module, REPORT_RELEASE_PRIVATE_KEY
    )
    result_key_transition = AuthorizedModeTransition(
        authority_module, FINAL_RELEASE_PRIVATE_KEY
    )
    composed_script, script_contract = compose_downstream_script(
        chain["runner_data"],
        primary_python_fd=int(chain["python_fd"]),
        qa_python_fd=int(chain["qa_python_fd"]),
        primary_python_record=chain["python_record"],
        qa_python_record=chain["qa_python_record"],
    )
    protected_paths = collect_protected_paths(
        chain["authorization"],
        chain["completion"],
        verification["recorded_payload"],
        verification["fresh_payload"],
        gate_validation["records"],
        replay,
    )
    protected_paths.add(FAILED_V16_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V17_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V18_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V19_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(REJECTED_V20_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V21_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(REJECTED_V22_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(REJECTED_V23_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(REJECTED_V24_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V25_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V26_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(REJECTED_V27_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V28_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V29_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(RECOVERY_CONTINUATION_PUBLIC_KEY)
    mutation_sentinel = PathMutationSentinel(
        protected_paths,
        all_event_directories={PORTABLE_ASSET_ROOT, TUMOR_REF},
    )

    reader.require_paths_still_bound()
    gate_inputs.require_paths_still_bound()
    require_bootstrap_stable()
    require_incident_absent()
    mutation_sentinel.assert_clean()
    child = run_downstream_leaf(reader, chain, composed_script)

    mutation_sentinel.assert_clean()
    reader.require_paths_still_bound()
    gate_inputs.require_paths_still_bound()
    require_bootstrap_stable()
    require_incident_absent()
    mutation_sentinel.assert_clean()
    terminal, require_terminal_release_stable = validate_terminal_releases(
        authority_module,
        reader,
        v7_validation,
        chain["reviewed_runtime_sources"],
        script_contract["runtime_fd_canonical_argv0_wrappers"]["primary_python"][
            "canonical_argv0"
        ],
        result_key_transition,
        report_key_transition,
    )
    if not ORIGINAL_RUNNER_CACHE.is_dir() or stat.S_IMODE(
        ORIGINAL_RUNNER_CACHE.stat().st_mode
    ) != 0o700:
        raise ContinuationError("Downstream Python cache root contract drift")

    reader.require_paths_still_bound()
    gate_inputs.require_paths_still_bound()
    require_bootstrap_stable()
    require_incident_absent()
    mutation_sentinel.assert_clean()
    gate_inputs.retain_until_process_exit()
    reader.retain_until_process_exit()
    mutation_sentinel.retain_until_process_exit()

    receipt = {
        "schema_name": CONTINUATION_SCHEMA,
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "task_type": "B_comprehensive_validation",
        "status": (
            "DOWNSTREAM_VALIDATED_AWAITING_SUPERVISOR_EXIT_ATTESTATION_SIGNATURE_"
            "AND_SUCCESS_WITNESS"
        ),
        "scope": "all_7_datasets_469849_sites_signed_task_b_release_chain",
        "canonical_command": supervised_child_command(),
        "clean_environment": EXPECTED_ENVIRONMENT,
        "supervisor_capability": supervisor_capability,
        "code": {
            "downstream_continuation": chain["self_record"],
            "promotion_tool": chain["promotion_record"],
            "continuation_verifier": chain["verifier_record"],
            "runner_gate_replay": chain["replay_record"],
            "completion_runner": chain["runner_record"],
            "release_source_authority_validator": chain["authority_record"],
            "schema_recovery_authority_validator": chain[
                "schema_recovery_validator_record"
            ],
        },
        "legacy_signed_code": {
            "downstream_continuation": chain["authorized_self_record"],
            "continuation_verifier": chain["authorized_verifier_record"],
            "runner_gate_replay": chain["authorized_replay_record"],
        },
        "schema_recovery_authority": chain["schema_recovery_authority"],
        "historical_identity_contract": {
            "historical_link_count_observation": "NOT_RECORDED_NOT_INFERRED",
            "scientific_payload_changed": False,
            "claim_ceiling_changed": False,
        },
        "self_binding": chain["self_binding"],
        "signed_promotion": {
            "authorization": chain["authorization_record"],
            "prepare_exit_attestation": chain["prepare_exit_attestation_record"],
            "authorization_signature": chain["authorization_signature_record"],
            "completion": chain["completion_record"],
            "promote_exit_attestation": chain["promote_exit_attestation_record"],
            "completion_signature": chain["completion_signature_record"],
            "historical_source_receipt": chain["source_record"],
            "canonical_source_receipt": chain["canonical_record"],
            "authorization_signature_verified": True,
            "completion_signature_verified": True,
            "bound_relation_count": chain["bound_relation_count"],
        },
        "fresh_promotion_verification": {
            "recorded_receipt": verification["recorded_record"],
            "fresh_stdout": verification["fresh_stdout"],
            "fresh_stderr": verification["fresh_stderr"],
            "fresh_actual_argv": verification["fresh_actual_argv"],
            "fresh_normalized_argv": verification["fresh_normalized_argv"],
            "fresh_pass_fds": verification["fresh_pass_fds"],
            "recorded_mode": "verify-and-record",
            "fresh_mode": "verify",
            "stable_fields_equal": True,
            "bound_relation_count": verification["bound_relation_count"],
            "pass": True,
        },
        "runner_only_replay": replay,
        "authority_semantics": dict(TERMINAL_AUTHORITY_SEMANTICS),
        "v7_validation": {
            "authority_id": v7_validation["authority_id"],
            "source_set_sha256": v7_validation["source_set_sha256"],
            "validation_binding_lease": v7_validation["validation_binding_lease"],
            "pass": True,
        },
        "parent_gate_validation": gate_validation,
        "mutation_monitor": {
            "mechanism": (
                "linux_inotify_parent_watch_before_leaf_snapshot_then_file_watch_"
                "with_mode000_parent_fallback_and_setup_identity_recheck"
            ),
            "all_event_directories": sorted(
                str(path) for path in (PORTABLE_ASSET_ROOT, TUMOR_REF)
            ),
            "protected_path_count": mutation_sentinel.protected_path_count,
            "protected_path_set_sha256": mutation_sentinel.protected_path_set_sha256,
            "direct_file_watch_count": mutation_sentinel.direct_file_watch_count,
            "parent_only_mode000_path_count": (
                mutation_sentinel.parent_only_mode000_path_count
            ),
            "parent_only_mode000_path_set_sha256": (
                mutation_sentinel.parent_only_mode000_path_set_sha256
            ),
            "parent_only_fallback_reason": (
                "inotify_eacces_on_exact_mode000_canonical_regular_file"
            ),
            "setup_order": (
                "named_parent_and_all_event_watches_before_leaf_snapshot_then_"
                "direct_file_watch_then_full_stat_identity_recheck"
            ),
            "setup_identity_recheck_count": (
                mutation_sentinel.setup_identity_recheck_count
            ),
            "setup_identity_recheck_path_set_sha256": (
                mutation_sentinel.setup_identity_recheck_path_set_sha256
            ),
            "transient_mutation_events_observed": 0,
            "retained_until_process_exit": True,
            "pass": True,
        },
        "composed_runner_stdin": script_contract,
        "child_execution": child,
        "terminal_releases": terminal,
        "terminal_signature_contract": {
            "algorithm": "Ed25519",
            "public_key": continuation_public_record,
            "private_key_pre_signature": with_stat(
                continuation_private_record, continuation_private_stat
            ),
            "execution_receipt": str(CONTINUATION_RECEIPT),
            "signed_artifact": str(CONTINUATION_EXIT_ATTESTATION),
            "signature": str(CONTINUATION_SIGNATURE),
            "supervisor_success_witness": str(CONTINUATION_SUCCESS_WITNESS),
            "required_private_key_post_signature_mode": "0o0",
            "unsigned_receipt_grants_release_authority": False,
            "success_witness_required_for_release_authority": True,
        },
        "governance": dict(TERMINAL_GOVERNANCE),
        "checks": dict(TERMINAL_CHECKS),
        "pass": True,
        "pass_semantics": TERMINAL_PASS_SEMANTICS,
    }
    def precommit() -> None:
        reader.require_paths_still_bound()
        gate_inputs.require_paths_still_bound()
        require_bootstrap_stable()
        require_incident_absent()
        mutation_sentinel.assert_clean()
        require_terminal_release_stable()

    commit_terminal_receipt(receipt, precommit, postcommit=precommit)


def validate_recorded_fresh_verification_claim(
    value: Any,
    verification: Mapping[str, Any],
    chain: Mapping[str, Any],
) -> None:
    expected_keys = {
        "recorded_receipt",
        "fresh_stdout",
        "fresh_stderr",
        "fresh_actual_argv",
        "fresh_normalized_argv",
        "fresh_pass_fds",
        "recorded_mode",
        "fresh_mode",
        "stable_fields_equal",
        "bound_relation_count",
        "pass",
    }
    if not isinstance(value, Mapping) or set(value) != expected_keys:
        raise ContinuationError("Fresh promotion verification claim schema drift")
    actual_argv = value.get("fresh_actual_argv")
    pass_fds = value.get("fresh_pass_fds")
    if (
        not isinstance(actual_argv, list)
        or len(actual_argv) != 5
        or actual_argv[1:3] != ["-I", "-B"]
        or actual_argv[4] != "--verify"
        or not isinstance(pass_fds, Mapping)
        or len(pass_fds) != 2
    ):
        raise ContinuationError("Fresh promotion verifier FD launch schema drift")
    try:
        python_fd = int(str(actual_argv[0]).removeprefix("/proc/self/fd/"))
        verifier_fd = int(str(actual_argv[3]).removeprefix("/proc/self/fd/"))
    except ValueError as error:
        raise ContinuationError("Fresh promotion verifier FD token drift") from error
    stdout_record = value.get("fresh_stdout")
    stderr_record = value.get("fresh_stderr")
    if (
        actual_argv[0] != f"/proc/self/fd/{python_fd}"
        or actual_argv[3] != f"/proc/self/fd/{verifier_fd}"
        or python_fd <= 2
        or verifier_fd <= 2
        or python_fd == verifier_fd
        or set(pass_fds) != {str(python_fd), str(verifier_fd)}
        or not strict_json_equal(
            pass_fds.get(str(python_fd)), chain["python_record"]
        )
        or not strict_json_equal(
            pass_fds.get(str(verifier_fd)), chain["verifier_record"]
        )
        or not strict_json_equal(
            value.get("recorded_receipt"), verification["recorded_record"]
        )
        or not strict_json_equal(
            value.get("fresh_normalized_argv"),
            [
                "/proc/self/fd/<bound-python-fd>",
                "-I",
                "-B",
                "/proc/self/fd/<bound-verifier-fd>",
                "--verify",
            ],
        )
        or not isinstance(stdout_record, Mapping)
        or set(stdout_record) != {"size_bytes", "sha256"}
        or type(stdout_record.get("size_bytes")) is not int
        or stdout_record["size_bytes"] <= 0
        or not is_canonical_sha256(stdout_record.get("sha256"))
        or not strict_json_equal(
            stderr_record, {"size_bytes": 0, "sha256": sha256_bytes(b"")}
        )
        or value.get("recorded_mode") != "verify-and-record"
        or value.get("fresh_mode") != "verify"
        or value.get("stable_fields_equal") is not True
        or type(value.get("bound_relation_count")) is not int
        or value["bound_relation_count"] <= 0
        or value.get("pass") is not True
    ):
        raise ContinuationError("Fresh promotion verification claim drift")


def validate_child_execution_claim(
    value: Any,
    reader: Any,
    chain: Mapping[str, Any],
) -> dict[str, int]:
    expected_keys = {"actual_argv", "normalized_argv", "pass_fds", "exit_code"}
    if not isinstance(value, Mapping) or set(value) != expected_keys:
        raise ContinuationError("Downstream child execution schema drift")
    _, _, bash_record, _ = open_bound(reader, BASH, expected_mode="0o755")
    actual_argv = value.get("actual_argv")
    pass_fds = value.get("pass_fds")
    if (
        not isinstance(actual_argv, list)
        or len(actual_argv) != 2
        or actual_argv[1] != "-s"
        or not isinstance(pass_fds, Mapping)
        or len(pass_fds) != 3
    ):
        raise ContinuationError("Downstream child FD launch schema drift")
    expected_records = {
        "bash": bash_record,
        "primary_python": chain["python_record"],
        "qa_python": chain["qa_python_record"],
    }
    bound_fds: dict[str, int] = {}
    for role, expected_record in expected_records.items():
        matches: list[int] = []
        for token, record in pass_fds.items():
            try:
                descriptor = int(token)
            except (TypeError, ValueError) as error:
                raise ContinuationError("Downstream child FD token drift") from error
            if str(descriptor) != token or descriptor <= 2:
                raise ContinuationError("Downstream child FD token is not canonical")
            if strict_json_equal(record, expected_record):
                matches.append(descriptor)
        if len(matches) != 1:
            raise ContinuationError(f"Downstream child {role} FD binding drift")
        bound_fds[role] = matches[0]
    if len(set(bound_fds.values())) != 3:
        raise ContinuationError("Downstream child runtime descriptors are not distinct")
    bash_fd = bound_fds["bash"]
    expected_pass_fds = {
        str(bound_fds["bash"]): bash_record,
        str(bound_fds["primary_python"]): chain["python_record"],
        str(bound_fds["qa_python"]): chain["qa_python_record"],
    }
    if (
        bash_fd <= 2
        or actual_argv[0] != f"/proc/self/fd/{bash_fd}"
        or not strict_json_equal(pass_fds, expected_pass_fds)
        or not strict_json_equal(
            value.get("normalized_argv"), ["/proc/self/fd/<bound-bash-fd>", "-s"]
        )
        or type(value.get("exit_code")) is not int
        or value.get("exit_code") != 0
    ):
        raise ContinuationError("Downstream child did not prove exit zero")
    return bound_fds


def validate_execution_receipt(
    *,
    authority_module: Any,
    bootstrap_record: Mapping[str, Any],
    v7_validation: Mapping[str, Any],
    reader: Any,
    chain: Mapping[str, Any],
    receipt_data: bytes,
    receipt_record: Mapping[str, Any],
    receipt_stat: os.stat_result,
    expected_supervisor_capability: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    payload = strict_json(receipt_data, "continuation execution receipt")
    require_exact_keys(payload, CONTINUATION_RECEIPT_KEYS, "continuation execution receipt")
    require_utc_timestamp(payload.get("created_at_utc"), "continuation execution timestamp")
    if (
        payload.get("schema_name") != CONTINUATION_SCHEMA
        or payload.get("schema_version") != "1.0.0"
        or payload.get("task_type") != "B_comprehensive_validation"
        or payload.get("status")
        != "DOWNSTREAM_VALIDATED_AWAITING_SUPERVISOR_EXIT_ATTESTATION_SIGNATURE_"
        "AND_SUCCESS_WITNESS"
        or payload.get("scope")
        != "all_7_datasets_469849_sites_signed_task_b_release_chain"
        or not strict_json_equal(
            payload.get("canonical_command"), supervised_child_command()
        )
        or not strict_json_equal(payload.get("clean_environment"), EXPECTED_ENVIRONMENT)
        or payload.get("pass") is not True
        or payload.get("pass_semantics") != TERMINAL_PASS_SEMANTICS
    ):
        raise ContinuationError("Continuation execution receipt top-level drift")
    supervisor_capability = validate_supervisor_capability_record(
        payload.get("supervisor_capability"), chain
    )
    if (
        expected_supervisor_capability is not None
        and not strict_json_equal(
            supervisor_capability, dict(expected_supervisor_capability)
        )
    ):
        raise ContinuationError("Execution receipt supervisor capability drift")

    expected_code = {
        "downstream_continuation": chain["self_record"],
        "promotion_tool": chain["promotion_record"],
        "continuation_verifier": chain["verifier_record"],
        "runner_gate_replay": chain["replay_record"],
        "completion_runner": chain["runner_record"],
        "release_source_authority_validator": chain["authority_record"],
        "schema_recovery_authority_validator": chain[
            "schema_recovery_validator_record"
        ],
    }
    expected_legacy_signed_code = {
        "downstream_continuation": chain["authorized_self_record"],
        "continuation_verifier": chain["authorized_verifier_record"],
        "runner_gate_replay": chain["authorized_replay_record"],
    }
    expected_historical_identity = {
        "historical_link_count_observation": "NOT_RECORDED_NOT_INFERRED",
        "scientific_payload_changed": False,
        "claim_ceiling_changed": False,
    }
    expected_signed_promotion = {
        "authorization": chain["authorization_record"],
        "prepare_exit_attestation": chain["prepare_exit_attestation_record"],
        "authorization_signature": chain["authorization_signature_record"],
        "completion": chain["completion_record"],
        "promote_exit_attestation": chain["promote_exit_attestation_record"],
        "completion_signature": chain["completion_signature_record"],
        "historical_source_receipt": chain["source_record"],
        "canonical_source_receipt": chain["canonical_record"],
        "authorization_signature_verified": True,
        "completion_signature_verified": True,
        "bound_relation_count": chain["bound_relation_count"],
    }
    if (
        not strict_json_equal(payload.get("code"), expected_code)
        or not strict_json_equal(
            payload.get("legacy_signed_code"), expected_legacy_signed_code
        )
        or not strict_json_equal(
            payload.get("schema_recovery_authority"),
            chain["schema_recovery_authority"],
        )
        or not strict_json_equal(
            payload.get("historical_identity_contract"),
            expected_historical_identity,
        )
        or not strict_json_equal(payload.get("self_binding"), chain["self_binding"])
        or not strict_json_equal(
            payload.get("signed_promotion"), expected_signed_promotion
        )
        or not strict_json_equal(
            payload.get("authority_semantics"), TERMINAL_AUTHORITY_SEMANTICS
        )
        or not strict_json_equal(payload.get("governance"), TERMINAL_GOVERNANCE)
        or not strict_json_equal(payload.get("checks"), TERMINAL_CHECKS)
    ):
        raise ContinuationError("Continuation execution receipt authority semantics drift")

    verification = load_recorded_verification(reader, chain)
    validate_recorded_fresh_verification_claim(
        payload.get("fresh_promotion_verification"), verification, chain
    )
    replay = validate_replay_evidence(reader, chain, verification)
    if not strict_json_equal(payload.get("runner_only_replay"), replay):
        raise ContinuationError("Continuation runner-only replay relation drift")

    gate_inputs = StreamBoundArtifactSet()
    gate_validation = validate_gate_inputs(
        gate_inputs, chain["reviewed_runtime_sources"]
    )
    if not strict_json_equal(payload.get("parent_gate_validation"), gate_validation):
        raise ContinuationError("Continuation parent gate validation drift")
    child_runtime_fds = validate_child_execution_claim(
        payload.get("child_execution"), reader, chain
    )
    _, script_contract = compose_downstream_script(
        chain["runner_data"],
        primary_python_fd=child_runtime_fds["primary_python"],
        qa_python_fd=child_runtime_fds["qa_python"],
        primary_python_record=chain["python_record"],
        qa_python_record=chain["qa_python_record"],
    )
    if not strict_json_equal(payload.get("composed_runner_stdin"), script_contract):
        raise ContinuationError("Continuation composed runner contract drift")

    terminal, require_terminal_release_stable = validate_terminal_releases(
        authority_module,
        reader,
        v7_validation,
        chain["reviewed_runtime_sources"],
        script_contract["runtime_fd_canonical_argv0_wrappers"]["primary_python"][
            "canonical_argv0"
        ],
        None,
        None,
    )
    if not strict_json_equal(payload.get("terminal_releases"), terminal):
        raise ContinuationError("Continuation terminal release inventory drift")
    expected_v7 = {
        "authority_id": v7_validation["authority_id"],
        "source_set_sha256": v7_validation["source_set_sha256"],
        "validation_binding_lease": v7_validation["validation_binding_lease"],
        "pass": True,
    }
    if not strict_json_equal(payload.get("v7_validation"), expected_v7):
        raise ContinuationError("Continuation v7 validation relation drift")

    protected_paths = collect_protected_paths(
        chain["authorization"],
        chain["completion"],
        verification["recorded_payload"],
        verification["recorded_payload"],
        gate_validation["records"],
    )
    protected_paths.add(FAILED_V16_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V17_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V18_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V19_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(REJECTED_V20_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V21_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(REJECTED_V22_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(REJECTED_V23_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(REJECTED_V24_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V25_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V26_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(REJECTED_V27_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V28_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(FAILED_V29_CONTINUATION_PUBLIC_KEY)
    protected_paths.add(RECOVERY_CONTINUATION_PUBLIC_KEY)
    mutation_sentinel = PathMutationSentinel(
        protected_paths,
        all_event_directories={PORTABLE_ASSET_ROOT, TUMOR_REF},
    )
    expected_mutation_monitor = {
        "mechanism": (
            "linux_inotify_parent_watch_before_leaf_snapshot_then_file_watch_"
            "with_mode000_parent_fallback_and_setup_identity_recheck"
        ),
        "all_event_directories": sorted(
            str(path) for path in (PORTABLE_ASSET_ROOT, TUMOR_REF)
        ),
        "protected_path_count": mutation_sentinel.protected_path_count,
        "protected_path_set_sha256": mutation_sentinel.protected_path_set_sha256,
        "direct_file_watch_count": mutation_sentinel.direct_file_watch_count,
        "parent_only_mode000_path_count": (
            mutation_sentinel.parent_only_mode000_path_count
        ),
        "parent_only_mode000_path_set_sha256": (
            mutation_sentinel.parent_only_mode000_path_set_sha256
        ),
        "parent_only_fallback_reason": (
            "inotify_eacces_on_exact_mode000_canonical_regular_file"
        ),
        "setup_order": (
            "named_parent_and_all_event_watches_before_leaf_snapshot_then_"
            "direct_file_watch_then_full_stat_identity_recheck"
        ),
        "setup_identity_recheck_count": (
            mutation_sentinel.setup_identity_recheck_count
        ),
        "setup_identity_recheck_path_set_sha256": (
            mutation_sentinel.setup_identity_recheck_path_set_sha256
        ),
        "transient_mutation_events_observed": 0,
        "retained_until_process_exit": True,
        "pass": True,
    }
    if not strict_json_equal(payload.get("mutation_monitor"), expected_mutation_monitor):
        raise ContinuationError("Continuation mutation-monitor contract drift")

    signature_contract = payload.get("terminal_signature_contract")
    expected_signature_keys = {
        "algorithm",
        "public_key",
        "private_key_pre_signature",
        "execution_receipt",
        "signed_artifact",
        "signature",
        "supervisor_success_witness",
        "required_private_key_post_signature_mode",
        "unsigned_receipt_grants_release_authority",
        "success_witness_required_for_release_authority",
    }
    private_pre = (
        signature_contract.get("private_key_pre_signature")
        if isinstance(signature_contract, Mapping)
        else None
    )
    private_live_stat = chain["recovery_continuation_private_key_stat"]
    if (
        not isinstance(signature_contract, Mapping)
        or set(signature_contract) != expected_signature_keys
        or signature_contract.get("algorithm") != "Ed25519"
        or not strict_json_equal(
            signature_contract.get("public_key"),
            chain["recovery_continuation_public_key_record"],
        )
        or signature_contract.get("execution_receipt") != str(CONTINUATION_RECEIPT)
        or signature_contract.get("signed_artifact")
        != str(CONTINUATION_EXIT_ATTESTATION)
        or signature_contract.get("signature") != str(CONTINUATION_SIGNATURE)
        or signature_contract.get("supervisor_success_witness")
        != str(CONTINUATION_SUCCESS_WITNESS)
        or signature_contract.get("required_private_key_post_signature_mode") != "0o0"
        or signature_contract.get("unsigned_receipt_grants_release_authority") is not False
        or signature_contract.get("success_witness_required_for_release_authority")
        is not True
        or not isinstance(private_pre, Mapping)
        or set(private_pre) != set(PRIVATE_KEY_PRE_SIGNATURE_KEYS)
        or private_pre.get("path") != str(RECOVERY_CONTINUATION_PRIVATE_KEY)
        or private_pre.get("mode") != "0o400"
        or type(private_pre.get("device")) is not int
        or type(private_pre.get("inode")) is not int
        or type(private_pre.get("link_count")) is not int
        or type(private_pre.get("mtime_ns")) is not int
        or type(private_pre.get("ctime_ns")) is not int
        or private_pre.get("device") != private_live_stat.st_dev
        or private_pre.get("inode") != private_live_stat.st_ino
        or private_pre.get("link_count") != private_live_stat.st_nlink
        or private_pre.get("mtime_ns") != private_live_stat.st_mtime_ns
        or private_live_stat.st_ctime_ns < private_pre["ctime_ns"]
    ):
        raise ContinuationError("Continuation execution-receipt signature contract drift")

    def require_execution_stable() -> None:
        reader.require_paths_still_bound()
        gate_inputs.require_paths_still_bound()
        require_terminal_release_stable()
        require_bootstrap_stable()
        require_incident_absent()
        mutation_sentinel.assert_clean()

    require_execution_stable()
    gate_inputs.retain_until_process_exit()
    mutation_sentinel.retain_until_process_exit()
    return {
        "payload": payload,
        "receipt": with_stat(receipt_record, receipt_stat),
        "complete_child_output_inventory": terminal[
            "complete_child_output_inventory"
        ],
        "private_key_pre_signature": dict(private_pre),
        "supervisor_capability": supervisor_capability,
        "terminal_releases": terminal,
        "require_execution_stable": require_execution_stable,
    }


def validate_exit_attestation(
    *,
    data: bytes,
    chain: Mapping[str, Any],
    execution_context: Mapping[str, Any],
) -> dict[str, Any]:
    attestation = strict_json(data, "continuation supervisor exit attestation")
    require_exact_keys(
        attestation,
        CONTINUATION_EXIT_ATTESTATION_KEYS,
        "continuation supervisor exit attestation",
    )
    require_utc_timestamp(
        attestation.get("created_at_utc"), "continuation exit-attestation timestamp"
    )
    child_wait = attestation.get("child_wait")
    signature_contract = attestation.get("terminal_signature_contract")
    supervisor_capability = validate_supervisor_capability_record(
        attestation.get("supervisor_capability"), chain
    )
    expected_wait_keys = {
        "pid",
        "returncode",
        "exited_normally",
        "terminating_signal",
        "stdout",
        "stderr",
    }
    expected_signature_keys = {
        "algorithm",
        "public_key",
        "private_key_pre_signature",
        "signed_artifact",
        "signature",
        "supervisor_success_witness",
        "required_private_key_post_signature_mode",
        "signature_alone_grants_release_authority",
        "success_witness_required_for_release_authority",
    }
    if (
        attestation.get("schema_name")
        != f"{CONTINUATION_SCHEMA}.supervisor_exit_attestation"
        or attestation.get("schema_version") != "1.0.0"
        or attestation.get("task_type") != "B_comprehensive_validation"
        or attestation.get("status")
        != "SUPERVISED_CHILD_NORMAL_EXIT_ZERO_ATTESTED_AWAITING_SIGNATURE_"
        "AND_SUCCESS_WITNESS"
        or attestation.get("scope")
        != "all_7_datasets_469849_sites_signed_task_b_release_chain"
        or not strict_json_equal(
            attestation.get("supervisor_command"), canonical_command()
        )
        or not strict_json_equal(
            attestation.get("supervised_child_command"), supervised_child_command()
        )
        or not strict_json_equal(
            attestation.get("clean_environment"), EXPECTED_ENVIRONMENT
        )
        or not strict_json_equal(
            attestation.get("supervisor_source"), chain["self_record"]
        )
        or not strict_json_equal(
            attestation.get("supervisor_source_binding"), chain["self_binding"]
        )
        or not strict_json_equal(
            attestation.get("python_runtime"), chain["python_record"]
        )
        or not strict_json_equal(
            supervisor_capability, execution_context["supervisor_capability"]
        )
        or not strict_json_equal(
            attestation.get("execution_receipt"), execution_context["receipt"]
        )
        or attestation.get("execution_receipt_pass_semantics")
        != TERMINAL_PASS_SEMANTICS
        or not strict_json_equal(attestation.get("checks"), EXIT_ATTESTATION_CHECKS)
        or attestation.get("pass") is not True
        or attestation.get("pass_semantics") != EXIT_ATTESTATION_PASS_SEMANTICS
        or not isinstance(child_wait, Mapping)
        or set(child_wait) != expected_wait_keys
        or type(child_wait.get("pid")) is not int
        or child_wait["pid"] <= 1
        or type(child_wait.get("returncode")) is not int
        or child_wait.get("returncode") != 0
        or child_wait.get("exited_normally") is not True
        or child_wait.get("terminating_signal") is not None
    ):
        raise ContinuationError("Continuation supervisor exit-attestation drift")
    for stream_name in ("stdout", "stderr"):
        stream = child_wait.get(stream_name)
        if (
            not isinstance(stream, Mapping)
            or set(stream) != {"size_bytes", "sha256"}
            or type(stream.get("size_bytes")) is not int
            or stream["size_bytes"] < 0
            or not is_canonical_sha256(stream.get("sha256"))
        ):
            raise ContinuationError(
                f"Continuation supervised-child {stream_name} record drift"
            )
    if (
        not isinstance(signature_contract, Mapping)
        or set(signature_contract) != expected_signature_keys
        or signature_contract.get("algorithm") != "Ed25519"
        or not strict_json_equal(
            signature_contract.get("public_key"),
            chain["recovery_continuation_public_key_record"],
        )
        or not strict_json_equal(
            signature_contract.get("private_key_pre_signature"),
            execution_context["private_key_pre_signature"],
        )
        or signature_contract.get("signed_artifact")
        != str(CONTINUATION_EXIT_ATTESTATION)
        or signature_contract.get("signature") != str(CONTINUATION_SIGNATURE)
        or signature_contract.get("supervisor_success_witness")
        != str(CONTINUATION_SUCCESS_WITNESS)
        or signature_contract.get("required_private_key_post_signature_mode") != "0o0"
        or signature_contract.get("signature_alone_grants_release_authority") is not False
        or signature_contract.get("success_witness_required_for_release_authority")
        is not True
    ):
        raise ContinuationError("Continuation exit-attestation signature contract drift")
    return attestation


def hash_open_descriptor(descriptor: int) -> dict[str, Any]:
    before = os.fstat(descriptor)
    digest = hashlib.sha256()
    offset = 0
    while offset < before.st_size:
        block = os.pread(
            descriptor,
            min(8 * 1024 * 1024, before.st_size - offset),
            offset,
        )
        if not block:
            raise ContinuationError("Short read while hashing supervised-child stream")
        digest.update(block)
        offset += len(block)
    after = os.fstat(descriptor)
    if not same_stat(before, after):
        raise ContinuationError("Supervised-child stream changed while hashing")
    return {"size_bytes": before.st_size, "sha256": digest.hexdigest()}


def stage_private_key_in_sealed_memfd(private_fd: int) -> int:
    source_before = os.fstat(private_fd)
    if (
        not stat.S_ISREG(source_before.st_mode)
        or stat.S_IMODE(source_before.st_mode) != 0o400
        or source_before.st_nlink != 1
        or source_before.st_size <= 0
        or source_before.st_size > 64 * 1024
    ):
        raise ContinuationError("Continuation private-key source contract drift")
    key_data = b"".join(
        os.pread(
            private_fd,
            min(64 * 1024, source_before.st_size - position),
            position,
        )
        for position in range(0, source_before.st_size, 64 * 1024)
    )
    source_after = os.fstat(private_fd)
    if len(key_data) != source_before.st_size or not same_stat(
        source_before, source_after
    ):
        raise ContinuationError("Continuation private key changed while staging")

    descriptor = linux_memfd_create("intersubmod-continuation-signing-key")
    try:
        offset = 0
        while offset < len(key_data):
            written = os.write(descriptor, key_data[offset:])
            if written <= 0:
                raise ContinuationError("Short write to continuation signing-key memfd")
            offset += written
        os.fchmod(descriptor, 0o400)
        os.fsync(descriptor)
        seal_mask = supervisor_capability_seal_mask()
        fcntl.fcntl(descriptor, LINUX_F_ADD_SEALS, seal_mask)
        staged = os.fstat(descriptor)
        if (
            not stat.S_ISREG(staged.st_mode)
            or stat.S_IMODE(staged.st_mode) != 0o400
            or staged.st_nlink != 0
            or staged.st_size != len(key_data)
            or fcntl.fcntl(descriptor, LINUX_F_GET_SEALS) != seal_mask
            or os.pread(descriptor, staged.st_size, 0) != key_data
        ):
            raise ContinuationError("Sealed continuation signing-key memfd drift")
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def sign_exit_attestation(
    authority_module: Any,
    chain: Mapping[str, Any],
    execution_context: Mapping[str, Any],
    committed_attestation: Mapping[str, Any],
) -> dict[str, Any]:
    require_absent(CONTINUATION_SIGNATURE, "Continuation terminal signature")
    require_absent(
        CONTINUATION_SUCCESS_WITNESS, "Continuation supervisor-success witness"
    )
    require_incident_absent()
    if (
        not isinstance(committed_attestation, Mapping)
        or set(committed_attestation) != {"descriptor", "data", "artifact"}
        or type(committed_attestation.get("descriptor")) is not int
        or committed_attestation["descriptor"] <= 2
        or not isinstance(committed_attestation.get("data"), bytes)
        or not isinstance(committed_attestation.get("artifact"), Mapping)
        or set(committed_attestation["artifact"]) != set(STAT_ARTIFACT_KEYS)
    ):
        raise ContinuationError("Committed exit-attestation handle contract drift")
    flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0)
    attestation_fd = os.dup(int(committed_attestation["descriptor"]))
    public_fd = os.open(RECOVERY_CONTINUATION_PUBLIC_KEY, flags)
    private_fd = os.open(RECOVERY_CONTINUATION_PRIVATE_KEY, flags)
    parent_flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0)
    if hasattr(os, "O_NOFOLLOW"):
        parent_flags |= os.O_NOFOLLOW
    signature_parent_fd = os.open(CONTINUATION_SIGNATURE.parent, parent_flags)
    key_parent_fd = os.open(RECOVERY_CONTINUATION_PRIVATE_KEY.parent, parent_flags)
    signature_fd = -1
    private_signing_fd = -1
    try:
        attestation_stat = os.fstat(attestation_fd)
        attestation_data = b"".join(
            os.pread(
                attestation_fd,
                min(1024 * 1024, attestation_stat.st_size - position),
                position,
            )
            for position in range(0, attestation_stat.st_size, 1024 * 1024)
        )
        attestation_live_fd = os.open(CONTINUATION_EXIT_ATTESTATION, flags)
        try:
            attestation_live_stat = os.fstat(attestation_live_fd)
            attestation_live_data = b"".join(
                os.pread(
                    attestation_live_fd,
                    min(1024 * 1024, attestation_live_stat.st_size - position),
                    position,
                )
                for position in range(
                    0, attestation_live_stat.st_size, 1024 * 1024
                )
            )
        finally:
            os.close(attestation_live_fd)
        public_stat = os.fstat(public_fd)
        private_stat = os.fstat(private_fd)
        if (
            stat.S_IMODE(attestation_stat.st_mode) != 0o444
            or attestation_stat.st_nlink != 1
            or not same_stat(attestation_stat, attestation_live_stat)
            or attestation_data != committed_attestation["data"]
            or attestation_live_data != committed_attestation["data"]
            or with_stat(
                {
                    "path": str(CONTINUATION_EXIT_ATTESTATION),
                    "size_bytes": len(attestation_data),
                    "sha256": sha256_bytes(attestation_data),
                    "mode": "0o444",
                },
                attestation_stat,
            )
            != committed_attestation["artifact"]
            or stat.S_IMODE(public_stat.st_mode) != 0o444
            or public_stat.st_nlink != 1
            or stat.S_IMODE(private_stat.st_mode) != 0o400
            or private_stat.st_nlink != 1
            or with_stat(
                {
                    "path": str(RECOVERY_CONTINUATION_PRIVATE_KEY),
                    "mode": "0o400",
                },
                private_stat,
            )
            != execution_context["private_key_pre_signature"]
        ):
            raise ContinuationError("Continuation key state drift before signing")

        private_signing_fd = stage_private_key_in_sealed_memfd(private_fd)

        # Retire the on-disk key before any signature bytes can exist. Signing
        # uses only the sealed, unlinked memfd and closes it before publication.
        os.fchmod(private_fd, 0o000)
        os.fsync(private_fd)
        os.fsync(key_parent_fd)
        retired = os.fstat(private_fd)
        private_live = os.lstat(RECOVERY_CONTINUATION_PRIVATE_KEY)
        if (
            stat.S_IMODE(retired.st_mode) != 0
            or retired.st_nlink != 1
            or not same_stat(retired, private_live)
            or RECOVERY_CONTINUATION_PRIVATE_KEY.resolve(strict=True)
            != RECOVERY_CONTINUATION_PRIVATE_KEY
        ):
            raise ContinuationError("Continuation private-key retirement drift")

        result = subprocess.run(
            [
                f"/proc/self/fd/{chain['openssl_fd']}",
                "pkeyutl",
                "-sign",
                "-rawin",
                "-inkey",
                f"/proc/self/fd/{private_signing_fd}",
                "-in",
                f"/proc/self/fd/{attestation_fd}",
            ],
            cwd=REPO_ROOT,
            env=EXPECTED_ENVIRONMENT,
            pass_fds=(
                int(chain["openssl_fd"]),
                private_signing_fd,
                attestation_fd,
            ),
            capture_output=True,
            check=False,
            timeout=120,
        )
        os.close(private_signing_fd)
        private_signing_fd = -1
        if result.returncode != 0 or result.stderr or len(result.stdout) != 64:
            raise ContinuationError(
                "Continuation attestation signing failed: "
                f"exit={result.returncode} stderr_sha256={sha256_bytes(result.stderr)}"
            )
        if not hasattr(os, "O_TMPFILE"):
            raise ContinuationError("O_TMPFILE is required for terminal signature")
        signature_fd = os.open(
            ".",
            os.O_RDWR | os.O_TMPFILE | os.O_CLOEXEC,
            0o400,
            dir_fd=signature_parent_fd,
        )
        if os.write(signature_fd, result.stdout) != 64:
            raise ContinuationError("Short write to continuation signature staging inode")
        os.fchmod(signature_fd, 0o444)
        os.fsync(signature_fd)
        staged = os.fstat(signature_fd)
        if (
            staged.st_size != 64
            or staged.st_nlink != 0
            or stat.S_IMODE(staged.st_mode) != 0o444
            or os.pread(signature_fd, 64, 0) != result.stdout
            or not authority_module.verify_ed25519_signature_fds(
                openssl_fd=int(chain["openssl_fd"]),
                data_fd=attestation_fd,
                public_key_fd=public_fd,
                signature_fd=signature_fd,
            )
        ):
            raise ContinuationError("Staged continuation signature verification failed")
        require_absent(CONTINUATION_SIGNATURE, "Continuation terminal signature")
        link_fd_no_replace(
            signature_fd,
            signature_parent_fd,
            CONTINUATION_SIGNATURE.name,
        )
        # The no-replace link is a provisional commit. Everything after it is
        # observational; any raised failure is incident-recorded by the lease holder.
        os.fsync(signature_parent_fd)
        signature_live_fd = os.open(
            CONTINUATION_SIGNATURE.name,
            flags,
            dir_fd=signature_parent_fd,
        )
        try:
            signature_live_stat = os.fstat(signature_live_fd)
            if (
                signature_live_stat.st_size != 64
                or signature_live_stat.st_nlink != 1
                or stat.S_IMODE(signature_live_stat.st_mode) != 0o444
                or not same_stat(os.fstat(signature_fd), signature_live_stat)
                or os.pread(signature_live_fd, 64, 0) != result.stdout
                or not authority_module.verify_ed25519_signature_fds(
                    openssl_fd=int(chain["openssl_fd"]),
                    data_fd=attestation_fd,
                    public_key_fd=public_fd,
                    signature_fd=signature_live_fd,
                )
            ):
                raise ContinuationError("Published continuation signature drift")
        finally:
            os.close(signature_live_fd)
        require_absent(
            CONTINUATION_SUCCESS_WITNESS,
            "Continuation supervisor-success witness before preverification",
        )

        return {
            "signature": {
                "path": str(CONTINUATION_SIGNATURE),
                "size_bytes": 64,
                "sha256": sha256_bytes(result.stdout),
                "mode": "0o444",
            },
            "retired_private_key": {
                "path": str(RECOVERY_CONTINUATION_PRIVATE_KEY),
                "mode": "0o0",
            },
            "signature_status": "PROVISIONAL_PENDING_INDEPENDENT_VERIFICATION",
            "release_authority_granted": False,
            "pass": True,
        }
    finally:
        if private_signing_fd >= 0:
            os.close(private_signing_fd)
        if signature_fd >= 0:
            os.close(signature_fd)
        os.close(key_parent_fd)
        os.close(signature_parent_fd)
        os.close(private_fd)
        os.close(public_fd)
        os.close(attestation_fd)


def supervise_and_sign() -> None:
    require_clean_runtime(mode="supervisor")
    capability_reservation = reserve_supervisor_capability_descriptor()
    require_incident_absent()
    for path in DOWNSTREAM_OUTPUT_SLOTS:
        require_absent(path, "Downstream output before supervised continuation")

    authority_module, bootstrap_record = bootstrap_release_authority()
    v7_validation = validate_v7(authority_module)
    reader = authority_module.BoundArtifactReader()
    chain = validate_promotion_chain(
        authority_module,
        reader,
        v7_validation,
        bootstrap_record,
        continuation_private_live_mode="0o400",
    )
    capability_fd, supervisor_capability = create_supervisor_capability(
        chain, capability_reservation
    )
    supervisor_sentinel = PathMutationSentinel(
        {
            CONTINUATION_RUNNER,
            PYTHON.resolve(strict=True),
            OPENSSL,
            RECOVERY_CONTINUATION_PUBLIC_KEY,
            *(contract["path"] for contract in REVIEWED_RUNTIME_SOURCE_CONTRACTS.values()),
        },
        all_event_directories={PORTABLE_ASSET_ROOT, TUMOR_REF},
    )
    capture_parent_flags = (
        os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0)
    )
    if hasattr(os, "O_NOFOLLOW"):
        capture_parent_flags |= os.O_NOFOLLOW
    capture_parent_fd = os.open(RESULT_ROOT, capture_parent_flags)
    stdout_fd = os.open(
        ".", os.O_RDWR | os.O_TMPFILE | os.O_CLOEXEC, 0o400, dir_fd=capture_parent_fd
    )
    stderr_fd = os.open(
        ".", os.O_RDWR | os.O_TMPFILE | os.O_CLOEXEC, 0o400, dir_fd=capture_parent_fd
    )
    try:
        reader.require_paths_still_bound()
        require_bootstrap_stable()
        supervisor_sentinel.assert_clean()
        try:
            process = subprocess.Popen(
                supervised_child_command(),
                executable=f"/proc/self/fd/{chain['python_fd']}",
                cwd=REPO_ROOT,
                env=EXPECTED_ENVIRONMENT,
                pass_fds=(int(chain["python_fd"]), capability_fd),
                stdout=stdout_fd,
                stderr=stderr_fd,
                close_fds=True,
            )
            returncode = process.wait()
        finally:
            os.close(capability_fd)
        os.fsync(stdout_fd)
        os.fsync(stderr_fd)
        stdout_record = hash_open_descriptor(stdout_fd)
        stderr_record = hash_open_descriptor(stderr_fd)
        if returncode != 0:
            raise ContinuationError(
                "Supervised continuation child did not exit zero: "
                f"returncode={returncode} stderr_sha256={stderr_record['sha256']}"
            )
        reader.require_paths_still_bound()
        require_bootstrap_stable()
        require_incident_absent()
        supervisor_sentinel.assert_clean()
        require_absent(
            CONTINUATION_EXIT_ATTESTATION,
            "Continuation supervisor exit attestation",
        )
        require_absent(CONTINUATION_SIGNATURE, "Continuation terminal signature")
        receipt_fd, receipt_data, receipt_record, receipt_stat = open_bound(
            reader, CONTINUATION_RECEIPT
        )
        del receipt_fd
        execution_context = validate_execution_receipt(
            authority_module=authority_module,
            bootstrap_record=bootstrap_record,
            v7_validation=v7_validation,
            reader=reader,
            chain=chain,
            receipt_data=receipt_data,
            receipt_record=receipt_record,
            receipt_stat=receipt_stat,
            expected_supervisor_capability=supervisor_capability,
        )
        child_wait = {
            "pid": process.pid,
            "returncode": returncode,
            "exited_normally": returncode == 0,
            "terminating_signal": None if returncode >= 0 else -returncode,
            "stdout": stdout_record,
            "stderr": stderr_record,
        }
        attestation = {
            "schema_name": f"{CONTINUATION_SCHEMA}.supervisor_exit_attestation",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "task_type": "B_comprehensive_validation",
            "status": (
                "SUPERVISED_CHILD_NORMAL_EXIT_ZERO_ATTESTED_AWAITING_SIGNATURE_"
                "AND_SUCCESS_WITNESS"
            ),
            "scope": "all_7_datasets_469849_sites_signed_task_b_release_chain",
            "supervisor_command": canonical_command(),
            "supervised_child_command": supervised_child_command(),
            "clean_environment": EXPECTED_ENVIRONMENT,
            "supervisor_source": chain["self_record"],
            "supervisor_source_binding": chain["self_binding"],
            "python_runtime": chain["python_record"],
            "supervisor_capability": supervisor_capability,
            "child_wait": child_wait,
            "execution_receipt": execution_context["receipt"],
            "execution_receipt_pass_semantics": TERMINAL_PASS_SEMANTICS,
            "terminal_signature_contract": {
                "algorithm": "Ed25519",
                "public_key": chain["recovery_continuation_public_key_record"],
                "private_key_pre_signature": execution_context[
                    "private_key_pre_signature"
                ],
                "signed_artifact": str(CONTINUATION_EXIT_ATTESTATION),
                "signature": str(CONTINUATION_SIGNATURE),
                "supervisor_success_witness": str(CONTINUATION_SUCCESS_WITNESS),
                "required_private_key_post_signature_mode": "0o0",
                "signature_alone_grants_release_authority": False,
                "success_witness_required_for_release_authority": True,
            },
            "checks": dict(EXIT_ATTESTATION_CHECKS),
            "pass": True,
            "pass_semantics": EXIT_ATTESTATION_PASS_SEMANTICS,
        }

        def attestation_commit_gate() -> None:
            execution_context["require_execution_stable"]()
            reader.require_paths_still_bound()
            require_bootstrap_stable()
            require_incident_absent()
            supervisor_sentinel.assert_clean()
            require_absent(CONTINUATION_SIGNATURE, "Continuation terminal signature")
            require_absent(
                CONTINUATION_SUCCESS_WITNESS,
                "Continuation supervisor-success witness",
            )

        committed_attestation = commit_terminal_receipt(
            attestation,
            attestation_commit_gate,
            path=CONTINUATION_EXIT_ATTESTATION,
            postcommit=attestation_commit_gate,
        )
        signature_result = sign_exit_attestation(
            authority_module,
            chain,
            execution_context,
            committed_attestation,
        )
    finally:
        os.close(stderr_fd)
        os.close(stdout_fd)
        os.close(capture_parent_fd)

    supervisor_sentinel.assert_clean()
    if not same_stat(
        os.fstat(int(chain["self_fd"])),
        os.stat(CONTINUATION_RUNNER, follow_symlinks=False),
    ):
        raise ContinuationError("Continuation source FD/path binding drift after signing")
    signed_verifier_argv = [
        str(PYTHON),
        "-I",
        "-B",
        f"/proc/self/fd/{chain['self_fd']}",
        "--verify-signed-terminal-prewitness",
    ]
    result = subprocess.run(
        signed_verifier_argv,
        executable=f"/proc/self/fd/{chain['python_fd']}",
        cwd=REPO_ROOT,
        env=EXPECTED_ENVIRONMENT,
        pass_fds=(int(chain["python_fd"]), int(chain["self_fd"])),
        capture_output=True,
        check=False,
        timeout=1800,
    )
    if result.returncode != 0 or result.stderr:
        raise ContinuationError(
            "Fresh signed-terminal prewitness verifier failed after signing: "
            f"exit={result.returncode} stderr_sha256={sha256_bytes(result.stderr)}"
        )
    supervisor_sentinel.assert_clean()
    verified = strict_json(
        result.stdout, "fresh signed-terminal prewitness verifier stdout"
    )
    expected_keys = {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "canonical_command",
        "continuation_source_binding",
        "continuation_receipt",
        "continuation_exit_attestation",
        "supervised_child_wait",
        "supervisor_capability",
        "threat_model",
        "continuation_signature",
        "continuation_public_key",
        "retired_continuation_private_key",
        "complete_child_output_inventory",
        "nested_dataset_release_signature_verified",
        "nested_report_release_signature_verified",
        "continuation_terminal_signature_verified",
        "supervised_child_normal_exit_zero_verified",
        "execution_receipt_full_nested_contract_verified",
        "success_witness_required",
        "ready_for_success_witness",
        "release_authority_granted",
        "pass",
        "pass_semantics",
    }
    require_exact_keys(
        verified,
        expected_keys,
        "fresh signed-terminal prewitness verification result",
    )
    require_utc_timestamp(
        verified.get("created_at_utc"),
        "fresh signed-terminal prewitness verification timestamp",
    )
    if (
        verified.get("schema_name")
        != f"{CONTINUATION_SCHEMA}.signed_terminal_prewitness_verification"
        or verified.get("schema_version") != "1.0.0"
        or not strict_json_equal(
            verified.get("canonical_command"),
            signed_terminal_prewitness_verification_command(),
        )
        or not strict_json_equal(
            verified.get("continuation_source_binding"), chain["self_binding"]
        )
        or not strict_json_equal(
            verified.get("continuation_receipt"), execution_context["receipt"]
        )
        or not strict_json_equal(
            verified.get("continuation_exit_attestation"),
            committed_attestation["artifact"],
        )
        or not strict_json_equal(verified.get("supervised_child_wait"), child_wait)
        or not strict_json_equal(
            verified.get("supervisor_capability"), supervisor_capability
        )
        or not strict_json_equal(
            verified.get("threat_model"), CONTINUATION_THREAT_MODEL
        )
        or not strict_json_equal(
            verified.get("continuation_signature"), signature_result["signature"]
        )
        or not strict_json_equal(
            verified.get("continuation_public_key"),
            chain["recovery_continuation_public_key_record"],
        )
        or not strict_json_equal(
            verified.get("retired_continuation_private_key"),
            signature_result["retired_private_key"],
        )
        or not strict_json_equal(
            verified.get("complete_child_output_inventory"),
            execution_context["complete_child_output_inventory"],
        )
        or verified.get("nested_dataset_release_signature_verified") is not True
        or verified.get("nested_report_release_signature_verified") is not True
        or verified.get("continuation_terminal_signature_verified") is not True
        or verified.get("success_witness_required") is not True
        or verified.get("ready_for_success_witness") is not True
        or verified.get("release_authority_granted") is not False
        or verified.get("supervised_child_normal_exit_zero_verified") is not True
        or verified.get("execution_receipt_full_nested_contract_verified") is not True
        or verified.get("pass") is not True
        or verified.get("pass_semantics")
        != PREWITNESS_VERIFICATION_PASS_SEMANTICS
    ):
        raise ContinuationError(
            "Fresh signed-terminal prewitness verifier result contract drift"
        )

    _, _, signature_record, signature_stat = open_bound(
        reader, CONTINUATION_SIGNATURE
    )
    signature_artifact = with_stat(signature_record, signature_stat)
    if not strict_json_equal(signature_record, signature_result["signature"]):
        raise ContinuationError("Supervisor signature result/path binding drift")
    preverification_record = {
        "command": signed_terminal_prewitness_verification_command(),
        "exit_code": 0,
        "stdout": {
            "size_bytes": len(result.stdout),
            "sha256": sha256_bytes(result.stdout),
        },
        "stderr": {
            "size_bytes": len(result.stderr),
            "sha256": sha256_bytes(result.stderr),
        },
        "schema_name": verified["schema_name"],
        "ready_for_success_witness": verified["ready_for_success_witness"],
        "release_authority_granted": verified["release_authority_granted"],
        "pass": verified["pass"],
        "pass_semantics": verified["pass_semantics"],
    }
    witness = {
        "schema_name": f"{CONTINUATION_SCHEMA}.supervisor_success_witness",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "task_type": "B_comprehensive_validation",
        "status": "SIGNED_TERMINAL_PREWITNESS_VERIFIED_DURABLE_SUCCESS_COMMITTED",
        "scope": "all_7_datasets_469849_sites_signed_task_b_release_chain",
        "supervisor_command": canonical_command(),
        "prewitness_verifier_command": (
            signed_terminal_prewitness_verification_command()
        ),
        "continuation_source": chain["self_record"],
        "continuation_source_binding": chain["self_binding"],
        "execution_receipt": execution_context["receipt"],
        "exit_attestation": committed_attestation["artifact"],
        "continuation_signature": signature_artifact,
        "continuation_public_key": chain["recovery_continuation_public_key_record"],
        "retired_continuation_private_key": signature_result[
            "retired_private_key"
        ],
        "supervisor_capability": supervisor_capability,
        "supervised_child_wait": child_wait,
        "complete_child_output_inventory": execution_context[
            "complete_child_output_inventory"
        ],
        "prewitness_verification": preverification_record,
        "checks": dict(SUCCESS_WITNESS_CHECKS),
        "release_authority_granted": False,
        "pass": True,
        "pass_semantics": SUCCESS_WITNESS_PASS_SEMANTICS,
    }

    def witness_inputs_stable() -> None:
        execution_context["require_execution_stable"]()
        reader.require_paths_still_bound()
        require_bootstrap_stable()
        require_incident_absent()
        supervisor_sentinel.assert_clean()
        if not same_stat(
            os.fstat(int(chain["self_fd"])),
            os.stat(CONTINUATION_RUNNER, follow_symlinks=False),
        ):
            raise ContinuationError(
                "Continuation source FD/path drift before success-witness commit"
            )

    def witness_precommit_gate() -> None:
        witness_inputs_stable()
        require_success_witness_state(required=False)

    commit_terminal_success_witness(witness, witness_precommit_gate)


def validate_signed_terminal_precommit_contract(
    *,
    authorization: Mapping[str, Any],
    completion: Mapping[str, Any],
    active_self_record: Mapping[str, Any],
    authorized_self_record: Mapping[str, Any],
    authorized_public_record: Mapping[str, Any],
    failed_v16_public_record: Mapping[str, Any],
    failed_v17_public_record: Mapping[str, Any],
    failed_v18_public_record: Mapping[str, Any],
    failed_v19_public_record: Mapping[str, Any],
    rejected_v20_public_record: Mapping[str, Any],
    failed_v21_public_record: Mapping[str, Any],
    rejected_v22_public_record: Mapping[str, Any],
    rejected_v23_public_record: Mapping[str, Any],
    rejected_v24_public_record: Mapping[str, Any],
    failed_v25_public_record: Mapping[str, Any],
    failed_v26_public_record: Mapping[str, Any],
    rejected_v27_public_record: Mapping[str, Any],
    failed_v28_public_record: Mapping[str, Any],
    failed_v29_public_record: Mapping[str, Any],
    recovery_public_record: Mapping[str, Any],
    canonical_record: Mapping[str, Any],
    source_record: Mapping[str, Any],
    recovery_authority: Mapping[str, Any],
) -> None:
    terminal_gate = authorization.get("continuation_gate")
    authorized_terminal_contract = (
        terminal_gate.get("terminal_signature_contract")
        if isinstance(terminal_gate, Mapping)
        else None
    )
    if (
        authorization.get("schema_name") != AUTHORIZATION_SCHEMA
        or authorization.get("schema_version") != SCHEMA_VERSION
        or authorization.get("authorization_id") != AUTHORIZATION_ID
        or authorization.get("pass") is not True
        or not strict_json_equal(
            authorization.get("reviewed_sources", {}).get("downstream_continuation"),
            authorized_self_record,
        )
        or not strict_json_equal(
            authorization.get("trusted_signing_keys", {}).get(
                "continuation_public_key"
            ),
            authorized_public_record,
        )
        or not isinstance(terminal_gate, Mapping)
        or not strict_json_equal(
            authorized_terminal_contract,
            {
                "algorithm": "Ed25519",
                "public_key": authorized_public_record,
                "private_key": {
                    "path": str(CONTINUATION_PRIVATE_KEY),
                    "mode": "0o400",
                },
                "execution_receipt": str(AUTHORIZED_CONTINUATION_RECEIPT),
                "signed_artifact": str(AUTHORIZED_CONTINUATION_EXIT_ATTESTATION),
                "signature": str(AUTHORIZED_CONTINUATION_SIGNATURE),
                "required_private_key_pre_signature_mode": "0o400",
                "required_private_key_post_signature_mode": "0o0",
                "supervisor_command": authorized_canonical_command(),
                "supervised_child_command": authorized_supervised_child_command(),
                "signed_terminal_verification_command": (
                    authorized_signed_terminal_verification_command()
                ),
            },
        )
        or not strict_json_equal(
            terminal_gate.get("downstream_continuation_command"),
            authorized_canonical_command(),
        )
        or not strict_json_equal(
            terminal_gate.get("supervised_child_command"),
            authorized_supervised_child_command(),
        )
        or terminal_gate.get("supervisor_capability_protocol")
        != SUPERVISOR_CAPABILITY_PROTOCOL
        or terminal_gate.get("threat_model") != CONTINUATION_THREAT_MODEL
        or not strict_json_equal(
            terminal_gate.get("policy"), AUTHORIZED_CONTINUATION_POLICY
        )
        or completion.get("schema_name") != COMPLETION_SCHEMA
        or completion.get("schema_version") != SCHEMA_VERSION
        or completion.get("status")
        != "PROVISIONAL_WRITE_AHEAD_CANONICAL_PUBLICATION_INTENT"
        or completion.get("pass") is not True
        or not strict_json_equal(completion.get("canonical_receipt"), canonical_record)
        or not strict_json_equal(completion.get("source_receipt"), source_record)
    ):
        raise ContinuationError("Legacy terminal key precommit contract drift")
    terminal_key_rotation = recovery_authority.get("terminal_key_rotation")
    terminal_key_state = recovery_authority.get("terminal_key_state")
    if (
        not isinstance(recovery_authority, Mapping)
        or recovery_authority.get("authority_id")
        != "20260724_tumor_ref_schema_recovery_v30"
        or recovery_authority.get("runtime_role") != "downstream_continuation"
        or not strict_json_equal(
            recovery_authority.get("runtime_source"), active_self_record
        )
        or not isinstance(terminal_key_rotation, Mapping)
        or set(terminal_key_rotation) != TERMINAL_KEY_ROTATION_KEYS
        or not strict_json_equal(
            terminal_key_rotation.get("legacy_signed_contract", {}).get(
                "public_key"
            ),
            authorized_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("failed_v16_contract", {}).get(
                "public_key"
            ),
            failed_v16_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("failed_v17_contract", {}).get(
                "public_key"
            ),
            failed_v17_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("failed_v18_contract", {}).get(
                "public_key"
            ),
            failed_v18_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("failed_v19_contract", {}).get(
                "public_key"
            ),
            failed_v19_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("rejected_v20_contract", {}).get(
                "public_key"
            ),
            rejected_v20_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("failed_v21_contract", {}).get(
                "public_key"
            ),
            failed_v21_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("rejected_v22_contract", {}).get(
                "public_key"
            ),
            rejected_v22_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("rejected_v23_contract", {}).get(
                "public_key"
            ),
            rejected_v23_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("rejected_v24_contract", {}).get(
                "public_key"
            ),
            rejected_v24_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("failed_v25_contract", {}).get(
                "public_key"
            ),
            failed_v25_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("failed_v26_contract", {}).get(
                "public_key"
            ),
            failed_v26_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("rejected_v27_contract", {}).get(
                "public_key"
            ),
            rejected_v27_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("failed_v28_contract", {}).get(
                "public_key"
            ),
            failed_v28_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("failed_v29_contract", {}).get(
                "public_key"
            ),
            failed_v29_public_record,
        )
        or not strict_json_equal(
            terminal_key_rotation.get("recovery_v30_contract", {}).get(
                "public_key"
            ),
            recovery_public_record,
        )
        or not isinstance(terminal_key_state, Mapping)
        or set(terminal_key_state) != TERMINAL_KEY_STATE_KEYS
        or not strict_json_equal(
            terminal_key_state.get("legacy_public_key"), authorized_public_record
        )
        or not strict_json_equal(
            terminal_key_state.get("failed_v16_public_key"),
            failed_v16_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("failed_v17_public_key"),
            failed_v17_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("failed_v18_public_key"),
            failed_v18_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("failed_v19_public_key"),
            failed_v19_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("rejected_v20_public_key"),
            rejected_v20_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("failed_v21_public_key"),
            failed_v21_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("rejected_v22_public_key"),
            rejected_v22_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("rejected_v23_public_key"),
            rejected_v23_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("rejected_v24_public_key"),
            rejected_v24_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("failed_v25_public_key"),
            failed_v25_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("failed_v26_public_key"),
            failed_v26_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("rejected_v27_public_key"),
            rejected_v27_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("failed_v28_public_key"),
            failed_v28_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("failed_v29_public_key"),
            failed_v29_public_record,
        )
        or not strict_json_equal(
            terminal_key_state.get("recovery_public_key"), recovery_public_record
        )
        or recovery_authority.get("pass") is not True
    ):
        raise ContinuationError("Recovery authority does not bind active terminal source")


def require_success_witness_state(*, required: bool) -> None:
    present = os.path.lexists(CONTINUATION_SUCCESS_WITNESS)
    if required and not present:
        raise ContinuationError(
            "Durable supervisor-success witness is absent; signature remains provisional"
        )
    if not required and present:
        raise ContinuationError(
            "Supervisor-success witness already exists before prewitness verification"
        )


def validate_supervisor_success_witness(
    *,
    data: bytes,
    witness_record: Mapping[str, Any],
    witness_stat: os.stat_result,
    self_record: Mapping[str, Any],
    self_binding: Mapping[str, Any],
    receipt_artifact: Mapping[str, Any],
    attestation_artifact: Mapping[str, Any],
    signature_artifact: Mapping[str, Any],
    public_record: Mapping[str, Any],
    private_record: Mapping[str, Any],
    execution_context: Mapping[str, Any],
    exit_attestation: Mapping[str, Any],
    child_inventory: Mapping[str, Any],
) -> dict[str, Any]:
    payload = strict_json(data, "durable supervisor-success witness")
    require_exact_keys(
        payload,
        CONTINUATION_SUCCESS_WITNESS_KEYS,
        "durable supervisor-success witness",
    )
    require_utc_timestamp(
        payload.get("created_at_utc"), "supervisor-success witness timestamp"
    )
    preverification = payload.get("prewitness_verification")
    expected_preverification_keys = {
        "command",
        "exit_code",
        "stdout",
        "stderr",
        "schema_name",
        "ready_for_success_witness",
        "release_authority_granted",
        "pass",
        "pass_semantics",
    }
    stdout = (
        preverification.get("stdout")
        if isinstance(preverification, Mapping)
        else None
    )
    stderr = (
        preverification.get("stderr")
        if isinstance(preverification, Mapping)
        else None
    )
    if (
        witness_record.get("mode") != "0o444"
        or witness_stat.st_nlink != 1
        or payload.get("schema_name")
        != f"{CONTINUATION_SCHEMA}.supervisor_success_witness"
        or payload.get("schema_version") != "1.0.0"
        or payload.get("task_type") != "B_comprehensive_validation"
        or payload.get("status")
        != "SIGNED_TERMINAL_PREWITNESS_VERIFIED_DURABLE_SUCCESS_COMMITTED"
        or payload.get("scope")
        != "all_7_datasets_469849_sites_signed_task_b_release_chain"
        or not strict_json_equal(payload.get("supervisor_command"), canonical_command())
        or not strict_json_equal(
            payload.get("prewitness_verifier_command"),
            signed_terminal_prewitness_verification_command(),
        )
        or not strict_json_equal(payload.get("continuation_source"), self_record)
        or not strict_json_equal(
            payload.get("continuation_source_binding"), self_binding
        )
        or not strict_json_equal(payload.get("execution_receipt"), receipt_artifact)
        or not strict_json_equal(payload.get("exit_attestation"), attestation_artifact)
        or not strict_json_equal(
            payload.get("continuation_signature"), signature_artifact
        )
        or not strict_json_equal(
            payload.get("continuation_public_key"), public_record
        )
        or not strict_json_equal(
            payload.get("retired_continuation_private_key"), private_record
        )
        or not strict_json_equal(
            payload.get("supervisor_capability"),
            execution_context["supervisor_capability"],
        )
        or not strict_json_equal(
            payload.get("supervised_child_wait"), exit_attestation["child_wait"]
        )
        or not strict_json_equal(
            payload.get("complete_child_output_inventory"), child_inventory
        )
        or not isinstance(preverification, Mapping)
        or set(preverification) != expected_preverification_keys
        or not strict_json_equal(
            preverification.get("command"),
            signed_terminal_prewitness_verification_command(),
        )
        or preverification.get("exit_code") != 0
        or preverification.get("schema_name")
        != f"{CONTINUATION_SCHEMA}.signed_terminal_prewitness_verification"
        or preverification.get("ready_for_success_witness") is not True
        or preverification.get("release_authority_granted") is not False
        or preverification.get("pass") is not True
        or preverification.get("pass_semantics")
        != PREWITNESS_VERIFICATION_PASS_SEMANTICS
        or not isinstance(stdout, Mapping)
        or set(stdout) != {"size_bytes", "sha256"}
        or type(stdout.get("size_bytes")) is not int
        or stdout["size_bytes"] <= 0
        or not is_canonical_sha256(stdout.get("sha256"))
        or not strict_json_equal(
            stderr, {"size_bytes": 0, "sha256": sha256_bytes(b"")}
        )
        or not strict_json_equal(payload.get("checks"), SUCCESS_WITNESS_CHECKS)
        or payload.get("release_authority_granted") is not False
        or payload.get("pass") is not True
        or payload.get("pass_semantics") != SUCCESS_WITNESS_PASS_SEMANTICS
    ):
        raise ContinuationError("Durable supervisor-success witness contract drift")
    return payload


def verify_signed_terminal(*, require_success_witness: bool = True) -> dict[str, Any]:
    require_clean_runtime(
        mode=(
            "signed_verifier"
            if require_success_witness
            else "signed_prewitness_verifier"
        )
    )
    if require_success_witness:
        _VERIFIER_LEASE_FDS.append(acquire_recovery_output_verifier_lease())
    require_incident_absent()
    require_success_witness_state(required=require_success_witness)
    authority_module, bootstrap_record = bootstrap_release_authority()
    v7_validation = validate_v7(authority_module)
    reader = authority_module.BoundArtifactReader()
    witness_data: bytes | None = None
    witness_record: Mapping[str, Any] | None = None
    witness_stat: os.stat_result | None = None
    if require_success_witness:
        _, witness_data, witness_record, witness_stat = open_bound(
            reader, CONTINUATION_SUCCESS_WITNESS
        )
    chain = validate_promotion_chain(
        authority_module,
        reader,
        v7_validation,
        bootstrap_record,
        continuation_private_live_mode="0o0",
    )

    _, self_data, self_record, self_stat = open_bound(reader, CONTINUATION_RUNNER)
    if self_record["mode"] != "0o444" or self_stat.st_nlink != 1:
        raise ContinuationError("Signed-terminal verifier source is not immutable")
    self_binding = verify_executing_source(self_data)

    openssl_fd, _, openssl_record, _ = open_bound(
        reader, OPENSSL, expected_mode="0o755"
    )
    if openssl_record["sha256"] != EXPECTED_OPENSSL_SHA256:
        raise ContinuationError("Signed-terminal OpenSSL identity drift")
    receipt_fd, receipt_data, receipt_record, receipt_stat = open_bound(
        reader, CONTINUATION_RECEIPT
    )
    attestation_fd, attestation_data, attestation_record, attestation_stat = open_bound(
        reader, CONTINUATION_EXIT_ATTESTATION
    )
    signature_fd, _, signature_record, signature_stat = open_bound(
        reader, CONTINUATION_SIGNATURE
    )
    public_fd, _, public_record, public_stat = open_bound(
        reader, RECOVERY_CONTINUATION_PUBLIC_KEY
    )
    _, private_record, private_stat = open_private_metadata(
        reader, RECOVERY_CONTINUATION_PRIVATE_KEY, "0o0"
    )
    if (
        signature_stat.st_size != 64
        or signature_stat.st_nlink != 1
        or public_stat.st_nlink != 1
        or private_stat.st_nlink != 1
        or public_record["sha256"]
        != EXPECTED_RECOVERY_CONTINUATION_PUBLIC_KEY_SHA256
        or not authority_module.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=attestation_fd,
            public_key_fd=public_fd,
            signature_fd=signature_fd,
        )
    ):
        raise ContinuationError("Detached continuation terminal signature failed")

    payload = strict_json(receipt_data, "signed continuation terminal receipt")
    require_exact_keys(payload, CONTINUATION_RECEIPT_KEYS, "continuation terminal receipt")
    require_utc_timestamp(payload.get("created_at_utc"), "continuation terminal timestamp")
    if (
        payload.get("schema_name") != CONTINUATION_SCHEMA
        or payload.get("schema_version") != "1.0.0"
        or payload.get("task_type") != "B_comprehensive_validation"
        or payload.get("status")
        != "DOWNSTREAM_VALIDATED_AWAITING_SUPERVISOR_EXIT_ATTESTATION_SIGNATURE_"
        "AND_SUCCESS_WITNESS"
        or payload.get("scope")
        != "all_7_datasets_469849_sites_signed_task_b_release_chain"
        or not strict_json_equal(
            payload.get("canonical_command"), supervised_child_command()
        )
        or not strict_json_equal(payload.get("clean_environment"), EXPECTED_ENVIRONMENT)
        or payload.get("pass") is not True
        or payload.get("pass_semantics") != TERMINAL_PASS_SEMANTICS
    ):
        raise ContinuationError("Continuation terminal receipt exact contract drift")

    signature_contract = payload.get("terminal_signature_contract")
    if not isinstance(signature_contract, Mapping) or set(signature_contract) != {
        "algorithm",
        "public_key",
        "private_key_pre_signature",
        "execution_receipt",
        "signed_artifact",
        "signature",
        "supervisor_success_witness",
        "required_private_key_post_signature_mode",
        "unsigned_receipt_grants_release_authority",
        "success_witness_required_for_release_authority",
    }:
        raise ContinuationError("Continuation signature contract schema drift")
    private_pre = signature_contract.get("private_key_pre_signature")
    if (
        signature_contract.get("algorithm") != "Ed25519"
        or not strict_json_equal(signature_contract.get("public_key"), public_record)
        or signature_contract.get("execution_receipt") != str(CONTINUATION_RECEIPT)
        or signature_contract.get("signed_artifact")
        != str(CONTINUATION_EXIT_ATTESTATION)
        or signature_contract.get("signature") != str(CONTINUATION_SIGNATURE)
        or signature_contract.get("supervisor_success_witness")
        != str(CONTINUATION_SUCCESS_WITNESS)
        or signature_contract.get("required_private_key_post_signature_mode") != "0o0"
        or signature_contract.get("unsigned_receipt_grants_release_authority") is not False
        or signature_contract.get("success_witness_required_for_release_authority")
        is not True
        or not isinstance(private_pre, Mapping)
        or set(private_pre) != set(PRIVATE_KEY_PRE_SIGNATURE_KEYS)
        or private_pre.get("path") != str(RECOVERY_CONTINUATION_PRIVATE_KEY)
        or private_pre.get("mode") != "0o400"
        or type(private_pre.get("device")) is not int
        or type(private_pre.get("inode")) is not int
        or type(private_pre.get("link_count")) is not int
        or type(private_pre.get("mtime_ns")) is not int
        or type(private_pre.get("ctime_ns")) is not int
        or private_pre.get("device") != private_stat.st_dev
        or private_pre.get("inode") != private_stat.st_ino
        or private_pre.get("link_count") != private_stat.st_nlink
        or private_pre.get("mtime_ns") != private_stat.st_mtime_ns
        or private_stat.st_ctime_ns < private_pre["ctime_ns"]
        or not strict_json_equal(
            private_record,
            {"path": str(RECOVERY_CONTINUATION_PRIVATE_KEY), "mode": "0o0"},
        )
    ):
        raise ContinuationError("Continuation one-time key lifecycle drift")

    auth_fd, auth_data, auth_basic, auth_stat = open_bound(reader, AUTHORIZATION)
    auth_record = with_stat(auth_basic, auth_stat)
    prepare_attestation_fd, _, prepare_attestation_basic, prepare_attestation_stat = open_bound(
        reader, PREPARE_EXIT_ATTESTATION
    )
    prepare_attestation_record = with_stat(
        prepare_attestation_basic, prepare_attestation_stat
    )
    auth_sig_fd, _, auth_sig_basic, auth_sig_stat = open_bound(
        reader, AUTHORIZATION_SIGNATURE
    )
    auth_sig_record = with_stat(auth_sig_basic, auth_sig_stat)
    auth_public_fd, _, auth_public_record, _ = open_bound(
        reader, AUTHORIZATION_PUBLIC_KEY
    )
    completion_fd, completion_data, completion_basic, completion_stat = open_bound(
        reader, COMPLETION
    )
    completion_record = with_stat(completion_basic, completion_stat)
    promote_attestation_fd, _, promote_attestation_basic, promote_attestation_stat = open_bound(
        reader, PROMOTE_EXIT_ATTESTATION
    )
    promote_attestation_record = with_stat(
        promote_attestation_basic, promote_attestation_stat
    )
    completion_sig_fd, _, completion_sig_basic, completion_sig_stat = open_bound(
        reader, COMPLETION_SIGNATURE
    )
    completion_sig_record = with_stat(completion_sig_basic, completion_sig_stat)
    completion_public_fd, _, completion_public_record, _ = open_bound(
        reader, COMPLETION_PUBLIC_KEY
    )
    _, source_data, source_basic, source_stat = open_bound(
        reader, HISTORICAL_SOURCE_RECEIPT
    )
    _, canonical_data, canonical_basic, canonical_stat = open_bound(
        reader, CANONICAL_SOURCE_RECEIPT
    )
    if (
        auth_sig_stat.st_size != 64
        or completion_sig_stat.st_size != 64
        or not authority_module.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=prepare_attestation_fd,
            public_key_fd=auth_public_fd,
            signature_fd=auth_sig_fd,
        )
        or not authority_module.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=promote_attestation_fd,
            public_key_fd=completion_public_fd,
            signature_fd=completion_sig_fd,
        )
        or source_data != canonical_data
        or source_basic["sha256"] != EXPECTED_SOURCE_RECEIPT_SHA256
        or canonical_basic["sha256"] != EXPECTED_SOURCE_RECEIPT_SHA256
        or (source_stat.st_dev, source_stat.st_ino)
        == (canonical_stat.st_dev, canonical_stat.st_ino)
    ):
        raise ContinuationError("Signed-terminal promotion chain signature drift")
    authorization = strict_json(auth_data, "signed-terminal promotion authorization")
    completion = strict_json(completion_data, "signed-terminal promotion completion")
    validate_signed_terminal_precommit_contract(
        authorization=authorization,
        completion=completion,
        active_self_record=self_record,
        authorized_self_record=chain["authorized_self_record"],
        authorized_public_record=chain["continuation_public_key_record"],
        failed_v16_public_record=chain[
            "failed_v16_continuation_public_key_record"
        ],
        failed_v17_public_record=chain[
            "failed_v17_continuation_public_key_record"
        ],
        failed_v18_public_record=chain[
            "failed_v18_continuation_public_key_record"
        ],
        failed_v19_public_record=chain[
            "failed_v19_continuation_public_key_record"
        ],
        rejected_v20_public_record=chain[
            "rejected_v20_continuation_public_key_record"
        ],
        failed_v21_public_record=chain[
            "failed_v21_continuation_public_key_record"
        ],
        rejected_v22_public_record=chain[
            "rejected_v22_continuation_public_key_record"
        ],
        rejected_v23_public_record=chain[
            "rejected_v23_continuation_public_key_record"
        ],
        rejected_v24_public_record=chain[
            "rejected_v24_continuation_public_key_record"
        ],
        failed_v25_public_record=chain[
            "failed_v25_continuation_public_key_record"
        ],
        failed_v26_public_record=chain[
            "failed_v26_continuation_public_key_record"
        ],
        rejected_v27_public_record=chain[
            "rejected_v27_continuation_public_key_record"
        ],
        failed_v28_public_record=chain[
            "failed_v28_continuation_public_key_record"
        ],
        failed_v29_public_record=chain[
            "failed_v29_continuation_public_key_record"
        ],
        recovery_public_record=public_record,
        canonical_record=canonical_basic,
        source_record=with_stat(source_basic, source_stat),
        recovery_authority=chain["schema_recovery_authority"],
    )
    signed_promotion = payload.get("signed_promotion")
    if (
        not isinstance(signed_promotion, Mapping)
        or not strict_json_equal(signed_promotion.get("authorization"), auth_record)
        or not strict_json_equal(
            signed_promotion.get("prepare_exit_attestation"),
            prepare_attestation_record,
        )
        or not strict_json_equal(
            signed_promotion.get("authorization_signature"), auth_sig_record
        )
        or not strict_json_equal(signed_promotion.get("completion"), completion_record)
        or not strict_json_equal(
            signed_promotion.get("promote_exit_attestation"),
            promote_attestation_record,
        )
        or not strict_json_equal(
            signed_promotion.get("completion_signature"), completion_sig_record
        )
        or not strict_json_equal(
            signed_promotion.get("historical_source_receipt"),
            with_stat(source_basic, source_stat),
        )
        or not strict_json_equal(
            signed_promotion.get("canonical_source_receipt"),
            with_stat(canonical_basic, canonical_stat),
        )
        or signed_promotion.get("authorization_signature_verified") is not True
        or signed_promotion.get("completion_signature_verified") is not True
    ):
        raise ContinuationError("Continuation receipt promotion cross-links drift")

    outputs = StreamBoundArtifactSet()
    (
        child_inventory,
        child_paths,
        child_directories,
        require_child_inventory_stable,
    ) = bind_complete_child_output_inventory(outputs)
    outputs.open(
        ORIGINAL_RUNNER_LOG,
        expected_mode="0o444",
        expected_link_count=1,
    )
    terminal_releases = payload.get("terminal_releases")
    if (
        not isinstance(terminal_releases, Mapping)
        or not strict_json_equal(
            terminal_releases.get("complete_child_output_inventory"), child_inventory
        )
        or terminal_releases.get(
            "all_declared_release_artifacts_mode_0444_link_count_one"
        )
        is not True
        or not strict_json_equal(
            terminal_releases.get("completion_log"),
            basic_projection(outputs.record_for(ORIGINAL_RUNNER_LOG)),
        )
        or terminal_releases.get("pass") is not True
    ):
        raise ContinuationError("Signed terminal child-output inventory drift")

    dataset_fd, _, dataset_record, dataset_stat = open_bound(
        reader, FINAL_DATASET_RELEASE_RECEIPT
    )
    dataset_sig_fd, _, dataset_sig_record, dataset_sig_stat = open_bound(
        reader, FINAL_DATASET_RELEASE_SIGNATURE
    )
    result_public_fd, _, result_public_record, _ = open_bound(
        reader, FINAL_RELEASE_PUBLIC_KEY
    )
    _, result_private_record, result_private_stat = open_private_metadata(
        reader, FINAL_RELEASE_PRIVATE_KEY, "0o0"
    )
    report_fd, _, report_record, report_stat = open_bound(
        reader, FINAL_REPORT_RELEASE_RECEIPT
    )
    report_sig_fd, _, report_sig_record, report_sig_stat = open_bound(
        reader, FINAL_REPORT_RELEASE_SIGNATURE
    )
    report_public_fd, _, report_public_record, _ = open_bound(
        reader, REPORT_RELEASE_PUBLIC_KEY
    )
    _, report_private_record, report_private_stat = open_private_metadata(
        reader, REPORT_RELEASE_PRIVATE_KEY, "0o0"
    )
    if (
        dataset_sig_stat.st_size != 64
        or report_sig_stat.st_size != 64
        or result_private_stat.st_nlink != 1
        or report_private_stat.st_nlink != 1
        or result_public_record["sha256"] != EXPECTED_FINAL_RELEASE_PUBLIC_KEY_SHA256
        or report_public_record["sha256"] != EXPECTED_REPORT_RELEASE_PUBLIC_KEY_SHA256
        or not authority_module.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=dataset_fd,
            public_key_fd=result_public_fd,
            signature_fd=dataset_sig_fd,
        )
        or not authority_module.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=report_fd,
            public_key_fd=report_public_fd,
            signature_fd=report_sig_fd,
        )
    ):
        raise ContinuationError("Nested dataset/report release signature drift")
    if (
        terminal_releases.get("dataset_release", {}).get("receipt")
        != with_stat(dataset_record, dataset_stat)
        or terminal_releases.get("dataset_release", {}).get("signature")
        != dataset_sig_record
        or terminal_releases.get("dataset_release", {}).get("public_key")
        != result_public_record
        or terminal_releases.get("dataset_release", {}).get("retired_private_key")
        != result_private_record
        or terminal_releases.get("report_release", {}).get("receipt")
        != with_stat(report_record, report_stat)
        or terminal_releases.get("report_release", {}).get("signature")
        != report_sig_record
        or terminal_releases.get("report_release", {}).get("public_key")
        != report_public_record
        or terminal_releases.get("report_release", {}).get("retired_private_key")
        != report_private_record
    ):
        raise ContinuationError("Nested terminal release artifact relation drift")

    execution_context = validate_execution_receipt(
        authority_module=authority_module,
        bootstrap_record=bootstrap_record,
        v7_validation=v7_validation,
        reader=reader,
        chain=chain,
        receipt_data=receipt_data,
        receipt_record=receipt_record,
        receipt_stat=receipt_stat,
    )
    exit_attestation = validate_exit_attestation(
        data=attestation_data,
        chain=chain,
        execution_context=execution_context,
    )
    success_witness_payload: Mapping[str, Any] | None = None
    success_witness_artifact: Mapping[str, Any] | None = None
    if require_success_witness:
        if (
            witness_data is None
            or witness_record is None
            or witness_stat is None
        ):
            raise ContinuationError("Durable supervisor-success witness was not bound")
        success_witness_payload = validate_supervisor_success_witness(
            data=witness_data,
            witness_record=witness_record,
            witness_stat=witness_stat,
            self_record=self_record,
            self_binding=self_binding,
            receipt_artifact=with_stat(receipt_record, receipt_stat),
            attestation_artifact=with_stat(attestation_record, attestation_stat),
            signature_artifact=with_stat(signature_record, signature_stat),
            public_record=public_record,
            private_record=private_record,
            execution_context=execution_context,
            exit_attestation=exit_attestation,
            child_inventory=child_inventory,
        )
        success_witness_artifact = with_stat(witness_record, witness_stat)

    protected_paths = {
        *child_paths,
        CONTINUATION_RUNNER,
        CONTINUATION_RECEIPT,
        CONTINUATION_EXIT_ATTESTATION,
        CONTINUATION_SIGNATURE,
        RECOVERY_CONTINUATION_PUBLIC_KEY,
        AUTHORIZATION,
        PREPARE_EXIT_ATTESTATION,
        AUTHORIZATION_SIGNATURE,
        AUTHORIZATION_PUBLIC_KEY,
        COMPLETION,
        PROMOTE_EXIT_ATTESTATION,
        COMPLETION_SIGNATURE,
        COMPLETION_PUBLIC_KEY,
        HISTORICAL_SOURCE_RECEIPT,
        CANONICAL_SOURCE_RECEIPT,
        FINAL_DATASET_RELEASE_RECEIPT,
        FINAL_DATASET_RELEASE_SIGNATURE,
        FINAL_RELEASE_PUBLIC_KEY,
        FINAL_REPORT_RELEASE_RECEIPT,
        FINAL_REPORT_RELEASE_SIGNATURE,
        REPORT_RELEASE_PUBLIC_KEY,
        ORIGINAL_RUNNER_LOG,
    }
    if require_success_witness:
        protected_paths.add(CONTINUATION_SUCCESS_WITNESS)
    sentinel = PathMutationSentinel(
        protected_paths,
        all_event_directories={*child_directories, PORTABLE_ASSET_ROOT},
    )
    require_child_inventory_stable()
    reader.require_paths_still_bound()
    require_bootstrap_stable()
    require_incident_absent()
    require_success_witness_state(required=require_success_witness)
    sentinel.assert_clean()
    outputs.retain_until_process_exit()
    reader.retain_until_process_exit()
    sentinel.retain_until_process_exit()
    result_payload = {
        "schema_name": (
            f"{CONTINUATION_SCHEMA}.signed_terminal_verification"
            if require_success_witness
            else f"{CONTINUATION_SCHEMA}.signed_terminal_prewitness_verification"
        ),
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "canonical_command": (
            signed_terminal_verification_command()
            if require_success_witness
            else signed_terminal_prewitness_verification_command()
        ),
        "continuation_source_binding": self_binding,
        "continuation_receipt": with_stat(receipt_record, receipt_stat),
        "continuation_exit_attestation": with_stat(
            attestation_record, attestation_stat
        ),
        "supervised_child_wait": exit_attestation["child_wait"],
        "supervisor_capability": execution_context["supervisor_capability"],
        "threat_model": CONTINUATION_THREAT_MODEL,
        "continuation_signature": signature_record,
        "continuation_public_key": public_record,
        "retired_continuation_private_key": private_record,
        "complete_child_output_inventory": child_inventory,
        "nested_dataset_release_signature_verified": True,
        "nested_report_release_signature_verified": True,
        "continuation_terminal_signature_verified": True,
        "supervised_child_normal_exit_zero_verified": True,
        "execution_receipt_full_nested_contract_verified": True,
        "release_authority_granted": require_success_witness,
        "pass": True,
        "pass_semantics": (
            "signed_terminal_release_chain_integrity_only_not_scientific_confirmation"
            if require_success_witness
            else PREWITNESS_VERIFICATION_PASS_SEMANTICS
        ),
    }
    if require_success_witness:
        if success_witness_payload is None or success_witness_artifact is None:
            raise ContinuationError("Validated supervisor-success witness was lost")
        result_payload.update(
            {
                "supervisor_success_witness": success_witness_artifact,
                "durable_supervisor_success_verified": True,
                "recovery_output_writer_lease_released_verified": True,
            }
        )
    else:
        result_payload.update(
            {
                "success_witness_required": True,
                "ready_for_success_witness": True,
            }
        )
    return result_payload


def record_continuation_incident(error: BaseException) -> None:
    if os.path.lexists(CONTINUATION_INCIDENT):
        return
    observed_paths = [
        str(path)
        for path in DOWNSTREAM_OUTPUT_SLOTS
        if os.path.lexists(path)
    ]
    if not observed_paths:
        return
    payload = {
        "schema_name": f"{CONTINUATION_SCHEMA}.incident",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "status": "DOWNSTREAM_ARTIFACTS_EXIST_WITHOUT_CLEAN_TERMINAL_RETURN",
        "error": {"type": type(error).__name__, "message": str(error)},
        "observed_paths": observed_paths,
        "continuation_receipt_exists": os.path.lexists(CONTINUATION_RECEIPT),
        "continuation_signature_exists": os.path.lexists(CONTINUATION_SIGNATURE),
        "continuation_signature_is_provisional": True,
        "supervisor_success_witness_exists": os.path.lexists(
            CONTINUATION_SUCCESS_WITNESS
        ),
        "success_witness_required_for_release_authority": True,
        "independent_signed_terminal_verification_required": True,
        "release_authority_granted": False,
        "required_action": (
            "Preserve and archive this incident and every observed output before any "
            "separately reviewed recovery; signed dataset/report receipts are provisional."
        ),
        "pass": False,
    }
    try:
        commit_terminal_receipt(payload, lambda: None, path=CONTINUATION_INCIDENT)
    except BaseException:
        return


def acquire_recovery_output_lease(*, lock_operation: int, conflict: str) -> int:
    path = RECOVERY_OUTPUT_LEASE_PUBLIC_KEY
    expected_sha256 = EXPECTED_RECOVERY_OUTPUT_LEASE_PUBLIC_KEY_SHA256
    if path.resolve(strict=True) != path:
        raise ContinuationError("Recovery writer-lock public-key path is not canonical")
    before = os.stat(path, follow_symlinks=False)
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    try:
        opened = os.fstat(descriptor)
        data = os.pread(descriptor, opened.st_size, 0)
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
        if (
            not stat.S_ISREG(opened.st_mode)
            or stat.S_IMODE(opened.st_mode) != 0o444
            or opened.st_nlink != 1
            or identity(before) != identity(opened)
            or identity(opened) != identity(after)
            or len(data) != opened.st_size
            or hashlib.sha256(data).hexdigest() != expected_sha256
        ):
            raise ContinuationError("Recovery output-lock public-key identity drift")
        try:
            fcntl.flock(descriptor, lock_operation | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ContinuationError(conflict) from error
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def acquire_recovery_output_writer_lease() -> int:
    return acquire_recovery_output_lease(
        lock_operation=fcntl.LOCK_EX,
        conflict="Another cooperating recovery output writer holds the v30 lease",
    )


def acquire_recovery_output_verifier_lease() -> int:
    return acquire_recovery_output_lease(
        lock_operation=fcntl.LOCK_SH,
        conflict=(
            "Recovery output writer is still active; terminal witness is not yet "
            "eligible for final release verification"
        ),
    )


def should_record_continuation_incident(
    *, signed_terminal_verification: bool, writer_lease_fd: int
) -> bool:
    return not signed_terminal_verification and writer_lease_fd >= 0


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) == 2 else ""
    signed_terminal_verification = mode == "--verify-signed-terminal"
    signed_prewitness_verification = (
        mode == "--verify-signed-terminal-prewitness"
    )
    read_only_verification = (
        signed_terminal_verification or signed_prewitness_verification
    )
    writer_lease_fd = -1
    try:
        if signed_terminal_verification:
            payload = verify_signed_terminal(require_success_witness=True)
            print(json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True))
            return 0
        if signed_prewitness_verification:
            payload = verify_signed_terminal(require_success_witness=False)
            print(json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True))
            return 0
        if mode == "--supervise-and-sign":
            writer_lease_fd = acquire_recovery_output_writer_lease()
            supervise_and_sign()
            os._exit(0)
        if mode == "--supervised-child":
            run()
            os._exit(0)
        raise ContinuationError(
            "Expected exactly --supervise-and-sign, --supervised-child, or "
            "--verify-signed-terminal/--verify-signed-terminal-prewitness"
        )
    except BaseException as error:
        failure = {
            "schema_name": f"{CONTINUATION_SCHEMA}.error",
            "schema_version": "1.0.0",
            "error": {"type": type(error).__name__, "message": str(error)},
            "pass": False,
        }
        try:
            os.write(
                2,
                (json.dumps(failure, ensure_ascii=True, sort_keys=True) + "\n").encode(
                    "ascii"
                ),
            )
        except OSError:
            pass
        if should_record_continuation_incident(
            signed_terminal_verification=read_only_verification,
            writer_lease_fd=writer_lease_fd,
        ):
            record_continuation_incident(error)
        os._exit(1)


if __name__ == "__main__":
    raise SystemExit(main())
