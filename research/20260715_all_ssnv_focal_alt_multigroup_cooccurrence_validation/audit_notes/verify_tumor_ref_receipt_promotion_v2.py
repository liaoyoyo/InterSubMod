#!/usr/bin/env python3
"""Verify the signed tumor-REF receipt promotion before Task-B continuation."""

from __future__ import annotations

from datetime import datetime, timezone
import ctypes
import hashlib
import json
import os
from pathlib import Path
import stat
import sys
from types import CodeType, FunctionType, ModuleType
from typing import Any, Callable, Mapping


AT_FDCWD = -100
AT_SYMLINK_FOLLOW = 0x400


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
SCRIPT_ROOT = TOPIC_ROOT / "scripts"
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
RESULT_ROOT = TOPIC_ROOT / "results"
WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)

PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")
PRIMARY_PYTHON_RUNTIME = Path(
    "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python3.11"
)
QA_PYTHON_RUNTIME = Path("/bip7_disk/liaoyoyo2001/miniconda3/bin/python3.9")
VERIFIER = AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_v2.py"
PROMOTION_TOOL = AUDIT_ROOT / "promote_tumor_ref_source_receipt_v2.py"
RUNNER_GATE_REPLAYER = AUDIT_ROOT / "replay_m2v5_runner_only_gates_v1.py"
DOWNSTREAM_CONTINUATION = (
    AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_v1.py"
)
RELEASE_SOURCE_AUTHORITY = SCRIPT_ROOT / "release_source_authority.py"
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
OPENSSL = Path("/usr/bin/openssl")

SOURCE_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json"
)
CANONICAL_RECEIPT = (
    WORKSPACE_ROOT
    / "tumor_ref_recovery_source_identity_v1"
    / "post_run_source_identity.receipt.json"
)
AUTHORIZATION = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_authorization.v3.json"
)
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
VERIFICATION_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.v3.json"
)
RUNNER_GATE_RECEIPT = RESULT_ROOT / "m2v5_runner_only_gate_replay.v1.json"
RUNNER_GATE_LOG = RESULT_ROOT / "m2v5_runner_only_gate_replay.v1.log"
PROMOTION_INCIDENT = RESULT_ROOT / "tumor_ref_source_receipt_promotion_incident.v3.json"
PROMOTION_CACHE_ROOTS = {
    "--prepare-authorization": (
        WORKSPACE_ROOT / ".python_cache_tumor_ref_promotion_v2_prepare"
    ),
    "--promote": WORKSPACE_ROOT / ".python_cache_tumor_ref_promotion_v2_promote",
}

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
CONTINUATION_TERMINAL_RECEIPT = RESULT_ROOT / "m2v5_downstream_continuation.v1.json"
CONTINUATION_EXIT_ATTESTATION = (
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.v1.json"
)
CONTINUATION_TERMINAL_SIGNATURE = Path(
    str(CONTINUATION_EXIT_ATTESTATION) + ".ed25519.sig"
)
SUPERVISOR_CAPABILITY_PROTOCOL = {
    "mechanism": "sealed_memfd_with_256_bit_nonce_and_live_parent_binding",
    "inherited_descriptor": 198,
    "required_memfd_mode": "0o400",
    "required_memfd_link_count": 0,
    "required_seals": ["F_SEAL_SEAL", "F_SEAL_SHRINK", "F_SEAL_GROW", "F_SEAL_WRITE"],
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
CONTINUATION_POLICY = (
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

INCIDENT_ARCHIVE_RECEIPT = RESULT_ROOT / "completion_attempt2_archive_receipt.v1.json"
REVIEW_SLOTS = (
    (
        "mendel",
        "Mendel",
        "mendel-v3-distinct-review-artifact",
        TOPIC_ROOT / "reviews" / "20260722_tumor_ref_receipt_promotion_mendel.v3.json",
    ),
    (
        "nash",
        "Nash",
        "nash-v3-distinct-review-artifact",
        TOPIC_ROOT / "reviews" / "20260722_tumor_ref_receipt_promotion_nash.v3.json",
    ),
    (
        "external_claude_opus",
        "External Claude Opus",
        "external-claude-opus-v3-distinct-review-artifact",
        TOPIC_ROOT
        / "reviews"
        / "20260722_tumor_ref_receipt_promotion_external_claude_opus.v3.json",
    ),
)

DOWNSTREAM_OUTPUTS = (
    WORKSPACE_ROOT / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT
    / "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5",
    RESULT_ROOT / "stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json",
    RESULT_ROOT
    / "frozen_input_immutability.post_m2v5_downstream_v3_source_authority_v5.json",
    WORKSPACE_ROOT / "all_ssnv_final_report_dataset_v5_m2v5_source_attested",
    WORKSPACE_ROOT / "all_ssnv_final_report_v5_m2v5_source_attested",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json.ed25519.sig",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json.ed25519.sig",
    WORKSPACE_ROOT / "m2v5_recovered_completion_chain_v2_source_authority_v5.log",
    WORKSPACE_ROOT / ".python_cache_m2v5_completion_v2_bound_bootstrap",
    RESULT_ROOT / "m2v5_downstream_continuation.v1.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.v1.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.v1.json.ed25519.sig",
)

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
EXPECTED_SOURCE_SIZE = 6_733
EXPECTED_SOURCE_SHA256 = (
    "d9fd6ecfc92b1e76ef6b68932faaab1e3d8e8540f521e5a7065f760bdd228e19"
)
EXPECTED_AUTHORIZATION_PUBLIC_KEY_SHA256 = (
    "e638b29a1d9c207dfa9849ff20a3867dacf702f7c1a899f1968ea1401a839677"
)
EXPECTED_COMPLETION_PUBLIC_KEY_SHA256 = (
    "6c55230f368db6b1fb0eedbc6c2362c0df7ee41a2699764fadb1f424e7036032"
)
EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "a98370415594ebc9e8eb166bb70cab3550067b8a83db7d73ad23fac0891cb8b5"
)
EXPECTED_RELEASE_SOURCE_AUTHORITY_SHA256 = (
    "36621dc735d17cf3322e34aeb5d72247869935e7b87306cdb600ff6d522c52af"
)
EXPECTED_OPENSSL_SHA256 = (
    "38c064f53a6619364c7947fc11cad09e1fc4218fc2c3e9016a7786b1341ecd08"
)
EXPECTED_PYTHON_SHA256 = (
    "777797a57eb75b28f530191628d26b14afada9ced2cb51c0ecae1eb62796062e"
)
EXPECTED_QA_PYTHON_SHA256 = (
    "cbd71e11d993aa3b7cc3ad553b97733f8054c1465ed7764dbbd386b9897bd5bc"
)
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
    "primary_python_runtime": (PRIMARY_PYTHON_RUNTIME, EXPECTED_PYTHON_SHA256, "0o775"),
    "qa_python_runtime": (QA_PYTHON_RUNTIME, EXPECTED_QA_PYTHON_SHA256, "0o755"),
    "portable_builder_module": (
        PORTABLE_BUILDER_MODULE,
        EXPECTED_PORTABLE_BUILDER_SHA256,
        "0o644",
    ),
    "portable_chart_extractor_module": (
        PORTABLE_CHART_EXTRACTOR_MODULE,
        EXPECTED_PORTABLE_CHART_EXTRACTOR_SHA256,
        "0o644",
    ),
    "portable_verifier_module": (
        PORTABLE_VERIFIER_MODULE,
        EXPECTED_PORTABLE_VERIFIER_SHA256,
        "0o644",
    ),
    "portable_browser_helpers_module": (
        PORTABLE_BROWSER_HELPERS_MODULE,
        EXPECTED_PORTABLE_BROWSER_HELPERS_SHA256,
        "0o644",
    ),
    "portable_browser_cli_module": (
        PORTABLE_BROWSER_CLI_MODULE,
        EXPECTED_PORTABLE_BROWSER_CLI_SHA256,
        "0o644",
    ),
    "portable_server_bundle": (
        PORTABLE_SERVER_BUNDLE,
        EXPECTED_PORTABLE_SERVER_BUNDLE_SHA256,
        "0o644",
    ),
    **{
        f"portable_reader_asset_part{index:03d}": (
            path,
            EXPECTED_PORTABLE_READER_ASSET_SHA256[index - 1],
            "0o644",
        )
        for index, path in enumerate(PORTABLE_READER_ASSET_PARTS, start=1)
    },
}
EXPECTED_INCIDENT_ARCHIVE_SHA256 = (
    "ad33b78379e387b68dab948cd60644e613c0387884fd4d6e35748fd1bb82dbed"
)
EXPECTED_V7_AUTHORITY_ID = (
    "20260722_all_ssnv_focal_alt_task_b_release_v7_strict_command_parity"
)
EXPECTED_V3_AUTHORITY_ID = "20260717_all_ssnv_focal_alt_task_b_release_v3"
EXPECTED_V3_PUBLIC_KEY_PATH = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
    "20260717_all_ssnv_v3/ed25519_public.pem"
)
EXPECTED_V3_PUBLIC_KEY_SHA256 = (
    "60ebac3ee2ebfbf69a80331b40410b365c0315d41363d59ee0b44f6dbf5040e4"
)
EXPECTED_V7_PUBLIC_KEY_PATH = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
    "20260722_all_ssnv_v10_strict_command_parity_bootstrap/ed25519_public.pem"
)
EXPECTED_V7_PUBLIC_KEY_SHA256 = (
    "cecb2287cf87f8bd948c7390bd5f7059578966e0fc917de420572ddd82d8f21b"
)
EXPECTED_V7_SOURCE_SET_SHA256 = (
    "f5ba8a9e786971c4261b51e283a0f9df6e807a8aca695bf5041271fee5420f58"
)

AUTHORIZATION_SCHEMA = (
    "intersubmod.tumor_ref_source_receipt_promotion_authorization"
)
COMPLETION_SCHEMA = "intersubmod.tumor_ref_source_receipt_promotion"
VERIFICATION_SCHEMA = "intersubmod.tumor_ref_source_receipt_promotion_verification"
PREPARE_EXIT_ATTESTATION_SCHEMA = (
    "intersubmod.tumor_ref_source_receipt_promotion.prepare_exit_attestation"
)
PROMOTE_EXIT_ATTESTATION_SCHEMA = (
    "intersubmod.tumor_ref_source_receipt_promotion.promote_exit_attestation"
)
REVIEW_SCHEMA = "intersubmod.tumor_ref_receipt_promotion_review"
SCHEMA_VERSION = "3.0.0"
EXIT_ATTESTATION_SCHEMA_VERSION = "1.0.0"
AUTHORIZATION_ID = "20260722_tumor_ref_receipt_promotion_v3"
AUTHORIZATION_STATUS = "AUTHORIZED_FOR_SINGLE_BYTE_IDENTICAL_PROMOTION"

BASIC_ARTIFACT_KEYS = frozenset({"path", "size_bytes", "sha256", "mode"})
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
AUTHORIZATION_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "authorization_id",
        "task_type",
        "status",
        "scope",
        "producer_exit_attestation",
        "commands",
        "reviewed_sources",
        "execution_binding",
        "producer_session",
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
        "producer_exit_attestation",
        "command",
        "authorization",
        "code",
        "execution_binding",
        "producer_session",
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
AUTHORIZATION_LINK_KEYS = frozenset(
    {
        "artifact",
        "producer_exit_attestation",
        "signature",
        "public_key",
        "retired_private_key",
        "producer_session",
        "signature_verified",
    }
)
SIGNATURE_CONTRACT_KEYS = frozenset(
    {
        "attestation_schema",
        "algorithm",
        "signed_artifact",
        "signature",
        "public_key",
        "private_key_lifecycle",
        "legacy_payload_signature",
    }
)
AUTHORIZATION_SIGNATURE_CONTRACT_KEYS = SIGNATURE_CONTRACT_KEYS
AUTHORIZED_COMPLETION_SIGNATURE_CONTRACT_KEYS = SIGNATURE_CONTRACT_KEYS
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
SIGNED_AUTHORITY_EVIDENCE_KEYS = frozenset(
    {
        "authority",
        "approval",
        "signature",
        "public_key",
        "private_key",
        "authority_id",
        "signature_verified",
    }
)

EXIT_ATTESTATION_KEYS = frozenset(
    {
        "created_at_utc",
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
CHILD_WAIT_KEYS = frozenset(
    {
        "waited_pid",
        "raw_wait_status",
        "wifexited",
        "exit_status",
        "wifsignaled",
        "terminating_signal",
        "core_dumped",
    }
)
PREPARE_EXIT_ATTESTATION_CHECKS = {
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
PROMOTE_EXIT_ATTESTATION_CHECKS = {
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
PREPARE_EXIT_ATTESTATION_PASS_SEMANTICS = (
    "parent_observed_fork_no_exec_child_normal_exit_zero_and_strictly_revalidated_"
    "the_authorization_output; this_attestation_requires_its_detached_authorization_"
    "key_signature_before_promotion"
)
PROMOTE_EXIT_ATTESTATION_PASS_SEMANTICS = (
    "parent_observed_fork_no_exec_child_normal_exit_zero_and_strictly_revalidated_"
    "the_completion_and_canonical_outputs_after_independently_verifying_the_signed_"
    "prepare_attestation; downstream_authority_still_requires_the_detached_completion_"
    "key_signature_and_all_continuation_gates"
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

_BOOTSTRAP_FD_LEASES: list[tuple[Path, int, os.stat_result]] = []
_DIRECT_PATH_LEASES: list[tuple[Path, os.stat_result]] = []
_INHERITED_VERIFIER_FD_LEASE: tuple[int, os.stat_result] | None = None
_INHERITED_PYTHON_FD_LEASE: tuple[int, os.stat_result] | None = None


class VerificationError(RuntimeError):
    """Raised when the signed continuation chain is not exact."""


def _same_stat(left: os.stat_result, right: os.stat_result) -> bool:
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


def _stable_code_payload(code: CodeType) -> dict[str, Any]:
    constants: list[Any] = []
    for value in code.co_consts:
        if isinstance(value, CodeType):
            constants.append({"code": _stable_code_payload(value)})
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


def _stable_code_digest(code: CodeType) -> str:
    return hashlib.sha256(
        json.dumps(_stable_code_payload(code), sort_keys=True).encode("utf-8")
    ).hexdigest()


def _code_inventory(namespace: Mapping[str, Any], filename: str) -> dict[str, str]:
    inventory: dict[str, str] = {}
    for name, value in namespace.items():
        if isinstance(value, FunctionType) and value.__code__.co_filename == filename:
            inventory[name] = _stable_code_digest(value.__code__)
        elif isinstance(value, type) and value.__module__ == namespace.get("__name__"):
            for method_name, method in vars(value).items():
                function = (
                    method.__func__
                    if isinstance(method, (staticmethod, classmethod))
                    else method
                )
                if (
                    isinstance(function, FunctionType)
                    and function.__code__.co_filename == filename
                ):
                    inventory[f"{name}.{method_name}"] = _stable_code_digest(
                        function.__code__
                    )
    return inventory


def _critical_globals_equal(
    live: Mapping[str, Any], shadow: Mapping[str, Any]
) -> bool:
    names = (
        "REPO_ROOT",
        "TOPIC_ROOT",
        "AUDIT_ROOT",
        "RESULT_ROOT",
        "WORKSPACE_ROOT",
        "PYTHON",
        "PRIMARY_PYTHON_RUNTIME",
        "QA_PYTHON_RUNTIME",
        "VERIFIER",
        "PROMOTION_TOOL",
        "RUNNER_GATE_REPLAYER",
        "DOWNSTREAM_CONTINUATION",
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
        "AUTHORIZATION_SIGNATURE",
        "PREPARE_EXIT_ATTESTATION",
        "LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE",
        "COMPLETION",
        "COMPLETION_SIGNATURE",
        "PROMOTE_EXIT_ATTESTATION",
        "LEGACY_COMPLETION_PAYLOAD_SIGNATURE",
        "CONTINUATION_PUBLIC_KEY",
        "CONTINUATION_PRIVATE_KEY",
        "CONTINUATION_TERMINAL_RECEIPT",
        "CONTINUATION_EXIT_ATTESTATION",
        "CONTINUATION_TERMINAL_SIGNATURE",
        "SUPERVISOR_CAPABILITY_PROTOCOL",
        "CONTINUATION_THREAT_MODEL",
        "CONTINUATION_POLICY",
        "VERIFICATION_RECEIPT",
        "RUNNER_GATE_RECEIPT",
        "RUNNER_GATE_LOG",
        "PROMOTION_INCIDENT",
        "REVIEW_SLOTS",
        "DOWNSTREAM_OUTPUTS",
        "EXPECTED_ENVIRONMENT",
        "AUTHORIZATION_SCHEMA",
        "COMPLETION_SCHEMA",
        "VERIFICATION_SCHEMA",
        "PREPARE_EXIT_ATTESTATION_SCHEMA",
        "PROMOTE_EXIT_ATTESTATION_SCHEMA",
        "SCHEMA_VERSION",
        "EXIT_ATTESTATION_SCHEMA_VERSION",
        "AUTHORIZATION_ID",
        "AUTHORIZATION_KEYS",
        "COMPLETION_KEYS",
        "EXIT_ATTESTATION_KEYS",
        "PRODUCER_SESSION_KEYS",
        "CHILD_WAIT_KEYS",
        "PREPARE_EXIT_ATTESTATION_CHECKS",
        "PROMOTE_EXIT_ATTESTATION_CHECKS",
        "PREPARE_EXIT_ATTESTATION_PASS_SEMANTICS",
        "PROMOTE_EXIT_ATTESTATION_PASS_SEMANTICS",
        "AUTHORIZATION_CHECKS",
        "COMPLETION_CHECKS",
        "EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256",
        "REVIEWED_RUNTIME_SOURCE_CONTRACTS",
    )
    return all(
        name in live
        and name in shadow
        and type(live[name]) is type(shadow[name])
        and live[name] == shadow[name]
        for name in names
    )


def _executing_module_code() -> CodeType:
    frame = sys._getframe()
    while frame is not None:
        if frame.f_code.co_name == "<module>" and frame.f_globals is globals():
            return frame.f_code
        frame = frame.f_back
    raise VerificationError("Executing verifier module frame is unavailable")


def _require_executing_verifier_matches_bound_source(data: bytes) -> None:
    module_code = _executing_module_code()
    filename = module_code.co_filename
    shadow: dict[str, Any] = {
        "__name__": "promotion_verifier_v2_bound_shadow",
        "__file__": filename,
        "__package__": "",
    }
    reviewed_module_code = compile(data, filename, "exec")
    exec(reviewed_module_code, shadow)
    live_inventory = _code_inventory(globals(), filename)
    shadow_inventory = _code_inventory(shadow, filename)
    if (
        not live_inventory
        or _stable_code_digest(module_code) != _stable_code_digest(reviewed_module_code)
        or live_inventory != shadow_inventory
        or not _critical_globals_equal(globals(), shadow)
    ):
        raise VerificationError("Executing verifier module differs from bound source bytes")


def _same_inode_and_mode(left: os.stat_result, right: os.stat_result) -> bool:
    return (
        left.st_dev == right.st_dev
        and left.st_ino == right.st_ino
        and left.st_mode == right.st_mode
    )


def _require_absolute_regular_path(path: Path) -> os.stat_result:
    if not path.is_absolute():
        raise VerificationError(f"Trust artifact path is not absolute: {path}")
    try:
        observed = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise VerificationError(f"Trust artifact is unavailable: {path}") from error
    if not stat.S_ISREG(observed.st_mode):
        raise VerificationError(f"Trust artifact is not a direct regular file: {path}")
    if path.resolve(strict=True) != path:
        raise VerificationError(f"Trust artifact path contains a symlink: {path}")
    return observed


def _bootstrap_open(path: Path) -> tuple[bytes, dict[str, Any]]:
    """Bind the authority module before its own reader class is available."""
    before_path = _require_absolute_regular_path(path)
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    try:
        first = os.fstat(descriptor)
        if not _same_stat(first, before_path):
            raise VerificationError(f"Bootstrap artifact path rebound: {path}")
        data = b"".join(
            os.pread(descriptor, min(8 * 1024 * 1024, first.st_size - offset), offset)
            for offset in range(0, first.st_size, 8 * 1024 * 1024)
        )
        second = os.fstat(descriptor)
        after_path = os.stat(path, follow_symlinks=False)
        if len(data) != first.st_size or not _same_stat(first, second) or not _same_stat(
            second, after_path
        ):
            raise VerificationError(f"Bootstrap artifact changed while reading: {path}")
    except Exception:
        os.close(descriptor)
        raise
    _BOOTSTRAP_FD_LEASES.append((path, descriptor, second))
    return data, {
        "path": str(path),
        "size_bytes": second.st_size,
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": oct(stat.S_IMODE(second.st_mode)),
    }


def _require_bootstrap_paths_still_bound() -> None:
    for path, descriptor, opened in _BOOTSTRAP_FD_LEASES:
        descriptor_stat = os.fstat(descriptor)
        live = os.stat(path, follow_symlinks=False)
        if not _same_stat(opened, descriptor_stat) or not _same_stat(opened, live):
            raise VerificationError(f"Bootstrap trust artifact binding changed: {path}")


def _require_direct_paths_still_bound() -> None:
    for path, opened in _DIRECT_PATH_LEASES:
        live = os.stat(path, follow_symlinks=False)
        if not _same_stat(opened, live):
            raise VerificationError(f"Direct trust artifact path binding changed: {path}")
    if _INHERITED_VERIFIER_FD_LEASE is not None:
        descriptor, opened = _INHERITED_VERIFIER_FD_LEASE
        descriptor_stat = os.fstat(descriptor)
        live = os.stat(VERIFIER, follow_symlinks=False)
        if not _same_stat(opened, descriptor_stat) or not _same_stat(opened, live):
            raise VerificationError("Inherited verifier descriptor binding changed")
    if _INHERITED_PYTHON_FD_LEASE is not None:
        descriptor, opened = _INHERITED_PYTHON_FD_LEASE
        descriptor_stat = os.fstat(descriptor)
        process_exe = os.stat("/proc/self/exe")
        if not _same_stat(opened, descriptor_stat) or not _same_stat(
            opened, process_exe
        ):
            raise VerificationError("Inherited Python descriptor binding changed")


def _load_authority_module() -> tuple[ModuleType, dict[str, Any]]:
    data, record = _bootstrap_open(RELEASE_SOURCE_AUTHORITY)
    if not _strict_json_equal(record, {
        "path": str(RELEASE_SOURCE_AUTHORITY),
        "size_bytes": len(data),
        "sha256": EXPECTED_RELEASE_SOURCE_AUTHORITY_SHA256,
        "mode": "0o444",
    }):
        raise VerificationError("release_source_authority.py identity drift")
    module = ModuleType("bound_release_source_authority_for_promotion_verifier_v2")
    module.__file__ = str(RELEASE_SOURCE_AUTHORITY)
    module.__package__ = ""
    exec(compile(data, str(RELEASE_SOURCE_AUTHORITY), "exec"), module.__dict__)
    if not hasattr(module, "BoundArtifactReader") or not hasattr(
        module, "verify_ed25519_signature_fds"
    ):
        raise VerificationError("Bound release authority API is incomplete")
    return module, record


def _observed_command() -> list[str]:
    raw = Path("/proc/self/cmdline").read_bytes()
    if not raw.endswith(b"\0"):
        raise VerificationError("Direct process command is unavailable")
    return [os.fsdecode(token) for token in raw[:-1].split(b"\0")]


def _expected_command(mode: str) -> list[str]:
    return [str(PYTHON), "-I", "-B", str(VERIFIER), mode]


def _promotion_command(mode: str) -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        "-X",
        f"pycache_prefix={PROMOTION_CACHE_ROOTS[mode]}",
        str(PROMOTION_TOOL),
        mode,
    ]


def _runner_gate_command() -> list[str]:
    return [str(PYTHON), "-I", "-B", str(RUNNER_GATE_REPLAYER)]


def _downstream_continuation_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(DOWNSTREAM_CONTINUATION),
        "--supervise-and-sign",
    ]


def _supervised_child_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(DOWNSTREAM_CONTINUATION),
        "--supervised-child",
    ]


def _require_clean_direct_runtime(mode: str) -> None:
    global _INHERITED_PYTHON_FD_LEASE, _INHERITED_VERIFIER_FD_LEASE

    if mode not in {"--verify", "--verify-and-record"}:
        raise VerificationError("Exactly one canonical verification mode is required")
    if dict(os.environ) != EXPECTED_ENVIRONMENT:
        raise VerificationError("Verification environment is not the exact clean environment")
    observed = _observed_command()
    canonical = _expected_command(mode)
    if observed == canonical and sys.argv == [str(VERIFIER), mode]:
        if mode != "--verify-and-record":
            raise VerificationError(
                "Fresh --verify must use inherited Python and verifier descriptors"
            )
        if Path(sys.executable) != PYTHON:
            raise VerificationError("Python launcher token drift")
    else:
        if (
            len(observed) != 5
            or observed[1:3] != canonical[1:3]
            or observed[4] != mode
            or len(sys.argv) != 2
            or sys.argv[0] != observed[3]
            or sys.argv[1] != mode
            or not observed[0].startswith("/proc/self/fd/")
            or not observed[3].startswith("/proc/self/fd/")
        ):
            raise VerificationError(
                "Verifier was not launched by the canonical or inherited-FD command"
            )
        inherited: list[tuple[int, os.stat_result]] = []
        for token, label, canonical_path in (
            (observed[0], "Python", PYTHON.resolve(strict=True)),
            (observed[3], "verifier", VERIFIER),
        ):
            descriptor_token = token.removeprefix("/proc/self/fd/")
            if (
                not descriptor_token.isdecimal()
                or str(int(descriptor_token)) != descriptor_token
                or int(descriptor_token) < 3
            ):
                raise VerificationError(
                    f"Inherited {label} descriptor token is not canonical"
                )
            descriptor = int(descriptor_token)
            try:
                descriptor_stat = os.fstat(descriptor)
            except OSError as error:
                raise VerificationError(
                    f"Inherited {label} descriptor is unavailable"
                ) from error
            canonical_stat = _require_absolute_regular_path(canonical_path)
            if not _same_stat(descriptor_stat, canonical_stat):
                raise VerificationError(
                    f"Inherited {label} descriptor is not the canonical inode"
                )
            inherited.append((descriptor, descriptor_stat))
        if sys.executable != observed[0]:
            raise VerificationError("Inherited Python launcher token drift")
        _INHERITED_PYTHON_FD_LEASE = inherited[0]
        _INHERITED_VERIFIER_FD_LEASE = inherited[1]
    if not (
        sys.flags.isolated == 1
        and sys.flags.ignore_environment == 1
        and sys.flags.no_user_site == 1
        and sys.flags.dont_write_bytecode == 1
    ):
        raise VerificationError("Python isolation flags are incomplete")


def _parse_json(data: bytes, label: str) -> dict[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        output: dict[str, Any] = {}
        for key, value in pairs:
            if key in output:
                raise VerificationError(f"Duplicate JSON key in {label}: {key}")
            output[key] = value
        return output

    def reject_nonfinite(value: str) -> None:
        raise VerificationError(f"Non-finite JSON value in {label}: {value}")

    try:
        payload = json.loads(
            data,
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_nonfinite,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise VerificationError(f"Invalid JSON in {label}") from error
    if not isinstance(payload, dict):
        raise VerificationError(f"{label} is not a JSON object")
    return payload


def _require_utc_timestamp(value: Any, label: str) -> datetime:
    if not isinstance(value, str):
        raise VerificationError(f"{label} is not a UTC timestamp string")
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise VerificationError(f"{label} is not ISO-8601") from error
    if (
        parsed.tzinfo is None
        or parsed.utcoffset() != timezone.utc.utcoffset(parsed)
        or parsed.microsecond != 0
        or value != parsed.isoformat(timespec="seconds")
    ):
        raise VerificationError(f"{label} is not UTC")
    return parsed


def _require_exact_keys(value: Any, expected: frozenset[str], label: str) -> None:
    if not isinstance(value, Mapping) or set(value) != set(expected):
        observed = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise VerificationError(
            f"{label} keys are not exact: observed={observed}, expected={sorted(expected)}"
        )


def _strict_json_equal(observed: Any, expected: Any) -> bool:
    if type(observed) is not type(expected):
        return False
    if isinstance(expected, dict):
        return set(observed) == set(expected) and all(
            _strict_json_equal(observed[key], expected[key]) for key in expected
        )
    if isinstance(expected, list):
        return len(observed) == len(expected) and all(
            _strict_json_equal(left, right)
            for left, right in zip(observed, expected, strict=True)
        )
    return observed == expected


def _require_json_equal(observed: Any, expected: Any, label: str) -> None:
    if not _strict_json_equal(observed, expected):
        raise VerificationError(f"{label} differs from the exact JSON contract")


def _require_absent(path: Path, label: str) -> None:
    if os.path.lexists(path):
        raise VerificationError(f"{label} must be absent: {path}")


def _require_lower_sha256(value: Any, label: str) -> None:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise VerificationError(f"{label} is not a lowercase SHA-256 string")


def _exit_attestation_schema(name: str) -> dict[str, str]:
    return {"name": name, "version": EXIT_ATTESTATION_SCHEMA_VERSION}


def _exit_attestation_announcement(
    path: Path, schema_name: str, public_key: Mapping[str, Any]
) -> dict[str, Any]:
    return {
        "path": str(path),
        "schema": _exit_attestation_schema(schema_name),
        "signing_key": dict(public_key),
    }


def _exit_attestation_signature_contract(
    *,
    path: Path,
    schema_name: str,
    signature: Path,
    public_key: Mapping[str, Any],
    private_key: Path,
    legacy_payload_signature: Path,
) -> dict[str, Any]:
    return {
        "attestation_schema": _exit_attestation_schema(schema_name),
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


def _validate_producer_session(
    value: Any, label: str, *, expected_mode: str
) -> dict[str, Any]:
    _require_exact_keys(value, PRODUCER_SESSION_KEYS, label)
    _require_lower_sha256(value.get("nonce_sha256"), f"{label} nonce_sha256")
    if (
        value.get("protocol") != "intersubmod.parent_fork_no_exec_child"
        or value.get("schema_version") != "1.0.0"
        or value.get("mode") != expected_mode
    ):
        raise VerificationError(f"{label} protocol/schema/mode contract drift")
    for field in (
        "parent_pid",
        "parent_start_time_ticks",
        "child_pid",
        "child_start_time_ticks",
    ):
        observed = value.get(field)
        if type(observed) is not int or observed <= 0:
            raise VerificationError(f"{label} {field} is not a positive JSON integer")
    if value["parent_pid"] <= 1 or value["child_pid"] <= 1:
        raise VerificationError(f"{label} producer PIDs are not process PIDs")
    if value["parent_pid"] == value["child_pid"]:
        raise VerificationError(f"{label} parent and child PIDs are identical")
    if value.get("fork_no_exec") is not True:
        raise VerificationError(f"{label} does not attest fork-without-exec")
    return dict(value)


def _validate_child_wait(
    value: Any, producer_session: Mapping[str, Any], label: str
) -> dict[str, Any]:
    _require_exact_keys(value, CHILD_WAIT_KEYS, label)
    expected = {
        "waited_pid": producer_session["child_pid"],
        "raw_wait_status": 0,
        "wifexited": True,
        "exit_status": 0,
        "wifsignaled": False,
        "terminating_signal": None,
        "core_dumped": False,
    }
    _require_json_equal(value, expected, label)
    raw_status = value["raw_wait_status"]
    if (
        os.WIFEXITED(raw_status) is not True
        or os.WEXITSTATUS(raw_status) != 0
        or os.WIFSIGNALED(raw_status) is not False
        or (hasattr(os, "WCOREDUMP") and os.WCOREDUMP(raw_status) is not False)
    ):
        raise VerificationError(f"{label} raw wait status is not normal exit zero")
    return dict(value)


def _validate_exit_attestation(
    payload: Mapping[str, Any],
    *,
    label: str,
    expected_schema: str,
    expected_phase: str,
    expected_command: list[str],
    expected_produced_artifacts: Mapping[str, Any],
    expected_producer_source: Mapping[str, Any],
    expected_python_runtime: Mapping[str, Any],
    expected_signature_contract: Mapping[str, Any],
    expected_checks: Mapping[str, Any],
    expected_pass_semantics: str,
) -> dict[str, Any]:
    _require_exact_keys(payload, EXIT_ATTESTATION_KEYS, label)
    created_at = _require_utc_timestamp(
        payload.get("created_at_utc"), f"{label} created_at_utc"
    )
    session = _validate_producer_session(
        payload.get("producer_session"),
        f"{label} producer_session",
        expected_mode=expected_command[-1],
    )
    child_wait = _validate_child_wait(
        payload.get("child_wait"), session, f"{label} child_wait"
    )
    if (
        payload.get("authorization_id") != AUTHORIZATION_ID
        or payload.get("phase") != expected_phase
        or not _strict_json_equal(payload.get("supervisor_command"), expected_command)
        or not _strict_json_equal(
            payload.get("produced_artifacts"), expected_produced_artifacts
        )
        or not _strict_json_equal(
            payload.get("producer_source"), expected_producer_source
        )
        or not _strict_json_equal(
            payload.get("python_runtime"), expected_python_runtime
        )
        or not _strict_json_equal(
            payload.get("signature_contract"), expected_signature_contract
        )
        or not _strict_json_equal(
            expected_signature_contract.get("attestation_schema"),
            _exit_attestation_schema(expected_schema),
        )
        or not _strict_json_equal(payload.get("checks"), expected_checks)
        or payload.get("pass") is not True
        or payload.get("pass_semantics") != expected_pass_semantics
    ):
        raise VerificationError(f"{label} exact envelope contract drift")
    return {
        "created_at_utc": created_at,
        "producer_session": session,
        "child_wait": child_wait,
    }


def _require_artifact_field_types(value: Any, label: str, *, with_stat: bool) -> None:
    expected_keys = STAT_ARTIFACT_KEYS if with_stat else BASIC_ARTIFACT_KEYS
    _require_exact_keys(value, expected_keys, label)
    if not isinstance(value.get("path"), str) or not value["path"].startswith("/"):
        raise VerificationError(f"{label} path is not an absolute string")
    if type(value.get("size_bytes")) is not int or value["size_bytes"] < 0:
        raise VerificationError(f"{label} size_bytes is not a non-negative JSON integer")
    _require_lower_sha256(value.get("sha256"), f"{label} sha256")
    mode = value.get("mode")
    if not isinstance(mode, str) or not mode.startswith("0o"):
        raise VerificationError(f"{label} mode is not an octal string")
    if with_stat:
        for field in ("device", "inode", "link_count", "mtime_ns", "ctime_ns"):
            observed = value.get(field)
            if type(observed) is not int or observed < 0:
                raise VerificationError(
                    f"{label} {field} is not a non-negative JSON integer"
                )


def _require_stat_artifact(value: Any, expected: Mapping[str, Any], label: str) -> None:
    _require_artifact_field_types(value, label, with_stat=True)
    _require_json_equal(value, dict(expected), label)


def _require_basic_artifact(value: Any, expected: Mapping[str, Any], label: str) -> None:
    _require_artifact_field_types(value, label, with_stat=False)
    _require_json_equal(value, _basic_projection(expected), label)


def _bind_artifact_record(
    reader: Any, value: Any, label: str, *, with_stat: bool = False
) -> tuple[int, bytes, dict[str, Any], os.stat_result]:
    _require_artifact_field_types(value, label, with_stat=with_stat)
    path_value = value.get("path")
    if not isinstance(path_value, str):
        raise VerificationError(f"{label} path is not a string")
    descriptor, data, basic, opened = _open_regular(reader, Path(path_value))
    observed = _stat_artifact(basic, opened) if with_stat else basic
    _require_json_equal(value, observed, label)
    return descriptor, data, observed, opened


def _open_regular(
    reader: Any, path: Path
) -> tuple[int, bytes, dict[str, Any], os.stat_result]:
    before = _require_absolute_regular_path(path)
    descriptor, data, record = reader.open(path, include_mode=True)
    opened = os.fstat(descriptor)
    if not _same_stat(before, opened) or record["path"] != str(path):
        raise VerificationError(f"Bound artifact path identity drift: {path}")
    _DIRECT_PATH_LEASES.append((path, opened))
    return descriptor, data, record, opened


def _open_private_metadata(
    reader: Any, path: Path, *, expected_mode: str
) -> tuple[dict[str, Any], os.stat_result]:
    before = _require_absolute_regular_path(path)
    descriptor, record = reader.open_metadata(path, include_mode=True)
    opened = os.fstat(descriptor)
    if (
        not _same_stat(before, opened)
        or opened.st_nlink != 1
        or not _strict_json_equal(
            record, {"path": str(path), "mode": expected_mode}
        )
    ):
        raise VerificationError(
            f"Private-key mode differs from {expected_mode}: {path}"
        )
    _DIRECT_PATH_LEASES.append((path, opened))
    return record, opened


def _stat_artifact(
    basic: Mapping[str, Any], opened: os.stat_result
) -> dict[str, Any]:
    return {
        **dict(basic),
        "device": opened.st_dev,
        "inode": opened.st_ino,
        "link_count": opened.st_nlink,
        "mtime_ns": opened.st_mtime_ns,
        "ctime_ns": opened.st_ctime_ns,
    }


def _basic_projection(record: Mapping[str, Any]) -> dict[str, Any]:
    return {key: record[key] for key in BASIC_ARTIFACT_KEYS}


def _require_mode_and_link_count(
    record: Mapping[str, Any], opened: os.stat_result, label: str, *, mode: str = "0o444"
) -> None:
    if record.get("mode") != mode or opened.st_nlink != 1:
        raise VerificationError(f"{label} mode/link-count contract drift")


def _validate_review(
    payload: Mapping[str, Any],
    *,
    expected_reviewer: str,
    expected_reviewer_agent_id: str,
    reviewed_sources: Mapping[str, Any],
    trusted_keys: Mapping[str, Any],
    artifact: Mapping[str, Any],
) -> dict[str, Any]:
    required = {
        "schema_name",
        "schema_version",
        "reviewer",
        "reviewer_agent_id",
        "verdict",
        "findings_closed",
        "reviewed_sources",
        "trusted_signing_keys",
        "conditions",
        "unresolved_findings",
        "unresolved_conditions",
        "pass",
    }
    if set(payload) != required:
        raise VerificationError(f"Promotion review fields are incomplete: {expected_reviewer}")
    if (
        payload.get("schema_name") != REVIEW_SCHEMA
        or payload.get("schema_version") != SCHEMA_VERSION
        or payload.get("reviewer") != expected_reviewer
        or payload.get("reviewer_agent_id") != expected_reviewer_agent_id
        or payload.get("verdict") != "APPROVE"
        or payload.get("findings_closed") is not True
        or not _strict_json_equal(payload.get("reviewed_sources"), reviewed_sources)
        or not _strict_json_equal(
            payload.get("trusted_signing_keys"), trusted_keys
        )
        or not _strict_json_equal(payload.get("conditions"), [])
        or not _strict_json_equal(payload.get("unresolved_findings"), [])
        or not _strict_json_equal(payload.get("unresolved_conditions"), [])
        or payload.get("pass") is not True
    ):
        raise VerificationError(f"Promotion review did not exactly approve v3: {expected_reviewer}")
    return {
        "artifact": dict(artifact),
        "reviewer": expected_reviewer,
        "reviewer_agent_id": expected_reviewer_agent_id,
        "verdict": "APPROVE",
        "findings_closed": True,
    }


def _validate_execution_binding(
    value: Any,
    promotion_record: Mapping[str, Any],
    promotion_data: bytes,
) -> None:
    _require_exact_keys(
        value,
        frozenset(
            {
                "source_sha256",
                "function_and_method_count",
                "inventory_sha256",
                "module_code_sha256",
                "critical_constants_equal",
                "canonical_commands_equal",
                "pass",
            }
        ),
        "promotion execution binding",
    )
    inventory_sha256 = value.get("inventory_sha256")
    module_code_sha256 = value.get("module_code_sha256")
    promotion_filename = str(PROMOTION_TOOL)
    shadow: dict[str, Any] = {
        "__name__": "promotion_v2_independent_binding_recount",
        "__file__": promotion_filename,
        "__package__": "",
    }
    promotion_module_code = compile(promotion_data, promotion_filename, "exec")
    exec(promotion_module_code, shadow)
    independent_inventory = _code_inventory(shadow, promotion_filename)
    independent_inventory_sha256 = hashlib.sha256(
        json.dumps(independent_inventory, sort_keys=True).encode("utf-8")
    ).hexdigest()
    independent_module_sha256 = _stable_code_digest(promotion_module_code)
    if (
        value.get("source_sha256") != promotion_record.get("sha256")
        or value.get("source_sha256") != hashlib.sha256(promotion_data).hexdigest()
        or type(value.get("function_and_method_count")) is not int
        or value["function_and_method_count"] <= 0
        or value["function_and_method_count"] != len(independent_inventory)
        or not isinstance(inventory_sha256, str)
        or len(inventory_sha256) != 64
        or any(character not in "0123456789abcdef" for character in inventory_sha256)
        or not isinstance(module_code_sha256, str)
        or len(module_code_sha256) != 64
        or any(character not in "0123456789abcdef" for character in module_code_sha256)
        or inventory_sha256 != independent_inventory_sha256
        or module_code_sha256 != independent_module_sha256
        or value.get("critical_constants_equal") is not True
        or value.get("canonical_commands_equal") is not True
        or value.get("pass") is not True
    ):
        raise VerificationError("Promotion execution/source binding contract drift")


def _historical_incident_disclosure(
    archive_record: Mapping[str, Any], archive_payload: Mapping[str, Any]
) -> dict[str, Any]:
    expected_check_keys = {
        "canonical_cache_path_absent_after_archive",
        "canonical_log_path_absent_after_archive",
        "cn_ccf_output_absent_after_failure",
        "completion_runner_bytes_unchanged",
        "failed_log_preserved",
        "failed_python_cache_preserved",
        "final_dataset_output_absent_after_failure",
        "final_report_output_absent_after_failure",
        "matched_normal_outputs_absent_after_failure",
        "no_failed_artifact_deleted",
        "strict_output_absent_after_failure",
    }
    checks = archive_payload.get("checks")
    if (
        archive_record.get("sha256") != EXPECTED_INCIDENT_ARCHIVE_SHA256
        or archive_record.get("mode") != "0o444"
        or archive_payload.get("schema_name")
        != "intersubmod.completion_attempt_archive_receipt"
        or archive_payload.get("schema_version") != "1.0.0"
        or archive_payload.get("task_type") != "B_comprehensive_validation"
        or archive_payload.get("attempt")
        != "m2v5_completion_attempt2_mode_tightening_conflict"
        or archive_payload.get("status")
        != "FAILED_BEFORE_DOWNSTREAM_OUTPUT_AND_ARCHIVED"
        or archive_payload.get("failure_class")
        != "post_execution_source_mode_tightening_conflicts_with_full_stat_identity_equality"
        or not isinstance(checks, Mapping)
        or set(checks) != expected_check_keys
        or any(value is not True for value in checks.values())
        or archive_payload.get("pass") is not True
    ):
        raise VerificationError("Historical completion-attempt incident evidence drift")
    return {
        "incident_id": "20260722_completion_attempt2_mode_tightening_conflict",
        "archive_receipt": dict(archive_record),
        "failure_class": archive_payload["failure_class"],
        "failed_before_downstream_output": True,
        "failed_evidence_preserved": True,
        "original_completion_runner_end_to_end_pass_claimed": False,
        "original_verifier_reexecuted_successfully_claimed": False,
        "scientific_payload_changed": False,
    }


def _validate_signed_authority_evidence(
    authority_module: ModuleType,
    reader: Any,
    openssl_fd: int,
    value: Any,
    *,
    expected_authority_id: str,
    label: str,
) -> None:
    _require_exact_keys(value, SIGNED_AUTHORITY_EVIDENCE_KEYS, label)
    authority_fd, authority_data, authority_record, authority_stat = _bind_artifact_record(
        reader, value["authority"], f"{label} authority"
    )
    approval_fd, approval_data, approval_record, approval_stat = _bind_artifact_record(
        reader, value["approval"], f"{label} approval"
    )
    signature_fd, _, signature_record, signature_stat = _bind_artifact_record(
        reader, value["signature"], f"{label} signature"
    )
    public_fd, _, public_record, public_stat = _bind_artifact_record(
        reader, value["public_key"], f"{label} public key"
    )
    expected_public_path, expected_public_sha256 = {
        EXPECTED_V3_AUTHORITY_ID: (
            EXPECTED_V3_PUBLIC_KEY_PATH,
            EXPECTED_V3_PUBLIC_KEY_SHA256,
        ),
        EXPECTED_V7_AUTHORITY_ID: (
            EXPECTED_V7_PUBLIC_KEY_PATH,
            EXPECTED_V7_PUBLIC_KEY_SHA256,
        ),
    }[expected_authority_id]
    for artifact_label, record, opened in (
        (f"{label} authority", authority_record, authority_stat),
        (f"{label} approval", approval_record, approval_stat),
        (f"{label} signature", signature_record, signature_stat),
        (f"{label} public key", public_record, public_stat),
    ):
        _require_mode_and_link_count(record, opened, artifact_label)
    if signature_record["size_bytes"] != 64:
        raise VerificationError(f"{label} detached signature size drift")

    private_value = value["private_key"]
    _require_exact_keys(private_value, frozenset({"path", "mode"}), f"{label} private key")
    if not isinstance(private_value.get("path"), str):
        raise VerificationError(f"{label} private key path is not a string")
    private_record, _ = _open_private_metadata(
        reader, Path(private_value["path"]), expected_mode="0o0"
    )
    _require_json_equal(private_value, private_record, f"{label} private key")

    authority_payload = _parse_json(authority_data, f"{label} authority")
    approval_payload = _parse_json(approval_data, f"{label} approval")
    if (
        value.get("authority_id") != expected_authority_id
        or value.get("signature_verified") is not True
        or public_record.get("path") != str(expected_public_path)
        or public_record.get("sha256") != expected_public_sha256
        or authority_payload.get("schema_name")
        != "intersubmod.release_source_authority"
        or authority_payload.get("schema_version") != "1.0.0"
        or authority_payload.get("authority_id") != expected_authority_id
        or authority_payload.get("approval_status")
        != "APPROVED_FOR_FULL_TASK_B_RUN"
        or approval_payload.get("schema_name")
        != "intersubmod.release_source_authority.approval"
        or approval_payload.get("schema_version") != "1.0.0"
        or approval_payload.get("authority_id") != expected_authority_id
        or approval_payload.get("approval_status")
        != "APPROVED_FOR_FULL_TASK_B_RUN"
        or not _strict_json_equal(
            approval_payload.get("authority_manifest"), authority_record
        )
        or not _strict_json_equal(approval_payload.get("public_key"), public_record)
        or approval_payload.get("private_key_retirement", {}).get("path")
        != private_record["path"]
        or private_record["path"]
        != str(
            expected_public_path.with_name(
                "ed25519_private_encrypted_unrecoverable_after_signing.pem"
            )
        )
        or approval_payload.get("private_key_retirement", {}).get(
            "required_post_signature_mode"
        )
        != "0o0"
    ):
        raise VerificationError(f"{label} signed authority contract drift")
    if not authority_module.verify_ed25519_signature_fds(
        openssl_fd=openssl_fd,
        data_fd=approval_fd,
        public_key_fd=public_fd,
        signature_fd=signature_fd,
    ):
        raise VerificationError(f"{label} detached signature verification failed")
    del authority_fd


def _validate_authorization_evidence(
    authority_module: ModuleType,
    reader: Any,
    openssl_fd: int,
    evidence: Any,
    *,
    archive_record: Mapping[str, Any],
    python_record: Mapping[str, Any],
    v7_validation: Mapping[str, Any],
    lease: Mapping[str, Any],
) -> dict[str, Any]:
    _require_exact_keys(evidence, AUTHORIZATION_EVIDENCE_KEYS, "authorization evidence")

    _, _, snapshot_record, _ = _bind_artifact_record(
        reader, evidence["source_snapshot"], "evidence source snapshot"
    )
    _, _, manifest_record, _ = _bind_artifact_record(
        reader, evidence["tumor_ref_manifest"], "evidence tumor-REF manifest"
    )
    _, prior_data, prior_record, _ = _bind_artifact_record(
        reader, evidence["prior_policy_review"], "evidence prior policy review"
    )
    _, _, builder_record, _ = _bind_artifact_record(
        reader, evidence["final_builder"], "evidence final builder"
    )
    _require_json_equal(
        evidence["python_runtime"], python_record, "evidence Python runtime"
    )
    if (
        "正式 receipt 仍依設計等待 completion runner 建立，不能以 pre-probe 取代".encode(
            "utf-8"
        )
        not in prior_data
        or "pre-probe 明確不可替代".encode("utf-8") not in prior_data
    ):
        raise VerificationError("Prior policy evidence no longer contains superseded text")

    transition = evidence["focal_source_identity_transition"]
    _require_exact_keys(
        transition,
        frozenset(
            {
                "during_execution",
                "current",
                "unchanged_fields",
                "observed_differences",
                "interpretation",
            }
        ),
        "focal source identity transition",
    )
    _require_exact_keys(
        transition["during_execution"],
        STAT_ARTIFACT_KEYS,
        "focal source during-execution identity",
    )
    _require_artifact_field_types(
        transition["during_execution"],
        "focal source during-execution identity",
        with_stat=True,
    )
    _, _, focal_current, _ = _bind_artifact_record(
        reader,
        transition["current"],
        "focal source current identity",
        with_stat=True,
    )
    unchanged_fields = ["path", "size_bytes", "sha256", "device", "inode", "mtime_ns"]
    during = transition["during_execution"]
    expected_differences = {
        "mode": {"during": "0o664", "current": "0o444"},
        "ctime_ns": {"during": during["ctime_ns"], "current": focal_current["ctime_ns"]},
    }
    if (
        not _strict_json_equal(transition["unchanged_fields"], unchanged_fields)
        or any(
            not _strict_json_equal(during[field], focal_current[field])
            for field in unchanged_fields
        )
        or during["mode"] != "0o664"
        or focal_current["mode"] != "0o444"
        or during["ctime_ns"] == focal_current["ctime_ns"]
        or not _strict_json_equal(
            transition["observed_differences"], expected_differences
        )
        or transition["interpretation"]
        != (
            "Observed bytes, path, device, inode, size, and mtime are unchanged; mode was "
            "later tightened from 0664 to 0444. This does not prove that no unobserved "
            "transient mutation ever occurred."
        )
    ):
        raise VerificationError("Focal source mode-only transition contract drift")

    archive = evidence["failed_attempt_archive"]
    _require_exact_keys(
        archive,
        frozenset({"receipt", "archived_log", "archived_cache"}),
        "failed-attempt archive evidence",
    )
    _require_json_equal(archive["receipt"], archive_record, "archive receipt cross-link")
    _bind_artifact_record(reader, archive["archived_log"], "failed-attempt archived log")
    _require_json_equal(
        archive["archived_cache"],
        {
            "path": str(
                WORKSPACE_ROOT
                / "audit_archives"
                / "20260722_completion_attempt2_mode_tightening_conflict"
                / "python_cache_m2v5_completion_v2_bound_bootstrap.failed"
            ),
            "regular_file_count": 45,
            "regular_file_total_bytes": 1_463_182,
            "sha256sum_line_manifest_sha256": (
                "0084670ad8f6440715a110c90aa399fc26573306c4e7a19ae40b9782d23264e2"
            ),
            "all_files_bound_until_process_exit": True,
        },
        "failed-attempt archived cache contract",
    )

    _validate_signed_authority_evidence(
        authority_module,
        reader,
        openssl_fd,
        evidence["historical_v3_authority"],
        expected_authority_id=EXPECTED_V3_AUTHORITY_ID,
        label="historical v3 authority",
    )
    expected_v7_authority = {
        "authority": v7_validation["authority_manifest"],
        "approval": v7_validation["detached_approval"],
        "signature": v7_validation["detached_approval_signature"],
        "public_key": v7_validation["external_public_key"],
        "private_key": {
            "path": v7_validation["retired_private_key"]["path"],
            "mode": v7_validation["retired_private_key"]["mode"],
        },
        "authority_id": EXPECTED_V7_AUTHORITY_ID,
        "signature_verified": True,
    }
    _require_json_equal(
        evidence["current_v7_authority"],
        expected_v7_authority,
        "current v7 authority evidence",
    )
    expected_runtime = {
        "authority_id": EXPECTED_V7_AUTHORITY_ID,
        "source_set_sha256": EXPECTED_V7_SOURCE_SET_SHA256,
        "required_roles_verified": v7_validation["required_roles_verified"],
        "binding_lease": dict(lease),
        "pass": True,
    }
    _require_json_equal(
        evidence["current_v7_runtime_validation"],
        expected_runtime,
        "current v7 runtime validation evidence",
    )

    v7_authority_path = Path(expected_v7_authority["authority"]["path"])
    _, v7_data, _, _ = _bind_artifact_record(
        reader, expected_v7_authority["authority"], "current v7 authority manifest"
    )
    v7_payload = _parse_json(v7_data, "current v7 authority manifest")
    del v7_authority_path
    v7_sources = v7_payload.get("sources")
    if (
        not isinstance(v7_sources, Mapping)
        or not _strict_json_equal(v7_sources.get("final_dataset_builder"), builder_record)
        or not _strict_json_equal(
            v7_sources.get("focal_alt_cluster_lib"), _basic_projection(focal_current)
        )
    ):
        raise VerificationError("Authorization evidence is not source-bound by v7")

    _require_json_equal(
        evidence["python_cache_contract"],
        {
            "paths": {key: str(value) for key, value in PROMOTION_CACHE_ROOTS.items()},
            "mode": "0o700",
            "active_mode_cache_fresh_and_empty_before_and_after_import": True,
            "bytecode_writes_disabled": True,
        },
        "promotion Python cache contract",
    )
    _require_json_equal(
        evidence["final_builder_gate_execution"],
        {
            "constants": [
                "REPO_ROOT",
                "SCHEMA_VERSION",
                "SCRIPT_DIR",
                "TUMOR_REF_SOURCE_IDENTITY_REQUIRED_ROLES",
                "TUMOR_REF_SOURCE_IDENTITY_SCHEMA_VERSION",
                "TUMOR_REF_SOURCE_IDENTITY_TRUE_CHECKS",
                "TUMOR_REF_SOURCE_IDENTITY_VERIFIER",
            ],
            "definitions": [
                "ContractError",
                "artifact",
                "load_json",
                "load_tumor_ref_source_identity_attestation",
                "parse_utc_timestamp",
                "require_nonempty_file",
                "require_schema_pass",
                "sha256",
                "source_script_token_matches_attested_source",
                "verify_declared_artifact",
                "verify_declared_path_artifact",
            ],
            "gate_first_line": 4729,
            "no_local_module_imports_executed": True,
        },
        "final builder bound-AST gate execution contract",
    )
    _require_json_equal(
        evidence["final_builder_source_path_gate"],
        {
            "status": "VERIFIED_BOUNDED_RETROSPECTIVE_SOURCE_IDENTITY",
            "release_gate_pass": True,
            "publishable_task_b_release": True,
        },
        "final builder source-path gate",
    )
    return {
        "source_snapshot": snapshot_record,
        "tumor_ref_manifest": manifest_record,
        "prior_policy_review": prior_record,
        "final_builder": builder_record,
        "current_v7_authority": expected_v7_authority,
        "current_v7_runtime_validation": expected_runtime,
    }


def _encode_json(payload: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def _link_fd_no_replace(source_fd: int, target_directory_fd: int, target_name: str) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    linkat = getattr(libc, "linkat", None)
    if linkat is None:
        raise VerificationError("linkat is required for terminal commit publication")
    linkat.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
    )
    linkat.restype = ctypes.c_int
    result = linkat(
        AT_FDCWD,
        os.fsencode(f"/proc/self/fd/{source_fd}"),
        target_directory_fd,
        os.fsencode(target_name),
        AT_SYMLINK_FOLLOW,
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise VerificationError(
            "Terminal commit publication failed: " + os.strerror(error_number)
        )


def _write_exclusive_readonly(
    path: Path,
    data: bytes,
    *,
    precommit: Callable[[], None] | None = None,
    terminal_commit: Callable[[], None] | None = None,
) -> dict[str, Any]:
    if path.resolve(strict=False) != path:
        raise VerificationError(f"Verification receipt path contains a symlink: {path}")
    if not hasattr(os, "O_TMPFILE"):
        raise VerificationError("O_TMPFILE is required for terminal commit publication")
    parent_stat = os.stat(path.parent, follow_symlinks=False)
    if not stat.S_ISDIR(parent_stat.st_mode):
        raise VerificationError("Verification receipt parent is not a directory")
    directory_flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0)
    if hasattr(os, "O_NOFOLLOW"):
        directory_flags |= os.O_NOFOLLOW
    parent_fd = os.open(path.parent, directory_flags)
    descriptor = -1
    reopened_fd = -1
    durable_fd = -1
    record: dict[str, Any] | None = None
    try:
        if not _same_stat(parent_stat, os.fstat(parent_fd)):
            raise VerificationError("Verification receipt parent path rebound before open")
        descriptor = os.open(
            ".", os.O_RDWR | os.O_TMPFILE | os.O_CLOEXEC, 0o400, dir_fd=parent_fd
        )
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise VerificationError("Short write while publishing verification receipt")
            offset += written
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
        observed = os.fstat(descriptor)
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 0
            or stat.S_IMODE(observed.st_mode) != 0o444
            or observed.st_size != len(data)
            or os.pread(descriptor, len(data), 0) != data
        ):
            raise VerificationError("Verification receipt publication check failed")
        if precommit is not None:
            precommit()
        if not _same_inode_and_mode(parent_stat, os.fstat(parent_fd)):
            raise VerificationError("Verification receipt parent path rebound before commit")
        if os.path.lexists(path):
            raise VerificationError("Verification receipt commit target is occupied")
        record = {
            "path": str(path),
            "size_bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
            "mode": "0o444",
        }
        _link_fd_no_replace(descriptor, parent_fd, path.name)
        committed = os.fstat(descriptor)
        linked = os.stat(path.name, dir_fd=parent_fd, follow_symlinks=False)
        reopened_flags = os.O_RDONLY | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            reopened_flags |= os.O_NOFOLLOW
        reopened_fd = os.open(path.name, reopened_flags, dir_fd=parent_fd)
        reopened = os.fstat(reopened_fd)
        if (
            committed.st_nlink != 1
            or not _same_stat(committed, linked)
            or not _same_stat(committed, reopened)
            or os.pread(reopened_fd, reopened.st_size, 0) != data
        ):
            raise VerificationError(
                "Verification receipt post-link path/inode binding failed"
            )
        os.fsync(parent_fd)
        parent_after = os.stat(path.parent, follow_symlinks=False)
        if (
            not _same_inode_and_mode(parent_stat, os.fstat(parent_fd))
            or not _same_inode_and_mode(parent_stat, parent_after)
            or path.resolve(strict=True) != path
        ):
            raise VerificationError(
                "Verification receipt parent/path rebound after durable commit"
            )
        durable_fd = os.open(path.name, reopened_flags, dir_fd=parent_fd)
        durable_source = os.fstat(descriptor)
        durable_opened = os.fstat(durable_fd)
        durable_linked = os.stat(path.name, dir_fd=parent_fd, follow_symlinks=False)
        if (
            durable_source.st_nlink != 1
            or not _same_stat(durable_source, durable_opened)
            or not _same_stat(durable_source, durable_linked)
            or os.pread(durable_fd, durable_opened.st_size, 0) != data
        ):
            raise VerificationError(
                "Verification receipt post-fsync path/inode binding failed"
            )
        if terminal_commit is not None:
            terminal_commit()
            raise VerificationError("Verification terminal commit callback returned")
    finally:
        if durable_fd >= 0:
            os.close(durable_fd)
        if reopened_fd >= 0:
            os.close(reopened_fd)
        if descriptor >= 0:
            os.close(descriptor)
        os.close(parent_fd)
    if record is None:
        raise VerificationError("Verification receipt publication did not complete")
    return record


def _verification_payload(
    *,
    mode: str,
    authorization_record: Mapping[str, Any],
    prepare_exit_attestation_record: Mapping[str, Any],
    authorization_signature_record: Mapping[str, Any],
    completion_record: Mapping[str, Any],
    promote_exit_attestation_record: Mapping[str, Any],
    completion_signature_record: Mapping[str, Any],
    source_record: Mapping[str, Any],
    canonical_record: Mapping[str, Any],
    trusted_keys: Mapping[str, Any],
    authorization_private_record: Mapping[str, Any],
    completion_private_record: Mapping[str, Any],
    incident_disclosure: Mapping[str, Any],
    release_authority_evidence: Mapping[str, Any],
    retained_descriptor_count: int,
) -> dict[str, Any]:
    return {
        "schema_name": VERIFICATION_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "task_type": "B_comprehensive_validation",
        "mode": mode.removeprefix("--"),
        "authorization": {
            "artifact": dict(authorization_record),
            "signature": dict(authorization_signature_record),
            "signature_verified": True,
        },
        "completion": {
            "artifact": dict(completion_record),
            "signature": dict(completion_signature_record),
            "signature_verified": True,
        },
        "producer_exit_attestations": {
            "prepare": {
                "artifact": dict(prepare_exit_attestation_record),
                "signature": dict(authorization_signature_record),
                "public_key": dict(trusted_keys["authorization_public_key"]),
                "signature_verified": True,
            },
            "promote": {
                "artifact": dict(promote_exit_attestation_record),
                "signature": dict(completion_signature_record),
                "public_key": dict(trusted_keys["completion_public_key"]),
                "signature_verified": True,
            },
        },
        "source_receipt": dict(source_record),
        "canonical_receipt": dict(canonical_record),
        "trusted_signing_keys": dict(trusted_keys),
        "private_key_retirement": {
            "authorization": dict(authorization_private_record),
            "completion": dict(completion_private_record),
        },
        "historical_incident_disclosure": dict(incident_disclosure),
        "release_source_authority_validator": dict(release_authority_evidence),
        "checks": {
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
        },
        "retained_critical_input_descriptor_count": retained_descriptor_count,
        "governance": {
            "original_completion_runner_end_to_end_pass_claimed": False,
            "original_verifier_reexecuted_successfully_claimed": False,
            "runner_only_gate_replay_evidence_required_before_downstream": True,
            "signed_downstream_continuation_revalidation_required": True,
            "verification_scope": "promotion_integrity_and_continuation_gate_only",
        },
        "pass": True,
        "pass_semantics": (
            "signed byte-identical promotion chain verified; this does not claim the original "
            "completion runner or original verifier completed successfully"
        ),
    }


def _verify(mode: str) -> tuple[dict[str, Any], dict[str, Any] | None]:
    _require_absent(
        LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
        "legacy authorization-payload signature",
    )
    _require_absent(
        LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
        "legacy completion-payload signature",
    )
    if os.path.lexists(PROMOTION_INCIDENT):
        raise VerificationError("Promotion incident exists; continuation is permanently refused")
    if any(os.path.lexists(path) for path in DOWNSTREAM_OUTPUTS):
        raise VerificationError("Downstream output exists before continuation verification")
    authority_module, bootstrap_authority_record = _load_authority_module()
    v7_validation = authority_module.validate_release_source_authority(
        authority_module.EXPECTED_SOURCE_PATHS
    )
    lease = v7_validation.get("validation_binding_lease")
    if (
        v7_validation.get("pass") is not True
        or v7_validation.get("authority_id") != EXPECTED_V7_AUTHORITY_ID
        or v7_validation.get("source_set_sha256")
        != EXPECTED_V7_SOURCE_SET_SHA256
        or not isinstance(v7_validation.get("all_source_roles_verified"), list)
        or len(v7_validation["all_source_roles_verified"]) != 23
        or not isinstance(lease, Mapping)
        or lease.get("policy")
        != "BOUND_FILE_DESCRIPTORS_RETAINED_UNTIL_PROCESS_EXIT"
        or type(lease.get("descriptor_count")) is not int
        or lease["descriptor_count"] <= 0
    ):
        raise VerificationError("v7 release source authority validation drift")
    with authority_module.BoundArtifactReader() as reader:
        self_fd, self_data, self_record, self_stat = _open_regular(reader, VERIFIER)
        del self_fd
        if (
            _INHERITED_VERIFIER_FD_LEASE is not None
            and not _same_stat(self_stat, _INHERITED_VERIFIER_FD_LEASE[1])
        ):
            raise VerificationError(
                "Bound verifier source differs from the inherited execution descriptor"
            )
        _require_executing_verifier_matches_bound_source(self_data)
        authority_fd, authority_bytes, authority_record, authority_stat = _open_regular(
            reader, RELEASE_SOURCE_AUTHORITY
        )
        del authority_fd
        if not _strict_json_equal(
            authority_record, bootstrap_authority_record
        ) or hashlib.sha256(authority_bytes).hexdigest() != (
            EXPECTED_RELEASE_SOURCE_AUTHORITY_SHA256
        ):
            raise VerificationError("Bound authority module differs from executed source bytes")
        _require_mode_and_link_count(self_record, self_stat, "continuation verifier")
        _require_mode_and_link_count(
            authority_record, authority_stat, "release source authority"
        )

        openssl_fd, _, openssl_record, openssl_stat = _open_regular(reader, OPENSSL)
        if (
            openssl_record["sha256"] != EXPECTED_OPENSSL_SHA256
            or openssl_record["mode"] != "0o755"
            or openssl_stat.st_nlink != 1
        ):
            raise VerificationError("OpenSSL identity drift")

        python_target = PYTHON.resolve(strict=True)
        if python_target != PRIMARY_PYTHON_RUNTIME:
            raise VerificationError("Primary Python launcher target drift")
        python_fd, _, python_record, python_stat = _open_regular(reader, python_target)
        del python_fd
        process_exe = os.stat("/proc/self/exe")
        if (
            python_stat.st_dev != process_exe.st_dev
            or python_stat.st_ino != process_exe.st_ino
            or python_record["sha256"] != EXPECTED_PYTHON_SHA256
            or python_record["mode"] != "0o775"
        ):
            raise VerificationError("Running interpreter is not the canonical Python target")

        source_fd, source_data, source_basic, source_stat = _open_regular(
            reader, SOURCE_RECEIPT
        )
        del source_fd
        canonical_fd, canonical_data, canonical_basic, canonical_stat = _open_regular(
            reader, CANONICAL_RECEIPT
        )
        del canonical_fd
        source_record = _stat_artifact(source_basic, source_stat)
        canonical_record = _stat_artifact(canonical_basic, canonical_stat)
        for label, data, record, opened in (
            ("source receipt", source_data, source_record, source_stat),
            ("canonical receipt", canonical_data, canonical_record, canonical_stat),
        ):
            if (
                len(data) != EXPECTED_SOURCE_SIZE
                or record["size_bytes"] != EXPECTED_SOURCE_SIZE
                or record["sha256"] != EXPECTED_SOURCE_SHA256
                or record["mode"] != "0o444"
                or opened.st_nlink != 1
            ):
                raise VerificationError(f"{label} fixed identity drift")
        if (
            source_data != canonical_data
            or (source_stat.st_dev, source_stat.st_ino)
            == (canonical_stat.st_dev, canonical_stat.st_ino)
        ):
            raise VerificationError(
                "Canonical receipt is not a byte-exact independent-inode source copy"
            )
        source_payload = _parse_json(source_data, "source preprobe receipt")
        if (
            source_payload.get("schema_name")
            != "intersubmod.retrospective_running_source_identity.receipt"
            or source_payload.get("schema_version") != "1.2.0"
            or source_payload.get("task_type") != "B_comprehensive_validation"
            or source_payload.get("pass") is not True
        ):
            raise VerificationError("Source preprobe receipt semantic contract drift")

        auth_fd, auth_data, auth_record, auth_stat = _open_regular(reader, AUTHORIZATION)
        prepare_exit_fd, prepare_exit_data, prepare_exit_record, prepare_exit_stat = (
            _open_regular(reader, PREPARE_EXIT_ATTESTATION)
        )
        auth_sig_fd, _, auth_sig_record, auth_sig_stat = _open_regular(
            reader, AUTHORIZATION_SIGNATURE
        )
        auth_public_fd, _, auth_public_record, auth_public_stat = _open_regular(
            reader, AUTHORIZATION_PUBLIC_KEY
        )
        auth_private_record, auth_private_stat = _open_private_metadata(
            reader, AUTHORIZATION_PRIVATE_KEY, expected_mode="0o0"
        )
        completion_fd, completion_data, completion_record, completion_stat = _open_regular(
            reader, COMPLETION
        )
        promote_exit_fd, promote_exit_data, promote_exit_record, promote_exit_stat = (
            _open_regular(reader, PROMOTE_EXIT_ATTESTATION)
        )
        completion_sig_fd, _, completion_sig_record, completion_sig_stat = _open_regular(
            reader, COMPLETION_SIGNATURE
        )
        completion_public_fd, _, completion_public_record, completion_public_stat = (
            _open_regular(reader, COMPLETION_PUBLIC_KEY)
        )
        completion_private_record, completion_private_stat = _open_private_metadata(
            reader, COMPLETION_PRIVATE_KEY, expected_mode="0o0"
        )
        continuation_public_fd, _, continuation_public_record, continuation_public_stat = (
            _open_regular(reader, CONTINUATION_PUBLIC_KEY)
        )
        del continuation_public_fd
        continuation_private_record, continuation_private_stat = _open_private_metadata(
            reader, CONTINUATION_PRIVATE_KEY, expected_mode="0o400"
        )

        for label, record, opened in (
            ("authorization", auth_record, auth_stat),
            ("prepare exit attestation", prepare_exit_record, prepare_exit_stat),
            ("authorization signature", auth_sig_record, auth_sig_stat),
            ("authorization public key", auth_public_record, auth_public_stat),
            ("completion", completion_record, completion_stat),
            ("promote exit attestation", promote_exit_record, promote_exit_stat),
            ("completion signature", completion_sig_record, completion_sig_stat),
            ("completion public key", completion_public_record, completion_public_stat),
            (
                "continuation public key",
                continuation_public_record,
                continuation_public_stat,
            ),
        ):
            _require_mode_and_link_count(record, opened, label)
        if auth_sig_record["size_bytes"] != 64 or completion_sig_record["size_bytes"] != 64:
            raise VerificationError("Ed25519 detached signature size drift")
        if (
            auth_public_record["size_bytes"] != 113
            or auth_public_record["sha256"]
            != EXPECTED_AUTHORIZATION_PUBLIC_KEY_SHA256
            or completion_public_record["size_bytes"] != 113
            or completion_public_record["sha256"]
            != EXPECTED_COMPLETION_PUBLIC_KEY_SHA256
            or continuation_public_record["size_bytes"] != 113
            or continuation_public_record["sha256"]
            != EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256
            or auth_public_record["sha256"] == completion_public_record["sha256"]
            or len(
                {
                    auth_public_record["sha256"],
                    completion_public_record["sha256"],
                    continuation_public_record["sha256"],
                }
            )
            != 3
            or (auth_public_stat.st_dev, auth_public_stat.st_ino)
            == (completion_public_stat.st_dev, completion_public_stat.st_ino)
            or len(
                {
                    (auth_public_stat.st_dev, auth_public_stat.st_ino),
                    (completion_public_stat.st_dev, completion_public_stat.st_ino),
                    (
                        continuation_public_stat.st_dev,
                        continuation_public_stat.st_ino,
                    ),
                }
            )
            != 3
            or stat.S_IMODE(auth_private_stat.st_mode) != 0
            or stat.S_IMODE(completion_private_stat.st_mode) != 0
            or stat.S_IMODE(continuation_private_stat.st_mode) != 0o400
            or continuation_private_stat.st_nlink != 1
        ):
            raise VerificationError("Pinned signing key or retirement contract drift")

        if not authority_module.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=prepare_exit_fd,
            public_key_fd=auth_public_fd,
            signature_fd=auth_sig_fd,
        ):
            raise VerificationError(
                "Prepare exit-attestation Ed25519 signature verification failed"
            )
        if not authority_module.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=promote_exit_fd,
            public_key_fd=completion_public_fd,
            signature_fd=completion_sig_fd,
        ):
            raise VerificationError(
                "Promote exit-attestation Ed25519 signature verification failed"
            )

        promotion_fd, promotion_data, promotion_record, promotion_stat = _open_regular(
            reader, PROMOTION_TOOL
        )
        del promotion_fd
        replayer_fd, _, replayer_record, replayer_stat = _open_regular(
            reader, RUNNER_GATE_REPLAYER
        )
        del replayer_fd
        continuation_fd, _, continuation_record, continuation_stat = _open_regular(
            reader, DOWNSTREAM_CONTINUATION
        )
        del continuation_fd
        for label, record, opened in (
            ("promotion tool", promotion_record, promotion_stat),
            ("runner gate replayer", replayer_record, replayer_stat),
            ("downstream continuation", continuation_record, continuation_stat),
        ):
            _require_mode_and_link_count(record, opened, label)
        reviewed_runtime_sources: dict[str, dict[str, Any]] = {}
        for role, (path, expected_sha256, expected_mode) in (
            REVIEWED_RUNTIME_SOURCE_CONTRACTS.items()
        ):
            _, _, record, opened = _open_regular(reader, path)
            _require_mode_and_link_count(
                record, opened, f"reviewed runtime source {role}", mode=expected_mode
            )
            if record.get("sha256") != expected_sha256:
                raise VerificationError(f"Reviewed runtime source SHA drift: {role}")
            reviewed_runtime_sources[role] = record
        observed_asset_parts = tuple(
            sorted(
                path
                for path in PORTABLE_ASSET_ROOT.iterdir()
                if path.name.startswith("portable-artifact-reader.html.gz.b64.part")
            )
        )
        if observed_asset_parts != PORTABLE_READER_ASSET_PARTS:
            raise VerificationError("Portable reader asset part membership drift")
        reviewed_sources = {
            "promotion_tool": promotion_record,
            "continuation_verifier": self_record,
            "runner_gate_replay": replayer_record,
            "downstream_continuation": continuation_record,
            **reviewed_runtime_sources,
        }
        trusted_keys = {
            "authorization_public_key": auth_public_record,
            "completion_public_key": completion_public_record,
            "continuation_public_key": continuation_public_record,
        }

        review_records: dict[str, dict[str, Any]] = {}
        reviewer_ids: set[str] = set()
        review_inodes: set[tuple[int, int]] = set()
        for role, expected_reviewer, expected_agent_id, review_path in REVIEW_SLOTS:
            _, review_data, review_record, review_stat = _open_regular(reader, review_path)
            _require_mode_and_link_count(
                review_record, review_stat, f"promotion review {expected_reviewer}"
            )
            review_payload = _parse_json(review_data, f"promotion review {expected_reviewer}")
            approval = _validate_review(
                review_payload,
                expected_reviewer=expected_reviewer,
                expected_reviewer_agent_id=expected_agent_id,
                reviewed_sources=reviewed_sources,
                trusted_keys=trusted_keys,
                artifact=review_record,
            )
            reviewer_id = approval["reviewer_agent_id"]
            inode_key = (review_stat.st_dev, review_stat.st_ino)
            if reviewer_id in reviewer_ids or inode_key in review_inodes:
                raise VerificationError("Promotion reviewer identities are not distinct")
            reviewer_ids.add(reviewer_id)
            review_inodes.add(inode_key)
            review_records[role] = approval

        _, archive_data, archive_record, archive_stat = _open_regular(
            reader, INCIDENT_ARCHIVE_RECEIPT
        )
        _require_mode_and_link_count(
            archive_record, archive_stat, "historical incident archive receipt"
        )
        archive_payload = _parse_json(archive_data, "historical incident archive receipt")
        incident_disclosure = _historical_incident_disclosure(
            archive_record, archive_payload
        )

        authorization_payload = _parse_json(auth_data, "promotion authorization")
        _require_exact_keys(
            authorization_payload, AUTHORIZATION_KEYS, "promotion authorization"
        )
        authorization_created_at = _require_utc_timestamp(
            authorization_payload.get("created_at_utc"),
            "promotion authorization created_at_utc",
        )
        validated_evidence = _validate_authorization_evidence(
            authority_module,
            reader,
            openssl_fd,
            authorization_payload.get("evidence"),
            archive_record=archive_record,
            python_record=python_record,
            v7_validation=v7_validation,
            lease=lease,
        )
        del validated_evidence
        _validate_execution_binding(
            authorization_payload.get("execution_binding"),
            promotion_record,
            promotion_data,
        )
        expected_commands = {
            "prepare_authorization": _promotion_command("--prepare-authorization"),
            "promote": _promotion_command("--promote"),
            "continuation_verify": _expected_command("--verify"),
            "runner_gate_replay": _runner_gate_command(),
            "downstream_continuation": _downstream_continuation_command(),
        }
        expected_authorization_signature_contract = (
            _exit_attestation_signature_contract(
                path=PREPARE_EXIT_ATTESTATION,
                schema_name=PREPARE_EXIT_ATTESTATION_SCHEMA,
                signature=AUTHORIZATION_SIGNATURE,
                public_key=auth_public_record,
                private_key=AUTHORIZATION_PRIVATE_KEY,
                legacy_payload_signature=LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
            )
        )
        expected_authorized_completion_contract = (
            _exit_attestation_signature_contract(
                path=PROMOTE_EXIT_ATTESTATION,
                schema_name=PROMOTE_EXIT_ATTESTATION_SCHEMA,
                signature=COMPLETION_SIGNATURE,
                public_key=completion_public_record,
                private_key=COMPLETION_PRIVATE_KEY,
                legacy_payload_signature=LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
            )
        )
        authorization_full_record = _stat_artifact(auth_record, auth_stat)
        prepare_exit_full_record = _stat_artifact(
            prepare_exit_record, prepare_exit_stat
        )
        authorization_signature_full_record = _stat_artifact(
            auth_sig_record, auth_sig_stat
        )
        promotion_source_full_record = _stat_artifact(
            promotion_record, promotion_stat
        )
        python_runtime_full_record = _stat_artifact(python_record, python_stat)
        prepare_exit_payload = _parse_json(
            prepare_exit_data, "prepare producer exit attestation"
        )
        prepare_exit_validation = _validate_exit_attestation(
            prepare_exit_payload,
            label="prepare producer exit attestation",
            expected_schema=PREPARE_EXIT_ATTESTATION_SCHEMA,
            expected_phase="prepare_authorization",
            expected_command=_promotion_command("--prepare-authorization"),
            expected_produced_artifacts={"authorization": authorization_full_record},
            expected_producer_source=promotion_source_full_record,
            expected_python_runtime=python_runtime_full_record,
            expected_signature_contract=expected_authorization_signature_contract,
            expected_checks=PREPARE_EXIT_ATTESTATION_CHECKS,
            expected_pass_semantics=PREPARE_EXIT_ATTESTATION_PASS_SEMANTICS,
        )
        if prepare_exit_validation["created_at_utc"] < authorization_created_at:
            raise VerificationError(
                "Prepare exit-attestation timestamp predates authorization"
            )
        expected_signed_prepare_chain = {
            "artifact": authorization_full_record,
            "producer_exit_attestation": prepare_exit_full_record,
            "signature": authorization_signature_full_record,
            "public_key": auth_public_record,
            "retired_private_key": auth_private_record,
            "producer_session": prepare_exit_validation["producer_session"],
            "signature_verified": True,
        }
        expected_continuation_gate = {
            "verifier": self_record,
            "verification_receipt": str(VERIFICATION_RECEIPT),
            "runner_gate_replay": replayer_record,
            "runner_gate_receipt": str(RUNNER_GATE_RECEIPT),
            "downstream_continuation": continuation_record,
            "downstream_continuation_command": _downstream_continuation_command(),
            "supervised_child_command": _supervised_child_command(),
            "supervisor_capability_protocol": SUPERVISOR_CAPABILITY_PROTOCOL,
            "threat_model": CONTINUATION_THREAT_MODEL,
            "terminal_signature_contract": {
                "algorithm": "Ed25519",
                "public_key": continuation_public_record,
                "private_key": continuation_private_record,
                "execution_receipt": str(CONTINUATION_TERMINAL_RECEIPT),
                "signed_artifact": str(CONTINUATION_EXIT_ATTESTATION),
                "signature": str(CONTINUATION_TERMINAL_SIGNATURE),
                "required_private_key_pre_signature_mode": "0o400",
                "required_private_key_post_signature_mode": "0o0",
                "supervisor_command": _downstream_continuation_command(),
                "supervised_child_command": _supervised_child_command(),
                "signed_terminal_verification_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--verify-signed-terminal",
                ],
            },
            "policy": CONTINUATION_POLICY,
        }
        expected_authorization_disclosure = {
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
        expected_prepare_exit_announcement = _exit_attestation_announcement(
            PREPARE_EXIT_ATTESTATION,
            PREPARE_EXIT_ATTESTATION_SCHEMA,
            auth_public_record,
        )
        if (
            authorization_payload.get("schema_name") != AUTHORIZATION_SCHEMA
            or authorization_payload.get("schema_version") != SCHEMA_VERSION
            or authorization_payload.get("authorization_id") != AUTHORIZATION_ID
            or authorization_payload.get("task_type") != "B_comprehensive_validation"
            or authorization_payload.get("status") != AUTHORIZATION_STATUS
            or authorization_payload.get("scope")
            != "canonical_path_publication_only_no_scientific_payload_change"
            or not _strict_json_equal(
                authorization_payload.get("producer_exit_attestation"),
                expected_prepare_exit_announcement,
            )
            or not _strict_json_equal(authorization_payload.get("commands"), expected_commands)
            or not _strict_json_equal(
                authorization_payload.get("reviewed_sources"), reviewed_sources
            )
            or not _strict_json_equal(
                authorization_payload.get("producer_session"),
                prepare_exit_validation["producer_session"],
            )
            or not _strict_json_equal(
                authorization_payload.get("review_artifacts"), review_records
            )
            or not _strict_json_equal(
                authorization_payload.get("trusted_signing_keys"), trusted_keys
            )
            or not _strict_json_equal(
                authorization_payload.get("source_receipt"), source_record
            )
            or authorization_payload.get("canonical_target_path")
            != str(CANONICAL_RECEIPT)
            or not _strict_json_equal(
                authorization_payload.get("authorization_signature_contract"),
                expected_authorization_signature_contract,
            )
            or not _strict_json_equal(
                authorization_payload.get("completion_signature_contract"),
                expected_authorized_completion_contract,
            )
            or not _strict_json_equal(
                authorization_payload.get("continuation_gate"), expected_continuation_gate
            )
            or not _strict_json_equal(
                authorization_payload.get("historical_incident_disclosure"),
                expected_authorization_disclosure,
            )
            or not _strict_json_equal(
                authorization_payload.get("downstream_output_absence"),
                [str(path) for path in DOWNSTREAM_OUTPUTS],
            )
            or not _strict_json_equal(
                authorization_payload.get("checks"), AUTHORIZATION_CHECKS
            )
            or authorization_payload.get("pass") is not True
            or authorization_payload.get("pass_semantics")
            != (
                "authorization_payload_only; parent exit-attestation publication and detached "
                "authorization-key signature remain mandatory; no canonical publication has "
                "occurred and this payload does not permit downstream execution"
            )
        ):
            raise VerificationError("Signed authorization exact schema/cross-link drift")

        completion_payload = _parse_json(completion_data, "promotion completion receipt")
        _require_exact_keys(completion_payload, COMPLETION_KEYS, "promotion completion")
        completion_created_at = _require_utc_timestamp(
            completion_payload.get("created_at_utc"),
            "promotion completion created_at_utc",
        )
        if completion_created_at < authorization_created_at:
            raise VerificationError(
                "Promotion completion timestamp predates authorization"
        )
        expected_authorization_link = expected_signed_prepare_chain
        expected_signature_contract = expected_authorized_completion_contract
        completion_full_record = _stat_artifact(completion_record, completion_stat)
        promote_exit_full_record = _stat_artifact(
            promote_exit_record, promote_exit_stat
        )
        promote_exit_payload = _parse_json(
            promote_exit_data, "promote producer exit attestation"
        )
        promote_exit_validation = _validate_exit_attestation(
            promote_exit_payload,
            label="promote producer exit attestation",
            expected_schema=PROMOTE_EXIT_ATTESTATION_SCHEMA,
            expected_phase="promote",
            expected_command=_promotion_command("--promote"),
            expected_produced_artifacts={
                "completion_receipt": completion_full_record,
                "canonical_receipt": canonical_record,
                "prepare_exit_attestation": prepare_exit_full_record,
                "prepare_exit_attestation_signature": (
                    authorization_signature_full_record
                ),
            },
            expected_producer_source=promotion_source_full_record,
            expected_python_runtime=python_runtime_full_record,
            expected_signature_contract=expected_signature_contract,
            expected_checks=PROMOTE_EXIT_ATTESTATION_CHECKS,
            expected_pass_semantics=PROMOTE_EXIT_ATTESTATION_PASS_SEMANTICS,
        )
        if (
            promote_exit_validation["created_at_utc"] < completion_created_at
            or promote_exit_validation["created_at_utc"]
            < prepare_exit_validation["created_at_utc"]
            or promote_exit_validation["producer_session"]["nonce_sha256"]
            == prepare_exit_validation["producer_session"]["nonce_sha256"]
        ):
            raise VerificationError(
                "Promote exit-attestation chronology/session separation drift"
            )
        _require_exact_keys(
            completion_payload.get("authorization"),
            AUTHORIZATION_LINK_KEYS,
            "completion authorization link",
        )
        _require_exact_keys(
            completion_payload.get("signature_contract"),
            SIGNATURE_CONTRACT_KEYS,
            "completion signature contract",
        )
        _require_stat_artifact(
            completion_payload.get("source_receipt"), source_record, "completion source receipt"
        )
        _require_basic_artifact(
            completion_payload.get("canonical_receipt"),
            canonical_record,
            "completion canonical receipt",
        )
        _validate_execution_binding(
            completion_payload.get("execution_binding"),
            promotion_record,
            promotion_data,
        )
        expected_completion_disclosure = {
            "original_completion_runner_end_to_end_pass_claimed": False,
            "original_post_run_verifier_success_claimed": False,
            "canonical_receipt_created_by": (
                "signed_v3_exit_attested_byte_identical_promotion"
            ),
            "failed_attempt_archive": authorization_payload["evidence"][
                "failed_attempt_archive"
            ],
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
        expected_governance = {
            "prior_review": authorization_payload["evidence"]["prior_policy_review"],
            "narrow_supersession": (
                "Only the canonical-emitter requirement is superseded; all content checks, "
                "source authorities, downstream gate replay, and release signatures remain."
            ),
            "continuation_verifier": self_record,
            "runner_gate_replay": replayer_record,
            "downstream_continuation": continuation_record,
            "downstream_continuation_command": _downstream_continuation_command(),
            "supervised_child_command": _supervised_child_command(),
            "continuation_verification_receipt": str(VERIFICATION_RECEIPT),
            "runner_gate_receipt": str(RUNNER_GATE_RECEIPT),
            "downstream_output_slots_absent_at_promotion": [
                str(path) for path in DOWNSTREAM_OUTPUTS
            ],
            "canonical_publication_order": (
                "provisional_completion_receipt_then_canonical_no_replace_commit"
            ),
            "canonical_alone_grants_downstream_authority": False,
            "write_ahead_receipt_marks_transaction_incomplete_until": (
                "detached_completion_signature_plus_live_canonical_verification"
            ),
        }
        expected_builder_gate = {
            "source_path_status": "VERIFIED_BOUNDED_RETROSPECTIVE_SOURCE_IDENTITY",
            "canonical_path_status": "VERIFIED_BOUNDED_RETROSPECTIVE_SOURCE_IDENTITY",
            "release_gate_pass": True,
            "publishable_task_b_release": True,
            "builder_loaded_from_v7_bound_bytes": True,
            "bytecode_disabled": True,
        }
        expected_promote_exit_announcement = _exit_attestation_announcement(
            PROMOTE_EXIT_ATTESTATION,
            PROMOTE_EXIT_ATTESTATION_SCHEMA,
            completion_public_record,
        )
        if (
            completion_payload.get("schema_name") != COMPLETION_SCHEMA
            or completion_payload.get("schema_version") != SCHEMA_VERSION
            or completion_payload.get("authorization_id") != AUTHORIZATION_ID
            or completion_payload.get("task_type") != "B_comprehensive_validation"
            or completion_payload.get("status")
            != "PROVISIONAL_WRITE_AHEAD_CANONICAL_PUBLICATION_INTENT"
            or completion_payload.get("scope")
            != "canonical_path_publication_only_no_scientific_payload_change"
            or not _strict_json_equal(
                completion_payload.get("producer_exit_attestation"),
                expected_promote_exit_announcement,
            )
            or not _strict_json_equal(
                completion_payload.get("command"), _promotion_command("--promote")
            )
            or not _strict_json_equal(
                completion_payload.get("authorization"), expected_authorization_link
            )
            or not _strict_json_equal(
                completion_payload.get("code"), reviewed_sources
            )
            or not _strict_json_equal(
                completion_payload.get("execution_binding"),
                authorization_payload.get("execution_binding"),
            )
            or not _strict_json_equal(
                completion_payload.get("producer_session"),
                promote_exit_validation["producer_session"],
            )
            or not _strict_json_equal(
                completion_payload.get("evidence"), authorization_payload.get("evidence")
            )
            or not _strict_json_equal(
                completion_payload.get("review_artifacts"), review_records
            )
            or not _strict_json_equal(
                completion_payload.get("builder_gate_replay"), expected_builder_gate
            )
            or not _strict_json_equal(
                completion_payload.get("historical_incident_disclosure"),
                expected_completion_disclosure,
            )
            or not _strict_json_equal(
                completion_payload.get("governance"), expected_governance
            )
            or not _strict_json_equal(
                completion_payload.get("signature_contract"),
                expected_signature_contract,
            )
            or not _strict_json_equal(
                completion_payload.get("checks"), COMPLETION_CHECKS
            )
            or completion_payload.get("pass") is not True
            or completion_payload.get("pass_semantics")
            != (
                "completion_payload_only; parent promote exit-attestation publication and "
                "detached completion-key signature remain mandatory, followed by continuation "
                "verification, runner-gate evidence, and the signed downstream terminal receipt"
            )
            or not _strict_json_equal(
                _basic_projection(completion_payload["source_receipt"]),
                _basic_projection(authorization_payload["source_receipt"]),
            )
            or completion_payload["canonical_receipt"]["path"]
            != authorization_payload["canonical_target_path"]
        ):
            raise VerificationError("Signed completion exact schema/cross-link drift")

        reader.require_paths_still_bound()
        _require_bootstrap_paths_still_bound()
        _require_direct_paths_still_bound()
        retained_count = reader.opened_descriptor_count + len(_BOOTSTRAP_FD_LEASES)
        payload = _verification_payload(
            mode=mode,
            authorization_record=auth_record,
            prepare_exit_attestation_record=prepare_exit_full_record,
            authorization_signature_record=auth_sig_record,
            completion_record=completion_record,
            promote_exit_attestation_record=promote_exit_full_record,
            completion_signature_record=completion_sig_record,
            source_record=source_record,
            canonical_record=canonical_record,
            trusted_keys=trusted_keys,
            authorization_private_record=auth_private_record,
            completion_private_record=completion_private_record,
            incident_disclosure=incident_disclosure,
            release_authority_evidence=authority_record,
            retained_descriptor_count=retained_count,
        )

        if os.path.lexists(PROMOTION_INCIDENT):
            raise VerificationError("Promotion incident appeared before verification commit")
        _require_absent(
            LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
            "legacy authorization-payload signature before verification commit",
        )
        _require_absent(
            LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
            "legacy completion-payload signature before verification commit",
        )
        if any(os.path.lexists(path) for path in DOWNSTREAM_OUTPUTS):
            raise VerificationError("Downstream output appeared before verification commit")
        reader.require_paths_still_bound()
        _require_bootstrap_paths_still_bound()
        _require_direct_paths_still_bound()
        reader.retain_until_process_exit()

        verification_record: dict[str, Any] | None = None
        if mode == "--verify-and-record":
            def verification_precommit() -> None:
                if os.path.lexists(PROMOTION_INCIDENT):
                    raise VerificationError(
                        "Promotion incident appeared during verification staging"
                    )
                _require_absent(
                    LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
                    "legacy authorization-payload signature during verification staging",
                )
                _require_absent(
                    LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
                    "legacy completion-payload signature during verification staging",
                )
                if any(os.path.lexists(path) for path in DOWNSTREAM_OUTPUTS):
                    raise VerificationError(
                        "Downstream output appeared during verification staging"
                    )
                reader.require_paths_still_bound()
                _require_bootstrap_paths_still_bound()
                _require_direct_paths_still_bound()

            _write_exclusive_readonly(
                VERIFICATION_RECEIPT,
                _encode_json(payload),
                precommit=verification_precommit,
                terminal_commit=lambda: os._exit(0),
            )
        return payload, verification_record


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) == 2 else ""
    try:
        _require_clean_direct_runtime(mode)
        payload, _ = _verify(mode)
        print(json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True))
        return 0
    except Exception as error:
        failure = {
            "schema_name": f"{VERIFICATION_SCHEMA}.error",
            "schema_version": SCHEMA_VERSION,
            "error": {
                "type": type(error).__name__,
                "message": str(error),
            },
            "pass": False,
        }
        print(json.dumps(failure, ensure_ascii=True, sort_keys=True), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
