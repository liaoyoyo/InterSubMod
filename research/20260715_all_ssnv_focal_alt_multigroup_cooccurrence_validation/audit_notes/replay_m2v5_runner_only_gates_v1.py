#!/usr/bin/env python3
"""Replay only the immutable completion runner's pre-downstream gates.

This is one-shot non-authoritative replay evidence. It verifies the signed tumor-REF receipt
promotion first, executes exactly physical lines 1-358 of the v7-authorized
completion runner through a bound Bash descriptor, and publishes immutable
evidence only when every check passes. It never executes line 359 or later and
does not itself authorize downstream execution.
"""

from __future__ import annotations

import base64
import ctypes
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import types
from typing import Any, Callable, Mapping, Sequence


AT_FDCWD = -100
AT_SYMLINK_FOLLOW = 0x400


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
REPLAYER = AUDIT_ROOT / "replay_m2v5_runner_only_gates_v1.py"

PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")
PYTHON_TARGET = Path(
    "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python3.11"
)
BASH = Path("/usr/bin/bash")
OPENSSL = Path("/usr/bin/openssl")
COMPLETION_RUNNER = SCRIPT_ROOT / "run_m2v5_recovered_completion_chain.sh"
AUTHORITY_VALIDATOR = SCRIPT_ROOT / "release_source_authority.py"
CONTINUATION_VERIFIER = AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_v2.py"
PROMOTION_TOOL = AUDIT_ROOT / "promote_tumor_ref_source_receipt_v2.py"
PROMOTION_CACHE_ROOTS = {
    "--prepare-authorization": (
        WORKSPACE_ROOT / ".python_cache_tumor_ref_promotion_v2_prepare"
    ),
    "--promote": WORKSPACE_ROOT / ".python_cache_tumor_ref_promotion_v2_promote",
}
DOWNSTREAM_CONTINUATION = (
    AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_v1.py"
)
QA_PYTHON_TARGET = Path("/bip7_disk/liaoyoyo2001/miniconda3/bin/python3.9")
PORTABLE_PLUGIN_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/"
    "data-analytics/0.2.8-13ceeea1f599"
)
PORTABLE_PLUGIN_SCRIPT_ROOT = (
    PORTABLE_PLUGIN_ROOT / "skills" / "build-report" / "scripts"
)
PORTABLE_ASSET_ROOT = PORTABLE_PLUGIN_ROOT / "assets"
PORTABLE_READER_ASSET_PARTS = tuple(
    PORTABLE_ASSET_ROOT / f"portable-artifact-reader.html.gz.b64.part{index:03d}"
    for index in range(1, 4)
)
REVIEWED_RUNTIME_SOURCE_CONTRACTS = {
    "primary_python_runtime": (
        PYTHON_TARGET,
        "777797a57eb75b28f530191628d26b14afada9ced2cb51c0ecae1eb62796062e",
        "0o775",
    ),
    "qa_python_runtime": (
        QA_PYTHON_TARGET,
        "cbd71e11d993aa3b7cc3ad553b97733f8054c1465ed7764dbbd386b9897bd5bc",
        "0o755",
    ),
    "portable_builder_module": (
        PORTABLE_PLUGIN_SCRIPT_ROOT / "build_portable_artifact.mjs",
        "0b86883810cf23c81f563a48a7708283a6a4cadf11228ad261efb64426ef728e",
        "0o644",
    ),
    "portable_chart_extractor_module": (
        PORTABLE_PLUGIN_SCRIPT_ROOT / "extract_portable_chart_svgs.mjs",
        "77e78593d8971563b6394334128dd3ee243f2182d3d087591d3a5499ece9283b",
        "0o644",
    ),
    "portable_verifier_module": (
        PORTABLE_PLUGIN_SCRIPT_ROOT / "verify_portable_artifact.mjs",
        "b495c4cc34113fb2918118eac302f4e7e2152c1b1b9b63e3646cea97ecbf9b3f",
        "0o644",
    ),
    "portable_browser_helpers_module": (
        PORTABLE_PLUGIN_SCRIPT_ROOT / "portable_browser_helpers.mjs",
        "84aa4f8a2a11376ebee6942d2d7e10a083d16c65dfdb73114d7b34e51c27f69d",
        "0o644",
    ),
    "portable_browser_cli_module": (
        PORTABLE_PLUGIN_SCRIPT_ROOT / "portable_browser_cli.mjs",
        "aac3b12fc12c7ad2e044533f881791dd3b23bd0ee31ddf6682dd7f6de99e6596",
        "0o644",
    ),
    "portable_server_bundle": (
        PORTABLE_PLUGIN_ROOT / "mcp" / "server.cjs",
        "eff59c6085d2ab6b6153c80a03749e764e160f8c6711da8433f7bd6762e1db66",
        "0o644",
    ),
    **{
        f"portable_reader_asset_part{index:03d}": (
            path,
            sha256,
            "0o644",
        )
        for index, (path, sha256) in enumerate(
            zip(
                PORTABLE_READER_ASSET_PARTS,
                (
                    "154f1d561bab28174f88d71ae710599709e5a5eda64dee08b025a699449dbbfc",
                    "9459ed4b76bd825daba9637564dd3122ca617a88dd16531d3e9bca5c7ced080e",
                    "cdfe2e6787faa61d37f043df6996d5f988c07360bb2b3d8af2aa3c18e8db3ac0",
                ),
                strict=True,
            ),
            start=1,
        )
    },
}

V7_AUTHORITY = (
    REPO_ROOT
    / "docs/provenance/source_authorities/"
    "20260722_all_ssnv_focal_alt_release_source_authority.v7.json"
)
V7_APPROVAL = Path(str(V7_AUTHORITY).removesuffix(".json") + ".approval.json")
V7_SIGNATURE = Path(str(V7_APPROVAL) + ".ed25519.sig")
V7_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
    "20260722_all_ssnv_v10_strict_command_parity_bootstrap/ed25519_public.pem"
)
V7_PRIVATE_KEY = V7_PUBLIC_KEY.with_name(
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)

AUTHORIZATION = RESULT_ROOT / "tumor_ref_source_receipt_promotion_authorization.v3.json"
PREPARE_EXIT_ATTESTATION = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_prepare_exit_attestation.v1.json"
)
AUTHORIZATION_SIGNATURE = Path(str(PREPARE_EXIT_ATTESTATION) + ".ed25519.sig")
LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE = Path(str(AUTHORIZATION) + ".ed25519.sig")
AUTHORIZATION_KEY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_promotion_authorization/"
    "20260722_tumor_ref_receipt_promotion_v4"
)
AUTHORIZATION_PUBLIC_KEY = AUTHORIZATION_KEY_ROOT / "ed25519_public.pem"
AUTHORIZATION_PRIVATE_KEY = AUTHORIZATION_KEY_ROOT / (
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)

PROMOTION_RECEIPT = RESULT_ROOT / "tumor_ref_source_receipt_promotion.v3.json"
PROMOTE_EXIT_ATTESTATION = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_promote_exit_attestation.v1.json"
)
PROMOTION_SIGNATURE = Path(str(PROMOTE_EXIT_ATTESTATION) + ".ed25519.sig")
LEGACY_COMPLETION_PAYLOAD_SIGNATURE = Path(str(PROMOTION_RECEIPT) + ".ed25519.sig")
PROMOTION_INCIDENT = RESULT_ROOT / "tumor_ref_source_receipt_promotion_incident.v3.json"
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
EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256 = (
    "a98370415594ebc9e8eb166bb70cab3550067b8a83db7d73ad23fac0891cb8b5"
)
EXPECTED_AUTHORIZATION_PUBLIC_KEY_SHA256 = (
    "e638b29a1d9c207dfa9849ff20a3867dacf702f7c1a899f1968ea1401a839677"
)
EXPECTED_COMPLETION_PUBLIC_KEY_SHA256 = (
    "6c55230f368db6b1fb0eedbc6c2362c0df7ee41a2699764fadb1f424e7036032"
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

HISTORICAL_SOURCE_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json"
)
CANONICAL_SOURCE_RECEIPT = (
    WORKSPACE_ROOT
    / "tumor_ref_recovery_source_identity_v1"
    / "post_run_source_identity.receipt.json"
)
PROMOTION_VERIFICATION_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.v3.json"
)

OUTPUT_LOG = RESULT_ROOT / "m2v5_runner_only_gate_replay.v1.log"
OUTPUT_RECEIPT = RESULT_ROOT / "m2v5_runner_only_gate_replay.v1.json"

EXPECTED_V7_AUTHORITY_ID = (
    "20260722_all_ssnv_focal_alt_task_b_release_v7_strict_command_parity"
)
EXPECTED_AUTHORIZATION_SCHEMA = (
    "intersubmod.tumor_ref_source_receipt_promotion_authorization"
)
EXPECTED_PROMOTION_SCHEMA = "intersubmod.tumor_ref_source_receipt_promotion"
EXPECTED_VERIFICATION_SCHEMA = (
    "intersubmod.tumor_ref_source_receipt_promotion_verification"
)
PROMOTION_SCHEMA_VERSION = "3.0.0"
EXIT_ATTESTATION_SCHEMA_VERSION = "1.0.0"
AUTHORIZATION_ID = "20260722_tumor_ref_receipt_promotion_v3"
PREPARE_EXIT_ATTESTATION_SCHEMA = (
    "intersubmod.tumor_ref_source_receipt_promotion.prepare_exit_attestation"
)
PROMOTE_EXIT_ATTESTATION_SCHEMA = (
    "intersubmod.tumor_ref_source_receipt_promotion.promote_exit_attestation"
)
EXPECTED_VERIFICATION_PASS_SEMANTICS = (
    "signed byte-identical promotion chain verified; this does not claim the original "
    "completion runner or original verifier completed successfully"
)
EXPECTED_RUNNER_LINE_COUNT_MINIMUM = 633
RUNNER_INCLUDED_LINE_START = 1
RUNNER_INCLUDED_LINE_END = 358
RUNNER_FIRST_EXCLUDED_LINE = 359
EXPECTED_RUNNER_PREFIX_SHA256 = (
    "a5409c1a25f11e7a85b1bc9305689e4a52cdd4e81281cef7fb6658780bb97814"
)
VERIFICATION_RECEIPT_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "task_type",
        "mode",
        "authorization",
        "completion",
        "producer_exit_attestations",
        "source_receipt",
        "canonical_receipt",
        "trusted_signing_keys",
        "private_key_retirement",
        "historical_incident_disclosure",
        "release_source_authority_validator",
        "checks",
        "retained_critical_input_descriptor_count",
        "governance",
        "pass",
        "pass_semantics",
    }
)
VERIFICATION_CHECK_KEYS = frozenset(
    {
        "clean_direct_runtime",
        "executing_verifier_definitions_match_bound_reviewed_source",
        "authority_module_executed_from_bound_source_bytes",
        "v7_release_source_authority_verified_with_process_lifetime_fd_lease",
        "authorization_schema_and_exact_keys_verified",
        "authorization_signature_verified_from_bound_fds",
        "prepare_exit_attestation_exact_schema_verified",
        "prepare_child_waitpid_normal_exit_zero_verified",
        "authorization_signature_verified_over_prepare_exit_attestation",
        "legacy_authorization_payload_signature_absent",
        "authorization_private_key_retired_mode_000",
        "completion_schema_and_exact_keys_verified",
        "completion_signature_verified_from_bound_fds",
        "promote_exit_attestation_exact_schema_verified",
        "promote_child_waitpid_normal_exit_zero_verified",
        "completion_signature_verified_over_promote_exit_attestation",
        "legacy_completion_payload_signature_absent",
        "promote_attestation_binds_signed_prepare_chain",
        "completion_private_key_retired_mode_000",
        "three_distinct_review_artifacts_verified",
        "four_core_plus_eleven_runtime_sources_verified",
        "portable_reader_asset_membership_exact",
        "source_and_canonical_bytes_size_sha256_equal",
        "source_and_canonical_independent_inodes",
        "canonical_regular_mode_0444_link_count_one",
        "historical_incident_disclosure_hash_bound",
        "exact_artifact_cross_links_verified",
        "critical_artifact_paths_rebound_before_success",
        "critical_descriptor_leases_retained_until_process_exit",
    }
)
VERIFICATION_GOVERNANCE_KEYS = frozenset(
    {
        "original_completion_runner_end_to_end_pass_claimed",
        "original_verifier_reexecuted_successfully_claimed",
        "runner_only_gate_replay_evidence_required_before_downstream",
        "signed_downstream_continuation_revalidation_required",
        "verification_scope",
    }
)

BASIC_ARTIFACT_KEYS = frozenset({"path", "size_bytes", "sha256", "mode"})
FULL_STAT_ARTIFACT_KEYS = frozenset(
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
PRODUCER_EXIT_RECEIPT_ENTRY_KEYS = frozenset(
    {"artifact", "signature", "public_key", "signature_verified"}
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
AUTHORIZATION_PASS_SEMANTICS = (
    "authorization_payload_only; parent exit-attestation publication and detached "
    "authorization-key signature remain mandatory; no canonical publication has "
    "occurred and this payload does not permit downstream execution"
)
COMPLETION_PASS_SEMANTICS = (
    "completion_payload_only; parent promote exit-attestation publication and "
    "detached completion-key signature remain mandatory, followed by continuation "
    "verification, runner-gate evidence, and the signed downstream terminal receipt"
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

# These runtime pins are transitively authorized by the signed exact-source
# identity of this gate. The completion runner itself is separately pinned by v7.
EXPECTED_RUNTIME_ARTIFACTS = {
    "python": {
        "path": str(PYTHON),
        "resolved_path": str(PYTHON_TARGET),
        "symlink_target": "python3.11",
        "size_bytes": 25_409_784,
        "sha256": "777797a57eb75b28f530191628d26b14afada9ced2cb51c0ecae1eb62796062e",
        "mode": "0o775",
    },
    "bash": {
        "path": str(BASH),
        "size_bytes": 1_396_520,
        "sha256": "59474588a312b6b6e73e5a42a59bf71e62b55416b6c9d5e4a6e1c630c2a9ecd4",
        "mode": "0o755",
    },
}

STRICT = WORKSPACE_ROOT / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5"
MATCHED_RUN = WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5"
MATCHED_ANALYSIS = (
    WORKSPACE_ROOT / "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5"
)
CN_CCF = WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5"
FINAL_DATASET = WORKSPACE_ROOT / "all_ssnv_final_report_dataset_v5_m2v5_source_attested"
FINAL_REPORT = WORKSPACE_ROOT / "all_ssnv_final_report_v5_m2v5_source_attested"
ORIGINAL_RUNNER_LOG = WORKSPACE_ROOT / (
    "m2v5_recovered_completion_chain_v2_source_authority_v5.log"
)
ORIGINAL_RUNNER_CACHE = WORKSPACE_ROOT / (
    ".python_cache_m2v5_completion_v2_bound_bootstrap"
)

DOWNSTREAM_OUTPUT_SLOTS = (
    STRICT,
    MATCHED_RUN,
    MATCHED_ANALYSIS,
    CN_CCF,
    RESULT_ROOT / "stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json",
    RESULT_ROOT / "frozen_input_immutability.post_m2v5_downstream_v3_source_authority_v5.json",
    FINAL_DATASET,
    FINAL_REPORT,
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json.ed25519.sig",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json.ed25519.sig",
    ORIGINAL_RUNNER_LOG,
    ORIGINAL_RUNNER_CACHE,
    RESULT_ROOT / "m2v5_downstream_continuation.v1.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.v1.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.v1.json.ed25519.sig",
)

_PROCESS_LIFETIME_READERS: list["BoundArtifactReader"] = []
_PROCESS_LIFETIME_MODULES: list[types.ModuleType] = []


class GateError(RuntimeError):
    """Raised when any continuation-gate invariant is not exact."""


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def parse_utc_timestamp(value: Any, label: str) -> datetime:
    if not isinstance(value, str):
        raise GateError(f"{label} is not a UTC timestamp string")
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise GateError(f"{label} is not ISO-8601") from error
    if (
        parsed.tzinfo is None
        or parsed.utcoffset() != timezone.utc.utcoffset(parsed)
        or parsed.microsecond != 0
        or value != parsed.isoformat(timespec="seconds")
    ):
        raise GateError(f"{label} is not canonical UTC")
    return parsed


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


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
        elif isinstance(value, type) and value.__module__ == namespace.get("__name__"):
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
    return inventory


def critical_globals_equal(
    live: Mapping[str, Any], shadow: Mapping[str, Any]
) -> bool:
    names = (
        "REPO_ROOT",
        "TOPIC_ROOT",
        "AUDIT_ROOT",
        "RESULT_ROOT",
        "WORKSPACE_ROOT",
        "REPLAYER",
        "PYTHON",
        "PYTHON_TARGET",
        "QA_PYTHON_TARGET",
        "PORTABLE_PLUGIN_ROOT",
        "PORTABLE_PLUGIN_SCRIPT_ROOT",
        "PORTABLE_ASSET_ROOT",
        "PORTABLE_READER_ASSET_PARTS",
        "REVIEWED_RUNTIME_SOURCE_CONTRACTS",
        "BASH",
        "OPENSSL",
        "COMPLETION_RUNNER",
        "AUTHORITY_VALIDATOR",
        "CONTINUATION_VERIFIER",
        "PROMOTION_TOOL",
        "PROMOTION_CACHE_ROOTS",
        "DOWNSTREAM_CONTINUATION",
        "AUTHORIZATION",
        "AUTHORIZATION_SIGNATURE",
        "PREPARE_EXIT_ATTESTATION",
        "LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE",
        "PROMOTION_RECEIPT",
        "PROMOTION_SIGNATURE",
        "PROMOTE_EXIT_ATTESTATION",
        "LEGACY_COMPLETION_PAYLOAD_SIGNATURE",
        "PROMOTION_INCIDENT",
        "PROMOTION_VERIFICATION_RECEIPT",
        "CONTINUATION_PUBLIC_KEY",
        "CONTINUATION_PRIVATE_KEY",
        "CONTINUATION_TERMINAL_RECEIPT",
        "CONTINUATION_EXIT_ATTESTATION",
        "CONTINUATION_TERMINAL_SIGNATURE",
        "SUPERVISOR_CAPABILITY_PROTOCOL",
        "CONTINUATION_THREAT_MODEL",
        "CONTINUATION_POLICY",
        "OUTPUT_LOG",
        "OUTPUT_RECEIPT",
        "DOWNSTREAM_OUTPUT_SLOTS",
        "EXPECTED_ENVIRONMENT",
        "EXPECTED_RUNNER_PREFIX_SHA256",
        "EXPECTED_VERIFICATION_PASS_SEMANTICS",
        "EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256",
        "EXPECTED_AUTHORIZATION_PUBLIC_KEY_SHA256",
        "EXPECTED_COMPLETION_PUBLIC_KEY_SHA256",
        "PROMOTION_SCHEMA_VERSION",
        "EXIT_ATTESTATION_SCHEMA_VERSION",
        "AUTHORIZATION_ID",
        "PREPARE_EXIT_ATTESTATION_SCHEMA",
        "PROMOTE_EXIT_ATTESTATION_SCHEMA",
        "AUTHORIZATION_KEYS",
        "COMPLETION_KEYS",
        "AUTHORIZATION_LINK_KEYS",
        "FULL_STAT_ARTIFACT_KEYS",
        "EXIT_ATTESTATION_KEYS",
        "PRODUCER_SESSION_KEYS",
        "CHILD_WAIT_KEYS",
        "PREPARE_EXIT_ATTESTATION_CHECKS",
        "PROMOTE_EXIT_ATTESTATION_CHECKS",
        "PREPARE_EXIT_ATTESTATION_PASS_SEMANTICS",
        "PROMOTE_EXIT_ATTESTATION_PASS_SEMANTICS",
        "AUTHORIZATION_CHECKS",
        "COMPLETION_CHECKS",
        "AUTHORIZATION_PASS_SEMANTICS",
        "COMPLETION_PASS_SEMANTICS",
        "VERIFICATION_RECEIPT_KEYS",
        "VERIFICATION_CHECK_KEYS",
    )
    return all(
        name in live
        and name in shadow
        and type(live[name]) is type(shadow[name])
        and live[name] == shadow[name]
        for name in names
    )


def executing_module_code() -> types.CodeType:
    frame = sys._getframe()
    while frame is not None:
        if frame.f_code.co_name == "<module>" and frame.f_globals is globals():
            return frame.f_code
        frame = frame.f_back
    raise GateError("Executing runner-gate module frame is unavailable")


def require_executing_replayer_matches_bound_source(data: bytes) -> None:
    module_code = executing_module_code()
    filename = module_code.co_filename
    shadow: dict[str, Any] = {
        "__name__": "runner_gate_replay_bound_shadow",
        "__file__": filename,
        "__package__": "",
    }
    reviewed_module_code = compile(data, filename, "exec")
    exec(reviewed_module_code, shadow)
    live_inventory = code_inventory(globals(), filename)
    shadow_inventory = code_inventory(shadow, filename)
    if (
        not live_inventory
        or stable_code_digest(module_code) != stable_code_digest(reviewed_module_code)
        or live_inventory != shadow_inventory
        or not critical_globals_equal(globals(), shadow)
    ):
        raise GateError("Executing runner-gate module differs from bound source bytes")


def parse_json(data: bytes, label: str) -> dict[str, Any]:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        output: dict[str, Any] = {}
        for key, value in pairs:
            if key in output:
                raise GateError(f"Duplicate JSON key in {label}: {key}")
            output[key] = value
        return output

    try:
        payload = json.loads(
            data,
            object_pairs_hook=no_duplicates,
            parse_constant=lambda value: (_ for _ in ()).throw(
                GateError(f"Non-finite JSON value in {label}: {value}")
            ),
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise GateError(f"Invalid JSON in {label}") from error
    if not isinstance(payload, dict):
        raise GateError(f"{label} is not a JSON object")
    return payload


def encode_json(payload: Mapping[str, Any]) -> bytes:
    return (json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode(
        "utf-8"
    )


def artifact_record(path: Path, data: bytes, opened: os.stat_result) -> dict[str, Any]:
    return {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": sha256_bytes(data),
        "mode": oct(stat.S_IMODE(opened.st_mode)),
    }


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


class BoundArtifactReader:
    """Bind regular trust artifacts to stable descriptors until process exit."""

    def __init__(self) -> None:
        self._opened: dict[Path, tuple[int, bytes | None, os.stat_result, dict[str, Any]]] = {}
        self._aliases: dict[Path, tuple[int, os.stat_result, str, Path]] = {}
        self._retained = False

    @property
    def descriptor_count(self) -> int:
        return len(self._opened) + len(self._aliases)

    def _canonical_regular_path(self, path: Path) -> tuple[Path, os.stat_result]:
        if not path.is_absolute():
            raise GateError(f"Trust path is not absolute: {path}")
        try:
            leaf = os.lstat(path)
        except FileNotFoundError as error:
            raise GateError(f"Trust artifact is missing: {path}") from error
        if stat.S_ISLNK(leaf.st_mode):
            raise GateError(f"Trust artifact is a symlink: {path}")
        resolved = path.resolve(strict=True)
        if resolved != path:
            raise GateError(f"Trust artifact path is not canonical: {path} -> {resolved}")
        if not stat.S_ISREG(leaf.st_mode):
            raise GateError(f"Trust artifact is not a regular file: {path}")
        return resolved, leaf

    def open(self, path: Path) -> tuple[int, bytes, dict[str, Any]]:
        resolved, leaf = self._canonical_regular_path(path)
        existing = self._opened.get(resolved)
        if existing is not None:
            descriptor, data, _, record = existing
            if data is None:
                raise GateError(f"Metadata-only artifact cannot be reopened for bytes: {resolved}")
            self.require_paths_still_bound()
            return descriptor, data, dict(record)

        flags = os.O_RDONLY | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(resolved, flags)
        try:
            before = os.fstat(descriptor)
            if not stat.S_ISREG(before.st_mode) or not _same_stat(before, leaf):
                raise GateError(f"Trust artifact changed before read: {resolved}")
            chunks: list[bytes] = []
            offset = 0
            while offset < before.st_size:
                block = os.pread(descriptor, min(8 * 1024 * 1024, before.st_size - offset), offset)
                if not block:
                    raise GateError(f"Short read from trust artifact: {resolved}")
                chunks.append(block)
                offset += len(block)
            data = b"".join(chunks)
            after = os.fstat(descriptor)
            live = os.lstat(resolved)
            if len(data) != before.st_size or not _same_stat(before, after) or not _same_stat(after, live):
                raise GateError(f"Trust artifact changed while reading: {resolved}")
            record = artifact_record(resolved, data, after)
            self._opened[resolved] = (descriptor, data, after, record)
            return descriptor, data, dict(record)
        except Exception:
            os.close(descriptor)
            raise

    def open_metadata(self, path: Path) -> tuple[int, dict[str, Any]]:
        if not hasattr(os, "O_PATH"):
            raise GateError("O_PATH is required for metadata-only key retirement binding")
        resolved, leaf = self._canonical_regular_path(path)
        existing = self._opened.get(resolved)
        if existing is not None:
            descriptor, _, _, record = existing
            self.require_paths_still_bound()
            return descriptor, dict(record)
        flags = os.O_PATH | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(resolved, flags)
        try:
            opened = os.fstat(descriptor)
            live = os.lstat(resolved)
            if not stat.S_ISREG(opened.st_mode) or not _same_stat(opened, leaf) or not _same_stat(opened, live):
                raise GateError(f"Metadata-only trust artifact changed: {resolved}")
            record = {
                "path": str(resolved),
                "mode": oct(stat.S_IMODE(opened.st_mode)),
            }
            self._opened[resolved] = (descriptor, None, opened, record)
            return descriptor, dict(record)
        except Exception:
            os.close(descriptor)
            raise

    def open_executable_alias(
        self,
        path: Path,
        *,
        expected_link_text: str,
        expected_resolved_path: Path,
    ) -> tuple[int, bytes, dict[str, Any]]:
        """Bind one reviewed executable symlink and the regular target it resolves to."""
        if not path.is_absolute() or not expected_resolved_path.is_absolute():
            raise GateError("Executable alias and target paths must be absolute")
        alias_stat = os.lstat(path)
        if not stat.S_ISLNK(alias_stat.st_mode) or os.readlink(path) != expected_link_text:
            raise GateError(f"Executable alias contract drift: {path}")
        resolved = path.resolve(strict=True)
        if resolved != expected_resolved_path:
            raise GateError(f"Executable alias resolved target drift: {path} -> {resolved}")
        if not hasattr(os, "O_PATH"):
            raise GateError("O_PATH is required for executable alias binding")
        alias_fd = os.open(path, os.O_PATH | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0))
        try:
            opened_alias = os.fstat(alias_fd)
            if not _same_stat(alias_stat, opened_alias):
                raise GateError(f"Executable alias changed while binding: {path}")
            target_fd, target_data, target_record = self.open(expected_resolved_path)
            self._aliases[path] = (
                alias_fd,
                opened_alias,
                expected_link_text,
                expected_resolved_path,
            )
            return target_fd, target_data, {
                **target_record,
                "path": str(path),
                "resolved_path": str(expected_resolved_path),
                "symlink_target": expected_link_text,
            }
        except Exception:
            os.close(alias_fd)
            raise

    def stat_for(self, path: Path) -> os.stat_result:
        resolved = path.resolve(strict=True)
        if resolved not in self._opened:
            raise GateError(f"Artifact was not bound: {resolved}")
        return self._opened[resolved][2]

    def require_paths_still_bound(self) -> None:
        for path, (descriptor, _, opened, _) in self._opened.items():
            descriptor_stat = os.fstat(descriptor)
            try:
                live = os.lstat(path)
            except FileNotFoundError as error:
                raise GateError(f"Bound artifact path disappeared: {path}") from error
            if stat.S_ISLNK(live.st_mode) or not _same_stat(descriptor_stat, opened) or not _same_stat(live, opened):
                raise GateError(f"Bound artifact path or descriptor drifted: {path}")
        for path, (descriptor, opened, link_text, resolved) in self._aliases.items():
            descriptor_stat = os.fstat(descriptor)
            live = os.lstat(path)
            if (
                not stat.S_ISLNK(live.st_mode)
                or not _same_stat(descriptor_stat, opened)
                or not _same_stat(live, opened)
                or os.readlink(path) != link_text
                or path.resolve(strict=True) != resolved
            ):
                raise GateError(f"Bound executable alias drifted: {path}")

    def retain_until_process_exit(self) -> None:
        if self._retained:
            raise GateError("Bound artifact reader was already retained")
        self.require_paths_still_bound()
        _PROCESS_LIFETIME_READERS.append(self)
        self._retained = True

    def close(self) -> None:
        if self._retained:
            return
        while self._aliases:
            _, (descriptor, _, _, _) = self._aliases.popitem()
            os.close(descriptor)
        while self._opened:
            _, (descriptor, _, _, _) = self._opened.popitem()
            os.close(descriptor)


def observed_command() -> list[str]:
    raw = Path("/proc/self/cmdline").read_bytes()
    if not raw.endswith(b"\0"):
        raise GateError("Direct process command is unavailable")
    return [os.fsdecode(token) for token in raw[:-1].split(b"\0")]


def expected_command() -> list[str]:
    return [str(PYTHON), "-I", "-B", str(Path(__file__).resolve(strict=True))]


def verifier_command() -> list[str]:
    return [str(PYTHON), "-I", "-B", str(CONTINUATION_VERIFIER), "--verify"]


def promotion_command(mode: str) -> list[str]:
    if mode not in PROMOTION_CACHE_ROOTS:
        raise GateError(f"Unsupported promotion mode: {mode}")
    return [
        str(PYTHON),
        "-I",
        "-B",
        "-X",
        f"pycache_prefix={PROMOTION_CACHE_ROOTS[mode]}",
        str(PROMOTION_TOOL),
        mode,
    ]


def fd_bound_verifier_command(python_descriptor: int, verifier_descriptor: int) -> list[str]:
    if python_descriptor < 3 or verifier_descriptor < 3:
        raise GateError("Python/verifier descriptor is not inheritable")
    return [
        f"/proc/self/fd/{python_descriptor}",
        "-I",
        "-B",
        f"/proc/self/fd/{verifier_descriptor}",
        "--verify",
    ]


def normalized_fd_bound_verifier_command() -> list[str]:
    return [
        "/proc/self/fd/<bound-python-fd>",
        "-I",
        "-B",
        "/proc/self/fd/<bound-verifier-fd>",
        "--verify",
    ]


def downstream_continuation_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(DOWNSTREAM_CONTINUATION),
        "--supervise-and-sign",
    ]


def supervised_child_command() -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        str(DOWNSTREAM_CONTINUATION),
        "--supervised-child",
    ]


def require_clean_runtime() -> None:
    if dict(os.environ) != EXPECTED_ENVIRONMENT:
        raise GateError("Runner gate environment is not the exact clean environment")
    if observed_command() != expected_command():
        raise GateError("Runner gate is not the canonical direct process")
    if sys.flags.isolated != 1 or sys.dont_write_bytecode is not True:
        raise GateError("Python isolated/no-bytecode flags are not active")
    if Path(sys.executable) != PYTHON:
        raise GateError("Python executable token differs from the canonical launcher alias")
    if len(sys.argv) != 1:
        raise GateError("Runner gate accepts no command-line arguments")


def require_absent(path: Path, label: str) -> None:
    if os.path.lexists(path):
        raise GateError(f"{label} output slot is occupied: {path}")


def absence_snapshot(paths: Sequence[Path], label: str) -> dict[str, bool]:
    observed: dict[str, bool] = {}
    for path in paths:
        absent = not os.path.lexists(path)
        observed[str(path)] = absent
        if not absent:
            raise GateError(f"{label} output slot is occupied: {path}")
    return observed


def strict_json_equal(observed: Any, expected: Any) -> bool:
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


def require_exact_keys(value: Any, expected: frozenset[str], label: str) -> None:
    if not isinstance(value, Mapping) or set(value) != set(expected):
        observed = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise GateError(
            f"{label} keys are not exact: observed={observed}, expected={sorted(expected)}"
        )


def require_lower_sha256(value: Any, label: str) -> None:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise GateError(f"{label} is not a lowercase SHA-256 string")


def exit_attestation_schema(name: str) -> dict[str, str]:
    return {"name": name, "version": EXIT_ATTESTATION_SCHEMA_VERSION}


def exit_attestation_announcement(
    path: Path, schema_name: str, public_key: Mapping[str, Any]
) -> dict[str, Any]:
    return {
        "path": str(path),
        "schema": exit_attestation_schema(schema_name),
        "signing_key": dict(public_key),
    }


def exit_attestation_signature_contract(
    *,
    path: Path,
    schema_name: str,
    signature: Path,
    public_key: Mapping[str, Any],
    private_key: Path,
    legacy_payload_signature: Path,
) -> dict[str, Any]:
    return {
        "attestation_schema": exit_attestation_schema(schema_name),
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


def require_basic_artifact_types(value: Any, label: str) -> None:
    if not isinstance(value, Mapping) or not BASIC_ARTIFACT_KEYS.issubset(value):
        raise GateError(f"{label} is not an artifact record")
    if not isinstance(value.get("path"), str) or not value["path"].startswith("/"):
        raise GateError(f"{label} path is not an absolute string")
    if type(value.get("size_bytes")) is not int or value["size_bytes"] < 0:
        raise GateError(f"{label} size_bytes is not a non-negative JSON integer")
    require_lower_sha256(value.get("sha256"), f"{label} sha256")
    mode = value.get("mode")
    if not isinstance(mode, str) or not mode.startswith("0o"):
        raise GateError(f"{label} mode is not an octal string")


def artifact_with_stat(
    record: Mapping[str, Any], opened: os.stat_result
) -> dict[str, Any]:
    return {
        **dict(record),
        "device": opened.st_dev,
        "inode": opened.st_ino,
        "link_count": opened.st_nlink,
        "mtime_ns": opened.st_mtime_ns,
        "ctime_ns": opened.st_ctime_ns,
    }


def require_full_stat_artifact_types(value: Any, label: str) -> None:
    require_exact_keys(value, FULL_STAT_ARTIFACT_KEYS, label)
    require_basic_artifact_types(value, label)
    for field in ("device", "inode", "link_count", "mtime_ns", "ctime_ns"):
        observed = value.get(field)
        if type(observed) is not int or observed < 0:
            raise GateError(f"{label} {field} is not a non-negative JSON integer")
    if value["inode"] <= 0 or value["link_count"] != 1:
        raise GateError(f"{label} is not a single-link regular artifact")


def require_exact_artifact_shape(value: Any, label: str, *, full_stat: bool) -> None:
    if full_stat:
        require_full_stat_artifact_types(value, label)
        return
    require_exact_keys(value, BASIC_ARTIFACT_KEYS, label)
    require_basic_artifact_types(value, label)


def validate_producer_session(
    value: Any, label: str, *, expected_mode: str
) -> dict[str, Any]:
    require_exact_keys(value, PRODUCER_SESSION_KEYS, label)
    require_lower_sha256(value.get("nonce_sha256"), f"{label} nonce_sha256")
    if (
        value.get("protocol") != "intersubmod.parent_fork_no_exec_child"
        or value.get("schema_version") != "1.0.0"
        or value.get("mode") != expected_mode
    ):
        raise GateError(f"{label} protocol/schema/mode contract drift")
    for field in (
        "parent_pid",
        "parent_start_time_ticks",
        "child_pid",
        "child_start_time_ticks",
    ):
        observed = value.get(field)
        if type(observed) is not int or observed <= 0:
            raise GateError(f"{label} {field} is not a positive JSON integer")
    if value["parent_pid"] <= 1 or value["child_pid"] <= 1:
        raise GateError(f"{label} producer PIDs are not process PIDs")
    if value["parent_pid"] == value["child_pid"]:
        raise GateError(f"{label} parent and child PIDs are identical")
    if value.get("fork_no_exec") is not True:
        raise GateError(f"{label} does not attest fork-without-exec")
    return dict(value)


def validate_child_wait(
    value: Any, producer_session: Mapping[str, Any], label: str
) -> dict[str, Any]:
    require_exact_keys(value, CHILD_WAIT_KEYS, label)
    expected = {
        "waited_pid": producer_session["child_pid"],
        "raw_wait_status": 0,
        "wifexited": True,
        "exit_status": 0,
        "wifsignaled": False,
        "terminating_signal": None,
        "core_dumped": False,
    }
    if not strict_json_equal(value, expected):
        raise GateError(f"{label} is not exact normal-exit-zero waitpid evidence")
    raw_status = value["raw_wait_status"]
    if (
        os.WIFEXITED(raw_status) is not True
        or os.WEXITSTATUS(raw_status) != 0
        or os.WIFSIGNALED(raw_status) is not False
        or (hasattr(os, "WCOREDUMP") and os.WCOREDUMP(raw_status) is not False)
    ):
        raise GateError(f"{label} raw wait status is not normal exit zero")
    return dict(value)


def validate_exit_attestation(
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
    require_exact_keys(payload, EXIT_ATTESTATION_KEYS, label)
    created_at = parse_utc_timestamp(
        payload.get("created_at_utc"), f"{label} created_at_utc"
    )
    session = validate_producer_session(
        payload.get("producer_session"),
        f"{label} producer_session",
        expected_mode=expected_command[-1],
    )
    child_wait = validate_child_wait(
        payload.get("child_wait"), session, f"{label} child_wait"
    )
    if (
        payload.get("authorization_id") != AUTHORIZATION_ID
        or payload.get("phase") != expected_phase
        or not strict_json_equal(payload.get("supervisor_command"), expected_command)
        or not strict_json_equal(
            payload.get("produced_artifacts"), expected_produced_artifacts
        )
        or not strict_json_equal(
            payload.get("producer_source"), expected_producer_source
        )
        or not strict_json_equal(
            payload.get("python_runtime"), expected_python_runtime
        )
        or not strict_json_equal(
            payload.get("signature_contract"), expected_signature_contract
        )
        or not strict_json_equal(
            expected_signature_contract.get("attestation_schema"),
            exit_attestation_schema(expected_schema),
        )
        or not strict_json_equal(payload.get("checks"), expected_checks)
        or payload.get("pass") is not True
        or payload.get("pass_semantics") != expected_pass_semantics
    ):
        raise GateError(f"{label} exact envelope contract drift")
    return {
        "created_at_utc": created_at,
        "producer_session": session,
        "child_wait": child_wait,
    }


def require_artifact_equal(observed: Any, expected: Mapping[str, Any], label: str) -> None:
    require_basic_artifact_types(observed, label)
    require_basic_artifact_types(expected, f"{label} expected")
    required = ("path", "size_bytes", "sha256", "mode")
    if any(
        not strict_json_equal(observed.get(key), expected.get(key)) for key in required
    ):
        raise GateError(f"{label} identity drift")


def require_metadata_equal(observed: Any, expected: Mapping[str, Any], label: str) -> None:
    if not isinstance(observed, Mapping) or any(
        not strict_json_equal(observed.get(key), expected.get(key))
        for key in ("path", "mode")
    ):
        raise GateError(f"{label} metadata identity drift")


def stream_evidence(data: bytes) -> dict[str, Any]:
    evidence: dict[str, Any] = {
        "size_bytes": len(data),
        "sha256": sha256_bytes(data),
    }
    try:
        evidence.update({"encoding": "utf-8", "text": data.decode("utf-8")})
    except UnicodeDecodeError:
        evidence.update(
            {
                "encoding": "base64",
                "base64": base64.b64encode(data).decode("ascii"),
            }
        )
    return evidence


def load_bound_authority_module(source_data: bytes, source_path: Path) -> types.ModuleType:
    module = types.ModuleType("_bound_release_source_authority_v7")
    module.__file__ = str(source_path)
    module.__package__ = ""
    code = compile(source_data, str(source_path), "exec", dont_inherit=True)
    exec(code, module.__dict__)
    _PROCESS_LIFETIME_MODULES.append(module)
    return module


def recursively_validate_artifact_relations(
    payload: Any,
    reader: BoundArtifactReader,
) -> dict[str, dict[str, Any]]:
    """Bind and validate every artifact-shaped absolute path in a receipt."""
    relations: dict[str, dict[str, Any]] = {}

    def visit(value: Any) -> None:
        if isinstance(value, Mapping):
            path_value = value.get("path")
            mode_value = value.get("mode")
            if isinstance(path_value, str) and path_value.startswith("/") and isinstance(mode_value, str):
                path = Path(path_value)
                if "sha256" in value or "size_bytes" in value:
                    _, _, live = reader.open(path)
                    for key in ("path", "size_bytes", "sha256", "mode"):
                        if key in value and not strict_json_equal(
                            value.get(key), live.get(key)
                        ):
                            raise GateError(f"Verification receipt relation drift: {path} field={key}")
                    opened = reader.stat_for(path)
                    stat_fields = {
                        "device": opened.st_dev,
                        "inode": opened.st_ino,
                        "link_count": opened.st_nlink,
                        "mtime_ns": opened.st_mtime_ns,
                        "ctime_ns": opened.st_ctime_ns,
                    }
                    for key, expected in stat_fields.items():
                        if key in value and not strict_json_equal(value.get(key), expected):
                            raise GateError(
                                f"Verification receipt stat relation drift: {path} field={key}"
                            )
                    relations[str(path)] = live
                elif mode_value == "0o0":
                    _, live = reader.open_metadata(path)
                    require_metadata_equal(value, live, f"verification relation {path}")
                    relations[str(path)] = live
            for child in value.values():
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)

    visit(payload)
    return relations


def validate_verification_receipt(
    payload: Mapping[str, Any],
    relations: Mapping[str, Mapping[str, Any]],
    *,
    expected_mode: str,
) -> None:
    checks = payload.get("checks")
    governance = payload.get("governance")
    if (
        set(payload) != set(VERIFICATION_RECEIPT_KEYS)
        or payload.get("schema_name") != EXPECTED_VERIFICATION_SCHEMA
        or payload.get("schema_version") != PROMOTION_SCHEMA_VERSION
        or payload.get("task_type") != "B_comprehensive_validation"
        or payload.get("mode") != expected_mode
        or not isinstance(checks, Mapping)
        or set(checks) != set(VERIFICATION_CHECK_KEYS)
        or not all(value is True for value in checks.values())
        or not isinstance(governance, Mapping)
        or set(governance) != set(VERIFICATION_GOVERNANCE_KEYS)
        or governance.get("original_completion_runner_end_to_end_pass_claimed") is not False
        or governance.get("original_verifier_reexecuted_successfully_claimed") is not False
        or governance.get("runner_only_gate_replay_evidence_required_before_downstream")
        is not True
        or governance.get("signed_downstream_continuation_revalidation_required")
        is not True
        or governance.get("verification_scope")
        != "promotion_integrity_and_continuation_gate_only"
        or payload.get("pass_semantics") != EXPECTED_VERIFICATION_PASS_SEMANTICS
        or type(payload.get("retained_critical_input_descriptor_count")) is not int
        or payload["retained_critical_input_descriptor_count"] <= 0
        or payload.get("pass") is not True
    ):
        raise GateError("Promotion verification receipt contract did not pass exactly")
    parse_utc_timestamp(
        payload.get("created_at_utc"), f"promotion verification {expected_mode} created_at_utc"
    )

    producer_exit_attestations = payload.get("producer_exit_attestations")
    require_exact_keys(
        producer_exit_attestations,
        frozenset({"prepare", "promote"}),
        "verification producer_exit_attestations",
    )
    for phase, expected_attestation, expected_signature, expected_public_key in (
        (
            "prepare",
            PREPARE_EXIT_ATTESTATION,
            AUTHORIZATION_SIGNATURE,
            AUTHORIZATION_PUBLIC_KEY,
        ),
        (
            "promote",
            PROMOTE_EXIT_ATTESTATION,
            PROMOTION_SIGNATURE,
            COMPLETION_PUBLIC_KEY,
        ),
    ):
        entry = producer_exit_attestations[phase]
        require_exact_keys(
            entry,
            PRODUCER_EXIT_RECEIPT_ENTRY_KEYS,
            f"verification {phase} exit-attestation entry",
        )
        for field, expected_path in (
            ("artifact", expected_attestation),
            ("signature", expected_signature),
            ("public_key", expected_public_key),
        ):
            relation = relations.get(str(expected_path))
            if relation is None:
                raise GateError(
                    f"Verification receipt omitted {phase} {field} relation: {expected_path}"
                )
            require_exact_artifact_shape(
                entry[field],
                f"verification {phase} {field}",
                full_stat=field == "artifact",
            )
            require_artifact_equal(
                entry[field], relation, f"verification {phase} {field}"
            )
        if entry.get("signature_verified") is not True:
            raise GateError(
                f"Verification receipt did not verify {phase} attestation signature"
            )

    authorization_summary = payload.get("authorization")
    completion_summary = payload.get("completion")
    if (
        not isinstance(authorization_summary, Mapping)
        or not isinstance(completion_summary, Mapping)
        or not strict_json_equal(
            authorization_summary.get("signature"),
            producer_exit_attestations["prepare"]["signature"],
        )
        or not strict_json_equal(
            completion_summary.get("signature"),
            producer_exit_attestations["promote"]["signature"],
        )
        or authorization_summary.get("signature_verified") is not True
        or completion_summary.get("signature_verified") is not True
    ):
        raise GateError("Verification receipt signature summaries are inconsistent")

    required_relation_paths = {
        str(AUTHORIZATION),
        str(PREPARE_EXIT_ATTESTATION),
        str(AUTHORIZATION_SIGNATURE),
        str(AUTHORIZATION_PUBLIC_KEY),
        str(AUTHORIZATION_PRIVATE_KEY),
        str(PROMOTION_RECEIPT),
        str(PROMOTE_EXIT_ATTESTATION),
        str(PROMOTION_SIGNATURE),
        str(COMPLETION_PUBLIC_KEY),
        str(COMPLETION_PRIVATE_KEY),
        str(HISTORICAL_SOURCE_RECEIPT),
        str(CANONICAL_SOURCE_RECEIPT),
        str(AUTHORITY_VALIDATOR),
    }
    missing = sorted(required_relation_paths - set(relations))
    if missing:
        raise GateError(f"Promotion verification receipt omitted trust relations: {missing}")


def require_fresh_and_recorded_verification_agree(
    fresh: Mapping[str, Any], recorded: Mapping[str, Any]
) -> None:
    comparable_fields = (
        "schema_name",
        "schema_version",
        "task_type",
        "authorization",
        "completion",
        "producer_exit_attestations",
        "source_receipt",
        "canonical_receipt",
        "trusted_signing_keys",
        "private_key_retirement",
        "historical_incident_disclosure",
        "release_source_authority_validator",
        "checks",
        "retained_critical_input_descriptor_count",
        "governance",
        "pass",
        "pass_semantics",
    )
    if any(
        not strict_json_equal(fresh.get(field), recorded.get(field))
        for field in comparable_fields
    ):
        raise GateError("Fresh --verify result differs from the recorded verification evidence")
    fresh_created = parse_utc_timestamp(
        fresh.get("created_at_utc"), "fresh promotion verification created_at_utc"
    )
    recorded_created = parse_utc_timestamp(
        recorded.get("created_at_utc"), "recorded promotion verification created_at_utc"
    )
    if fresh_created < recorded_created:
        raise GateError("Fresh promotion verification predates the recorded receipt")


def validate_authorized_gate_contract(
    authorization_payload: Mapping[str, Any],
    *,
    self_record: Mapping[str, Any],
    verifier_record: Mapping[str, Any],
    continuation_record: Mapping[str, Any],
    continuation_public_record: Mapping[str, Any],
    continuation_private_record: Mapping[str, Any],
    python_record: Mapping[str, Any],
    bash_record: Mapping[str, Any],
    runner_record: Mapping[str, Any],
    runner_prefix_sha256: str,
    runner_line_count: int,
) -> Mapping[str, Any]:
    reviewed = authorization_payload.get("reviewed_sources")
    if not isinstance(reviewed, Mapping):
        raise GateError("Signed authorization has no reviewed_sources mapping")
    require_artifact_equal(
        reviewed.get("runner_gate_replay"), self_record, "authorization runner gate source"
    )
    require_artifact_equal(
        reviewed.get("continuation_verifier"), verifier_record, "authorization verifier source"
    )
    require_artifact_equal(
        reviewed.get("downstream_continuation"),
        continuation_record,
        "authorization downstream continuation source",
    )

    commands = authorization_payload.get("commands")
    continuation = authorization_payload.get("continuation_gate")
    terminal_signature_contract = (
        continuation.get("terminal_signature_contract")
        if isinstance(continuation, Mapping)
        else None
    )
    authorized_absence = authorization_payload.get("downstream_output_absence")
    expected_terminal_signature_contract = {
        "algorithm": "Ed25519",
        "public_key": dict(continuation_public_record),
        "private_key": dict(continuation_private_record),
        "execution_receipt": str(CONTINUATION_TERMINAL_RECEIPT),
        "signed_artifact": str(CONTINUATION_EXIT_ATTESTATION),
        "signature": str(CONTINUATION_TERMINAL_SIGNATURE),
        "required_private_key_pre_signature_mode": "0o400",
        "required_private_key_post_signature_mode": "0o0",
        "supervisor_command": downstream_continuation_command(),
        "supervised_child_command": supervised_child_command(),
        "signed_terminal_verification_command": [
            str(PYTHON),
            "-I",
            "-B",
            str(DOWNSTREAM_CONTINUATION),
            "--verify-signed-terminal",
        ],
    }
    if (
        not isinstance(commands, Mapping)
        or not strict_json_equal(commands.get("runner_gate_replay"), expected_command())
        or not strict_json_equal(commands.get("continuation_verify"), verifier_command())
        or not strict_json_equal(
            commands.get("downstream_continuation"), downstream_continuation_command()
        )
        or not isinstance(continuation, Mapping)
        or not strict_json_equal(continuation.get("verifier"), verifier_record)
        or not strict_json_equal(continuation.get("runner_gate_replay"), self_record)
        or not strict_json_equal(
            continuation.get("downstream_continuation"), continuation_record
        )
        or not strict_json_equal(
            continuation.get("downstream_continuation_command"),
            downstream_continuation_command(),
        )
        or not strict_json_equal(
            continuation.get("supervised_child_command"), supervised_child_command()
        )
        or not strict_json_equal(
            continuation.get("supervisor_capability_protocol"),
            SUPERVISOR_CAPABILITY_PROTOCOL,
        )
        or not strict_json_equal(
            continuation.get("threat_model"), CONTINUATION_THREAT_MODEL
        )
        or continuation.get("verification_receipt") != str(PROMOTION_VERIFICATION_RECEIPT)
        or continuation.get("runner_gate_receipt") != str(OUTPUT_RECEIPT)
        or not strict_json_equal(
            terminal_signature_contract, expected_terminal_signature_contract
        )
        or continuation.get("policy") != CONTINUATION_POLICY
        or not strict_json_equal(
            authorized_absence, [str(path) for path in DOWNSTREAM_OUTPUT_SLOTS]
        )
    ):
        raise GateError("Signed authorization continuation-gate contract drift")
    if not strict_json_equal(dict(python_record), EXPECTED_RUNTIME_ARTIFACTS["python"]):
        raise GateError("Pinned Python launcher alias/target identity drift")
    if not strict_json_equal(dict(bash_record), EXPECTED_RUNTIME_ARTIFACTS["bash"]):
        raise GateError("Pinned Bash runtime identity drift")
    if (
        runner_prefix_sha256 != EXPECTED_RUNNER_PREFIX_SHA256
        or runner_line_count < EXPECTED_RUNNER_LINE_COUNT_MINIMUM
    ):
        raise GateError("Authorization-bound runner slice constant drift")
    return {
        "identity_anchor": (
            "runtime and slice constants are transitively bound by the signed exact SHA-256 "
            "of runner_gate_replay; completion_runner is independently bound by v7"
        ),
        "authorization_commands": {
            "continuation_verify": verifier_command(),
            "runner_gate_replay": expected_command(),
            "downstream_continuation": downstream_continuation_command(),
        },
        "runtime_artifacts": {
            "python": dict(python_record),
            "bash": dict(bash_record),
        },
        "completion_runner": dict(runner_record),
        "runner_slice": {
            "included_line_start": RUNNER_INCLUDED_LINE_START,
            "included_line_end": RUNNER_INCLUDED_LINE_END,
            "first_excluded_line": RUNNER_FIRST_EXCLUDED_LINE,
            "minimum_total_lines": EXPECTED_RUNNER_LINE_COUNT_MINIMUM,
            "observed_total_lines": runner_line_count,
            "sha256": runner_prefix_sha256,
        },
        "outputs": {"log": str(OUTPUT_LOG), "receipt": str(OUTPUT_RECEIPT)},
        "downstream_output_slots": [str(path) for path in DOWNSTREAM_OUTPUT_SLOTS],
        "downstream_not_executed": True,
    }


def link_fd_no_replace(source_fd: int, target_directory_fd: int, target_name: str) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    linkat = getattr(libc, "linkat", None)
    if linkat is None:
        raise GateError("linkat is required for terminal commit publication")
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
        raise GateError("No-replace publication failed: " + os.strerror(error_number))


def write_exclusive_readonly(
    path: Path,
    data: bytes,
    *,
    precommit: Callable[[], None] | None = None,
    terminal_commit: Callable[[], None] | None = None,
) -> dict[str, Any]:
    if not path.parent.is_absolute() or path.parent.resolve(strict=True) != path.parent:
        raise GateError(f"Output parent is not canonical: {path.parent}")
    if stat.S_ISLNK(os.lstat(path.parent).st_mode):
        raise GateError(f"Output parent is a symlink: {path.parent}")
    if not hasattr(os, "O_TMPFILE"):
        raise GateError("O_TMPFILE is required for staged publication")
    parent_flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0)
    if hasattr(os, "O_NOFOLLOW"):
        parent_flags |= os.O_NOFOLLOW
    parent_fd = os.open(path.parent, parent_flags)
    parent_opened = os.fstat(parent_fd)
    descriptor = -1
    reopened_fd = -1
    durable_fd = -1
    try:
        descriptor = os.open(
            ".", os.O_RDWR | os.O_TMPFILE | os.O_CLOEXEC, 0o400, dir_fd=parent_fd
        )
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise GateError(f"Short write while publishing {path}")
            offset += written
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 0
            or opened.st_size != len(data)
            or stat.S_IMODE(opened.st_mode) != 0o444
            or os.pread(descriptor, len(data), 0) != data
        ):
            raise GateError(f"Staged immutable publication verification failed: {path}")
        if precommit is not None:
            precommit()
        parent_live = os.stat(path.parent, follow_symlinks=False)
        if (
            parent_opened.st_dev,
            parent_opened.st_ino,
            parent_opened.st_mode,
        ) != (parent_live.st_dev, parent_live.st_ino, parent_live.st_mode):
            raise GateError(f"Staged publication parent path rebound: {path.parent}")
        record = artifact_record(path, data, opened)
        require_absent(path, "staged publication target")
        link_fd_no_replace(descriptor, parent_fd, path.name)
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
            raise GateError(f"Committed publication path/inode drift: {path}")
        os.fsync(parent_fd)
        parent_after = os.stat(path.parent, follow_symlinks=False)
        if (
            parent_opened.st_dev,
            parent_opened.st_ino,
            parent_opened.st_mode,
        ) != (parent_after.st_dev, parent_after.st_ino, parent_after.st_mode):
            raise GateError(f"Publication parent changed after commit: {path.parent}")
        if path.resolve(strict=True) != path:
            raise GateError(f"Committed publication path is not canonical: {path}")
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
            raise GateError(f"Post-fsync publication path/inode drift: {path}")
        if terminal_commit is not None:
            terminal_commit()
            raise GateError("Terminal publication callback returned")
        return record
    finally:
        if durable_fd >= 0:
            os.close(durable_fd)
        if reopened_fd >= 0:
            os.close(reopened_fd)
        if descriptor >= 0:
            os.close(descriptor)
        os.close(parent_fd)


def build_log(
    verifier_result: subprocess.CompletedProcess[bytes],
    replay_result: subprocess.CompletedProcess[bytes],
    runner_prefix_sha256: str,
) -> bytes:
    def decoded(data: bytes) -> str:
        return data.decode("utf-8", errors="backslashreplace")

    sections = [
        "InterSubMod M2v5 runner-only gate replay v1\n",
        "CLAIM: This does not assert that the original completion runner passed end to end.\n",
        "CLAIM: This does not assert that the original tumor-REF verifier was reexecuted successfully.\n",
        "CLAIM: This pathname replay is non-authoritative and does not authorize downstream execution.\n",
        f"Authorized promotion verifier command: {json.dumps(verifier_command())}\n",
        f"Executed promotion verifier command: {json.dumps(normalized_fd_bound_verifier_command())}\n",
        f"Promotion verifier exit code: {verifier_result.returncode}\n",
        "--- promotion verifier stdout ---\n",
        decoded(verifier_result.stdout),
        "\n--- promotion verifier stderr ---\n",
        decoded(verifier_result.stderr),
        "\n--- runner-only replay ---\n",
        f"Runner: {COMPLETION_RUNNER}\n",
        f"Included physical lines: 1-358\n",
        f"First excluded physical line: 359\n",
        f"Runner prefix SHA-256: {runner_prefix_sha256}\n",
        f"Bash command: {json.dumps([f'/proc/self/fd/<bound-bash-fd>', '-s'])}\n",
        f"Runner-only exit code: {replay_result.returncode}\n",
        "--- runner-only stdout ---\n",
        decoded(replay_result.stdout),
        "\n--- runner-only stderr ---\n",
        decoded(replay_result.stderr),
        "\n--- terminal status ---\nPASS EVIDENCE ONLY: runner lines 1-358; downstream not executed or authorized.\n",
    ]
    return "".join(sections).encode("utf-8")


def main() -> int:
    require_clean_runtime()
    require_absent(
        LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
        "legacy authorization-payload signature",
    )
    require_absent(
        LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
        "legacy completion-payload signature",
    )
    require_absent(OUTPUT_LOG, "runner gate log")
    require_absent(OUTPUT_RECEIPT, "runner gate receipt")
    require_absent(PROMOTION_INCIDENT, "promotion failure incident")
    downstream_before = absence_snapshot(DOWNSTREAM_OUTPUT_SLOTS, "pre-verifier downstream")

    reader = BoundArtifactReader()
    try:
        self_fd, self_data, self_record = reader.open(Path(__file__).resolve(strict=True))
        require_executing_replayer_matches_bound_source(self_data)
        verifier_fd, verifier_data, verifier_record = reader.open(CONTINUATION_VERIFIER)
        continuation_fd, continuation_data, continuation_record = reader.open(
            DOWNSTREAM_CONTINUATION
        )
        del continuation_fd, continuation_data
        python_fd, _, python_record = reader.open_executable_alias(
            PYTHON,
            expected_link_text="python3.11",
            expected_resolved_path=PYTHON_TARGET,
        )
        bash_fd, _, bash_record = reader.open(BASH)
        openssl_fd, _, openssl_record = reader.open(OPENSSL)
        reviewed_runtime_sources: dict[str, dict[str, Any]] = {}
        for role, (path, expected_sha256, expected_mode) in (
            REVIEWED_RUNTIME_SOURCE_CONTRACTS.items()
        ):
            _, _, record = reader.open(path)
            if (
                record.get("sha256") != expected_sha256
                or record.get("mode") != expected_mode
            ):
                raise GateError(f"Reviewed runtime source identity drift: {role}")
            reviewed_runtime_sources[role] = record
        observed_asset_parts = tuple(
            sorted(
                path
                for path in PORTABLE_ASSET_ROOT.iterdir()
                if path.name.startswith("portable-artifact-reader.html.gz.b64.part")
            )
        )
        if observed_asset_parts != PORTABLE_READER_ASSET_PARTS:
            raise GateError("Portable reader asset part membership drift")

        authorization_fd, authorization_data, authorization_record = reader.open(AUTHORIZATION)
        prepare_exit_fd, prepare_exit_data, prepare_exit_record = reader.open(
            PREPARE_EXIT_ATTESTATION
        )
        authorization_signature_fd, _, authorization_signature_record = reader.open(
            AUTHORIZATION_SIGNATURE
        )
        authorization_public_fd, _, authorization_public_record = reader.open(
            AUTHORIZATION_PUBLIC_KEY
        )
        _, authorization_private_record = reader.open_metadata(AUTHORIZATION_PRIVATE_KEY)
        promotion_fd, promotion_data, promotion_record = reader.open(PROMOTION_RECEIPT)
        promote_exit_fd, promote_exit_data, promote_exit_record = reader.open(
            PROMOTE_EXIT_ATTESTATION
        )
        promotion_signature_fd, _, promotion_signature_record = reader.open(PROMOTION_SIGNATURE)
        completion_public_fd, _, completion_public_record = reader.open(COMPLETION_PUBLIC_KEY)
        _, completion_private_record = reader.open_metadata(COMPLETION_PRIVATE_KEY)
        _, _, continuation_public_record = reader.open(CONTINUATION_PUBLIC_KEY)
        _, continuation_private_record = reader.open_metadata(CONTINUATION_PRIVATE_KEY)
        _, source_data, source_record = reader.open(HISTORICAL_SOURCE_RECEIPT)
        _, canonical_data, canonical_record = reader.open(CANONICAL_SOURCE_RECEIPT)
        verification_fd, verification_data, verification_record = reader.open(
            PROMOTION_VERIFICATION_RECEIPT
        )

        v7_fd, v7_data, v7_record = reader.open(V7_AUTHORITY)
        _, v7_approval_data, v7_approval_record = reader.open(V7_APPROVAL)
        _, _, v7_signature_record = reader.open(V7_SIGNATURE)
        _, _, v7_public_record = reader.open(V7_PUBLIC_KEY)
        _, v7_private_record = reader.open_metadata(V7_PRIVATE_KEY)
        _, authority_source_data, authority_source_record = reader.open(AUTHORITY_VALIDATOR)
        runner_fd, runner_data, runner_record = reader.open(COMPLETION_RUNNER)

        if any(
            record.get("mode") != "0o444"
            for record in (
                self_record,
                verifier_record,
                continuation_record,
                authorization_record,
                prepare_exit_record,
                authorization_signature_record,
                authorization_public_record,
                promotion_record,
                promote_exit_record,
                promotion_signature_record,
                completion_public_record,
                continuation_public_record,
                source_record,
                canonical_record,
                verification_record,
                v7_record,
                v7_approval_record,
                v7_signature_record,
                v7_public_record,
                authority_source_record,
                runner_record,
            )
        ):
            raise GateError("A required immutable artifact is not mode 0444")
        if authorization_private_record.get("mode") != "0o0":
            raise GateError("Promotion authorization private key is not retired to mode 0000")
        if completion_private_record.get("mode") != "0o0":
            raise GateError("Promotion completion private key is not retired to mode 0000")
        if continuation_private_record.get("mode") != "0o400":
            raise GateError("Continuation terminal private key is not pre-signing mode 0400")
        if (
            authorization_signature_record.get("size_bytes") != 64
            or promotion_signature_record.get("size_bytes") != 64
            or authorization_public_record.get("sha256")
            != EXPECTED_AUTHORIZATION_PUBLIC_KEY_SHA256
            or completion_public_record.get("sha256")
            != EXPECTED_COMPLETION_PUBLIC_KEY_SHA256
            or continuation_public_record.get("sha256")
            != EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256
        ):
            raise GateError("Promotion signature or pinned public-key identity drift")
        if v7_private_record.get("mode") != "0o0":
            raise GateError("v7 authority private key is not retired to mode 0000")

        # The continuation verifier is the first executed gate. Its source and
        # interpreter paths remain descriptor-bound across the direct process.
        reader.require_paths_still_bound()
        verifier_result = subprocess.run(
            fd_bound_verifier_command(python_fd, verifier_fd),
            cwd=REPO_ROOT,
            env=EXPECTED_ENVIRONMENT,
            pass_fds=(python_fd, verifier_fd),
            capture_output=True,
            check=False,
            timeout=300,
        )
        reader.require_paths_still_bound()
        if verifier_result.returncode != 0:
            raise GateError(
                "Promotion continuation verifier failed: "
                f"exit={verifier_result.returncode} stderr_sha256={sha256_bytes(verifier_result.stderr)}"
            )
        verifier_stdout = parse_json(verifier_result.stdout, "promotion verifier stdout")
        if verifier_stdout.get("pass") is not True:
            raise GateError("Promotion continuation verifier stdout did not declare pass=true")
        verification_payload = parse_json(
            verification_data, "promotion verification receipt"
        )
        fresh_verification_relations = recursively_validate_artifact_relations(
            verifier_stdout, reader
        )
        verification_relations = recursively_validate_artifact_relations(
            verification_payload, reader
        )
        validate_verification_receipt(
            verifier_stdout,
            fresh_verification_relations,
            expected_mode="verify",
        )
        validate_verification_receipt(
            verification_payload,
            verification_relations,
            expected_mode="verify-and-record",
        )
        require_fresh_and_recorded_verification_agree(
            verifier_stdout, verification_payload
        )

        authority_payload = parse_json(v7_data, "v7 source authority")
        approval_payload = parse_json(v7_approval_data, "v7 source authority approval")
        if (
            authority_payload.get("authority_id") != EXPECTED_V7_AUTHORITY_ID
            or authority_payload.get("approval_status") != "APPROVED_FOR_FULL_TASK_B_RUN"
            or approval_payload.get("authority_id") != EXPECTED_V7_AUTHORITY_ID
            or approval_payload.get("approval_status") != "APPROVED_FOR_FULL_TASK_B_RUN"
            or not strict_json_equal(
                approval_payload.get("authority_manifest"), v7_record
            )
            or not strict_json_equal(approval_payload.get("public_key"), v7_public_record)
        ):
            raise GateError("Bound v7 authority/approval relation drift")
        require_artifact_equal(
            authority_payload.get("sources", {}).get("release_source_authority_validator"),
            authority_source_record,
            "v7 authority validator",
        )
        require_artifact_equal(
            authority_payload.get("sources", {}).get("completion_runner"),
            runner_record,
            "v7 completion runner",
        )

        authority_module = load_bound_authority_module(
            authority_source_data, AUTHORITY_VALIDATOR
        )
        v7_validation = authority_module.validate_release_source_authority(
            authority_module.EXPECTED_SOURCE_PATHS
        )
        if (
            v7_validation.get("pass") is not True
            or v7_validation.get("authority_id") != EXPECTED_V7_AUTHORITY_ID
            or not strict_json_equal(
                v7_validation.get("authority_manifest"), v7_record
            )
            or not strict_json_equal(
                v7_validation.get("detached_approval"), v7_approval_record
            )
            or not strict_json_equal(
                v7_validation.get("detached_approval_signature"), v7_signature_record
            )
            or not strict_json_equal(
                v7_validation.get("external_public_key"), v7_public_record
            )
            or v7_validation.get("retired_private_key", {}).get("path")
            != str(V7_PRIVATE_KEY)
            or v7_validation.get("retired_private_key", {}).get("mode") != "0o0"
            or not strict_json_equal(
                v7_validation.get("signature_verifier"), openssl_record
            )
        ):
            raise GateError("v7 source authority process-lifetime validation failed")

        authorization_payload = parse_json(authorization_data, "promotion authorization")
        promotion_payload = parse_json(promotion_data, "promotion completion receipt")
        prepare_exit_payload = parse_json(
            prepare_exit_data, "prepare producer exit attestation"
        )
        promote_exit_payload = parse_json(
            promote_exit_data, "promote producer exit attestation"
        )
        require_exact_keys(
            authorization_payload, AUTHORIZATION_KEYS, "promotion authorization"
        )
        require_exact_keys(
            promotion_payload, COMPLETION_KEYS, "promotion completion receipt"
        )
        promotion_tool_record = reader.open(PROMOTION_TOOL)[2]
        expected_reviewed_sources = {
            "promotion_tool": promotion_tool_record,
            "continuation_verifier": verifier_record,
            "runner_gate_replay": self_record,
            "downstream_continuation": continuation_record,
            **reviewed_runtime_sources,
        }
        if not strict_json_equal(
            authorization_payload.get("reviewed_sources"), expected_reviewed_sources
        ):
            raise GateError("Signed reviewed source/runtime authority drift")
        authorization_created = parse_utc_timestamp(
            authorization_payload.get("created_at_utc"),
            "promotion authorization created_at_utc",
        )
        promotion_created = parse_utc_timestamp(
            promotion_payload.get("created_at_utc"),
            "promotion completion created_at_utc",
        )
        if promotion_created < authorization_created:
            raise GateError("Promotion completion predates signed authorization")
        expected_authorization_signature_contract = (
            exit_attestation_signature_contract(
                path=PREPARE_EXIT_ATTESTATION,
                schema_name=PREPARE_EXIT_ATTESTATION_SCHEMA,
                signature=AUTHORIZATION_SIGNATURE,
                public_key=authorization_public_record,
                private_key=AUTHORIZATION_PRIVATE_KEY,
                legacy_payload_signature=LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
            )
        )
        expected_authorized_completion_signature_contract = (
            exit_attestation_signature_contract(
                path=PROMOTE_EXIT_ATTESTATION,
                schema_name=PROMOTE_EXIT_ATTESTATION_SCHEMA,
                signature=PROMOTION_SIGNATURE,
                public_key=completion_public_record,
                private_key=COMPLETION_PRIVATE_KEY,
                legacy_payload_signature=LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
            )
        )
        expected_completion_signature_contract = (
            expected_authorized_completion_signature_contract
        )
        authorization_full_record = artifact_with_stat(
            authorization_record, reader.stat_for(AUTHORIZATION)
        )
        prepare_exit_full_record = artifact_with_stat(
            prepare_exit_record, reader.stat_for(PREPARE_EXIT_ATTESTATION)
        )
        authorization_signature_full_record = artifact_with_stat(
            authorization_signature_record, reader.stat_for(AUTHORIZATION_SIGNATURE)
        )
        promotion_full_record = artifact_with_stat(
            promotion_record, reader.stat_for(PROMOTION_RECEIPT)
        )
        promote_exit_full_record = artifact_with_stat(
            promote_exit_record, reader.stat_for(PROMOTE_EXIT_ATTESTATION)
        )
        promotion_signature_full_record = artifact_with_stat(
            promotion_signature_record, reader.stat_for(PROMOTION_SIGNATURE)
        )
        canonical_full_record = artifact_with_stat(
            canonical_record, reader.stat_for(CANONICAL_SOURCE_RECEIPT)
        )
        promotion_tool_full_record = artifact_with_stat(
            promotion_tool_record, reader.stat_for(PROMOTION_TOOL)
        )
        primary_python_full_record = artifact_with_stat(
            reviewed_runtime_sources["primary_python_runtime"],
            reader.stat_for(PYTHON_TARGET),
        )
        prepare_exit_validation = validate_exit_attestation(
            prepare_exit_payload,
            label="prepare producer exit attestation",
            expected_schema=PREPARE_EXIT_ATTESTATION_SCHEMA,
            expected_phase="prepare_authorization",
            expected_command=promotion_command("--prepare-authorization"),
            expected_produced_artifacts={"authorization": authorization_full_record},
            expected_producer_source=promotion_tool_full_record,
            expected_python_runtime=primary_python_full_record,
            expected_signature_contract=expected_authorization_signature_contract,
            expected_checks=PREPARE_EXIT_ATTESTATION_CHECKS,
            expected_pass_semantics=PREPARE_EXIT_ATTESTATION_PASS_SEMANTICS,
        )
        expected_signed_prepare_chain = {
            "artifact": authorization_full_record,
            "producer_exit_attestation": prepare_exit_full_record,
            "signature": authorization_signature_full_record,
            "public_key": authorization_public_record,
            "retired_private_key": authorization_private_record,
            "producer_session": prepare_exit_validation["producer_session"],
            "signature_verified": True,
        }
        promote_exit_validation = validate_exit_attestation(
            promote_exit_payload,
            label="promote producer exit attestation",
            expected_schema=PROMOTE_EXIT_ATTESTATION_SCHEMA,
            expected_phase="promote",
            expected_command=promotion_command("--promote"),
            expected_produced_artifacts={
                "completion_receipt": promotion_full_record,
                "canonical_receipt": canonical_full_record,
                "prepare_exit_attestation": prepare_exit_full_record,
                "prepare_exit_attestation_signature": (
                    authorization_signature_full_record
                ),
            },
            expected_producer_source=promotion_tool_full_record,
            expected_python_runtime=primary_python_full_record,
            expected_signature_contract=expected_completion_signature_contract,
            expected_checks=PROMOTE_EXIT_ATTESTATION_CHECKS,
            expected_pass_semantics=PROMOTE_EXIT_ATTESTATION_PASS_SEMANTICS,
        )
        if (
            prepare_exit_validation["created_at_utc"] < authorization_created
            or promote_exit_validation["created_at_utc"] < promotion_created
            or promote_exit_validation["created_at_utc"]
            < prepare_exit_validation["created_at_utc"]
            or promote_exit_validation["producer_session"]["nonce_sha256"]
            == prepare_exit_validation["producer_session"]["nonce_sha256"]
        ):
            raise GateError("Promotion producer attestation chronology/session drift")
        recorded_verification_created = parse_utc_timestamp(
            verification_payload.get("created_at_utc"),
            "recorded promotion verification created_at_utc",
        )
        if recorded_verification_created < promotion_created:
            raise GateError("Recorded promotion verification predates completion")
        authorization_relations = recursively_validate_artifact_relations(
            authorization_payload, reader
        )
        promotion_relations = recursively_validate_artifact_relations(
            promotion_payload, reader
        )
        expected_prepare_exit_announcement = exit_attestation_announcement(
            PREPARE_EXIT_ATTESTATION,
            PREPARE_EXIT_ATTESTATION_SCHEMA,
            authorization_public_record,
        )
        expected_promote_exit_announcement = exit_attestation_announcement(
            PROMOTE_EXIT_ATTESTATION,
            PROMOTE_EXIT_ATTESTATION_SCHEMA,
            completion_public_record,
        )
        if (
            authorization_payload.get("schema_name") != EXPECTED_AUTHORIZATION_SCHEMA
            or authorization_payload.get("schema_version") != PROMOTION_SCHEMA_VERSION
            or authorization_payload.get("authorization_id") != AUTHORIZATION_ID
            or authorization_payload.get("task_type") != "B_comprehensive_validation"
            or authorization_payload.get("status")
            != "AUTHORIZED_FOR_SINGLE_BYTE_IDENTICAL_PROMOTION"
            or authorization_payload.get("scope")
            != "canonical_path_publication_only_no_scientific_payload_change"
            or not strict_json_equal(
                authorization_payload.get("producer_exit_attestation"),
                expected_prepare_exit_announcement,
            )
            or authorization_payload.get("pass") is not True
            or not strict_json_equal(
                authorization_payload.get("producer_session"),
                prepare_exit_validation["producer_session"],
            )
            or not strict_json_equal(
                authorization_payload.get("authorization_signature_contract"),
                expected_authorization_signature_contract,
            )
            or not strict_json_equal(
                authorization_payload.get("completion_signature_contract"),
                expected_authorized_completion_signature_contract,
            )
            or not strict_json_equal(
                authorization_payload.get("checks"), AUTHORIZATION_CHECKS
            )
            or authorization_payload.get("pass_semantics")
            != AUTHORIZATION_PASS_SEMANTICS
            or promotion_payload.get("schema_name") != EXPECTED_PROMOTION_SCHEMA
            or promotion_payload.get("schema_version") != PROMOTION_SCHEMA_VERSION
            or promotion_payload.get("authorization_id") != AUTHORIZATION_ID
            or promotion_payload.get("task_type") != "B_comprehensive_validation"
            or promotion_payload.get("status")
            != "PROVISIONAL_WRITE_AHEAD_CANONICAL_PUBLICATION_INTENT"
            or promotion_payload.get("scope")
            != "canonical_path_publication_only_no_scientific_payload_change"
            or not strict_json_equal(
                promotion_payload.get("producer_exit_attestation"),
                expected_promote_exit_announcement,
            )
            or not strict_json_equal(
                promotion_payload.get("command"), promotion_command("--promote")
            )
            or not strict_json_equal(
                promotion_payload.get("producer_session"),
                promote_exit_validation["producer_session"],
            )
            or not strict_json_equal(
                promotion_payload.get("authorization"), expected_signed_prepare_chain
            )
            or not strict_json_equal(
                promotion_payload.get("code"), expected_reviewed_sources
            )
            or not strict_json_equal(
                promotion_payload.get("signature_contract"),
                expected_completion_signature_contract,
            )
            or not strict_json_equal(
                promotion_payload.get("checks"), COMPLETION_CHECKS
            )
            or promotion_payload.get("pass_semantics") != COMPLETION_PASS_SEMANTICS
            or promotion_payload.get("governance", {}).get(
                "write_ahead_receipt_marks_transaction_incomplete_until"
            )
            != "detached_completion_signature_plus_live_canonical_verification"
            or promotion_payload.get("governance", {}).get(
                "canonical_alone_grants_downstream_authority"
            )
            is not False
            or promotion_payload.get("pass") is not True
        ):
            raise GateError("Authorization or promotion completion contract drift")
        require_absent(
            LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
            "legacy authorization-payload signature before detached verification",
        )
        require_absent(
            LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
            "legacy completion-payload signature before detached verification",
        )
        if not authority_module.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=prepare_exit_fd,
            public_key_fd=authorization_public_fd,
            signature_fd=authorization_signature_fd,
        ):
            raise GateError("Prepare exit-attestation detached signature failed")
        if not authority_module.verify_ed25519_signature_fds(
            openssl_fd=openssl_fd,
            data_fd=promote_exit_fd,
            public_key_fd=completion_public_fd,
            signature_fd=promotion_signature_fd,
        ):
            raise GateError("Promote exit-attestation detached signature failed")
        require_artifact_equal(
            authorization_payload.get("trusted_signing_keys", {}).get(
                "authorization_public_key"
            ),
            authorization_public_record,
            "authorization trusted authorization key",
        )
        require_artifact_equal(
            authorization_payload.get("trusted_signing_keys", {}).get(
                "completion_public_key"
            ),
            completion_public_record,
            "authorization trusted completion key",
        )
        require_artifact_equal(
            authorization_payload.get("trusted_signing_keys", {}).get(
                "continuation_public_key"
            ),
            continuation_public_record,
            "authorization trusted continuation key",
        )
        require_artifact_equal(
            promotion_payload.get("authorization", {}).get("artifact"),
            authorization_record,
            "promotion authorization evidence",
        )
        require_artifact_equal(
            promotion_payload.get("authorization", {}).get("signature"),
            authorization_signature_record,
            "promotion authorization signature evidence",
        )
        require_artifact_equal(
            promotion_payload.get("authorization", {}).get("public_key"),
            authorization_public_record,
            "promotion authorization key evidence",
        )
        require_artifact_equal(
            promotion_payload.get("signature_contract", {}).get("public_key"),
            completion_public_record,
            "promotion completion key contract",
        )
        require_artifact_equal(
            promotion_payload.get("source_receipt"),
            source_record,
            "promotion historical source receipt",
        )
        require_artifact_equal(
            promotion_payload.get("canonical_receipt"),
            canonical_record,
            "promotion canonical source receipt",
        )
        if (
            source_data != canonical_data
            or source_record["sha256"] != canonical_record["sha256"]
            or reader.stat_for(HISTORICAL_SOURCE_RECEIPT).st_ino
            == reader.stat_for(CANONICAL_SOURCE_RECEIPT).st_ino
            or promotion_payload.get("historical_incident_disclosure", {}).get(
                "original_completion_runner_end_to_end_pass_claimed"
            )
            is not False
            or promotion_payload.get("historical_incident_disclosure", {}).get(
                "original_post_run_verifier_success_claimed"
            )
            is not False
            or promotion_payload.get("governance", {}).get(
                "continuation_verification_receipt"
            )
            != str(PROMOTION_VERIFICATION_RECEIPT)
            or promotion_payload.get("governance", {}).get("runner_gate_receipt")
            != str(OUTPUT_RECEIPT)
        ):
            raise GateError("Canonical promotion bytes or claim ceiling drift")

        runner_lines = runner_data.splitlines(keepends=True)
        if len(runner_lines) < EXPECTED_RUNNER_LINE_COUNT_MINIMUM:
            raise GateError(f"Completion runner has only {len(runner_lines)} physical lines")
        if not all(line.endswith((b"\n", b"\r\n")) for line in runner_lines[:-1]):
            raise GateError("Completion runner contains a non-terminal unterminated physical line")
        runner_prefix = b"".join(runner_lines[:RUNNER_INCLUDED_LINE_END])
        runner_prefix_sha256 = sha256_bytes(runner_prefix)
        if runner_lines[RUNNER_FIRST_EXCLUDED_LINE - 1] not in (b"\n", b"\r\n"):
            raise GateError("Physical line 359 is not the reviewed blank downstream boundary")
        if not runner_lines[RUNNER_FIRST_EXCLUDED_LINE].startswith(b"for path in \\"):
            raise GateError("Physical line 360 is not the reviewed downstream output-slot boundary")
        if any(
            marker in runner_prefix
            for marker in (
                b"/usr/bin/mkdir --mode=0700 -- \"${PYTHON_CACHE_ROOT}\"",
                b"exec > >(tee \"${LOG_PATH}\")",
                b"Verify bounded tumor-REF recovery source identity",
                b"Run strict multi-seed and column-null robustness",
                b"Run matched-normal controls for formal G2 selection",
                b"Attach authority-locked CN and fit-local CCF annotations",
                b"Build source-attested final all-sSNV dataset",
            )
        ):
            raise GateError("Runner prefix contains a downstream write/execution marker")

        runtime_contract = validate_authorized_gate_contract(
            authorization_payload,
            self_record=self_record,
            verifier_record=verifier_record,
            continuation_record=continuation_record,
            continuation_public_record=continuation_public_record,
            continuation_private_record=continuation_private_record,
            python_record=python_record,
            bash_record=bash_record,
            runner_record=runner_record,
            runner_prefix_sha256=runner_prefix_sha256,
            runner_line_count=len(runner_lines),
        )
        require_artifact_equal(
            promotion_payload.get("code", {}).get("runner_gate_replay"),
            self_record,
            "promotion runner gate source",
        )
        require_artifact_equal(
            promotion_payload.get("code", {}).get("continuation_verifier"),
            verifier_record,
            "promotion continuation verifier source",
        )
        require_artifact_equal(
            promotion_payload.get("code", {}).get("downstream_continuation"),
            continuation_record,
            "promotion downstream continuation source",
        )
        require_artifact_equal(
            authorization_relations.get(str(V7_AUTHORITY)),
            v7_record,
            "authorization v7 authority",
        )
        require_artifact_equal(
            authorization_relations.get(str(V7_APPROVAL)),
            v7_approval_record,
            "authorization v7 approval",
        )
        require_artifact_equal(
            authorization_relations.get(str(V7_SIGNATURE)),
            v7_signature_record,
            "authorization v7 signature",
        )
        if (
            authorization_payload.get("evidence", {})
            .get("current_v7_runtime_validation", {})
            .get("authority_id")
            != EXPECTED_V7_AUTHORITY_ID
            or authorization_payload.get("evidence", {})
            .get("current_v7_runtime_validation", {})
            .get("source_set_sha256")
            != v7_validation.get("source_set_sha256")
            or str(V7_AUTHORITY) not in promotion_relations
            or str(V7_APPROVAL) not in promotion_relations
            or str(V7_SIGNATURE) not in promotion_relations
        ):
            raise GateError("Authorization v7 transitive identity drift")

        downstream_after_verifier = absence_snapshot(
            DOWNSTREAM_OUTPUT_SLOTS, "post-verifier downstream"
        )
        reader.require_paths_still_bound()
        replay_result = subprocess.run(
            [f"/proc/self/fd/{bash_fd}", "-s"],
            input=runner_prefix,
            cwd=REPO_ROOT,
            env=EXPECTED_ENVIRONMENT,
            pass_fds=(bash_fd,),
            capture_output=True,
            check=False,
            timeout=900,
        )
        reader.require_paths_still_bound()
        downstream_after_replay = absence_snapshot(
            DOWNSTREAM_OUTPUT_SLOTS, "post-runner-replay downstream"
        )
        if replay_result.returncode != 0:
            raise GateError(
                "Runner-only preflight replay failed: "
                f"exit={replay_result.returncode} stderr_sha256={sha256_bytes(replay_result.stderr)}"
            )

        shell_observation_checks = {
            "executing_replayer_definitions_match_bound_reviewed_source": True,
            "exact_runner_prefix_bytes_lines_1_358": True,
            "shell_pathname_checks_exited_zero_non_authoritative": True,
            "parent_did_not_treat_shell_result_as_downstream_authority": True,
            "signed_downstream_continuation_must_revalidate_bound_inputs": True,
            "runner_prefix_exit_zero_under_set_Eeuo_pipefail": True,
        }
        log_data = build_log(verifier_result, replay_result, runner_prefix_sha256)
        reader.require_paths_still_bound()
        def log_precommit() -> None:
            require_absent(PROMOTION_INCIDENT, "promotion failure incident before log commit")
            require_absent(
                LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
                "legacy authorization-payload signature before log commit",
            )
            require_absent(
                LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
                "legacy completion-payload signature before log commit",
            )
            absence_snapshot(DOWNSTREAM_OUTPUT_SLOTS, "pre-runner-log commit downstream")
            reader.require_paths_still_bound()

        log_record = write_exclusive_readonly(
            OUTPUT_LOG, log_data, precommit=log_precommit
        )
        _, observed_log_data, bound_log_record = reader.open(OUTPUT_LOG)
        if observed_log_data != log_data or bound_log_record != log_record:
            raise GateError("Published runner gate log identity drift")

        replay_created_at = now_utc()
        replay_created = parse_utc_timestamp(
            replay_created_at, "runner-only replay created_at_utc"
        )
        fresh_verification_created = parse_utc_timestamp(
            verifier_stdout.get("created_at_utc"), "fresh verification created_at_utc"
        )
        if replay_created < fresh_verification_created:
            raise GateError("Runner-only replay receipt predates fresh verification")
        receipt_payload = {
            "schema_name": "intersubmod.m2v5_runner_only_gate_replay",
            "schema_version": "1.0.0",
            "created_at_utc": replay_created_at,
            "task_type": "B_comprehensive_validation",
            "status": "RUNNER_PREFIX_REPLAY_EVIDENCE_PASS_NON_AUTHORITATIVE",
            "scope": "completion_runner_physical_lines_1_358_only",
            "canonical_command": expected_command(),
            "clean_environment": EXPECTED_ENVIRONMENT,
            "code": {
                "runner_gate_replay": self_record,
                "continuation_verifier": verifier_record,
                "downstream_continuation": continuation_record,
                "release_source_authority_validator": authority_source_record,
            },
            "runtime": {
                "python": python_record,
                "bash": bash_record,
                "openssl": openssl_record,
            },
            "promotion_verifier": {
                "authorized_command": verifier_command(),
                "actual_argv": fd_bound_verifier_command(python_fd, verifier_fd),
                "normalized_argv": normalized_fd_bound_verifier_command(),
                "pass_fds": [python_fd, verifier_fd],
                "fd_artifact_binding": {
                    "python": {"fd": python_fd, "artifact": python_record},
                    "verifier": {"fd": verifier_fd, "artifact": verifier_record},
                },
                "source_descriptor_bound": True,
                "exit_code": verifier_result.returncode,
                "stdout": stream_evidence(verifier_result.stdout),
                "stderr": stream_evidence(verifier_result.stderr),
                "verification_receipt": verification_record,
                "fresh_verification_mode": verifier_stdout["mode"],
                "recorded_verification_mode": verification_payload["mode"],
                "fresh_and_recorded_evidence_agree": True,
                "pass": True,
            },
            "promotion_trust_chain": {
                "authorization": authorization_full_record,
                "prepare_exit_attestation": prepare_exit_full_record,
                "authorization_signature": authorization_signature_full_record,
                "authorization_signature_target": str(PREPARE_EXIT_ATTESTATION),
                "authorization_public_key": authorization_public_record,
                "retired_authorization_private_key": authorization_private_record,
                "promotion_receipt": promotion_full_record,
                "promote_exit_attestation": promote_exit_full_record,
                "promotion_signature": promotion_signature_full_record,
                "promotion_signature_target": str(PROMOTE_EXIT_ATTESTATION),
                "completion_public_key": completion_public_record,
                "retired_completion_private_key": completion_private_record,
                "continuation_public_key": continuation_public_record,
                "pre_signature_continuation_private_key": continuation_private_record,
                "historical_source_receipt": source_record,
                "canonical_source_receipt": canonical_record,
                "reviewed_sources": expected_reviewed_sources,
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
            },
            "v7_validation": {
                "authority_id": v7_validation["authority_id"],
                "authority": v7_record,
                "approval": v7_approval_record,
                "signature": v7_signature_record,
                "public_key": v7_public_record,
                "retired_private_key": v7_private_record,
                "source_set_sha256": v7_validation["source_set_sha256"],
                "required_roles_verified": v7_validation["required_roles_verified"],
                "all_source_roles_verified": v7_validation["all_source_roles_verified"],
                "process_lifetime_binding": v7_validation["validation_binding_lease"],
                "pass": True,
            },
            "completion_runner": {
                "artifact": runner_record,
                "observed_total_physical_lines": len(runner_lines),
                "included_line_start": RUNNER_INCLUDED_LINE_START,
                "included_line_end": RUNNER_INCLUDED_LINE_END,
                "first_excluded_line": RUNNER_FIRST_EXCLUDED_LINE,
                "first_excluded_line_is_blank": True,
                "first_downstream_executable_line": 360,
                "first_downstream_executable_line_prefix": "for path in \\",
                "stdin_size_bytes": len(runner_prefix),
                "stdin_sha256": runner_prefix_sha256,
                "authorization_bound_runtime_contract": runtime_contract,
            },
            "runner_only_replay": {
                "command": ["/proc/self/fd/<bound-bash-fd>", "-s"],
                "stdin": {
                    "source": str(COMPLETION_RUNNER),
                    "physical_lines": "1-358",
                    "size_bytes": len(runner_prefix),
                    "sha256": runner_prefix_sha256,
                },
                "exit_code": replay_result.returncode,
                "stdout": stream_evidence(replay_result.stdout),
                "stderr": stream_evidence(replay_result.stderr),
                "shell_observation_checks": shell_observation_checks,
                "authoritative_gate_validation": False,
                "pass": True,
            },
            "output_slot_checks": {
                "before_verifier": downstream_before,
                "after_verifier": downstream_after_verifier,
                "after_runner_replay": downstream_after_replay,
                "all_absent": True,
            },
            "log": log_record,
            "claims": {
                "original_completion_runner_end_to_end_pass_claimed": False,
                "original_tumor_ref_verifier_successfully_reexecuted_claimed": False,
                "promotion_integrity_verified": True,
                "runner_lines_1_358_replayed_successfully": True,
                "shell_pathname_gate_result_is_authoritative": False,
                "runner_line_359_or_later_executed": False,
                "strict_executed": False,
                "matched_normal_executed": False,
                "cn_ccf_executed": False,
                "final_dataset_or_report_executed": False,
                "downstream_authorized_by_this_receipt": False,
            },
            "downstream_not_executed": True,
            "downstream_authorized_after_this_gate": False,
            "pass": True,
            "pass_semantics": (
                "non_authoritative_pathname_replay_evidence_only; signed downstream "
                "continuation must independently revalidate all gate inputs from bound "
                "descriptors before any downstream command"
            ),
        }
        receipt_data = encode_json(receipt_payload)
        require_absent(PROMOTION_INCIDENT, "promotion failure incident before replay commit")
        require_absent(
            LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
            "legacy authorization-payload signature before replay commit",
        )
        require_absent(
            LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
            "legacy completion-payload signature before replay commit",
        )
        absence_snapshot(DOWNSTREAM_OUTPUT_SLOTS, "pre-replay-receipt commit downstream")
        reader.require_paths_still_bound()
        reader.retain_until_process_exit()

        def receipt_precommit() -> None:
            require_absent(
                PROMOTION_INCIDENT, "promotion failure incident during replay staging"
            )
            require_absent(
                LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
                "legacy authorization-payload signature during replay staging",
            )
            require_absent(
                LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
                "legacy completion-payload signature during replay staging",
            )
            if (
                absence_snapshot(
                    DOWNSTREAM_OUTPUT_SLOTS, "terminal replay-receipt commit downstream"
                )
                != downstream_after_replay
            ):
                raise GateError("Downstream absence contract changed during receipt staging")
            reader.require_paths_still_bound()

        write_exclusive_readonly(
            OUTPUT_RECEIPT,
            receipt_data,
            precommit=receipt_precommit,
            terminal_commit=lambda: os._exit(0),
        )
    except Exception:
        reader.close()
        raise


if __name__ == "__main__":
    raise SystemExit(main())
