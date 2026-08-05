#!/usr/bin/env python3
"""Authorize and publish the historical tumor-REF source receipt exactly once."""

from __future__ import annotations

import ast
import ctypes
import hashlib
import json
import os
from pathlib import Path
import stat
import sys
import types
from datetime import datetime, timezone
from typing import Any, Callable, Mapping


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
PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")
PRIMARY_PYTHON_RUNTIME = Path(
    "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python3.11"
)
QA_PYTHON_RUNTIME = Path("/bip7_disk/liaoyoyo2001/miniconda3/bin/python3.9")
OPENSSL = Path("/usr/bin/openssl")
PYTHON_CACHE_ROOTS = {
    "--prepare-authorization": WORKSPACE_ROOT / ".python_cache_tumor_ref_promotion_v2_prepare",
    "--promote": WORKSPACE_ROOT / ".python_cache_tumor_ref_promotion_v2_promote",
}

PROMOTION_TOOL = AUDIT_ROOT / "promote_tumor_ref_source_receipt_v2.py"
CONTINUATION_VERIFIER = AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_v2.py"
RUNNER_GATE_REPLAY = AUDIT_ROOT / "replay_m2v5_runner_only_gates_v1.py"
DOWNSTREAM_CONTINUATION = (
    AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_v1.py"
)

SOURCE_RECEIPT = (
    RESULT_ROOT / "tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json"
)
CANONICAL_RECEIPT = (
    WORKSPACE_ROOT
    / "tumor_ref_recovery_source_identity_v1"
    / "post_run_source_identity.receipt.json"
)
AUTHORIZATION = RESULT_ROOT / "tumor_ref_source_receipt_promotion_authorization.v3.json"
PREPARE_EXIT_ATTESTATION = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_prepare_exit_attestation.v1.json"
)
AUTHORIZATION_SIGNATURE = Path(str(PREPARE_EXIT_ATTESTATION) + ".ed25519.sig")
LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE = Path(str(AUTHORIZATION) + ".ed25519.sig")
COMPLETION_RECEIPT = RESULT_ROOT / "tumor_ref_source_receipt_promotion.v3.json"
PROMOTE_EXIT_ATTESTATION = (
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_promote_exit_attestation.v1.json"
)
COMPLETION_SIGNATURE = Path(str(PROMOTE_EXIT_ATTESTATION) + ".ed25519.sig")
LEGACY_COMPLETION_PAYLOAD_SIGNATURE = Path(str(COMPLETION_RECEIPT) + ".ed25519.sig")
PROMOTION_INCIDENT = RESULT_ROOT / "tumor_ref_source_receipt_promotion_incident.v3.json"
VERIFICATION_RECEIPT = RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.v3.json"
RUNNER_GATE_RECEIPT = RESULT_ROOT / "m2v5_runner_only_gate_replay.v1.json"
RUNNER_GATE_LOG = RESULT_ROOT / "m2v5_runner_only_gate_replay.v1.log"

PREPARE_EXIT_ATTESTATION_SCHEMA = {
    "name": (
        "intersubmod.tumor_ref_source_receipt_promotion."
        "prepare_exit_attestation"
    ),
    "version": "1.0.0",
}
PROMOTE_EXIT_ATTESTATION_SCHEMA = {
    "name": (
        "intersubmod.tumor_ref_source_receipt_promotion."
        "promote_exit_attestation"
    ),
    "version": "1.0.0",
}
EXIT_ATTESTATION_KEYS = frozenset(
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

AUTHORIZATION_KEY_DIR = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_promotion_authorization/"
    "20260722_tumor_ref_receipt_promotion_v4"
)
AUTHORIZATION_PUBLIC_KEY = AUTHORIZATION_KEY_DIR / "ed25519_public.pem"
AUTHORIZATION_PRIVATE_KEY = (
    AUTHORIZATION_KEY_DIR / "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
AUTHORIZATION_PUBLIC_KEY_SHA256 = (
    "e638b29a1d9c207dfa9849ff20a3867dacf702f7c1a899f1968ea1401a839677"
)
COMPLETION_KEY_DIR = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_promotion_completion/"
    "20260722_tumor_ref_receipt_promotion_v4"
)
COMPLETION_PUBLIC_KEY = COMPLETION_KEY_DIR / "ed25519_public.pem"
COMPLETION_PRIVATE_KEY = (
    COMPLETION_KEY_DIR / "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
COMPLETION_PUBLIC_KEY_SHA256 = (
    "6c55230f368db6b1fb0eedbc6c2362c0df7ee41a2699764fadb1f424e7036032"
)
CONTINUATION_KEY_DIR = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260722_m2v5_terminal_v2"
)
CONTINUATION_PUBLIC_KEY = CONTINUATION_KEY_DIR / "ed25519_public.pem"
CONTINUATION_PRIVATE_KEY = (
    CONTINUATION_KEY_DIR / "ed25519_private_one_time_resident.pem"
)
CONTINUATION_PUBLIC_KEY_SHA256 = (
    "a98370415594ebc9e8eb166bb70cab3550067b8a83db7d73ad23fac0891cb8b5"
)
PROTECTED_PRIVATE_KEY_PATHS = frozenset(
    {
        AUTHORIZATION_PRIVATE_KEY,
        COMPLETION_PRIVATE_KEY,
        CONTINUATION_PRIVATE_KEY,
    }
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

SNAPSHOT = (
    WORKSPACE_ROOT
    / "tumor_ref_recovery_source_identity_v1"
    / "observed_during_execution.snapshot.json"
)
TUMOR_REF_MANIFEST = (
    WORKSPACE_ROOT
    / "all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel"
    / "run_manifest.json"
)
FOCAL_ALT_CLUSTER_LIB = (
    REPO_ROOT
    / "research"
    / "20260715_single_fp_alt_multicluster_subclone_validation"
    / "scripts"
    / "focal_alt_cluster_lib.py"
)
ARCHIVE_RECEIPT = RESULT_ROOT / "completion_attempt2_archive_receipt.v1.json"
PRIOR_REVIEW = (
    TOPIC_ROOT
    / "reviews"
    / "20260717_tumor_REF來源attestation與release接縫獨立稽核_01.md"
)
V3_AUTHORITY = (
    REPO_ROOT
    / "docs"
    / "provenance"
    / "source_authorities"
    / "20260717_all_ssnv_focal_alt_release_source_authority.v3.json"
)
V3_APPROVAL = V3_AUTHORITY.with_name(V3_AUTHORITY.name.removesuffix(".json") + ".approval.json")
V3_SIGNATURE = Path(str(V3_APPROVAL) + ".ed25519.sig")
V3_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
    "20260717_all_ssnv_v3/ed25519_public.pem"
)
V3_PUBLIC_KEY_SHA256 = "60ebac3ee2ebfbf69a80331b40410b365c0315d41363d59ee0b44f6dbf5040e4"
V3_PRIVATE_KEY = V3_PUBLIC_KEY.with_name(
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
V7_AUTHORITY = (
    REPO_ROOT
    / "docs"
    / "provenance"
    / "source_authorities"
    / "20260722_all_ssnv_focal_alt_release_source_authority.v7.json"
)
V7_APPROVAL = V7_AUTHORITY.with_name(V7_AUTHORITY.name.removesuffix(".json") + ".approval.json")
V7_SIGNATURE = Path(str(V7_APPROVAL) + ".ed25519.sig")
V7_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
    "20260722_all_ssnv_v10_strict_command_parity_bootstrap/ed25519_public.pem"
)
V7_PUBLIC_KEY_SHA256 = "cecb2287cf87f8bd948c7390bd5f7059578966e0fc917de420572ddd82d8f21b"
V7_PRIVATE_KEY = V7_PUBLIC_KEY.with_name(
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
SOURCE_AUTHORITY_MODULE = SCRIPT_ROOT / "release_source_authority.py"
FINAL_DATASET_BUILDER = SCRIPT_ROOT / "build_all_ssnv_final_report_dataset.py"
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

REVIEW_SLOTS = {
    "mendel": {
        "path": TOPIC_ROOT / "reviews" / "20260722_tumor_ref_receipt_promotion_mendel.v3.json",
        "reviewer": "Mendel",
        "reviewer_agent_id": "mendel-v3-distinct-review-artifact",
    },
    "nash": {
        "path": TOPIC_ROOT / "reviews" / "20260722_tumor_ref_receipt_promotion_nash.v3.json",
        "reviewer": "Nash",
        "reviewer_agent_id": "nash-v3-distinct-review-artifact",
    },
    "external_claude_opus": {
        "path": TOPIC_ROOT
        / "reviews"
        / "20260722_tumor_ref_receipt_promotion_external_claude_opus.v3.json",
        "reviewer": "External Claude Opus",
        "reviewer_agent_id": "external-claude-opus-v3-distinct-review-artifact",
    },
}

EXPECTED_SOURCE_SHA256 = "d9fd6ecfc92b1e76ef6b68932faaab1e3d8e8540f521e5a7065f760bdd228e19"
EXPECTED_SOURCE_SIZE = 6733
EXPECTED_SNAPSHOT_SHA256 = "1c960fbe8fc512b8288ae3e67519043d35fb26173abbf63fee79624673ef52bc"
EXPECTED_MANIFEST_SHA256 = "7aec5c8a94a9bf95d7b21af9515dfd5add7b4a214b5dfaf48d4bb98315d1e615"
EXPECTED_FOCAL_SHA256 = "ed5f3c99461248248b20a9f49597ab5de7340a1e2055d77a1d83dbcc2799b72a"
EXPECTED_FOCAL_SIZE = 20253
EXPECTED_ARCHIVE_RECEIPT_SHA256 = "ad33b78379e387b68dab948cd60644e613c0387884fd4d6e35748fd1bb82dbed"
EXPECTED_PRIOR_REVIEW_SHA256 = "c4ca5d9538127af67f79545ca7f84a4f227a167e4bec99eb3bcff7acc6de9aa4"
EXPECTED_AUTHORITY_MODULE_SHA256 = "36621dc735d17cf3322e34aeb5d72247869935e7b87306cdb600ff6d522c52af"
EXPECTED_PYTHON_SHA256 = "777797a57eb75b28f530191628d26b14afada9ced2cb51c0ecae1eb62796062e"
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
EXPECTED_V3_ID = "20260717_all_ssnv_focal_alt_task_b_release_v3"
EXPECTED_V7_ID = "20260722_all_ssnv_focal_alt_task_b_release_v7_strict_command_parity"

AUTHORIZATION_ID = "20260722_tumor_ref_receipt_promotion_v3"
REVIEWED_SOURCE_ROLES = {
    "promotion_tool": PROMOTION_TOOL,
    "continuation_verifier": CONTINUATION_VERIFIER,
    "runner_gate_replay": RUNNER_GATE_REPLAY,
    "downstream_continuation": DOWNSTREAM_CONTINUATION,
    "primary_python_runtime": PRIMARY_PYTHON_RUNTIME,
    "qa_python_runtime": QA_PYTHON_RUNTIME,
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
REVIEWED_SOURCE_RUNTIME_CONTRACTS = {
    "primary_python_runtime": (EXPECTED_PYTHON_SHA256, "0o775"),
    "qa_python_runtime": (EXPECTED_QA_PYTHON_SHA256, "0o755"),
    "portable_builder_module": (EXPECTED_PORTABLE_BUILDER_SHA256, "0o644"),
    "portable_chart_extractor_module": (
        EXPECTED_PORTABLE_CHART_EXTRACTOR_SHA256,
        "0o644",
    ),
    "portable_verifier_module": (EXPECTED_PORTABLE_VERIFIER_SHA256, "0o644"),
    "portable_browser_helpers_module": (
        EXPECTED_PORTABLE_BROWSER_HELPERS_SHA256,
        "0o644",
    ),
    "portable_browser_cli_module": (EXPECTED_PORTABLE_BROWSER_CLI_SHA256, "0o644"),
    "portable_server_bundle": (EXPECTED_PORTABLE_SERVER_BUNDLE_SHA256, "0o644"),
    **{
        f"portable_reader_asset_part{index:03d}": (
            EXPECTED_PORTABLE_READER_ASSET_SHA256[index - 1],
            "0o644",
        )
        for index in range(1, 4)
    },
}
TRUSTED_KEY_ROLES = {
    "authorization_public_key": AUTHORIZATION_PUBLIC_KEY,
    "completion_public_key": COMPLETION_PUBLIC_KEY,
    "continuation_public_key": CONTINUATION_PUBLIC_KEY,
}

ARCHIVE_DIR = (
    WORKSPACE_ROOT
    / "audit_archives"
    / "20260722_completion_attempt2_mode_tightening_conflict"
)
ARCHIVED_LOG = ARCHIVE_DIR / "m2v5_recovered_completion_chain_v2_source_authority_v5.failed.log"
ARCHIVED_CACHE = ARCHIVE_DIR / "python_cache_m2v5_completion_v2_bound_bootstrap.failed"

DOWNSTREAM_OUTPUTS = (
    WORKSPACE_ROOT / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5",
    WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5",
    RESULT_ROOT / "stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json",
    RESULT_ROOT / "frozen_input_immutability.post_m2v5_downstream_v3_source_authority_v5.json",
    WORKSPACE_ROOT / "all_ssnv_final_report_dataset_v5_m2v5_source_attested",
    WORKSPACE_ROOT / "all_ssnv_final_report_v5_m2v5_source_attested",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json",
    Path(str(RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json") + ".ed25519.sig"),
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json",
    Path(str(RESULT_ROOT / "task_b_final_report_release_receipt.v1.json") + ".ed25519.sig"),
    WORKSPACE_ROOT / "m2v5_recovered_completion_chain_v2_source_authority_v5.log",
    WORKSPACE_ROOT / ".python_cache_m2v5_completion_v2_bound_bootstrap",
    RESULT_ROOT / "m2v5_downstream_continuation.v1.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.v1.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.v1.json.ed25519.sig",
)

CLEAN_ENV = {
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
FORBIDDEN_ENV = {
    "BASH_ENV",
    "ENV",
    "LD_AUDIT",
    "LD_DEBUG",
    "LD_LIBRARY_PATH",
    "LD_PRELOAD",
    "OPENSSL_CONF",
    "OPENSSL_MODULES",
    "PYTHONHOME",
    "PYTHONINSPECT",
    "PYTHONPATH",
    "PYTHONSTARTUP",
}
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


class PromotionError(RuntimeError):
    """Raised when authorization or promotion cannot be proven exactly."""


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_json(data: bytes, label: str) -> dict[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        output: dict[str, Any] = {}
        for key, value in pairs:
            if key in output:
                raise PromotionError(f"Duplicate JSON key in {label}: {key}")
            output[key] = value
        return output

    def reject_nonfinite(value: str) -> None:
        raise PromotionError(f"Non-finite JSON value in {label}: {value}")

    try:
        payload = json.loads(
            data,
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_nonfinite,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PromotionError(f"Invalid JSON in {label}") from error
    if not isinstance(payload, dict):
        raise PromotionError(f"{label} is not a JSON object")
    return payload


def strict_json_equal(left: Any, right: Any) -> bool:
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(
            strict_json_equal(left[key], right[key]) for key in left
        )
    if isinstance(left, list):
        return len(left) == len(right) and all(
            strict_json_equal(left_item, right_item)
            for left_item, right_item in zip(left, right)
        )
    return left == right


def process_identity(pid: int) -> tuple[int, int]:
    if type(pid) is not int or pid <= 0:
        raise PromotionError("Process PID must be a positive JSON integer")
    try:
        raw = Path(f"/proc/{pid}/stat").read_bytes()
    except OSError as error:
        raise PromotionError(f"Unable to bind process identity for PID {pid}") from error
    separator = raw.rfind(b") ")
    if separator < 0:
        raise PromotionError(f"Malformed /proc stat record for PID {pid}")
    fields = raw[separator + 2 :].split()
    if len(fields) <= 19 or not fields[1].isdigit() or not fields[19].isdigit():
        raise PromotionError(f"Incomplete /proc stat identity for PID {pid}")
    parent_pid = int(fields[1])
    start_time_ticks = int(fields[19])
    if parent_pid < 0 or start_time_ticks <= 0:
        raise PromotionError(f"Invalid /proc stat identity for PID {pid}")
    return parent_pid, start_time_ticks


def build_producer_session(
    mode: str,
    nonce: bytes,
    *,
    parent_pid: int,
    parent_start_time_ticks: int,
    child_pid: int,
) -> dict[str, Any]:
    if mode not in PYTHON_CACHE_ROOTS:
        raise PromotionError(f"Unsupported producer mode: {mode}")
    if type(nonce) is not bytes or len(nonce) != 32:
        raise PromotionError("Producer session requires an exact 256-bit nonce")
    if (
        type(parent_pid) is not int
        or parent_pid <= 0
        or type(parent_start_time_ticks) is not int
        or parent_start_time_ticks <= 0
        or type(child_pid) is not int
        or child_pid <= 0
        or child_pid == parent_pid
    ):
        raise PromotionError("Producer process identity is invalid")
    _, observed_parent_start = process_identity(parent_pid)
    if observed_parent_start != parent_start_time_ticks:
        raise PromotionError("Producer parent PID/start-time binding changed")
    observed_child_parent, child_start_time_ticks = process_identity(child_pid)
    if observed_child_parent != parent_pid:
        raise PromotionError("Forked producer child is not bound to the live parent")
    session = {
        "protocol": "intersubmod.parent_fork_no_exec_child",
        "schema_version": "1.0.0",
        "mode": mode,
        "nonce_sha256": sha256_bytes(nonce),
        "parent_pid": parent_pid,
        "parent_start_time_ticks": parent_start_time_ticks,
        "child_pid": child_pid,
        "child_start_time_ticks": child_start_time_ticks,
        "fork_no_exec": True,
    }
    validate_producer_session(session, expected_mode=mode)
    return session


def validate_producer_session(
    value: Any, *, expected_mode: str | None = None
) -> dict[str, Any]:
    if type(value) is not dict or set(value) != PRODUCER_SESSION_KEYS:
        raise PromotionError("Producer session does not have the exact JSON contract")
    if (
        type(value["protocol"]) is not str
        or value["protocol"] != "intersubmod.parent_fork_no_exec_child"
        or type(value["schema_version"]) is not str
        or value["schema_version"] != "1.0.0"
        or type(value["mode"]) is not str
        or value["mode"] not in PYTHON_CACHE_ROOTS
        or (expected_mode is not None and value["mode"] != expected_mode)
        or type(value["nonce_sha256"]) is not str
        or len(value["nonce_sha256"]) != 64
        or any(character not in "0123456789abcdef" for character in value["nonce_sha256"])
        or type(value["fork_no_exec"]) is not bool
        or value["fork_no_exec"] is not True
    ):
        raise PromotionError("Producer session string or boolean field differs")
    for key in (
        "parent_pid",
        "parent_start_time_ticks",
        "child_pid",
        "child_start_time_ticks",
    ):
        if type(value[key]) is not int or value[key] <= 0:
            raise PromotionError(f"Producer session {key} is not a positive JSON integer")
    if value["parent_pid"] == value["child_pid"]:
        raise PromotionError("Producer parent and child PID must differ")
    return dict(value)


def decode_child_wait(waited_pid: int, raw_wait_status: int) -> dict[str, Any]:
    if type(waited_pid) is not int or waited_pid <= 0:
        raise PromotionError("waitpid returned an invalid child PID")
    if type(raw_wait_status) is not int or raw_wait_status < 0:
        raise PromotionError("waitpid returned an invalid raw status")
    exited = os.WIFEXITED(raw_wait_status)
    signaled = os.WIFSIGNALED(raw_wait_status)
    value = {
        "waited_pid": waited_pid,
        "raw_wait_status": raw_wait_status,
        "wifexited": bool(exited),
        "exit_status": os.WEXITSTATUS(raw_wait_status) if exited else None,
        "wifsignaled": bool(signaled),
        "terminating_signal": os.WTERMSIG(raw_wait_status) if signaled else None,
        "core_dumped": bool(os.WCOREDUMP(raw_wait_status)) if signaled else False,
    }
    validate_child_wait(value, expected_pid=waited_pid, require_success=False)
    return value


def validate_child_wait(
    value: Any,
    *,
    expected_pid: int | None = None,
    require_success: bool = True,
) -> dict[str, Any]:
    if type(value) is not dict or set(value) != CHILD_WAIT_KEYS:
        raise PromotionError("Child wait status does not have the exact JSON contract")
    if (
        type(value["waited_pid"]) is not int
        or value["waited_pid"] <= 0
        or type(value["raw_wait_status"]) is not int
        or value["raw_wait_status"] < 0
        or type(value["wifexited"]) is not bool
        or type(value["wifsignaled"]) is not bool
        or type(value["core_dumped"]) is not bool
        or (value["exit_status"] is not None and type(value["exit_status"]) is not int)
        or (
            value["terminating_signal"] is not None
            and type(value["terminating_signal"]) is not int
        )
        or (expected_pid is not None and value["waited_pid"] != expected_pid)
    ):
        raise PromotionError("Child wait field has an invalid strict JSON type or value")
    decoded = decode_child_wait_fields(value["waited_pid"], value["raw_wait_status"])
    if not strict_json_equal(value, decoded):
        raise PromotionError("Child wait decoded fields differ from the raw waitpid status")
    if require_success and not (
        value["raw_wait_status"] == 0
        and value["wifexited"] is True
        and value["exit_status"] == 0
        and value["wifsignaled"] is False
        and value["terminating_signal"] is None
        and value["core_dumped"] is False
    ):
        raise PromotionError("Producer child did not terminate by normal exit zero")
    return dict(value)


def decode_child_wait_fields(waited_pid: int, raw_wait_status: int) -> dict[str, Any]:
    exited = os.WIFEXITED(raw_wait_status)
    signaled = os.WIFSIGNALED(raw_wait_status)
    return {
        "waited_pid": waited_pid,
        "raw_wait_status": raw_wait_status,
        "wifexited": bool(exited),
        "exit_status": os.WEXITSTATUS(raw_wait_status) if exited else None,
        "wifsignaled": bool(signaled),
        "terminating_signal": os.WTERMSIG(raw_wait_status) if signaled else None,
        "core_dumped": bool(os.WCOREDUMP(raw_wait_status)) if signaled else False,
    }


def waitpid_exact(child_pid: int) -> dict[str, Any]:
    while True:
        try:
            waited_pid, raw_wait_status = os.waitpid(child_pid, 0)
            break
        except InterruptedError:
            continue
    if waited_pid != child_pid:
        raise PromotionError("waitpid reaped a process other than the producer child")
    return decode_child_wait(waited_pid, raw_wait_status)


def fork_no_exec_and_wait(
    mode: str, child_action: Callable[[dict[str, Any]], None]
) -> tuple[dict[str, Any], dict[str, Any]]:
    if mode not in PYTHON_CACHE_ROOTS or not callable(child_action):
        raise PromotionError("Invalid fork/no-exec producer request")
    nonce = os.urandom(32)
    parent_pid = os.getpid()
    observed_parent_pid, parent_start_time_ticks = process_identity(parent_pid)
    if observed_parent_pid < 0:
        raise PromotionError("Live producer parent identity is invalid")
    child_pid = os.fork()
    if child_pid == 0:
        try:
            child_session = build_producer_session(
                mode,
                nonce,
                parent_pid=parent_pid,
                parent_start_time_ticks=parent_start_time_ticks,
                child_pid=os.getpid(),
            )
            child_action(child_session)
        except BaseException as error:
            message = (
                f"fork/no-exec producer child failed: {type(error).__name__}: {error}\n"
            ).encode("utf-8", errors="replace")
            try:
                os.write(2, message[:8192])
            finally:
                os._exit(1)
        os._exit(0)
    try:
        parent_session = build_producer_session(
            mode,
            nonce,
            parent_pid=parent_pid,
            parent_start_time_ticks=parent_start_time_ticks,
            child_pid=child_pid,
        )
    except BaseException:
        waitpid_exact(child_pid)
        raise
    child_wait = waitpid_exact(child_pid)
    validate_producer_session(parent_session, expected_mode=mode)
    validate_child_wait(child_wait, expected_pid=child_pid, require_success=False)
    return parent_session, child_wait


def require_utc_timestamp(value: Any, label: str) -> datetime:
    if not isinstance(value, str):
        raise PromotionError(f"{label} is not a UTC timestamp string")
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise PromotionError(f"{label} is not ISO-8601") from error
    if (
        parsed.tzinfo is None
        or parsed.utcoffset() != timezone.utc.utcoffset(parsed)
        or parsed.microsecond != 0
        or value != parsed.isoformat(timespec="seconds")
    ):
        raise PromotionError(f"{label} is not canonical UTC")
    return parsed


def artifact_with_stat(fd: int, record: Mapping[str, Any]) -> dict[str, Any]:
    opened = os.fstat(fd)
    return {
        **dict(record),
        "device": opened.st_dev,
        "inode": opened.st_ino,
        "link_count": opened.st_nlink,
        "mtime_ns": opened.st_mtime_ns,
        "ctime_ns": opened.st_ctime_ns,
    }


def validate_full_stat_artifact(value: Any, label: str) -> dict[str, Any]:
    if type(value) is not dict or set(value) != FULL_STAT_ARTIFACT_KEYS:
        raise PromotionError(f"{label} does not have the exact full-stat artifact contract")
    if (
        type(value["path"]) is not str
        or not value["path"].startswith("/")
        or type(value["sha256"]) is not str
        or len(value["sha256"]) != 64
        or any(character not in "0123456789abcdef" for character in value["sha256"])
        or type(value["mode"]) is not str
        or not value["mode"].startswith("0o")
    ):
        raise PromotionError(f"{label} has an invalid string identity field")
    for key in (
        "size_bytes",
        "device",
        "inode",
        "link_count",
        "mtime_ns",
        "ctime_ns",
    ):
        if type(value[key]) is not int or value[key] < 0:
            raise PromotionError(f"{label} {key} is not a non-negative JSON integer")
    if value["inode"] <= 0 or value["link_count"] != 1:
        raise PromotionError(f"{label} is not a single-link regular artifact")
    return dict(value)


def reopen_full_stat_artifact(
    reader: Any,
    path: Path,
    *,
    expected_basic: Mapping[str, Any] | None = None,
    expected_sha256: str | None = None,
    expected_size: int | None = None,
    expected_mode: str = "0o444",
) -> tuple[int, bytes, dict[str, Any]]:
    descriptor, data, record = open_artifact(
        reader,
        path,
        expected_sha256=expected_sha256,
        expected_size=expected_size,
        expected_mode=expected_mode,
        with_stat=True,
    )
    validate_full_stat_artifact(record, str(path))
    if expected_basic is not None:
        basic = {
            key: record[key] for key in ("path", "size_bytes", "sha256", "mode")
        }
        if not strict_json_equal(basic, dict(expected_basic)):
            raise PromotionError(f"Full-stat artifact differs from the prior binding: {path}")
    return descriptor, data, record


def exit_attestation_announcement(
    path: Path,
    schema: Mapping[str, str],
    public_key: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "path": str(path),
        "schema": dict(schema),
        "signing_key": dict(public_key),
    }


def exit_attestation_signature_contract(
    *,
    path: Path,
    schema: Mapping[str, str],
    signature: Path,
    public_key: Mapping[str, Any],
    private_key_pre_signature: Mapping[str, Any],
    legacy_payload_signature: Path,
) -> dict[str, Any]:
    return {
        "attestation_schema": dict(schema),
        "algorithm": "Ed25519",
        "signed_artifact": str(path),
        "signature": str(signature),
        "public_key": dict(public_key),
        "private_key_lifecycle": {
            "pre_signature": dict(private_key_pre_signature),
            "required_post_signature_mode": "0o0",
            "signing_must_follow_attestation_publication": True,
        },
        "legacy_payload_signature": {
            "path": str(legacy_payload_signature),
            "required_absent": True,
        },
    }


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


def source_code_inventory(namespace: Mapping[str, Any]) -> dict[str, str]:

    inventory: dict[str, str] = {}
    for name, value in namespace.items():
        if isinstance(value, types.FunctionType) and value.__code__.co_filename == str(PROMOTION_TOOL):
            inventory[name] = stable_code_digest(value.__code__)
        elif isinstance(value, type) and value.__module__ == namespace.get("__name__"):
            for method_name, method in vars(value).items():
                function = method.__func__ if isinstance(method, (staticmethod, classmethod)) else method
                if isinstance(function, types.FunctionType):
                    inventory[f"{name}.{method_name}"] = stable_code_digest(function.__code__)
    return inventory


def executing_module_code() -> types.CodeType:
    frame = sys._getframe()
    while frame is not None:
        if frame.f_code.co_name == "<module>" and frame.f_code.co_filename == str(PROMOTION_TOOL):
            return frame.f_code
        frame = frame.f_back
    raise PromotionError("Executing promotion module frame is unavailable")


def verify_executing_source_matches_bound_bytes(data: bytes) -> dict[str, Any]:
    shadow_name = "promotion_v2_bound_source_shadow"
    shadow: dict[str, Any] = {
        "__name__": shadow_name,
        "__file__": str(PROMOTION_TOOL),
        "__package__": "",
    }
    reviewed_module_code = compile(data, str(PROMOTION_TOOL), "exec")
    exec(reviewed_module_code, shadow)
    module_code_sha256 = stable_code_digest(executing_module_code())
    reviewed_module_code_sha256 = stable_code_digest(reviewed_module_code)
    live_inventory = source_code_inventory(globals())
    shadow_inventory = source_code_inventory(shadow)
    critical_names = (
        "EXPECTED_SOURCE_SHA256",
        "EXPECTED_SOURCE_SIZE",
        "EXPECTED_SNAPSHOT_SHA256",
        "EXPECTED_MANIFEST_SHA256",
        "EXPECTED_FOCAL_SHA256",
        "EXPECTED_ARCHIVE_RECEIPT_SHA256",
        "EXPECTED_PRIOR_REVIEW_SHA256",
        "EXPECTED_AUTHORITY_MODULE_SHA256",
        "EXPECTED_PYTHON_SHA256",
        "EXPECTED_QA_PYTHON_SHA256",
        "EXPECTED_PORTABLE_BUILDER_SHA256",
        "EXPECTED_PORTABLE_CHART_EXTRACTOR_SHA256",
        "EXPECTED_PORTABLE_VERIFIER_SHA256",
        "EXPECTED_PORTABLE_BROWSER_HELPERS_SHA256",
        "EXPECTED_PORTABLE_BROWSER_CLI_SHA256",
        "EXPECTED_PORTABLE_SERVER_BUNDLE_SHA256",
        "EXPECTED_PORTABLE_READER_ASSET_SHA256",
        "AUTHORIZATION_PUBLIC_KEY_SHA256",
        "COMPLETION_PUBLIC_KEY_SHA256",
        "CONTINUATION_PUBLIC_KEY_SHA256",
        "AUTHORIZATION_PUBLIC_KEY",
        "AUTHORIZATION_PRIVATE_KEY",
        "COMPLETION_PUBLIC_KEY",
        "COMPLETION_PRIVATE_KEY",
        "CONTINUATION_PUBLIC_KEY",
        "AUTHORIZATION",
        "AUTHORIZATION_SIGNATURE",
        "COMPLETION_RECEIPT",
        "COMPLETION_SIGNATURE",
        "CONTINUATION_PRIVATE_KEY",
        "PROTECTED_PRIVATE_KEY_PATHS",
        "PREPARE_EXIT_ATTESTATION",
        "PROMOTE_EXIT_ATTESTATION",
        "PREPARE_EXIT_ATTESTATION_SCHEMA",
        "PROMOTE_EXIT_ATTESTATION_SCHEMA",
        "EXIT_ATTESTATION_KEYS",
        "PRODUCER_SESSION_KEYS",
        "CHILD_WAIT_KEYS",
        "FULL_STAT_ARTIFACT_KEYS",
        "PREPARE_EXIT_CHECKS",
        "PROMOTE_EXIT_CHECKS",
        "PREPARE_EXIT_PASS_SEMANTICS",
        "PROMOTE_EXIT_PASS_SEMANTICS",
        "LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE",
        "LEGACY_COMPLETION_PAYLOAD_SIGNATURE",
        "CONTINUATION_TERMINAL_RECEIPT",
        "CONTINUATION_EXIT_ATTESTATION",
        "CONTINUATION_TERMINAL_SIGNATURE",
        "SUPERVISOR_CAPABILITY_PROTOCOL",
        "CONTINUATION_THREAT_MODEL",
        "CONTINUATION_POLICY",
        "V3_PUBLIC_KEY_SHA256",
        "V7_PUBLIC_KEY_SHA256",
        "EXPECTED_V3_ID",
        "EXPECTED_V7_ID",
        "AUTHORIZATION_ID",
        "REVIEWED_SOURCE_ROLES",
        "REVIEWED_SOURCE_RUNTIME_CONTRACTS",
        "TRUSTED_KEY_ROLES",
        "REVIEW_SLOTS",
        "DOWNSTREAM_OUTPUTS",
        "DOWNSTREAM_CONTINUATION",
        "AUTHORIZATION_CHECKS",
        "COMPLETION_CHECKS",
        "CLEAN_ENV",
        "FORBIDDEN_ENV",
    )
    constant_match = all(
        type(globals()[name]) is type(shadow[name])
        and globals()[name] == shadow[name]
        for name in critical_names
    )
    command_match = all(
        expected_command(mode) == shadow["expected_command"](mode)
        for mode in PYTHON_CACHE_ROOTS
    )
    binding_checks = {
        "live_inventory_nonempty": bool(live_inventory),
        "module_code_equal": module_code_sha256 == reviewed_module_code_sha256,
        "inventory_equal": live_inventory == shadow_inventory,
        "critical_constants_equal": constant_match,
        "canonical_commands_equal": command_match,
    }
    if not all(binding_checks.values()):
        raise PromotionError(
            "Executing promotion definitions differ from bound reviewed source bytes: "
            + json.dumps(binding_checks, sort_keys=True)
        )
    return {
        "source_sha256": sha256_bytes(data),
        "function_and_method_count": len(live_inventory),
        "inventory_sha256": sha256_bytes(
            json.dumps(live_inventory, sort_keys=True).encode("utf-8")
        ),
        "module_code_sha256": module_code_sha256,
        "critical_constants_equal": True,
        "canonical_commands_equal": True,
        "pass": True,
    }


def expected_command(mode: str) -> list[str]:
    return [
        str(PYTHON),
        "-I",
        "-B",
        "-X",
        f"pycache_prefix={PYTHON_CACHE_ROOTS[mode]}",
        str(PROMOTION_TOOL),
        mode,
    ]


def observed_command() -> list[str]:
    raw = Path("/proc/self/cmdline").read_bytes()
    if not raw.endswith(b"\0"):
        raise PromotionError("Unable to bind process command")
    return [os.fsdecode(token) for token in raw[:-1].split(b"\0")]


def require_clean_invocation(mode: str) -> None:
    if observed_command() != expected_command(mode):
        raise PromotionError("Promotion tool was not invoked by the canonical direct command")
    for key, value in CLEAN_ENV.items():
        if os.environ.get(key) != value:
            raise PromotionError(f"Clean-environment value differs: {key}")
    unexpected = set(os.environ) - set(CLEAN_ENV)
    if unexpected:
        raise PromotionError(f"Unexpected environment variables: {sorted(unexpected)}")
    present_forbidden = FORBIDDEN_ENV.intersection(os.environ)
    if present_forbidden:
        raise PromotionError(f"Forbidden environment variables are present: {sorted(present_forbidden)}")
    cache_root = PYTHON_CACHE_ROOTS[mode]
    if (
        cache_root.is_symlink()
        or not cache_root.is_dir()
        or stat.S_IMODE(cache_root.stat().st_mode) != 0o700
        or any(cache_root.iterdir())
    ):
        raise PromotionError("Python cache prefix must be a fresh empty mode-0700 directory")


def require_active_cache_empty() -> None:
    cache_root = PYTHON_CACHE_ROOTS[sys.argv[1]]
    if any(cache_root.iterdir()):
        raise PromotionError("Fresh Python cache prefix was populated during execution")


def bootstrap_authority_module() -> tuple[int, os.stat_result, Any]:
    if SOURCE_AUTHORITY_MODULE.is_symlink():
        raise PromotionError("Source-authority validator must not be a symlink")
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(SOURCE_AUTHORITY_MODULE, flags)
    opened = os.fstat(descriptor)
    data = os.pread(descriptor, opened.st_size, 0)
    live = SOURCE_AUTHORITY_MODULE.stat()
    if (
        not stat.S_ISREG(opened.st_mode)
        or stat.S_IMODE(opened.st_mode) != 0o444
        or opened.st_size != len(data)
        or sha256_bytes(data) != EXPECTED_AUTHORITY_MODULE_SHA256
        or (opened.st_dev, opened.st_ino, opened.st_mode, opened.st_size, opened.st_mtime_ns, opened.st_ctime_ns)
        != (live.st_dev, live.st_ino, live.st_mode, live.st_size, live.st_mtime_ns, live.st_ctime_ns)
    ):
        os.close(descriptor)
        raise PromotionError("Source-authority bootstrap identity differs from v7")
    module = types.ModuleType("release_source_authority")
    module.__file__ = str(SOURCE_AUTHORITY_MODULE)
    module.__package__ = ""
    exec(compile(data, str(SOURCE_AUTHORITY_MODULE), "exec"), module.__dict__)
    sys.modules["release_source_authority"] = module
    return descriptor, opened, module


def require_bootstrap_stable(descriptor: int, opened: os.stat_result) -> None:
    current = os.fstat(descriptor)
    live = SOURCE_AUTHORITY_MODULE.stat()
    fields = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns")
    if any(getattr(current, key) != getattr(opened, key) for key in fields) or any(
        getattr(live, key) != getattr(opened, key) for key in fields
    ):
        raise PromotionError("Source-authority bootstrap binding changed")


def link_fd_no_replace(source_fd: int, target_directory_fd: int, target_name: str) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    linkat = getattr(libc, "linkat", None)
    if linkat is None:
        raise PromotionError("linkat is required for atomic no-replace publication")
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
        raise PromotionError(
            "Atomic no-replace publication failed: " + os.strerror(error_number)
        )


def write_exclusive_readonly(
    path: Path,
    payload: Mapping[str, Any],
    *,
    precommit: Callable[[], None] | None = None,
) -> dict[str, Any]:
    if path in PROTECTED_PRIVATE_KEY_PATHS:
        raise PromotionError("Artifact publication cannot target a protected private key")
    encoded = (json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if not hasattr(os, "O_TMPFILE"):
        raise PromotionError("O_TMPFILE is required for staged publication")
    if path.parent.resolve(strict=True) != path.parent or stat.S_ISLNK(
        os.lstat(path.parent).st_mode
    ):
        raise PromotionError(f"Publication parent is not canonical: {path.parent}")
    parent_flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0)
    if hasattr(os, "O_NOFOLLOW"):
        parent_flags |= os.O_NOFOLLOW
    parent_fd = os.open(path.parent, parent_flags)
    parent_opened = os.fstat(parent_fd)
    descriptor = -1
    record: dict[str, Any] | None = None
    try:
        descriptor = os.open(
            ".", os.O_RDWR | os.O_TMPFILE | os.O_CLOEXEC, 0o400, dir_fd=parent_fd
        )
        offset = 0
        while offset < len(encoded):
            written = os.write(descriptor, encoded[offset:])
            if written <= 0:
                raise PromotionError(f"Short write while publishing {path}")
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
            or opened.st_size != len(encoded)
            or observed != encoded
        ):
            raise PromotionError(f"Staged artifact validation failed: {path}")
        record = {
            "path": str(path),
            "size_bytes": opened.st_size,
            "sha256": sha256_bytes(observed),
            "mode": oct(stat.S_IMODE(opened.st_mode)),
        }
        if precommit is not None:
            precommit()
        parent_live = os.stat(path.parent, follow_symlinks=False)
        if (
            parent_opened.st_dev,
            parent_opened.st_ino,
            parent_opened.st_mode,
        ) != (parent_live.st_dev, parent_live.st_ino, parent_live.st_mode):
            raise PromotionError(f"Publication parent path rebound: {path.parent}")
        require_path_absent(path, "Atomic publication target")
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
            committed_identity = (
                committed.st_dev,
                committed.st_ino,
                committed.st_mode,
                committed.st_nlink,
                committed.st_size,
                committed.st_mtime_ns,
                committed.st_ctime_ns,
            )
            if (
                committed.st_nlink != 1
                or committed_identity
                != (
                    linked.st_dev,
                    linked.st_ino,
                    linked.st_mode,
                    linked.st_nlink,
                    linked.st_size,
                    linked.st_mtime_ns,
                    linked.st_ctime_ns,
                )
                or committed_identity
                != (
                    reopened.st_dev,
                    reopened.st_ino,
                    reopened.st_mode,
                    reopened.st_nlink,
                    reopened.st_size,
                    reopened.st_mtime_ns,
                    reopened.st_ctime_ns,
                )
                or reopened_data != encoded
            ):
                raise PromotionError(f"Committed artifact path/inode drift: {path}")
        finally:
            os.close(reopened_fd)
        os.fsync(parent_fd)
        parent_after = os.stat(path.parent, follow_symlinks=False)
        if (
            parent_opened.st_dev,
            parent_opened.st_ino,
            parent_opened.st_mode,
        ) != (parent_after.st_dev, parent_after.st_ino, parent_after.st_mode):
            raise PromotionError(f"Publication parent changed after commit: {path.parent}")
        if path.resolve(strict=True) != path:
            raise PromotionError(f"Committed artifact path is not canonical: {path}")
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
            durable_identity = (
                durable_source.st_dev,
                durable_source.st_ino,
                durable_source.st_mode,
                durable_source.st_nlink,
                durable_source.st_size,
                durable_source.st_mtime_ns,
                durable_source.st_ctime_ns,
            )
            if (
                durable_source.st_nlink != 1
                or durable_identity
                != tuple(
                    getattr(durable_opened, field)
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
                or durable_identity
                != tuple(
                    getattr(durable_linked, field)
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
                or durable_data != encoded
            ):
                raise PromotionError(f"Post-fsync artifact path/inode drift: {path}")
        finally:
            os.close(durable_fd)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        os.close(parent_fd)
    if record is None or record["mode"] != "0o444":
        raise PromotionError(f"Staged artifact did not reach commit: {path}")
    return record


def require_path_absent(path: Path, label: str) -> None:
    if os.path.lexists(path):
        raise PromotionError(f"{label} output slot is occupied: {path}")


def require_downstream_absent() -> list[str]:
    for path in DOWNSTREAM_OUTPUTS:
        require_path_absent(path, "Downstream")
    return [str(path) for path in DOWNSTREAM_OUTPUTS]


def open_artifact(
    reader: Any,
    path: Path,
    *,
    expected_sha256: str | None = None,
    expected_size: int | None = None,
    expected_mode: str = "0o444",
    with_stat: bool = False,
) -> tuple[int, bytes, dict[str, Any]]:
    descriptor, data, record = reader.open(path, include_mode=True)
    if record["mode"] != expected_mode:
        raise PromotionError(f"Unexpected mode for {path}: {record['mode']}")
    if expected_sha256 is not None and record["sha256"] != expected_sha256:
        raise PromotionError(f"Unexpected SHA-256 for {path}")
    if expected_size is not None and record["size_bytes"] != expected_size:
        raise PromotionError(f"Unexpected size for {path}")
    if os.fstat(descriptor).st_nlink != 1:
        raise PromotionError(f"Trust artifact must have link count one: {path}")
    return descriptor, data, artifact_with_stat(descriptor, record) if with_stat else record


def open_private_key_metadata(
    reader: Any, path: Path, *, expected_mode: str
) -> tuple[int, dict[str, Any]]:
    descriptor, record = reader.open_metadata(path, include_mode=True)
    if record != {"path": str(path.resolve(strict=True)), "mode": expected_mode}:
        raise PromotionError(f"Private-key retirement state differs: {path}")
    return descriptor, record


def verify_signature(
    authority_module: Any,
    reader: Any,
    *,
    data_path: Path,
    signature_path: Path,
    public_key_path: Path,
    public_key_sha256: str,
) -> dict[str, Any]:
    openssl_fd, _, openssl_record = open_artifact(
        reader,
        OPENSSL,
        expected_sha256=authority_module.OPENSSL_SHA256,
        expected_mode="0o755",
    )
    data_fd, _, data_record = open_artifact(reader, data_path)
    public_fd, _, public_record = open_artifact(
        reader, public_key_path, expected_sha256=public_key_sha256
    )
    signature_fd, _, signature_record = open_artifact(reader, signature_path)
    if not authority_module.verify_ed25519_signature_fds(
        openssl_fd=openssl_fd,
        data_fd=data_fd,
        public_key_fd=public_fd,
        signature_fd=signature_fd,
    ):
        raise PromotionError(f"Detached Ed25519 signature failed: {signature_path}")
    return {
        "data": data_record,
        "signature": signature_record,
        "public_key": public_record,
        "openssl": openssl_record,
        "verified": True,
    }


def validate_signed_release_authority(
    authority_module: Any,
    reader: Any,
    *,
    authority_path: Path,
    approval_path: Path,
    signature_path: Path,
    expected_authority_id: str,
    expected_public_key_path: Path,
    expected_public_key_sha256: str,
    expected_private_key_path: Path,
) -> dict[str, Any]:
    authority_fd, authority_data, authority_record = open_artifact(reader, authority_path)
    approval_fd, approval_data, approval_record = open_artifact(reader, approval_path)
    signature_fd, _, signature_record = open_artifact(reader, signature_path)
    authority_payload = load_json(authority_data, authority_path.name)
    approval_payload = load_json(approval_data, approval_path.name)
    if (
        authority_payload.get("schema_name") != "intersubmod.release_source_authority"
        or authority_payload.get("schema_version") != "1.0.0"
        or authority_payload.get("authority_id") != expected_authority_id
        or authority_payload.get("approval_status") != "APPROVED_FOR_FULL_TASK_B_RUN"
        or approval_payload.get("schema_name") != "intersubmod.release_source_authority.approval"
        or approval_payload.get("schema_version") != "1.0.0"
        or approval_payload.get("authority_id") != expected_authority_id
        or approval_payload.get("approval_status") != "APPROVED_FOR_FULL_TASK_B_RUN"
        or approval_payload.get("authority_manifest") != authority_record
    ):
        raise PromotionError(f"Signed source-authority contract differs: {authority_path}")
    public_key_field = approval_payload.get("public_key")
    private_key_field = approval_payload.get("private_key_retirement")
    if not isinstance(public_key_field, Mapping) or not isinstance(private_key_field, Mapping):
        raise PromotionError(f"Source authority lacks key contract: {authority_path}")
    public_key_path = Path(str(public_key_field.get("path")))
    private_key_path = Path(str(private_key_field.get("path")))
    if (
        public_key_path != expected_public_key_path
        or private_key_path != expected_private_key_path
    ):
        raise PromotionError(f"Source-authority key path lacks independent pin: {authority_path}")
    public_fd, _, public_record = open_artifact(
        reader, public_key_path, expected_sha256=expected_public_key_sha256
    )
    private_fd, private_record = open_private_key_metadata(
        reader, private_key_path, expected_mode="0o0"
    )
    openssl_fd, _, openssl_record = open_artifact(
        reader,
        OPENSSL,
        expected_sha256=authority_module.OPENSSL_SHA256,
        expected_mode="0o755",
    )
    if public_key_field != public_record:
        raise PromotionError(f"Source-authority public key differs: {authority_path}")
    if (
        private_key_field.get("path") != private_record["path"]
        or private_key_field.get("required_post_signature_mode") != "0o0"
    ):
        raise PromotionError(f"Source-authority private-key contract differs: {authority_path}")
    if not authority_module.verify_ed25519_signature_fds(
        openssl_fd=openssl_fd,
        data_fd=approval_fd,
        public_key_fd=public_fd,
        signature_fd=signature_fd,
    ):
        raise PromotionError(f"Source-authority detached signature failed: {authority_path}")
    return {
        "authority": authority_record,
        "approval": approval_record,
        "signature": signature_record,
        "public_key": public_record,
        "private_key": private_record,
        "openssl": openssl_record,
        "authority_id": expected_authority_id,
        "signature_verified": True,
        "bound_descriptors": [authority_fd, approval_fd, signature_fd, public_fd, private_fd],
    }


def validate_review(
    reader: Any,
    *,
    role: str,
    slot: Mapping[str, Any],
    reviewed_sources: Mapping[str, Mapping[str, Any]],
    trusted_keys: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any], tuple[int, int]]:
    path = Path(slot["path"])
    descriptor, data, record = open_artifact(reader, path)
    payload = load_json(data, f"promotion review {role}")
    required_keys = {
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
    exact_contract = (
        set(payload) == required_keys
        and payload.get("schema_name") == "intersubmod.tumor_ref_receipt_promotion_review"
        and payload.get("schema_version") == "3.0.0"
        and payload.get("reviewer") == slot["reviewer"]
        and payload.get("reviewer_agent_id") == slot["reviewer_agent_id"]
        and payload.get("verdict") == "APPROVE"
        and payload.get("pass") is True
        and payload.get("findings_closed") is True
        and payload.get("conditions") == []
        and payload.get("unresolved_findings") == []
        and payload.get("unresolved_conditions") == []
        and payload.get("reviewed_sources") == reviewed_sources
        and payload.get("trusted_signing_keys") == trusted_keys
    )
    if not exact_contract:
        raise PromotionError(f"Promotion review does not satisfy exact v3 contract: {path}")
    opened = os.fstat(descriptor)
    return record, payload, (opened.st_dev, opened.st_ino)


def validate_archive_evidence(reader: Any) -> dict[str, Any]:
    _, archive_data, archive_record = open_artifact(
        reader, ARCHIVE_RECEIPT, expected_sha256=EXPECTED_ARCHIVE_RECEIPT_SHA256
    )
    payload = load_json(archive_data, "completion attempt-2 archive receipt")
    required_checks = {
        "failed_log_preserved",
        "failed_python_cache_preserved",
        "no_failed_artifact_deleted",
        "canonical_log_path_absent_after_archive",
        "canonical_cache_path_absent_after_archive",
        "strict_output_absent_after_failure",
        "matched_normal_outputs_absent_after_failure",
        "cn_ccf_output_absent_after_failure",
        "final_dataset_output_absent_after_failure",
        "final_report_output_absent_after_failure",
        "completion_runner_bytes_unchanged",
    }
    if (
        payload.get("schema_name") != "intersubmod.completion_attempt_archive_receipt"
        or payload.get("schema_version") != "1.0.0"
        or payload.get("status") != "FAILED_BEFORE_DOWNSTREAM_OUTPUT_AND_ARCHIVED"
        or payload.get("pass") is not True
        or set(payload.get("checks", {})) != required_checks
        or not all(payload["checks"].values())
    ):
        raise PromotionError("Attempt-2 archive receipt contract differs")
    _, _, log_record = open_artifact(
        reader,
        ARCHIVED_LOG,
        expected_sha256="870f4d1530491508d194ca49d36012651dfdd94f8a9ac072bfa712ac52b4e218",
        expected_size=3802,
        expected_mode="0o644",
    )
    if ARCHIVED_CACHE.is_symlink() or not ARCHIVED_CACHE.is_dir():
        raise PromotionError("Archived Python cache root is not a direct directory")
    files = sorted(path.resolve(strict=True) for path in ARCHIVED_CACHE.rglob("*") if path.is_file())
    if len(files) != 45:
        raise PromotionError("Archived Python-cache file count differs")
    cache_lines: list[str] = []
    total_bytes = 0
    for path in files:
        _, data, record = open_artifact(reader, path, expected_mode=oct(path.stat().st_mode & 0o777))
        total_bytes += len(data)
        cache_lines.append(f"{record['sha256']}  {path}\n")
    replayed_files = sorted(
        path.resolve(strict=True) for path in ARCHIVED_CACHE.rglob("*") if path.is_file()
    )
    if replayed_files != files:
        raise PromotionError("Archived Python-cache file set changed during replay")
    cache_manifest_sha256 = sha256_bytes("".join(cache_lines).encode("utf-8"))
    if (
        total_bytes != 1_463_182
        or cache_manifest_sha256
        != "0084670ad8f6440715a110c90aa399fc26573306c4e7a19ae40b9782d23264e2"
    ):
        raise PromotionError("Archived Python-cache byte manifest differs")
    return {
        "receipt": archive_record,
        "archived_log": log_record,
        "archived_cache": {
            "path": str(ARCHIVED_CACHE.resolve(strict=True)),
            "regular_file_count": len(files),
            "regular_file_total_bytes": total_bytes,
            "sha256sum_line_manifest_sha256": cache_manifest_sha256,
            "all_files_bound_until_process_exit": True,
        },
    }


BUILDER_GATE_CONSTANTS = {
    "SCRIPT_DIR",
    "REPO_ROOT",
    "TUMOR_REF_SOURCE_IDENTITY_VERIFIER",
    "SCHEMA_VERSION",
    "TUMOR_REF_SOURCE_IDENTITY_SCHEMA_VERSION",
    "TUMOR_REF_SOURCE_IDENTITY_REQUIRED_ROLES",
    "TUMOR_REF_SOURCE_IDENTITY_TRUE_CHECKS",
}
BUILDER_GATE_DEFINITIONS = {
    "ContractError",
    "sha256",
    "artifact",
    "source_script_token_matches_attested_source",
    "require_nonempty_file",
    "load_json",
    "require_schema_pass",
    "verify_declared_artifact",
    "verify_declared_path_artifact",
    "parse_utc_timestamp",
    "load_tumor_ref_source_identity_attestation",
}


def assigned_names(node: ast.Assign | ast.AnnAssign) -> set[str]:
    targets = node.targets if isinstance(node, ast.Assign) else [node.target]
    names: set[str] = set()
    for target in targets:
        if isinstance(target, ast.Name):
            names.add(target.id)
    return names


def load_builder_gate_from_bound_ast(data: bytes) -> Any:
    """Execute only the exact builder gate and its local helpers from bound v7 bytes."""
    tree = ast.parse(data, filename=str(FINAL_DATASET_BUILDER))
    selected: list[ast.stmt] = []
    observed: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module == "__future__":
            selected.append(node)
            continue
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            names = assigned_names(node)
            if names.intersection(BUILDER_GATE_CONSTANTS):
                selected.append(node)
                observed.update(names.intersection(BUILDER_GATE_CONSTANTS))
            continue
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if node.name in BUILDER_GATE_DEFINITIONS:
                selected.append(node)
                observed.add(node.name)
    expected = BUILDER_GATE_CONSTANTS | BUILDER_GATE_DEFINITIONS
    if observed != expected:
        raise PromotionError(
            f"Unable to extract exact final-builder source gate: missing={sorted(expected - observed)} "
            f"extra={sorted(observed - expected)}"
        )
    module_name = "build_all_ssnv_final_report_dataset_promotion_gate_bound"
    module = types.ModuleType(module_name)
    module.__file__ = str(FINAL_DATASET_BUILDER)
    module.__package__ = ""
    module.__dict__.update(
        {
            "Path": Path,
            "Any": Any,
            "Mapping": Mapping,
            "datetime": datetime,
            "timezone": timezone,
            "hashlib": hashlib,
            "json": json,
            "frozenset": frozenset,
        }
    )
    gate_tree = ast.Module(body=selected, type_ignores=[])
    ast.fix_missing_locations(gate_tree)
    exec(compile(gate_tree, str(FINAL_DATASET_BUILDER), "exec"), module.__dict__)
    if (
        module.load_tumor_ref_source_identity_attestation.__code__.co_firstlineno
        != 4729
    ):
        raise PromotionError("Extracted final-builder gate line identity differs")
    return module


def validate_core_evidence(authority_module: Any, reader: Any) -> dict[str, Any]:
    if PYTHON.resolve(strict=True) != PRIMARY_PYTHON_RUNTIME:
        raise PromotionError("Primary Python launcher target drift")
    python_fd, _, python_record = open_artifact(
        reader,
        PRIMARY_PYTHON_RUNTIME,
        expected_sha256=EXPECTED_PYTHON_SHA256,
        expected_mode="0o775",
    )
    process_executable = Path("/proc/self/exe").stat()
    python_opened = os.fstat(python_fd)
    if (process_executable.st_dev, process_executable.st_ino) != (
        python_opened.st_dev,
        python_opened.st_ino,
    ):
        raise PromotionError("Running interpreter differs from the reviewed canonical Python")
    source_fd, source_data, source_record = open_artifact(
        reader,
        SOURCE_RECEIPT,
        expected_sha256=EXPECTED_SOURCE_SHA256,
        expected_size=EXPECTED_SOURCE_SIZE,
        with_stat=True,
    )
    source_payload = load_json(source_data, "immutable preprobe source receipt")
    if (
        source_payload.get("schema_name")
        != "intersubmod.retrospective_running_source_identity.receipt"
        or source_payload.get("schema_version") != "1.2.0"
        or source_payload.get("pass") is not True
    ):
        raise PromotionError("Immutable preprobe source receipt contract differs")

    _, snapshot_data, snapshot_record = open_artifact(
        reader, SNAPSHOT, expected_sha256=EXPECTED_SNAPSHOT_SHA256
    )
    _, manifest_data, manifest_record = open_artifact(
        reader,
        TUMOR_REF_MANIFEST,
        expected_sha256=EXPECTED_MANIFEST_SHA256,
        expected_mode="0o664",
    )
    snapshot_payload = load_json(snapshot_data, "source snapshot")
    manifest_payload = load_json(manifest_data, "tumor-REF run manifest")
    if (
        snapshot_payload.get("pass") is not True
        or manifest_payload.get("pass") is not True
        or manifest_payload.get("counts", {}).get("processed_sites") != 102_842
    ):
        raise PromotionError("Snapshot or tumor-REF producer manifest is not full-scope PASS")

    focal_fd, _, focal_record = open_artifact(
        reader,
        FOCAL_ALT_CLUSTER_LIB,
        expected_sha256=EXPECTED_FOCAL_SHA256,
        expected_size=EXPECTED_FOCAL_SIZE,
        with_stat=True,
    )
    during = snapshot_payload.get("source_identity_during_execution", {}).get(
        "focal_alt_cluster_lib"
    )
    if not isinstance(during, Mapping):
        raise PromotionError("Source snapshot lacks focal library identity")
    unchanged_fields = ("path", "size_bytes", "sha256", "device", "inode", "mtime_ns")
    if any(focal_record.get(field) != during.get(field) for field in unchanged_fields):
        raise PromotionError("Focal source differs beyond the authorized mode-only transition")
    if (
        during.get("mode") != "0o664"
        or focal_record.get("mode") != "0o444"
        or focal_record.get("ctime_ns") == during.get("ctime_ns")
    ):
        raise PromotionError("Focal source mode/ctime transition contract differs")

    _, prior_data, prior_record = open_artifact(
        reader,
        PRIOR_REVIEW,
        expected_sha256=EXPECTED_PRIOR_REVIEW_SHA256,
        expected_mode="0o664",
    )
    prior_text = prior_data.decode("utf-8")
    if (
        "正式 receipt 仍依設計等待 completion runner 建立，不能以 pre-probe 取代" not in prior_text
        or "pre-probe 明確不可替代" not in prior_text
    ):
        raise PromotionError("Prior review no longer contains the superseded policy text")

    archive_evidence = validate_archive_evidence(reader)
    v3 = validate_signed_release_authority(
        authority_module,
        reader,
        authority_path=V3_AUTHORITY,
        approval_path=V3_APPROVAL,
        signature_path=V3_SIGNATURE,
        expected_authority_id=EXPECTED_V3_ID,
        expected_public_key_path=V3_PUBLIC_KEY,
        expected_public_key_sha256=V3_PUBLIC_KEY_SHA256,
        expected_private_key_path=V3_PRIVATE_KEY,
    )
    v7 = validate_signed_release_authority(
        authority_module,
        reader,
        authority_path=V7_AUTHORITY,
        approval_path=V7_APPROVAL,
        signature_path=V7_SIGNATURE,
        expected_authority_id=EXPECTED_V7_ID,
        expected_public_key_path=V7_PUBLIC_KEY,
        expected_public_key_sha256=V7_PUBLIC_KEY_SHA256,
        expected_private_key_path=V7_PRIVATE_KEY,
    )
    v3_payload = load_json(reader.open(V3_AUTHORITY, include_mode=True)[1], "v3 authority replay")
    v3_focal = v3_payload.get("sources", {}).get("focal_alt_cluster_lib")
    focal_basic = {
        key: focal_record[key] for key in ("path", "size_bytes", "sha256", "mode")
    }
    if v3_focal != focal_basic:
        raise PromotionError("v3 authority does not bind the current focal source")

    v7_validation = authority_module.validate_release_source_authority(
        authority_module.EXPECTED_SOURCE_PATHS
    )
    if (
        v7_validation.get("pass") is not True
        or v7_validation.get("authority_id") != EXPECTED_V7_ID
        or v7_validation.get("source_set_sha256")
        != "f5ba8a9e786971c4261b51e283a0f9df6e807a8aca695bf5041271fee5420f58"
    ):
        raise PromotionError("Current v7 source authority validation failed")
    builder_expected = load_json(
        reader.open(V7_AUTHORITY, include_mode=True)[1], "v7 authority replay"
    ).get("sources", {}).get("final_dataset_builder")
    builder_fd, builder_data, builder_record = open_artifact(reader, FINAL_DATASET_BUILDER)
    if builder_expected != builder_record:
        raise PromotionError("Final builder differs from v7 authority")
    require_active_cache_empty()
    builder = load_builder_gate_from_bound_ast(builder_data)
    require_active_cache_empty()
    pre_status, _ = builder.load_tumor_ref_source_identity_attestation(
        SOURCE_RECEIPT,
        TUMOR_REF_MANIFEST,
        manifest_payload,
    )
    if (
        pre_status.get("release_gate_pass") is not True
        or pre_status.get("publishable_task_b_release") is not True
    ):
        raise PromotionError("Final builder rejected the immutable source receipt")

    return {
        "source_fd": source_fd,
        "source_data": source_data,
        "source_record": source_record,
        "source_payload": source_payload,
        "snapshot": snapshot_record,
        "manifest": manifest_record,
        "manifest_payload": manifest_payload,
        "focal_transition": {
            "during_execution": dict(during),
            "current": focal_record,
            "unchanged_fields": list(unchanged_fields),
            "observed_differences": {
                "mode": {"during": "0o664", "current": "0o444"},
                "ctime_ns": {
                    "during": during.get("ctime_ns"),
                    "current": focal_record.get("ctime_ns"),
                },
            },
            "interpretation": (
                "Observed bytes, path, device, inode, size, and mtime are unchanged; "
                "mode was later tightened from 0664 to 0444. This does not prove that "
                "no unobserved transient mutation ever occurred."
            ),
            "bound_fd": focal_fd,
        },
        "prior_review": prior_record,
        "archive_evidence": archive_evidence,
        "v3": v3,
        "v7": v7,
        "v7_validation": {
            "authority_id": v7_validation["authority_id"],
            "source_set_sha256": v7_validation["source_set_sha256"],
            "required_roles_verified": v7_validation["required_roles_verified"],
            "binding_lease": v7_validation["validation_binding_lease"],
            "pass": True,
        },
        "builder": builder,
        "builder_fd": builder_fd,
        "builder_record": builder_record,
        "builder_gate_ast_contract": {
            "constants": sorted(BUILDER_GATE_CONSTANTS),
            "definitions": sorted(BUILDER_GATE_DEFINITIONS),
            "gate_first_line": 4729,
            "no_local_module_imports_executed": True,
        },
        "python_runtime": python_record,
        "builder_pre_status": pre_status,
    }


def collect_reviewed_inputs(reader: Any) -> dict[str, Any]:
    reviewed_sources: dict[str, dict[str, Any]] = {}
    source_inodes: set[tuple[int, int]] = set()
    for role, path in REVIEWED_SOURCE_ROLES.items():
        runtime_contract = REVIEWED_SOURCE_RUNTIME_CONTRACTS.get(role)
        descriptor, data, record = open_artifact(
            reader,
            path,
            expected_sha256=(runtime_contract[0] if runtime_contract else None),
            expected_mode=(runtime_contract[1] if runtime_contract else "0o444"),
        )
        reviewed_sources[role] = record
        opened = os.fstat(descriptor)
        source_inodes.add((opened.st_dev, opened.st_ino))
        if role == "promotion_tool":
            execution_binding = verify_executing_source_matches_bound_bytes(data)
    observed_asset_parts = tuple(
        sorted(
            path
            for path in PORTABLE_ASSET_ROOT.iterdir()
            if path.name.startswith("portable-artifact-reader.html.gz.b64.part")
        )
    )
    if observed_asset_parts != PORTABLE_READER_ASSET_PARTS:
        raise PromotionError("Portable reader asset part membership drift")
    if len(source_inodes) != len(REVIEWED_SOURCE_ROLES):
        raise PromotionError("Reviewed tools must be distinct regular files")
    trusted_keys: dict[str, dict[str, Any]] = {}
    expected_key_sha = {
        "authorization_public_key": AUTHORIZATION_PUBLIC_KEY_SHA256,
        "completion_public_key": COMPLETION_PUBLIC_KEY_SHA256,
        "continuation_public_key": CONTINUATION_PUBLIC_KEY_SHA256,
    }
    for role, path in TRUSTED_KEY_ROLES.items():
        _, _, record = open_artifact(
            reader, path, expected_sha256=expected_key_sha[role]
        )
        trusted_keys[role] = record
    return {
        "reviewed_sources": reviewed_sources,
        "trusted_keys": trusted_keys,
        "execution_binding": execution_binding,
    }


def validate_reviews(
    reader: Any,
    reviewed_sources: Mapping[str, Mapping[str, Any]],
    trusted_keys: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    reviewer_ids: set[str] = set()
    inodes: set[tuple[int, int]] = set()
    for role, slot in REVIEW_SLOTS.items():
        record, payload, inode = validate_review(
            reader,
            role=role,
            slot=slot,
            reviewed_sources=reviewed_sources,
            trusted_keys=trusted_keys,
        )
        reviewer_id = str(payload["reviewer_agent_id"])
        if reviewer_id in reviewer_ids or inode in inodes:
            raise PromotionError("Promotion reviews must be distinct")
        reviewer_ids.add(reviewer_id)
        inodes.add(inode)
        records[role] = {
            "artifact": record,
            "reviewer": payload["reviewer"],
            "reviewer_agent_id": reviewer_id,
            "verdict": "APPROVE",
            "findings_closed": True,
        }
    return records


def authority_for_receipt(value: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: item
        for key, item in value.items()
        if key not in {"bound_descriptors", "openssl"}
    }


def core_evidence_for_receipt(core: Mapping[str, Any]) -> dict[str, Any]:
    focal_transition = {
        key: value
        for key, value in core["focal_transition"].items()
        if key != "bound_fd"
    }
    return {
        "source_snapshot": core["snapshot"],
        "tumor_ref_manifest": core["manifest"],
        "focal_source_identity_transition": focal_transition,
        "prior_policy_review": core["prior_review"],
        "failed_attempt_archive": core["archive_evidence"],
        "historical_v3_authority": authority_for_receipt(core["v3"]),
        "current_v7_authority": authority_for_receipt(core["v7"]),
        "current_v7_runtime_validation": core["v7_validation"],
        "final_builder": core["builder_record"],
        "final_builder_gate_execution": core["builder_gate_ast_contract"],
        "python_runtime": core["python_runtime"],
        "python_cache_contract": {
            "paths": {key: str(value) for key, value in PYTHON_CACHE_ROOTS.items()},
            "mode": "0o700",
            "active_mode_cache_fresh_and_empty_before_and_after_import": True,
            "bytecode_writes_disabled": True,
        },
        "final_builder_source_path_gate": {
            "status": core["builder_pre_status"].get("status"),
            "release_gate_pass": core["builder_pre_status"].get("release_gate_pass"),
            "publishable_task_b_release": core["builder_pre_status"].get(
                "publishable_task_b_release"
            ),
        },
    }


def publish_exact_copy(
    reader: Any,
    core: Mapping[str, Any],
    *,
    publish_provisional_completion: Callable[
        [dict[str, Any], dict[str, Any]], dict[str, Any]
    ],
    precommit: Callable[[], None],
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    require_path_absent(CANONICAL_RECEIPT, "Canonical receipt")
    if not hasattr(os, "O_TMPFILE"):
        raise PromotionError("O_TMPFILE is required for staged canonical publication")
    source_fd = int(core["source_fd"])
    source_data = bytes(core["source_data"])
    source_opened = os.fstat(source_fd)
    if CANONICAL_RECEIPT.parent.resolve(strict=True) != CANONICAL_RECEIPT.parent or stat.S_ISLNK(
        os.lstat(CANONICAL_RECEIPT.parent).st_mode
    ):
        raise PromotionError("Canonical receipt parent is not canonical")
    parent_flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0)
    if hasattr(os, "O_NOFOLLOW"):
        parent_flags |= os.O_NOFOLLOW
    parent_fd = os.open(CANONICAL_RECEIPT.parent, parent_flags)
    parent_opened = os.fstat(parent_fd)
    target_fd = -1
    completion_record: dict[str, Any] | None = None
    projected_record: dict[str, Any] | None = None
    try:
        target_fd = os.open(
            ".", os.O_RDWR | os.O_TMPFILE | os.O_CLOEXEC, 0o400, dir_fd=parent_fd
        )
        offset = 0
        while offset < len(source_data):
            written = os.write(target_fd, source_data[offset:])
            if written <= 0:
                raise PromotionError("Short write during canonical receipt publication")
            offset += written
        os.fchmod(target_fd, 0o444)
        os.fsync(target_fd)
        opened = os.fstat(target_fd)
        target_data = b"".join(
            os.pread(target_fd, min(1024 * 1024, opened.st_size - position), position)
            for position in range(0, opened.st_size, 1024 * 1024)
        )
        if (
            not stat.S_ISREG(opened.st_mode)
            or stat.S_IMODE(opened.st_mode) != 0o444
            or opened.st_nlink != 0
            or opened.st_ino == source_opened.st_ino
            or target_data != source_data
        ):
            raise PromotionError("Staged canonical receipt metadata or bytes differ")
        staged_path = Path(f"/proc/self/fd/{target_fd}")
        post_status, _ = core["builder"].load_tumor_ref_source_identity_attestation(
            staged_path,
            TUMOR_REF_MANIFEST,
            core["manifest_payload"],
        )
        if (
            post_status.get("release_gate_pass") is not True
            or post_status.get("publishable_task_b_release") is not True
        ):
            raise PromotionError("Final builder rejected the staged canonical receipt")
        projected_record = {
            "path": str(CANONICAL_RECEIPT),
            "size_bytes": opened.st_size,
            "sha256": sha256_bytes(target_data),
            "mode": "0o444",
        }
        completion_record = publish_provisional_completion(
            projected_record, post_status
        )
        _, _, rebound_completion = open_artifact(reader, COMPLETION_RECEIPT)
        if rebound_completion != completion_record:
            raise PromotionError("Provisional completion receipt identity rebound")
        require_downstream_absent()
        require_path_absent(PROMOTION_INCIDENT, "Promotion incident")
        reader.require_paths_still_bound()
        precommit()
        parent_live = os.stat(CANONICAL_RECEIPT.parent, follow_symlinks=False)
        if (
            parent_opened.st_dev,
            parent_opened.st_ino,
            parent_opened.st_mode,
        ) != (parent_live.st_dev, parent_live.st_ino, parent_live.st_mode):
            raise PromotionError("Canonical receipt parent path rebound")
        link_fd_no_replace(target_fd, parent_fd, CANONICAL_RECEIPT.name)
        committed = os.fstat(target_fd)
        linked = os.stat(
            CANONICAL_RECEIPT.name, dir_fd=parent_fd, follow_symlinks=False
        )
        reopened_flags = os.O_RDONLY | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            reopened_flags |= os.O_NOFOLLOW
        reopened_fd = os.open(
            CANONICAL_RECEIPT.name, reopened_flags, dir_fd=parent_fd
        )
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
            committed_identity = (
                committed.st_dev,
                committed.st_ino,
                committed.st_mode,
                committed.st_nlink,
                committed.st_size,
                committed.st_mtime_ns,
                committed.st_ctime_ns,
            )
            if (
                committed.st_nlink != 1
                or committed_identity
                != (
                    linked.st_dev,
                    linked.st_ino,
                    linked.st_mode,
                    linked.st_nlink,
                    linked.st_size,
                    linked.st_mtime_ns,
                    linked.st_ctime_ns,
                )
                or committed_identity
                != (
                    reopened.st_dev,
                    reopened.st_ino,
                    reopened.st_mode,
                    reopened.st_nlink,
                    reopened.st_size,
                    reopened.st_mtime_ns,
                    reopened.st_ctime_ns,
                )
                or reopened_data != source_data
            ):
                raise PromotionError("Canonical receipt post-link identity drift")
        finally:
            os.close(reopened_fd)
        os.fsync(parent_fd)
        parent_after = os.stat(CANONICAL_RECEIPT.parent, follow_symlinks=False)
        if (
            parent_opened.st_dev,
            parent_opened.st_ino,
            parent_opened.st_mode,
        ) != (parent_after.st_dev, parent_after.st_ino, parent_after.st_mode):
            raise PromotionError("Canonical parent changed after commit")
        if CANONICAL_RECEIPT.resolve(strict=True) != CANONICAL_RECEIPT:
            raise PromotionError("Canonical receipt path is not canonical after commit")
        durable_fd = os.open(
            CANONICAL_RECEIPT.name, reopened_flags, dir_fd=parent_fd
        )
        try:
            durable_source = os.fstat(target_fd)
            durable_opened = os.fstat(durable_fd)
            durable_linked = os.stat(
                CANONICAL_RECEIPT.name, dir_fd=parent_fd, follow_symlinks=False
            )
            durable_data = b"".join(
                os.pread(
                    durable_fd,
                    min(1024 * 1024, durable_opened.st_size - position),
                    position,
                )
                for position in range(0, durable_opened.st_size, 1024 * 1024)
            )
            durable_identity = (
                durable_source.st_dev,
                durable_source.st_ino,
                durable_source.st_mode,
                durable_source.st_nlink,
                durable_source.st_size,
                durable_source.st_mtime_ns,
                durable_source.st_ctime_ns,
            )
            if (
                durable_source.st_nlink != 1
                or durable_identity
                != tuple(
                    getattr(durable_opened, field)
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
                or durable_identity
                != tuple(
                    getattr(durable_linked, field)
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
                or durable_data != source_data
            ):
                raise PromotionError("Canonical receipt post-fsync identity drift")
        finally:
            os.close(durable_fd)
    finally:
        if target_fd >= 0:
            os.close(target_fd)
        os.close(parent_fd)
    _, reopened_data, target_record = open_artifact(
        reader,
        CANONICAL_RECEIPT,
        expected_sha256=EXPECTED_SOURCE_SHA256,
        expected_size=EXPECTED_SOURCE_SIZE,
        with_stat=True,
    )
    if completion_record is None or projected_record is None:
        raise PromotionError("Canonical transaction did not publish its provisional receipt")
    if (
        reopened_data != source_data
        or target_record["inode"] == core["source_record"]["inode"]
        or target_record["device"] != core["source_record"]["device"]
        or {
            key: target_record[key]
            for key in ("path", "size_bytes", "sha256", "mode")
        }
        != projected_record
    ):
        raise PromotionError("Canonical receipt is not an independent byte-exact local copy")
    return target_record, post_status, completion_record


def reopen_producer_identities(
    reader: Any,
    inputs: Mapping[str, Any],
    core: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    _, _, producer_source = reopen_full_stat_artifact(
        reader,
        PROMOTION_TOOL,
        expected_basic=inputs["reviewed_sources"]["promotion_tool"],
    )
    _, _, python_runtime = reopen_full_stat_artifact(
        reader,
        PRIMARY_PYTHON_RUNTIME,
        expected_basic=core["python_runtime"],
        expected_sha256=EXPECTED_PYTHON_SHA256,
        expected_mode="0o775",
    )
    return producer_source, python_runtime


def build_prepare_exit_attestation(
    *,
    created_at_utc: str,
    producer_session: Mapping[str, Any],
    child_wait: Mapping[str, Any],
    authorization_record: Mapping[str, Any],
    producer_source: Mapping[str, Any],
    python_runtime: Mapping[str, Any],
    public_key: Mapping[str, Any],
    private_key_pre_signature: Mapping[str, Any],
) -> dict[str, Any]:
    validate_producer_session(producer_session, expected_mode="--prepare-authorization")
    validate_child_wait(
        child_wait,
        expected_pid=producer_session["child_pid"],
        require_success=True,
    )
    validate_full_stat_artifact(dict(authorization_record), "prepare authorization")
    validate_full_stat_artifact(dict(producer_source), "prepare producer source")
    validate_full_stat_artifact(dict(python_runtime), "prepare Python runtime")
    return {
        "authorization_id": AUTHORIZATION_ID,
        "phase": "prepare_authorization",
        "supervisor_command": expected_command("--prepare-authorization"),
        "producer_session": dict(producer_session),
        "child_wait": dict(child_wait),
        "produced_artifacts": {
            "authorization": dict(authorization_record),
        },
        "producer_source": dict(producer_source),
        "python_runtime": dict(python_runtime),
        "signature_contract": exit_attestation_signature_contract(
            path=PREPARE_EXIT_ATTESTATION,
            schema=PREPARE_EXIT_ATTESTATION_SCHEMA,
            signature=AUTHORIZATION_SIGNATURE,
            public_key=public_key,
            private_key_pre_signature=private_key_pre_signature,
            legacy_payload_signature=LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
        ),
        "checks": dict(PREPARE_EXIT_CHECKS),
        "pass": True,
        "pass_semantics": PREPARE_EXIT_PASS_SEMANTICS,
        "created_at_utc": created_at_utc,
    }


def build_promote_exit_attestation(
    *,
    created_at_utc: str,
    producer_session: Mapping[str, Any],
    child_wait: Mapping[str, Any],
    completion_record: Mapping[str, Any],
    canonical_record: Mapping[str, Any],
    prepare_attestation_record: Mapping[str, Any],
    prepare_signature_record: Mapping[str, Any],
    producer_source: Mapping[str, Any],
    python_runtime: Mapping[str, Any],
    public_key: Mapping[str, Any],
    private_key_pre_signature: Mapping[str, Any],
) -> dict[str, Any]:
    validate_producer_session(producer_session, expected_mode="--promote")
    validate_child_wait(
        child_wait,
        expected_pid=producer_session["child_pid"],
        require_success=True,
    )
    for label, record in (
        ("promote completion receipt", completion_record),
        ("promote canonical receipt", canonical_record),
        ("promote prepare exit attestation", prepare_attestation_record),
        ("promote prepare attestation signature", prepare_signature_record),
        ("promote producer source", producer_source),
        ("promote Python runtime", python_runtime),
    ):
        validate_full_stat_artifact(dict(record), label)
    return {
        "authorization_id": AUTHORIZATION_ID,
        "phase": "promote",
        "supervisor_command": expected_command("--promote"),
        "producer_session": dict(producer_session),
        "child_wait": dict(child_wait),
        "produced_artifacts": {
            "completion_receipt": dict(completion_record),
            "canonical_receipt": dict(canonical_record),
            "prepare_exit_attestation": dict(prepare_attestation_record),
            "prepare_exit_attestation_signature": dict(prepare_signature_record),
        },
        "producer_source": dict(producer_source),
        "python_runtime": dict(python_runtime),
        "signature_contract": exit_attestation_signature_contract(
            path=PROMOTE_EXIT_ATTESTATION,
            schema=PROMOTE_EXIT_ATTESTATION_SCHEMA,
            signature=COMPLETION_SIGNATURE,
            public_key=public_key,
            private_key_pre_signature=private_key_pre_signature,
            legacy_payload_signature=LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
        ),
        "checks": dict(PROMOTE_EXIT_CHECKS),
        "pass": True,
        "pass_semantics": PROMOTE_EXIT_PASS_SEMANTICS,
        "created_at_utc": created_at_utc,
    }


def prepare_authorization(
    authority_module: Any,
    reader: Any,
    precommit: Callable[[], None],
    producer_session: Mapping[str, Any],
) -> dict[str, Any]:
    validate_producer_session(
        producer_session, expected_mode="--prepare-authorization"
    )
    for path, label in (
        (AUTHORIZATION, "Authorization"),
        (PREPARE_EXIT_ATTESTATION, "Prepare exit attestation"),
        (AUTHORIZATION_SIGNATURE, "Authorization signature"),
        (LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE, "Legacy authorization payload signature"),
        (COMPLETION_RECEIPT, "Completion receipt"),
        (PROMOTE_EXIT_ATTESTATION, "Promote exit attestation"),
        (COMPLETION_SIGNATURE, "Completion signature"),
        (LEGACY_COMPLETION_PAYLOAD_SIGNATURE, "Legacy completion payload signature"),
        (PROMOTION_INCIDENT, "Promotion incident"),
        (CANONICAL_RECEIPT, "Canonical receipt"),
        (VERIFICATION_RECEIPT, "Continuation verification"),
        (RUNNER_GATE_RECEIPT, "Runner gate receipt"),
        (RUNNER_GATE_LOG, "Runner gate log"),
    ):
        require_path_absent(path, label)
    absent_outputs = require_downstream_absent()
    inputs = collect_reviewed_inputs(reader)
    core = validate_core_evidence(authority_module, reader)
    reviews = validate_reviews(
        reader, inputs["reviewed_sources"], inputs["trusted_keys"]
    )
    _, authorization_private = open_private_key_metadata(
        reader, AUTHORIZATION_PRIVATE_KEY, expected_mode="0o400"
    )
    _, completion_private = open_private_key_metadata(
        reader, COMPLETION_PRIVATE_KEY, expected_mode="0o400"
    )
    _, continuation_private = open_private_key_metadata(
        reader, CONTINUATION_PRIVATE_KEY, expected_mode="0o400"
    )
    payload = {
        "schema_name": "intersubmod.tumor_ref_source_receipt_promotion_authorization",
        "schema_version": "3.0.0",
        "created_at_utc": now_utc(),
        "authorization_id": AUTHORIZATION_ID,
        "task_type": "B_comprehensive_validation",
        "status": "AUTHORIZED_FOR_SINGLE_BYTE_IDENTICAL_PROMOTION",
        "scope": "canonical_path_publication_only_no_scientific_payload_change",
        "producer_session": dict(producer_session),
        "producer_exit_attestation": exit_attestation_announcement(
            PREPARE_EXIT_ATTESTATION,
            PREPARE_EXIT_ATTESTATION_SCHEMA,
            inputs["trusted_keys"]["authorization_public_key"],
        ),
        "commands": {
            "prepare_authorization": expected_command("--prepare-authorization"),
            "promote": expected_command("--promote"),
            "continuation_verify": [
                str(PYTHON),
                "-I",
                "-B",
                str(CONTINUATION_VERIFIER),
                "--verify",
            ],
            "runner_gate_replay": [str(PYTHON), "-I", "-B", str(RUNNER_GATE_REPLAY)],
            "downstream_continuation": [
                str(PYTHON),
                "-I",
                "-B",
                str(DOWNSTREAM_CONTINUATION),
                "--supervise-and-sign",
            ],
        },
        "reviewed_sources": inputs["reviewed_sources"],
        "execution_binding": inputs["execution_binding"],
        "review_artifacts": reviews,
        "trusted_signing_keys": inputs["trusted_keys"],
        "source_receipt": core["source_record"],
        "canonical_target_path": str(CANONICAL_RECEIPT),
        "evidence": core_evidence_for_receipt(core),
        "authorization_signature_contract": exit_attestation_signature_contract(
            path=PREPARE_EXIT_ATTESTATION,
            schema=PREPARE_EXIT_ATTESTATION_SCHEMA,
            signature=AUTHORIZATION_SIGNATURE,
            public_key=inputs["trusted_keys"]["authorization_public_key"],
            private_key_pre_signature=authorization_private,
            legacy_payload_signature=LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
        ),
        "completion_signature_contract": exit_attestation_signature_contract(
            path=PROMOTE_EXIT_ATTESTATION,
            schema=PROMOTE_EXIT_ATTESTATION_SCHEMA,
            signature=COMPLETION_SIGNATURE,
            public_key=inputs["trusted_keys"]["completion_public_key"],
            private_key_pre_signature=completion_private,
            legacy_payload_signature=LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
        ),
        "continuation_gate": {
            "verifier": inputs["reviewed_sources"]["continuation_verifier"],
            "verification_receipt": str(VERIFICATION_RECEIPT),
            "runner_gate_replay": inputs["reviewed_sources"]["runner_gate_replay"],
            "runner_gate_receipt": str(RUNNER_GATE_RECEIPT),
            "downstream_continuation": inputs["reviewed_sources"][
                "downstream_continuation"
            ],
            "downstream_continuation_command": [
                str(PYTHON),
                "-I",
                "-B",
                str(DOWNSTREAM_CONTINUATION),
                "--supervise-and-sign",
            ],
            "supervised_child_command": [
                str(PYTHON),
                "-I",
                "-B",
                str(DOWNSTREAM_CONTINUATION),
                "--supervised-child",
            ],
            "supervisor_capability_protocol": SUPERVISOR_CAPABILITY_PROTOCOL,
            "threat_model": CONTINUATION_THREAT_MODEL,
            "terminal_signature_contract": {
                "algorithm": "Ed25519",
                "public_key": inputs["trusted_keys"]["continuation_public_key"],
                "private_key": continuation_private,
                "execution_receipt": str(CONTINUATION_TERMINAL_RECEIPT),
                "signed_artifact": str(CONTINUATION_EXIT_ATTESTATION),
                "signature": str(CONTINUATION_TERMINAL_SIGNATURE),
                "required_private_key_pre_signature_mode": "0o400",
                "required_private_key_post_signature_mode": "0o0",
                "supervisor_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--supervise-and-sign",
                ],
                "supervised_child_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--supervised-child",
                ],
                "signed_terminal_verification_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--verify-signed-terminal",
                ],
            },
            "policy": CONTINUATION_POLICY,
        },
        "historical_incident_disclosure": {
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
        },
        "downstream_output_absence": absent_outputs,
        "checks": dict(AUTHORIZATION_CHECKS),
        "pass": True,
        "pass_semantics": (
            "authorization_payload_only; parent exit-attestation publication and detached "
            "authorization-key signature remain mandatory; no canonical publication has "
            "occurred and this payload does not permit downstream execution"
        ),
    }
    canonical_payload = expected_authorization_payload(
        created_at_utc=payload["created_at_utc"],
        producer_session=producer_session,
        inputs=inputs,
        reviews=reviews,
        core=core,
        completion_private=completion_private,
        continuation_private=continuation_private,
        downstream_absence=absent_outputs,
    )
    if not strict_json_equal(payload, canonical_payload):
        raise PromotionError("Authorization producer payload differs from its exact v3 contract")
    reader.require_paths_still_bound()
    reader.retain_until_process_exit()
    authorization_record = write_exclusive_readonly(
        AUTHORIZATION, payload, precommit=precommit
    )
    return {
        "authorization": authorization_record,
        "prepare_exit_attestation_pending": str(PREPARE_EXIT_ATTESTATION),
        "attestation_signature_pending": str(AUTHORIZATION_SIGNATURE),
        "canonical_receipt_created": False,
        "pass": True,
    }


def expected_authorization_payload(
    *,
    created_at_utc: str,
    producer_session: Mapping[str, Any],
    inputs: Mapping[str, Any],
    reviews: Mapping[str, Any],
    core: Mapping[str, Any],
    completion_private: Mapping[str, Any],
    continuation_private: Mapping[str, Any],
    downstream_absence: list[str],
) -> dict[str, Any]:
    validate_producer_session(
        producer_session, expected_mode="--prepare-authorization"
    )
    authorization_private_pre = {
        "path": str(AUTHORIZATION_PRIVATE_KEY),
        "mode": "0o400",
    }
    return {
        "schema_name": "intersubmod.tumor_ref_source_receipt_promotion_authorization",
        "schema_version": "3.0.0",
        "created_at_utc": created_at_utc,
        "authorization_id": AUTHORIZATION_ID,
        "task_type": "B_comprehensive_validation",
        "status": "AUTHORIZED_FOR_SINGLE_BYTE_IDENTICAL_PROMOTION",
        "scope": "canonical_path_publication_only_no_scientific_payload_change",
        "producer_session": dict(producer_session),
        "producer_exit_attestation": exit_attestation_announcement(
            PREPARE_EXIT_ATTESTATION,
            PREPARE_EXIT_ATTESTATION_SCHEMA,
            inputs["trusted_keys"]["authorization_public_key"],
        ),
        "commands": {
            "prepare_authorization": expected_command("--prepare-authorization"),
            "promote": expected_command("--promote"),
            "continuation_verify": [
                str(PYTHON), "-I", "-B", str(CONTINUATION_VERIFIER), "--verify"
            ],
            "runner_gate_replay": [
                str(PYTHON), "-I", "-B", str(RUNNER_GATE_REPLAY)
            ],
            "downstream_continuation": [
                str(PYTHON),
                "-I",
                "-B",
                str(DOWNSTREAM_CONTINUATION),
                "--supervise-and-sign",
            ],
        },
        "reviewed_sources": inputs["reviewed_sources"],
        "execution_binding": inputs["execution_binding"],
        "review_artifacts": reviews,
        "trusted_signing_keys": inputs["trusted_keys"],
        "source_receipt": core["source_record"],
        "canonical_target_path": str(CANONICAL_RECEIPT),
        "evidence": core_evidence_for_receipt(core),
        "authorization_signature_contract": exit_attestation_signature_contract(
            path=PREPARE_EXIT_ATTESTATION,
            schema=PREPARE_EXIT_ATTESTATION_SCHEMA,
            signature=AUTHORIZATION_SIGNATURE,
            public_key=inputs["trusted_keys"]["authorization_public_key"],
            private_key_pre_signature=authorization_private_pre,
            legacy_payload_signature=LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
        ),
        "completion_signature_contract": exit_attestation_signature_contract(
            path=PROMOTE_EXIT_ATTESTATION,
            schema=PROMOTE_EXIT_ATTESTATION_SCHEMA,
            signature=COMPLETION_SIGNATURE,
            public_key=inputs["trusted_keys"]["completion_public_key"],
            private_key_pre_signature=completion_private,
            legacy_payload_signature=LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
        ),
        "continuation_gate": {
            "verifier": inputs["reviewed_sources"]["continuation_verifier"],
            "verification_receipt": str(VERIFICATION_RECEIPT),
            "runner_gate_replay": inputs["reviewed_sources"]["runner_gate_replay"],
            "runner_gate_receipt": str(RUNNER_GATE_RECEIPT),
            "downstream_continuation": inputs["reviewed_sources"][
                "downstream_continuation"
            ],
            "downstream_continuation_command": [
                str(PYTHON),
                "-I",
                "-B",
                str(DOWNSTREAM_CONTINUATION),
                "--supervise-and-sign",
            ],
            "supervised_child_command": [
                str(PYTHON),
                "-I",
                "-B",
                str(DOWNSTREAM_CONTINUATION),
                "--supervised-child",
            ],
            "supervisor_capability_protocol": SUPERVISOR_CAPABILITY_PROTOCOL,
            "threat_model": CONTINUATION_THREAT_MODEL,
            "terminal_signature_contract": {
                "algorithm": "Ed25519",
                "public_key": inputs["trusted_keys"]["continuation_public_key"],
                "private_key": continuation_private,
                "execution_receipt": str(CONTINUATION_TERMINAL_RECEIPT),
                "signed_artifact": str(CONTINUATION_EXIT_ATTESTATION),
                "signature": str(CONTINUATION_TERMINAL_SIGNATURE),
                "required_private_key_pre_signature_mode": "0o400",
                "required_private_key_post_signature_mode": "0o0",
                "supervisor_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--supervise-and-sign",
                ],
                "supervised_child_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--supervised-child",
                ],
                "signed_terminal_verification_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--verify-signed-terminal",
                ],
            },
            "policy": CONTINUATION_POLICY,
        },
        "historical_incident_disclosure": {
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
        },
        "downstream_output_absence": downstream_absence,
        "checks": dict(AUTHORIZATION_CHECKS),
        "pass": True,
        "pass_semantics": (
            "authorization_payload_only; parent exit-attestation publication and detached "
            "authorization-key signature remain mandatory; no canonical publication has "
            "occurred and this payload does not permit downstream execution"
        ),
    }


def validate_signed_authorization(
    authority_module: Any,
    reader: Any,
    *,
    inputs: Mapping[str, Any],
    reviews: Mapping[str, Any],
    core: Mapping[str, Any],
) -> dict[str, Any]:
    require_path_absent(
        LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
        "Legacy authorization payload signature",
    )
    require_path_absent(
        LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
        "Legacy completion payload signature",
    )
    authorization_fd, authorization_data, authorization_record = reopen_full_stat_artifact(
        reader, AUTHORIZATION
    )
    _, prepare_attestation_data, prepare_attestation_record = reopen_full_stat_artifact(
        reader, PREPARE_EXIT_ATTESTATION
    )
    signature_check = verify_signature(
        authority_module,
        reader,
        data_path=PREPARE_EXIT_ATTESTATION,
        signature_path=AUTHORIZATION_SIGNATURE,
        public_key_path=AUTHORIZATION_PUBLIC_KEY,
        public_key_sha256=AUTHORIZATION_PUBLIC_KEY_SHA256,
    )
    _, _, prepare_signature_record = reopen_full_stat_artifact(
        reader, AUTHORIZATION_SIGNATURE
    )
    _, retired_private = open_private_key_metadata(
        reader, AUTHORIZATION_PRIVATE_KEY, expected_mode="0o0"
    )
    _, completion_private = open_private_key_metadata(
        reader, COMPLETION_PRIVATE_KEY, expected_mode="0o400"
    )
    _, continuation_private = open_private_key_metadata(
        reader, CONTINUATION_PRIVATE_KEY, expected_mode="0o400"
    )
    payload = load_json(authorization_data, "signed promotion authorization")
    created_at = payload.get("created_at_utc")
    require_utc_timestamp(created_at, "signed authorization created_at_utc")
    producer_session = validate_producer_session(
        payload.get("producer_session"), expected_mode="--prepare-authorization"
    )
    current_absence = require_downstream_absent()
    expected_payload = {
        "schema_name": "intersubmod.tumor_ref_source_receipt_promotion_authorization",
        "schema_version": "3.0.0",
        "created_at_utc": created_at,
        "authorization_id": AUTHORIZATION_ID,
        "task_type": "B_comprehensive_validation",
        "status": "AUTHORIZED_FOR_SINGLE_BYTE_IDENTICAL_PROMOTION",
        "scope": "canonical_path_publication_only_no_scientific_payload_change",
        "producer_session": producer_session,
        "producer_exit_attestation": exit_attestation_announcement(
            PREPARE_EXIT_ATTESTATION,
            PREPARE_EXIT_ATTESTATION_SCHEMA,
            inputs["trusted_keys"]["authorization_public_key"],
        ),
        "commands": {
            "prepare_authorization": expected_command("--prepare-authorization"),
            "promote": expected_command("--promote"),
            "continuation_verify": [
                str(PYTHON), "-I", "-B", str(CONTINUATION_VERIFIER), "--verify"
            ],
            "runner_gate_replay": [
                str(PYTHON), "-I", "-B", str(RUNNER_GATE_REPLAY)
            ],
            "downstream_continuation": [
                str(PYTHON),
                "-I",
                "-B",
                str(DOWNSTREAM_CONTINUATION),
                "--supervise-and-sign",
            ],
        },
        "reviewed_sources": inputs["reviewed_sources"],
        "execution_binding": inputs["execution_binding"],
        "review_artifacts": reviews,
        "trusted_signing_keys": inputs["trusted_keys"],
        "source_receipt": core["source_record"],
        "canonical_target_path": str(CANONICAL_RECEIPT),
        "evidence": core_evidence_for_receipt(core),
        "authorization_signature_contract": exit_attestation_signature_contract(
            path=PREPARE_EXIT_ATTESTATION,
            schema=PREPARE_EXIT_ATTESTATION_SCHEMA,
            signature=AUTHORIZATION_SIGNATURE,
            public_key=inputs["trusted_keys"]["authorization_public_key"],
            private_key_pre_signature={
                "path": str(AUTHORIZATION_PRIVATE_KEY),
                "mode": "0o400",
            },
            legacy_payload_signature=LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
        ),
        "completion_signature_contract": exit_attestation_signature_contract(
            path=PROMOTE_EXIT_ATTESTATION,
            schema=PROMOTE_EXIT_ATTESTATION_SCHEMA,
            signature=COMPLETION_SIGNATURE,
            public_key=inputs["trusted_keys"]["completion_public_key"],
            private_key_pre_signature=completion_private,
            legacy_payload_signature=LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
        ),
        "continuation_gate": {
            "verifier": inputs["reviewed_sources"]["continuation_verifier"],
            "verification_receipt": str(VERIFICATION_RECEIPT),
            "runner_gate_replay": inputs["reviewed_sources"]["runner_gate_replay"],
            "runner_gate_receipt": str(RUNNER_GATE_RECEIPT),
            "downstream_continuation": inputs["reviewed_sources"][
                "downstream_continuation"
            ],
            "downstream_continuation_command": [
                str(PYTHON),
                "-I",
                "-B",
                str(DOWNSTREAM_CONTINUATION),
                "--supervise-and-sign",
            ],
            "supervised_child_command": [
                str(PYTHON),
                "-I",
                "-B",
                str(DOWNSTREAM_CONTINUATION),
                "--supervised-child",
            ],
            "supervisor_capability_protocol": SUPERVISOR_CAPABILITY_PROTOCOL,
            "threat_model": CONTINUATION_THREAT_MODEL,
            "terminal_signature_contract": {
                "algorithm": "Ed25519",
                "public_key": inputs["trusted_keys"]["continuation_public_key"],
                "private_key": continuation_private,
                "execution_receipt": str(CONTINUATION_TERMINAL_RECEIPT),
                "signed_artifact": str(CONTINUATION_EXIT_ATTESTATION),
                "signature": str(CONTINUATION_TERMINAL_SIGNATURE),
                "required_private_key_pre_signature_mode": "0o400",
                "required_private_key_post_signature_mode": "0o0",
                "supervisor_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--supervise-and-sign",
                ],
                "supervised_child_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--supervised-child",
                ],
                "signed_terminal_verification_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--verify-signed-terminal",
                ],
            },
            "policy": CONTINUATION_POLICY,
        },
        "historical_incident_disclosure": {
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
        },
        "downstream_output_absence": current_absence,
        "checks": dict(AUTHORIZATION_CHECKS),
        "pass": True,
        "pass_semantics": (
            "authorization_payload_only; parent exit-attestation publication and detached "
            "authorization-key signature remain mandatory; no canonical publication has "
            "occurred and this payload does not permit downstream execution"
        ),
    }
    canonical_expected_payload = expected_authorization_payload(
        created_at_utc=created_at,
        producer_session=producer_session,
        inputs=inputs,
        reviews=reviews,
        core=core,
        completion_private=completion_private,
        continuation_private=continuation_private,
        downstream_absence=current_absence,
    )
    if not strict_json_equal(expected_payload, canonical_expected_payload):
        raise PromotionError("Authorization validator contract differs from its v3 builder")
    if not strict_json_equal(payload, expected_payload):
        raise PromotionError(
            "Signed authorization is not the exact pre-publication contract"
        )
    prepare_payload = load_json(
        prepare_attestation_data, "signed prepare exit attestation"
    )
    if set(prepare_payload) != EXIT_ATTESTATION_KEYS:
        raise PromotionError("Prepare exit attestation top-level keys differ")
    prepare_created_at = prepare_payload.get("created_at_utc")
    if require_utc_timestamp(
        prepare_created_at, "prepare exit attestation created_at_utc"
    ) < require_utc_timestamp(created_at, "authorization created_at_utc"):
        raise PromotionError("Prepare exit attestation predates its authorization")
    producer_source, python_runtime = reopen_producer_identities(reader, inputs, core)
    expected_prepare_payload = build_prepare_exit_attestation(
        created_at_utc=prepare_created_at,
        producer_session=producer_session,
        child_wait=prepare_payload.get("child_wait"),
        authorization_record=authorization_record,
        producer_source=producer_source,
        python_runtime=python_runtime,
        public_key=inputs["trusted_keys"]["authorization_public_key"],
        private_key_pre_signature={
            "path": str(AUTHORIZATION_PRIVATE_KEY),
            "mode": "0o400",
        },
    )
    if not strict_json_equal(prepare_payload, expected_prepare_payload):
        raise PromotionError("Prepare exit attestation differs from its exact v1 contract")
    prepare_attestation_basic = {
        key: prepare_attestation_record[key]
        for key in ("path", "size_bytes", "sha256", "mode")
    }
    prepare_signature_basic = {
        key: prepare_signature_record[key]
        for key in ("path", "size_bytes", "sha256", "mode")
    }
    if (
        not strict_json_equal(signature_check["data"], prepare_attestation_basic)
        or not strict_json_equal(signature_check["signature"], prepare_signature_basic)
    ):
        raise PromotionError(
            "Authorization signature was not checked over the bound prepare attestation"
        )
    return {
        "artifact": authorization_record,
        "producer_exit_attestation": prepare_attestation_record,
        "signature": prepare_signature_record,
        "public_key": signature_check["public_key"],
        "retired_private_key": retired_private,
        "producer_session": producer_session,
        "signature_verified": True,
        "bound_fd": authorization_fd,
    }


def attest_prepare_child_output(
    authority_module: Any,
    reader: Any,
    precommit: Callable[[], None],
    producer_session: Mapping[str, Any],
    child_wait: Mapping[str, Any],
) -> dict[str, Any]:
    validate_producer_session(
        producer_session, expected_mode="--prepare-authorization"
    )
    validate_child_wait(
        child_wait,
        expected_pid=producer_session["child_pid"],
        require_success=True,
    )
    for path, label in (
        (PREPARE_EXIT_ATTESTATION, "Prepare exit attestation"),
        (AUTHORIZATION_SIGNATURE, "Authorization signature"),
        (LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE, "Legacy authorization payload signature"),
        (COMPLETION_RECEIPT, "Completion receipt"),
        (PROMOTE_EXIT_ATTESTATION, "Promote exit attestation"),
        (COMPLETION_SIGNATURE, "Completion signature"),
        (LEGACY_COMPLETION_PAYLOAD_SIGNATURE, "Legacy completion payload signature"),
        (PROMOTION_INCIDENT, "Promotion incident"),
        (CANONICAL_RECEIPT, "Canonical receipt"),
        (VERIFICATION_RECEIPT, "Continuation verification"),
        (RUNNER_GATE_RECEIPT, "Runner gate receipt"),
        (RUNNER_GATE_LOG, "Runner gate log"),
    ):
        require_path_absent(path, label)
    downstream_absence = require_downstream_absent()
    inputs = collect_reviewed_inputs(reader)
    core = validate_core_evidence(authority_module, reader)
    reviews = validate_reviews(
        reader, inputs["reviewed_sources"], inputs["trusted_keys"]
    )
    _, authorization_private = open_private_key_metadata(
        reader, AUTHORIZATION_PRIVATE_KEY, expected_mode="0o400"
    )
    _, completion_private = open_private_key_metadata(
        reader, COMPLETION_PRIVATE_KEY, expected_mode="0o400"
    )
    _, continuation_private = open_private_key_metadata(
        reader, CONTINUATION_PRIVATE_KEY, expected_mode="0o400"
    )
    _, authorization_data, authorization_record = reopen_full_stat_artifact(
        reader, AUTHORIZATION
    )
    authorization_payload = load_json(
        authorization_data, "parent-reopened promotion authorization"
    )
    authorization_created_at = authorization_payload.get("created_at_utc")
    require_utc_timestamp(
        authorization_created_at, "parent-reopened authorization created_at_utc"
    )
    if not strict_json_equal(
        authorization_payload.get("producer_session"), dict(producer_session)
    ):
        raise PromotionError("Authorization producer session differs in the parent")
    expected_payload = expected_authorization_payload(
        created_at_utc=authorization_created_at,
        producer_session=producer_session,
        inputs=inputs,
        reviews=reviews,
        core=core,
        completion_private=completion_private,
        continuation_private=continuation_private,
        downstream_absence=downstream_absence,
    )
    if not strict_json_equal(authorization_payload, expected_payload):
        raise PromotionError("Parent rejected the child authorization payload")
    producer_source, python_runtime = reopen_producer_identities(reader, inputs, core)
    attestation_created_at = now_utc()
    if require_utc_timestamp(
        attestation_created_at, "prepare exit attestation created_at_utc"
    ) < require_utc_timestamp(
        authorization_created_at, "authorization created_at_utc"
    ):
        raise PromotionError("Prepare exit attestation timestamp ordering failed")
    attestation = build_prepare_exit_attestation(
        created_at_utc=attestation_created_at,
        producer_session=producer_session,
        child_wait=child_wait,
        authorization_record=authorization_record,
        producer_source=producer_source,
        python_runtime=python_runtime,
        public_key=inputs["trusted_keys"]["authorization_public_key"],
        private_key_pre_signature=authorization_private,
    )
    if set(attestation) != EXIT_ATTESTATION_KEYS:
        raise PromotionError("Prepare exit attestation builder emitted unexpected keys")
    reader.require_paths_still_bound()
    reader.retain_until_process_exit()

    def attestation_precommit() -> None:
        require_path_absent(
            LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
            "Legacy authorization payload signature",
        )
        require_path_absent(
            LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
            "Legacy completion payload signature",
        )
        reader.require_paths_still_bound()
        precommit()

    attestation_record = write_exclusive_readonly(
        PREPARE_EXIT_ATTESTATION,
        attestation,
        precommit=attestation_precommit,
    )
    return {
        "authorization": authorization_record,
        "prepare_exit_attestation": attestation_record,
        "authorization_attestation_signature_pending": str(AUTHORIZATION_SIGNATURE),
        "canonical_receipt_created": False,
        "pass": True,
    }


def promote(
    authority_module: Any,
    reader: Any,
    precommit: Callable[[], None],
    producer_session: Mapping[str, Any],
) -> dict[str, Any]:
    validate_producer_session(producer_session, expected_mode="--promote")
    for path, label in (
        (COMPLETION_RECEIPT, "Completion receipt"),
        (PROMOTE_EXIT_ATTESTATION, "Promote exit attestation"),
        (COMPLETION_SIGNATURE, "Completion signature"),
        (LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE, "Legacy authorization payload signature"),
        (LEGACY_COMPLETION_PAYLOAD_SIGNATURE, "Legacy completion payload signature"),
        (PROMOTION_INCIDENT, "Promotion incident"),
        (CANONICAL_RECEIPT, "Canonical receipt"),
        (VERIFICATION_RECEIPT, "Continuation verification"),
        (RUNNER_GATE_RECEIPT, "Runner gate receipt"),
        (RUNNER_GATE_LOG, "Runner gate log"),
    ):
        require_path_absent(path, label)
    absent_outputs = require_downstream_absent()
    inputs = collect_reviewed_inputs(reader)
    reviews = validate_reviews(
        reader, inputs["reviewed_sources"], inputs["trusted_keys"]
    )
    core = validate_core_evidence(authority_module, reader)
    authorization = validate_signed_authorization(
        authority_module,
        reader,
        inputs=inputs,
        reviews=reviews,
        core=core,
    )
    _, completion_private = open_private_key_metadata(
        reader, COMPLETION_PRIVATE_KEY, expected_mode="0o400"
    )
    authorization_public = {
        key: value for key, value in authorization.items() if key != "bound_fd"
    }

    def publish_completion(
        projected_target: dict[str, Any], staged_post_status: dict[str, Any]
    ) -> dict[str, Any]:
        terminal_absence = require_downstream_absent()
        if terminal_absence != absent_outputs:
            raise PromotionError("Downstream absence contract changed during promotion")
        payload = {
            "schema_name": "intersubmod.tumor_ref_source_receipt_promotion",
            "schema_version": "3.0.0",
            "created_at_utc": now_utc(),
            "authorization_id": AUTHORIZATION_ID,
            "task_type": "B_comprehensive_validation",
            "status": "PROVISIONAL_WRITE_AHEAD_CANONICAL_PUBLICATION_INTENT",
            "scope": "canonical_path_publication_only_no_scientific_payload_change",
            "producer_session": dict(producer_session),
            "producer_exit_attestation": exit_attestation_announcement(
                PROMOTE_EXIT_ATTESTATION,
                PROMOTE_EXIT_ATTESTATION_SCHEMA,
                inputs["trusted_keys"]["completion_public_key"],
            ),
            "command": expected_command("--promote"),
            "authorization": authorization_public,
            "code": inputs["reviewed_sources"],
            "execution_binding": inputs["execution_binding"],
            "source_receipt": core["source_record"],
            "canonical_receipt": projected_target,
            "evidence": core_evidence_for_receipt(core),
            "review_artifacts": reviews,
            "builder_gate_replay": {
                "source_path_status": core["builder_pre_status"].get("status"),
                "canonical_path_status": staged_post_status.get("status"),
                "release_gate_pass": True,
                "publishable_task_b_release": True,
                "builder_loaded_from_v7_bound_bytes": True,
                "bytecode_disabled": True,
            },
            "historical_incident_disclosure": {
                "original_completion_runner_end_to_end_pass_claimed": False,
                "original_post_run_verifier_success_claimed": False,
                "canonical_receipt_created_by": (
                    "signed_v3_exit_attested_byte_identical_promotion"
                ),
                "failed_attempt_archive": core["archive_evidence"],
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
            },
            "governance": {
                "prior_review": core["prior_review"],
                "narrow_supersession": (
                    "Only the canonical-emitter requirement is superseded; all content checks, "
                    "source authorities, downstream gate replay, and release signatures remain."
                ),
                "continuation_verifier": inputs["reviewed_sources"][
                    "continuation_verifier"
                ],
                "runner_gate_replay": inputs["reviewed_sources"][
                    "runner_gate_replay"
                ],
                "downstream_continuation": inputs["reviewed_sources"][
                    "downstream_continuation"
                ],
                "downstream_continuation_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--supervise-and-sign",
                ],
                "supervised_child_command": [
                    str(PYTHON),
                    "-I",
                    "-B",
                    str(DOWNSTREAM_CONTINUATION),
                    "--supervised-child",
                ],
                "continuation_verification_receipt": str(VERIFICATION_RECEIPT),
                "runner_gate_receipt": str(RUNNER_GATE_RECEIPT),
                "downstream_output_slots_absent_at_promotion": terminal_absence,
                "canonical_publication_order": (
                    "provisional_completion_receipt_then_canonical_no_replace_commit"
                ),
                "canonical_alone_grants_downstream_authority": False,
                "write_ahead_receipt_marks_transaction_incomplete_until": (
                    "detached_completion_signature_plus_live_canonical_verification"
                ),
            },
            "signature_contract": exit_attestation_signature_contract(
                path=PROMOTE_EXIT_ATTESTATION,
                schema=PROMOTE_EXIT_ATTESTATION_SCHEMA,
                signature=COMPLETION_SIGNATURE,
                public_key=inputs["trusted_keys"]["completion_public_key"],
                private_key_pre_signature=completion_private,
                legacy_payload_signature=LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
            ),
            "checks": dict(COMPLETION_CHECKS),
            "pass": True,
            "pass_semantics": (
                "completion_payload_only; parent promote exit-attestation publication and "
                "detached completion-key signature remain mandatory, followed by continuation "
                "verification, runner-gate evidence, and the signed downstream terminal receipt"
            ),
        }
        canonical_payload = expected_completion_payload(
            created_at_utc=payload["created_at_utc"],
            producer_session=producer_session,
            authorization=authorization_public,
            inputs=inputs,
            core=core,
            reviews=reviews,
            canonical_receipt=projected_target,
            canonical_path_status=staged_post_status,
            completion_private=completion_private,
            downstream_absence=terminal_absence,
        )
        if not strict_json_equal(payload, canonical_payload):
            raise PromotionError("Completion producer payload differs from its exact v3 contract")
        reader.require_paths_still_bound()
        def completion_precommit() -> None:
            reader.require_paths_still_bound()
            precommit()

        return write_exclusive_readonly(
            COMPLETION_RECEIPT, payload, precommit=completion_precommit
        )

    reader.require_paths_still_bound()
    reader.retain_until_process_exit()
    target_record, _, completion_record = publish_exact_copy(
        reader,
        core,
        publish_provisional_completion=publish_completion,
        precommit=precommit,
    )
    return {
        "canonical_receipt": target_record,
        "completion_receipt": completion_record,
        "promote_exit_attestation_pending": str(PROMOTE_EXIT_ATTESTATION),
        "completion_attestation_signature_pending": str(COMPLETION_SIGNATURE),
        "continuation_verification_pending": str(VERIFICATION_RECEIPT),
        "downstream_started": False,
        "pass": True,
    }


def expected_completion_payload(
    *,
    created_at_utc: str,
    producer_session: Mapping[str, Any],
    authorization: Mapping[str, Any],
    inputs: Mapping[str, Any],
    core: Mapping[str, Any],
    reviews: Mapping[str, Any],
    canonical_receipt: Mapping[str, Any],
    canonical_path_status: Mapping[str, Any],
    completion_private: Mapping[str, Any],
    downstream_absence: list[str],
) -> dict[str, Any]:
    validate_producer_session(producer_session, expected_mode="--promote")
    return {
        "schema_name": "intersubmod.tumor_ref_source_receipt_promotion",
        "schema_version": "3.0.0",
        "created_at_utc": created_at_utc,
        "authorization_id": AUTHORIZATION_ID,
        "task_type": "B_comprehensive_validation",
        "status": "PROVISIONAL_WRITE_AHEAD_CANONICAL_PUBLICATION_INTENT",
        "scope": "canonical_path_publication_only_no_scientific_payload_change",
        "producer_session": dict(producer_session),
        "producer_exit_attestation": exit_attestation_announcement(
            PROMOTE_EXIT_ATTESTATION,
            PROMOTE_EXIT_ATTESTATION_SCHEMA,
            inputs["trusted_keys"]["completion_public_key"],
        ),
        "command": expected_command("--promote"),
        "authorization": dict(authorization),
        "code": inputs["reviewed_sources"],
        "execution_binding": inputs["execution_binding"],
        "source_receipt": core["source_record"],
        "canonical_receipt": dict(canonical_receipt),
        "evidence": core_evidence_for_receipt(core),
        "review_artifacts": reviews,
        "builder_gate_replay": {
            "source_path_status": core["builder_pre_status"].get("status"),
            "canonical_path_status": canonical_path_status.get("status"),
            "release_gate_pass": True,
            "publishable_task_b_release": True,
            "builder_loaded_from_v7_bound_bytes": True,
            "bytecode_disabled": True,
        },
        "historical_incident_disclosure": {
            "original_completion_runner_end_to_end_pass_claimed": False,
            "original_post_run_verifier_success_claimed": False,
            "canonical_receipt_created_by": (
                "signed_v3_exit_attested_byte_identical_promotion"
            ),
            "failed_attempt_archive": core["archive_evidence"],
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
        },
        "governance": {
            "prior_review": core["prior_review"],
            "narrow_supersession": (
                "Only the canonical-emitter requirement is superseded; all content checks, "
                "source authorities, downstream gate replay, and release signatures remain."
            ),
            "continuation_verifier": inputs["reviewed_sources"][
                "continuation_verifier"
            ],
            "runner_gate_replay": inputs["reviewed_sources"]["runner_gate_replay"],
            "downstream_continuation": inputs["reviewed_sources"][
                "downstream_continuation"
            ],
            "downstream_continuation_command": [
                str(PYTHON),
                "-I",
                "-B",
                str(DOWNSTREAM_CONTINUATION),
                "--supervise-and-sign",
            ],
            "supervised_child_command": [
                str(PYTHON),
                "-I",
                "-B",
                str(DOWNSTREAM_CONTINUATION),
                "--supervised-child",
            ],
            "continuation_verification_receipt": str(VERIFICATION_RECEIPT),
            "runner_gate_receipt": str(RUNNER_GATE_RECEIPT),
            "downstream_output_slots_absent_at_promotion": downstream_absence,
            "canonical_publication_order": (
                "provisional_completion_receipt_then_canonical_no_replace_commit"
            ),
            "canonical_alone_grants_downstream_authority": False,
            "write_ahead_receipt_marks_transaction_incomplete_until": (
                "detached_completion_signature_plus_live_canonical_verification"
            ),
        },
        "signature_contract": exit_attestation_signature_contract(
            path=PROMOTE_EXIT_ATTESTATION,
            schema=PROMOTE_EXIT_ATTESTATION_SCHEMA,
            signature=COMPLETION_SIGNATURE,
            public_key=inputs["trusted_keys"]["completion_public_key"],
            private_key_pre_signature=completion_private,
            legacy_payload_signature=LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
        ),
        "checks": dict(COMPLETION_CHECKS),
        "pass": True,
        "pass_semantics": (
            "completion_payload_only; parent promote exit-attestation publication and "
            "detached completion-key signature remain mandatory, followed by continuation "
            "verification, runner-gate evidence, and the signed downstream terminal receipt"
        ),
    }


def attest_promote_child_output(
    authority_module: Any,
    reader: Any,
    precommit: Callable[[], None],
    producer_session: Mapping[str, Any],
    child_wait: Mapping[str, Any],
) -> dict[str, Any]:
    validate_producer_session(producer_session, expected_mode="--promote")
    validate_child_wait(
        child_wait,
        expected_pid=producer_session["child_pid"],
        require_success=True,
    )
    for path, label in (
        (PROMOTE_EXIT_ATTESTATION, "Promote exit attestation"),
        (COMPLETION_SIGNATURE, "Completion signature"),
        (LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE, "Legacy authorization payload signature"),
        (LEGACY_COMPLETION_PAYLOAD_SIGNATURE, "Legacy completion payload signature"),
        (PROMOTION_INCIDENT, "Promotion incident"),
        (VERIFICATION_RECEIPT, "Continuation verification"),
        (RUNNER_GATE_RECEIPT, "Runner gate receipt"),
        (RUNNER_GATE_LOG, "Runner gate log"),
    ):
        require_path_absent(path, label)
    downstream_absence = require_downstream_absent()
    inputs = collect_reviewed_inputs(reader)
    reviews = validate_reviews(
        reader, inputs["reviewed_sources"], inputs["trusted_keys"]
    )
    core = validate_core_evidence(authority_module, reader)
    authorization = validate_signed_authorization(
        authority_module,
        reader,
        inputs=inputs,
        reviews=reviews,
        core=core,
    )
    _, completion_private = open_private_key_metadata(
        reader, COMPLETION_PRIVATE_KEY, expected_mode="0o400"
    )
    authorization_public = {
        key: value for key, value in authorization.items() if key != "bound_fd"
    }
    _, completion_data, completion_record = reopen_full_stat_artifact(
        reader, COMPLETION_RECEIPT
    )
    _, _, canonical_record = reopen_full_stat_artifact(
        reader,
        CANONICAL_RECEIPT,
        expected_sha256=EXPECTED_SOURCE_SHA256,
        expected_size=EXPECTED_SOURCE_SIZE,
    )
    completion_payload = load_json(
        completion_data, "parent-reopened promotion completion receipt"
    )
    completion_created_at = completion_payload.get("created_at_utc")
    require_utc_timestamp(
        completion_created_at, "parent-reopened completion created_at_utc"
    )
    if not strict_json_equal(
        completion_payload.get("producer_session"), dict(producer_session)
    ):
        raise PromotionError("Completion producer session differs in the parent")
    canonical_status, _ = core["builder"].load_tumor_ref_source_identity_attestation(
        CANONICAL_RECEIPT,
        TUMOR_REF_MANIFEST,
        core["manifest_payload"],
    )
    if (
        canonical_status.get("release_gate_pass") is not True
        or canonical_status.get("publishable_task_b_release") is not True
    ):
        raise PromotionError("Parent final-builder replay rejected the canonical receipt")
    canonical_basic = {
        key: canonical_record[key]
        for key in ("path", "size_bytes", "sha256", "mode")
    }
    expected_payload = expected_completion_payload(
        created_at_utc=completion_created_at,
        producer_session=producer_session,
        authorization=authorization_public,
        inputs=inputs,
        core=core,
        reviews=reviews,
        canonical_receipt=canonical_basic,
        canonical_path_status=canonical_status,
        completion_private=completion_private,
        downstream_absence=downstream_absence,
    )
    if not strict_json_equal(completion_payload, expected_payload):
        raise PromotionError("Parent rejected the child completion payload")
    _, prepare_attestation_data, prepare_attestation_record = reopen_full_stat_artifact(
        reader,
        PREPARE_EXIT_ATTESTATION,
        expected_basic={
            key: authorization["producer_exit_attestation"][key]
            for key in ("path", "size_bytes", "sha256", "mode")
        },
    )
    _, _, prepare_signature_record = reopen_full_stat_artifact(
        reader,
        AUTHORIZATION_SIGNATURE,
        expected_basic={
            key: authorization["signature"][key]
            for key in ("path", "size_bytes", "sha256", "mode")
        },
    )
    prepare_attestation_payload = load_json(
        prepare_attestation_data, "parent-reopened prepare exit attestation"
    )
    prepare_created_at = prepare_attestation_payload.get("created_at_utc")
    if require_utc_timestamp(
        completion_created_at, "completion created_at_utc"
    ) < require_utc_timestamp(
        prepare_created_at, "prepare exit attestation created_at_utc"
    ):
        raise PromotionError("Completion receipt predates the signed prepare attestation")
    producer_source, python_runtime = reopen_producer_identities(reader, inputs, core)
    attestation_created_at = now_utc()
    if require_utc_timestamp(
        attestation_created_at, "promote exit attestation created_at_utc"
    ) < require_utc_timestamp(completion_created_at, "completion created_at_utc"):
        raise PromotionError("Promote exit attestation timestamp ordering failed")
    attestation = build_promote_exit_attestation(
        created_at_utc=attestation_created_at,
        producer_session=producer_session,
        child_wait=child_wait,
        completion_record=completion_record,
        canonical_record=canonical_record,
        prepare_attestation_record=prepare_attestation_record,
        prepare_signature_record=prepare_signature_record,
        producer_source=producer_source,
        python_runtime=python_runtime,
        public_key=inputs["trusted_keys"]["completion_public_key"],
        private_key_pre_signature=completion_private,
    )
    if set(attestation) != EXIT_ATTESTATION_KEYS:
        raise PromotionError("Promote exit attestation builder emitted unexpected keys")
    reader.require_paths_still_bound()
    reader.retain_until_process_exit()

    def attestation_precommit() -> None:
        require_path_absent(
            LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
            "Legacy authorization payload signature",
        )
        require_path_absent(
            LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
            "Legacy completion payload signature",
        )
        reader.require_paths_still_bound()
        precommit()

    attestation_record = write_exclusive_readonly(
        PROMOTE_EXIT_ATTESTATION,
        attestation,
        precommit=attestation_precommit,
    )
    return {
        "completion_receipt": completion_record,
        "canonical_receipt": canonical_record,
        "promote_exit_attestation": attestation_record,
        "completion_attestation_signature_pending": str(COMPLETION_SIGNATURE),
        "downstream_started": False,
        "pass": True,
    }


def record_partial_publication_incident(error: BaseException) -> None:
    if (
        not (
            os.path.lexists(CANONICAL_RECEIPT)
            or os.path.lexists(COMPLETION_RECEIPT)
        )
        or os.path.lexists(PROMOTION_INCIDENT)
    ):
        return
    canonical: dict[str, Any]
    try:
        if CANONICAL_RECEIPT.is_symlink() or not CANONICAL_RECEIPT.is_file():
            canonical = {"path": str(CANONICAL_RECEIPT), "valid_regular_file": False}
        else:
            data = CANONICAL_RECEIPT.read_bytes()
            opened = CANONICAL_RECEIPT.stat()
            canonical = {
                "path": str(CANONICAL_RECEIPT.resolve(strict=True)),
                "size_bytes": opened.st_size,
                "sha256": sha256_bytes(data),
                "mode": oct(stat.S_IMODE(opened.st_mode)),
                "device": opened.st_dev,
                "inode": opened.st_ino,
                "link_count": opened.st_nlink,
            }
    except OSError as inspection_error:
        canonical = {
            "path": str(CANONICAL_RECEIPT),
            "inspection_error": type(inspection_error).__name__,
        }
    completion: dict[str, Any] = {
        "path": str(COMPLETION_RECEIPT),
        "exists": os.path.lexists(COMPLETION_RECEIPT),
    }
    if completion["exists"]:
        try:
            if COMPLETION_RECEIPT.is_symlink() or not COMPLETION_RECEIPT.is_file():
                completion["valid_regular_file"] = False
            else:
                completion_data = COMPLETION_RECEIPT.read_bytes()
                completion_stat = COMPLETION_RECEIPT.stat()
                completion.update(
                    {
                        "size_bytes": completion_stat.st_size,
                        "sha256": sha256_bytes(completion_data),
                        "mode": oct(stat.S_IMODE(completion_stat.st_mode)),
                        "device": completion_stat.st_dev,
                        "inode": completion_stat.st_ino,
                        "link_count": completion_stat.st_nlink,
                    }
                )
        except OSError as inspection_error:
            completion["inspection_error"] = type(inspection_error).__name__
    incident = {
        "schema_name": "intersubmod.tumor_ref_source_receipt_promotion_incident",
        "schema_version": "3.0.0",
        "created_at_utc": now_utc(),
        "status": "PROMOTION_ARTIFACT_PUBLISHED_WITH_INCOMPLETE_TRANSACTION",
        "error_class": type(error).__name__,
        "error": str(error),
        "canonical_receipt": canonical,
        "completion_receipt": completion,
        "authorization": str(AUTHORIZATION),
        "required_action": (
            "Do not run downstream. Preserve and archive this incident and canonical artifact "
            "before any separately reviewed recovery."
        ),
        "downstream_permitted": False,
        "pass": False,
    }
    write_exclusive_readonly(PROMOTION_INCIDENT, incident)


def execute_phase_in_forked_child(
    mode: str, producer_session: Mapping[str, Any]
) -> None:
    require_clean_invocation(mode)
    validate_producer_session(producer_session, expected_mode=mode)
    if producer_session["child_pid"] != os.getpid():
        raise PromotionError("Producer session is not executing in its declared child")
    observed_parent_pid, observed_child_start = process_identity(os.getpid())
    if (
        observed_parent_pid != producer_session["parent_pid"]
        or observed_child_start != producer_session["child_start_time_ticks"]
    ):
        raise PromotionError("Forked producer child identity changed before execution")
    bootstrap_fd, bootstrap_stat, authority_module = bootstrap_authority_module()
    reader = authority_module.BoundArtifactReader()
    reader.__enter__()
    try:
        precommit = lambda: require_bootstrap_stable(bootstrap_fd, bootstrap_stat)
        if mode == "--prepare-authorization":
            prepare_authorization(
                authority_module,
                reader,
                precommit,
                producer_session,
            )
        else:
            promote(
                authority_module,
                reader,
                precommit,
                producer_session,
            )
        require_active_cache_empty()
    except BaseException as error:
        if mode == "--promote":
            record_partial_publication_incident(error)
        raise
    finally:
        reader.__exit__(None, None, None)
        os.close(bootstrap_fd)


def attest_phase_in_parent(
    mode: str,
    producer_session: Mapping[str, Any],
    child_wait: Mapping[str, Any],
) -> dict[str, Any]:
    validate_producer_session(producer_session, expected_mode=mode)
    validate_child_wait(
        child_wait,
        expected_pid=producer_session["child_pid"],
        require_success=True,
    )
    if producer_session["parent_pid"] != os.getpid():
        raise PromotionError("Exit attestation is not executing in the declared parent")
    _, observed_parent_start = process_identity(os.getpid())
    if observed_parent_start != producer_session["parent_start_time_ticks"]:
        raise PromotionError("Exit-attestation parent identity changed")
    bootstrap_fd, bootstrap_stat, authority_module = bootstrap_authority_module()
    reader = authority_module.BoundArtifactReader()
    reader.__enter__()
    try:
        precommit = lambda: require_bootstrap_stable(bootstrap_fd, bootstrap_stat)
        result = (
            attest_prepare_child_output(
                authority_module,
                reader,
                precommit,
                producer_session,
                child_wait,
            )
            if mode == "--prepare-authorization"
            else attest_promote_child_output(
                authority_module,
                reader,
                precommit,
                producer_session,
                child_wait,
            )
        )
        require_active_cache_empty()
        return result
    finally:
        reader.__exit__(None, None, None)
        os.close(bootstrap_fd)


def main() -> int:
    if len(sys.argv) != 2 or sys.argv[1] not in {"--prepare-authorization", "--promote"}:
        raise PromotionError("Expected exactly --prepare-authorization or --promote")
    mode = sys.argv[1]
    require_clean_invocation(mode)
    producer_session, child_wait = fork_no_exec_and_wait(
        mode,
        lambda session: execute_phase_in_forked_child(mode, session),
    )
    try:
        validate_child_wait(
            child_wait,
            expected_pid=producer_session["child_pid"],
            require_success=True,
        )
    except BaseException as error:
        if mode == "--promote":
            record_partial_publication_incident(error)
        raise
    try:
        attest_phase_in_parent(mode, producer_session, child_wait)
    except BaseException as error:
        if mode == "--promote":
            record_partial_publication_incident(error)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
