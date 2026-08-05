#!/usr/bin/env python3
"""Read-only source/runtime probe for transition-relation recovery v30."""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
from typing import Any


sys.dont_write_bytecode = True

REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
RESULT_ROOT = TOPIC_ROOT / "results"
WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")
PYTHON_TARGET = Path(
    "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python3.11"
)
PILOT_VERIFIER = AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_v3_schema_recovery.py"
EXPECTED_REGRESSION_SUMMARY = "733 passed"

EXPECTED_SOURCES = {
    "authority_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v30.py",
        None,
        None,
        None,
    ),
    "continuation_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v30.py",
        None,
        None,
        "_require_executing_verifier_matches_bound_source",
    ),
    "runner_gate_replay": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v30.py",
        None,
        None,
        "require_executing_replayer_matches_bound_source",
    ),
    "downstream_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v30.py",
        None,
        None,
        "verify_executing_source",
    ),
    "regression_tests": (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v30.py",
        None,
        None,
        None,
    ),
    "ceremony_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v30.py",
        None,
        None,
        None,
    ),
    "four_role_bootstrap_builder": (
        AUDIT_ROOT / "bootstrap_v30_four_role_keys.py",
        63_736,
        "e971d5db502ebeadca3df28c71be42875de957c16dea2fbc2caeb6579bca8f4a",
        None,
    ),
    "recovery_final_dataset_builder": (
        AUDIT_ROOT / "build_all_ssnv_final_report_dataset_schema_recovery_v29.py",
        351_388,
        "cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f",
        None,
    ),
    "recovery_result_report_finalizer": (
        AUDIT_ROOT / "finalize_task_b_result_release_schema_recovery_v30.py",
        33_422,
        "b9645cf4f57653f078357421592cef19a168db8389ee68ae90898b9c8c63d318",
        None,
    ),
    "recovery_report_builder": (
        AUDIT_ROOT / "build_all_ssnv_report_artifact_schema_recovery_v29.py",
        238_719,
        "fe19be151bfd72978f87f4a003e8fb8732f54c21e49f702ea5c575b92f4ee9ae",
        None,
    ),
    "readonly_probe": (
        AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v30.py",
        None,
        None,
        None,
    ),
}

FORBIDDEN_OUTPUT_SLOTS = (
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
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v2.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v2.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v2.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v2.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v2.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v2.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v2.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v3.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v3.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v3.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v3.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v3.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v3.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v3.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v4.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v4.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v4.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v4.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v4.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v4.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v4.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v5.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v5.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v5.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v5.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v5.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v5.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v5.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v6.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v6.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v6.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v6.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v6.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v6.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v6.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v7.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v7.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v7.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v7.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v7.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v7.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v7.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v8.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v8.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v8.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v8.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v8.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v8.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v8.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v9.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v9.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v9.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v9.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v9.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v9.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v9.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v10.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v10.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v10.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v10.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v10.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v10.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v10.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v11.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v11.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v11.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v11.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v11.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v11.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v11.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v12.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v12.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v12.log",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v12.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v12.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v12.json.ed25519.sig",
    RESULT_ROOT
    / "m2v5_downstream_continuation_supervisor_success.recovery.v12.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v12.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v13.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v13.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v13.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v13.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v13.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v13.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v13.json.ed25519.sig",
    RESULT_ROOT
    / "m2v5_downstream_continuation_supervisor_success.recovery.v13.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v13.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v15.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v15.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v15.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v15.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v15.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v15.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v15.json.ed25519.sig",
    RESULT_ROOT
    / "m2v5_downstream_continuation_supervisor_success.recovery.v15.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v15.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v15.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v15.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v15.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v15.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v15.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v15.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v15.json.ed25519.sig",
    RESULT_ROOT
    / "m2v5_downstream_continuation_supervisor_success.recovery.v15.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v15.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v29.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v29.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v29.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v29.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v29.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v29.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v29.json.ed25519.sig",
    RESULT_ROOT
    / "m2v5_downstream_continuation_supervisor_success.recovery.v29.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v29.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_incident.v3.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_authorization.v3.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion.v3.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v2.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v2.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v3.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v4.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v5.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v6.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v7.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v8.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v9.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v10.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v11.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v12.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v13.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v15.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v15.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v29.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v3.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v3.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v4.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v4.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v5.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v5.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v6.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v6.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v7.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v7.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v8.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v8.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v9.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v9.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v10.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v10.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v11.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v11.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v12.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v12.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v13.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v13.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v15.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v15.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v15.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v15.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v29.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v29.json.ed25519.sig",
    WORKSPACE_ROOT / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v6",
    WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v6",
    WORKSPACE_ROOT
    / "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v6",
    WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v6",
    WORKSPACE_ROOT / "all_ssnv_final_report_dataset_v5_m2v5_source_attested",
    WORKSPACE_ROOT / "all_ssnv_final_report_v5_m2v5_source_attested",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_mendel.v2.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_nash.v2.json",
    TOPIC_ROOT
    / "reviews"
    / "20260722_tumor_ref_schema_recovery_external_claude_opus.v2.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_mendel.v3.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_nash.v3.json",
    TOPIC_ROOT
    / "reviews"
    / "20260722_tumor_ref_schema_recovery_external_claude_opus.v3.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_mendel.v4.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_nash.v4.json",
    TOPIC_ROOT
    / "reviews"
    / "20260722_tumor_ref_schema_recovery_external_claude_opus.v4.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_mendel.v5.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_nash.v5.json",
    TOPIC_ROOT
    / "reviews"
    / "20260722_tumor_ref_schema_recovery_external_claude_opus.v5.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_mendel.v6.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_nash.v6.json",
    TOPIC_ROOT
    / "reviews"
    / "20260722_tumor_ref_schema_recovery_external_claude_opus.v6.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_mendel.v7.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_nash.v7.json",
    TOPIC_ROOT
    / "reviews"
    / "20260722_tumor_ref_schema_recovery_external_claude_opus.v7.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_mendel.v8.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_nash.v8.json",
    TOPIC_ROOT
    / "reviews"
    / "20260722_tumor_ref_schema_recovery_external_claude_opus.v8.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_mendel.v9.json",
    TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_nash.v9.json",
    TOPIC_ROOT
    / "reviews"
    / "20260722_tumor_ref_schema_recovery_external_claude_opus.v9.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v10.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v10.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v10.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v11.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v11.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v11.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v12.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v12.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v12.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v13.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v13.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v13.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v15.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v15.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v15.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v15.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v15.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v15.json",
)

FAILED_V16_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v16.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v16.json",
    RESULT_ROOT
    / "tumor_ref_promotion_schema_recovery_authority.v16.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v16.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v16.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v16.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v16.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v16.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v16.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v16.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v16.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v16.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v16.json.ed25519.sig",
    RESULT_ROOT
    / "m2v5_downstream_continuation_supervisor_success.recovery.v16.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v16.json",
)
FAILED_V17_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v17.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v17.json",
    RESULT_ROOT
    / "tumor_ref_promotion_schema_recovery_authority.v17.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v17.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v17.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v17.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v17.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v17.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v17.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v17.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v17.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v17.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v17.json.ed25519.sig",
    RESULT_ROOT
    / "m2v5_downstream_continuation_supervisor_success.recovery.v17.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v17.json",
)
FAILED_V18_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v18.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v18.json",
    RESULT_ROOT
    / "tumor_ref_promotion_schema_recovery_authority.v18.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v18.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v18.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v18.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v18.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v18.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v18.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v18.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v18.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v18.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v18.json.ed25519.sig",
    RESULT_ROOT
    / "m2v5_downstream_continuation_supervisor_success.recovery.v18.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v18.json",
)
FAILED_V19_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v19.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v19.json",
    RESULT_ROOT
    / "tumor_ref_promotion_schema_recovery_authority.v19.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v19.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v19.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v19.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v19.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v19.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v19.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v19.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v19.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v19.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v19.json.ed25519.sig",
    RESULT_ROOT
    / "m2v5_downstream_continuation_supervisor_success.recovery.v19.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v19.json",
)
REJECTED_V20_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v20.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v20.json",
    RESULT_ROOT
    / "tumor_ref_promotion_schema_recovery_authority.v20.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v20.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v20.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v20.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v20.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v20.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v20.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v20.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v20.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v20.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v20.json.ed25519.sig",
    RESULT_ROOT
    / "m2v5_downstream_continuation_supervisor_success.recovery.v20.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v20.json",
)
FAILED_V21_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v21.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v21.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v21.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v21.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v21.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v21.json",
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
REJECTED_V22_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v22.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v22.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v22.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v22.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v22.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v22.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v22.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v22.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v22.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v22.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v22.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v22.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v22.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v22.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v22.json",
)
REJECTED_V23_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v23.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v23.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v23.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v23.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v23.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v23.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v23.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v23.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v23.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v23.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v23.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v23.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v23.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v23.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v23.json",
)
FAILED_V14_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v14.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v14.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v14.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v14.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v14.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v14.json",
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
REJECTED_V24_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v24.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v24.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v24.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v24.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v24.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v24.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v24.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v24.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v24.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v24.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v24.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v24.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v24.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v24.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v24.json",
)
FAILED_V25_OUTPUT_SLOTS = (
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
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_mendel.v25.json",
    TOPIC_ROOT / "reviews" / "20260723_tumor_ref_schema_recovery_nash.v25.json",
    TOPIC_ROOT
    / "reviews"
    / "20260723_tumor_ref_schema_recovery_external_claude_opus.v25.json",
)
FAILED_V26_OUTPUT_SLOTS = (
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
    TOPIC_ROOT / "reviews" / "20260724_tumor_ref_schema_recovery_mendel.v26.json",
    TOPIC_ROOT / "reviews" / "20260724_tumor_ref_schema_recovery_nash.v26.json",
    TOPIC_ROOT
    / "reviews"
    / "20260724_tumor_ref_schema_recovery_external_claude_opus.v26.json",
)
REJECTED_V27_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v27.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v27.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v27.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260724_tumor_ref_schema_recovery_mendel.v27.json",
    TOPIC_ROOT / "reviews" / "20260724_tumor_ref_schema_recovery_nash.v27.json",
    TOPIC_ROOT
    / "reviews"
    / "20260724_tumor_ref_schema_recovery_external_claude_opus.v27.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v27.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v27.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v27.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v27.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v27.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v27.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v27.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v27.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v27.json",
)
FAILED_V28_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v28.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v28.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v28.json.ed25519.sig",
    TOPIC_ROOT / "reviews" / "20260724_tumor_ref_schema_recovery_mendel.v28.json",
    TOPIC_ROOT / "reviews" / "20260724_tumor_ref_schema_recovery_nash.v28.json",
    TOPIC_ROOT
    / "reviews"
    / "20260724_tumor_ref_schema_recovery_external_claude_opus.v28.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v28.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v28.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v28.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v28.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v28.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v28.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v28.json",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json.ed25519.sig",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json.ed25519.sig",
)
FAILED_V29_REVIEW_OUTPUT_SLOTS = (
    TOPIC_ROOT / "reviews" / "20260724_tumor_ref_schema_recovery_mendel.v29.json",
    TOPIC_ROOT / "reviews" / "20260724_tumor_ref_schema_recovery_nash.v29.json",
    TOPIC_ROOT
    / "reviews"
    / "20260724_tumor_ref_schema_recovery_external_claude_opus.v29.json",
)
CURRENT_V30_OUTPUT_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v30.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v30.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v30.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v30.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v30.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v30.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v30.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v30.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v30.json",
)
FORBIDDEN_OUTPUT_SLOTS = tuple(
    dict.fromkeys(
        (
            *FORBIDDEN_OUTPUT_SLOTS,
            *FAILED_V16_OUTPUT_SLOTS,
            *FAILED_V17_OUTPUT_SLOTS,
            *FAILED_V18_OUTPUT_SLOTS,
            *FAILED_V19_OUTPUT_SLOTS,
            *REJECTED_V20_OUTPUT_SLOTS,
            *FAILED_V21_OUTPUT_SLOTS,
            *FAILED_V21_ORIGINAL_OBSERVED_OUTPUT_SLOTS,
            *REJECTED_V22_OUTPUT_SLOTS,
            *REJECTED_V23_OUTPUT_SLOTS,
            *FAILED_V14_OUTPUT_SLOTS,
            *REJECTED_V24_OUTPUT_SLOTS,
            *FAILED_V25_OUTPUT_SLOTS,
            *FAILED_V26_OUTPUT_SLOTS,
            *REJECTED_V27_OUTPUT_SLOTS,
            *FAILED_V28_OUTPUT_SLOTS,
            *FAILED_V29_REVIEW_OUTPUT_SLOTS,
            *CURRENT_V30_OUTPUT_SLOTS,
        )
    )
)

REVIEW_EVIDENCE_PATHS = (
    TOPIC_ROOT / "reviews" / "20260724_tumor_ref_schema_recovery_mendel.v30.json",
    TOPIC_ROOT / "reviews" / "20260724_tumor_ref_schema_recovery_nash.v30.json",
    TOPIC_ROOT
    / "reviews"
    / "20260724_tumor_ref_schema_recovery_external_claude_opus.v30.json",
)

CEREMONY_STAGING_PATTERNS = tuple(
    f".tumor_ref_promotion_schema_recovery_authority.v{generation}.bundle.staging.*"
    for generation in range(2, 31)
)

FAILED_V12_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v12_signed_authority_c_python_executable_alias_relation_mismatch"
)
FAILED_V14_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v14_signed_authority_r_metadata_only_relation_schema_mismatch"
)
REJECTED_V15_ROOT = (
    AUDIT_ROOT
    / "rejected_pre_authority_reviews"
    / "20260723_v15_candidate_runtime_schema_and_provenance_findings"
)
FAILED_V16_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v16_signed_authority_c_legacy_eight_key_stat_schema_mismatch"
)
FAILED_V17_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v17_signed_authority_c_metadata_enrichment_relation_conflict"
)
FAILED_V18_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v18_signed_authority_c_tumor_ref_summary_alias_noncanonical"
)
FAILED_V19_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v19_signed_authority_c_mode000_inotify_permission_denied"
)
FAILED_V25_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260724_v25_signed_authority_c_tumor_ref_pre_audit_path_relation_mismatch"
)
FAILED_V26_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260724_v26_signed_authority_c_signed_promotion_runtime_role_projection_mismatch"
)
REJECTED_V27_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260724_v27_pre_authority_runtime_contract_and_review_transport_rejection"
)
FAILED_V28_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260724_v28_signed_dataset_c_report_metric_key_order_mismatch"
)
REJECTED_V29_ROUND1_ROOT = (
    AUDIT_ROOT
    / "rejected_pre_authority_reviews"
    / "20260724_v29_round1_failed_v28_terminal_binding_gap"
)
FAILED_V29_ROOT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260724_v29_signed_authority_v_pass_r_archived_key_live_path_mismatch"
)
FAILED_V29_KEY_ARCHIVES = {
    "authority_v29": Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
        "archive/20260724_v29_signed_authority_v_pass_r_archived_key_live_path_failure_01"
    ),
    "terminal_v19": Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
        "archive/20260724_m2v5_terminal_v19_unused_after_signed_v29_r_failure_01"
    ),
    "result_v6": Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
        "archive/20260724_all_ssnv_result_v6_unused_v29_r_failure_01"
    ),
    "report_v6": Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
        "archive/20260724_all_ssnv_report_v6_unused_v29_r_failure_01"
    ),
}
V30_BOOTSTRAP_PREPARE = (
    AUDIT_ROOT / "20260724_v30_four_role_key_bootstrap_prepare.json"
)
V30_BOOTSTRAP_PROGRESS = (
    AUDIT_ROOT / "20260724_v30_four_role_key_bootstrap_progress.jsonl"
)
V30_BOOTSTRAP_RECEIPT = (
    AUDIT_ROOT / "20260724_v30_four_role_key_bootstrap_receipt.json"
)
V30_BOOTSTRAP_SUCCESS = (
    AUDIT_ROOT / "20260724_v30_four_role_key_bootstrap_success.json"
)

REQUIRED_PRIOR_INPUTS = {
    "rejected_v2_generation": (
        AUDIT_ROOT / "20260722_recovery_v2_formal_rejection_and_key_retirement.v1.json",
        4_482,
        "f0c64e46ab34872149ca0cae830e79f48d818b591f8100160ac0cce44ed5ba84",
    ),
    "rejected_v3_generation": (
        AUDIT_ROOT / "20260722_recovery_v3_formal_rejection_and_key_retirement.v1.json",
        4_473,
        "c31009328e2130449422e01a5fd766446f6673e6f2625be0f4f380e3f41e4ef5",
    ),
    "rejected_v4_generation": (
        AUDIT_ROOT / "20260722_recovery_v4_formal_rejection_and_key_retirement.v1.json",
        4_676,
        "a1044ae1a0580b9a6587e30ddfeea22afeac84128cae0d8a28ee3a64619b9fb3",
    ),
    "rejected_v5_generation": (
        AUDIT_ROOT / "20260722_recovery_v5_formal_rejection_and_key_retirement.v1.json",
        4_322,
        "0495af623ba463b822f4f823ce28d49745c11c056ae391a36deeda06e0e78047",
    ),
    "rejected_v6_generation": (
        AUDIT_ROOT / "20260722_recovery_v6_formal_rejection_and_key_retirement.v1.json",
        6_506,
        "2a8a4a779d1e6df8a31e90adb64ae37e89d075ce78a3f0cba8942aea4359c9ab",
    ),
    "rejected_v7_generation": (
        AUDIT_ROOT / "20260722_recovery_v7_formal_rejection_and_key_retirement.v1.json",
        5_732,
        "07858f6f990cb9d72b9fca421d3e064b68d8641ec056c96a25c8adacb98c5956",
    ),
    "rejected_v8_generation": (
        AUDIT_ROOT / "20260722_recovery_v8_formal_rejection_and_key_retirement.v1.json",
        5_695,
        "fe1724c88f2e168322fa16858db7450afd7145d8a0cbda097f019bfa052f6cd9",
    ),
    "rejected_v13_pre_authority_superseded_evidence": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v13_candidate_nash_medium_findings"
        / "rejection_evidence.json",
        5_830,
        "a6dad5b9c6ec7d26a6390751b8296ce0d7eea1f9dfc419633458a88a36dc05fa",
    ),
    "rejected_v13_pre_authority_evidence": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v13_candidate_nash_medium_findings"
        / "rejection_evidence.v2.json",
        6_694,
        "b42de4167900069fbe563244b9058a8291fe1baf097fb0a5f5e77dee74ee02eb",
    ),
    "rejected_v13_key_archive_ledger": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v13_candidate_nash_medium_findings"
        / "key_archive_ledger.v1.jsonl",
        2_783,
        "0f4b9e0d6d48187cd7b63bafa46472efb83d8f1f3d790d4503e8423102799321",
    ),
    "authority_v1": (
        RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v1.json",
        15_076,
        "af7c2cf68b7edbfe3883cd65b273213985bb9058e1799984860c3a09e1725fda",
    ),
    "authority_v1_signature": (
        RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v1.json.ed25519.sig",
        64,
        "ab2026b3ae7fce9c2e602e06dd95760f1d386a5ec90525471600f50ba5b60243",
    ),
    "review_v1_mendel": (
        TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_mendel.v1.json",
        4_788,
        "7341d17279611b2c1a6a97815b68bd2f73fa4fc9d9a56e6f6b802cb1efba9502",
    ),
    "review_v1_nash": (
        TOPIC_ROOT / "reviews" / "20260722_tumor_ref_schema_recovery_nash.v1.json",
        4_786,
        "3631d4a9c087a64cddf80e9184949fb5ac0b7cf5fd6bb780248c08dc38c1c1e2",
    ),
    "review_v1_external": (
        TOPIC_ROOT
        / "reviews"
        / "20260722_tumor_ref_schema_recovery_external_claude_opus.v1.json",
        4_802,
        "f51fe4042a40244ae866b78a64438fe77bda19642d8a93a7b4e9f15caad0c556",
    ),
    "verification_receipt_v1": (
        RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v1.json",
        13_422,
        "b1f9b22d570ca68f91fd3039c2d2f9c1956c1bd769f4b30c983ddf6a122a13a5",
    ),
    "runner_failure_evidence_v1": (
        AUDIT_ROOT / "20260722_rrec_v1_historical_transition_live_binding_failure.v1.json",
        4_686,
        "7af8f793c1a8ca534445e598608ca3edf0d01f237d4130f43537b6f8c3fb7cb0",
    ),
    "public_key_v1": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260722_v1/ed25519_public.pem"
        ),
        113,
        "1614f15414b3719187f2fcf281794f206e6f12b7c32fd9ab7d532ef167ba1f34",
    ),
    "failed_v9_evidence": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v9_signed_authority_runtime_17_23_slot_mismatch"
        / "failure_evidence.json",
        4_741,
        "024d58b94c61c7b94ff3e842d47039fb228cb69e09e29806a000caef7f041bcc",
    ),
    "failed_v9_authority": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v9_signed_authority_runtime_17_23_slot_mismatch"
        / "tumor_ref_promotion_schema_recovery_authority.v9.bundle"
        / "authority.json",
        33_050,
        "28ab25f3866a37ee1df42dba20c9bae0633311b0761e0cfabaff442cf079e58a",
    ),
    "failed_v9_signature": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v9_signed_authority_runtime_17_23_slot_mismatch"
        / "tumor_ref_promotion_schema_recovery_authority.v9.bundle"
        / "authority.ed25519.sig",
        64,
        "bdcf725c65758d8d89ddd5706ffd27f766a1c7c0fe14af8f205c620a040f9ece",
    ),
    "failed_v9_commit": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v9_signed_authority_runtime_17_23_slot_mismatch"
        / "tumor_ref_promotion_schema_recovery_authority.v9.bundle"
        / "commit.json",
        802,
        "11a17c5942178423517aac9055398801c31d966e8e23c9d2088b565dad498a7c",
    ),
    "failed_v9_review_mendel": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v9_signed_authority_runtime_17_23_slot_mismatch"
        / "20260722_tumor_ref_schema_recovery_mendel.v9.json",
        2_320,
        "73db88ac230d03316f7065e6d5ce11fbf080dc90d65ecc269161792a53967b3e",
    ),
    "failed_v9_review_nash": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v9_signed_authority_runtime_17_23_slot_mismatch"
        / "20260722_tumor_ref_schema_recovery_nash.v9.json",
        2_293,
        "af22f7d609e89fd7a4a4e1a729bbd184d5d7b2f42553bfee08ca344a6b985564",
    ),
    "failed_v9_review_external": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v9_signed_authority_runtime_17_23_slot_mismatch"
        / "20260722_tumor_ref_schema_recovery_external_claude_opus.v9.json",
        3_064,
        "aca38fcfcd6dfee6dbe6312c7997d48f2c4e5eea8c988b988e5acc65a2b2c864",
    ),
    "failed_v9_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260722_v9/ed25519_public.pem"
        ),
        113,
        "1aaa9f2fc71fa854f042fe8a8c6e99cc0a1d72d4224f4e85953e5d747c348b7b",
    ),
    "failed_v10_evidence": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
        / "failure_evidence.json",
        6_190,
        "25768002648008ad0ae6c459731f33ece68d0c74f484b4084de0fa89a6227701",
    ),
    "failed_v10_authority": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
        / "tumor_ref_promotion_schema_recovery_authority.v10.bundle"
        / "authority.json",
        37_868,
        "3e1c3851c9784ec3ace3135283c9a839902b3ce404de7157e96b707aefb9fd11",
    ),
    "failed_v10_signature": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
        / "tumor_ref_promotion_schema_recovery_authority.v10.bundle"
        / "authority.ed25519.sig",
        64,
        "3927b9deee2ab87e21bf7d44c5ca9c622fcc77c87dfeccb98fd8d96959856bba",
    ),
    "failed_v10_commit": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
        / "tumor_ref_promotion_schema_recovery_authority.v10.bundle"
        / "commit.json",
        803,
        "dcd77c2fc958f17ca8d4e1c80e566880bbb7db924c41aa6f5e8a7073e759af2e",
    ),
    "failed_v10_review_mendel": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
        / "20260723_tumor_ref_schema_recovery_mendel.v10.json",
        2_281,
        "a8ec6f3da0c68641f9d477a801a6cdcb7205a299f717fc6262f7b4d8df0c9a36",
    ),
    "failed_v10_review_nash": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
        / "20260723_tumor_ref_schema_recovery_nash.v10.json",
        2_767,
        "ae56105735feb27a00a0dcb17d5a7f1acfbf888b68beab9f480757e78c55cf41",
    ),
    "failed_v10_review_external": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v10.json",
        2_692,
        "c5bb30530e348209d1538ecb30964f9f9edfbbf90b7f63630efd8ae961737680",
    ),
    "failed_v10_verification_receipt": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
        / "tumor_ref_source_receipt_promotion_verification.recovery.v10.json",
        14_015,
        "b005ba27638a290e9fd3c6e1ee07b4a47824b27dc7d9eaa1f7c1f4cfcef16d9c",
    ),
    "failed_v10_replay_receipt": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
        / "m2v5_runner_only_gate_replay.recovery.v10.json",
        63_837,
        "25011cbd0e55fd8aabb07b78907f438570767791185dfa9b818e48045fe50d16",
    ),
    "failed_v10_replay_log": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
        / "m2v5_runner_only_gate_replay.recovery.v10.log",
        15_741,
        "9a76c28ac30e82d8963a5d4b7d4ecfda51744aaf6fc1a7f7e8db6b3e118341aa",
    ),
    "failed_v10_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v10/ed25519_public.pem"
        ),
        113,
        "ac17fe53932a426216b12ed875b1ef84c946e0362f8a34359bff2d23bb51bc4e",
    ),
    "failed_v11_evidence": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
        / "failure_evidence.json",
        9_812,
        "bc76d6026a20857145feaeb188a4cc5f933ed8b262e2f3f7d107d4bde13f7619",
    ),
    "failed_v11_authority": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
        / "tumor_ref_promotion_schema_recovery_authority.v11.bundle"
        / "authority.json",
        44_364,
        "2fe51e9f3a746781324dfb70156d96481d59f0c673b8623f7da28df20cf924ad",
    ),
    "failed_v11_signature": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
        / "tumor_ref_promotion_schema_recovery_authority.v11.bundle"
        / "authority.ed25519.sig",
        64,
        "8e95a60f23cd75bf4b8b39490c5142b7cf3fac04873beffa924f76cdc3fee69c",
    ),
    "failed_v11_commit": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
        / "tumor_ref_promotion_schema_recovery_authority.v11.bundle"
        / "commit.json",
        803,
        "f4900248e51d64e961454b433558fca8bb1d7f074fd257bc03bde73bbef493c8",
    ),
    "failed_v11_review_mendel": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
        / "20260723_tumor_ref_schema_recovery_mendel.v11.json",
        2_713,
        "ade00ba448525616bdc2b5953d3720a44f996d763ea498f594461dcfd2514007",
    ),
    "failed_v11_review_nash": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
        / "20260723_tumor_ref_schema_recovery_nash.v11.json",
        2_534,
        "ee23983fe4fbc5336ef063dc4effaffb08e621b392c7b4fa1f490a67f469c42d",
    ),
    "failed_v11_review_external": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v11.json",
        2_990,
        "340a580ebc849d3b88ddbf3a50de9742b5357c4d6dbadd1d4e310fc11a459be3",
    ),
    "failed_v11_verification_receipt": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
        / "tumor_ref_source_receipt_promotion_verification.recovery.v11.json",
        14_017,
        "011eb6afef7e0e5d20c02e241c90e19b53f7116bfc800801bf503d116f4164c1",
    ),
    "failed_v11_replay_receipt": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
        / "m2v5_runner_only_gate_replay.recovery.v11.json",
        63_842,
        "b3ecee416e260461c19f2130269220cf064889c9850ae54e4bcd94cbf2489d29",
    ),
    "failed_v11_replay_log": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
        / "m2v5_runner_only_gate_replay.recovery.v11.log",
        15_743,
        "46f5c44d1b7974484f9d8becaa8c75b43b2aba1fddad4e2eba8caf5b830b4b1f",
    ),
    "failed_v11_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v11/ed25519_public.pem"
        ),
        113,
        "2c1e3b609cb0f58dd23233d257ef66f3950ee5953dcda18474ce209525180b7b",
    ),
    "failed_v12_evidence": (
        FAILED_V12_ROOT / "failure_evidence.json",
        11_521,
        "a5655e8c19fa0ef52e008bf045aefa09a6be1172a5adcfe8e153e62a0cbac62b",
    ),
    "failed_v12_authority": (
        FAILED_V12_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v12.bundle"
        / "authority.json",
        53_842,
        "df7db2a31b7f7a417363344f85296d7c3d9bf470e084fd1d81e871be082dc3ca",
    ),
    "failed_v12_signature": (
        FAILED_V12_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v12.bundle"
        / "authority.ed25519.sig",
        64,
        "0e2d4a54757e78f396b5fd1f9f7aae7100e556547d35683160e76d925f72a5fc",
    ),
    "failed_v12_commit": (
        FAILED_V12_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v12.bundle"
        / "commit.json",
        803,
        "eb40f8f7bf0e564e4a261dcfd3163d1688378479d31d08012281a309fbaddaa5",
    ),
    "failed_v12_review_mendel": (
        FAILED_V12_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v12.json",
        2_249,
        "b2c22ce612bac10d607d082375b7f53b2d3c231033a68cd6b24f96399d7a9aad",
    ),
    "failed_v12_review_nash": (
        FAILED_V12_ROOT / "20260723_tumor_ref_schema_recovery_nash.v12.json",
        2_241,
        "a3e85b961fa48fa5dbd21c3e4016ec6cdd3fb6394649b5ea9599a63a94cdd585",
    ),
    "failed_v12_review_external": (
        FAILED_V12_ROOT
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v12.json",
        2_562,
        "9ef0c7e1d9e37cb1819ab10faaf43da05bb887c0ea4a78a17edfcbd5734fd009",
    ),
    "failed_v12_verification_receipt": (
        FAILED_V12_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v12.json",
        14_017,
        "b2ffe1da3e99f34e33480ec48873f72199d4e262051833d2d80edcc79a0d9e98",
    ),
    "failed_v12_replay_receipt": (
        FAILED_V12_ROOT / "m2v5_runner_only_gate_replay.recovery.v12.json",
        64_924,
        "4c38a5887b057987f589ef86578f30bbd7629e2b97de11d64d7f0b472217b443",
    ),
    "failed_v12_replay_log": (
        FAILED_V12_ROOT / "m2v5_runner_only_gate_replay.recovery.v12.log",
        15_743,
        "52e2804978801b9986d8a72e155e4210c228e07adcbeb274f7f6f092914b0d01",
    ),
    "failed_v12_source_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v12.py",
        134_395,
        "6096b5823144b65af926e5d438f5480236e4bbc1eccec708bcfa6f92e5dddcc9",
    ),
    "failed_v12_source_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v12.py",
        127_760,
        "42fc979b01f7e8cbe68f2f11321f6a384d3ef769a2e8d4bdf460ba7cf847bde2",
    ),
    "failed_v12_source_replayer": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v12.py",
        132_252,
        "10ca48e72f76cb445e67298860ab290ab66a2521db6587e34a04ca238b9b6b65",
    ),
    "failed_v12_source_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v12.py",
        319_106,
        "ab423dced07620b667745122f911a6c90f83c0a672f0cf945cefc6ce7a265300",
    ),
    "failed_v12_source_probe": (
        AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v12.py",
        45_396,
        "7b4853700b27303fb45da1a26dd5bf16b867050393902a536011a61ada440ca7",
    ),
    "failed_v12_source_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v12.py",
        51_664,
        "b29bff632e360d88dee5496ce488dda0bc6fa2e39c61b27231094f2504c593d9",
    ),
    "failed_v12_source_tests": (
        AUDIT_ROOT / "schema_recovery_tests" / "test_tumor_ref_schema_recovery_v12.py",
        83_583,
        "508530fe569ce0aa48dcfefd45221f2f9bad81f476f331d78deff5f11f6622e7",
    ),
    "failed_v12_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v12/ed25519_public.pem"
        ),
        113,
        "b6c0734543c9608f8af830c89bdde071a36b3cc38ab9d5c53399fd63a2d46d3d",
    ),
    "failed_v14_evidence": (
        FAILED_V14_ROOT / "failure_evidence.json",
        12_657,
        "8e3f5c13be270279a11f4e3e11d708a5815d787e5b0eccaf524885616166af9c",
    ),
    "failed_v14_authority": (
        FAILED_V14_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v14.bundle"
        / "authority.json",
        73_225,
        "d534c226a32da7a5b08701aa6abdf693e6aac7fb30e0e78927854c86ac486c0a",
    ),
    "failed_v14_signature": (
        FAILED_V14_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v14.bundle"
        / "authority.ed25519.sig",
        64,
        "e25f55a711e109997a9c250ad80eb4031b49ee93c9d54a7ad93da121c818a274",
    ),
    "failed_v14_commit": (
        FAILED_V14_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v14.bundle"
        / "commit.json",
        803,
        "b4f901d3218408d2ec47cb38620108859f5d75ef212865d6eca3688d31d58fc2",
    ),
    "failed_v14_verification_receipt": (
        FAILED_V14_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v14.json",
        18_477,
        "579bec69e9c59bf77ca5c6444180d0c5596203368983100e6d32a3d82cc72f3f",
    ),
    "failed_v14_formal_replay_stderr": (
        FAILED_V14_ROOT / "formal_replay_stderr.log",
        1_493,
        "7c5391e3eb128b72736986720fb482ad36995e70deb97b5208b1639cf0fc59a7",
    ),
    "failed_v14_review_mendel": (
        FAILED_V14_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v14.json",
        2_483,
        "c81627a0dc9a9cd339325336ba5c4b403a6e1736b4484bcc2b2eaa05d2783d76",
    ),
    "failed_v14_review_nash": (
        FAILED_V14_ROOT / "20260723_tumor_ref_schema_recovery_nash.v14.json",
        2_450,
        "6579099ecdbab100d64777769fb833279091df9abe8a991d1b4f25d55fb67091",
    ),
    "failed_v14_review_external": (
        FAILED_V14_ROOT
        / "20260723_tumor_ref_schema_recovery_external_claude_opus.v14.json",
        2_636,
        "882aaf0000c0ef5eac89daaaa28da6e70d5bf5671e611df1ac001ebe69008b4e",
    ),
    "failed_v14_source_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v14.py",
        179_003,
        "1d9480b0268fe17c1d6bfb5467671c76a1623d197835821d453ae75a99e9410c",
    ),
    "failed_v14_source_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v14.py",
        127_906,
        "13b159e3f1fcd961fd693331cd8f6f3a3a0862a663cceb648d755f22cb91bd3a",
    ),
    "failed_v14_source_replayer": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v14.py",
        148_758,
        "14a012a3407bbe8d8b0855545f7060763107400bc4206571d2612b11a61adcb7",
    ),
    "failed_v14_source_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v14.py",
        344_261,
        "e254a21b4dc95d6c25085256080c1d66f95dd3b474f47aa07096bd8367ca90c1",
    ),
    "failed_v14_source_probe": (
        AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v14.py",
        56_196,
        "69742abbdb033f7a004e0a780364cb301b96a7ef95019d64392e853e1e521dd8",
    ),
    "failed_v14_source_tests": (
        AUDIT_ROOT / "schema_recovery_tests" / "test_tumor_ref_schema_recovery_v14.py",
        109_952,
        "933ad2162b62cf58da9d924aed2ad081a6b9e1d5b5e896a79eac497b3481e860",
    ),
    "failed_v14_source_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v14.py",
        53_376,
        "5d9bb98d0192ef01cf3e85ab239c519958b9ec5667b9887adeb1949dea50c326",
    ),
    "failed_v14_archive_script": (
        AUDIT_ROOT / "archive_v14_signed_r_metadata_relation_failure.py",
        19_863,
        "181d8eac8f3fea95b6ef70a0be6436eeca42d5a7bbb240eb04cd8496fffd1358",
    ),
    "failed_v14_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v14/ed25519_public.pem"
        ),
        113,
        "91f3b81a0dfab1911b492269dc40ef150eed76bf3b16c4143ba541d16ffdc8a3",
    ),
    "failed_v14_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v4/ed25519_public.pem"
        ),
        113,
        "091785646a7cadf2295f97668838b8279a2bb3b8a798e317dcf1ec6aba33d427",
    ),
    "rejected_v15_evidence": (
        REJECTED_V15_ROOT / "rejection_evidence.json",
        6_290,
        "40b329951a199fdf33e2efa7c73b2d7779d12d5c707d237ed7cc9e1734930297",
    ),
    "rejected_v15_key_archive_ledger": (
        REJECTED_V15_ROOT / "key_archive_ledger.v1.jsonl",
        2_801,
        "6409605672db27706ca1666225a8dae20ae31a61f4893414f2cf4fb284b9f9bc",
    ),
    "rejected_v15_review_mendel": (
        REJECTED_V15_ROOT / "mendel_request_changes.json",
        1_494,
        "6bbcb70233be88a2306abe42fc170af097c3a6a17835c2f97ca6456633e1d4d4",
    ),
    "rejected_v15_review_nash": (
        REJECTED_V15_ROOT / "nash_request_changes.json",
        1_342,
        "07dd41f18969fc1a05a72ea7f1f84c732f2c2a98518cf6a96d45dbcea3528e60",
    ),
    "rejected_v15_review_external": (
        REJECTED_V15_ROOT / "external_claude_opus_approve.json",
        1_497,
        "653ff58b82b41d2dab5ccdd7a40925ee2b6bc418e5acf08e53a7f8f06b8199aa",
    ),
    "rejected_v15_external_envelope": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v15.attempt2.claude_cli.envelope.json",
        18_191,
        "86671481a742a61b899eb2fe59a278d27e754de9c261d9084349fd2a320bd332",
    ),
    "rejected_v15_external_stderr": (
        AUDIT_ROOT
        / "20260723_external_claude_schema_recovery_review_v15.attempt2.claude_cli.stderr.txt",
        302,
        "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a",
    ),
    "rejected_v15_source_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v15.py",
        198_269,
        "e09220da97cffcdf8edcdd44f447fac88292eb2b989543ea774a8547f85fe090",
    ),
    "rejected_v15_source_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v15.py",
        127_906,
        "aa8c5ed29a46ba4a4618e9b40bf7d627df7cf15ffc6daea240ad16aafe3bd043",
    ),
    "rejected_v15_source_replayer": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v15.py",
        150_445,
        "0f0811a9171653eea2f59e81b1e5283f281457f40b6c64960c9f32f34d52a4d9",
    ),
    "rejected_v15_source_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v15.py",
        344_261,
        "df5e93539b7fe33d4b8ccc992656ed5688f265e1287d84176564231c30db989f",
    ),
    "rejected_v15_source_probe": (
        AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v15.py",
        62_434,
        "737b1d31fbbe59c589c5cda275d5d219ce5a76d5163f38b7afef97093ce03bdf",
    ),
    "rejected_v15_source_tests": (
        AUDIT_ROOT / "schema_recovery_tests" / "test_tumor_ref_schema_recovery_v15.py",
        113_042,
        "9e379b4f1ba2f9d27228b9bba5c8ea03ad9f28d6ca94c8a93d8774a0db754bbf",
    ),
    "rejected_v15_source_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v15.py",
        53_376,
        "bfc4de402d3f1d6fac836d179ed3a080b8008c3874f1f48b634bc4b621778582",
    ),
    "rejected_v15_authority_archive_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "archive/20260723_v15_unused_pre_authority_review_rejection_01/"
            "ed25519_public.pem"
        ),
        113,
        "797c9c174b72ca588d6a63b16c4ab1f0b1bf465763a4ee030e367e1e807aaf4d",
    ),
    "rejected_v15_authority_archive_record": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "archive/20260723_v15_unused_pre_authority_review_rejection_01/"
            "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
        ),
        909,
        "b15ea9c3790bd25b864a5120514077515c0fc7f3e6d39f9cab296aab73887b1a",
    ),
    "rejected_v15_terminal_archive_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260723_m2v5_terminal_v5_unused_pre_authority_review_rejection_01/"
            "ed25519_public.pem"
        ),
        113,
        "17a0b67f22b706da120ca46ce37e9104cc53ba2fd524fd37c92b261668e00f84",
    ),
    "rejected_v15_terminal_archive_record": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260723_m2v5_terminal_v5_unused_pre_authority_review_rejection_01/"
            "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
        ),
        959,
        "1ab26be741aa8634d7cadf26ab7bd33e9d9a13eba232f1003a7f4b466abddbfc",
    ),
    "failed_v16_evidence": (
        FAILED_V16_ROOT / "failure_evidence.json",
        17_151,
        "c6ba5ced3fa16b06bb2e4dbd974c8db3e8f4c9f1c9304f1021da061835bf0c25",
    ),
    "failed_v16_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v16/ed25519_public.pem"
        ),
        113,
        "540b64ed3615618efed89637069f772787fdd025acadbbd27b6e334423d2345e",
    ),
    "failed_v16_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v6/ed25519_public.pem"
        ),
        113,
        "066949c0c36be413cd2fb60670e5a2fbc583ab9a3a4264ebf4d3766aba39867f",
    ),
    "failed_v17_evidence": (
        FAILED_V17_ROOT / "failure_evidence.json",
        40_696,
        "b237a8d03b4a483abd47b3239d5d8d13442794835c090dfe4122b60834b025a9",
    ),
    "failed_v17_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v17/ed25519_public.pem"
        ),
        113,
        "ef3413162e1cb9850fb966a1505d71cf185fe79812415b27a93f4d5887154eac",
    ),
    "failed_v17_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v7/ed25519_public.pem"
        ),
        113,
        "4b83e655f1a7a778691e27dc3df2257a230001c702afdf6e703e578c706e0b03",
    ),
    "failed_v18_evidence": (
        FAILED_V18_ROOT / "failure_evidence.json",
        40_663,
        "d4cca03adc3382f48177cf8cec705189c66ffbc063c5533bcf794629ccaead2b",
    ),
    "failed_v18_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v18/ed25519_public.pem"
        ),
        113,
        "57062811b5b67f514f1454cac3e96fe55f7ab94be2fbbbace43431f68fbb65e4",
    ),
    "failed_v18_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v8/ed25519_public.pem"
        ),
        113,
        "9fa7667e8076cee90a93fee44cfe08bc45a68e9ce28d01507b5c057463102a93",
    ),
    "failed_v19_evidence": (
        FAILED_V19_ROOT / "failure_evidence.json",
        18_407,
        "45cb20e1f2346a278fef33ec83baeca527ad726e8c96526a7977215a1437268b",
    ),
    "failed_v19_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v19/ed25519_public.pem"
        ),
        113,
        "d494f66e8aea206e37e0e803ccb4e0ceb9cf2b244e71eccba6b370f96ddee2e0",
    ),
    "failed_v19_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v9/ed25519_public.pem"
        ),
        113,
        "34d11ff6a699b96aaa8624b30f45fbc8845e6958dc04122b19c413a1a2c12c2d",
    ),
    "rejected_v19_round1_evidence": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v19_candidate_mutation_monitor_provenance_finding"
        / "rejection_evidence.json",
        4_015,
        "d99750a48971e1e872661b6f6435a8b1adfb0a7f0dbdeeb6e47fe938e0c02b2f",
    ),
    "rejected_v20_evidence": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v20_candidate_scope_gate_and_setup_toctou_findings"
        / "rejection_evidence.json",
        10_996,
        "8b851d38304df5aa6f350bd934d5d82c7eacc287faaa3aee9c28cd73a75b0a34",
    ),
    "rejected_v20_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v20/ed25519_public.pem"
        ),
        113,
        "122fb38c395c37d1a4a2786c385110397db30f4b2db0ae3e4944f55355656fa9",
    ),
    "rejected_v20_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v10/ed25519_public.pem"
        ),
        113,
        "09794c2ce162af3bf2f3117f6d11dea0c4bd626cbe50946267609058c6c0c291",
    ),
    "failed_v21_evidence": (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v21_signed_authority_c_pre_audit_canonical_launcher_relation_mismatch"
        / "failure_evidence.json",
        9_865,
        "1180a6db140929298e97ebfd84aa3e1a1efece07b1ceab3d123056c99dcd7109",
    ),
    "failed_v21_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v21/ed25519_public.pem"
        ),
        113,
        "b66aa240810122b7d9be0ff89840339a7a58563163a55b5d2a837178dffebafb",
    ),
    "failed_v21_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v11/ed25519_public.pem"
        ),
        113,
        "825111aa7dcd25e60c6357243e228422d2cf07855d682aae30d90ebdcd5559d2",
    ),
    "rejected_v22_evidence": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v22_candidate_builder_self_binding_review_provenance_v5_watch_findings"
        / "rejection_evidence.json",
        11_084,
        "62d1690831d8279aecb717d4bd6f7058708ea6cbcae23dfed556064ae7cf5cf4",
    ),
    "rejected_v22_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v22/ed25519_public.pem"
        ),
        113,
        "a09adf97ee4b4afb2f90a0d12faf3db4ef1a3a0218676db333a6af6591c76051",
    ),
    "rejected_v22_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v12/ed25519_public.pem"
        ),
        113,
        "94ebaec5e5fc994dc75ffcf189f124726154fccc89ff1ea32add8398f06bf5b2",
    ),
    "rejected_v23_evidence": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v23_nonexistent_node_runtime_path"
        / "rejection_evidence.json",
        11_042,
        "2c6625c0bc3c711e6fba690f91bbffd50c30eb5a885c6b6b6f0a7f1b44a01858",
    ),
    "rejected_v23_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v23/ed25519_public.pem"
        ),
        113,
        "468fb1526e08b4e2cee44ee27d3cd20438beec2298769fca2a413e164759fbe1",
    ),
    "rejected_v23_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v13/ed25519_public.pem"
        ),
        113,
        "d050c0dfea29b469e271967ea7759d4cbb36c3cbec493010163be5938e81e54a",
    ),
    "rejected_v24_evidence": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "rejection_evidence.json",
        13_158,
        "d8ec62cd282190737fd951eafb1eec89d54e9dec36b0af70584fba243e52857b",
    ),
    "rejected_v24_authority_validator": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "validate_tumor_ref_schema_recovery_authority_v24.py",
        442_595,
        "8e2c5b5d4363583082a95cb742dff838e4ee9ac3fb9ff800bfe4dcb4b5387ba2",
    ),
    "rejected_v24_ceremony_builder": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "build_tumor_ref_schema_recovery_authority_v24.py",
        60_648,
        "703623d52452bde8f8ae80cc86dbbb9cce6133e6b4cdf870efdcafedb5f02e9a",
    ),
    "rejected_v24_continuation_verifier": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "verify_tumor_ref_receipt_promotion_recovery_v24.py",
        129_383,
        "eea1927f92d09217e9cb648a791c36f0cb0cb1c9ac4d4cf94c0be4b506fad02d",
    ),
    "rejected_v24_downstream_continuation": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "continue_m2v5_after_tumor_ref_promotion_recovery_v24.py",
        388_261,
        "6968b819e3a2a43f7d3d7aee3d75268279ff31af2ebd1d19a3ea8b3854ab30d2",
    ),
    "rejected_v24_readonly_probe": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "probe_tumor_ref_schema_recovery_sources_v24.py",
        109_466,
        "127bf610fe30a4a34fd76c079a0085a263bfed304f446c77f3543425e9a7174c",
    ),
    "rejected_v24_regression_tests": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "test_tumor_ref_schema_recovery_v24.py",
        169_789,
        "6451daf6bf0dfb74834fffe0bc7cea2ce45e53fb2ab834ee4fb0eccf033df9b3",
    ),
    "rejected_v24_runner_gate_replay": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "replay_m2v5_runner_only_gates_recovery_v24.py",
        153_508,
        "4640cccade1e8166d13132cb75347f2bb77db857cedaeafe62d917c62b07a3d8",
    ),
    "rejected_v24_mendel_review": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "mendel_request_changes.json",
        3_979,
        "c3e11dc27f36bb1c42ea1c27010edaa5aa8361f6ebc3a36339485771a2cb8e3a",
    ),
    "rejected_v24_nash_review": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "nash_request_changes.json",
        4_826,
        "404f40e7111b6b46f3c1292ccd8cf78c057312d250872ddf38c8b3eef80d672c",
    ),
    "rejected_v24_external_review": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "20260723_external_claude_schema_recovery_review_v24.claude_cli.envelope.json",
        19_758,
        "e4b67862fd33d2b37c00a6281e804f563e307a66cbcc27a48f0078b22619f340",
    ),
    "rejected_v24_external_prompt": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "20260723_external_claude_schema_recovery_review_v24_prompt.md",
        7_217,
        "89683e8a831780d24a0e28c85307701aeb6b5bf369a2de768b1f0622307cdcfa",
    ),
    "rejected_v24_external_schema": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "20260723_external_claude_schema_recovery_review_v24.schema.json",
        4_332,
        "e47ddd8f41086038a9b64bafb3f0c0fee11db2dbc03e1738ecfa6cd4a93ad053",
    ),
    "rejected_v24_external_runner": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "run_external_claude_schema_recovery_review_v24.py",
        4_748,
        "4a0b1c2960141798b44486c5b6b3eb9bd89ea024fae7792f64a2bfe942d265e3",
    ),
    "rejected_v24_external_stderr": (
        AUDIT_ROOT
        / "rejected_pre_authority_reviews"
        / "20260723_v24_runtime_execution_and_inventory_scope_findings"
        / "20260723_external_claude_schema_recovery_review_v24.claude_cli.stderr.txt",
        302,
        "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a",
    ),
    "rejected_v24_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260723_v24/ed25519_public.pem"
        ),
        113,
        "c15734463d661be0ead9a8292a3c0c23bd3391233882606df0e930eabb933879",
    ),
    "rejected_v24_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v14/ed25519_public.pem"
        ),
        113,
        "dcff9fbb753ebdda8525597354be18be7b3fa3d22a4f40a020cd9f677736cc8a",
    ),
    "failed_v25_evidence": (
        FAILED_V25_ROOT / "failure_evidence.json",
        588_727,
        "b49d2dac5b610ec34f214a54bfd0c159f795b08ab693761a874e09b284e393be",
    ),
    "failed_v25_authority": (
        FAILED_V25_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v25.bundle"
        / "authority.json",
        263_274,
        "490097e6839ec847efadf1087bdc6eaf162e2d01d66be91c4237da83fe365bfa",
    ),
    "failed_v25_signature": (
        FAILED_V25_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v25.bundle"
        / "authority.ed25519.sig",
        64,
        "912d806fe591778685d9beb951511fd6bf9b2911c80c155824567a0cae623729",
    ),
    "failed_v25_commit": (
        FAILED_V25_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v25.bundle"
        / "commit.json",
        804,
        "51d1ec8e0555bf5ff62cc068db98fab84ccf7e6f13143bd671fac579c67c44ff",
    ),
    "failed_v25_source_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v25.py",
        461_398,
        "58d19165741dfec488fa0603dcbb01118ed0f9addddf982abbdeb16b609ed2cc",
    ),
    "failed_v25_source_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v25.py",
        129_383,
        "f5e6f00b936b7bd315ce8652c441d121f6ca0d1e5c1d679f8e946efc86f062f0",
    ),
    "failed_v25_source_replayer": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v25.py",
        153_508,
        "871c675a6c5967b0593212ffb6b5db5749879a7ef7b821f94eaa0f09b2e1ecfb",
    ),
    "failed_v25_source_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v25.py",
        392_379,
        "151d6810a4f3dc27a50ba5e99aa23ad1f297ce9fb9de2053b816f6b55f2e8e66",
    ),
    "failed_v25_source_probe": (
        AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v25.py",
        119_153,
        "c73f67e2dfae33770165b691702332e2fed73fb497fbfdd381aae056d798e1fb",
    ),
    "failed_v25_source_tests": (
        AUDIT_ROOT / "schema_recovery_tests" / "test_tumor_ref_schema_recovery_v25.py",
        174_699,
        "d5f753e85269f55477a8680a3b15581219560b41e8ba1e9a0aacb86feb2284e4",
    ),
    "failed_v25_source_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v25.py",
        61_065,
        "9da7c2a2b221450908685270924a1ab65c14a0dda6df096664c9550dd40fde7c",
    ),
    "failed_v25_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "archive/20260723_v25_signed_c_tumor_ref_pre_audit_path_relation_failure_01/"
            "ed25519_public.pem"
        ),
        113,
        "8b5dc9f1715a3b8b0dafee138a66135f1a29a8bcc324ae15f28a0fd7f3ea54a6",
    ),
    "failed_v25_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260723_m2v5_terminal_v15_unused_after_signed_v25_c_failure_01/"
            "ed25519_public.pem"
        ),
        113,
        "b0056c3f60d7a7204d782ac1cea31e1e9411200b371295e38b9afc3d5f67a1d1",
    ),
    "failed_v26_evidence": (
        FAILED_V26_ROOT / "failure_evidence.json",
        2_552,
        "0ee3e2923e6cad144b3decff0c8005fb845af3c1b934ba3cc2256b0751961732",
    ),
    "failed_v26_authority": (
        FAILED_V26_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v26.bundle"
        / "authority.json",
        273_966,
        "26e632c9334a1438d754e03397760aa3c59624a12649762d9752eeb13c6f033e",
    ),
    "failed_v26_signature": (
        FAILED_V26_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v26.bundle"
        / "authority.ed25519.sig",
        64,
        "ce9078a9280b81b4ae89637394a6e25b23c2a15ab973ca90559817cfafc2b51a",
    ),
    "failed_v26_commit": (
        FAILED_V26_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v26.bundle"
        / "commit.json",
        804,
        "f581ceedbd5ed02221e3171f826c048c2b3f8811a2de2a76a87f01935c906de4",
    ),
    "failed_v26_source_validator": (
        AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v26.py",
        476_922,
        "ead51f4b3134268f7ef163ede38827d38d1f07383e2d7c44f949d971b2c315ab",
    ),
    "failed_v26_source_verifier": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v26.py",
        129_377,
        "e32a0d57d49af05a123ae0225f8b6a39fd221e9efccfadf96445a62d08e9b885",
    ),
    "failed_v26_source_replayer": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v26.py",
        153_502,
        "2f19d2cdcc96c65ed96c2d0aacbe7fcc24cc303adc01f70d77b9e3150b72355b",
    ),
    "failed_v26_source_continuation": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v26.py",
        398_596,
        "a79d3630570f523068b1220b20aa78369b4ba7b665672a6b01fcbcdb47eedd54",
    ),
    "failed_v26_source_probe": (
        AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v26.py",
        126_050,
        "4fcf5df25224d2e18a663bc846957cc4c43913162e70d792cef85d7fbe25b187",
    ),
    "failed_v26_source_tests": (
        AUDIT_ROOT / "schema_recovery_tests" / "test_tumor_ref_schema_recovery_v26.py",
        188_158,
        "1670c508bb0ba2be3b5e8439b0465b54c8415a240f1d4a5b928a443a76a87f99",
    ),
    "failed_v26_source_builder": (
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v26.py",
        61_065,
        "4ab964b97b63e2903c9c3e431545942743768b5485fda0fbeaabd686f29cf840",
    ),
    "failed_v26_source_dataset_builder": (
        AUDIT_ROOT / "build_all_ssnv_final_report_dataset_schema_recovery_v26.py",
        351_388,
        "cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f",
    ),
    "failed_v26_source_finalizer": (
        AUDIT_ROOT / "finalize_task_b_result_release_schema_recovery_v26.py",
        33_500,
        "ab758306b80c26c17c169285273f99e235ba2c16496cd1b9a2a8fdce974f73d1",
    ),
    "failed_v26_source_report_builder": (
        AUDIT_ROOT / "build_all_ssnv_report_artifact_schema_recovery_v26.py",
        238_722,
        "2bf11448db06c7729e94bd28bfee2f2d83269728eedd09dac185a87b10488527",
    ),
    "failed_v26_review_mendel": (
        FAILED_V26_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v26.json",
        2_866,
        "bf5a7690fccbc7d82c46c792f3447ffead507ab3ab4269ddd7b105e3ddd9da78",
    ),
    "failed_v26_review_nash": (
        FAILED_V26_ROOT / "20260724_tumor_ref_schema_recovery_nash.v26.json",
        3_129,
        "9b3caa622f146df192261ab6613d6fd35edc13707d72c8c5eb7836c7177c56d3",
    ),
    "failed_v26_review_external": (
        FAILED_V26_ROOT
        / "20260724_tumor_ref_schema_recovery_external_claude_opus.v26.json",
        6_452,
        "8881218b0595ed370c9e018ff05263e1e6e4456d49cdaf09a87639ec16fa7603",
    ),
    "failed_v26_verification_receipt": (
        FAILED_V26_ROOT
        / "tumor_ref_source_receipt_promotion_verification.recovery.v26.json",
        40_915,
        "b0658144a819c6d19cc1f446f15909155d533886b3a51a7ee06be8b01c59f28e",
    ),
    "failed_v26_replay_receipt": (
        FAILED_V26_ROOT / "m2v5_runner_only_gate_replay.recovery.v26.json",
        120_133,
        "dccabc1cd95a40bea1a7ca52fbfa7ef0496d6fb19341e514eb8e5de50442939e",
    ),
    "failed_v26_replay_log": (
        FAILED_V26_ROOT / "m2v5_runner_only_gate_replay.recovery.v26.log",
        42_641,
        "f0ae55837a7d1f402de7678c598694131d7b1c1f1c160758d2df7dddb08f8c41",
    ),
    "failed_v26_replay_witness": (
        FAILED_V26_ROOT
        / "m2v5_runner_only_gate_replay.recovery.v26.success_witness.json",
        2_579,
        "6a2764f739c27669152b946a7725c67481b4e8cc609b56a1f33375b4966421b0",
    ),
    "failed_v26_authority_key_archive_record": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "archive/20260724_v26_signed_c_runtime_role_projection_failure_01/"
            "FAILED_KEY_ARCHIVE_RECORD.v1.json"
        ),
        1_078,
        "a3957acd555181a4a4b575db7dc6f0190db6cb882f020f2f30a445aa6eb33372",
    ),
    "failed_v26_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "archive/20260724_v26_signed_c_runtime_role_projection_failure_01/"
            "ed25519_public.pem"
        ),
        113,
        "36fde24fb209c62465d74b4f42a70ca2361bab3c0e326cac25539e0f42049a41",
    ),
    "failed_v26_terminal_key_archive_record": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260724_m2v5_terminal_v16_unused_after_signed_v26_c_failure_01/"
            "FAILED_KEY_ARCHIVE_RECORD.v1.json"
        ),
        1_143,
        "ba9b0d191fb49aa1c5fe897009fd91ce1592c119b96133e70880354f9f73477e",
    ),
    "failed_v26_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260724_m2v5_terminal_v16_unused_after_signed_v26_c_failure_01/"
            "ed25519_public.pem"
        ),
        113,
        "b61e7f75ba3e418098c3c7f9fb19da0b261380d7985833eb6a59b99d6e1aeaee",
    ),
    "rejected_v27_evidence": (
        REJECTED_V27_ROOT / "rejection_evidence.json",
        4_778,
        "3c73542487bfa9a5bc5a107c8b6bc50d8e1af19062ab1925f7821b00a0cce604",
    ),
    "rejected_v27_source_validator": (
        REJECTED_V27_ROOT
        / "sources/validate_tumor_ref_schema_recovery_authority_v27.py",
        497_952,
        "12d194b97f6dac302d3e5234c858ceff3ff00170e5d99344bc9af8f29f601801",
    ),
    "rejected_v27_source_verifier": (
        REJECTED_V27_ROOT
        / "sources/verify_tumor_ref_receipt_promotion_recovery_v27.py",
        129_377,
        "5e83ed0add71ba832bdc6c5b3729caba198d1c1bc9d831426526d1059016af25",
    ),
    "rejected_v27_source_replayer": (
        REJECTED_V27_ROOT / "sources/replay_m2v5_runner_only_gates_recovery_v27.py",
        153_502,
        "5302d3307048fd24df3bd82015b4729248c3226a057c1cf873e6487569ff4aa9",
    ),
    "rejected_v27_source_continuation": (
        REJECTED_V27_ROOT
        / "sources/continue_m2v5_after_tumor_ref_promotion_recovery_v27.py",
        403_619,
        "1005bde75fcbc35fbe53457195bda6e55e01a8df52a3974ab61776c69a170b95",
    ),
    "rejected_v27_source_probe": (
        REJECTED_V27_ROOT / "sources/probe_tumor_ref_schema_recovery_sources_v27.py",
        134_418,
        "8be86c609904d4d891db872a1339b34d33f3dfab3a61002b60d47f8360ea2a75",
    ),
    "rejected_v27_source_tests": (
        REJECTED_V27_ROOT / "sources/test_tumor_ref_schema_recovery_v27.py",
        194_893,
        "a9bb4e654643e82758a6687a6c37f0fb323cc480e0d44c1a9f9c8ee1fbd05891",
    ),
    "rejected_v27_source_builder": (
        REJECTED_V27_ROOT
        / "sources/build_tumor_ref_schema_recovery_authority_v27.py",
        61_065,
        "23137d0f174e435c7cdfcd3c9b4b35efd91063104ff06c17ff0df75c6aa8e5c0",
    ),
    "rejected_v27_source_dataset_builder": (
        REJECTED_V27_ROOT
        / "sources/build_all_ssnv_final_report_dataset_schema_recovery_v27.py",
        351_388,
        "cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f",
    ),
    "rejected_v27_source_finalizer": (
        REJECTED_V27_ROOT
        / "sources/finalize_task_b_result_release_schema_recovery_v27.py",
        33_500,
        "35508a5049e5bba2b1e72e8b08a2e9d0f35665cfb9208a082549923d6a2b3fc7",
    ),
    "rejected_v27_source_report_builder": (
        REJECTED_V27_ROOT
        / "sources/build_all_ssnv_report_artifact_schema_recovery_v27.py",
        238_722,
        "f16cc51bbec3ee23e7598f8ded3eb50f85e13e1634d5c795c8335024cbfeac32",
    ),
    "rejected_v27_review_mendel": (
        REJECTED_V27_ROOT / "reviews/mendel_request_changes.json",
        3_949,
        "50bbe7cfb2b2ca32e402240c698666d4a45462527083e081cc65c8e2525f95fb",
    ),
    "rejected_v27_review_nash": (
        REJECTED_V27_ROOT / "reviews/nash_request_changes.json",
        4_623,
        "cede19b5565cf64f5569eda2a7319236d9da3e0cdaa109cc87ec6feb853f74b7",
    ),
    "rejected_v27_review_external": (
        REJECTED_V27_ROOT / "reviews/external_claude_approve.json",
        6_413,
        "6a97e4a693a6b2b3b79b15281b09276d2ae58cf1fb4f87fedadcde32e7d4a54e",
    ),
    "rejected_v27_transport_prompt": (
        REJECTED_V27_ROOT
        / "review_transport/20260724_external_claude_schema_recovery_review_v27_prompt.md",
        8_652,
        "4df1602a33fcd8f5842c5b50361a632027514db2332711332b6cf9774a1d4983",
    ),
    "rejected_v27_transport_schema": (
        REJECTED_V27_ROOT
        / "review_transport/20260724_external_claude_schema_recovery_review_v27.schema.json",
        4_332,
        "a03238bd84fed5008a5429865e154753066d8fb1f303bea500a20fc53442cc54",
    ),
    "rejected_v27_transport_runner": (
        REJECTED_V27_ROOT
        / "review_transport/run_external_claude_schema_recovery_review_v27.py",
        4_761,
        "4bbc029342a767bfe1a61d53d2c6314627c32a9e996c0160f818b54f389fe1ca",
    ),
    "rejected_v27_transport_publisher": (
        REJECTED_V27_ROOT
        / "review_transport/publish_tumor_ref_schema_recovery_reviews_v27.py",
        8_542,
        "882d3ac0f061381920236b7ab17db09319fda4550c1f1626cf2ff274a9d36a1f",
    ),
    "rejected_v27_transport_envelope": (
        REJECTED_V27_ROOT
        / "review_transport/20260724_external_claude_schema_recovery_review_v27.claude_cli.envelope.json",
        15_658,
        "3a9ca228970554f506420e5ca4fccc0724a19c8c6e2ae0a42a3454b0db5a2f5d",
    ),
    "rejected_v27_transport_stderr": (
        REJECTED_V27_ROOT
        / "review_transport/20260724_external_claude_schema_recovery_review_v27.claude_cli.stderr.txt",
        302,
        "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a",
    ),
    "rejected_v27_authority_key_archive_record": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "archive/20260724_v27_unused_pre_authority_review_rejection_01/"
            "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
        ),
        967,
        "1e653a445197726e8ffa2d66361708814a3a9b6e27bb7c129d9c2c314f7d139a",
    ),
    "rejected_v27_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "archive/20260724_v27_unused_pre_authority_review_rejection_01/"
            "ed25519_public.pem"
        ),
        113,
        "db8c8b5c5ea1a7e4235efb59966c6d41981ed9b6085af6226c07226b56ab3cb3",
    ),
    "rejected_v27_terminal_key_archive_record": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260724_m2v5_terminal_v17_unused_pre_authority_review_rejection_01/"
            "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
        ),
        1_022,
        "59ee02af13e34d0c2394cae45f915972d86e0078223c9c35d5af454a28a0c3cc",
    ),
    "rejected_v27_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260724_m2v5_terminal_v17_unused_pre_authority_review_rejection_01/"
            "ed25519_public.pem"
        ),
        113,
        "355979601138f8dca29534db58e3862f30f30b49526c01d10a97b01ba91c26f9",
    ),
    "failed_v28_evidence": (
        FAILED_V28_ROOT / "failure_evidence.json",
        588_302,
        "a49fd9339448f75834e18b71f6b43a257ab99aaf26d833f8e9a2644d4e65fbcb",
    ),
    "failed_v28_summary": (
        FAILED_V28_ROOT / "SUMMARY.md",
        1_116,
        "a4b98ed12d2f7566e8e9eed293d0090a3fb9d269bd986f02e87670053c7453c8",
    ),
    "failed_v28_authority": (
        FAILED_V28_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v28.bundle/authority.json",
        303_975,
        "01e75a6eb54b0510001abbd72907da1b269db4e80d5842e9f5b35163d18454b6",
    ),
    "failed_v28_authority_signature": (
        FAILED_V28_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v28.bundle/authority.ed25519.sig",
        64,
        "069dd80d3c97bf7e940ef9d9d85490e993cfdad1c3cde632f7e76b9302856ce5",
    ),
    "failed_v28_authority_commit": (
        FAILED_V28_ROOT
        / "tumor_ref_promotion_schema_recovery_authority.v28.bundle/commit.json",
        804,
        "a307a4dabac00457aca3c607eecacf7a2091eb484a6db5d562987b0985b5392e",
    ),
    "failed_v28_provisional_final_dataset": (
        FAILED_V28_ROOT
        / "observed_downstream_outputs/all_ssnv_final_report_dataset_v5_m2v5_source_attested/final_report_dataset.json",
        196_804,
        "fc246a0cde487894137e416507b5026230237d4b6681400a5ab515784eea8feb",
    ),
    "failed_v28_dataset_release_receipt": (
        FAILED_V28_ROOT
        / "observed_downstream_outputs/task_b_final_dataset_release_receipt.v1.json",
        22_217,
        "e2e6896aee1b6e82756a9dabd8e743dc86df7281720990d9a98a7a7496466543",
    ),
    "failed_v28_dataset_release_signature": (
        FAILED_V28_ROOT
        / "observed_downstream_outputs/task_b_final_dataset_release_receipt.v1.json.ed25519.sig",
        64,
        "4f0f2c7c910aa32fc7fc78d52f576a5c4dbea36b88584f87ca7a7cc0730b8425",
    ),
    "failed_v28_authority_key_archive_record": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "archive/20260724_v28_signed_dataset_c_report_metric_key_order_failure_01/"
            "FAILED_KEY_ARCHIVE_RECORD.v1.json"
        ),
        1_112,
        "cdfb48e11b5adaf1a76340e6246d1edb73cf96ad9c7591e1401ece368a4a8c49",
    ),
    "failed_v28_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "archive/20260724_v28_signed_dataset_c_report_metric_key_order_failure_01/"
            "ed25519_public.pem"
        ),
        113,
        "0a4d3a99befa388255c506e5a2a77c2acfd831a8b586d7a5b140585015584e3e",
    ),
    "failed_v28_terminal_key_archive_record": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260724_m2v5_terminal_v18_unused_after_signed_v28_c_failure_01/"
            "FAILED_KEY_ARCHIVE_RECORD.v1.json"
        ),
        1_152,
        "456076174d9322335eba48b3f070601b2e80ef8d77566fc0d84c1810a55c3036",
    ),
    "failed_v28_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260724_m2v5_terminal_v18_unused_after_signed_v28_c_failure_01/"
            "ed25519_public.pem"
        ),
        113,
        "ded65f60581476b08c530c04d795c55bcb393558f393acb0f171629dce47add7",
    ),
    "failed_v28_result_key_archive_record": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
            "archive/20260724_all_ssnv_result_v5_consumed_v28_c_failure_01/"
            "FAILED_KEY_ARCHIVE_RECORD.v1.json"
        ),
        1_237,
        "a25b458c0118248dff43fb7fbc77d14aea20b8df76b28305646aee7b917d9d26",
    ),
    "failed_v28_result_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
            "archive/20260724_all_ssnv_result_v5_consumed_v28_c_failure_01/"
            "ed25519_public.pem"
        ),
        113,
        "5b7d5d026835ec6ec677bcd886bc16ac71444117dabc44a084ce3ede3a4db5a9",
    ),
    "failed_v28_report_key_archive_record": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
            "archive/20260724_all_ssnv_report_v5_unused_v28_c_failure_01/"
            "FAILED_KEY_ARCHIVE_RECORD.v1.json"
        ),
        1_231,
        "5e9d9114233b742626a1d300aaf721bfe3c4a4527717d39db1a930a4fe8f3576",
    ),
    "failed_v28_report_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
            "archive/20260724_all_ssnv_report_v5_unused_v28_c_failure_01/"
            "ed25519_public.pem"
        ),
        113,
        "572e27167e1eea4c39ca53546ba33868b2874ec7fa3ed1682db821b2e50fa439",
    ),
    "rejected_v29_round1_evidence": (
        REJECTED_V29_ROUND1_ROOT / "rejection_evidence.json",
        6_363,
        "6d44a75b32737c0e16a7c05c3801ade734977336982a7f08260b70610442cafa",
    ),
    "rejected_v29_round1_summary": (
        REJECTED_V29_ROUND1_ROOT / "SUMMARY.md",
        3_045,
        "56122d97ee422d9d80d5508ed0655cd7c7fb9be75887b2b4bc38ac0be8c4fa92",
    ),
    "failed_v29_failure_evidence": (
        FAILED_V29_ROOT / "failure_evidence.json",
        91_001,
        "4765ee8326b7d209a5b09b9f6676bd5945bd842e558f36eb1464d62e9474bf4c",
    ),
    "failed_v29_summary": (
        FAILED_V29_ROOT / "SUMMARY.md",
        1_233,
        "a6013d073738c6b891fcc8bda531ebc76df6c0af745b00f0ce1135d23c76cc39",
    ),
    "failed_v29_authority_key_archive_record": (
        FAILED_V29_KEY_ARCHIVES["authority_v29"] / "FAILED_KEY_ARCHIVE_RECORD.v1.json",
        1_139,
        "c7e2c9cce9cd12643fb104c37524beccc99556155ea15b824dccb90ff983ca72",
    ),
    "failed_v29_authority_public_key": (
        FAILED_V29_KEY_ARCHIVES["authority_v29"] / "ed25519_public.pem",
        113,
        "819c0ee9d6c6729ad41f3ea1e8071b90ddf506981c68b6927b8abc98e78e9cda",
    ),
    "authorized_v2_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260722_m2v5_terminal_v2/ed25519_public.pem"
        ),
        113,
        "a98370415594ebc9e8eb166bb70cab3550067b8a83db7d73ad23fac0891cb8b5",
    ),
    "failed_v29_terminal_key_archive_record": (
        FAILED_V29_KEY_ARCHIVES["terminal_v19"] / "FAILED_KEY_ARCHIVE_RECORD.v1.json",
        1_182,
        "66dd4099a57ce9aa94bf0d465fd4dd49dc6e648e5f9287d7737bd2fc1a5311c9",
    ),
    "failed_v29_terminal_public_key": (
        FAILED_V29_KEY_ARCHIVES["terminal_v19"] / "ed25519_public.pem",
        113,
        "04d6acab01d56b0bfe25726a242904afd38bbc3ee47d4e3db29e9eb154e23e8b",
    ),
    "failed_v29_result_key_archive_record": (
        FAILED_V29_KEY_ARCHIVES["result_v6"] / "FAILED_KEY_ARCHIVE_RECORD.v1.json",
        1_206,
        "22fe8a7bdc18839c425c7b8abe881262ee450a7f6b2284e7e44a164484854e07",
    ),
    "failed_v29_result_public_key": (
        FAILED_V29_KEY_ARCHIVES["result_v6"] / "ed25519_public.pem",
        113,
        "84438d0a91200108ee06ad7600a3c5428804f37567c351ba33036843ae837864",
    ),
    "failed_v29_report_key_archive_record": (
        FAILED_V29_KEY_ARCHIVES["report_v6"] / "FAILED_KEY_ARCHIVE_RECORD.v1.json",
        1_206,
        "6d02b12e986ce5b41d2dd18d0051f3425540e4b642ee427e442582fdf147ecfe",
    ),
    "failed_v29_report_public_key": (
        FAILED_V29_KEY_ARCHIVES["report_v6"] / "ed25519_public.pem",
        113,
        "79a684d855ee2d0010691c2a42439389d0e9148d0b84157e04b322c188df6c59",
    ),
    "v30_bootstrap_prepare": (
        V30_BOOTSTRAP_PREPARE,
        19_584,
        "3fb0232160a517d1b202462d64a73643e5932ed541845d58cca778240460b91b",
    ),
    "v30_bootstrap_progress": (
        V30_BOOTSTRAP_PROGRESS,
        2_238,
        "328776bf5124c6e11d15236891e8a3e393bc05b4d90a3ac7f6b7f5d0d1c46343",
    ),
    "v30_bootstrap_receipt": (
        V30_BOOTSTRAP_RECEIPT,
        27_517,
        "6c86539fd8c15bb3e98299474fcd27ebe426bf06d99366d4ea36f0d9212fe4be",
    ),
    "v30_bootstrap_success": (
        V30_BOOTSTRAP_SUCCESS,
        26_802,
        "3d4d09eb64f561fa6e2c1d5bd544f8e371aa44eb44a5daf47ac53b69843dd898",
    ),
    "active_v30_authority_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260724_v30_four_role_rebootstrap/ed25519_public.pem"
        ),
        113,
        "a5b0e0b2c2a9f220d988597b47c8eb1d5446de401a932102948d829ffd0611ed",
    ),
    "recovery_v30_terminal_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260724_m2v5_terminal_v21/ed25519_public.pem"
        ),
        113,
        "3ea7624ed42caba9bd51ade25a4c9a037f0b84689b4e2a3563c8205bbb136fcd",
    ),
    "recovery_v30_result_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
            "20260724_all_ssnv_result_v9_v30_recovery/ed25519_public.pem"
        ),
        113,
        "0d985d9afce029c06275b6932d51db050f807db50ed27c3bf66cd6f9e201f267",
    ),
    "recovery_v30_report_public_key": (
        Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
            "20260724_all_ssnv_report_v9_v30_recovery/ed25519_public.pem"
        ),
        113,
        "8eaf44b95e216320b45ca4109cf34b2c86b8fd4dcf06bb58025c1e27126bf5b5",
    ),
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


def artifact(path: Path) -> dict[str, Any]:
    if path.resolve(strict=True) != path:
        raise RuntimeError(f"Non-canonical source path: {path}")
    data = path.read_bytes()
    opened = os.stat(path, follow_symlinks=False)
    if (
        not stat.S_ISREG(opened.st_mode)
        or stat.S_IMODE(opened.st_mode) != 0o444
        or opened.st_nlink != 1
    ):
        raise RuntimeError(f"Frozen source mode/link drift: {path}")
    return {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": "0o444",
    }


def validate_rejected_v29_round1_archive() -> dict[str, Any]:
    evidence_path = REJECTED_V29_ROUND1_ROOT / "rejection_evidence.json"
    summary_path = REJECTED_V29_ROUND1_ROOT / "SUMMARY.md"
    evidence_record = artifact(evidence_path)
    summary_record = artifact(summary_path)
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    inventory = evidence.get("inventory")
    if (
        evidence_record["size_bytes"] != 6_363
        or evidence_record["sha256"]
        != "6d44a75b32737c0e16a7c05c3801ade734977336982a7f08260b70610442cafa"
        or summary_record["size_bytes"] != 3_045
        or summary_record["sha256"]
        != "56122d97ee422d9d80d5508ed0655cd7c7fb9be75887b2b4bc38ac0be8c4fa92"
        or evidence.get("candidate") != "v29_round1"
        or evidence.get("status")
        != "REJECTED_BEFORE_AUTHORITY_BY_STRICTEST_REPRODUCIBLE_REVIEW_WINS"
        or evidence.get("reviewed_source_set_sha256")
        != "0087b195ce2b7bfb0495f4f3e2c879b11484a645218a318402d14b47879a019e"
        or evidence.get("inventory_file_count") != 20
        or not isinstance(inventory, dict)
        or len(inventory) != 20
        or evidence.get("pass") is not True
    ):
        raise RuntimeError("Rejected v29 round1 evidence drift")
    expected_files = set(inventory) | {"rejection_evidence.json", "SUMMARY.md"}
    observed_files = {
        path.relative_to(REJECTED_V29_ROUND1_ROOT).as_posix()
        for path in REJECTED_V29_ROUND1_ROOT.rglob("*")
        if path.is_file()
    }
    if observed_files != expected_files:
        raise RuntimeError("Rejected v29 round1 archive file-set drift")
    payload_records: dict[str, dict[str, Any]] = {}
    root = REJECTED_V29_ROUND1_ROOT.resolve(strict=True)
    for relative_path, expected in inventory.items():
        relative = Path(relative_path)
        path = (REJECTED_V29_ROUND1_ROOT / relative).resolve(strict=True)
        if relative.is_absolute() or not path.is_relative_to(root):
            raise RuntimeError("Rejected v29 round1 inventory path escaped archive")
        record = artifact(path)
        if (
            not isinstance(expected, list)
            or len(expected) != 2
            or record["size_bytes"] != expected[0]
            or record["sha256"] != expected[1]
        ):
            raise RuntimeError(
                f"Rejected v29 round1 payload identity drift: {relative_path}"
            )
        payload_records[relative_path] = record
    return {
        "evidence": evidence_record,
        "summary": summary_record,
        "payload_inventory": payload_records,
        "payload_inventory_count": len(payload_records),
        "strictest_reproducible_review_wins": True,
        "pass": True,
    }


def bound_fd_artifact(descriptor: int, path: Path, *, expected_mode: int) -> dict[str, Any]:
    opened = os.fstat(descriptor)
    live = os.stat(path, follow_symlinks=False)
    identity = lambda value: (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )
    data = b"".join(
        os.pread(descriptor, min(8 * 1024 * 1024, opened.st_size - offset), offset)
        for offset in range(0, opened.st_size, 8 * 1024 * 1024)
    )
    if (
        identity(opened) != identity(live)
        or not stat.S_ISREG(opened.st_mode)
        or stat.S_IMODE(opened.st_mode) != expected_mode
        or opened.st_nlink != 1
        or len(data) != opened.st_size
    ):
        raise RuntimeError(f"Bound execution source identity drift: {path}")
    return {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": oct(expected_mode),
    }


def ensure_bound_self_execution() -> None:
    if "INTERSUBMOD_PROBE_SOURCE_FD" in os.environ:
        return
    canonical = AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v30.py"
    if (
        Path(__file__).resolve(strict=True) != canonical
        or Path(sys.argv[0]).resolve(strict=True) != canonical
        or sys.argv[0] != str(canonical)
    ):
        raise RuntimeError("Read-only probe bootstrap requires its canonical path")
    python_fd = os.open(PYTHON_TARGET, os.O_RDONLY)
    source_fd = os.open(canonical, os.O_RDONLY)
    os.set_inheritable(python_fd, True)
    os.set_inheritable(source_fd, True)
    environment = {
        **EXPECTED_ENVIRONMENT,
        "INTERSUBMOD_PROBE_SOURCE_FD": str(source_fd),
        "INTERSUBMOD_PROBE_CANONICAL_PATH": str(canonical),
        "INTERSUBMOD_PROBE_PYTHON_FD": str(python_fd),
        "INTERSUBMOD_PROBE_PYTHON_ALIAS": str(PYTHON),
    }
    os.execve(
        f"/proc/self/fd/{python_fd}",
        [str(PYTHON), "-I", "-B", f"/proc/self/fd/{source_fd}"],
        environment,
    )


def verify_self_execution_binding() -> dict[str, Any]:
    source_fd = int(os.environ["INTERSUBMOD_PROBE_SOURCE_FD"])
    python_fd = int(os.environ["INTERSUBMOD_PROBE_PYTHON_FD"])
    canonical_path = Path(os.environ["INTERSUBMOD_PROBE_CANONICAL_PATH"])
    python_alias = os.environ["INTERSUBMOD_PROBE_PYTHON_ALIAS"]
    source_token = f"/proc/self/fd/{source_fd}"
    python_token = f"/proc/self/fd/{python_fd}"
    cmdline = [os.fsdecode(token) for token in Path("/proc/self/cmdline").read_bytes().split(b"\0") if token]
    source_record = bound_fd_artifact(source_fd, canonical_path, expected_mode=0o444)
    python_record = bound_fd_artifact(python_fd, PYTHON_TARGET, expected_mode=0o775)
    if (
        canonical_path != AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v30.py"
        or sys.executable != python_alias
        or not cmdline
        or cmdline[0] != python_alias
        or source_token not in cmdline
        or python_token == cmdline[0]
        or Path(__file__).as_posix() != source_token
    ):
        raise RuntimeError("Read-only probe did not execute its bound source FD")
    return {
        "source": source_record,
        "python_runtime": python_record,
        "canonical_python_argv0": python_alias,
        "script_argv0": sys.argv[0],
        "proc_cmdline_source_fd_present": True,
        "pass": True,
    }


def validate_rejected_v15_key_archives() -> dict[str, Any]:
    contracts = (
        (
            "schema_recovery_authority_v15",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/"
                "intersubmod_verifier_schema_recovery/20260723_v15"
            ),
            Path(
                "/bip7_disk/liaoyoyo2001/.config/"
                "intersubmod_verifier_schema_recovery/archive/"
                "20260723_v15_unused_pre_authority_review_rejection_01"
            ),
            "ed25519_private_one_time.pem",
        ),
        (
            "downstream_terminal_v5",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/"
                "intersubmod_downstream_continuation/20260723_m2v5_terminal_v5"
            ),
            Path(
                "/bip7_disk/liaoyoyo2001/.config/"
                "intersubmod_downstream_continuation/archive/"
                "20260723_m2v5_terminal_v5_unused_pre_authority_review_rejection_01"
            ),
            "ed25519_private_one_time_resident.pem",
        ),
    )
    observed: dict[str, Any] = {}
    for role, original_root, archive_root, private_name in contracts:
        root_stat = os.lstat(archive_root)
        if (
            os.path.lexists(original_root)
            or archive_root.resolve(strict=True) != archive_root
            or not stat.S_ISDIR(root_stat.st_mode)
            or stat.S_IMODE(root_stat.st_mode) != 0o700
        ):
            raise RuntimeError(f"Rejected v15 key archive root drift: {role}")
        private_path = archive_root / private_name
        before = os.lstat(private_path)
        descriptor = os.open(
            private_path,
            os.O_PATH | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0),
        )
        try:
            opened = os.fstat(descriptor)
            if (
                not stat.S_ISREG(opened.st_mode)
                or stat.S_IMODE(opened.st_mode) != 0o400
                or opened.st_size != 119
                or opened.st_nlink != 1
                or (before.st_dev, before.st_ino, before.st_mode, before.st_size)
                != (opened.st_dev, opened.st_ino, opened.st_mode, opened.st_size)
            ):
                raise RuntimeError(f"Rejected v15 private-key metadata drift: {role}")
            observed[role] = {
                "archive": str(archive_root),
                "private_key": {
                    "path": str(private_path),
                    "mode": "0o400",
                    "size_bytes": 119,
                    "link_count": 1,
                },
                "private_key_bytes_read": False,
                "pass": True,
            }
        finally:
            os.close(descriptor)
    return observed


def classify_review_evidence(paths: tuple[Path, ...]) -> tuple[dict[str, bool], str]:
    state = {str(path): os.path.lexists(path) for path in paths}
    if any(state.values()) and not all(state.values()):
        raise RuntimeError("Review evidence is only partially present")
    return state, "all_present" if all(state.values()) else "all_absent"


def recursive_self_check(path: Path, function_name: str) -> Any:
    code = """import json,sys
from pathlib import Path
p=Path(sys.argv[1]).resolve()
function_name=sys.argv[2]
data=p.read_bytes()
namespace={"__name__":"recovery_probe_shadow","__file__":str(p),"__package__":""}
result={}
def trace(frame,event,arg):
    if (event=="return" and frame.f_code.co_name=="<module>"
            and frame.f_code.co_filename==str(p)
            and function_name in namespace and not result):
        observed=namespace[function_name](data)
        result["value"]="PASS" if observed is None else observed
    return trace
sys.settrace(trace)
exec(compile(data,str(p),"exec"),namespace)
sys.settrace(None)
print(json.dumps(result,sort_keys=True))
"""
    result = subprocess.run(
        [str(PYTHON), "-I", "-B", "-c", code, str(path), function_name],
        cwd=REPO_ROOT,
        env=EXPECTED_ENVIRONMENT,
        capture_output=True,
        check=False,
        timeout=180,
    )
    if result.returncode != 0 or result.stderr:
        raise RuntimeError(
            f"Recursive source self-check failed: {path.name}; exit={result.returncode}"
        )
    payload = json.loads(result.stdout)
    if "value" not in payload:
        raise RuntimeError(f"Recursive source self-check returned no value: {path.name}")
    return payload["value"]


def run_archived_state_runner_replay() -> dict[str, Any]:
    replayer_path = EXPECTED_SOURCES["runner_gate_replay"][0]
    spec = importlib.util.spec_from_file_location(
        "schema_recovery_v30_archive_rebase_probe", replayer_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load v30 replayer for archive-state canary")
    replayer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(replayer)

    runner_path = replayer.COMPLETION_RUNNER
    runner_data = runner_path.read_bytes()
    runner_lines = runner_data.splitlines(keepends=True)
    if (
        len(runner_lines) < replayer.EXPECTED_RUNNER_LINE_COUNT_MINIMUM
        or replayer.RUNNER_INCLUDED_LINE_START != 1
        or replayer.RUNNER_INCLUDED_LINE_END != 358
        or replayer.RUNNER_FIRST_EXCLUDED_LINE != 359
    ):
        raise RuntimeError("Reviewed runner prefix boundary drift")
    runner_prefix = b"".join(runner_lines[: replayer.RUNNER_INCLUDED_LINE_END])
    if (
        len(runner_prefix.splitlines()) != 358
        or hashlib.sha256(runner_prefix).hexdigest()
        != replayer.EXPECTED_RUNNER_PREFIX_SHA256
    ):
        raise RuntimeError("Original runner prefix identity drift")

    result_archive_path = replayer.FAILED_V28_RESULT_KEY_ARCHIVE_RECORD
    result_public_path = replayer.FAILED_V28_RESULT_PUBLIC_KEY
    report_archive_path = replayer.FAILED_V28_REPORT_KEY_ARCHIVE_RECORD
    report_public_path = replayer.FAILED_V28_REPORT_PUBLIC_KEY
    report_private_path = replayer.FAILED_V28_REPORT_PRIVATE_KEY
    report_private_fd = os.open(
        report_private_path,
        os.O_PATH | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        report_private_stat = os.fstat(report_private_fd)
        if (
            not stat.S_ISREG(report_private_stat.st_mode)
            or stat.S_IMODE(report_private_stat.st_mode) != 0o400
            or report_private_stat.st_size != 290
            or report_private_stat.st_nlink != 1
        ):
            raise RuntimeError("Archived report private-key metadata drift")
        archive_inputs = replayer.validate_archive_rebase_inputs(
            result_archive_data=result_archive_path.read_bytes(),
            result_archive_record=artifact(result_archive_path),
            result_public_record=artifact(result_public_path),
            report_archive_data=report_archive_path.read_bytes(),
            report_archive_record=artifact(report_archive_path),
            report_public_record=artifact(report_public_path),
            report_private_record={
                "path": str(report_private_path),
                "mode": "0o400",
                "size_bytes": report_private_stat.st_size,
            },
            report_private_link_count=report_private_stat.st_nlink,
        )
    finally:
        os.close(report_private_fd)

    old_live_paths = tuple(replayer.FAILED_V28_LIVE_KEY_PATHS)
    downstream_paths = tuple(replayer.DOWNSTREAM_OUTPUT_SLOTS)
    protected_paths = tuple(
        dict.fromkeys((*FORBIDDEN_OUTPUT_SLOTS, *old_live_paths, *downstream_paths))
    )
    before = {str(path): os.path.lexists(path) for path in protected_paths}
    if any(before.values()):
        raise RuntimeError("Archive-state canary started with a protected path present")

    raw_result = subprocess.run(
        [str(replayer.BASH), "-s"],
        input=runner_prefix,
        cwd=REPO_ROOT,
        env=EXPECTED_ENVIRONMENT,
        capture_output=True,
        check=False,
        timeout=300,
    )
    expected_raw_stderr = (
        b"Required file is missing or empty: /bip7_disk/liaoyoyo2001/.config/"
        b"intersubmod_result_release_authority/"
        b"20260719_all_ssnv_result_v5_post_reboot_bootstrap/"
        b"ed25519_public.pem\n"
    )
    if (
        raw_result.returncode != 1
        or raw_result.stdout != b""
        or raw_result.stderr != expected_raw_stderr
        or hashlib.sha256(raw_result.stderr).hexdigest()
        != "18477fe5a82535b402f6f04e72ebcea147f5c5cee0475dba6e250ae75078c4a0"
    ):
        raise RuntimeError("Raw archived-state runner prefix did not fail as expected")

    transformed, archive_rebase = replayer.archive_rebase_runner_prefix(runner_prefix)
    if (
        archive_rebase.get("mapping_count") != 3
        or archive_rebase.get("private_key_bytes_read") is not False
        or archive_rebase.get("host_live_key_paths_created") is not False
        or archive_rebase.get("pass") is not True
    ):
        raise RuntimeError("Archive-rebase contract summary drift")
    for original, archived, role in replayer.RUNNER_ARCHIVE_REBASE:
        if (
            runner_prefix.count(original) != 1
            or runner_prefix.count(archived) != 0
            or transformed.count(original) != 0
            or transformed.count(archived) != 1
            or not any(
                mapping.get("role") == role
                and mapping.get("replacement_count") == 1
                for mapping in archive_rebase["mappings"]
            )
        ):
            raise RuntimeError(f"Archive-rebase mapping cardinality drift: {role}")
    if (
        len(transformed.splitlines()) != 358
        or hashlib.sha256(transformed).hexdigest()
        != replayer.EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SHA256
    ):
        raise RuntimeError("Archive-rebased runner prefix identity drift")

    transformed_result = subprocess.run(
        [str(replayer.BASH), "-s"],
        input=transformed,
        cwd=REPO_ROOT,
        env=EXPECTED_ENVIRONMENT,
        capture_output=True,
        check=False,
        timeout=300,
    )
    after = {str(path): os.path.lexists(path) for path in protected_paths}
    if (
        transformed_result.returncode != 0
        or transformed_result.stdout != b""
        or transformed_result.stderr != b""
        or before != after
        or any(after.values())
    ):
        raise RuntimeError("Archive-rebased runner canary failed or changed protected state")
    return {
        "contract": "archived_state_runner_prefix_ab_canary_v1",
        "runner_source": artifact(runner_path),
        "line_range": {"start": 1, "end": 358, "physical_line_count": 358},
        "raw_execution": {
            "exit_code": raw_result.returncode,
            "stdin_size_bytes": len(runner_prefix),
            "stdin_sha256": hashlib.sha256(runner_prefix).hexdigest(),
            "stdout_sha256": hashlib.sha256(raw_result.stdout).hexdigest(),
            "stderr_sha256": hashlib.sha256(raw_result.stderr).hexdigest(),
            "expected_archived_live_path_failure_verified": True,
        },
        "archive_rebase_inputs": archive_inputs,
        "archive_rebase": archive_rebase,
        "transformed_execution": {
            "exit_code": transformed_result.returncode,
            "stdin_size_bytes": len(transformed),
            "stdin_sha256": hashlib.sha256(transformed).hexdigest(),
            "stdout_sha256": hashlib.sha256(transformed_result.stdout).hexdigest(),
            "stderr_sha256": hashlib.sha256(transformed_result.stderr).hexdigest(),
            "exact_transformed_bytes_executed": True,
        },
        "old_live_key_path_count": len(old_live_paths),
        "downstream_output_slot_count": len(downstream_paths),
        "protected_path_count": len(protected_paths),
        "protected_paths_absent_before_and_after": True,
        "report_private_key_bytes_read": False,
        "non_authoritative": True,
        "pass": True,
    }


def run_pilot_verifier() -> dict[str, Any]:
    pilot = artifact(PILOT_VERIFIER)
    if (
        pilot["size_bytes"] != 121_770
        or pilot["sha256"]
        != "05753bd833fdef06e95af05ec8f582df1d38ed5e06c2a2dc824d1d65d5fc4372"
    ):
        raise RuntimeError("Frozen eight-key pilot verifier identity drift")
    python_fd = os.open(PYTHON_TARGET, os.O_RDONLY | os.O_CLOEXEC)
    verifier_fd = os.open(PILOT_VERIFIER, os.O_RDONLY | os.O_CLOEXEC)
    try:
        command = [
            f"/proc/self/fd/{python_fd}",
            "-I",
            "-B",
            f"/proc/self/fd/{verifier_fd}",
            "--verify",
        ]
        result = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=EXPECTED_ENVIRONMENT,
            pass_fds=(python_fd, verifier_fd),
            capture_output=True,
            check=False,
            timeout=600,
        )
    finally:
        os.close(verifier_fd)
        os.close(python_fd)
    if result.returncode != 0:
        raise RuntimeError(
            "Eight-key pilot verifier failed: "
            f"exit={result.returncode} stderr_sha256={hashlib.sha256(result.stderr).hexdigest()}"
        )
    payload = json.loads(result.stdout)
    if payload.get("pass") is not True or payload.get("mode") != "verify":
        raise RuntimeError("Eight-key pilot verifier did not return pass=true")
    return {
        "source": pilot,
        "normalized_command": [
            "/proc/self/fd/<bound-python-fd>",
            "-I",
            "-B",
            "/proc/self/fd/<bound-verifier-fd>",
            "--verify",
        ],
        "exit_code": result.returncode,
        "stdout_sha256": hashlib.sha256(result.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(result.stderr).hexdigest(),
        "pass": True,
    }


def parse_exact_pytest_summary(stdout: bytes) -> str | None:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    expected = EXPECTED_REGRESSION_SUMMARY.encode("ascii")
    if (
        not lines
        or re.fullmatch(
            re.escape(expected) + rb" in [0-9]+(?:\.[0-9]+)?s", lines[-1]
        )
        is None
    ):
        return None
    return lines[-1].split(b" in ", 1)[0].decode("ascii")


def run_regression_tests() -> dict[str, Any]:
    test_path = (
        AUDIT_ROOT
        / "schema_recovery_tests"
        / "test_tumor_ref_schema_recovery_v30.py"
    )
    python_fd = os.open(PYTHON_TARGET, os.O_RDONLY | os.O_CLOEXEC)
    test_fd = os.open(test_path, os.O_RDONLY | os.O_CLOEXEC)
    try:
        source_record = bound_fd_artifact(test_fd, test_path, expected_mode=0o444)
        python_record = bound_fd_artifact(python_fd, PYTHON_TARGET, expected_mode=0o775)
        bootstrap = """import hashlib,os,stat,sys,types
source_token=sys.argv[1]
fd=int(source_token.rsplit('/',1)[-1])
canonical=sys.argv[2]
repo_root=sys.argv[3]
pytest_args=sys.argv[4:]
opened=os.fstat(fd)
live=os.stat(canonical,follow_symlinks=False)
data=b''.join(os.pread(fd,min(8*1024*1024,opened.st_size-o),o) for o in range(0,opened.st_size,8*1024*1024))
if ((opened.st_dev,opened.st_ino,opened.st_size)!=(live.st_dev,live.st_ino,live.st_size)
        or not stat.S_ISREG(opened.st_mode) or stat.S_IMODE(opened.st_mode)!=0o444
        or opened.st_nlink!=1 or len(data)!=opened.st_size
        or hashlib.sha256(data).hexdigest()!=os.environ['INTERSUBMOD_REGRESSION_TEST_SHA256']):
    raise SystemExit('bound regression source identity drift')
name=os.path.relpath(canonical,repo_root)[:-3].replace(os.sep,'.')
module=types.ModuleType(name)
module.__file__=canonical
module.__package__=name.rpartition('.')[0]
module.__dict__['_INTERSUBMOD_BOUND_SOURCE_FD']=fd
sys.modules[name]=module
exec(compile(data,canonical,'exec',dont_inherit=True),module.__dict__)
import pytest
raise SystemExit(pytest.main(pytest_args))
"""
        command = [
            str(PYTHON),
            "-I",
            "-B",
            "-c",
            bootstrap,
            f"/proc/self/fd/{test_fd}",
            str(test_path),
            str(REPO_ROOT),
            "-q",
            "-p",
            "no:cacheprovider",
            str(test_path),
        ]
        environment = {
            **EXPECTED_ENVIRONMENT,
            "INTERSUBMOD_REGRESSION_TEST_FD": str(test_fd),
            "INTERSUBMOD_REGRESSION_TEST_CANONICAL_PATH": str(test_path),
            "INTERSUBMOD_REGRESSION_TEST_SHA256": source_record["sha256"],
            "INTERSUBMOD_REGRESSION_PYTHON_FD": str(python_fd),
            "INTERSUBMOD_REGRESSION_PYTHON_ALIAS": str(PYTHON),
        }
        result = subprocess.run(
            command,
            executable=f"/proc/self/fd/{python_fd}",
            cwd=REPO_ROOT,
            env=environment,
            pass_fds=(python_fd, test_fd),
            capture_output=True,
            check=False,
            timeout=300,
        )
        if (
            source_record
            != bound_fd_artifact(test_fd, test_path, expected_mode=0o444)
            or python_record
            != bound_fd_artifact(python_fd, PYTHON_TARGET, expected_mode=0o775)
        ):
            raise RuntimeError("Regression execution FD identity changed")
    finally:
        os.close(test_fd)
        os.close(python_fd)
    summary = parse_exact_pytest_summary(result.stdout)
    if result.returncode != 0 or summary is None:
        raise RuntimeError(
            f"Schema-recovery regression tests failed: exit={result.returncode}"
        )
    return {
        "command": command,
        "exit_code": result.returncode,
        "stdout_sha256": hashlib.sha256(result.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(result.stderr).hexdigest(),
        "summary": summary,
        "source_execution_binding": {
            "source": source_record,
            "python_runtime": python_record,
            "execution": (
                "pytest_preloaded_module_compiled_from_bound_source_fd_"
                "via_bound_python_fd"
            ),
            "canonical_python_argv0": str(PYTHON),
            "pass": True,
        },
        "pass": True,
    }


def run_prior_recorded_verifier_recheck() -> dict[str, Any]:
    verifier = AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v1.py"
    recorded_path = REQUIRED_PRIOR_INPUTS["verification_receipt_v1"][0]
    python_fd = os.open(PYTHON_TARGET, os.O_RDONLY | os.O_CLOEXEC)
    verifier_fd = os.open(verifier, os.O_RDONLY | os.O_CLOEXEC)
    try:
        command = [
            f"/proc/self/fd/{python_fd}",
            "-I",
            "-B",
            f"/proc/self/fd/{verifier_fd}",
            "--verify",
        ]
        result = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=EXPECTED_ENVIRONMENT,
            pass_fds=(python_fd, verifier_fd),
            capture_output=True,
            check=False,
            timeout=600,
        )
    finally:
        os.close(verifier_fd)
        os.close(python_fd)
    if result.returncode != 0 or result.stderr:
        raise RuntimeError("Fresh inherited-FD verifier recheck failed")
    fresh = json.loads(result.stdout)
    recorded = json.loads(recorded_path.read_bytes())
    replayer_path = AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v1.py"
    spec = importlib.util.spec_from_file_location(
        "schema_recovery_v30_probe_replayer", replayer_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load v2 replayer for receipt comparison")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.require_fresh_and_recorded_verification_agree(fresh, recorded)
    return {
        "normalized_command": [
            "/proc/self/fd/<bound-python-fd>",
            "-I",
            "-B",
            "/proc/self/fd/<bound-verifier-fd>",
            "--verify",
        ],
        "exit_code": result.returncode,
        "fresh_mode": fresh.get("mode"),
        "recorded_mode": recorded.get("mode"),
        "recorded_receipt_sha256": hashlib.sha256(recorded_path.read_bytes()).hexdigest(),
        "pass": True,
    }


def run_static_contract() -> dict[str, Any]:
    validator = EXPECTED_SOURCES["authority_validator"][0]
    python_fd = os.open(PYTHON_TARGET, os.O_RDONLY | os.O_CLOEXEC)
    validator_fd = os.open(validator, os.O_RDONLY | os.O_CLOEXEC)
    try:
        validator_record = bound_fd_artifact(
            validator_fd, validator, expected_mode=0o444
        )
        result = subprocess.run(
            [
                str(PYTHON),
                "-I",
                "-B",
                f"/proc/self/fd/{validator_fd}",
                "--static-contract",
            ],
            executable=f"/proc/self/fd/{python_fd}",
            cwd=REPO_ROOT,
            env=EXPECTED_ENVIRONMENT,
            pass_fds=(python_fd, validator_fd),
            capture_output=True,
            check=False,
            timeout=300,
        )
        if validator_record != bound_fd_artifact(
            validator_fd, validator, expected_mode=0o444
        ):
            raise RuntimeError("Static validator execution FD identity changed")
    finally:
        os.close(validator_fd)
        os.close(python_fd)
    if result.returncode != 0 or result.stderr:
        raise RuntimeError("v30 recovery static contract failed")
    payload = json.loads(result.stdout)
    absence_contract = payload.get("ceremony_absence_contract", {})
    failed_recoveries = payload.get("prior_failed_signed_recovery", {})
    failed_v9 = failed_recoveries.get("v9", {})
    failed_v10 = failed_recoveries.get("v10", {})
    failed_v11 = failed_recoveries.get("v11", {})
    failed_v12 = failed_recoveries.get("v12", {})
    failed_v14 = failed_recoveries.get("v14", {})
    failed_v16 = failed_recoveries.get("v16", {})
    failed_v17 = failed_recoveries.get("v17", {})
    failed_v18 = failed_recoveries.get("v18", {})
    failed_v19 = failed_recoveries.get("v19", {})
    failed_v21 = failed_recoveries.get("v21", {})
    failed_v25 = failed_recoveries.get("v25", {})
    failed_v26 = failed_recoveries.get("v26", {})
    failed_v28 = failed_recoveries.get("v28", {})
    failed_v29 = failed_recoveries.get("v29", {})
    v30_key_bootstrap = payload.get("fresh_key_bootstrap", {})
    rejected_v13 = payload.get("rejected_unsigned_generations", {}).get("v13", {})
    rejected_v15 = payload.get("rejected_unsigned_generations", {}).get("v15", {})
    rejected_v19_round1 = payload.get("rejected_unsigned_generations", {}).get(
        "v19_round1", {}
    )
    rejected_v20 = payload.get("rejected_unsigned_generations", {}).get("v20", {})
    rejected_v22 = payload.get("rejected_unsigned_generations", {}).get("v22", {})
    rejected_v23 = payload.get("rejected_unsigned_generations", {}).get("v23", {})
    rejected_v24 = payload.get("rejected_unsigned_generations", {}).get("v24", {})
    rejected_v27 = payload.get("rejected_unsigned_generations", {}).get("v27", {})
    rejected_v29_round1 = payload.get("rejected_unsigned_generations", {}).get(
        "v29_round1", {}
    )
    terminal_key_state = payload.get("terminal_key_state", {})
    expected_paths_sha256 = hashlib.sha256(
        json.dumps(
            sorted(str(path) for path in FORBIDDEN_OUTPUT_SLOTS),
            ensure_ascii=True,
            separators=(",", ":"),
        ).encode("ascii")
    ).hexdigest()
    expected_patterns_sha256 = hashlib.sha256(
        json.dumps(
            list(CEREMONY_STAGING_PATTERNS),
            ensure_ascii=True,
            separators=(",", ":"),
        ).encode("ascii")
    ).hexdigest()
    if (
        payload.get("pass") is not True
        or payload.get("authority_id") != "20260724_tumor_ref_schema_recovery_v30"
        or payload.get("prior_recovery_chain", {}).get("authority_signature_verified")
        is not True
        or failed_v9.get("authority_signature_verified") is not True
        or failed_v9.get("atomic_commit_verified") is not True
        or failed_v9.get("historical_authorized_output_slot_count") != 17
        or failed_v9.get("current_recovery_output_slot_count") != 23
        or failed_v9.get("runtime_outputs_created") is not False
        or failed_v10.get("authority_signature_verified") is not True
        or failed_v10.get("atomic_commit_verified") is not True
        or failed_v10.get("verification_receipt_created") is not True
        or failed_v10.get("replay_receipt_and_log_created") is not True
        or failed_v10.get("continuation_outputs_created") is not False
        or failed_v10.get("canonical_downstream_outputs_created") is not False
        or failed_v11.get("authority_signature_verified") is not True
        or failed_v11.get("basic_to_full_stat_mismatch_reproduced") is not True
        or failed_v12.get("authority_signature_verified") is not True
        or failed_v12.get("atomic_commit_verified") is not True
        or failed_v12.get("executable_alias_relation_mismatch_reproduced") is not True
        or failed_v12.get("continuation_outputs_created") is not False
        or failed_v14.get("authority_signature_verified") is not True
        or failed_v14.get("atomic_commit_verified") is not True
        or failed_v14.get("verification_receipt_created") is not True
        or failed_v14.get("metadata_only_relation_schema_mismatch_reproduced")
        is not True
        or failed_v14.get("replay_outputs_created") is not False
        or failed_v14.get("continuation_outputs_created") is not False
        or failed_v16.get("authority_signature_verified") is not True
        or failed_v16.get("atomic_commit_verified") is not True
        or failed_v16.get("verification_receipt_created") is not True
        or failed_v16.get("replay_receipt_and_log_created") is not True
        or failed_v16.get("replay_success_witness_created") is not True
        or failed_v16.get("legacy_eight_key_stat_schema_mismatch_reproduced")
        is not True
        or failed_v16.get("continuation_outputs_created") is not False
        or failed_v16.get("canonical_downstream_outputs_created") is not False
        or failed_v17.get("authority_signature_verified") is not True
        or failed_v17.get("atomic_commit_verified") is not True
        or failed_v17.get("verification_receipt_created") is not True
        or failed_v17.get("replay_receipt_and_log_created") is not True
        or failed_v17.get("replay_success_witness_created") is not True
        or failed_v17.get("continuation_child_started") is not True
        or failed_v17.get("fresh_verifier_passed") is not True
        or failed_v17.get("metadata_enrichment_conflict_reproduced") is not True
        or failed_v17.get("syscall_trace_evidence_bound") is not True
        or failed_v17.get("continuation_outputs_created") is not False
        or failed_v17.get("canonical_downstream_outputs_created") is not False
        or failed_v18.get("authority_signature_verified") is not True
        or failed_v18.get("atomic_commit_verified") is not True
        or failed_v18.get("verification_receipt_created") is not True
        or failed_v18.get("replay_receipt_and_log_created") is not True
        or failed_v18.get("replay_success_witness_created") is not True
        or failed_v18.get("continuation_child_started") is not True
        or failed_v18.get("fresh_verifier_passed") is not True
        or failed_v18.get("tumor_ref_summary_alias_noncanonical_reproduced")
        is not True
        or failed_v18.get("syscall_trace_evidence_bound") is not True
        or failed_v18.get("continuation_outputs_created") is not False
        or failed_v18.get("canonical_downstream_outputs_created") is not False
        or failed_v19.get("authority_signature_verified") is not True
        or failed_v19.get("atomic_commit_verified") is not True
        or failed_v19.get("verification_receipt_created") is not True
        or failed_v19.get("replay_receipt_and_log_created") is not True
        or failed_v19.get("replay_success_witness_created") is not True
        or failed_v19.get("continuation_child_started") is not True
        or failed_v19.get("fresh_verifier_passed") is not True
        or failed_v19.get("mode000_inotify_eacces_reproduced") is not True
        or failed_v19.get("syscall_trace_evidence_bound") is not True
        or failed_v19.get("continuation_outputs_created") is not False
        or failed_v19.get("canonical_downstream_outputs_created") is not False
        or failed_v21.get("authority_signature_verified") is not True
        or failed_v21.get("atomic_commit_verified") is not True
        or failed_v21.get("verification_receipt_created") is not True
        or failed_v21.get("replay_receipt_log_and_success_witness_created")
        is not True
        or failed_v21.get("continuation_incident_created") is not True
        or failed_v21.get("canonical_launcher_relation_mismatch_reproduced")
        is not True
        or failed_v21.get("intermediate_downstream_outputs_archived") is not True
        or failed_v21.get("original_observed_downstream_slots_absent") is not True
        or failed_v21.get("continuation_terminal_outputs_created") is not False
        or failed_v21.get("final_dataset_created") is not False
        or failed_v25.get("pass") is not True
        or failed_v25.get("failure", {}).get("stage")
        != "C_FINAL_DATASET_INTEGRATION"
        or failed_v25.get("failure", {}).get("contract")
        != "tumor_ref_v1_declared_vs_v6_canonical_audit_path"
        or failed_v25.get("failure", {}).get("scientific_payload_changed") is not False
        or failed_v25.get("failure", {}).get("terminal_signature_created") is not False
        or failed_v25.get("authority_retired_private_key", {}).get("mode") != "0o0"
        or failed_v25.get("terminal_unused_private_key", {}).get("mode") != "0o400"
        or failed_v26.get("pass") is not True
        or failed_v26.get("failure", {}).get("stage")
        != "C_PRE_DOWNSTREAM_VALIDATE_PROMOTION_CHAIN"
        or failed_v26.get("failure", {}).get("contract")
        != "historical_signed_runtime_projection_exact_11_roles"
        or failed_v26.get("failure", {}).get("current_reviewed_runtime_role_count")
        != 14
        or failed_v26.get("failure", {}).get("scientific_payload_changed") is not False
        or failed_v26.get("failure", {}).get("terminal_signature_created") is not False
        or failed_v26.get("authority_retired_private_key", {}).get("mode") != "0o0"
        or failed_v26.get("terminal_unused_private_key", {}).get("mode") != "0o400"
        or failed_v28.get("pass") is not True
        or failed_v29.get("pass") is not True
        or failed_v29.get("authority_signature_verified") is not True
        or failed_v29.get("failure", {}).get("stage")
        != "R_RUNNER_PREFIX_BEFORE_FIRST_R_OUTPUT"
        or failed_v29.get("failure", {}).get("root_cause")
        != "historical_live_key_paths_absent_after_required_quarantine"
        or failed_v29.get("failure", {}).get("scientific_payload_changed") is not False
        or failed_v29.get("failure", {}).get("claim_ceiling_changed") is not False
        or len(failed_v29.get("archived_artifacts", {})) != 15
        or len(failed_v29.get("key_archives", {})) != 4
        or "v30_key_bootstrap" in failed_recoveries
        or v30_key_bootstrap.get("pass") is not True
        or v30_key_bootstrap.get("release_authority_granted") is not False
        or v30_key_bootstrap.get("four_role_public_keys_pairwise_distinct") is not True
        or len(v30_key_bootstrap.get("public_keys", {})) != 4
        or v30_key_bootstrap.get("pre_authority_archive", {}).get("pass") is not True
        or v30_key_bootstrap.get("failed_v8_pre_ready_archive", {}).get("pass")
        is not True
        or rejected_v13.get("authority_created") is not False
        or rejected_v13.get("formal_reviews_published") is not False
        or rejected_v13.get("signatures_created") is not False
        or rejected_v13.get("outputs_remain_absent") is not True
        or rejected_v13.get("blocking_finding_count") != 2
        or rejected_v13.get("pass") is not True
        or rejected_v15.get("authority_created") is not False
        or rejected_v15.get("formal_reviews_published") is not False
        or rejected_v15.get("signatures_created") is not False
        or rejected_v15.get("outputs_remain_absent") is not True
        or rejected_v15.get("strictest_review_wins") is not True
        or rejected_v15.get("pass") is not True
        or rejected_v19_round1.get("authority_created") is not False
        or rejected_v19_round1.get("signature_created") is not False
        or rejected_v19_round1.get("strictest_review_wins") is not True
        or rejected_v19_round1.get("findings_corrected_in_active_candidate")
        is not True
        or rejected_v19_round1.get("pass") is not True
        or rejected_v20.get("authority_created") is not False
        or rejected_v20.get("signature_created") is not False
        or rejected_v20.get("formal_reviews_published") is not False
        or rejected_v20.get("outputs_remain_absent") is not True
        or rejected_v20.get("strictest_review_wins") is not True
        or rejected_v20.get("blocking_high_count") != 1
        or rejected_v20.get("blocking_medium_count") != 2
        or rejected_v20.get("findings_corrected_in_active_candidate") is not True
        or rejected_v20.get("pass") is not True
        or rejected_v22.get("authority_created") is not False
        or rejected_v22.get("signature_created") is not False
        or rejected_v22.get("formal_reviews_published") is not False
        or rejected_v22.get("outputs_remain_absent") is not True
        or rejected_v22.get("strictest_review_wins") is not True
        or rejected_v22.get("review_authorship_cryptographically_proven") is not False
        or rejected_v22.get("findings_corrected_in_active_candidate") is not True
        or rejected_v22.get("pass") is not True
        or rejected_v23.get("authority_created") is not False
        or rejected_v23.get("signature_created") is not False
        or rejected_v23.get("formal_reviews_published") is not False
        or rejected_v23.get("outputs_remain_absent") is not True
        or rejected_v23.get("strictest_review_wins") is not True
        or rejected_v23.get("review_authorship_cryptographically_proven") is not False
        or rejected_v23.get("findings_corrected_in_active_candidate") is not True
        or rejected_v23.get("pass") is not True
        or rejected_v24.get("authority_created") is not False
        or rejected_v24.get("signature_created") is not False
        or rejected_v24.get("formal_reviews_published") is not False
        or rejected_v24.get("outputs_remain_absent") is not True
        or rejected_v24.get("strictest_review_wins") is not True
        or rejected_v24.get("review_authorship_cryptographically_proven") is not False
        or rejected_v24.get("findings_corrected_in_active_candidate") is not True
        or rejected_v24.get("pass") is not True
        or rejected_v27.get("authority_created") is not False
        or rejected_v27.get("signature_created") is not False
        or rejected_v27.get("formal_reviews_published") is not False
        or rejected_v27.get("outputs_remain_absent") is not True
        or rejected_v27.get("strictest_reproducible_finding_wins") is not True
        or rejected_v27.get("review_authorship_cryptographically_proven") is not False
        or rejected_v27.get("findings_corrected_in_active_candidate") is not True
        or rejected_v27.get("pass") is not True
        or rejected_v29_round1.get("payload_inventory_count") != 20
        or rejected_v29_round1.get("mendel_initial_verdict")
        != "APPROVE_SUPERSEDED"
        or rejected_v29_round1.get("mendel_corrected_verdict")
        != "REQUEST_CHANGES"
        or rejected_v29_round1.get("nash_verdict") != "REQUEST_CHANGES"
        or rejected_v29_round1.get("strictest_reproducible_review_wins")
        is not True
        or rejected_v29_round1.get("authority_created_before_correction")
        is not False
        or rejected_v29_round1.get("any_key_consumed_before_correction")
        is not False
        or rejected_v29_round1.get("findings_corrected_in_active_candidate")
        is not True
        or rejected_v29_round1.get("pass") is not True
        or terminal_key_state.get("legacy_private_key", {}).get("mode") != "0o400"
        or terminal_key_state.get("failed_v16_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("failed_v17_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("failed_v18_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("failed_v19_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("rejected_v20_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("failed_v21_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("rejected_v22_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("rejected_v23_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("rejected_v24_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("failed_v25_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("failed_v26_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("rejected_v27_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("failed_v28_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("failed_v29_private_key", {}).get("mode")
        != "0o400"
        or terminal_key_state.get("recovery_private_key", {}).get("mode") != "0o400"
        or len(
            {
                terminal_key_state.get("legacy_public_key", {}).get("sha256"),
                terminal_key_state.get("failed_v16_public_key", {}).get("sha256"),
                terminal_key_state.get("failed_v17_public_key", {}).get("sha256"),
                terminal_key_state.get("failed_v18_public_key", {}).get("sha256"),
                terminal_key_state.get("failed_v19_public_key", {}).get("sha256"),
                terminal_key_state.get("rejected_v20_public_key", {}).get("sha256"),
                terminal_key_state.get("failed_v21_public_key", {}).get("sha256"),
                terminal_key_state.get("rejected_v22_public_key", {}).get("sha256"),
                terminal_key_state.get("rejected_v23_public_key", {}).get("sha256"),
                terminal_key_state.get("rejected_v24_public_key", {}).get("sha256"),
                terminal_key_state.get("failed_v25_public_key", {}).get("sha256"),
                terminal_key_state.get("failed_v26_public_key", {}).get("sha256"),
                terminal_key_state.get("rejected_v27_public_key", {}).get("sha256"),
                terminal_key_state.get("failed_v28_public_key", {}).get("sha256"),
                terminal_key_state.get("failed_v29_public_key", {}).get("sha256"),
                terminal_key_state.get("recovery_public_key", {}).get("sha256"),
            }
        )
        != 16
        or set(payload.get("recovery_sources", {})) != set(EXPECTED_SOURCES)
        or absence_contract.get("forbidden_output_slot_count")
        != len(FORBIDDEN_OUTPUT_SLOTS)
        or absence_contract.get("forbidden_output_slots_sha256")
        != expected_paths_sha256
        or absence_contract.get("staging_pattern_count")
        != len(CEREMONY_STAGING_PATTERNS)
        or absence_contract.get("staging_patterns_sha256")
        != expected_patterns_sha256
        or absence_contract.get("pass") is not True
    ):
        raise RuntimeError("v30 recovery static contract payload drift")
    return {
        "authority_id": payload["authority_id"],
        "prior_authority_signature_verified": True,
        "prior_failed_signed_generations_verified": [
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
        "rejected_pre_authority_generations_verified": [
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
        "terminal_key_rotation_verified": True,
        "v30_key_bootstrap_verified": True,
        "reviewed_source_count": len(payload["recovery_sources"]),
        "ceremony_absence_contract": absence_contract,
        "pass": True,
    }


def main() -> None:
    ensure_bound_self_execution()
    self_execution_binding = verify_self_execution_binding()
    before = {str(path): os.path.lexists(path) for path in FORBIDDEN_OUTPUT_SLOTS}
    staging_patterns = CEREMONY_STAGING_PATTERNS
    staging_before = {
        pattern: sorted(str(path) for path in RESULT_ROOT.glob(pattern))
        for pattern in staging_patterns
    }
    if any(before.values()):
        raise RuntimeError("A recovery/downstream output exists before read-only probe")
    if any(staging_before.values()):
        raise RuntimeError("A recovery authority staging directory exists before probe")
    reviews_before, review_state_before = classify_review_evidence(
        REVIEW_EVIDENCE_PATHS
    )
    prior_records = {}
    for role, (path, expected_size, expected_sha256) in REQUIRED_PRIOR_INPUTS.items():
        record = artifact(path)
        if (
            record["size_bytes"] != expected_size
            or record["sha256"] != expected_sha256
        ):
            raise RuntimeError(f"Required prior recovery input drift: {role}")
        prior_records[role] = record
    rejected_v15_key_archives = validate_rejected_v15_key_archives()
    rejected_v29_round1_before = validate_rejected_v29_round1_archive()
    for generation in (
        "20260722_v1",
        "20260722_v2",
        "20260722_v3",
        "20260722_v4",
        "20260722_v5",
        "20260722_v6",
        "20260722_v7",
        "20260722_v8",
        "20260722_v9",
        "20260723_v10",
        "20260723_v11",
        "20260723_v12",
        "20260723_v14",
        "20260723_v16",
        "20260723_v17",
        "20260723_v18",
        "20260723_v19",
        "20260723_v21",
    ):
        prior_private = Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery"
        ) / generation / "ed25519_private_one_time.pem"
        prior_private_stat = os.stat(prior_private, follow_symlinks=False)
        if (
            not stat.S_ISREG(prior_private_stat.st_mode)
            or stat.S_IMODE(prior_private_stat.st_mode) != 0
            or prior_private_stat.st_nlink != 1
        ):
            raise RuntimeError(f"Prior recovery private key is not retired: {generation}")
    failed_v25_private = Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
        "archive/20260723_v25_signed_c_tumor_ref_pre_audit_path_relation_failure_01/"
        "ed25519_private_one_time.pem"
    )
    failed_v25_private_stat = os.stat(failed_v25_private, follow_symlinks=False)
    if (
        not stat.S_ISREG(failed_v25_private_stat.st_mode)
        or stat.S_IMODE(failed_v25_private_stat.st_mode) != 0
        or failed_v25_private_stat.st_nlink != 1
    ):
        raise RuntimeError("Failed v25 recovery private key is not retired")
    failed_v26_private = Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
        "archive/20260724_v26_signed_c_runtime_role_projection_failure_01/"
        "ed25519_private_one_time.pem"
    )
    failed_v26_private_stat = os.stat(failed_v26_private, follow_symlinks=False)
    if (
        not stat.S_ISREG(failed_v26_private_stat.st_mode)
        or stat.S_IMODE(failed_v26_private_stat.st_mode) != 0
        or failed_v26_private_stat.st_nlink != 1
    ):
        raise RuntimeError("Failed v26 recovery private key is not retired")
    failed_v29_authority_private = (
        FAILED_V29_KEY_ARCHIVES["authority_v29"] / "ed25519_private_one_time.pem"
    )
    failed_v29_authority_private_stat = os.stat(
        failed_v29_authority_private, follow_symlinks=False
    )
    if (
        not stat.S_ISREG(failed_v29_authority_private_stat.st_mode)
        or stat.S_IMODE(failed_v29_authority_private_stat.st_mode) != 0
        or failed_v29_authority_private_stat.st_nlink != 1
        or failed_v29_authority_private_stat.st_size != 119
    ):
        raise RuntimeError("Failed v29 recovery authority private key is not retired")
    for label, path, expected_size in (
        (
            "failed v28 authority",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
                "archive/20260724_v28_signed_dataset_c_report_metric_key_order_failure_01/"
                "ed25519_private_one_time.pem"
            ),
            119,
        ),
        (
            "failed v28 result",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
                "archive/20260724_all_ssnv_result_v5_consumed_v28_c_failure_01/"
                "ed25519_private_encrypted_unrecoverable_after_signing.pem"
            ),
            290,
        ),
    ):
        opened = os.stat(path, follow_symlinks=False)
        if (
            not stat.S_ISREG(opened.st_mode)
            or stat.S_IMODE(opened.st_mode) != 0
            or opened.st_nlink != 1
            or opened.st_size != expected_size
        ):
            raise RuntimeError(f"{label} private key is not retired")
    for label, path in (
        (
            "rejected v20 recovery",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
                "20260723_v20/ed25519_private_one_time.pem"
            ),
        ),
        (
            "rejected v22 recovery",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
                "20260723_v22/ed25519_private_one_time.pem"
            ),
        ),
        (
            "rejected v23 recovery",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
                "20260723_v23/ed25519_private_one_time.pem"
            ),
        ),
        (
            "rejected v24 recovery",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
                "20260723_v24/ed25519_private_one_time.pem"
            ),
        ),
        (
            "rejected v27 archived unused recovery",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
                "archive/20260724_v27_unused_pre_authority_review_rejection_01/"
                "ed25519_private_one_time.pem"
            ),
        ),
        (
            "active v30 recovery",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
                "20260724_v30_four_role_rebootstrap/ed25519_private_one_time.pem"
            ),
        ),
    ):
        opened = os.stat(path, follow_symlinks=False)
        if (
            not stat.S_ISREG(opened.st_mode)
            or stat.S_IMODE(opened.st_mode) != 0o400
            or opened.st_nlink != 1
            or opened.st_size != 119
        ):
            raise RuntimeError(f"{label} private key is not mode 0400")
    for label, path in (
        (
            "authorized v2 terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "20260722_m2v5_terminal_v2/ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "failed v16 quarantined terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "20260723_m2v5_terminal_v6/ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "failed v17 quarantined terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "20260723_m2v5_terminal_v7/ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "failed v18 quarantined terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "20260723_m2v5_terminal_v8/ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "failed v19 quarantined terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "20260723_m2v5_terminal_v9/ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "rejected v20 quarantined terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "20260723_m2v5_terminal_v10/ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "failed v21 quarantined terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "20260723_m2v5_terminal_v11/ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "rejected v22 quarantined terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "20260723_m2v5_terminal_v12/ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "rejected v23 quarantined terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "20260723_m2v5_terminal_v13/ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "rejected v24 quarantined terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "20260723_m2v5_terminal_v14/ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "failed v25 archived unused terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "archive/20260723_m2v5_terminal_v15_unused_after_signed_v25_c_failure_01/"
                "ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "failed v26 archived unused terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "archive/20260724_m2v5_terminal_v16_unused_after_signed_v26_c_failure_01/"
                "ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "rejected v27 archived unused terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "archive/20260724_m2v5_terminal_v17_unused_pre_authority_review_rejection_01/"
                "ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "failed v28 archived unused terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "archive/20260724_m2v5_terminal_v18_unused_after_signed_v28_c_failure_01/"
                "ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "failed v29 archived terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "archive/20260724_m2v5_terminal_v19_unused_after_signed_v29_r_failure_01/"
                "ed25519_private_one_time_resident.pem"
            ),
        ),
        (
            "recovery v30 terminal",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
                "20260724_m2v5_terminal_v21/ed25519_private_one_time_resident.pem"
            ),
        ),
    ):
        opened = os.stat(path, follow_symlinks=False)
        if (
            not stat.S_ISREG(opened.st_mode)
            or stat.S_IMODE(opened.st_mode) != 0o400
            or opened.st_nlink != 1
            or opened.st_size != 119
        ):
            raise RuntimeError(f"{label} private key is not mode 0400")
    for label, path in (
        (
            "failed v28 archived unused report",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
                "archive/20260724_all_ssnv_report_v5_unused_v28_c_failure_01/"
                "ed25519_private_encrypted_unrecoverable_after_signing.pem"
            ),
        ),
        (
            "failed v29 archived result",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
                "archive/20260724_all_ssnv_result_v6_unused_v29_r_failure_01/"
                "ed25519_private_one_time.pem"
            ),
        ),
        (
            "failed v29 archived report",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
                "archive/20260724_all_ssnv_report_v6_unused_v29_r_failure_01/"
                "ed25519_private_one_time.pem"
            ),
        ),
        (
            "recovery v30 result",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
                "20260724_all_ssnv_result_v9_v30_recovery/ed25519_private_one_time.pem"
            ),
        ),
        (
            "recovery v30 report",
            Path(
                "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
                "20260724_all_ssnv_report_v9_v30_recovery/ed25519_private_one_time.pem"
            ),
        ),
    ):
        opened = os.stat(path, follow_symlinks=False)
        if (
            not stat.S_ISREG(opened.st_mode)
            or stat.S_IMODE(opened.st_mode) != 0o400
            or opened.st_nlink != 1
            or opened.st_size != 290
        ):
            raise RuntimeError(f"{label} private key is not mode 0400")
    records = {}
    self_checks = {}
    for role, (path, expected_size, expected_sha256, self_check) in EXPECTED_SOURCES.items():
        record = artifact(path)
        data = path.read_bytes()
        ast.parse(data, filename=str(path))
        compile(data, str(path), "exec")
        if (
            expected_size is not None
            and (
                record["size_bytes"] != expected_size
                or record["sha256"] != expected_sha256
            )
        ):
            raise RuntimeError(f"Frozen recovery source identity drift: {role}")
        records[role] = record
        if self_check is not None:
            self_checks[role] = recursive_self_check(path, self_check)
    if self_checks["continuation_verifier"] != "PASS":
        raise RuntimeError("Recovery verifier recursive source binding failed")
    if self_checks["runner_gate_replay"] != "PASS":
        raise RuntimeError("Recovery replayer recursive source binding failed")
    if self_checks["downstream_continuation"].get("pass") is not True:
        raise RuntimeError("Recovery continuation recursive source binding failed")

    archived_state_runner_replay = run_archived_state_runner_replay()
    tests = run_regression_tests()
    pilot = run_pilot_verifier()
    prior_recorded_verifier = run_prior_recorded_verifier_recheck()
    static_contract = run_static_contract()
    rejected_v29_round1_after = validate_rejected_v29_round1_archive()
    after = {str(path): os.path.lexists(path) for path in FORBIDDEN_OUTPUT_SLOTS}
    staging_after = {
        pattern: sorted(str(path) for path in RESULT_ROOT.glob(pattern))
        for pattern in staging_patterns
    }
    reviews_after, review_state_after = classify_review_evidence(REVIEW_EVIDENCE_PATHS)
    if before != after or any(after.values()) or staging_before != staging_after:
        raise RuntimeError("Read-only probe created a recovery/downstream output")
    if reviews_before != reviews_after:
        raise RuntimeError("Read-only probe changed review evidence state")
    if review_state_before != review_state_after:
        raise RuntimeError("Read-only probe changed review evidence classification")
    if rejected_v29_round1_before != rejected_v29_round1_after:
        raise RuntimeError("Read-only probe changed rejected v29 round1 archive")
    print(
        json.dumps(
            {
                "schema_name": "intersubmod.tumor_ref_schema_recovery_source_probe",
                "schema_version": "1.0.0",
                "task_type": "B_comprehensive_validation",
                "sources": records,
                "self_execution_binding": self_execution_binding,
                "prior_inputs": prior_records,
                "prior_inputs_verified": True,
                "rejected_v15_key_archives": rejected_v15_key_archives,
                "rejected_v29_round1_archive": rejected_v29_round1_after,
                "rejected_v29_round1_verified": True,
                "recursive_self_checks": self_checks,
                "archived_state_runner_replay": archived_state_runner_replay,
                "regression_tests": tests,
                "historical_eight_key_pilot": pilot,
                "prior_recorded_verifier_recheck": prior_recorded_verifier,
                "static_contract": static_contract,
                "forbidden_output_slots_checked": len(FORBIDDEN_OUTPUT_SLOTS),
                "review_evidence_state": review_state_after,
                "no_output_writes": True,
                "no_output_writes_scope": (
                    "protected_recovery_and_downstream_namespaces_only"
                ),
                "pass": True,
            },
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
