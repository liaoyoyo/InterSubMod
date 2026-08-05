#!/usr/bin/env python3
"""Build the fail-closed machine dataset for the final all-sSNV report.

This program reconciles the full screen, stable read assignments,
cooccurrence results, strict confirmation, tumor-REF controls, and matched
normal controls. It emits data products only; it does not write report prose
or figures.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence, TextIO

from scipy.stats import beta, binomtest


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
TUMOR_REF_SOURCE_IDENTITY_VERIFIER = (
    SCRIPT_DIR / "verify_retrospective_running_source_identity_v2.py"
)
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import m2_screen_gate as M2_GATE  # noqa: E402
import audit_stable_primary_artifacts as PRIMARY_AUDITOR  # noqa: E402
import finalize_cooccurrence_release_receipt as RELEASE_FINALIZER  # noqa: E402
import release_source_authority as SOURCE_AUTHORITY  # noqa: E402
import run_matched_normal_candidate_controls as MATCHED_NORMAL_RUNNER  # noqa: E402
import analyze_matched_normal_candidate_controls as MATCHED_NORMAL_ANALYZER  # noqa: E402
import annotate_candidate_cn_ccf as CN_CCF_ANNOTATOR  # noqa: E402
import run_strict_methyl_candidate_confirmation as STRICT_PRODUCER  # noqa: E402


SCHEMA_VERSION = "2.0.0"
INPUT_MANIFEST_SCHEMA_VERSION = "1.0.0"
SCREEN_OUTPUT_SCHEMA_VERSION = "1.2.0"
ASSIGNMENT_SCHEMA_VERSION = "1.0.0"
TUMOR_REF_SCHEMA_VERSION = "1.0.0"
TUMOR_REF_SOURCE_IDENTITY_SCHEMA_VERSION = "1.2.0"
MATCHED_NORMAL_ANALYSIS_SCHEMA_VERSION = "2.0.0"
MATCHED_NORMAL_PAIRED_RUNNER_SCHEMA_VERSION = "1.0.0"
COOCCURRENCE_SCHEMA_VERSION = "2.0.0"
COOCCURRENCE_PREFLIGHT_SCHEMA_VERSION = "2.0.0"
COOCCURRENCE_RELEASE_SCHEMA_VERSION = "1.0.0"
COOCCURRENCE_PARTNER_WINDOW_BP = RELEASE_FINALIZER.ANALYZER.PAIR_WINDOW_BP
ANALYSIS_SCOPE_IDENTITY_CONTRACT = (
    "selected_focal_ref_alt_alignments_full_sam_core_all_typed_aux_except_RG_"
    "latest_HP_PS_assignment_label_and_partner_calls_v1"
)
RAW_IDENTITY_RELEASE_CONTRACT = (
    "typed_aux_dependency_attested_sparse_audit_analysis_payload_v3"
)
RAW_DUPLICATE_EQUIVALENCE_CONTRACT = "sam_core_and_all_aux_tags_except_RG_exact_v1"
RAW_IDENTITY_MISSING_POLICY = "hard_fail_before_site_result"
RAW_IDENTITY_CONFLICT_POLICY = "hard_fail_before_site_result"
STRICT_FORMAL_SCHEMA_VERSION = "3.1.0"
STRICT_SCHEMA_VERSIONS = {STRICT_FORMAL_SCHEMA_VERSION}
CN_CCF_ANNOTATION_SCHEMA_VERSION = "1.0.0"
PRIMARY_ARTIFACT_AUDIT_SCHEMA_VERSION = "2.0.0"
INDEPENDENT_PRIMARY_ARTIFACT_PATHS = {
    "reads": Path("reads") / "reads.tsv",
    "distance_matrix": Path("distance") / "BERNOULLI" / "matrix.csv",
    "methylation_matrix": Path("methylation") / "methylation.csv",
}
INDEPENDENT_PRIMARY_RECOUNT_WORKERS = 40
INDEPENDENT_PRIMARY_RECOUNT_CHUNK_SIZE = 64
JOINT_SIGNATURE_P_MAX = 0.05
COMPLETE_MARKER_CRAMERS_V_MIN = 0.30
COMPLETE_MARKER_DELTA_ALT_MIN = 0.50
PAIR_GLOBAL_Q_MAX = 0.05
PAIR_CONDITIONAL_P_MAX = 0.05
PAIR_CONDITIONAL_PERMUTATIONS = 999
FOUR_STATE_ERROR_CEILING = 0.02
FOUR_STATE_FAMILYWISE_CONFIDENCE = 0.95
FOUR_STATE_RELATION_FAMILY_SIZE = 3
FOUR_STATE_PER_RELATION_CONFIDENCE = 1.0 - (
    1.0 - FOUR_STATE_FAMILYWISE_CONFIDENCE
) / FOUR_STATE_RELATION_FAMILY_SIZE
FOUR_STATE_MULTIPLICITY_METHOD = "bonferroni_three_relation_models"
FOUR_STATE_MIN_SPLIT_COUNT = 3
EXPECTED_SCREEN_SITES = 469_849
BACKGROUND_CONTROL_REPLICATION_GATE_CONTRACT = (
    "lenient_coarse_modal_K2_without_membership_ARI_requirement_v1"
)
BACKGROUND_CONTROL_RELATION_TO_PRIMARY_M1 = (
    "lenient_predicate_superset_of_ARI_qualified_predicate_on_same_background_payload"
)
TUMOR_REF_SOURCE_IDENTITY_REQUIRED_ROLES = frozenset(
    {"analyzer", "focal_alt_cluster_lib"}
)
COOCCURRENCE_CODE_PATHS = {
    "analyzer": SCRIPT_DIR / "analyze_methyl_ssnv_cooccurrence.py",
    "ssnv_cooccurrence_lib": REPO_ROOT / "scripts" / "ssnv_cooccurrence_lib.py",
    "latest_tag_join": SCRIPT_DIR / "latest_tag_join.py",
    "m2_screen_gate": SCRIPT_DIR / "m2_screen_gate.py",
    "release_finalizer": SCRIPT_DIR / "finalize_cooccurrence_release_receipt.py",
    "release_runner": SCRIPT_DIR / "run_cooccurrence_v6_source_locked.sh",
    "source_authority_validator": SCRIPT_DIR / "release_source_authority.py",
}
COOCCURRENCE_PREFLIGHT_SOURCE = SCRIPT_DIR / "audit_cooccurrence_task_contract_preflight.py"
COOCCURRENCE_RELEASE_CODE_PATHS = {
    "release_finalizer": SCRIPT_DIR / "finalize_cooccurrence_release_receipt.py",
    "release_runner": SCRIPT_DIR / "run_cooccurrence_v6_source_locked.sh",
}
RAW_IDENTITY_DUPLICATE_COLUMNS = (
    "sample",
    "chrom",
    "pos",
    "projection_read_name",
    "projection_chrom",
    "projection_start",
    "projection_end",
    "projection_mapq",
    "projection_strand",
    "record_count",
    "classification",
    "differing_auxiliary_tags",
)
COOCCURRENCE_ORACLE_FOCAL_IDS = frozenset(
    {
        "H2009_chr3_193395128",
        "HCC1395_DORADO_chr5_750311",
        "COLO829_chr4_92183865",
    }
)
COOCCURRENCE_SHARED_READSET_ORACLE_ID = (
    "COLO829_chr20_52761564_52761565_shared_readset"
)
TUMOR_REF_SOURCE_IDENTITY_TRUE_CHECKS = frozenset(
    {
        "snapshot_pass",
        "producer_manifest_pass",
        "producer_processed_all_102842_sites",
        "source_role_sets_exact",
        "snapshot_creator_source_unchanged",
        "post_run_verifier_source_identity_recorded",
        "process_start_within_manifest_start_clock_tolerance",
        "snapshot_observed_between_manifest_start_and_finish",
        "all_sources_predate_manifest_start",
        "all_source_identities_unchanged_after_snapshot",
        "manifest_source_artifacts_match_snapshot",
        "manifest_command_exactly_matches_live_process_command",
        "manifest_script_token_matches_attested_source",
        "snapshot_and_manifest_parsed_from_hash_bound_bytes",
        "manifest_hash_bound_to_receipt",
    }
)
MANIFEST_HASH_REQUIRED_FIELDS = frozenset(
    {
        "all_ssnv_vcf",
        "all_ssnv_vcf_index",
        "latest_read_tag_sidecar_index",
        "site_ledger",
        "site_ledger_index",
        "layered_region_view",
        "truth_fp_vcf",
        "truth_fp_vcf_index",
        "truth_tp_vcf",
        "truth_tp_vcf_index",
    }
)
MANIFEST_LARGE_METADATA_IDENTITY_FIELDS = frozenset(
    {"raw_alignment", "raw_alignment_index", "latest_read_tag_sidecar"}
)
DATASETS = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
TRUTH_LABELS = ("TP", "FP", "UNASSESSED")
BIOLOGICAL_IDS = {
    "HCC1395": "HCC1395",
    "HCC1395_DORADO": "HCC1395",
    "COLO829": "COLO829",
    "H1437": "H1437",
    "H2009": "H2009",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
}
FORMAL_SELECTION_COLUMN = "multi_marker_molecular_haplotype_base_candidate"
DEFAULT_SELECTION_COLUMN = "genetically_anchored_multi_marker_candidate"
BY_SELECTION_COLUMN = "genetically_anchored_multi_marker_candidate_by_sensitivity"
STRICT_SELECTION_COLUMN = FORMAL_SELECTION_COLUMN
ASSIGNMENT_SCHEMA = "intersubmod.all_ssnv_stable_multigroup_read_assignment"
PASS_SEMANTICS = "integration_integrity_only_not_scientific_confirmation"
STRICT_PASS_SEMANTICS = "execution_integrity_only_not_scientific_confirmation"
CANONICAL_CLAIM_CONTRACT = SOURCE_AUTHORITY.CLAIM_CONTRACT_PATH
CANONICAL_WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
CANONICAL_RESULT_ROOT = SCRIPT_DIR.parent / "results"
CANONICAL_TASK_B_PATHS = {
    "manifest": CANONICAL_RESULT_ROOT / "all_ssnv_input_manifest.json",
    "screen_dir": (
        CANONICAL_WORKSPACE_ROOT
        / "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
    ),
    "cooccurrence_dir": (
        CANONICAL_WORKSPACE_ROOT
        / "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity"
    ),
    "strict_dir": (
        CANONICAL_WORKSPACE_ROOT
        / "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5"
    ),
    "tumor_ref_dir": (
        CANONICAL_WORKSPACE_ROOT
        / "all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel"
    ),
    "primary_artifact_audit_pre": (
        CANONICAL_RESULT_ROOT
        / "stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream.json"
    ),
    "primary_artifact_audit_post": (
        CANONICAL_RESULT_ROOT
        / "stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json"
    ),
    "cooccurrence_preflight": (
        CANONICAL_RESULT_ROOT
        / "cooccurrence_task_contract_preflight.v9_command_parity_full_runtime.json"
    ),
    "tumor_ref_source_identity_receipt": (
        CANONICAL_WORKSPACE_ROOT
        / "tumor_ref_recovery_source_identity_v1"
        / "post_run_source_identity.receipt.json"
    ),
    "independent_m2_audit": CANONICAL_RESULT_ROOT / "independent_m2_gate_recount.v3.json",
    "matched_normal_dir": (
        CANONICAL_WORKSPACE_ROOT
        / "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5"
    ),
    "cn_ccf_annotations": (
        CANONICAL_WORKSPACE_ROOT
        / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5"
    ),
}
CANONICAL_MATCHED_NORMAL_DIRS = (
    CANONICAL_TASK_B_PATHS["matched_normal_dir"],
    CANONICAL_WORKSPACE_ROOT
    / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5",
)
CANONICAL_FINAL_DATASET_DIR = (
    CANONICAL_WORKSPACE_ROOT / "all_ssnv_final_report_dataset_v5_m2v5_source_attested"
)
CLAIM_IDS = ("M1", "M2", "G1", "G2", "R1", "B1", "C1", "L1", "L2")
CLAIM_STATUSES = {"PASS", "FAIL", "NOT_EVALUABLE", "NOT_RUN"}
FOUR_STATE_RELATIONS = frozenset(
    {
        "NOT_IDENTIFIABLE_INSUFFICIENT_FOUR_STATE_DEPTH",
        "FOCAL_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
        "PARTNER_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
        "BRANCHING_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
        "MULTIPLE_MUTATION_ORDER_MODELS_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
        "NOT_IDENTIFIABLE_FIXED_ERROR_CEILING",
        "INCOMPATIBLE_OR_COMPLEX_UNDER_FIXED_ERROR_MODEL",
    }
)
COMPATIBLE_RELATIONS = frozenset(
    {
        "FOCAL_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
        "PARTNER_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
        "BRANCHING_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
        "MULTIPLE_MUTATION_ORDER_MODELS_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
    }
)
UNRESOLVED_FOUR_STATE_RELATIONS = frozenset(
    {
        "NOT_IDENTIFIABLE_INSUFFICIENT_FOUR_STATE_DEPTH",
        "NOT_IDENTIFIABLE_FIXED_ERROR_CEILING",
    }
)

CN_CCF_OUTPUT_COLUMNS = (
    "sample",
    "chrom",
    "pos",
    "ref",
    "alt",
    "mutation_id",
    "callset_status",
    "cn_status",
    "cn_authority_sample",
    "cn_transfer_policy",
    "independent_cn",
    "cn_segment_id",
    "cn_segment_start",
    "cn_segment_end",
    "savana_total_cn_raw",
    "savana_minor_cn_raw",
    "savana_total_cn_discrete",
    "savana_major_cn_discrete",
    "savana_minor_cn_discrete",
    "purity_model_value",
    "pyclone_status",
    "pyclone_primary_bundle_id",
    "pyclone_primary_model_role",
    "pyclone_fit_local_cluster_id",
    "pyclone_vi_cellular_prevalence",
    "pyclone_vi_cellular_prevalence_std",
    "pyclone_vi_assignment_probability",
    "pyclone_sensitivity_status",
    "pyclone_sensitivity_bundle_id",
    "pyclone_sensitivity_fit_local_cluster_id",
    "pyclone_sensitivity_cellular_prevalence",
    "pyclone_sensitivity_cellular_prevalence_std",
    "pyclone_sensitivity_assignment_probability",
    "authority_config_sha256",
    "input_provenance_sha256",
    "analysis_summary_sha256",
    "cn_source_sha256",
    "pyclone_primary_metadata_sha256",
    "pyclone_primary_results_sha256",
    "pyclone_primary_status_sha256",
    "pyclone_sensitivity_metadata_sha256",
    "pyclone_sensitivity_results_sha256",
    "pyclone_sensitivity_status_sha256",
    "claim_ceiling",
)
CN_CCF_NATIVE_DETAIL_FIELDS = tuple(
    field_name
    for field_name in CN_CCF_OUTPUT_COLUMNS
    if field_name not in {"sample", "chrom", "pos", "ref", "alt"}
)
CN_CCF_CLAIM_CEILING = "conditional_annotation_only_not_orthogonal_confirmation"
CN_CCF_PASS_SEMANTICS = "execution_integrity_only_not_orthogonal_confirmation"
CN_STATUS_ENUM = frozenset(
    {
        "AVAILABLE_EXACT_SEGMENT",
        "SHARED_CN_SENSITIVITY",
        "BLOCKED_CN_MISFIT",
        "NO_CN_SEGMENT",
        "RECEIPT_FAIL",
    }
)
PYCLONE_STATUS_ENUM = frozenset(
    {
        "MATCHED_PRIMARY",
        "MATCHED_SENSITIVITY_ONLY",
        "NOT_IN_FIT_UNIVERSE",
        "BLOCKED_CN_MISFIT",
        "RECEIPT_FAIL",
    }
)

CandidateKey = tuple[str, str, int, str, str]
BiologicalKey = tuple[str, str, int, str, str]
PairKey = tuple[str, str, int, str, str, str, int, str, str]
CrossPlatformPairKey = tuple[str, int, str, str, str, int, str, str]


class ContractError(RuntimeError):
    """Raised when independently produced machine outputs cannot be reconciled."""


@dataclass(frozen=True)
class InputBundle:
    manifest: Path
    screen_dir: Path
    cooccurrence_dir: Path
    strict_dir: Path
    tumor_ref_dir: Path
    primary_artifact_audit_pre: Path
    primary_artifact_audit_post: Path
    cooccurrence_preflight: Path
    tumor_ref_source_identity_receipt: Path | None = None
    independent_m2_audit: Path | None = None
    matched_normal_dir: Path | None = None
    cn_ccf_annotations: Path | None = None

    @property
    def screen_sites(self) -> Path:
        return self.screen_dir / "all_ssnv_site_results.tsv.gz"

    @property
    def screen_assignments(self) -> Path:
        return self.screen_dir / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"

    @property
    def screen_summary(self) -> Path:
        return self.screen_dir / "all_ssnv_summary.json"

    @property
    def screen_receipt(self) -> Path:
        return self.screen_dir / "run_manifest.json"

    @property
    def cooccurrence_sites(self) -> Path:
        return self.cooccurrence_dir / "methyl_ssnv_site_results.tsv.gz"

    @property
    def cooccurrence_pairs(self) -> Path:
        return self.cooccurrence_dir / "methyl_ssnv_pair_results.tsv.gz"

    @property
    def cooccurrence_summary(self) -> Path:
        return self.cooccurrence_dir / "summary.json"

    @property
    def cooccurrence_receipt(self) -> Path:
        return self.cooccurrence_dir / "run_receipt.json"

    @property
    def cooccurrence_release_receipt(self) -> Path:
        return self.cooccurrence_dir / "release_receipt.json"

    @property
    def cooccurrence_raw_identity_duplicates(self) -> Path:
        return self.cooccurrence_dir / "raw_identity_duplicate_audit.tsv.gz"

    @property
    def cooccurrence_oracle_cases(self) -> Path:
        return self.cooccurrence_dir / "oracle_cases.json"

    @property
    def strict_sites(self) -> Path:
        return self.strict_dir / "strict_methyl_candidate_confirmation_sites.tsv.gz"

    @property
    def strict_summary(self) -> Path:
        return self.strict_dir / "strict_methyl_candidate_confirmation_summary.json"

    @property
    def strict_receipt(self) -> Path:
        return self.strict_dir / "run_manifest.json"

    @property
    def tumor_ref_sites(self) -> Path:
        return self.tumor_ref_dir / "all_ssnv_tumor_ref_control_site_results.tsv.gz"

    @property
    def tumor_ref_summary(self) -> Path:
        return self.tumor_ref_dir / "all_ssnv_tumor_ref_control_summary.json"

    @property
    def tumor_ref_receipt(self) -> Path:
        return self.tumor_ref_dir / "run_manifest.json"

    @property
    def matched_normal_sites(self) -> Path | None:
        if self.matched_normal_dir is None:
            return None
        return self.matched_normal_dir / "matched_normal_candidate_controls.tsv.gz"

    @property
    def matched_normal_summary(self) -> Path | None:
        if self.matched_normal_dir is None:
            return None
        return self.matched_normal_dir / "matched_normal_candidate_controls_summary.json"

    @property
    def matched_normal_receipt(self) -> Path | None:
        if self.matched_normal_dir is None:
            return None
        return self.matched_normal_dir / "run_receipt.json"


@dataclass
class ScreenStageCounts:
    screened: int = 0
    evaluable: int = 0
    stable: int = 0
    robust: int = 0


@dataclass
class LatestTagAudit:
    n_sites: int = 0
    join_status_counts: Counter[str] = field(default_factory=Counter)
    n_reads_tsv_site_rows: int = 0
    n_exact_hp_ps_site_read_joins: int = 0
    n_ps_present_site_read_joins: int = 0
    n_source_hp_replaced_site_read_joins: int = 0
    n_sidecar_rows_fetched_site_observations: int = 0
    n_sidecar_rows_eligible_site_observations: int = 0
    n_projection_multimatch_site_reads: int = 0
    n_sites_with_zero_reads_tsv_rows: int = 0

    def add(self, row: Mapping[str, Any], *, source: str) -> None:
        status = str(row["latest_tag_join_status"])
        field_map = {
            "latest_tag_rows_fetched": "n_sidecar_rows_fetched_site_observations",
            "latest_tag_rows_eligible": "n_sidecar_rows_eligible_site_observations",
            "latest_tag_reads_joined": "n_exact_hp_ps_site_read_joins",
            "latest_tag_ps_present": "n_ps_present_site_read_joins",
            "source_hp_replaced_reads": "n_source_hp_replaced_site_read_joins",
            "latest_tag_projection_multimatch_reads": "n_projection_multimatch_site_reads",
            "n_reads_total": "n_reads_tsv_site_rows",
        }
        counts: dict[str, int] = {}
        for input_field in field_map:
            value = optional_int(row[input_field], field_name=input_field)
            if value is None or value < 0:
                raise ContractError(
                    f"Screen latest HP/PS count must be nonnegative in {source}: "
                    f"{input_field}={row[input_field]!r}"
                )
            counts[input_field] = value
        joined = counts["latest_tag_reads_joined"]
        if status != "PASS":
            raise ContractError(f"Screen latest HP/PS join status is not PASS in {source}: {status!r}")
        if joined != counts["n_reads_total"]:
            raise ContractError(
                f"Screen latest HP/PS joined/read-row mismatch in {source}: "
                f"joined={joined} reads={counts['n_reads_total']}"
            )
        if not (
            counts["latest_tag_rows_fetched"]
            >= counts["latest_tag_rows_eligible"]
            >= joined
        ):
            raise ContractError(f"Screen latest HP/PS fetched/eligible/joined drift in {source}")
        for field_name in (
            "latest_tag_ps_present",
            "source_hp_replaced_reads",
            "latest_tag_projection_multimatch_reads",
        ):
            if counts[field_name] > joined:
                raise ContractError(
                    f"Screen latest HP/PS count exceeds joined reads in {source}: {field_name}"
                )
        self.n_sites += 1
        self.join_status_counts[status] += 1
        for input_field, attribute in field_map.items():
            setattr(self, attribute, getattr(self, attribute) + counts[input_field])
        self.n_sites_with_zero_reads_tsv_rows += counts["n_reads_total"] == 0

    def payload(self) -> dict[str, Any]:
        all_sites_pass = self.join_status_counts == Counter({"PASS": self.n_sites})
        every_row_joined = (
            self.n_exact_hp_ps_site_read_joins == self.n_reads_tsv_site_rows
        )
        projection_unique = self.n_projection_multimatch_site_reads == 0
        return {
            "authoritative_tag_source": "same_run_LongPhase_S_external_HP_PS_sidecar",
            "embedded_reads_tsv_hp_used_for_analysis": False,
            "join_occurs_before_focal_ALT_selection": True,
            "counting_unit": "site_read_observation_not_globally_unique_read",
            "n_sites": self.n_sites,
            "join_status_counts": dict(self.join_status_counts),
            "all_sites_pass": all_sites_pass,
            "n_reads_tsv_site_rows": self.n_reads_tsv_site_rows,
            "n_exact_hp_ps_site_read_joins": self.n_exact_hp_ps_site_read_joins,
            "every_reads_tsv_row_joined": every_row_joined,
            "n_ps_present_site_read_joins": self.n_ps_present_site_read_joins,
            "n_source_hp_replaced_site_read_joins": (
                self.n_source_hp_replaced_site_read_joins
            ),
            "n_sidecar_rows_fetched_site_observations": (
                self.n_sidecar_rows_fetched_site_observations
            ),
            "n_sidecar_rows_eligible_site_observations": (
                self.n_sidecar_rows_eligible_site_observations
            ),
            "n_projection_multimatch_site_reads": self.n_projection_multimatch_site_reads,
            "all_projection_identities_unique": projection_unique,
            "n_sites_with_zero_reads_tsv_rows": self.n_sites_with_zero_reads_tsv_rows,
            "pass": all_sites_pass and every_row_joined and projection_unique,
        }


@dataclass
class ScreenData:
    all_keys: set[CandidateKey]
    all_rows: dict[CandidateKey, dict[str, str]]
    stable_rows: dict[CandidateKey, dict[str, str]]
    stage_counts: dict[tuple[str, str], ScreenStageCounts]
    per_sample_truth: dict[str, Counter[str]]
    claim_dataset_sets: dict[str, set[CandidateKey]]
    claim_biological_sets: dict[str, set[BiologicalKey]]
    latest_tag_audits: dict[tuple[str, str], LatestTagAudit]


@dataclass
class PairAggregate:
    n_pairs: int = 0
    n_testable_pairs: int = 0
    n_confirmed_pairs: int = 0
    n_confirmed_pairs_by: int = 0
    n_confirmed_complete_four_state_pairs: int = 0
    n_confirmed_compatible_relation_pairs: int = 0
    n_confirmed_topology_inferable_pairs: int = 0
    n_cross_platform_exact_pair_rows: int = 0
    n_cross_platform_both_confirmed_pair_rows: int = 0
    n_exact_testable_pairs: int = 0
    n_formal_pair_by_confirmed: int = 0
    n_same_pair_four_state_witnesses: int = 0


@dataclass
class PairData:
    aggregates: dict[CandidateKey, PairAggregate]
    n_rows: int
    n_testable: int
    n_confirmed: int
    n_confirmed_by: int
    exact_cross_platform_pairs: set[CrossPlatformPairKey]
    both_confirmed_cross_platform_pairs: set[CrossPlatformPairKey]
    replicated_cross_platform_pairs: set[CrossPlatformPairKey]
    rows: dict[PairKey, dict[str, str]]
    by_focal: dict[CandidateKey, list[dict[str, str]]]


@dataclass
class IntegrationData:
    manifest: dict[str, Any]
    screen: ScreenData
    cooccurrence_rows: dict[CandidateKey, dict[str, str]]
    pair_data: PairData
    strict_rows: dict[CandidateKey, dict[str, str]]
    strict_status: dict[str, Any]
    tumor_ref_rows: dict[CandidateKey, dict[str, str]]
    matched_normal_rows: dict[CandidateKey, dict[str, str]]
    matched_normal_status: dict[str, Any]
    cn_ccf_rows: dict[CandidateKey, dict[str, Any]] = field(default_factory=dict)
    cn_ccf_status: dict[str, Any] = field(default_factory=dict)
    candidate_pair_details: list[dict[str, Any]] = field(default_factory=list)
    candidates: list[dict[str, Any]] = field(default_factory=list)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def current_release_source_authority() -> dict[str, Any]:
    return SOURCE_AUTHORITY.validate_release_source_authority(
        SOURCE_AUTHORITY.EXPECTED_SOURCE_PATHS
    )


def validate_release_bound_producer(
    receipt: Mapping[str, Any],
    *,
    receipt_path: Path,
    producer_module: Any,
    producer_role: str,
    expected_command: Sequence[str],
    label: str,
) -> Path:
    source_authority = current_release_source_authority()
    expected_code = producer_module.capture_source_identity()
    expected_modes = producer_module.capture_source_modes()
    source_lock = receipt.get("source_lock")
    if receipt.get("source_authority") != source_authority:
        raise ContractError(f"{label} source authority drift")
    if receipt.get("code") != expected_code:
        raise ContractError(f"{label} producer source identity drift")
    if set(expected_code) != {producer_role}:
        raise ContractError(f"{label} producer role drift")
    if list(receipt.get("command") or []) != list(expected_command):
        raise ContractError(f"{label} producer command is not canonical")
    if (
        not isinstance(source_lock, Mapping)
        or source_lock.get("source_identity_before") != expected_code
        or source_lock.get("source_identity_after_compute") != expected_code
        or source_lock.get("source_modes_before") != expected_modes
        or source_lock.get("source_modes_after_compute") != expected_modes
        or source_lock.get("all_sources_read_only_and_unchanged") is not True
        or expected_modes != {producer_role: "0o444"}
    ):
        raise ContractError(f"{label} producer source lock drift")
    if oct(receipt_path.resolve(strict=True).stat().st_mode & 0o777) != "0o444":
        raise ContractError(f"{label} terminal receipt is not mode 0444")
    return Path(producer_module.__file__).resolve()


def validate_canonical_task_b_paths(bundle: InputBundle, output_dir: Path) -> None:
    for field, expected in CANONICAL_TASK_B_PATHS.items():
        observed = getattr(bundle, field)
        if field == "matched_normal_dir":
            allowed = {path.resolve() for path in CANONICAL_MATCHED_NORMAL_DIRS}
            if observed is None or Path(observed).resolve() not in allowed:
                raise ContractError(
                    "Task-B canonical path drift for matched_normal_dir: "
                    f"observed={observed} "
                    f"expected_one_of={sorted(str(path) for path in allowed)}"
                )
            continue
        if observed is None or Path(observed).resolve() != expected.resolve():
            raise ContractError(
                f"Task-B canonical path drift for {field}: "
                f"observed={observed} expected={expected}"
            )
    if output_dir.resolve() != CANONICAL_FINAL_DATASET_DIR.resolve():
        raise ContractError(
            "Task-B final dataset output path is not canonical: "
            f"observed={output_dir.resolve()} expected={CANONICAL_FINAL_DATASET_DIR.resolve()}"
        )


def canonical_task_b_final_builder_command(
    bundle: InputBundle, output_dir: Path
) -> list[str]:
    optional_paths = {
        "tumor_ref_source_identity_receipt": bundle.tumor_ref_source_identity_receipt,
        "independent_m2_audit": bundle.independent_m2_audit,
        "matched_normal_dir": bundle.matched_normal_dir,
        "cn_ccf_annotations": bundle.cn_ccf_annotations,
    }
    if any(path is None for path in optional_paths.values()):
        missing = sorted(name for name, path in optional_paths.items() if path is None)
        raise ContractError(f"Canonical Task-B builder command lacks inputs: {missing}")
    return [
        *STRICT_PRODUCER.canonical_python_prefix(bundle.strict_dir),
        str(Path(__file__).resolve()),
        "--manifest",
        str(bundle.manifest.resolve()),
        "--screen-dir",
        str(bundle.screen_dir.resolve()),
        "--cooccurrence-dir",
        str(bundle.cooccurrence_dir.resolve()),
        "--strict-dir",
        str(bundle.strict_dir.resolve()),
        "--tumor-ref-dir",
        str(bundle.tumor_ref_dir.resolve()),
        "--primary-artifact-audit-pre",
        str(bundle.primary_artifact_audit_pre.resolve()),
        "--primary-artifact-audit-post",
        str(bundle.primary_artifact_audit_post.resolve()),
        "--cooccurrence-preflight",
        str(bundle.cooccurrence_preflight.resolve()),
        "--tumor-ref-source-identity-receipt",
        str(bundle.tumor_ref_source_identity_receipt.resolve()),
        "--independent-m2-audit",
        str(bundle.independent_m2_audit.resolve()),
        "--matched-normal-dir",
        str(bundle.matched_normal_dir.resolve()),
        "--cn-ccf-annotations",
        str(bundle.cn_ccf_annotations.resolve()),
        "--output-dir",
        str(output_dir.resolve()),
    ]


def canonical_task_b_final_builder_commands(output_dir: Path) -> list[list[str]]:
    commands: list[list[str]] = []
    for matched_normal_dir in CANONICAL_MATCHED_NORMAL_DIRS:
        paths = dict(CANONICAL_TASK_B_PATHS)
        paths["matched_normal_dir"] = matched_normal_dir
        commands.append(
            canonical_task_b_final_builder_command(
                InputBundle(**paths),
                output_dir,
            )
        )
    return commands


def observed_process_command() -> list[str]:
    raw = Path("/proc/self/cmdline").read_bytes()
    if not raw.endswith(b"\0"):
        raise ContractError("Process command line is unavailable or malformed")
    return [os.fsdecode(token) for token in raw[:-1].split(b"\0")]


def source_script_token_matches_attested_source(token: str, source_path: str) -> bool:
    """Validate an absolute or exact repository-relative producer script token."""
    if not token or "\x00" in token:
        return False
    token_path = Path(token)
    source = Path(source_path)
    if token_path.is_absolute():
        try:
            return token_path.resolve(strict=True) == source.resolve(strict=True)
        except FileNotFoundError:
            return False
    raw_parts = token.split("/")
    if any(part in {"", ".", ".."} for part in raw_parts) or not token_path.parts:
        return False
    try:
        expected_relative = source.resolve(strict=True).relative_to(REPO_ROOT)
    except (FileNotFoundError, ValueError):
        return False
    return token_path.as_posix() == expected_relative.as_posix()


def require_nonempty_file(path: Path, label: str) -> Path:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ContractError(f"Missing or empty {label}: {path}")
    return path


def open_text(path: Path) -> TextIO:
    require_nonempty_file(path, "input table")
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def load_json(path: Path, label: str) -> dict[str, Any]:
    require_nonempty_file(path, label)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise ContractError(f"Invalid JSON in {label}: {path}") from error
    if not isinstance(payload, dict):
        raise ContractError(f"{label} is not a JSON object: {path}")
    return payload


def require_schema_pass(
    payload: Mapping[str, Any],
    schema: str,
    label: str,
    *,
    pass_required: bool = True,
    schema_version: str = SCHEMA_VERSION,
) -> None:
    if payload.get("schema_name") != schema:
        raise ContractError(
            f"Unexpected {label} schema: {payload.get('schema_name')!r} != {schema!r}"
        )
    if payload.get("schema_version") != schema_version:
        raise ContractError(
            f"Unexpected {label} schema version: {payload.get('schema_version')!r}"
        )
    if pass_required and payload.get("pass") is not True:
        raise ContractError(f"{label} is not passing")


def require_schema_pass_one_of(
    payload: Mapping[str, Any],
    schema: str,
    label: str,
    *,
    schema_versions: set[str],
    pass_required: bool = True,
) -> str:
    if payload.get("schema_name") != schema:
        raise ContractError(
            f"Unexpected {label} schema: {payload.get('schema_name')!r} != {schema!r}"
        )
    observed = str(payload.get("schema_version", ""))
    if observed not in schema_versions:
        raise ContractError(
            f"Unexpected {label} schema version: {observed!r}; expected one of {sorted(schema_versions)}"
        )
    if pass_required and payload.get("pass") is not True:
        raise ContractError(f"{label} is not passing")
    return observed


def verify_declared_artifact(
    reference: Any,
    path: Path,
    label: str,
    *,
    require_sha256: bool = True,
) -> None:
    if not isinstance(reference, Mapping):
        raise ContractError(f"{label} receipt lacks an artifact object")
    declared_path = reference.get("path")
    if not declared_path or Path(str(declared_path)).resolve() != path.resolve():
        raise ContractError(
            f"{label} path mismatch: declared={declared_path!r} observed={str(path.resolve())!r}"
        )
    if "size_bytes" in reference and int(reference["size_bytes"]) != path.stat().st_size:
        raise ContractError(f"{label} size mismatch")
    if require_sha256 and not reference.get("sha256"):
        raise ContractError(f"{label} receipt lacks SHA-256")
    if "sha256" in reference and str(reference["sha256"]) != sha256(path):
        raise ContractError(f"{label} SHA-256 mismatch")


def verify_declared_path_artifact(reference: Any, label: str) -> Path:
    """Verify a receipt artifact when its path is not otherwise fixed by the bundle."""
    if not isinstance(reference, Mapping):
        raise ContractError(f"{label} receipt lacks an artifact object")
    declared_path = reference.get("path")
    if not declared_path:
        raise ContractError(f"{label} receipt lacks path")
    path = Path(str(declared_path)).expanduser().resolve()
    require_nonempty_file(path, label)
    verify_declared_artifact(reference, path, label)
    return path


def parse_bool(value: Any, *, field_name: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes", "y", "t"}:
        return True
    if normalized in {"false", "0", "no", "n", "f"}:
        return False
    raise ContractError(f"Invalid boolean for {field_name}: {value!r}")


def optional_bool(value: Any, *, field_name: str) -> bool | None:
    if value is None or str(value).strip() in {"", "NA", "N/A", "None", "null"}:
        return None
    return parse_bool(value, field_name=field_name)


def optional_int(value: Any, *, field_name: str) -> int | None:
    if value is None or str(value).strip() in {"", "NA", "N/A", "None", "null"}:
        return None
    try:
        return int(value)
    except (TypeError, ValueError) as error:
        raise ContractError(f"Invalid integer for {field_name}: {value!r}") from error


def optional_float(value: Any, *, field_name: str) -> float | None:
    if value is None or str(value).strip() in {"", "NA", "N/A", "None", "null"}:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ContractError(f"Invalid float for {field_name}: {value!r}") from error
    if not math.isfinite(parsed):
        raise ContractError(f"Non-finite float for {field_name}: {value!r}")
    return parsed


def row_value(
    row: Mapping[str, Any],
    *field_names: str,
    default: Any = None,
) -> Any:
    for field_name in field_names:
        if field_name in row and str(row[field_name]).strip() not in {"", "NA", "N/A", "None", "null"}:
            return row[field_name]
    return default


def required_row_value(row: Mapping[str, Any], *field_names: str, source: str) -> Any:
    value = row_value(row, *field_names)
    if value is None:
        raise ContractError(f"{source} lacks required compatible field: {field_names}")
    return value


def optional_json(value: Any, *, field_name: str) -> Any:
    if value is None or str(value).strip() in {"", "NA", "N/A", "None", "null"}:
        return None
    if isinstance(value, (dict, list)):
        return value
    try:
        return json.loads(str(value))
    except json.JSONDecodeError as error:
        raise ContractError(f"Invalid JSON value for {field_name}: {value!r}") from error


def require_optional_float_equal(
    observed: Any,
    expected: float | None,
    *,
    field_name: str,
    source: str,
) -> None:
    parsed = optional_float(observed, field_name=field_name)
    if parsed is None or expected is None:
        if parsed is not None or expected is not None:
            raise ContractError(
                f"{source} {field_name} nullability drift: observed={parsed} expected={expected}"
            )
        return
    if not math.isclose(parsed, expected, rel_tol=1e-12, abs_tol=1e-12):
        raise ContractError(
            f"{source} {field_name} drift: observed={parsed} expected={expected}"
        )


def benjamini_hochberg(p_values: Iterable[float | None]) -> list[float | None]:
    values = list(p_values)
    valid = [
        (index, float(value))
        for index, value in enumerate(values)
        if value is not None and math.isfinite(float(value))
    ]
    valid.sort(key=lambda item: item[1])
    adjusted: list[float | None] = [None] * len(values)
    running = 1.0
    total = len(valid)
    for rank_index in range(total - 1, -1, -1):
        original_index, p_value = valid[rank_index]
        rank = rank_index + 1
        running = min(running, p_value * total / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def benjamini_yekutieli(p_values: Iterable[float | None]) -> list[float | None]:
    values = list(p_values)
    valid_count = sum(value is not None for value in values)
    if valid_count == 0:
        return [None] * len(values)
    harmonic = sum(1.0 / rank for rank in range(1, valid_count + 1))
    bh = benjamini_hochberg(values)
    return [None if value is None else min(1.0, value * harmonic) for value in bh]


def recompute_binomial_violation_gate(
    violations: int,
    denominator: int,
    *,
    error_ceiling: float = FOUR_STATE_ERROR_CEILING,
    confidence: float = FOUR_STATE_PER_RELATION_CONFIDENCE,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "violations": violations,
        "denominator": denominator,
        "rate": violations / denominator if denominator else None,
        "p_exact_greater": None,
        "upper_bound": None,
        "threshold": error_ceiling,
        "confidence": confidence,
        "status": "NOT_IDENTIFIABLE_NO_DENOMINATOR",
    }
    if denominator <= 0:
        return result
    p_value = float(
        binomtest(violations, denominator, error_ceiling, alternative="greater").pvalue
    )
    upper_bound = (
        1.0
        if violations == denominator
        else float(beta.ppf(confidence, violations + 1, denominator - violations))
    )
    if upper_bound <= error_ceiling + 1e-15:
        status = "COMPATIBLE_WITH_FIXED_ERROR_CEILING"
    elif p_value <= 1.0 - confidence + 1e-15:
        status = "VIOLATES_FIXED_ERROR_CEILING"
    else:
        status = "NOT_IDENTIFIABLE_LOW_PRECISION"
    return {
        **result,
        "p_exact_greater": p_value,
        "upper_bound": upper_bound,
        "status": status,
    }


def recompute_four_state_from_counts(state_counts: Mapping[str, Any]) -> dict[str, Any]:
    expected_keys = {"RR", "AR", "RA", "AA", "O", "X"}
    if set(state_counts) != expected_keys:
        raise ContractError(
            "Four-state count keys drift: "
            f"observed={sorted(state_counts)} expected={sorted(expected_keys)}"
        )
    parsed: dict[str, int] = {}
    for state in sorted(expected_keys):
        try:
            count = int(state_counts[state])
        except (TypeError, ValueError) as error:
            raise ContractError(f"Invalid four-state count for {state}: {state_counts[state]!r}") from error
        if count < 0 or str(state_counts[state]).strip() != str(count):
            raise ContractError(f"Invalid four-state count for {state}: {state_counts[state]!r}")
        parsed[state] = count
    rr, ar, ra, aa = (parsed[state] for state in ("RR", "AR", "RA", "AA"))
    called = rr + ar + ra + aa
    focal_ref = rr + ra
    focal_alt = ar + aa
    partner_ref = rr + ar
    partner_alt = ra + aa
    focal_gate = recompute_binomial_violation_gate(ra, partner_alt)
    partner_gate = recompute_binomial_violation_gate(ar, focal_alt)
    branching_gate = recompute_binomial_violation_gate(aa, ar + ra + aa)
    complete = called >= 10 and min(focal_ref, focal_alt, partner_ref, partner_alt) >= 3
    focal_split = ar >= FOUR_STATE_MIN_SPLIT_COUNT and aa >= FOUR_STATE_MIN_SPLIT_COUNT
    partner_split = ra >= FOUR_STATE_MIN_SPLIT_COUNT and aa >= FOUR_STATE_MIN_SPLIT_COUNT
    branching_split = ar >= FOUR_STATE_MIN_SPLIT_COUNT and ra >= FOUR_STATE_MIN_SPLIT_COUNT
    models: list[str] = []
    if complete and focal_split and focal_gate["status"] == "COMPATIBLE_WITH_FIXED_ERROR_CEILING":
        models.append("FOCAL_ANCESTOR")
    if complete and partner_split and partner_gate["status"] == "COMPATIBLE_WITH_FIXED_ERROR_CEILING":
        models.append("PARTNER_ANCESTOR")
    if complete and branching_split and branching_gate["status"] == "COMPATIBLE_WITH_FIXED_ERROR_CEILING":
        models.append("BRANCHING")
    relevant_statuses = []
    if focal_split:
        relevant_statuses.append(str(focal_gate["status"]))
    if partner_split:
        relevant_statuses.append(str(partner_gate["status"]))
    if branching_split:
        relevant_statuses.append(str(branching_gate["status"]))
    if not complete:
        relation = "NOT_IDENTIFIABLE_INSUFFICIENT_FOUR_STATE_DEPTH"
    elif len(models) > 1:
        relation = "MULTIPLE_MUTATION_ORDER_MODELS_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    elif models == ["FOCAL_ANCESTOR"]:
        relation = "FOCAL_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    elif models == ["PARTNER_ANCESTOR"]:
        relation = "PARTNER_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    elif models == ["BRANCHING"]:
        relation = "BRANCHING_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    elif not relevant_statuses or any(
        status.startswith("NOT_IDENTIFIABLE") for status in relevant_statuses
    ):
        relation = "NOT_IDENTIFIABLE_FIXED_ERROR_CEILING"
    else:
        relation = "INCOMPATIBLE_OR_COMPLEX_UNDER_FIXED_ERROR_MODEL"
    minimum_zero_depth = math.ceil(
        math.log(1.0 - FOUR_STATE_PER_RELATION_CONFIDENCE)
        / math.log(1.0 - FOUR_STATE_ERROR_CEILING)
    )
    status = (
        "NOT_IDENTIFIABLE_NO_FOCAL_REF"
        if focal_ref == 0
        else "NOT_IDENTIFIABLE_INSUFFICIENT_FOUR_STATE"
        if not complete
        else relation
    )
    return {
        "n_joint": called,
        "n_called_depth": called,
        "n_focal_ref": focal_ref,
        "n_focal_alt": focal_alt,
        "n_partner_ref": partner_ref,
        "n_partner_alt": partner_alt,
        "focal_ancestor_violation_rate": focal_gate["rate"],
        "partner_ancestor_violation_rate": partner_gate["rate"],
        "error_ceiling": FOUR_STATE_ERROR_CEILING,
        "confidence": FOUR_STATE_PER_RELATION_CONFIDENCE,
        "familywise_confidence": FOUR_STATE_FAMILYWISE_CONFIDENCE,
        "relation_family_size": FOUR_STATE_RELATION_FAMILY_SIZE,
        "multiplicity_method": FOUR_STATE_MULTIPLICITY_METHOD,
        "minimum_zero_violation_depth": minimum_zero_depth,
        "focal_ancestor_violation": focal_gate,
        "partner_ancestor_violation": partner_gate,
        "branching_violation": branching_gate,
        "complete_four_state_testable": complete,
        "relation_compatibility": relation,
        "compatible_relation_models": models,
        "n_compatible_relation_models": len(models),
        "status": status,
    }


def claim_status(passed: bool, *, evaluable: bool, ran: bool = True) -> str:
    if not ran:
        return "NOT_RUN"
    if not evaluable:
        return "NOT_EVALUABLE"
    return "PASS" if passed else "FAIL"


def validate_claim_status(value: str, *, field_name: str) -> str:
    if value not in CLAIM_STATUSES:
        raise ContractError(f"Invalid claim status for {field_name}: {value!r}")
    return value


def candidate_key(row: Mapping[str, Any], *, source: str) -> CandidateKey:
    try:
        sample = str(row["sample"]).strip()
        chrom = str(row["chrom"]).strip()
        pos = int(row["pos"])
        ref = str(row["ref"]).strip().upper()
        alt = str(row["alt"]).strip().upper()
    except (KeyError, TypeError, ValueError) as error:
        raise ContractError(f"Malformed candidate key in {source}: {row!r}") from error
    if sample not in DATASETS or not chrom or pos <= 0 or not ref or not alt:
        raise ContractError(f"Invalid candidate key in {source}: {(sample, chrom, pos, ref, alt)!r}")
    return sample, chrom, pos, ref, alt


def assignment_key(row: Mapping[str, Any], *, source: str) -> CandidateKey:
    posthoc = row.get("posthoc")
    if not isinstance(posthoc, Mapping):
        raise ContractError(f"Assignment lacks posthoc identity in {source}")
    merged = {
        "sample": row.get("sample"),
        "chrom": row.get("chrom"),
        "pos": row.get("pos"),
        "ref": posthoc.get("ref"),
        "alt": posthoc.get("alt"),
    }
    return candidate_key(merged, source=source)


def biological_key(key: CandidateKey, biological_id: str | None = None) -> BiologicalKey:
    sample, chrom, pos, ref, alt = key
    observed = biological_id or BIOLOGICAL_IDS[sample]
    return observed, chrom, pos, ref, alt


def compact_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def spatially_separated_positions(
    positions: Iterable[int], minimum_separation: int = 20
) -> list[int]:
    if minimum_separation < 1:
        raise ValueError("minimum_separation must be positive")
    selected: list[int] = []
    for position in sorted(set(int(value) for value in positions)):
        if not selected or position - selected[-1] >= minimum_separation:
            selected.append(position)
    return selected


def rate(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def require_columns(reader: csv.DictReader, required: Iterable[str], label: str) -> None:
    fields = list(reader.fieldnames or [])
    if not fields or len(fields) != len(set(fields)):
        raise ContractError(f"{label} has missing or duplicate header fields")
    missing = sorted(set(required).difference(fields))
    if missing:
        raise ContractError(f"{label} missing columns: {missing}")


def load_manifest(path: Path, expected_screen_sites: int) -> dict[str, Any]:
    manifest = load_json(path, "all-sSNV input manifest")
    require_schema_pass(
        manifest,
        "intersubmod.all_ssnv_focal_alt_input_manifest",
        "all-sSNV input manifest",
        schema_version=INPUT_MANIFEST_SCHEMA_VERSION,
    )
    samples = manifest.get("samples")
    if not isinstance(samples, list) or tuple(row.get("sample") for row in samples) != DATASETS:
        raise ContractError("Manifest must contain the canonical 7/7 dataset order")
    totals = manifest.get("totals")
    if not isinstance(totals, Mapping):
        raise ContractError("Manifest lacks totals")
    if int(totals.get("all_ssnv", -1)) != expected_screen_sites:
        raise ContractError(
            f"Manifest all-sSNV total mismatch: {totals.get('all_ssnv')} != {expected_screen_sites}"
        )
    truth_total = sum(int(totals.get(f"truth_{label.lower()}", -1)) for label in TRUTH_LABELS)
    if truth_total != expected_screen_sites:
        raise ContractError(
            f"Manifest TP+FP+UNASSESSED mismatch: {truth_total} != {expected_screen_sites}"
        )
    sample_total = 0
    sample_truth = Counter()
    for entry in samples:
        sample = str(entry["sample"])
        expected_biological = BIOLOGICAL_IDS[sample]
        if entry.get("biological_id") != expected_biological:
            raise ContractError(
                f"Manifest biological_id drift for {sample}: {entry.get('biological_id')!r}"
            )
        counts = entry.get("counts")
        if not isinstance(counts, Mapping):
            raise ContractError(f"Manifest sample lacks counts: {sample}")
        n_sample = int(counts.get("all_ssnv", -1))
        sample_truth_total = 0
        for label in TRUTH_LABELS:
            value = int(counts.get(f"truth_{label.lower()}", -1))
            if value < 0:
                raise ContractError(f"Negative manifest truth count for {sample}/{label}")
            sample_truth[label] += value
            sample_truth_total += value
        if sample_truth_total != n_sample:
            raise ContractError(
                f"Manifest sample TP+FP+UNASSESSED mismatch for {sample}: "
                f"{sample_truth_total} != {n_sample}"
            )
        sample_total += n_sample
    if sample_total != expected_screen_sites:
        raise ContractError(f"Manifest sample totals mismatch: {sample_total} != {expected_screen_sites}")
    for label in TRUTH_LABELS:
        if sample_truth[label] != int(totals[f"truth_{label.lower()}"]):
            raise ContractError(f"Manifest pooled/sample truth mismatch for {label}")
    return manifest


def _stage_keys(sample: str, truth_label: str) -> tuple[tuple[str, str], ...]:
    return (
        ("pooled", "ALL"),
        ("sample", sample),
        ("truth", truth_label),
        ("sample_truth", f"{sample}|{truth_label}"),
    )


def load_screen(path: Path, manifest: Mapping[str, Any], expected_screen_sites: int) -> ScreenData:
    all_keys: set[CandidateKey] = set()
    all_rows: dict[CandidateKey, dict[str, str]] = {}
    stable_rows: dict[CandidateKey, dict[str, str]] = {}
    stage_counts: dict[tuple[str, str], ScreenStageCounts] = defaultdict(ScreenStageCounts)
    per_sample_truth: dict[str, Counter[str]] = {sample: Counter() for sample in DATASETS}
    claim_dataset_sets = {name: set() for name in ("evaluable", "stable", "robust")}
    claim_biological_sets = {name: set() for name in ("evaluable", "stable", "robust")}
    latest_tag_audits: dict[tuple[str, str], LatestTagAudit] = defaultdict(LatestTagAudit)
    required = {
        "sample",
        "biological_id",
        "truth_label",
        "chrom",
        "pos",
        "ref",
        "alt",
        "analysis_status",
        "stable_null_multigroup",
        "modal_assignment_ari_min",
        "hp_axis_confound",
        "technical_axis_confound",
        "residual_unexplained_multigroup",
        "phase_anchored_robust_epigenetic_candidate",
        "ssnv_branch",
        "latest_tag_join_status",
        "latest_tag_rows_fetched",
        "latest_tag_rows_eligible",
        "latest_tag_reads_joined",
        "latest_tag_ps_present",
        "latest_tag_projection_multimatch_reads",
        "source_hp_replaced_reads",
        "n_reads_total",
        *M2_GATE.AXIS_SOURCE_FIELDS,
    }
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader, required, "screen site table")
        for row_number, row in enumerate(reader, 2):
            key = candidate_key(row, source=f"screen row {row_number}")
            if key in all_keys:
                raise ContractError(f"Duplicate screen candidate key: {key}")
            all_keys.add(key)
            all_rows[key] = dict(row)
            sample = key[0]
            biological_id = str(row["biological_id"])
            if biological_id != BIOLOGICAL_IDS[sample]:
                raise ContractError(
                    f"Screen biological_id drift at {key}: {biological_id!r}"
                )
            truth_label = str(row["truth_label"])
            if truth_label not in TRUTH_LABELS:
                raise ContractError(f"Invalid truth label at {key}: {truth_label!r}")
            per_sample_truth[sample][truth_label] += 1
            branch = str(row.get("ssnv_branch") or "NA")
            for audit_key in (
                ("pooled", "ALL"),
                ("sample", sample),
                ("truth", truth_label),
                ("biological", biological_id),
                ("branch", branch),
            ):
                latest_tag_audits[audit_key].add(
                    row, source=f"screen row {row_number} {key}"
                )
            evaluable = row["analysis_status"] == "evaluable"
            stable = parse_bool(
                row["stable_null_multigroup"], field_name="stable_null_multigroup"
            )
            modal_ari = optional_float(
                row["modal_assignment_ari_min"], field_name="modal_assignment_ari_min"
            )
            hp_confound = parse_bool(
                row["hp_axis_confound"], field_name="hp_axis_confound"
            )
            technical_confound = parse_bool(
                row["technical_axis_confound"], field_name="technical_axis_confound"
            )
            residual = parse_bool(
                row["residual_unexplained_multigroup"],
                field_name="residual_unexplained_multigroup",
            )
            robust = parse_bool(
                row["phase_anchored_robust_epigenetic_candidate"],
                field_name="phase_anchored_robust_epigenetic_candidate",
            )
            if stable and not evaluable:
                raise ContractError(f"Stable screen row is not evaluable: {key}")
            if stable and (modal_ari is None or modal_ari < 0.8):
                raise ContractError(f"M1 stable row fails modal assignment ARI>=0.8: {key}")
            if residual != (stable and not hp_confound and not technical_confound):
                raise ContractError(f"Residual-confound gate drift at {key}")
            if robust and not stable:
                raise ContractError(f"Robust epigenetic screen row is not stable: {key}")
            if robust and not residual:
                raise ContractError(f"M2 robust row is not residual-unexplained: {key}")
            bio_key = biological_key(key, biological_id)
            if evaluable:
                claim_dataset_sets["evaluable"].add(key)
                claim_biological_sets["evaluable"].add(bio_key)
            if stable:
                claim_dataset_sets["stable"].add(key)
                claim_biological_sets["stable"].add(bio_key)
                stable_rows[key] = dict(row)
            if robust:
                claim_dataset_sets["robust"].add(key)
                claim_biological_sets["robust"].add(bio_key)
            for stage_key in _stage_keys(sample, truth_label):
                counts = stage_counts[stage_key]
                counts.screened += 1
                counts.evaluable += int(evaluable)
                counts.stable += int(stable)
                counts.robust += int(robust)
    if len(all_keys) != expected_screen_sites:
        raise ContractError(
            f"Screen row count mismatch: {len(all_keys)} != {expected_screen_sites}"
        )
    pooled_tag_audit = latest_tag_audits[("pooled", "ALL")].payload()
    if pooled_tag_audit["n_sites"] != expected_screen_sites or pooled_tag_audit["pass"] is not True:
        raise ContractError(
            "Screen latest HP/PS terminal join audit does not prove every screen row was joined"
        )
    manifest_by_sample = {entry["sample"]: entry for entry in manifest["samples"]}
    for sample in DATASETS:
        expected_counts = manifest_by_sample[sample]["counts"]
        observed = stage_counts[("sample", sample)].screened
        if observed != int(expected_counts["all_ssnv"]):
            raise ContractError(
                f"Screen/manifest sample count mismatch for {sample}: "
                f"{observed} != {expected_counts['all_ssnv']}"
            )
        for label in TRUTH_LABELS:
            expected_truth = int(expected_counts[f"truth_{label.lower()}"])
            if per_sample_truth[sample][label] != expected_truth:
                raise ContractError(
                    f"Screen/manifest truth count mismatch for {sample}/{label}: "
                    f"{per_sample_truth[sample][label]} != {expected_truth}"
                )
    pooled_truth_total = sum(stage_counts[("truth", label)].screened for label in TRUTH_LABELS)
    if pooled_truth_total != expected_screen_sites:
        raise ContractError(
            f"Screen TP+FP+UNASSESSED mismatch: {pooled_truth_total} != {expected_screen_sites}"
        )
    for stage_key in (
        ("pooled", "ALL"),
        *(("sample", sample) for sample in DATASETS),
        *(("truth", label) for label in TRUTH_LABELS),
        *(
            ("sample_truth", f"{sample}|{label}")
            for sample in DATASETS
            for label in TRUTH_LABELS
        ),
    ):
        stage_counts[stage_key]
    return ScreenData(
        all_keys=all_keys,
        all_rows=all_rows,
        stable_rows=stable_rows,
        stage_counts=dict(stage_counts),
        per_sample_truth=per_sample_truth,
        claim_dataset_sets=claim_dataset_sets,
        claim_biological_sets=claim_biological_sets,
        latest_tag_audits=dict(latest_tag_audits),
    )


def load_assignment_keys(path: Path, screen: ScreenData) -> set[CandidateKey]:
    keys: set[CandidateKey] = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as error:
                raise ContractError(f"Invalid assignment JSON at line {line_number}") from error
            if not isinstance(row, Mapping):
                raise ContractError(f"Assignment line {line_number} is not an object")
            if (
                row.get("schema_name") != ASSIGNMENT_SCHEMA
                or row.get("schema_version") != ASSIGNMENT_SCHEMA_VERSION
            ):
                raise ContractError(f"Unexpected assignment schema at line {line_number}")
            if row.get("strict_confirm_candidate") is not True:
                raise ContractError(f"Assignment line {line_number} is not a stable candidate")
            key = assignment_key(row, source=f"assignment line {line_number}")
            if key in keys:
                raise ContractError(f"Duplicate stable assignment key: {key}")
            keys.add(key)
            if key not in screen.stable_rows:
                raise ContractError(f"Stable assignment absent from screen stable set: {key}")
            posthoc = row["posthoc"]
            screen_row = screen.stable_rows[key]
            for field_name in ("biological_id", "truth_label"):
                value = posthoc.get(field_name)
                if value not in {None, ""} and str(value) != str(screen_row[field_name]):
                    raise ContractError(
                        f"Assignment/screen {field_name} mismatch at {key}: "
                        f"{value!r} != {screen_row[field_name]!r}"
                    )
    expected = set(screen.stable_rows)
    if keys != expected:
        raise ContractError(
            "Screen stable/assignment key-set mismatch: "
            f"screen={len(expected)} assignments={len(keys)} "
            f"missing={sorted(expected - keys)[:3]} extra={sorted(keys - expected)[:3]}"
        )
    return keys


def independently_verify_primary_artifact_task(
    task: tuple[CandidateKey, str, dict[str, Any]],
) -> list[str]:
    key, region_raw, contracts = task
    region_dir = Path(region_raw).expanduser().resolve()
    key_text = "|".join(str(value) for value in key)
    lines: list[str] = []
    for role, relative_path in INDEPENDENT_PRIMARY_ARTIFACT_PATHS.items():
        record = contracts.get(role)
        if not isinstance(record, Mapping):
            raise ContractError(f"Independent recount missing {role} at {key}")
        if not {"path", "size_bytes", "sha256"}.issubset(record):
            raise ContractError(f"Independent recount malformed {role} identity at {key}")
        expected_path = (region_dir / relative_path).resolve()
        declared_path = Path(str(record["path"])).expanduser().resolve()
        if declared_path != expected_path:
            raise ContractError(f"Independent recount {role} path mismatch at {key}")
        if not expected_path.is_file():
            raise ContractError(f"Independent recount missing {role} file at {key}")
        before = expected_path.stat()
        observed_sha = sha256(expected_path)
        after = expected_path.stat()
        if before.st_size != after.st_size or before.st_mtime_ns != after.st_mtime_ns:
            raise ContractError(f"Independent recount {role} changed while hashing at {key}")
        if int(record["size_bytes"]) != after.st_size:
            raise ContractError(f"Independent recount {role} size mismatch at {key}")
        if str(record["sha256"]) != observed_sha:
            raise ContractError(f"Independent recount {role} SHA-256 mismatch at {key}")
        lines.append(
            "|".join(
                (key_text, role, str(expected_path), str(after.st_size), observed_sha)
            )
        )
    return lines


def independently_recount_primary_artifacts(
    assignments_path: Path,
    expected_keys: set[CandidateKey],
    expected_digest: str,
) -> dict[str, Any]:
    """Rehash all primary artifacts without using the primary-auditor implementation."""
    assignments_before = artifact(assignments_path)
    tasks: dict[CandidateKey, tuple[CandidateKey, str, dict[str, Any]]] = {}
    with open_text(assignments_path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                assignment = json.loads(line)
            except json.JSONDecodeError as error:
                raise ContractError(
                    f"Independent recount malformed assignment JSON at line {line_number}"
                ) from error
            if not isinstance(assignment, dict):
                raise ContractError(
                    f"Independent recount assignment line {line_number} is not an object"
                )
            if (
                assignment.get("schema_name") != ASSIGNMENT_SCHEMA
                or assignment.get("schema_version") != ASSIGNMENT_SCHEMA_VERSION
                or assignment.get("screen_contract")
                != "phylo-v4.1_column_null95_modal_K10_RNULL40_min_group3"
                or assignment.get("artifact_identity_contract") != "sha256_size_path_v1"
                or assignment.get("strict_confirm_candidate") is not True
            ):
                raise ContractError(
                    f"Independent recount assignment contract drift at line {line_number}"
                )
            posthoc = assignment.get("posthoc")
            if not isinstance(posthoc, dict):
                raise ContractError(
                    f"Independent recount assignment posthoc missing at line {line_number}"
                )
            key = candidate_key(
                {
                    "sample": assignment.get("sample") or assignment.get("dataset"),
                    "chrom": assignment.get("chrom"),
                    "pos": assignment.get("pos"),
                    "ref": posthoc.get("ref"),
                    "alt": posthoc.get("alt"),
                },
                source=f"independent primary recount line {line_number}",
            )
            if key in tasks:
                raise ContractError(f"Independent recount duplicate assignment key: {key}")
            region_raw = assignment.get("region_dir")
            contracts = assignment.get("primary_artifacts")
            if not isinstance(region_raw, str) or not region_raw:
                raise ContractError(f"Independent recount region_dir missing at {key}")
            if not isinstance(contracts, dict) or set(contracts) != set(
                INDEPENDENT_PRIMARY_ARTIFACT_PATHS
            ):
                raise ContractError(f"Independent recount artifact role drift at {key}")
            tasks[key] = (key, region_raw, contracts)
    if set(tasks) != expected_keys:
        raise ContractError(
            "Independent primary recount assignment key-set mismatch: "
            f"expected={len(expected_keys)} observed={len(tasks)}"
        )
    all_lines: list[str] = []
    ordered_tasks = [tasks[key] for key in sorted(tasks)]
    with ProcessPoolExecutor(max_workers=INDEPENDENT_PRIMARY_RECOUNT_WORKERS) as executor:
        for lines in executor.map(
            independently_verify_primary_artifact_task,
            ordered_tasks,
            chunksize=INDEPENDENT_PRIMARY_RECOUNT_CHUNK_SIZE,
        ):
            all_lines.extend(lines)
    observed_digest = hashlib.sha256(
        "\n".join(sorted(all_lines)).encode("utf-8")
    ).hexdigest()
    assignments_after = artifact(assignments_path)
    if assignments_after != assignments_before:
        raise ContractError("Stable assignments changed during independent primary recount")
    if observed_digest != expected_digest:
        raise ContractError(
            "Independent primary artifact-set digest does not match pre/post audit: "
            f"{observed_digest} != {expected_digest}"
        )
    return {
        "contract": "independent_all_primary_artifact_path_size_sha256_recount_v1",
        "implementation_independence": True,
        "primary_auditor_functions_called": False,
        "assignments": assignments_after,
        "stable_sites": len(tasks),
        "primary_artifacts_verified": len(all_lines),
        "artifact_roles": sorted(INDEPENDENT_PRIMARY_ARTIFACT_PATHS),
        "artifact_set_sha256": observed_digest,
        "workers": INDEPENDENT_PRIMARY_RECOUNT_WORKERS,
        "chunk_size": INDEPENDENT_PRIMARY_RECOUNT_CHUNK_SIZE,
        "pass": True,
    }


def load_keyed_tsv(
    path: Path,
    *,
    label: str,
    required_fields: Iterable[str],
) -> dict[CandidateKey, dict[str, str]]:
    rows: dict[CandidateKey, dict[str, str]] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(
            reader,
            {"sample", "chrom", "pos", "ref", "alt", *required_fields},
            label,
        )
        for row_number, row in enumerate(reader, 2):
            key = candidate_key(row, source=f"{label} row {row_number}")
            if key in rows:
                raise ContractError(f"Duplicate {label} candidate key: {key}")
            rows[key] = dict(row)
    return rows


def validate_recovery_screen_contract(
    summary: Mapping[str, Any],
    receipt: Mapping[str, Any],
) -> dict[str, Any] | None:
    execution = receipt.get("execution")
    if not isinstance(execution, Mapping):
        raise ContractError("Screen receipt lacks execution contract")
    if execution.get("recovery_merge") is not True:
        if "recovery_merge" in summary or "recovery" in receipt:
            raise ContractError("Screen recovery metadata exists without recovery execution flag")
        return None
    summary_recovery = summary.get("recovery_merge")
    recovery = receipt.get("recovery")
    if not isinstance(summary_recovery, Mapping) or not isinstance(recovery, Mapping):
        raise ContractError("Recovery screen lacks summary/run-manifest recovery metadata")
    required_summary = {
        "replacement_sample": "HCC1954",
        "method_change": False,
        "prefix_source_locked": True,
        "replacement_read_only_pinned_source": True,
        "serial_parallel_exact_equivalence_pass": True,
    }
    for field_name, expected in required_summary.items():
        if summary_recovery.get(field_name) != expected:
            raise ContractError(f"Recovery summary contract drift: {field_name}")
    required_recovery = {
        "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_recovery_merge",
        "replacement_sample": "HCC1954",
        "prefix_samples": list(DATASETS[:-1]),
        "method_change": False,
        "serial_parallel_exact_equivalence_required": True,
        "prefix_replacement_source_dependencies_exact": True,
        "recovery_source_identity_unchanged_during_merge": True,
    }
    for field_name, expected in required_recovery.items():
        if recovery.get(field_name) != expected:
            raise ContractError(f"Recovery run-manifest contract drift: {field_name}")
    prefix_artifacts = recovery.get("prefix_source_artifacts")
    if not isinstance(prefix_artifacts, Mapping):
        raise ContractError("Recovery receipt lacks prefix source artifacts")
    prefix_paths = {
        role: verify_declared_path_artifact(reference, f"recovery prefix {role}")
        for role, reference in prefix_artifacts.items()
    }
    required_prefix_roles = {
        "sites",
        "assignments",
        "summary",
        "receipt",
        "source_lock_receipt",
    }
    if set(prefix_paths) != required_prefix_roles:
        raise ContractError("Recovery prefix artifact roles drift")
    replacement_summary_path = verify_declared_path_artifact(
        recovery.get("replacement_summary"), "recovery replacement summary"
    )
    replacement_receipt_path = verify_declared_path_artifact(
        recovery.get("replacement_receipt"), "recovery replacement receipt"
    )
    equivalence_path = verify_declared_path_artifact(
        recovery.get("serial_parallel_exact_equivalence_receipt"),
        "recovery serial/parallel equivalence receipt",
    )
    prefix_receipt = load_json(prefix_paths["receipt"], "recovery prefix receipt")
    prefix_source_lock = load_json(
        prefix_paths["source_lock_receipt"], "recovery prefix source-lock receipt"
    )
    replacement_receipt = load_json(
        replacement_receipt_path, "recovery replacement receipt"
    )
    replacement_summary = load_json(
        replacement_summary_path, "recovery replacement summary"
    )
    equivalence = load_json(equivalence_path, "recovery equivalence receipt")
    if (
        prefix_source_lock.get("schema_name") != "intersubmod.source_locked_screen_run"
        or prefix_source_lock.get("pass") is not True
        or prefix_source_lock.get("source_identity_before")
        != prefix_source_lock.get("source_identity_after")
        or prefix_source_lock.get("child_run_manifest") != artifact(prefix_paths["receipt"])
    ):
        raise ContractError("Recovery prefix source-lock receipt failed independent validation")
    source_lock_checks = prefix_source_lock.get("checks") or {}
    if not source_lock_checks or not all(value is True for value in source_lock_checks.values()):
        raise ContractError("Recovery prefix source-lock checks are incomplete")
    if replacement_receipt.get("pass") is not True or replacement_summary.get("pass") is not True:
        raise ContractError("Recovery replacement output is not passing")
    if equivalence.get("schema_name") != "intersubmod.phylo_parallel_exact_equivalence" or equivalence.get("pass") is not True:
        raise ContractError("Recovery equivalence receipt is not passing")
    equivalence_checks = equivalence.get("checks") or {}
    expected_equivalence_checks = {
        "pinned_analyzer_sha256_exact",
        "source_identity_unchanged",
        "synthetic_full_payload_exact",
        "real_nested_full_payload_exact",
        "real_nested_parallel_triggered",
    }
    if set(equivalence_checks) != expected_equivalence_checks or not all(
        equivalence_checks.values()
    ):
        raise ContractError("Recovery equivalence checks are incomplete")
    real_fixture = equivalence.get("real_fixture") or {}
    if (
        real_fixture.get("full_row_and_assignment_payload_exact") is not True
        or int(real_fixture.get("n_alt_after_peel", 0))
        < int(real_fixture.get("parallel_min_reads", 1))
    ):
        raise ContractError("Recovery real nested equivalence fixture is invalid")
    pinned_sha = str(recovery.get("pinned_analyzer_sha256", ""))
    if len(pinned_sha) != 64 or summary_recovery.get("pinned_analyzer_sha256") != pinned_sha:
        raise ContractError("Recovery pinned analyzer SHA-256 is invalid")
    prefix_sources = prefix_receipt.get("source_code") or {}
    replacement_sources = replacement_receipt.get("source_code") or {}
    equivalence_analyzer = (
        ((equivalence.get("inputs") or {}).get("source_identity_before") or {}).get(
            "analyzer"
        )
        or {}
    )
    for label, reference in (
        ("prefix analyzer", prefix_sources.get("analyzer")),
        ("replacement analyzer", replacement_sources.get("analyzer")),
        ("equivalence analyzer", equivalence_analyzer),
    ):
        path = verify_declared_path_artifact(reference, f"recovery {label}")
        if sha256(path) != pinned_sha:
            raise ContractError(f"Recovery {label} does not match pinned SHA-256")
    for role in ("focal_alt_cluster_lib", "latest_tag_join", "claim_contract_v2"):
        if prefix_sources.get(role) != replacement_sources.get(role):
            raise ContractError(f"Recovery prefix/replacement dependency drift: {role}")
        verify_declared_path_artifact(prefix_sources.get(role), f"recovery dependency {role}")
    merger_sources = receipt.get("source_code")
    if not isinstance(merger_sources, Mapping):
        raise ContractError("Recovery receipt lacks merger source identities")
    for role, reference in merger_sources.items():
        verify_declared_path_artifact(reference, f"recovery merger source {role}")
    return {
        "mode": "source_locked_prefix_plus_seed_parallel_replacement",
        "pinned_analyzer_sha256": pinned_sha,
        "prefix_source_lock_pass": True,
        "serial_parallel_exact_equivalence_pass": True,
        "real_nested_equivalence_fixture": dict(real_fixture),
    }


def validate_screen_receipts(
    bundle: InputBundle,
    manifest: Mapping[str, Any],
    screen: ScreenData,
    expected_screen_sites: int,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any] | None]:
    summary = load_json(bundle.screen_summary, "screen summary")
    receipt = load_json(bundle.screen_receipt, "screen run manifest")
    require_schema_pass(
        summary,
        "intersubmod.all_ssnv_focal_alt_multigroup_screen",
        "screen summary",
        schema_version=SCREEN_OUTPUT_SCHEMA_VERSION,
    )
    require_schema_pass(
        receipt,
        "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest",
        "screen run manifest",
        schema_version=SCREEN_OUTPUT_SCHEMA_VERSION,
    )
    for payload, label in ((summary, "summary"), (receipt, "run manifest")):
        if payload.get("status") != "EXECUTION_PASS":
            raise ContractError(f"Screen {label} status is not EXECUTION_PASS")
        if payload.get("pass_semantics") != "execution_integrity_only_not_scientific_confirmation":
            raise ContractError(f"Screen {label} pass semantics drift")
    scope = summary.get("scope")
    if not isinstance(scope, Mapping):
        raise ContractError("Screen summary lacks scope")
    if scope.get("full_469849") is not True:
        raise ContractError("Screen summary is not marked full_469849")
    if tuple(scope.get("selected_samples") or []) != DATASETS or tuple(
        scope.get("selected_datasets") or []
    ) != DATASETS:
        raise ContractError("Screen summary does not cover canonical 7/7 datasets")
    for name in ("expected_sites", "processed_sites"):
        if int(scope.get(name, -1)) != expected_screen_sites:
            raise ContractError(f"Screen summary {name} mismatch")
    pooled = summary.get("pooled_site_weighted")
    if not isinstance(pooled, Mapping):
        raise ContractError("Screen summary lacks pooled_site_weighted")
    if int(pooled.get("n_sites", -1)) != expected_screen_sites:
        raise ContractError("Screen summary pooled site count mismatch")
    stable_count = len(screen.stable_rows)
    if int(pooled.get("n_stable_null_multigroup", -1)) != stable_count:
        raise ContractError("Screen summary stable count mismatch")
    if int(summary.get("n_stable_assignment_records", -1)) != stable_count:
        raise ContractError("Screen summary assignment count mismatch")
    pooled_tag_audit = screen.latest_tag_audits[("pooled", "ALL")].payload()
    if pooled_tag_audit["n_sites"] != expected_screen_sites or pooled_tag_audit["pass"] is not True:
        raise ContractError("Recomputed pooled screen latest HP/PS audit is not passing")

    def require_tag_audit(container: Any, expected: Mapping[str, Any], label: str) -> None:
        if not isinstance(container, Mapping):
            raise ContractError(f"{label} is missing")
        observed = container.get("latest_hp_ps_terminal_join_audit")
        if observed != expected:
            raise ContractError(f"{label} latest HP/PS terminal join aggregate mismatch")

    require_tag_audit(summary, pooled_tag_audit, "Screen summary")
    require_tag_audit(pooled, pooled_tag_audit, "Screen pooled summary")
    summary_strata = (
        ("per_dataset", "sample", DATASETS),
        ("per_sample", "sample", DATASETS),
        ("posthoc_truth_strata", "truth", TRUTH_LABELS),
        (
            "posthoc_biological_id_strata",
            "biological",
            tuple(sorted(set(BIOLOGICAL_IDS.values()))),
        ),
        (
            "posthoc_ledger_branch_strata",
            "branch",
            tuple(
                sorted(
                    key_value
                    for key_type, key_value in screen.latest_tag_audits
                    if key_type == "branch"
                )
            ),
        ),
    )
    for container_name, audit_type, keys in summary_strata:
        container = summary.get(container_name)
        if not isinstance(container, Mapping) or set(container) != set(keys):
            raise ContractError(f"Screen summary {container_name} strata drift")
        for key_value in keys:
            require_tag_audit(
                container[key_value],
                screen.latest_tag_audits[(audit_type, key_value)].payload(),
                f"Screen summary {container_name}[{key_value}]",
            )
    receipt_counts = receipt.get("counts") or {}
    expected_receipt_counts = {
        "expected_sites": expected_screen_sites,
        "processed_sites": expected_screen_sites,
        "stable_assignment_records": stable_count,
        "reads_tsv_site_rows": pooled_tag_audit["n_reads_tsv_site_rows"],
        "exact_hp_ps_site_read_joins": pooled_tag_audit[
            "n_exact_hp_ps_site_read_joins"
        ],
        "ps_present_site_read_joins": pooled_tag_audit[
            "n_ps_present_site_read_joins"
        ],
        "source_hp_replaced_site_read_joins": pooled_tag_audit[
            "n_source_hp_replaced_site_read_joins"
        ],
    }
    for field_name, expected in expected_receipt_counts.items():
        if int(receipt_counts.get(field_name, -1)) != expected:
            raise ContractError(f"Screen receipt count mismatch for {field_name}")
    execution = receipt.get("execution")
    if not isinstance(execution, Mapping) or tuple(
        execution.get("selected_samples") or []
    ) != DATASETS or tuple(execution.get("selected_datasets") or []) != DATASETS:
        raise ContractError("Screen receipt execution scope does not cover canonical 7/7 datasets")
    recovery_source_validation = validate_recovery_screen_contract(summary, receipt)
    contracts = receipt.get("contracts")
    if not isinstance(contracts, Mapping):
        raise ContractError("Screen receipt lacks contracts")
    if contracts.get("latest_hp_ps_terminal_join") != pooled_tag_audit:
        raise ContractError("Screen receipt latest HP/PS terminal join contract mismatch")
    required_contracts = {
        "truth_and_cooccurrence_enter_clustering": False,
        "screen_global_fdr_calibrated": False,
        "strict_methyl_partition_robustness_status": "NOT_EVALUATED_AT_M1_SCREEN",
        "strict_confirm_status_legacy_alias": "NOT_RUN",
        "strict_confirm_candidate_is_formal_r1_claim": False,
        "existing_results_overwritten": False,
    }
    for field_name, expected in required_contracts.items():
        if contracts.get(field_name) != expected:
            raise ContractError(f"Screen receipt contract mismatch for {field_name}")
    clustering = summary.get("clustering_contract")
    if not isinstance(clustering, Mapping):
        raise ContractError("Screen summary lacks clustering contract")
    for field_name in (
        "m1_stability_gate_contract",
        "prior_screen_thresholds",
        "stable_null_multigroup_basis",
    ):
        if contracts.get(field_name) != clustering.get(field_name):
            raise ContractError(f"Screen summary/run-manifest contract drift for {field_name}")
    outputs = receipt.get("outputs") or {}
    verify_declared_artifact(outputs.get("site_results"), bundle.screen_sites, "screen sites")
    verify_declared_artifact(
        outputs.get("stable_assignments"), bundle.screen_assignments, "screen assignments"
    )
    verify_declared_artifact(outputs.get("summary"), bundle.screen_summary, "screen summary")
    verify_declared_artifact(
        receipt.get("input_manifest"), bundle.manifest, "screen input manifest"
    )
    for label in TRUTH_LABELS:
        expected_truth = int(manifest["totals"][f"truth_{label.lower()}"])
        observed_truth = screen.stage_counts[("truth", label)].screened
        if observed_truth != expected_truth:
            raise ContractError(f"Screen truth-stratum denominator mismatch for {label}")
    return summary, receipt, recovery_source_validation


def load_cooccurrence_sites(
    path: Path, screen: ScreenData
) -> tuple[
    dict[CandidateKey, dict[str, str]], set[CandidateKey], set[CandidateKey]
]:
    required = {
        "truth_label",
        "biological_id",
        "stable_null_multigroup",
        "modal_assignment_ari_min",
        "hp_axis_confound",
        "technical_axis_confound",
        "residual_unexplained_multigroup",
        "phase_anchored_robust_epigenetic_candidate",
        "m2_screen_gate_contract",
        "m2_screen_evaluable",
        "m2_screen_eligible",
        "m2_screen_eligibility_status",
        "m2_categorical_level_counts",
        "m2_axis_statuses",
        "m2_indeterminate_axes",
        "m2_low_power_axes",
        "m2_aligned_axes",
        "m2_constant_axes",
        "m2_aligned_below_negative_evaluability_power_axes",
        "n_partner_markers",
        "n_pair_rows_reconciled",
        "pair_row_count_reconciliation_pass",
        "n_endpoint_a_testable_markers",
        "n_endpoint_a_exact_identifiable_markers",
        "n_pair_by_confirmed",
        "pair_by_confirmed_positions",
        "n_spatially_separated_pair_by_20bp",
        "spatially_separated_pair_by_positions_20bp",
        "top_marker_positions",
        "n_top_marker_pair_by_confirmed",
        "top_marker_pair_by_confirmed_positions",
        "joint_signature_complete_marker_support",
        "joint_signature_n_complete_reads",
        "joint_signature_n_complete_marker_effect_supported",
        "joint_signature_complete_marker_effect_supported_positions",
        "joint_signature_testable",
        "joint_signature_p_conditional_perm",
        "joint_signature_permutations",
        "joint_signature_permutable",
        "joint_signature_conditional_status",
        "joint_signature_sensitivity_pass",
        "joint_signature_global_fdr_family_status",
        "joint_signature_q_global_bh",
        "joint_signature_q_global_by",
        "joint_signature_global_bh_discovery",
        "joint_signature_global_by_discovery",
        "n_same_complete_read_effect_supported_top_pair_by",
        "same_complete_read_effect_supported_top_pair_by_positions",
        "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp",
        "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp",
        FORMAL_SELECTION_COLUMN,
    }
    rows = load_keyed_tsv(
        path, label="cooccurrence site table", required_fields=required
    )
    if set(rows) != set(screen.stable_rows):
        expected = set(screen.stable_rows)
        observed = set(rows)
        raise ContractError(
            "Cooccurrence/screen stable key-set mismatch: "
            f"stable={len(expected)} cooccurrence={len(observed)} "
            f"missing={sorted(expected - observed)[:3]} extra={sorted(observed - expected)[:3]}"
        )
    pre_m2_formal_pair_sites: set[CandidateKey] = set()
    g2_sites: set[CandidateKey] = set()
    for key, row in rows.items():
        screen_row = screen.stable_rows[key]
        for field_name in ("truth_label", "biological_id"):
            if str(row[field_name]) != str(screen_row[field_name]):
                raise ContractError(f"Cooccurrence/screen {field_name} mismatch at {key}")
        for field_name in (
            "stable_null_multigroup",
            "hp_axis_confound",
            "technical_axis_confound",
            "residual_unexplained_multigroup",
            "phase_anchored_robust_epigenetic_candidate",
        ):
            if parse_bool(row[field_name], field_name=field_name) != parse_bool(
                screen_row[field_name], field_name=field_name
            ):
                raise ContractError(f"Cooccurrence/screen {field_name} mismatch at {key}")
        if optional_float(
            row["modal_assignment_ari_min"], field_name="modal_assignment_ari_min"
        ) != optional_float(
            screen_row["modal_assignment_ari_min"], field_name="modal_assignment_ari_min"
        ):
            raise ContractError(f"Cooccurrence/screen modal ARI mismatch at {key}")
        try:
            categorical_level_counts = optional_json(
                row["m2_categorical_level_counts"],
                field_name="m2_categorical_level_counts",
            )
            expected_gate = M2_GATE.evaluate_m2_screen(
                screen_row,
                categorical_level_counts=categorical_level_counts,
            )
        except M2_GATE.M2GateError as error:
            raise ContractError(f"M2 source evidence is internally inconsistent at {key}") from error
        m2 = parse_bool(row["m2_screen_eligible"], field_name="m2_screen_eligible")
        m2_evaluable = parse_bool(
            row["m2_screen_evaluable"], field_name="m2_screen_evaluable"
        )
        observed_axis_statuses = optional_json(
            row["m2_axis_statuses"], field_name="m2_axis_statuses"
        )
        observed_indeterminate_axes = optional_json(
            row["m2_indeterminate_axes"], field_name="m2_indeterminate_axes"
        )
        observed_low_power_axes = optional_json(
            row["m2_low_power_axes"], field_name="m2_low_power_axes"
        )
        observed_aligned_axes = optional_json(
            row["m2_aligned_axes"], field_name="m2_aligned_axes"
        )
        observed_constant_axes = optional_json(
            row["m2_constant_axes"], field_name="m2_constant_axes"
        )
        observed_aligned_below_power_axes = optional_json(
            row["m2_aligned_below_negative_evaluability_power_axes"],
            field_name="m2_aligned_below_negative_evaluability_power_axes",
        )
        gate_checks = {
            "contract": row["m2_screen_gate_contract"] == expected_gate["contract"],
            "eligible": m2 == expected_gate["eligible"],
            "evaluable": m2_evaluable == expected_gate["evaluable"],
            "status": row["m2_screen_eligibility_status"] == expected_gate["status"],
            "axis_statuses": observed_axis_statuses == expected_gate["axis_statuses"],
            "indeterminate_axes": observed_indeterminate_axes
            == expected_gate["indeterminate_axes"],
            "low_power_axes": observed_low_power_axes == expected_gate["low_power_axes"],
            "aligned_axes": observed_aligned_axes == expected_gate["aligned_axes"],
            "constant_axes": observed_constant_axes == expected_gate["constant_axes"],
            "categorical_level_counts": (
                categorical_level_counts == expected_gate["categorical_level_counts"]
            ),
            "aligned_below_negative_evaluability_power_axes": (
                observed_aligned_below_power_axes
                == expected_gate["aligned_below_negative_evaluability_power_axes"]
            ),
        }
        if not all(gate_checks.values()):
            failed = sorted(name for name, passed in gate_checks.items() if not passed)
            raise ContractError(f"M2 screen gate drift at {key}: {failed}")
        if not parse_bool(
            row["pair_row_count_reconciliation_pass"],
            field_name="pair_row_count_reconciliation_pass",
        ):
            raise ContractError(f"Pair-row reconciliation is not passing at {key}")
        n_partner = optional_int(row["n_partner_markers"], field_name="n_partner_markers")
        n_reconciled = optional_int(
            row["n_pair_rows_reconciled"], field_name="n_pair_rows_reconciled"
        )
        n_confirmed = optional_int(
            row["n_pair_by_confirmed"], field_name="n_pair_by_confirmed"
        )
        n_spaced = optional_int(
            row["n_spatially_separated_pair_by_20bp"],
            field_name="n_spatially_separated_pair_by_20bp",
        )
        n_top_confirmed = optional_int(
            row["n_top_marker_pair_by_confirmed"],
            field_name="n_top_marker_pair_by_confirmed",
        )
        if None in {n_partner, n_reconciled, n_confirmed, n_spaced, n_top_confirmed}:
            raise ContractError(f"Cooccurrence formal marker count is missing at {key}")
        if n_partner != n_reconciled or n_spaced > n_confirmed or n_top_confirmed > n_confirmed:
            raise ContractError(f"Cooccurrence site pair-count drift at {key}")
        confirmed_positions = optional_json(
            row["pair_by_confirmed_positions"], field_name="pair_by_confirmed_positions"
        )
        spaced_positions = optional_json(
            row["spatially_separated_pair_by_positions_20bp"],
            field_name="spatially_separated_pair_by_positions_20bp",
        )
        top_confirmed_positions = optional_json(
            row["top_marker_pair_by_confirmed_positions"],
            field_name="top_marker_pair_by_confirmed_positions",
        )
        if not all(
            isinstance(value, list)
            for value in (confirmed_positions, spaced_positions, top_confirmed_positions)
        ):
            raise ContractError(f"Cooccurrence formal marker positions are not lists at {key}")
        if (
            len(confirmed_positions) != n_confirmed
            or len(spaced_positions) != n_spaced
            or len(top_confirmed_positions) != n_top_confirmed
        ):
            raise ContractError(f"Cooccurrence formal marker position-count drift at {key}")
        if spaced_positions != spatially_separated_positions(confirmed_positions):
            raise ContractError(f"Cooccurrence spatially separated marker drift at {key}")
        top_marker_positions = optional_json(
            row["top_marker_positions"], field_name="top_marker_positions"
        )
        if not isinstance(top_marker_positions, list):
            raise ContractError(f"Cooccurrence top marker positions are not a list at {key}")
        top_marker_positions = [int(value) for value in top_marker_positions]
        if len(top_marker_positions) != len(set(top_marker_positions)):
            raise ContractError(f"Cooccurrence top marker positions contain duplicates at {key}")
        expected_top_confirmed = sorted(
            set(int(value) for value in top_marker_positions).intersection(
                int(value) for value in confirmed_positions
            )
        )
        if top_confirmed_positions != expected_top_confirmed:
            raise ContractError(f"Cooccurrence top/confirmed marker intersection drift at {key}")
        top_spaced_positions = spatially_separated_positions(
            top_confirmed_positions, minimum_separation=20
        )
        optional_top_spaced_fields = (
            (
                "n_top_spatially_separated_pair_by_20bp",
                len(top_spaced_positions),
            ),
            (
                "top_spatially_separated_pair_by_positions_20bp",
                top_spaced_positions,
            ),
        )
        for field_name, expected_value in optional_top_spaced_fields:
            if field_name not in row:
                continue
            observed_value = (
                optional_int(row[field_name], field_name=field_name)
                if field_name.startswith("n_")
                else optional_json(row[field_name], field_name=field_name)
            )
            if observed_value != expected_value:
                raise ContractError(f"Cooccurrence top-spaced marker drift at {key}: {field_name}")
        complete_marker_support = optional_json(
            row["joint_signature_complete_marker_support"],
            field_name="joint_signature_complete_marker_support",
        )
        if not isinstance(complete_marker_support, list):
            raise ContractError(
                f"Cooccurrence complete-read marker support is not a list at {key}"
            )
        support_positions: set[int] = set()
        recomputed_effect_positions: set[int] = set()
        for support_index, support in enumerate(complete_marker_support):
            if not isinstance(support, Mapping):
                raise ContractError(
                    f"Cooccurrence complete-read marker support record is invalid at "
                    f"{key}/{support_index}"
                )
            position = optional_int(
                support.get("position"), field_name="complete marker support position"
            )
            if position is None or position in support_positions:
                raise ContractError(
                    f"Cooccurrence complete-read marker support position is missing/duplicate "
                    f"at {key}/{support_index}"
                )
            support_positions.add(position)
            testable = parse_bool(
                support.get("testable"), field_name="complete marker support testable"
            )
            cramers_v = optional_float(
                support.get("cramers_v"), field_name="complete marker support cramers_v"
            )
            delta_alt = optional_float(
                support.get("delta_alt_fraction"),
                field_name="complete marker support delta_alt_fraction",
            )
            expected_effect_pass = bool(
                testable
                and cramers_v is not None
                and cramers_v >= COMPLETE_MARKER_CRAMERS_V_MIN
                and delta_alt is not None
                and delta_alt >= COMPLETE_MARKER_DELTA_ALT_MIN
            )
            declared_effect_pass = parse_bool(
                support.get("effect_gate_pass"),
                field_name="complete marker support effect_gate_pass",
            )
            if declared_effect_pass != expected_effect_pass:
                raise ContractError(
                    f"Cooccurrence complete-read marker effect gate drift at {key}/{position}"
                )
            if expected_effect_pass:
                recomputed_effect_positions.add(position)
        if support_positions != set(top_marker_positions):
            raise ContractError(
                f"Cooccurrence complete-read support/top-marker set drift at {key}"
            )
        complete_effect_positions_raw = optional_json(
            row["joint_signature_complete_marker_effect_supported_positions"],
            field_name="joint_signature_complete_marker_effect_supported_positions",
        )
        if not isinstance(complete_effect_positions_raw, list):
            raise ContractError(
                f"Cooccurrence complete-read effect-supported positions are not a list at {key}"
            )
        complete_effect_positions = [int(value) for value in complete_effect_positions_raw]
        if len(complete_effect_positions) != len(set(complete_effect_positions)):
            raise ContractError(
                f"Cooccurrence complete-read effect-supported positions contain duplicates at {key}"
            )
        if set(complete_effect_positions) != recomputed_effect_positions:
            raise ContractError(
                f"Cooccurrence complete-read effect-supported positions drift at {key}"
            )
        declared_complete_effect_count = optional_int(
            row["joint_signature_n_complete_marker_effect_supported"],
            field_name="joint_signature_n_complete_marker_effect_supported",
        )
        if declared_complete_effect_count != len(complete_effect_positions):
            raise ContractError(
                f"Cooccurrence complete-read effect-supported position-count drift at {key}"
            )
        same_complete_supported_positions = sorted(
            set(top_confirmed_positions).intersection(complete_effect_positions)
        )
        same_complete_spaced_positions = spatially_separated_positions(
            same_complete_supported_positions, minimum_separation=20
        )
        same_complete_fields = (
            (
                "n_same_complete_read_effect_supported_top_pair_by",
                len(same_complete_supported_positions),
            ),
            (
                "same_complete_read_effect_supported_top_pair_by_positions",
                same_complete_supported_positions,
            ),
            (
                "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp",
                len(same_complete_spaced_positions),
            ),
            (
                "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp",
                same_complete_spaced_positions,
            ),
        )
        for field_name, expected_value in same_complete_fields:
            observed_value = (
                optional_int(row[field_name], field_name=field_name)
                if field_name.startswith("n_")
                else optional_json(row[field_name], field_name=field_name)
            )
            if observed_value != expected_value:
                raise ContractError(
                    f"Cooccurrence same-complete-read marker drift at {key}: {field_name}"
                )
        joint_testable = parse_bool(
            row["joint_signature_testable"], field_name="joint_signature_testable"
        )
        joint_permutable = parse_bool(
            row["joint_signature_permutable"], field_name="joint_signature_permutable"
        )
        joint_pass = parse_bool(
            row["joint_signature_sensitivity_pass"],
            field_name="joint_signature_sensitivity_pass",
        )
        joint_global_by = parse_bool(
            row["joint_signature_global_by_discovery"],
            field_name="joint_signature_global_by_discovery",
        )
        joint_p = optional_float(
            row["joint_signature_p_conditional_perm"],
            field_name="joint_signature_p_conditional_perm",
        )
        joint_permutations = optional_int(
            row["joint_signature_permutations"],
            field_name="joint_signature_permutations",
        )
        if joint_p is not None and not 0.0 <= joint_p <= 1.0:
            raise ContractError(f"Joint-signature p-value is outside [0,1] at {key}")
        if joint_permutable:
            if joint_permutations != PAIR_CONDITIONAL_PERMUTATIONS:
                raise ContractError(f"Joint-signature permutation-count drift at {key}")
            if str(row["joint_signature_conditional_status"]) != "PERMUTABLE":
                raise ContractError(f"Joint-signature conditional-status drift at {key}")
        elif joint_permutations != 0 or joint_p is not None:
            raise ContractError(f"Nonpermutable joint signature carries conditional evidence at {key}")
        expected_joint_pass = bool(
            joint_testable
            and joint_permutable
            and joint_permutations == PAIR_CONDITIONAL_PERMUTATIONS
            and joint_p is not None
            and joint_p <= JOINT_SIGNATURE_P_MAX
        )
        if joint_pass != expected_joint_pass:
            raise ContractError(f"Joint-signature sensitivity-pass drift at {key}")
        selected = parse_bool(
            row[FORMAL_SELECTION_COLUMN], field_name=FORMAL_SELECTION_COLUMN
        )
        expected_selected = bool(
            m2
            and n_spaced >= 2
            and len(top_spaced_positions) >= 2
            and len(same_complete_spaced_positions) >= 2
            and joint_global_by
        )
        if selected != expected_selected:
            raise ContractError(f"G2 formal selection drift at {key}")
        if n_confirmed >= 1:
            pre_m2_formal_pair_sites.add(key)
        if selected:
            g2_sites.add(key)
        for legacy_field in (DEFAULT_SELECTION_COLUMN, BY_SELECTION_COLUMN):
            if legacy_field in row and parse_bool(
                row[legacy_field], field_name=legacy_field
            ) != selected:
                raise ContractError(f"Deprecated selection alias drift at {key}: {legacy_field}")
    validate_joint_signature_global_gates(rows)
    return rows, pre_m2_formal_pair_sites, g2_sites


def validate_joint_signature_global_gates(
    rows: Mapping[CandidateKey, Mapping[str, Any]],
) -> None:
    ordered = list(rows.items())
    p_values: list[float | None] = []
    for key, row in ordered:
        m2_eligible = parse_bool(row["m2_screen_eligible"], field_name="m2_screen_eligible")
        testable = parse_bool(
            row["joint_signature_testable"], field_name="joint_signature_testable"
        )
        permutable = parse_bool(
            row["joint_signature_permutable"], field_name="joint_signature_permutable"
        )
        permutations = optional_int(
            row["joint_signature_permutations"], field_name="joint_signature_permutations"
        )
        p_value = optional_float(
            row["joint_signature_p_conditional_perm"],
            field_name="joint_signature_p_conditional_perm",
        )
        if not m2_eligible:
            status = "INELIGIBLE_M2_SCREEN"
        elif not testable:
            status = "ELIGIBLE_M2_JOINT_SIGNATURE_NOT_TESTABLE"
        elif not permutable or permutations != PAIR_CONDITIONAL_PERMUTATIONS or p_value is None:
            status = "ELIGIBLE_M2_JOINT_SIGNATURE_NOT_PERMUTABLE"
        else:
            status = "ELIGIBLE_M2_JOINT_SIGNATURE_GLOBAL_FDR_FAMILY"
        if str(row["joint_signature_global_fdr_family_status"]) != status:
            raise ContractError(f"Joint-signature global FDR family drift at {key}")
        p_values.append(
            p_value
            if status == "ELIGIBLE_M2_JOINT_SIGNATURE_GLOBAL_FDR_FAMILY"
            else None
        )
    expected_bh = benjamini_hochberg(p_values)
    expected_by = benjamini_yekutieli(p_values)
    for index, (key, row) in enumerate(ordered):
        for field_name, expected_value in (
            ("joint_signature_q_global_bh", expected_bh[index]),
            ("joint_signature_q_global_by", expected_by[index]),
        ):
            require_optional_float_equal(
                row[field_name],
                expected_value,
                field_name=field_name,
                source=f"joint-signature global adjustment {key}",
            )
        expected_bh_discovery = bool(
            expected_bh[index] is not None and expected_bh[index] <= JOINT_SIGNATURE_P_MAX
        )
        expected_by_discovery = bool(
            expected_by[index] is not None and expected_by[index] <= JOINT_SIGNATURE_P_MAX
        )
        for field_name, expected_value in (
            ("joint_signature_global_bh_discovery", expected_bh_discovery),
            ("joint_signature_global_by_discovery", expected_by_discovery),
        ):
            if parse_bool(row[field_name], field_name=field_name) != expected_value:
                raise ContractError(f"Joint-signature global discovery drift at {key}: {field_name}")


def load_primary_artifact_audit(
    path: Path,
    bundle: InputBundle,
    assignment_keys: set[CandidateKey],
    *,
    label: str,
    expected_consumer_receipts: Sequence[Path],
) -> dict[str, Any]:
    payload = load_json(path, label)
    require_schema_pass(
        payload,
        "intersubmod.stable_primary_artifact_audit",
        label,
        schema_version=PRIMARY_ARTIFACT_AUDIT_SCHEMA_VERSION,
    )
    if payload.get("task_type") != "B_comprehensive_validation":
        raise ContractError(f"{label} task_type is not full scope")
    if payload.get("scope") != "complete_primary_stable_null_multigroup_set":
        raise ContractError(f"{label} scope drift")
    if payload.get("pass_semantics") != "execution_integrity_and_primary_artifact_identity_only":
        raise ContractError(f"{label} pass semantics drift")
    source_authority = current_release_source_authority()
    if payload.get("source_authority") != source_authority:
        raise ContractError(f"{label} source authority drift")
    expected_code = PRIMARY_AUDITOR.capture_source_identity()
    expected_modes = PRIMARY_AUDITOR.capture_source_modes()
    if payload.get("code") != expected_code:
        raise ContractError(f"{label} producer source identity drift")
    source_lock = payload.get("source_lock")
    if (
        not isinstance(source_lock, Mapping)
        or source_lock.get("source_identity_before") != expected_code
        or source_lock.get("source_identity_after_compute") != expected_code
        or source_lock.get("source_modes_before") != expected_modes
        or source_lock.get("source_modes_after_compute") != expected_modes
        or source_lock.get("all_sources_read_only_and_unchanged") is not True
        or set(expected_modes.values()) != {"0o444"}
    ):
        raise ContractError(f"{label} producer source lock drift")
    expected_command = PRIMARY_AUDITOR.canonical_command(
        site_results=bundle.screen_sites,
        assignments=bundle.screen_assignments,
        consumer_receipts=expected_consumer_receipts,
        output=path,
        **PRIMARY_AUDITOR.CANONICAL_PARAMETERS,
    )
    if payload.get("command") != expected_command:
        raise ContractError(f"{label} command is not canonical")
    if payload.get("execution") != PRIMARY_AUDITOR.CANONICAL_PARAMETERS:
        raise ContractError(f"{label} execution parameters are not canonical")
    if oct(path.stat().st_mode & 0o777) != "0o444":
        raise ContractError(f"{label} receipt is not mode 0444")
    interpretation = payload.get("scientific_interpretation")
    if not isinstance(interpretation, Mapping) or interpretation.get("is_negative_result") is not False:
        raise ContractError(f"{label} must declare is_negative_result=false")
    inputs = payload.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ContractError(f"{label} lacks inputs")
    verify_declared_artifact(inputs.get("site_results"), bundle.screen_sites, f"{label} site results")
    verify_declared_artifact(
        inputs.get("stable_assignments"),
        bundle.screen_assignments,
        f"{label} stable assignments",
    )
    consumer_receipts = inputs.get("consumer_receipts")
    if not isinstance(consumer_receipts, list):
        raise ContractError(f"{label} lacks consumer_receipts list")
    if len(consumer_receipts) != len(expected_consumer_receipts):
        raise ContractError(
            f"{label} consumer receipt count drift: "
            f"{len(consumer_receipts)} != {len(expected_consumer_receipts)}"
        )
    for index, (reference, expected_path) in enumerate(
        zip(consumer_receipts, expected_consumer_receipts, strict=True)
    ):
        verify_declared_artifact(
            reference,
            expected_path,
            f"{label} consumer receipt {index + 1}",
        )
    counts = payload.get("counts")
    if not isinstance(counts, Mapping):
        raise ContractError(f"{label} lacks counts")
    expected_sites = len(assignment_keys)
    expected_artifacts = expected_sites * 3
    expected_counts = {
        "stable_sites": expected_sites,
        "assignment_records": expected_sites,
        "primary_artifacts_expected": expected_artifacts,
        "primary_artifacts_verified": expected_artifacts,
    }
    for field_name, expected in expected_counts.items():
        if int(counts.get(field_name, -1)) != expected:
            raise ContractError(f"{label} count drift: {field_name}")
    verification = payload.get("verification")
    if not isinstance(verification, Mapping):
        raise ContractError(f"{label} lacks verification")
    expected_roles = ["distance_matrix", "methylation_matrix", "reads"]
    gates = {
        "stable_site_assignment_key_sets_exact": verification.get(
            "stable_site_assignment_key_sets_exact"
        )
        is True,
        "artifact_roles_exact": verification.get("artifact_roles_exact") == expected_roles,
        "path_size_sha256_verified": verification.get("path_size_sha256_verified") is True,
        "errors": int(verification.get("errors", -1)) == 0,
    }
    failed = sorted(name for name, passed in gates.items() if not passed)
    if failed:
        raise ContractError(f"{label} failed gates: {failed}")
    artifact_set_sha256 = str(verification.get("artifact_set_sha256", ""))
    if len(artifact_set_sha256) != 64 or any(
        character not in "0123456789abcdef" for character in artifact_set_sha256.lower()
    ):
        raise ContractError(f"{label} artifact_set_sha256 is invalid")
    return payload


def parse_utc_timestamp(value: Any, *, label: str) -> datetime:
    if not isinstance(value, str) or not value.strip():
        raise ContractError(f"{label} timestamp is missing")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise ContractError(f"{label} timestamp is invalid: {value!r}") from error
    if parsed.tzinfo is None:
        raise ContractError(f"{label} timestamp lacks timezone")
    return parsed.astimezone(timezone.utc)


def validate_primary_artifact_audit_window(
    pre_audit: Mapping[str, Any],
    post_audit: Mapping[str, Any],
    producer_receipts: Mapping[str, Mapping[str, Any]],
) -> None:
    def audit_window(
        audit: Mapping[str, Any], *, label: str
    ) -> tuple[datetime, datetime]:
        started = parse_utc_timestamp(
            audit.get("started_at_utc"), label=f"{label} start"
        )
        finished = parse_utc_timestamp(
            audit.get("finished_at_utc"), label=f"{label} finish"
        )
        created = parse_utc_timestamp(
            audit.get("created_at_utc"), label=f"{label} creation"
        )
        if started > finished:
            raise ContractError(f"{label} execution time window is reversed")
        if created != finished:
            raise ContractError(f"{label} created_at_utc must equal finished_at_utc")
        return started, finished

    _, pre_completed = audit_window(
        pre_audit, label="pre-downstream primary artifact audit"
    )
    post_started, _ = audit_window(
        post_audit, label="post-downstream primary artifact audit"
    )
    if pre_completed > post_started:
        raise ContractError("Primary artifact audit time window is reversed")
    for producer, receipt in producer_receipts.items():
        started = parse_utc_timestamp(
            receipt.get("started_at_utc"), label=f"{producer} start"
        )
        finished = parse_utc_timestamp(
            receipt.get("finished_at_utc"), label=f"{producer} finish"
        )
        if started > finished:
            raise ContractError(f"{producer} execution time window is reversed")
        if started < pre_completed or finished > post_started:
            raise ContractError(
                f"{producer} did not execute inside the primary artifact pre/post audit window"
            )


def pair_key(row: Mapping[str, Any], *, source: str) -> PairKey:
    try:
        focal = candidate_key(
            {
                "sample": row["sample"],
                "chrom": row["focal_chrom"],
                "pos": row["focal_pos"],
                "ref": row["focal_ref"],
                "alt": row["focal_alt"],
            },
            source=source,
        )
        partner_chrom = str(row["partner_chrom"]).strip()
        partner_pos = int(row["partner_pos"])
        partner_ref = str(row["partner_ref"]).strip().upper()
        partner_alt = str(row["partner_alt"]).strip().upper()
    except (KeyError, TypeError, ValueError) as error:
        raise ContractError(f"Malformed pair key in {source}") from error
    if not partner_chrom or partner_pos <= 0 or not partner_ref or not partner_alt:
        raise ContractError(f"Invalid partner identity in {source}")
    return (*focal, partner_chrom, partner_pos, partner_ref, partner_alt)


def cross_platform_key(key: PairKey) -> CrossPlatformPairKey:
    return key[1:]


def validate_four_state_pair_row(row: Mapping[str, Any], key: PairKey) -> None:
    counts = optional_json(row["endpoint_b_state_counts"], field_name="endpoint_b_state_counts")
    if not isinstance(counts, Mapping):
        raise ContractError(f"Four-state counts are not an object at {key}")
    expected = recompute_four_state_from_counts(counts)
    integer_fields = {
        "endpoint_b_n_joint": "n_joint",
        "endpoint_b_n_called_depth": "n_called_depth",
        "endpoint_b_n_focal_ref": "n_focal_ref",
        "endpoint_b_n_focal_alt": "n_focal_alt",
        "endpoint_b_n_partner_ref": "n_partner_ref",
        "endpoint_b_n_partner_alt": "n_partner_alt",
        "endpoint_b_minimum_zero_violation_depth": "minimum_zero_violation_depth",
        "endpoint_b_n_compatible_relation_models": "n_compatible_relation_models",
    }
    for field_name, expected_name in integer_fields.items():
        observed = optional_int(row[field_name], field_name=field_name)
        if observed != expected[expected_name]:
            raise ContractError(
                f"Four-state integer drift at {key}: {field_name}={observed}, "
                f"expected={expected[expected_name]}"
            )
    float_fields = {
        "endpoint_b_focal_ancestor_violation_rate": expected[
            "focal_ancestor_violation_rate"
        ],
        "endpoint_b_partner_ancestor_violation_rate": expected[
            "partner_ancestor_violation_rate"
        ],
        "endpoint_b_error_ceiling": expected["error_ceiling"],
        "endpoint_b_error_model_confidence": expected["confidence"],
        "endpoint_b_familywise_confidence": expected["familywise_confidence"],
        "endpoint_b_focal_ancestor_violation_p_exact": expected[
            "focal_ancestor_violation"
        ]["p_exact_greater"],
        "endpoint_b_focal_ancestor_violation_upper_bound": expected[
            "focal_ancestor_violation"
        ]["upper_bound"],
        "endpoint_b_focal_ancestor_violation_threshold": expected[
            "focal_ancestor_violation"
        ]["threshold"],
        "endpoint_b_partner_ancestor_violation_p_exact": expected[
            "partner_ancestor_violation"
        ]["p_exact_greater"],
        "endpoint_b_partner_ancestor_violation_upper_bound": expected[
            "partner_ancestor_violation"
        ]["upper_bound"],
        "endpoint_b_partner_ancestor_violation_threshold": expected[
            "partner_ancestor_violation"
        ]["threshold"],
        "endpoint_b_branching_violation_p_exact": expected["branching_violation"][
            "p_exact_greater"
        ],
        "endpoint_b_branching_violation_upper_bound": expected["branching_violation"][
            "upper_bound"
        ],
        "endpoint_b_branching_violation_threshold": expected["branching_violation"][
            "threshold"
        ],
    }
    for field_name, expected_value in float_fields.items():
        require_optional_float_equal(
            row[field_name], expected_value, field_name=field_name, source=f"four-state {key}"
        )
    if optional_int(
        row["endpoint_b_relation_family_size"],
        field_name="endpoint_b_relation_family_size",
    ) != expected["relation_family_size"]:
        raise ContractError(f"Four-state relation-family size drift at {key}")
    if str(row["endpoint_b_multiplicity_method"]) != expected["multiplicity_method"]:
        raise ContractError(f"Four-state multiplicity-method drift at {key}")
    status_fields = {
        "endpoint_b_status": expected["status"],
        "endpoint_b_focal_ancestor_violation_status": expected[
            "focal_ancestor_violation"
        ]["status"],
        "endpoint_b_partner_ancestor_violation_status": expected[
            "partner_ancestor_violation"
        ]["status"],
        "endpoint_b_branching_violation_status": expected["branching_violation"]["status"],
        "endpoint_b_relation_compatibility": expected["relation_compatibility"],
    }
    for field_name, expected_value in status_fields.items():
        if str(row[field_name]) != expected_value:
            raise ContractError(
                f"Four-state status drift at {key}: {field_name}={row[field_name]!r}, "
                f"expected={expected_value!r}"
            )
    observed_complete = parse_bool(
        row["endpoint_b_complete_four_state_testable"],
        field_name="endpoint_b_complete_four_state_testable",
    )
    if observed_complete != expected["complete_four_state_testable"]:
        raise ContractError(f"Four-state testability drift at {key}")
    observed_models = optional_json(
        row["endpoint_b_compatible_relation_models"],
        field_name="endpoint_b_compatible_relation_models",
    )
    if observed_models != expected["compatible_relation_models"]:
        raise ContractError(
            f"Four-state compatible-model drift at {key}: "
            f"observed={observed_models!r} expected={expected['compatible_relation_models']!r}"
        )


def validate_pair_global_gates(
    rows: Mapping[PairKey, Mapping[str, Any]],
    cooccurrence_rows: Mapping[CandidateKey, Mapping[str, Any]],
) -> None:
    ordered = list(rows.items())
    endpoint_p: list[float | None] = []
    expected_family_statuses: list[str] = []
    for key, row in ordered:
        focal_key: CandidateKey = key[:5]
        site = cooccurrence_rows.get(focal_key)
        if site is None:
            raise ContractError(f"Pair focal site is absent from cooccurrence site table: {key}")
        m2_eligible = parse_bool(
            site["m2_screen_eligible"], field_name="m2_screen_eligible"
        )
        testable = parse_bool(row["endpoint_a_testable"], field_name="endpoint_a_testable")
        exact_status = str(row["endpoint_a_exact_state_space_status"])
        exact_p = optional_float(
            row["endpoint_a_p_fixed_margins_exact"],
            field_name="endpoint_a_p_fixed_margins_exact",
        )
        if not m2_eligible:
            expected_status = "INELIGIBLE_M2_SCREEN"
        elif not testable:
            expected_status = "ELIGIBLE_M2_ENDPOINT_A_NOT_TESTABLE"
        elif exact_status != "EXACT_ENUMERATED" or exact_p is None:
            expected_status = "ELIGIBLE_M2_EXACT_NOT_IDENTIFIABLE"
        else:
            expected_status = "ELIGIBLE_M2_EXACT_FAMILY"
        if str(row["endpoint_a_global_fdr_family_status"]) != expected_status:
            raise ContractError(f"Endpoint-A global FDR family membership drift at {key}")
        expected_family_statuses.append(expected_status)
        endpoint_p.append(
            exact_p if expected_status == "ELIGIBLE_M2_EXACT_FAMILY" else None
        )
    callability_p = [
        optional_float(row["callability_p_analytic"], field_name="callability_p_analytic")
        for _, row in ordered
    ]
    expected_endpoint_bh = benjamini_hochberg(endpoint_p)
    expected_endpoint_by = benjamini_yekutieli(endpoint_p)
    expected_callability_bh = benjamini_hochberg(callability_p)
    expected_callability_by = benjamini_yekutieli(callability_p)
    for index, (key, row) in enumerate(ordered):
        for field_name, expected_value in (
            ("endpoint_a_q_global_bh", expected_endpoint_bh[index]),
            ("endpoint_a_q_global_by", expected_endpoint_by[index]),
            ("callability_q_global_bh", expected_callability_bh[index]),
            ("callability_q_global_by", expected_callability_by[index]),
        ):
            require_optional_float_equal(
                row[field_name],
                expected_value,
                field_name=field_name,
                source=f"pair global adjustment {key}",
            )
        testable = parse_bool(row["endpoint_a_testable"], field_name="endpoint_a_testable")
        exact_family = expected_family_statuses[index] == "ELIGIBLE_M2_EXACT_FAMILY"
        expected_exact_bh = bool(
            testable
            and exact_family
            and expected_endpoint_bh[index] is not None
            and expected_endpoint_bh[index] <= PAIR_GLOBAL_Q_MAX
        )
        expected_exact_by = bool(
            testable
            and exact_family
            and expected_endpoint_by[index] is not None
            and expected_endpoint_by[index] <= PAIR_GLOBAL_Q_MAX
        )
        effect = optional_float(row["endpoint_a_cramers_v"], field_name="endpoint_a_cramers_v")
        delta = optional_float(
            row["endpoint_a_delta_alt_fraction"], field_name="endpoint_a_delta_alt_fraction"
        )
        expected_effect = bool(
            testable
            and effect is not None
            and effect >= COMPLETE_MARKER_CRAMERS_V_MIN
            and delta is not None
            and delta >= COMPLETE_MARKER_DELTA_ALT_MIN
        )
        core_counts = optional_json(
            row["core_partner_call_counts"], field_name="core_partner_call_counts"
        )
        if not isinstance(core_counts, Mapping):
            raise ContractError(f"Core partner call counts are not an object at {key}")
        try:
            noncallable = int(core_counts.get("O", 0)) + int(core_counts.get("X", 0))
        except (TypeError, ValueError) as error:
            raise ContractError(f"Invalid core callability counts at {key}") from error
        declared_noncallable = optional_int(
            row["callability_noncallable_core_reads"],
            field_name="callability_noncallable_core_reads",
        )
        if declared_noncallable != noncallable:
            raise ContractError(f"Callability noncallable-read count drift at {key}")
        callability_testable = parse_bool(
            row["callability_testable"], field_name="callability_testable"
        )
        callability_v = optional_float(
            row["callability_cramers_v"], field_name="callability_cramers_v"
        )
        callability_by = expected_callability_by[index]
        if noncallable == 0:
            expected_callability_status = "PASS_ALL_CORE_READS_CALLABLE"
            expected_callability_pass = True
        elif not callability_testable:
            expected_callability_status = "NOT_IDENTIFIABLE_DIFFERENTIAL_CALLABILITY"
            expected_callability_pass = False
        elif callability_by is None or callability_v is None:
            expected_callability_status = "NOT_IDENTIFIABLE_MISSING_CALLABILITY_STATISTIC"
            expected_callability_pass = False
        elif callability_by <= PAIR_GLOBAL_Q_MAX or callability_v >= COMPLETE_MARKER_CRAMERS_V_MIN:
            expected_callability_status = "FAIL_DIFFERENTIAL_CALLABILITY_SIGNAL"
            expected_callability_pass = False
        else:
            expected_callability_status = "PASS_NO_STRONG_DIFFERENTIAL_CALLABILITY_DETECTED"
            expected_callability_pass = True
        if str(row["callability_gate_status"]) != expected_callability_status:
            raise ContractError(f"Callability status drift at {key}")
        if parse_bool(row["callability_gate_pass"], field_name="callability_gate_pass") != expected_callability_pass:
            raise ContractError(f"Callability gate drift at {key}")
        conditional_p = optional_float(
            row["endpoint_a_p_conditional_perm"],
            field_name="endpoint_a_p_conditional_perm",
        )
        conditional_permutations = optional_int(
            row["endpoint_a_permutations"], field_name="endpoint_a_permutations"
        )
        conditional_permutable = parse_bool(
            row["endpoint_a_permutable"], field_name="endpoint_a_permutable"
        )
        if conditional_permutable:
            if conditional_permutations != PAIR_CONDITIONAL_PERMUTATIONS:
                raise ContractError(f"Conditional permutation count drift at {key}")
            if str(row["endpoint_a_conditional_status"]) != "PERMUTABLE":
                raise ContractError(f"Conditional status drift at {key}")
        elif conditional_permutations != 0 or conditional_p is not None:
            raise ContractError(f"Nonpermutable pair carries conditional evidence at {key}")
        expected_conditional = bool(
            conditional_permutable
            and conditional_permutations == PAIR_CONDITIONAL_PERMUTATIONS
            and conditional_p is not None
            and conditional_p <= PAIR_CONDITIONAL_P_MAX
        )
        expected_formal = bool(
            expected_exact_by
            and expected_effect
            and expected_conditional
            and expected_callability_pass
        )
        boolean_expectations = {
            "endpoint_a_effect_gate_pass": expected_effect,
            "endpoint_a_exact_bh_discovery": expected_exact_bh,
            "endpoint_a_exact_by_discovery": expected_exact_by,
            "endpoint_a_conditional_sensitivity_pass": expected_conditional,
            "endpoint_a_formal_pair_by_confirmed": expected_formal,
            "endpoint_a_confirmed_association": expected_formal,
            "endpoint_a_confirmed_by_sensitivity": expected_formal,
        }
        for field_name, expected_value in boolean_expectations.items():
            if parse_bool(row[field_name], field_name=field_name) != expected_value:
                raise ContractError(f"Formal G1 gate drift at {key}: {field_name}")


def load_pair_data(
    path: Path,
    stable_keys: set[CandidateKey],
    cooccurrence_rows: Mapping[CandidateKey, Mapping[str, str]],
) -> PairData:
    required = {
        "sample",
        "focal_chrom",
        "focal_pos",
        "focal_ref",
        "focal_alt",
        "partner_chrom",
        "partner_pos",
        "partner_ref",
        "partner_alt",
        "distance_bp",
        "core_partner_call_counts",
        "all_partner_call_counts",
        "endpoint_a_testable",
        "endpoint_a_p_fixed_margins_exact",
        "endpoint_a_exact_state_space_status",
        "endpoint_a_global_fdr_family_status",
        "endpoint_a_q_global_bh",
        "endpoint_a_q_global_by",
        "endpoint_a_cramers_v",
        "endpoint_a_delta_alt_fraction",
        "endpoint_a_effect_gate_pass",
        "endpoint_a_exact_bh_discovery",
        "endpoint_a_exact_by_discovery",
        "endpoint_a_p_conditional_perm",
        "endpoint_a_permutations",
        "endpoint_a_permutable",
        "endpoint_a_conditional_status",
        "endpoint_a_conditional_sensitivity_pass",
        "endpoint_a_formal_pair_by_confirmed",
        "endpoint_a_confirmed_association",
        "endpoint_a_confirmed_by_sensitivity",
        "callability_testable",
        "callability_p_analytic",
        "callability_q_global_bh",
        "callability_q_global_by",
        "callability_cramers_v",
        "callability_noncallable_core_reads",
        "callability_gate_status",
        "callability_gate_pass",
        "endpoint_b_status",
        "endpoint_b_n_joint",
        "endpoint_b_n_called_depth",
        "endpoint_b_state_counts",
        "endpoint_b_n_focal_ref",
        "endpoint_b_n_focal_alt",
        "endpoint_b_n_partner_ref",
        "endpoint_b_n_partner_alt",
        "endpoint_b_focal_ancestor_violation_rate",
        "endpoint_b_partner_ancestor_violation_rate",
        "endpoint_b_error_ceiling",
        "endpoint_b_error_model_confidence",
        "endpoint_b_familywise_confidence",
        "endpoint_b_relation_family_size",
        "endpoint_b_multiplicity_method",
        "endpoint_b_minimum_zero_violation_depth",
        "endpoint_b_focal_ancestor_violation_p_exact",
        "endpoint_b_focal_ancestor_violation_upper_bound",
        "endpoint_b_focal_ancestor_violation_threshold",
        "endpoint_b_focal_ancestor_violation_status",
        "endpoint_b_partner_ancestor_violation_p_exact",
        "endpoint_b_partner_ancestor_violation_upper_bound",
        "endpoint_b_partner_ancestor_violation_threshold",
        "endpoint_b_partner_ancestor_violation_status",
        "endpoint_b_branching_violation_p_exact",
        "endpoint_b_branching_violation_upper_bound",
        "endpoint_b_branching_violation_threshold",
        "endpoint_b_branching_violation_status",
        "endpoint_b_complete_four_state_testable",
        "endpoint_b_relation_compatibility",
        "endpoint_b_compatible_relation_models",
        "endpoint_b_n_compatible_relation_models",
        "top_coverage_marker",
        "focal_truth_label",
        "partner_truth_label",
        "focal_ssnv_branch",
        "partner_ssnv_branch",
        "focal_component_id",
        "partner_component_id",
        "topology_scope",
        "topology_region",
        "topology_order_status",
        "cross_platform_exact_pair_present",
        "cross_platform_biological_n",
        "cross_platform_conditional_by_effect_direction_replication_pass",
        "cross_platform_replication_status",
    }
    aggregates: dict[CandidateKey, PairAggregate] = defaultdict(PairAggregate)
    seen: set[PairKey] = set()
    rows: dict[PairKey, dict[str, str]] = {}
    by_focal: dict[CandidateKey, list[dict[str, str]]] = defaultdict(list)
    hcc_rows: dict[
        CrossPlatformPairKey, dict[str, tuple[bool, bool, bool, str]]
    ] = defaultdict(dict)
    n_testable = 0
    n_confirmed = 0
    n_confirmed_by = 0
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader, required, "cooccurrence pair table")
        for row_number, row in enumerate(reader, 2):
            key = pair_key(row, source=f"cooccurrence pair row {row_number}")
            if key in seen:
                raise ContractError(f"Duplicate cooccurrence pair key: {key}")
            seen.add(key)
            rows[key] = dict(row)
            focal_key = key[:5]
            if focal_key not in stable_keys:
                raise ContractError(f"Pair focal is absent from stable set: {focal_key}")
            testable = parse_bool(row["endpoint_a_testable"], field_name="endpoint_a_testable")
            formal = parse_bool(
                row["endpoint_a_formal_pair_by_confirmed"],
                field_name="endpoint_a_formal_pair_by_confirmed",
            )
            confirmed = parse_bool(
                row["endpoint_a_confirmed_association"],
                field_name="endpoint_a_confirmed_association",
            )
            confirmed_by = parse_bool(
                row["endpoint_a_confirmed_by_sensitivity"],
                field_name="endpoint_a_confirmed_by_sensitivity",
            )
            if confirmed != formal or confirmed_by != formal:
                raise ContractError(f"Deprecated pair-confirmation alias drift: {key}")
            exact_by = parse_bool(
                row["endpoint_a_exact_by_discovery"],
                field_name="endpoint_a_exact_by_discovery",
            )
            conditional_pass = parse_bool(
                row["endpoint_a_conditional_sensitivity_pass"],
                field_name="endpoint_a_conditional_sensitivity_pass",
            )
            effect = optional_float(
                row["endpoint_a_cramers_v"], field_name="endpoint_a_cramers_v"
            )
            delta = optional_float(
                row["endpoint_a_delta_alt_fraction"],
                field_name="endpoint_a_delta_alt_fraction",
            )
            if formal and not (
                exact_by
                and conditional_pass
                and effect is not None
                and effect >= 0.30
                and delta is not None
                and delta >= 0.50
            ):
                raise ContractError(f"Formal G1 pair gate drift: {key}")
            complete = parse_bool(
                row["endpoint_b_complete_four_state_testable"],
                field_name="endpoint_b_complete_four_state_testable",
            )
            if formal and not testable:
                raise ContractError(f"Confirmed pair is not testable: {key}")
            exact_p = optional_float(
                row["endpoint_a_p_fixed_margins_exact"],
                field_name="endpoint_a_p_fixed_margins_exact",
            )
            exact_testable = bool(
                testable
                and exact_p is not None
                and row["endpoint_a_global_fdr_family_status"]
                == "ELIGIBLE_M2_EXACT_FAMILY"
            )
            aggregate = aggregates[focal_key]
            aggregate.n_pairs += 1
            aggregate.n_testable_pairs += int(testable)
            aggregate.n_exact_testable_pairs += int(exact_testable)
            aggregate.n_confirmed_pairs += int(formal)
            aggregate.n_confirmed_pairs_by += int(formal)
            aggregate.n_formal_pair_by_confirmed += int(formal)
            aggregate.n_confirmed_complete_four_state_pairs += int(formal and complete)
            relation = str(row["endpoint_b_relation_compatibility"])
            if relation not in FOUR_STATE_RELATIONS:
                raise ContractError(
                    f"Four-state relation is outside the producer enum at {key}: {relation!r}"
                )
            validate_four_state_pair_row(row, key)
            same_pair_witness = bool(formal and complete and relation in COMPATIBLE_RELATIONS)
            aggregate.n_confirmed_compatible_relation_pairs += int(same_pair_witness)
            aggregate.n_same_pair_four_state_witnesses += int(same_pair_witness)
            topology = str(row["topology_order_status"])
            aggregate.n_confirmed_topology_inferable_pairs += int(
                formal
                and topology
                and not topology.startswith("NOT_INFERABLE")
                and topology != "None"
            )
            n_testable += int(testable)
            n_confirmed += int(formal)
            n_confirmed_by += int(formal)
            by_focal[focal_key].append(dict(row))
            if key[0] in {"HCC1395", "HCC1395_DORADO"}:
                cross_key = cross_platform_key(key)
                if key[0] in hcc_rows[cross_key]:
                    raise ContractError(f"Duplicate cross-platform pair for {key[0]}: {cross_key}")
                declared_exact = parse_bool(
                    row["cross_platform_exact_pair_present"],
                    field_name="cross_platform_exact_pair_present",
                )
                replicated = parse_bool(
                    row[
                        "cross_platform_conditional_by_effect_direction_replication_pass"
                    ],
                    field_name=(
                        "cross_platform_conditional_by_effect_direction_replication_pass"
                    ),
                )
                biological_n = optional_int(
                    row["cross_platform_biological_n"],
                    field_name="cross_platform_biological_n",
                )
                if declared_exact and biological_n != 1:
                    raise ContractError(f"HCC technical pair must declare biological_n=1: {key}")
                hcc_rows[cross_key][key[0]] = (
                    formal,
                    declared_exact,
                    replicated,
                    str(row["cross_platform_replication_status"]),
                )
    validate_pair_global_gates(rows, cooccurrence_rows)
    exact_pairs: set[CrossPlatformPairKey] = set()
    both_confirmed_pairs: set[CrossPlatformPairKey] = set()
    replicated_pairs: set[CrossPlatformPairKey] = set()
    for cross_key, by_sample in hcc_rows.items():
        exact = set(by_sample) == {"HCC1395", "HCC1395_DORADO"}
        both_confirmed = exact and all(values[0] for values in by_sample.values())
        replicated = exact and all(values[2] for values in by_sample.values())
        if replicated and not both_confirmed:
            raise ContractError(
                f"Cross-platform replication passes without the same exact G1 pair: {cross_key}"
            )
        if exact:
            exact_pairs.add(cross_key)
        if both_confirmed:
            both_confirmed_pairs.add(cross_key)
        if replicated:
            replicated_pairs.add(cross_key)
        for sample, (formal, declared_exact, declared_replicated, declared_status) in by_sample.items():
            if declared_exact != exact:
                raise ContractError(
                    f"Cross-platform exact-pair flag drift for {sample}/{cross_key}"
                )
            if not declared_status:
                raise ContractError(f"Cross-platform status is empty for {sample}/{cross_key}")
            focal_key = (sample, *cross_key[:4])
            aggregate = aggregates[focal_key]
            aggregate.n_cross_platform_exact_pair_rows += int(exact)
            aggregate.n_cross_platform_both_confirmed_pair_rows += int(both_confirmed)
            if declared_replicated != replicated:
                raise ContractError(
                    f"Cross-platform replication flag drift for {sample}/{cross_key}"
                )
    for focal_key in stable_keys:
        site = cooccurrence_rows[focal_key]
        focal_rows = by_focal.get(focal_key, [])
        if int(site["n_pair_rows_reconciled"]) != len(focal_rows):
            raise ContractError(f"Pair-site exact row-count reconciliation failed at {focal_key}")
        observed_positions = sorted(
            int(row["partner_pos"])
            for row in focal_rows
            if parse_bool(
                row["endpoint_a_formal_pair_by_confirmed"],
                field_name="endpoint_a_formal_pair_by_confirmed",
            )
        )
        declared_positions = optional_json(
            site["pair_by_confirmed_positions"], field_name="pair_by_confirmed_positions"
        )
        if declared_positions != observed_positions:
            raise ContractError(f"Pair-site formal position reconciliation failed at {focal_key}")
        observed_top_positions = sorted(
            int(row["partner_pos"])
            for row in focal_rows
            if parse_bool(
                row["top_coverage_marker"], field_name="top_coverage_marker"
            )
        )
        declared_top_positions = optional_json(
            site["top_marker_positions"], field_name="top_marker_positions"
        )
        if not isinstance(declared_top_positions, list):
            raise ContractError(f"Site top-marker positions are not a list at {focal_key}")
        if sorted(int(value) for value in declared_top_positions) != observed_top_positions:
            raise ContractError(f"Pair-site top-marker reconciliation failed at {focal_key}")
    return PairData(
        aggregates=dict(aggregates),
        n_rows=len(seen),
        n_testable=n_testable,
        n_confirmed=n_confirmed,
        n_confirmed_by=n_confirmed_by,
        exact_cross_platform_pairs=exact_pairs,
        both_confirmed_cross_platform_pairs=both_confirmed_pairs,
        replicated_cross_platform_pairs=replicated_pairs,
        rows=rows,
        by_focal={key: list(value) for key, value in by_focal.items()},
    )


def validate_cooccurrence_input_identity_and_m2_provenance(
    summary: Mapping[str, Any], receipt: Mapping[str, Any]
) -> None:
    policy = summary.get("frozen_input_identity_policy")
    if not isinstance(policy, Mapping) or receipt.get(
        "frozen_manifest_input_identity_policy"
    ) != policy:
        raise ContractError("Cooccurrence frozen-input identity policy is missing or drifted")
    if set(policy.get("hash_required_fields") or []) != MANIFEST_HASH_REQUIRED_FIELDS:
        raise ContractError("Cooccurrence hash-required manifest roles drift")
    if set(policy.get("explicit_large_metadata_identity_fields") or []) != (
        MANIFEST_LARGE_METADATA_IDENTITY_FIELDS
    ):
        raise ContractError("Cooccurrence explicit large-artifact identity roles drift")
    artifacts = policy.get("artifacts")
    if not isinstance(artifacts, list):
        raise ContractError("Cooccurrence frozen-input identity inventory is malformed")
    expected_keys = {
        (sample, field_name)
        for sample in DATASETS
        for field_name in (
            MANIFEST_HASH_REQUIRED_FIELDS | MANIFEST_LARGE_METADATA_IDENTITY_FIELDS
        )
    }
    observed_keys: set[tuple[str, str]] = set()
    n_sha = 0
    n_large_metadata = 0
    for row in artifacts:
        if not isinstance(row, Mapping):
            raise ContractError("Cooccurrence frozen-input identity row is malformed")
        key = (str(row.get("sample")), str(row.get("field")))
        if key in observed_keys:
            raise ContractError(f"Duplicate cooccurrence frozen-input identity row: {key}")
        observed_keys.add(key)
        mode = str(row.get("identity_mode"))
        if key[1] in MANIFEST_HASH_REQUIRED_FIELDS:
            if mode != "sha256_size_path_v1" or not isinstance(row.get("sha256"), str):
                raise ContractError(f"Required SHA-256 identity missing for {key}")
            n_sha += 1
        elif key[1] in MANIFEST_LARGE_METADATA_IDENTITY_FIELDS:
            if mode == "sha256_size_path_v1":
                n_sha += 1
            elif mode == "explicit_large_artifact_size_mtime_path_v1":
                if row.get("sha256") is not None:
                    raise ContractError(f"Metadata-only identity falsely declares SHA-256 for {key}")
                n_large_metadata += 1
            else:
                raise ContractError(f"Invalid large-artifact identity mode for {key}: {mode}")
        if int(row.get("size_bytes", 0)) <= 0 or int(row.get("mtime_ns", 0)) <= 0:
            raise ContractError(f"Incomplete frozen-input size/mtime identity for {key}")
    if observed_keys != expected_keys:
        raise ContractError("Cooccurrence frozen-input identity inventory coverage drift")
    if int(policy.get("n_sha256_identity_artifacts", -1)) != n_sha or int(
        policy.get("n_explicit_large_size_mtime_identity_artifacts", -1)
    ) != n_large_metadata:
        raise ContractError("Cooccurrence frozen-input identity mode counts drift")
    if "weaker than full SHA-256" not in str(policy.get("claim_guardrail", "")):
        raise ContractError("Cooccurrence large-artifact identity limitation is missing")

    provenance = summary.get("m2_axis_statistic_provenance")
    if not isinstance(provenance, Mapping) or receipt.get(
        "m2_axis_statistic_provenance"
    ) != provenance:
        raise ContractError("Cooccurrence M2 axis-statistic provenance is missing or drifted")
    if (
        provenance.get("downstream_raw_read_axis_statistic_recomputed") is not False
        or provenance.get("downstream_axis_classification_recomputed") is not True
        or set(provenance.get("downstream_validation") or [])
        != {
            "axis sample-size reconciliation",
            "499-permutation add-one p-value grid",
            "effect threshold classification",
            "80-percent planning-power evaluability",
        }
        or "not an independent raw-read remeasurement"
        not in str(provenance.get("claim_guardrail", ""))
    ):
        raise ContractError("Cooccurrence M2 axis-statistic provenance contract drift")


def projection_digest(projections: Iterable[tuple[Any, ...]]) -> str:
    digest = hashlib.sha256()
    for projection in sorted(projections):
        digest.update(compact_json(list(projection)).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def validate_cooccurrence_source_and_preflight_contract(
    bundle: InputBundle,
    summary: Mapping[str, Any],
    receipt: Mapping[str, Any],
    preflight: Mapping[str, Any],
) -> None:
    source_authority = current_release_source_authority()
    require_schema_pass(
        preflight,
        "intersubmod.cooccurrence_task_contract_preflight",
        "cooccurrence raw-identity preflight",
        schema_version=COOCCURRENCE_PREFLIGHT_SCHEMA_VERSION,
    )
    if (
        preflight.get("task_type") != "B_comprehensive_validation"
        or preflight.get("raw_identity_release_contract")
        != RAW_IDENTITY_RELEASE_CONTRACT
        or receipt.get("raw_identity_release_contract")
        != RAW_IDENTITY_RELEASE_CONTRACT
        or summary.get("raw_identity_release_contract")
        != RAW_IDENTITY_RELEASE_CONTRACT
    ):
        raise ContractError("Cooccurrence raw-identity release contract drift")
    if (
        preflight.get("source_authority") != source_authority
        or receipt.get("source_authority") != source_authority
    ):
        raise ContractError("Cooccurrence preflight/producer source authority drift")
    checks = preflight.get("checks")
    if not isinstance(checks, Mapping) or not checks or not all(
        value is True for value in checks.values()
    ):
        raise ContractError("Cooccurrence preflight contains a failed or malformed check")

    preflight_inputs = preflight.get("inputs")
    if not isinstance(preflight_inputs, Mapping):
        raise ContractError("Cooccurrence preflight lacks input identities")
    for name, path in (
        ("manifest", bundle.manifest),
        ("assignments", bundle.screen_assignments),
        ("sites", bundle.screen_sites),
        ("primary_artifact_audit_pre", bundle.primary_artifact_audit_pre),
    ):
        verify_declared_artifact(
            preflight_inputs.get(name), path, f"cooccurrence preflight {name}"
        )
    if bundle.independent_m2_audit is not None:
        verify_declared_artifact(
            preflight_inputs.get("independent_m2_audit"),
            bundle.independent_m2_audit,
            "cooccurrence preflight independent M2 audit",
        )
    preflight_input_lock = preflight.get("input_lock")
    if (
        not isinstance(preflight_input_lock, Mapping)
        or preflight_input_lock.get("identity_before") != preflight_inputs
        or preflight_input_lock.get("identity_after") != preflight_inputs
        or preflight_input_lock.get("all_primary_inputs_unchanged") is not True
    ):
        raise ContractError("Cooccurrence preflight input lock drift")

    receipt_inputs = receipt.get("inputs")
    if not isinstance(receipt_inputs, Mapping):
        raise ContractError("Cooccurrence receipt lacks input identities")
    producer_locked_inputs = {
        name: value
        for name, value in receipt_inputs.items()
        if name != "intersubmod_runs"
    }
    producer_input_lock = receipt.get("input_lock")
    if (
        not isinstance(producer_input_lock, Mapping)
        or receipt.get("release_status")
        != "PRODUCER_PASS_PENDING_RUNNER_RECONCILIATION"
        or producer_input_lock.get("identity_before") != producer_locked_inputs
        or producer_input_lock.get("identity_after_compute") != producer_locked_inputs
        or producer_input_lock.get("identity_after_output") != producer_locked_inputs
        or producer_input_lock.get("all_primary_inputs_unchanged") is not True
    ):
        raise ContractError("Cooccurrence producer input lock/release status drift")
    verify_declared_artifact(
        receipt_inputs.get("preflight_receipt"),
        bundle.cooccurrence_preflight,
        "cooccurrence preflight receipt binding",
    )
    gate = receipt.get("preflight_gate")
    if (
        not isinstance(gate, Mapping)
        or gate.get("schema_name")
        != "intersubmod.cooccurrence_task_contract_preflight"
        or gate.get("schema_version") != COOCCURRENCE_PREFLIGHT_SCHEMA_VERSION
        or gate.get("pass") is not True
        or int(gate.get("task_count", -1)) != 102_842
        or gate.get("raw_identity_release_contract")
        != RAW_IDENTITY_RELEASE_CONTRACT
    ):
        raise ContractError("Cooccurrence receipt preflight gate binding drift")

    code = receipt.get("code")
    source_lock = receipt.get("source_lock")
    if not isinstance(code, Mapping) or set(code) != set(COOCCURRENCE_CODE_PATHS):
        raise ContractError("Cooccurrence receipt code role set is not exact")
    for role, path in COOCCURRENCE_CODE_PATHS.items():
        verify_declared_artifact(code.get(role), path, f"cooccurrence code {role}")
    if not isinstance(source_lock, Mapping):
        raise ContractError("Cooccurrence receipt lacks source-lock evidence")
    for identity_field in (
        "source_identity_before",
        "source_identity_after_compute",
        "source_identity_after_output",
    ):
        identity = source_lock.get(identity_field)
        if not isinstance(identity, Mapping) or dict(identity) != dict(code):
            raise ContractError(f"Cooccurrence source-lock identity drift: {identity_field}")
    for mode_field in (
        "source_modes_before",
        "source_modes_after_compute",
        "source_modes_after_output",
    ):
        modes = source_lock.get(mode_field)
        if (
            not isinstance(modes, Mapping)
            or set(modes) != set(COOCCURRENCE_CODE_PATHS)
            or set(modes.values()) != {"0o444"}
        ):
            raise ContractError(f"Cooccurrence source-lock mode drift: {mode_field}")
    if source_lock.get("all_sources_read_only_and_unchanged") is not True:
        raise ContractError("Cooccurrence source-lock pass marker is absent")

    preflight_code = preflight.get("code")
    if not isinstance(preflight_code, Mapping):
        raise ContractError("Cooccurrence preflight lacks source lock")
    pre_before = preflight_code.get("source_identity_before")
    pre_after = preflight_code.get("source_identity_after")
    expected_preflight_roles = {"preflight", *COOCCURRENCE_CODE_PATHS}
    if (
        not isinstance(pre_before, Mapping)
        or not isinstance(pre_after, Mapping)
        or set(pre_before) != expected_preflight_roles
        or dict(pre_before) != dict(pre_after)
    ):
        raise ContractError("Cooccurrence preflight source identity map drift")
    verify_declared_artifact(
        pre_before.get("preflight"),
        COOCCURRENCE_PREFLIGHT_SOURCE,
        "cooccurrence preflight source",
    )
    for role in COOCCURRENCE_CODE_PATHS:
        if pre_before.get(role) != code.get(role):
            raise ContractError(f"Preflight/analyzer source identity differs for {role}")
    for mode_field in ("source_modes_before", "source_modes_after"):
        modes = preflight_code.get(mode_field)
        if (
            not isinstance(modes, Mapping)
            or set(modes) != expected_preflight_roles
            or set(modes.values()) != {"0o444"}
        ):
            raise ContractError(f"Cooccurrence preflight source mode drift: {mode_field}")

    raw_summary = summary.get("raw_bam_identity_recovery_audit")
    raw_preflight = (preflight.get("observed") or {}).get(
        "raw_bam_identity_recovery"
    )
    counts = receipt.get("counts")
    if not all(isinstance(value, Mapping) for value in (raw_summary, raw_preflight, counts)):
        raise ContractError("Cooccurrence raw-identity policy evidence is incomplete")
    for raw_payload, label in (
        (raw_summary, "summary"),
        (raw_preflight, "preflight"),
    ):
        if (
            raw_payload.get("equivalence_contract")
            != RAW_DUPLICATE_EQUIVALENCE_CONTRACT
            or raw_payload.get("allowed_differing_auxiliary_tags") != ["RG"]
            or raw_payload.get("missing_projection_policy")
            != RAW_IDENTITY_MISSING_POLICY
            or raw_payload.get("conflicting_analysis_payload_policy")
            != RAW_IDENTITY_CONFLICT_POLICY
            or raw_payload.get("failure_counts_materialized") is not False
        ):
            raise ContractError(f"Cooccurrence {label} fail-closed policy drift")
    if (
        counts.get("raw_identity_missing_projection_policy")
        != RAW_IDENTITY_MISSING_POLICY
        or counts.get("raw_identity_conflicting_analysis_payload_policy")
        != RAW_IDENTITY_CONFLICT_POLICY
        or counts.get("raw_identity_failure_counts_materialized") is not False
        or counts.get("raw_identity_all_site_results_passed_invariant_validation")
        is not True
    ):
        raise ContractError("Cooccurrence receipt fail-closed policy drift")


def validate_raw_identity_duplicate_artifact(
    bundle: InputBundle,
    cooccurrence_rows: Mapping[CandidateKey, Mapping[str, str]],
    summary: Mapping[str, Any],
    receipt: Mapping[str, Any],
    preflight: Mapping[str, Any],
) -> None:
    by_locus: dict[tuple[str, str, int], Mapping[str, str]] = {}
    for key, site_row in cooccurrence_rows.items():
        locus = (key[0], key[1], key[2])
        if locus in by_locus:
            raise ContractError(f"Multiple cooccurrence alleles share raw-audit locus {locus}")
        by_locus[locus] = site_row

    sparse_by_locus: dict[tuple[str, str, int], list[dict[str, Any]]] = defaultdict(list)
    seen_sparse_keys: set[tuple[Any, ...]] = set()
    total_extra_records = 0
    classifications: Counter[str] = Counter()
    with open_text(bundle.cooccurrence_raw_identity_duplicates) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != RAW_IDENTITY_DUPLICATE_COLUMNS:
            raise ContractError("Raw identity duplicate audit header drift")
        for line_number, row in enumerate(reader, 2):
            try:
                locus = (row["sample"], row["chrom"], int(row["pos"]))
                projection = (
                    row["projection_read_name"],
                    row["projection_chrom"],
                    int(row["projection_start"]),
                    int(row["projection_end"]),
                    int(row["projection_mapq"]),
                    row["projection_strand"],
                )
                record_count = int(row["record_count"])
                differing_tags = json.loads(row["differing_auxiliary_tags"])
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
                raise ContractError(
                    f"Malformed raw identity duplicate audit row {line_number}"
                ) from error
            if locus not in by_locus or projection[1] != locus[1] or record_count < 2:
                raise ContractError(f"Invalid raw identity duplicate audit row {line_number}")
            classification = row["classification"]
            if (
                classification == "EXACT_DUPLICATE"
                and differing_tags != []
            ) or (
                classification == "RG_ONLY_DUPLICATE"
                and differing_tags != ["RG"]
            ):
                raise ContractError(f"Raw identity duplicate class/tag drift at row {line_number}")
            if classification not in {"EXACT_DUPLICATE", "RG_ONLY_DUPLICATE"}:
                raise ContractError(f"Unexpected raw identity duplicate class at row {line_number}")
            sparse_key = (*locus, *projection)
            if sparse_key in seen_sparse_keys:
                raise ContractError(f"Duplicate sparse raw identity projection at row {line_number}")
            seen_sparse_keys.add(sparse_key)
            normalized = {
                "projection": projection,
                "record_count": record_count,
                "classification": classification,
            }
            sparse_by_locus[locus].append(normalized)
            total_extra_records += record_count - 1
            classifications[classification] += 1

    site_weighted_digest = hashlib.sha256()
    aggregate: Counter[str] = Counter()
    for locus, row in sorted(by_locus.items()):
        try:
            expected = int(row["raw_identity_expected_projections"])
            matched = int(row["raw_identity_matched_records"])
            duplicates = int(row["raw_identity_duplicate_projections_collapsed"])
            extras = int(row["raw_identity_duplicate_extra_records_collapsed"])
            exact = int(row["raw_identity_exact_duplicate_projections_collapsed"])
            rg_only = int(row["raw_identity_rg_only_duplicate_projections_collapsed"])
            digest = row["raw_identity_duplicate_projection_sha256"]
            alignment_digest = row["raw_identity_alignment_payload_sha256"]
            recovered_digest = row["raw_identity_recovered_payload_sha256"]
            analysis_site_digest = row["raw_identity_analysis_site_payload_sha256"]
            focal_reads = int(row["n_all_focal_ref_alt_reads"])
        except (KeyError, TypeError, ValueError) as error:
            raise ContractError(f"Malformed raw identity site fields at {locus}") from error
        if (
            min(expected, matched, duplicates, extras, exact, rg_only) < 0
            or expected < 1
            or expected != focal_reads
            or matched != expected + extras
            or duplicates > expected
            or extras < duplicates
            or exact + rg_only != duplicates
            or (duplicates == 0) != (extras == 0)
            or re.fullmatch(r"[0-9a-f]{64}", digest or "") is None
            or re.fullmatch(r"[0-9a-f]{64}", alignment_digest or "") is None
            or re.fullmatch(r"[0-9a-f]{64}", recovered_digest or "") is None
            or re.fullmatch(r"[0-9a-f]{64}", analysis_site_digest or "") is None
        ):
            raise ContractError(f"Raw identity site invariant failed at {locus}")
        sparse_rows = sorted(
            sparse_by_locus.get(locus, []), key=lambda value: value["projection"]
        )
        if (
            len(sparse_rows) != duplicates
            or sum(value["record_count"] - 1 for value in sparse_rows) != extras
            or sum(value["classification"] == "EXACT_DUPLICATE" for value in sparse_rows)
            != exact
            or sum(value["classification"] == "RG_ONLY_DUPLICATE" for value in sparse_rows)
            != rg_only
            or projection_digest(value["projection"] for value in sparse_rows) != digest
        ):
            raise ContractError(f"Raw identity sparse/site reconciliation failed at {locus}")
        payload = [
            *locus,
            expected,
            matched,
            duplicates,
            extras,
            exact,
            rg_only,
            digest,
            alignment_digest,
            recovered_digest,
            analysis_site_digest,
        ]
        site_weighted_digest.update(compact_json(payload).encode("utf-8"))
        site_weighted_digest.update(b"\n")
        aggregate.update(
            {
                "site_tasks": 1,
                "expected_projection_occurrences": expected,
                "matched_record_occurrences": matched,
                "sites_with_collapsed_duplicates": int(duplicates > 0),
                "duplicate_projection_occurrences_collapsed": duplicates,
                "duplicate_extra_record_occurrences_collapsed": extras,
                "exact_duplicate_projection_occurrences_collapsed": exact,
                "rg_only_duplicate_projection_occurrences_collapsed": rg_only,
            }
        )

    raw_summary = summary["raw_bam_identity_recovery_audit"]
    raw_preflight = preflight["observed"]["raw_bam_identity_recovery"]
    if (
        raw_summary.get("analysis_scope_identity_contract")
        != ANALYSIS_SCOPE_IDENTITY_CONTRACT
        or raw_preflight.get("analysis_scope_identity_contract")
        != ANALYSIS_SCOPE_IDENTITY_CONTRACT
    ):
        raise ContractError("Cooccurrence analysis-scope identity contract drift")
    raw_preflight_aggregate = raw_preflight["aggregate"]
    receipt_counts = receipt["counts"]
    mapping = {
        "site_tasks": "n_site_tasks",
        "expected_projection_occurrences": "n_expected_projection_occurrences",
        "matched_record_occurrences": "n_matched_record_occurrences",
        "sites_with_collapsed_duplicates": "n_sites_with_collapsed_duplicates",
        "duplicate_projection_occurrences_collapsed": (
            "n_duplicate_projection_occurrences_collapsed"
        ),
        "duplicate_extra_record_occurrences_collapsed": (
            "n_duplicate_extra_record_occurrences_collapsed"
        ),
        "exact_duplicate_projection_occurrences_collapsed": (
            "n_exact_duplicate_projection_occurrences_collapsed"
        ),
        "rg_only_duplicate_projection_occurrences_collapsed": (
            "n_rg_only_duplicate_projection_occurrences_collapsed"
        ),
    }
    for aggregate_field, summary_field in mapping.items():
        observed = int(aggregate.get(aggregate_field, 0))
        if (
            int(raw_summary.get(summary_field, -1)) != observed
            or int(raw_preflight_aggregate.get(aggregate_field, -1)) != observed
        ):
            raise ContractError(f"Raw identity aggregate drift: {aggregate_field}")
    digest = site_weighted_digest.hexdigest()
    if (
        raw_summary.get("site_weighted_audit_sha256") != digest
        or raw_preflight.get("site_weighted_audit_sha256") != digest
        or int(receipt_counts.get("raw_identity_sparse_duplicate_rows", -1))
        != len(seen_sparse_keys)
        or len(seen_sparse_keys)
        != int(aggregate["duplicate_projection_occurrences_collapsed"])
        or total_extra_records
        != int(aggregate["duplicate_extra_record_occurrences_collapsed"])
        or classifications["EXACT_DUPLICATE"]
        != int(aggregate["exact_duplicate_projection_occurrences_collapsed"])
        or classifications["RG_ONLY_DUPLICATE"]
        != int(aggregate["rg_only_duplicate_projection_occurrences_collapsed"])
    ):
        raise ContractError("Raw identity sparse artifact aggregate/digest drift")


def validate_cooccurrence_oracle_artifact(bundle: InputBundle) -> None:
    payload = load_json(bundle.cooccurrence_oracle_cases, "cooccurrence oracle cases")
    if (
        payload.get("schema_name")
        != "intersubmod.methyl_ssnv_cooccurrence_validation.oracle_cases"
        or payload.get("schema_version") != COOCCURRENCE_SCHEMA_VERSION
        or int(payload.get("partner_window_bp", -1))
        != COOCCURRENCE_PARTNER_WINDOW_BP
    ):
        raise ContractError("Cooccurrence oracle schema or window drift")
    focal_cases = payload.get("focal_cases")
    if not isinstance(focal_cases, list) or len(focal_cases) != len(
        COOCCURRENCE_ORACLE_FOCAL_IDS
    ):
        raise ContractError("Cooccurrence focal oracle inventory drift")
    observed_ids: set[str] = set()
    for case in focal_cases:
        if not isinstance(case, Mapping):
            raise ContractError("Cooccurrence focal oracle case is malformed")
        case_id = str(case.get("case_id", ""))
        expected_positions = case.get("expected_partner_positions")
        observed_positions = case.get("observed_partner_positions")
        if (
            case_id in observed_ids
            or case.get("case_type") != "focal_partner_window"
            or case.get("present") is not True
            or case.get("partner_window_oracle_pass") is not True
            or not isinstance(expected_positions, list)
            or observed_positions != expected_positions
        ):
            raise ContractError(f"Cooccurrence focal oracle failed: {case_id!r}")
        observed_ids.add(case_id)
    if observed_ids != COOCCURRENCE_ORACLE_FOCAL_IDS:
        raise ContractError("Cooccurrence focal oracle case IDs drift")
    shared = payload.get("shared_readset_case")
    if (
        not isinstance(shared, Mapping)
        or shared.get("case_id") != COOCCURRENCE_SHARED_READSET_ORACLE_ID
        or shared.get("case_type") != "shared_alt_readset_dedup"
        or shared.get("same_alt_readset") is not True
        or shared.get("one_alt_readset_representative") is not True
        or shared.get("present_positions") != shared.get("positions")
    ):
        raise ContractError("Cooccurrence shared-readset oracle failed")


def validate_cooccurrence_release_receipt(
    bundle: InputBundle,
    producer_receipt: Mapping[str, Any],
    *,
    expected_stable_sites: int,
) -> dict[str, Any]:
    release = load_json(
        bundle.cooccurrence_release_receipt,
        "cooccurrence release receipt",
    )
    require_schema_pass(
        release,
        "intersubmod.cooccurrence_release_receipt",
        "cooccurrence release receipt",
        schema_version=COOCCURRENCE_RELEASE_SCHEMA_VERSION,
    )
    checks = release.get("checks")
    if (
        release.get("task_type") != "B_comprehensive_validation"
        or release.get("scope") != "all_manifest_samples"
        or release.get("raw_identity_release_contract")
        != RAW_IDENTITY_RELEASE_CONTRACT
        or release.get("release_status") != "RELEASE_RECONCILIATION_PASS"
        or release.get("producer_status") != "EXECUTION_PASS"
        or release.get("pass_semantics")
        != "runner_reconciled_release_integrity_only_not_scientific_claim"
        or not isinstance(checks, Mapping)
        or set(checks) != RELEASE_FINALIZER.RELEASE_CHECK_KEYS
        or not all(value is True for value in checks.values())
    ):
        raise ContractError("Cooccurrence release reconciliation did not pass")
    if release.get("source_authority") != current_release_source_authority():
        raise ContractError("Cooccurrence release source authority drift")
    if oct(bundle.cooccurrence_release_receipt.stat().st_mode & 0o777) != "0o444":
        raise ContractError("Cooccurrence release receipt is not mode 0444")
    inputs = release.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ContractError("Cooccurrence release receipt lacks input identities")
    for role, path in (
        ("preflight", bundle.cooccurrence_preflight),
        ("producer_receipt", bundle.cooccurrence_receipt),
        ("summary", bundle.cooccurrence_summary),
        ("sites", bundle.cooccurrence_sites),
        ("pairs", bundle.cooccurrence_pairs),
        ("duplicates", bundle.cooccurrence_raw_identity_duplicates),
        ("oracle", bundle.cooccurrence_oracle_cases),
    ):
        verify_declared_artifact(inputs.get(role), path, f"cooccurrence release {role}")
    expected_command = RELEASE_FINALIZER.canonical_command(
        argparse.Namespace(
            preflight=bundle.cooccurrence_preflight,
            producer_receipt=bundle.cooccurrence_receipt,
            summary=bundle.cooccurrence_summary,
            sites=bundle.cooccurrence_sites,
            pairs=bundle.cooccurrence_pairs,
            duplicates=bundle.cooccurrence_raw_identity_duplicates,
            oracle=bundle.cooccurrence_oracle_cases,
            runner_source=COOCCURRENCE_RELEASE_CODE_PATHS["release_runner"],
            output=bundle.cooccurrence_release_receipt,
        )
    )
    if release.get("command") != expected_command:
        raise ContractError("Cooccurrence release finalizer command is not canonical")
    code = release.get("code")
    if not isinstance(code, Mapping) or set(code) != set(COOCCURRENCE_RELEASE_CODE_PATHS):
        raise ContractError("Cooccurrence release code role set is not exact")
    for role, path in COOCCURRENCE_RELEASE_CODE_PATHS.items():
        verify_declared_artifact(code.get(role), path, f"cooccurrence release code {role}")
    if release.get("source_modes") != {
        role: "0o444" for role in COOCCURRENCE_RELEASE_CODE_PATHS
    }:
        raise ContractError("Cooccurrence release source-mode contract drift")
    recomputed = RELEASE_FINALIZER.recompute_output_contract(
        sites_path=bundle.cooccurrence_sites,
        pairs_path=bundle.cooccurrence_pairs,
        duplicates_path=bundle.cooccurrence_raw_identity_duplicates,
        oracle_path=bundle.cooccurrence_oracle_cases,
        expected_stable_sites=expected_stable_sites,
    )
    if release.get("recomputed") != recomputed:
        raise ContractError(
            "Cooccurrence release independently recomputed payload drift"
        )
    if producer_receipt.get("release_status") != (
        "PRODUCER_PASS_PENDING_RUNNER_RECONCILIATION"
    ):
        raise ContractError("Cooccurrence producer status bypassed release reconciliation")
    try:
        producer_finished = datetime.fromisoformat(
            str(producer_receipt["finished_at_utc"])
        )
        release_started = datetime.fromisoformat(str(release["started_at_utc"]))
        release_finished = datetime.fromisoformat(str(release["finished_at_utc"]))
    except (KeyError, TypeError, ValueError) as error:
        raise ContractError("Cooccurrence release chronology is malformed") from error
    if not producer_finished <= release_started <= release_finished:
        raise ContractError("Cooccurrence release chronology drifted")
    return release


def validate_cooccurrence_receipts(
    bundle: InputBundle,
    cooccurrence_rows: Mapping[CandidateKey, Mapping[str, str]],
    pair_data: PairData,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    summary = load_json(bundle.cooccurrence_summary, "cooccurrence summary")
    receipt = load_json(bundle.cooccurrence_receipt, "cooccurrence run receipt")
    preflight = load_json(
        bundle.cooccurrence_preflight, "cooccurrence raw-identity preflight"
    )
    require_schema_pass(
        summary,
        "intersubmod.methyl_ssnv_cooccurrence_validation.summary",
        "cooccurrence summary",
        schema_version=COOCCURRENCE_SCHEMA_VERSION,
    )
    require_schema_pass(
        receipt,
        "intersubmod.methyl_ssnv_cooccurrence_validation.run_receipt",
        "cooccurrence run receipt",
        schema_version=COOCCURRENCE_SCHEMA_VERSION,
    )
    if summary.get("task_type") != "B_comprehensive_validation":
        raise ContractError("Cooccurrence summary is not Task Type B")
    if tuple(summary.get("samples") or []) != DATASETS:
        raise ContractError("Cooccurrence summary does not cover canonical 7/7 datasets")
    validate_cooccurrence_input_identity_and_m2_provenance(summary, receipt)
    expected_counts = {
        "n_stable_sites": len(cooccurrence_rows),
        "n_focal_partner_pairs": pair_data.n_rows,
        "n_endpoint_a_testable_pairs": pair_data.n_testable,
        "n_pair_by_confirmed": pair_data.n_confirmed,
        "n_multi_marker_molecular_haplotype_base_candidates": sum(
            parse_bool(
                row[FORMAL_SELECTION_COLUMN], field_name=FORMAL_SELECTION_COLUMN
            )
            for row in cooccurrence_rows.values()
        ),
    }
    for field_name, expected in expected_counts.items():
        if int(summary.get(field_name, -1)) != expected:
            raise ContractError(f"Cooccurrence summary count mismatch for {field_name}")
    receipt_counts = receipt.get("counts") or {}
    if int(receipt_counts.get("stable_sites", -1)) != len(cooccurrence_rows):
        raise ContractError("Cooccurrence receipt stable count mismatch")
    if int(receipt_counts.get("focal_partner_pairs", -1)) != pair_data.n_rows:
        raise ContractError("Cooccurrence receipt pair count mismatch")
    if int(receipt_counts.get("pair_by_confirmed", -1)) != pair_data.n_confirmed:
        raise ContractError("Cooccurrence receipt formal G1 pair count mismatch")
    replication = summary.get("cross_platform_replication") or {}
    if int(replication.get("n_exact_pairs_present_both", -1)) != len(
        pair_data.exact_cross_platform_pairs
    ):
        raise ContractError("Cooccurrence cross-platform exact-pair count mismatch")
    if int(replication.get("n_exact_pairs_confirmed_both", -1)) != len(
        pair_data.both_confirmed_cross_platform_pairs
    ):
        raise ContractError("Cooccurrence cross-platform confirmed-pair count mismatch")
    if int(
        replication.get(
            "n_exact_pairs_conditional_by_effect_direction_replicated", -1
        )
    ) != len(pair_data.replicated_cross_platform_pairs):
        raise ContractError("Cooccurrence technical-replication count mismatch")
    reconciliation = summary.get("site_pair_count_reconciliation")
    if not isinstance(reconciliation, Mapping):
        raise ContractError("Cooccurrence summary lacks site-pair reconciliation")
    if int(reconciliation.get("n_sites_reconciled", -1)) != len(
        cooccurrence_rows
    ):
        raise ContractError("Cooccurrence site-pair reconciliation site count mismatch")
    if int(reconciliation.get("n_pair_rows_reconciled", -1)) != pair_data.n_rows:
        raise ContractError("Cooccurrence site-pair reconciliation pair count mismatch")
    outputs = receipt.get("outputs") or {}
    verify_declared_artifact(outputs.get("site_tsv"), bundle.cooccurrence_sites, "cooccurrence sites")
    verify_declared_artifact(outputs.get("pair_tsv"), bundle.cooccurrence_pairs, "cooccurrence pairs")
    verify_declared_artifact(
        outputs.get("summary_json"), bundle.cooccurrence_summary, "cooccurrence summary"
    )
    verify_declared_artifact(
        outputs.get("raw_identity_duplicate_audit_tsv"),
        bundle.cooccurrence_raw_identity_duplicates,
        "cooccurrence raw identity duplicate audit",
    )
    verify_declared_artifact(
        outputs.get("case_json"),
        bundle.cooccurrence_oracle_cases,
        "cooccurrence oracle cases",
    )
    validate_cooccurrence_oracle_artifact(bundle)
    inputs = receipt.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ContractError("Cooccurrence receipt lacks SHA-locked inputs")
    verify_declared_artifact(inputs.get("manifest"), bundle.manifest, "cooccurrence manifest")
    verify_declared_artifact(
        inputs.get("assignments"), bundle.screen_assignments, "cooccurrence assignments"
    )
    verify_declared_artifact(inputs.get("sites"), bundle.screen_sites, "cooccurrence screen sites")
    verify_declared_artifact(
        inputs.get("primary_artifact_audit_pre"),
        bundle.primary_artifact_audit_pre,
        "cooccurrence pre-downstream primary artifact audit",
    )
    if bundle.independent_m2_audit is not None:
        verify_declared_artifact(
            inputs.get("independent_m2_audit"),
            bundle.independent_m2_audit,
            "cooccurrence independent M2 audit",
        )
    code = receipt.get("code")
    if not isinstance(code, Mapping):
        raise ContractError("Cooccurrence receipt lacks SHA-locked code artifacts")
    verify_declared_artifact(
        code.get("analyzer"),
        SCRIPT_DIR / "analyze_methyl_ssnv_cooccurrence.py",
        "cooccurrence analyzer code",
    )
    verify_declared_artifact(
        code.get("m2_screen_gate"),
        Path(M2_GATE.__file__).resolve(),
        "M2 screen gate code",
    )
    validate_cooccurrence_source_and_preflight_contract(
        bundle, summary, receipt, preflight
    )
    validate_raw_identity_duplicate_artifact(
        bundle, cooccurrence_rows, summary, receipt, preflight
    )
    validate_cooccurrence_release_receipt(
        bundle,
        receipt,
        expected_stable_sites=len(cooccurrence_rows),
    )
    return summary, receipt, preflight


def _nested_count(payload: Mapping[str, Any], field_name: str) -> int | None:
    candidates: list[Any] = [payload.get(field_name)]
    for container_name in ("counts", "selection_contract"):
        container = payload.get(container_name)
        if isinstance(container, Mapping):
            candidates.append(container.get(field_name))
    observed = [int(value) for value in candidates if value is not None]
    if not observed:
        return None
    if len(set(observed)) != 1:
        raise ContractError(f"Conflicting {field_name} values in N/A receipt: {observed}")
    return observed[0]


def find_explicit_na_receipt(
    directory: Path,
    *,
    prefix: str,
    schema_versions: set[str] | None = None,
    expected_selection_column: str = STRICT_SELECTION_COLUMN,
    require_execution_status: bool = False,
) -> tuple[Path, dict[str, Any]]:
    candidates = (
        directory / "run_receipt.json",
        directory / "run_manifest.json",
        directory / "not_applicable_receipt.json",
    )
    existing = [path for path in candidates if path.is_file() and path.stat().st_size > 0]
    if len(existing) != 1:
        raise ContractError(
            f"Expected exactly one explicit N/A receipt under {directory}, found {existing}"
        )
    path = existing[0]
    payload = load_json(path, "not-applicable receipt")
    schema_name = str(payload.get("schema_name", ""))
    if not schema_name.startswith(prefix):
        raise ContractError(f"Unexpected N/A receipt schema: {schema_name!r}")
    accepted_versions = schema_versions or {SCHEMA_VERSION}
    if payload.get("schema_version") not in accepted_versions or payload.get("pass") is not True:
        raise ContractError(
            f"N/A receipt must have schema_version in {sorted(accepted_versions)} and pass=true"
        )
    if payload.get("status") != "NOT_APPLICABLE":
        raise ContractError("N/A receipt status must be NOT_APPLICABLE")
    if require_execution_status and payload.get("execution_status") != "NOT_APPLICABLE":
        raise ContractError("Strict N/A receipt execution_status must be NOT_APPLICABLE")
    reason = str(
        payload.get("reason")
        or payload.get("not_applicable_reason")
        or payload.get("status_reason")
        or ""
    ).upper()
    if reason not in {"ZERO_SELECTED_CANDIDATES", "NO_SELECTED_CANDIDATES"}:
        raise ContractError(
            "N/A receipt reason must be ZERO_SELECTED_CANDIDATES or NO_SELECTED_CANDIDATES"
        )
    selected_count = _nested_count(payload, "n_selected_candidates")
    if selected_count != 0:
        raise ContractError("N/A receipt must explicitly declare n_selected_candidates=0")
    observed_selection_column = payload.get("selection_column")
    contract = payload.get("selection_contract")
    if observed_selection_column is None and isinstance(contract, Mapping):
        observed_selection_column = contract.get("selection_column")
    if observed_selection_column is None:
        raise ContractError("N/A receipt must explicitly declare its selection column")
    if observed_selection_column != expected_selection_column:
        raise ContractError(
            f"N/A receipt selection-column drift: {observed_selection_column!r} "
            f"!= {expected_selection_column!r}"
        )
    return path, payload


def validate_matched_normal_paired_na_receipt(
    path: Path,
    receipt: Mapping[str, Any],
    bundle: InputBundle,
) -> list[Path]:
    require_schema_pass(
        receipt,
        "intersubmod.matched_normal_candidate_control.not_applicable",
        "matched-normal paired-runner N/A receipt",
        schema_version=MATCHED_NORMAL_PAIRED_RUNNER_SCHEMA_VERSION,
    )
    if receipt.get("status") != "NOT_APPLICABLE" or receipt.get(
        "execution_status"
    ) != "NOT_APPLICABLE" or receipt.get("reason") != "ZERO_SELECTED_CANDIDATES":
        raise ContractError("Matched-normal N/A status/reason drift")
    if receipt.get("task_type") != "B_comprehensive_validation":
        raise ContractError("Matched-normal N/A receipt is not Task Type B")
    if receipt.get("selection_column") != STRICT_SELECTION_COLUMN or _nested_count(
        receipt, "n_selected_candidates"
    ) != 0:
        raise ContractError("Matched-normal N/A selection contract drift")
    if receipt.get("pass_semantics") != (
        "receipt_integrity_only_not_normal_control_execution_or_negative_evidence"
    ):
        raise ContractError("Matched-normal N/A pass semantics drift")
    producer_path = validate_release_bound_producer(
        receipt,
        receipt_path=path,
        producer_module=MATCHED_NORMAL_RUNNER,
        producer_role="matched_normal_runner",
        expected_command=MATCHED_NORMAL_RUNNER.canonical_task_b_command(),
        label="matched-normal N/A runner",
    )
    for field_name in (
        "sample_outputs_created",
        "cpp_executed",
        "normal_control_executed",
        "not_evaluable_is_negative",
    ):
        if receipt.get(field_name) is not False:
            raise ContractError(f"Matched-normal N/A contract drift for {field_name}")
    inputs = receipt.get("inputs")
    expected_provenance = {
        "candidate_table",
        "all_ssnv_manifest",
        "normal_audit",
        "binary",
        "reference",
        "runner_script",
    }
    if not isinstance(inputs, Mapping) or set(inputs) != expected_provenance:
        raise ContractError("Matched-normal N/A six-item provenance set drift")
    verify_declared_artifact(
        inputs.get("candidate_table"),
        bundle.cooccurrence_sites,
        "matched-normal N/A candidate table",
    )
    verify_declared_artifact(
        inputs.get("all_ssnv_manifest"),
        bundle.manifest,
        "matched-normal N/A all-sSNV manifest",
    )
    paths = [path, producer_path]
    for field_name in ("normal_audit", "binary", "reference", "runner_script"):
        paths.append(
            verify_declared_path_artifact(
                inputs.get(field_name), f"matched-normal N/A {field_name}"
            )
        )
    outputs = receipt.get("outputs")
    if not isinstance(outputs, Mapping) or Path(
        str(outputs.get("output_root", ""))
    ).resolve() != path.parent.resolve() or Path(
        str(outputs.get("not_applicable_receipt", ""))
    ).resolve() != path.resolve() or outputs.get("sample_outputs") != []:
        raise ContractError("Matched-normal N/A output contract drift")
    return paths


def require_strict_execution_pass(payload: Mapping[str, Any], label: str) -> None:
    execution_status = payload.get("execution_status")
    if execution_status not in {"PASS", "EXECUTION_PASS"}:
        raise ContractError(
            f"{label} must explicitly declare execution_status=PASS/EXECUTION_PASS; "
            f"observed={execution_status!r}"
        )


def validate_strict_execution_authority(
    receipt: Mapping[str, Any],
    *,
    receipt_path: Path,
    bundle: InputBundle,
) -> None:
    source_authority = current_release_source_authority()
    if receipt.get("source_authority") != source_authority:
        raise ContractError("Strict receipt source authority drift")
    expected_code = STRICT_PRODUCER.capture_source_identity()
    expected_modes = STRICT_PRODUCER.capture_source_modes()
    if receipt.get("code") != expected_code:
        raise ContractError("Strict receipt producer source identity drift")
    current_runtime = STRICT_PRODUCER.capture_runtime_environment_identity()
    formal_cache_required = (
        bundle.strict_dir.resolve()
        == CANONICAL_TASK_B_PATHS["strict_dir"].resolve()
    )
    try:
        STRICT_PRODUCER.require_canonical_runtime(
            current_runtime,
            output_dir=bundle.strict_dir,
            formal_cache_required=formal_cache_required,
        )
    except STRICT_PRODUCER.ConfirmationContractError as error:
        raise ContractError("Final builder is not running in the canonical strict runtime") from error
    runtime_lock = receipt.get("runtime_environment_lock")
    if (
        not isinstance(runtime_lock, Mapping)
        or runtime_lock.get("identity_before") != current_runtime
        or runtime_lock.get("identity_after_compute") != current_runtime
        or runtime_lock.get("isolated_runtime_unchanged") is not True
    ):
        raise ContractError("Strict receipt runtime environment lock drift")
    source_lock = receipt.get("source_lock")
    if (
        not isinstance(source_lock, Mapping)
        or source_lock.get("source_identity_before") != expected_code
        or source_lock.get("source_identity_after_compute") != expected_code
        or source_lock.get("source_modes_before") != expected_modes
        or source_lock.get("source_modes_after_compute") != expected_modes
        or source_lock.get("all_sources_read_only_and_unchanged") is not True
        or set(expected_modes.values()) != {"0o444"}
    ):
        raise ContractError("Strict receipt producer source lock drift")
    config = STRICT_PRODUCER.ConfirmationConfig()
    expected_command = STRICT_PRODUCER.canonical_command(
        candidate_table=bundle.cooccurrence_sites,
        assignments=bundle.screen_assignments,
        cooccurrence_receipt=bundle.cooccurrence_receipt,
        cooccurrence_release_receipt=bundle.cooccurrence_release_receipt,
        output_dir=bundle.strict_dir,
        config=config,
    )
    if receipt.get("command") != expected_command:
        raise ContractError("Strict receipt command is not canonical")
    parameters = receipt.get("parameters")
    expected_parameters = {
        "permutations_per_seed_per_null": STRICT_PRODUCER.DEFAULT_PERMUTATIONS,
        "seeds": STRICT_PRODUCER.DEFAULT_SEEDS,
        "formal_minimum_permutations": STRICT_PRODUCER.FORMAL_MIN_PERMUTATIONS,
        "formal_minimum_seeds": STRICT_PRODUCER.FORMAL_MIN_SEEDS,
        "formal_parameter_gate": True,
    }
    if not isinstance(parameters, Mapping) or any(
        parameters.get(field) != value
        for field, value in expected_parameters.items()
    ):
        raise ContractError("Strict receipt formal parameters drift")
    if oct(receipt_path.stat().st_mode & 0o777) != "0o444":
        raise ContractError("Strict receipt is not mode 0444")


def load_strict(
    bundle: InputBundle, candidate_keys: set[CandidateKey]
) -> tuple[dict[CandidateKey, dict[str, str]], dict[str, Any], list[Path]]:
    if not candidate_keys:
        path, receipt = find_explicit_na_receipt(
            bundle.strict_dir,
            prefix="intersubmod.strict_methyl_candidate_confirmation",
            schema_versions=STRICT_SCHEMA_VERSIONS,
            expected_selection_column=STRICT_SELECTION_COLUMN,
            require_execution_status=True,
        )
        if bundle.strict_sites.exists() or bundle.strict_summary.exists():
            raise ContractError("Zero-candidate strict N/A run must not include result table/summary")
        inputs = receipt.get("inputs")
        if not isinstance(inputs, Mapping):
            raise ContractError("Strict N/A receipt lacks SHA-locked inputs")
        if receipt.get("pass_semantics") != STRICT_PASS_SEMANTICS:
            raise ContractError("Strict N/A pass semantics drift")
        validate_strict_execution_authority(
            receipt,
            receipt_path=path,
            bundle=bundle,
        )
        interpretation = receipt.get("scientific_interpretation")
        if receipt.get("is_negative_result") is not False or not isinstance(
            interpretation, Mapping
        ) or interpretation.get("is_negative_result") is not False:
            raise ContractError(
                "Strict N/A receipt must explicitly declare is_negative_result=false"
            )
        verify_declared_artifact(
            inputs.get("candidate_table"), bundle.cooccurrence_sites, "strict N/A candidate table"
        )
        verify_declared_artifact(
            inputs.get("assignments"), bundle.screen_assignments, "strict N/A assignments"
        )
        verify_declared_artifact(
            inputs.get("cooccurrence_receipt"),
            bundle.cooccurrence_receipt,
            "strict N/A cooccurrence receipt",
        )
        verify_declared_artifact(
            inputs.get("cooccurrence_release_receipt"),
            bundle.cooccurrence_release_receipt,
            "strict N/A cooccurrence release receipt",
        )
        return (
            {},
            {
                "status": "NOT_APPLICABLE_VALIDATED_RECEIPT",
                "execution_status": "NOT_APPLICABLE",
                "is_negative_result": False,
                "receipt": receipt,
                "receipt_path": path,
            },
            [path],
        )
    rows = load_keyed_tsv(
        bundle.strict_sites,
        label="strict confirmation table",
        required_fields={
            "strict_analysis_class",
            "strict_formal_parameter_gate",
            "strict_formal_selection_contract_gate",
            "strict_formal_selection_column",
            "strict_selection_column",
            "strict_candidate_selection_gate",
            "strict_cooccurrence_receipt_gate",
            "strict_artifact_identity_gate",
            "strict_screen_artifact_hash_lock_gate",
            "strict_postselection_fdr_calibrated",
            "strict_combined_empirical_p_postselection_descriptive",
            "strict_postselection_bh_q_descriptive",
            "strict_postselection_by_q_descriptive",
            "strict_null_robustness_pass",
            "strict_methyl_partition_robustness_status",
            "strict_scientific_status",
            "strict_not_evaluable_reason",
            "strict_failure_reason",
        },
    )
    if set(rows) != candidate_keys:
        raise ContractError(
            "Strict/cooccurrence selected-candidate key-set mismatch: "
            f"selected={len(candidate_keys)} strict={len(rows)} "
            f"missing={sorted(candidate_keys - set(rows))[:3]} "
            f"extra={sorted(set(rows) - candidate_keys)[:3]}"
    )
    for key, row in rows.items():
        if (
            row["strict_selection_column"] != STRICT_SELECTION_COLUMN
            or row["strict_formal_selection_column"] != STRICT_SELECTION_COLUMN
        ):
            raise ContractError(f"Strict selection-column drift at {key}")
        if not parse_bool(
            row["strict_candidate_selection_gate"],
            field_name="strict_candidate_selection_gate",
        ):
            raise ContractError(f"Strict candidate selection gate is false at {key}")
        if parse_bool(
            row["strict_postselection_fdr_calibrated"],
            field_name="strict_postselection_fdr_calibrated",
        ):
            raise ContractError(f"Strict post-selection diagnostics claim calibrated FDR at {key}")
        for field_name in (
            "strict_formal_parameter_gate",
            "strict_formal_selection_contract_gate",
            "strict_cooccurrence_receipt_gate",
            "strict_artifact_identity_gate",
            "strict_screen_artifact_hash_lock_gate",
        ):
            if not parse_bool(row[field_name], field_name=field_name):
                raise ContractError(f"Strict provenance/formal gate is false at {key}: {field_name}")
        robustness_pass = parse_bool(
            row["strict_null_robustness_pass"],
            field_name="strict_null_robustness_pass",
        )
        analysis_class = row["strict_analysis_class"]
        robustness_status = row["strict_methyl_partition_robustness_status"]
        if analysis_class == "EXPLORATORY_ONLY":
            expected_status = "EXPLORATORY_ONLY"
            if robustness_pass:
                raise ContractError(f"Exploratory strict row cannot pass formal R1 at {key}")
        elif robustness_status == "NOT_EVALUABLE":
            expected_status = "NOT_EVALUABLE"
            if robustness_pass:
                raise ContractError(f"Not-evaluable strict row cannot pass R1 at {key}")
        else:
            expected_status = (
                "ROBUSTNESS_PASS_NOT_FDR_CALIBRATED"
                if robustness_pass
                else "ROBUSTNESS_FAIL"
            )
        if row["strict_scientific_status"] != expected_status:
            raise ContractError(
                f"Strict scientific-status/robustness mismatch at {key}: "
                f"{row['strict_scientific_status']!r} != {expected_status!r}"
            )
    summary = load_json(bundle.strict_summary, "strict confirmation summary")
    receipt = load_json(bundle.strict_receipt, "strict confirmation run manifest")
    summary_version = require_schema_pass_one_of(
        summary,
        "intersubmod.strict_methyl_candidate_confirmation_summary",
        "strict confirmation summary",
        schema_versions=STRICT_SCHEMA_VERSIONS,
    )
    receipt_version = require_schema_pass_one_of(
        receipt,
        "intersubmod.strict_methyl_candidate_confirmation_run_manifest",
        "strict confirmation run manifest",
        schema_versions=STRICT_SCHEMA_VERSIONS,
    )
    if summary_version != receipt_version:
        raise ContractError("Strict summary/receipt schema-version mismatch")
    require_strict_execution_pass(summary, "Strict confirmation summary")
    require_strict_execution_pass(receipt, "Strict confirmation run manifest")
    validate_strict_execution_authority(
        receipt,
        receipt_path=bundle.strict_receipt,
        bundle=bundle,
    )
    if summary.get("pass_semantics") != STRICT_PASS_SEMANTICS or receipt.get(
        "pass_semantics"
    ) != STRICT_PASS_SEMANTICS:
        raise ContractError("Strict summary/receipt pass semantics are not execution-only")
    diagnostic_contract = summary.get("postselection_diagnostic_contract")
    if not isinstance(diagnostic_contract, Mapping):
        raise ContractError("Strict summary lacks post-selection diagnostic contract")
    if diagnostic_contract.get("fdr_calibrated") is not False:
        raise ContractError("Strict post-selection diagnostics must declare fdr_calibrated=false")
    if diagnostic_contract.get("bh_by_values_are_descriptive_only") is not True:
        raise ContractError("Strict post-selection BH/BY must be explicitly descriptive only")
    selected_count = len(candidate_keys)
    robustness_count = sum(
        parse_bool(
            row["strict_null_robustness_pass"],
            field_name="strict_null_robustness_pass",
        )
        for row in rows.values()
    )
    robustness_evaluable_count = sum(
        row["strict_methyl_partition_robustness_status"]
        not in {"NOT_EVALUABLE", "EXPLORATORY_ONLY"}
        for row in rows.values()
    )
    for payload, label in ((summary, "summary"), (receipt, "receipt")):
        counts = payload.get("counts") or {}
        if int(counts.get("n_selected_candidates", -1)) != selected_count:
            raise ContractError(f"Strict {label} selected count mismatch")
        if int(counts.get("n_null_robustness_pass", -1)) != robustness_count:
            raise ContractError(f"Strict {label} robustness-pass count mismatch")
        if int(counts.get("n_null_robustness_fail_retained", -1)) != (
            robustness_evaluable_count - robustness_count
        ):
            raise ContractError(f"Strict {label} robustness-fail count mismatch")
    selection = summary.get("selection_contract") or {}
    if selection.get("selection_column") != STRICT_SELECTION_COLUMN:
        raise ContractError("Strict summary selection-column drift")
    if int(selection.get("n_selected_candidates", -1)) != selected_count:
        raise ContractError("Strict summary selection denominator mismatch")
    outputs = receipt.get("outputs") or {}
    verify_declared_artifact(outputs.get("site_results"), bundle.strict_sites, "strict sites")
    verify_declared_artifact(outputs.get("summary"), bundle.strict_summary, "strict summary")
    inputs = receipt.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ContractError("Strict receipt lacks SHA-locked inputs")
    verify_declared_artifact(
        inputs.get("candidate_table"), bundle.cooccurrence_sites, "strict candidate table"
    )
    verify_declared_artifact(
        inputs.get("assignments"), bundle.screen_assignments, "strict assignments"
    )
    verify_declared_artifact(
        inputs.get("cooccurrence_receipt"),
        bundle.cooccurrence_receipt,
        "strict cooccurrence receipt",
    )
    verify_declared_artifact(
        inputs.get("cooccurrence_release_receipt"),
        bundle.cooccurrence_release_receipt,
        "strict cooccurrence release receipt",
    )
    return (
        rows,
        {
            "status": "COMPLETE",
            "execution_status": summary["execution_status"],
            "schema_version": summary_version,
            "summary": summary,
            "receipt": receipt,
            "receipt_path": bundle.strict_receipt,
        },
        [bundle.strict_sites, bundle.strict_summary, bundle.strict_receipt],
    )


def replay_strict_statistics(
    bundle: InputBundle,
    observed_rows: Mapping[CandidateKey, Mapping[str, Any]],
    strict_status: Mapping[str, Any],
) -> dict[str, Any]:
    """Reexecute every strict statistic and compare all emitted strict fields."""
    config = STRICT_PRODUCER.ConfirmationConfig()
    input_fields = STRICT_PRODUCER.candidate_table_fields(bundle.cooccurrence_sites)
    selection = STRICT_PRODUCER.resolve_selection_column(input_fields)
    _, selected, n_input_rows = STRICT_PRODUCER.load_candidates(
        bundle.cooccurrence_sites,
        selection.resolved_column,
    )
    candidate_keys = STRICT_PRODUCER.candidate_table_keys(bundle.cooccurrence_sites)
    assignments = STRICT_PRODUCER.load_assignments(bundle.screen_assignments)
    selected_keys = [
        STRICT_PRODUCER.site_key(row, source="strict replay candidate")
        for row in selected
    ]
    cooccurrence_contract = STRICT_PRODUCER.validate_cooccurrence_receipt(
        bundle.cooccurrence_receipt,
        release_receipt_path=bundle.cooccurrence_release_receipt,
        candidate_table=bundle.cooccurrence_sites,
        assignments_path=bundle.screen_assignments,
        n_candidate_rows=n_input_rows,
        candidate_keys=candidate_keys,
        selected_keys=selected_keys,
        selection=selection,
    )
    missing_assignments = [key for key in selected_keys if key not in assignments]
    if missing_assignments:
        raise ContractError(
            "Strict deterministic replay is missing assignment records: "
            f"{missing_assignments[:8]}"
        )
    prepared = [
        STRICT_PRODUCER.prepare_site(row, assignments[key])
        for row, key in zip(selected, selected_keys, strict=True)
    ]
    results = [
        STRICT_PRODUCER.analyze_selected_site(site, config) for site in prepared
    ]
    STRICT_PRODUCER.apply_postselection_diagnostics(results)
    replay_rows = [
        STRICT_PRODUCER.flatten_result(
            result,
            selection,
            cooccurrence_contract,
            len(results),
            config,
        )
        for result in results
    ]
    replay_rows.sort(
        key=lambda row: (str(row["sample"]), str(row["chrom"]), int(row["pos"]))
    )
    receipt = strict_status.get("receipt")
    if not isinstance(receipt, Mapping):
        raise ContractError("Strict deterministic replay lacks receipt payload")
    declared = receipt.get("analysis_replay")
    observed = STRICT_PRODUCER.strict_analysis_replay_record(
        list(observed_rows.values())
    )
    recomputed = STRICT_PRODUCER.strict_analysis_replay_record(replay_rows)
    if declared != observed:
        raise ContractError("Strict receipt analysis digest does not match emitted rows")
    if recomputed != observed:
        raise ContractError(
            "Strict statistics deterministic replay does not match emitted rows"
        )
    return {
        "contract": "authorized_strict_statistics_full_deterministic_reexecution_v1",
        "implementation_independence": False,
        "execution_independence": True,
        "n_input_rows": n_input_rows,
        "n_selected_rows": len(replay_rows),
        "all_strict_output_fields_replayed": True,
        "analysis_replay": recomputed,
        "pass": True,
    }


def load_tumor_ref(
    bundle: InputBundle, stable_keys: set[CandidateKey]
) -> tuple[dict[CandidateKey, dict[str, str]], dict[str, Any], dict[str, Any]]:
    rows = load_keyed_tsv(
        bundle.tumor_ref_sites,
        label="tumor-REF control table",
        required_fields={"ref_evaluable", "ref_stable_null_multigroup"},
    )
    if set(rows) != stable_keys:
        raise ContractError(
            "Tumor-REF/stable key-set mismatch: "
            f"stable={len(stable_keys)} tumor_ref={len(rows)} "
            f"missing={sorted(stable_keys - set(rows))[:3]} extra={sorted(set(rows) - stable_keys)[:3]}"
        )
    for row in rows.values():
        parse_bool(row["ref_evaluable"], field_name="ref_evaluable")
        parse_bool(row["ref_stable_null_multigroup"], field_name="ref_stable_null_multigroup")
    summary = load_json(bundle.tumor_ref_summary, "tumor-REF summary")
    receipt = load_json(bundle.tumor_ref_receipt, "tumor-REF run manifest")
    require_schema_pass(
        summary,
        "intersubmod.all_ssnv_tumor_ref_controls.summary",
        "tumor-REF summary",
        schema_version=TUMOR_REF_SCHEMA_VERSION,
    )
    require_schema_pass(
        receipt,
        "intersubmod.all_ssnv_tumor_ref_controls.run_manifest",
        "tumor-REF run manifest",
        schema_version=TUMOR_REF_SCHEMA_VERSION,
    )
    if summary.get("task_type") != "B_comprehensive_validation":
        raise ContractError("Tumor-REF summary is not Task Type B")
    counts = receipt.get("counts") or {}
    for field_name in ("primary_stable_sites", "processed_sites"):
        if int(counts.get(field_name, -1)) != len(stable_keys):
            raise ContractError(f"Tumor-REF receipt count mismatch for {field_name}")
    outputs = receipt.get("outputs") or {}
    verify_declared_artifact(outputs.get("site_results"), bundle.tumor_ref_sites, "tumor-REF sites")
    verify_declared_artifact(outputs.get("summary"), bundle.tumor_ref_summary, "tumor-REF summary")
    inputs = receipt.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ContractError("Tumor-REF receipt lacks SHA-locked inputs")
    verify_declared_artifact(
        inputs.get("site_results"), bundle.screen_sites, "tumor-REF screen sites"
    )
    verify_declared_artifact(
        inputs.get("stable_assignments"),
        bundle.screen_assignments,
        "tumor-REF assignments",
    )
    verify_declared_artifact(
        inputs.get("primary_artifact_audit_pre"),
        bundle.primary_artifact_audit_pre,
        "tumor-REF pre-downstream primary artifact audit",
    )
    return rows, summary, receipt


def load_independent_m2_audit(
    path: Path | None,
    bundle: InputBundle,
    screen: ScreenData,
    cooccurrence_rows: Mapping[CandidateKey, Mapping[str, Any]],
) -> tuple[dict[str, Any], list[Path]]:
    if path is None:
        return {"status": "NOT_INCLUDED", "pass": None}, []
    payload = load_json(path, "logic-independent M2 audit")
    require_schema_pass(
        payload,
        "intersubmod.independent_m2_gate_recount",
        "logic-independent M2 audit",
        schema_version="2.0.0",
    )
    if payload.get("task_type") != "B_comprehensive_validation":
        raise ContractError("Logic-independent M2 audit is not Task Type B")
    if payload.get("contract") != M2_GATE.GATE_CONTRACT:
        raise ContractError("Logic-independent M2 gate contract drift")
    independence = payload.get("logic_independence")
    expected_independence = {
        "production_gate_imported": False,
        "production_gate_functions_called": False,
        "screen_effect_and_p_values_reused_as_frozen_inputs": True,
        "assignment_categories_and_coarse_group_counts_recomputed": True,
    }
    if independence != expected_independence:
        raise ContractError("Logic-independent M2 audit independence declaration drift")
    inputs = payload.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ContractError("Logic-independent M2 audit lacks input identities")
    verify_declared_artifact(
        inputs.get("site_results"), bundle.screen_sites, "independent M2 screen sites"
    )
    verify_declared_artifact(
        inputs.get("stable_assignments"),
        bundle.screen_assignments,
        "independent M2 stable assignments",
    )
    verify_declared_artifact(
        inputs.get("production_gate_source_reference_only"),
        Path(M2_GATE.__file__).resolve(),
        "independent M2 production-gate reference",
    )
    claim_contract = verify_declared_path_artifact(
        inputs.get("claim_contract"), "independent M2 claim contract"
    )
    if claim_contract.resolve() != CANONICAL_CLAIM_CONTRACT.resolve():
        raise ContractError("Logic-independent M2 claim-contract-v5 path is not canonical")
    if (
        claim_contract.stat().st_size != SOURCE_AUTHORITY.CLAIM_CONTRACT_SIZE
        or sha256(claim_contract) != SOURCE_AUTHORITY.CLAIM_CONTRACT_SHA256
        or oct(claim_contract.stat().st_mode & 0o777) != "0o444"
    ):
        raise ContractError("Logic-independent M2 claim-contract-v5 identity drift")
    code = payload.get("code")
    if not isinstance(code, Mapping):
        raise ContractError("Logic-independent M2 audit lacks code identity")
    audit_code = verify_declared_path_artifact(
        code.get("independent_recount"), "logic-independent M2 audit code"
    )
    if audit_code.resolve() != SOURCE_AUTHORITY.EXPECTED_SOURCE_PATHS[
        "independent_m2_auditor"
    ].resolve():
        raise ContractError("Logic-independent M2 audit source path is not canonical")
    checks = payload.get("checks")
    if not isinstance(checks, Mapping) or not checks or any(value is not True for value in checks.values()):
        raise ContractError("Logic-independent M2 audit checks are incomplete or nonpassing")
    counts = payload.get("counts")
    if not isinstance(counts, Mapping):
        raise ContractError("Logic-independent M2 audit lacks counts")
    expected_counts = {
        "all_rows": len(screen.all_keys),
        "m1_stable_rows": len(screen.stable_rows),
        "eligible": sum(
            parse_bool(row["m2_screen_eligible"], field_name="m2_screen_eligible")
            for row in cooccurrence_rows.values()
        ),
        "evaluable_ineligible": sum(
            parse_bool(row["m2_screen_evaluable"], field_name="m2_screen_evaluable")
            and not parse_bool(row["m2_screen_eligible"], field_name="m2_screen_eligible")
            for row in cooccurrence_rows.values()
        ),
        "not_evaluable_axis_indeterminate": sum(
            str(row["m2_screen_eligibility_status"]).startswith(
                "NOT_EVALUABLE_M2_AXIS_INDETERMINATE"
            )
            for row in cooccurrence_rows.values()
        ),
        "not_evaluable_group_count_gt10": sum(
            str(row["m2_screen_eligibility_status"]).startswith(
                "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM"
            )
            for row in cooccurrence_rows.values()
        ),
    }
    observed_counts = {name: int(counts.get(name, -1)) for name in expected_counts}
    if observed_counts != expected_counts:
        raise ContractError(
            f"Logic-independent M2/cooccurrence count drift: {observed_counts} != {expected_counts}"
        )
    aligned_below_power = int(
        payload.get(
            "n_evaluable_sites_with_aligned_axis_below_negative_evaluability_power",
            -1,
        )
    )
    if aligned_below_power < 0:
        raise ContractError("Logic-independent M2 aligned-below-power count is missing")
    return (
        {
            "status": "PASS_LOGIC_INDEPENDENT_RECOUNT",
            "schema_version": payload["schema_version"],
            "contract": payload["contract"],
            "counts": observed_counts,
            "n_evaluable_sites_with_aligned_axis_below_negative_evaluability_power": (
                aligned_below_power
            ),
            "constant_axis_assignment_proof_counts": payload.get(
                "constant_axis_assignment_proof_counts"
            ),
            "production_gate_imported": False,
            "production_gate_functions_called": False,
            "receipt": artifact(path),
        },
        [path, claim_contract, audit_code],
    )


def load_tumor_ref_source_identity_attestation(
    receipt_path: Path | None,
    tumor_ref_manifest_path: Path,
    tumor_ref_manifest: Mapping[str, Any],
) -> tuple[dict[str, Any], list[Path]]:
    if receipt_path is None:
        return (
            {
                "status": "NOT_INCLUDED_INTERMEDIATE_TERMINAL_BUILD",
                "release_gate_pass": False,
                "publishable_task_b_release": False,
                "interpretation": (
                    "Computational intermediate only; a passing bounded retrospective source-file "
                    "identity receipt is required for the final release bundle."
                ),
            },
            [],
        )
    require_nonempty_file(receipt_path, "tumor-REF source identity receipt")
    receipt = load_json(receipt_path, "tumor-REF source identity receipt")
    require_schema_pass(
        receipt,
        "intersubmod.retrospective_running_source_identity.receipt",
        "tumor-REF source identity receipt",
        schema_version=TUMOR_REF_SOURCE_IDENTITY_SCHEMA_VERSION,
    )
    if (
        receipt.get("task_type") != "B_comprehensive_validation"
        or receipt.get("audit_class") != "bounded_retrospective_source_file_identity"
    ):
        raise ContractError("Tumor-REF source identity receipt scope/class drift")
    checks = receipt.get("checks")
    expected_check_keys = set(TUMOR_REF_SOURCE_IDENTITY_TRUE_CHECKS) | {
        "process_start_clock_tolerance_seconds"
    }
    if not isinstance(checks, Mapping) or set(checks) != expected_check_keys:
        raise ContractError("Tumor-REF source identity receipt check set drift")
    if any(checks.get(key) is not True for key in TUMOR_REF_SOURCE_IDENTITY_TRUE_CHECKS):
        raise ContractError("Tumor-REF source identity receipt has a nonpassing check")
    if float(checks.get("process_start_clock_tolerance_seconds", -1)) != 2.0:
        raise ContractError("Tumor-REF source identity clock-tolerance drift")
    verify_declared_artifact(
        receipt.get("tumor_ref_run_manifest"),
        tumor_ref_manifest_path,
        "source identity tumor-REF run manifest",
    )
    snapshot_path = verify_declared_path_artifact(
        receipt.get("snapshot"), "source identity during-execution snapshot"
    )
    snapshot = load_json(snapshot_path, "source identity during-execution snapshot")
    during = receipt.get("source_identity_during_execution")
    after = receipt.get("source_identity_after_execution")
    if (
        not isinstance(during, Mapping)
        or set(during) != TUMOR_REF_SOURCE_IDENTITY_REQUIRED_ROLES
        or after != during
    ):
        raise ContractError("Tumor-REF source identities are missing, extra, or changed")
    manifest_sources = tumor_ref_manifest.get("source_code")
    if (
        not isinstance(manifest_sources, Mapping)
        or set(manifest_sources) != TUMOR_REF_SOURCE_IDENTITY_REQUIRED_ROLES
    ):
        raise ContractError("Tumor-REF run manifest source role set drift")
    for role in sorted(TUMOR_REF_SOURCE_IDENTITY_REQUIRED_ROLES):
        source = during[role]
        if not isinstance(source, Mapping):
            raise ContractError(f"Tumor-REF source identity is malformed for {role}")
        expected = {
            key: source.get(key) for key in ("path", "size_bytes", "sha256")
        }
        if manifest_sources.get(role) != expected:
            raise ContractError(f"Tumor-REF manifest/source receipt mismatch for {role}")
    creator_before = receipt.get("snapshot_creator_source_identity")
    creator_after = receipt.get("snapshot_creator_source_identity_after_execution")
    if not isinstance(creator_before, Mapping) or creator_after != creator_before:
        raise ContractError("Tumor-REF snapshot creator source changed across the audit")
    if (
        snapshot.get("source_identity_during_execution") != during
        or snapshot.get("snapshot_creator_source_identity") != creator_before
    ):
        raise ContractError("Tumor-REF receipt identities do not match the hash-bound snapshot")
    verifier = receipt.get("post_run_verifier_source_identity")
    if not isinstance(verifier, Mapping):
        raise ContractError("Tumor-REF post-run verifier source identity is missing")
    expected_verifier = {
        key: verifier.get(key) for key in ("path", "size_bytes", "sha256")
    }
    trusted_verifier = artifact(TUMOR_REF_SOURCE_IDENTITY_VERIFIER)
    if expected_verifier != trusted_verifier:
        raise ContractError("Tumor-REF receipt was not produced by the trusted v2 verifier")
    if receipt.get("auditor_source_identity_after_execution") != verifier:
        raise ContractError("Tumor-REF post-run auditor/verifier identity mismatch")
    manifest_command = tumor_ref_manifest.get("command")
    process = snapshot.get("process")
    live_cmdline = process.get("cmdline") if isinstance(process, Mapping) else None
    if (
        not isinstance(manifest_command, list)
        or not manifest_command
        or not all(isinstance(token, str) for token in manifest_command)
        or not isinstance(live_cmdline, list)
        or len(live_cmdline) < 2
        or not all(isinstance(token, str) for token in live_cmdline)
        or live_cmdline[1:] != manifest_command
    ):
        raise ContractError("Tumor-REF snapshot and manifest command binding mismatch")
    analyzer_path = str(during["analyzer"].get("path", ""))
    if (
        not source_script_token_matches_attested_source(
            manifest_command[0], analyzer_path
        )
        or snapshot.get("expected_command_fragment") != analyzer_path
    ):
        raise ContractError("Tumor-REF manifest script token is not bound to the analyzer")
    expected_command_binding = {
        "live_python_launcher_token": live_cmdline[0],
        "manifest_script_token": manifest_command[0],
        "manifest_script_token_mode": (
            "absolute"
            if Path(manifest_command[0]).is_absolute()
            else "repo_relative_exact"
        ),
        "attested_analyzer_path": analyzer_path,
        "live_after_launcher_exactly_equals_manifest": True,
        "relative_token_rejects_dot_and_parent_segments": True,
        "repo_relative_token_must_equal_attested_source_relative_to_repo_root": True,
    }
    if receipt.get("command_binding") != expected_command_binding:
        raise ContractError("Tumor-REF source identity command-binding receipt drift")
    if receipt.get("pass_semantics") != (
        "Named producer source files were observed during execution, predated the run, "
        "remained identity-equal afterward, and match the passing producer manifest."
    ):
        raise ContractError("Tumor-REF source identity pass semantics drift")
    limitation = str(receipt.get("limitation", ""))
    if "not a prelaunch lock" not in limitation or "environment attestation" not in limitation:
        raise ContractError("Tumor-REF source identity limitation disclosure is incomplete")
    manifest_finished = parse_utc_timestamp(
        tumor_ref_manifest.get("finished_at_utc"), label="tumor-REF manifest finished"
    )
    receipt_created = parse_utc_timestamp(
        receipt.get("created_at_utc"), label="tumor-REF source receipt created"
    )
    if receipt_created < manifest_finished:
        raise ContractError("Tumor-REF source identity receipt predates producer completion")
    return (
        {
            "status": "VERIFIED_BOUNDED_RETROSPECTIVE_SOURCE_IDENTITY",
            "release_gate_pass": True,
            "publishable_task_b_release": True,
            "audit_class": receipt["audit_class"],
            "receipt": artifact(receipt_path),
            "snapshot": artifact(snapshot_path),
            "source_roles": sorted(TUMOR_REF_SOURCE_IDENTITY_REQUIRED_ROLES),
            "source_sha256": {
                role: str(during[role]["sha256"])
                for role in sorted(TUMOR_REF_SOURCE_IDENTITY_REQUIRED_ROLES)
            },
            "limitation": limitation,
        },
        [receipt_path, snapshot_path],
    )


def load_matched_normal(
    bundle: InputBundle, candidate_keys: set[CandidateKey]
) -> tuple[dict[CandidateKey, dict[str, str]], dict[str, Any], list[Path]]:
    if bundle.matched_normal_dir is None:
        if candidate_keys:
            raise ContractError("Matched-normal output directory is required for nonzero candidates")
        return {}, {"status": "NOT_APPLICABLE_ZERO_SELECTED_CANDIDATES_OMITTED"}, []
    assert bundle.matched_normal_sites is not None
    assert bundle.matched_normal_summary is not None
    assert bundle.matched_normal_receipt is not None
    native_presence = {
        "site table": bundle.matched_normal_sites.is_file(),
        "summary": bundle.matched_normal_summary.is_file(),
        "run receipt": bundle.matched_normal_receipt.is_file(),
    }
    if not candidate_keys and not native_presence["site table"] and not native_presence["summary"]:
        path, receipt = find_explicit_na_receipt(
            bundle.matched_normal_dir,
            prefix="intersubmod.matched_normal_candidate_control",
            schema_versions={MATCHED_NORMAL_PAIRED_RUNNER_SCHEMA_VERSION},
            expected_selection_column=STRICT_SELECTION_COLUMN,
        )
        paths = validate_matched_normal_paired_na_receipt(path, receipt, bundle)
        return (
            {},
            {
                "status": "NOT_APPLICABLE_VALIDATED_RECEIPT",
                "receipt": receipt,
                "paired_receipt": receipt,
                "paired_receipt_path": path,
            },
            paths,
        )
    if not all(native_presence.values()):
        raise ContractError(f"Matched-normal native output set is incomplete: {native_presence}")
    rows = load_keyed_tsv(
        bundle.matched_normal_sites,
        label="matched-normal control table",
        required_fields={
            "normal_called_reads",
            "normal_alt_reads",
            "normal_ref_reads",
            "normal_unknown_reads",
            "normal_focal_callability_gate",
            "normal_control_evaluable",
            "normal_control_not_evaluable_reason",
            "normal_ref_methyl_stable_multigroup",
            "normal_ref_methyl_nonreplication_gate",
            "normal_control_status",
            "normal_evaluable",
            "normal_stable_multigroup",
            "normal_genetic_alt_support_present",
            "tumor_ref_evaluable",
            "tumor_ref_stable_multigroup",
            "primary_group_assignment_coverage",
            "primary_identity_collision_count",
            "primary_identity_missing_count",
        },
    )
    if set(rows) != candidate_keys:
        raise ContractError(
            "Matched-normal/cooccurrence candidate key-set mismatch: "
            f"selected={len(candidate_keys)} normal={len(rows)} "
            f"missing={sorted(candidate_keys - set(rows))[:3]} "
            f"extra={sorted(set(rows) - candidate_keys)[:3]}"
        )
    for key, row in rows.items():
        for field_name in (
            "normal_evaluable",
            "normal_stable_multigroup",
            "normal_genetic_alt_support_present",
            "tumor_ref_evaluable",
            "tumor_ref_stable_multigroup",
        ):
            parse_bool(row[field_name], field_name=field_name)
        called = optional_int(row["normal_called_reads"], field_name="normal_called_reads")
        alt_reads = optional_int(row["normal_alt_reads"], field_name="normal_alt_reads")
        ref_reads = optional_int(row["normal_ref_reads"], field_name="normal_ref_reads")
        unknown_reads = optional_int(
            row["normal_unknown_reads"], field_name="normal_unknown_reads"
        )
        if None in {called, alt_reads, ref_reads, unknown_reads}:
            raise ContractError(f"Matched-normal read counts are missing at {key}")
        if called != alt_reads + ref_reads:
            raise ContractError(f"Matched-normal called-depth reconciliation failed at {key}")
        callability = parse_bool(
            row["normal_focal_callability_gate"],
            field_name="normal_focal_callability_gate",
        )
        if callability != (called >= 10 and ref_reads >= 5):
            raise ContractError(f"Matched-normal focal callability gate drift at {key}")
        control_evaluable = parse_bool(
            row["normal_control_evaluable"], field_name="normal_control_evaluable"
        )
        ref_methyl_evaluable = parse_bool(
            row["normal_evaluable"], field_name="normal_evaluable"
        )
        if control_evaluable != (callability and ref_methyl_evaluable):
            raise ContractError(f"Matched-normal REF-only evaluability drift at {key}")
        if parse_bool(
            row["normal_stable_multigroup"], field_name="normal_stable_multigroup"
        ) != parse_bool(
            row["normal_ref_methyl_stable_multigroup"],
            field_name="normal_ref_methyl_stable_multigroup",
        ):
            raise ContractError(f"Matched-normal deprecated methyl alias drift at {key}")
        coverage = optional_float(
            row["primary_group_assignment_coverage"],
            field_name="primary_group_assignment_coverage",
        )
        collisions = optional_int(
            row["primary_identity_collision_count"],
            field_name="primary_identity_collision_count",
        )
        missing = optional_int(
            row["primary_identity_missing_count"],
            field_name="primary_identity_missing_count",
        )
        if coverage != 1.0 or collisions != 0 or missing != 0:
            raise ContractError(f"Matched-normal primary identity join is not exact at {key}")
    summary = load_json(bundle.matched_normal_summary, "matched-normal summary")
    receipt = load_json(bundle.matched_normal_receipt, "matched-normal run receipt")
    require_schema_pass(
        summary,
        "intersubmod.matched_normal_candidate_control_analysis",
        "matched-normal summary",
        schema_version=MATCHED_NORMAL_ANALYSIS_SCHEMA_VERSION,
    )
    require_schema_pass(
        receipt,
        "intersubmod.matched_normal_candidate_control_analysis_run",
        "matched-normal run receipt",
        schema_version=MATCHED_NORMAL_ANALYSIS_SCHEMA_VERSION,
    )
    analyzer_producer_path = validate_release_bound_producer(
        receipt,
        receipt_path=bundle.matched_normal_receipt,
        producer_module=MATCHED_NORMAL_ANALYZER,
        producer_role="matched_normal_analyzer",
        expected_command=MATCHED_NORMAL_ANALYZER.canonical_task_b_command(),
        label="matched-normal analyzer",
    )
    for field_name in ("source_authority", "code", "source_lock"):
        if summary.get(field_name) != receipt.get(field_name):
            raise ContractError(
                f"Matched-normal summary/analyzer provenance drift for {field_name}"
            )
    for output_path, label in (
        (bundle.matched_normal_sites, "matched-normal site table"),
        (bundle.matched_normal_summary, "matched-normal summary"),
    ):
        if oct(output_path.resolve(strict=True).stat().st_mode & 0o777) != "0o444":
            raise ContractError(f"{label} is not mode 0444")
    expected_pass_semantics = "execution_and_identity_integrity_only_not_background_negativity"
    if summary.get("pass_semantics") != expected_pass_semantics or receipt.get(
        "pass_semantics"
    ) != expected_pass_semantics:
        raise ContractError("Matched-normal analysis pass semantics drift")
    if summary.get("not_evaluable_is_negative_result") is not False or receipt.get(
        "not_evaluable_is_negative_result"
    ) is not False:
        raise ContractError(
            "Matched-normal analysis must declare not_evaluable_is_negative_result=false"
        )
    pooled = summary.get("pooled") or {}
    if int(pooled.get("n_candidates", -1)) != len(candidate_keys):
        raise ContractError("Matched-normal summary candidate count mismatch")
    if pooled.get("all_primary_group_assignments_exact") is not True:
        raise ContractError("Matched-normal summary identity join is not exact")
    counts = receipt.get("counts") or {}
    if counts != pooled:
        raise ContractError("Matched-normal summary/receipt pooled counts mismatch")
    if int(counts.get("n_candidates", -1)) != len(candidate_keys):
        raise ContractError("Matched-normal receipt candidate count mismatch")
    if counts.get("all_primary_group_assignments_exact") is not True:
        raise ContractError("Matched-normal receipt identity join is not exact")
    outputs = receipt.get("outputs") or {}
    verify_declared_artifact(outputs.get("site_table"), bundle.matched_normal_sites, "matched-normal sites")
    verify_declared_artifact(outputs.get("summary"), bundle.matched_normal_summary, "matched-normal summary")
    inputs = receipt.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ContractError("Matched-normal receipt lacks SHA-locked inputs")
    verify_declared_artifact(
        inputs.get("primary_assignments"),
        bundle.screen_assignments,
        "matched-normal primary assignments",
    )
    paired_receipt_path = verify_declared_path_artifact(
        inputs.get("paired_run_receipt"), "matched-normal paired-run receipt"
    )
    paired_receipt = load_json(paired_receipt_path, "matched-normal paired-run receipt")
    paired_paths = [paired_receipt_path]
    if not candidate_keys and paired_receipt.get("schema_name") == (
        "intersubmod.matched_normal_candidate_control.not_applicable"
    ):
        paired_paths = validate_matched_normal_paired_na_receipt(
            paired_receipt_path, paired_receipt, bundle
        )
    else:
        require_schema_pass(
            paired_receipt,
            "intersubmod.matched_normal_candidate_control_run",
            "matched-normal paired-run receipt",
            schema_version=MATCHED_NORMAL_PAIRED_RUNNER_SCHEMA_VERSION,
        )
        paired_producer_path = validate_release_bound_producer(
            paired_receipt,
            receipt_path=paired_receipt_path,
            producer_module=MATCHED_NORMAL_RUNNER,
            producer_role="matched_normal_runner",
            expected_command=MATCHED_NORMAL_RUNNER.canonical_task_b_command(),
            label="matched-normal paired runner",
        )
        paired_paths.append(paired_producer_path)
        summary_artifact_validation = summary.get("paired_artifact_identity_validation")
        receipt_artifact_validation = receipt.get("paired_artifact_identity_validation")
        if summary_artifact_validation != receipt_artifact_validation or not isinstance(
            summary_artifact_validation, Mapping
        ):
            raise ContractError("Matched-normal paired artifact validation drift")
        if summary_artifact_validation.get(
            "all_runner_artifact_set_sha256_recomputed"
        ) is not True:
            raise ContractError("Matched-normal paired artifact SHA-256 was not recomputed")
        validated_samples = summary_artifact_validation.get("samples")
        if not isinstance(validated_samples, Mapping) or not validated_samples:
            raise ContractError("Matched-normal paired artifact validation lacks samples")
        paired_sample_receipts = paired_receipt.get("receipts")
        if not isinstance(paired_sample_receipts, list) or not paired_sample_receipts:
            raise ContractError("Matched-normal paired receipt lacks sample receipts")
        declared_digests: dict[str, str] = {}
        for sample_receipt in paired_sample_receipts:
            if not isinstance(sample_receipt, Mapping):
                raise ContractError("Malformed matched-normal paired sample receipt")
            sample = str(sample_receipt.get("sample", ""))
            validation = sample_receipt.get("validation")
            if not sample or not isinstance(validation, Mapping):
                raise ContractError("Matched-normal paired sample validation is missing")
            digest = str(validation.get("artifact_set_sha256", ""))
            if len(digest) != 64:
                raise ContractError("Matched-normal paired sample artifact digest is invalid")
            if sample in declared_digests:
                raise ContractError("Duplicate matched-normal paired sample receipt")
            declared_digests[sample] = digest
        if set(validated_samples) != set(declared_digests):
            raise ContractError("Matched-normal paired artifact validation sample-set drift")
        for sample, declared_digest in declared_digests.items():
            observed = validated_samples[sample]
            if (
                not isinstance(observed, Mapping)
                or observed.get("pass") is not True
                or observed.get("artifact_set_sha256") != declared_digest
                or int(observed.get("artifacts_verified", 0)) <= 0
            ):
                raise ContractError(
                    f"Matched-normal paired artifact validation failed for {sample}"
                )
    return (
        rows,
        {
            "status": (
                "COMPLETE"
                if candidate_keys
                else "NOT_APPLICABLE_VALIDATED_NATIVE_ZERO_ROWS"
            ),
            "summary": summary,
            "receipt": receipt,
            "paired_receipt": paired_receipt,
            "receipt_path": bundle.matched_normal_receipt,
            "paired_receipt_path": paired_receipt_path,
        },
        [
            bundle.matched_normal_sites,
            bundle.matched_normal_summary,
            bundle.matched_normal_receipt,
            analyzer_producer_path,
            *paired_paths,
        ],
    )


def load_cn_ccf_annotations(
    path: Path | None,
    candidate_keys: set[CandidateKey],
) -> tuple[dict[CandidateKey, dict[str, Any]], dict[str, Any], list[Path]]:
    if path is None:
        if not candidate_keys:
            return (
                {},
                {
                    "status": "NOT_APPLICABLE_ZERO_SELECTED_CANDIDATES_OMITTED",
                    "reason": "ZERO_SELECTED_CANDIDATES",
                    "is_negative_result": False,
                    "c1_formed": False,
                    "default_diploid_imputation_allowed": False,
                },
                [],
            )
        return (
            {},
            {
                "status": "NOT_RUN",
                "reason": "NO_AUTHORITY_REVIEWED_CN_CCF_ANNOTATION_SOURCE",
                "default_diploid_imputation_allowed": False,
            },
            [],
        )
    resolved = path.expanduser().resolve()
    if resolved.is_dir():
        annotation_path = resolved / "candidate_cn_ccf_annotations.tsv.gz"
        receipt_path = resolved / "receipt.json"
    elif resolved.name == "candidate_cn_ccf_annotations.tsv.gz":
        annotation_path = resolved
        receipt_path = resolved.parent / "receipt.json"
    else:
        raise ContractError(
            "CN/CCF input must be the native annotation directory or "
            "candidate_cn_ccf_annotations.tsv.gz"
        )
    require_nonempty_file(annotation_path, "native CN/CCF annotation")
    receipt = load_json(receipt_path, "native CN/CCF receipt")
    require_schema_pass(
        receipt,
        "intersubmod.candidate_cn_ccf_annotation_receipt",
        "native CN/CCF receipt",
        schema_version=CN_CCF_ANNOTATION_SCHEMA_VERSION,
    )
    cn_producer_path = validate_release_bound_producer(
        receipt,
        receipt_path=receipt_path,
        producer_module=CN_CCF_ANNOTATOR,
        producer_role="cn_ccf_annotator",
        expected_command=CN_CCF_ANNOTATOR.canonical_task_b_command(),
        label="native CN/CCF annotator",
    )
    for output_path, label in (
        (annotation_path, "native CN/CCF annotation"),
        (receipt_path, "native CN/CCF receipt"),
    ):
        if oct(output_path.resolve(strict=True).stat().st_mode & 0o777) != "0o444":
            raise ContractError(f"{label} is not mode 0444")
    expected_zero = not candidate_keys
    expected_status = "NOT_APPLICABLE" if expected_zero else "PASS"
    expected_execution = "NOT_APPLICABLE" if expected_zero else "EXECUTION_PASS"
    if receipt.get("status") != expected_status or receipt.get(
        "execution_status"
    ) != expected_execution:
        raise ContractError("Native CN/CCF receipt status/execution_status drift")
    if int(receipt.get("n_selected_candidates", -1)) != len(candidate_keys):
        raise ContractError("Native CN/CCF selected-candidate count mismatch")
    if expected_zero and receipt.get("reason") != "ZERO_SELECTED_CANDIDATES":
        raise ContractError("Native CN/CCF zero-row receipt reason drift")
    if not expected_zero and receipt.get("reason") not in {None, ""}:
        raise ContractError("Native CN/CCF nonzero receipt has an unexpected reason")
    if receipt.get("claim_ceiling") != CN_CCF_CLAIM_CEILING:
        raise ContractError("Native CN/CCF receipt claim ceiling drift")
    if receipt.get("pass_semantics") != CN_CCF_PASS_SEMANTICS:
        raise ContractError("Native CN/CCF receipt pass semantics drift")
    interpretation = receipt.get("scientific_interpretation")
    if not isinstance(interpretation, Mapping) or interpretation.get(
        "is_negative_result"
    ) is not False or interpretation.get("c1_formed") is not False:
        raise ContractError("Native CN/CCF receipt scientific interpretation drift")

    output_reference = receipt.get("output")
    verify_declared_artifact(
        output_reference, annotation_path, "native CN/CCF annotation output"
    )
    if not isinstance(output_reference, Mapping) or tuple(
        output_reference.get("columns") or []
    ) != CN_CCF_OUTPUT_COLUMNS:
        raise ContractError("Native CN/CCF receipt output columns drift")
    with open_text(annotation_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != CN_CCF_OUTPUT_COLUMNS:
            raise ContractError("Native CN/CCF annotation header drift")
        rows = [dict(row) for row in reader]

    receipt_input = receipt.get("input")
    input_path = verify_declared_path_artifact(receipt_input, "native CN/CCF candidate input")
    if not isinstance(receipt_input, Mapping):
        raise ContractError("Native CN/CCF receipt input is malformed")
    rows_in = int(receipt_input.get("rows_in", -1))
    rows_read_total = int(receipt_input.get("rows_read_total", -1))
    if rows_in != len(candidate_keys) or rows_read_total < rows_in:
        raise ContractError("Native CN/CCF input row counts drift")
    conservation = receipt.get("conservation")
    expected_conservation = {
        "rows_in": len(candidate_keys),
        "rows_out": len(rows),
        "rows_in_equals_rows_out": len(candidate_keys) == len(rows),
    }
    if conservation != expected_conservation or len(rows) != len(candidate_keys):
        raise ContractError("Native CN/CCF row conservation failed")
    if int(output_reference.get("rows_out", -1)) != len(rows):
        raise ContractError("Native CN/CCF output row count mismatch")

    authority = receipt.get("authority")
    if not isinstance(authority, Mapping):
        raise ContractError("Native CN/CCF receipt lacks authority chain")
    source_paths: list[Path] = [
        annotation_path,
        receipt_path,
        input_path,
        cn_producer_path,
    ]
    for field_name in ("config", "input_provenance", "analysis_summary"):
        authority_path = verify_declared_path_artifact(
            authority.get(field_name), f"native CN/CCF authority {field_name}"
        )
        source_paths.append(authority_path)
    all_source_hashes = authority.get("all_source_hashes")
    if not isinstance(all_source_hashes, Mapping) or not all_source_hashes:
        raise ContractError("Native CN/CCF authority all_source_hashes is empty")
    verified_source_hashes: dict[str, str] = {}
    for source_path_text, declared_sha in all_source_hashes.items():
        source_path = require_nonempty_file(
            Path(str(source_path_text)).expanduser().resolve(),
            "native CN/CCF authority source",
        )
        observed_sha = sha256(source_path)
        if str(declared_sha) != observed_sha:
            raise ContractError(f"Native CN/CCF authority source SHA-256 mismatch: {source_path}")
        verified_source_hashes[str(source_path)] = observed_sha
        source_paths.append(source_path)
    by_key: dict[CandidateKey, dict[str, Any]] = {}
    for row_number, row in enumerate(rows, 2):
        key = candidate_key(row, source=f"native CN/CCF annotation row {row_number}")
        if key in by_key:
            raise ContractError(f"Duplicate native CN/CCF annotation key: {key}")
        if key not in candidate_keys:
            raise ContractError(f"Native CN/CCF annotation key is not a G2 candidate: {key}")
        if row["mutation_id"] != f"{key[1]}:{key[2]}:{key[3]}>{key[4]}":
            raise ContractError(f"Native CN/CCF mutation_id drift at {key}")
        if row["callset_status"] != "INPUT_CANDIDATE_KEY_VALIDATED_NOT_RECHECKED":
            raise ContractError(f"Native CN/CCF callset status drift at {key}")
        cn_status = row["cn_status"]
        pyclone_status = row["pyclone_status"]
        sensitivity_status = row["pyclone_sensitivity_status"]
        if cn_status not in CN_STATUS_ENUM:
            raise ContractError(f"Invalid native CN status at {key}: {cn_status!r}")
        if pyclone_status not in PYCLONE_STATUS_ENUM or sensitivity_status not in PYCLONE_STATUS_ENUM:
            raise ContractError(f"Invalid native PyClone status at {key}")
        if "RECEIPT_FAIL" in {cn_status, pyclone_status, sensitivity_status}:
            raise ContractError(f"Passing native CN/CCF receipt contains RECEIPT_FAIL at {key}")
        independent_cn = parse_bool(row["independent_cn"], field_name="independent_cn")
        if key[0] == "HCC1395_DORADO":
            if independent_cn or row["cn_authority_sample"] != "HCC1395" or row[
                "cn_transfer_policy"
            ] != "SAME_CELL_LINE_SHARED_CN_SENSITIVITY" or cn_status not in {
                "SHARED_CN_SENSITIVITY",
                "NO_CN_SEGMENT",
            }:
                raise ContractError(f"DORADO shared-CN sensitivity contract drift at {key}")
        elif key[0] in {"COLO829", "HCC1937"}:
            if independent_cn or cn_status != "BLOCKED_CN_MISFIT" or pyclone_status != (
                "BLOCKED_CN_MISFIT"
            ) or sensitivity_status != "BLOCKED_CN_MISFIT" or row[
                "cn_transfer_policy"
            ] != "BLOCKED_CN_MISFIT_NO_CN2_IMPUTATION":
                raise ContractError(f"Blocked native CN contract drift at {key}")
            for field_name in (
                "savana_total_cn_raw",
                "savana_minor_cn_raw",
                "savana_total_cn_discrete",
                "savana_major_cn_discrete",
                "savana_minor_cn_discrete",
                "purity_model_value",
            ):
                if row[field_name] not in {None, ""}:
                    raise ContractError(f"Blocked sample received imputed native CN at {key}")
        elif not independent_cn or cn_status == "SHARED_CN_SENSITIVITY":
            raise ContractError(f"Sample-specific native CN independence drift at {key}")
        if cn_status in {"AVAILABLE_EXACT_SEGMENT", "SHARED_CN_SENSITIVITY"}:
            for field_name in (
                "cn_segment_id",
                "cn_segment_start",
                "cn_segment_end",
                "savana_total_cn_raw",
                "savana_minor_cn_raw",
                "savana_total_cn_discrete",
                "savana_major_cn_discrete",
                "savana_minor_cn_discrete",
                "purity_model_value",
                "cn_source_sha256",
            ):
                if row[field_name] in {None, ""}:
                    raise ContractError(f"Native exact-segment CN field is missing at {key}: {field_name}")
        if pyclone_status == "MATCHED_PRIMARY":
            for field_name in (
                "pyclone_fit_local_cluster_id",
                "pyclone_vi_cellular_prevalence",
                "pyclone_vi_cellular_prevalence_std",
                "pyclone_vi_assignment_probability",
            ):
                if row[field_name] in {None, ""}:
                    raise ContractError(f"Native matched PyClone field is missing at {key}: {field_name}")
        authority_hashes = {
            "authority_config_sha256": str(authority["config"]["sha256"]),
            "input_provenance_sha256": str(authority["input_provenance"]["sha256"]),
            "analysis_summary_sha256": str(authority["analysis_summary"]["sha256"]),
        }
        for field_name, expected_sha in authority_hashes.items():
            if row[field_name] != expected_sha:
                raise ContractError(f"Native CN/CCF row authority hash drift at {key}: {field_name}")
        source_hash_values = set(verified_source_hashes.values())
        for field_name in (
            "cn_source_sha256",
            "pyclone_primary_metadata_sha256",
            "pyclone_primary_results_sha256",
            "pyclone_primary_status_sha256",
            "pyclone_sensitivity_metadata_sha256",
            "pyclone_sensitivity_results_sha256",
            "pyclone_sensitivity_status_sha256",
        ):
            value = str(row[field_name]).strip()
            if value and value not in source_hash_values:
                raise ContractError(f"Native CN/CCF row source hash is not authority-locked at {key}")
        if row["claim_ceiling"] != CN_CCF_CLAIM_CEILING:
            raise ContractError(f"Native CN/CCF row claim ceiling drift at {key}")
        matching_cn_paths = sorted(
            source_path_text
            for source_path_text, source_sha in verified_source_hashes.items()
            if source_sha == row["cn_source_sha256"] and row["cn_source_sha256"]
        )
        normalized = dict(row)
        normalized.update(
            {
                "independent_cn": independent_cn,
                "_cn_source_path": matching_cn_paths[0] if matching_cn_paths else None,
            }
        )
        by_key[key] = normalized
    if set(by_key) != candidate_keys:
        raise ContractError(
            "Native CN/CCF candidate key-set mismatch: "
            f"expected={len(candidate_keys)} observed={len(by_key)}"
        )
    expected_status_counts = {
        field_name: dict(sorted(Counter(row[field_name] for row in rows).items()))
        for field_name in ("cn_status", "pyclone_status", "pyclone_sensitivity_status")
    }
    if receipt.get("status_counts") != expected_status_counts:
        raise ContractError("Native CN/CCF receipt status counts drift")
    status_enums = receipt.get("status_enums")
    if not isinstance(status_enums, Mapping) or set(status_enums.get("cn_status") or []) != (
        CN_STATUS_ENUM
    ) or set(status_enums.get("pyclone_status") or []) != PYCLONE_STATUS_ENUM:
        raise ContractError("Native CN/CCF receipt status enum drift")
    return (
        by_key,
        {
            "status": (
                "NOT_APPLICABLE_VALIDATED_NATIVE_ZERO_ROWS" if expected_zero else "COMPLETE"
            ),
            "annotation_rows": len(by_key),
            "candidate_rows_without_annotation": 0,
            "default_diploid_imputation_allowed": False,
            "annotation_artifact": artifact(annotation_path),
            "receipt_artifact": artifact(receipt_path),
            "receipt_status": receipt["status"],
            "execution_status": receipt["execution_status"],
            "c1_formed": False,
        },
        list(dict.fromkeys(source_paths)),
    )


def dataset_role(sample: str) -> str:
    if sample == "HCC1395":
        return "technical_pair_primary"
    if sample == "HCC1395_DORADO":
        return "technical_replicate"
    return "independent_biological_dataset"


def biological_n_contribution(sample: str) -> int:
    return 0 if sample == "HCC1395_DORADO" else 1


def _typed_value(row: Mapping[str, Any], field_name: str, kind: str) -> Any:
    value = row.get(field_name)
    if kind == "bool":
        return optional_bool(value, field_name=field_name)
    if kind == "int":
        return optional_int(value, field_name=field_name)
    if kind == "float":
        return optional_float(value, field_name=field_name)
    return None if value is None or value == "" else value


CANDIDATE_FIELDS = (
    "candidate_key",
    "sample",
    "biological_id",
    "dataset_role",
    "biological_n_contribution",
    "biological_site_key",
    "chrom",
    "pos",
    "ref",
    "alt",
    "truth_label",
    "ssnv_branch",
    "component_id",
    "selected_group_id",
    "m1_status",
    "m1_status_semantics",
    "modal_assignment_ari_min",
    "m2_status",
    "m2_screen_gate_contract",
    "m2_screen_evaluable",
    "m2_screen_eligibility_status",
    "m2_axis_statuses",
    "m2_indeterminate_axes",
    "m2_low_power_axes",
    "hp_axis_confound",
    "technical_axis_confound",
    "residual_unexplained_multigroup",
    "g1_status",
    "n_pair_by_confirmed",
    "pair_by_confirmed_positions",
    "g2_status",
    "n_spatially_separated_pair_by_20bp",
    "spatially_separated_pair_by_positions_20bp",
    "top_marker_positions",
    "n_top_marker_pair_by_confirmed",
    "n_top_spatially_separated_pair_by_20bp",
    "top_spatially_separated_pair_by_positions_20bp",
    "joint_signature_n_complete_reads",
    "joint_signature_testable",
    "joint_signature_p_conditional_perm",
    "joint_signature_permutations",
    "joint_signature_permutable",
    "joint_signature_conditional_status",
    "joint_signature_sensitivity_pass",
    "joint_signature_global_fdr_family_status",
    "joint_signature_q_global_bh",
    "joint_signature_q_global_by",
    "joint_signature_global_bh_discovery",
    "joint_signature_global_by_discovery",
    "r1_status",
    "r1_reason",
    "strict_analysis_class",
    "strict_assignment_concordance_ari_min",
    "tumor_ref_status",
    "tumor_ref_evaluable",
    "tumor_ref_stable_null_multigroup",
    "tumor_ref_read_count",
    "normal_control_status",
    "normal_called_reads",
    "normal_ref_reads",
    "normal_alt_reads",
    "normal_unknown_reads",
    "normal_focal_callability_gate",
    "normal_ref_methyl_evaluable",
    "normal_ref_methyl_stable_multigroup",
    "normal_ref_methyl_nonreplication_gate",
    "n_same_pair_four_state_witnesses",
    "same_pair_four_state_witness_keys",
    "n_four_state_compatible_formal_pair_opportunities",
    "four_state_compatible_formal_pair_opportunity_keys",
    "b1_prespecified_pair_key",
    "b1_prespecified_pair_relation_status",
    "b1_prespecified_pair_is_witness",
    "b1_pair_selection_contract",
    "b1_uses_posthoc_compatible_pair_search",
    "b1_status",
    "b1_reason",
    "technical_replication_status",
    "technical_exact_pair_opportunities",
    "technical_replicated_exact_pairs",
    "technical_biological_n",
    "cn_ccf_status",
    "cn_ccf_source_path",
    "cn_ccf_source_sha256",
    "cn_ccf_authority_reviewed",
    "cn_ccf_exact_locus_covered",
    "cn_total",
    "cn_minor",
    "purity",
    "mutation_multiplicity",
    "ccf",
    *CN_CCF_NATIVE_DETAIL_FIELDS,
    "c1_status",
    "c1_reason",
    "l1_status",
    "l2_status",
    "deprecated_aliases",
)


PAIR_DETAIL_FIELDS = (
    "witness_pair_key",
    "sample",
    "biological_id",
    "focal_chrom",
    "focal_pos",
    "focal_ref",
    "focal_alt",
    "partner_chrom",
    "partner_pos",
    "partner_ref",
    "partner_alt",
    "focal_truth_label",
    "partner_truth_label",
    "focal_ssnv_branch",
    "partner_ssnv_branch",
    "focal_component_id",
    "partner_component_id",
    "distance_bp",
    "g1_same_pair_status",
    "exact_state_space_status",
    "exact_p",
    "global_bh_q",
    "global_by_q",
    "conditional_p",
    "conditional_permutations",
    "conditional_status",
    "conditional_sensitivity_pass",
    "callability_noncallable_core_reads",
    "callability_gate_status",
    "callability_gate_pass",
    "cramers_v",
    "delta_alt_fraction",
    "four_state_testable",
    "four_state_rr",
    "four_state_ar",
    "four_state_ra",
    "four_state_aa",
    "four_state_o",
    "four_state_x",
    "four_state_called_depth",
    "four_state_focal_ref_depth",
    "four_state_focal_alt_depth",
    "four_state_partner_ref_depth",
    "four_state_partner_alt_depth",
    "four_state_error_ceiling",
    "four_state_error_model_confidence",
    "four_state_familywise_confidence",
    "four_state_relation_family_size",
    "four_state_multiplicity_method",
    "four_state_minimum_zero_violation_depth",
    "focal_ancestor_violation_p_exact",
    "focal_ancestor_violation_upper_bound",
    "focal_ancestor_violation_status",
    "partner_ancestor_violation_p_exact",
    "partner_ancestor_violation_upper_bound",
    "partner_ancestor_violation_status",
    "branching_violation_p_exact",
    "branching_violation_upper_bound",
    "branching_violation_status",
    "four_state_relation_status",
    "four_state_compatible_relation_models",
    "four_state_n_compatible_relation_models",
    "same_pair_four_state_witness",
    "b1_prespecified_pair",
    "topology_scope",
    "topology_region",
    "topology_order_status",
    "tumor_ref_status",
    "tumor_ref_read_count",
    "normal_control_status",
    "normal_called_reads",
    "normal_ref_reads",
    "normal_alt_reads",
    "normal_unknown_reads",
    "technical_exact_pair_present",
    "technical_replication_pass",
    "technical_replication_status",
    "technical_core_read_name_overlap_present",
    "technical_molecule_independence_status",
    "technical_biological_n",
)


CLAIM_RULES = {
    "M1": {
        "name": "focal-ALT operational stable methyl-multigroup screen flag",
        "denominator": (
            "all LongPhase-S recalibrated FILTER=PASS chr1-22 biallelic focal "
            "dataset-sites in scope"
        ),
        "guardrail": (
            "M1 FAIL means not flagged by the operational screen, including technically "
            "non-evaluable sites; it is not evidence of homogeneity. Global null-validity "
            "state was not exported for nonstable sites."
        ),
    },
    "M2": {
        "name": "residual robust epigenetic partition",
        "denominator": (
            "M1 PASS sites with 2-10 observed methyl groups where each measured axis is "
            "either an observed aligned confound, adequately powered for non-alignment, "
            "or assignment-proven constant"
        ),
        "guardrail": "Residual measured-axis structure is not a somatic clone.",
    },
    "G1": {
        "name": "LongPhase-S callset-anchored local co-segregation",
        "denominator": "M2 sites with at least one exact-testable partner",
        "guardrail": "A LongPhase-S PASS partner is not independently confirmed somatic truth.",
    },
    "G2": {
        "name": "multi-marker molecular-haplotype base candidate",
        "denominator": "M2 sites with >=2 spatially separated exact-identifiable top markers and an executable 999-permutation joint conditional test",
        "guardrail": "Spatial separation does not imply statistical independence or a cellular clone.",
    },
    "R1": {
        "name": "strict methyl-partition robustness",
        "denominator": "G2 PASS sites with formal 999-permutation x 10-seed strict evaluation",
        "guardrail": "This is a robustness audit, not post-selection FDR confirmation.",
    },
    "B1": {
        "name": "background-controlled molecular-haplotype candidate",
        "denominator": "G2 sites with evaluable strict, tumor-REF, normal REF-only, and same-pair four-state controls",
        "guardrail": "Bulk molecular-haplotype evidence does not establish cellular identity or clone count.",
    },
    "C1": {
        "name": "CN/CCF-conditioned candidate",
        "denominator": "B1 PASS sites covered by authority-reviewed exact-locus CN/purity/CCF models",
        "guardrail": "The result is conditional on the reviewed CN/CCF model.",
    },
    "L1": {
        "name": "cellular subclone supported",
        "denominator": "C1 PASS sites with orthogonal same-cell-population evidence",
        "guardrail": "No single-cell, colony, spatial, or multi-region identity evidence is integrated here.",
    },
    "L2": {
        "name": "lineage/order supported",
        "denominator": "L1 PASS sites with identifiable >=3-state order evidence",
        "guardrail": "Four-state bulk compatibility alone cannot identify a unique lineage tree.",
    },
}


def spatially_separated_count(positions: Iterable[int], minimum_bp: int = 20) -> int:
    selected: list[int] = []
    for position in sorted(set(int(value) for value in positions)):
        if not selected or position - selected[-1] >= minimum_bp:
            selected.append(position)
    return len(selected)


def strict_r1_status(row: Mapping[str, Any]) -> tuple[str, str]:
    analysis_class = str(row.get("strict_analysis_class", ""))
    robustness_status = str(row.get("strict_methyl_partition_robustness_status", ""))
    if analysis_class == "EXPLORATORY_ONLY":
        return "NOT_EVALUABLE", "STRICT_CONFIGURATION_EXPLORATORY_ONLY"
    if robustness_status == "NOT_EVALUABLE":
        return (
            "NOT_EVALUABLE",
            str(row.get("strict_not_evaluable_reason") or "STRICT_NOT_EVALUABLE"),
        )
    passed = parse_bool(
        row["strict_null_robustness_pass"], field_name="strict_null_robustness_pass"
    )
    return ("PASS", "PASS") if passed else (
        "FAIL",
        str(row.get("strict_failure_reason") or "STRICT_ROBUSTNESS_FAIL"),
    )


def tumor_ref_claim_status(row: Mapping[str, Any]) -> tuple[str, str]:
    evaluable = parse_bool(row["ref_evaluable"], field_name="ref_evaluable")
    if not evaluable:
        return "NOT_EVALUABLE", str(row.get("ref_not_testable_reason") or "TUMOR_REF_NOT_EVALUABLE")
    stable = parse_bool(
        row["ref_stable_null_multigroup"], field_name="ref_stable_null_multigroup"
    )
    return (
        "FAIL",
        "TUMOR_REF_LENIENT_COARSE_MODAL_PARTITION_REPRODUCED",
    ) if stable else (
        "PASS",
        "TUMOR_REF_LENIENT_COARSE_MODAL_PARTITION_NOT_REPRODUCED",
    )


def normal_claim_status(row: Mapping[str, Any]) -> tuple[str, str]:
    evaluable = parse_bool(
        row["normal_control_evaluable"], field_name="normal_control_evaluable"
    )
    if not evaluable:
        return (
            "NOT_EVALUABLE",
            str(row.get("normal_control_not_evaluable_reason") or "NORMAL_NOT_EVALUABLE"),
        )
    if int(row["normal_alt_reads"]) != 0:
        return "FAIL", "NORMAL_ALT_SUPPORT_PRESENT"
    stable = parse_bool(
        row["normal_ref_methyl_stable_multigroup"],
        field_name="normal_ref_methyl_stable_multigroup",
    )
    return (
        "FAIL",
        "NORMAL_REF_LENIENT_COARSE_MODAL_PARTITION_REPRODUCED",
    ) if stable else (
        "PASS",
        "NORMAL_REF_LENIENT_COARSE_MODAL_PARTITION_NOT_REPRODUCED_AND_ALT_ZERO",
    )


def pair_is_four_state_witness(row: Mapping[str, Any]) -> bool:
    return bool(
        parse_bool(
            row["endpoint_a_formal_pair_by_confirmed"],
            field_name="endpoint_a_formal_pair_by_confirmed",
        )
        and parse_bool(
            row["endpoint_b_complete_four_state_testable"],
            field_name="endpoint_b_complete_four_state_testable",
        )
        and str(row["endpoint_b_relation_compatibility"]) in COMPATIBLE_RELATIONS
    )


def select_b1_prespecified_pair(
    formal_pairs: Sequence[Mapping[str, Any]],
) -> Mapping[str, Any] | None:
    """Select one G1-formal pair without inspecting its four-state result."""
    if not formal_pairs:
        return None
    return min(
        formal_pairs,
        key=lambda row: (
            -int(row["endpoint_a_n_informative"]),
            abs(int(row["distance_bp"])),
            str(row["partner_chrom"]),
            int(row["partner_pos"]),
            str(row["partner_ref"]),
            str(row["partner_alt"]),
        ),
    )


def four_state_claim_decision(row: Mapping[str, Any]) -> tuple[str, str]:
    complete = parse_bool(
        row["endpoint_b_complete_four_state_testable"],
        field_name="endpoint_b_complete_four_state_testable",
    )
    relation = str(row["endpoint_b_relation_compatibility"])
    if not complete:
        return "NOT_EVALUABLE", "PRESPECIFIED_G1_PAIR_HAS_INSUFFICIENT_FOUR_STATE_DEPTH"
    if relation in COMPATIBLE_RELATIONS:
        return "PASS", "PRESPECIFIED_G1_PAIR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    if relation in UNRESOLVED_FOUR_STATE_RELATIONS:
        return "NOT_EVALUABLE", f"PRESPECIFIED_G1_PAIR_{relation}"
    if relation == "INCOMPATIBLE_OR_COMPLEX_UNDER_FIXED_ERROR_MODEL":
        return "FAIL", "PRESPECIFIED_G1_PAIR_INCOMPATIBLE_OR_COMPLEX_UNDER_FIXED_ERROR_MODEL"
    raise ContractError(f"Unexpected four-state relation in B1 decision: {relation!r}")


def cn_ccf_claim_status(
    row: Mapping[str, Any] | None,
) -> tuple[str, str]:
    if row is None:
        return "NOT_EVALUABLE", "NO_NATIVE_AUTHORITY_LOCKED_CN_PYCLONE_ANNOTATION"
    sample = str(row.get("sample", ""))
    cn_status = str(row.get("cn_status", ""))
    if sample in {"COLO829", "HCC1937"} or cn_status == "BLOCKED_CN_MISFIT":
        return "NOT_EVALUABLE", "BLOCKED_CN_MISFIT_NO_CN2_IMPUTATION"
    if sample == "HCC1395_DORADO" or cn_status == "SHARED_CN_SENSITIVITY":
        return (
            "NOT_EVALUABLE",
            "SHARED_CN_SENSITIVITY_NOT_INDEPENDENT_AND_NO_PRESPECIFIED_JOINT_WITNESS_MODEL",
        )
    if cn_status == "NO_CN_SEGMENT":
        return "NOT_EVALUABLE", "NO_AUTHORITY_LOCKED_CN_SEGMENT_AT_FOCAL_LOCUS"
    return (
        "NOT_EVALUABLE",
        "FOCAL_ONLY_CN_PYCLONE_ANNOTATION_NO_PRESPECIFIED_JOINT_WITNESS_PAIR_MODEL",
    )


def build_pair_detail(
    data: IntegrationData,
    focal_key: CandidateKey,
    pair: Mapping[str, Any],
    tumor_status: str,
    normal_status: str,
    b1_prespecified_pair_identity: PairKey | None,
) -> dict[str, Any]:
    state_counts = optional_json(
        pair.get("endpoint_b_state_counts"), field_name="endpoint_b_state_counts"
    ) or {}
    partner_counts = optional_json(
        pair.get("all_partner_call_counts"), field_name="all_partner_call_counts"
    ) or {}
    pair_identity = pair_key(pair, source="candidate pair detail")
    return {
        "witness_pair_key": compact_json(list(pair_identity)),
        "sample": focal_key[0],
        "biological_id": BIOLOGICAL_IDS[focal_key[0]],
        "focal_chrom": pair["focal_chrom"],
        "focal_pos": int(pair["focal_pos"]),
        "focal_ref": pair["focal_ref"],
        "focal_alt": pair["focal_alt"],
        "partner_chrom": pair["partner_chrom"],
        "partner_pos": int(pair["partner_pos"]),
        "partner_ref": pair["partner_ref"],
        "partner_alt": pair["partner_alt"],
        "focal_truth_label": pair["focal_truth_label"],
        "partner_truth_label": pair["partner_truth_label"],
        "focal_ssnv_branch": pair["focal_ssnv_branch"],
        "partner_ssnv_branch": pair["partner_ssnv_branch"],
        "focal_component_id": row_value(pair, "focal_component_id"),
        "partner_component_id": row_value(pair, "partner_component_id"),
        "distance_bp": optional_int(pair.get("distance_bp"), field_name="distance_bp"),
        "g1_same_pair_status": claim_status(
            parse_bool(
                pair["endpoint_a_formal_pair_by_confirmed"],
                field_name="endpoint_a_formal_pair_by_confirmed",
            ),
            evaluable=optional_float(
                pair.get("endpoint_a_p_fixed_margins_exact"),
                field_name="endpoint_a_p_fixed_margins_exact",
            )
            is not None,
        ),
        "exact_state_space_status": pair["endpoint_a_exact_state_space_status"],
        "exact_p": optional_float(
            pair.get("endpoint_a_p_fixed_margins_exact"),
            field_name="endpoint_a_p_fixed_margins_exact",
        ),
        "global_bh_q": optional_float(
            pair.get("endpoint_a_q_global_bh"), field_name="endpoint_a_q_global_bh"
        ),
        "global_by_q": optional_float(
            pair.get("endpoint_a_q_global_by"), field_name="endpoint_a_q_global_by"
        ),
        "conditional_p": optional_float(
            pair.get("endpoint_a_p_conditional_perm"),
            field_name="endpoint_a_p_conditional_perm",
        ),
        "conditional_permutations": optional_int(
            pair.get("endpoint_a_permutations"), field_name="endpoint_a_permutations"
        ),
        "conditional_status": pair["endpoint_a_conditional_status"],
        "conditional_sensitivity_pass": parse_bool(
            pair["endpoint_a_conditional_sensitivity_pass"],
            field_name="endpoint_a_conditional_sensitivity_pass",
        ),
        "callability_noncallable_core_reads": optional_int(
            pair["callability_noncallable_core_reads"],
            field_name="callability_noncallable_core_reads",
        ),
        "callability_gate_status": pair["callability_gate_status"],
        "callability_gate_pass": parse_bool(
            pair["callability_gate_pass"], field_name="callability_gate_pass"
        ),
        "cramers_v": optional_float(
            pair.get("endpoint_a_cramers_v"), field_name="endpoint_a_cramers_v"
        ),
        "delta_alt_fraction": optional_float(
            pair.get("endpoint_a_delta_alt_fraction"),
            field_name="endpoint_a_delta_alt_fraction",
        ),
        "four_state_testable": parse_bool(
            pair["endpoint_b_complete_four_state_testable"],
            field_name="endpoint_b_complete_four_state_testable",
        ),
        "four_state_rr": int(state_counts.get("RR", 0)),
        "four_state_ar": int(state_counts.get("AR", 0)),
        "four_state_ra": int(state_counts.get("RA", 0)),
        "four_state_aa": int(state_counts.get("AA", 0)),
        "four_state_o": int(state_counts.get("O", partner_counts.get("O", 0))),
        "four_state_x": int(state_counts.get("X", partner_counts.get("X", 0))),
        "four_state_called_depth": optional_int(
            pair.get("endpoint_b_n_called_depth"), field_name="endpoint_b_n_called_depth"
        ),
        "four_state_focal_ref_depth": optional_int(
            pair.get("endpoint_b_n_focal_ref"), field_name="endpoint_b_n_focal_ref"
        ),
        "four_state_focal_alt_depth": optional_int(
            pair.get("endpoint_b_n_focal_alt"), field_name="endpoint_b_n_focal_alt"
        ),
        "four_state_partner_ref_depth": optional_int(
            pair.get("endpoint_b_n_partner_ref"), field_name="endpoint_b_n_partner_ref"
        ),
        "four_state_partner_alt_depth": optional_int(
            pair.get("endpoint_b_n_partner_alt"), field_name="endpoint_b_n_partner_alt"
        ),
        "four_state_error_ceiling": optional_float(
            pair.get("endpoint_b_error_ceiling"), field_name="endpoint_b_error_ceiling"
        ),
        "four_state_error_model_confidence": optional_float(
            pair.get("endpoint_b_error_model_confidence"),
            field_name="endpoint_b_error_model_confidence",
        ),
        "four_state_familywise_confidence": optional_float(
            pair.get("endpoint_b_familywise_confidence"),
            field_name="endpoint_b_familywise_confidence",
        ),
        "four_state_relation_family_size": optional_int(
            pair.get("endpoint_b_relation_family_size"),
            field_name="endpoint_b_relation_family_size",
        ),
        "four_state_multiplicity_method": row_value(
            pair, "endpoint_b_multiplicity_method"
        ),
        "four_state_minimum_zero_violation_depth": optional_int(
            pair.get("endpoint_b_minimum_zero_violation_depth"),
            field_name="endpoint_b_minimum_zero_violation_depth",
        ),
        "focal_ancestor_violation_p_exact": optional_float(
            pair.get("endpoint_b_focal_ancestor_violation_p_exact"),
            field_name="endpoint_b_focal_ancestor_violation_p_exact",
        ),
        "focal_ancestor_violation_upper_bound": optional_float(
            pair.get("endpoint_b_focal_ancestor_violation_upper_bound"),
            field_name="endpoint_b_focal_ancestor_violation_upper_bound",
        ),
        "focal_ancestor_violation_status": pair[
            "endpoint_b_focal_ancestor_violation_status"
        ],
        "partner_ancestor_violation_p_exact": optional_float(
            pair.get("endpoint_b_partner_ancestor_violation_p_exact"),
            field_name="endpoint_b_partner_ancestor_violation_p_exact",
        ),
        "partner_ancestor_violation_upper_bound": optional_float(
            pair.get("endpoint_b_partner_ancestor_violation_upper_bound"),
            field_name="endpoint_b_partner_ancestor_violation_upper_bound",
        ),
        "partner_ancestor_violation_status": pair[
            "endpoint_b_partner_ancestor_violation_status"
        ],
        "branching_violation_p_exact": optional_float(
            pair.get("endpoint_b_branching_violation_p_exact"),
            field_name="endpoint_b_branching_violation_p_exact",
        ),
        "branching_violation_upper_bound": optional_float(
            pair.get("endpoint_b_branching_violation_upper_bound"),
            field_name="endpoint_b_branching_violation_upper_bound",
        ),
        "branching_violation_status": pair["endpoint_b_branching_violation_status"],
        "four_state_relation_status": pair["endpoint_b_relation_compatibility"],
        "four_state_compatible_relation_models": optional_json(
            pair["endpoint_b_compatible_relation_models"],
            field_name="endpoint_b_compatible_relation_models",
        ),
        "four_state_n_compatible_relation_models": optional_int(
            pair["endpoint_b_n_compatible_relation_models"],
            field_name="endpoint_b_n_compatible_relation_models",
        ),
        "same_pair_four_state_witness": pair_is_four_state_witness(pair),
        "b1_prespecified_pair": pair_identity == b1_prespecified_pair_identity,
        "topology_scope": pair["topology_scope"],
        "topology_region": row_value(pair, "topology_region"),
        "topology_order_status": pair["topology_order_status"],
        "tumor_ref_status": tumor_status,
        "tumor_ref_read_count": optional_int(
            data.tumor_ref_rows[focal_key].get("n_tumor_ref"), field_name="n_tumor_ref"
        ),
        "normal_control_status": normal_status,
        "normal_called_reads": int(data.matched_normal_rows[focal_key]["normal_called_reads"]),
        "normal_ref_reads": int(data.matched_normal_rows[focal_key]["normal_ref_reads"]),
        "normal_alt_reads": int(data.matched_normal_rows[focal_key]["normal_alt_reads"]),
        "normal_unknown_reads": int(data.matched_normal_rows[focal_key]["normal_unknown_reads"]),
        "technical_exact_pair_present": parse_bool(
            pair["cross_platform_exact_pair_present"],
            field_name="cross_platform_exact_pair_present",
        ),
        "technical_replication_pass": parse_bool(
            pair["cross_platform_conditional_by_effect_direction_replication_pass"],
            field_name="cross_platform_conditional_by_effect_direction_replication_pass",
        ),
        "technical_replication_status": pair["cross_platform_replication_status"],
        "technical_core_read_name_overlap_present": optional_bool(
            pair.get("cross_platform_core_read_name_overlap_present"),
            field_name="cross_platform_core_read_name_overlap_present",
        ),
        "technical_molecule_independence_status": row_value(
            pair, "cross_platform_molecule_independence_status"
        ),
        "technical_biological_n": optional_int(
            pair.get("cross_platform_biological_n"),
            field_name="cross_platform_biological_n",
        ),
    }


def build_candidate_catalog_v2(
    data: IntegrationData,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    candidates: list[dict[str, Any]] = []
    pair_details: list[dict[str, Any]] = []
    candidate_keys = {
        key
        for key, row in data.cooccurrence_rows.items()
        if parse_bool(row[FORMAL_SELECTION_COLUMN], field_name=FORMAL_SELECTION_COLUMN)
    }
    for key in sorted(candidate_keys):
        screen = data.screen.stable_rows[key]
        cooccurrence = data.cooccurrence_rows[key]
        strict = data.strict_rows[key]
        tumor_ref = data.tumor_ref_rows[key]
        normal = data.matched_normal_rows[key]
        pair_rows = data.pair_data.by_focal.get(key, [])
        sample, chrom, pos, ref, alt = key
        r1_status, r1_reason = strict_r1_status(strict)
        tumor_status, tumor_reason = tumor_ref_claim_status(tumor_ref)
        normal_status, normal_reason = normal_claim_status(normal)
        formal_pairs = [
            row
            for row in pair_rows
            if parse_bool(
                row["endpoint_a_formal_pair_by_confirmed"],
                field_name="endpoint_a_formal_pair_by_confirmed",
            )
        ]
        compatible_opportunity_rows = [
            row for row in formal_pairs if pair_is_four_state_witness(row)
        ]
        compatible_opportunity_keys = [
            compact_json(list(pair_key(row, source="same-pair witness")))
            for row in compatible_opportunity_rows
        ]
        prespecified_pair = select_b1_prespecified_pair(formal_pairs)
        prespecified_pair_identity = (
            pair_key(prespecified_pair, source="B1 prespecified pair")
            if prespecified_pair is not None
            else None
        )
        prespecified_pair_key = (
            compact_json(list(prespecified_pair_identity))
            if prespecified_pair is not None
            else None
        )
        prespecified_four_state_status, prespecified_four_state_reason = (
            four_state_claim_decision(prespecified_pair)
            if prespecified_pair is not None
            else ("NOT_EVALUABLE", "NO_G1_FORMAL_PAIR_AVAILABLE_FOR_B1")
        )
        prespecified_witness_rows = (
            [prespecified_pair]
            if prespecified_pair is not None and prespecified_four_state_status == "PASS"
            else []
        )
        prespecified_witness_keys = (
            [prespecified_pair_key] if prespecified_witness_rows else []
        )
        if r1_status == "NOT_EVALUABLE":
            b1_status, b1_reason = "NOT_EVALUABLE", r1_reason
        elif r1_status == "FAIL":
            b1_status, b1_reason = "FAIL", r1_reason
        elif tumor_status == "NOT_EVALUABLE":
            b1_status, b1_reason = "NOT_EVALUABLE", tumor_reason
        elif tumor_status == "FAIL":
            b1_status, b1_reason = "FAIL", tumor_reason
        elif normal_status == "NOT_EVALUABLE":
            b1_status, b1_reason = "NOT_EVALUABLE", normal_reason
        elif normal_status == "FAIL":
            b1_status, b1_reason = "FAIL", normal_reason
        else:
            b1_status, b1_reason = prespecified_four_state_status, prespecified_four_state_reason
        cn_row = data.cn_ccf_rows.get(key)
        if b1_status == "PASS":
            c1_status, c1_reason = cn_ccf_claim_status(cn_row)
        else:
            c1_status, c1_reason = "NOT_RUN", "UPSTREAM_B1_NOT_PASS"
        top_confirmed_positions = optional_json(
            cooccurrence["top_marker_pair_by_confirmed_positions"],
            field_name="top_marker_pair_by_confirmed_positions",
        )
        if not isinstance(top_confirmed_positions, list):
            raise ContractError(f"Candidate top confirmed positions are not a list at {key}")
        top_spaced_positions = spatially_separated_positions(top_confirmed_positions)
        exact_technical = sum(
            parse_bool(
                row["cross_platform_exact_pair_present"],
                field_name="cross_platform_exact_pair_present",
            )
            for row in pair_rows
        )
        replicated_technical = sum(
            parse_bool(
                row["cross_platform_conditional_by_effect_direction_replication_pass"],
                field_name="cross_platform_conditional_by_effect_direction_replication_pass",
            )
            for row in pair_rows
        )
        technical_status = (
            "NOT_EVALUABLE_NON_HCC_TECHNICAL_PAIR"
            if sample not in {"HCC1395", "HCC1395_DORADO"}
            else "NOT_EVALUABLE_NO_EXACT_SHARED_PAIR"
            if exact_technical == 0
            else "ANY_CONCORDANT_EXACT_PAIR_OBSERVED"
            if replicated_technical > 0
            else "NO_CONCORDANT_EXACT_PAIR_OBSERVED"
        )
        row = {
            "candidate_key": compact_json(list(key)),
            "sample": sample,
            "biological_id": BIOLOGICAL_IDS[sample],
            "dataset_role": dataset_role(sample),
            "biological_n_contribution": biological_n_contribution(sample),
            "biological_site_key": compact_json(list(biological_key(key))),
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "truth_label": screen["truth_label"],
            "ssnv_branch": row_value(screen, "ssnv_branch"),
            "component_id": row_value(screen, "component_id"),
            "selected_group_id": row_value(screen, "selected_group_id"),
            "m1_status": "PASS",
            "m1_status_semantics": "OPERATIONAL_STABLE_MULTIGROUP_SCREEN_FLAG",
            "modal_assignment_ari_min": optional_float(
                screen["modal_assignment_ari_min"], field_name="modal_assignment_ari_min"
            ),
            "m2_status": "PASS",
            "m2_screen_gate_contract": cooccurrence["m2_screen_gate_contract"],
            "m2_screen_evaluable": parse_bool(
                cooccurrence["m2_screen_evaluable"], field_name="m2_screen_evaluable"
            ),
            "m2_screen_eligibility_status": cooccurrence[
                "m2_screen_eligibility_status"
            ],
            "m2_axis_statuses": optional_json(
                cooccurrence["m2_axis_statuses"], field_name="m2_axis_statuses"
            ),
            "m2_indeterminate_axes": optional_json(
                cooccurrence["m2_indeterminate_axes"], field_name="m2_indeterminate_axes"
            ),
            "m2_low_power_axes": optional_json(
                cooccurrence["m2_low_power_axes"], field_name="m2_low_power_axes"
            ),
            "hp_axis_confound": parse_bool(screen["hp_axis_confound"], field_name="hp_axis_confound"),
            "technical_axis_confound": parse_bool(
                screen["technical_axis_confound"], field_name="technical_axis_confound"
            ),
            "residual_unexplained_multigroup": parse_bool(
                screen["residual_unexplained_multigroup"],
                field_name="residual_unexplained_multigroup",
            ),
            "g1_status": "PASS",
            "n_pair_by_confirmed": int(cooccurrence["n_pair_by_confirmed"]),
            "pair_by_confirmed_positions": optional_json(
                cooccurrence["pair_by_confirmed_positions"],
                field_name="pair_by_confirmed_positions",
            ),
            "g2_status": "PASS",
            "n_spatially_separated_pair_by_20bp": int(
                cooccurrence["n_spatially_separated_pair_by_20bp"]
            ),
            "spatially_separated_pair_by_positions_20bp": optional_json(
                cooccurrence["spatially_separated_pair_by_positions_20bp"],
                field_name="spatially_separated_pair_by_positions_20bp",
            ),
            "top_marker_positions": optional_json(
                cooccurrence["top_marker_positions"], field_name="top_marker_positions"
            ),
            "n_top_marker_pair_by_confirmed": int(
                cooccurrence["n_top_marker_pair_by_confirmed"]
            ),
            "n_top_spatially_separated_pair_by_20bp": len(top_spaced_positions),
            "top_spatially_separated_pair_by_positions_20bp": top_spaced_positions,
            "joint_signature_n_complete_reads": int(
                cooccurrence["joint_signature_n_complete_reads"]
            ),
            "joint_signature_testable": parse_bool(
                cooccurrence["joint_signature_testable"], field_name="joint_signature_testable"
            ),
            "joint_signature_p_conditional_perm": optional_float(
                cooccurrence["joint_signature_p_conditional_perm"],
                field_name="joint_signature_p_conditional_perm",
            ),
            "joint_signature_permutations": optional_int(
                cooccurrence["joint_signature_permutations"],
                field_name="joint_signature_permutations",
            ),
            "joint_signature_permutable": parse_bool(
                cooccurrence["joint_signature_permutable"],
                field_name="joint_signature_permutable",
            ),
            "joint_signature_conditional_status": cooccurrence[
                "joint_signature_conditional_status"
            ],
            "joint_signature_sensitivity_pass": parse_bool(
                cooccurrence["joint_signature_sensitivity_pass"],
                field_name="joint_signature_sensitivity_pass",
            ),
            "joint_signature_global_fdr_family_status": cooccurrence[
                "joint_signature_global_fdr_family_status"
            ],
            "joint_signature_q_global_bh": optional_float(
                cooccurrence["joint_signature_q_global_bh"],
                field_name="joint_signature_q_global_bh",
            ),
            "joint_signature_q_global_by": optional_float(
                cooccurrence["joint_signature_q_global_by"],
                field_name="joint_signature_q_global_by",
            ),
            "joint_signature_global_bh_discovery": parse_bool(
                cooccurrence["joint_signature_global_bh_discovery"],
                field_name="joint_signature_global_bh_discovery",
            ),
            "joint_signature_global_by_discovery": parse_bool(
                cooccurrence["joint_signature_global_by_discovery"],
                field_name="joint_signature_global_by_discovery",
            ),
            "r1_status": r1_status,
            "r1_reason": r1_reason,
            "strict_analysis_class": strict["strict_analysis_class"],
            "strict_assignment_concordance_ari_min": optional_float(
                strict.get("strict_assignment_concordance_ari_min"),
                field_name="strict_assignment_concordance_ari_min",
            ),
            "tumor_ref_status": tumor_status,
            "tumor_ref_evaluable": parse_bool(tumor_ref["ref_evaluable"], field_name="ref_evaluable"),
            "tumor_ref_stable_null_multigroup": parse_bool(
                tumor_ref["ref_stable_null_multigroup"],
                field_name="ref_stable_null_multigroup",
            ),
            "tumor_ref_read_count": optional_int(
                tumor_ref.get("n_tumor_ref"), field_name="n_tumor_ref"
            ),
            "normal_control_status": normal_status,
            "normal_called_reads": int(normal["normal_called_reads"]),
            "normal_ref_reads": int(normal["normal_ref_reads"]),
            "normal_alt_reads": int(normal["normal_alt_reads"]),
            "normal_unknown_reads": int(normal["normal_unknown_reads"]),
            "normal_focal_callability_gate": parse_bool(
                normal["normal_focal_callability_gate"],
                field_name="normal_focal_callability_gate",
            ),
            "normal_ref_methyl_evaluable": parse_bool(
                normal["normal_evaluable"], field_name="normal_evaluable"
            ),
            "normal_ref_methyl_stable_multigroup": parse_bool(
                normal["normal_ref_methyl_stable_multigroup"],
                field_name="normal_ref_methyl_stable_multigroup",
            ),
            "normal_ref_methyl_nonreplication_gate": optional_bool(
                normal.get("normal_ref_methyl_nonreplication_gate"),
                field_name="normal_ref_methyl_nonreplication_gate",
            ),
            "n_same_pair_four_state_witnesses": len(prespecified_witness_rows),
            "same_pair_four_state_witness_keys": prespecified_witness_keys,
            "n_four_state_compatible_formal_pair_opportunities": len(
                compatible_opportunity_rows
            ),
            "four_state_compatible_formal_pair_opportunity_keys": (
                compatible_opportunity_keys
            ),
            "b1_prespecified_pair_key": prespecified_pair_key,
            "b1_prespecified_pair_relation_status": (
                prespecified_pair["endpoint_b_relation_compatibility"]
                if prespecified_pair is not None
                else None
            ),
            "b1_prespecified_pair_is_witness": bool(prespecified_witness_rows),
            "b1_pair_selection_contract": (
                "among_G1_formal_pairs_max_endpoint_a_n_informative_then_abs_distance_"
                "then_partner_identity_without_four_state_result"
            ),
            "b1_uses_posthoc_compatible_pair_search": False,
            "b1_status": b1_status,
            "b1_reason": b1_reason,
            "technical_replication_status": technical_status,
            "technical_exact_pair_opportunities": exact_technical,
            "technical_replicated_exact_pairs": replicated_technical,
            "technical_biological_n": 1 if sample in {"HCC1395", "HCC1395_DORADO"} else None,
            "cn_ccf_status": cn_row.get("cn_status") if cn_row else "N/A",
            "cn_ccf_source_path": cn_row.get("_cn_source_path") if cn_row else None,
            "cn_ccf_source_sha256": cn_row.get("cn_source_sha256") if cn_row else None,
            "cn_ccf_authority_reviewed": cn_row is not None,
            "cn_ccf_exact_locus_covered": (
                cn_row is not None
                and cn_row.get("cn_status")
                in {"AVAILABLE_EXACT_SEGMENT", "SHARED_CN_SENSITIVITY"}
            ),
            "cn_total": optional_float(
                cn_row.get("savana_total_cn_raw"), field_name="savana_total_cn_raw"
            ) if cn_row else None,
            "cn_minor": optional_float(
                cn_row.get("savana_minor_cn_raw"), field_name="savana_minor_cn_raw"
            ) if cn_row else None,
            "purity": optional_float(
                cn_row.get("purity_model_value"), field_name="purity_model_value"
            ) if cn_row else None,
            "mutation_multiplicity": None,
            "ccf": optional_float(
                cn_row.get("pyclone_vi_cellular_prevalence"),
                field_name="pyclone_vi_cellular_prevalence",
            ) if cn_row else None,
            "c1_status": c1_status,
            "c1_reason": c1_reason,
            "l1_status": "NOT_EVALUABLE" if c1_status == "PASS" else "NOT_RUN",
            "l2_status": "NOT_RUN",
            "deprecated_aliases": {
                DEFAULT_SELECTION_COLUMN: row_value(cooccurrence, DEFAULT_SELECTION_COLUMN),
                BY_SELECTION_COLUMN: row_value(cooccurrence, BY_SELECTION_COLUMN),
                "formal_claim": False,
            },
        }
        for field_name in CN_CCF_NATIVE_DETAIL_FIELDS:
            row[field_name] = cn_row.get(field_name) if cn_row else None
        for status_field in (
            "m1_status",
            "m2_status",
            "g1_status",
            "g2_status",
            "r1_status",
            "b1_status",
            "c1_status",
            "l1_status",
            "l2_status",
        ):
            validate_claim_status(row[status_field], field_name=status_field)
        candidates.append(row)
        pair_details.extend(
            build_pair_detail(
                data,
                key,
                pair,
                tumor_status,
                normal_status,
                prespecified_pair_identity,
            )
            for pair in pair_rows
        )
    if len({row["candidate_key"] for row in candidates}) != len(candidates):
        raise ContractError("Candidate catalog key uniqueness check failed")
    if len({row["witness_pair_key"] for row in pair_details}) != len(pair_details):
        raise ContractError("Candidate pair-detail key uniqueness check failed")
    return candidates, pair_details


def g2_evaluable(data: IntegrationData, key: CandidateKey) -> bool:
    cooccurrence = data.cooccurrence_rows[key]
    top_exact_identifiable_positions = [
        int(row["partner_pos"])
        for row in data.pair_data.by_focal.get(key, [])
        if parse_bool(row.get("top_coverage_marker"), field_name="top_coverage_marker")
        and parse_bool(row["endpoint_a_testable"], field_name="endpoint_a_testable")
        and row["endpoint_a_global_fdr_family_status"] == "ELIGIBLE_M2_EXACT_FAMILY"
        and optional_float(
            row["endpoint_a_p_fixed_margins_exact"],
            field_name="endpoint_a_p_fixed_margins_exact",
        )
        is not None
    ]
    joint_p = optional_float(
        cooccurrence["joint_signature_p_conditional_perm"],
        field_name="joint_signature_p_conditional_perm",
    )
    joint_permutations = optional_int(
        cooccurrence["joint_signature_permutations"],
        field_name="joint_signature_permutations",
    )
    return bool(
        spatially_separated_count(top_exact_identifiable_positions) >= 2
        and parse_bool(
            cooccurrence["joint_signature_testable"], field_name="joint_signature_testable"
        )
        and parse_bool(
            cooccurrence["joint_signature_permutable"],
            field_name="joint_signature_permutable",
        )
        and joint_p is not None
        and joint_permutations == PAIR_CONDITIONAL_PERMUTATIONS
    )


def build_site_claim_statuses(data: IntegrationData) -> dict[CandidateKey, dict[str, str]]:
    candidate_by_key = {
        candidate_key(row, source="candidate catalog"): row for row in data.candidates
    }
    statuses: dict[CandidateKey, dict[str, str]] = {}
    for key in data.screen.all_keys:
        screen_row = data.screen.stable_rows.get(key)
        stage: dict[str, str] = {}
        stable = key in data.screen.claim_dataset_sets["stable"]
        # M1 is an all-site operational flag yield. A site that cannot enter the
        # methyl screen is still "not flagged", not an excluded denominator row.
        # Its technical evaluability is preserved separately in
        # m1_operational_screen rather than converted into biological evidence.
        stage["M1"] = "PASS" if stable else "FAIL"
        cooccurrence = data.cooccurrence_rows.get(key)
        m2_eligible = bool(
            stable
            and cooccurrence is not None
            and parse_bool(
                cooccurrence["m2_screen_eligible"], field_name="m2_screen_eligible"
            )
        )
        m2_evaluable = bool(
            stable
            and cooccurrence is not None
            and parse_bool(
                cooccurrence["m2_screen_evaluable"], field_name="m2_screen_evaluable"
            )
        )
        stage["M2"] = claim_status(m2_eligible, evaluable=m2_evaluable) if stable else "NOT_RUN"
        if m2_eligible and screen_row is not None:
            pairs = data.pair_data.by_focal.get(key, [])
            exact_testable = any(
                optional_float(
                    row.get("endpoint_a_p_fixed_margins_exact"),
                    field_name="endpoint_a_p_fixed_margins_exact",
                )
                is not None
                and row["endpoint_a_global_fdr_family_status"]
                == "ELIGIBLE_M2_EXACT_FAMILY"
                for row in pairs
            )
            g1_pass = any(
                parse_bool(
                    row["endpoint_a_formal_pair_by_confirmed"],
                    field_name="endpoint_a_formal_pair_by_confirmed",
                )
                for row in pairs
            )
            stage["G1"] = claim_status(g1_pass, evaluable=exact_testable)
            stage["G2"] = claim_status(
                parse_bool(
                    data.cooccurrence_rows[key][FORMAL_SELECTION_COLUMN],
                    field_name=FORMAL_SELECTION_COLUMN,
                ),
                evaluable=g2_evaluable(data, key),
            )
        else:
            stage["G1"] = "NOT_RUN"
            stage["G2"] = "NOT_RUN"
        candidate = candidate_by_key.get(key)
        if candidate is None:
            stage.update({claim_id: "NOT_RUN" for claim_id in ("R1", "B1", "C1", "L1", "L2")})
        else:
            stage.update(
                {
                    "R1": candidate["r1_status"],
                    "B1": candidate["b1_status"],
                    "C1": candidate["c1_status"],
                    "L1": candidate["l1_status"],
                    "L2": candidate["l2_status"],
                }
            )
        for claim_id, status in stage.items():
            validate_claim_status(status, field_name=f"{claim_id}@{key}")
        statuses[key] = stage
    return statuses


def row_in_stratum_v2(
    key: CandidateKey,
    truth_label: str,
    stratum_type: str,
    stratum: str,
) -> bool:
    if stratum_type == "pooled":
        return True
    if stratum_type == "sample":
        return key[0] == stratum
    if stratum_type == "truth":
        return truth_label == stratum
    if stratum_type == "sample_truth":
        sample, truth = stratum.split("|", 1)
        return key[0] == sample and truth_label == truth
    raise ValueError(f"Unknown stratum type: {stratum_type}")


def build_metrics_v2(
    data: IntegrationData,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[CandidateKey, dict[str, str]]]:
    statuses = build_site_claim_statuses(data)
    strata = [("pooled", "ALL")]
    strata.extend(("sample", sample) for sample in DATASETS)
    strata.extend(("truth", truth) for truth in TRUTH_LABELS)
    strata.extend(
        ("sample_truth", f"{sample}|{truth}")
        for sample in DATASETS
        for truth in TRUTH_LABELS
    )
    nested: dict[str, Any] = {
        "pooled": {},
        "per_sample": {},
        "truth_strata": {},
        "sample_by_truth": {},
    }
    rows: list[dict[str, Any]] = []
    for stratum_type, stratum in strata:
        keys = [
            key
            for key, screen_row in data.screen.all_rows.items()
            if row_in_stratum_v2(
                key,
                screen_row["truth_label"],
                stratum_type,
                stratum,
            )
        ]
        metric_payload: dict[str, Any] = {}
        for claim_id in CLAIM_IDS:
            counter = Counter(statuses[key][claim_id] for key in keys)
            numerator = counter["PASS"]
            denominator = counter["PASS"] + counter["FAIL"]
            payload = {
                "numerator": numerator,
                "denominator": denominator,
                "ratio": rate(numerator, denominator),
                "not_evaluable": counter["NOT_EVALUABLE"],
                "not_run": counter["NOT_RUN"],
                "population": len(keys),
                "denominator_definition": CLAIM_RULES[claim_id]["denominator"],
            }
            metric_payload[claim_id] = payload
            sample_value = stratum if stratum_type == "sample" else None
            rows.append(
                {
                    "stratum_type": stratum_type,
                    "stratum": stratum,
                    "sample": sample_value,
                    "biological_id": BIOLOGICAL_IDS.get(sample_value) if sample_value else None,
                    "dataset_role": dataset_role(sample_value) if sample_value else None,
                    "biological_n_contribution": biological_n_contribution(sample_value) if sample_value else None,
                    "claim_id": claim_id,
                    **payload,
                }
            )
        if stratum_type == "pooled":
            nested["pooled"] = metric_payload
        elif stratum_type == "sample":
            nested["per_sample"][stratum] = metric_payload
        elif stratum_type == "truth":
            nested["truth_strata"][stratum] = metric_payload
        else:
            nested["sample_by_truth"][stratum] = metric_payload
    return nested, rows, statuses


def build_claim_ladder_v2(
    data: IntegrationData,
    statuses: Mapping[CandidateKey, Mapping[str, str]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for claim_id in CLAIM_IDS:
        counter = Counter(stage[claim_id] for stage in statuses.values())
        numerator = counter["PASS"]
        denominator = counter["PASS"] + counter["FAIL"]
        biological_primary_statuses: dict[BiologicalKey, str] = {}
        dorado_sensitivity_statuses: dict[BiologicalKey, str] = {}
        for key, stage in statuses.items():
            target = (
                dorado_sensitivity_statuses
                if key[0] == "HCC1395_DORADO"
                else biological_primary_statuses
            )
            bkey = biological_key(key)
            if bkey in target:
                raise ContractError(
                    f"Duplicate fixed-primary biological site for {claim_id}: {bkey}"
                )
            target[bkey] = stage[claim_id]
        biological_pass = sum(
            value == "PASS" for value in biological_primary_statuses.values()
        )
        biological_evaluable = sum(
            value in {"PASS", "FAIL"} for value in biological_primary_statuses.values()
        )
        hcc_primary = {
            key: value
            for key, value in biological_primary_statuses.items()
            if key[0] == "HCC1395"
        }
        shared_hcc_keys = set(hcc_primary).intersection(dorado_sensitivity_statuses)
        joint_evaluable_hcc_keys = {
            key
            for key in shared_hcc_keys
            if hcc_primary[key] in {"PASS", "FAIL"}
            and dorado_sensitivity_statuses[key] in {"PASS", "FAIL"}
        }
        hcc_pass_fail_concordant = sum(
            hcc_primary[key] == dorado_sensitivity_statuses[key]
            for key in joint_evaluable_hcc_keys
        )
        aggregate_status = (
            "PASS"
            if numerator > 0
            else "FAIL"
            if denominator > 0
            else "NOT_EVALUABLE"
            if counter["NOT_EVALUABLE"] > 0
            else "NOT_RUN"
        )
        rows.append(
            {
                "claim_id": claim_id,
                "claim_name": CLAIM_RULES[claim_id]["name"],
                "dataset_numerator": numerator,
                "dataset_denominator": denominator,
                "dataset_ratio": rate(numerator, denominator),
                "dataset_not_evaluable": counter["NOT_EVALUABLE"],
                "dataset_not_run": counter["NOT_RUN"],
                "biological_numerator": biological_pass,
                "biological_denominator": biological_evaluable,
                "biological_ratio": rate(biological_pass, biological_evaluable),
                "biological_primary_population": len(biological_primary_statuses),
                "biological_aggregation_policy": (
                    "fixed_primary_HCC1395_with_DORADO_as_technical_sensitivity_only"
                ),
                "hcc1395_primary_numerator": sum(
                    value == "PASS" for value in hcc_primary.values()
                ),
                "hcc1395_primary_denominator": sum(
                    value in {"PASS", "FAIL"} for value in hcc_primary.values()
                ),
                "dorado_sensitivity_numerator": sum(
                    value == "PASS" for value in dorado_sensitivity_statuses.values()
                ),
                "dorado_sensitivity_denominator": sum(
                    value in {"PASS", "FAIL"}
                    for value in dorado_sensitivity_statuses.values()
                ),
                "hcc1395_dorado_exact_site_overlap": len(shared_hcc_keys),
                "hcc1395_dorado_joint_evaluable": len(joint_evaluable_hcc_keys),
                "hcc1395_dorado_pass_fail_concordant": hcc_pass_fail_concordant,
                "hcc1395_dorado_pass_fail_discordant": (
                    len(joint_evaluable_hcc_keys) - hcc_pass_fail_concordant
                ),
                "status": aggregate_status,
                "denominator_definition": CLAIM_RULES[claim_id]["denominator"],
                "guardrail": CLAIM_RULES[claim_id]["guardrail"],
                "technical_replication_is_required_gate": False,
                "automatic_upgrade_prohibited": True,
            }
        )
    return rows


def input_artifacts(bundle: InputBundle, extra_paths: Iterable[Path]) -> dict[str, Any]:
    paths = {
        "manifest": bundle.manifest,
        "screen_sites": bundle.screen_sites,
        "screen_assignments": bundle.screen_assignments,
        "screen_summary": bundle.screen_summary,
        "screen_receipt": bundle.screen_receipt,
        "cooccurrence_sites": bundle.cooccurrence_sites,
        "cooccurrence_pairs": bundle.cooccurrence_pairs,
        "cooccurrence_summary": bundle.cooccurrence_summary,
        "cooccurrence_receipt": bundle.cooccurrence_receipt,
        "cooccurrence_release_receipt": bundle.cooccurrence_release_receipt,
        "cooccurrence_preflight": bundle.cooccurrence_preflight,
        "cooccurrence_raw_identity_duplicates": (
            bundle.cooccurrence_raw_identity_duplicates
        ),
        "cooccurrence_oracle_cases": bundle.cooccurrence_oracle_cases,
        "tumor_ref_sites": bundle.tumor_ref_sites,
        "tumor_ref_summary": bundle.tumor_ref_summary,
        "tumor_ref_receipt": bundle.tumor_ref_receipt,
        "primary_artifact_audit_pre": bundle.primary_artifact_audit_pre,
        "primary_artifact_audit_post": bundle.primary_artifact_audit_post,
    }
    if bundle.tumor_ref_source_identity_receipt is not None:
        paths["tumor_ref_source_identity_receipt"] = (
            bundle.tumor_ref_source_identity_receipt
        )
    if bundle.independent_m2_audit is not None:
        paths["independent_m2_audit"] = bundle.independent_m2_audit
    known = {path.resolve() for path in paths.values()}
    for index, path in enumerate(extra_paths, 1):
        if path.resolve() not in known:
            paths[f"downstream_extra_{index}"] = path
            known.add(path.resolve())
    return {name: artifact(path) for name, path in paths.items()}


def serialize_tsv(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> str:
    from io import StringIO

    handle = StringIO(newline="")
    writer = csv.DictWriter(
        handle,
        fieldnames=list(fields),
        delimiter="\t",
        extrasaction="raise",
        lineterminator="\n",
    )
    writer.writeheader()
    for row in rows:
        encoded = {}
        for field_name in fields:
            value = row.get(field_name)
            if value is None:
                encoded[field_name] = ""
            elif isinstance(value, bool):
                encoded[field_name] = "true" if value else "false"
            elif isinstance(value, (dict, list, tuple)):
                encoded[field_name] = compact_json(value)
            else:
                encoded[field_name] = value
        writer.writerow(encoded)
    return handle.getvalue()


def assert_output_available(path: Path) -> None:
    if os.path.lexists(path):
        raise FileExistsError(f"Refusing to overwrite existing output directory: {path}")


def downstream_primary_consumer_receipts(
    bundle: InputBundle,
    strict_status: Mapping[str, Any],
    matched_normal_status: Mapping[str, Any],
) -> list[Path]:
    paths = [
        bundle.cooccurrence_preflight,
        bundle.cooccurrence_receipt,
        bundle.cooccurrence_release_receipt,
        bundle.tumor_ref_receipt,
    ]
    if bundle.independent_m2_audit is not None:
        paths.append(bundle.independent_m2_audit)
    strict_path = strict_status.get("receipt_path")
    if not isinstance(strict_path, Path):
        raise ContractError("Strict consumer receipt path is missing")
    paths.append(strict_path)
    paired_path = matched_normal_status.get("paired_receipt_path")
    if paired_path is not None:
        if not isinstance(paired_path, Path):
            raise ContractError("Matched-normal paired receipt path is invalid")
        paths.append(paired_path)
    analysis_path = matched_normal_status.get("receipt_path")
    if analysis_path is not None:
        if not isinstance(analysis_path, Path):
            raise ContractError("Matched-normal analysis receipt path is invalid")
        paths.append(analysis_path)
    if len({path.resolve() for path in paths}) != len(paths):
        raise ContractError("Primary consumer receipt list contains duplicates")
    return paths


def downstream_primary_consumer_payloads(
    cooccurrence_preflight: Mapping[str, Any],
    cooccurrence_receipt: Mapping[str, Any],
    tumor_ref_receipt: Mapping[str, Any],
    strict_status: Mapping[str, Any],
    matched_normal_status: Mapping[str, Any],
) -> dict[str, Mapping[str, Any]]:
    payloads: dict[str, Mapping[str, Any]] = {
        "cooccurrence-preflight": cooccurrence_preflight,
        "cooccurrence": cooccurrence_receipt,
        "tumor-REF": tumor_ref_receipt,
    }
    strict_receipt = strict_status.get("receipt")
    if not isinstance(strict_receipt, Mapping):
        raise ContractError("Strict consumer receipt payload is missing")
    payloads["strict"] = strict_receipt
    paired = matched_normal_status.get("paired_receipt")
    if paired is not None:
        if not isinstance(paired, Mapping):
            raise ContractError("Matched-normal paired receipt payload is invalid")
        payloads["matched-normal-paired"] = paired
    analysis = matched_normal_status.get("receipt")
    if analysis is not None and analysis is not paired:
        if not isinstance(analysis, Mapping):
            raise ContractError("Matched-normal analysis receipt payload is invalid")
        payloads["matched-normal-analysis"] = analysis
    return payloads


def integrate_inputs(
    bundle: InputBundle,
    *,
    expected_screen_sites: int = EXPECTED_SCREEN_SITES,
) -> tuple[IntegrationData, dict[str, Any], list[Path]]:
    release_source_authority = current_release_source_authority()
    manifest = load_manifest(bundle.manifest, expected_screen_sites)
    screen = load_screen(bundle.screen_sites, manifest, expected_screen_sites)
    assignment_keys = load_assignment_keys(bundle.screen_assignments, screen)
    primary_artifact_audit_pre = load_primary_artifact_audit(
        bundle.primary_artifact_audit_pre,
        bundle,
        assignment_keys,
        label="pre-downstream primary artifact audit",
        expected_consumer_receipts=(),
    )
    screen_summary, screen_receipt, screen_recovery_source_validation = validate_screen_receipts(
        bundle, manifest, screen, expected_screen_sites
    )
    cooccurrence_rows, pre_m2_formal_pair_keys, strict_candidate_keys = load_cooccurrence_sites(
        bundle.cooccurrence_sites, screen
    )
    if set(cooccurrence_rows) != assignment_keys:
        raise ContractError("Cooccurrence stable sites do not exactly match stable assignments")
    pair_data = load_pair_data(
        bundle.cooccurrence_pairs, assignment_keys, cooccurrence_rows
    )
    (
        cooccurrence_summary,
        cooccurrence_receipt,
        cooccurrence_preflight,
    ) = validate_cooccurrence_receipts(bundle, cooccurrence_rows, pair_data)
    independent_m2_audit, independent_m2_paths = load_independent_m2_audit(
        bundle.independent_m2_audit,
        bundle,
        screen,
        cooccurrence_rows,
    )
    strict_rows, strict_status, strict_paths = load_strict(bundle, strict_candidate_keys)
    tumor_ref_rows, tumor_ref_summary, tumor_ref_receipt = load_tumor_ref(
        bundle, assignment_keys
    )
    tumor_ref_source_identity, source_identity_paths = (
        load_tumor_ref_source_identity_attestation(
            bundle.tumor_ref_source_identity_receipt,
            bundle.tumor_ref_receipt,
            tumor_ref_receipt,
        )
    )
    if not strict_candidate_keys.issubset(tumor_ref_rows):
        raise ContractError("Tumor-REF candidate join coverage is not exact")
    matched_normal_rows, matched_normal_status, normal_paths = load_matched_normal(
        bundle, strict_candidate_keys
    )
    cn_ccf_rows, cn_ccf_status, cn_ccf_paths = load_cn_ccf_annotations(
        bundle.cn_ccf_annotations, strict_candidate_keys
    )
    primary_consumer_receipts = downstream_primary_consumer_receipts(
        bundle, strict_status, matched_normal_status
    )
    primary_artifact_audit_post = load_primary_artifact_audit(
        bundle.primary_artifact_audit_post,
        bundle,
        assignment_keys,
        label="post-downstream primary artifact audit",
        expected_consumer_receipts=primary_consumer_receipts,
    )
    pre_artifact_set = primary_artifact_audit_pre["verification"]["artifact_set_sha256"]
    post_artifact_set = primary_artifact_audit_post["verification"]["artifact_set_sha256"]
    if pre_artifact_set != post_artifact_set:
        raise ContractError("Stable primary artifact set changed across downstream analyses")
    validate_primary_artifact_audit_window(
        primary_artifact_audit_pre,
        primary_artifact_audit_post,
        downstream_primary_consumer_payloads(
            cooccurrence_preflight,
            cooccurrence_receipt,
            tumor_ref_receipt,
            strict_status,
            matched_normal_status,
        ),
    )
    data = IntegrationData(
        manifest=manifest,
        screen=screen,
        cooccurrence_rows=cooccurrence_rows,
        pair_data=pair_data,
        strict_rows=strict_rows,
        strict_status=strict_status,
        tumor_ref_rows=tumor_ref_rows,
        matched_normal_rows=matched_normal_rows,
        matched_normal_status=matched_normal_status,
        cn_ccf_rows=cn_ccf_rows,
        cn_ccf_status=cn_ccf_status,
    )
    if strict_candidate_keys:
        data.candidates, data.candidate_pair_details = build_candidate_catalog_v2(data)
    validations = {
        "task_type_b_full_scope": True,
        "release_source_authority": release_source_authority,
        "screen_expected_sites": expected_screen_sites,
        "screen_observed_sites": len(screen.all_keys),
        "screen_dataset_count": len(DATASETS),
        "screen_biological_sample_count": len(set(BIOLOGICAL_IDS.values())),
        "truth_tp_fp_unassessed_equals_total": True,
        "stable_screen_assignment_keys_exact": True,
        "stable_primary_artifacts_pre_downstream_verified": True,
        "stable_primary_artifacts_post_downstream_verified": True,
        "stable_primary_artifacts_unchanged_across_downstream": True,
        "stable_primary_artifact_audit_window_contains_data_plane_consumers": True,
        "screen_recovery_source_validation": screen_recovery_source_validation,
        "stable_primary_artifact_set_sha256": pre_artifact_set,
        "cooccurrence_stable_keys_exact": True,
        "independent_m2_audit": independent_m2_audit,
        "pre_m2_sites_with_at_least_one_formal_pair_count": len(pre_m2_formal_pair_keys),
        "g2_candidate_count": len(strict_candidate_keys),
        "candidate_keys_unique": len(strict_candidate_keys) == len(data.candidates),
        "strict_candidate_keys_exact": set(strict_rows) == strict_candidate_keys,
        "strict_selection_column": STRICT_SELECTION_COLUMN,
        "strict_legacy_inference_fields_are_formal": False,
        "strict_postselection_bh_by_calibrated": False,
        "strict_zero_candidate_status": strict_status["status"],
        "tumor_ref_stable_keys_exact": set(tumor_ref_rows) == assignment_keys,
        "tumor_ref_candidate_join_exact": strict_candidate_keys.issubset(tumor_ref_rows),
        "tumor_ref_source_identity_attestation": tumor_ref_source_identity,
        "matched_normal_candidate_keys_exact": (
            set(matched_normal_rows) == strict_candidate_keys
        ),
        "matched_normal_status": matched_normal_status["status"],
        "many_to_many_candidate_joins": False,
        "technical_replicate_counted_as_biological_n": False,
        "technical_replication_is_b1_gate": False,
        "same_pair_witness_gate": all(
            row["n_same_pair_four_state_witnesses"]
            == len(row["same_pair_four_state_witness_keys"])
            for row in data.candidates
        ),
        "b1_single_prespecified_pair_gate": all(
            row["n_same_pair_four_state_witnesses"] in {0, 1}
            and row["b1_prespecified_pair_is_witness"]
            == (row["n_same_pair_four_state_witnesses"] == 1)
            and row["b1_uses_posthoc_compatible_pair_search"] is False
            for row in data.candidates
        ),
        "cn_ccf_default_diploid_imputation_allowed": False,
        "claim_automatic_upgrade": False,
        "producer_receipts": {
            "screen": screen_receipt["schema_name"],
            "cooccurrence": cooccurrence_receipt["schema_name"],
            "tumor_ref": tumor_ref_receipt["schema_name"],
        },
        "producer_summaries": {
            "screen": screen_summary["schema_name"],
            "cooccurrence": cooccurrence_summary["schema_name"],
            "tumor_ref": tumor_ref_summary["schema_name"],
        },
    }
    return data, validations, [
        *strict_paths,
        *normal_paths,
        *cn_ccf_paths,
        *source_identity_paths,
        *independent_m2_paths,
        SOURCE_AUTHORITY.AUTHORITY_PATH,
        SOURCE_AUTHORITY.CLAIM_CONTRACT_PATH,
    ]


def build_outputs(
    bundle: InputBundle,
    output_dir: Path,
    *,
    expected_screen_sites: int = EXPECTED_SCREEN_SITES,
    command: Sequence[str] | None = None,
) -> dict[str, Path]:
    output_dir = output_dir.expanduser().resolve()
    production_release = expected_screen_sites == EXPECTED_SCREEN_SITES
    if production_release:
        validate_canonical_task_b_paths(bundle, output_dir)
        expected_command = canonical_task_b_final_builder_command(bundle, output_dir)
        if list(command or []) != expected_command:
            raise ContractError("Task-B final dataset builder command is not canonical")
        if observed_process_command() != expected_command:
            raise ContractError(
                "Task-B final dataset builder was not executed as the direct canonical process"
            )
    assert_output_available(output_dir)
    data, validations, downstream_paths = integrate_inputs(
        bundle, expected_screen_sites=expected_screen_sites
    )
    validations["formal_task_b_release_eligible"] = production_release
    if production_release:
        validations["strict_statistics_deterministic_replay"] = (
            replay_strict_statistics(bundle, data.strict_rows, data.strict_status)
        )
        validations["independent_primary_artifact_recount"] = (
            independently_recount_primary_artifacts(
                bundle.screen_assignments,
                set(data.screen.stable_rows),
                str(validations["stable_primary_artifact_set_sha256"]),
            )
        )
    else:
        validations["task_type_b_full_scope"] = False
        validations["strict_statistics_deterministic_replay"] = {
            "status": "NOT_RUN_TEST_FIXTURE_ONLY",
            "pass": None,
        }
        validations["independent_primary_artifact_recount"] = {
            "status": "NOT_RUN_TEST_FIXTURE_ONLY",
            "pass": None,
        }
    metrics, metric_rows, site_statuses = build_metrics_v2(data)
    claim_rows = build_claim_ladder_v2(data, site_statuses)
    inputs = input_artifacts(bundle, downstream_paths)
    hcc_pair_union = {
        cross_platform_key(key)
        for key in data.pair_data.rows
        if key[0] in {"HCC1395", "HCC1395_DORADO"}
    }
    technical_denominator = len(data.pair_data.exact_cross_platform_pairs)
    technical_numerator = len(data.pair_data.replicated_cross_platform_pairs)
    technical_not_evaluable = len(
        hcc_pair_union - data.pair_data.exact_cross_platform_pairs
    )
    focal_partner_truth_counts = Counter(
        (str(row["focal_truth_label"]), str(row["partner_truth_label"]))
        for row in data.pair_data.rows.values()
    )
    focal_partner_formal_truth_counts = Counter(
        (str(row["focal_truth_label"]), str(row["partner_truth_label"]))
        for row in data.pair_data.rows.values()
        if parse_bool(
            row["endpoint_a_formal_pair_by_confirmed"],
            field_name="endpoint_a_formal_pair_by_confirmed",
        )
    )
    focal_partner_truth_matrix = [
        {
            "focal_truth_label": focal_truth,
            "partner_truth_label": partner_truth,
            "n_all_pair_rows": focal_partner_truth_counts[(focal_truth, partner_truth)],
            "n_g1_formal_pair_rows": focal_partner_formal_truth_counts[
                (focal_truth, partner_truth)
            ],
        }
        for focal_truth in TRUTH_LABELS
        for partner_truth in TRUTH_LABELS
    ]
    technical_status = (
        "ANY_CONCORDANT_EXACT_PAIR_OBSERVED"
        if technical_numerator > 0
        else "NO_CONCORDANT_EXACT_PAIR_OBSERVED"
        if technical_denominator > 0
        else "NOT_EVALUABLE"
    )
    claim_counts = {
        claim_id: Counter(stage[claim_id] for stage in site_statuses.values())
        for claim_id in CLAIM_IDS
    }
    m2_not_evaluable_reason_counts: Counter[str] = Counter()
    m2_group_limit_examples: list[dict[str, Any]] = []
    for key, row in data.cooccurrence_rows.items():
        if site_statuses[key]["M2"] != "NOT_EVALUABLE":
            continue
        status = str(row["m2_screen_eligibility_status"])
        reason = status.split(":", 1)[0]
        m2_not_evaluable_reason_counts[reason] += 1
        if reason != "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM":
            continue
        screen_row = data.screen.stable_rows[key]
        observed_groups = len(
            M2_GATE.parse_cluster_sizes(screen_row["cluster_sizes"], max_groups=None)
        )
        m2_group_limit_examples.append(
            {
                "dataset": key[0],
                "chrom": key[1],
                "pos": key[2],
                "ref": key[3],
                "alt": key[4],
                "observed_methyl_groups": observed_groups,
                "maximum_supported_methyl_groups": M2_GATE.MAX_METHYL_GROUPS,
                "m1_status": site_statuses[key]["M1"],
                "m2_status": site_statuses[key]["M2"],
                "g1_status": site_statuses[key]["G1"],
                "g2_status": site_statuses[key]["G2"],
                "b1_status": site_statuses[key]["B1"],
                "reason": status,
            }
        )
    m2_group_limit_examples.sort(
        key=lambda row: (
            row["dataset"],
            int(str(row["chrom"]).removeprefix("chr")),
            int(row["pos"]),
            row["ref"],
            row["alt"],
        )
    )
    m2_not_evaluable_total = claim_counts["M2"]["NOT_EVALUABLE"]
    if sum(m2_not_evaluable_reason_counts.values()) != m2_not_evaluable_total:
        raise RuntimeError(
            "M2 NOT_EVALUABLE reason census does not reconcile: "
            f"reasons={sum(m2_not_evaluable_reason_counts.values())} "
            f"metric={m2_not_evaluable_total}"
        )
    final_dataset = {
        "schema_name": "intersubmod.all_ssnv_final_report_dataset",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": now_utc(),
        "task_type": (
            "B_comprehensive_validation"
            if production_release
            else "NON_RELEASE_TEST_FIXTURE_VALIDATION"
        ),
        "formal_task_b_release_eligible": production_release,
        "pass_semantics": PASS_SEMANTICS,
        "scope": {
            "datasets": list(DATASETS),
            "dataset_count": len(DATASETS),
            "biological_samples": sorted(set(BIOLOGICAL_IDS.values())),
            "biological_sample_count": len(set(BIOLOGICAL_IDS.values())),
            "expected_screen_sites": expected_screen_sites,
            "observed_screen_sites": len(data.screen.all_keys),
            "technical_replicate": {
                "primary": "HCC1395",
                "replicate": "HCC1395_DORADO",
                "biological_id": "HCC1395",
                "counts_as_independent_biological_n": False,
            },
        },
        "denominator_contract": {
            "tables_may_have_different_row_counts": True,
            "every_ratio_has_numerator_and_denominator": True,
            "every_claim_metric_has_not_evaluable_and_not_run": True,
            "site_weighted_dataset_rows_are_not_independent_biological_replicates": True,
            "missing_values_are_preserved_as_json_null_and_blank_tsv_cells": True,
            "truth_is_posthoc_stratification_only": True,
            "site_metric_truth_axis": "focal_site_only",
            "pair_truth_axes": "focal_x_partner_reported_separately",
        },
        "m1_operational_screen": {
            "status_semantics": "FLAGGED_VS_NOT_FLAGGED_OPERATIONAL_SCREEN_ONLY",
            "denominator_definition": CLAIM_RULES["M1"]["denominator"],
            "n_all_dataset_sites": len(data.screen.all_keys),
            "n_screen_evaluable": len(data.screen.claim_dataset_sets["evaluable"]),
            "n_screen_not_evaluable": (
                len(data.screen.all_keys)
                - len(data.screen.claim_dataset_sets["evaluable"])
            ),
            "n_flagged_stable_multigroup": len(data.screen.claim_dataset_sets["stable"]),
            "n_not_flagged_all": (
                len(data.screen.all_keys)
                - len(data.screen.claim_dataset_sets["stable"])
            ),
            "flag_yield": rate(
                len(data.screen.claim_dataset_sets["stable"]),
                len(data.screen.all_keys),
            ),
            "flag_yield_among_screen_evaluable": rate(
                len(data.screen.claim_dataset_sets["stable"]),
                len(data.screen.claim_dataset_sets["evaluable"]),
            ),
            "n_evaluable_not_flagged": (
                len(data.screen.claim_dataset_sets["evaluable"])
                - len(data.screen.claim_dataset_sets["stable"])
            ),
            "global_null_validity_exported_for_nonstable_sites": False,
            "nonflagged_scientific_interpretation": (
                "NOT_IDENTIFIABLE_AS_TRUE_NEGATIVE_VS_NULL_INVALID_FROM_SITE_TSV"
            ),
            "biological_prevalence_estimate": None,
        },
        "m2_evaluability_contract": {
            "gate_contract": M2_GATE.GATE_CONTRACT,
            "denominator_definition": CLAIM_RULES["M2"]["denominator"],
            "minimum_supported_methyl_groups": 2,
            "maximum_supported_methyl_groups": M2_GATE.MAX_METHYL_GROUPS,
            "categorical_planning_level_ceilings": dict(
                M2_GATE.CATEGORICAL_LEVEL_CEILINGS
            ),
            "assignment_observed_levels_role": (
                "constant-axis proof only; observed levels do not replace the "
                "source-locked planning level ceilings"
            ),
            "n_m1_pass": claim_counts["M1"]["PASS"],
            "n_m2_evaluable": (
                claim_counts["M2"]["PASS"] + claim_counts["M2"]["FAIL"]
            ),
            "n_m2_not_evaluable": m2_not_evaluable_total,
            "not_evaluable_reason_counts": dict(
                sorted(m2_not_evaluable_reason_counts.items())
            ),
            "n_group_count_exceeds_planning_model_maximum": len(
                m2_group_limit_examples
            ),
            "group_count_exceeds_planning_model_examples": m2_group_limit_examples[:50],
            "group_count_exceeds_examples_complete": len(m2_group_limit_examples) <= 50,
            "group_count_exceeds_claim_behavior": {
                "M1": "PASS retained",
                "M2": "NOT_EVALUABLE excluded from PASS/FAIL denominator",
                "G1": "NOT_RUN",
                "G2": "NOT_RUN",
                "B1": "NOT_RUN",
            },
            "independent_logic_audit": validations["independent_m2_audit"],
        },
        "tumor_ref_source_identity_attestation": validations[
            "tumor_ref_source_identity_attestation"
        ],
        "background_control_replication_gate": {
            "contract": BACKGROUND_CONTROL_REPLICATION_GATE_CONTRACT,
            "applies_to": ["tumor_REF", "matched_normal_REF"],
            "required_conditions": [
                "coarse_ng>=2",
                "modal_fraction>=0.7_via_unstable_false",
            ],
            "membership_ari_minimum_required": False,
            "relation_to_primary_m1_replication_flags": (
                BACKGROUND_CONTROL_RELATION_TO_PRIMARY_M1
            ),
            "b1_pass_direction": "requires_no_lenient_background_replication",
            "false_positive_direction": (
                "cannot_increase_B1_passes_vs_ARI_qualified_predicate_on_same_background_payload"
            ),
            "false_negative_direction": (
                "may_conservatively_reduce_B1_passes_when_K_is_stable_but_membership_is_not"
            ),
            "scientific_interpretation": (
                "background nonreplication guardrail, not an exact primary-M1 replay"
            ),
        },
        "m2_axis_statistic_provenance": {
            "axis_effect_and_permutation_p_source": (
                "source_locked_focal_alt_multigroup_screen_producer"
            ),
            "screen_recovery_source_validation": validations[
                "screen_recovery_source_validation"
            ],
            "downstream_raw_read_axis_statistic_recomputed": False,
            "downstream_axis_classification_recomputed": True,
            "downstream_checks": [
                "axis sample-size reconciliation",
                "499-permutation add-one p-value grid",
                "effect threshold classification",
                "80-percent planning-power evaluability",
                "assignment-derived categorical constant-axis proof",
                "asymmetric positive-confound versus negative-evaluability power decision",
            ],
            "categorical_planning_level_ceilings": dict(
                M2_GATE.CATEGORICAL_LEVEL_CEILINGS
            ),
            "assignment_observed_levels_role": (
                "constant-axis proof only; observed levels do not replace the "
                "source-locked planning level ceilings"
            ),
            "independent_logic_audit": validations["independent_m2_audit"],
            "claim_guardrail": (
                "M2 effect/p are producer-derived; terminal validation independently reclassifies "
                "the frozen statistics but is not an independent raw-read remeasurement."
            ),
        },
        "focal_partner_truth_matrix": focal_partner_truth_matrix,
        "validations": validations,
        "input_artifacts": inputs,
        "provenance_chain": {
            "screen_receipt_to_manifest_and_screen_outputs": "VERIFIED_SHA256",
            "cooccurrence_receipt_to_manifest_screen_assignments_and_sites": "VERIFIED_SHA256",
            "strict_receipt_to_cooccurrence_receipt_candidate_table_and_assignments": (
                "NOT_APPLICABLE_ZERO_G2" if not data.candidates else "VERIFIED_SHA256"
            ),
            "tumor_ref_receipt_to_screen_sites_and_assignments": "VERIFIED_SHA256",
            "tumor_ref_bounded_source_identity_release_gate": (
                "VERIFIED_SHA256"
                if validations["tumor_ref_source_identity_attestation"][
                    "release_gate_pass"
                ]
                else "NOT_INCLUDED_INTERMEDIATE_TERMINAL_BUILD"
            ),
            "normal_receipt_to_primary_assignments_and_paired_run": (
                "NOT_APPLICABLE_ZERO_G2" if not data.candidates else "VERIFIED_SHA256"
            ),
            "hash_mismatch_policy": "HARD_FAIL_BEFORE_OUTPUT",
        },
        "funnel_metrics": metrics,
        "technical_replication": {
            "endpoint": "HCC1395_vs_HCC1395_DORADO_exact_pair_technical_concordance",
            "status": technical_status,
            "numerator": technical_numerator,
            "denominator": technical_denominator,
            "ratio": rate(technical_numerator, technical_denominator),
            "not_evaluable_one_platform_only": technical_not_evaluable,
            "denominator_definition": (
                "exact shared focal+partner CHROM/POS/REF/ALT opportunities in both HCC1395 pipelines"
            ),
            "biological_n": 1,
            "independent_biological_replication_n": 0,
            "replication_claim_status": "NOT_EVALUABLE_BIOLOGICAL_N1",
            "inferential_confidence_interval": None,
            "pair_independence_assumption_met": False,
            "required_for_b1": False,
            "same_exact_pair_evidence_required": True,
        },
        "biological_replication": {
            "status": "NOT_RUN",
            "independent_biological_replication_n": 0,
            "hcc1395_and_dorado_biological_n": 1,
            "guardrail": "A technical reprocessing pair is not independent biological replication.",
        },
        "claim_ladder": claim_rows,
        "candidate_catalog": data.candidates,
        "candidate_witness_pairs": data.candidate_pair_details,
        "cn_ccf_annotation": data.cn_ccf_status,
        "counts": {
            "screen_sites": len(data.screen.all_keys),
            "stable_sites": len(data.screen.stable_rows),
            "cooccurrence_pair_rows": data.pair_data.n_rows,
            "pre_m2_sites_with_at_least_one_formal_pair": validations[
                "pre_m2_sites_with_at_least_one_formal_pair_count"
            ],
            "g2_candidates": validations["g2_candidate_count"],
            "candidate_catalog_rows": len(data.candidates),
            "candidate_witness_pair_rows": len(data.candidate_pair_details),
            "same_pair_four_state_witnesses": sum(
                row["n_same_pair_four_state_witnesses"] for row in data.candidates
            ),
            "claim_status_counts": {
                claim_id: {status: claim_counts[claim_id][status] for status in sorted(CLAIM_STATUSES)}
                for claim_id in CLAIM_IDS
            },
        },
        "guardrails": {
            "automatic_claim_upgrade_prohibited": True,
            "strict_postselection_bh_by_are_uncalibrated_descriptions": True,
            "formal_pair_global_fdr_uses_fixed_margins_exact_p": True,
            "conditional_permutation_is_secondary_sensitivity": True,
            "spatially_separated_markers_are_not_statistically_independent": True,
            "same_pair_witness_required": True,
            "technical_replication_is_separate_from_b1": True,
            "cn_two_default_imputation_prohibited": True,
            "l1_requires_orthogonal_cellular_identity": True,
            "l2_requires_identifiable_order_evidence": True,
        },
        "pass": True,
    }
    metric_fields = (
        "stratum_type",
        "stratum",
        "sample",
        "biological_id",
        "dataset_role",
        "biological_n_contribution",
        "claim_id",
        "numerator",
        "denominator",
        "ratio",
        "not_evaluable",
        "not_run",
        "population",
        "denominator_definition",
    )
    candidate_fields = CANDIDATE_FIELDS
    claim_fields = tuple(claim_rows[0])
    serialized = {
        "final_report_dataset.json": json.dumps(
            final_dataset, ensure_ascii=False, indent=2, sort_keys=True
        )
        + "\n",
        "per_sample_metrics.tsv": serialize_tsv(metric_rows, metric_fields),
        "candidate_catalog.tsv": serialize_tsv(data.candidates, candidate_fields),
        "candidate_witness_pairs.tsv": serialize_tsv(
            data.candidate_pair_details, PAIR_DETAIL_FIELDS
        ),
        "claim_ladder.tsv": serialize_tsv(claim_rows, claim_fields),
    }
    assert_output_available(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=False, exist_ok=False)
    output_paths: dict[str, Path] = {}
    for filename, content in serialized.items():
        path = output_dir / filename
        path.write_text(content, encoding="utf-8")
        path.chmod(0o444)
        output_paths[filename] = path
    receipt_path = output_dir / "run_receipt.json"
    receipt = {
        "schema_name": "intersubmod.all_ssnv_final_report_dataset_run_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": now_utc(),
        "task_type": (
            "B_comprehensive_validation"
            if production_release
            else "NON_RELEASE_TEST_FIXTURE_VALIDATION"
        ),
        "formal_task_b_release_eligible": production_release,
        "pass_semantics": PASS_SEMANTICS,
        "command": list(command or []),
        "source_authority": validations["release_source_authority"],
        "inputs": inputs,
        "code": {
            "final_report_dataset_builder": artifact(Path(__file__).resolve()),
        },
        "source_modes": {
            "final_report_dataset_builder": oct(
                Path(__file__).resolve().stat().st_mode & 0o777
            ),
        },
        "validations": validations,
        "contracts": {
            "output_overwrite_allowed": False,
            "all_candidate_joins": "exact one-to-one on sample,chrom,pos,ref,alt",
            "different_source_row_counts_allowed": True,
            "ratios_include_explicit_numerator_denominator": True,
            "claim_metrics_include_not_evaluable_and_not_run": True,
            "technical_replicate_is_not_biological_n": True,
            "technical_replication_is_b1_gate": False,
            "same_pair_witness_required": True,
            "cn_two_default_imputation_allowed": False,
            "claim_automatic_upgrade": False,
            "report_text_or_figures_generated": False,
            "tumor_ref_source_identity_release_gate": (
                "PASS"
                if validations["tumor_ref_source_identity_attestation"][
                    "release_gate_pass"
                ]
                else "NOT_INCLUDED_INTERMEDIATE_ONLY"
            ),
        },
        "outputs": {
            name.removesuffix(".json").removesuffix(".tsv"): artifact(path)
            for name, path in output_paths.items()
        },
        "counts": final_dataset["counts"],
        "pass": True,
    }
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    receipt_path.chmod(0o444)
    output_paths["run_receipt.json"] = receipt_path
    return output_paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", "--input-manifest", dest="manifest", type=Path, required=True)
    parser.add_argument(
        "--screen-dir", "--screen-output-dir", dest="screen_dir", type=Path, required=True
    )
    parser.add_argument(
        "--cooccurrence-dir",
        "--cooccurrence-output-dir",
        dest="cooccurrence_dir",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--strict-dir", "--strict-output-dir", dest="strict_dir", type=Path, required=True
    )
    parser.add_argument(
        "--tumor-ref-dir",
        "--tumor-ref-output-dir",
        dest="tumor_ref_dir",
        type=Path,
        required=True,
    )
    parser.add_argument("--primary-artifact-audit-pre", type=Path, required=True)
    parser.add_argument("--primary-artifact-audit-post", type=Path, required=True)
    parser.add_argument("--cooccurrence-preflight", type=Path, required=True)
    parser.add_argument(
        "--tumor-ref-source-identity-receipt",
        type=Path,
        default=None,
        help=(
            "Passing bounded retrospective source identity receipt; required for the final "
            "Task Type B release, optional only for the intermediate terminal build"
        ),
    )
    parser.add_argument(
        "--independent-m2-audit",
        type=Path,
        default=None,
        help="Logic-independent M2 v2 recount receipt for final provenance and count reconciliation",
    )
    parser.add_argument(
        "--matched-normal-dir",
        "--matched-normal-output-dir",
        "--normal-dir",
        dest="matched_normal_dir",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--cn-ccf-annotations",
        type=Path,
        default=None,
        help=(
            "Optional native annotate_candidate_cn_ccf.py output directory (or its TSV.GZ)"
        ),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    bundle = InputBundle(
        manifest=args.manifest,
        screen_dir=args.screen_dir,
        cooccurrence_dir=args.cooccurrence_dir,
        strict_dir=args.strict_dir,
        tumor_ref_dir=args.tumor_ref_dir,
        primary_artifact_audit_pre=args.primary_artifact_audit_pre,
        primary_artifact_audit_post=args.primary_artifact_audit_post,
        cooccurrence_preflight=args.cooccurrence_preflight,
        tumor_ref_source_identity_receipt=args.tumor_ref_source_identity_receipt,
        independent_m2_audit=args.independent_m2_audit,
        matched_normal_dir=args.matched_normal_dir,
        cn_ccf_annotations=args.cn_ccf_annotations,
    )
    command_args = list(argv) if argv is not None else sys.argv[1:]
    outputs = build_outputs(
        bundle,
        args.output_dir,
        command=[
            *STRICT_PRODUCER.canonical_python_prefix(bundle.strict_dir),
            str(Path(__file__).resolve()),
            *command_args,
        ],
    )
    print(
        json.dumps(
            {
                "input_manifest": str(bundle.manifest.resolve()),
                "screen_dir": str(bundle.screen_dir.resolve()),
                "cooccurrence_dir": str(bundle.cooccurrence_dir.resolve()),
                "strict_dir": str(bundle.strict_dir.resolve()),
                "tumor_ref_dir": str(bundle.tumor_ref_dir.resolve()),
                "primary_artifact_audit_pre": str(
                    bundle.primary_artifact_audit_pre.resolve()
                ),
                "primary_artifact_audit_post": str(
                    bundle.primary_artifact_audit_post.resolve()
                ),
                "cooccurrence_preflight": str(
                    bundle.cooccurrence_preflight.resolve()
                ),
                "tumor_ref_source_identity_receipt": (
                    str(bundle.tumor_ref_source_identity_receipt.resolve())
                    if bundle.tumor_ref_source_identity_receipt is not None
                    else None
                ),
                "independent_m2_audit": (
                    str(bundle.independent_m2_audit.resolve())
                    if bundle.independent_m2_audit is not None
                    else None
                ),
                "matched_normal_dir": (
                    str(bundle.matched_normal_dir.resolve())
                    if bundle.matched_normal_dir is not None
                    else None
                ),
                "cn_ccf_annotations": (
                    str(bundle.cn_ccf_annotations.resolve())
                    if bundle.cn_ccf_annotations is not None
                    else None
                ),
                "output_dir": str(args.output_dir.resolve()),
                "outputs": {name: str(path.resolve()) for name, path in outputs.items()},
                "pass": True,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
