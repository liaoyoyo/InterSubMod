#!/usr/bin/env python3
"""Validate methyl-group associations with nearby frozen LongPhase-S sSNVs.

Selection is deliberately separated from annotation: stable methyl-group sites
and nearby PASS biallelic sSNVs are selected first. Truth, ledger, and topology
are joined only after all read-level statistics have been computed.
"""

from __future__ import annotations

import argparse
import atexit
from bisect import bisect_left, bisect_right
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np
import pysam
from scipy.stats import chi2_contingency, fisher_exact


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import latest_tag_join as TAGS  # noqa: E402
import m2_screen_gate as M2_GATE  # noqa: E402

REPO_ROOT = SCRIPT_DIR.parents[2]
LIB_PATH = REPO_ROOT / "scripts" / "ssnv_cooccurrence_lib.py"
LIB_SPEC = importlib.util.spec_from_file_location("intersubmod_ssnv_cooccurrence_lib", LIB_PATH)
if LIB_SPEC is None or LIB_SPEC.loader is None:
    raise RuntimeError(f"Cannot load shared cooccurrence library: {LIB_PATH}")
LIB = importlib.util.module_from_spec(LIB_SPEC)
sys.modules[LIB_SPEC.name] = LIB
LIB_SPEC.loader.exec_module(LIB)


SCHEMA_NAME = "intersubmod.methyl_ssnv_cooccurrence_validation"
SCHEMA_VERSION = "2.0.0"
ASSIGNMENT_SCHEMA = "intersubmod.all_ssnv_stable_multigroup_read_assignment"
MANIFEST_SCHEMA = "intersubmod.all_ssnv_focal_alt_input_manifest"
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
    {
        "raw_alignment",
        "raw_alignment_index",
        "latest_read_tag_sidecar",
    }
)
MANIFEST_FROZEN_FIELDS = (
    "all_ssnv_vcf",
    "all_ssnv_vcf_index",
    "raw_alignment",
    "raw_alignment_index",
    "latest_read_tag_sidecar",
    "latest_read_tag_sidecar_index",
    "site_ledger",
    "site_ledger_index",
    "layered_region_view",
    "truth_fp_vcf",
    "truth_fp_vcf_index",
    "truth_tp_vcf",
    "truth_tp_vcf_index",
)
SCREEN_CONTRACT = "phylo-v4.1_column_null95_modal_K10_RNULL40_min_group3"
PAIR_WINDOW_BP = 5_000
MIN_BASE_QUALITY = 20
MIN_MAPQ = 20
MIN_QUERY_LENGTH = 1_000
MIN_PERMUTATIONS = 999
RAW_DUPLICATE_EQUIVALENCE_CONTRACT = "sam_core_and_all_aux_tags_except_RG_exact_v1"
RAW_DUPLICATE_ALLOWED_DIFFERING_TAGS = frozenset({"RG"})
RAW_IDENTITY_MISSING_POLICY = "hard_fail_before_site_result"
RAW_IDENTITY_CONFLICT_POLICY = "hard_fail_before_site_result"
ENDPOINT_A_Q_MAX = 0.05
ENDPOINT_A_V_MIN = 0.30
ENDPOINT_A_DELTA_MIN = 0.50
CONDITIONAL_P_MAX = 0.05
EXACT_STATE_SPACE_CEILING = LIB.DEFAULT_EXACT_STATE_SPACE_CEILING
SPACED_MARKER_MIN_SEPARATION_BP = 20
# Deprecated compatibility alias. Genomic spacing does not imply statistical independence.
INDEPENDENT_MARKER_SEPARATION_BP = SPACED_MARKER_MIN_SEPARATION_BP
DEFAULT_TOP_MARKERS = 3
CROSS_PLATFORM_FOUR_STATE_COMPATIBLE_RELATIONS = frozenset(
    {
        "FOCAL_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
        "PARTNER_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
        "BRANCHING_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
        "MULTIPLE_MUTATION_ORDER_MODELS_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
    }
)

EXPECTED_GEOMETRY = {
    "n_all_ssnv_sites": 469_849,
    "n_sites_with_at_least_one_partner": 238_664,
    "n_unordered_pairs": 864_255,
    "unordered_pairs_by_sample": {
        "HCC1395": 287_727,
        "HCC1395_DORADO": 285_217,
        "COLO829": 6_975,
        "H1437": 142_185,
        "H2009": 122_004,
        "HCC1937": 10_572,
        "HCC1954": 9_575,
    },
}

PAIR_OUTPUT_NAME = "methyl_ssnv_pair_results.tsv.gz"
SITE_OUTPUT_NAME = "methyl_ssnv_site_results.tsv.gz"
CASE_OUTPUT_NAME = "oracle_cases.json"
SUMMARY_OUTPUT_NAME = "summary.json"
RECEIPT_OUTPUT_NAME = "run_receipt.json"
RAW_IDENTITY_DUPLICATE_OUTPUT_NAME = "raw_identity_duplicate_audit.tsv.gz"

ORACLE_FOCALS = (
    {
        "case_id": "H2009_chr3_193395128",
        "sample": "H2009",
        "chrom": "chr3",
        "pos": 193_395_128,
        "expected_partner_positions": [],
    },
    {
        "case_id": "HCC1395_DORADO_chr5_750311",
        "sample": "HCC1395_DORADO",
        "chrom": "chr5",
        "pos": 750_311,
        "expected_partner_positions": [748_455, 750_429, 750_926, 751_297],
    },
    {
        "case_id": "COLO829_chr4_92183865",
        "sample": "COLO829",
        "chrom": "chr4",
        "pos": 92_183_865,
        "expected_partner_positions": [92_182_834],
    },
)

SHARED_READSET_ORACLE = {
    "case_id": "COLO829_chr20_52761564_52761565_shared_readset",
    "sample": "COLO829",
    "chrom": "chr20",
    "positions": [52_761_564, 52_761_565],
}

PAIR_COLUMNS = [
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
    "distance_bp",
    "partner_universe_contract",
    "n_all_focal_ref_alt_reads",
    "n_core_focal_alt_reads",
    "core_partner_call_counts",
    "all_partner_call_counts",
    "endpoint_a_testable",
    "endpoint_a_reason",
    "endpoint_a_groups",
    "endpoint_a_table",
    "endpoint_a_n_informative",
    "endpoint_a_analytic_test",
    "endpoint_a_p_analytic",
    "endpoint_a_p_fixed_margins_exact",
    "endpoint_a_exact_method",
    "endpoint_a_exact_state_space_status",
    "endpoint_a_exact_state_space_size",
    "endpoint_a_exact_state_space_lower_bound",
    "endpoint_a_exact_state_space_ceiling",
    "endpoint_a_exact_fallback",
    "endpoint_a_global_fdr_family_status",
    "endpoint_a_q_global_bh",
    "endpoint_a_q_global_by",
    "endpoint_a_cramers_v",
    "endpoint_a_delta_alt_fraction",
    "endpoint_a_dominant_group",
    "endpoint_a_direction_contrast",
    "endpoint_a_effect_direction",
    "endpoint_a_effect_gate_pass",
    "endpoint_a_exact_bh_discovery",
    "endpoint_a_exact_by_discovery",
    "endpoint_a_pre_candidate",
    "endpoint_a_pre_candidate_by_sensitivity",
    "endpoint_a_conditional_strata",
    "endpoint_a_p_conditional_perm",
    "endpoint_a_permutations",
    "endpoint_a_permutable",
    "endpoint_a_conditional_status",
    "endpoint_a_conditional_sensitivity_pass",
    "endpoint_a_formal_pair_by_confirmed",
    "endpoint_a_confirmed_association",
    "endpoint_a_confirmed_by_sensitivity",
    "callability_testable",
    "callability_reason",
    "callability_table",
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
    "endpoint_b_claim_guardrail",
    "top_coverage_rank",
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
    "topology_claim_guardrail",
    "cross_platform_exact_pair_present",
    "cross_platform_both_conditional_by",
    "cross_platform_direction_compatible",
    "cross_platform_effect_compatible",
    "cross_platform_four_state_compatible",
    "cross_platform_cramers_v_abs_difference",
    "cross_platform_delta_alt_fraction_abs_difference",
    "cross_platform_core_read_name_intersection_n",
    "cross_platform_core_read_name_union_n",
    "cross_platform_core_read_name_jaccard",
    "cross_platform_core_read_name_overlap_present",
    "cross_platform_molecule_independence_status",
    "cross_platform_biological_n",
    "cross_platform_conditional_by_effect_direction_replication_pass",
    "cross_platform_replication_status",
]

SITE_COLUMNS = [
    "sample",
    "biological_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "truth_label",
    "ssnv_branch",
    "component_id",
    "selected_group_id",
    "alt_readset_sha256",
    "region_dir",
    "screen_contract",
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
    "n_all_focal_ref_alt_reads",
    "n_focal_ref_reads",
    "n_focal_alt_reads",
    "n_core_focal_alt_reads",
    "raw_identity_expected_projections",
    "raw_identity_matched_records",
    "raw_identity_duplicate_projections_collapsed",
    "raw_identity_duplicate_extra_records_collapsed",
    "raw_identity_exact_duplicate_projections_collapsed",
    "raw_identity_rg_only_duplicate_projections_collapsed",
    "raw_identity_duplicate_projection_sha256",
    "n_partner_markers",
    "n_pair_rows_reconciled",
    "pair_row_count_reconciliation_pass",
    "n_endpoint_a_testable_markers",
    "n_endpoint_a_exact_identifiable_markers",
    "n_endpoint_a_exact_not_identifiable_markers",
    "n_endpoint_a_conditional_permutable_markers",
    "n_pair_bh_discoveries",
    "pair_bh_discovery_positions",
    "n_pair_by_discoveries",
    "pair_by_discovery_positions",
    "n_pair_by_confirmed",
    "pair_by_confirmed_positions",
    "n_spatially_separated_pair_by_20bp",
    "spatially_separated_pair_by_positions_20bp",
    "spaced_marker_20bp_contract",
    "n_endpoint_a_pre_candidates",
    "n_confirmed_markers",
    "confirmed_marker_positions",
    "n_independent_confirmed_markers_20bp",
    "genetically_anchored_multi_marker_candidate",
    "n_confirmed_markers_by_sensitivity",
    "confirmed_marker_positions_by_sensitivity",
    "n_independent_confirmed_markers_20bp_by_sensitivity",
    "genetically_anchored_multi_marker_candidate_by_sensitivity",
    "top_marker_positions",
    "n_top_marker_pair_by_confirmed",
    "top_marker_pair_by_confirmed_positions",
    "joint_signature_complete_marker_support",
    "joint_signature_n_complete_marker_effect_supported",
    "joint_signature_complete_marker_effect_supported_positions",
    "joint_signature_complete_marker_support_contract",
    "n_same_complete_read_effect_supported_top_pair_by",
    "same_complete_read_effect_supported_top_pair_by_positions",
    "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp",
    "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp",
    "top_marker_signatures",
    "joint_signature_n_complete_reads",
    "joint_signature_testable",
    "joint_signature_reason",
    "joint_signature_groups",
    "joint_signature_categories",
    "joint_signature_table",
    "joint_signature_cramers_v",
    "joint_signature_conditional_strata",
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
    "joint_signature_postselection_fdr_calibrated",
    "joint_signature_claim_guardrail",
    "multi_marker_molecular_haplotype_base_candidate",
    "signature_n_complete_reads",
    "signature_testable",
    "signature_reason",
    "signature_table",
    "signature_cramers_v",
    "signature_p_analytic",
    "signature_q_global_bh",
    "signature_q_global_by",
    "signature_by_sensitivity_pass",
    "biological_site_key",
    "biological_site_representative",
    "component_dedup_key",
    "component_representative",
    "alt_readset_dedup_key",
    "alt_readset_representative",
]

RAW_IDENTITY_DUPLICATE_COLUMNS = [
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
]


_BAM_HANDLES: dict[tuple[str, str], pysam.AlignmentFile] = {}
_TABIX_HANDLES: dict[tuple[str, str], pysam.TabixFile] = {}
_VCF_HANDLES: dict[tuple[str, str], pysam.VariantFile] = {}
_LEDGER_HEADERS: dict[str, list[str]] = {}
_LEDGER_ROWS: dict[tuple[str, str, int, str, str], dict[str, str]] = {}
_TRUTH_LABELS: dict[tuple[str, str, str, int, str, str], str] = {}


@dataclass(frozen=True)
class RecoveredRead:
    read_id: str
    read_name: str
    expected_focal_call: str
    focal_call: str
    latest_hp: str
    latest_ps: int | None
    strand: str
    label: str | None
    partner_calls: dict[int, str]
    readset_identity: str = ""


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def json_text(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def open_text(path: Path) -> Any:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path, include_hash: bool = True) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
    }
    if include_hash:
        result["sha256"] = sha256(path)
    return result


def code_source_paths() -> dict[str, Path]:
    return {
        "analyzer": Path(__file__).resolve(),
        "ssnv_cooccurrence_lib": Path(LIB.__file__).resolve(),
        "latest_tag_join": Path(TAGS.__file__).resolve(),
        "m2_screen_gate": Path(M2_GATE.__file__).resolve(),
    }


def capture_code_source_identity() -> dict[str, dict[str, Any]]:
    return {name: artifact(path) for name, path in code_source_paths().items()}


def capture_code_source_modes() -> dict[str, str]:
    return {
        name: oct(path.stat().st_mode & 0o777)
        for name, path in code_source_paths().items()
    }


def require_read_only_source_modes(modes: Mapping[str, str], stage: str) -> None:
    unlocked = {name: mode for name, mode in modes.items() if mode != "0o444"}
    if unlocked:
        raise RuntimeError(f"Code sources are not mode 0444 at {stage}: {unlocked}")


def parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes"}:
        return True
    if normalized in {"0", "false", "no", ""}:
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def is_tumor(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "tumor", "t"}


def hp_family(tag: Any) -> str:
    value = str(tag).strip()
    if value in {"1", "HP1", "1-1", "1-2"}:
        return "HP1-side"
    if value in {"2", "HP2", "2-1", "2-2"}:
        return "HP2-side"
    if value == "3":
        return "HP3-ambiguous"
    if value == "4":
        return "HP4-both"
    return "untagged"


def conditional_stratum(read: RecoveredRead) -> str:
    phase_set = "." if read.latest_ps is None else str(read.latest_ps)
    return f"{hp_family(read.latest_hp)}|PS={phase_set}|strand={read.strand}"


def m2_categorical_level_counts(
    assignment: Mapping[str, Any], screen_row: Mapping[str, Any]
) -> dict[str, int]:
    core_reads = assignment.get("core_reads")
    if not isinstance(core_reads, list) or not core_reads:
        raise RuntimeError("Stable assignment has no core_reads for M2 category proof")
    required = ("label", "latest_hp", "strand")
    for index, read in enumerate(core_reads):
        if not isinstance(read, Mapping):
            raise RuntimeError(f"Stable assignment core read {index} is not an object")
        missing = [field for field in required if field not in read or read[field] is None]
        if missing:
            raise RuntimeError(
                f"Stable assignment core read {index} lacks M2 category fields: {missing}"
            )
    assignment_groups = Counter(str(read["label"]) for read in core_reads)
    screen_groups = M2_GATE.parse_cluster_sizes(screen_row["cluster_sizes"], max_groups=None)
    if dict(assignment_groups) != screen_groups:
        raise RuntimeError(
            "Stable assignment coarse-label counts conflict with screen cluster_sizes: "
            f"assignment={dict(assignment_groups)} screen={screen_groups}"
        )
    return {
        "hp_exact": len({str(read["latest_hp"]) for read in core_reads}),
        "hp_family": len({hp_family(read["latest_hp"]) for read in core_reads}),
        "strand": len({str(read["strand"]) for read in core_reads}),
    }


def m2_screen_eligibility(
    screen_row: dict[str, Any],
    categorical_level_counts: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    return M2_GATE.evaluate_m2_screen(
        screen_row, categorical_level_counts=categorical_level_counts
    )


def deterministic_seed(*parts: Any) -> int:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big") % (2**32)


def focal_alt_readset_digest(reads: Sequence[RecoveredRead]) -> str:
    identities = sorted(read.readset_identity for read in reads if read.focal_call == "A")
    if not identities or any(not identity for identity in identities):
        raise RuntimeError("Cannot compute focal ALT readset digest from missing read identities")
    return hashlib.sha256("\n".join(identities).encode("utf-8")).hexdigest()


def variant_key(chrom: str, pos: int | str, ref: str, alt: str) -> tuple[str, int, str, str]:
    return chrom, int(pos), ref.upper(), alt.upper()


def site_key(sample: str, chrom: str, pos: int | str) -> tuple[str, str, int]:
    return sample, chrom, int(pos)


def verify_frozen_artifact(record: dict[str, Any], require_hash: bool = False) -> Path:
    path = Path(record["path"])
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"Frozen artifact missing or empty: {path}")
    if "size_bytes" in record and path.stat().st_size != int(record["size_bytes"]):
        raise RuntimeError(f"Frozen artifact size drift: {path}")
    if require_hash and "sha256" not in record:
        raise RuntimeError(f"Frozen artifact lacks required SHA-256: {path}")
    if "sha256" not in record and not {"size_bytes", "mtime_ns"}.issubset(record):
        raise RuntimeError(
            f"Frozen unhashed artifact lacks required size+mtime identity: {path}"
        )
    if "sha256" in record and sha256(path) != record["sha256"]:
        raise RuntimeError(f"Frozen artifact hash drift: {path}")
    if "sha256" not in record and "mtime_ns" in record:
        if path.stat().st_mtime_ns != int(record["mtime_ns"]):
            raise RuntimeError(f"Frozen unhashed artifact mtime drift: {path}")
    return path


def verify_manifest_frozen_field(field: str, record: dict[str, Any]) -> dict[str, Any]:
    if field in MANIFEST_HASH_REQUIRED_FIELDS:
        path = verify_frozen_artifact(record, require_hash=True)
        identity_mode = "sha256_size_path_v1"
    elif field in MANIFEST_LARGE_METADATA_IDENTITY_FIELDS:
        path = verify_frozen_artifact(record, require_hash=False)
        identity_mode = (
            "sha256_size_path_v1"
            if "sha256" in record
            else "explicit_large_artifact_size_mtime_path_v1"
        )
    else:
        raise RuntimeError(f"Frozen manifest field lacks an identity policy: {field}")
    return {
        "field": field,
        "path": str(path.resolve()),
        "identity_mode": identity_mode,
        "size_bytes": path.stat().st_size,
        "mtime_ns": path.stat().st_mtime_ns,
        "sha256": record.get("sha256"),
    }


def close_persistent_handles() -> None:
    for collection in (_BAM_HANDLES, _TABIX_HANDLES, _VCF_HANDLES):
        for handle in collection.values():
            try:
                handle.close()
            except Exception:
                pass
        collection.clear()


atexit.register(close_persistent_handles)


def get_bam(path: str, index_path: str) -> pysam.AlignmentFile:
    key = (path, index_path)
    if key not in _BAM_HANDLES:
        _BAM_HANDLES[key] = pysam.AlignmentFile(path, "rb", index_filename=index_path)
    return _BAM_HANDLES[key]


def get_tabix(path: str, index_path: str) -> pysam.TabixFile:
    key = (path, index_path)
    if key not in _TABIX_HANDLES:
        _TABIX_HANDLES[key] = pysam.TabixFile(path, index=index_path)
    return _TABIX_HANDLES[key]


def get_vcf(path: str, index_path: str) -> pysam.VariantFile:
    key = (path, index_path)
    if key not in _VCF_HANDLES:
        _VCF_HANDLES[key] = pysam.VariantFile(path, index_filename=index_path)
    return _VCF_HANDLES[key]


def cigar_query_length(alignment: pysam.AlignedSegment) -> int:
    """Match htslib bam_cigar2qlen: M/I/S/=/X consume query bases."""
    return sum(
        int(length)
        for operation, length in (alignment.cigartuples or [])
        if operation in {0, 1, 4, 7, 8}
    )


def eligible_alignment(alignment: pysam.AlignedSegment) -> bool:
    if alignment.flag & TAGS.EXCLUDED_BY_INTERSUBMOD_FLAGS:
        return False
    if int(alignment.mapping_quality) < MIN_MAPQ:
        return False
    if cigar_query_length(alignment) < MIN_QUERY_LENGTH:
        return False
    if not alignment.has_tag("MM") or not alignment.has_tag("ML"):
        return False
    return True


def alignment_projection(alignment: pysam.AlignedSegment) -> TAGS.ProjectionKey:
    if alignment.reference_name is None or alignment.reference_end is None:
        raise RuntimeError(f"Eligible alignment lacks reference coordinates: {alignment.query_name}")
    return TAGS.projection_key(
        TAGS.cxx_read_name(alignment.query_name, alignment.flag),
        alignment.reference_name,
        int(alignment.reference_start),
        int(alignment.reference_end),
        int(alignment.mapping_quality),
        TAGS.strand_from_flag(alignment.flag),
    )


def _normalized_aux_value(value: Any) -> Any:
    if hasattr(value, "typecode") and hasattr(value, "tobytes"):
        return (
            "typed_array",
            str(value.typecode),
            len(value),
            value.tobytes(),
        )
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        return tuple(_normalized_aux_value(item) for item in value)
    if isinstance(value, bytearray):
        return bytes(value)
    return value


def raw_alignment_core_payload(alignment: pysam.AlignedSegment) -> dict[str, Any]:
    """Return every SAM core field that can distinguish two raw records."""
    qualities = alignment.query_qualities
    return {
        "query_name": alignment.query_name,
        "flag": int(alignment.flag),
        "reference_id": int(alignment.reference_id),
        "reference_name": alignment.reference_name,
        "reference_start": int(alignment.reference_start),
        "reference_end": int(alignment.reference_end or -1),
        "mapping_quality": int(alignment.mapping_quality),
        "cigartuples": tuple(alignment.cigartuples or ()),
        "next_reference_id": int(alignment.next_reference_id),
        "next_reference_start": int(alignment.next_reference_start),
        "template_length": int(alignment.template_length),
        "query_sequence": alignment.query_sequence,
        "query_qualities": None if qualities is None else bytes(qualities),
    }


def raw_alignment_aux_payload(
    alignment: pysam.AlignedSegment,
    *,
    excluded_tags: frozenset[str] = frozenset(),
) -> dict[str, tuple[tuple[str, Any], ...]]:
    """Return typed auxiliary values, preserving duplicate tag occurrences."""
    grouped: dict[str, list[tuple[str, Any]]] = defaultdict(list)
    for tag, value, value_type in alignment.get_tags(with_value_type=True):
        if tag not in excluded_tags:
            grouped[str(tag)].append((str(value_type), _normalized_aux_value(value)))
    return {tag: tuple(values) for tag, values in sorted(grouped.items())}


def raw_alignment_equivalence_differences(
    first: pysam.AlignedSegment,
    other: pysam.AlignedSegment,
) -> dict[str, list[str]]:
    first_core = raw_alignment_core_payload(first)
    other_core = raw_alignment_core_payload(other)
    core_fields = sorted(field for field in first_core if first_core[field] != other_core[field])
    first_tags = raw_alignment_aux_payload(first)
    other_tags = raw_alignment_aux_payload(other)
    auxiliary_tags = sorted(
        tag
        for tag in set(first_tags).union(other_tags)
        if first_tags.get(tag) != other_tags.get(tag)
    )
    return {"core_fields": core_fields, "auxiliary_tags": auxiliary_tags}


def raw_alignments_are_equivalent_duplicates(
    alignments: Sequence[pysam.AlignedSegment],
) -> tuple[bool, dict[str, Any]]:
    if len(alignments) < 2:
        raise ValueError("Duplicate equivalence requires at least two alignments")
    first = alignments[0]
    core_conflicts: set[str] = set()
    tag_conflicts: set[str] = set()
    for other in alignments[1:]:
        differences = raw_alignment_equivalence_differences(first, other)
        core_conflicts.update(differences["core_fields"])
        tag_conflicts.update(differences["auxiliary_tags"])
    disallowed_tags = tag_conflicts.difference(RAW_DUPLICATE_ALLOWED_DIFFERING_TAGS)
    equivalent = not core_conflicts and not disallowed_tags
    return equivalent, {
        "contract": RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
        "core_conflicts": sorted(core_conflicts),
        "differing_auxiliary_tags": sorted(tag_conflicts),
        "disallowed_auxiliary_tags": sorted(disallowed_tags),
        "classification": (
            "RG_ONLY_DUPLICATE"
            if equivalent and tag_conflicts
            else "EXACT_DUPLICATE"
            if equivalent
            else "CONFLICTING_DUPLICATE"
        ),
    }


def projection_digest(projections: Iterable[TAGS.ProjectionKey]) -> str:
    digest = hashlib.sha256()
    for projection in sorted(projections):
        digest.update(json_text(list(projection)).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _nonnegative_int(value: Any, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise RuntimeError(f"Raw identity audit field {field} is not an integer: {value!r}")
    normalized = int(value)
    if normalized < 0:
        raise RuntimeError(f"Raw identity audit field {field} is negative: {normalized}")
    return normalized


def _require_sha256(value: Any, field: str) -> str:
    normalized = str(value)
    if re.fullmatch(r"[0-9a-f]{64}", normalized) is None:
        raise RuntimeError(f"Raw identity audit field {field} is not lowercase SHA-256")
    return normalized


def validate_raw_identity_audit(
    audit: Mapping[str, Any],
    *,
    expected_recovered_count: int | None = None,
) -> None:
    required = {
        "equivalence_contract",
        "allowed_differing_auxiliary_tags",
        "missing_projection_policy",
        "conflicting_analysis_payload_policy",
        "expected_projections",
        "matched_records",
        "duplicate_projections_collapsed",
        "duplicate_extra_records_collapsed",
        "exact_duplicate_projections_collapsed",
        "rg_only_duplicate_projections_collapsed",
        "duplicate_projection_sha256",
        "duplicate_projection_examples",
        "duplicate_projection_records",
    }
    missing = sorted(required.difference(audit))
    if missing:
        raise RuntimeError(f"Raw identity resolver omitted audit fields: {missing}")
    if audit["equivalence_contract"] != RAW_DUPLICATE_EQUIVALENCE_CONTRACT:
        raise RuntimeError("Raw identity audit equivalence contract drifted")
    if list(audit["allowed_differing_auxiliary_tags"]) != sorted(
        RAW_DUPLICATE_ALLOWED_DIFFERING_TAGS
    ):
        raise RuntimeError("Raw identity audit allowed-tag contract drifted")
    if audit["missing_projection_policy"] != RAW_IDENTITY_MISSING_POLICY:
        raise RuntimeError("Raw identity missing-projection policy drifted")
    if audit["conflicting_analysis_payload_policy"] != RAW_IDENTITY_CONFLICT_POLICY:
        raise RuntimeError("Raw identity conflicting-payload policy drifted")

    expected = _nonnegative_int(audit["expected_projections"], "expected_projections")
    matched = _nonnegative_int(audit["matched_records"], "matched_records")
    duplicates = _nonnegative_int(
        audit["duplicate_projections_collapsed"], "duplicate_projections_collapsed"
    )
    extras = _nonnegative_int(
        audit["duplicate_extra_records_collapsed"], "duplicate_extra_records_collapsed"
    )
    exact = _nonnegative_int(
        audit["exact_duplicate_projections_collapsed"],
        "exact_duplicate_projections_collapsed",
    )
    rg_only = _nonnegative_int(
        audit["rg_only_duplicate_projections_collapsed"],
        "rg_only_duplicate_projections_collapsed",
    )
    if expected < 1:
        raise RuntimeError("Raw identity audit has no expected focal REF/ALT projection")
    if expected_recovered_count is not None and expected != int(expected_recovered_count):
        raise RuntimeError(
            "Raw identity expected projection count does not match recovered reads: "
            f"{expected} != {expected_recovered_count}"
        )
    if matched != expected + extras:
        raise RuntimeError(
            f"Raw identity matched-record invariant failed: {matched} != {expected} + {extras}"
        )
    if duplicates > expected:
        raise RuntimeError("Raw identity duplicate projection count exceeds expected projections")
    if extras < duplicates:
        raise RuntimeError("Raw identity duplicate extra-record count is smaller than projections")
    if exact + rg_only != duplicates:
        raise RuntimeError("Raw identity duplicate classification counts do not reconcile")
    if (duplicates == 0) != (extras == 0):
        raise RuntimeError("Raw identity duplicate and extra-record zero states disagree")

    digest = _require_sha256(
        audit["duplicate_projection_sha256"], "duplicate_projection_sha256"
    )
    records = audit["duplicate_projection_records"]
    examples = audit["duplicate_projection_examples"]
    if not isinstance(records, list) or not isinstance(examples, list):
        raise RuntimeError("Raw identity duplicate records/examples must be lists")
    if len(records) != duplicates:
        raise RuntimeError("Raw identity duplicate sparse-record count does not reconcile")
    if len(examples) != min(duplicates, 10):
        raise RuntimeError("Raw identity bounded example count does not reconcile")

    projections: list[Any] = []
    observed_exact = 0
    observed_rg_only = 0
    observed_extras = 0
    for record in records:
        if not isinstance(record, Mapping):
            raise RuntimeError("Raw identity duplicate sparse record is not an object")
        projection = record.get("projection")
        if not isinstance(projection, (list, tuple)) or len(projection) != 6:
            raise RuntimeError("Raw identity duplicate projection is malformed")
        record_count = _nonnegative_int(record.get("record_count"), "record_count")
        if record_count < 2:
            raise RuntimeError("Raw identity duplicate sparse record has fewer than two records")
        classification = record.get("classification")
        differing_tags = record.get("differing_auxiliary_tags")
        if not isinstance(differing_tags, list):
            raise RuntimeError("Raw identity differing auxiliary tags are not a list")
        if classification == "EXACT_DUPLICATE":
            if differing_tags:
                raise RuntimeError("EXACT duplicate unexpectedly reports differing tags")
            observed_exact += 1
        elif classification == "RG_ONLY_DUPLICATE":
            if differing_tags != ["RG"]:
                raise RuntimeError("RG-only duplicate has a non-RG difference set")
            observed_rg_only += 1
        else:
            raise RuntimeError(f"Unexpected raw identity duplicate classification: {classification}")
        observed_extras += record_count - 1
        projections.append(tuple(projection))
    if len(set(projections)) != len(projections):
        raise RuntimeError("Raw identity duplicate sparse records contain repeated projections")
    if observed_exact != exact or observed_rg_only != rg_only or observed_extras != extras:
        raise RuntimeError("Raw identity duplicate sparse records do not reconcile with counts")
    if examples != [list(projection) for projection in projections[:10]]:
        raise RuntimeError("Raw identity duplicate examples do not match sparse records")
    if projection_digest(projections) != digest:
        raise RuntimeError("Raw identity duplicate projection digest does not match sparse records")
    if duplicates == 0 and digest != projection_digest([]):
        raise RuntimeError("Empty raw identity duplicate set has a non-empty-set digest")
    if duplicates > 0 and digest == projection_digest([]):
        raise RuntimeError("Non-empty raw identity duplicate set has the empty-set digest")


def validate_raw_identity_site_fields(
    row: Mapping[str, Any],
    *,
    expected_recovered_count: int | None = None,
) -> None:
    required = {
        "raw_identity_expected_projections",
        "raw_identity_matched_records",
        "raw_identity_duplicate_projections_collapsed",
        "raw_identity_duplicate_extra_records_collapsed",
        "raw_identity_exact_duplicate_projections_collapsed",
        "raw_identity_rg_only_duplicate_projections_collapsed",
        "raw_identity_duplicate_projection_sha256",
    }
    missing = sorted(required.difference(row))
    if missing:
        raise RuntimeError(f"Raw identity site row omitted audit fields: {missing}")
    expected = _nonnegative_int(
        row["raw_identity_expected_projections"], "raw_identity_expected_projections"
    )
    matched = _nonnegative_int(
        row["raw_identity_matched_records"], "raw_identity_matched_records"
    )
    duplicates = _nonnegative_int(
        row["raw_identity_duplicate_projections_collapsed"],
        "raw_identity_duplicate_projections_collapsed",
    )
    extras = _nonnegative_int(
        row["raw_identity_duplicate_extra_records_collapsed"],
        "raw_identity_duplicate_extra_records_collapsed",
    )
    exact = _nonnegative_int(
        row["raw_identity_exact_duplicate_projections_collapsed"],
        "raw_identity_exact_duplicate_projections_collapsed",
    )
    rg_only = _nonnegative_int(
        row["raw_identity_rg_only_duplicate_projections_collapsed"],
        "raw_identity_rg_only_duplicate_projections_collapsed",
    )
    digest = _require_sha256(
        row["raw_identity_duplicate_projection_sha256"],
        "raw_identity_duplicate_projection_sha256",
    )
    if expected < 1:
        raise RuntimeError("Raw identity site row has no expected projection")
    if expected_recovered_count is not None and expected != int(expected_recovered_count):
        raise RuntimeError("Raw identity site expected count differs from recovered-read count")
    if matched != expected + extras:
        raise RuntimeError("Raw identity site matched-record invariant failed")
    if duplicates > expected or extras < duplicates:
        raise RuntimeError("Raw identity site duplicate-count bounds failed")
    if exact + rg_only != duplicates:
        raise RuntimeError("Raw identity site classification invariant failed")
    if (duplicates == 0) != (extras == 0):
        raise RuntimeError("Raw identity site duplicate zero-state invariant failed")
    if duplicates == 0 and digest != projection_digest([]):
        raise RuntimeError("Raw identity site empty duplicate digest invariant failed")
    if duplicates > 0 and digest == projection_digest([]):
        raise RuntimeError("Raw identity site non-empty duplicate digest invariant failed")


def raw_identity_site_digest_payload(row: Mapping[str, Any]) -> list[Any]:
    validate_raw_identity_site_fields(
        row,
        expected_recovered_count=(
            int(row["n_all_focal_ref_alt_reads"])
            if "n_all_focal_ref_alt_reads" in row
            else None
        ),
    )
    return [
        row["sample"],
        row["chrom"],
        int(row["pos"]),
        int(row["raw_identity_expected_projections"]),
        int(row["raw_identity_matched_records"]),
        int(row["raw_identity_duplicate_projections_collapsed"]),
        int(row["raw_identity_duplicate_extra_records_collapsed"]),
        int(row["raw_identity_exact_duplicate_projections_collapsed"]),
        int(row["raw_identity_rg_only_duplicate_projections_collapsed"]),
        row["raw_identity_duplicate_projection_sha256"],
    ]


def resolve_unique_eligible_alignments(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos1: int,
    expected_projections: Iterable[TAGS.ProjectionKey],
    audit: dict[str, Any] | None = None,
) -> dict[TAGS.ProjectionKey, pysam.AlignedSegment]:
    """Resolve reads.tsv projections, collapsing only exact or RG-only raw duplicates."""
    expected = set(expected_projections)
    if not expected:
        raise RuntimeError(f"No focal REF/ALT reads to resolve at {chrom}:{pos1}")
    matches: dict[TAGS.ProjectionKey, list[pysam.AlignedSegment]] = defaultdict(list)
    try:
        fetched = bam.fetch(chrom, pos1 - 1, pos1)
    except (ValueError, OSError) as error:
        raise RuntimeError(f"Cannot fetch raw BAM at {chrom}:{pos1}: {error}") from error
    for alignment in fetched:
        if not eligible_alignment(alignment):
            continue
        projection = alignment_projection(alignment)
        if projection in expected:
            matches[projection].append(alignment)
    missing = sorted(expected.difference(matches))
    conflicting: list[dict[str, Any]] = []
    duplicate_projections: list[TAGS.ProjectionKey] = []
    exact_duplicate_projections: list[TAGS.ProjectionKey] = []
    rg_only_duplicate_projections: list[TAGS.ProjectionKey] = []
    duplicate_projection_records: list[dict[str, Any]] = []
    for key, values in sorted(matches.items()):
        if len(values) == 1:
            continue
        duplicate_projections.append(key)
        equivalent, details = raw_alignments_are_equivalent_duplicates(values)
        if not equivalent:
            conflicting.append({"projection": key, "record_count": len(values), **details})
        elif details["classification"] == "EXACT_DUPLICATE":
            exact_duplicate_projections.append(key)
        else:
            rg_only_duplicate_projections.append(key)
        if equivalent:
            duplicate_projection_records.append(
                {
                    "projection": list(key),
                    "record_count": len(values),
                    "classification": details["classification"],
                    "differing_auxiliary_tags": details[
                        "differing_auxiliary_tags"
                    ],
                }
            )
    if missing or conflicting:
        raise RuntimeError(
            f"Raw BAM identity recovery failed at {chrom}:{pos1}: "
            f"missing={missing[:3]} conflicting_analysis_payload={conflicting[:3]}"
        )
    resolved_audit = {
        "equivalence_contract": RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
        "allowed_differing_auxiliary_tags": sorted(
            RAW_DUPLICATE_ALLOWED_DIFFERING_TAGS
        ),
        "missing_projection_policy": RAW_IDENTITY_MISSING_POLICY,
        "conflicting_analysis_payload_policy": RAW_IDENTITY_CONFLICT_POLICY,
        "expected_projections": len(expected),
        "matched_records": sum(len(values) for values in matches.values()),
        "duplicate_projections_collapsed": len(duplicate_projections),
        "duplicate_extra_records_collapsed": sum(
            len(matches[key]) - 1 for key in duplicate_projections
        ),
        "exact_duplicate_projections_collapsed": len(
            exact_duplicate_projections
        ),
        "rg_only_duplicate_projections_collapsed": len(
            rg_only_duplicate_projections
        ),
        "duplicate_projection_sha256": projection_digest(duplicate_projections),
        "duplicate_projection_examples": [
            list(key) for key in duplicate_projections[:10]
        ],
        "duplicate_projection_records": duplicate_projection_records,
    }
    validate_raw_identity_audit(
        resolved_audit,
        expected_recovered_count=len(expected),
    )
    if audit is not None:
        audit.clear()
        audit.update(resolved_audit)
    return {key: values[0] for key, values in matches.items()}


def _record_is_pass_biallelic_snv(record: pysam.VariantRecord) -> bool:
    filters = set(record.filter.keys())
    pass_filter = not filters or filters == {"PASS"}
    alleles = record.alleles or ()
    return pass_filter and len(alleles) == 2 and all(len(str(allele)) == 1 for allele in alleles)


def discover_focal_and_partners(
    vcf: pysam.VariantFile,
    chrom: str,
    focal_pos: int,
) -> tuple[LIB.Variant, list[LIB.Variant]]:
    low = max(1, focal_pos - PAIR_WINDOW_BP)
    high = focal_pos + PAIR_WINDOW_BP
    records: list[LIB.Variant] = []
    try:
        fetched = vcf.fetch(chrom, low - 1, high)
    except (ValueError, OSError) as error:
        raise RuntimeError(f"Cannot fetch frozen LongPhase-S VCF at {chrom}:{focal_pos}: {error}") from error
    for record in fetched:
        if not _record_is_pass_biallelic_snv(record):
            raise RuntimeError(
                f"Frozen all-sSNV VCF contains non-PASS/non-biallelic-SNV record: "
                f"{record.chrom}:{record.pos}"
            )
        if abs(int(record.pos) - focal_pos) <= PAIR_WINDOW_BP:
            records.append(LIB.Variant(record.chrom, int(record.pos), record.ref, record.alts[0]))
    by_position: dict[int, LIB.Variant] = {}
    for variant in records:
        if variant.pos in by_position:
            raise RuntimeError(f"Multiple frozen sSNVs share {chrom}:{variant.pos}")
        by_position[variant.pos] = variant
    if focal_pos not in by_position:
        raise RuntimeError(f"Stable focal site absent from frozen all-sSNV VCF: {chrom}:{focal_pos}")
    focal = by_position.pop(focal_pos)
    partners = sorted(by_position.values(), key=lambda marker: marker.pos)
    return focal, partners


def count_partner_geometry(
    positions_by_chrom: dict[str, Sequence[int]],
    window_bp: int = PAIR_WINDOW_BP,
) -> dict[str, int]:
    """Count every unordered within-window pair and every site incident to one.

    The right pointer is monotone, but every focal i contributes j-i-1 edges.
    Counting only newly advanced j values would undercount dense windows.
    """
    n_sites = 0
    n_with_partner = 0
    n_pairs = 0
    for chrom in sorted(positions_by_chrom):
        positions = sorted(int(value) for value in positions_by_chrom[chrom])
        if any(right <= left for left, right in zip(positions, positions[1:])):
            raise RuntimeError(f"sSNV positions must be unique and increasing after sort: {chrom}")
        n_sites += len(positions)
        right = 1
        for index, position in enumerate(positions):
            right = max(right, index + 1)
            while right < len(positions) and positions[right] - position <= window_bp:
                right += 1
            n_pairs += right - index - 1
            has_left = index > 0 and position - positions[index - 1] <= window_bp
            has_right = index + 1 < len(positions) and positions[index + 1] - position <= window_bp
            n_with_partner += has_left or has_right
    return {
        "n_sites": n_sites,
        "n_sites_with_at_least_one_partner": n_with_partner,
        "n_unordered_pairs": n_pairs,
    }


def audit_partner_geometry(entries: Sequence[dict[str, Any]], full_scope: bool) -> dict[str, Any]:
    by_sample: dict[str, dict[str, int]] = {}
    for entry in entries:
        positions_by_chrom: dict[str, list[int]] = defaultdict(list)
        with pysam.VariantFile(
            entry["all_ssnv_vcf"]["path"],
            index_filename=entry["all_ssnv_vcf_index"]["path"],
        ) as vcf:
            for record in vcf:
                if not _record_is_pass_biallelic_snv(record):
                    raise RuntimeError(
                        f"Geometry audit encountered non-PASS/non-biallelic-SNV: "
                        f"{entry['sample']}:{record.chrom}:{record.pos}"
                    )
                positions_by_chrom[str(record.chrom)].append(int(record.pos))
        observed = count_partner_geometry(positions_by_chrom)
        expected_sites = int(entry["counts"]["all_ssnv"])
        expected_pairs = EXPECTED_GEOMETRY["unordered_pairs_by_sample"].get(entry["sample"])
        if observed["n_sites"] != expected_sites:
            raise RuntimeError(
                f"{entry['sample']} geometry site count mismatch: "
                f"{observed['n_sites']} != {expected_sites}"
            )
        if expected_pairs is None or observed["n_unordered_pairs"] != expected_pairs:
            raise RuntimeError(
                f"{entry['sample']} geometry unordered-pair mismatch: "
                f"{observed['n_unordered_pairs']} != {expected_pairs}"
            )
        by_sample[entry["sample"]] = observed
    aggregate = {
        "n_sites": sum(row["n_sites"] for row in by_sample.values()),
        "n_sites_with_at_least_one_partner": sum(
            row["n_sites_with_at_least_one_partner"] for row in by_sample.values()
        ),
        "n_unordered_pairs": sum(row["n_unordered_pairs"] for row in by_sample.values()),
    }
    if full_scope:
        expected = {
            "n_sites": EXPECTED_GEOMETRY["n_all_ssnv_sites"],
            "n_sites_with_at_least_one_partner": EXPECTED_GEOMETRY[
                "n_sites_with_at_least_one_partner"
            ],
            "n_unordered_pairs": EXPECTED_GEOMETRY["n_unordered_pairs"],
        }
        if aggregate != expected:
            raise RuntimeError(f"Full-scope geometry oracle mismatch: observed={aggregate} expected={expected}")
    return {
        "contract": "inclusive_abs_distance_le_5000; each_i_adds_j_minus_i_minus_1",
        "full_scope": full_scope,
        "by_sample": by_sample,
        "aggregate": aggregate,
        "expected_full_scope": EXPECTED_GEOMETRY,
        "pass": True,
    }


def build_stable_partner_oracle(
    entries: Sequence[dict[str, Any]],
    assignments: dict[tuple[str, str, int], dict[str, Any]],
) -> dict[tuple[str, str, int], set[tuple[str, int, str, str]]]:
    """Independently enumerate each stable focal's partners by full-VCF scan and bisect."""
    focals_by_sample: dict[str, list[tuple[str, str, int]]] = defaultdict(list)
    for key in assignments:
        focals_by_sample[key[0]].append(key)
    oracle: dict[tuple[str, str, int], set[tuple[str, int, str, str]]] = {}
    for entry in entries:
        sample = entry["sample"]
        variants_by_chrom: dict[str, list[LIB.Variant]] = defaultdict(list)
        with pysam.VariantFile(
            entry["all_ssnv_vcf"]["path"],
            index_filename=entry["all_ssnv_vcf_index"]["path"],
        ) as vcf:
            for record in vcf:
                if not _record_is_pass_biallelic_snv(record):
                    raise RuntimeError(
                        f"Stable-partner oracle encountered invalid record: "
                        f"{sample}:{record.chrom}:{record.pos}"
                    )
                variants_by_chrom[str(record.chrom)].append(
                    LIB.Variant(str(record.chrom), int(record.pos), record.ref, record.alts[0])
                )
        positions_by_chrom: dict[str, list[int]] = {}
        for chrom, variants in variants_by_chrom.items():
            variants.sort(key=lambda variant: variant.pos)
            positions = [variant.pos for variant in variants]
            if any(right <= left for left, right in zip(positions, positions[1:])):
                raise RuntimeError(f"Stable-partner oracle found duplicate/non-increasing positions: {sample}:{chrom}")
            positions_by_chrom[chrom] = positions
        for key in focals_by_sample.get(sample, []):
            _, chrom, focal_pos = key
            variants = variants_by_chrom.get(chrom, [])
            positions = positions_by_chrom.get(chrom, [])
            focal_index = bisect_left(positions, focal_pos)
            if focal_index >= len(positions) or positions[focal_index] != focal_pos:
                raise RuntimeError(f"Stable focal absent from sequential VCF oracle: {key}")
            left = bisect_left(positions, focal_pos - PAIR_WINDOW_BP)
            right = bisect_right(positions, focal_pos + PAIR_WINDOW_BP)
            oracle[key] = {
                variant_key(variant.chrom, variant.pos, variant.ref, variant.alt)
                for variant in variants[left:right]
                if variant.pos != focal_pos
            }
    if set(oracle) != set(assignments):
        raise RuntimeError(
            f"Stable-partner oracle focal denominator mismatch: "
            f"oracle={len(oracle)} assignments={len(assignments)}"
        )
    return oracle


def validate_stable_partner_rows(
    oracle: dict[tuple[str, str, int], set[tuple[str, int, str, str]]],
    pair_rows: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    observed: dict[tuple[str, str, int], set[tuple[str, int, str, str]]] = {
        key: set() for key in oracle
    }
    for row in pair_rows:
        key = site_key(row["sample"], row["focal_chrom"], row["focal_pos"])
        if key not in observed:
            raise RuntimeError(f"Analyzed pair row has non-stable focal: {key}")
        partner = variant_key(
            row["partner_chrom"], row["partner_pos"], row["partner_ref"], row["partner_alt"]
        )
        if partner in observed[key]:
            raise RuntimeError(f"Duplicate analyzed focal-partner row: {key} -> {partner}")
        observed[key].add(partner)
    mismatches = []
    for key in sorted(oracle):
        if observed[key] != oracle[key]:
            mismatches.append(
                {
                    "focal": key,
                    "missing": sorted(oracle[key].difference(observed[key]))[:5],
                    "extra": sorted(observed[key].difference(oracle[key]))[:5],
                }
            )
    if mismatches:
        raise RuntimeError(f"Stable focal partner enumeration mismatch: {mismatches[:3]}")
    expected_pairs = sum(len(values) for values in oracle.values())
    return {
        "contract": "sequential_full_VCF_plus_bisect_oracle_vs_indexed_worker_fetch",
        "n_stable_focals": len(oracle),
        "n_stable_focals_with_partner": sum(bool(values) for values in oracle.values()),
        "n_expected_directed_pair_rows": expected_pairs,
        "n_observed_directed_pair_rows": len(pair_rows),
        "pass": True,
    }


def categorical_association(
    labels: Sequence[str],
    categories: Sequence[str],
    min_total: int = 10,
    min_group: int = 3,
    min_category: int = 3,
) -> dict[str, Any]:
    if len(labels) != len(categories):
        raise ValueError("labels/categories length mismatch")
    groups = sorted(set(labels))
    states = sorted(set(categories))
    table = np.zeros((len(groups), len(states)), dtype=int)
    group_index = {value: index for index, value in enumerate(groups)}
    state_index = {value: index for index, value in enumerate(states)}
    for label, category in zip(labels, categories):
        table[group_index[label], state_index[category]] += 1
    result: dict[str, Any] = {
        "testable": False,
        "reason": "insufficient_joint_information",
        "groups": groups,
        "categories": states,
        "table": table.tolist(),
        "n": int(table.sum()),
        "cramers_v": None,
        "p_analytic": None,
        "analytic_test": None,
    }
    if len(groups) < 2 or len(states) < 2 or int(table.sum()) < min_total:
        return result
    if int(table.sum(axis=1).min()) < min_group:
        result["reason"] = "group_below_minimum"
        return result
    if int(table.sum(axis=0).min()) < min_category:
        result["reason"] = "category_below_minimum"
        return result
    if table.shape == (2, 2):
        p_value = float(fisher_exact(table)[1])
        test = "fisher_exact_2x2"
    else:
        p_value = float(chi2_contingency(table, correction=False)[1])
        test = "pearson_chi_square"
    result.update(
        {
            "testable": True,
            "reason": "ok",
            "cramers_v": LIB.cramer_v(table),
            "p_analytic": p_value,
            "analytic_test": test,
        }
    )
    return result


def benjamini_yekutieli(p_values: Iterable[float | None]) -> list[float | None]:
    return LIB.benjamini_yekutieli(p_values)


def endpoint_a_gate(row: dict[str, Any], q_field: str = "endpoint_a_q_global_bh") -> bool:
    q_value = row.get(q_field)
    return bool(
        row.get("endpoint_a_testable")
        and row.get("endpoint_a_global_fdr_family_status") == "ELIGIBLE_M2_EXACT_FAMILY"
        and q_value is not None
        and float(q_value) <= ENDPOINT_A_Q_MAX
    )


def endpoint_a_effect_gate(row: dict[str, Any]) -> bool:
    effect = row.get("endpoint_a_cramers_v")
    delta = row.get("endpoint_a_delta_alt_fraction")
    return bool(
        row.get("endpoint_a_testable")
        and effect is not None
        and delta is not None
        and float(effect) >= ENDPOINT_A_V_MIN
        and float(delta) >= ENDPOINT_A_DELTA_MIN
    )


def apply_endpoint_a_bh_and_gate(rows: list[dict[str, Any]]) -> None:
    endpoint_p = [
        row.get("endpoint_a_p_fixed_margins_exact")
        if row.get("endpoint_a_global_fdr_family_status") == "ELIGIBLE_M2_EXACT_FAMILY"
        else None
        for row in rows
    ]
    callability_p = [row.get("callability_p_analytic") for row in rows]
    q_values = LIB.benjamini_hochberg(endpoint_p)
    by_values = benjamini_yekutieli(endpoint_p)
    callability_q = LIB.benjamini_hochberg(callability_p)
    callability_by = benjamini_yekutieli(callability_p)
    for row, q_value, by_value, call_q, call_by in zip(
        rows, q_values, by_values, callability_q, callability_by
    ):
        row["endpoint_a_q_global_bh"] = q_value
        row["endpoint_a_q_global_by"] = by_value
        row["callability_q_global_bh"] = call_q
        row["callability_q_global_by"] = call_by
        core_counts = row.get("core_partner_call_counts") or {}
        noncallable = int(core_counts.get("O", 0)) + int(core_counts.get("X", 0))
        row["callability_noncallable_core_reads"] = noncallable
        callability_v = row.get("callability_cramers_v")
        if noncallable == 0:
            callability_status = "PASS_ALL_CORE_READS_CALLABLE"
            callability_pass = True
        elif not row.get("callability_testable"):
            callability_status = "NOT_IDENTIFIABLE_DIFFERENTIAL_CALLABILITY"
            callability_pass = False
        elif call_by is None or callability_v is None:
            callability_status = "NOT_IDENTIFIABLE_MISSING_CALLABILITY_STATISTIC"
            callability_pass = False
        elif float(call_by) <= ENDPOINT_A_Q_MAX or float(callability_v) >= ENDPOINT_A_V_MIN:
            callability_status = "FAIL_DIFFERENTIAL_CALLABILITY_SIGNAL"
            callability_pass = False
        else:
            callability_status = "PASS_NO_STRONG_DIFFERENTIAL_CALLABILITY_DETECTED"
            callability_pass = True
        row["callability_gate_status"] = callability_status
        row["callability_gate_pass"] = callability_pass
        row["endpoint_a_effect_gate_pass"] = endpoint_a_effect_gate(row)
        row["endpoint_a_exact_bh_discovery"] = endpoint_a_gate(row)
        row["endpoint_a_exact_by_discovery"] = endpoint_a_gate(
            row, q_field="endpoint_a_q_global_by"
        )
        # Deprecated aliases retained for old consumers; formal claims use exact_* fields.
        row["endpoint_a_pre_candidate"] = bool(
            row["endpoint_a_exact_bh_discovery"] and row["endpoint_a_effect_gate_pass"]
        )
        row["endpoint_a_pre_candidate_by_sensitivity"] = bool(
            row["endpoint_a_exact_by_discovery"] and row["endpoint_a_effect_gate_pass"]
        )


def apply_joint_signature_global_fdr(rows: list[dict[str, Any]]) -> dict[str, Any]:
    p_values: list[float | None] = []
    statuses: list[str] = []
    for row in rows:
        p_value = row.get("joint_signature_p_conditional_perm")
        if not bool(row.get("m2_screen_eligible", False)):
            status = "INELIGIBLE_M2_SCREEN"
        elif not bool(row.get("joint_signature_testable", False)):
            status = "ELIGIBLE_M2_JOINT_SIGNATURE_NOT_TESTABLE"
        elif (
            not bool(row.get("joint_signature_permutable", False))
            or int(row.get("joint_signature_permutations", 0)) != MIN_PERMUTATIONS
            or p_value is None
        ):
            status = "ELIGIBLE_M2_JOINT_SIGNATURE_NOT_PERMUTABLE"
        else:
            status = "ELIGIBLE_M2_JOINT_SIGNATURE_GLOBAL_FDR_FAMILY"
        statuses.append(status)
        p_values.append(
            float(p_value)
            if status == "ELIGIBLE_M2_JOINT_SIGNATURE_GLOBAL_FDR_FAMILY"
            else None
        )
    bh_values = LIB.benjamini_hochberg(p_values)
    by_values = benjamini_yekutieli(p_values)
    for row, status, bh_value, by_value in zip(rows, statuses, bh_values, by_values):
        row["joint_signature_global_fdr_family_status"] = status
        row["joint_signature_q_global_bh"] = bh_value
        row["joint_signature_q_global_by"] = by_value
        row["joint_signature_global_bh_discovery"] = bool(
            bh_value is not None and float(bh_value) <= ENDPOINT_A_Q_MAX
        )
        row["joint_signature_global_by_discovery"] = bool(
            by_value is not None and float(by_value) <= ENDPOINT_A_Q_MAX
        )
    return {
        "family_contract": (
            "all_M2_eligible_joint_signature_testable_permutable_sites_global_adjustment"
        ),
        "n_family_sites": sum(value is not None for value in p_values),
        "n_global_bh_discoveries": sum(
            bool(row["joint_signature_global_bh_discovery"]) for row in rows
        ),
        "n_global_by_discoveries": sum(
            bool(row["joint_signature_global_by_discovery"]) for row in rows
        ),
        "formal_gate": "global_BY_q_le_0.05",
        "pass": True,
    }


def endpoint_b_status(summary: dict[str, Any]) -> str:
    if int(summary["n_focal_ref"]) == 0:
        return "NOT_IDENTIFIABLE_NO_FOCAL_REF"
    if not summary["complete_four_state_testable"]:
        return "NOT_IDENTIFIABLE_INSUFFICIENT_FOUR_STATE"
    return str(summary["relation_compatibility"])


def _validate_assignment_against_rows(
    assignment: dict[str, Any],
    joined_by_read_id: dict[str, dict[str, Any]],
) -> dict[str, str]:
    core_reads = assignment.get("core_reads")
    if not isinstance(core_reads, list) or not core_reads:
        raise RuntimeError("Stable assignment has no core_reads")
    labels_by_read_id: dict[str, str] = {}
    for core in core_reads:
        read_id = str(core["read_id"])
        if read_id in labels_by_read_id:
            raise RuntimeError(f"Duplicate core read_id in stable assignment: {read_id}")
        if read_id not in joined_by_read_id:
            raise RuntimeError(f"Stable core read absent from focal REF/ALT reads.tsv rows: {read_id}")
        row = joined_by_read_id[read_id]
        if row["expected_focal_call"] != "A":
            raise RuntimeError(f"Stable core read is not focal ALT in reads.tsv: {read_id}")
        expected_projection = TAGS.projection_key(
            str(core["read_name"]),
            str(core["chrom"]),
            int(core["start"]),
            int(core["end"]),
            int(core["mapq"]),
            str(core["strand"]),
        )
        if expected_projection != row["projection"]:
            raise RuntimeError(f"Stable assignment/read projection conflict: {read_id}")
        if str(core.get("latest_hp")) != row["latest_hp"]:
            raise RuntimeError(f"Stable assignment/latest HP conflict: {read_id}")
        expected_ps = core.get("latest_ps")
        if expected_ps in {"", "."}:
            expected_ps = None
        if expected_ps is not None:
            expected_ps = int(expected_ps)
        if expected_ps != row["latest_ps"]:
            raise RuntimeError(f"Stable assignment/latest PS conflict: {read_id}")
        if int(core.get("latest_tag_full_identity_count", -1)) != int(
            row["latest_tag_full_identity_count"]
        ):
            raise RuntimeError(f"Stable assignment/latest sidecar identity-count conflict: {read_id}")
        labels_by_read_id[read_id] = str(core["label"])
    parallel_fields = {
        "read_ids": [str(value) for value in assignment.get("read_ids", [])],
        "read_names": [str(value) for value in assignment.get("read_names", [])],
        "labels": [str(value) for value in assignment.get("labels", [])],
        "latest_hp": [str(value) for value in assignment.get("latest_hp", [])],
        "latest_ps": [
            None if value in {None, "", "."} else int(value)
            for value in assignment.get("latest_ps", [])
        ],
    }
    expected_parallel = {
        "read_ids": [str(row["read_id"]) for row in core_reads],
        "read_names": [str(row["read_name"]) for row in core_reads],
        "labels": [str(row["label"]) for row in core_reads],
        "latest_hp": [str(row["latest_hp"]) for row in core_reads],
        "latest_ps": [
            None if row.get("latest_ps") in {None, "", "."} else int(row["latest_ps"])
            for row in core_reads
        ],
    }
    if parallel_fields != expected_parallel:
        raise RuntimeError("Stable assignment parallel read/label arrays conflict with core_reads")
    return labels_by_read_id


def recover_site_reads(
    assignment: dict[str, Any],
    focal: LIB.Variant,
    partners: Sequence[LIB.Variant],
    entry: dict[str, Any],
    raw_identity_audit: dict[str, Any] | None = None,
) -> list[RecoveredRead]:
    reads_path = Path(assignment["region_dir"]) / "reads" / "reads.tsv"
    if not reads_path.is_file() or reads_path.stat().st_size <= 0:
        raise FileNotFoundError(f"InterSubMod reads.tsv missing: {reads_path}")
    sidecar = get_tabix(
        entry["latest_read_tag_sidecar"]["path"],
        entry["latest_read_tag_sidecar_index"]["path"],
    )
    lookup = TAGS.fetch_site_lookup(sidecar, focal.chrom, focal.pos)
    joined = TAGS.join_read_rows(reads_path, lookup)
    joined_by_read_id: dict[str, dict[str, Any]] = {}
    by_projection: dict[TAGS.ProjectionKey, dict[str, Any]] = {}
    for row, tags, full_identity_count in joined:
        if int(full_identity_count) != 1:
            raise RuntimeError(
                "Latest LongPhase-S sidecar projection is not a unique full alignment identity: "
                f"{row.get('read_name')} {row.get('chr')}:{row.get('start')}-{row.get('end')} "
                f"count={full_identity_count}"
            )
        if not is_tumor(row.get("is_tumor")):
            continue
        support = str(row.get("alt_support", "")).strip().upper()
        if support not in {"REF", "ALT"}:
            continue
        if "read_id" not in row or not str(row["read_id"]).strip():
            raise RuntimeError(f"reads.tsv has focal REF/ALT row without read_id: {reads_path}")
        read_id = str(row["read_id"])
        projection = TAGS.read_projection(row)
        if read_id in joined_by_read_id:
            raise RuntimeError(f"Duplicate focal REF/ALT read_id in reads.tsv: {read_id}")
        if projection in by_projection:
            raise RuntimeError(f"Duplicate focal REF/ALT reads.tsv identity projection: {projection}")
        normalized = {
            "read_id": read_id,
            "read_name": str(row["read_name"]),
            "expected_focal_call": "A" if support == "ALT" else "R",
            "latest_hp": str(tags.hp),
            "latest_ps": tags.ps,
            "strand": str(row["strand"]),
            "projection": projection,
            "latest_tag_full_identity_count": int(full_identity_count),
        }
        joined_by_read_id[read_id] = normalized
        by_projection[projection] = normalized
    labels_by_read_id = _validate_assignment_against_rows(assignment, joined_by_read_id)
    bam = get_bam(entry["raw_alignment"]["path"], entry["raw_alignment_index"]["path"])
    alignments = resolve_unique_eligible_alignments(
        bam,
        focal.chrom,
        focal.pos,
        by_projection,
        raw_identity_audit,
    )
    recovered: list[RecoveredRead] = []
    for projection, row in sorted(by_projection.items(), key=lambda item: item[1]["read_id"]):
        alignment = alignments[projection]
        focal_call = LIB.call_snv(alignment, focal, min_base_quality=MIN_BASE_QUALITY)
        if focal_call != row["expected_focal_call"]:
            raise RuntimeError(
                f"Focal call/read.tsv conflict for {row['read_id']} at {focal.chrom}:{focal.pos}: "
                f"raw={focal_call} reads.tsv={row['expected_focal_call']}"
            )
        recovered.append(
            RecoveredRead(
                read_id=row["read_id"],
                read_name=row["read_name"],
                expected_focal_call=row["expected_focal_call"],
                focal_call=focal_call,
                latest_hp=row["latest_hp"],
                latest_ps=row["latest_ps"],
                strand=row["strand"],
                label=labels_by_read_id.get(row["read_id"]),
                partner_calls=LIB.call_snvs(alignment, partners, min_base_quality=MIN_BASE_QUALITY),
                readset_identity=(
                    f"{row['read_name']}|{projection[1]}|{projection[2]}|{projection[3]}|{row['strand']}"
                ),
            )
        )
    if set(labels_by_read_id) != {read.read_id for read in recovered if read.label is not None}:
        raise RuntimeError(f"Not all stable core reads were recovered at {focal.chrom}:{focal.pos}")
    return recovered


def _pair_result(
    entry: dict[str, Any],
    focal: LIB.Variant,
    partner: LIB.Variant,
    recovered: Sequence[RecoveredRead],
    *,
    m2_eligible: bool = True,
    exact_state_space_ceiling: int = EXACT_STATE_SPACE_CEILING,
) -> dict[str, Any]:
    core = [read for read in recovered if read.label is not None]
    labels = [str(read.label) for read in core]
    core_calls = [read.partner_calls[partner.pos] for read in core]
    all_calls = [read.partner_calls[partner.pos] for read in recovered]
    endpoint_a = LIB.methyl_group_marker_association(labels, core_calls)
    exact = (
        LIB.fisher_freeman_halton_kx2(
            endpoint_a["table"], state_space_ceiling=exact_state_space_ceiling
        )
        if endpoint_a["testable"]
        else {
            "method": "fisher_freeman_halton_fixed_margins_probability_ordering",
            "p_value": None,
            "identifiable": False,
            "state_space_status": "NOT_IDENTIFIABLE_ENDPOINT_A_NOT_TESTABLE",
            "state_space_size": None,
            "state_space_lower_bound": None,
            "state_space_ceiling": exact_state_space_ceiling,
            "fallback": "fail_closed_not_identifiable",
        }
    )
    if not m2_eligible:
        family_status = "INELIGIBLE_M2_SCREEN"
    elif not endpoint_a["testable"]:
        family_status = "ELIGIBLE_M2_ENDPOINT_A_NOT_TESTABLE"
    elif not exact["identifiable"]:
        family_status = "ELIGIBLE_M2_EXACT_NOT_IDENTIFIABLE"
    else:
        family_status = "ELIGIBLE_M2_EXACT_FAMILY"
    callability = categorical_association(
        labels,
        ["CALLABLE_RA" if call in {"R", "A"} else "NONCALLABLE_OX" for call in core_calls],
    )
    endpoint_b = LIB.four_state_summary([read.focal_call for read in recovered], all_calls)
    conditional_core = [
        (read, call) for read, call in zip(core, core_calls) if call in {"R", "A"}
    ]
    row: dict[str, Any] = {
        "sample": entry["sample"],
        "biological_id": entry["biological_id"],
        "focal_chrom": focal.chrom,
        "focal_pos": focal.pos,
        "focal_ref": focal.ref,
        "focal_alt": focal.alt,
        "partner_chrom": partner.chrom,
        "partner_pos": partner.pos,
        "partner_ref": partner.ref,
        "partner_alt": partner.alt,
        "distance_bp": partner.pos - focal.pos,
        "partner_universe_contract": (
            "all_frozen_latest_LongPhase-S_PASS_biallelic_sSNVs_within_focal_plus_or_minus_5000bp"
        ),
        "n_all_focal_ref_alt_reads": len(recovered),
        "n_core_focal_alt_reads": len(core),
        "core_partner_call_counts": dict(sorted(Counter(core_calls).items())),
        "all_partner_call_counts": dict(sorted(Counter(all_calls).items())),
        "endpoint_a_testable": endpoint_a["testable"],
        "endpoint_a_reason": endpoint_a["reason"],
        "endpoint_a_groups": endpoint_a["groups"],
        "endpoint_a_table": endpoint_a["table"],
        "endpoint_a_n_informative": endpoint_a["n_informative"],
        "endpoint_a_analytic_test": endpoint_a["analytic_test"],
        "endpoint_a_p_analytic": endpoint_a["p_analytic"],
        "endpoint_a_p_fixed_margins_exact": exact["p_value"],
        "endpoint_a_exact_method": exact["method"],
        "endpoint_a_exact_state_space_status": exact["state_space_status"],
        "endpoint_a_exact_state_space_size": exact["state_space_size"],
        "endpoint_a_exact_state_space_lower_bound": exact["state_space_lower_bound"],
        "endpoint_a_exact_state_space_ceiling": exact["state_space_ceiling"],
        "endpoint_a_exact_fallback": exact["fallback"],
        "endpoint_a_global_fdr_family_status": family_status,
        "endpoint_a_q_global_bh": None,
        "endpoint_a_q_global_by": None,
        "endpoint_a_cramers_v": endpoint_a["cramers_v"],
        "endpoint_a_delta_alt_fraction": endpoint_a["delta_alt_fraction"],
        "endpoint_a_dominant_group": endpoint_a["dominant_group"],
        "endpoint_a_direction_contrast": endpoint_a["direction_contrast"],
        "endpoint_a_effect_direction": endpoint_a["direction"],
        "endpoint_a_effect_gate_pass": False,
        "endpoint_a_exact_bh_discovery": False,
        "endpoint_a_exact_by_discovery": False,
        "endpoint_a_pre_candidate": False,
        "endpoint_a_pre_candidate_by_sensitivity": False,
        "endpoint_a_conditional_strata": "latest_HP_family_x_PS_x_strand",
        "endpoint_a_p_conditional_perm": None,
        "endpoint_a_permutations": 0,
        "endpoint_a_permutable": False,
        "endpoint_a_conditional_status": "NOT_RUN_NOT_EXACT_BY_DISCOVERY",
        "endpoint_a_conditional_sensitivity_pass": False,
        "endpoint_a_formal_pair_by_confirmed": False,
        "endpoint_a_confirmed_association": False,
        "endpoint_a_confirmed_by_sensitivity": False,
        "callability_testable": callability["testable"],
        "callability_reason": callability["reason"],
        "callability_table": callability["table"],
        "callability_p_analytic": callability["p_analytic"],
        "callability_q_global_bh": None,
        "callability_q_global_by": None,
        "callability_cramers_v": callability["cramers_v"],
        "callability_noncallable_core_reads": sum(
            call not in {"R", "A"} for call in core_calls
        ),
        "callability_gate_status": "NOT_RUN_GLOBAL_ADJUSTMENT_PENDING",
        "callability_gate_pass": False,
        "endpoint_b_status": endpoint_b_status(endpoint_b),
        "endpoint_b_n_joint": endpoint_b["n_joint"],
        "endpoint_b_n_called_depth": endpoint_b["n_called_depth"],
        "endpoint_b_state_counts": endpoint_b["state_counts_with_noncallable"],
        "endpoint_b_n_focal_ref": endpoint_b["n_focal_ref"],
        "endpoint_b_n_focal_alt": endpoint_b["n_focal_alt"],
        "endpoint_b_n_partner_ref": endpoint_b["n_partner_ref"],
        "endpoint_b_n_partner_alt": endpoint_b["n_partner_alt"],
        "endpoint_b_focal_ancestor_violation_rate": endpoint_b["focal_ancestor_violation_rate"],
        "endpoint_b_partner_ancestor_violation_rate": endpoint_b["partner_ancestor_violation_rate"],
        "endpoint_b_error_ceiling": endpoint_b["error_ceiling"],
        "endpoint_b_error_model_confidence": endpoint_b["confidence"],
        "endpoint_b_familywise_confidence": endpoint_b["familywise_confidence"],
        "endpoint_b_relation_family_size": endpoint_b["relation_family_size"],
        "endpoint_b_multiplicity_method": endpoint_b["multiplicity_method"],
        "endpoint_b_minimum_zero_violation_depth": endpoint_b["minimum_zero_violation_depth"],
        "endpoint_b_focal_ancestor_violation_p_exact": endpoint_b["focal_ancestor_violation"][
            "p_exact_greater"
        ],
        "endpoint_b_focal_ancestor_violation_upper_bound": endpoint_b[
            "focal_ancestor_violation"
        ]["upper_bound"],
        "endpoint_b_focal_ancestor_violation_threshold": endpoint_b["focal_ancestor_violation"][
            "threshold"
        ],
        "endpoint_b_focal_ancestor_violation_status": endpoint_b["focal_ancestor_violation"][
            "status"
        ],
        "endpoint_b_partner_ancestor_violation_p_exact": endpoint_b["partner_ancestor_violation"][
            "p_exact_greater"
        ],
        "endpoint_b_partner_ancestor_violation_upper_bound": endpoint_b[
            "partner_ancestor_violation"
        ]["upper_bound"],
        "endpoint_b_partner_ancestor_violation_threshold": endpoint_b[
            "partner_ancestor_violation"
        ]["threshold"],
        "endpoint_b_partner_ancestor_violation_status": endpoint_b["partner_ancestor_violation"][
            "status"
        ],
        "endpoint_b_branching_violation_p_exact": endpoint_b["branching_violation"][
            "p_exact_greater"
        ],
        "endpoint_b_branching_violation_upper_bound": endpoint_b["branching_violation"][
            "upper_bound"
        ],
        "endpoint_b_branching_violation_threshold": endpoint_b["branching_violation"][
            "threshold"
        ],
        "endpoint_b_branching_violation_status": endpoint_b["branching_violation"]["status"],
        "endpoint_b_complete_four_state_testable": endpoint_b["complete_four_state_testable"],
        "endpoint_b_relation_compatibility": endpoint_b["relation_compatibility"],
        "endpoint_b_compatible_relation_models": endpoint_b["compatible_relation_models"],
        "endpoint_b_n_compatible_relation_models": endpoint_b["n_compatible_relation_models"],
        "endpoint_b_claim_guardrail": endpoint_b["claim_guardrail"],
        "top_coverage_rank": None,
        "top_coverage_marker": False,
        "focal_truth_label": None,
        "partner_truth_label": None,
        "focal_ssnv_branch": None,
        "partner_ssnv_branch": None,
        "focal_component_id": None,
        "partner_component_id": None,
        "topology_scope": None,
        "topology_region": None,
        "topology_order_status": None,
        "topology_claim_guardrail": (
            "Topology is posthoc context; compatibility or tree context is not proof of cellular ancestry."
        ),
        "cross_platform_exact_pair_present": False,
        "cross_platform_both_conditional_by": False,
        "cross_platform_direction_compatible": None,
        "cross_platform_effect_compatible": None,
        "cross_platform_four_state_compatible": None,
        "cross_platform_cramers_v_abs_difference": None,
        "cross_platform_delta_alt_fraction_abs_difference": None,
        "cross_platform_core_read_name_intersection_n": None,
        "cross_platform_core_read_name_union_n": None,
        "cross_platform_core_read_name_jaccard": None,
        "cross_platform_core_read_name_overlap_present": None,
        "cross_platform_molecule_independence_status": "NOT_APPLICABLE",
        "cross_platform_biological_n": None,
        "cross_platform_conditional_by_effect_direction_replication_pass": False,
        "cross_platform_replication_status": "NOT_APPLICABLE",
        "_permutation_payload": {
            "labels": [str(read.label) for read, _ in conditional_core],
            "categories": [call for _, call in conditional_core],
            "strata": [conditional_stratum(read) for read, _ in conditional_core],
            "core_raw_read_names": sorted({read.read_name for read in core}),
        },
    }
    return row


def _signature_result(
    recovered: Sequence[RecoveredRead],
    ranked_partners: Sequence[LIB.Variant],
    *,
    seed: int,
    permutations: int = MIN_PERMUTATIONS,
) -> dict[str, Any]:
    core = [read for read in recovered if read.label is not None]
    labels: list[str] = []
    signatures: list[str] = []
    strata: list[str] = []
    complete_reads: list[RecoveredRead] = []
    for read in core:
        calls = [read.partner_calls[partner.pos] for partner in ranked_partners]
        if all(call in {"R", "A"} for call in calls):
            complete_reads.append(read)
            labels.append(str(read.label))
            signatures.append("".join(calls))
            strata.append(conditional_stratum(read))
    if len(ranked_partners) < 2:
        association = {
            "testable": False,
            "reason": "fewer_than_two_top_markers",
            "groups": sorted(set(labels)),
            "categories": sorted(set(signatures)),
            "table": [],
            "cramers_v": None,
        }
    else:
        association = LIB.categorical_association_table(labels, signatures)
    conditional = (
        LIB.conditional_label_permutation(
            labels,
            signatures,
            strata,
            seed=seed,
            permutations=permutations,
        )
        if association["testable"]
        else {
            "p_conditional_perm": None,
            "permutations": 0,
            "permutable": False,
            "status": "NOT_IDENTIFIABLE_JOINT_SIGNATURE_NOT_TESTABLE",
        }
    )
    sensitivity_pass = bool(
        conditional["permutable"]
        and conditional["p_conditional_perm"] is not None
        and float(conditional["p_conditional_perm"]) <= CONDITIONAL_P_MAX
    )
    marker_support = []
    effect_supported_positions = []
    for partner in ranked_partners:
        marker_calls = [read.partner_calls[partner.pos] for read in complete_reads]
        marker_association = LIB.methyl_group_marker_association(labels, marker_calls)
        effect_pass = bool(
            marker_association["testable"]
            and marker_association["cramers_v"] is not None
            and marker_association["delta_alt_fraction"] is not None
            and float(marker_association["cramers_v"]) >= ENDPOINT_A_V_MIN
            and float(marker_association["delta_alt_fraction"]) >= ENDPOINT_A_DELTA_MIN
        )
        marker_support.append(
            {
                "position": partner.pos,
                "testable": marker_association["testable"],
                "reason": marker_association["reason"],
                "n_informative": marker_association["n_informative"],
                "cramers_v": marker_association["cramers_v"],
                "delta_alt_fraction": marker_association["delta_alt_fraction"],
                "effect_gate_pass": effect_pass,
            }
        )
        if effect_pass:
            effect_supported_positions.append(partner.pos)
    return {
        "top_marker_positions": [partner.pos for partner in ranked_partners],
        "top_marker_signatures": dict(sorted(Counter(signatures).items())),
        "joint_signature_n_complete_reads": len(signatures),
        "joint_signature_testable": association["testable"],
        "joint_signature_reason": association["reason"],
        "joint_signature_groups": association["groups"],
        "joint_signature_categories": association["categories"],
        "joint_signature_table": association["table"],
        "joint_signature_cramers_v": association["cramers_v"],
        "joint_signature_conditional_strata": "latest_HP_family_x_PS_x_strand",
        "joint_signature_p_conditional_perm": conditional["p_conditional_perm"],
        "joint_signature_permutations": conditional["permutations"],
        "joint_signature_permutable": conditional["permutable"],
        "joint_signature_conditional_status": conditional["status"],
        "joint_signature_sensitivity_pass": sensitivity_pass,
        "joint_signature_postselection_fdr_calibrated": False,
        "joint_signature_claim_guardrail": (
            "Coverage-selected joint-signature conditional permutation is a postselection robustness gate, "
            "not a global-FDR-calibrated discovery test."
        ),
        "joint_signature_complete_marker_support": marker_support,
        "joint_signature_n_complete_marker_effect_supported": len(effect_supported_positions),
        "joint_signature_complete_marker_effect_supported_positions": effect_supported_positions,
        "joint_signature_complete_marker_support_contract": (
            "each contributing marker is R/A-testable and retains Cramers_V>=0.30 plus "
            "delta_ALT_fraction>=0.50 in the same complete-read set; the joint signature separately "
            "passes HP-family_x_PS_x_strand conditional sensitivity"
        ),
        "multi_marker_molecular_haplotype_base_candidate": False,
        # Deprecated aliases retained without analytic or global-FDR semantics.
        "signature_n_complete_reads": len(signatures),
        "signature_testable": association["testable"],
        "signature_reason": association["reason"],
        "signature_table": association["table"],
        "signature_cramers_v": association["cramers_v"],
        "signature_p_analytic": None,
        "signature_q_global_bh": None,
        "signature_q_global_by": None,
        "signature_by_sensitivity_pass": sensitivity_pass,
    }


def analyze_site_task(task: dict[str, Any]) -> dict[str, Any]:
    entry = task["entry"]
    assignment = task["assignment"]
    screen_row = task["screen_row"]
    vcf = get_vcf(entry["all_ssnv_vcf"]["path"], entry["all_ssnv_vcf_index"]["path"])
    focal, partners = discover_focal_and_partners(vcf, assignment["chrom"], int(assignment["pos"]))
    for field, observed in (("ref", focal.ref), ("alt", focal.alt)):
        screen_value = str(screen_row.get(field, observed)).upper()
        if screen_value != observed.upper():
            raise RuntimeError(
                f"Stable site/frozen VCF {field} conflict at {focal.chrom}:{focal.pos}: "
                f"{screen_value} != {observed}"
            )
    raw_identity_audit: dict[str, Any] = {}
    recovered = recover_site_reads(
        assignment,
        focal,
        partners,
        entry,
        raw_identity_audit,
    )
    validate_raw_identity_audit(
        raw_identity_audit,
        expected_recovered_count=len(recovered),
    )
    raw_identity_duplicate_rows = []
    for duplicate_record in raw_identity_audit["duplicate_projection_records"]:
        projection = duplicate_record["projection"]
        raw_identity_duplicate_rows.append(
            {
                "sample": entry["sample"],
                "chrom": focal.chrom,
                "pos": focal.pos,
                "projection_read_name": projection[0],
                "projection_chrom": projection[1],
                "projection_start": projection[2],
                "projection_end": projection[3],
                "projection_mapq": projection[4],
                "projection_strand": projection[5],
                "record_count": duplicate_record["record_count"],
                "classification": duplicate_record["classification"],
                "differing_auxiliary_tags": duplicate_record[
                    "differing_auxiliary_tags"
                ],
            }
        )
    observed_readset_digest = focal_alt_readset_digest(recovered)
    expected_readset_digest = str(screen_row.get("alt_readset_sha256") or "")
    if observed_readset_digest != expected_readset_digest:
        raise RuntimeError(
            f"Focal ALT readset digest mismatch at {entry['sample']}:{focal.chrom}:{focal.pos}: "
            f"raw={observed_readset_digest} site_tsv={expected_readset_digest}"
        )
    categorical_level_counts = m2_categorical_level_counts(assignment, screen_row)
    m2_gate = m2_screen_eligibility(screen_row, categorical_level_counts)
    m2_eligible = bool(m2_gate["eligible"])
    pair_rows = [
        _pair_result(
            entry,
            focal,
            partner,
            recovered,
            m2_eligible=m2_eligible,
            exact_state_space_ceiling=int(
                task.get("exact_state_space_ceiling", EXACT_STATE_SPACE_CEILING)
            ),
        )
        for partner in partners
    ]
    pair_by_position = {int(row["partner_pos"]): row for row in pair_rows}
    ranked = sorted(
        partners,
        key=lambda partner: (
            -int(pair_by_position[partner.pos]["endpoint_a_n_informative"]),
            abs(partner.pos - focal.pos),
            partner.pos,
            partner.ref,
            partner.alt,
        ),
    )[: int(task["top_markers"])]
    for rank, partner in enumerate(ranked, start=1):
        pair_by_position[partner.pos]["top_coverage_rank"] = rank
        pair_by_position[partner.pos]["top_coverage_marker"] = True
    n_ref = sum(read.focal_call == "R" for read in recovered)
    n_alt = sum(read.focal_call == "A" for read in recovered)
    core_count = sum(read.label is not None for read in recovered)
    site_row: dict[str, Any] = {
        "sample": entry["sample"],
        "biological_id": entry["biological_id"],
        "chrom": focal.chrom,
        "pos": focal.pos,
        "ref": focal.ref,
        "alt": focal.alt,
        "truth_label": None,
        "ssnv_branch": None,
        "component_id": None,
        "selected_group_id": None,
        "alt_readset_sha256": screen_row.get("alt_readset_sha256"),
        "region_dir": assignment["region_dir"],
        "screen_contract": screen_row["screen_contract"],
        "stable_null_multigroup": parse_bool(screen_row.get("stable_null_multigroup", False)),
        "modal_assignment_ari_min": (
            float(screen_row["modal_assignment_ari_min"])
            if screen_row.get("modal_assignment_ari_min") not in {None, ""}
            else None
        ),
        "hp_axis_confound": parse_bool(screen_row.get("hp_axis_confound", True)),
        "technical_axis_confound": parse_bool(screen_row.get("technical_axis_confound", True)),
        "residual_unexplained_multigroup": parse_bool(
            screen_row.get("residual_unexplained_multigroup", False)
        ),
        "phase_anchored_robust_epigenetic_candidate": parse_bool(
            screen_row.get("phase_anchored_robust_epigenetic_candidate", False)
        ),
        "m2_screen_gate_contract": m2_gate["contract"],
        "m2_screen_evaluable": bool(m2_gate["evaluable"]),
        "m2_screen_eligible": m2_eligible,
        "m2_screen_eligibility_status": m2_gate["status"],
        "m2_categorical_level_counts": categorical_level_counts,
        "m2_axis_statuses": m2_gate["axis_statuses"],
        "m2_indeterminate_axes": m2_gate["indeterminate_axes"],
        "m2_low_power_axes": m2_gate["low_power_axes"],
        "m2_aligned_axes": m2_gate["aligned_axes"],
        "m2_constant_axes": m2_gate["constant_axes"],
        "m2_aligned_below_negative_evaluability_power_axes": m2_gate[
            "aligned_below_negative_evaluability_power_axes"
        ],
        "n_all_focal_ref_alt_reads": len(recovered),
        "n_focal_ref_reads": n_ref,
        "n_focal_alt_reads": n_alt,
        "n_core_focal_alt_reads": core_count,
        "raw_identity_expected_projections": int(raw_identity_audit["expected_projections"]),
        "raw_identity_matched_records": int(raw_identity_audit["matched_records"]),
        "raw_identity_duplicate_projections_collapsed": int(
            raw_identity_audit["duplicate_projections_collapsed"]
        ),
        "raw_identity_duplicate_extra_records_collapsed": int(
            raw_identity_audit["duplicate_extra_records_collapsed"]
        ),
        "raw_identity_exact_duplicate_projections_collapsed": int(
            raw_identity_audit["exact_duplicate_projections_collapsed"]
        ),
        "raw_identity_rg_only_duplicate_projections_collapsed": int(
            raw_identity_audit["rg_only_duplicate_projections_collapsed"]
        ),
        "raw_identity_duplicate_projection_sha256": str(
            raw_identity_audit["duplicate_projection_sha256"]
        ),
        "n_partner_markers": len(partners),
        "n_pair_rows_reconciled": 0,
        "pair_row_count_reconciliation_pass": False,
        "n_endpoint_a_testable_markers": 0,
        "n_endpoint_a_exact_identifiable_markers": 0,
        "n_endpoint_a_exact_not_identifiable_markers": 0,
        "n_endpoint_a_conditional_permutable_markers": 0,
        "n_pair_bh_discoveries": 0,
        "pair_bh_discovery_positions": [],
        "n_pair_by_discoveries": 0,
        "pair_by_discovery_positions": [],
        "n_pair_by_confirmed": 0,
        "pair_by_confirmed_positions": [],
        "n_spatially_separated_pair_by_20bp": 0,
        "spatially_separated_pair_by_positions_20bp": [],
        "spaced_marker_20bp_contract": (
            "genomic_distance_at_least_20bp_only_not_statistical_independence"
        ),
        "n_endpoint_a_pre_candidates": 0,
        "n_confirmed_markers": 0,
        "confirmed_marker_positions": [],
        "n_independent_confirmed_markers_20bp": 0,
        "genetically_anchored_multi_marker_candidate": False,
        "n_confirmed_markers_by_sensitivity": 0,
        "confirmed_marker_positions_by_sensitivity": [],
        "n_independent_confirmed_markers_20bp_by_sensitivity": 0,
        "genetically_anchored_multi_marker_candidate_by_sensitivity": False,
        "n_top_marker_pair_by_confirmed": 0,
        "top_marker_pair_by_confirmed_positions": [],
        "n_same_complete_read_effect_supported_top_pair_by": 0,
        "same_complete_read_effect_supported_top_pair_by_positions": [],
        "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp": 0,
        "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp": [],
        **_signature_result(
            recovered,
            ranked,
            seed=deterministic_seed(
                entry["sample"], focal.chrom, focal.pos, "joint_signature_conditional"
            ),
            permutations=MIN_PERMUTATIONS,
        ),
        "joint_signature_global_fdr_family_status": "NOT_RUN_GLOBAL_ADJUSTMENT_PENDING",
        "joint_signature_q_global_bh": None,
        "joint_signature_q_global_by": None,
        "joint_signature_global_bh_discovery": False,
        "joint_signature_global_by_discovery": False,
        "biological_site_key": None,
        "biological_site_representative": False,
        "component_dedup_key": None,
        "component_representative": False,
        "alt_readset_dedup_key": None,
        "alt_readset_representative": False,
    }
    validate_raw_identity_site_fields(
        site_row,
        expected_recovered_count=len(recovered),
    )
    return {
        "site": site_row,
        "pairs": pair_rows,
        "raw_identity_duplicates": raw_identity_duplicate_rows,
    }


def analyze_task_chunk(tasks: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    return [analyze_site_task(task) for task in tasks]


def chunked(values: Sequence[Any], size: int) -> Iterator[list[Any]]:
    if size < 1:
        raise ValueError("chunk size must be positive")
    for offset in range(0, len(values), size):
        yield list(values[offset : offset + size])


def bounded_process_map(
    tasks: Sequence[dict[str, Any]],
    workers: int,
    chunk_size: int,
    max_pending_chunks: int,
) -> Iterator[dict[str, Any]]:
    chunks = iter(chunked(tasks, chunk_size))
    if workers == 1:
        for task_chunk in chunks:
            yield from analyze_task_chunk(task_chunk)
        return
    if max_pending_chunks < workers:
        raise ValueError("max_pending_chunks must be at least workers")
    with ProcessPoolExecutor(max_workers=workers) as executor:
        pending: dict[Any, int] = {}
        next_index = 0

        def submit_one() -> bool:
            nonlocal next_index
            try:
                task_chunk = next(chunks)
            except StopIteration:
                return False
            pending[executor.submit(analyze_task_chunk, task_chunk)] = next_index
            next_index += 1
            return True

        while len(pending) < max_pending_chunks and submit_one():
            pass
        completed_by_index: dict[int, list[dict[str, Any]]] = {}
        emit_index = 0
        try:
            while pending:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    index = pending.pop(future)
                    completed_by_index[index] = future.result()
                    submit_one()
                while emit_index in completed_by_index:
                    yield from completed_by_index.pop(emit_index)
                    emit_index += 1
        except Exception:
            for future in pending:
                future.cancel()
            raise


def validate_completed_sample_run(sample_root: Path, entry: dict[str, Any]) -> dict[str, Any]:
    receipt_path = sample_root / "run_receipt.json"
    if not receipt_path.is_file() or receipt_path.stat().st_size <= 0:
        raise FileNotFoundError(f"Completed InterSubMod run receipt missing: {receipt_path}")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    validation = receipt.get("validation") or {}
    expected = int(entry["counts"]["all_ssnv"])
    exact_counts = {
        "expected_vcf_sites": expected,
        "reads_files": expected,
        "bernoulli_matrix_files": expected,
        "methylation_files": expected,
    }
    for field, expected_value in exact_counts.items():
        if int(validation.get(field, -1)) != expected_value:
            raise RuntimeError(f"{entry['sample']} run receipt {field} mismatch")
    summary_rows = int(validation.get("summary_rows", -1))
    if not 0 < summary_rows <= expected:
        raise RuntimeError(f"{entry['sample']} invalid summary_rows in run receipt")
    gates = {
        "schema": receipt.get("schema_name") == "intersubmod.all_ssnv_site_run",
        "sample": receipt.get("sample") == entry["sample"],
        "pass": receipt.get("pass") is True,
        "exit_code": receipt.get("exit_code") == 0,
        "validation_pass": validation.get("pass") is True,
        "vcf_sha256": receipt.get("input_vcf_sha256") == entry["all_ssnv_vcf"].get("sha256"),
        "output_dir": Path(receipt.get("output_dir", "")).resolve() == sample_root.resolve(),
        "sidecar": Path((receipt.get("latest_read_tag_sidecar") or {}).get("path", "")).resolve()
        == Path(entry["latest_read_tag_sidecar"]["path"]).resolve(),
    }
    failed = sorted(name for name, passed in gates.items() if not passed)
    if failed:
        raise RuntimeError(f"{entry['sample']} completed-run receipt failed gates: {failed}")
    return {"receipt": artifact(receipt_path), "validation": validation}


def load_manifest(
    path: Path,
    intersubmod_root: Path,
    requested_samples: set[str] | None,
    allow_partial_scope: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if (
        manifest.get("schema_name") != MANIFEST_SCHEMA
        or manifest.get("schema_version") != "1.0.0"
        or manifest.get("pass") is not True
    ):
        raise RuntimeError("Input manifest schema/pass gate failed")
    all_entries = manifest.get("samples")
    if not isinstance(all_entries, list) or not all_entries:
        raise RuntimeError("Input manifest has no samples")
    known = {entry["sample"] for entry in all_entries}
    if requested_samples:
        unknown = sorted(requested_samples.difference(known))
        if unknown:
            raise ValueError(f"Unknown requested samples: {unknown}")
        if requested_samples != known and not allow_partial_scope:
            raise RuntimeError("Subset execution requires --allow-partial-scope")
        entries = [entry for entry in all_entries if entry["sample"] in requested_samples]
    else:
        entries = list(all_entries)
    run_validation: dict[str, Any] = {}
    for entry in entries:
        identity_rows = [
            verify_manifest_frozen_field(field, entry[field])
            for field in MANIFEST_FROZEN_FIELDS
        ]
        entry["_frozen_artifact_identity_audit"] = identity_rows
        sample_root = (intersubmod_root / entry["sample"]).resolve()
        run_validation[entry["sample"]] = validate_completed_sample_run(sample_root, entry)
        entry["_validated_intersubmod_sample_root"] = str(sample_root)
    return manifest, entries, run_validation


def load_assignments(path: Path, selected_samples: set[str]) -> dict[tuple[str, str, int], dict[str, Any]]:
    assignments: dict[tuple[str, str, int], dict[str, Any]] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("schema_name") != ASSIGNMENT_SCHEMA or row.get("schema_version") != "1.0.0":
                raise RuntimeError(f"Unexpected assignment schema at line {line_number}")
            if row.get("screen_contract") != SCREEN_CONTRACT:
                raise RuntimeError(f"Unexpected assignment screen contract at line {line_number}")
            if row.get("artifact_identity_contract") != "sha256_size_path_v1":
                raise RuntimeError(f"Unexpected assignment artifact identity contract at line {line_number}")
            if row.get("strict_confirm_candidate") is not True:
                raise RuntimeError(f"Non-stable assignment emitted at line {line_number}")
            if row["sample"] not in selected_samples:
                continue
            key = site_key(row["sample"], row["chrom"], row["pos"])
            if key in assignments:
                raise RuntimeError(f"Duplicate stable assignment: {key}")
            assignments[key] = row
    if not assignments:
        raise RuntimeError("No stable assignments selected")
    return assignments


def load_stable_site_rows(
    path: Path,
    assignments: dict[tuple[str, str, int], dict[str, Any]],
    selected_samples: set[str] | None = None,
    require_candidate_screen_fields: bool = False,
) -> dict[tuple[str, str, int], dict[str, str]]:
    if selected_samples is None:
        selected_samples = {key[0] for key in assignments}
    rows: dict[tuple[str, str, int], dict[str, str]] = {}
    missing_assignments: list[tuple[str, str, int]] = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "sample",
            "chrom",
            "pos",
            "stable_null_multigroup",
            "ref",
            "alt",
            "alt_readset_sha256",
            "screen_contract",
        }
        if require_candidate_screen_fields:
            required.update(
                {
                    "cluster_sizes",
                    "modal_assignment_ari_min",
                    "hp_axis_confound",
                    "technical_axis_confound",
                    "residual_unexplained_multigroup",
                    "phase_anchored_robust_epigenetic_candidate",
                    *M2_GATE.AXIS_SOURCE_FIELDS,
                }
            )
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise RuntimeError(f"Stable site TSV missing columns: {missing}")
        for row in reader:
            if row["sample"] not in selected_samples:
                continue
            key = site_key(row["sample"], row["chrom"], row["pos"])
            stable = parse_bool(row["stable_null_multigroup"])
            if stable and key not in assignments:
                missing_assignments.append(key)
                continue
            if key not in assignments:
                continue
            if not stable:
                raise RuntimeError(f"Assignment/site stable flag conflict: {key}")
            if row["screen_contract"] != SCREEN_CONTRACT:
                raise RuntimeError(f"Stable site screen contract conflict: {key}")
            if assignment_contract := assignments[key].get("screen_contract"):
                if assignment_contract != row["screen_contract"]:
                    raise RuntimeError(f"Assignment/site screen contract conflict: {key}")
            if key in rows:
                raise RuntimeError(f"Duplicate stable site TSV row: {key}")
            rows[key] = row
    missing_rows = sorted(set(assignments).difference(rows))
    if missing_rows:
        raise RuntimeError(f"Stable assignments missing site TSV rows: {missing_rows[:5]}")
    if missing_assignments:
        raise RuntimeError(
            f"Stable site TSV rows missing read assignments: {sorted(missing_assignments)[:5]}"
        )
    return rows


def build_tasks(
    entries: Sequence[dict[str, Any]],
    assignments: dict[tuple[str, str, int], dict[str, Any]],
    site_rows: dict[tuple[str, str, int], dict[str, str]],
    top_markers: int,
    exact_state_space_ceiling: int = EXACT_STATE_SPACE_CEILING,
) -> list[dict[str, Any]]:
    entries_by_sample = {entry["sample"]: entry for entry in entries}
    tasks = []
    for key in sorted(assignments):
        assignment = assignments[key]
        if assignment["sample"] not in entries_by_sample:
            raise RuntimeError(f"Assignment sample absent from selected manifest: {assignment['sample']}")
        screen_region = site_rows[key].get("region_dir")
        if screen_region and Path(screen_region).resolve() != Path(assignment["region_dir"]).resolve():
            raise RuntimeError(f"Stable assignment/site region_dir conflict: {key}")
        validated_root = entries_by_sample[assignment["sample"]].get(
            "_validated_intersubmod_sample_root"
        )
        if validated_root:
            try:
                Path(assignment["region_dir"]).resolve().relative_to(Path(validated_root))
            except ValueError as error:
                raise RuntimeError(
                    f"Stable assignment region_dir is outside validated v2 sample root: {key}"
                ) from error
        worker_assignment = dict(assignment)
        worker_assignment.pop("posthoc", None)
        worker_screen_row = {
            "region_dir": site_rows[key].get("region_dir"),
            "alt_readset_sha256": site_rows[key]["alt_readset_sha256"],
            "ref": site_rows[key]["ref"],
            "alt": site_rows[key]["alt"],
            "screen_contract": site_rows[key]["screen_contract"],
        }
        for field in (
            "stable_null_multigroup",
            "cluster_sizes",
            "modal_assignment_ari_min",
            "hp_axis_confound",
            "technical_axis_confound",
            "residual_unexplained_multigroup",
            "phase_anchored_robust_epigenetic_candidate",
            *M2_GATE.AXIS_SOURCE_FIELDS,
        ):
            if field in site_rows[key]:
                worker_screen_row[field] = site_rows[key][field]
        source_entry = entries_by_sample[assignment["sample"]]
        worker_entry = {
            field: source_entry[field]
            for field in (
                "sample",
                "biological_id",
                "all_ssnv_vcf",
                "all_ssnv_vcf_index",
                "latest_read_tag_sidecar",
                "latest_read_tag_sidecar_index",
                "raw_alignment",
                "raw_alignment_index",
            )
            if field in source_entry
        }
        tasks.append(
            {
                "entry": worker_entry,
                "assignment": worker_assignment,
                "screen_row": worker_screen_row,
                "top_markers": top_markers,
                "exact_state_space_ceiling": exact_state_space_ceiling,
            }
        )
    return tasks


def _conditional_sensitivity_task(task: dict[str, Any]) -> dict[str, Any]:
    payload = task["payload"]
    return {
        "index": task["index"],
        "result": LIB.conditional_label_permutation(
            payload["labels"],
            payload["categories"],
            payload["strata"],
            seed=int(task["seed"]),
            permutations=int(task["permutations"]),
        ),
    }


def _conditional_sensitivity_chunk(tasks: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    return [_conditional_sensitivity_task(task) for task in tasks]


def run_conditional_permutations(
    pair_rows: list[dict[str, Any]],
    permutations: int,
    workers: int = 1,
    chunk_size: int = 16,
) -> dict[str, Any]:
    """Run secondary sensitivity only for global exact-BY effect discoveries."""
    tasks: list[dict[str, Any]] = []
    for index, row in enumerate(pair_rows):
        payload = row["_permutation_payload"]
        row["_permutation_payload"] = {
            "core_raw_read_names": payload["core_raw_read_names"]
        }
        if not row["endpoint_a_exact_by_discovery"]:
            continue
        tasks.append(
            {
                "index": index,
                "payload": payload,
                "permutations": permutations,
                "seed": deterministic_seed(
                    row["sample"],
                    row["focal_chrom"],
                    row["focal_pos"],
                    row["focal_ref"],
                    row["focal_alt"],
                    row["partner_chrom"],
                    row["partner_pos"],
                    row["partner_ref"],
                    row["partner_alt"],
                    "conditional_sensitivity",
                ),
            }
        )
    results: list[dict[str, Any]] = []
    task_chunks = list(chunked(tasks, chunk_size))
    if workers <= 1 or len(task_chunks) <= 1:
        for task_chunk in task_chunks:
            results.extend(_conditional_sensitivity_chunk(task_chunk))
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            for chunk_results in executor.map(_conditional_sensitivity_chunk, task_chunks):
                results.extend(chunk_results)
    seen: set[int] = set()
    for payload in results:
        index = int(payload["index"])
        if index in seen:
            raise RuntimeError(f"Duplicate conditional sensitivity result index: {index}")
        seen.add(index)
        row = pair_rows[index]
        result = payload["result"]
        row["endpoint_a_p_conditional_perm"] = result["p_conditional_perm"]
        row["endpoint_a_permutations"] = result["permutations"]
        row["endpoint_a_permutable"] = result["permutable"]
        row["endpoint_a_conditional_status"] = result["status"]
        row["endpoint_a_conditional_sensitivity_pass"] = bool(
            result["permutable"]
            and result["p_conditional_perm"] is not None
            and float(result["p_conditional_perm"]) <= CONDITIONAL_P_MAX
        )
        row["endpoint_a_formal_pair_by_confirmed"] = bool(
            row["endpoint_a_exact_by_discovery"]
            and row["endpoint_a_effect_gate_pass"]
            and row["endpoint_a_conditional_sensitivity_pass"]
            and row["callability_gate_pass"]
        )
        # Deprecated aliases. Formal G1 uses endpoint_a_formal_pair_by_confirmed.
        row["endpoint_a_confirmed_association"] = row["endpoint_a_formal_pair_by_confirmed"]
        row["endpoint_a_confirmed_by_sensitivity"] = row[
            "endpoint_a_formal_pair_by_confirmed"
        ]
    expected = {int(task["index"]) for task in tasks}
    if seen != expected:
        raise RuntimeError("Conditional sensitivity worker result reconciliation failed")
    return {
        "contract": "secondary_only_after_global_fixed_margins_exact_BY_discovery",
        "n_exact_by_discoveries_submitted": len(tasks),
        "n_conditional_permutable": sum(
            row["endpoint_a_permutable"] for row in pair_rows if row["endpoint_a_exact_by_discovery"]
        ),
        "n_conditional_sensitivity_pass": sum(
            row["endpoint_a_conditional_sensitivity_pass"] for row in pair_rows
        ),
        "permutations": permutations,
        "global_fdr_calibrated": False,
    }


def finalize_site_inference(
    site_rows: list[dict[str, Any]], pair_rows: list[dict[str, Any]]
) -> dict[str, Any]:
    by_site: dict[tuple[str, str, int], list[dict[str, Any]]] = defaultdict(list)
    for pair in pair_rows:
        by_site[site_key(pair["sample"], pair["focal_chrom"], pair["focal_pos"])].append(pair)
    reconciled_pair_rows = 0
    eligible_sites = 0
    eligible_family_pairs = 0
    for site in site_rows:
        pairs = by_site.get(site_key(site["sample"], site["chrom"], site["pos"]), [])
        formal_contract = "n_partner_markers" in site
        declared_pairs = int(site.get("n_partner_markers", len(pairs)))
        if len(pairs) != declared_pairs:
            raise RuntimeError(
                f"Per-focal pair reconciliation failed for {site['sample']}:{site['chrom']}:{site['pos']}: "
                f"site={declared_pairs} pair_rows={len(pairs)}"
            )
        reconciled_pair_rows += len(pairs)
        m2_eligible = bool(site.get("m2_screen_eligible", True))
        eligible_sites += m2_eligible
        eligible_family_pairs += sum(
            pair.get("endpoint_a_global_fdr_family_status") == "ELIGIBLE_M2_EXACT_FAMILY"
            for pair in pairs
        )
        bh_positions = sorted(
            int(pair["partner_pos"])
            for pair in pairs
            if pair.get("endpoint_a_exact_bh_discovery", pair.get("endpoint_a_pre_candidate", False))
        )
        by_positions = sorted(
            int(pair["partner_pos"])
            for pair in pairs
            if pair.get(
                "endpoint_a_exact_by_discovery",
                pair.get("endpoint_a_pre_candidate_by_sensitivity", False),
            )
        )
        confirmed_positions = sorted(
            int(pair["partner_pos"])
            for pair in pairs
            if pair.get(
                "endpoint_a_formal_pair_by_confirmed",
                pair.get("endpoint_a_confirmed_association", False),
            )
        )
        spaced_positions = LIB.spatially_separated_positions(
            confirmed_positions, minimum_separation=SPACED_MARKER_MIN_SEPARATION_BP
        )
        top_positions = {
            int(value) for value in site.get("top_marker_positions", confirmed_positions)
        }
        top_confirmed_positions = sorted(top_positions.intersection(confirmed_positions))
        top_spaced_positions = LIB.spatially_separated_positions(
            top_confirmed_positions, minimum_separation=SPACED_MARKER_MIN_SEPARATION_BP
        )
        complete_effect_positions = {
            int(value)
            for value in site.get(
                "joint_signature_complete_marker_effect_supported_positions", []
            )
        }
        same_complete_supported_positions = sorted(
            set(top_confirmed_positions).intersection(complete_effect_positions)
        )
        same_complete_spaced_positions = LIB.spatially_separated_positions(
            same_complete_supported_positions,
            minimum_separation=SPACED_MARKER_MIN_SEPARATION_BP,
        )
        candidate = bool(
            formal_contract
            and m2_eligible
            and len(spaced_positions) >= 2
            and len(top_spaced_positions) >= 2
            and len(same_complete_spaced_positions) >= 2
            and site.get("joint_signature_global_by_discovery", False)
        )
        legacy_candidate = candidate if formal_contract else len(spaced_positions) >= 2
        site.update(
            {
                "n_pair_rows_reconciled": len(pairs),
                "pair_row_count_reconciliation_pass": True,
                "n_endpoint_a_testable_markers": sum(
                    pair.get("endpoint_a_testable", False) for pair in pairs
                ),
                "n_endpoint_a_exact_identifiable_markers": sum(
                    pair.get("endpoint_a_p_fixed_margins_exact") is not None for pair in pairs
                ),
                "n_endpoint_a_exact_not_identifiable_markers": sum(
                    pair.get("endpoint_a_testable", False)
                    and pair.get("endpoint_a_p_fixed_margins_exact") is None
                    for pair in pairs
                ),
                "n_endpoint_a_conditional_permutable_markers": sum(
                    pair.get("endpoint_a_permutable", False) for pair in pairs
                ),
                "n_pair_bh_discoveries": len(bh_positions),
                "pair_bh_discovery_positions": bh_positions,
                "n_pair_by_discoveries": len(by_positions),
                "pair_by_discovery_positions": by_positions,
                "n_pair_by_confirmed": len(confirmed_positions),
                "pair_by_confirmed_positions": confirmed_positions,
                "n_spatially_separated_pair_by_20bp": len(spaced_positions),
                "spatially_separated_pair_by_positions_20bp": spaced_positions,
                "spaced_marker_20bp_contract": (
                    "genomic_distance_at_least_20bp_only_not_statistical_independence"
                ),
                "n_top_marker_pair_by_confirmed": len(top_confirmed_positions),
                "top_marker_pair_by_confirmed_positions": top_confirmed_positions,
                "n_same_complete_read_effect_supported_top_pair_by": len(
                    same_complete_supported_positions
                ),
                "same_complete_read_effect_supported_top_pair_by_positions": (
                    same_complete_supported_positions
                ),
                "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp": len(
                    same_complete_spaced_positions
                ),
                "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp": (
                    same_complete_spaced_positions
                ),
                "multi_marker_molecular_haplotype_base_candidate": candidate,
                "n_endpoint_a_pre_candidates": len(bh_positions),
                # Deprecated aliases; never describe these markers as statistically independent.
                "n_confirmed_markers": len(confirmed_positions),
                "confirmed_marker_positions": confirmed_positions,
                "n_independent_confirmed_markers_20bp": len(spaced_positions),
                "genetically_anchored_multi_marker_candidate": legacy_candidate,
                "n_confirmed_markers_by_sensitivity": len(confirmed_positions),
                "confirmed_marker_positions_by_sensitivity": confirmed_positions,
                "n_independent_confirmed_markers_20bp_by_sensitivity": len(spaced_positions),
                "genetically_anchored_multi_marker_candidate_by_sensitivity": legacy_candidate,
            }
        )
    if reconciled_pair_rows != len(pair_rows):
        raise RuntimeError(
            f"Aggregate pair reconciliation failed: sites={reconciled_pair_rows} rows={len(pair_rows)}"
        )
    return {
        "contract": "every_site_count_recomputed_from_exact_directed_pair_rows",
        "pass": True,
        "n_sites_reconciled": len(site_rows),
        "n_sites_with_zero_pairs": sum(int(site.get("n_partner_markers", 0)) == 0 for site in site_rows),
        "n_pair_rows_reconciled": reconciled_pair_rows,
        "n_m2_global_fdr_eligible_sites": eligible_sites,
        "n_m2_exact_global_fdr_family_pairs": eligible_family_pairs,
    }


def _ledger_header(path: str) -> list[str]:
    if path not in _LEDGER_HEADERS:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            _LEDGER_HEADERS[path] = handle.readline().rstrip("\n").split("\t")
    return _LEDGER_HEADERS[path]


def fetch_ledger_row(entry: dict[str, Any], variant: LIB.Variant) -> dict[str, str]:
    path = entry["site_ledger"]["path"]
    cache_key = (path, variant.chrom, variant.pos, variant.ref.upper(), variant.alt.upper())
    if cache_key in _LEDGER_ROWS:
        return _LEDGER_ROWS[cache_key]
    tabix = get_tabix(path, entry["site_ledger_index"]["path"])
    header = _ledger_header(path)
    matches = []
    for line in tabix.fetch(variant.chrom, variant.pos - 1, variant.pos):
        values = line.rstrip("\n").split("\t")
        if len(values) != len(header):
            raise RuntimeError(f"Malformed ledger row at {variant.chrom}:{variant.pos}")
        row = dict(zip(header, values))
        if variant_key(row["chrom"], row["pos"], row["ref"], row["alt"]) == variant_key(
            variant.chrom, variant.pos, variant.ref, variant.alt
        ):
            matches.append(row)
    if len(matches) != 1:
        raise RuntimeError(
            f"Ledger exact join must be unique for {variant.chrom}:{variant.pos}:{variant.ref}>{variant.alt}; "
            f"found={len(matches)}"
        )
    _LEDGER_ROWS[cache_key] = matches[0]
    return matches[0]


def _vcf_has_exact(vcf: pysam.VariantFile, variant: LIB.Variant) -> bool:
    matches = 0
    for record in vcf.fetch(variant.chrom, variant.pos - 1, variant.pos):
        for alt in record.alts or ():
            if variant_key(record.chrom, record.pos, record.ref, alt) == variant_key(
                variant.chrom, variant.pos, variant.ref, variant.alt
            ):
                matches += 1
    if matches > 1:
        raise RuntimeError(f"Duplicate exact truth record: {variant}")
    return matches == 1


def truth_label(entry: dict[str, Any], variant: LIB.Variant) -> str:
    key = (
        entry["sample"],
        entry["truth_tp_vcf"]["path"],
        variant.chrom,
        variant.pos,
        variant.ref.upper(),
        variant.alt.upper(),
    )
    if key in _TRUTH_LABELS:
        return _TRUTH_LABELS[key]
    tp = _vcf_has_exact(
        get_vcf(entry["truth_tp_vcf"]["path"], entry["truth_tp_vcf_index"]["path"]), variant
    )
    fp = _vcf_has_exact(
        get_vcf(entry["truth_fp_vcf"]["path"], entry["truth_fp_vcf_index"]["path"]), variant
    )
    if tp and fp:
        raise RuntimeError(f"Truth TP/FP overlap at {variant}")
    label = "TP" if tp else "FP" if fp else "UNASSESSED"
    _TRUTH_LABELS[key] = label
    return label


def topology_direct_scope(
    focal_branch: str,
    partner_branch: str,
    focal_pos: int,
    partner_pos: int,
    region: dict[str, Any] | None,
) -> str:
    if "max_snv_excluded" in {focal_branch, partner_branch}:
        return "INDIRECT_MAX_SNV_EXCLUDED"
    if focal_branch != "retained" or partner_branch != "retained":
        return "INDIRECT_NON_RETAINED"
    if region is None:
        return "OUTSIDE_SHARED_LAYERED_REGION"
    positions = {int(value) for value in region.get("positions", [])}
    if focal_pos not in positions or partner_pos not in positions:
        return "OUTSIDE_REGION_POSITIONS"
    return "DIRECT_RETAINED_SHARED_REGION"


def _node_state(node: str, n_positions: int) -> str | None:
    if node == "ROOT":
        return "R" * n_positions
    state = node.rsplit("_", 1)[-1]
    if len(state) != n_positions or any(value not in {"R", "A"} for value in state):
        return None
    return state


def _tree_pair_order(tree: dict[str, Any], focal_index: int, partner_index: int, n_positions: int) -> str | None:
    if tree.get("recurrence"):
        return None
    parents: dict[str, str] = {}
    focal_transitions: list[str] = []
    partner_transitions: list[str] = []
    for edge in tree.get("edges", []):
        if not isinstance(edge, list) or len(edge) != 2:
            return None
        parent, child = str(edge[0]), str(edge[1])
        if child in parents and parents[child] != parent:
            return None
        parents[child] = parent
        parent_state = _node_state(parent, n_positions)
        child_state = _node_state(child, n_positions)
        if parent_state is None or child_state is None:
            return None
        if parent_state[focal_index] == "R" and child_state[focal_index] == "A":
            focal_transitions.append(child)
        if parent_state[partner_index] == "R" and child_state[partner_index] == "A":
            partner_transitions.append(child)
    if len(focal_transitions) != 1 or len(partner_transitions) != 1:
        return None
    focal_node = focal_transitions[0]
    partner_node = partner_transitions[0]
    if focal_node == partner_node:
        return "COINCIDENT_TRANSITION"

    def ancestors(node: str) -> set[str]:
        result = set()
        while node in parents:
            node = parents[node]
            if node in result:
                return set()
            result.add(node)
        return result

    if focal_node in ancestors(partner_node):
        return "FOCAL_BEFORE_PARTNER"
    if partner_node in ancestors(focal_node):
        return "PARTNER_BEFORE_FOCAL"
    return "BRANCHING_RELATION"


def infer_topology_order(region: dict[str, Any], focal_pos: int, partner_pos: int) -> str:
    positions = [int(value) for value in region.get("positions", [])]
    if focal_pos not in positions or partner_pos not in positions:
        return "NOT_INFERABLE_OUTSIDE_REGION_POSITIONS"
    focal_index = positions.index(focal_pos)
    partner_index = positions.index(partner_pos)
    relevant = []
    for lineage in region.get("lineages", []):
        coverage = lineage.get("obs_col_coverage") or {}
        if not lineage.get("is_primary_lineage"):
            continue
        if str(focal_pos) in coverage and str(partner_pos) in coverage:
            relevant.append(lineage)
    if not relevant:
        return "NOT_INFERABLE_NO_SHARED_PRIMARY_LINEAGE"
    orders = []
    for lineage in relevant:
        if lineage.get("capped"):
            return "NOT_INFERABLE_CAPPED"
        if lineage.get("analysis_candidate_set_complete") is not True:
            return "NOT_INFERABLE_INCOMPLETE_CANDIDATE_SET"
        if int(lineage.get("n_trees") or 0) != 1:
            return "NOT_INFERABLE_MULTI_TREE"
        trees = lineage.get("trees") or []
        if len(trees) != 1:
            return "NOT_INFERABLE_TREE_STORAGE_MISMATCH"
        order = _tree_pair_order(trees[0], focal_index, partner_index, len(positions))
        if order is None:
            return "NOT_INFERABLE_TREE_ENCODING_OR_RECURRENCE"
        orders.append(order)
    if len(set(orders)) != 1:
        return "NOT_INFERABLE_LINEAGE_DISAGREEMENT"
    return orders[0]


def _region_index(path: Path, required: set[tuple[str, int]]) -> dict[tuple[str, int], dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    index: dict[tuple[str, int], dict[str, Any]] = {}
    for region in payload.get("regions", []):
        chrom = str(region["chrom"])
        for position in region.get("positions", []):
            key = (chrom, int(position))
            if key not in required:
                continue
            if key in index and index[key].get("region") != region.get("region"):
                raise RuntimeError(f"Topology position belongs to multiple regions: {key}")
            index[key] = region
    return index


def attach_posthoc_annotations(
    entries: Sequence[dict[str, Any]],
    site_rows: list[dict[str, Any]],
    pair_rows: list[dict[str, Any]],
    source_site_rows: dict[tuple[str, str, int], dict[str, str]] | None = None,
    source_assignments: dict[tuple[str, str, int], dict[str, Any]] | None = None,
) -> None:
    entries_by_sample = {entry["sample"]: entry for entry in entries}
    required_by_sample: dict[str, set[tuple[str, int]]] = defaultdict(set)
    for pair in pair_rows:
        required_by_sample[pair["sample"]].update(
            {
                (pair["focal_chrom"], int(pair["focal_pos"])),
                (pair["partner_chrom"], int(pair["partner_pos"])),
            }
        )
    region_indexes = {
        sample: _region_index(Path(entries_by_sample[sample]["layered_region_view"]["path"]), required)
        for sample, required in required_by_sample.items()
    }
    site_by_key = {site_key(row["sample"], row["chrom"], row["pos"]): row for row in site_rows}
    for pair in pair_rows:
        entry = entries_by_sample[pair["sample"]]
        focal = LIB.Variant(pair["focal_chrom"], int(pair["focal_pos"]), pair["focal_ref"], pair["focal_alt"])
        partner = LIB.Variant(
            pair["partner_chrom"], int(pair["partner_pos"]), pair["partner_ref"], pair["partner_alt"]
        )
        focal_ledger = fetch_ledger_row(entry, focal)
        partner_ledger = fetch_ledger_row(entry, partner)
        pair.update(
            {
                "focal_truth_label": truth_label(entry, focal),
                "partner_truth_label": truth_label(entry, partner),
                "focal_ssnv_branch": focal_ledger["ssnv_branch"],
                "partner_ssnv_branch": partner_ledger["ssnv_branch"],
                "focal_component_id": focal_ledger.get("component_id"),
                "partner_component_id": partner_ledger.get("component_id"),
            }
        )
        screen_site = site_by_key[site_key(pair["sample"], pair["focal_chrom"], pair["focal_pos"])]
        for field, value in (
            ("truth_label", pair["focal_truth_label"]),
            ("ssnv_branch", pair["focal_ssnv_branch"]),
            ("component_id", pair["focal_component_id"]),
        ):
            existing = screen_site.get(field)
            if existing not in {None, ""} and str(existing) != str(value):
                raise RuntimeError(f"Stable site/posthoc {field} conflict for {pair['sample']}:{focal}")
            screen_site[field] = value
        index = region_indexes.get(pair["sample"], {})
        focal_region = index.get((focal.chrom, focal.pos))
        partner_region = index.get((partner.chrom, partner.pos))
        shared_region = (
            focal_region
            if focal_region is not None
            and partner_region is not None
            and focal_region.get("region") == partner_region.get("region")
            else None
        )
        scope = topology_direct_scope(
            pair["focal_ssnv_branch"],
            pair["partner_ssnv_branch"],
            focal.pos,
            partner.pos,
            shared_region,
        )
        pair["topology_scope"] = scope
        pair["topology_region"] = shared_region.get("region") if shared_region else None
        pair["topology_order_status"] = (
            infer_topology_order(shared_region, focal.pos, partner.pos)
            if scope == "DIRECT_RETAINED_SHARED_REGION"
            else "NOT_INFERABLE_INDIRECT_SCOPE"
        )
    for site in site_rows:
        entry = entries_by_sample[site["sample"]]
        focal = LIB.Variant(site["chrom"], int(site["pos"]), site["ref"], site["alt"])
        ledger = fetch_ledger_row(entry, focal)
        observed = {
            "truth_label": truth_label(entry, focal),
            "ssnv_branch": ledger["ssnv_branch"],
            "component_id": ledger.get("component_id"),
            "selected_group_id": ledger.get("selected_group_id"),
        }
        for field, value in observed.items():
            if site.get(field) not in {None, ""} and str(site[field]) != str(value):
                raise RuntimeError(f"Stable site/posthoc {field} conflict for {site['sample']}:{focal}")
            site[field] = value
        key = site_key(site["sample"], site["chrom"], site["pos"])
        if source_site_rows is not None:
            source = source_site_rows[key]
            authoritative = {
                "biological_id": entry["biological_id"],
                "ref": focal.ref,
                "alt": focal.alt,
                **observed,
            }
            for field, value in authoritative.items():
                source_value = source.get(field)
                if source_value not in {None, ""} and str(source_value) != str(value):
                    raise RuntimeError(
                        f"Source stable site/authoritative posthoc {field} conflict: {key}"
                    )
        if source_assignments is not None:
            posthoc = source_assignments[key].get("posthoc") or {}
            authoritative = {
                "biological_id": entry["biological_id"],
                "ref": focal.ref,
                "alt": focal.alt,
                **observed,
            }
            for field, value in authoritative.items():
                assignment_value = posthoc.get(field)
                if assignment_value not in {None, ""} and str(assignment_value) != str(value):
                    raise RuntimeError(
                        f"Source stable assignment/authoritative posthoc {field} conflict: {key}"
                    )


def apply_dedup_statistics(site_rows: list[dict[str, Any]]) -> dict[str, Any]:
    key_fields = {
        "biological_site": lambda row: (
            row["biological_id"], row["chrom"], int(row["pos"]), row["ref"], row["alt"]
        ),
        "component": lambda row: (
            row["biological_id"],
            row.get("component_id") or f"SITE:{row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}",
        ),
        "alt_readset": lambda row: (
            row["biological_id"],
            row.get("alt_readset_sha256")
            or f"MISSING:{row['sample']}:{row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}",
        ),
    }
    representatives: dict[str, set[Any]] = {name: set() for name in key_fields}
    for row in sorted(site_rows, key=lambda value: (value["sample"], value["chrom"], int(value["pos"]))):
        for name, key_function in key_fields.items():
            key = key_function(row)
            row[f"{name}_key" if name == "biological_site" else f"{name}_dedup_key"] = json_text(key)
            representative_field = f"{name}_representative"
            row[representative_field] = key not in representatives[name]
            representatives[name].add(key)
    return {
        "n_dataset_sites": len(site_rows),
        "n_unique_biological_sites": len(representatives["biological_site"]),
        "n_unique_components": len(representatives["component"]),
        "n_unique_alt_readsets": len(representatives["alt_readset"]),
    }


def apply_cross_platform_replication(pair_rows: list[dict[str, Any]]) -> dict[str, Any]:
    def exact_key(row: dict[str, Any]) -> tuple[Any, ...]:
        return (
            row["focal_chrom"],
            int(row["focal_pos"]),
            row["focal_ref"],
            row["focal_alt"],
            row["partner_chrom"],
            int(row["partner_pos"]),
            row["partner_ref"],
            row["partner_alt"],
        )

    by_sample: dict[str, dict[tuple[Any, ...], dict[str, Any]]] = {
        "HCC1395": {},
        "HCC1395_DORADO": {},
    }
    for row in pair_rows:
        if row["sample"] in by_sample:
            key = exact_key(row)
            if key in by_sample[row["sample"]]:
                raise RuntimeError(f"Duplicate focal-partner pair for cross-platform join: {row['sample']} {key}")
            by_sample[row["sample"]][key] = row
    shared = set(by_sample["HCC1395"]).intersection(by_sample["HCC1395_DORADO"])
    for sample, rows in by_sample.items():
        other = "HCC1395_DORADO" if sample == "HCC1395" else "HCC1395"
        for key, row in rows.items():
            row["cross_platform_biological_n"] = 1
            if key not in by_sample[other]:
                row["cross_platform_replication_status"] = "ONE_PLATFORM_ONLY"
                row["cross_platform_molecule_independence_status"] = "NOT_EVALUABLE_ONE_PLATFORM_ONLY"
                continue
            row["cross_platform_exact_pair_present"] = True
            peer = by_sample[other][key]
            row_confirmed = bool(row.get("endpoint_a_formal_pair_by_confirmed"))
            peer_confirmed = bool(peer.get("endpoint_a_formal_pair_by_confirmed"))
            both_conditional_by = row_confirmed and peer_confirmed
            direction = row.get("endpoint_a_effect_direction")
            peer_direction = peer.get("endpoint_a_effect_direction")
            direction_compatible = bool(
                direction
                and peer_direction
                and direction == peer_direction
                and not str(direction).startswith("NOT_IDENTIFIABLE")
                and direction != "NO_DIRECTIONAL_CONTRAST"
            )
            row_v = row.get("endpoint_a_cramers_v")
            peer_v = peer.get("endpoint_a_cramers_v")
            row_delta = row.get("endpoint_a_delta_alt_fraction")
            peer_delta = peer.get("endpoint_a_delta_alt_fraction")
            effect_compatible = bool(
                row_v is not None
                and peer_v is not None
                and row_delta is not None
                and peer_delta is not None
                and float(row_v) >= ENDPOINT_A_V_MIN
                and float(peer_v) >= ENDPOINT_A_V_MIN
                and float(row_delta) >= ENDPOINT_A_DELTA_MIN
                and float(peer_delta) >= ENDPOINT_A_DELTA_MIN
            )
            relation = row.get("endpoint_b_relation_compatibility")
            peer_relation = peer.get("endpoint_b_relation_compatibility")
            four_state_compatible = bool(
                relation == peer_relation
                and relation in CROSS_PLATFORM_FOUR_STATE_COMPATIBLE_RELATIONS
            )
            row_names = set((row.get("_permutation_payload") or {}).get("core_raw_read_names") or [])
            peer_names = set((peer.get("_permutation_payload") or {}).get("core_raw_read_names") or [])
            intersection = row_names.intersection(peer_names)
            union = row_names.union(peer_names)
            overlap_present = bool(intersection) if row_names and peer_names else None
            molecule_status = (
                "SAME_RAW_MOLECULES_OVERLAP_NOT_INDEPENDENT_TECHNICAL_LIBRARY_REPLICATION"
                if overlap_present
                else "DISTINCT_CORE_READ_NAME_SETS_TECHNICAL_ONLY_BIOLOGICAL_N1"
                if overlap_present is False
                else "NOT_IDENTIFIABLE_CORE_READ_NAMES_UNAVAILABLE"
            )
            formal_pass = both_conditional_by and direction_compatible and effect_compatible
            for target in (row, peer):
                target["cross_platform_exact_pair_present"] = True
                target["cross_platform_both_conditional_by"] = both_conditional_by
                target["cross_platform_direction_compatible"] = direction_compatible
                target["cross_platform_effect_compatible"] = effect_compatible
                target["cross_platform_four_state_compatible"] = four_state_compatible
                target["cross_platform_cramers_v_abs_difference"] = (
                    abs(float(row_v) - float(peer_v)) if row_v is not None and peer_v is not None else None
                )
                target["cross_platform_delta_alt_fraction_abs_difference"] = (
                    abs(float(row_delta) - float(peer_delta))
                    if row_delta is not None and peer_delta is not None
                    else None
                )
                target["cross_platform_core_read_name_intersection_n"] = len(intersection)
                target["cross_platform_core_read_name_union_n"] = len(union)
                target["cross_platform_core_read_name_jaccard"] = (
                    len(intersection) / len(union) if union else None
                )
                target["cross_platform_core_read_name_overlap_present"] = overlap_present
                target["cross_platform_molecule_independence_status"] = molecule_status
                target["cross_platform_biological_n"] = 1
                target["cross_platform_conditional_by_effect_direction_replication_pass"] = formal_pass
            if not both_conditional_by:
                row["cross_platform_replication_status"] = "EXACT_PAIR_NOT_CONDITIONAL_BY_BOTH"
            elif not direction_compatible:
                row["cross_platform_replication_status"] = "EXACT_PAIR_DIRECTION_NOT_COMPATIBLE"
            elif not effect_compatible:
                row["cross_platform_replication_status"] = "EXACT_PAIR_EFFECT_NOT_COMPATIBLE"
            elif overlap_present:
                row["cross_platform_replication_status"] = (
                    "SAME_MOLECULE_REPROCESSING_CONDITIONAL_BY_CONCORDANCE"
                )
            else:
                row["cross_platform_replication_status"] = (
                    "DISTINCT_MOLECULE_SET_TECHNICAL_CONDITIONAL_BY_CONCORDANCE"
                )
    replicated = sum(
        by_sample["HCC1395"][key].get(
            "cross_platform_conditional_by_effect_direction_replication_pass", False
        )
        for key in shared
    )
    overlap_pairs = sum(
        by_sample["HCC1395"][key].get("cross_platform_core_read_name_overlap_present") is True
        for key in shared
    )
    both_conditional_by_pairs = sum(
        by_sample["HCC1395"][key].get("cross_platform_both_conditional_by", False)
        for key in shared
    )
    for row in pair_rows:
        row.pop("_permutation_payload", None)
    return {
        "contract": "exact_focal_and_partner_chrom_pos_ref_alt",
        "inference_contract": (
            "both platforms fixed-margins exact global-BY plus HP-family_x_PS_x_strand "
            "conditional sensitivity and compatible effect direction/floors"
        ),
        "replication_scope": "technical_only_same_biological_sample",
        "biological_n": 1,
        "independent_technical_library_claim_allowed": False,
        "n_hcc1395_pairs": len(by_sample["HCC1395"]),
        "n_hcc1395_dorado_pairs": len(by_sample["HCC1395_DORADO"]),
        "n_exact_pairs_present_both": len(shared),
        "n_exact_pairs_conditional_by_effect_direction_replicated": replicated,
        "n_exact_pairs_with_core_raw_read_name_overlap": overlap_pairs,
        "n_exact_pairs_confirmed_both": both_conditional_by_pairs,
        "n_exact_pairs_confirmed_both_semantics": "deprecated_alias_for_both_conditional_BY",
    }


def extract_oracle_cases(
    site_rows: Sequence[dict[str, Any]],
    pair_rows: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    site_index: dict[tuple[str, str, int], list[dict[str, Any]]] = defaultdict(list)
    pair_index: dict[tuple[str, str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in site_rows:
        site_index[site_key(row["sample"], row["chrom"], row["pos"])].append(row)
    for row in pair_rows:
        pair_index[site_key(row["sample"], row["focal_chrom"], row["focal_pos"])].append(row)
    cases = []
    for oracle in ORACLE_FOCALS:
        key = site_key(oracle["sample"], oracle["chrom"], oracle["pos"])
        sites = site_index.get(key, [])
        pairs = sorted(pair_index.get(key, []), key=lambda row: int(row["partner_pos"]))
        if len(sites) > 1:
            raise RuntimeError(f"Oracle site is not unique: {key}")
        observed_positions = [int(row["partner_pos"]) for row in pairs]
        cases.append(
            {
                "case_id": oracle["case_id"],
                "case_type": "focal_partner_window",
                "sample": oracle["sample"],
                "focal": {"chrom": oracle["chrom"], "pos": oracle["pos"]},
                "present": len(sites) == 1,
                "site": sites[0] if sites else None,
                "pairs": pairs,
                "expected_partner_positions": list(oracle["expected_partner_positions"]),
                "observed_partner_positions": observed_positions,
                "partner_window_oracle_pass": observed_positions
                == list(oracle["expected_partner_positions"]),
            }
        )
    shared = SHARED_READSET_ORACLE
    shared_sites = []
    for position in shared["positions"]:
        matches = site_index.get(site_key(shared["sample"], shared["chrom"], position), [])
        if len(matches) > 1:
            raise RuntimeError(f"Shared-readset oracle site is not unique: {position}")
        if matches:
            shared_sites.append(matches[0])
    digests = [row.get("alt_readset_sha256") for row in shared_sites]
    shared_case = {
        "case_id": shared["case_id"],
        "case_type": "shared_alt_readset_dedup",
        "sample": shared["sample"],
        "chrom": shared["chrom"],
        "positions": list(shared["positions"]),
        "present_positions": [int(row["pos"]) for row in shared_sites],
        "sites": shared_sites,
        "same_alt_readset": len(shared_sites) == 2 and len(set(digests)) == 1 and bool(digests[0]),
        "one_alt_readset_representative": len(shared_sites) == 2
        and sum(bool(row.get("alt_readset_representative")) for row in shared_sites) == 1,
    }
    return {
        "schema_name": f"{SCHEMA_NAME}.oracle_cases",
        "schema_version": SCHEMA_VERSION,
        "partner_window_bp": PAIR_WINDOW_BP,
        "focal_cases": cases,
        "shared_readset_case": shared_case,
    }


def validate_oracle_cases(payload: dict[str, Any], full_scope: bool) -> None:
    required = {"schema_name", "schema_version", "partner_window_bp", "focal_cases", "shared_readset_case"}
    missing = sorted(required.difference(payload))
    if missing:
        raise RuntimeError(f"Oracle extraction schema missing fields: {missing}")
    if not full_scope:
        return
    failed = [
        case["case_id"]
        for case in payload["focal_cases"]
        if not case["present"] or not case["partner_window_oracle_pass"]
    ]
    shared = payload["shared_readset_case"]
    if not shared["same_alt_readset"] or not shared["one_alt_readset_representative"]:
        failed.append(shared["case_id"])
    if failed:
        raise RuntimeError(f"Full-scope oracle validation failed: {failed}")


def _tsv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (dict, list, tuple)):
        return json_text(value)
    if isinstance(value, float) and not math.isfinite(value):
        raise RuntimeError(f"Non-finite value cannot be emitted: {value}")
    return value


def write_tsv_gz(path: Path, rows: Sequence[dict[str, Any]], columns: Sequence[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            if set(row) != set(columns):
                missing = sorted(set(columns).difference(row))
                extra = sorted(set(row).difference(columns))
                raise RuntimeError(f"TSV schema mismatch for {path.name}: missing={missing} extra={extra}")
            writer.writerow({column: _tsv_value(row[column]) for column in columns})


def create_output_dir(path: Path) -> None:
    if os.path.lexists(path):
        raise FileExistsError(f"Refusing to overwrite existing output path: {path}")
    if not path.parent.is_dir():
        raise FileNotFoundError(f"Output parent directory must already exist: {path.parent}")
    path.mkdir()


def build_summary(
    entries: Sequence[dict[str, Any]],
    site_rows: Sequence[dict[str, Any]],
    pair_rows: Sequence[dict[str, Any]],
    dedup: dict[str, Any],
    replication: dict[str, Any],
    geometry: dict[str, Any],
    stable_partner_audit: dict[str, Any],
    conditional_sensitivity_audit: dict[str, Any],
    joint_signature_global_fdr_audit: dict[str, Any],
    site_pair_reconciliation: dict[str, Any],
    full_scope: bool,
) -> dict[str, Any]:
    pairs_per_focal = Counter(
        site_key(row["sample"], row["focal_chrom"], row["focal_pos"]) for row in pair_rows
    )
    candidates = [
        row for row in site_rows if row["multi_marker_molecular_haplotype_base_candidate"]
    ]
    identity_rows = [
        {"sample": entry["sample"], **row}
        for entry in entries
        for row in entry["_frozen_artifact_identity_audit"]
    ]
    raw_identity_site_digest = hashlib.sha256()
    for row in sorted(
        site_rows,
        key=lambda value: (value["sample"], value["chrom"], int(value["pos"])),
    ):
        raw_identity_site_digest.update(
            json_text(raw_identity_site_digest_payload(row)).encode("utf-8")
        )
        raw_identity_site_digest.update(b"\n")
    raw_identity_summary = {
        "equivalence_contract": RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
        "allowed_differing_auxiliary_tags": sorted(
            RAW_DUPLICATE_ALLOWED_DIFFERING_TAGS
        ),
        "deterministic_representative": "first_BAM_fetch_record_after_equivalence_proof",
        "n_site_tasks": len(site_rows),
        "n_expected_projection_occurrences": sum(
            int(row["raw_identity_expected_projections"]) for row in site_rows
        ),
        "n_matched_record_occurrences": sum(
            int(row["raw_identity_matched_records"]) for row in site_rows
        ),
        "n_sites_with_collapsed_duplicates": sum(
            int(row["raw_identity_duplicate_projections_collapsed"]) > 0
            for row in site_rows
        ),
        "n_duplicate_projection_occurrences_collapsed": sum(
            int(row["raw_identity_duplicate_projections_collapsed"])
            for row in site_rows
        ),
        "n_sparse_duplicate_rows": sum(
            int(row["raw_identity_duplicate_projections_collapsed"])
            for row in site_rows
        ),
        "n_duplicate_extra_record_occurrences_collapsed": sum(
            int(row["raw_identity_duplicate_extra_records_collapsed"])
            for row in site_rows
        ),
        "n_exact_duplicate_projection_occurrences_collapsed": sum(
            int(row["raw_identity_exact_duplicate_projections_collapsed"])
            for row in site_rows
        ),
        "n_rg_only_duplicate_projection_occurrences_collapsed": sum(
            int(row["raw_identity_rg_only_duplicate_projections_collapsed"])
            for row in site_rows
        ),
        "site_weighted_audit_sha256": raw_identity_site_digest.hexdigest(),
        "all_site_results_passed_invariant_validation": True,
        "missing_projection_policy": RAW_IDENTITY_MISSING_POLICY,
        "conflicting_analysis_payload_policy": RAW_IDENTITY_CONFLICT_POLICY,
        "failure_counts_materialized": False,
        "count_semantics": (
            "site-weighted projection occurrences; one long read can occur at multiple focal sites"
        ),
    }

    def candidate_dedup_counts(rows: Sequence[dict[str, Any]]) -> dict[str, int]:
        return {
            "n_dataset_sites": len(rows),
            "n_unique_biological_sites": len({row["biological_site_key"] for row in rows}),
            "n_unique_components": len({row["component_dedup_key"] for row in rows}),
            "n_unique_alt_readsets": len({row["alt_readset_dedup_key"] for row in rows}),
        }
    return {
        "schema_name": f"{SCHEMA_NAME}.summary",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": now_utc(),
        "status": "EXECUTION_PASS",
        "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
        "task_type": "B_comprehensive_validation" if full_scope else "A_partial_explicit_scope",
        "scope": "all_manifest_samples" if full_scope else "explicit_sample_subset",
        "samples": [entry["sample"] for entry in entries],
        "selection_contract": {
            "focal_sites": "stable_null_multigroup_only",
            "partners": "all frozen latest LongPhase-S PASS biallelic sSNVs within inclusive +/-5kb",
            "analysis_pair_rows": (
                "directed stable-focal-to-partner rows; reciprocal rows are retained when both sites are stable"
            ),
            "region_dir_binding": "inside validated v2 InterSubMod sample receipt output_dir",
            "selected_group_restriction": False,
            "truth_ledger_topology_role": "posthoc_annotation_only",
        },
        "frozen_input_identity_policy": {
            "hash_required_fields": sorted(MANIFEST_HASH_REQUIRED_FIELDS),
            "explicit_large_metadata_identity_fields": sorted(
                MANIFEST_LARGE_METADATA_IDENTITY_FIELDS
            ),
            "n_sha256_identity_artifacts": sum(
                row["identity_mode"] == "sha256_size_path_v1" for row in identity_rows
            ),
            "n_explicit_large_size_mtime_identity_artifacts": sum(
                row["identity_mode"]
                == "explicit_large_artifact_size_mtime_path_v1"
                for row in identity_rows
            ),
            "artifacts": identity_rows,
            "claim_guardrail": (
                "Only the explicitly listed large BAM/index/tag-sidecar roles may use "
                "size+mtime+path identity; this is weaker than full SHA-256 and is not "
                "reported as bytewise identity."
            ),
        },
        "m2_axis_statistic_provenance": {
            "effect_and_permutation_p_source": (
                "frozen source-locked screen-producer site rows"
            ),
            "downstream_raw_read_axis_statistic_recomputed": False,
            "downstream_axis_classification_recomputed": True,
            "downstream_validation": [
                "axis sample-size reconciliation",
                "499-permutation add-one p-value grid",
                "effect threshold classification",
                "80-percent planning-power evaluability",
            ],
            "claim_guardrail": (
                "The cooccurrence stage independently reclassifies producer-derived axis "
                "statistics; it is not an independent raw-read remeasurement."
            ),
        },
        "read_contract": {
            "raw_bam_identity_recovery": (
                "one record, or multiple records proven equivalent under the declared contract"
            ),
            "raw_bam_duplicate_equivalence_contract": RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
            "raw_bam_duplicate_allowed_differing_auxiliary_tags": sorted(
                RAW_DUPLICATE_ALLOWED_DIFFERING_TAGS
            ),
            "raw_bam_duplicate_conflict_policy": "hard_fail",
            "excluded_flags": ["unmapped", "secondary", "supplementary", "duplicate"],
            "minimum_mapq": MIN_MAPQ,
            "minimum_cigar_query_length": MIN_QUERY_LENGTH,
            "require_MM_ML": True,
            "minimum_partner_base_quality": MIN_BASE_QUALITY,
            "partner_calls": ["R", "A", "O", "X"],
            "focal_reads": "all tumor REF+ALT reads.tsv rows",
            "focal_call_conflict_policy": "hard_fail",
        },
        "raw_bam_identity_recovery_audit": raw_identity_summary,
        "endpoint_a_contract": {
            "population": "core focal-ALT stable methyl-group reads",
            "descriptive_only": "analytic Fisher 2x2 or Pearson chi-square p; never a gate",
            "primary_test": (
                "Fisher-Freeman-Halton Kx2 fixed-margins probability-ordered two-sided exact p"
            ),
            "exact_state_space_ceiling": sorted(
                {row["endpoint_a_exact_state_space_ceiling"] for row in pair_rows}
            ),
            "state_space_fallback": (
                "fail closed as NOT_IDENTIFIABLE_STATE_SPACE_LIMIT; never substitute chi-square"
            ),
            "global_fdr_family": (
                "pairs with exact-identifiable Kx2 tables from pre-screened M2 residual-unexplained "
                "sites only; stable/confounded pairs retain descriptive exact p but q is null"
            ),
            "formal_pair_gate": {
                "q_max": ENDPOINT_A_Q_MAX,
                "cramers_v_min": ENDPOINT_A_V_MIN,
                "delta_alt_fraction_min": ENDPOINT_A_DELTA_MIN,
            },
            "effect_direction": (
                "partner ALT fraction in the unique dominant methyl group minus the pooled other "
                "groups; tied dominant group is NOT_IDENTIFIABLE"
            ),
            "conditional_permutation": (
                "secondary sensitivity only for exact-BY discoveries; methyl labels within "
                "latest HP-family x PS x strand strata"
            ),
            "conditional_p_max_inclusive": CONDITIONAL_P_MAX,
            "conditional_permutation_global_fdr_calibrated": False,
        },
        "multiplicity_dependence_contract": {
            "primary_descriptive": "global BH on fixed-margins exact p within the locked M2 family",
            "formal": "global BY on the same fixed-margins exact-p family",
            "dependence_sources": [
                "multiple partners share one focal read set and methyl labels",
                "dense windows share partner markers",
                "biological/component/alt-readset duplicate focal units",
            ],
            "claim_guardrail": (
                "Spatially separated pair hits are not statistically independent or independent "
                "biological replicates."
            ),
        },
        "global_exact_fdr_family": {
            "n_m2_eligible_sites": sum(row["m2_screen_eligible"] for row in site_rows),
            "n_m2_evaluable_sites": sum(row["m2_screen_evaluable"] for row in site_rows),
            "n_m2_axis_indeterminate_sites": sum(
                bool(row["m2_indeterminate_axes"]) for row in site_rows
            ),
            "n_m2_axis_low_power_sites": sum(
                bool(row["m2_low_power_axes"]) for row in site_rows
            ),
            "n_m2_aligned_below_negative_evaluability_power_sites": sum(
                bool(row["m2_aligned_below_negative_evaluability_power_axes"])
                for row in site_rows
            ),
            "m2_screen_gate_contract": M2_GATE.GATE_CONTRACT,
            "m2_axis_status_counts": {
                axis: dict(
                    sorted(
                        Counter(
                            row["m2_axis_statuses"][axis]["status"] for row in site_rows
                        ).items()
                    )
                )
                for axis, _kind, _effect_suffix, _threshold in M2_GATE.AXIS_SPECS
            },
            "n_exact_family_pairs": sum(
                row["endpoint_a_global_fdr_family_status"] == "ELIGIBLE_M2_EXACT_FAMILY"
                for row in pair_rows
            ),
            "n_ineligible_m2_screen_pairs": sum(
                row["endpoint_a_global_fdr_family_status"] == "INELIGIBLE_M2_SCREEN"
                for row in pair_rows
            ),
            "state_space_status_counts": dict(
                sorted(Counter(row["endpoint_a_exact_state_space_status"] for row in pair_rows).items())
            ),
        },
        "conditional_sensitivity_audit": conditional_sensitivity_audit,
        "joint_signature_global_fdr_audit": joint_signature_global_fdr_audit,
        "joint_signature_contract": {
            "marker_selection": "coverage descending with deterministic positional ties; association-blind",
            "molecules": "one joint category per complete read spanning every selected top marker",
            "per_marker_same_complete_read_gate": (
                "at least two spatially separated exact-BY/effect/conditional pairs must each remain "
                "R/A-testable with Cramers_V>=0.30 and delta_ALT_fraction>=0.50 in the same "
                "complete-read set"
            ),
            "test": "fixed 999 label permutations within HP-family x PS x strand",
            "role": (
                "global-BY formal candidate gate over all M2-eligible testable/permutable sites; "
                "raw p<=0.05 is retained as sensitivity only"
            ),
            "global_fdr_calibrated": True,
            "formal_q_method": "Benjamini-Yekutieli",
            "formal_q_max": ENDPOINT_A_Q_MAX,
            "postselection_fdr_calibrated": False,
            "analytic_chi_square_used": False,
        },
        "partner_geometry_audit": geometry,
        "stable_focal_partner_enumeration_audit": stable_partner_audit,
        "site_pair_count_reconciliation": site_pair_reconciliation,
        "endpoint_b_guardrail": (
            "Mutation-order compatibility under a fixed error model is not proof of cellular ancestry."
        ),
        "top_marker_contract": (
            "per focal: core focal-ALT R/A callable coverage descending, then absolute distance, "
            "position, ref, alt; never p-value"
        ),
        "n_stable_sites": len(site_rows),
        "n_focal_partner_pairs": len(pair_rows),
        "n_endpoint_a_testable_pairs": sum(row["endpoint_a_testable"] for row in pair_rows),
        "n_endpoint_a_exact_bh_discoveries": sum(
            row["endpoint_a_exact_bh_discovery"] for row in pair_rows
        ),
        "n_endpoint_a_exact_by_discoveries": sum(
            row["endpoint_a_exact_by_discovery"] for row in pair_rows
        ),
        "n_pair_by_confirmed": sum(
            row["endpoint_a_formal_pair_by_confirmed"] for row in pair_rows
        ),
        "n_pair_callability_gate_pass": sum(row["callability_gate_pass"] for row in pair_rows),
        "callability_gate_status_counts": dict(
            sorted(Counter(row["callability_gate_status"] for row in pair_rows).items())
        ),
        "n_multi_marker_molecular_haplotype_base_candidates": len(candidates),
        "deprecated_legacy_alias_counts": {
            "n_endpoint_a_pre_candidates": sum(
                row["endpoint_a_pre_candidate"] for row in pair_rows
            ),
            "n_genetically_anchored_multi_marker_candidates": sum(
                row["genetically_anchored_multi_marker_candidate"] for row in site_rows
            ),
            "guardrail": "Legacy names are aliases only and must not be used as formal claims.",
        },
        "dedup_denominators": dedup,
        "candidate_dedup_denominators": candidate_dedup_counts(candidates),
        "stable_focal_pair_density": {
            "n_focals_with_at_least_one_partner": len(pairs_per_focal),
            "median_partners_per_nonempty_focal": (
                float(np.median(list(pairs_per_focal.values()))) if pairs_per_focal else 0.0
            ),
            "max_partners_per_focal": max(pairs_per_focal.values(), default=0),
        },
        "cross_platform_replication": replication,
        "pass": True,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--primary-artifact-audit-pre", type=Path, required=True)
    parser.add_argument("--intersubmod-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--allow-partial-scope", action="store_true")
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--chunk-size", type=int, default=8)
    parser.add_argument("--max-pending-chunks", type=int, default=16)
    parser.add_argument("--permutations", type=int, default=MIN_PERMUTATIONS)
    parser.add_argument("--top-markers", type=int, default=DEFAULT_TOP_MARKERS)
    parser.add_argument(
        "--exact-state-space-ceiling", type=int, default=EXACT_STATE_SPACE_CEILING
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    started_at_utc = now_utc()
    if args.permutations != MIN_PERMUTATIONS:
        raise ValueError(f"--permutations must equal canonical value {MIN_PERMUTATIONS}")
    if args.workers < 1 or args.chunk_size < 1 or args.max_pending_chunks < args.workers:
        raise ValueError("Invalid bounded-process settings")
    if args.top_markers != DEFAULT_TOP_MARKERS:
        raise ValueError(f"--top-markers must equal canonical value {DEFAULT_TOP_MARKERS}")
    if args.exact_state_space_ceiling != EXACT_STATE_SPACE_CEILING:
        raise ValueError(
            "--exact-state-space-ceiling must equal canonical value "
            f"{EXACT_STATE_SPACE_CEILING}"
        )
    for path in (
        args.manifest,
        args.assignments,
        args.sites,
        args.primary_artifact_audit_pre,
    ):
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(path)
    if os.path.lexists(args.output_dir):
        raise FileExistsError(f"Refusing to overwrite existing output path: {args.output_dir}")
    source_identity_before = capture_code_source_identity()
    source_modes_before = capture_code_source_modes()
    require_read_only_source_modes(source_modes_before, "start")

    requested = set(args.sample) or None
    manifest, entries, run_validation = load_manifest(
        args.manifest,
        args.intersubmod_root,
        requested,
        args.allow_partial_scope,
    )
    selected_samples = {entry["sample"] for entry in entries}
    full_scope = selected_samples == {entry["sample"] for entry in manifest["samples"]}
    geometry = audit_partner_geometry(entries, full_scope=full_scope)
    assignments = load_assignments(args.assignments, selected_samples)
    stable_sites = load_stable_site_rows(
        args.sites,
        assignments,
        selected_samples,
        require_candidate_screen_fields=True,
    )
    stable_partner_oracle = build_stable_partner_oracle(entries, assignments)
    tasks = build_tasks(
        entries,
        assignments,
        stable_sites,
        args.top_markers,
        exact_state_space_ceiling=args.exact_state_space_ceiling,
    )

    site_rows: list[dict[str, Any]] = []
    pair_rows: list[dict[str, Any]] = []
    raw_identity_duplicate_rows: list[dict[str, Any]] = []
    for result in bounded_process_map(
        tasks,
        workers=args.workers,
        chunk_size=args.chunk_size,
        max_pending_chunks=args.max_pending_chunks,
    ):
        site_rows.append(result["site"])
        pair_rows.extend(result["pairs"])
        raw_identity_duplicate_rows.extend(result["raw_identity_duplicates"])
    site_rows.sort(key=lambda row: (row["sample"], row["chrom"], int(row["pos"])))
    pair_rows.sort(
        key=lambda row: (
            row["sample"],
            row["focal_chrom"],
            int(row["focal_pos"]),
            int(row["partner_pos"]),
            row["partner_ref"],
            row["partner_alt"],
        )
    )
    raw_identity_duplicate_rows.sort(
        key=lambda row: (
            row["sample"],
            row["chrom"],
            int(row["pos"]),
            row["projection_read_name"],
            row["projection_chrom"],
            int(row["projection_start"]),
            int(row["projection_end"]),
            int(row["projection_mapq"]),
            row["projection_strand"],
        )
    )
    if len(raw_identity_duplicate_rows) != sum(
        int(row["raw_identity_duplicate_projections_collapsed"])
        for row in site_rows
    ):
        raise RuntimeError("Raw identity sparse duplicate output count does not reconcile")
    stable_partner_audit = validate_stable_partner_rows(stable_partner_oracle, pair_rows)

    apply_endpoint_a_bh_and_gate(pair_rows)
    conditional_sensitivity_audit = run_conditional_permutations(
        pair_rows,
        args.permutations,
        workers=args.workers,
        chunk_size=max(1, args.chunk_size * 2),
    )
    joint_signature_global_fdr_audit = apply_joint_signature_global_fdr(site_rows)
    site_pair_reconciliation = finalize_site_inference(site_rows, pair_rows)
    attach_posthoc_annotations(
        entries,
        site_rows,
        pair_rows,
        source_site_rows=stable_sites,
        source_assignments=assignments,
    )
    dedup = apply_dedup_statistics(site_rows)
    replication = apply_cross_platform_replication(pair_rows)
    oracle_cases = extract_oracle_cases(site_rows, pair_rows)
    validate_oracle_cases(oracle_cases, full_scope=full_scope)
    summary = build_summary(
        entries,
        site_rows,
        pair_rows,
        dedup,
        replication,
        geometry,
        stable_partner_audit,
        conditional_sensitivity_audit,
        joint_signature_global_fdr_audit,
        site_pair_reconciliation,
        full_scope,
    )

    source_identity_after_compute = capture_code_source_identity()
    source_modes_after_compute = capture_code_source_modes()
    require_read_only_source_modes(source_modes_after_compute, "after_compute")
    if source_identity_before != source_identity_after_compute:
        raise RuntimeError("Code source identity changed during analysis")

    create_output_dir(args.output_dir)
    pair_path = args.output_dir / PAIR_OUTPUT_NAME
    site_path = args.output_dir / SITE_OUTPUT_NAME
    case_path = args.output_dir / CASE_OUTPUT_NAME
    summary_path = args.output_dir / SUMMARY_OUTPUT_NAME
    receipt_path = args.output_dir / RECEIPT_OUTPUT_NAME
    raw_identity_duplicate_path = args.output_dir / RAW_IDENTITY_DUPLICATE_OUTPUT_NAME
    write_tsv_gz(pair_path, pair_rows, PAIR_COLUMNS)
    write_tsv_gz(site_path, site_rows, SITE_COLUMNS)
    write_tsv_gz(
        raw_identity_duplicate_path,
        raw_identity_duplicate_rows,
        RAW_IDENTITY_DUPLICATE_COLUMNS,
    )
    case_path.write_text(json.dumps(oracle_cases, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    source_identity_after_output = capture_code_source_identity()
    source_modes_after_output = capture_code_source_modes()
    require_read_only_source_modes(source_modes_after_output, "after_output")
    if source_identity_before != source_identity_after_output:
        raise RuntimeError("Code source identity changed while outputs were written")
    code_artifacts = source_identity_after_output
    input_artifacts = {
        "manifest": artifact(args.manifest),
        "assignments": artifact(args.assignments),
        "sites": artifact(args.sites),
        "primary_artifact_audit_pre": artifact(args.primary_artifact_audit_pre),
        "intersubmod_runs": run_validation,
    }
    output_artifacts = {
        "pair_tsv": artifact(pair_path),
        "site_tsv": artifact(site_path),
        "case_json": artifact(case_path),
        "summary_json": artifact(summary_path),
        "raw_identity_duplicate_audit_tsv": artifact(raw_identity_duplicate_path),
    }
    finished_at_utc = now_utc()
    if datetime.fromisoformat(finished_at_utc) < datetime.fromisoformat(started_at_utc):
        raise RuntimeError("Run completion timestamp precedes start timestamp")
    receipt = {
        "schema_name": f"{SCHEMA_NAME}.run_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": finished_at_utc,
        "started_at_utc": started_at_utc,
        "finished_at_utc": finished_at_utc,
        "pass": True,
        "command": [sys.executable, str(Path(__file__).resolve()), *(argv if argv is not None else sys.argv[1:])],
        "code": code_artifacts,
        "source_lock": {
            "source_identity_before": source_identity_before,
            "source_identity_after_compute": source_identity_after_compute,
            "source_identity_after_output": source_identity_after_output,
            "source_modes_before": source_modes_before,
            "source_modes_after_compute": source_modes_after_compute,
            "source_modes_after_output": source_modes_after_output,
            "all_sources_read_only_and_unchanged": True,
        },
        "inputs": input_artifacts,
        "parameters": {
            "workers": args.workers,
            "chunk_size": args.chunk_size,
            "max_pending_chunks": args.max_pending_chunks,
            "permutations": args.permutations,
            "conditional_permutation_role": "secondary_sensitivity_not_global_FDR",
            "joint_signature_permutations": MIN_PERMUTATIONS,
            "top_markers": args.top_markers,
            "pair_window_bp": PAIR_WINDOW_BP,
            "exact_state_space_ceiling": args.exact_state_space_ceiling,
            "exact_state_space_fallback": "fail_closed_NOT_IDENTIFIABLE_STATE_SPACE_LIMIT",
            "canonical_formal_parameter_gate": True,
        },
        "frozen_manifest_input_identity_policy": summary[
            "frozen_input_identity_policy"
        ],
        "m2_axis_statistic_provenance": summary[
            "m2_axis_statistic_provenance"
        ],
        "counts": {
            "stable_sites": len(site_rows),
            "focal_partner_pairs": len(pair_rows),
            "m2_global_fdr_eligible_sites": site_pair_reconciliation[
                "n_m2_global_fdr_eligible_sites"
            ],
            "m2_exact_global_fdr_family_pairs": site_pair_reconciliation[
                "n_m2_exact_global_fdr_family_pairs"
            ],
            "pair_by_confirmed": sum(
                row["endpoint_a_formal_pair_by_confirmed"] for row in pair_rows
            ),
            "multi_marker_molecular_haplotype_base_candidates": sum(
                row["multi_marker_molecular_haplotype_base_candidate"] for row in site_rows
            ),
            "stable_partner_oracle_pairs": stable_partner_audit[
                "n_expected_directed_pair_rows"
            ],
            "all_ssnv_unordered_geometry_pairs": geometry["aggregate"]["n_unordered_pairs"],
            "raw_identity_expected_projection_occurrences": summary[
                "raw_bam_identity_recovery_audit"
            ]["n_expected_projection_occurrences"],
            "raw_identity_duplicate_projection_occurrences_collapsed": summary[
                "raw_bam_identity_recovery_audit"
            ]["n_duplicate_projection_occurrences_collapsed"],
            "raw_identity_duplicate_extra_record_occurrences_collapsed": summary[
                "raw_bam_identity_recovery_audit"
            ]["n_duplicate_extra_record_occurrences_collapsed"],
            "raw_identity_sparse_duplicate_rows": len(raw_identity_duplicate_rows),
            "raw_identity_all_site_results_passed_invariant_validation": summary[
                "raw_bam_identity_recovery_audit"
            ]["all_site_results_passed_invariant_validation"],
            "raw_identity_missing_projection_policy": RAW_IDENTITY_MISSING_POLICY,
            "raw_identity_conflicting_analysis_payload_policy": (
                RAW_IDENTITY_CONFLICT_POLICY
            ),
            "raw_identity_failure_counts_materialized": False,
        },
        "outputs": output_artifacts,
    }
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(args.output_dir.resolve()), "pass": True, **receipt["counts"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
