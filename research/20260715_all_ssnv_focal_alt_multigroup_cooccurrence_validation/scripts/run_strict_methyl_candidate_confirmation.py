#!/usr/bin/env python3
"""Strict methyl-partition robustness audit for preselected molecular-haplotype candidates.

The candidate selection column is an immutable upstream decision. This sidecar
never reruns or relaxes molecular-haplotype discovery, cooccurrence, effect-size,
or spatially-separated-marker gates. It subjects selected sites only to stricter
methylation null and fixed-K partition-stability checks using hash-locked primary
InterSubMod artifacts.

Because the methylation nulls are evaluated only after methyl-genetic candidate
selection, their p/q values are post-selection diagnostics, not calibrated FDR
confirmation. Scientific output is therefore a robustness status, not a claim
that a cellular subclone has been confirmed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib
import inspect
import io
import json
import math
import os
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence, TextIO

import numpy as np
import scipy
import sklearn
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import adjusted_rand_score


TOPIC_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = TOPIC_ROOT.parents[1]
SCRIPT_ROOT = TOPIC_ROOT / "scripts"
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))
FP_SCRIPT_ROOT = (
    TOPIC_ROOT.parent / "20260715_single_fp_alt_multicluster_subclone_validation" / "scripts"
)
if str(FP_SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(FP_SCRIPT_ROOT))

import focal_alt_cluster_lib as F  # noqa: E402
import release_source_authority as SOURCE_AUTHORITY  # noqa: E402


SCHEMA_VERSION = "3.1.0"
GENOME_WIDE_SITE_REFERENCE = 469_849
ASSIGNMENT_SCHEMA = "intersubmod.all_ssnv_stable_multigroup_read_assignment"
ASSIGNMENT_SCHEMA_VERSION = "1.0.0"
COOCCURRENCE_RECEIPT_SCHEMA = "intersubmod.methyl_ssnv_cooccurrence_validation.run_receipt"
COOCCURRENCE_SUMMARY_SCHEMA = "intersubmod.methyl_ssnv_cooccurrence_validation.summary"
COOCCURRENCE_FORMAL_SCHEMA_VERSION = "2.0.0"
COOCCURRENCE_RELEASE_SCHEMA = "intersubmod.cooccurrence_release_receipt"
COOCCURRENCE_RELEASE_SCHEMA_VERSION = "1.0.0"
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
COOCCURRENCE_CODE_PATHS = {
    "analyzer": TOPIC_ROOT / "scripts" / "analyze_methyl_ssnv_cooccurrence.py",
    "ssnv_cooccurrence_lib": REPO_ROOT / "scripts" / "ssnv_cooccurrence_lib.py",
    "latest_tag_join": TOPIC_ROOT / "scripts" / "latest_tag_join.py",
    "m2_screen_gate": TOPIC_ROOT / "scripts" / "m2_screen_gate.py",
    "release_finalizer": (
        TOPIC_ROOT / "scripts" / "finalize_cooccurrence_release_receipt.py"
    ),
    "release_runner": TOPIC_ROOT / "scripts" / "run_cooccurrence_v6_source_locked.sh",
    "source_authority_validator": (
        TOPIC_ROOT / "scripts" / "release_source_authority.py"
    ),
}
COOCCURRENCE_PREFLIGHT_SOURCE = (
    TOPIC_ROOT / "scripts" / "audit_cooccurrence_task_contract_preflight.py"
)
COOCCURRENCE_RELEASE_CODE_PATHS = {
    "release_finalizer": (
        TOPIC_ROOT / "scripts" / "finalize_cooccurrence_release_receipt.py"
    ),
    "release_runner": TOPIC_ROOT / "scripts" / "run_cooccurrence_v6_source_locked.sh",
}
SCREEN_CONTRACT = "phylo-v4.1_column_null95_modal_K10_RNULL40_min_group3"
ARTIFACT_IDENTITY_CONTRACT = "sha256_size_path_v1"
FORMAL_SELECTION_COLUMN = "multi_marker_molecular_haplotype_base_candidate"
LEGACY_SELECTION_COLUMNS = (
    "genetically_anchored_multi_marker_candidate_by_sensitivity",
)
DEFAULT_SELECTION_COLUMN = FORMAL_SELECTION_COLUMN
DEFAULT_SEED = 20260715
DEFAULT_PERMUTATIONS = 999
DEFAULT_SEEDS = 10
FORMAL_MIN_PERMUTATIONS = 999
FORMAL_MIN_SEEDS = 10
DEFAULT_ASSIGNMENT_ARI = 0.80
DEFAULT_MIN_VALID_NULL_FRACTION = 0.80
CANONICAL_PYTHON_CACHE_DIRNAME = ".python_cache_m2v5_completion_v2_bound_bootstrap"
NULL_MODES = ("column", "row_circular")
MODE_SEED_OFFSETS = {"column": 0, "row_circular": 1_000_000}
POSTSELECTION_SCOPE = (
    "upstream_multi_marker_molecular_haplotype_base_candidate_"
    "methyl_partition_diagnostic_only_not_fdr_calibrated"
)
GUARDRAIL = (
    "This strict sidecar evaluates methyl-partition robustness only. It does not rerun or confirm "
    "molecular-haplotype discovery, genetic cooccurrence, a cellular subclone, clone count, or lineage. "
    "Matched-normal evidence, copy-number/purity-aware cellular-fraction evidence, and orthogonal "
    "cellular truth remain separate requirements."
)
STRICT_CODE_PATHS = {
    "strict_producer": Path(__file__).resolve(),
    "focal_alt_cluster_lib": Path(F.__file__).resolve(),
    "source_authority_validator": Path(SOURCE_AUTHORITY.__file__).resolve(),
}


class ConfirmationContractError(RuntimeError):
    """Raised when an input or primary-artifact identity contract is violated."""


@dataclass(frozen=True)
class ConfirmationConfig:
    permutations: int = DEFAULT_PERMUTATIONS
    seeds: int = DEFAULT_SEEDS
    seed: int = DEFAULT_SEED
    assignment_ari_threshold: float = DEFAULT_ASSIGNMENT_ARI
    min_valid_null_fraction: float = DEFAULT_MIN_VALID_NULL_FRACTION


@dataclass(frozen=True)
class SelectionColumnContract:
    formal_column: str
    resolved_column: str
    legacy_fallback_used: bool


@dataclass
class PreparedSite:
    input_row: dict[str, str]
    source_row_number: int
    sample: str
    chrom: str
    pos: int
    region_dir: Path
    distance: np.ndarray
    methylation: np.ndarray
    kept_ids: list[str]
    screen_labels: list[str]
    n_reads_total: int
    n_alt_raw: int
    n_alt_after_peel: int


MODE_OUTPUT_SUFFIXES = [
    "null_mode",
    "base_seed",
    "observed_between_within",
    "null_threshold",
    "empirical_p_postselection_descriptive",
    "postselection_bh_q_descriptive",
    "postselection_by_q_descriptive",
    "postselection_fdr_calibrated",
    "n_valid_null",
    "minimum_valid_null",
    "null_valid_gate",
    "effect_gate",
    "null_threshold_gate",
    "multigroup_gate",
    "modal_groups",
    "modal_fraction",
    "unstable",
    "group_sizes",
    "minimum_group_size",
    "group_size_gate",
    "multiseed_pair_count",
    "multiseed_ari_median",
    "multiseed_ari_min",
    "multiseed_assignment_gate",
    "multiseed_exact_partition_gate",
    "screen_assignment_ari",
    "screen_assignment_concordance_gate",
    "screen_assignment_exact_partition_gate",
    "fixed_k_loo_partition_stability_feasible",
    "fixed_k_loo_partition_stability_evaluated",
    "fixed_k_loo_partition_stability_invalid",
    "fixed_k_loo_partition_stability_ari_median",
    "fixed_k_loo_partition_stability_ari_min",
    "fixed_k_loo_partition_stability_all_group_size_gate",
    "fixed_k_loo_partition_stability_ari_gate",
    "fixed_k_loo_partition_stability_exact_partition_count",
    "fixed_k_loo_partition_stability_exact_partition_gate",
    "failure",
]

STRICT_BASE_FIELDS = [
    "strict_source_candidate_row",
    "strict_analysis_class",
    "strict_robustness_target",
    "strict_formal_parameter_gate",
    "strict_formal_selection_contract_gate",
    "strict_formal_selection_column",
    "strict_selection_column",
    "strict_selection_source_column",
    "strict_selection_legacy_fallback_used",
    "strict_candidate_selection_gate",
    "strict_cooccurrence_receipt_gate",
    "strict_artifact_identity_gate",
    "strict_screen_artifact_hash_lock_gate",
    "strict_candidate_table_key_count",
    "strict_postselection_family_size",
    "strict_postselection_scope",
    "strict_postselection_fdr_calibrated",
    "strict_n_reads_total",
    "strict_n_alt_raw",
    "strict_n_alt_after_peel",
    "strict_minimum_reads_gate",
    "strict_combined_empirical_p_postselection_descriptive",
    "strict_postselection_bh_q_descriptive",
    "strict_postselection_by_q_descriptive",
    "strict_cross_null_assignment_ari",
    "strict_cross_null_exact_partition_gate",
    "strict_assignment_concordance_ari_min",
    "strict_methyl_partition_robustness_evaluable",
    "strict_null_robustness_pass",
    "strict_methyl_partition_robustness_status",
    "strict_scientific_status",
    "strict_not_evaluable_reason",
    "strict_failure_reason",
]

STRICT_MODE_FIELDS = [
    f"strict_{mode}_{suffix}" for mode in NULL_MODES for suffix in MODE_OUTPUT_SUFFIXES
]
STRICT_OUTPUT_FIELDS = STRICT_BASE_FIELDS + STRICT_MODE_FIELDS

STRICT_GATE_FIELDS = [
    "strict_formal_parameter_gate",
    "strict_formal_selection_contract_gate",
    "strict_candidate_selection_gate",
    "strict_cooccurrence_receipt_gate",
    "strict_artifact_identity_gate",
    "strict_screen_artifact_hash_lock_gate",
    "strict_minimum_reads_gate",
    "strict_column_null_valid_gate",
    "strict_column_effect_gate",
    "strict_column_null_threshold_gate",
    "strict_column_multigroup_gate",
    "strict_column_group_size_gate",
    "strict_column_multiseed_assignment_gate",
    "strict_column_multiseed_exact_partition_gate",
    "strict_column_screen_assignment_exact_partition_gate",
    "strict_column_fixed_k_loo_partition_stability_exact_partition_gate",
    "strict_row_circular_null_valid_gate",
    "strict_row_circular_effect_gate",
    "strict_row_circular_null_threshold_gate",
    "strict_row_circular_multigroup_gate",
    "strict_row_circular_group_size_gate",
    "strict_row_circular_multiseed_assignment_gate",
    "strict_row_circular_multiseed_exact_partition_gate",
    "strict_row_circular_screen_assignment_exact_partition_gate",
    "strict_row_circular_fixed_k_loo_partition_stability_exact_partition_gate",
    "strict_cross_null_exact_partition_gate",
]


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


def capture_source_identity() -> dict[str, dict[str, Any]]:
    return {role: artifact(path) for role, path in STRICT_CODE_PATHS.items()}


def capture_source_modes() -> dict[str, str]:
    return {
        role: oct(path.stat().st_mode & 0o777)
        for role, path in STRICT_CODE_PATHS.items()
    }


def capture_runtime_environment_identity() -> dict[str, Any]:
    implementation_modules = {
        "numpy_core_multiarray": importlib.import_module("numpy.core._multiarray_umath"),
        "numpy_random_generator": importlib.import_module("numpy.random._generator"),
        "scipy_hierarchy_extension": importlib.import_module(
            "scipy.cluster._hierarchy"
        ),
        "sklearn_pairwise_fast": importlib.import_module(
            "sklearn.metrics._pairwise_fast"
        ),
    }
    implementation_paths = {
        "scipy_fcluster_callable": Path(inspect.getfile(fcluster)).resolve(),
        "sklearn_adjusted_rand_score_callable": Path(
            inspect.getfile(adjusted_rand_score)
        ).resolve(),
        **{
            name: Path(module.__file__).resolve()
            for name, module in implementation_modules.items()
        },
    }
    return {
        "contract": "strict_python_executable_packages_and_direct_extensions_v1",
        "python_version": sys.version,
        "python_executable": artifact(Path(sys.executable)),
        "interpreter_flags": {
            "isolated": int(sys.flags.isolated),
            "no_user_site": int(sys.flags.no_user_site),
            "safe_path": bool(getattr(sys.flags, "safe_path", False)),
        },
        "python_cache": {
            "configured_prefix": sys.pycache_prefix,
            "xoption_prefix": sys._xoptions.get("pycache_prefix"),
        },
        "environment": {
            "PYTHONPATH": os.environ.get("PYTHONPATH"),
            "PYTHONNOUSERSITE": os.environ.get("PYTHONNOUSERSITE"),
            "PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED"),
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
            "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
            "NUMEXPR_NUM_THREADS": os.environ.get("NUMEXPR_NUM_THREADS"),
            "BLIS_NUM_THREADS": os.environ.get("BLIS_NUM_THREADS"),
        },
        "packages": {
            "numpy": {
                "version": np.__version__,
                "module": artifact(Path(np.__file__).resolve()),
            },
            "scipy": {
                "version": scipy.__version__,
                "module": artifact(Path(scipy.__file__).resolve()),
            },
            "scikit_learn": {
                "version": sklearn.__version__,
                "module": artifact(Path(sklearn.__file__).resolve()),
            },
        },
        "implementation_artifacts": {
            role: artifact(path)
            for role, path in sorted(implementation_paths.items())
        },
    }


def canonical_python_cache_root(output_dir: Path) -> Path:
    return output_dir.resolve().parent / CANONICAL_PYTHON_CACHE_DIRNAME


def canonical_python_prefix(output_dir: Path) -> list[str]:
    return [
        sys.executable,
        "-I",
        "-X",
        f"pycache_prefix={canonical_python_cache_root(output_dir)}",
    ]


def require_canonical_runtime(
    identity: Mapping[str, Any],
    *,
    output_dir: Path,
    formal_cache_required: bool = True,
) -> None:
    flags = identity["interpreter_flags"]
    environment = identity["environment"]
    python_cache = identity["python_cache"]
    expected_python_cache = str(canonical_python_cache_root(output_dir))
    if not (
        flags["isolated"] == 1
        and flags["no_user_site"] == 1
        and flags["safe_path"] is True
        and (
            not formal_cache_required
            or (
                python_cache["configured_prefix"] == expected_python_cache
                and python_cache["xoption_prefix"] == expected_python_cache
            )
        )
        and environment["PYTHONNOUSERSITE"] == "1"
        and environment["PYTHONPATH"] is None
        and environment["PYTHONHASHSEED"] == "0"
        and all(
            environment[name] == "1"
            for name in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
                "BLIS_NUM_THREADS",
            )
        )
    ):
        raise ConfirmationContractError(
            "Strict producer requires the canonical python -I -X pycache_prefix, "
            "no PYTHONPATH, PYTHONHASHSEED=0, PYTHONNOUSERSITE=1, and "
            "single-threaded numeric-library settings"
        )


def canonical_command(
    *,
    candidate_table: Path,
    assignments: Path,
    cooccurrence_receipt: Path,
    cooccurrence_release_receipt: Path,
    output_dir: Path,
    config: ConfirmationConfig,
) -> list[str]:
    return [
        *canonical_python_prefix(output_dir),
        str(Path(__file__).resolve()),
        "--candidate-table",
        str(candidate_table.resolve()),
        "--assignments",
        str(assignments.resolve()),
        "--cooccurrence-receipt",
        str(cooccurrence_receipt.resolve()),
        "--cooccurrence-release-receipt",
        str(cooccurrence_release_receipt.resolve()),
        "--output-dir",
        str(output_dir.resolve()),
        "--permutations",
        str(config.permutations),
        "--seeds",
        str(config.seeds),
        "--seed",
        str(config.seed),
        "--assignment-ari-threshold",
        str(config.assignment_ari_threshold),
        "--min-valid-null-fraction",
        str(config.min_valid_null_fraction),
    ]


def observed_process_command() -> list[str]:
    raw = Path("/proc/self/cmdline").read_bytes()
    if not raw.endswith(b"\0"):
        raise ConfirmationContractError("Process command line is unavailable or malformed")
    return [os.fsdecode(token) for token in raw[:-1].split(b"\0")]


def compact_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def open_text_auto(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def parse_boolean(value: Any, *, field: str, row_number: int | None = None) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "t", "yes", "y"}:
        return True
    if normalized in {"0", "false", "f", "no", "n", ""}:
        return False
    location = f" at row {row_number}" if row_number is not None else ""
    raise ConfirmationContractError(f"Invalid boolean in {field}{location}: {value!r}")


def validate_config(config: ConfirmationConfig) -> None:
    if config.permutations < 1:
        raise ConfirmationContractError("--permutations must be >= 1")
    if config.seeds < 2:
        raise ConfirmationContractError("--seeds must be >= 2 for multi-seed ARI")
    if not 0.0 <= config.assignment_ari_threshold <= 1.0:
        raise ConfirmationContractError("--assignment-ari-threshold must be in [0, 1]")
    if not 0.0 < config.min_valid_null_fraction <= 1.0:
        raise ConfirmationContractError("--min-valid-null-fraction must be in (0, 1]")


def formal_parameter_gate(config: ConfirmationConfig) -> bool:
    return bool(
        config.permutations == DEFAULT_PERMUTATIONS
        and config.seeds == DEFAULT_SEEDS
        and config.seed == DEFAULT_SEED
        and math.isclose(
            config.assignment_ari_threshold,
            DEFAULT_ASSIGNMENT_ARI,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        and math.isclose(
            config.min_valid_null_fraction,
            DEFAULT_MIN_VALID_NULL_FRACTION,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    )


def resolve_selection_column(fields: Sequence[str]) -> SelectionColumnContract:
    available = set(fields)
    if FORMAL_SELECTION_COLUMN in available:
        return SelectionColumnContract(
            formal_column=FORMAL_SELECTION_COLUMN,
            resolved_column=FORMAL_SELECTION_COLUMN,
            legacy_fallback_used=False,
        )
    for legacy_column in LEGACY_SELECTION_COLUMNS:
        if legacy_column in available:
            return SelectionColumnContract(
                formal_column=FORMAL_SELECTION_COLUMN,
                resolved_column=legacy_column,
                legacy_fallback_used=True,
            )
    raise ConfirmationContractError(
        "Candidate table lacks the formal selection column and supported legacy aliases: "
        f"formal={FORMAL_SELECTION_COLUMN!r} legacy={list(LEGACY_SELECTION_COLUMNS)!r}"
    )


def candidate_table_fields(path: Path) -> list[str]:
    if not path.is_file():
        raise ConfirmationContractError(f"Candidate table is not a file: {path}")
    with open_text_auto(path) as handle:
        fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    if not fields or len(fields) != len(set(fields)):
        raise ConfirmationContractError("Candidate table has missing or duplicate header fields")
    return fields


def candidate_table_keys(path: Path) -> list[tuple[str, str, int]]:
    keys: list[tuple[str, str, int]] = []
    with open_text_auto(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=2):
            keys.append(site_key(row, source=f"candidate row {row_number}"))
    if len(keys) != len(set(keys)):
        raise ConfirmationContractError("Candidate table contains duplicate site keys")
    return keys


def scientific_analysis_class(
    config: ConfirmationConfig, selection: SelectionColumnContract
) -> str:
    if formal_parameter_gate(config) and not selection.legacy_fallback_used:
        return "FORMAL"
    return "EXPLORATORY_ONLY"


def site_key(record: dict[str, Any], *, source: str) -> tuple[str, str, int]:
    try:
        sample = str(record["sample"]).strip()
        chrom = str(record["chrom"]).strip()
        pos = int(record["pos"])
    except (KeyError, TypeError, ValueError) as error:
        raise ConfirmationContractError(f"Invalid site key in {source}: {record!r}") from error
    if not sample or not chrom or pos <= 0:
        raise ConfirmationContractError(f"Invalid site key in {source}: {(sample, chrom, pos)!r}")
    return sample, chrom, pos


def load_candidates(
    path: Path, selection_column: str
) -> tuple[list[str], list[dict[str, str]], int]:
    if not path.is_file():
        raise ConfirmationContractError(f"Candidate table is not a file: {path}")
    with open_text_auto(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if not fields or len(fields) != len(set(fields)):
            raise ConfirmationContractError("Candidate table has missing or duplicate header fields")
        required = {
            "sample",
            "chrom",
            "pos",
            "ref",
            "alt",
            "region_dir",
            "alt_readset_sha256",
            "screen_contract",
            selection_column,
        }
        missing = sorted(required.difference(fields))
        if missing:
            raise ConfirmationContractError(f"Candidate table missing required fields: {missing}")
        collisions = sorted(set(fields).intersection(STRICT_OUTPUT_FIELDS))
        if collisions:
            raise ConfirmationContractError(
                f"Candidate table collides with reserved strict output fields: {collisions}"
            )

        selected: list[dict[str, str]] = []
        seen: set[tuple[str, str, int]] = set()
        n_input = 0
        for row_number, row in enumerate(reader, start=2):
            n_input += 1
            key = site_key(row, source=f"candidate row {row_number}")
            if key in seen:
                raise ConfirmationContractError(f"Duplicate candidate site key: {key}")
            seen.add(key)
            if not str(row.get("region_dir", "")).strip():
                raise ConfirmationContractError(f"Empty region_dir at candidate row {row_number}")
            if str(row.get("screen_contract", "")).strip() != SCREEN_CONTRACT:
                raise ConfirmationContractError(
                    f"Unexpected screen_contract at candidate row {row_number}"
                )
            if not str(row.get("alt_readset_sha256", "")).strip():
                raise ConfirmationContractError(
                    f"Empty alt_readset_sha256 at candidate row {row_number}"
                )
            for allele in ("ref", "alt"):
                if not str(row.get(allele, "")).strip():
                    raise ConfirmationContractError(f"Empty {allele} at candidate row {row_number}")
            if parse_boolean(row[selection_column], field=selection_column, row_number=row_number):
                copy = dict(row)
                copy["__source_row_number"] = str(row_number)
                selected.append(copy)
    return fields, selected, n_input


def load_assignments(path: Path) -> dict[tuple[str, str, int], dict[str, Any]]:
    if not path.is_file():
        raise ConfirmationContractError(f"Assignments input is not a file: {path}")
    records: dict[tuple[str, str, int], dict[str, Any]] = {}
    with open_text_auto(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                raise ConfirmationContractError(
                    f"Malformed assignments JSON at line {line_number}: {error}"
                ) from error
            if not isinstance(record, dict):
                raise ConfirmationContractError(
                    f"Assignments line {line_number} is not a JSON object"
                )
            if (
                record.get("schema_name") != ASSIGNMENT_SCHEMA
                or record.get("schema_version") != ASSIGNMENT_SCHEMA_VERSION
            ):
                raise ConfirmationContractError(
                    f"Unexpected assignment schema at line {line_number}"
                )
            if record.get("screen_contract") != SCREEN_CONTRACT:
                raise ConfirmationContractError(
                    f"Unexpected assignment screen contract at line {line_number}"
                )
            if record.get("artifact_identity_contract") != ARTIFACT_IDENTITY_CONTRACT:
                raise ConfirmationContractError(
                    f"Unexpected assignment artifact identity contract at line {line_number}"
                )
            key = site_key(record, source=f"assignments line {line_number}")
            if key in records:
                raise ConfirmationContractError(f"Duplicate assignments site key: {key}")
            records[key] = record
    return records


def load_reads_strict(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    required = {
        "read_id",
        "read_name",
        "chr",
        "start",
        "end",
        "mapq",
        "alt_support",
        "is_tumor",
        "strand",
    }
    if not path.is_file():
        raise ConfirmationContractError(f"Missing primary reads artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = sorted(required.difference(fields))
        if missing:
            raise ConfirmationContractError(f"reads.tsv missing fields {missing}: {path}")
        order: list[str] = []
        rows: dict[str, dict[str, str]] = {}
        for line_number, row in enumerate(reader, start=2):
            read_id = str(row["read_id"])
            if read_id in rows:
                raise ConfirmationContractError(
                    f"Duplicate read_id {read_id!r} in {path} line {line_number}"
                )
            try:
                int(row["start"])
                int(row["end"])
                int(row["mapq"])
            except ValueError as error:
                raise ConfirmationContractError(
                    f"Invalid numeric read identity in {path} line {line_number}"
                ) from error
            order.append(read_id)
            rows[read_id] = dict(row)
    if not order:
        raise ConfirmationContractError(f"Empty primary reads artifact: {path}")
    return order, rows


def load_csv_matrix_strict(
    path: Path, *, square: bool
) -> tuple[list[str], list[str], np.ndarray]:
    if not path.is_file():
        raise ConfirmationContractError(f"Missing primary matrix artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows or len(rows[0]) < 2 or rows[0][0] != "read_id":
        raise ConfirmationContractError(f"Invalid matrix header: {path}")
    columns = rows[0][1:]
    if len(columns) != len(set(columns)):
        raise ConfirmationContractError(f"Duplicate matrix columns: {path}")
    identifiers: list[str] = []
    values: list[list[float]] = []
    for line_number, row in enumerate(rows[1:], start=2):
        if not row:
            continue
        if len(row) != len(rows[0]):
            raise ConfirmationContractError(
                f"Matrix width mismatch in {path} line {line_number}: {len(row)} != {len(rows[0])}"
            )
        read_id = row[0]
        if read_id in identifiers:
            raise ConfirmationContractError(f"Duplicate matrix read_id {read_id!r}: {path}")
        try:
            parsed = [F.parse_float(value) for value in row[1:]]
        except ValueError as error:
            raise ConfirmationContractError(
                f"Invalid matrix value in {path} line {line_number}"
            ) from error
        identifiers.append(read_id)
        values.append(parsed)
    matrix = np.asarray(values, dtype=float)
    expected_shape = (len(identifiers), len(columns))
    if matrix.shape != expected_shape:
        raise ConfirmationContractError(
            f"Matrix shape mismatch in {path}: {matrix.shape} != {expected_shape}"
        )
    if square:
        if identifiers != columns or matrix.shape[0] != matrix.shape[1]:
            raise ConfirmationContractError(f"BERNOULLI row/column read identity mismatch: {path}")
        if np.isinf(matrix).any() or np.nanmin(matrix) < -1.0:
            raise ConfirmationContractError(f"Invalid BERNOULLI values: {path}")
        if not np.allclose(matrix, matrix.T, atol=1e-6, rtol=0.0, equal_nan=True):
            raise ConfirmationContractError(f"Asymmetric BERNOULLI matrix: {path}")
        if not np.allclose(np.diag(matrix), 0.0, atol=1e-6, rtol=0.0):
            raise ConfirmationContractError(f"Non-zero BERNOULLI diagonal: {path}")
    else:
        finite = matrix[np.isfinite(matrix)]
        if finite.size and (float(finite.min()) < 0.0 or float(finite.max()) > 1.0):
            raise ConfirmationContractError(f"Methylation values outside [0, 1]: {path}")
    return identifiers, columns, matrix


def readset_digest(rows: Iterable[dict[str, str]]) -> str:
    identities = sorted(
        f"{row['read_name']}|{row['chr']}|{row['start']}|{row['end']}|{row['strand']}"
        for row in rows
    )
    return hashlib.sha256("\n".join(identities).encode()).hexdigest()


def verify_artifact_identity(record: Any, expected_path: Path, *, label: str) -> None:
    if not isinstance(record, dict):
        raise ConfirmationContractError(f"Missing {label} artifact identity record")
    required = {"path", "size_bytes", "sha256"}
    missing = sorted(required.difference(record))
    if missing:
        raise ConfirmationContractError(f"{label} artifact identity missing fields: {missing}")
    observed_path = Path(str(record["path"])).expanduser().resolve()
    expected_path = expected_path.expanduser().resolve()
    if observed_path != expected_path:
        raise ConfirmationContractError(
            f"{label} artifact path mismatch: {observed_path} != {expected_path}"
        )
    if not expected_path.is_file():
        raise ConfirmationContractError(f"Missing {label} artifact: {expected_path}")
    if int(record["size_bytes"]) != expected_path.stat().st_size:
        raise ConfirmationContractError(f"{label} artifact size mismatch: {expected_path}")
    if str(record["sha256"]) != sha256(expected_path):
        raise ConfirmationContractError(f"{label} artifact sha256 mismatch: {expected_path}")


def load_json_object(path: Path, *, label: str) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ConfirmationContractError(f"Missing or empty {label}: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ConfirmationContractError(f"Malformed {label}: {path}: {error}") from error
    if not isinstance(payload, dict):
        raise ConfirmationContractError(f"{label} must be a JSON object: {path}")
    return payload


def require_count(container: Any, names: Sequence[str], *, label: str) -> tuple[int, str]:
    if not isinstance(container, dict):
        raise ConfirmationContractError(f"Missing count container for {label}")
    for name in names:
        if name not in container:
            continue
        try:
            value = int(container[name])
        except (TypeError, ValueError) as error:
            raise ConfirmationContractError(f"Invalid {label} count in {name}") from error
        if value < 0:
            raise ConfirmationContractError(f"Negative {label} count in {name}")
        return value, name
    raise ConfirmationContractError(
        f"Missing {label} count; expected one of {list(names)!r}"
    )


def selection_family_count(
    receipt: dict[str, Any],
    summary: dict[str, Any],
    selection: SelectionColumnContract,
) -> tuple[int, str]:
    if selection.legacy_fallback_used:
        names = (
            "n_genetically_anchored_multi_marker_candidates_by_sensitivity",
            "candidate_family_size",
            "postselection_family_size",
        )
    else:
        names = (
            "n_multi_marker_molecular_haplotype_base_candidates",
            "multi_marker_molecular_haplotype_base_candidates",
            "multi_marker_molecular_haplotype_base_candidate_family_size",
            "candidate_family_size",
            "postselection_family_size",
        )
    observed: list[tuple[str, int]] = []
    for container_name, container in (
        ("receipt.counts", receipt.get("counts")),
        ("summary", summary),
    ):
        if not isinstance(container, dict):
            continue
        for name in names:
            if name not in container:
                continue
            value, _ = require_count(container, (name,), label="selected candidate family")
            observed.append((f"{container_name}.{name}", value))
    if not observed:
        raise ConfirmationContractError(
            "Cooccurrence receipt/summary does not record the selected candidate family size "
            f"for source column {selection.resolved_column!r}"
        )
    if len({value for _, value in observed}) != 1:
        raise ConfirmationContractError(
            f"Conflicting cooccurrence selected candidate family counts: {observed!r}"
        )
    return observed[0][1], ",".join(name for name, _ in observed)


def validate_cooccurrence_receipt(
    receipt_path: Path,
    *,
    release_receipt_path: Path | None = None,
    candidate_table: Path,
    assignments_path: Path,
    n_candidate_rows: int,
    candidate_keys: Sequence[tuple[str, str, int]],
    selected_keys: Sequence[tuple[str, str, int]],
    selection: SelectionColumnContract,
) -> dict[str, Any]:
    current_source_authority = SOURCE_AUTHORITY.validate_release_source_authority(
        {
            "release_source_authority_validator",
            "cooccurrence_analyzer",
            "cooccurrence_preflight",
            "cooccurrence_release_finalizer",
            "cooccurrence_release_runner",
            "ssnv_cooccurrence_lib",
            "latest_tag_join",
            "m2_screen_gate",
            "independent_m2_auditor",
            "primary_artifact_auditor",
            "strict_producer",
            "focal_alt_cluster_lib",
        }
    )
    receipt = load_json_object(receipt_path, label="cooccurrence receipt")
    if release_receipt_path is None:
        release_receipt_path = (
            receipt_path.parent / "cooccurrence_release_receipt.json"
        )
    release_receipt = load_json_object(
        release_receipt_path,
        label="cooccurrence release receipt",
    )
    release_inputs = release_receipt.get("inputs")
    release_checks = release_receipt.get("checks")
    if (
        release_receipt.get("schema_name") != COOCCURRENCE_RELEASE_SCHEMA
        or release_receipt.get("schema_version") != COOCCURRENCE_RELEASE_SCHEMA_VERSION
        or release_receipt.get("task_type") != "B_comprehensive_validation"
        or release_receipt.get("scope") != "all_manifest_samples"
        or release_receipt.get("raw_identity_release_contract")
        != RAW_IDENTITY_RELEASE_CONTRACT
        or release_receipt.get("release_status")
        != "RELEASE_RECONCILIATION_PASS"
        or release_receipt.get("pass") is not True
        or not isinstance(release_checks, dict)
        or not release_checks
        or not all(value is True for value in release_checks.values())
        or not isinstance(release_inputs, dict)
    ):
        raise ConfirmationContractError("Cooccurrence release receipt did not pass")
    if release_receipt.get("source_authority") != current_source_authority:
        raise ConfirmationContractError("Cooccurrence release source authority drift")
    verify_artifact_identity(
        release_inputs.get("producer_receipt"),
        receipt_path,
        label="cooccurrence released producer receipt",
    )
    if oct(release_receipt_path.stat().st_mode & 0o777) != "0o444":
        raise ConfirmationContractError("Cooccurrence release receipt is not mode 0444")
    release_code = release_receipt.get("code")
    if (
        not isinstance(release_code, dict)
        or set(release_code) != set(COOCCURRENCE_RELEASE_CODE_PATHS)
        or release_receipt.get("source_modes")
        != {role: "0o444" for role in COOCCURRENCE_RELEASE_CODE_PATHS}
    ):
        raise ConfirmationContractError("Cooccurrence release code contract drift")
    for role, source_path in COOCCURRENCE_RELEASE_CODE_PATHS.items():
        verify_artifact_identity(
            release_code.get(role),
            source_path,
            label=f"cooccurrence release code {role}",
        )
    if receipt.get("schema_name") != COOCCURRENCE_RECEIPT_SCHEMA:
        raise ConfirmationContractError("Unexpected cooccurrence receipt schema")
    if (
        not selection.legacy_fallback_used
        and receipt.get("schema_version") != COOCCURRENCE_FORMAL_SCHEMA_VERSION
    ):
        raise ConfirmationContractError(
            "Formal molecular-haplotype selection requires cooccurrence schema "
            f"{COOCCURRENCE_FORMAL_SCHEMA_VERSION}"
        )
    if receipt.get("pass") is not True:
        raise ConfirmationContractError("Cooccurrence receipt is not passing")
    if receipt.get("source_authority") != current_source_authority:
        raise ConfirmationContractError("Cooccurrence producer source authority drift")
    receipt_scope_sources: list[str] = []
    if receipt.get("task_type") != "B_comprehensive_validation":
        raise ConfirmationContractError("Cooccurrence receipt task_type is not full scope")
    receipt_scope_sources.append("receipt.task_type")
    if receipt.get("full_scope") is not True:
        raise ConfirmationContractError("Cooccurrence receipt full_scope is not true")
    receipt_scope_sources.append("receipt.full_scope")
    receipt_scope = receipt.get("scope")
    if isinstance(receipt_scope, str):
        if receipt_scope != "all_manifest_samples":
            raise ConfirmationContractError("Cooccurrence receipt scope is not all samples")
        receipt_scope_sources.append("receipt.scope")
    elif isinstance(receipt_scope, dict):
        if receipt_scope.get("full_scope") is not True:
            raise ConfirmationContractError("Cooccurrence receipt scope.full_scope is not true")
        receipt_scope_sources.append("receipt.scope.full_scope")
    else:
        raise ConfirmationContractError("Cooccurrence receipt scope is missing")

    outputs = receipt.get("outputs")
    inputs = receipt.get("inputs")
    if not isinstance(outputs, dict) or not isinstance(inputs, dict):
        raise ConfirmationContractError("Cooccurrence receipt lacks inputs/outputs")
    verify_artifact_identity(outputs.get("site_tsv"), candidate_table, label="cooccurrence site_tsv")
    verify_artifact_identity(
        inputs.get("assignments"), assignments_path, label="cooccurrence assignments"
    )
    if receipt.get("raw_identity_release_contract") != RAW_IDENTITY_RELEASE_CONTRACT:
        raise ConfirmationContractError("Cooccurrence raw-identity release contract drift")
    for output_name, label in (
        ("raw_identity_duplicate_audit_tsv", "raw identity duplicate audit"),
        ("case_json", "oracle cases"),
    ):
        record = outputs.get(output_name)
        if not isinstance(record, dict) or not record.get("path"):
            raise ConfirmationContractError(f"Cooccurrence receipt lacks {label}")
        verify_artifact_identity(
            record,
            Path(str(record["path"])),
            label=f"cooccurrence {label}",
        )
    preflight_record = inputs.get("preflight_receipt")
    if not isinstance(preflight_record, dict) or not preflight_record.get("path"):
        raise ConfirmationContractError("Cooccurrence receipt lacks preflight identity")
    verify_artifact_identity(
        preflight_record,
        Path(str(preflight_record["path"])),
        label="cooccurrence preflight",
    )
    preflight = load_json_object(
        Path(str(preflight_record["path"])), label="cooccurrence preflight"
    )
    preflight_raw_identity = (preflight.get("observed") or {}).get(
        "raw_bam_identity_recovery"
    )
    if (
        preflight.get("schema_name")
        != "intersubmod.cooccurrence_task_contract_preflight"
        or preflight.get("schema_version") != "2.0.0"
        or preflight.get("pass") is not True
        or preflight.get("raw_identity_release_contract")
        != RAW_IDENTITY_RELEASE_CONTRACT
        or int((preflight.get("observed") or {}).get("task_count", -1)) != 102_842
        or not isinstance(preflight_raw_identity, dict)
        or int((preflight_raw_identity.get("aggregate") or {}).get("site_tasks", -1))
        != 102_842
        or preflight_raw_identity.get("equivalence_contract")
        != RAW_DUPLICATE_EQUIVALENCE_CONTRACT
        or preflight_raw_identity.get("allowed_differing_auxiliary_tags") != ["RG"]
        or preflight_raw_identity.get("analysis_scope_identity_contract")
        != ANALYSIS_SCOPE_IDENTITY_CONTRACT
        or preflight_raw_identity.get("missing_projection_policy")
        != RAW_IDENTITY_MISSING_POLICY
        or preflight_raw_identity.get("conflicting_analysis_payload_policy")
        != RAW_IDENTITY_CONFLICT_POLICY
        or preflight_raw_identity.get("failure_counts_materialized") is not False
        or preflight_raw_identity.get(
            "all_result_rows_passed_invariant_validation"
        )
        is not True
    ):
        raise ConfirmationContractError("Cooccurrence preflight gate drift")
    if preflight.get("source_authority") != current_source_authority:
        raise ConfirmationContractError("Cooccurrence preflight source authority drift")
    preflight_gate = receipt.get("preflight_gate")
    if (
        not isinstance(preflight_gate, dict)
        or preflight_gate.get("schema_name")
        != "intersubmod.cooccurrence_task_contract_preflight"
        or preflight_gate.get("schema_version") != "2.0.0"
        or int(preflight_gate.get("task_count", -1)) != 102_842
        or preflight_gate.get("pass") is not True
        or preflight_gate.get("raw_identity_release_contract")
        != RAW_IDENTITY_RELEASE_CONTRACT
    ):
        raise ConfirmationContractError("Cooccurrence receipt preflight gate drift")

    code = receipt.get("code")
    source_lock = receipt.get("source_lock")
    if not isinstance(code, dict) or set(code) != set(COOCCURRENCE_CODE_PATHS):
        raise ConfirmationContractError("Cooccurrence code role set is not exact")
    for role, source_path in COOCCURRENCE_CODE_PATHS.items():
        verify_artifact_identity(
            code.get(role), source_path, label=f"cooccurrence code {role}"
        )
    if not isinstance(source_lock, dict):
        raise ConfirmationContractError("Cooccurrence receipt lacks source lock")
    for field in (
        "source_identity_before",
        "source_identity_after_compute",
        "source_identity_after_output",
    ):
        if source_lock.get(field) != code:
            raise ConfirmationContractError(f"Cooccurrence source identity drift: {field}")
    for field in (
        "source_modes_before",
        "source_modes_after_compute",
        "source_modes_after_output",
    ):
        modes = source_lock.get(field)
        if (
            not isinstance(modes, dict)
            or set(modes) != set(COOCCURRENCE_CODE_PATHS)
            or set(modes.values()) != {"0o444"}
        ):
            raise ConfirmationContractError(f"Cooccurrence source mode drift: {field}")
    if source_lock.get("all_sources_read_only_and_unchanged") is not True:
        raise ConfirmationContractError("Cooccurrence source-lock pass marker is absent")

    preflight_code = preflight.get("code")
    if not isinstance(preflight_code, dict):
        raise ConfirmationContractError("Cooccurrence preflight lacks source lock")
    preflight_before = preflight_code.get("source_identity_before")
    preflight_after = preflight_code.get("source_identity_after")
    expected_preflight_roles = {"preflight", *COOCCURRENCE_CODE_PATHS}
    if (
        not isinstance(preflight_before, dict)
        or not isinstance(preflight_after, dict)
        or set(preflight_before) != expected_preflight_roles
        or preflight_before != preflight_after
    ):
        raise ConfirmationContractError("Cooccurrence preflight source identity drift")
    verify_artifact_identity(
        preflight_before.get("preflight"),
        COOCCURRENCE_PREFLIGHT_SOURCE,
        label="cooccurrence preflight source",
    )
    for role in COOCCURRENCE_CODE_PATHS:
        if preflight_before.get(role) != code.get(role):
            raise ConfirmationContractError(
                f"Cooccurrence preflight/producer source drift: {role}"
            )
    for field in ("source_modes_before", "source_modes_after"):
        modes = preflight_code.get(field)
        if (
            not isinstance(modes, dict)
            or set(modes) != expected_preflight_roles
            or set(modes.values()) != {"0o444"}
        ):
            raise ConfirmationContractError(
                f"Cooccurrence preflight source mode drift: {field}"
            )
    counts = receipt.get("counts")
    if (
        not isinstance(counts, dict)
        or counts.get("raw_identity_missing_projection_policy")
        != RAW_IDENTITY_MISSING_POLICY
        or counts.get("raw_identity_conflicting_analysis_payload_policy")
        != RAW_IDENTITY_CONFLICT_POLICY
        or counts.get("raw_identity_failure_counts_materialized") is not False
        or counts.get("raw_identity_all_site_results_passed_invariant_validation")
        is not True
    ):
        raise ConfirmationContractError("Cooccurrence raw-identity receipt policy drift")

    summary_record = outputs.get("summary_json")
    if not isinstance(summary_record, dict) or not summary_record.get("path"):
        raise ConfirmationContractError("Cooccurrence receipt lacks summary_json identity")
    summary_path = Path(str(summary_record["path"])).expanduser().resolve()
    verify_artifact_identity(summary_record, summary_path, label="cooccurrence summary_json")
    summary = load_json_object(summary_path, label="cooccurrence summary")
    if summary.get("schema_name") != COOCCURRENCE_SUMMARY_SCHEMA:
        raise ConfirmationContractError("Unexpected cooccurrence summary schema")
    if summary.get("schema_version") != receipt.get("schema_version"):
        raise ConfirmationContractError("Cooccurrence receipt/summary schema-version mismatch")
    raw_summary = summary.get("raw_bam_identity_recovery_audit")
    if (
        summary.get("raw_identity_release_contract")
        != RAW_IDENTITY_RELEASE_CONTRACT
        or not isinstance(raw_summary, dict)
        or raw_summary.get("equivalence_contract")
        != RAW_DUPLICATE_EQUIVALENCE_CONTRACT
        or raw_summary.get("allowed_differing_auxiliary_tags") != ["RG"]
        or raw_summary.get("analysis_scope_identity_contract")
        != ANALYSIS_SCOPE_IDENTITY_CONTRACT
        or raw_summary.get("missing_projection_policy")
        != RAW_IDENTITY_MISSING_POLICY
        or raw_summary.get("conflicting_analysis_payload_policy")
        != RAW_IDENTITY_CONFLICT_POLICY
        or raw_summary.get("failure_counts_materialized") is not False
        or raw_summary.get("all_site_results_passed_invariant_validation") is not True
    ):
        raise ConfirmationContractError("Cooccurrence raw-identity summary policy drift")
    geometry = summary.get("partner_geometry_audit")
    full_scope_gate = bool(
        summary.get("task_type") == "B_comprehensive_validation"
        and summary.get("scope") == "all_manifest_samples"
        and isinstance(geometry, dict)
        and geometry.get("full_scope") is True
    )
    if not full_scope_gate:
        raise ConfirmationContractError(
            "Cooccurrence receipt/summary is not Task-B full scope"
        )

    unique_candidate_keys = set(candidate_keys)
    unique_selected_keys = set(selected_keys)
    if len(candidate_keys) != n_candidate_rows or len(unique_candidate_keys) != n_candidate_rows:
        raise ConfirmationContractError(
            "Candidate row/key count mismatch or duplicate site keys"
        )
    if len(unique_selected_keys) != len(selected_keys):
        raise ConfirmationContractError("Selected candidate key count is not unique")

    receipt_rows, receipt_rows_field = require_count(
        receipt.get("counts"),
        (
            "candidate_table_rows",
            "candidate_row_count",
            "stable_site_rows",
            "stable_sites",
        ),
        label="candidate table rows",
    )
    summary_rows, summary_rows_field = require_count(
        summary,
        ("candidate_table_rows", "candidate_row_count", "n_stable_sites"),
        label="candidate table rows",
    )
    if receipt_rows != n_candidate_rows or summary_rows != n_candidate_rows:
        raise ConfirmationContractError(
            "Cooccurrence candidate row count mismatch: "
            f"receipt={receipt_rows} summary={summary_rows} observed={n_candidate_rows}"
        )

    explicit_key_counts: list[tuple[str, int]] = []
    for container_name, container in (
        ("receipt.counts", receipt.get("counts")),
        ("summary", summary),
    ):
        if not isinstance(container, dict):
            continue
        for field in (
            "candidate_unique_key_count",
            "candidate_key_count",
            "stable_site_unique_key_count",
        ):
            if field in container:
                try:
                    explicit_key_counts.append((f"{container_name}.{field}", int(container[field])))
                except (TypeError, ValueError) as error:
                    raise ConfirmationContractError(
                        f"Invalid candidate key count in {container_name}.{field}"
                    ) from error
    if any(value != len(unique_candidate_keys) for _, value in explicit_key_counts):
        raise ConfirmationContractError(
            "Cooccurrence candidate unique-key count mismatch: "
            f"recorded={explicit_key_counts!r} observed={len(unique_candidate_keys)}"
        )

    expected_family_size, family_count_source = selection_family_count(
        receipt, summary, selection
    )
    if expected_family_size != len(selected_keys) or len(unique_selected_keys) != len(selected_keys):
        raise ConfirmationContractError(
            "Cooccurrence selected candidate family mismatch: "
            f"receipt={expected_family_size} rows={len(selected_keys)} "
            f"keys={len(unique_selected_keys)}"
        )

    return {
        "receipt": artifact(receipt_path),
        "release_receipt": artifact(release_receipt_path),
        "schema_name": receipt["schema_name"],
        "schema_version": receipt.get("schema_version"),
        "full_scope_gate": True,
        "full_scope_sources": [
            "summary.task_type",
            "summary.scope",
            "summary.partner_geometry_audit.full_scope",
            *receipt_scope_sources,
        ],
        "site_tsv": artifact(candidate_table),
        "assignments": artifact(assignments_path),
        "summary_json": artifact(summary_path),
        "candidate_table_row_count": n_candidate_rows,
        "candidate_table_unique_key_count": len(unique_candidate_keys),
        "candidate_row_count_sources": [
            f"receipt.counts.{receipt_rows_field}",
            f"summary.{summary_rows_field}",
        ],
        "explicit_candidate_key_count_sources": [name for name, _ in explicit_key_counts],
        "selected_candidate_row_count": len(selected_keys),
        "selected_candidate_unique_key_count": len(unique_selected_keys),
        "postselection_family_size": expected_family_size,
        "family_count_source": family_count_source,
        "screen_assignment_artifact_hash_lock": True,
    }


def verify_assignment_artifacts_and_site_identity(
    assignment: dict[str, Any],
    candidate: dict[str, str],
    region_dir: Path,
    reads_path: Path,
    distance_path: Path,
    methylation_path: Path,
) -> None:
    posthoc = assignment.get("posthoc")
    if not isinstance(posthoc, dict):
        raise ConfirmationContractError("Assignment posthoc site identity is missing")
    for field in ("ref", "alt"):
        if str(posthoc.get(field, "")).upper() != str(candidate[field]).upper():
            raise ConfirmationContractError(
                f"Assignment/candidate {field} mismatch for {site_key(candidate, source='candidate')}"
            )
    if assignment.get("screen_contract") != candidate.get("screen_contract"):
        raise ConfirmationContractError("Assignment/candidate screen_contract mismatch")
    assignment_region = assignment.get("region_dir")
    if not assignment_region or Path(str(assignment_region)).resolve() != region_dir.resolve():
        raise ConfirmationContractError(
            f"Assignment region_dir mismatch for {site_key(candidate, source='candidate')}"
        )
    artifacts = assignment.get("primary_artifacts")
    if not isinstance(artifacts, dict):
        raise ConfirmationContractError("Assignment primary_artifacts is missing")
    verify_artifact_identity(artifacts.get("reads"), reads_path, label="reads")
    verify_artifact_identity(
        artifacts.get("distance_matrix"), distance_path, label="distance_matrix"
    )
    verify_artifact_identity(
        artifacts.get("methylation_matrix"), methylation_path, label="methylation_matrix"
    )


def require_exact_list(actual: Any, expected: list[str], *, label: str) -> list[str]:
    if not isinstance(actual, list) or any(not isinstance(value, str) for value in actual):
        raise ConfirmationContractError(f"{label} must be a list of strings")
    if actual != expected:
        raise ConfirmationContractError(
            f"{label} identity/order mismatch: assignment={actual[:8]!r} primary={expected[:8]!r}"
        )
    return actual


def validate_assignment_identity(
    assignment: dict[str, Any],
    candidate: dict[str, str],
    region_dir: Path,
    reads: dict[str, dict[str, str]],
    kept_ids: list[str],
) -> list[str]:
    assignment_region = assignment.get("region_dir")
    if not assignment_region or Path(str(assignment_region)).resolve() != region_dir.resolve():
        raise ConfirmationContractError(
            f"Assignment region_dir mismatch for {site_key(candidate, source='candidate')}"
        )
    if "strict_confirm_candidate" in assignment and not parse_boolean(
        assignment["strict_confirm_candidate"], field="assignments.strict_confirm_candidate"
    ):
        raise ConfirmationContractError("Selected candidate is not screen-positive in assignments")

    require_exact_list(
        assignment.get("all_after_peel_read_ids"),
        kept_ids,
        label="assignments.all_after_peel_read_ids",
    )
    screen_labels = assignment.get("coarse_labels_all_after_peel")
    if not isinstance(screen_labels, list) or len(screen_labels) != len(kept_ids):
        raise ConfirmationContractError(
            "assignments.coarse_labels_all_after_peel length/type mismatch"
        )
    screen_labels = [str(value) for value in screen_labels]
    expected_core_ids = [
        read_id
        for read_id, label in zip(kept_ids, screen_labels)
        if label not in {"other", "outlier"}
    ]
    require_exact_list(assignment.get("read_ids"), expected_core_ids, label="assignments.read_ids")
    expected_core_labels = [
        label for label in screen_labels if label not in {"other", "outlier"}
    ]
    require_exact_list(
        assignment.get("labels"), expected_core_labels, label="assignments.labels"
    )
    expected_names = [reads[read_id]["read_name"] for read_id in expected_core_ids]
    require_exact_list(
        assignment.get("read_names"), expected_names, label="assignments.read_names"
    )

    core_reads = assignment.get("core_reads")
    if not isinstance(core_reads, list) or len(core_reads) != len(expected_core_ids):
        raise ConfirmationContractError("assignments.core_reads length/type mismatch")
    for index, (record, read_id, label) in enumerate(
        zip(core_reads, expected_core_ids, expected_core_labels)
    ):
        if not isinstance(record, dict):
            raise ConfirmationContractError(f"assignments.core_reads[{index}] is not an object")
        primary = reads[read_id]
        expected = {
            "read_id": read_id,
            "read_name": primary["read_name"],
            "chrom": primary["chr"],
            "start": int(primary["start"]),
            "end": int(primary["end"]),
            "mapq": int(primary["mapq"]),
            "strand": primary["strand"],
            "label": label,
        }
        for field, value in expected.items():
            if record.get(field) != value:
                raise ConfirmationContractError(
                    f"assignments.core_reads[{index}].{field} identity mismatch: "
                    f"{record.get(field)!r} != {value!r}"
                )
    return screen_labels


def prepare_site(
    candidate: dict[str, str], assignment: dict[str, Any]
) -> PreparedSite:
    sample, chrom, pos = site_key(candidate, source="selected candidate")
    region_dir = Path(candidate["region_dir"]).expanduser().resolve()
    reads_path = region_dir / "reads/reads.tsv"
    distance_path = region_dir / "distance/BERNOULLI/matrix.csv"
    methylation_path = region_dir / "methylation/methylation.csv"

    verify_assignment_artifacts_and_site_identity(
        assignment,
        candidate,
        region_dir,
        reads_path,
        distance_path,
        methylation_path,
    )
    reads_order, reads = load_reads_strict(reads_path)
    distance_ids, _, distance = load_csv_matrix_strict(distance_path, square=True)
    methylation_ids, _, methylation = load_csv_matrix_strict(methylation_path, square=False)
    if reads_order != distance_ids or reads_order != methylation_ids:
        raise ConfirmationContractError(
            f"Primary reads/BERNOULLI/methylation read identity mismatch at {(sample, chrom, pos)}"
        )
    wrong_chrom = [read_id for read_id in reads_order if reads[read_id]["chr"] != chrom]
    if wrong_chrom:
        raise ConfirmationContractError(
            f"Primary read chromosome mismatch at {(sample, chrom, pos)}: {wrong_chrom[:5]}"
        )

    alt_ids = [
        read_id
        for read_id in distance_ids
        if F.is_tumor(reads[read_id]["is_tumor"])
        and reads[read_id]["alt_support"] == "ALT"
    ]
    alt_id_set = set(alt_ids)
    alt_rows = [reads[read_id] for read_id in reads_order if read_id in alt_id_set]
    candidate_digest = str(candidate.get("alt_readset_sha256", "")).strip()
    if candidate_digest != readset_digest(alt_rows):
        raise ConfirmationContractError(
            f"Candidate alt_readset_sha256 mismatch at {(sample, chrom, pos)}"
        )

    distance_index = {read_id: index for index, read_id in enumerate(distance_ids)}
    alt_indices = [distance_index[read_id] for read_id in alt_ids]
    sub_distance = distance[np.ix_(alt_indices, alt_indices)]
    kept = F.peel_complete(sub_distance)
    kept_ids = [alt_ids[index] for index in kept]
    sub_distance = sub_distance[np.ix_(kept, kept)]
    methylation_index = {read_id: index for index, read_id in enumerate(methylation_ids)}
    sub_methylation = methylation[[methylation_index[read_id] for read_id in kept_ids]]

    if candidate.get("n_alt_after_peel") not in {None, ""}:
        try:
            expected_n = int(candidate["n_alt_after_peel"])
        except ValueError as error:
            raise ConfirmationContractError("Invalid candidate n_alt_after_peel") from error
        if expected_n != len(kept_ids):
            raise ConfirmationContractError(
                f"Candidate n_alt_after_peel mismatch at {(sample, chrom, pos)}: "
                f"{expected_n} != {len(kept_ids)}"
            )

    screen_labels = validate_assignment_identity(
        assignment, candidate, region_dir, reads, kept_ids
    )
    return PreparedSite(
        input_row={key: value for key, value in candidate.items() if not key.startswith("__")},
        source_row_number=int(candidate["__source_row_number"]),
        sample=sample,
        chrom=chrom,
        pos=pos,
        region_dir=region_dir,
        distance=sub_distance,
        methylation=sub_methylation,
        kept_ids=kept_ids,
        screen_labels=screen_labels,
        n_reads_total=len(reads_order),
        n_alt_raw=len(alt_ids),
        n_alt_after_peel=len(kept_ids),
    )


def empirical_p_value(observed: float, null_values: Sequence[float]) -> float | None:
    """Return the plus-one corrected upper-tail empirical p-value."""
    finite = [float(value) for value in null_values if math.isfinite(float(value))]
    if not finite or not math.isfinite(float(observed)):
        return None
    exceed = sum(value >= float(observed) for value in finite)
    return (exceed + 1) / (len(finite) + 1)


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    """BH-adjust finite p-values while preserving their original order."""
    values = [float(value) for value in p_values]
    if any(not math.isfinite(value) or value < 0.0 or value > 1.0 for value in values):
        raise ValueError("BH p-values must be finite and in [0, 1]")
    if not values:
        return []
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    adjusted = [1.0] * len(values)
    running = 1.0
    total = len(values)
    for rank_index in range(total - 1, -1, -1):
        original = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[original] * total / rank)
        adjusted[original] = min(1.0, running)
    return adjusted


def benjamini_yekutieli(p_values: Sequence[float]) -> list[float]:
    """BY-adjust finite p-values while preserving their original order."""
    bh = benjamini_hochberg(p_values)
    harmonic = sum(1.0 / rank for rank in range(1, len(bh) + 1))
    return [min(1.0, value * harmonic) for value in bh]


def partition_signature(labels: Sequence[Any]) -> tuple[tuple[int, ...], ...]:
    """Canonicalize an unlabeled partition as sorted index sets."""
    groups: dict[str, list[int]] = {}
    for index, label in enumerate(labels):
        groups.setdefault(str(label), []).append(index)
    return tuple(sorted((tuple(indices) for indices in groups.values()), key=lambda value: value))


def same_partition(first: Sequence[Any], second: Sequence[Any]) -> bool:
    return len(first) == len(second) and partition_signature(first) == partition_signature(second)


def fixed_k_labels(distance: np.ndarray, groups: int) -> list[int] | None:
    if groups < 2 or distance.shape[0] <= groups:
        return None
    try:
        tree, _ = F.link_tree(distance)
    except (ValueError, FloatingPointError):
        return None
    labels = fcluster(tree, groups, criterion="maxclust").astype(int).tolist()
    return labels if len(set(labels)) == groups else None


def fixed_k_loo_partition_stability(
    distance: np.ndarray,
    reference_labels: list[str],
    modal_groups: int,
    ari_threshold: float,
) -> dict[str, Any]:
    """Audit partition stability after one-read omission while keeping K fixed.

    This diagnostic never reruns candidate discovery, chooses a new K, or
    reevaluates the upstream molecular-haplotype gate.
    """
    n_reads = len(reference_labels)
    if n_reads != distance.shape[0] or n_reads < 2 * F.MIN_SIZE + 1 or modal_groups < 2:
        return {
            "feasible": False,
            "evaluated": 0,
            "invalid": n_reads,
            "ari_median": None,
            "ari_min": None,
            "all_group_size_gate": False,
            "stability_gate": False,
            "exact_partition_count": 0,
            "exact_partition_gate": False,
        }

    aris: list[float] = []
    all_group_sizes_valid = True
    exact_partition_count = 0
    invalid = 0
    for omitted in range(n_reads):
        retained = [index for index in range(n_reads) if index != omitted]
        expected = [reference_labels[index] for index in retained]
        expected_sizes = Counter(expected)
        if len(expected_sizes) != modal_groups or min(expected_sizes.values()) < F.MIN_SIZE:
            invalid += 1
            all_group_sizes_valid = False
            continue
        labels = fixed_k_labels(distance[np.ix_(retained, retained)], modal_groups)
        if labels is None:
            invalid += 1
            all_group_sizes_valid = False
            continue
        sizes = Counter(labels)
        if min(sizes.values()) < F.MIN_SIZE:
            all_group_sizes_valid = False
        aris.append(float(adjusted_rand_score(expected, labels)))
        if same_partition(expected, labels):
            exact_partition_count += 1
    feasible = invalid == 0 and len(aris) == n_reads
    ari_median = float(np.median(aris)) if aris else None
    ari_min = min(aris) if aris else None
    gate = bool(
        feasible
        and all_group_sizes_valid
        and ari_min is not None
        and ari_min >= ari_threshold
    )
    return {
        "feasible": feasible,
        "evaluated": len(aris),
        "invalid": invalid,
        "ari_median": ari_median,
        "ari_min": ari_min,
        "all_group_size_gate": all_group_sizes_valid and feasible,
        "stability_gate": gate,
        "exact_partition_count": exact_partition_count,
        "exact_partition_gate": bool(feasible and exact_partition_count == n_reads),
    }


def leave_one_out_stability(
    distance: np.ndarray,
    reference_labels: list[str],
    modal_groups: int,
    ari_threshold: float,
) -> dict[str, Any]:
    """Deprecated compatibility alias for fixed-K LOO partition stability."""
    return fixed_k_loo_partition_stability(
        distance, reference_labels, modal_groups, ari_threshold
    )


def validate_trace_empirical_p(empirical_p: Any, n_valid_null: int) -> None:
    if empirical_p is None:
        return
    value = float(empirical_p)
    scaled = value * (n_valid_null + 1)
    if value <= 0.0 or value > 1.0 or not math.isclose(scaled, round(scaled), abs_tol=1e-8):
        raise ConfirmationContractError(
            f"Invalid empirical p-value/grid from null trace: p={value} n={n_valid_null}"
        )


def run_null_mode(
    site: PreparedSite, mode: str, config: ConfirmationConfig
) -> dict[str, Any]:
    if mode not in NULL_MODES:
        raise ValueError(f"Unsupported null mode: {mode}")
    base_seed = F.stable_seed(
        site.sample,
        site.chrom,
        site.pos,
        offset=config.seed + MODE_SEED_OFFSETS[mode],
    )
    inference = F.analyze_phylo(
        site.distance,
        site.methylation,
        base_seed=base_seed,
        seeds=1,
        rnull=config.permutations,
        null_mode=mode,
        empirical_alpha=None,
        min_valid_null_fraction=config.min_valid_null_fraction,
    )
    stability = F.analyze_phylo(
        site.distance,
        site.methylation,
        base_seed=base_seed,
        seeds=config.seeds,
        rnull=config.permutations,
        null_mode=mode,
        empirical_alpha=None,
        min_valid_null_fraction=config.min_valid_null_fraction,
    )
    labels = [str(value) for value in inference["coarse_labels"]]
    stability_labels = [str(value) for value in stability["coarse_labels"]]
    if len(labels) != len(site.kept_ids):
        raise ConfirmationContractError(
            f"Strict {mode} assignment length mismatch at {(site.sample, site.chrom, site.pos)}"
        )
    root = inference["coarse_split_trace"][0] if inference["coarse_split_trace"] else {}
    empirical_p = root.get("empirical_p")
    n_valid_null = int(root.get("n_valid_null") or 0)
    validate_trace_empirical_p(empirical_p, n_valid_null)
    minimum_valid = int(math.ceil(config.permutations * config.min_valid_null_fraction))
    observed = root.get("observed_between_within")
    group_sizes = Counter(label for label in labels if label not in {"other", "outlier"})
    minimum_group_size = min(group_sizes.values()) if group_sizes else 0
    screen_ari = float(adjusted_rand_score(site.screen_labels, labels))
    pair_count = int(stability["modal_assignment_pair_count"])
    expected_pair_count = config.seeds * (config.seeds - 1) // 2
    multiseed_ari_min = float(stability["modal_assignment_ari_min"])
    multiseed_exact_partition_gate = bool(
        float(stability["modal_fraction"]) == 1.0
        and not bool(stability["unstable"])
        and pair_count == expected_pair_count
        and multiseed_ari_min == 1.0
        and same_partition(labels, stability_labels)
    )
    screen_exact = same_partition(site.screen_labels, labels)
    fixed_k_loo = leave_one_out_stability(
        site.distance,
        labels,
        int(inference["coarse_ng"]),
        config.assignment_ari_threshold,
    )
    return {
        "null_mode": mode,
        "base_seed": base_seed,
        "observed_between_within": observed,
        "null_threshold": root.get("null_threshold"),
        "empirical_p_postselection_descriptive": empirical_p,
        "postselection_bh_q_descriptive": None,
        "postselection_by_q_descriptive": None,
        "postselection_fdr_calibrated": False,
        "n_valid_null": n_valid_null,
        "minimum_valid_null": minimum_valid,
        "null_valid_gate": bool(empirical_p is not None and n_valid_null >= minimum_valid),
        "effect_gate": bool(observed is not None and float(observed) >= F.SEP_MIN),
        "null_threshold_gate": bool(root.get("passed", False)),
        "multigroup_gate": bool(inference["coarse_ng"] >= 2 and not inference["unstable"]),
        "modal_groups": int(inference["coarse_ng"]),
        "modal_fraction": float(stability["modal_fraction"]),
        "unstable": bool(stability["unstable"]),
        "group_sizes": compact_json(group_sizes),
        "minimum_group_size": minimum_group_size,
        "group_size_gate": bool(len(group_sizes) >= 2 and minimum_group_size >= F.MIN_SIZE),
        "multiseed_pair_count": pair_count,
        "multiseed_ari_median": float(stability["modal_assignment_ari_median"]),
        "multiseed_ari_min": multiseed_ari_min,
        "multiseed_assignment_gate": bool(
            pair_count >= 1 and multiseed_ari_min >= config.assignment_ari_threshold
        ),
        "multiseed_exact_partition_gate": multiseed_exact_partition_gate,
        "screen_assignment_ari": screen_ari,
        "screen_assignment_concordance_gate": bool(
            screen_ari >= config.assignment_ari_threshold
        ),
        "screen_assignment_exact_partition_gate": screen_exact,
        "fixed_k_loo_partition_stability_feasible": fixed_k_loo["feasible"],
        "fixed_k_loo_partition_stability_evaluated": fixed_k_loo["evaluated"],
        "fixed_k_loo_partition_stability_invalid": fixed_k_loo["invalid"],
        "fixed_k_loo_partition_stability_ari_median": fixed_k_loo["ari_median"],
        "fixed_k_loo_partition_stability_ari_min": fixed_k_loo["ari_min"],
        "fixed_k_loo_partition_stability_all_group_size_gate": fixed_k_loo[
            "all_group_size_gate"
        ],
        "fixed_k_loo_partition_stability_ari_gate": fixed_k_loo["stability_gate"],
        "fixed_k_loo_partition_stability_exact_partition_count": fixed_k_loo[
            "exact_partition_count"
        ],
        "fixed_k_loo_partition_stability_exact_partition_gate": fixed_k_loo[
            "exact_partition_gate"
        ],
        "failure": root.get("failure") or "",
        "labels": labels,
    }


def analyze_selected_site(site: PreparedSite, config: ConfirmationConfig) -> dict[str, Any]:
    modes = {mode: run_null_mode(site, mode, config) for mode in NULL_MODES}
    cross_null_ari = float(
        adjusted_rand_score(modes["column"]["labels"], modes["row_circular"]["labels"])
    )
    concordance_values = [
        cross_null_ari,
        *[float(modes[mode]["screen_assignment_ari"]) for mode in NULL_MODES],
        *[float(modes[mode]["multiseed_ari_min"]) for mode in NULL_MODES],
    ]
    return {
        "site": site,
        "modes": modes,
        "strict_cross_null_assignment_ari": cross_null_ari,
        "strict_cross_null_exact_partition_gate": same_partition(
            modes["column"]["labels"], modes["row_circular"]["labels"]
        ),
        "strict_assignment_concordance_ari_min": min(concordance_values),
    }


def finite_p_or_one(value: Any) -> float:
    if value is None:
        return 1.0
    parsed = float(value)
    return parsed if math.isfinite(parsed) and 0.0 <= parsed <= 1.0 else 1.0


def apply_postselection_diagnostics(results: list[dict[str, Any]]) -> None:
    for mode in NULL_MODES:
        values = [
            finite_p_or_one(
                result["modes"][mode]["empirical_p_postselection_descriptive"]
            )
            for result in results
        ]
        bh_values = benjamini_hochberg(values)
        by_values = benjamini_yekutieli(values)
        for result, bh_value, by_value in zip(results, bh_values, by_values):
            result["modes"][mode]["postselection_bh_q_descriptive"] = bh_value
            result["modes"][mode]["postselection_by_q_descriptive"] = by_value

    combined_values: list[float] = []
    for result in results:
        p_values = [
            result["modes"][mode]["empirical_p_postselection_descriptive"]
            for mode in NULL_MODES
        ]
        if any(value is None for value in p_values):
            result["strict_combined_empirical_p_postselection_descriptive"] = None
            combined_values.append(1.0)
        else:
            combined = max(float(value) for value in p_values)
            result["strict_combined_empirical_p_postselection_descriptive"] = combined
            combined_values.append(combined)
    combined_bh = benjamini_hochberg(combined_values)
    combined_by = benjamini_yekutieli(combined_values)
    for result, bh_value, by_value in zip(results, combined_bh, combined_by):
        result["strict_postselection_bh_q_descriptive"] = bh_value
        result["strict_postselection_by_q_descriptive"] = by_value


def flatten_result(
    result: dict[str, Any],
    selection: SelectionColumnContract,
    cooccurrence_contract: dict[str, Any],
    candidate_set_size: int,
    config: ConfirmationConfig,
) -> dict[str, Any]:
    site: PreparedSite = result["site"]
    analysis_class = scientific_analysis_class(config, selection)
    output: dict[str, Any] = dict(site.input_row)
    output.update(
        {
            "strict_source_candidate_row": site.source_row_number,
            "strict_analysis_class": analysis_class,
            "strict_robustness_target": "R1_STRICT_METHYL_PARTITION_ROBUSTNESS_ONLY",
            "strict_formal_parameter_gate": formal_parameter_gate(config),
            "strict_formal_selection_contract_gate": not selection.legacy_fallback_used,
            "strict_formal_selection_column": selection.formal_column,
            "strict_selection_column": selection.formal_column,
            "strict_selection_source_column": selection.resolved_column,
            "strict_selection_legacy_fallback_used": selection.legacy_fallback_used,
            "strict_candidate_selection_gate": True,
            "strict_cooccurrence_receipt_gate": cooccurrence_contract["full_scope_gate"],
            "strict_artifact_identity_gate": True,
            "strict_screen_artifact_hash_lock_gate": True,
            "strict_candidate_table_key_count": cooccurrence_contract[
                "candidate_table_unique_key_count"
            ],
            "strict_postselection_family_size": candidate_set_size,
            "strict_postselection_scope": POSTSELECTION_SCOPE,
            "strict_postselection_fdr_calibrated": False,
            "strict_n_reads_total": site.n_reads_total,
            "strict_n_alt_raw": site.n_alt_raw,
            "strict_n_alt_after_peel": site.n_alt_after_peel,
            "strict_minimum_reads_gate": site.n_alt_after_peel >= 2 * F.MIN_SIZE,
            "strict_combined_empirical_p_postselection_descriptive": result[
                "strict_combined_empirical_p_postselection_descriptive"
            ],
            "strict_postselection_bh_q_descriptive": result[
                "strict_postselection_bh_q_descriptive"
            ],
            "strict_postselection_by_q_descriptive": result[
                "strict_postselection_by_q_descriptive"
            ],
            "strict_cross_null_assignment_ari": result[
                "strict_cross_null_assignment_ari"
            ],
            "strict_cross_null_exact_partition_gate": result[
                "strict_cross_null_exact_partition_gate"
            ],
            "strict_assignment_concordance_ari_min": result[
                "strict_assignment_concordance_ari_min"
            ],
        }
    )
    for mode in NULL_MODES:
        metrics = result["modes"][mode]
        for suffix in MODE_OUTPUT_SUFFIXES:
            output[f"strict_{mode}_{suffix}"] = metrics[suffix]
    failed = [gate for gate in STRICT_GATE_FIELDS if not bool(output[gate])]
    if not output["strict_formal_parameter_gate"]:
        failed.insert(0, "strict_formal_parameter_gate")
    if not output["strict_formal_selection_contract_gate"]:
        failed.insert(0, "strict_formal_selection_contract_gate")
    not_evaluable: list[str] = []
    if not bool(output["strict_minimum_reads_gate"]):
        not_evaluable.append("strict_minimum_reads_gate")
    for mode in NULL_MODES:
        observed = output[f"strict_{mode}_observed_between_within"]
        observed_finite = observed is not None and math.isfinite(float(observed))
        if not observed_finite:
            not_evaluable.append(f"strict_{mode}_observed_effect_undefined")
        elif float(observed) >= F.SEP_MIN and not bool(
            output[f"strict_{mode}_null_valid_gate"]
        ):
            not_evaluable.append(f"strict_{mode}_null_valid_gate")
    evaluable = bool(analysis_class == "FORMAL" and not not_evaluable)
    output["strict_methyl_partition_robustness_evaluable"] = evaluable
    output["strict_null_robustness_pass"] = bool(evaluable and not failed)
    if analysis_class == "EXPLORATORY_ONLY":
        output["strict_methyl_partition_robustness_status"] = "EXPLORATORY_ONLY"
        output["strict_scientific_status"] = "EXPLORATORY_ONLY"
        output["strict_not_evaluable_reason"] = ""
    elif not evaluable:
        output["strict_methyl_partition_robustness_status"] = "NOT_EVALUABLE"
        output["strict_scientific_status"] = "NOT_EVALUABLE"
        output["strict_not_evaluable_reason"] = ";".join(not_evaluable)
    elif output["strict_null_robustness_pass"]:
        output["strict_methyl_partition_robustness_status"] = (
            "R1_STRICT_METHYL_PARTITION_ROBUSTNESS_PASS"
        )
        output["strict_scientific_status"] = "ROBUSTNESS_PASS_NOT_FDR_CALIBRATED"
        output["strict_not_evaluable_reason"] = ""
    else:
        output["strict_methyl_partition_robustness_status"] = (
            "R1_STRICT_METHYL_PARTITION_ROBUSTNESS_FAIL"
        )
        output["strict_scientific_status"] = "ROBUSTNESS_FAIL"
        output["strict_not_evaluable_reason"] = ""
    output["strict_failure_reason"] = ";".join(failed) if failed else "PASS"
    return output


def strict_analysis_replay_digest(rows: Sequence[Mapping[str, Any]]) -> str:
    fields = ["sample", "chrom", "pos", "ref", "alt", *STRICT_OUTPUT_FIELDS]
    canonical_rows = [
        ["" if row.get(field) is None else str(row.get(field)) for field in fields]
        for row in rows
    ]
    canonical_rows.sort(key=lambda row: (row[0], row[1], int(row[2]), row[3], row[4]))
    payload = json.dumps(
        {"fields": fields, "rows": canonical_rows},
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def strict_analysis_replay_record(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    return {
        "contract": "all_strict_output_fields_canonical_string_rows_sha256_v1",
        "fields": ["sample", "chrom", "pos", "ref", "alt", *STRICT_OUTPUT_FIELDS],
        "n_rows": len(rows),
        "sha256": strict_analysis_replay_digest(rows),
    }


def write_deterministic_tsv_gz(
    path: Path, rows: list[dict[str, Any]], fields: list[str]
) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(fileobj=raw_handle, mode="wb", filename="", mtime=0) as gzip_handle:
            with io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="") as text_handle:
                writer = csv.DictWriter(
                    text_handle,
                    fieldnames=fields,
                    delimiter="\t",
                    extrasaction="raise",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)


def summarize_rows(
    rows: list[dict[str, Any]],
    *,
    n_input_rows: int,
    n_assignment_records: int,
    selection: SelectionColumnContract,
    cooccurrence_contract: dict[str, Any],
    config: ConfirmationConfig,
) -> dict[str, Any]:
    gate_counts = {
        gate: sum(bool(row[gate]) for row in rows) for gate in STRICT_GATE_FIELDS
    }
    failed_gate_counts = {gate: len(rows) - count for gate, count in gate_counts.items()}
    by_sample: dict[str, dict[str, int]] = {}
    for sample in sorted({str(row["sample"]) for row in rows}):
        sample_rows = [row for row in rows if row["sample"] == sample]
        sample_evaluable = sum(
            bool(row["strict_methyl_partition_robustness_evaluable"])
            for row in sample_rows
        )
        sample_pass = sum(
            bool(row["strict_null_robustness_pass"]) for row in sample_rows
        )
        by_sample[sample] = {
            "n_selected_candidates": len(sample_rows),
            "n_methyl_partition_robustness_evaluable": sample_evaluable,
            "n_methyl_partition_robustness_not_evaluable": sum(
                row["strict_methyl_partition_robustness_status"] == "NOT_EVALUABLE"
                for row in sample_rows
            ),
            "n_exploratory_only": sum(
                row["strict_analysis_class"] == "EXPLORATORY_ONLY"
                for row in sample_rows
            ),
            "n_null_robustness_pass": sample_pass,
            "n_null_robustness_fail": sample_evaluable - sample_pass,
        }
    n_robust = sum(bool(row["strict_null_robustness_pass"]) for row in rows)
    n_evaluable = sum(
        bool(row["strict_methyl_partition_robustness_evaluable"]) for row in rows
    )
    n_not_evaluable = sum(
        row["strict_methyl_partition_robustness_status"] == "NOT_EVALUABLE"
        for row in rows
    )
    n_exploratory = sum(
        row["strict_analysis_class"] == "EXPLORATORY_ONLY" for row in rows
    )
    analysis_class = scientific_analysis_class(config, selection)
    return {
        "schema_name": "intersubmod.strict_methyl_candidate_confirmation_summary",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": now_utc(),
        "status": "EXECUTION_PASS",
        "execution_status": "EXECUTION_PASS",
        "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
        "scientific_output_class": analysis_class,
        "robustness_target": "R1_strict_methyl_partition_robustness_only",
        "selection_contract": {
            "formal_selection_column": selection.formal_column,
            "selection_column": selection.formal_column,
            "source_selection_column": selection.resolved_column,
            "legacy_alias_consumed": selection.legacy_fallback_used,
            "formal_selection_contract_gate": not selection.legacy_fallback_used,
            "legacy_alias_is_formal_claim": False,
            "immutable_upstream_gate": True,
            "false_rows_processed": False,
            "upstream_thresholds_recomputed_or_relaxed": False,
            "n_input_rows": n_input_rows,
            "n_selected_candidates": len(rows),
            "n_not_selected": n_input_rows - len(rows),
            "n_assignment_records_available": n_assignment_records,
        },
        "cooccurrence_receipt_contract": cooccurrence_contract,
        "denominator_contract": {
            "formal_family": "canonical G2 multi-marker molecular-haplotype base candidates",
            "robustness_denominator": "formal rows with valid minimum-read and both-null evaluations",
            "not_evaluable_is_robustness_failure": False,
            "exploratory_rows_enter_formal_denominator": False,
        },
        "postselection_diagnostic_contract": {
            "scope": POSTSELECTION_SCOPE,
            "family_size": len(rows),
            "genome_wide_site_reference": GENOME_WIDE_SITE_REFERENCE,
            "fdr_calibrated": False,
            "bh_by_values_are_descriptive_only": True,
            "reason": "candidate family was selected using methylation and genetic association before this audit",
            "combined_site_p_descriptive": "max(column empirical p, row-circular empirical p)",
            "missing_or_non_testable_p_in_diagnostics": 1.0,
        },
        "parameters": {
            "permutations_per_seed_per_null": config.permutations,
            "seeds": config.seeds,
            "formal_minimum_permutations": FORMAL_MIN_PERMUTATIONS,
            "formal_minimum_seeds": FORMAL_MIN_SEEDS,
            "formal_parameter_gate": formal_parameter_gate(config),
            "fixed_seed": config.seed,
            "minimum_valid_null_fraction": config.min_valid_null_fraction,
            "minimum_group_size": F.MIN_SIZE,
            "minimum_between_within_ratio": F.SEP_MIN,
            "assignment_ari_threshold": config.assignment_ari_threshold,
            "null_modes": {
                "column": "independent within-CpG column permutation among observed reads",
                "row_circular": "independent circular CpG shift within each read, including identity",
            },
            "inference_seed_contract": "one pre-fixed seed; no modal-seed representative selection",
            "multiseed_contract": "separate stability audit requiring exact partition across every seed",
            "fixed_k_loo_partition_stability": (
                "omit every read once while keeping the already observed K fixed; average-linkage "
                "reclustering on the primary BERNOULLI submatrix; require original groups remain >= "
                "minimum size and exact partition; this does not rerun discovery or select K"
            ),
        },
        "counts": {
            "n_selected_candidates": len(rows),
            "n_methyl_partition_robustness_evaluable": n_evaluable,
            "n_methyl_partition_robustness_not_evaluable": n_not_evaluable,
            "n_null_robustness_pass": n_robust,
            "n_null_robustness_fail_retained": n_evaluable - n_robust,
            "n_exploratory_only": n_exploratory,
        },
        "gate_pass_counts": gate_counts,
        "gate_fail_counts": failed_gate_counts,
        "failure_reason_counts": dict(Counter(row["strict_failure_reason"] for row in rows)),
        "per_sample": by_sample,
        "guardrail": GUARDRAIL,
        "pass": True,
    }


def assert_output_available(path: Path) -> None:
    if os.path.lexists(path):
        raise FileExistsError(f"Refusing to overwrite existing output path: {path}")


def write_not_applicable_receipt(
    output_dir: Path,
    *,
    candidate_table: Path,
    assignments_path: Path,
    cooccurrence_receipt_path: Path,
    cooccurrence_release_receipt_path: Path,
    n_input_rows: int,
    n_assignment_records: int,
    selection: SelectionColumnContract,
    cooccurrence_contract: dict[str, Any],
    config: ConfirmationConfig,
    started_at_utc: str,
    command: Sequence[str],
    source_authority: Mapping[str, Any],
    source_identity_before: Mapping[str, Any],
    source_identity_after_compute: Mapping[str, Any],
    source_modes_before: Mapping[str, str],
    source_modes_after_compute: Mapping[str, str],
    runtime_environment_before: Mapping[str, Any],
    runtime_environment_after_compute: Mapping[str, Any],
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=False)
    receipt_path = output_dir / "not_applicable_receipt.json"
    finished_at_utc = now_utc()
    payload = {
        "schema_name": "intersubmod.strict_methyl_candidate_confirmation.not_applicable",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": finished_at_utc,
        "started_at_utc": started_at_utc,
        "finished_at_utc": finished_at_utc,
        "status": "NOT_APPLICABLE",
        "execution_status": "NOT_APPLICABLE",
        "scientific_output_class": scientific_analysis_class(config, selection),
        "robustness_target": "R1_strict_methyl_partition_robustness_only",
        "reason": "ZERO_SELECTED_CANDIDATES",
        "command": list(command),
        "source_authority": dict(source_authority),
        "selection_column": selection.formal_column,
        "selection_contract": {
            "formal_selection_column": selection.formal_column,
            "selection_column": selection.formal_column,
            "source_selection_column": selection.resolved_column,
            "legacy_alias_consumed": selection.legacy_fallback_used,
            "formal_selection_contract_gate": not selection.legacy_fallback_used,
            "legacy_alias_is_formal_claim": False,
            "immutable_upstream_gate": True,
            "n_input_rows": n_input_rows,
            "n_selected_candidates": 0,
            "n_not_selected": n_input_rows,
            "n_assignment_records_available": n_assignment_records,
        },
        "counts": {
            "n_selected_candidates": 0,
            "n_methyl_partition_robustness_evaluable": 0,
            "n_methyl_partition_robustness_not_evaluable": 0,
            "n_null_robustness_pass": 0,
            "n_null_robustness_fail_retained": 0,
            "n_exploratory_only": 0,
        },
        "inputs": {
            "candidate_table": artifact(candidate_table),
            "assignments": artifact(assignments_path),
            "cooccurrence_receipt": artifact(cooccurrence_receipt_path),
            "cooccurrence_release_receipt": artifact(
                cooccurrence_release_receipt_path
            ),
        },
        "cooccurrence_receipt_contract": cooccurrence_contract,
        "denominator_contract": {
            "not_evaluable_is_robustness_failure": False,
            "exploratory_rows_enter_formal_denominator": False,
        },
        "parameters": {
            "permutations_per_seed_per_null": config.permutations,
            "seeds": config.seeds,
            "formal_minimum_permutations": FORMAL_MIN_PERMUTATIONS,
            "formal_minimum_seeds": FORMAL_MIN_SEEDS,
            "formal_parameter_gate": formal_parameter_gate(config),
        },
        "code": dict(source_identity_after_compute),
        "source_lock": {
            "source_identity_before": dict(source_identity_before),
            "source_identity_after_compute": dict(source_identity_after_compute),
            "source_modes_before": dict(source_modes_before),
            "source_modes_after_compute": dict(source_modes_after_compute),
            "all_sources_read_only_and_unchanged": True,
        },
        "runtime_environment_lock": {
            "identity_before": dict(runtime_environment_before),
            "identity_after_compute": dict(runtime_environment_after_compute),
            "isolated_runtime_unchanged": True,
        },
        "analysis_replay": strict_analysis_replay_record([]),
        "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
        "is_negative_result": False,
        "scientific_interpretation": {
            "is_negative_result": False,
            "statement": (
                "ZERO_SELECTED_CANDIDATES means R1 was not applicable; it is not "
                "evidence that methyl substructure or subclones are absent."
            ),
        },
        "guardrail": GUARDRAIL,
        "pass": True,
    }
    receipt_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    receipt_path.chmod(0o444)
    return receipt_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run an R1 strict methyl-partition robustness audit only for sites selected by a "
            "frozen upstream G2 multi-marker molecular-haplotype base-candidate gate."
        )
    )
    parser.add_argument("--candidate-table", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--cooccurrence-receipt", type=Path, required=True)
    parser.add_argument(
        "--cooccurrence-release-receipt", type=Path, required=True
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    parser.add_argument("--seeds", type=int, default=DEFAULT_SEEDS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--assignment-ari-threshold", type=float, default=DEFAULT_ASSIGNMENT_ARI
    )
    parser.add_argument(
        "--min-valid-null-fraction", type=float, default=DEFAULT_MIN_VALID_NULL_FRACTION
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    started_at_utc = now_utc()
    args = build_parser().parse_args(argv)
    output_dir = args.output_dir.expanduser().resolve()
    assert_output_available(output_dir)
    config = ConfirmationConfig(
        permutations=args.permutations,
        seeds=args.seeds,
        seed=args.seed,
        assignment_ari_threshold=args.assignment_ari_threshold,
        min_valid_null_fraction=args.min_valid_null_fraction,
    )
    validate_config(config)
    formal_release = formal_parameter_gate(config)
    actual_command = [
        *canonical_python_prefix(output_dir),
        str(Path(__file__).resolve()),
        *(list(argv) if argv is not None else sys.argv[1:]),
    ]
    expected_command = canonical_command(
        candidate_table=args.candidate_table,
        assignments=args.assignments,
        cooccurrence_receipt=args.cooccurrence_receipt,
        cooccurrence_release_receipt=args.cooccurrence_release_receipt,
        output_dir=output_dir,
        config=config,
    )
    if actual_command != expected_command:
        raise ConfirmationContractError("Strict producer command is not canonical")
    if formal_release:
        if argv is not None or observed_process_command() != expected_command:
            raise ConfirmationContractError(
                "Formal strict producer is direct-CLI canonical-process only"
            )
    runtime_environment_before = capture_runtime_environment_identity()
    require_canonical_runtime(
        runtime_environment_before,
        output_dir=output_dir,
        formal_cache_required=formal_release,
    )
    source_authority = SOURCE_AUTHORITY.validate_release_source_authority(
        {
            "release_source_authority_validator",
            "strict_producer",
            "focal_alt_cluster_lib",
        }
    )
    source_identity_before = capture_source_identity()
    source_modes_before = capture_source_modes()
    if set(source_modes_before.values()) != {"0o444"}:
        raise ConfirmationContractError(
            f"Strict producer sources are not mode 0444: {source_modes_before}"
        )

    input_fields = candidate_table_fields(args.candidate_table)
    selection = resolve_selection_column(input_fields)
    input_fields, selected, n_input_rows = load_candidates(
        args.candidate_table, selection.resolved_column
    )
    candidate_keys = candidate_table_keys(args.candidate_table)
    assignments = load_assignments(args.assignments)
    selected_keys = [site_key(row, source="selected candidate") for row in selected]
    cooccurrence_contract = validate_cooccurrence_receipt(
        args.cooccurrence_receipt,
        release_receipt_path=args.cooccurrence_release_receipt,
        candidate_table=args.candidate_table,
        assignments_path=args.assignments,
        n_candidate_rows=n_input_rows,
        candidate_keys=candidate_keys,
        selected_keys=selected_keys,
        selection=selection,
    )
    missing_assignments = [key for key in selected_keys if key not in assignments]
    if missing_assignments:
        raise ConfirmationContractError(
            f"Selected candidates missing assignment records: {missing_assignments[:8]}"
        )
    if not selected:
        runtime_environment_after_compute = capture_runtime_environment_identity()
        require_canonical_runtime(
            runtime_environment_after_compute,
            output_dir=output_dir,
            formal_cache_required=formal_release,
        )
        if runtime_environment_after_compute != runtime_environment_before:
            raise ConfirmationContractError(
                "Strict producer runtime environment identity changed during execution"
            )
        source_identity_after_compute = capture_source_identity()
        source_modes_after_compute = capture_source_modes()
        if (
            source_identity_after_compute != source_identity_before
            or source_modes_after_compute != source_modes_before
        ):
            raise ConfirmationContractError(
                "Strict producer source identity changed during execution"
            )
        receipt_path = write_not_applicable_receipt(
            output_dir,
            candidate_table=args.candidate_table,
            assignments_path=args.assignments,
            cooccurrence_receipt_path=args.cooccurrence_receipt,
            cooccurrence_release_receipt_path=(
                args.cooccurrence_release_receipt
            ),
            n_input_rows=n_input_rows,
            n_assignment_records=len(assignments),
            selection=selection,
            cooccurrence_contract=cooccurrence_contract,
            config=config,
            started_at_utc=started_at_utc,
            command=actual_command,
            source_authority=source_authority,
            source_identity_before=source_identity_before,
            source_identity_after_compute=source_identity_after_compute,
            source_modes_before=source_modes_before,
            source_modes_after_compute=source_modes_after_compute,
            runtime_environment_before=runtime_environment_before,
            runtime_environment_after_compute=runtime_environment_after_compute,
        )
        print(
            json.dumps(
                {
                    "input_candidate_table": str(args.candidate_table.resolve()),
                    "input_assignments": str(args.assignments.resolve()),
                    "input_cooccurrence_receipt": str(args.cooccurrence_receipt.resolve()),
                    "output_dir": str(output_dir),
                    "not_applicable_receipt": str(receipt_path),
                    "n_selected_candidates": 0,
                    "execution_status": "NOT_APPLICABLE",
                    "scientific_output_class": scientific_analysis_class(config, selection),
                    "reason": "ZERO_SELECTED_CANDIDATES",
                    "pass": True,
                },
                indent=2,
            )
        )
        return 0

    prepared = [
        prepare_site(row, assignments[key]) for row, key in zip(selected, selected_keys)
    ]
    results = [analyze_selected_site(site, config) for site in prepared]
    apply_postselection_diagnostics(results)
    rows = [
        flatten_result(
            result,
            selection,
            cooccurrence_contract,
            len(results),
            config,
        )
        for result in results
    ]
    rows.sort(key=lambda row: (str(row["sample"]), str(row["chrom"]), int(row["pos"])))

    runtime_environment_after_compute = capture_runtime_environment_identity()
    require_canonical_runtime(
        runtime_environment_after_compute,
        output_dir=output_dir,
        formal_cache_required=formal_release,
    )
    if runtime_environment_after_compute != runtime_environment_before:
        raise ConfirmationContractError(
            "Strict producer runtime environment identity changed during execution"
        )
    source_identity_after_compute = capture_source_identity()
    source_modes_after_compute = capture_source_modes()
    if (
        source_identity_after_compute != source_identity_before
        or source_modes_after_compute != source_modes_before
    ):
        raise ConfirmationContractError(
            "Strict producer source identity changed during execution"
        )

    output_dir.mkdir(parents=True, exist_ok=False)
    site_path = output_dir / "strict_methyl_candidate_confirmation_sites.tsv.gz"
    summary_path = output_dir / "strict_methyl_candidate_confirmation_summary.json"
    manifest_path = output_dir / "run_manifest.json"
    output_fields = input_fields + STRICT_OUTPUT_FIELDS
    write_deterministic_tsv_gz(site_path, rows, output_fields)
    summary = summarize_rows(
        rows,
        n_input_rows=n_input_rows,
        n_assignment_records=len(assignments),
        selection=selection,
        cooccurrence_contract=cooccurrence_contract,
        config=config,
    )
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    finished_at_utc = now_utc()
    manifest = {
        "schema_name": "intersubmod.strict_methyl_candidate_confirmation_run_manifest",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": finished_at_utc,
        "started_at_utc": started_at_utc,
        "finished_at_utc": finished_at_utc,
        "status": "EXECUTION_PASS",
        "execution_status": "EXECUTION_PASS",
        "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
        "scientific_output_class": scientific_analysis_class(config, selection),
        "robustness_target": "R1_strict_methyl_partition_robustness_only",
        "command": actual_command,
        "source_authority": source_authority,
        "inputs": {
            "candidate_table": artifact(args.candidate_table),
            "assignments": artifact(args.assignments),
            "cooccurrence_receipt": artifact(args.cooccurrence_receipt),
            "cooccurrence_release_receipt": artifact(
                args.cooccurrence_release_receipt
            ),
        },
        "code": source_identity_after_compute,
        "source_lock": {
            "source_identity_before": source_identity_before,
            "source_identity_after_compute": source_identity_after_compute,
            "source_modes_before": source_modes_before,
            "source_modes_after_compute": source_modes_after_compute,
            "all_sources_read_only_and_unchanged": True,
        },
        "runtime_environment_lock": {
            "identity_before": runtime_environment_before,
            "identity_after_compute": runtime_environment_after_compute,
            "isolated_runtime_unchanged": True,
        },
        "analysis_replay": strict_analysis_replay_record(rows),
        "selection_contract": summary["selection_contract"],
        "cooccurrence_receipt_contract": cooccurrence_contract,
        "denominator_contract": summary["denominator_contract"],
        "postselection_diagnostic_contract": summary["postselection_diagnostic_contract"],
        "parameters": summary["parameters"],
        "contracts": {
            "output_overwrite_allowed": False,
            "primary_artifact_identity_mismatch": "HARD_FAIL_NO_OUTPUT",
            "assignment_read_identity_mismatch": "HARD_FAIL_NO_OUTPUT",
            "screen_positive_strict_fail_rows_retained": True,
            "screen_positive_all_status_rows_retained": True,
            "strict_gate_fields": STRICT_GATE_FIELDS,
            "candidate_selection_relaxed": False,
            "postselection_fdr_claim_allowed": False,
            "legacy_selection_alias_is_formal_claim": False,
            "formal_parameter_gate": formal_parameter_gate(config),
            "formal_selection_contract_gate": not selection.legacy_fallback_used,
            "fixed_k_loo_reruns_discovery": False,
            "scientific_claim_ceiling": (
                "R1_methyl_partition_robustness_only_no_molecular_haplotype_or_clone_upgrade"
            ),
        },
        "outputs": {
            "site_results": artifact(site_path),
            "summary": artifact(summary_path),
        },
        "counts": summary["counts"],
        "guardrail": GUARDRAIL,
        "pass": True,
    }
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    for path in (site_path, summary_path, manifest_path):
        path.chmod(0o444)
    print(
        json.dumps(
            {
                "input_candidate_table": str(args.candidate_table.resolve()),
                "input_assignments": str(args.assignments.resolve()),
                "input_cooccurrence_receipt": str(args.cooccurrence_receipt.resolve()),
                "output_dir": str(output_dir),
                "site_results": str(site_path),
                "summary": str(summary_path),
                "run_manifest": str(manifest_path),
                "n_selected_candidates": len(rows),
                "n_null_robustness_pass": summary["counts"]["n_null_robustness_pass"],
                "scientific_output_class": summary["scientific_output_class"],
                "pass_semantics": summary["pass_semantics"],
                "pass": True,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
