#!/usr/bin/env python3
"""Analyze exact-raw-HP multi-sSNV patterns against read-level methylation.

This module consumes the frozen pattern census and the hash-bound 2026-07-15
all-sSNV artifact catalog.  It rebuilds a union read x CpG matrix from every
marker's +/-5 kb tile, masks marker-proximal CpGs, computes the canonical
BERNOULLI distance, and runs within-(read-group x strand) restricted
PERMANOVA/PERMDISP.  Methylation is an association-only overlay and never
changes topology selection or direction.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import hashlib
import io
import json
import math
import os
import statistics
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np
from statsmodels.stats.multitest import multipletests


SCHEMA_NAME = "intersubmod.pattern_methyl_evidence"
SCHEMA_VERSION = "1.0.0"
DEFAULT_SEED = 20260727
DEFAULT_PERMUTATIONS = 999
DEFAULT_SENSITIVITY_PERMUTATIONS = 199
ADAPTIVE_SELECTION_SCHEMA_NAME = (
    "intersubmod.pattern_methyl_adaptive_refinement_selection"
)
ADAPTIVE_SELECTION_SCHEMA_VERSION = "1.0.0"
ADAPTIVE_SCREEN_PERMUTATIONS = 999
ADAPTIVE_REFINEMENT_PERMUTATIONS = 49999
MIN_COMPLETE_N = 40
MIN_STATE_N = 5
MIN_COMMON_CPG = 10
MIN_PAIR_CPG = 3
STATE_CPG_COVERAGE = 0.80
READ_CPG_COVERAGE = 0.80
MARKER_MASK_BP = 100
DISTAL_BP = 1000

RESULT_FIELDS = (
    "schema_version",
    "dataset",
    "chrom",
    "region_id",
    "unit_id",
    "phase_set",
    "hp_family",
    "hp_raw",
    "active_positions",
    "n_active_bits",
    "pair_full4",
    "k_ge_3",
    "input_n_complete",
    "input_state_counts_json",
    "analysis_n",
    "analysis_state_counts_json",
    "n_common_cpg",
    "n_distal_cpg",
    "qname_join_fraction_min",
    "tile_overlap_conflicts",
    "exchangeable_strata",
    "exchangeable_n",
    "permanova_pseudo_f",
    "permanova_r2",
    "permanova_p",
    "permanova_permutations_requested",
    "permanova_permutations_realized",
    "permdisp_f",
    "permdisp_p",
    "best_pair",
    "best_pair_hamming",
    "best_pair_between_mean",
    "best_pair_pooled_within_mean",
    "best_pair_distance_contrast",
    "best_pair_standardized_effect",
    "best_pair_topology_relation",
    "max_geometry_smd",
    "geometry_feature",
    "all_states_n8",
    "all_states_n10",
    "equal_n_r2",
    "equal_n_retention",
    "rarefaction_median_r2",
    "rarefaction_retention",
    "distal_r2",
    "distal_p",
    "distal_permutations_realized",
    "distal_retention",
    "multiplicity_family",
    "q_by",
    "p_holm",
    "assessment",
    "evaluation_status",
    "invalid_reason",
)


class AnalysisContractError(RuntimeError):
    """Raised when a frozen input or identity contract is violated."""


@dataclass(frozen=True)
class AnalysisConfig:
    permutations: int = DEFAULT_PERMUTATIONS
    sensitivity_permutations: int = DEFAULT_SENSITIVITY_PERMUTATIONS
    seed: int = DEFAULT_SEED
    marker_mask_bp: int = MARKER_MASK_BP
    distal_bp: int = DISTAL_BP
    min_common_cpg: int = MIN_COMMON_CPG
    min_pair_cpg: int = MIN_PAIR_CPG
    state_cpg_coverage: float = STATE_CPG_COVERAGE
    read_cpg_coverage: float = READ_CPG_COVERAGE


_WORKER_CATALOG: dict[tuple[str, str, int], dict[str, str]] = {}
_WORKER_TOPOLOGY: dict[str, dict[str, Any]] = {}
_WORKER_CONFIG = AnalysisConfig()


def bool_value(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def parse_int_list(value: object) -> tuple[int, ...]:
    text = str(value or "").strip()
    if not text:
        return ()
    if text.startswith("["):
        payload = json.loads(text)
        return tuple(int(item) for item in payload)
    return tuple(int(item) for item in text.split(",") if item)


def stable_seed(seed: int, *parts: object) -> int:
    payload = "\x1f".join(str(part) for part in (seed,) + parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big") & 0x7FFFFFFF


def sha256_read_name(read_name: str) -> str:
    return hashlib.sha256(read_name.encode("utf-8")).hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size <= 0:
        raise AnalysisContractError(f"bound file missing or empty: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def assert_file_identity(
    path: Path,
    expected: Mapping[str, Any],
    *,
    label: str,
) -> dict[str, Any]:
    """Re-attest a file against a frozen path/size/SHA-256 identity."""
    observed = file_identity(path)
    normalized_expected = {
        "path": str(Path(str(expected.get("path", ""))).expanduser().resolve()),
        "size_bytes": int(expected.get("size_bytes", -1)),
        "sha256": str(expected.get("sha256", "")),
    }
    if observed != normalized_expected:
        raise AnalysisContractError(
            f"{label} source identity drift: "
            f"observed={observed} expected={normalized_expected}"
        )
    return observed


def catalog_file_identity(
    record: Mapping[str, str], prefix: str
) -> dict[str, Any]:
    return {
        "path": str(Path(record[f"{prefix}_path"]).expanduser().resolve()),
        "size_bytes": int(record[f"{prefix}_size_bytes"]),
        "sha256": record[f"{prefix}_sha256"],
    }


def attest_catalog_record(
    record: Mapping[str, str],
    *,
    marker_key: tuple[str, str, int],
) -> None:
    for prefix in ("reads", "methylation"):
        expected = catalog_file_identity(record, prefix)
        assert_file_identity(
            Path(expected["path"]),
            expected,
            label=f"catalog {marker_key} {prefix}",
        )


def attest_catalog_artifacts(
    catalog: Mapping[tuple[str, str, int], Mapping[str, str]],
) -> None:
    seen: set[tuple[str, int, str]] = set()
    for marker_key, record in sorted(catalog.items()):
        for prefix in ("reads", "methylation"):
            expected = catalog_file_identity(record, prefix)
            identity_key = (
                expected["path"],
                expected["size_bytes"],
                expected["sha256"],
            )
            if identity_key in seen:
                continue
            seen.add(identity_key)
            assert_file_identity(
                Path(expected["path"]),
                expected,
                label=f"catalog {marker_key} {prefix}",
            )


def _resolve_receipt_binding(receipt_path: Path, binding: Mapping[str, Any]) -> Path:
    raw_path = binding.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        raise AnalysisContractError(
            f"adaptive receipt file binding lacks path: {binding!r}"
        )
    bound = Path(raw_path).expanduser()
    if not bound.is_absolute():
        bound = receipt_path.parent / bound
    return bound.resolve()


def _validate_receipt_file_binding(
    receipt_path: Path,
    binding: object,
    actual_path: Path | None,
    *,
    label: str,
) -> dict[str, Any]:
    if not isinstance(binding, dict):
        raise AnalysisContractError(
            f"adaptive receipt {label} binding is not an object"
        )
    bound_path = _resolve_receipt_binding(receipt_path, binding)
    if actual_path is not None and bound_path != actual_path.expanduser().resolve():
        raise AnalysisContractError(
            f"adaptive receipt {label} path does not bind supplied input"
        )
    identity = file_identity(bound_path)
    if binding.get("sha256") != identity["sha256"]:
        raise AnalysisContractError(
            f"adaptive receipt {label} SHA-256 mismatch"
        )
    try:
        bound_size = int(binding.get("size_bytes", -1))
    except (TypeError, ValueError) as exc:
        raise AnalysisContractError(
            f"adaptive receipt {label} size is invalid"
        ) from exc
    if bound_size != identity["size_bytes"]:
        raise AnalysisContractError(f"adaptive receipt {label} size mismatch")
    return identity


def validate_adaptive_refinement_binding(
    unit_key_file: Path,
    unit_key_receipt: Path,
    *,
    analyzer_permutations: int,
    selected_key_count: int,
) -> dict[str, dict[str, Any]]:
    """Validate the frozen selector receipt and return summary source bindings."""
    receipt_identity = file_identity(unit_key_receipt)
    receipt_path = Path(receipt_identity["path"])
    try:
        payload = json.loads(receipt_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise AnalysisContractError(
            f"adaptive unit-key receipt is invalid JSON: {receipt_path}"
        ) from exc
    if not isinstance(payload, dict):
        raise AnalysisContractError("adaptive unit-key receipt root is not an object")
    if payload.get("schema_name") != ADAPTIVE_SELECTION_SCHEMA_NAME:
        raise AnalysisContractError("adaptive unit-key receipt schema_name mismatch")
    if payload.get("schema_version") != ADAPTIVE_SELECTION_SCHEMA_VERSION:
        raise AnalysisContractError("adaptive unit-key receipt schema_version mismatch")
    if payload.get("all_pass") is not True:
        raise AnalysisContractError("adaptive unit-key receipt is not all_pass")

    contract = payload.get("contract")
    if not isinstance(contract, dict):
        raise AnalysisContractError("adaptive unit-key receipt contract missing")
    if contract.get("family") != "CONFIRMATORY_FULL4_OR_LONG":
        raise AnalysisContractError("adaptive unit-key receipt family mismatch")
    try:
        screen_permutations = int(contract.get("screen_permutations", -1))
        refinement_permutations = int(
            contract.get("refinement_permutations", -1)
        )
    except (TypeError, ValueError) as exc:
        raise AnalysisContractError(
            "adaptive unit-key receipt permutation budget is invalid"
        ) from exc
    if screen_permutations != ADAPTIVE_SCREEN_PERMUTATIONS:
        raise AnalysisContractError(
            "adaptive unit-key receipt screen budget mismatch: "
            f"{screen_permutations}"
        )
    if refinement_permutations != ADAPTIVE_REFINEMENT_PERMUTATIONS:
        raise AnalysisContractError(
            "adaptive unit-key receipt refinement budget mismatch: "
            f"{refinement_permutations}"
        )
    try:
        screen_floor = float(contract.get("screen_floor", float("nan")))
    except (TypeError, ValueError) as exc:
        raise AnalysisContractError(
            "adaptive unit-key receipt screen floor is invalid"
        ) from exc
    expected_floor = 1.0 / (ADAPTIVE_SCREEN_PERMUTATIONS + 1)
    if not math.isclose(screen_floor, expected_floor, rel_tol=0.0, abs_tol=1e-15):
        raise AnalysisContractError(
            f"adaptive unit-key receipt screen floor mismatch: {screen_floor}"
        )
    if analyzer_permutations != refinement_permutations:
        raise AnalysisContractError(
            "analyzer permutation budget does not match adaptive receipt: "
            f"{analyzer_permutations}!={refinement_permutations}"
        )

    counts = payload.get("counts")
    if not isinstance(counts, dict):
        raise AnalysisContractError("adaptive unit-key receipt counts missing")
    try:
        receipt_selected = int(counts.get("selected_for_refinement", -1))
    except (TypeError, ValueError) as exc:
        raise AnalysisContractError(
            "adaptive receipt selected count is invalid"
        ) from exc
    if receipt_selected != selected_key_count:
        raise AnalysisContractError(
            "adaptive receipt selected count does not match unit-key TSV: "
            f"{receipt_selected}!={selected_key_count}"
        )

    outputs = payload.get("outputs")
    inputs = payload.get("inputs")
    if not isinstance(outputs, dict) or not isinstance(inputs, dict):
        raise AnalysisContractError("adaptive unit-key receipt bindings missing")
    unit_key_identity = _validate_receipt_file_binding(
        receipt_path,
        outputs.get("unit_keys"),
        unit_key_file,
        label="unit_keys",
    )
    screen_identity = _validate_receipt_file_binding(
        receipt_path,
        inputs.get("screen_evidence"),
        None,
        label="screen_evidence",
    )
    return {
        "unit_key_file": unit_key_identity,
        "unit_key_receipt": receipt_identity,
        "screen_evidence": screen_identity,
    }


def natural_chrom_key(chrom: str) -> tuple[int, str]:
    if chrom.startswith("chr") and chrom[3:].isdigit():
        return int(chrom[3:]), chrom
    return 10**9, chrom


def bernoulli_distance_matrix(matrix: np.ndarray, min_common: int = MIN_PAIR_CPG) -> np.ndarray:
    """Vectorized parity implementation of DistanceMatrix.cpp BERNOULLI."""
    if matrix.ndim != 2:
        raise ValueError("methylation matrix must be two-dimensional")
    finite = np.isfinite(matrix).astype(float)
    values = np.where(np.isfinite(matrix), matrix, 0.0)
    confidence = 2.0 * np.abs(values - 0.5) * finite
    weighted_probability = confidence * values
    weight_sum = confidence @ confidence.T
    weighted_disagreement = (
        weighted_probability @ confidence.T
        + confidence @ weighted_probability.T
        - 2.0 * (weighted_probability @ weighted_probability.T)
    )
    common = finite @ finite.T
    with np.errstate(divide="ignore", invalid="ignore"):
        distance = weighted_disagreement / weight_sum
    distance[(common < min_common) | (weight_sum < 1e-9)] = np.nan
    np.fill_diagonal(distance, 0.0)
    return distance


def compute_ss(distance: np.ndarray, labels: Sequence[str]) -> tuple[float, float, float]:
    n = distance.shape[0]
    squared = distance * distance
    total = float(np.triu(squared, 1).sum()) / float(n)
    within = 0.0
    labels_array = np.asarray(labels, dtype=object)
    for label in sorted(set(labels)):
        indices = np.flatnonzero(labels_array == label)
        within += float(np.triu(squared[np.ix_(indices, indices)], 1).sum()) / float(
            len(indices)
        )
    between = max(0.0, total - within)
    return total, within, between


def pseudo_f_r2(distance: np.ndarray, labels: Sequence[str]) -> tuple[float, float]:
    n = distance.shape[0]
    k = len(set(labels))
    total, within, between = compute_ss(distance, labels)
    if k < 2 or n <= k:
        return 0.0, 0.0
    if within <= 0.0:
        pseudo_f = 1e9 if between > 0.0 else 0.0
    else:
        pseudo_f = (between / (k - 1)) / (within / (n - k))
    r_squared = min(1.0, max(0.0, between / total)) if total > 0.0 else 0.0
    return float(pseudo_f), float(r_squared)


def distances_to_centroid(distance: np.ndarray, labels: Sequence[str]) -> np.ndarray:
    labels_array = np.asarray(labels, dtype=object)
    output = np.zeros(distance.shape[0], dtype=float)
    squared = distance * distance
    for label in sorted(set(labels)):
        indices = np.flatnonzero(labels_array == label)
        size = len(indices)
        if size <= 1:
            continue
        block = squared[np.ix_(indices, indices)]
        correction = float(block.sum()) / (2.0 * size * size)
        values = block.sum(axis=1) / float(size) - correction
        output[indices] = np.sqrt(np.maximum(0.0, values))
    return output


def anova_f(values: np.ndarray, labels: Sequence[str]) -> float:
    labels_array = np.asarray(labels, dtype=object)
    groups = [values[labels_array == label] for label in sorted(set(labels))]
    n = len(values)
    k = len(groups)
    if k < 2 or n <= k:
        return 0.0
    grand = float(values.mean())
    between = sum(len(group) * (float(group.mean()) - grand) ** 2 for group in groups)
    within = sum(float(((group - group.mean()) ** 2).sum()) for group in groups)
    tolerance = 128.0 * np.finfo(float).eps * max(1.0, float((values * values).sum()))
    if within <= tolerance:
        return 1e9 if between > tolerance else 0.0
    return float((between / (k - 1)) / (within / (n - k)))


def exchangeable_indices(labels: Sequence[str], strata: Sequence[str]) -> np.ndarray:
    labels_by_stratum: dict[str, set[str]] = defaultdict(set)
    for label, stratum in zip(labels, strata):
        labels_by_stratum[stratum].add(label)
    exchangeable = {
        stratum for stratum, values in labels_by_stratum.items() if len(values) >= 2
    }
    return np.asarray(
        [index for index, stratum in enumerate(strata) if stratum in exchangeable],
        dtype=int,
    )


def restricted_permanova_permdisp(
    distance: np.ndarray,
    labels: Sequence[str],
    strata: Sequence[str],
    permutations: int,
    seed: int,
) -> dict[str, float | int]:
    if distance.shape[0] != len(labels) or len(labels) != len(strata):
        raise ValueError("distance, labels, and strata lengths differ")
    if len(set(labels)) < 2:
        raise ValueError("insufficient groups")
    stratum_indices: dict[str, np.ndarray] = {}
    strata_array = np.asarray(strata, dtype=object)
    labels_array = np.asarray(labels, dtype=object)
    for stratum in sorted(set(strata)):
        indices = np.flatnonzero(strata_array == stratum)
        if len(set(labels_array[indices].tolist())) < 2:
            raise ValueError("no_exchangeable_labels")
        stratum_indices[stratum] = indices

    observed_f, observed_r2 = pseudo_f_r2(distance, labels)
    observed_dispersion = anova_f(distances_to_centroid(distance, labels), labels)
    permanova_extreme = 1
    permdisp_extreme = 1
    rng = np.random.default_rng(seed)
    for _ in range(permutations):
        permuted = labels_array.copy()
        for indices in stratum_indices.values():
            permuted[indices] = rng.permutation(permuted[indices])
        permuted_list = permuted.tolist()
        permuted_f, _ = pseudo_f_r2(distance, permuted_list)
        if permuted_f >= observed_f:
            permanova_extreme += 1
        permuted_dispersion = anova_f(
            distances_to_centroid(distance, permuted_list), permuted_list
        )
        if permuted_dispersion >= observed_dispersion:
            permdisp_extreme += 1
    denominator = permutations + 1
    return {
        "pseudo_f": observed_f,
        "r_squared": observed_r2,
        "p_value": permanova_extreme / denominator,
        "permdisp_f": observed_dispersion,
        "permdisp_p": permdisp_extreme / denominator,
        "requested": permutations,
        "realized": permutations,
    }


def hamming_distance(first: str, second: str) -> int:
    if len(first) != len(second):
        raise ValueError("patterns have different lengths")
    return sum(left != right for left, right in zip(first, second))


def pattern_vertex(pattern: str) -> int:
    if any(code not in {"R", "A"} for code in pattern):
        raise ValueError("partial pattern cannot map to a topology vertex")
    return sum((1 << index) for index, code in enumerate(pattern) if code == "A")


def topology_pair_relation(
    first: str, second: str, topology: Mapping[str, Any] | None
) -> str:
    hamming = hamming_distance(first, second)
    if hamming > 1:
        return "PAIR_BAND_ONLY_HAMMING_GT1"
    if topology is None:
        return "HAMMING1_TOPOLOGY_UNAVAILABLE"
    if not bool(topology.get("best_tree_unique")):
        return "HAMMING1_NOT_UNANIMOUS"
    vertices = {pattern_vertex(first), pattern_vertex(second)}
    for edge in topology.get("representative_best_edges", []):
        if {int(edge["parent_vertex"]), int(edge["child_vertex"])} == vertices:
            return "HAMMING1_GLOBAL_BEST_UNANIMOUS"
    return "HAMMING1_NOT_IN_GLOBAL_BEST"


def geometry_max_smd(rows: Sequence[Mapping[str, Any]], labels: Sequence[str]) -> tuple[float, str]:
    features: dict[str, np.ndarray] = {
        "mapq": np.asarray([float(row["mapq"]) for row in rows]),
        "read_length": np.asarray(
            [float(row["end0"]) - float(row["start0"]) for row in rows]
        ),
        "start0": np.asarray([float(row["start0"]) for row in rows]),
        "end0": np.asarray([float(row["end0"]) for row in rows]),
    }
    label_array = np.asarray(labels, dtype=object)
    unique = sorted(set(labels))
    best = (0.0, "")
    for feature_name, values in features.items():
        for left_index, left in enumerate(unique):
            first = values[label_array == left]
            for right in unique[left_index + 1 :]:
                second = values[label_array == right]
                denominator_df = len(first) + len(second) - 2
                if denominator_df <= 0:
                    continue
                pooled_variance = (
                    (len(first) - 1) * float(first.var(ddof=1))
                    + (len(second) - 1) * float(second.var(ddof=1))
                ) / denominator_df
                mean_difference = abs(float(first.mean()) - float(second.mean()))
                if pooled_variance <= 0.0:
                    effect = math.inf if mean_difference > 0.0 else 0.0
                else:
                    effect = mean_difference / math.sqrt(pooled_variance)
                if effect > best[0]:
                    best = (float(effect), f"{feature_name}:{left}|{right}")
    return best


def pairwise_effects(
    distance: np.ndarray,
    labels: Sequence[str],
    topology: Mapping[str, Any] | None,
) -> list[dict[str, Any]]:
    label_array = np.asarray(labels, dtype=object)
    output: list[dict[str, Any]] = []
    states = sorted(set(labels))
    for first_index, first in enumerate(states):
        first_indices = np.flatnonzero(label_array == first)
        first_within = distance[np.ix_(first_indices, first_indices)]
        first_values = first_within[np.triu_indices(len(first_indices), 1)]
        for second in states[first_index + 1 :]:
            second_indices = np.flatnonzero(label_array == second)
            second_within = distance[np.ix_(second_indices, second_indices)]
            second_values = second_within[np.triu_indices(len(second_indices), 1)]
            between_values = distance[np.ix_(first_indices, second_indices)].ravel()
            within_values = np.concatenate([first_values, second_values])
            between_mean = float(between_values.mean())
            within_mean = float(within_values.mean()) if len(within_values) else math.nan
            contrast = between_mean - within_mean
            pooled = np.concatenate([between_values, within_values])
            standard_deviation = float(pooled.std(ddof=1)) if len(pooled) > 1 else 0.0
            standardized = contrast / standard_deviation if standard_deviation > 0 else 0.0
            output.append(
                {
                    "first": first,
                    "second": second,
                    "hamming": hamming_distance(first, second),
                    "between_mean": between_mean,
                    "pooled_within_mean": within_mean,
                    "distance_contrast": contrast,
                    "standardized_effect": standardized,
                    "topology_relation": topology_pair_relation(first, second, topology),
                }
            )
    return sorted(
        output,
        key=lambda row: (
            -float(row["distance_contrast"]),
            -float(row["standardized_effect"]),
            row["first"],
            row["second"],
        ),
    )


def load_formal_counts(path: Path) -> dict[tuple[str, str, str, str], dict[str, str]]:
    output: dict[tuple[str, str, str, str], dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "dataset",
            "chrom",
            "region_id",
            "unit_id",
            "phase_set",
            "hp_family",
            "hp_raw",
            "n_total",
            "n_complete",
            "n_partial",
            "formal_n5",
            "state_count_json",
            "partial_state_count_json",
            "active_positions",
            "n_active_bits",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise AnalysisContractError(
                f"pattern counts missing columns: {sorted(required - set(reader.fieldnames or []))}"
            )
        for row in reader:
            if not bool_value(row["formal_n5"]):
                continue
            key = (row["dataset"], row["chrom"], row["region_id"], row["hp_raw"])
            if key in output:
                raise AnalysisContractError(f"duplicate formal count key: {key}")
            try:
                complete_counts = {
                    str(state): int(count)
                    for state, count in json.loads(row["state_count_json"]).items()
                }
                partial_counts = {
                    str(state): int(count)
                    for state, count in json.loads(
                        row["partial_state_count_json"]
                    ).items()
                }
            except (AttributeError, TypeError, ValueError, json.JSONDecodeError) as exc:
                raise AnalysisContractError(
                    f"invalid state-count JSON for formal count key: {key}"
                ) from exc
            n_active_bits = int(row["n_active_bits"])
            if (
                n_active_bits < 2
                or len(parse_int_list(row["active_positions"])) != n_active_bits
                or any(
                    len(state) != n_active_bits
                    or not set(state).issubset({"R", "A"})
                    or count <= 0
                    for state, count in complete_counts.items()
                )
                or any(
                    len(state) != n_active_bits
                    or "X" not in state
                    or not set(state).issubset({"R", "A", "X"})
                    or count <= 0
                    for state, count in partial_counts.items()
                )
                or sum(complete_counts.values()) != int(row["n_complete"])
                or sum(partial_counts.values()) != int(row["n_partial"])
                or int(row["n_total"])
                != int(row["n_complete"]) + int(row["n_partial"])
                or int(row["n_complete"]) < MIN_COMPLETE_N
                or sum(count >= MIN_STATE_N for count in complete_counts.values()) < 2
            ):
                raise AnalysisContractError(
                    f"formal count invariants failed for key: {key}"
                )
            output[key] = dict(row)
    if not output:
        raise AnalysisContractError("pattern counts contain no formal_n5 units")
    return output


def load_unit_keys(path: Path) -> set[tuple[str, str, str, str]]:
    output: set[tuple[str, str, str, str]] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"dataset", "chrom", "region_id", "hp_raw"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise AnalysisContractError(
                f"unit-key file missing columns: "
                f"{sorted(required - set(reader.fieldnames or []))}"
            )
        for line_number, row in enumerate(reader, start=2):
            key = tuple(row[field].strip() for field in (
                "dataset",
                "chrom",
                "region_id",
                "hp_raw",
            ))
            if any(not value for value in key):
                raise AnalysisContractError(
                    f"unit-key file contains an empty key at line {line_number}"
                )
            if key in output:
                raise AnalysisContractError(
                    f"unit-key file contains duplicate key: {key}"
                )
            output.add(key)
    if not output:
        raise AnalysisContractError("unit-key file contains no units")
    return output


def load_catalog(path: Path) -> dict[tuple[str, str, int], dict[str, str]]:
    output: dict[tuple[str, str, int], dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "dataset",
            "chrom",
            "position1",
            "status",
            "reads_path",
            "reads_size_bytes",
            "reads_sha256",
            "methylation_path",
            "methylation_size_bytes",
            "methylation_sha256",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise AnalysisContractError(
                f"artifact catalog missing columns: {sorted(required - set(reader.fieldnames or []))}"
            )
        for row in reader:
            key = (row["dataset"], row["chrom"], int(row["position1"]))
            if key in output:
                raise AnalysisContractError(f"duplicate artifact marker: {key}")
            if row["status"] != "PASS":
                raise AnalysisContractError(f"artifact marker not PASS: {key} {row['status']}")
            for field, size_field, sha_field in (
                ("reads_path", "reads_size_bytes", "reads_sha256"),
                (
                    "methylation_path",
                    "methylation_size_bytes",
                    "methylation_sha256",
                ),
            ):
                artifact = Path(row[field]).resolve()
                if not artifact.is_file() or artifact.stat().st_size <= 0:
                    raise AnalysisContractError(f"catalog path unavailable: {artifact}")
                try:
                    expected_size = int(row[size_field])
                except ValueError as exc:
                    raise AnalysisContractError(
                        f"catalog has invalid {size_field}: {key}"
                    ) from exc
                if artifact.stat().st_size != expected_size:
                    raise AnalysisContractError(
                        f"catalog artifact size mismatch: {artifact}"
                    )
                observed_sha = sha256_file(artifact)
                if observed_sha != row[sha_field]:
                    raise AnalysisContractError(
                        f"catalog artifact SHA-256 mismatch: {artifact}"
                    )
            output[key] = dict(row)
    if not output:
        raise AnalysisContractError("artifact catalog is empty")
    return output


def load_topology(
    topology_root: Path, datasets: Iterable[str]
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    output: dict[str, dict[str, Any]] = {}
    sources: list[dict[str, Any]] = []
    for dataset in sorted(set(datasets)):
        path = topology_root / "samples" / dataset / f"{dataset}.topology.jsonl"
        if not path.is_file():
            raise AnalysisContractError(f"topology JSONL missing: {path}")
        before = {
            "dataset": dataset,
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                try:
                    row = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise AnalysisContractError(
                        f"invalid topology JSON {path}:{line_number}: {exc}"
                    ) from exc
                if (
                    row.get("schema_name")
                    != "intersubmod.exact_ps_cpp_topology_af.unit"
                    or row.get("schema_version") != "1.0.0"
                    or row.get("sample") != dataset
                ):
                    raise AnalysisContractError(
                        f"topology authority mismatch: {path}:{line_number}"
                    )
                region_id = str(row.get("region_id", ""))
                if not region_id:
                    raise AnalysisContractError(
                        f"topology row missing region_id: {path}:{line_number}"
                    )
                topology_key = f"{dataset}\x1f{region_id}"
                if topology_key in output:
                    raise AnalysisContractError(f"duplicate topology region_id: {topology_key}")
                output[topology_key] = {
                    "sample": row.get("sample"),
                    "chrom": row.get("chrom"),
                    "unit_id": row.get("unit_id"),
                    "phase_set": row.get("phase_set"),
                    "hp_family": row.get("hp_family"),
                    "active_positions": row.get("active_positions"),
                    "best_tree_unique": bool(row.get("best_tree_unique", False)),
                    "representative_best_edges": row.get("representative_best_edges", []),
                    "representative_best_vertices": row.get(
                        "representative_best_vertices", []
                    ),
                    "best_tree_tie_count": row.get("best_tree_tie_count"),
                }
        after_size = path.stat().st_size
        after_sha = sha256_file(path)
        if (
            after_size != before["size_bytes"]
            or after_sha != before["sha256"]
        ):
            raise AnalysisContractError(
                f"topology JSONL changed while loading: {path}"
            )
        sources.append(before)
    return output, sources


def validate_formal_topology_bindings(
    formal_counts: Mapping[tuple[str, str, str, str], Mapping[str, str]],
    topology: Mapping[str, Mapping[str, Any]],
) -> None:
    for (dataset, chrom, region_id, _hp_raw), row in formal_counts.items():
        key = f"{dataset}\x1f{region_id}"
        topology_row = topology.get(key)
        if topology_row is None:
            raise AnalysisContractError(
                f"formal unit absent from topology authority: {dataset} {region_id}"
            )
        expected = {
            "sample": dataset,
            "chrom": chrom,
            "unit_id": row.get("unit_id"),
            "phase_set": row.get("phase_set"),
            "hp_family": row.get("hp_family"),
            "active_positions": list(parse_int_list(row.get("active_positions"))),
        }
        observed = {field: topology_row.get(field) for field in expected}
        if observed != expected:
            raise AnalysisContractError(
                f"formal/topology identity mismatch for {dataset} {region_id}: "
                f"observed={observed} expected={expected}"
            )


def iter_candidate_groups(
    shard_paths: Sequence[Path],
    formal_counts: Mapping[tuple[str, str, str, str], Mapping[str, str]],
    expected_identities: Mapping[Path, Mapping[str, Any]] | None = None,
) -> Iterator[tuple[tuple[str, str, str, str], list[dict[str, str]]]]:
    """Stream census shards and retain only pre-methyl formal units."""
    for path in shard_paths:
        expected = (
            expected_identities.get(path.resolve())
            if expected_identities is not None
            else None
        )
        if expected_identities is not None and expected is None:
            raise AnalysisContractError(
                f"candidate shard lacks frozen identity: {path}"
            )
        if expected is not None:
            assert_file_identity(
                path, expected, label=f"candidate shard pre-read {path.name}"
            )
        opener = gzip.open if path.suffix == ".gz" else open
        current_key: tuple[str, str, str, str] | None = None
        rows: list[dict[str, str]] = []
        seen_keys: set[tuple[str, str, str, str]] = set()
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "dataset",
                "chrom",
                "region_id",
                "hp_raw",
                "qname_sha256",
                "pattern",
                "complete_pattern",
            }
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise AnalysisContractError(
                    f"candidate shard {path} missing columns: "
                    f"{sorted(required - set(reader.fieldnames or []))}"
                )
            for row in reader:
                key = (row["dataset"], row["chrom"], row["region_id"], row["hp_raw"])
                if key not in formal_counts:
                    continue
                if current_key is None:
                    current_key = key
                if key != current_key:
                    if key in seen_keys:
                        raise AnalysisContractError(
                            f"candidate shard is not grouped by unit key: {path} {key}"
                        )
                    seen_keys.add(current_key)
                    yield current_key, rows
                    current_key = key
                    rows = []
                rows.append(dict(row))
        if current_key is not None:
            yield current_key, rows
        if expected is not None:
            assert_file_identity(
                path, expected, label=f"candidate shard post-read {path.name}"
            )


def validate_candidate_unit(
    key: tuple[str, str, str, str],
    rows: Sequence[Mapping[str, str]],
    count_row: Mapping[str, str],
) -> None:
    """Bind a candidate-read unit exactly to its frozen census count row."""
    if not rows:
        raise AnalysisContractError(f"candidate unit contains no rows: {key}")
    identity_fields = (
        "dataset",
        "chrom",
        "region_id",
        "unit_id",
        "phase_set",
        "hp_family",
        "hp_raw",
        "active_positions",
        "n_active_bits",
    )
    expected = {field: str(count_row.get(field, "")) for field in identity_fields}
    complete_counts: Counter[str] = Counter()
    partial_counts: Counter[str] = Counter()
    seen_qnames: set[str] = set()
    n_active_bits = int(expected["n_active_bits"])
    for row in rows:
        observed = {field: str(row.get(field, "")) for field in identity_fields}
        if observed != expected:
            raise AnalysisContractError(
                f"candidate/count identity mismatch for {key}: "
                f"observed={observed} expected={expected}"
            )
        qname = str(row.get("qname_sha256", ""))
        if (
            len(qname) != 64
            or any(character not in "0123456789abcdef" for character in qname)
            or qname in seen_qnames
        ):
            raise AnalysisContractError(
                f"candidate qname identity invalid or duplicated for {key}: {qname}"
            )
        seen_qnames.add(qname)
        pattern = str(row.get("pattern", ""))
        if (
            len(pattern) != n_active_bits
            or not set(pattern).issubset({"R", "A", "X"})
        ):
            raise AnalysisContractError(
                f"candidate pattern invalid for {key}: {pattern!r}"
            )
        complete = bool_value(row.get("complete_pattern"))
        if complete != ("X" not in pattern):
            raise AnalysisContractError(
                f"candidate complete-pattern flag mismatch for {key}: {pattern!r}"
            )
        (complete_counts if complete else partial_counts)[pattern] += 1

    try:
        expected_complete = {
            str(state): int(count)
            for state, count in json.loads(
                str(count_row["state_count_json"])
            ).items()
        }
        expected_partial = {
            str(state): int(count)
            for state, count in json.loads(
                str(count_row["partial_state_count_json"])
            ).items()
        }
    except (AttributeError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise AnalysisContractError(
            f"invalid count JSON while binding candidate unit: {key}"
        ) from exc
    observed_counts = {
        "n_total": len(rows),
        "n_complete": sum(complete_counts.values()),
        "n_partial": sum(partial_counts.values()),
        "state_count_json": dict(sorted(complete_counts.items())),
        "partial_state_count_json": dict(sorted(partial_counts.items())),
    }
    expected_counts = {
        "n_total": int(count_row["n_total"]),
        "n_complete": int(count_row["n_complete"]),
        "n_partial": int(count_row["n_partial"]),
        "state_count_json": dict(sorted(expected_complete.items())),
        "partial_state_count_json": dict(sorted(expected_partial.items())),
    }
    if observed_counts != expected_counts:
        raise AnalysisContractError(
            f"candidate/count conservation mismatch for {key}: "
            f"observed={observed_counts} expected={expected_counts}"
        )


def load_methyl_union(
    rows: Sequence[Mapping[str, str]],
    active_positions: Sequence[int],
    catalog: Mapping[tuple[str, str, int], Mapping[str, str]],
) -> tuple[dict[str, dict[int, float]], dict[str, Any]]:
    dataset = rows[0]["dataset"]
    chrom = rows[0]["chrom"]
    targets = {row["qname_sha256"] for row in rows}
    methyl_by_qname: dict[str, dict[int, float]] = {target: {} for target in targets}
    marker_join_fractions: list[float] = []
    conflicts = 0

    for position in active_positions:
        marker_key = (dataset, chrom, int(position))
        if marker_key not in catalog:
            raise AnalysisContractError(f"formal marker absent from artifact catalog: {marker_key}")
        record = catalog[marker_key]
        reads_path = Path(record["reads_path"])
        methyl_path = Path(record["methylation_path"])
        reads_identity = catalog_file_identity(record, "reads")
        methyl_identity = catalog_file_identity(record, "methylation")
        assert_file_identity(
            reads_path,
            reads_identity,
            label=f"catalog {marker_key} reads pre-read",
        )
        read_id_to_qname: dict[str, str] = {}
        qname_to_read_id: dict[str, str] = {}
        with reads_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or not {"read_id", "read_name"}.issubset(
                reader.fieldnames
            ):
                raise AnalysisContractError(f"invalid reads header: {reads_path}")
            for read_row in reader:
                digest = sha256_read_name(read_row["read_name"])
                if digest not in targets:
                    continue
                read_id = read_row["read_id"]
                if digest in qname_to_read_id or read_id in read_id_to_qname:
                    raise AnalysisContractError(
                        f"ambiguous qname/read_id join at {reads_path}: {digest}"
                    )
                qname_to_read_id[digest] = read_id
                read_id_to_qname[read_id] = digest
        assert_file_identity(
            reads_path,
            reads_identity,
            label=f"catalog {marker_key} reads post-read",
        )
        marker_join_fractions.append(len(qname_to_read_id) / max(1, len(targets)))

        assert_file_identity(
            methyl_path,
            methyl_identity,
            label=f"catalog {marker_key} methylation pre-read",
        )
        seen_methyl_ids: set[str] = set()
        with methyl_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            if not header or header[0] != "read_id":
                raise AnalysisContractError(f"invalid methylation header: {methyl_path}")
            cpg_positions = [int(value) for value in header[1:]]
            for methyl_row in reader:
                if not methyl_row:
                    continue
                read_id = methyl_row[0]
                qname = read_id_to_qname.get(read_id)
                if qname is None:
                    continue
                if read_id in seen_methyl_ids:
                    raise AnalysisContractError(
                        f"duplicate methylation read_id at {methyl_path}: {read_id}"
                    )
                seen_methyl_ids.add(read_id)
                if len(methyl_row) != len(header):
                    raise AnalysisContractError(
                        f"methylation row width mismatch at {methyl_path}: {read_id}"
                    )
                destination = methyl_by_qname[qname]
                for cpg_position, value in zip(cpg_positions, methyl_row[1:]):
                    if value in {"", "NA", "NaN", "nan", "."}:
                        continue
                    probability = float(value)
                    if not 0.0 <= probability <= 1.0:
                        raise AnalysisContractError(
                            f"methyl probability outside [0,1] at {methyl_path}"
                        )
                    previous = destination.get(cpg_position)
                    if previous is not None and not math.isclose(
                        previous, probability, rel_tol=0.0, abs_tol=5e-5
                    ):
                        conflicts += 1
                    else:
                        destination[cpg_position] = probability
        assert_file_identity(
            methyl_path,
            methyl_identity,
            label=f"catalog {marker_key} methylation post-read",
        )
        if set(read_id_to_qname) != seen_methyl_ids:
            raise AnalysisContractError(
                f"reads/methylation selected ID mismatch at {methyl_path}"
            )
    if conflicts:
        raise AnalysisContractError(
            f"overlapping marker tiles have {conflicts} methylation value conflicts"
        )
    return methyl_by_qname, {
        "join_fraction_min": min(marker_join_fractions) if marker_join_fractions else 0.0,
        "overlap_conflicts": conflicts,
    }


def build_common_basis(
    rows: Sequence[Mapping[str, str]],
    methyl_by_qname: Mapping[str, Mapping[int, float]],
    active_positions: Sequence[int],
    config: AnalysisConfig,
) -> tuple[list[dict[str, str]], np.ndarray, tuple[int, ...], dict[str, int]]:
    complete_rows = [
        dict(row)
        for row in rows
        if bool_value(row["complete_pattern"])
        and set(row["pattern"]).issubset({"R", "A"})
    ]
    input_counts = Counter(row["pattern"] for row in complete_rows)
    eligible_states = sorted(
        state for state, count in input_counts.items() if count >= MIN_STATE_N
    )
    complete_rows = [row for row in complete_rows if row["pattern"] in eligible_states]
    if len(eligible_states) < 2 or len(complete_rows) < MIN_COMPLETE_N:
        raise ValueError("post_join_insufficient_complete_support")

    all_cpgs = sorted(
        {
            position
            for row in complete_rows
            for position in methyl_by_qname.get(row["qname_sha256"], {})
            if all(
                abs(position - marker) > config.marker_mask_bp
                for marker in active_positions
            )
        }
    )
    # Exchangeability, read coverage, and state-wise CpG coverage are mutually
    # dependent. Apply a conservative monotone fixed point: rows and CpGs may
    # only be removed, and every pre-methyl formal state stays frozen.
    retained_rows = complete_rows
    retained_basis = all_cpgs
    for _ in range(len(complete_rows) + len(all_cpgs) + 2):
        labels = [row["pattern"] for row in retained_rows]
        strata = [
            f"{row['read_group']}|{row['strand']}" for row in retained_rows
        ]
        exchangeable = exchangeable_indices(labels, strata)
        if len(exchangeable) == 0:
            raise ValueError("no_exchangeable_labels")
        exchangeable_rows = [
            retained_rows[index] for index in exchangeable.tolist()
        ]
        exchangeable_counts = Counter(
            row["pattern"] for row in exchangeable_rows
        )
        if (
            len(exchangeable_rows) < MIN_COMPLETE_N
            or any(
                exchangeable_counts[state] < MIN_STATE_N
                for state in eligible_states
            )
        ):
            raise ValueError("post_exchangeability_formal_state_support_lost")

        next_basis = [
            cpg
            for cpg in retained_basis
            if all(
                sum(
                    cpg in methyl_by_qname.get(row["qname_sha256"], {})
                    for row in exchangeable_rows
                    if row["pattern"] == state
                )
                / exchangeable_counts[state]
                >= config.state_cpg_coverage
                for state in eligible_states
            )
        ]
        if len(next_basis) < config.min_common_cpg:
            reason = (
                "insufficient_common_cpg"
                if retained_basis == all_cpgs
                else "post_read_filter_insufficient_common_cpg"
            )
            raise ValueError(reason)

        next_rows = []
        for row in exchangeable_rows:
            methyl = methyl_by_qname.get(row["qname_sha256"], {})
            finite_fraction = (
                sum(cpg in methyl for cpg in next_basis)
                / len(next_basis)
            )
            if finite_fraction >= config.read_cpg_coverage:
                next_rows.append(row)
        next_counts = Counter(row["pattern"] for row in next_rows)
        if (
            len(next_rows) < MIN_COMPLETE_N
            or any(next_counts[state] < MIN_STATE_N for state in eligible_states)
        ):
            raise ValueError("post_cpg_filter_formal_state_support_lost")

        if (
            [row["qname_sha256"] for row in next_rows]
            == [row["qname_sha256"] for row in retained_rows]
            and next_basis == retained_basis
        ):
            matrix = np.asarray(
                [
                    [
                        methyl_by_qname.get(row["qname_sha256"], {}).get(
                            cpg, math.nan
                        )
                        for cpg in next_basis
                    ]
                    for row in next_rows
                ],
                dtype=float,
            )
            return (
                next_rows,
                matrix,
                tuple(next_basis),
                dict(sorted(next_counts.items())),
            )
        retained_rows = next_rows
        retained_basis = next_basis
    raise ValueError("common_basis_filter_did_not_converge")


def restrict_to_exchangeable(
    rows: list[dict[str, str]],
    matrix: np.ndarray,
    distance: np.ndarray,
) -> tuple[list[dict[str, str]], np.ndarray, np.ndarray, list[str]]:
    """Retain exchangeable strata without silently dropping a formal state."""
    formal_states = set(row["pattern"] for row in rows)
    labels = [row["pattern"] for row in rows]
    strata = [f"{row['read_group']}|{row['strand']}" for row in rows]
    keep = exchangeable_indices(labels, strata)
    if len(keep) == 0:
        raise ValueError("no_exchangeable_labels")
    retained_rows = [rows[index] for index in keep.tolist()]
    retained_matrix = matrix[keep]
    retained_distance = distance[np.ix_(keep, keep)]
    counts = Counter(row["pattern"] for row in retained_rows)
    if (
        len(retained_rows) < MIN_COMPLETE_N
        or any(counts[state] < MIN_STATE_N for state in formal_states)
    ):
        raise ValueError("post_exchangeability_formal_state_support_lost")
    retained_strata = [
        f"{row['read_group']}|{row['strand']}" for row in retained_rows
    ]
    if len(exchangeable_indices(
        [row["pattern"] for row in retained_rows], retained_strata
    )) != len(retained_rows):
        raise ValueError("exchangeability_filter_did_not_converge")
    return retained_rows, retained_matrix, retained_distance, retained_strata


def equal_n_r2(
    distance: np.ndarray, labels: Sequence[str], seed: int
) -> tuple[float, int]:
    labels_array = np.asarray(labels, dtype=object)
    states = sorted(set(labels))
    minimum = min(int((labels_array == state).sum()) for state in states)
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    for state in states:
        candidates = np.flatnonzero(labels_array == state)
        selected.extend(sorted(rng.choice(candidates, size=minimum, replace=False).tolist()))
    selected_array = np.asarray(sorted(selected), dtype=int)
    _, r_squared = pseudo_f_r2(
        distance[np.ix_(selected_array, selected_array)],
        labels_array[selected_array].tolist(),
    )
    return r_squared, minimum


def rarefaction_r2(
    matrix: np.ndarray,
    labels: Sequence[str],
    repetitions: int,
    seed: int,
    min_pair_cpg: int,
) -> list[float]:
    rng = np.random.default_rng(seed)
    n_cpg = matrix.shape[1]
    subset_size = max(MIN_COMMON_CPG, int(round(0.80 * n_cpg)))
    if subset_size > n_cpg:
        return []
    values: list[float] = []
    for _ in range(repetitions):
        columns = np.sort(rng.choice(n_cpg, size=subset_size, replace=False))
        distance = bernoulli_distance_matrix(matrix[:, columns], min_pair_cpg)
        if not np.isfinite(distance).all():
            continue
        _, r_squared = pseudo_f_r2(distance, labels)
        values.append(r_squared)
    return values


def empty_result(
    count_row: Mapping[str, str],
    *,
    reason: str,
    evaluation_status: str = "NOT_EVALUABLE",
) -> dict[str, Any]:
    result = {field: "" for field in RESULT_FIELDS}
    for field in (
        "dataset",
        "chrom",
        "region_id",
        "unit_id",
        "phase_set",
        "hp_family",
        "hp_raw",
        "active_positions",
        "n_active_bits",
        "pair_full4",
        "k_ge_3",
    ):
        result[field] = count_row.get(field, "")
    result.update(
        {
            "schema_version": SCHEMA_VERSION,
            "input_n_complete": int(count_row.get("n_complete", 0)),
            "input_state_counts_json": count_row.get("state_count_json", "{}"),
            "assessment": "NOT_EVALUABLE",
            "evaluation_status": evaluation_status,
            "invalid_reason": reason,
            "permanova_permutations_requested": 0,
            "permanova_permutations_realized": 0,
        }
    )
    return result


def analyze_unit(
    key: tuple[str, str, str, str],
    rows: list[dict[str, str]],
    count_row: Mapping[str, str],
    catalog: Mapping[tuple[str, str, int], Mapping[str, str]],
    topology: Mapping[str, Mapping[str, Any]],
    config: AnalysisConfig,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    dataset, chrom, region_id, hp_raw = key
    if not rows:
        return empty_result(count_row, reason="no_candidate_rows"), None
    validate_candidate_unit(key, rows, count_row)
    active_positions = parse_int_list(count_row["active_positions"])
    if len(active_positions) < 2:
        return empty_result(count_row, reason="insufficient_active_positions"), None
    if any(
        (
            row["dataset"],
            row["chrom"],
            row["region_id"],
            row["hp_raw"],
        )
        != key
        for row in rows
    ):
        raise AnalysisContractError(f"candidate row escaped group key: {key}")
    base = empty_result(count_row, reason="")
    try:
        methyl_by_qname, join_audit = load_methyl_union(rows, active_positions, catalog)
        common_rows, matrix, cpg_positions, common_counts = build_common_basis(
            rows, methyl_by_qname, active_positions, config
        )
        complete_distance = bernoulli_distance_matrix(
            matrix, config.min_pair_cpg
        )
        invalid = ~np.isfinite(complete_distance)
        np.fill_diagonal(invalid, False)
        if invalid.any():
            raise ValueError("invalid_pair_distance")
        labels = [row["pattern"] for row in common_rows]
        strata = [
            f"{row['read_group']}|{row['strand']}" for row in common_rows
        ]
        if len(exchangeable_indices(labels, strata)) != len(common_rows):
            raise ValueError("exchangeability_filter_did_not_converge")
        counts = Counter(labels)
        unit_seed = stable_seed(config.seed, dataset, chrom, region_id, hp_raw)
        primary = restricted_permanova_permdisp(
            complete_distance, labels, strata, config.permutations, unit_seed
        )
        primary_r2 = float(primary["r_squared"])
        geometry_smd, geometry_feature = geometry_max_smd(common_rows, labels)
        topology_row = topology.get(f"{dataset}\x1f{region_id}")
        effects = pairwise_effects(complete_distance, labels, topology_row)
        best = effects[0] if effects else None
        if best is None:
            raise ValueError("no_pairwise_effect")

        balanced_r2, balanced_n = equal_n_r2(
            complete_distance, labels, stable_seed(unit_seed, "equal_n")
        )
        rarefaction_values = rarefaction_r2(
            matrix,
            labels,
            repetitions=10,
            seed=stable_seed(unit_seed, "rarefaction"),
            min_pair_cpg=config.min_pair_cpg,
        )
        rarefaction_median = (
            float(statistics.median(rarefaction_values))
            if rarefaction_values
            else math.nan
        )

        distal_columns = np.asarray(
            [
                index
                for index, position in enumerate(cpg_positions)
                if all(abs(position - marker) > config.distal_bp for marker in active_positions)
            ],
            dtype=int,
        )
        distal_r2 = math.nan
        distal_p = math.nan
        distal_permutations_realized = 0
        if len(distal_columns) >= config.min_common_cpg:
            distal_distance = bernoulli_distance_matrix(
                matrix[:, distal_columns], config.min_pair_cpg
            )
            if np.isfinite(distal_distance).all():
                _, distal_r2 = pseudo_f_r2(distal_distance, labels)
                if float(primary["p_value"]) <= 0.05 or primary_r2 >= 0.10:
                    distal_test = restricted_permanova_permdisp(
                        distal_distance,
                        labels,
                        strata,
                        config.sensitivity_permutations,
                        stable_seed(unit_seed, "distal"),
                    )
                    distal_p = float(distal_test["p_value"])
                    distal_permutations_realized = int(distal_test["realized"])

        base.update(
            {
                "analysis_n": len(common_rows),
                "analysis_state_counts_json": json.dumps(
                    dict(sorted(counts.items())), separators=(",", ":")
                ),
                "n_common_cpg": len(cpg_positions),
                "n_distal_cpg": len(distal_columns),
                "qname_join_fraction_min": join_audit["join_fraction_min"],
                "tile_overlap_conflicts": join_audit["overlap_conflicts"],
                "exchangeable_strata": len(set(strata)),
                "exchangeable_n": len(common_rows),
                "permanova_pseudo_f": primary["pseudo_f"],
                "permanova_r2": primary_r2,
                "permanova_p": primary["p_value"],
                "permanova_permutations_requested": primary["requested"],
                "permanova_permutations_realized": primary["realized"],
                "permdisp_f": primary["permdisp_f"],
                "permdisp_p": primary["permdisp_p"],
                "best_pair": f"{best['first']}|{best['second']}",
                "best_pair_hamming": best["hamming"],
                "best_pair_between_mean": best["between_mean"],
                "best_pair_pooled_within_mean": best["pooled_within_mean"],
                "best_pair_distance_contrast": best["distance_contrast"],
                "best_pair_standardized_effect": best["standardized_effect"],
                "best_pair_topology_relation": best["topology_relation"],
                "max_geometry_smd": geometry_smd,
                "geometry_feature": geometry_feature,
                "all_states_n8": all(count >= 8 for count in counts.values()),
                "all_states_n10": all(count >= 10 for count in counts.values()),
                "equal_n_r2": balanced_r2,
                "equal_n_retention": balanced_r2 / primary_r2
                if primary_r2 > 0
                else math.nan,
                "rarefaction_median_r2": rarefaction_median,
                "rarefaction_retention": rarefaction_median / primary_r2
                if primary_r2 > 0 and math.isfinite(rarefaction_median)
                else math.nan,
                "distal_r2": distal_r2,
                "distal_p": distal_p,
                "distal_permutations_realized": distal_permutations_realized,
                "distal_retention": distal_r2 / primary_r2
                if primary_r2 > 0 and math.isfinite(distal_r2)
                else math.nan,
                "assessment": "PENDING_MULTIPLICITY",
                "evaluation_status": "EVALUABLE",
                "invalid_reason": "",
            }
        )
        detail: dict[str, Any] | None = None
        if float(primary["p_value"]) <= 0.05 or primary_r2 >= 0.05:
            label_array = np.asarray(labels, dtype=object)
            profiles = {
                state: [
                    None if not math.isfinite(value) else float(value)
                    for value in np.nanmean(matrix[label_array == state], axis=0)
                ]
                for state in sorted(counts)
            }
            detail = {
                "schema_name": f"{SCHEMA_NAME}.detail",
                "schema_version": SCHEMA_VERSION,
                "dataset": dataset,
                "chrom": chrom,
                "region_id": region_id,
                "hp_raw": hp_raw,
                "active_positions": list(active_positions),
                "cpg_positions": list(cpg_positions),
                "state_counts": dict(sorted(counts.items())),
                "state_mean_profiles": profiles,
                "pairwise_effects": effects,
                "topology": topology_row,
                "balanced_n_per_state": balanced_n,
                "rarefaction_r2": rarefaction_values,
                "read_order": [
                    {
                        "qname_sha256": row["qname_sha256"],
                        "pattern": row["pattern"],
                        "read_group": row["read_group"],
                        "strand": row["strand"],
                    }
                    for row in common_rows
                ],
                "distance_matrix": complete_distance.round(6).tolist()
                if len(common_rows) <= 160
                else None,
            }
        return base, detail
    except ValueError as exc:
        return empty_result(count_row, reason=str(exc)), None
    except AnalysisContractError:
        raise
    except Exception as exc:
        raise AnalysisContractError(
            f"unexpected analysis failure for {key}: {type(exc).__name__}: {exc}"
        ) from exc


def initialize_worker(
    catalog: dict[tuple[str, str, int], dict[str, str]],
    topology: dict[str, dict[str, Any]],
    config: AnalysisConfig,
) -> None:
    global _WORKER_CATALOG, _WORKER_TOPOLOGY, _WORKER_CONFIG
    _WORKER_CATALOG = catalog
    _WORKER_TOPOLOGY = topology
    _WORKER_CONFIG = config


def worker_analyze(
    job: tuple[
        tuple[str, str, str, str],
        list[dict[str, str]],
        dict[str, str],
    ]
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    key, rows, count_row = job
    return analyze_unit(
        key,
        rows,
        count_row,
        _WORKER_CATALOG,
        _WORKER_TOPOLOGY,
        _WORKER_CONFIG,
    )


def assign_multiplicity_family(row: Mapping[str, Any]) -> str:
    if bool_value(row.get("pair_full4")) or bool_value(row.get("k_ge_3")):
        return "CONFIRMATORY_FULL4_OR_LONG"
    return "SECONDARY_PAIR_CONTRAST"


def adjust_p_values(rows: list[dict[str, Any]]) -> None:
    families: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        family = assign_multiplicity_family(row)
        row["multiplicity_family"] = family
        if row["evaluation_status"] == "EVALUABLE":
            families[family].append(index)
    for indices in families.values():
        p_values = [float(rows[index]["permanova_p"]) for index in indices]
        by_values = multipletests(p_values, alpha=0.05, method="fdr_by")[1]
        holm_values = multipletests(p_values, alpha=0.05, method="holm")[1]
        for index, by_value, holm_value in zip(indices, by_values, holm_values):
            rows[index]["q_by"] = float(by_value)
            rows[index]["p_holm"] = float(holm_value)


def finite_at_least(value: object, threshold: float) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(number) and number >= threshold


def classify_rows(rows: list[dict[str, Any]]) -> None:
    """Apply frozen gates after multiplicity without cross-tag causal inference."""
    for row in rows:
        if row["evaluation_status"] != "EVALUABLE":
            row["assessment"] = "NOT_EVALUABLE"
            continue
        significant_effect = (
            float(row.get("q_by", 1.0)) <= 0.05
            and float(row["permanova_r2"]) >= 0.10
            and float(row["best_pair_distance_contrast"]) >= 0.10
            and float(row["best_pair_standardized_effect"]) >= 0.50
        )
        raw_effect = (
            float(row["permanova_p"]) <= 0.05
            and float(row["permanova_r2"]) >= 0.10
        )
        confound_ok = (
            float(row["permdisp_p"]) >= 0.05
            and float(row["max_geometry_smd"]) < 0.50
        )
        robustness_ok = (
            bool_value(row["all_states_n8"])
            and finite_at_least(row["equal_n_retention"], 0.50)
            and finite_at_least(row["rarefaction_retention"], 0.50)
        )
        confirmatory = row["multiplicity_family"] == "CONFIRMATORY_FULL4_OR_LONG"
        if raw_effect and not confound_ok:
            row["assessment"] = "CONFOUNDED"
        elif significant_effect and confound_ok and robustness_ok and confirmatory:
            if not finite_at_least(row["distal_retention"], 0.50):
                row["assessment"] = "LOCAL_CIS_COMPATIBLE"
            else:
                row["assessment"] = "ROBUST_ASSOCIATION"
        else:
            row["assessment"] = "EVALUABLE_NO_ROBUST_ASSOCIATION"


def output_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.10g}"
    return str(value)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def load_candidate_shards(
    explicit_paths: Sequence[Path],
    manifest_path: Path | None,
) -> list[Path]:
    """Resolve candidate shards in deterministic, hash-bound order."""
    shards: list[Path] = []
    if manifest_path is not None:
        if not manifest_path.is_file():
            raise AnalysisContractError(f"candidate manifest missing: {manifest_path}")
        with manifest_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "sort_ordinal",
                "relative_path",
                "sha256",
                "size_bytes",
            }
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise AnalysisContractError(
                    "candidate manifest missing columns: "
                    f"{sorted(required - set(reader.fieldnames or []))}"
                )
            rows = sorted(
                (dict(row) for row in reader),
                key=lambda row: int(row["sort_ordinal"]),
            )
        ordinals = [int(row["sort_ordinal"]) for row in rows]
        if ordinals != list(range(len(rows))):
            raise AnalysisContractError(
                "candidate manifest sort_ordinal must be contiguous from zero"
            )
        for row in rows:
            relative_path = Path(row["relative_path"])
            if relative_path.is_absolute() or ".." in relative_path.parts:
                raise AnalysisContractError(
                    f"unsafe candidate manifest relative_path: {relative_path}"
                )
            shard = manifest_path.parent / relative_path
            if not shard.is_file():
                raise AnalysisContractError(f"candidate shard missing: {shard}")
            if shard.stat().st_size != int(row["size_bytes"]):
                raise AnalysisContractError(
                    f"candidate shard size mismatch: {shard}"
                )
            if sha256_file(shard) != row["sha256"]:
                raise AnalysisContractError(
                    f"candidate shard SHA-256 mismatch: {shard}"
                )
            shards.append(shard)

    shards.extend(explicit_paths)
    if not shards:
        raise AnalysisContractError(
            "provide --candidate-manifest and/or --candidate-shard"
        )
    resolved = [path.resolve() for path in shards]
    if len(set(resolved)) != len(resolved):
        raise AnalysisContractError("duplicate candidate shard path")
    for shard in resolved:
        if not shard.is_file():
            raise AnalysisContractError(f"candidate shard missing: {shard}")
    return resolved


def archive_failed_staging(staging: Path, output_parent: Path) -> Path:
    """Retain a failed analyzer staging directory without deleting evidence."""
    archive_root = output_parent / "_failed_staging_archive"
    archive_root.mkdir(parents=True, exist_ok=True)
    destination = archive_root / staging.name
    suffix = 1
    while os.path.lexists(destination):
        destination = archive_root / f"{staging.name}.{suffix}"
        suffix += 1
    os.rename(staging, destination)
    return destination


@contextlib.contextmanager
def deterministic_gzip_text_writer(
    path: Path, *, newline: str | None = None
) -> Iterator[io.TextIOWrapper]:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, compresslevel=6, mtime=0
        ) as gzip_handle:
            with io.TextIOWrapper(
                gzip_handle, encoding="utf-8", newline=newline, write_through=True
            ) as text_handle:
                yield text_handle


def write_outputs(
    output_dir: Path,
    rows: list[dict[str, Any]],
    details: list[dict[str, Any]],
    config: AnalysisConfig,
    sources: Mapping[str, Any],
    *,
    provisional_refinement: bool = False,
    published_output_dir: Path | None = None,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "pattern_methyl_evidence.v1.tsv.gz"
    detail_path = output_dir / "pattern_methyl_details.v1.jsonl.gz"
    summary_path = output_dir / "analysis_summary.json"
    for path in (evidence_path, detail_path, summary_path):
        if path.exists():
            raise AnalysisContractError(f"refusing to overwrite output: {path}")

    sorted_rows = sorted(
        rows,
        key=lambda row: (
            row["dataset"],
            natural_chrom_key(str(row["chrom"])),
            row["region_id"],
            row["hp_raw"],
        ),
    )
    with deterministic_gzip_text_writer(evidence_path, newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=RESULT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in sorted_rows:
            writer.writerow({field: output_value(row.get(field, "")) for field in RESULT_FIELDS})

    assessment_by_key = {
        (row["dataset"], row["chrom"], row["region_id"], row["hp_raw"]): row[
            "assessment"
        ]
        for row in sorted_rows
    }
    with deterministic_gzip_text_writer(detail_path) as handle:
        for detail in sorted(
            details,
            key=lambda row: (
                row["dataset"],
                natural_chrom_key(str(row["chrom"])),
                row["region_id"],
                row["hp_raw"],
            ),
        ):
            key = (
                detail["dataset"],
                detail["chrom"],
                detail["region_id"],
                detail["hp_raw"],
            )
            detail["assessment"] = assessment_by_key[key]
            handle.write(json.dumps(detail, ensure_ascii=False, separators=(",", ":")))
            handle.write("\n")

    assessment_counts = Counter(row["assessment"] for row in sorted_rows)
    invalid_counts = Counter(
        row["invalid_reason"]
        for row in sorted_rows
        if row["evaluation_status"] != "EVALUABLE"
    )
    family_counts = Counter(row["multiplicity_family"] for row in sorted_rows)
    summary = {
        "schema_name": f"{SCHEMA_NAME}.summary",
        "schema_version": SCHEMA_VERSION,
        "claim_ceiling": (
            "pattern-conditioned regional methylation association only; "
            "not ancestry, clone, causality, or topology rescoring"
        ),
        "result_status": (
            "PROVISIONAL_SUBSET_REFINEMENT_REQUIRES_FULL_FAMILY_MERGE"
            if provisional_refinement
            else "AUTHORITATIVE_WITHIN_DECLARED_FAMILY"
        ),
        "config": {
            "permutations": config.permutations,
            "sensitivity_permutations": config.sensitivity_permutations,
            "seed": config.seed,
            "marker_mask_bp": config.marker_mask_bp,
            "distal_bp": config.distal_bp,
            "min_common_cpg": config.min_common_cpg,
            "min_pair_cpg": config.min_pair_cpg,
            "state_cpg_coverage": config.state_cpg_coverage,
            "read_cpg_coverage": config.read_cpg_coverage,
        },
        "counts": {
            "units_total": len(sorted_rows),
            "units_evaluable": sum(
                row["evaluation_status"] == "EVALUABLE" for row in sorted_rows
            ),
            "detail_records": len(details),
            "assessment": dict(sorted(assessment_counts.items())),
            "multiplicity_family": dict(sorted(family_counts.items())),
            "invalid_reason": dict(sorted(invalid_counts.items())),
        },
        "sources": sources,
        "outputs": {},
    }
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    binding_root = (
        published_output_dir.resolve()
        if published_output_dir is not None
        else output_dir.resolve()
    )
    for label, path in {
        "evidence": evidence_path,
        "details": detail_path,
    }.items():
        summary["outputs"][label] = {
            "path": str(binding_root / path.name),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pattern-counts", type=Path, required=True)
    parser.add_argument("--artifact-catalog", type=Path, required=True)
    parser.add_argument(
        "--candidate-shard", type=Path, action="append", default=[]
    )
    parser.add_argument("--candidate-manifest", type=Path)
    parser.add_argument("--unit-key-file", type=Path)
    parser.add_argument("--unit-key-receipt", type=Path)
    parser.add_argument("--topology-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    parser.add_argument(
        "--sensitivity-permutations",
        type=int,
        default=DEFAULT_SENSITIVITY_PERMUTATIONS,
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--only-family",
        choices=("all", "confirmatory", "secondary"),
        default="all",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.workers <= 0 or args.permutations <= 0:
        raise AnalysisContractError("workers and permutations must be positive")
    if (args.unit_key_file is None) != (args.unit_key_receipt is None):
        raise AnalysisContractError(
            "--unit-key-file and --unit-key-receipt must be provided together"
        )
    config = AnalysisConfig(
        permutations=args.permutations,
        sensitivity_permutations=args.sensitivity_permutations,
        seed=args.seed,
    )
    output_dir = args.output_dir.expanduser().absolute()
    if os.path.lexists(output_dir):
        raise AnalysisContractError(f"refusing to overwrite output: {output_dir}")

    pattern_counts_identity = file_identity(args.pattern_counts)
    formal_counts = load_formal_counts(args.pattern_counts)
    assert_file_identity(
        args.pattern_counts,
        pattern_counts_identity,
        label="pattern_counts post-load",
    )
    if args.only_family != "all":
        want_confirmatory = args.only_family == "confirmatory"
        formal_counts = {
            key: row
            for key, row in formal_counts.items()
            if (
                bool_value(row.get("pair_full4")) or bool_value(row.get("k_ge_3"))
            )
            == want_confirmatory
        }
    family_datasets = {key[0] for key in formal_counts}
    adaptive_sources: dict[str, dict[str, Any] | None] = {
        "unit_key_file": None,
        "unit_key_receipt": None,
        "screen_evidence": None,
    }
    if args.unit_key_file is not None:
        selected_keys = load_unit_keys(args.unit_key_file)
        adaptive_sources = validate_adaptive_refinement_binding(
            args.unit_key_file,
            args.unit_key_receipt,
            analyzer_permutations=args.permutations,
            selected_key_count=len(selected_keys),
        )
        missing_keys = selected_keys - set(formal_counts)
        if missing_keys:
            raise AnalysisContractError(
                "unit-key selection contains units absent after family filtering: "
                f"{sorted(missing_keys)[:3]}"
            )
        formal_counts = {
            key: row for key, row in formal_counts.items() if key in selected_keys
        }
    artifact_catalog_identity = file_identity(args.artifact_catalog)
    catalog = load_catalog(args.artifact_catalog)
    assert_file_identity(
        args.artifact_catalog,
        artifact_catalog_identity,
        label="artifact_catalog post-load",
    )
    topology, topology_sources = load_topology(
        args.topology_root, family_datasets
    )
    validate_formal_topology_bindings(formal_counts, topology)
    candidate_manifest_identity = (
        file_identity(args.candidate_manifest)
        if args.candidate_manifest is not None
        else None
    )
    shards = load_candidate_shards(args.candidate_shard, args.candidate_manifest)
    if candidate_manifest_identity is not None:
        assert_file_identity(
            args.candidate_manifest,
            candidate_manifest_identity,
            label="candidate_manifest post-load",
        )
    shard_identities = {
        path.resolve(): file_identity(path)
        for path in shards
    }

    jobs = (
        (key, rows, dict(formal_counts[key]))
        for key, rows in iter_candidate_groups(
            shards, formal_counts, shard_identities
        )
    )
    results: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str, str]] = set()
    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=initialize_worker,
        initargs=(catalog, topology, config),
    ) as executor:
        for result, detail in executor.map(worker_analyze, jobs, chunksize=1):
            key = (
                result["dataset"],
                result["chrom"],
                result["region_id"],
                result["hp_raw"],
            )
            if key in seen:
                raise AnalysisContractError(f"duplicate analyzed unit: {key}")
            seen.add(key)
            results.append(result)
            if detail is not None:
                details.append(detail)
            if len(results) % 100 == 0:
                print(
                    f"analyzed={len(results)}/{len(formal_counts)} "
                    f"evaluable={sum(row['evaluation_status'] == 'EVALUABLE' for row in results)}",
                    file=sys.stderr,
                    flush=True,
                )
    missing = set(formal_counts) - seen
    if missing:
        preview = sorted(missing)[:3]
        raise AnalysisContractError(
            f"formal units missing from candidate shards: {len(missing)} preview={preview}"
        )
    provisional_refinement = args.unit_key_file is not None
    if provisional_refinement:
        for row in results:
            row["multiplicity_family"] = assign_multiplicity_family(row)
            row["q_by"] = ""
            row["p_holm"] = ""
            row["assessment"] = (
                "PROVISIONAL_REFINEMENT"
                if row["evaluation_status"] == "EVALUABLE"
                else "NOT_EVALUABLE"
            )
    else:
        adjust_p_values(results)
        classify_rows(results)
    sources = {
        "pattern_counts": pattern_counts_identity,
        "artifact_catalog": artifact_catalog_identity,
        "candidate_shards": [
            shard_identities[path.resolve()] for path in shards
        ],
        "candidate_manifest": candidate_manifest_identity,
        **adaptive_sources,
        "topology_root": str(args.topology_root.resolve()),
        "topology_jsonl": topology_sources,
    }
    def reattest_all_sources() -> None:
        assert_file_identity(
            args.pattern_counts,
            pattern_counts_identity,
            label="pattern_counts publish gate",
        )
        assert_file_identity(
            args.artifact_catalog,
            artifact_catalog_identity,
            label="artifact_catalog publish gate",
        )
        if candidate_manifest_identity is not None:
            assert_file_identity(
                args.candidate_manifest,
                candidate_manifest_identity,
                label="candidate_manifest publish gate",
            )
        for path in shards:
            assert_file_identity(
                path,
                shard_identities[path.resolve()],
                label=f"candidate shard publish gate {path.name}",
            )
        for source in topology_sources:
            assert_file_identity(
                Path(source["path"]),
                source,
                label=f"topology publish gate {source['dataset']}",
            )
        for label in ("unit_key_file", "unit_key_receipt", "screen_evidence"):
            identity = adaptive_sources.get(label)
            if identity is not None:
                assert_file_identity(
                    Path(identity["path"]),
                    identity,
                    label=f"adaptive {label} publish gate",
                )
        attest_catalog_artifacts(catalog)

    reattest_all_sources()
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.staging.",
            dir=output_dir.parent,
        )
    )
    try:
        summary = write_outputs(
            staging,
            results,
            details,
            config,
            sources,
            provisional_refinement=provisional_refinement,
            published_output_dir=output_dir,
        )
        reattest_all_sources()
        if os.path.lexists(output_dir):
            raise AnalysisContractError(
                f"output appeared during staging; refusing publish: {output_dir}"
            )
        os.rename(staging, output_dir)
    except Exception as exc:
        if os.path.lexists(staging):
            try:
                archived = archive_failed_staging(staging, output_dir.parent)
            except Exception as archive_exc:
                raise AnalysisContractError(
                    "analysis failed and staging archive failed; "
                    f"staging retained at {staging}; "
                    f"archive_error={type(archive_exc).__name__}: {archive_exc}"
                ) from exc
            raise AnalysisContractError(
                f"analysis failed; staging archived at {archived}: "
                f"{type(exc).__name__}: {exc}"
            ) from exc
        raise
    print(json.dumps(summary["counts"], ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AnalysisContractError as exc:
        print(f"FAIL_CLOSED: {exc}", file=sys.stderr)
        raise SystemExit(2)
