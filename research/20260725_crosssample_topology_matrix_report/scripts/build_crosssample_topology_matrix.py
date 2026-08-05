#!/usr/bin/env python3
"""Build fail-closed cross-dataset topology and VAF comparison matrices."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


DATASETS = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)

# KB-verified project metadata. HCC1937 is breast, not lung.
SAMPLE_METADATA = {
    "HCC1395": {"biological_id": "HCC1395", "cancer_type": "breast"},
    "HCC1395_DORADO": {
        "biological_id": "HCC1395",
        "cancer_type": "breast",
    },
    "COLO829": {"biological_id": "COLO829", "cancer_type": "melanoma"},
    "H1437": {
        "biological_id": "H1437",
        "cancer_type": "lung_adenocarcinoma",
    },
    "H2009": {
        "biological_id": "H2009",
        "cancer_type": "lung_adenocarcinoma",
    },
    "HCC1937": {"biological_id": "HCC1937", "cancer_type": "breast"},
    "HCC1954": {"biological_id": "HCC1954", "cancer_type": "breast"},
}

STRICT_METRICS = (
    "candidate_sites",
    "active_sites",
    "w_sites",
    "primary_edges",
    "exact_components",
)

COARSE_CATEGORIES = (
    "Single-only",
    "Sister-only",
    "Direct-only",
    "Sister+direct",
    "Unresolved-cross-coarse",
)

VAF_SOURCES = ("latest_lps_pass", "baseline_caller_pass")
VAF_METRICS = (
    "js_divergence_50bin_nats",
    "ks_statistic",
    "wasserstein_vaf",
)


class AnalysisError(RuntimeError):
    """Raised when a required input or invariant fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AnalysisError(message)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build 7x7 strict-linkage, exact-topology, coarse-composition, and "
            "truth-confirmed raw-VAF matrices."
        )
    )
    parser.add_argument("--strict-pairwise-json", required=True, type=Path)
    parser.add_argument(
        "--topology-root",
        required=True,
        type=Path,
        help="Root containing samples/<sample>/<sample>.topology.jsonl.",
    )
    parser.add_argument(
        "--census-root",
        required=True,
        type=Path,
        help="Root containing <sample>.census.jsonl and receipt.v2.json.",
    )
    parser.add_argument("--vaf-distance-tsv", required=True, type=Path)
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="A new, non-existing directory. Existing outputs are never overwritten.",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"required file is missing: {path}")
    return {
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def canonical_pair(left: str, right: str) -> tuple[str, str]:
    require(left in DATASETS and right in DATASETS, f"unknown dataset pair: {left}, {right}")
    require(left != right, f"self pair is not allowed in a long pair table: {left}")
    order = {sample: index for index, sample in enumerate(DATASETS)}
    return (left, right) if order[left] < order[right] else (right, left)


def biological_pair(left: str, right: str) -> tuple[str, str]:
    return tuple(sorted((left, right)))


def pair_class(left: str, right: str) -> str:
    left_meta = SAMPLE_METADATA[left]
    right_meta = SAMPLE_METADATA[right]
    if left_meta["biological_id"] == right_meta["biological_id"]:
        return "same_biological_id"
    if left_meta["cancer_type"] == right_meta["cancer_type"]:
        return "same_cancer_different_id"
    return "cross_cancer"


def safe_rate(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def parse_int_string(value: Any, label: str) -> int:
    require(
        isinstance(value, (str, int)) and str(value).isdigit(),
        f"{label} must be a non-negative decimal integer",
    )
    return int(value)


def float_or_none(value: Any) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def serialize(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return format(value, ".15g")
    return str(value)


def write_tsv_exclusive(path: Path, rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> None:
    require(rows, f"refusing to write an empty TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: serialize(row.get(column)) for column in columns})


def write_json_exclusive(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")


def write_matrix_exclusive(
    path: Path,
    matrix: Mapping[str, Mapping[str, float | int | None]],
) -> None:
    rows: list[dict[str, Any]] = []
    for row_sample in DATASETS:
        row: dict[str, Any] = {"dataset": row_sample}
        row.update({column_sample: matrix[row_sample][column_sample] for column_sample in DATASETS})
        rows.append(row)
    write_tsv_exclusive(path, rows, ("dataset", *DATASETS))


def blank_matrix(diagonal: float | int | None) -> dict[str, dict[str, float | int | None]]:
    return {
        left: {right: diagonal if left == right else None for right in DATASETS}
        for left in DATASETS
    }


def assert_symmetric_matrix(
    matrix: Mapping[str, Mapping[str, float | int | None]],
    label: str,
    *,
    diagonal: float | int | None,
    tolerance: float = 1e-12,
) -> None:
    for left in DATASETS:
        observed = matrix[left][left]
        if diagonal is None:
            require(observed is None, f"{label} diagonal must be null for {left}")
        elif isinstance(diagonal, float):
            require(
                observed is not None and math.isclose(float(observed), diagonal, abs_tol=tolerance),
                f"{label} diagonal mismatch for {left}: {observed}",
            )
        else:
            require(observed == diagonal, f"{label} diagonal mismatch for {left}: {observed}")
        for right in DATASETS:
            a = matrix[left][right]
            b = matrix[right][left]
            if a is None or b is None:
                require(a is None and b is None, f"{label} null asymmetry: {left}, {right}")
            elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
                require(
                    math.isclose(float(a), float(b), abs_tol=tolerance),
                    f"{label} asymmetry: {left}, {right}: {a} != {b}",
                )
            else:
                require(a == b, f"{label} asymmetry: {left}, {right}: {a} != {b}")


def assert_symmetric_count_matrix(
    matrix: Mapping[str, Mapping[str, float | int | None]],
    label: str,
) -> None:
    """Validate a symmetric count matrix whose self-count varies by dataset."""

    for left in DATASETS:
        observed = matrix[left][left]
        require(
            isinstance(observed, int) and not isinstance(observed, bool) and observed >= 0,
            f"{label} diagonal must be a non-negative integer for {left}: {observed}",
        )
        for right in DATASETS:
            a = matrix[left][right]
            b = matrix[right][left]
            require(a == b, f"{label} asymmetry: {left}, {right}: {a} != {b}")
            if a is not None:
                require(
                    isinstance(a, int) and not isinstance(a, bool) and a >= 0,
                    f"{label} must contain only non-negative integers or null: "
                    f"{left}, {right}: {a}",
                )


def assert_directional_diagonal(
    matrix: Mapping[str, Mapping[str, float | int | None]],
    label: str,
) -> None:
    for sample in DATASETS:
        observed = matrix[sample][sample]
        require(
            observed is not None and math.isclose(float(observed), 1.0, abs_tol=1e-12),
            f"{label} directional diagonal mismatch for {sample}: {observed}",
        )


@dataclass(frozen=True)
class JoinedTopology:
    sample: str
    chrom: str
    coordinate_signature: tuple[str, tuple[int, ...]]
    active_positions: tuple[int, ...]
    active_bit_count: int
    topology_signatures: frozenset[str]
    coarse_classes: frozenset[str]
    resolution_class: str
    recurrence_required: bool
    ancestor_relations: frozenset[tuple[int, int]] | None


@dataclass
class SampleTopologyData:
    sample: str
    coordinate_index: dict[tuple[str, tuple[int, ...]], list[JoinedTopology]]
    coarse_counts: Counter[str]
    topology_rows: int
    ranked_rows: int
    census_rows: int
    duplicate_coordinate_keys: int
    duplicate_coordinate_rows: int
    topology_join_key_duplicates: int
    census_join_key_duplicates: int


def verified_ancestor_relations(row: Mapping[str, Any]) -> frozenset[tuple[int, int]] | None:
    """Return strict mutation-position ancestor pairs for a verified unique tree."""

    if not bool(row.get("best_tree_unique")):
        return None
    if bool(row.get("recurrence_required")) or not bool(row.get("has_recfree")):
        return None
    positions = tuple(int(value) for value in row.get("active_positions") or ())
    active_k = int(row.get("active_bit_count", -1))
    require(active_k == len(positions), "active_bit_count/active_positions mismatch")
    require(len(set(positions)) == len(positions), "active_positions contain duplicates")
    edges = row.get("representative_best_edges")
    require(isinstance(edges, list), "unique-tree representative edges are missing")
    if active_k == 0:
        return frozenset()

    parent_by_child: dict[int, int] = {}
    acquisition_node_by_bit: dict[int, int] = {}
    vertices: set[int] = {0}
    for edge in edges:
        require(isinstance(edge, Mapping), "representative edge must be an object")
        parent = int(edge["parent_vertex"])
        child = int(edge["child_vertex"])
        bit = int(edge["acquired_active_bit"])
        acquired_position = int(edge["acquired_position"])
        require(child not in parent_by_child, "unique tree has a child with multiple parents")
        require(bit not in acquisition_node_by_bit, "recurrence-free tree reacquires an active bit")
        require(0 <= bit < active_k, f"acquired active bit is out of range: {bit}")
        require(
            acquired_position == positions[bit],
            "acquired position does not match active-bit coordinate",
        )
        parent_by_child[child] = parent
        acquisition_node_by_bit[bit] = child
        vertices.update((parent, child))

    require(
        set(acquisition_node_by_bit) == set(range(active_k)),
        "unique recurrence-free tree does not acquire every active bit exactly once",
    )
    require(len(edges) == len(vertices) - 1, "representative edges do not form a tree")

    ancestor_nodes_by_node: dict[int, set[int]] = {}
    for vertex in vertices:
        current = vertex
        ancestors: set[int] = set()
        seen: set[int] = set()
        while current != 0:
            require(current not in seen, "cycle detected in representative tree")
            seen.add(current)
            require(current in parent_by_child, f"representative tree is disconnected at {current}")
            current = parent_by_child[current]
            ancestors.add(current)
        ancestor_nodes_by_node[vertex] = ancestors

    relations: set[tuple[int, int]] = set()
    for left_bit, right_bit in itertools.permutations(range(active_k), 2):
        left_node = acquisition_node_by_bit[left_bit]
        right_node = acquisition_node_by_bit[right_bit]
        if left_node in ancestor_nodes_by_node[right_node]:
            relations.add((positions[left_bit], positions[right_bit]))
    return frozenset(relations)


def load_topology_sample(
    sample: str,
    topology_path: Path,
    census_path: Path,
) -> SampleTopologyData:
    topology_by_key: dict[tuple[str, int], Mapping[str, Any]] = {}
    topology_rows = 0
    ranked_rows = 0
    topology_duplicates = 0
    with topology_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            row = json.loads(line)
            topology_rows += 1
            require(row.get("sample") == sample, f"{topology_path}:{line_number} sample mismatch")
            if row.get("unit_status") != "ranked":
                continue
            ranked_rows += 1
            key = (str(row["unit_id"]), int(row["group_index"]))
            if key in topology_by_key:
                topology_duplicates += 1
            else:
                topology_by_key[key] = row

    census_by_key: dict[tuple[str, int], Mapping[str, Any]] = {}
    census_rows = 0
    census_duplicates = 0
    with census_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            row = json.loads(line)
            census_rows += 1
            require(row.get("sample") == sample, f"{census_path}:{line_number} sample mismatch")
            require(
                row.get("canonical_reproduction_pass") is True,
                f"{census_path}:{line_number} failed canonical reproduction",
            )
            key = (str(row["unit_id"]), int(row["group_index"]))
            if key in census_by_key:
                census_duplicates += 1
            else:
                census_by_key[key] = row

    require(topology_duplicates == 0, f"{sample} ranked topology join-key duplicates")
    require(census_duplicates == 0, f"{sample} census join-key duplicates")
    require(
        set(topology_by_key) == set(census_by_key),
        f"{sample} topology/census join-key sets differ",
    )
    require(ranked_rows == census_rows, f"{sample} ranked topology/census row mismatch")

    coordinate_index: dict[tuple[str, tuple[int, ...]], list[JoinedTopology]] = defaultdict(list)
    coarse_counts: Counter[str] = Counter()
    for key in sorted(topology_by_key):
        topology = topology_by_key[key]
        census = census_by_key[key]
        require(
            topology.get("chrom") == census.get("chrom"),
            f"{sample} topology/census chromosome mismatch at {key}",
        )
        require(
            int(topology["active_bit_count"]) == int(census["active_bit_count"]),
            f"{sample} topology/census active-k mismatch at {key}",
        )
        require(
            str(topology["best_tree_tie_count"]) == str(census["best_tree_tie_count"]),
            f"{sample} topology/census tie-count mismatch at {key}",
        )
        positions = tuple(sorted(int(value) for value in topology.get("active_positions") or ()))
        require(
            len(positions) == int(topology["active_bit_count"]),
            f"{sample} coordinate signature active-k mismatch at {key}",
        )
        require(len(set(positions)) == len(positions), f"{sample} duplicate active coordinate at {key}")

        signatures_payload = census.get("topology_signatures")
        require(isinstance(signatures_payload, list) and signatures_payload, "missing signatures")
        signature_set = frozenset(str(entry["shape_signature"]) for entry in signatures_payload)
        coarse_set = frozenset(str(entry["coarse_class"]) for entry in signatures_payload)
        require(
            len(signature_set) == int(census["topology_signature_count"]),
            f"{sample} topology-signature set count mismatch at {key}",
        )
        require(
            len(coarse_set) == int(census["coarse_class_count"]),
            f"{sample} coarse-class set count mismatch at {key}",
        )
        resolution = str(census["resolution_class"])
        if resolution == "UNIQUE_TREE":
            require(len(signature_set) == 1, f"{sample} UNIQUE_TREE has multiple signatures at {key}")
            require(
                topology.get("best_tree_unique") is True,
                f"{sample} census/topology unique-tree mismatch at {key}",
            )

        if len(coarse_set) == 1:
            coarse_class = next(iter(coarse_set))
            require(coarse_class in COARSE_CATEGORIES[:-1], f"unknown coarse class: {coarse_class}")
        else:
            coarse_class = "Unresolved-cross-coarse"
        coarse_counts[coarse_class] += 1

        ancestor_relations = (
            verified_ancestor_relations(topology) if resolution == "UNIQUE_TREE" else None
        )
        coordinate_signature = (str(topology["chrom"]), positions)
        coordinate_index[coordinate_signature].append(
            JoinedTopology(
                sample=sample,
                chrom=str(topology["chrom"]),
                coordinate_signature=coordinate_signature,
                active_positions=positions,
                active_bit_count=len(positions),
                topology_signatures=signature_set,
                coarse_classes=coarse_set,
                resolution_class=resolution,
                recurrence_required=bool(topology.get("recurrence_required")),
                ancestor_relations=ancestor_relations,
            )
        )

    require(sum(coarse_counts.values()) == ranked_rows, f"{sample} coarse census not conserved")
    duplicate_keys = sum(len(records) > 1 for records in coordinate_index.values())
    duplicate_rows = sum(
        len(records) for records in coordinate_index.values() if len(records) > 1
    )
    return SampleTopologyData(
        sample=sample,
        coordinate_index=dict(coordinate_index),
        coarse_counts=coarse_counts,
        topology_rows=topology_rows,
        ranked_rows=ranked_rows,
        census_rows=census_rows,
        duplicate_coordinate_keys=duplicate_keys,
        duplicate_coordinate_rows=duplicate_rows,
        topology_join_key_duplicates=topology_duplicates,
        census_join_key_duplicates=census_duplicates,
    )


def exact_coordinate_diagnostics(
    left: SampleTopologyData,
    right: SampleTopologyData,
) -> dict[str, Any]:
    shared_keys = set(left.coordinate_index) & set(right.coordinate_index)
    one_to_one_keys = {
        key
        for key in shared_keys
        if len(left.coordinate_index[key]) == 1 and len(right.coordinate_index[key]) == 1
    }
    signature_equal = 0
    signature_compatible = 0
    coarse_equal = 0
    coarse_compatible = 0
    both_unique = 0
    both_unique_shape_equal = 0
    both_unique_k_ge_2 = 0
    both_unique_k_ge_2_shape_equal = 0
    both_unique_ancestor_evaluable = 0
    both_unique_ancestor_equal = 0

    for key in one_to_one_keys:
        left_row = left.coordinate_index[key][0]
        right_row = right.coordinate_index[key][0]
        if left_row.topology_signatures == right_row.topology_signatures:
            signature_equal += 1
        if left_row.topology_signatures & right_row.topology_signatures:
            signature_compatible += 1
        if left_row.coarse_classes == right_row.coarse_classes:
            coarse_equal += 1
        if left_row.coarse_classes & right_row.coarse_classes:
            coarse_compatible += 1

        unique = (
            left_row.resolution_class == "UNIQUE_TREE"
            and right_row.resolution_class == "UNIQUE_TREE"
        )
        if not unique:
            continue
        both_unique += 1
        shape_equal = left_row.topology_signatures == right_row.topology_signatures
        if shape_equal:
            both_unique_shape_equal += 1
        if left_row.active_bit_count < 2:
            continue
        both_unique_k_ge_2 += 1
        if shape_equal:
            both_unique_k_ge_2_shape_equal += 1
        if (
            left_row.ancestor_relations is not None
            and right_row.ancestor_relations is not None
        ):
            both_unique_ancestor_evaluable += 1
            if left_row.ancestor_relations == right_row.ancestor_relations:
                both_unique_ancestor_equal += 1

    one_to_one = len(one_to_one_keys)
    return {
        "exact_coordinate_shared_signature_n": len(shared_keys),
        "exact_coordinate_one_to_one_n": one_to_one,
        "exact_coordinate_duplicate_excluded_signature_n": len(shared_keys) - one_to_one,
        "topology_signature_set_equal_n": signature_equal,
        "topology_signature_set_equal_rate": safe_rate(signature_equal, one_to_one),
        "topology_signature_set_compatible_n": signature_compatible,
        "topology_signature_set_compatible_rate": safe_rate(signature_compatible, one_to_one),
        "coarse_class_set_equal_n": coarse_equal,
        "coarse_class_set_equal_rate": safe_rate(coarse_equal, one_to_one),
        "coarse_class_set_compatible_n": coarse_compatible,
        "coarse_class_set_compatible_rate": safe_rate(coarse_compatible, one_to_one),
        "both_unique_n": both_unique,
        "both_unique_shape_equal_n": both_unique_shape_equal,
        "both_unique_shape_equal_rate": safe_rate(both_unique_shape_equal, both_unique),
        "both_unique_k_ge_2_n": both_unique_k_ge_2,
        "both_unique_k_ge_2_shape_equal_n": both_unique_k_ge_2_shape_equal,
        "both_unique_k_ge_2_shape_equal_rate": safe_rate(
            both_unique_k_ge_2_shape_equal, both_unique_k_ge_2
        ),
        "both_unique_ancestor_evaluable_n": both_unique_ancestor_evaluable,
        "both_unique_ancestor_equal_n": both_unique_ancestor_equal,
        "both_unique_ancestor_equal_rate": safe_rate(
            both_unique_ancestor_equal, both_unique_ancestor_evaluable
        ),
    }


def composition_metrics(
    left: SampleTopologyData,
    right: SampleTopologyData,
) -> dict[str, float]:
    left_total = sum(left.coarse_counts.values())
    right_total = sum(right.coarse_counts.values())
    require(left_total > 0 and right_total > 0, "coarse composition has a zero denominator")
    left_p = [left.coarse_counts[category] / left_total for category in COARSE_CATEGORIES]
    right_p = [right.coarse_counts[category] / right_total for category in COARSE_CATEGORIES]
    tvd = 0.5 * sum(abs(a - b) for a, b in zip(left_p, right_p))
    midpoint = [(a + b) / 2.0 for a, b in zip(left_p, right_p)]

    def kl_divergence(probabilities: Sequence[float], reference: Sequence[float]) -> float:
        return sum(
            p * math.log(p / q)
            for p, q in zip(probabilities, reference)
            if p > 0.0
        )

    jsd_nats = 0.5 * kl_divergence(left_p, midpoint) + 0.5 * kl_divergence(
        right_p, midpoint
    )
    return {
        "coarse_topology_tvd": tvd,
        "coarse_topology_similarity": 1.0 - tvd,
        "coarse_topology_jsd_nats": jsd_nats,
        "coarse_topology_js_distance": math.sqrt(jsd_nats / math.log(2.0)),
        "coarse_topology_js_similarity": 1.0
        - math.sqrt(jsd_nats / math.log(2.0)),
        "coarse_topology_js_divergence_complement": 1.0
        - jsd_nats / math.log(2.0),
    }


def load_vaf_distances(path: Path) -> dict[tuple[str, str], dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(rows, "VAF distribution-distance TSV is empty")
    expected_pairs = {
        canonical_pair(left, right) for left, right in itertools.combinations(DATASETS, 2)
    }
    indexed: dict[tuple[str, str, str, str], Mapping[str, str]] = {}
    for row in rows:
        source = row["source"]
        metric = row["metric"]
        require(source in VAF_SOURCES, f"unexpected VAF source: {source}")
        require(metric in VAF_METRICS, f"unexpected VAF metric: {metric}")
        pair = canonical_pair(row["sample_a"], row["sample_b"])
        key = (source, metric, *pair)
        require(key not in indexed, f"duplicate VAF distance row: {key}")
        require(
            row["distribution_scope"] == "truth_confirmed_marginal_raw_vaf",
            f"unexpected VAF distribution scope: {row['distribution_scope']}",
        )
        indexed[key] = row

    for source in VAF_SOURCES:
        for metric in VAF_METRICS:
            observed = {
                (left, right)
                for row_source, row_metric, left, right in indexed
                if row_source == source and row_metric == metric
            }
            require(
                observed == expected_pairs,
                f"VAF source/metric pair coverage mismatch: {source}/{metric}",
            )

    result: dict[tuple[str, str], dict[str, Any]] = {}
    for pair in sorted(expected_pairs, key=lambda value: (DATASETS.index(value[0]), DATASETS.index(value[1]))):
        record: dict[str, Any] = {}
        for source in VAF_SOURCES:
            source_prefix = "vaf_latest" if source == "latest_lps_pass" else "vaf_baseline"
            for metric in VAF_METRICS:
                row = indexed[(source, metric, *pair)]
                distance = float(row["distance"])
                short_metric = {
                    "js_divergence_50bin_nats": "jsd_nats",
                    "ks_statistic": "ks",
                    "wasserstein_vaf": "wasserstein",
                }[metric]
                record[f"{source_prefix}_{short_metric}"] = distance
                if metric == "js_divergence_50bin_nats":
                    require(
                        -1e-15 <= distance <= math.log(2.0) + 1e-12,
                        f"JSD is outside [0, ln(2)]: {source}/{pair}/{distance}",
                    )
                    normalized_js_distance = math.sqrt(distance / math.log(2.0))
                    record[f"{source_prefix}_js_distance"] = normalized_js_distance
                    record[f"{source_prefix}_js_similarity"] = (
                        1.0 - normalized_js_distance
                    )
                    record[f"{source_prefix}_js_divergence_complement"] = (
                        1.0 - distance / math.log(2.0)
                    )
                    if row["sample_a"] == pair[0]:
                        record[f"{source_prefix}_n_a"] = int(row["n_a"])
                        record[f"{source_prefix}_n_b"] = int(row["n_b"])
                    else:
                        record[f"{source_prefix}_n_a"] = int(row["n_b"])
                        record[f"{source_prefix}_n_b"] = int(row["n_a"])
        result[pair] = record
    return result


def mean_min_max(values: Iterable[float | int | None]) -> tuple[float | None, float | None, float | None]:
    clean = [float(value) for value in values if value is not None]
    if not clean:
        return None, None, None
    return statistics.fmean(clean), min(clean), max(clean)


def aggregate_rows(
    rows: Sequence[Mapping[str, Any]],
    group_fields: Sequence[str],
    numeric_fields: Sequence[str],
) -> list[dict[str, Any]]:
    groups: dict[tuple[Any, ...], list[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[tuple(row[field] for field in group_fields)].append(row)
    output: list[dict[str, Any]] = []
    for group_key in sorted(groups, key=lambda value: tuple(str(item) for item in value)):
        members = groups[group_key]
        result = {field: value for field, value in zip(group_fields, group_key)}
        result["pair_n"] = len(members)
        for field in numeric_fields:
            mean, minimum, maximum = mean_min_max(float_or_none(row.get(field)) for row in members)
            result[f"{field}_mean"] = mean
            result[f"{field}_min"] = minimum
            result[f"{field}_max"] = maximum
        output.append(result)
    return output


def exact_group_label_permutation(
    biological_pair_rows: Sequence[Mapping[str, Any]],
    metric_fields: Sequence[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Enumerate all 6!/(3!2!1!) group assignments and test within-group enrichment."""

    biological_ids = sorted(
        {
            str(row[field])
            for row in biological_pair_rows
            for field in ("biological_id_a", "biological_id_b")
        }
    )
    require(len(biological_ids) == 6, "exact permutation requires six biological IDs")
    expected_pairs = {
        biological_pair(left, right)
        for left, right in itertools.combinations(biological_ids, 2)
    }
    indexed = {
        biological_pair(
            str(row["biological_id_a"]),
            str(row["biological_id_b"]),
        ): row
        for row in biological_pair_rows
    }
    require(
        set(indexed) == expected_pairs and len(indexed) == 15,
        "biological pair table is incomplete for exact permutation",
    )

    observed_groups = {
        3: frozenset(("HCC1395", "HCC1937", "HCC1954")),
        2: frozenset(("H1437", "H2009")),
        1: frozenset(("COLO829",)),
    }
    permutation_rows: list[dict[str, Any]] = []
    assignment_id = 0
    all_ids = frozenset(biological_ids)
    for size3_tuple in itertools.combinations(biological_ids, 3):
        size3 = frozenset(size3_tuple)
        remaining_after_size3 = sorted(all_ids - size3)
        for size2_tuple in itertools.combinations(remaining_after_size3, 2):
            size2 = frozenset(size2_tuple)
            size1 = all_ids - size3 - size2
            require(len(size1) == 1, "permutation singleton group has invalid size")
            assignment_id += 1
            group_by_id = {
                biological_id: group_size
                for group_size, members in ((3, size3), (2, size2), (1, size1))
                for biological_id in members
            }
            within_pairs = [
                pair
                for pair in sorted(expected_pairs)
                if group_by_id[pair[0]] == group_by_id[pair[1]]
            ]
            between_pairs = sorted(expected_pairs - set(within_pairs))
            require(
                len(within_pairs) == 4 and len(between_pairs) == 11,
                "permutation within/between denominator mismatch",
            )
            is_observed = (
                size3 == observed_groups[3]
                and size2 == observed_groups[2]
                and size1 == observed_groups[1]
            )
            for metric in metric_fields:
                source_field = f"{metric}_mean"
                within_values = [
                    float_or_none(indexed[pair].get(source_field))
                    for pair in within_pairs
                ]
                between_values = [
                    float_or_none(indexed[pair].get(source_field))
                    for pair in between_pairs
                ]
                require(
                    all(value is not None for value in within_values + between_values),
                    f"exact permutation metric is incomplete: {metric}",
                )
                within_mean = statistics.fmean(
                    float(value) for value in within_values if value is not None
                )
                between_mean = statistics.fmean(
                    float(value) for value in between_values if value is not None
                )
                permutation_rows.append(
                    {
                        "assignment_id": assignment_id,
                        "size3_ids": ";".join(sorted(size3)),
                        "size2_ids": ";".join(sorted(size2)),
                        "size1_id": next(iter(size1)),
                        "is_observed_cancer_assignment": is_observed,
                        "metric": metric,
                        "within_pair_n": len(within_pairs),
                        "between_pair_n": len(between_pairs),
                        "within_mean": within_mean,
                        "between_mean": between_mean,
                        "within_minus_between": within_mean - between_mean,
                    }
                )

    require(assignment_id == 60, "exact permutation assignment count is not 60")
    tests: list[dict[str, Any]] = []
    for metric in metric_fields:
        metric_rows = [row for row in permutation_rows if row["metric"] == metric]
        observed = [
            row for row in metric_rows if row["is_observed_cancer_assignment"]
        ]
        require(
            len(metric_rows) == 60 and len(observed) == 1,
            f"exact permutation coverage mismatch: {metric}",
        )
        observed_difference = float(observed[0]["within_minus_between"])
        greater_or_equal_n = sum(
            float(row["within_minus_between"]) >= observed_difference - 1e-15
            for row in metric_rows
        )
        tests.append(
            {
                "metric": metric,
                "technical_source_handling": (
                    "HCC1395 and HCC1395_DORADO technical pair values are "
                    "macro-averaged before the 6-biological-ID test"
                ),
                "group_sizes": "3/2/1",
                "biological_id_n": 6,
                "within_pair_n": observed[0]["within_pair_n"],
                "between_pair_n": observed[0]["between_pair_n"],
                "observed_within_mean": observed[0]["within_mean"],
                "observed_between_mean": observed[0]["between_mean"],
                "observed_within_minus_between": observed_difference,
                "permutation_n": len(metric_rows),
                "greater_or_equal_n": greater_or_equal_n,
                "one_sided_exact_p": greater_or_equal_n / len(metric_rows),
                "alternative": "within-group similarity > between-group similarity",
            }
        )
    return tests, permutation_rows


def main() -> int:
    args = parse_args()
    require(not args.output_dir.exists(), f"output directory already exists: {args.output_dir}")
    for path in (
        args.strict_pairwise_json,
        args.vaf_distance_tsv,
        args.census_root / "summary.json",
        args.census_root / "receipt.v2.json",
        args.topology_root / "summary" / "all7_summary.json",
        args.topology_root / "cohort_receipt.json",
    ):
        require(path.is_file(), f"required input is missing: {path}")

    source_paths: dict[str, Path] = {
        "strict_pairwise_json": args.strict_pairwise_json,
        "vaf_distance_tsv": args.vaf_distance_tsv,
        "census_summary": args.census_root / "summary.json",
        "census_receipt": args.census_root / "receipt.v2.json",
        "topology_summary": args.topology_root / "summary" / "all7_summary.json",
        "topology_cohort_receipt": args.topology_root / "cohort_receipt.json",
        "analysis_script": Path(__file__).resolve(),
    }
    for sample in DATASETS:
        source_paths[f"topology_{sample}"] = (
            args.topology_root / "samples" / sample / f"{sample}.topology.jsonl"
        )
        source_paths[f"census_{sample}"] = args.census_root / f"{sample}.census.jsonl"
    source_identities_before = {label: identity(path) for label, path in source_paths.items()}

    strict = json.loads(args.strict_pairwise_json.read_text(encoding="utf-8"))
    require(strict.get("all_validations_pass") is True, "strict pairwise validations did not pass")
    require(tuple(strict["scope"]["datasets"]) == DATASETS, "strict dataset order/scope mismatch")
    require(int(strict["scope"]["pair_count"]) == 21, "strict pair count is not 21")
    require(
        strict["definitions"]["phase_invariance"]
        == "PS and HP labels are intentionally excluded across datasets",
        "strict phase-invariance definition changed",
    )

    expected_pairs = {
        canonical_pair(left, right) for left, right in itertools.combinations(DATASETS, 2)
    }
    strict_pairs: dict[tuple[str, str], Mapping[str, Any]] = {}
    for row in strict["pairwise"]:
        pair = canonical_pair(row["left"], row["right"])
        require(pair not in strict_pairs, f"duplicate strict pair: {pair}")
        strict_pairs[pair] = row
        for metric in STRICT_METRICS:
            payload = row[metric]
            require(
                payload["union_n"]
                == payload["left_n"] + payload["right_n"] - payload["intersection_n"],
                f"strict union conservation failed: {pair}/{metric}",
            )
            expected_jaccard = payload["intersection_n"] / payload["union_n"]
            require(
                math.isclose(payload["jaccard"], expected_jaccard, abs_tol=1e-15),
                f"strict Jaccard mismatch: {pair}/{metric}",
            )
    require(set(strict_pairs) == expected_pairs, "strict pair coverage is not the complete 7 choose 2")

    topology_summary = json.loads(
        (args.topology_root / "summary" / "all7_summary.json").read_text(encoding="utf-8")
    )
    topology_receipt_sha: dict[str, str] = {}
    for sample_row in topology_summary["samples"]:
        sample = sample_row["sample"]
        for source_identity in sample_row["source_identities"]:
            if source_identity["path"].endswith(f"/{sample}.topology.jsonl"):
                topology_receipt_sha[sample] = source_identity["sha256"]
    require(set(topology_receipt_sha) == set(DATASETS), "topology summary lacks JSONL identities")

    census_receipt = json.loads((args.census_root / "receipt.v2.json").read_text(encoding="utf-8"))
    require(census_receipt["checks"]["all_pass"] is True, "census receipt is not all-pass")
    require(
        tuple(census_receipt["scope"]["datasets"]) == DATASETS,
        "census receipt dataset order/scope mismatch",
    )

    samples: dict[str, SampleTopologyData] = {}
    for sample in DATASETS:
        topology_path = source_paths[f"topology_{sample}"]
        census_path = source_paths[f"census_{sample}"]
        require(
            source_identities_before[f"topology_{sample}"]["sha256"]
            == topology_receipt_sha[sample],
            f"{sample} topology SHA does not match all7 summary",
        )
        expected_census = census_receipt["sample_outputs"][sample]["census"]
        require(
            source_identities_before[f"census_{sample}"]["sha256"]
            == expected_census["sha256"],
            f"{sample} census SHA does not match receipt.v2",
        )
        samples[sample] = load_topology_sample(sample, topology_path, census_path)
        require(
            samples[sample].census_rows
            == int(census_receipt["sample_outputs"][sample]["row_count"]),
            f"{sample} census row count does not match receipt.v2",
        )

    vaf_by_pair = load_vaf_distances(args.vaf_distance_tsv)

    sample_composition_rows: list[dict[str, Any]] = []
    for sample in DATASETS:
        data = samples[sample]
        total = sum(data.coarse_counts.values())
        row: dict[str, Any] = {
            "dataset": sample,
            "biological_id": SAMPLE_METADATA[sample]["biological_id"],
            "cancer_type": SAMPLE_METADATA[sample]["cancer_type"],
            "ranked_topology_units": total,
            "topology_jsonl_rows": data.topology_rows,
            "census_jsonl_rows": data.census_rows,
            "duplicate_coordinate_keys": data.duplicate_coordinate_keys,
            "rows_in_duplicate_coordinate_keys": data.duplicate_coordinate_rows,
        }
        for category in COARSE_CATEGORIES:
            slug = category.lower().replace("+", "_plus_").replace("-", "_")
            row[f"coarse_{slug}_n"] = data.coarse_counts[category]
            row[f"coarse_{slug}_fraction"] = data.coarse_counts[category] / total
        sample_composition_rows.append(row)

    diagnostics_by_pair: dict[tuple[str, str], dict[str, Any]] = {}
    composition_by_pair: dict[tuple[str, str], dict[str, float]] = {}
    for pair in expected_pairs:
        diagnostics_by_pair[pair] = exact_coordinate_diagnostics(
            samples[pair[0]], samples[pair[1]]
        )
        composition_by_pair[pair] = composition_metrics(samples[pair[0]], samples[pair[1]])

    pair_rows: list[dict[str, Any]] = []
    for left, right in itertools.combinations(DATASETS, 2):
        pair = canonical_pair(left, right)
        strict_row = strict_pairs[pair]
        row: dict[str, Any] = {
            "dataset_a": pair[0],
            "dataset_b": pair[1],
            "biological_id_a": SAMPLE_METADATA[pair[0]]["biological_id"],
            "biological_id_b": SAMPLE_METADATA[pair[1]]["biological_id"],
            "cancer_type_a": SAMPLE_METADATA[pair[0]]["cancer_type"],
            "cancer_type_b": SAMPLE_METADATA[pair[1]]["cancer_type"],
            "pair_class": pair_class(*pair),
            "is_hcc1395_technical_pair": set(pair)
            == {"HCC1395", "HCC1395_DORADO"},
        }
        for metric in STRICT_METRICS:
            metric_row = strict_row[metric]
            row[f"{metric}_intersection_n"] = metric_row["intersection_n"]
            row[f"{metric}_union_n"] = metric_row["union_n"]
            row[f"{metric}_jaccard"] = metric_row["jaccard"]
            row[f"{metric}_overlap_coefficient"] = metric_row["overlap_coefficient"]
            if strict_row["left"] == pair[0]:
                row[f"{metric}_a_contained_in_b"] = metric_row["left_recall"]
                row[f"{metric}_b_contained_in_a"] = metric_row["right_recall"]
            else:
                row[f"{metric}_a_contained_in_b"] = metric_row["right_recall"]
                row[f"{metric}_b_contained_in_a"] = metric_row["left_recall"]
        row["joint_w_coverage_of_shared_candidates"] = strict_row[
            "joint_w_coverage_of_shared_candidates"
        ]
        row.update(composition_by_pair[pair])
        row.update(diagnostics_by_pair[pair])
        row.update(vaf_by_pair[pair])

        row["shared_candidate_projection_available"] = False
        for field in (
            "shared_candidate_sites_n",
            "shared_candidate_active_dor_to_hcc",
            "shared_candidate_w_dor_to_hcc",
            "shared_candidate_primary_edges_dor_to_hcc",
            "shared_candidate_alt_informative_edges_dor_to_hcc",
            "shared_candidate_aa_edges_dor_to_hcc",
            "shared_candidate_co_membership_dor_to_hcc",
        ):
            row[field] = None
        if row["is_hcc1395_technical_pair"]:
            target = strict["target_pair"]
            projection = target["shared_candidate_projection"]
            require(
                target["left"] == "HCC1395" and target["right"] == "HCC1395_DORADO",
                "target shared-candidate direction changed",
            )
            row["shared_candidate_projection_available"] = True
            row["shared_candidate_sites_n"] = projection["shared_candidate_sites_n"]
            row["shared_candidate_active_dor_to_hcc"] = projection["active_sites"][
                "right_recall"
            ]
            row["shared_candidate_w_dor_to_hcc"] = projection["w_sites"]["right_recall"]
            row["shared_candidate_primary_edges_dor_to_hcc"] = projection[
                "primary_edges"
            ]["right_recall"]
            row["shared_candidate_alt_informative_edges_dor_to_hcc"] = projection[
                "alt_informative_edges"
            ]["right_recall"]
            row["shared_candidate_aa_edges_dor_to_hcc"] = projection["aa_edges"][
                "right_recall"
            ]
            row["shared_candidate_co_membership_dor_to_hcc"] = projection[
                "co_membership_pairs"
            ]["right_recall"]
        pair_rows.append(row)

    require(len(pair_rows) == 21, "technical long pair table does not contain 21 rows")
    require(
        Counter(row["pair_class"] for row in pair_rows)
        == {
            "same_biological_id": 1,
            "same_cancer_different_id": 6,
            "cross_cancer": 14,
        },
        "technical pair-class census mismatch",
    )

    symmetric_matrices: dict[str, dict[str, dict[str, float | int | None]]] = {}
    directional_matrices: dict[str, dict[str, dict[str, float | int | None]]] = {}
    for metric in STRICT_METRICS:
        jaccard_matrix = blank_matrix(1.0)
        containment_matrix = blank_matrix(1.0)
        for row in pair_rows:
            left = row["dataset_a"]
            right = row["dataset_b"]
            jaccard = row[f"{metric}_jaccard"]
            jaccard_matrix[left][right] = jaccard
            jaccard_matrix[right][left] = jaccard
            containment_matrix[left][right] = row[f"{metric}_a_contained_in_b"]
            containment_matrix[right][left] = row[f"{metric}_b_contained_in_a"]
        symmetric_matrices[f"{metric}_jaccard"] = jaccard_matrix
        directional_matrices[f"{metric}_containment_directional"] = containment_matrix

    symmetric_pair_fields = (
        "coarse_topology_tvd",
        "coarse_topology_similarity",
        "coarse_topology_jsd_nats",
        "coarse_topology_js_distance",
        "coarse_topology_js_similarity",
        "coarse_topology_js_divergence_complement",
        "topology_signature_set_equal_rate",
        "topology_signature_set_compatible_rate",
        "both_unique_shape_equal_rate",
        "both_unique_k_ge_2_shape_equal_rate",
        "both_unique_ancestor_equal_rate",
        "exact_coordinate_one_to_one_n",
        "vaf_latest_jsd_nats",
        "vaf_latest_js_distance",
        "vaf_latest_js_similarity",
        "vaf_latest_js_divergence_complement",
        "vaf_baseline_jsd_nats",
        "vaf_baseline_js_distance",
        "vaf_baseline_js_similarity",
        "vaf_baseline_js_divergence_complement",
    )
    diagonal_by_field: dict[str, float | int | None] = {
        "coarse_topology_tvd": 0.0,
        "coarse_topology_similarity": 1.0,
        "coarse_topology_jsd_nats": 0.0,
        "coarse_topology_js_distance": 0.0,
        "coarse_topology_js_similarity": 1.0,
        "coarse_topology_js_divergence_complement": 1.0,
        "topology_signature_set_equal_rate": 1.0,
        "topology_signature_set_compatible_rate": 1.0,
        "both_unique_shape_equal_rate": 1.0,
        "both_unique_k_ge_2_shape_equal_rate": 1.0,
        "both_unique_ancestor_equal_rate": 1.0,
        "exact_coordinate_one_to_one_n": None,
        "vaf_latest_jsd_nats": 0.0,
        "vaf_latest_js_distance": 0.0,
        "vaf_latest_js_similarity": 1.0,
        "vaf_latest_js_divergence_complement": 1.0,
        "vaf_baseline_jsd_nats": 0.0,
        "vaf_baseline_js_distance": 0.0,
        "vaf_baseline_js_similarity": 1.0,
        "vaf_baseline_js_divergence_complement": 1.0,
    }
    for field in symmetric_pair_fields:
        matrix = blank_matrix(diagonal_by_field[field])
        if field == "exact_coordinate_one_to_one_n":
            for sample in DATASETS:
                self_diag = exact_coordinate_diagnostics(samples[sample], samples[sample])
                matrix[sample][sample] = self_diag[field]
        for row in pair_rows:
            left = row["dataset_a"]
            right = row["dataset_b"]
            matrix[left][right] = row[field]
            matrix[right][left] = row[field]
        symmetric_matrices[field] = matrix

    for label, matrix in symmetric_matrices.items():
        if label == "exact_coordinate_one_to_one_n":
            assert_symmetric_count_matrix(matrix, label)
        else:
            assert_symmetric_matrix(
                matrix,
                label,
                diagonal=diagonal_by_field.get(label, 1.0),
            )
    for label, matrix in directional_matrices.items():
        assert_directional_diagonal(matrix, label)

    numeric_summary_fields = [
        *(f"{metric}_jaccard" for metric in STRICT_METRICS),
        "coarse_topology_similarity",
        "coarse_topology_tvd",
        "coarse_topology_jsd_nats",
        "coarse_topology_js_distance",
        "coarse_topology_js_similarity",
        "coarse_topology_js_divergence_complement",
        "topology_signature_set_equal_rate",
        "topology_signature_set_compatible_rate",
        "both_unique_shape_equal_rate",
        "both_unique_k_ge_2_shape_equal_rate",
        "both_unique_ancestor_equal_rate",
        "vaf_latest_jsd_nats",
        "vaf_latest_js_distance",
        "vaf_latest_js_similarity",
        "vaf_latest_js_divergence_complement",
        "vaf_latest_ks",
        "vaf_latest_wasserstein",
        "vaf_baseline_jsd_nats",
        "vaf_baseline_js_distance",
        "vaf_baseline_js_similarity",
        "vaf_baseline_js_divergence_complement",
        "vaf_baseline_ks",
        "vaf_baseline_wasserstein",
    ]
    technical_class_summary = aggregate_rows(pair_rows, ("pair_class",), numeric_summary_fields)

    biological_groups: dict[tuple[str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for row in pair_rows:
        left_bio = row["biological_id_a"]
        right_bio = row["biological_id_b"]
        if left_bio == right_bio:
            continue
        biological_groups[biological_pair(left_bio, right_bio)].append(row)

    biological_pair_rows: list[dict[str, Any]] = []
    for bio_pair in sorted(biological_groups):
        members = biological_groups[bio_pair]
        left_bio, right_bio = bio_pair
        representative_left = next(
            sample
            for sample in DATASETS
            if SAMPLE_METADATA[sample]["biological_id"] == left_bio
        )
        representative_right = next(
            sample
            for sample in DATASETS
            if SAMPLE_METADATA[sample]["biological_id"] == right_bio
        )
        record: dict[str, Any] = {
            "biological_id_a": left_bio,
            "biological_id_b": right_bio,
            "cancer_type_a": SAMPLE_METADATA[representative_left]["cancer_type"],
            "cancer_type_b": SAMPLE_METADATA[representative_right]["cancer_type"],
            "pair_class": pair_class(representative_left, representative_right),
            "technical_pair_n": len(members),
            "technical_pairs": ";".join(
                f"{row['dataset_a']}--{row['dataset_b']}" for row in members
            ),
        }
        for field in numeric_summary_fields:
            mean, minimum, maximum = mean_min_max(float_or_none(row.get(field)) for row in members)
            record[f"{field}_mean"] = mean
            record[f"{field}_min"] = minimum
            record[f"{field}_max"] = maximum
        biological_pair_rows.append(record)
    require(len(biological_pair_rows) == 15, "biological pair table is not 6 choose 2")
    require(
        Counter(row["pair_class"] for row in biological_pair_rows)
        == {"same_cancer_different_id": 4, "cross_cancer": 11},
        "biological pair-class census mismatch",
    )

    biological_summary_fields = [f"{field}_mean" for field in numeric_summary_fields]
    biological_class_summary = aggregate_rows(
        biological_pair_rows,
        ("pair_class",),
        biological_summary_fields,
    )
    same_cancer_biological_pairs = [
        row
        for row in biological_pair_rows
        if row["pair_class"] == "same_cancer_different_id"
    ]
    same_cancer_summary = aggregate_rows(
        same_cancer_biological_pairs,
        (),
        biological_summary_fields,
    )
    require(len(same_cancer_biological_pairs) == 4, "deduplicated same-cancer pairs != 4")
    require(
        len(same_cancer_summary) == 1 and same_cancer_summary[0]["pair_n"] == 4,
        "deduplicated same-cancer summary denominator mismatch",
    )

    exact_permutation_metric_fields = (
        "candidate_sites_jaccard",
        "active_sites_jaccard",
        "w_sites_jaccard",
        "primary_edges_jaccard",
        "exact_components_jaccard",
        "coarse_topology_similarity",
        "coarse_topology_js_similarity",
        "vaf_latest_js_similarity",
        "vaf_baseline_js_similarity",
    )
    exact_permutation_tests, exact_permutation_distribution = (
        exact_group_label_permutation(
            biological_pair_rows,
            exact_permutation_metric_fields,
        )
    )

    target_rows = [row for row in pair_rows if row["is_hcc1395_technical_pair"]]
    require(len(target_rows) == 1, "HCC1395 technical target row is not unique")
    target = target_rows[0]
    target_containment_rows = [
        {
            "source_dataset": "HCC1395_DORADO",
            "target_dataset": "HCC1395",
            "projection_universe": "shared exact candidate allele keys",
            "metric": metric,
            "source_n": strict["target_pair"]["shared_candidate_projection"][metric][
                "right_n"
            ],
            "target_n": strict["target_pair"]["shared_candidate_projection"][metric][
                "left_n"
            ],
            "intersection_n": strict["target_pair"]["shared_candidate_projection"][metric][
                "intersection_n"
            ],
            "source_contained_in_target": strict["target_pair"][
                "shared_candidate_projection"
            ][metric]["right_recall"],
        }
        for metric in (
            "active_sites",
            "w_sites",
            "primary_edges",
            "alt_informative_edges",
            "aa_edges",
            "co_membership_pairs",
        )
    ]

    measure_definitions = {
        "schema_name": "intersubmod.crosssample_topology_matrix.measure_definitions",
        "schema_version": "1.1.0",
        "dataset_grain": "7 technical datasets / 6 biological IDs / GRCh38 chr1-22",
        "pair_class": {
            "same_biological_id": "same biological ID; technical-source comparison",
            "same_cancer_different_id": "different biological IDs with the same cancer type",
            "cross_cancer": "different cancer types",
        },
        "strict_set_metrics": strict["definitions"],
        "directional_containment": (
            "Matrix[row, column] = |row set intersect column set| / |row set|. "
            "It is directional and is not required to be symmetric."
        ),
        "coarse_topology_composition": {
            "categories": list(COARSE_CATEGORIES),
            "resolved_rule": (
                "A census unit is assigned to a coarse class only when all AF-optimal "
                "tree signatures occupy exactly one coarse class; otherwise it is "
                "Unresolved-cross-coarse."
            ),
            "tvd": "0.5 * sum absolute category-proportion differences",
            "similarity": "1 - TVD",
            "jsd_nats": "Jensen-Shannon divergence on the five category proportions",
            "js_distance": "sqrt(JSD/ln(2)); normalized Jensen-Shannon distance",
            "js_similarity": "1 - sqrt(JSD/ln(2))",
            "js_divergence_complement": (
                "1 - JSD/ln(2); retained under an explicit name and not "
                "interpreted as 1 minus Jensen-Shannon distance"
            ),
        },
        "exact_coordinate_diagnostics": {
            "join_key_within_sample": "(unit_id, group_index)",
            "coordinate_signature": "(chrom, sorted active_positions); coordinate-only, not REF/ALT",
            "fail_closed_1_to_1": (
                "A shared coordinate signature is evaluable only if exactly one ranked "
                "unit exists in each dataset. Duplicate signature keys are excluded."
            ),
            "topology_signature_equality": "the two complete unlabeled rooted-shape sets are equal",
            "topology_signature_compatibility": "the two complete rooted-shape sets have a non-empty intersection",
            "both_unique_shape": "both census rows are UNIQUE_TREE; shape equality compares their sole shape",
            "ancestor_relation": (
                "Strict mutation-position ancestry derived from the verified unique, "
                "recurrence-free representative arborescence. Primary denominator is k>=2; "
                "k=1 empty relations are excluded."
            ),
        },
        "vaf_distribution": {
            "scope": "truth_confirmed_marginal_raw_vaf",
            "sources": {
                "latest_lps_pass": "latest LongPhase-S PASS universe",
                "baseline_caller_pass": "caller-PASS baseline universe",
            },
            "js_divergence_50bin_nats": (
                "Jensen-Shannon divergence between 50-bin marginal raw-VAF histograms, "
                "using natural logarithms; bounded by ln(2)."
            ),
            "js_distance": (
                "sqrt(js_divergence_50bin_nats / ln(2)); normalized "
                "Jensen-Shannon distance"
            ),
            "js_similarity": (
                "1 - sqrt(js_divergence_50bin_nats / ln(2))"
            ),
            "js_divergence_complement": (
                "1 - js_divergence_50bin_nats / ln(2); explicitly not "
                "1 minus Jensen-Shannon distance"
            ),
            "ks_statistic": "two-sample Kolmogorov-Smirnov statistic",
            "wasserstein_vaf": "1D Wasserstein distance in raw-VAF units",
        },
        "exact_group_label_permutation": {
            "unit": "six biological IDs after technical-source macro-averaging",
            "observed_group_sizes": "breast=3, lung_adenocarcinoma=2, melanoma=1",
            "assignment_count": "6!/(3!*2!*1!) = 60 exhaustive assignments",
            "statistic": (
                "mean similarity over four within-group biological-ID pairs "
                "minus mean similarity over eleven between-group pairs"
            ),
            "one_sided_p": (
                "fraction of all 60 assignments with statistic greater than "
                "or equal to the observed statistic"
            ),
            "technical_source_handling": (
                "For every biological pair involving HCC1395, HCC1395 and "
                "HCC1395_DORADO technical-pair values are macro-averaged first."
            ),
        },
        "claim_ceiling": (
            "Exact-coordinate comparisons evaluate mathematical mutation-state tree shapes "
            "and unique-tree ancestor relations on 1:1 coordinate signatures. They do not "
            "validate cellular clones, biological ancestry, or a global clone tree."
        ),
    }

    output_dir = args.output_dir
    matrix_dir = output_dir / "matrices"
    output_dir.mkdir(parents=True, exist_ok=False)
    matrix_dir.mkdir(parents=True, exist_ok=False)

    pair_columns = list(pair_rows[0].keys())
    write_tsv_exclusive(output_dir / "technical_pair_metrics.tsv", pair_rows, pair_columns)
    write_tsv_exclusive(
        output_dir / "sample_topology_composition.tsv",
        sample_composition_rows,
        list(sample_composition_rows[0].keys()),
    )
    write_tsv_exclusive(
        output_dir / "technical_pair_class_summary.tsv",
        technical_class_summary,
        list(technical_class_summary[0].keys()),
    )
    write_tsv_exclusive(
        output_dir / "biological_pair_aggregates.tsv",
        biological_pair_rows,
        list(biological_pair_rows[0].keys()),
    )
    write_tsv_exclusive(
        output_dir / "biological_pair_class_summary.tsv",
        biological_class_summary,
        list(biological_class_summary[0].keys()),
    )
    write_tsv_exclusive(
        output_dir / "same_cancer_biological_pairs.tsv",
        same_cancer_biological_pairs,
        list(same_cancer_biological_pairs[0].keys()),
    )
    write_tsv_exclusive(
        output_dir / "same_cancer_biological_summary.tsv",
        same_cancer_summary,
        list(same_cancer_summary[0].keys()),
    )
    write_tsv_exclusive(
        output_dir / "target_shared_candidate_containment.tsv",
        target_containment_rows,
        list(target_containment_rows[0].keys()),
    )
    write_tsv_exclusive(
        output_dir / "biological_id_exact_permutation_tests.tsv",
        exact_permutation_tests,
        list(exact_permutation_tests[0].keys()),
    )
    write_tsv_exclusive(
        output_dir / "biological_id_exact_permutation_distribution.tsv",
        exact_permutation_distribution,
        list(exact_permutation_distribution[0].keys()),
    )
    write_json_exclusive(output_dir / "measure_definitions.json", measure_definitions)

    for label, matrix in sorted(symmetric_matrices.items()):
        write_matrix_exclusive(matrix_dir / f"{label}.tsv", matrix)
    for label, matrix in sorted(directional_matrices.items()):
        write_matrix_exclusive(matrix_dir / f"{label}.tsv", matrix)

    summary_payload = {
        "schema_name": "intersubmod.crosssample_topology_matrix.summary",
        "schema_version": "1.1.0",
        "scope": {
            "datasets": list(DATASETS),
            "technical_pairs": 21,
            "biological_ids": 6,
            "biological_pairs": 15,
            "chromosomes": "GRCh38 chr1-22",
        },
        "pair_class_census": dict(Counter(row["pair_class"] for row in pair_rows)),
        "biological_pair_class_census": dict(
            Counter(row["pair_class"] for row in biological_pair_rows)
        ),
        "hcc1395_target": target,
        "same_cancer_biological_id_deduplicated_pair_n": len(
            same_cancer_biological_pairs
        ),
        "exact_group_label_permutation": {
            "assignment_n": 60,
            "tests": exact_permutation_tests,
        },
        "claim_ceiling": measure_definitions["claim_ceiling"],
    }
    write_json_exclusive(output_dir / "analysis_summary.json", summary_payload)

    source_identities_after = {label: identity(path) for label, path in source_paths.items()}
    require(
        source_identities_before == source_identities_after,
        "one or more source files changed during analysis",
    )

    output_paths = sorted(
        path
        for path in output_dir.rglob("*")
        if path.is_file() and path.name != "validation_receipt.json"
    )
    output_identities = {
        str(path.relative_to(output_dir)): identity(path) for path in output_paths
    }
    checks = {
        "all_pass": True,
        "dataset_count_7": len(DATASETS) == 7,
        "technical_pair_count_21": len(pair_rows) == 21,
        "biological_id_count_6": len(
            {metadata["biological_id"] for metadata in SAMPLE_METADATA.values()}
        )
        == 6,
        "biological_pair_count_15": len(biological_pair_rows) == 15,
        "strict_pair_keys_complete_and_unique": set(strict_pairs) == expected_pairs,
        "topology_census_join_1_to_1": all(
            data.topology_join_key_duplicates == 0
            and data.census_join_key_duplicates == 0
            and data.ranked_rows == data.census_rows
            for data in samples.values()
        ),
        "duplicate_coordinate_signatures_excluded_fail_closed": all(
            row["exact_coordinate_one_to_one_n"]
            + row["exact_coordinate_duplicate_excluded_signature_n"]
            == row["exact_coordinate_shared_signature_n"]
            for row in pair_rows
        ),
        "all_symmetric_matrices_pass": True,
        "all_directional_diagonals_equal_1": True,
        "technical_pair_class_census_1_6_14": Counter(
            row["pair_class"] for row in pair_rows
        )
        == {
            "same_biological_id": 1,
            "same_cancer_different_id": 6,
            "cross_cancer": 14,
        },
        "biological_pair_class_census_4_11": Counter(
            row["pair_class"] for row in biological_pair_rows
        )
        == {"same_cancer_different_id": 4, "cross_cancer": 11},
        "same_cancer_biological_pairs_deduplicated_4": len(
            same_cancer_biological_pairs
        )
        == 4,
        "exact_label_permutation_assignment_count_60": len(
            {row["assignment_id"] for row in exact_permutation_distribution}
        )
        == 60,
        "exact_label_permutation_observed_assignment_unique": sum(
            row["is_observed_cancer_assignment"]
            for row in exact_permutation_distribution
            if row["metric"] == exact_permutation_metric_fields[0]
        )
        == 1,
        "exact_label_permutation_metric_coverage": len(exact_permutation_tests)
        == len(exact_permutation_metric_fields)
        and all(
            row["permutation_n"] == 60
            and 0.0 < row["one_sided_exact_p"] <= 1.0
            for row in exact_permutation_tests
        ),
        "cancer_mapping_hcc1937_breast": SAMPLE_METADATA["HCC1937"]["cancer_type"]
        == "breast",
        "source_sha256_stable_before_after": source_identities_before
        == source_identities_after,
        "hcc_shared_projection_unique": len(target_rows) == 1
        and bool(target["shared_candidate_projection_available"]),
    }
    require(all(checks.values()), f"validation checks failed: {checks}")
    receipt = {
        "schema_name": "intersubmod.crosssample_topology_matrix.validation_receipt",
        "schema_version": "1.1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "serves_goals": ["G4", "G5"],
        "checks": checks,
        "scope": summary_payload["scope"],
        "pair_class_census": summary_payload["pair_class_census"],
        "biological_pair_class_census": summary_payload[
            "biological_pair_class_census"
        ],
        "sample_metadata": SAMPLE_METADATA,
        "sources": source_identities_before,
        "outputs": output_identities,
        "claim_ceiling": measure_definitions["claim_ceiling"],
    }
    write_json_exclusive(output_dir / "validation_receipt.json", receipt)

    print(
        json.dumps(
            {
                "status": "PASS",
                "technical_pairs": len(pair_rows),
                "biological_pairs": len(biological_pair_rows),
                "matrices": len(symmetric_matrices) + len(directional_matrices),
                "output_dir": str(output_dir.resolve()),
                "exact_label_permutation": {
                    row["metric"]: {
                        "within_minus_between": row[
                            "observed_within_minus_between"
                        ],
                        "one_sided_exact_p": row["one_sided_exact_p"],
                    }
                    for row in exact_permutation_tests
                },
                "hcc_target": {
                    "candidate_jaccard": target["candidate_sites_jaccard"],
                    "primary_edge_jaccard": target["primary_edges_jaccard"],
                    "shared_primary_edge_dor_to_hcc": target[
                        "shared_candidate_primary_edges_dor_to_hcc"
                    ],
                    "coarse_topology_similarity": target[
                        "coarse_topology_similarity"
                    ],
                    "one_to_one_coordinate_signatures": target[
                        "exact_coordinate_one_to_one_n"
                    ],
                    "topology_signature_compatible_rate": target[
                        "topology_signature_set_compatible_rate"
                    ],
                    "both_unique_ancestor_equal_rate": target[
                        "both_unique_ancestor_equal_rate"
                    ],
                    "vaf_latest_js_similarity": target[
                        "vaf_latest_js_similarity"
                    ],
                },
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AnalysisError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
