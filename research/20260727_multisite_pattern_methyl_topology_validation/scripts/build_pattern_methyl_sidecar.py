#!/usr/bin/env python3
"""Build a hash-bound, overlay-only pattern-methylation sidecar.

The sidecar joins the merged exact-PS x exact-raw-HP methylation evidence to
seven frozen topology JSONL authorities.  It emits references for complete
R/A states and association overlays, but never rewrites topology counts, AF,
edge incidence, or the selected representative tree.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import sys
import tempfile
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA_NAME = "intersubmod.pattern_methyl_sidecar"
SCHEMA_VERSION = "1.0.0"
RECEIPT_SCHEMA_NAME = f"{SCHEMA_NAME}.receipt"
EVIDENCE_SCHEMA_NAME = "intersubmod.pattern_methyl_evidence"
DETAIL_SCHEMA_NAME = f"{EVIDENCE_SCHEMA_NAME}.detail"
TOPOLOGY_SCHEMA_NAME = "intersubmod.exact_ps_cpp_topology_af.unit"
OUTPUT_FILENAME = "pattern_methyl_sidecar.v1.json.gz"
RECEIPT_FILENAME = "pattern_methyl_sidecar.v1.receipt.json"
FAILED_STAGING_DIRECTORY = "_failed_staging_archive"
TOPOLOGY_FILE_COUNT = 7
DEFAULT_MATRIX_TOP_N = 24
SHA_CHUNK_BYTES = 8 * 1024 * 1024

HP_FAMILY_BY_RAW = {
    ".": "none",
    "1": "1",
    "1-1": "1",
    "1-2": "1",
    "2": "2",
    "2-1": "2",
    "2-2": "2",
    "3": "3",
    "4": "4",
}
HP_ORDER = {
    value: index
    for index, value in enumerate(
        (".", "1", "1-1", "1-2", "2", "2-1", "2-2", "3", "4")
    )
}
ASSESSMENT_ORDER = {
    "ROBUST_ASSOCIATION": 0,
    "LOCAL_CIS_COMPATIBLE": 1,
    "TAG_DEPENDENT": 2,
    "CONFOUNDED": 3,
    "EVALUABLE_NO_ROBUST_ASSOCIATION": 4,
    "NOT_EVALUABLE": 5,
}
ASSESSMENTS = frozenset(ASSESSMENT_ORDER)
EVALUATION_STATUSES = frozenset({"EVALUABLE", "NOT_EVALUABLE"})
HALO_RELATION = "HAMMING1_GLOBAL_BEST_UNANIMOUS"
HAMMING_GT1_RELATION = "PAIR_BAND_ONLY_HAMMING_GT1"
PARTIAL_X_SOURCE_STATUS = "UNAVAILABLE_NO_HASH_BOUND_PATTERN_COUNTS"
PAIR_EFFECT_FIELDS = (
    "hamming",
    "between_mean",
    "pooled_within_mean",
    "distance_contrast",
    "standardized_effect",
    "topology_relation",
)

EVIDENCE_FIELDS = (
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

INTEGER_FIELDS = frozenset(
    {
        "n_active_bits",
        "input_n_complete",
        "analysis_n",
        "n_common_cpg",
        "n_distal_cpg",
        "tile_overlap_conflicts",
        "exchangeable_strata",
        "exchangeable_n",
        "permanova_permutations_requested",
        "permanova_permutations_realized",
        "best_pair_hamming",
        "distal_permutations_realized",
    }
)
FLOAT_FIELDS = frozenset(
    {
        "qname_join_fraction_min",
        "permanova_pseudo_f",
        "permanova_r2",
        "permanova_p",
        "permdisp_f",
        "permdisp_p",
        "best_pair_between_mean",
        "best_pair_pooled_within_mean",
        "best_pair_distance_contrast",
        "best_pair_standardized_effect",
        "max_geometry_smd",
        "equal_n_r2",
        "equal_n_retention",
        "rarefaction_median_r2",
        "rarefaction_retention",
        "distal_r2",
        "distal_p",
        "distal_retention",
        "q_by",
        "p_holm",
    }
)
BOOLEAN_FIELDS = frozenset({"pair_full4", "k_ge_3", "all_states_n8", "all_states_n10"})
IDENTITY_FIELDS = frozenset(
    {
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
        "input_state_counts_json",
        "analysis_state_counts_json",
    }
)
TOPOLOGY_SNAPSHOT_FIELDS = (
    "sample",
    "chrom",
    "unit_id",
    "phase_set",
    "hp_family",
    "active_positions",
    "best_tree_unique",
    "representative_best_edges",
    "representative_best_vertices",
    "best_tree_tie_count",
)


class SidecarContractError(RuntimeError):
    """Raised when an input cannot be safely promoted into the sidecar."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SidecarContractError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(SHA_CHUNK_BYTES)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def file_binding(path: Path | str) -> dict[str, Any]:
    candidate = Path(path)
    try:
        resolved = candidate.resolve(strict=True)
    except FileNotFoundError as exc:
        raise SidecarContractError(f"missing input: {candidate}") from exc
    require(resolved.is_file(), f"input is not a regular file: {resolved}")
    size = resolved.stat().st_size
    require(size > 0, f"input is empty: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": size,
        "sha256": sha256_file(resolved),
    }


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def reject_json_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant: {value}")


def load_json_line(text: str, label: str) -> Any:
    try:
        return json.loads(text, parse_constant=reject_json_constant)
    except (json.JSONDecodeError, ValueError) as exc:
        raise SidecarContractError(f"invalid JSON at {label}: {exc}") from exc


def validate_json_value(value: Any, label: str) -> Any:
    """Validate JSON compatibility while preserving explicit null values."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        require(math.isfinite(value), f"non-finite value at {label}")
        return value
    if isinstance(value, list):
        return [
            validate_json_value(item, f"{label}[{index}]")
            for index, item in enumerate(value)
        ]
    if isinstance(value, dict):
        require(
            all(isinstance(key, str) for key in value),
            f"non-string JSON key at {label}",
        )
        return {
            key: validate_json_value(item, f"{label}.{key}")
            for key, item in value.items()
        }
    raise SidecarContractError(f"unsupported JSON value at {label}: {type(value).__name__}")


def natural_chrom_key(chrom: str) -> tuple[int, str]:
    if chrom.startswith("chr") and chrom[3:].isdigit():
        return int(chrom[3:]), chrom
    return 10**9, chrom


def evidence_key(row: Mapping[str, Any]) -> tuple[str, str, str, str]:
    return (
        str(row["dataset"]),
        str(row["chrom"]),
        str(row["region_id"]),
        str(row["hp_raw"]),
    )


def topology_key(row: Mapping[str, Any]) -> tuple[str, str]:
    return str(row["sample"]), str(row["region_id"])


def optional_integer(value: Any, label: str) -> int | None:
    if value in {"", None}:
        return None
    try:
        result = int(str(value))
    except ValueError as exc:
        raise SidecarContractError(f"invalid integer {label}: {value!r}") from exc
    require(result >= 0, f"negative integer {label}: {result}")
    return result


def optional_float(value: Any, label: str) -> float | None:
    if value in {"", None, "NA", "NaN", "nan"}:
        return None
    try:
        result = float(str(value))
    except ValueError as exc:
        raise SidecarContractError(f"invalid float {label}: {value!r}") from exc
    require(math.isfinite(result), f"non-finite float {label}: {value!r}")
    return result


def optional_boolean(value: Any, label: str) -> bool | None:
    if value in {"", None}:
        return None
    text = str(value).strip().lower()
    require(text in {"true", "false"}, f"invalid boolean {label}: {value!r}")
    return text == "true"


def parse_positions(value: Any, label: str) -> tuple[int, ...]:
    text = str(value or "").strip()
    require(text != "", f"missing positions at {label}")
    try:
        if text.startswith("["):
            raw = load_json_line(text, label)
            require(isinstance(raw, list), f"positions are not a list at {label}")
            positions = tuple(int(item) for item in raw)
        else:
            positions = tuple(int(item) for item in text.split(",") if item)
    except (TypeError, ValueError) as exc:
        raise SidecarContractError(f"invalid positions at {label}: {text!r}") from exc
    require(positions, f"empty positions at {label}")
    require(
        all(position > 0 for position in positions),
        f"non-positive position at {label}",
    )
    require(
        tuple(sorted(set(positions))) == positions,
        f"positions are not strictly increasing at {label}",
    )
    return positions


def validate_pattern(pattern: str, n_active_bits: int, label: str) -> None:
    require(
        len(pattern) == n_active_bits,
        f"pattern length mismatch at {label}: {pattern!r}",
    )
    require(
        set(pattern) <= {"R", "A", "X"},
        f"pattern contains non-R/A/X code at {label}: {pattern!r}",
    )


def parse_state_counts(
    value: Any,
    *,
    n_active_bits: int,
    label: str,
) -> dict[str, int]:
    if value in {"", None}:
        return {}
    payload = load_json_line(str(value), label)
    require(isinstance(payload, dict), f"state counts are not an object at {label}")
    counts: dict[str, int] = {}
    for raw_pattern, raw_count in payload.items():
        require(isinstance(raw_pattern, str), f"non-string pattern at {label}")
        validate_pattern(raw_pattern, n_active_bits, label)
        require(
            isinstance(raw_count, int)
            and not isinstance(raw_count, bool)
            and raw_count >= 0,
            f"invalid state count at {label}.{raw_pattern}: {raw_count!r}",
        )
        counts[raw_pattern] = raw_count
    return dict(sorted(counts.items()))


def complete_count_sum(counts: Mapping[str, int]) -> int:
    return sum(count for pattern, count in counts.items() if "X" not in pattern)


def parse_evidence_row(
    raw: Mapping[str, Any], path: Path, line_number: int
) -> dict[str, Any]:
    label = f"{path}:{line_number}"
    require(None not in raw, f"extra TSV cells at {label}")
    require(
        raw.get("schema_version") == SCHEMA_VERSION,
        f"unsupported evidence schema_version at {label}: {raw.get('schema_version')!r}",
    )
    output: dict[str, Any] = {}
    for field in EVIDENCE_FIELDS:
        value = raw.get(field)
        field_label = f"{label}.{field}"
        if field in INTEGER_FIELDS:
            output[field] = optional_integer(value, field_label)
        elif field in FLOAT_FIELDS:
            output[field] = optional_float(value, field_label)
        elif field in BOOLEAN_FIELDS:
            output[field] = optional_boolean(value, field_label)
        else:
            output[field] = None if value in {"", None} else str(value)

    for field in (
        "dataset",
        "chrom",
        "region_id",
        "unit_id",
        "phase_set",
        "hp_family",
        "hp_raw",
        "assessment",
        "evaluation_status",
    ):
        require(output[field] is not None, f"missing required field {field} at {label}")
    require(
        output["chrom"].startswith("chr") and output["chrom"][3:].isdigit(),
        f"invalid chromosome at {label}: {output['chrom']!r}",
    )
    hp_raw = str(output["hp_raw"])
    require(hp_raw in HP_FAMILY_BY_RAW, f"unknown exact raw HP at {label}: {hp_raw!r}")
    require(
        output["hp_family"] == HP_FAMILY_BY_RAW[hp_raw],
        f"raw HP/family mismatch at {label}: {hp_raw!r}/{output['hp_family']!r}",
    )
    require(
        output["assessment"] in ASSESSMENTS,
        f"invalid assessment at {label}: {output['assessment']!r}",
    )
    require(
        output["evaluation_status"] in EVALUATION_STATUSES,
        f"invalid evaluation_status at {label}: {output['evaluation_status']!r}",
    )

    positions = parse_positions(raw.get("active_positions"), f"{label}.active_positions")
    output["active_positions"] = list(positions)
    require(
        output["n_active_bits"] == len(positions),
        f"active-bit count mismatch at {label}",
    )
    require(
        output["n_active_bits"] is not None and output["n_active_bits"] >= 2,
        f"evidence unit has fewer than two active bits at {label}",
    )
    input_counts = parse_state_counts(
        raw.get("input_state_counts_json"),
        n_active_bits=int(output["n_active_bits"]),
        label=f"{label}.input_state_counts_json",
    )
    analysis_counts = parse_state_counts(
        raw.get("analysis_state_counts_json"),
        n_active_bits=int(output["n_active_bits"]),
        label=f"{label}.analysis_state_counts_json",
    )
    output["input_state_counts"] = input_counts
    output["analysis_state_counts"] = analysis_counts
    require(
        all("X" not in pattern for pattern in input_counts),
        f"input_state_counts contains X but analyzer evidence is complete-R/A-only at {label}",
    )
    require(
        all("X" not in pattern for pattern in analysis_counts),
        f"analysis_state_counts contains X at {label}",
    )
    require(
        output["input_n_complete"] is not None
        and complete_count_sum(input_counts) == output["input_n_complete"],
        f"input complete-count conservation failed at {label}",
    )
    if output["analysis_n"] is None:
        require(
            not analysis_counts,
            f"analysis state counts exist without analysis_n at {label}",
        )
    else:
        require(
            complete_count_sum(analysis_counts) == output["analysis_n"],
            f"analysis count conservation failed at {label}",
        )
    if output["evaluation_status"] == "EVALUABLE":
        require(
            output["analysis_n"] is not None and bool(analysis_counts),
            f"evaluable row lacks analysis counts at {label}",
        )
        require(
            len(analysis_counts) >= 2,
            f"evaluable row has fewer than two analysis states at {label}",
        )

    if output["best_pair"] is None:
        require(
            all(
                output[field] is None
                for field in (
                    "best_pair_hamming",
                    "best_pair_between_mean",
                    "best_pair_pooled_within_mean",
                    "best_pair_distance_contrast",
                    "best_pair_standardized_effect",
                    "best_pair_topology_relation",
                )
            ),
            f"best-pair metadata exists without best_pair at {label}",
        )
        require(
            output["evaluation_status"] != "EVALUABLE",
            f"evaluable row lacks best_pair at {label}",
        )
    else:
        parts = str(output["best_pair"]).split("|")
        require(len(parts) == 2 and parts[0] != parts[1], f"invalid best_pair at {label}")
        for pattern in parts:
            validate_pattern(pattern, int(output["n_active_bits"]), f"{label}.best_pair")
            require("X" not in pattern, f"best_pair contains X at {label}")
            require(
                pattern in analysis_counts,
                f"best_pair references absent analysis state at {label}: {pattern}",
            )
        require(
            output["best_pair_hamming"] == hamming_distance(parts[0], parts[1]),
            f"best_pair Hamming mismatch at {label}",
        )
        require(
            all(
                output[field] is not None
                for field in (
                    "best_pair_between_mean",
                    "best_pair_pooled_within_mean",
                    "best_pair_distance_contrast",
                    "best_pair_standardized_effect",
                    "best_pair_topology_relation",
                )
            ),
            f"best_pair has incomplete metrics at {label}",
        )
    return output


def load_evidence(
    path: Path | str,
) -> tuple[dict[tuple[str, str, str, str], dict[str, Any]], dict[str, Any]]:
    binding = file_binding(path)
    resolved = Path(binding["path"])
    rows: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    try:
        with gzip.open(resolved, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            require(
                tuple(reader.fieldnames or ()) == EVIDENCE_FIELDS,
                f"evidence TSV header drift: {resolved}",
            )
            for line_number, raw in enumerate(reader, start=2):
                parsed = parse_evidence_row(raw, resolved, line_number)
                key = evidence_key(parsed)
                require(key not in rows, f"duplicate exact-raw-HP evidence key: {key}")
                rows[key] = parsed
    except (OSError, UnicodeError, csv.Error) as exc:
        raise SidecarContractError(f"failed to read evidence TSV {resolved}: {exc}") from exc
    require(rows, f"evidence TSV contains no rows: {resolved}")
    binding["row_count"] = len(rows)
    binding["schema_name"] = EVIDENCE_SCHEMA_NAME
    binding["schema_version"] = SCHEMA_VERSION
    return rows, binding


def validate_topology_row(
    row: Mapping[str, Any],
    label: str,
    *,
    require_evidence_linked_ranked: bool = False,
) -> None:
    require(
        row.get("schema_name") == TOPOLOGY_SCHEMA_NAME
        and row.get("schema_version") == SCHEMA_VERSION,
        f"unexpected topology schema at {label}",
    )
    for field in (
        "sample",
        "chrom",
        "region_id",
        "unit_id",
        "phase_set",
        "hp_family",
    ):
        require(
            isinstance(row.get(field), str) and bool(row.get(field)),
            f"missing topology field {field} at {label}",
        )
    positions = row.get("active_positions")
    require(
        isinstance(positions, list)
        and all(isinstance(value, int) and value > 0 for value in positions)
        and positions == sorted(set(positions)),
        f"invalid topology active_positions at {label}",
    )
    require(
        isinstance(row.get("active_bit_count"), int)
        and not isinstance(row.get("active_bit_count"), bool)
        and row.get("active_bit_count") == len(positions),
        f"topology active_bit_count mismatch at {label}",
    )
    if not require_evidence_linked_ranked:
        return
    require(
        len(positions) >= 2,
        f"evidence-linked topology requires at least two active positions at {label}",
    )
    best_tree_unique = row.get("best_tree_unique")
    if best_tree_unique is None:
        require(
            isinstance(row.get("unit_status"), str)
            and row.get("unit_status") not in {"", "ranked"},
            f"evidence-linked topology lacks a selected tree without an unranked status at {label}",
        )
        require(
            row.get("representative_best_vertices") in (None, []),
            f"unranked topology unexpectedly contains selected vertices at {label}",
        )
        require(
            row.get("representative_best_edges") == [],
            f"unranked topology unexpectedly contains selected edges at {label}",
        )
        return
    require(
        isinstance(best_tree_unique, bool),
        f"topology best_tree_unique is not boolean at {label}",
    )
    vertices = row.get("representative_best_vertices")
    edges = row.get("representative_best_edges")
    require(isinstance(vertices, list), f"topology vertices are not a list at {label}")
    require(isinstance(edges, list), f"topology edges are not a list at {label}")
    vertex_ids: set[int] = set()
    for index, vertex in enumerate(vertices):
        require(isinstance(vertex, dict), f"invalid topology vertex at {label}[{index}]")
        value = vertex.get("vertex")
        require(
            isinstance(value, int) and not isinstance(value, bool) and value >= 0,
            f"invalid topology vertex id at {label}[{index}]",
        )
        require(value not in vertex_ids, f"duplicate topology vertex at {label}: {value}")
        vertex_ids.add(value)
    edge_ids: set[tuple[int, int]] = set()
    for index, edge in enumerate(edges):
        require(isinstance(edge, dict), f"invalid topology edge at {label}[{index}]")
        parent = edge.get("parent_vertex")
        child = edge.get("child_vertex")
        require(
            isinstance(parent, int)
            and not isinstance(parent, bool)
            and isinstance(child, int)
            and not isinstance(child, bool),
            f"non-integer topology edge endpoint at {label}[{index}]",
        )
        delta = parent ^ child
        require(
            parent >= 0
            and child >= 0
            and parent & ~child == 0
            and delta > 0
            and delta & (delta - 1) == 0,
            f"topology edge is not one-bit acquisition at {label}: {parent}->{child}",
        )
        require(
            parent in vertex_ids and child in vertex_ids,
            f"topology edge endpoint absent from selected vertices at {label}",
        )
        require(
            (parent, child) not in edge_ids,
            f"duplicate topology edge at {label}: {parent}->{child}",
        )
        edge_ids.add((parent, child))


def topology_protected_hashes(row: Mapping[str, Any]) -> dict[str, str]:
    counts = {
        key: value
        for key, value in row.items()
        if "count" in key
        or key
        in {
            "active_bit_count",
            "original_bit_count",
            "objective_h",
            "search_nodes",
            "total_tree_count",
        }
    }
    af = {
        "af_basis": row.get("af_basis"),
        "af_coverage": row.get("af_coverage"),
        "best_score_fraction": row.get("best_score_fraction"),
        "input_vaf_eligible": row.get("input_vaf_eligible"),
    }
    selected = {
        "best_tree_unique": row.get("best_tree_unique"),
        "representative_best_edges": row.get("representative_best_edges"),
        "representative_best_vertices": row.get("representative_best_vertices"),
        "representative_best_morphology": row.get("representative_best_morphology"),
    }
    protected = {"counts": counts, "af": af, "selected_exemplar": selected}
    return {
        "counts_sha256": canonical_sha256(counts),
        "af_sha256": canonical_sha256(af),
        "selected_exemplar_sha256": canonical_sha256(selected),
        "protected_topology_sha256": canonical_sha256(protected),
    }


def load_topologies(
    paths: Sequence[Path | str],
    needed_regions: set[tuple[str, str]],
) -> tuple[dict[tuple[str, str], dict[str, Any]], list[dict[str, Any]]]:
    require(
        len(paths) == TOPOLOGY_FILE_COUNT,
        f"exactly {TOPOLOGY_FILE_COUNT} topology JSONL files are required",
    )
    resolved_paths = [Path(path).resolve(strict=True) for path in paths]
    require(
        len(set(resolved_paths)) == TOPOLOGY_FILE_COUNT,
        "topology JSONL paths are not unique",
    )
    records: dict[tuple[str, str], dict[str, Any]] = {}
    bindings: list[dict[str, Any]] = []
    samples_seen: set[str] = set()
    for path in sorted(resolved_paths, key=str):
        require(path.is_file() and path.stat().st_size > 0, f"invalid topology file: {path}")
        digest = hashlib.sha256()
        file_sample: str | None = None
        row_count = 0
        region_ids: set[str] = set()
        try:
            with path.open("rb") as handle:
                for line_number, raw_line in enumerate(handle, start=1):
                    digest.update(raw_line)
                    row_count += 1
                    text = raw_line.decode("utf-8").strip()
                    require(text != "", f"blank topology JSONL row: {path}:{line_number}")
                    row = validate_json_value(
                        load_json_line(text, f"{path}:{line_number}"),
                        f"{path}:{line_number}",
                    )
                    require(
                        isinstance(row, dict),
                        f"topology row is not an object: {path}:{line_number}",
                    )
                    validate_topology_row(row, f"{path}:{line_number}")
                    sample = str(row["sample"])
                    if file_sample is None:
                        file_sample = sample
                    require(
                        sample == file_sample,
                        f"multiple samples in topology file {path}: {file_sample}/{sample}",
                    )
                    region_id = str(row["region_id"])
                    require(
                        region_id not in region_ids,
                        f"duplicate topology region in {path}: {region_id}",
                    )
                    region_ids.add(region_id)
                    key = (sample, region_id)
                    if key in needed_regions:
                        validate_topology_row(
                            row,
                            f"{path}:{line_number}",
                            require_evidence_linked_ranked=True,
                        )
                        require(key not in records, f"duplicate topology authority key: {key}")
                        stored = dict(row)
                        stored["_source_line"] = line_number
                        stored["_row_sha256"] = hashlib.sha256(
                            raw_line.rstrip(b"\r\n")
                        ).hexdigest()
                        stored["_protected_hashes"] = topology_protected_hashes(row)
                        records[key] = stored
        except (OSError, UnicodeError) as exc:
            raise SidecarContractError(f"failed to read topology JSONL {path}: {exc}") from exc
        require(file_sample is not None and row_count > 0, f"empty topology JSONL: {path}")
        require(
            file_sample not in samples_seen,
            f"duplicate topology sample authority: {file_sample}",
        )
        samples_seen.add(file_sample)
        bindings.append(
            {
                "dataset": file_sample,
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": digest.hexdigest(),
                "row_count": row_count,
                "schema_name": TOPOLOGY_SCHEMA_NAME,
                "schema_version": SCHEMA_VERSION,
            }
        )
    missing = sorted(needed_regions - set(records))
    require(
        not missing,
        f"evidence regions missing from topology authorities: {len(missing)} preview={missing[:3]}",
    )
    return records, sorted(bindings, key=lambda item: item["dataset"])


def hamming_distance(first: str, second: str) -> int:
    require(len(first) == len(second), "patterns have different lengths")
    return sum(left != right for left, right in zip(first, second))


def pattern_vertex(pattern: str) -> int:
    require(
        set(pattern) <= {"R", "A"},
        f"topology vertex requires a complete R/A pattern: {pattern!r}",
    )
    return sum(1 << index for index, code in enumerate(pattern) if code == "A")


def selected_edge_for_pair(
    first: str,
    second: str,
    topology: Mapping[str, Any],
) -> Mapping[str, Any] | None:
    first_vertex = pattern_vertex(first)
    second_vertex = pattern_vertex(second)
    vertices = {first_vertex, second_vertex}
    matches = [
        edge
        for edge in topology.get("representative_best_edges", [])
        if {int(edge["parent_vertex"]), int(edge["child_vertex"])} == vertices
    ]
    require(len(matches) <= 1, f"multiple selected edges for state pair {first}|{second}")
    return matches[0] if matches else None


def expected_topology_relation(
    first: str,
    second: str,
    topology: Mapping[str, Any],
) -> str:
    require(
        "X" not in first and "X" not in second,
        "topology relation requires complete R/A states",
    )
    hamming = hamming_distance(first, second)
    if hamming > 1:
        return HAMMING_GT1_RELATION
    if not bool(topology.get("best_tree_unique")):
        return "HAMMING1_NOT_UNANIMOUS"
    if selected_edge_for_pair(first, second, topology) is not None:
        return HALO_RELATION
    return "HAMMING1_NOT_IN_GLOBAL_BEST"


def normalized_topology_snapshot(topology: Mapping[str, Any]) -> dict[str, Any]:
    """Match the analyzer's immutable topology projection exactly."""
    return {
        "sample": topology.get("sample"),
        "chrom": topology.get("chrom"),
        "unit_id": topology.get("unit_id"),
        "phase_set": topology.get("phase_set"),
        "hp_family": topology.get("hp_family"),
        "active_positions": topology.get("active_positions"),
        "best_tree_unique": bool(topology.get("best_tree_unique", False)),
        "representative_best_edges": topology.get("representative_best_edges", []),
        "representative_best_vertices": topology.get(
            "representative_best_vertices", []
        ),
        "best_tree_tie_count": topology.get("best_tree_tie_count"),
    }


def validate_detail(
    payload: Mapping[str, Any],
    *,
    evidence: Mapping[str, Any],
    topology: Mapping[str, Any],
    label: str,
) -> None:
    require(
        payload.get("schema_name") == DETAIL_SCHEMA_NAME
        and payload.get("schema_version") == SCHEMA_VERSION,
        f"unexpected detail schema at {label}",
    )
    for field in ("dataset", "chrom", "region_id", "hp_raw"):
        require(
            payload.get(field) == evidence.get(field),
            f"detail/evidence identity mismatch at {label}.{field}",
        )
    require(
        payload.get("assessment") == evidence.get("assessment"),
        f"detail/evidence assessment mismatch at {label}",
    )
    positions = payload.get("active_positions")
    require(
        positions == evidence.get("active_positions"),
        f"detail/evidence active_positions mismatch at {label}",
    )
    state_counts = payload.get("state_counts")
    require(isinstance(state_counts, dict), f"detail state_counts missing at {label}")
    normalized_counts: dict[str, int] = {}
    for pattern, count in state_counts.items():
        require(isinstance(pattern, str), f"detail state pattern is not a string at {label}")
        validate_pattern(pattern, int(evidence["n_active_bits"]), f"{label}.state_counts")
        require("X" not in pattern, f"detail analysis state contains X at {label}")
        require(
            isinstance(count, int) and not isinstance(count, bool) and count >= 0,
            f"detail state count invalid at {label}.{pattern}",
        )
        normalized_counts[pattern] = count
    require(
        dict(sorted(normalized_counts.items())) == evidence["analysis_state_counts"],
        f"detail/evidence state counts mismatch at {label}",
    )
    topology_snapshot = payload.get("topology")
    require(isinstance(topology_snapshot, dict), f"detail topology snapshot missing at {label}")
    expected_snapshot = normalized_topology_snapshot(topology)
    require(
        topology_snapshot == expected_snapshot,
        f"detail topology snapshot drift at {label}",
    )

    cpg_positions = payload.get("cpg_positions")
    require(
        isinstance(cpg_positions, list)
        and all(isinstance(item, int) and item > 0 for item in cpg_positions)
        and cpg_positions == sorted(set(cpg_positions)),
        f"invalid detail cpg_positions at {label}",
    )
    profiles = payload.get("state_mean_profiles")
    require(isinstance(profiles, dict), f"detail state_mean_profiles missing at {label}")
    require(
        set(profiles) == set(normalized_counts),
        f"detail profile/state key mismatch at {label}",
    )
    for pattern, profile in profiles.items():
        validate_pattern(str(pattern), int(evidence["n_active_bits"]), f"{label}.profiles")
        require(
            isinstance(profile, list) and len(profile) == len(cpg_positions),
            f"detail profile width mismatch at {label}.{pattern}",
        )
        validate_json_value(profile, f"{label}.profiles.{pattern}")

    read_order = payload.get("read_order")
    require(isinstance(read_order, list), f"detail read_order missing at {label}")
    require(
        len(read_order) == evidence.get("analysis_n"),
        f"detail read_order length mismatch at {label}",
    )
    qnames: set[str] = set()
    read_pattern_counts: Counter[str] = Counter()
    for index, read in enumerate(read_order):
        require(isinstance(read, dict), f"invalid read_order row at {label}[{index}]")
        qname = read.get("qname_sha256")
        pattern = read.get("pattern")
        require(
            isinstance(qname, str) and len(qname) == 64,
            f"invalid qname digest at {label}[{index}]",
        )
        require(qname not in qnames, f"duplicate qname digest at {label}: {qname}")
        qnames.add(qname)
        require(isinstance(pattern, str), f"missing read pattern at {label}[{index}]")
        validate_pattern(pattern, int(evidence["n_active_bits"]), f"{label}.read_order")
        require("X" not in pattern, f"detail read_order contains X at {label}[{index}]")
        read_pattern_counts[pattern] += 1
    require(
        dict(sorted(read_pattern_counts.items())) == dict(sorted(normalized_counts.items())),
        f"detail read_order/state count mismatch at {label}",
    )

    matrix = payload.get("distance_matrix")
    require(matrix is None or isinstance(matrix, list), f"invalid distance_matrix at {label}")
    if isinstance(matrix, list):
        require(len(matrix) == len(read_order), f"distance matrix height mismatch at {label}")
        for index, row in enumerate(matrix):
            require(
                isinstance(row, list) and len(row) == len(read_order),
                f"distance matrix width mismatch at {label}[{index}]",
            )
            validate_json_value(row, f"{label}.distance_matrix[{index}]")
    effects = payload.get("pairwise_effects")
    require(isinstance(effects, list), f"detail pairwise_effects missing at {label}")
    expected_pairs = set(combinations(sorted(normalized_counts), 2))
    normalized_effects: dict[tuple[str, str], dict[str, Any]] = {}
    for index, effect in enumerate(effects):
        require(
            isinstance(effect, dict),
            f"detail pair effect is not an object at {label}[{index}]",
        )
        first = effect.get("first")
        second = effect.get("second")
        require(
            isinstance(first, str)
            and isinstance(second, str)
            and "X" not in first
            and "X" not in second,
            f"detail pair effect contains X or invalid state at {label}[{index}]",
        )
        require(
            first in normalized_counts and second in normalized_counts,
            f"detail pair effect references absent state at {label}[{index}]",
        )
        association, _ = normalize_pair_effect(
            effect,
            topology=topology,
            n_active_bits=int(evidence["n_active_bits"]),
            label=f"{label}.pairwise_effects[{index}]",
        )
        require(
            all(
                association[field] is not None
                for field in (
                    "between_mean",
                    "pooled_within_mean",
                    "distance_contrast",
                    "standardized_effect",
                )
            ),
            f"detail pair effect has null numeric metrics at {label}[{index}]",
        )
        pair_key = tuple(sorted((first, second)))
        require(
            pair_key not in normalized_effects,
            f"duplicate unordered detail pair at {label}: {pair_key}",
        )
        normalized_effects[pair_key] = association
    actual_pairs = set(normalized_effects)
    require(
        actual_pairs == expected_pairs,
        f"detail pairwise_effects is not the exact unordered analysis-state pair set "
        f"at {label}: missing={sorted(expected_pairs - actual_pairs)} "
        f"unexpected={sorted(actual_pairs - expected_pairs)}",
    )

    ranked_effects = sorted(
        normalized_effects.items(),
        key=lambda item: (
            -float(item[1]["distance_contrast"]),
            -float(item[1]["standardized_effect"]),
            item[0],
        ),
    )
    require(ranked_effects, f"detail pairwise_effects is empty at {label}")
    detail_best_key, detail_best = ranked_effects[0]
    evidence_best = evidence_best_pair_effect(evidence)
    require(evidence_best is not None, f"detail exists without evidence best_pair at {label}")
    evidence_best_key = tuple(
        sorted((str(evidence_best["first"]), str(evidence_best["second"])))
    )
    require(
        evidence_best_key == detail_best_key,
        f"evidence/detail best_pair mismatch at {label}: "
        f"{evidence_best_key}!={detail_best_key}",
    )
    for field in PAIR_EFFECT_FIELDS:
        evidence_value = evidence_best.get(field)
        detail_value = detail_best.get(field)
        values_match = (
            matches_evidence_float_serialization(evidence_value, detail_value)
            if field
            in {
                "between_mean",
                "pooled_within_mean",
                "distance_contrast",
                "standardized_effect",
            }
            else evidence_value == detail_value
        )
        require(
            values_match,
            f"evidence/detail best_pair metric mismatch at {label}.{field}: "
            f"{evidence_value!r}!={detail_value!r}",
        )


def load_details(
    path: Path | str,
    evidence_rows: Mapping[tuple[str, str, str, str], Mapping[str, Any]],
    topologies: Mapping[tuple[str, str], Mapping[str, Any]],
) -> tuple[dict[tuple[str, str, str, str], dict[str, Any]], dict[str, Any]]:
    binding = file_binding(path)
    resolved = Path(binding["path"])
    details: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    record_count = 0
    try:
        with gzip.open(resolved, "rt", encoding="utf-8", newline="") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                text = raw_line.strip()
                require(text != "", f"blank detail JSONL row: {resolved}:{line_number}")
                payload = validate_json_value(
                    load_json_line(text, f"{resolved}:{line_number}"),
                    f"{resolved}:{line_number}",
                )
                require(
                    isinstance(payload, dict),
                    f"detail row is not an object: {resolved}:{line_number}",
                )
                key = evidence_key(payload)
                require(
                    key in evidence_rows,
                    f"orphan detail row at {resolved}:{line_number}: {key}",
                )
                require(key not in details, f"duplicate detail key: {key}")
                evidence = evidence_rows[key]
                topology = topologies[(key[0], key[2])]
                validate_detail(
                    payload,
                    evidence=evidence,
                    topology=topology,
                    label=f"{resolved}:{line_number}",
                )
                stored = dict(payload)
                stored["_source_line"] = line_number
                stored["_record_sha256"] = hashlib.sha256(
                    text.encode("utf-8")
                ).hexdigest()
                details[key] = stored
                record_count += 1
    except (OSError, UnicodeError) as exc:
        raise SidecarContractError(f"failed to read details JSONL {resolved}: {exc}") from exc
    binding["row_count"] = record_count
    binding["schema_name"] = DETAIL_SCHEMA_NAME
    binding["schema_version"] = SCHEMA_VERSION
    return details, binding


def verify_evidence_topology_identity(
    evidence: Mapping[str, Any],
    topology: Mapping[str, Any],
    key: tuple[str, str, str, str],
) -> None:
    expected = {
        "sample": evidence["dataset"],
        "chrom": evidence["chrom"],
        "region_id": evidence["region_id"],
        "unit_id": evidence["unit_id"],
        "phase_set": evidence["phase_set"],
        "hp_family": evidence["hp_family"],
        "active_positions": evidence["active_positions"],
        "active_bit_count": evidence["n_active_bits"],
    }
    for field, value in expected.items():
        require(
            topology.get(field) == value,
            f"evidence/topology identity mismatch {key}.{field}: "
            f"{value!r}!={topology.get(field)!r}",
        )


def node_projection(
    evidence: Mapping[str, Any],
    topology: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    input_counts = evidence["input_state_counts"]
    analysis_counts = evidence["analysis_state_counts"]
    selected_labels = {
        int(item["vertex"]): item.get("label")
        for item in topology.get("representative_best_vertices", [])
    }
    selected_exemplar_available = isinstance(topology.get("best_tree_unique"), bool)
    nodes: list[dict[str, Any]] = []
    for pattern in sorted(set(input_counts) | set(analysis_counts)):
        item = {
            "pattern": pattern,
            "input_n": input_counts.get(pattern),
            "analysis_n": analysis_counts.get(pattern),
        }
        vertex = pattern_vertex(pattern)
        item.update(
            {
                "vertex": vertex,
                "projection": (
                    "TOPOLOGY_VERTEX_REFERENCE_ONLY"
                    if selected_exemplar_available
                    else "TOPOLOGY_VERTEX_REFERENCE_ONLY_SELECTED_EXEMPLAR_UNAVAILABLE"
                ),
                "selected_exemplar_member": vertex in selected_labels,
                "selected_exemplar_label": selected_labels.get(vertex),
            }
        )
        nodes.append(item)
    return sorted(nodes, key=lambda item: (item["vertex"], item["pattern"])), []


def numeric_effect_value(effect: Mapping[str, Any], field: str, label: str) -> float | None:
    value = effect.get(field)
    if value is None:
        return None
    require(
        isinstance(value, (int, float)) and not isinstance(value, bool),
        f"non-numeric pair effect {label}.{field}: {value!r}",
    )
    number = float(value)
    require(math.isfinite(number), f"non-finite pair effect {label}.{field}")
    return number


def matches_evidence_float_serialization(
    evidence_value: Any,
    detail_value: Any,
) -> bool:
    """Reproduce the analyzer's ``.10g`` TSV serialization before comparison."""
    try:
        evidence_number = float(evidence_value)
        detail_number = float(detail_value)
    except (TypeError, ValueError):
        return False
    if not math.isfinite(evidence_number) or not math.isfinite(detail_number):
        return False
    return evidence_number == float(f"{detail_number:.10g}")


def normalize_pair_effect(
    effect: Mapping[str, Any],
    *,
    topology: Mapping[str, Any],
    n_active_bits: int,
    label: str,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    require(isinstance(effect, dict), f"pair effect is not an object at {label}")
    first = effect.get("first")
    second = effect.get("second")
    require(
        isinstance(first, str) and isinstance(second, str) and first != second,
        f"invalid pair states at {label}",
    )
    validate_pattern(first, n_active_bits, label)
    validate_pattern(second, n_active_bits, label)
    require(
        "X" not in first and "X" not in second,
        f"pair effect contains X at {label}",
    )
    hamming = hamming_distance(first, second)
    source_hamming = effect.get("hamming")
    require(
        isinstance(source_hamming, int)
        and not isinstance(source_hamming, bool)
        and source_hamming == hamming,
        f"pair Hamming mismatch at {label}",
    )
    source_relation = effect.get("topology_relation")
    require(isinstance(source_relation, str), f"missing pair topology_relation at {label}")
    expected_relation = expected_topology_relation(first, second, topology)
    require(
        source_relation == expected_relation,
        f"pair/topology relation mismatch at {label}: "
        f"{source_relation!r}!={expected_relation!r}",
    )
    association = {
        "first": first,
        "second": second,
        "hamming": hamming,
        "topology_relation": expected_relation,
        "projection_kind": (
            "PAIR_BAND_ONLY"
            if hamming > 1
            else "EDGE_HALO_ELIGIBLE"
            if expected_relation == HALO_RELATION
            else "PAIR_ASSOCIATION_ONLY"
        ),
        "between_mean": numeric_effect_value(effect, "between_mean", label),
        "pooled_within_mean": numeric_effect_value(
            effect, "pooled_within_mean", label
        ),
        "distance_contrast": numeric_effect_value(effect, "distance_contrast", label),
        "standardized_effect": numeric_effect_value(
            effect, "standardized_effect", label
        ),
    }
    halo: dict[str, Any] | None = None
    if expected_relation == HALO_RELATION:
        require(hamming == 1, f"invalid halo projection at {label}")
        edge = selected_edge_for_pair(first, second, topology)
        require(edge is not None, f"halo lacks selected topology edge at {label}")
        halo = {
            "relation": HALO_RELATION,
            "first": first,
            "second": second,
            "parent_vertex": int(edge["parent_vertex"]),
            "child_vertex": int(edge["child_vertex"]),
            "distance_contrast": association["distance_contrast"],
            "standardized_effect": association["standardized_effect"],
        }
    return association, halo


def evidence_best_pair_effect(evidence: Mapping[str, Any]) -> dict[str, Any] | None:
    best_pair = evidence.get("best_pair")
    if best_pair is None:
        return None
    first, second = str(best_pair).split("|")
    return {
        "first": first,
        "second": second,
        "hamming": evidence["best_pair_hamming"],
        "between_mean": evidence["best_pair_between_mean"],
        "pooled_within_mean": evidence["best_pair_pooled_within_mean"],
        "distance_contrast": evidence["best_pair_distance_contrast"],
        "standardized_effect": evidence["best_pair_standardized_effect"],
        "topology_relation": evidence["best_pair_topology_relation"],
    }


def matrix_rank_key(
    key: tuple[str, str, str, str],
    evidence: Mapping[str, Any],
) -> tuple[Any, ...]:
    q_by = evidence.get("q_by")
    p_value = evidence.get("permanova_p")
    r_squared = evidence.get("permanova_r2")
    effect = evidence.get("best_pair_distance_contrast")
    return (
        ASSESSMENT_ORDER[str(evidence["assessment"])],
        math.inf if q_by is None else float(q_by),
        math.inf if p_value is None else float(p_value),
        -(float(r_squared) if r_squared is not None else -math.inf),
        -(float(effect) if effect is not None else -math.inf),
        key,
    )


def choose_matrix_cases(
    evidence_rows: Mapping[tuple[str, str, str, str], Mapping[str, Any]],
    details: Mapping[tuple[str, str, str, str], Mapping[str, Any]],
    top_n: int,
) -> tuple[dict[tuple[str, str, str, str], int], int]:
    require(top_n >= 0, "matrix_top_n must be non-negative")
    candidates = [
        key
        for key, detail in details.items()
        if evidence_rows[key]["evaluation_status"] == "EVALUABLE"
        and isinstance(detail.get("distance_matrix"), list)
    ]
    candidates.sort(key=lambda key: matrix_rank_key(key, evidence_rows[key]))
    selected = {
        key: rank
        for rank, key in enumerate(candidates[:top_n], start=1)
    }
    return selected, len(candidates)


def build_bundle(
    key: tuple[str, str, str, str],
    evidence: Mapping[str, Any],
    topology: Mapping[str, Any],
    detail: Mapping[str, Any] | None,
    matrix_rank: int | None,
) -> dict[str, Any]:
    verify_evidence_topology_identity(evidence, topology, key)
    nodes, _ = node_projection(evidence, topology)
    bundle_id = hashlib.sha256("\x1f".join(key).encode("utf-8")).hexdigest()
    compact_evidence = {
        field: evidence.get(field)
        for field in EVIDENCE_FIELDS
        if field not in IDENTITY_FIELDS
        and field not in {"input_state_counts_json", "analysis_state_counts_json"}
    }
    pair_effects: list[Mapping[str, Any]] = []
    if detail is not None:
        raw_effects = detail.get("pairwise_effects")
        require(isinstance(raw_effects, list), f"pairwise effects missing for {key}")
        pair_effects = raw_effects
    else:
        fallback = evidence_best_pair_effect(evidence)
        if fallback is not None:
            pair_effects = [fallback]

    associations: list[dict[str, Any]] = []
    edge_halos: list[dict[str, Any]] = []
    seen_pairs: set[tuple[str, str]] = set()
    for index, effect in enumerate(pair_effects):
        association, halo = normalize_pair_effect(
            effect,
            topology=topology,
            n_active_bits=int(evidence["n_active_bits"]),
            label=f"{key}.pairwise_effects[{index}]",
        )
        pair_key = tuple(sorted((association["first"], association["second"])))
        require(pair_key not in seen_pairs, f"duplicate pair association for {key}: {pair_key}")
        seen_pairs.add(pair_key)
        associations.append(association)
        if halo is not None:
            edge_halos.append(halo)
    associations.sort(key=lambda item: (item["first"], item["second"]))
    edge_halos.sort(
        key=lambda item: (
            item["parent_vertex"],
            item["child_vertex"],
            item["first"],
            item["second"],
        )
    )
    require(
        all("X" not in item["first"] + item["second"] for item in associations),
        f"X-containing pair escaped complete-R/A source contract for {key}",
    )
    complete_associations = associations
    pair_bands = [
        item for item in complete_associations if item["hamming"] > 1
    ]
    require(
        all(item["topology_relation"] == HAMMING_GT1_RELATION for item in pair_bands),
        f"Hamming>1 association escaped pair-band contract for {key}",
    )
    require(
        all(item["relation"] == HALO_RELATION for item in edge_halos),
        f"non-unanimous edge halo for {key}",
    )

    detail_summary: dict[str, Any] = {
        "available": detail is not None,
        "large_matrix_embedded": False,
        "matrix_rank": None,
        "matrix_source_status": "NOT_AVAILABLE",
        "source_line": None,
        "record_sha256": None,
    }
    if detail is not None:
        matrix = detail.get("distance_matrix")
        detail_summary.update(
            {
                "source_line": detail["_source_line"],
                "record_sha256": detail["_record_sha256"],
                "matrix_source_status": (
                    "AVAILABLE_NOT_EMBEDDED"
                    if isinstance(matrix, list)
                    else "SOURCE_NULL"
                ),
            }
        )
        if matrix_rank is not None:
            require(
                evidence["evaluation_status"] == "EVALUABLE"
                and isinstance(matrix, list),
                f"non-evaluable or null matrix selected for embedding: {key}",
            )
            detail_summary.update(
                {
                    "large_matrix_embedded": True,
                    "matrix_rank": matrix_rank,
                    "matrix_source_status": "EMBEDDED_TOP_CASE",
                    "top_case_payload": {
                        "cpg_positions": detail["cpg_positions"],
                        "state_mean_profiles": detail["state_mean_profiles"],
                        "balanced_n_per_state": detail.get("balanced_n_per_state"),
                        "rarefaction_r2": detail.get("rarefaction_r2"),
                        "read_order": detail["read_order"],
                        "distance_matrix": matrix,
                    },
                }
            )

    return {
        "bundle_id": bundle_id,
        "grain": {
            "dataset": evidence["dataset"],
            "chrom": evidence["chrom"],
            "region_id": evidence["region_id"],
            "phase_set": evidence["phase_set"],
            "hp_family": evidence["hp_family"],
            "hp_raw": evidence["hp_raw"],
        },
        "evidence": compact_evidence,
        "input_state_counts": evidence["input_state_counts"],
        "analysis_state_counts": evidence["analysis_state_counts"],
        "complete_node_projection": nodes,
        "partial_subcube_evidence": {
            "states": [],
            "pairs": [],
            "source_status": PARTIAL_X_SOURCE_STATUS,
            "topology_projection_forbidden": True,
        },
        "pair_associations": complete_associations,
        "pair_bands": pair_bands,
        "edge_halos": edge_halos,
        "detail": detail_summary,
    }


def topology_anchor(
    topology: Mapping[str, Any],
) -> dict[str, Any]:
    selected_exemplar_available = isinstance(topology.get("best_tree_unique"), bool)
    return {
        "authority_dataset": topology["sample"],
        "source_line": topology["_source_line"],
        "row_sha256": topology["_row_sha256"],
        **topology["_protected_hashes"],
        "overlay_only": True,
        "selected_exemplar_available": selected_exemplar_available,
        "authority_unit_status": topology.get("unit_status"),
        "authority_read_af_status": topology.get("read_af_status"),
        "topology_counts_af_selected_tree_embedded": False,
    }


def build_sidecar(
    evidence_tsv: Path | str,
    details_jsonl: Path | str,
    topology_jsonls: Sequence[Path | str],
    *,
    matrix_top_n: int = DEFAULT_MATRIX_TOP_N,
) -> dict[str, Any]:
    """Load, validate, and assemble the complete deterministic sidecar payload."""
    evidence_rows, evidence_binding = load_evidence(evidence_tsv)
    needed_regions = {(key[0], key[2]) for key in evidence_rows}
    topologies, topology_bindings = load_topologies(topology_jsonls, needed_regions)
    topology_binding_by_dataset = {
        item["dataset"]: item for item in topology_bindings
    }
    evidence_datasets = {key[0] for key in evidence_rows}
    require(
        evidence_datasets <= set(topology_binding_by_dataset),
        "evidence dataset lacks topology authority",
    )
    for key, evidence in evidence_rows.items():
        verify_evidence_topology_identity(
            evidence,
            topologies[(key[0], key[2])],
            key,
        )
        best = evidence_best_pair_effect(evidence)
        if best is not None:
            normalize_pair_effect(
                best,
                topology=topologies[(key[0], key[2])],
                n_active_bits=int(evidence["n_active_bits"]),
                label=f"{key}.best_pair",
            )

    details, details_binding = load_details(
        details_jsonl, evidence_rows, topologies
    )
    matrix_ranks, matrix_candidate_count = choose_matrix_cases(
        evidence_rows, details, matrix_top_n
    )

    region_bundles: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    top_cases: list[dict[str, Any]] = []
    for key in sorted(
        evidence_rows,
        key=lambda item: (
            item[0],
            natural_chrom_key(item[1]),
            item[2],
            HP_ORDER.get(item[3], 10**9),
            item[3],
        ),
    ):
        evidence = evidence_rows[key]
        topology = topologies[(key[0], key[2])]
        bundle = build_bundle(
            key,
            evidence,
            topology,
            details.get(key),
            matrix_ranks.get(key),
        )
        region_bundles[(key[0], key[2])].append(bundle)
        if key in matrix_ranks:
            top_cases.append(
                {
                    "rank": matrix_ranks[key],
                    "bundle_id": bundle["bundle_id"],
                    "dataset": key[0],
                    "chrom": key[1],
                    "region_id": key[2],
                    "phase_set": evidence["phase_set"],
                    "hp_raw": key[3],
                    "assessment": evidence["assessment"],
                }
            )
    top_cases.sort(key=lambda item: item["rank"])

    datasets: list[dict[str, Any]] = []
    for dataset in sorted(evidence_datasets):
        regions: list[dict[str, Any]] = []
        dataset_region_keys = sorted(
            (key for key in region_bundles if key[0] == dataset),
            key=lambda item: (
                natural_chrom_key(str(topologies[item]["chrom"])),
                item[1],
            ),
        )
        for region_key in dataset_region_keys:
            topology = topologies[region_key]
            bundles = sorted(
                region_bundles[region_key],
                key=lambda item: (
                    HP_ORDER.get(item["grain"]["hp_raw"], 10**9),
                    item["grain"]["hp_raw"],
                ),
            )
            regions.append(
                {
                    "region_id": topology["region_id"],
                    "chrom": topology["chrom"],
                    "unit_id": topology["unit_id"],
                    "phase_set": topology["phase_set"],
                    "hp_family": topology["hp_family"],
                    "active_positions": topology["active_positions"],
                    "active_bit_count": topology["active_bit_count"],
                    "topology_anchor": topology_anchor(
                        topology
                    ),
                    "exact_raw_hp_bundles": bundles,
                }
            )
        datasets.append(
            {
                "dataset": dataset,
                "topology_authority": topology_binding_by_dataset[dataset],
                "regions": regions,
            }
        )

    all_bundles = [
        bundle
        for dataset in datasets
        for region in dataset["regions"]
        for bundle in region["exact_raw_hp_bundles"]
    ]
    summary = {
        "topology_authority_files": len(topology_bindings),
        "topology_authority_datasets": len(topology_binding_by_dataset),
        "evidence_datasets": len(datasets),
        "regions": sum(len(dataset["regions"]) for dataset in datasets),
        "exact_raw_hp_bundles": len(all_bundles),
        "evaluable_bundles": sum(
            bundle["evidence"]["evaluation_status"] == "EVALUABLE"
            for bundle in all_bundles
        ),
        "partial_subcube_states": sum(
            len(bundle["partial_subcube_evidence"]["states"])
            for bundle in all_bundles
        ),
        "partial_x_source_status": PARTIAL_X_SOURCE_STATUS,
        "pair_bands_hamming_gt1": sum(
            len(bundle["pair_bands"]) for bundle in all_bundles
        ),
        "edge_halos": sum(len(bundle["edge_halos"]) for bundle in all_bundles),
        "matrix_candidates": matrix_candidate_count,
        "matrices_embedded": len(matrix_ranks),
        "matrix_top_n": matrix_top_n,
        "assessment": dict(
            sorted(
                Counter(
                    str(bundle["evidence"]["assessment"]) for bundle in all_bundles
                ).items()
            )
        ),
    }
    checks = {
        "seven_unique_topology_authorities": len(topology_bindings)
        == TOPOLOGY_FILE_COUNT
        == len(topology_binding_by_dataset),
        "all_evidence_exact_ps_topology_bound": True,
        "exact_raw_hp_bundles_not_collapsed": True,
        "complete_ra_states_reference_vertices_only": True,
        "x_states_never_projected": all(
            "X" not in node["pattern"]
            for bundle in all_bundles
            for node in bundle["complete_node_projection"]
        ),
        "analyzer_state_counts_reads_effects_are_complete_ra_only": all(
            all("X" not in pattern for pattern in bundle["input_state_counts"])
            and all("X" not in pattern for pattern in bundle["analysis_state_counts"])
            and all(
                "X" not in pair["first"] + pair["second"]
                for pair in bundle["pair_associations"]
            )
            for bundle in all_bundles
        ),
        "partial_x_not_claimed_without_hash_bound_source": (
            all(
                not bundle["partial_subcube_evidence"]["states"]
                and not bundle["partial_subcube_evidence"]["pairs"]
                and bundle["partial_subcube_evidence"]["source_status"]
                == PARTIAL_X_SOURCE_STATUS
                for bundle in all_bundles
            )
        ),
        "hamming_gt1_is_pair_band_only": all(
            pair["hamming"] > 1
            and pair["topology_relation"] == HAMMING_GT1_RELATION
            for bundle in all_bundles
            for pair in bundle["pair_bands"]
        ),
        "edge_halo_only_global_best_unanimous": all(
            halo["relation"] == HALO_RELATION
            for bundle in all_bundles
            for halo in bundle["edge_halos"]
        ),
        "topology_counts_af_selected_tree_not_embedded_or_rescored": True,
        "large_matrices_only_evaluable_top_cases": all(
            not bundle["detail"]["large_matrix_embedded"]
            or bundle["evidence"]["evaluation_status"] == "EVALUABLE"
            for bundle in all_bundles
        ),
        "missing_numeric_values_preserved_as_json_null": True,
    }
    require(all(checks.values()), "derived sidecar checks did not all pass")
    payload = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "status": "PASS",
        "claim_ceiling": (
            "pattern-conditioned regional methylation association only; "
            "not ancestry, clone identity, causality, or topology rescoring"
        ),
        "contract": {
            "primary_grain": (
                "dataset x chrom x exact phase_set x exact raw hp x topology region"
            ),
            "complete_ra_projection": "TOPOLOGY_VERTEX_REFERENCE_ONLY",
            "partial_x_projection": "NOT_EMITTED_WITHOUT_HASH_BOUND_SOURCE",
            "partial_x_source_status": PARTIAL_X_SOURCE_STATUS,
            "hamming_gt1_projection": "PAIR_BAND_ONLY",
            "edge_halo_gate": HALO_RELATION,
            "topology_authority": (
                "external immutable JSONL; counts, AF, edge incidence and selected "
                "representative tree are neither embedded nor recomputed"
            ),
            "missing_value_encoding": "JSON_NULL",
            "large_matrix_policy": (
                "only deterministic top-N evaluable cases with an available source matrix"
            ),
        },
        "sources": {
            "evidence_tsv": evidence_binding,
            "details_jsonl": details_binding,
            "topology_jsonl": topology_bindings,
        },
        "external_tables": {
            "full_evidence_tsv": evidence_binding,
            "full_details_jsonl": details_binding,
        },
        "summary": summary,
        "checks": checks,
        "top_cases": top_cases,
        "datasets": datasets,
    }
    return validate_json_value(payload, "sidecar")


def archive_failed_staging(temporary: Path, target: Path) -> Path | None:
    """Preserve a failed staging file beside the intended output for audit."""
    if not temporary.exists():
        return None
    archive_root = target.parent / FAILED_STAGING_DIRECTORY
    archive_root.mkdir(parents=True, exist_ok=True)
    archived = archive_root / temporary.name
    suffix = 1
    while archived.exists():
        archived = archive_root / f"{temporary.name}.{suffix}"
        suffix += 1
    temporary.rename(archived)
    return archived


def deterministic_gzip_json(path: Path, payload: Mapping[str, Any]) -> None:
    temporary_handle = tempfile.NamedTemporaryFile(
        mode="wb",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    )
    temporary = Path(temporary_handle.name)
    temporary_handle.close()
    try:
        with temporary.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=raw_handle,
                compresslevel=6,
                mtime=0,
            ) as gzip_handle:
                with io.TextIOWrapper(
                    gzip_handle,
                    encoding="utf-8",
                    newline="\n",
                    write_through=True,
                ) as text_handle:
                    json.dump(
                        payload,
                        text_handle,
                        ensure_ascii=False,
                        sort_keys=True,
                        separators=(",", ":"),
                        allow_nan=False,
                    )
                    text_handle.write("\n")
            raw_handle.flush()
            os.fsync(raw_handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        archive_failed_staging(temporary, path)
        raise


def atomic_write_text(path: Path, text: str) -> None:
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    )
    temporary = Path(handle.name)
    try:
        with handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        archive_failed_staging(temporary, path)
        raise


def write_outputs(
    payload: Mapping[str, Any],
    output_dir: Path | str,
) -> dict[str, Any]:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    sidecar_path = output_root / OUTPUT_FILENAME
    receipt_path = output_root / RECEIPT_FILENAME
    require(not sidecar_path.exists(), f"refusing to overwrite output: {sidecar_path}")
    require(not receipt_path.exists(), f"refusing to overwrite output: {receipt_path}")

    deterministic_gzip_json(sidecar_path, payload)
    output_binding = {
        "relative_path": OUTPUT_FILENAME,
        "size_bytes": sidecar_path.stat().st_size,
        "sha256": sha256_file(sidecar_path),
    }
    with gzip.open(sidecar_path, "rt", encoding="utf-8") as handle:
        round_trip = load_json_line(handle.read(), str(sidecar_path))
    require(
        round_trip.get("schema_name") == SCHEMA_NAME
        and round_trip.get("schema_version") == SCHEMA_VERSION
        and round_trip.get("status") == "PASS",
        "sidecar gzip round-trip validation failed",
    )
    receipt = {
        "schema_name": RECEIPT_SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "status": "PASS",
        "all_pass": True,
        "overlay_only": True,
        "topology_authority_unchanged": {
            "counts": True,
            "af": True,
            "selected_tree": True,
        },
        "claim_ceiling": payload["claim_ceiling"],
        "parameters": {
            "matrix_top_n": payload["summary"]["matrix_top_n"],
            "output_filename": OUTPUT_FILENAME,
        },
        "inputs": payload["sources"],
        "output": output_binding,
        "summary": payload["summary"],
        "checks": {
            **payload["checks"],
            "sidecar_gzip_round_trip": True,
            "deterministic_gzip_mtime_zero": True,
            "receipt_hash_binds_all_inputs_and_output": True,
        },
    }
    receipt_text = json.dumps(
        receipt,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ) + "\n"
    atomic_write_text(receipt_path, receipt_text)
    return receipt


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-tsv", type=Path, required=True)
    parser.add_argument("--details-jsonl", type=Path, required=True)
    parser.add_argument(
        "--topology-jsonl",
        type=Path,
        action="append",
        required=True,
        help="Repeat exactly seven times, once per technical dataset.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--matrix-top-n", type=int, default=DEFAULT_MATRIX_TOP_N)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    try:
        args = parse_args(argv)
        payload = build_sidecar(
            args.evidence_tsv,
            args.details_jsonl,
            args.topology_jsonl,
            matrix_top_n=args.matrix_top_n,
        )
        receipt = write_outputs(payload, args.output_dir)
    except (SidecarContractError, OSError, ValueError) as exc:
        print(f"FAIL_CLOSED: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": receipt["status"],
                **receipt["summary"],
                "output": str((args.output_dir / OUTPUT_FILENAME).resolve()),
                "receipt": str((args.output_dir / RECEIPT_FILENAME).resolve()),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
