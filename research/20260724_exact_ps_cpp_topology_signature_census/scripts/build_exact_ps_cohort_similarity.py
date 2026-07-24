#!/usr/bin/env python3
"""Build a provenance-bound exact-PS cohort similarity sidecar.

The sidecar keeps two deliberately separate comparison grains:

1. all-seven distribution profiles, compared with equal-dimension-weighted
   Jensen-Shannon distance; and
2. HCC1395 versus HCC1395_DORADO locus-matched technical concordance.

It never pairs a repeated ordered-locus key by choosing an HP label.  Such keys
are compared only as multisets because HP1/HP2 labels can flip across runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import statistics
import sys
import time
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


HERE = Path(__file__).resolve().parent
TOPOLOGY_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1"
)
CENSUS_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260724_exact_ps_cpp_topology_signature_census/all7_v1"
)

SCHEMA_NAME = "intersubmod.exact_ps_cohort_similarity"
SCHEMA_VERSION = "1.0.0"
SOURCE_PROFILE = "20260724_exact_ps_strict_guard1000_signature_census_all7_v1"
TASK_TYPE = "B_COMPREHENSIVE_VALIDATION"

EXPECTED_SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
EXPECTED_CHROMS = tuple(f"chr{index}" for index in range(1, 23))
HCC_PAIR = ("HCC1395", "HCC1395_DORADO")

STATUS_BINS = (
    "ranked",
    "family_incomplete",
    "no_active_alt",
    "zero_denominator",
    "recurrence_required",
)
RESOLUTION_BINS = (
    "UNIQUE_TREE",
    "TIED_SAME_TOPOLOGY",
    "TIED_CROSS_TOPOLOGY",
)
MORPHOLOGY_BINS = (
    "Single-only",
    "Sister-only",
    "Direct-only",
    "Sister+direct",
    "Unresolved-cross-coarse",
)
ACTIVE_K_BINS = tuple(str(index) for index in range(13))
HP_BINS = ("1", "2")
DIMENSION_ORDER = (
    "status",
    "resolution",
    "morphology",
    "active_k",
    "hp_symmetric",
)
DIMENSION_WEIGHT = 1.0 / len(DIMENSION_ORDER)

BOOTSTRAP_SEED = 20260725
BOOTSTRAP_REPLICATES = 2000


class SimilarityError(RuntimeError):
    """Raised when a frozen source or arithmetic contract fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SimilarityError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing source artifact: {path}")
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def verify_identity(spec: Mapping[str, Any], label: str) -> Path:
    path = Path(str(spec.get("path", "")))
    require(path.is_file(), f"missing bound artifact {label}: {path}")
    require(
        isinstance(spec.get("size_bytes"), int)
        and path.stat().st_size == int(spec["size_bytes"]),
        f"size drift for {label}: {path}",
    )
    require(
        isinstance(spec.get("sha256"), str)
        and sha256_file(path) == spec["sha256"],
        f"SHA-256 drift for {label}: {path}",
    )
    return path


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing JSON: {path}")
    try:
        with path.open(encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise SimilarityError(f"cannot parse {path}: {exc}") from exc
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def iter_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    require(path.is_file(), f"missing JSONL: {path}")
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            try:
                row = json.loads(raw)
            except json.JSONDecodeError as exc:
                raise SimilarityError(
                    f"invalid JSONL {path}:{line_number}: {exc}"
                ) from exc
            require(isinstance(row, dict), f"non-object JSONL row {path}:{line_number}")
            yield row


def all_true(value: Any, label: str) -> None:
    if isinstance(value, Mapping):
        for key, child in value.items():
            all_true(child, f"{label}.{key}")
    elif isinstance(value, bool):
        require(value, f"false authority check: {label}")


def normalized_counts(
    counts: Mapping[str, int], bins: Sequence[str], denominator: int, label: str
) -> list[float]:
    require(denominator > 0, f"{label}: denominator must be positive")
    observed = sum(int(counts.get(key, 0)) for key in bins)
    require(observed == denominator, f"{label}: {observed} != {denominator}")
    return [int(counts.get(key, 0)) / denominator for key in bins]


def dimension_record(
    counts: Mapping[str, int],
    bins: Sequence[str],
    denominator_name: str,
    denominator: int,
    label: str,
) -> dict[str, Any]:
    proportions = normalized_counts(counts, bins, denominator, label)
    return {
        "denominator_name": denominator_name,
        "denominator": denominator,
        "bins": list(bins),
        "counts": {key: int(counts.get(key, 0)) for key in bins},
        "proportions": {
            key: proportions[index] for index, key in enumerate(bins)
        },
    }


def jensen_shannon_similarity(
    left: Sequence[float], right: Sequence[float]
) -> float:
    """Return 1 - Jensen-Shannon distance using base-2 logarithms."""

    require(len(left) == len(right) and len(left) > 0, "JSD vector shape mismatch")
    require(
        math.isclose(sum(left), 1.0, abs_tol=1e-12)
        and math.isclose(sum(right), 1.0, abs_tol=1e-12),
        "JSD vectors must each sum to one",
    )
    midpoint = [(a + b) / 2.0 for a, b in zip(left, right)]

    def kl_divergence(values: Sequence[float]) -> float:
        return sum(
            value * math.log2(value / reference)
            for value, reference in zip(values, midpoint)
            if value > 0.0
        )

    divergence = (kl_divergence(left) + kl_divergence(right)) / 2.0
    distance = math.sqrt(max(0.0, divergence))
    require(-1e-12 <= distance <= 1.0 + 1e-12, "JSD distance outside [0,1]")
    return 1.0 - min(1.0, max(0.0, distance))


def composite_vector(profile: Mapping[str, Any]) -> list[float]:
    vector: list[float] = []
    for dimension in DIMENSION_ORDER:
        if dimension == "hp_symmetric":
            values = list(profile["dimensions"]["hp"]["comparison_vector"])
        else:
            record = profile["dimensions"][dimension]
            values = [record["proportions"][key] for key in record["bins"]]
        vector.extend(DIMENSION_WEIGHT * float(value) for value in values)
    require(math.isclose(sum(vector), 1.0, abs_tol=1e-12), "composite != 1")
    return vector


def chromosome_sort_key(key: tuple[str, tuple[int, ...]]) -> tuple[Any, ...]:
    chrom, positions = key
    return (int(chrom.removeprefix("chr")), positions)


def ordered_locus_key(group: Mapping[str, Any]) -> tuple[str, tuple[int, ...]]:
    positions = tuple(int(value) for value in group.get("positions") or ())
    require(positions and tuple(sorted(positions)) == positions, "positions not ordered")
    return str(group.get("chrom")), positions


def partition_locus_keys(
    left: Mapping[tuple[str, tuple[int, ...]], Sequence[str]],
    right: Mapping[tuple[str, tuple[int, ...]], Sequence[str]],
) -> tuple[
    list[tuple[tuple[str, tuple[int, ...]], str, str]],
    list[dict[str, Any]],
]:
    """Partition common keys without ever selecting among repeated HP units."""

    unambiguous: list[tuple[tuple[str, tuple[int, ...]], str, str]] = []
    ambiguous: list[dict[str, Any]] = []
    for key in sorted(set(left) & set(right), key=chromosome_sort_key):
        left_ids = list(left[key])
        right_ids = list(right[key])
        if len(left_ids) == 1 and len(right_ids) == 1:
            unambiguous.append((key, left_ids[0], right_ids[0]))
        else:
            ambiguous.append(
                {
                    "chrom": key[0],
                    "positions": list(key[1]),
                    "left_multiplicity": len(left_ids),
                    "right_multiplicity": len(right_ids),
                    "reason": "repeated_ordered_locus_key_no_hp_pairing",
                }
            )
    return unambiguous, ambiguous


def normalized_selected_edges(row: Mapping[str, Any]) -> tuple[tuple[int, int, int], ...]:
    edges = row.get("representative_best_edges")
    require(isinstance(edges, list), "unique selected tree lacks representative edges")
    normalized = tuple(
        sorted(
            (
                int(edge["parent_vertex"]),
                int(edge["child_vertex"]),
                int(edge["acquired_position"]),
            )
            for edge in edges
        )
    )
    require(len(normalized) == len(set(normalized)), "duplicate selected edge")
    return normalized


def signature_set(row: Mapping[str, Any]) -> tuple[str, ...]:
    return tuple(
        sorted(str(item["shape_sha256"]) for item in row["topology_signatures"])
    )


def weighted_signature_set(row: Mapping[str, Any]) -> tuple[tuple[str, str], ...]:
    return tuple(
        sorted(
            (str(item["shape_sha256"]), str(item["tree_count"]))
            for item in row["topology_signatures"]
        )
    )


def coarse_set(row: Mapping[str, Any]) -> tuple[str, ...]:
    return tuple(sorted(str(key) for key in row["coarse_class_tree_counts"]))


def weighted_coarse_set(row: Mapping[str, Any]) -> tuple[tuple[str, str], ...]:
    return tuple(
        sorted(
            (str(key), str(value))
            for key, value in row["coarse_class_tree_counts"].items()
        )
    )


def counter_equal(left: Iterable[Any], right: Iterable[Any]) -> bool:
    return Counter(left) == Counter(right)


def rate_record(numerator: int, denominator: int) -> dict[str, Any]:
    return {
        "numerator": numerator,
        "denominator": denominator,
        "rate": None if denominator == 0 else numerator / denominator,
    }


def source_paths(sample: str) -> dict[str, Path]:
    sample_root = TOPOLOGY_ROOT / "samples" / sample
    return {
        "pipeline_receipt": sample_root / "pipeline_receipt.json",
        "mlhp": sample_root / f"{sample}.exact_ps_mlhp.json",
        "mlhp_receipt": sample_root / f"{sample}.exact_ps_mlhp.json.receipt.json",
        "topology": sample_root / f"{sample}.topology.jsonl",
        "topology_receipt": sample_root / f"{sample}.topology.receipt.json",
        "census": CENSUS_ROOT / f"{sample}.census.jsonl",
    }


def validate_authority() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    cohort_receipt_path = TOPOLOGY_ROOT / "cohort_receipt.json"
    cohort_receipt = load_json(cohort_receipt_path)
    require(
        cohort_receipt.get("schema_name")
        == "intersubmod.exact_ps_cpp_topology_cohort.cohort_receipt"
        and cohort_receipt.get("schema_version") == "1.0.0",
        "unexpected primary cohort receipt schema",
    )
    require(
        cohort_receipt.get("all_pass") is True
        and cohort_receipt.get("technical_all_pass") is True,
        "primary cohort technical gate did not pass",
    )
    primary_scope = cohort_receipt.get("scope") or {}
    require(primary_scope.get("autosomes_complete") is True, "autosomes incomplete")
    require(primary_scope.get("canonical_all7") is True, "not canonical all-seven")
    require(
        tuple(primary_scope.get("samples") or ()) == EXPECTED_SAMPLES,
        "primary sample scope/order mismatch",
    )
    require(
        tuple(primary_scope.get("chromosomes") or ()) == EXPECTED_CHROMS,
        "primary chromosome scope mismatch",
    )

    primary_summary_path = TOPOLOGY_ROOT / "summary" / "all7_summary.json"
    primary_summary = load_json(primary_summary_path)
    require(
        primary_summary.get("schema_name")
        == "intersubmod.exact_ps_cpp_topology_af.cohort_summary"
        and primary_summary.get("schema_version") == "1.0.0",
        "unexpected primary summary schema",
    )
    require(
        primary_summary.get("all_pass") is True
        and primary_summary.get("technical_all_pass") is True,
        "primary summary technical gate did not pass",
    )
    all_true(primary_summary.get("checks") or {}, "primary_summary.checks")
    verify_identity(
        primary_summary.get("cohort_receipt_identity") or {},
        "primary_summary.cohort_receipt",
    )

    census_receipt_path = CENSUS_ROOT / "receipt.v2.json"
    census_receipt = load_json(census_receipt_path)
    require(
        census_receipt.get("schema_name")
        == "intersubmod.exact_ps_cpp_topology_signature_census.receipt"
        and census_receipt.get("schema_version") == "1.0.0",
        "unexpected census receipt schema",
    )
    require(
        census_receipt.get("task_type") == TASK_TYPE,
        "census receipt is not Task B",
    )
    all_true(census_receipt.get("checks") or {}, "census_receipt.checks")
    census_scope = census_receipt.get("scope") or {}
    require(
        tuple(census_scope.get("datasets") or ()) == EXPECTED_SAMPLES,
        "census sample scope/order mismatch",
    )
    require(
        census_scope.get("chromosomes") == "chr1-22 via frozen MLHP inputs",
        "census chromosome scope mismatch",
    )

    census_summary_path = verify_identity(
        census_receipt.get("summary") or {}, "census_receipt.summary"
    )
    census_summary = load_json(census_summary_path)
    require(
        census_summary.get("schema_name")
        == "intersubmod.exact_ps_cpp_topology_signature_census.summary"
        and census_summary.get("schema_version") == "1.0.0",
        "unexpected census summary schema",
    )
    all_true(
        (census_summary.get("cohort") or {}).get("checks") or {},
        "census_summary.cohort.checks",
    )
    return primary_summary, census_summary, {
        "cohort_receipt": cohort_receipt,
        "cohort_receipt_path": cohort_receipt_path,
        "primary_summary_path": primary_summary_path,
        "census_receipt": census_receipt,
        "census_receipt_path": census_receipt_path,
        "census_summary_path": census_summary_path,
    }


def validate_sample_source_identities(
    sample: str, authority: Mapping[str, Any]
) -> dict[str, Any]:
    paths = source_paths(sample)
    cohort_receipt = authority["cohort_receipt"]
    receipt_specs = (cohort_receipt.get("sample_receipts") or {}).get(sample)
    require(isinstance(receipt_specs, Mapping), f"{sample}: missing receipt specs")
    verify_identity(receipt_specs.get("partition") or {}, f"{sample}.partition")
    pipeline_path = verify_identity(
        receipt_specs.get("pipeline") or {}, f"{sample}.pipeline"
    )
    require(
        pipeline_path.resolve() == paths["pipeline_receipt"].resolve(),
        f"{sample}: unexpected pipeline path",
    )
    pipeline = load_json(pipeline_path)
    require(
        pipeline.get("schema_name")
        == "intersubmod.exact_ps_cpp_topology_cohort.sample_pipeline_receipt"
        and pipeline.get("schema_version") == "1.0.0",
        f"{sample}: pipeline schema drift",
    )
    require(
        pipeline.get("sample") == sample
        and pipeline.get("all_pass") is True
        and pipeline.get("technical_all_pass") is True,
        f"{sample}: pipeline identity/PASS failure",
    )
    bound_by_path = {
        Path(str(item.get("path", ""))).resolve(): item
        for item in pipeline.get("bound_artifacts") or []
        if isinstance(item, Mapping)
    }
    for role in ("mlhp", "mlhp_receipt", "topology", "topology_receipt"):
        expected = paths[role].resolve()
        require(expected in bound_by_path, f"{sample}: pipeline does not bind {role}")
        verify_identity(bound_by_path[expected], f"{sample}.{role}")

    mlhp_receipt = load_json(paths["mlhp_receipt"])
    require(
        mlhp_receipt.get("schema_name")
        == "intersubmod.exact_ps_layered_tree_input.receipt"
        and mlhp_receipt.get("schema_version") == "1.0.0"
        and mlhp_receipt.get("all_pass") is True,
        f"{sample}: MLHP receipt schema/PASS failure",
    )
    funnel = mlhp_receipt.get("funnel") or {}
    require(
        funnel.get("check_cross_ps_zero") is True
        and funnel.get("check_cross_hp_zero") is True,
        f"{sample}: MLHP cross-PS/HP gate failed",
    )
    require(
        int((mlhp_receipt.get("counts") or {}).get("python_cpp_mismatch_count", -1))
        == 0,
        f"{sample}: MLHP Python/C++ mismatch",
    )
    topology_receipt = load_json(paths["topology_receipt"])
    require(
        topology_receipt.get("schema_name")
        == "intersubmod.exact_ps_cpp_topology_af.receipt"
        and topology_receipt.get("schema_version") == "1.0.0"
        and topology_receipt.get("all_pass") is True
        and topology_receipt.get("all_units_well_formed") is True,
        f"{sample}: topology receipt schema/PASS failure",
    )

    census_spec = (
        (authority["census_receipt"].get("sample_outputs") or {})
        .get(sample, {})
        .get("census")
        or {}
    )
    census_path = verify_identity(census_spec, f"{sample}.census")
    require(
        census_path.resolve() == paths["census"].resolve(),
        f"{sample}: census path drift",
    )

    def frozen(spec: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "path": str(Path(str(spec["path"])).resolve()),
            "size_bytes": int(spec["size_bytes"]),
            "sha256": str(spec["sha256"]),
        }

    return {
        "pipeline_receipt": frozen(receipt_specs["pipeline"]),
        "mlhp": frozen(bound_by_path[paths["mlhp"].resolve()]),
        "mlhp_receipt": frozen(bound_by_path[paths["mlhp_receipt"].resolve()]),
        "topology": frozen(bound_by_path[paths["topology"].resolve()]),
        "topology_receipt": frozen(
            bound_by_path[paths["topology_receipt"].resolve()]
        ),
        "census": frozen(census_spec),
    }


def raw_sample_profile(
    sample: str,
    primary_record: Mapping[str, Any],
    census_record: Mapping[str, Any],
) -> tuple[
    dict[str, Any],
    dict[str, dict[str, Counter[str]]],
    dict[str, Any],
]:
    """Read and validate one sample, returning profile and HCC-only details."""

    paths = source_paths(sample)
    topology_rows: dict[str, dict[str, Any]] = {}
    per_chrom = {
        chrom: {
            "status": Counter(),
            "resolution": Counter(),
            "morphology": Counter(),
            "active_k": Counter(),
            "hp": Counter(),
        }
        for chrom in EXPECTED_CHROMS
    }
    status_counts: Counter[str] = Counter()
    active_k_counts: Counter[str] = Counter()
    hp_counts: Counter[str] = Counter()

    for expected_index, row in enumerate(iter_jsonl(paths["topology"])):
        require(
            row.get("schema_name") == "intersubmod.exact_ps_cpp_topology_af.unit"
            and row.get("schema_version") == "1.0.0",
            f"{sample}: topology row schema drift",
        )
        require(row.get("sample") == sample, f"{sample}: topology sample mismatch")
        require(
            row.get("group_index") == expected_index,
            f"{sample}: non-contiguous topology group_index {expected_index}",
        )
        chrom = str(row.get("chrom"))
        require(chrom in per_chrom, f"{sample}: out-of-scope chrom {chrom}")
        region_id = str(row.get("region_id") or "")
        require(region_id and region_id not in topology_rows, f"{sample}: duplicate region")
        status = str(row.get("unit_status"))
        require(status in STATUS_BINS, f"{sample}: unknown status {status}")
        active_k = int(row.get("active_bit_count", -1))
        require(0 <= active_k <= 12, f"{sample}: active_k outside [0,12]")
        hp_family = str(row.get("hp_family"))
        require(hp_family in HP_BINS, f"{sample}: non-primary HP {hp_family}")
        active_positions = tuple(int(value) for value in row.get("active_positions") or ())
        require(
            tuple(sorted(active_positions)) == active_positions
            and len(active_positions) == active_k,
            f"{sample}: active-position contract drift",
        )
        compact = {
            "sample": sample,
            "group_index": expected_index,
            "region_id": region_id,
            "unit_id": str(row.get("unit_id")),
            "block_id": str(row.get("block_id")),
            "chrom": chrom,
            "phase_set": str(row.get("phase_set")),
            "hp_family": hp_family,
            "unit_status": status,
            "active_positions": active_positions,
            "active_bit_count": active_k,
            "representative_best_edges": row.get("representative_best_edges"),
        }
        topology_rows[region_id] = compact
        status_counts[status] += 1
        active_k_counts[str(active_k)] += 1
        hp_counts[hp_family] += 1
        per_chrom[chrom]["status"][status] += 1
        per_chrom[chrom]["active_k"][str(active_k)] += 1
        per_chrom[chrom]["hp"][hp_family] += 1

    expected_groups = int(primary_record.get("groups_total", -1))
    require(len(topology_rows) == expected_groups, f"{sample}: topology row count")
    require(
        dict(status_counts)
        == {
            key: int(value)
            for key, value in (primary_record.get("status_census") or {})
            .get("unit_status", {})
            .items()
        },
        f"{sample}: raw/summary status mismatch",
    )
    require(
        {key: int(active_k_counts.get(key, 0)) for key in ACTIVE_K_BINS}
        == {
            key: int((primary_record.get("active_k_distribution") or {}).get(key, 0))
            for key in ACTIVE_K_BINS
        },
        f"{sample}: raw/summary active-k mismatch",
    )

    census_rows: dict[str, dict[str, Any]] = {}
    resolution_counts: Counter[str] = Counter()
    morphology_counts: Counter[str] = Counter()
    for row in iter_jsonl(paths["census"]):
        require(
            row.get("schema_name")
            == "intersubmod.exact_ps_cpp_topology_signature_census.unit"
            and row.get("schema_version") == "1.0.0",
            f"{sample}: census row schema drift",
        )
        require(row.get("sample") == sample, f"{sample}: census sample mismatch")
        region_id = str(row.get("region_id") or "")
        require(
            region_id in topology_rows
            and topology_rows[region_id]["unit_status"] == "ranked",
            f"{sample}: census does not join one ranked topology row",
        )
        require(region_id not in census_rows, f"{sample}: duplicate census region")
        require(
            row.get("canonical_reproduction_pass") is True,
            f"{sample}: census canonical reproduction failure",
        )
        resolution = str(row.get("resolution_class"))
        require(resolution in RESOLUTION_BINS, f"{sample}: unknown resolution")
        signatures = row.get("topology_signatures") or []
        coarse_counts = row.get("coarse_class_tree_counts") or {}
        require(
            sum(int(item["tree_count"]) for item in signatures)
            == int(row["best_tree_tie_count"])
            == sum(int(value) for value in coarse_counts.values()),
            f"{sample}: census tree-count conservation failed",
        )
        if int(row.get("coarse_class_count", -1)) == 1:
            require(len(coarse_counts) == 1, f"{sample}: one-coarse mismatch")
            morphology = str(next(iter(coarse_counts)))
            require(morphology in MORPHOLOGY_BINS[:-1], f"{sample}: coarse class")
        else:
            require(
                int(row.get("coarse_class_count", -1)) > 1,
                f"{sample}: invalid coarse_class_count",
            )
            morphology = "Unresolved-cross-coarse"
        chrom = topology_rows[region_id]["chrom"]
        resolution_counts[resolution] += 1
        morphology_counts[morphology] += 1
        per_chrom[chrom]["resolution"][resolution] += 1
        per_chrom[chrom]["morphology"][morphology] += 1
        census_rows[region_id] = {
            "resolution_class": resolution,
            "topology_signatures": signatures,
            "coarse_class_tree_counts": coarse_counts,
        }

    ranked_regions = {
        region_id
        for region_id, row in topology_rows.items()
        if row["unit_status"] == "ranked"
    }
    require(set(census_rows) == ranked_regions, f"{sample}: ranked/census key mismatch")
    ranked = len(census_rows)
    require(ranked == int(census_record["ranked_units"]), f"{sample}: ranked count")
    expected_resolution = {
        key: int((census_record.get("resolution") or {}).get(key, {}).get("n", 0))
        for key in RESOLUTION_BINS
    }
    require(
        {key: int(resolution_counts.get(key, 0)) for key in RESOLUTION_BINS}
        == expected_resolution,
        f"{sample}: raw/summary resolution mismatch",
    )
    expected_morphology = {
        key: int(
            (census_record.get("resolved_coarse_class") or {})
            .get(key, {})
            .get("n", 0)
        )
        for key in MORPHOLOGY_BINS[:-1]
    }
    expected_morphology["Unresolved-cross-coarse"] = int(
        (census_record.get("cross_coarse_class") or {}).get("n", 0)
    )
    require(
        {key: int(morphology_counts.get(key, 0)) for key in MORPHOLOGY_BINS}
        == expected_morphology,
        f"{sample}: raw/summary morphology mismatch",
    )

    mlhp = load_json(paths["mlhp"])
    require(
        mlhp.get("schema_name") == "intersubmod.exact_ps_layered_tree_input"
        and mlhp.get("schema_version") == "1.0.0"
        and mlhp.get("sample") == sample,
        f"{sample}: MLHP schema/sample mismatch",
    )
    require(
        tuple((mlhp.get("scope") or {}).get("chromosomes") or ()) == EXPECTED_CHROMS,
        f"{sample}: MLHP chromosome scope mismatch",
    )
    groups = mlhp.get("groups") or []
    require(len(groups) == expected_groups, f"{sample}: MLHP/topology row count")
    loci: defaultdict[tuple[str, tuple[int, ...]], list[str]] = defaultdict(list)
    tree_input_sites: set[tuple[str, int]] = set()
    for index, group in enumerate(groups):
        region_id = str(group.get("region_id") or "")
        require(region_id in topology_rows, f"{sample}: MLHP region missing topology")
        topology = topology_rows[region_id]
        for field in ("unit_id", "block_id", "chrom", "phase_set", "hp_family"):
            require(
                str(group.get(field)) == topology[field],
                f"{sample}: MLHP/topology {field} mismatch at {index}",
            )
        require(
            topology["group_index"] == index,
            f"{sample}: MLHP array ordinal/group_index mismatch",
        )
        positions = tuple(int(value) for value in group.get("positions") or ())
        require(
            len(positions) == int(group.get("n_sSNV", -1)),
            f"{sample}: MLHP positions/n_sSNV mismatch",
        )
        key = ordered_locus_key(group)
        loci[key].append(region_id)
        tree_input_sites.update((key[0], position) for position in key[1])

    group_count = len(topology_rows)
    status_dimension = dimension_record(
        status_counts, STATUS_BINS, "all_groups", group_count, f"{sample}.status"
    )
    resolution_dimension = dimension_record(
        resolution_counts,
        RESOLUTION_BINS,
        "ranked_units",
        ranked,
        f"{sample}.resolution",
    )
    morphology_dimension = dimension_record(
        morphology_counts,
        MORPHOLOGY_BINS,
        "ranked_units",
        ranked,
        f"{sample}.morphology",
    )
    active_k_dimension = dimension_record(
        active_k_counts,
        ACTIVE_K_BINS,
        "all_groups",
        group_count,
        f"{sample}.active_k",
    )
    hp_dimension = dimension_record(
        hp_counts, HP_BINS, "all_groups", group_count, f"{sample}.hp"
    )
    raw_hp = [hp_dimension["proportions"][key] for key in HP_BINS]
    hp_dimension["comparison_transform"] = (
        "sort HP1/HP2 proportions as minor/major; HP labels are not aligned "
        "across datasets"
    )
    hp_dimension["comparison_bins"] = ["minor_hp_share", "major_hp_share"]
    hp_dimension["comparison_vector"] = sorted(raw_hp)

    output_profile = {
        "sample": sample,
        "denominators": {
            "all_groups": group_count,
            "ranked_units": ranked,
        },
        "dimensions": {
            "status": status_dimension,
            "resolution": resolution_dimension,
            "morphology": morphology_dimension,
            "active_k": active_k_dimension,
            "hp": hp_dimension,
        },
    }
    details = {
        "topology": topology_rows,
        "census": census_rows,
        "loci": loci,
        "tree_input_sites": tree_input_sites,
        "active_sites": {
            (row["chrom"], position)
            for row in topology_rows.values()
            for position in row["active_positions"]
        },
    }
    return output_profile, per_chrom, details


def build_pairwise_profiles(
    profiles: Mapping[str, Mapping[str, Any]]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    vectors = {sample: composite_vector(profiles[sample]) for sample in EXPECTED_SAMPLES}
    rows: list[dict[str, Any]] = []
    for left, right in combinations(EXPECTED_SAMPLES, 2):
        components: dict[str, Any] = {}
        for dimension in DIMENSION_ORDER:
            if dimension == "hp_symmetric":
                left_vector = profiles[left]["dimensions"]["hp"]["comparison_vector"]
                right_vector = profiles[right]["dimensions"]["hp"]["comparison_vector"]
            else:
                left_record = profiles[left]["dimensions"][dimension]
                right_record = profiles[right]["dimensions"][dimension]
                left_vector = [
                    left_record["proportions"][key] for key in left_record["bins"]
                ]
                right_vector = [
                    right_record["proportions"][key] for key in right_record["bins"]
                ]
            similarity = jensen_shannon_similarity(left_vector, right_vector)
            components[dimension] = {
                "similarity": similarity,
                "distance": 1.0 - similarity,
            }
        similarity = jensen_shannon_similarity(vectors[left], vectors[right])
        rows.append(
            {
                "pair_id": f"{left}__{right}",
                "sample_a": left,
                "sample_b": right,
                "relationship": (
                    "same_biological_sample_technical_pair"
                    if (left, right) == HCC_PAIR
                    else "cross_biological_sample_distribution_only"
                ),
                "profile_similarity": similarity,
                "profile_distance": 1.0 - similarity,
                "components": components,
            }
        )
    ranked = sorted(
        rows, key=lambda row: (-row["profile_similarity"], row["pair_id"])
    )
    for rank, row in enumerate(ranked, 1):
        row["similarity_rank_of_21"] = rank

    by_pair = {
        frozenset((row["sample_a"], row["sample_b"])): row for row in rows
    }
    similarity_matrix: list[list[float]] = []
    distance_matrix: list[list[float]] = []
    for left in EXPECTED_SAMPLES:
        similarity_line: list[float] = []
        distance_line: list[float] = []
        for right in EXPECTED_SAMPLES:
            if left == right:
                similarity = 1.0
            else:
                similarity = by_pair[frozenset((left, right))]["profile_similarity"]
            similarity_line.append(similarity)
            distance_line.append(1.0 - similarity)
        similarity_matrix.append(similarity_line)
        distance_matrix.append(distance_line)
    matrix = {
        "sample_order": list(EXPECTED_SAMPLES),
        "similarity": similarity_matrix,
        "distance": distance_matrix,
    }
    return rows, matrix


def aggregate_bootstrap_profile(
    sample: str,
    draw: Sequence[str],
    per_chrom: Mapping[str, Mapping[str, Mapping[str, Counter[str]]]],
) -> list[float]:
    aggregate = {
        dimension: Counter()
        for dimension in ("status", "resolution", "morphology", "active_k", "hp")
    }
    for chrom in draw:
        for dimension in aggregate:
            aggregate[dimension].update(per_chrom[sample][chrom][dimension])

    dimension_specs = (
        ("status", STATUS_BINS),
        ("resolution", RESOLUTION_BINS),
        ("morphology", MORPHOLOGY_BINS),
        ("active_k", ACTIVE_K_BINS),
    )
    values: list[list[float]] = []
    for dimension, bins in dimension_specs:
        denominator = sum(aggregate[dimension].values())
        values.append(
            normalized_counts(
                aggregate[dimension], bins, denominator, f"bootstrap.{sample}.{dimension}"
            )
        )
    hp_denominator = sum(aggregate["hp"].values())
    values.append(
        sorted(
            normalized_counts(
                aggregate["hp"], HP_BINS, hp_denominator, f"bootstrap.{sample}.hp"
            )
        )
    )
    vector = [
        DIMENSION_WEIGHT * value for dimension in values for value in dimension
    ]
    require(math.isclose(sum(vector), 1.0, abs_tol=1e-12), "bootstrap vector")
    return vector


def quantile_nearest_round(values: Sequence[float | int], probability: float) -> Any:
    require(values and 0.0 <= probability <= 1.0, "invalid quantile request")
    ordered = sorted(values)
    return ordered[round((len(ordered) - 1) * probability)]


def chromosome_block_bootstrap(
    pairwise: Sequence[Mapping[str, Any]],
    per_chrom: Mapping[str, Mapping[str, Mapping[str, Counter[str]]]],
) -> dict[str, Any]:
    rng = random.Random(BOOTSTRAP_SEED)
    similarities: list[float] = []
    ranks: list[int] = []
    gaps: list[float] = []
    pair_keys = list(combinations(EXPECTED_SAMPLES, 2))
    for _ in range(BOOTSTRAP_REPLICATES):
        draw = [rng.choice(EXPECTED_CHROMS) for _ in EXPECTED_CHROMS]
        vectors = {
            sample: aggregate_bootstrap_profile(sample, draw, per_chrom)
            for sample in EXPECTED_SAMPLES
        }
        replicate_rows = sorted(
            (
                (jensen_shannon_similarity(vectors[left], vectors[right]), left, right)
                for left, right in pair_keys
            ),
            reverse=True,
        )
        hcc_index = next(
            index
            for index, (_, left, right) in enumerate(replicate_rows)
            if (left, right) == HCC_PAIR
        )
        hcc_similarity = replicate_rows[hcc_index][0]
        best_other = max(
            value
            for value, left, right in replicate_rows
            if (left, right) != HCC_PAIR
        )
        similarities.append(hcc_similarity)
        ranks.append(hcc_index + 1)
        gaps.append(hcc_similarity - best_other)

    point = next(
        row for row in pairwise if (row["sample_a"], row["sample_b"]) == HCC_PAIR
    )
    rank_counts = Counter(ranks)
    return {
        "pair_id": point["pair_id"],
        "method": (
            "resample 22 autosomes with replacement; aggregate raw category "
            "counts; recompute five-dimension equal-weight similarity"
        ),
        "seed": BOOTSTRAP_SEED,
        "replicates": BOOTSTRAP_REPLICATES,
        "quantile_contract": (
            "median=statistics.median; p2.5/p97.5="
            "sorted_values[round((n-1)*p)]"
        ),
        "point_similarity": point["profile_similarity"],
        "similarity": {
            "median": statistics.median(similarities),
            "p2_5": quantile_nearest_round(similarities, 0.025),
            "p97_5": quantile_nearest_round(similarities, 0.975),
        },
        "rank_of_21": {
            "point": point["similarity_rank_of_21"],
            "median": statistics.median(ranks),
            "p2_5": quantile_nearest_round(ranks, 0.025),
            "p97_5": quantile_nearest_round(ranks, 0.975),
            "rank_counts": {str(rank): rank_counts[rank] for rank in sorted(rank_counts)},
            "rank_1_rate": rank_counts[1] / BOOTSTRAP_REPLICATES,
            "top_3_rate": sum(value for rank, value in rank_counts.items() if rank <= 3)
            / BOOTSTRAP_REPLICATES,
        },
        "gap_vs_best_other_pair": {
            "median": statistics.median(gaps),
            "p2_5": quantile_nearest_round(gaps, 0.025),
            "p97_5": quantile_nearest_round(gaps, 0.975),
            "positive_rate": sum(value > 0.0 for value in gaps)
            / BOOTSTRAP_REPLICATES,
        },
    }


def hcc_fingerprint(
    region_id: str,
    topology: Mapping[str, Mapping[str, Any]],
    census: Mapping[str, Mapping[str, Any]],
    level: str,
) -> Any:
    row = topology[region_id]
    status = row["unit_status"]
    active = row["active_positions"]
    if level == "status":
        return status
    if level == "active_status":
        return status, active
    require(status == "ranked", f"{region_id}: ranked fingerprint on {status}")
    signature = census[region_id]
    if level == "active":
        return active
    if level == "resolution":
        return active, signature["resolution_class"]
    if level == "shape":
        return active, signature_set(signature)
    if level == "weighted_shape":
        return active, weighted_signature_set(signature)
    if level == "coarse":
        return active, coarse_set(signature)
    if level == "weighted_coarse":
        return active, weighted_coarse_set(signature)
    if level == "selected":
        require(
            signature["resolution_class"] == "UNIQUE_TREE",
            f"{region_id}: selected fingerprint on tied unit",
        )
        return active, normalized_selected_edges(row)
    raise SimilarityError(f"unknown HCC fingerprint level: {level}")


def build_hcc_comparison(
    details: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    left_name, right_name = HCC_PAIR
    left = details[left_name]
    right = details[right_name]
    left_loci = left["loci"]
    right_loci = right["loci"]
    common_keys = set(left_loci) & set(right_loci)
    union_keys = set(left_loci) | set(right_loci)
    unambiguous, ambiguous = partition_locus_keys(left_loci, right_loci)
    require(
        not {
            (item["chrom"], tuple(item["positions"])) for item in ambiguous
        }
        & {key for key, _, _ in unambiguous},
        "ambiguous HCC key leaked into pair records",
    )

    strict = Counter()
    status_pairs: Counter[str] = Counter()
    pair_records: list[dict[str, Any]] = []
    for key, left_id, right_id in unambiguous:
        left_topology = left["topology"][left_id]
        right_topology = right["topology"][right_id]
        left_status = left_topology["unit_status"]
        right_status = right_topology["unit_status"]
        strict["paired_units"] += 1
        status_pairs[f"{left_status}|{right_status}"] += 1
        same_status = left_status == right_status
        strict["same_status"] += int(same_status)
        record: dict[str, Any] = {
            "chrom": key[0],
            "ordered_positions": list(key[1]),
            "sample_a_region_id": left_id,
            "sample_b_region_id": right_id,
            "sample_a_status": left_status,
            "sample_b_status": right_status,
            "same_status": same_status,
            "both_ranked": False,
            "same_active_positions": None,
            "same_resolution_class": None,
            "same_shape_signature_set": None,
            "same_weighted_shape_signature_set": None,
            "same_coarse_class_set": None,
            "same_weighted_coarse_counts": None,
            "both_unique_tree": None,
            "same_selected_labeled_tree": None,
            "same_selected_unlabeled_shape": None,
            "not_comparable_reason": None,
        }
        if left_status != "ranked" or right_status != "ranked":
            record["not_comparable_reason"] = "one_or_both_units_not_ranked"
            pair_records.append(record)
            continue
        strict["both_ranked"] += 1
        record["both_ranked"] = True
        same_active = (
            left_topology["active_positions"] == right_topology["active_positions"]
        )
        record["same_active_positions"] = same_active
        if not same_active:
            strict["both_ranked_different_active_loci"] += 1
            record["not_comparable_reason"] = "ranked_active_locus_sets_differ"
            pair_records.append(record)
            continue
        strict["both_ranked_same_active_loci"] += 1
        left_census = left["census"][left_id]
        right_census = right["census"][right_id]
        comparisons = {
            "same_resolution_class": (
                left_census["resolution_class"] == right_census["resolution_class"]
            ),
            "same_shape_signature_set": (
                signature_set(left_census) == signature_set(right_census)
            ),
            "same_weighted_shape_signature_set": (
                weighted_signature_set(left_census)
                == weighted_signature_set(right_census)
            ),
            "same_coarse_class_set": (
                coarse_set(left_census) == coarse_set(right_census)
            ),
            "same_weighted_coarse_counts": (
                weighted_coarse_set(left_census)
                == weighted_coarse_set(right_census)
            ),
        }
        for name, matches in comparisons.items():
            record[name] = matches
            strict[name] += int(matches)
        both_unique = (
            left_census["resolution_class"] == "UNIQUE_TREE"
            and right_census["resolution_class"] == "UNIQUE_TREE"
        )
        record["both_unique_tree"] = both_unique
        if both_unique:
            strict["both_unique_tree"] += 1
            same_tree = normalized_selected_edges(
                left_topology
            ) == normalized_selected_edges(right_topology)
            same_shape = signature_set(left_census) == signature_set(right_census)
            record["same_selected_labeled_tree"] = same_tree
            record["same_selected_unlabeled_shape"] = same_shape
            strict["same_selected_labeled_tree"] += int(same_tree)
            strict["same_selected_unlabeled_shape"] += int(same_shape)
        else:
            record["not_comparable_reason"] = (
                "selected_tree_requires_both_unique_tree"
            )
        pair_records.append(record)

    multiplicity_pairs: Counter[str] = Counter()
    multiset = Counter()
    for key in common_keys:
        left_ids = left_loci[key]
        right_ids = right_loci[key]
        multiplicity_pairs[f"{len(left_ids)}|{len(right_ids)}"] += 1
        if len(left_ids) != len(right_ids):
            continue
        multiset["same_multiplicity_keys"] += 1
        if counter_equal(
            (
                hcc_fingerprint(region_id, left["topology"], left["census"], "status")
                for region_id in left_ids
            ),
            (
                hcc_fingerprint(
                    region_id, right["topology"], right["census"], "status"
                )
                for region_id in right_ids
            ),
        ):
            multiset["same_status_multiset"] += 1
        if counter_equal(
            (
                hcc_fingerprint(
                    region_id, left["topology"], left["census"], "active_status"
                )
                for region_id in left_ids
            ),
            (
                hcc_fingerprint(
                    region_id, right["topology"], right["census"], "active_status"
                )
                for region_id in right_ids
            ),
        ):
            multiset["same_active_status_multiset"] += 1

        all_ranked = all(
            left["topology"][region_id]["unit_status"] == "ranked"
            for region_id in left_ids
        ) and all(
            right["topology"][region_id]["unit_status"] == "ranked"
            for region_id in right_ids
        )
        if all_ranked:
            multiset["all_ranked_same_multiplicity_keys"] += 1
            for level, output_key in (
                ("active", "same_active_locus_multiset"),
                ("resolution", "same_resolution_multiset"),
                ("shape", "same_shape_signature_multiset"),
                ("weighted_shape", "same_weighted_shape_signature_multiset"),
                ("coarse", "same_coarse_class_multiset"),
                ("weighted_coarse", "same_weighted_coarse_count_multiset"),
            ):
                if counter_equal(
                    (
                        hcc_fingerprint(
                            region_id, left["topology"], left["census"], level
                        )
                        for region_id in left_ids
                    ),
                    (
                        hcc_fingerprint(
                            region_id, right["topology"], right["census"], level
                        )
                        for region_id in right_ids
                    ),
                ):
                    multiset[output_key] += 1

        all_unique = all_ranked and all(
            left["census"][region_id]["resolution_class"] == "UNIQUE_TREE"
            for region_id in left_ids
        ) and all(
            right["census"][region_id]["resolution_class"] == "UNIQUE_TREE"
            for region_id in right_ids
        )
        if all_unique:
            multiset["all_unique_same_multiplicity_keys"] += 1
            if counter_equal(
                (
                    hcc_fingerprint(
                        region_id, left["topology"], left["census"], "selected"
                    )
                    for region_id in left_ids
                ),
                (
                    hcc_fingerprint(
                        region_id, right["topology"], right["census"], "selected"
                    )
                    for region_id in right_ids
                ),
            ):
                multiset["same_selected_tree_multiset"] += 1

    left_sites = left["tree_input_sites"]
    right_sites = right["tree_input_sites"]
    left_active = left["active_sites"]
    right_active = right["active_sites"]
    same_active_denominator = strict["both_ranked_same_active_loci"]
    unique_denominator = strict["both_unique_tree"]
    all_ranked_multiset_denominator = multiset[
        "all_ranked_same_multiplicity_keys"
    ]
    all_unique_multiset_denominator = multiset[
        "all_unique_same_multiplicity_keys"
    ]
    return {
        "samples": list(HCC_PAIR),
        "relationship": "same_biological_sample_technical_pair",
        "claim_boundary": (
            "Coordinate-matched exact-PS technical concordance; not proof of "
            "identical cellular clones, biological ancestry, or interchangeable runs."
        ),
        "pairing_contract": {
            "ordered_locus_key": ["chrom", "ordered_original_positions"],
            "unambiguous_unit_pair": (
                "key multiplicity must equal one in each dataset"
            ),
            "repeated_key_policy": (
                "never align HP1/HP2 labels; compare the complete unit multiset"
            ),
            "selected_tree_eligibility": (
                "both units ranked, identical active ordered positions, and "
                "both resolution_class=UNIQUE_TREE"
            ),
            "selected_tree_fingerprint": [
                "parent_vertex",
                "child_vertex",
                "acquired_position",
            ],
        },
        "tree_input_site_overlap": {
            "sample_a": len(left_sites),
            "sample_b": len(right_sites),
            "intersection": len(left_sites & right_sites),
            "union": len(left_sites | right_sites),
            "jaccard": len(left_sites & right_sites) / len(left_sites | right_sites),
        },
        "active_site_overlap": {
            "sample_a": len(left_active),
            "sample_b": len(right_active),
            "intersection": len(left_active & right_active),
            "union": len(left_active | right_active),
            "jaccard": len(left_active & right_active) / len(left_active | right_active),
        },
        "ordered_locus_key_overlap": {
            "sample_a_unique_keys": len(left_loci),
            "sample_b_unique_keys": len(right_loci),
            "intersection": len(common_keys),
            "union": len(union_keys),
            "jaccard": len(common_keys) / len(union_keys),
            "unambiguous_1_to_1_keys": len(unambiguous),
            "ambiguous_common_keys_not_paired": len(ambiguous),
            "same_multiplicity_keys": multiset["same_multiplicity_keys"],
            "multiplicity_pair_counts": {
                key: multiplicity_pairs[key] for key in sorted(multiplicity_pairs)
            },
        },
        "unambiguous_comparison": {
            "pair_count": strict["paired_units"],
            "status_pair_counts": {
                key: status_pairs[key] for key in sorted(status_pairs)
            },
            "same_status": rate_record(strict["same_status"], strict["paired_units"]),
            "both_ranked": rate_record(strict["both_ranked"], strict["paired_units"]),
            "same_active_loci_given_both_ranked": rate_record(
                strict["both_ranked_same_active_loci"], strict["both_ranked"]
            ),
            "different_active_loci_given_both_ranked": rate_record(
                strict["both_ranked_different_active_loci"], strict["both_ranked"]
            ),
            "same_resolution_class": rate_record(
                strict["same_resolution_class"], same_active_denominator
            ),
            "same_shape_signature_set": rate_record(
                strict["same_shape_signature_set"], same_active_denominator
            ),
            "same_weighted_shape_signature_set": rate_record(
                strict["same_weighted_shape_signature_set"],
                same_active_denominator,
            ),
            "same_coarse_class_set": rate_record(
                strict["same_coarse_class_set"], same_active_denominator
            ),
            "same_weighted_coarse_counts": rate_record(
                strict["same_weighted_coarse_counts"], same_active_denominator
            ),
            "both_unique_tree": rate_record(
                strict["both_unique_tree"], same_active_denominator
            ),
            "same_selected_labeled_tree": rate_record(
                strict["same_selected_labeled_tree"], unique_denominator
            ),
            "same_selected_unlabeled_shape": rate_record(
                strict["same_selected_unlabeled_shape"], unique_denominator
            ),
        },
        "unambiguous_pair_records": pair_records,
        "multiset_comparison": {
            "same_multiplicity_keys": multiset["same_multiplicity_keys"],
            "same_status_multiset": rate_record(
                multiset["same_status_multiset"],
                multiset["same_multiplicity_keys"],
            ),
            "same_active_status_multiset": rate_record(
                multiset["same_active_status_multiset"],
                multiset["same_multiplicity_keys"],
            ),
            "all_ranked_same_multiplicity_keys": all_ranked_multiset_denominator,
            "same_active_locus_multiset": rate_record(
                multiset["same_active_locus_multiset"],
                all_ranked_multiset_denominator,
            ),
            "same_resolution_multiset": rate_record(
                multiset["same_resolution_multiset"],
                all_ranked_multiset_denominator,
            ),
            "same_shape_signature_multiset": rate_record(
                multiset["same_shape_signature_multiset"],
                all_ranked_multiset_denominator,
            ),
            "same_weighted_shape_signature_multiset": rate_record(
                multiset["same_weighted_shape_signature_multiset"],
                all_ranked_multiset_denominator,
            ),
            "same_coarse_class_multiset": rate_record(
                multiset["same_coarse_class_multiset"],
                all_ranked_multiset_denominator,
            ),
            "same_weighted_coarse_count_multiset": rate_record(
                multiset["same_weighted_coarse_count_multiset"],
                all_ranked_multiset_denominator,
            ),
            "all_unique_same_multiplicity_keys": all_unique_multiset_denominator,
            "same_selected_tree_multiset": rate_record(
                multiset["same_selected_tree_multiset"],
                all_unique_multiset_denominator,
            ),
        },
        "ambiguous_common_key_examples": ambiguous[:20],
    }


def build() -> dict[str, Any]:
    primary_summary, census_summary, authority = validate_authority()
    primary_records = {
        row["sample"]: row for row in primary_summary.get("samples") or []
    }
    census_records = {
        row["sample"]: row for row in census_summary.get("samples") or []
    }
    require(
        tuple(primary_records) == EXPECTED_SAMPLES,
        "primary summary sample order mismatch",
    )
    require(
        tuple(census_records) == EXPECTED_SAMPLES,
        "census summary sample order mismatch",
    )

    profiles: dict[str, dict[str, Any]] = {}
    per_chrom: dict[str, dict[str, dict[str, Counter[str]]]] = {}
    hcc_details: dict[str, dict[str, Any]] = {}
    sample_identities: dict[str, dict[str, Any]] = {}
    for sample in EXPECTED_SAMPLES:
        sample_identities[sample] = validate_sample_source_identities(
            sample, authority
        )
        profile, chrom_counts, details = raw_sample_profile(
            sample, primary_records[sample], census_records[sample]
        )
        expected_census_rows = int(
            (
                (authority["census_receipt"].get("sample_outputs") or {})
                .get(sample, {})
                .get("row_count", -1)
            )
        )
        require(
            profile["denominators"]["ranked_units"] == expected_census_rows,
            f"{sample}: raw census row count differs from receipt",
        )
        profiles[sample] = profile
        per_chrom[sample] = chrom_counts
        if sample in HCC_PAIR:
            hcc_details[sample] = details

    require(
        sum(profile["denominators"]["all_groups"] for profile in profiles.values())
        == 98955,
        "all-seven group total drift",
    )
    require(
        sum(profile["denominators"]["ranked_units"] for profile in profiles.values())
        == 71955,
        "all-seven ranked total drift",
    )

    pairwise, matrix = build_pairwise_profiles(profiles)
    require(len(pairwise) == 21, "expected exactly 21 unordered pairs")
    bootstrap = chromosome_block_bootstrap(pairwise, per_chrom)
    hcc = build_hcc_comparison(hcc_details)

    for index in range(len(EXPECTED_SAMPLES)):
        require(
            math.isclose(matrix["similarity"][index][index], 1.0)
            and math.isclose(matrix["distance"][index][index], 0.0),
            "matrix diagonal drift",
        )
        for other in range(len(EXPECTED_SAMPLES)):
            require(
                math.isclose(
                    matrix["similarity"][index][other],
                    matrix["similarity"][other][index],
                    abs_tol=1e-15,
                ),
                "similarity matrix is not symmetric",
            )

    source_identities = {
        "producer": identity(Path(__file__)),
        "primary_cohort_receipt": identity(authority["cohort_receipt_path"]),
        "primary_summary": identity(authority["primary_summary_path"]),
        "census_receipt_v2": identity(authority["census_receipt_path"]),
        "census_summary": identity(authority["census_summary_path"]),
        "samples": sample_identities,
    }
    return {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "task_type": TASK_TYPE,
        "source_profile": SOURCE_PROFILE,
        "scope": {
            "datasets": list(EXPECTED_SAMPLES),
            "dataset_count": 7,
            "biological_sample_count": 6,
            "chromosomes": list(EXPECTED_CHROMS),
            "autosomes_complete": True,
            "partial": False,
        },
        "metric_contract": {
            "name": "one_minus_jensen_shannon_distance",
            "formula": "similarity = 1 - sqrt(JS_divergence_base2)",
            "range": [0.0, 1.0],
            "dimension_order": list(DIMENSION_ORDER),
            "dimension_weights": {
                dimension: DIMENSION_WEIGHT for dimension in DIMENSION_ORDER
            },
            "dimensions": {
                "status": {
                    "denominator": "all_groups",
                    "bins": list(STATUS_BINS),
                },
                "resolution": {
                    "denominator": "ranked_units",
                    "bins": list(RESOLUTION_BINS),
                },
                "morphology": {
                    "denominator": "ranked_units",
                    "bins": list(MORPHOLOGY_BINS),
                },
                "active_k": {
                    "denominator": "all_groups",
                    "bins": list(ACTIVE_K_BINS),
                    "definition": "active ALT bit count; not original n_sSNV",
                },
                "hp_symmetric": {
                    "denominator": "all_groups",
                    "bins": ["minor_hp_share", "major_hp_share"],
                    "definition": (
                        "sorted HP1/HP2 shares; direct HP labels are not compared "
                        "between datasets"
                    ),
                },
            },
            "composite_construction": (
                "concatenate each normalized dimension multiplied by 0.2; "
                "category count therefore does not change dimension weight"
            ),
            "interpretation": (
                "distribution-profile similarity only; no locus identity, clone "
                "identity, or biological ancestry claim"
            ),
        },
        "sample_profiles": profiles,
        "pairwise_profiles": pairwise,
        "matrix": matrix,
        "hcc1395_dorado_chromosome_block_bootstrap": bootstrap,
        "hcc1395_dorado_strict_comparison": hcc,
        "source_identities": source_identities,
        "checks": {
            "all_source_identities_verified": True,
            "scope_is_exactly_all7_chr1_22": True,
            "profile_denominators_conserved": True,
            "pair_count_equals_21": len(pairwise) == 21,
            "matrix_symmetric_with_unit_diagonal": True,
            "hcc_repeated_ordered_locus_keys_not_arbitrarily_paired": (
                hcc["ordered_locus_key_overlap"][
                    "ambiguous_common_keys_not_paired"
                ]
                > 0
                and len(hcc["unambiguous_pair_records"])
                == hcc["ordered_locus_key_overlap"]["unambiguous_1_to_1_keys"]
            ),
            "bootstrap_contract_fixed": (
                bootstrap["seed"] == BOOTSTRAP_SEED
                and bootstrap["replicates"] == BOOTSTRAP_REPLICATES
            ),
        },
        "claim_boundary": (
            "All-seven Jensen-Shannon results compare dataset-level category "
            "distributions. Only the HCC1395 technical pair receives coordinate-"
            "matched analysis. Neither result validates identical cellular clones, "
            "biological ancestry, CN/LOH correctness, or calibrated probability."
        ),
    }


def write_new_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Write once and fail closed if the target already exists."""

    path = path.resolve()
    if path.exists():
        raise FileExistsError(f"output already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    document = json.dumps(
        payload, ensure_ascii=False, indent=2, sort_keys=True
    ) + "\n"
    try:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(document)
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError:
        raise FileExistsError(f"output already exists: {path}") from None


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="New JSON sidecar path. Existing targets are never overwritten.",
    )
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if output.exists():
        print(f"ERROR output already exists: {output}", file=sys.stderr)
        return 2

    started = time.perf_counter()
    print(f"INPUT topology_root={TOPOLOGY_ROOT}", flush=True)
    print(f"INPUT census_root={CENSUS_ROOT}", flush=True)
    print(f"COMMAND {Path(sys.executable).resolve()} {Path(__file__).resolve()} "
          f"--output {output}", flush=True)
    try:
        payload = build()
        write_new_json(output, payload)
    except (SimilarityError, FileExistsError, OSError, ValueError) as exc:
        print(f"ERROR {exc}", file=sys.stderr)
        return 2
    elapsed = time.perf_counter() - started
    hcc_pair = next(
        row
        for row in payload["pairwise_profiles"]
        if row["pair_id"] == "HCC1395__HCC1395_DORADO"
    )
    strict = payload["hcc1395_dorado_strict_comparison"][
        "unambiguous_comparison"
    ]
    print(f"OUTPUT {output}", flush=True)
    print(
        "RESULT datasets=7 pairs=21 groups=98955 ranked=71955 "
        f"hcc_similarity={hcc_pair['profile_similarity']:.6f} "
        f"hcc_strict_pairs={strict['pair_count']} "
        f"selected_tree={strict['same_selected_labeled_tree']['numerator']}/"
        f"{strict['same_selected_labeled_tree']['denominator']} "
        f"elapsed_seconds={elapsed:.3f}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
