#!/usr/bin/env python3
"""Run a within-region B-side site-label permutation null.

The input is the compact, machine-readable output produced by
``build_site_topology_containment.py``.  No solver or read data are reopened.
For every region, alignment, and mapped HP component, a uniformly random
permutation of the B-side shared-site labels is applied to every frozen
F/R/P relation signature.  The same component permutation is shared by the
read-full, VAF-official, and VAF-normalized endpoints within an iteration.

The VAF endpoints are conditional on their frozen selected top sets.  VAF
scores are deliberately not recomputed under the permutation.

This null preserves, by construction, component k, projected Q-set size,
pre-projection candidate multiplicity, and HP mapping.  Identity and global
HP-swap alignments are recomputed separately.  The swap-tolerant endpoint
then applies the same whole-region best-of-two rule as the observed analysis:
Jaccard, minimum directional coverage, intersection size, category rank, and
identity as the final tie-break.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import itertools
import json
import math
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


ENDPOINTS = ("read_full", "vaf_official", "vaf_normalized")
MAPPINGS = ("identity", "global_hp_swap")
COMPARISONS = ("identity", "global_hp_swap", "swap_tolerant")
METRICS = ("exact", "robust_nested", "any_compatible", "unique_exact")
K_BINS = ("k=2", "k=3", "k=4", "k=5+")
CATEGORY_RANK = {
    "disjoint": 0,
    "overlap": 1,
    "a_proper_subset_b": 2,
    "b_proper_subset_a": 2,
    "exact": 3,
}


@dataclass
class ComponentPlan:
    k: int
    permutation_count: int
    intersections: dict[str, np.ndarray]


@dataclass
class MappingPlan:
    valid_mapping: bool
    components: list[ComponentPlan]
    endpoint_evaluable: dict[str, bool]
    set_size_a: dict[str, int]
    set_size_b: dict[str, int]
    observed: dict[str, dict]
    k_max: int


@dataclass
class RegionPlan:
    match_id: str
    chrom: str
    mappings: dict[str, MappingPlan]
    observed_swap: dict[str, dict]
    swap_k_bin: dict[str, str]


@dataclass
class GroupStats:
    endpoint: str
    comparison: str
    k_bin: str
    denominator: int = 0
    observed_success: dict[str, int] = field(
        default_factory=lambda: {metric: 0 for metric in METRICS}
    )
    chromosome_denominator: collections.Counter = field(default_factory=collections.Counter)
    chromosome_observed: dict[str, collections.Counter] = field(
        default_factory=lambda: {metric: collections.Counter() for metric in METRICS}
    )
    null_counts: dict[str, np.ndarray] = field(default_factory=dict)
    chromosome_null_success_sum: dict[str, collections.Counter] = field(
        default_factory=lambda: {metric: collections.Counter() for metric in METRICS}
    )


def parse_args() -> argparse.Namespace:
    topic_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=topic_dir / "data" / "hcc1395_site_topology_signature_sets.jsonl.gz",
        help="Compact JSONL.gz signature-set input.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=topic_dir / "data",
        help="Directory for null TSV/JSON/check artifacts.",
    )
    parser.add_argument("--permutations", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260712)
    parser.add_argument("--bootstrap-replicates", type=int, default=5000)
    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=None,
        help="Defaults to seed+1; permutation seed remains exactly --seed.",
    )
    parser.add_argument("--progress-every", type=int, default=1000)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )


def write_tsv(path: Path, rows: Sequence[dict], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_jsonl_gz(path: Path) -> list[dict]:
    rows: list[dict] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.strip():
                try:
                    rows.append(json.loads(line))
                except json.JSONDecodeError as exc:
                    raise RuntimeError(f"invalid JSON at {path}:{line_number}: {exc}") from exc
    return rows


@lru_cache(maxsize=None)
def pair_order(k: int) -> tuple[tuple[int, int], ...]:
    return tuple(itertools.combinations(range(k), 2))


@lru_cache(maxsize=None)
def pair_index(k: int) -> dict[tuple[int, int], int]:
    return {pair: index for index, pair in enumerate(pair_order(k))}


@lru_cache(maxsize=None)
def permutations_for_k(k: int) -> tuple[tuple[int, ...], ...]:
    return tuple(itertools.permutations(range(k)))


def flip_direction(code: str) -> str:
    if code == "F":
        return "R"
    if code == "R":
        return "F"
    if code == "P":
        return "P"
    raise RuntimeError(f"unknown relation code: {code!r}")


@lru_cache(maxsize=None)
def transform_signature(signature: str, k: int, permutation: tuple[int, ...]) -> str:
    """Relabel a relation vector by a vertex permutation.

    ``permutation[new_rank]`` is the old vertex assigned to that new genomic
    rank.  The emitted vector is again ordered by lexicographic new-rank pairs.
    F/R is flipped whenever the queried old-vertex order is reversed.
    """
    expected = math.comb(k, 2)
    if len(signature) != expected:
        raise RuntimeError(
            f"relation vector length mismatch: k={k}, expected={expected}, observed={len(signature)}"
        )
    if tuple(sorted(permutation)) != tuple(range(k)):
        raise RuntimeError(f"not a permutation for k={k}: {permutation}")
    old_index = pair_index(k)
    transformed: list[str] = []
    for new_left, new_right in pair_order(k):
        old_left = permutation[new_left]
        old_right = permutation[new_right]
        if old_left < old_right:
            transformed.append(signature[old_index[(old_left, old_right)]])
        else:
            transformed.append(flip_direction(signature[old_index[(old_right, old_left)]]))
    return "".join(transformed)


@lru_cache(maxsize=None)
def intersection_profile(
    k: int,
    q_a: tuple[str, ...],
    q_b: tuple[str, ...],
) -> tuple[int, ...]:
    """Return |Q_A intersect perm(Q_B)| for every lexicographic permutation."""
    set_a = set(q_a)
    set_b = set(q_b)
    if not set_a or not set_b:
        raise RuntimeError("evaluable component cannot have an empty Q set")
    output: list[int] = []
    for permutation in permutations_for_k(k):
        transformed_b = {transform_signature(value, k, permutation) for value in set_b}
        if len(transformed_b) != len(set_b):
            raise RuntimeError("site relabeling changed Q-set cardinality")
        output.append(len(set_a & transformed_b))
    return tuple(output)


def category_from_counts(size_a: int, size_b: int, intersection: int) -> str:
    if size_a == size_b == intersection:
        return "exact"
    if intersection == size_a and size_a < size_b:
        return "a_proper_subset_b"
    if intersection == size_b and size_b < size_a:
        return "b_proper_subset_a"
    if intersection:
        return "overlap"
    return "disjoint"


def metric_values(category: str, unique_exact: bool) -> dict[str, bool]:
    return {
        "exact": category == "exact",
        "robust_nested": category
        in {"exact", "a_proper_subset_b", "b_proper_subset_a"},
        "any_compatible": category != "disjoint",
        "unique_exact": unique_exact,
    }


def metric_arrays(size_a: int | np.ndarray, size_b: int | np.ndarray, intersection: np.ndarray):
    exact = (intersection == size_a) & (intersection == size_b)
    a_subset = (intersection == size_a) & (np.asarray(size_a) < np.asarray(size_b))
    b_subset = (intersection == size_b) & (np.asarray(size_b) < np.asarray(size_a))
    return {
        "exact": exact,
        "robust_nested": exact | a_subset | b_subset,
        "any_compatible": intersection > 0,
        "unique_exact": exact & (np.asarray(size_a) == 1) & (np.asarray(size_b) == 1),
    }


def k_bin(k: int) -> str:
    if k <= 2:
        return "k=2"
    if k == 3:
        return "k=3"
    if k == 4:
        return "k=4"
    return "k=5+"


def validate_component(component: Mapping, endpoint: str, side: str) -> None:
    k = int(component["shared_k"])
    projection = component[endpoint][side]
    q = projection.get("q", [])
    if len(q) != len(set(q)):
        raise RuntimeError("compact Q set contains duplicates")
    if int(projection.get("projected_topology_set_size", len(q))) != len(q):
        raise RuntimeError("compact projected_topology_set_size does not match Q cardinality")
    for signature in q:
        if len(signature) != math.comb(k, 2):
            raise RuntimeError(f"compact relation length mismatch: k={k}, q={signature!r}")
        if set(signature) - {"F", "R", "P"}:
            raise RuntimeError(f"compact relation alphabet mismatch: {signature!r}")
    if q and int(projection.get("candidate_tree_count", 0)) < len(q):
        raise RuntimeError("candidate multiplicity is smaller than projected Q-set size")


def make_mapping_plan(alignment: Mapping) -> MappingPlan:
    valid = bool(alignment.get("valid_mapping"))
    components: list[ComponentPlan] = []
    if valid:
        for component in alignment.get("components", []):
            k = int(component["shared_k"])
            intersections: dict[str, np.ndarray] = {}
            for endpoint in ENDPOINTS:
                for side in ("a", "b"):
                    validate_component(component, endpoint, side)
                left = component[endpoint]["a"]
                right = component[endpoint]["b"]
                if left.get("status") == "not_evaluable" or right.get("status") == "not_evaluable":
                    continue
                q_a = tuple(left["q"])
                q_b = tuple(right["q"])
                if not q_a or not q_b:
                    continue
                intersections[endpoint] = np.asarray(
                    intersection_profile(k, q_a, q_b), dtype=np.int32
                )
            components.append(
                ComponentPlan(
                    k=k,
                    permutation_count=math.factorial(k),
                    intersections=intersections,
                )
            )

    endpoint_evaluable: dict[str, bool] = {}
    set_size_a: dict[str, int] = {}
    set_size_b: dict[str, int] = {}
    observed: dict[str, dict] = {}
    for endpoint in ENDPOINTS:
        source = dict(alignment.get(endpoint, {}))
        observed[endpoint] = source
        evaluable = valid and source.get("status") == "evaluable"
        endpoint_evaluable[endpoint] = evaluable
        if evaluable:
            set_size_a[endpoint] = int(source["topology_set_size_a"])
            set_size_b[endpoint] = int(source["topology_set_size_b"])
            if any(endpoint not in component.intersections for component in components):
                raise RuntimeError("evaluable mapping endpoint lacks a component intersection profile")
        else:
            set_size_a[endpoint] = 0
            set_size_b[endpoint] = 0
    return MappingPlan(
        valid_mapping=valid,
        components=components,
        endpoint_evaluable=endpoint_evaluable,
        set_size_a=set_size_a,
        set_size_b=set_size_b,
        observed=observed,
        k_max=max((component.k for component in components), default=0),
    )


def scalar_rank(summary: Mapping | None, mapping_name: str) -> tuple:
    if not summary or summary.get("status") != "evaluable":
        return (-1.0, -1.0, -1, -1, -1)
    return (
        float(summary["jaccard"]),
        min(float(summary["coverage_a"]), float(summary["coverage_b"])),
        int(summary["intersection_size"]),
        CATEGORY_RANK[str(summary["category"])],
        1 if mapping_name == "identity" else 0,
    )


def observed_best_mapping(mappings: Mapping[str, MappingPlan], endpoint: str) -> tuple[str, dict] | None:
    candidates = [
        (scalar_rank(plan.observed[endpoint], mapping_name), mapping_name, plan.observed[endpoint])
        for mapping_name, plan in mappings.items()
        if plan.valid_mapping
    ]
    if not candidates:
        return None
    _, mapping_name, summary = max(candidates, key=lambda item: item[0])
    return mapping_name, summary


def build_region_plans(records: Sequence[dict], checks: list[dict]) -> list[RegionPlan]:
    plans: list[RegionPlan] = []
    schema_versions = {str(record.get("schema_version")) for record in records}
    add_check(checks, "input_schema_version", "['1.0']", str(sorted(schema_versions)), schema_versions == {"1.0"})
    observed_mapping_checks = 0
    identity_reconstruction_checks = 0

    for record in records:
        mappings = {
            mapping: make_mapping_plan(record["alignments"][mapping]) for mapping in MAPPINGS
        }
        observed_swap = {endpoint: dict(record["swap_tolerant"][endpoint]) for endpoint in ENDPOINTS}
        swap_k_bins: dict[str, str] = {}

        for mapping_name, mapping in mappings.items():
            for endpoint in ENDPOINTS:
                if not mapping.endpoint_evaluable[endpoint]:
                    continue
                intersection = 1
                for component in mapping.components:
                    intersection *= int(component.intersections[endpoint][0])
                size_a = mapping.set_size_a[endpoint]
                size_b = mapping.set_size_b[endpoint]
                category = category_from_counts(size_a, size_b, intersection)
                source = mapping.observed[endpoint]
                identity_reconstruction_checks += 1
                if (
                    category != source.get("category")
                    or intersection != int(source.get("intersection_size", -1))
                    or (size_a == size_b == intersection == 1) != bool(source.get("unique_exact"))
                ):
                    raise RuntimeError(
                        f"identity relabeling failed to reconstruct observed mapping: "
                        f"{record['match_id']} {mapping_name} {endpoint}"
                    )

        for endpoint in ENDPOINTS:
            best = observed_best_mapping(mappings, endpoint)
            stored = observed_swap[endpoint]
            if stored.get("status") == "evaluable":
                if best is None:
                    raise RuntimeError("stored swap-tolerant result is evaluable without a valid mapping")
                mapping_name, summary = best
                observed_mapping_checks += 1
                if (
                    mapping_name != stored.get("selected_mapping")
                    or summary.get("category") != stored.get("category")
                    or int(summary.get("intersection_size", -1))
                    != int(stored.get("intersection_size", -2))
                ):
                    raise RuntimeError(
                        f"best-of-two reconstruction mismatch: {record['match_id']} {endpoint}"
                    )
                swap_k_bins[endpoint] = k_bin(mappings[mapping_name].k_max)
            else:
                swap_k_bins[endpoint] = ""

        plans.append(
            RegionPlan(
                match_id=str(record["match_id"]),
                chrom=str(record["chrom"]),
                mappings=mappings,
                observed_swap=observed_swap,
                swap_k_bin=swap_k_bins,
            )
        )

    add_check(
        checks,
        "identity_permutation_reconstructs_observed_alignments",
        ">0 exact reconstructions",
        str(identity_reconstruction_checks),
        identity_reconstruction_checks > 0,
    )
    add_check(
        checks,
        "observed_best_of_two_reconstructed",
        ">0 exact reconstructions",
        str(observed_mapping_checks),
        observed_mapping_checks > 0,
    )
    return plans


def add_check(
    checks: list[dict], check: str, expected: str, observed: str, passed: bool
) -> None:
    checks.append(
        {
            "check": check,
            "expected": expected,
            "observed": observed,
            "status": "PASS" if passed else "FAIL",
        }
    )


def get_observed_summary(region: RegionPlan, endpoint: str, comparison: str) -> dict:
    if comparison in MAPPINGS:
        return region.mappings[comparison].observed[endpoint]
    return region.observed_swap[endpoint]


def get_fixed_k_bin(region: RegionPlan, endpoint: str, comparison: str) -> str:
    if comparison in MAPPINGS:
        return k_bin(region.mappings[comparison].k_max)
    return region.swap_k_bin[endpoint]


def initialize_groups(
    regions: Sequence[RegionPlan], permutations: int, checks: list[dict]
) -> dict[tuple[str, str, str], GroupStats]:
    groups: dict[tuple[str, str, str], GroupStats] = {}
    for endpoint in ENDPOINTS:
        for comparison in COMPARISONS:
            for label in ("all", *K_BINS):
                groups[(endpoint, comparison, label)] = GroupStats(
                    endpoint=endpoint,
                    comparison=comparison,
                    k_bin=label,
                    null_counts={
                        metric: np.zeros(permutations, dtype=np.int32) for metric in METRICS
                    },
                )

    for region in regions:
        for endpoint in ENDPOINTS:
            for comparison in COMPARISONS:
                summary = get_observed_summary(region, endpoint, comparison)
                if summary.get("status") != "evaluable":
                    continue
                label = get_fixed_k_bin(region, endpoint, comparison)
                if label not in K_BINS:
                    raise RuntimeError(
                        f"evaluable result lacks a fixed complexity bin: "
                        f"{region.match_id} {endpoint} {comparison}"
                    )
                category = str(summary["category"])
                values = metric_values(category, bool(summary.get("unique_exact")))
                for group_label in ("all", label):
                    group = groups[(endpoint, comparison, group_label)]
                    group.denominator += 1
                    group.chromosome_denominator[region.chrom] += 1
                    for metric, success in values.items():
                        group.observed_success[metric] += int(success)
                        group.chromosome_observed[metric][region.chrom] += int(success)

    for endpoint in ENDPOINTS:
        for comparison in COMPARISONS:
            all_group = groups[(endpoint, comparison, "all")]
            stratum_total = sum(groups[(endpoint, comparison, label)].denominator for label in K_BINS)
            add_check(
                checks,
                f"complexity_denominator_conservation:{endpoint}:{comparison}",
                str(all_group.denominator),
                str(stratum_total),
                all_group.denominator == stratum_total,
            )
    return groups


def mapping_null_arrays(
    plan: MappingPlan,
    endpoint: str,
    permutation_indices: Sequence[np.ndarray],
) -> dict | None:
    if not plan.endpoint_evaluable[endpoint]:
        return None
    if len(permutation_indices) != len(plan.components):
        raise RuntimeError("component permutation index count mismatch")
    permutations = len(permutation_indices[0]) if permutation_indices else 0
    intersection = np.ones(permutations, dtype=np.int64)
    for component, indices in zip(plan.components, permutation_indices):
        intersection *= component.intersections[endpoint][indices]
    size_a = plan.set_size_a[endpoint]
    size_b = plan.set_size_b[endpoint]
    union = size_a + size_b - intersection
    max_size = max(size_a, size_b)
    category_rank = np.zeros(permutations, dtype=np.int8)
    exact = (intersection == size_a) & (intersection == size_b)
    subset = ((intersection == size_a) & (size_a < size_b)) | (
        (intersection == size_b) & (size_b < size_a)
    )
    overlap = (intersection > 0) & ~exact & ~subset
    category_rank[overlap] = 1
    category_rank[subset] = 2
    category_rank[exact] = 3
    return {
        "intersection": intersection,
        "size_a": size_a,
        "size_b": size_b,
        "union": union,
        "max_size": max_size,
        "category_rank": category_rank,
    }


def select_best_arrays(identity: dict | None, swapped: dict | None) -> dict | None:
    if identity is None:
        return swapped
    if swapped is None:
        return identity

    i_inter, s_inter = identity["intersection"], swapped["intersection"]
    i_union, s_union = identity["union"], swapped["union"]
    # Compare Jaccard exactly by cross multiplication.
    i_better = i_inter * s_union > s_inter * i_union
    tied = i_inter * s_union == s_inter * i_union

    # min(coverage_a, coverage_b) = intersection / max(size_a, size_b).
    i_cov_better = i_inter * swapped["max_size"] > s_inter * identity["max_size"]
    cov_tied = i_inter * swapped["max_size"] == s_inter * identity["max_size"]
    i_better |= tied & i_cov_better
    tied &= cov_tied

    i_better |= tied & (i_inter > s_inter)
    tied &= i_inter == s_inter
    i_better |= tied & (identity["category_rank"] > swapped["category_rank"])
    tied &= identity["category_rank"] == swapped["category_rank"]
    # Identity is the final observed tie-break.
    i_better |= tied

    return {
        "intersection": np.where(i_better, i_inter, s_inter),
        "size_a": np.where(i_better, identity["size_a"], swapped["size_a"]),
        "size_b": np.where(i_better, identity["size_b"], swapped["size_b"]),
    }


def add_null_result(
    group: GroupStats,
    region: RegionPlan,
    arrays: dict,
) -> None:
    masks = metric_arrays(arrays["size_a"], arrays["size_b"], arrays["intersection"])
    for metric, mask in masks.items():
        values = np.asarray(mask, dtype=np.int32)
        group.null_counts[metric] += values
        group.chromosome_null_success_sum[metric][region.chrom] += int(values.sum())


def run_null(
    regions: Sequence[RegionPlan],
    groups: dict[tuple[str, str, str], GroupStats],
    permutations: int,
    seed: int,
    progress_every: int,
) -> None:
    rng = np.random.default_rng(seed)
    for region_index, region in enumerate(regions, start=1):
        mapping_indices: dict[str, list[np.ndarray]] = {}
        for mapping_name in MAPPINGS:
            plan = region.mappings[mapping_name]
            mapping_indices[mapping_name] = [
                rng.integers(0, component.permutation_count, size=permutations, dtype=np.int32)
                for component in plan.components
            ]

        for endpoint in ENDPOINTS:
            mapping_arrays = {
                mapping_name: mapping_null_arrays(
                    region.mappings[mapping_name],
                    endpoint,
                    mapping_indices[mapping_name],
                )
                for mapping_name in MAPPINGS
            }
            comparison_arrays = {
                "identity": mapping_arrays["identity"],
                "global_hp_swap": mapping_arrays["global_hp_swap"],
                "swap_tolerant": select_best_arrays(
                    mapping_arrays["identity"], mapping_arrays["global_hp_swap"]
                ),
            }
            for comparison, arrays in comparison_arrays.items():
                observed = get_observed_summary(region, endpoint, comparison)
                if observed.get("status") != "evaluable":
                    continue
                if arrays is None:
                    raise RuntimeError("fixed evaluable population became non-evaluable under null")
                label = get_fixed_k_bin(region, endpoint, comparison)
                for group_label in ("all", label):
                    add_null_result(groups[(endpoint, comparison, group_label)], region, arrays)

        if progress_every and region_index % progress_every == 0:
            print(f"processed_regions={region_index}/{len(regions)}", flush=True)


def quantile(values: np.ndarray, probability: float) -> float:
    if not len(values):
        return float("nan")
    ordered = np.sort(np.asarray(values, dtype=float))
    position = (len(ordered) - 1) * probability
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return float(ordered[lower])
    weight = position - lower
    return float(ordered[lower] * (1.0 - weight) + ordered[upper] * weight)


def bootstrap_and_loco(
    group: GroupStats,
    metric: str,
    chromosomes: Sequence[str],
    bootstrap_indices: np.ndarray,
    permutations: int,
) -> tuple[dict, list[dict]]:
    denominator = np.asarray(
        [group.chromosome_denominator[chrom] for chrom in chromosomes], dtype=float
    )
    observed = np.asarray(
        [group.chromosome_observed[metric][chrom] for chrom in chromosomes], dtype=float
    )
    null_expected = np.asarray(
        [
            group.chromosome_null_success_sum[metric][chrom] / permutations
            for chrom in chromosomes
        ],
        dtype=float,
    )
    sampled_denominator = denominator[bootstrap_indices].sum(axis=1)
    sampled_excess_numerator = (observed - null_expected)[bootstrap_indices].sum(axis=1)
    valid = sampled_denominator > 0
    bootstrap_excess = sampled_excess_numerator[valid] / sampled_denominator[valid]

    total_denominator = denominator.sum()
    total_observed = observed.sum()
    total_null = null_expected.sum()
    loco_rows: list[dict] = []
    for index, chromosome in enumerate(chromosomes):
        held_denominator = total_denominator - denominator[index]
        if held_denominator:
            held_observed_rate = (total_observed - observed[index]) / held_denominator
            held_null_rate = (total_null - null_expected[index]) / held_denominator
            held_excess = held_observed_rate - held_null_rate
        else:
            held_observed_rate = held_null_rate = held_excess = float("nan")
        loco_rows.append(
            {
                "endpoint": group.endpoint,
                "comparison": group.comparison,
                "k_bin": group.k_bin,
                "metric": metric,
                "held_out_chromosome": chromosome,
                "denominator": int(held_denominator),
                "observed_rate": held_observed_rate,
                "null_mean_rate": held_null_rate,
                "observed_minus_null_mean": held_excess,
            }
        )

    finite_loco = np.asarray(
        [row["observed_minus_null_mean"] for row in loco_rows if math.isfinite(row["observed_minus_null_mean"])],
        dtype=float,
    )
    summary = {
        "block_bootstrap_replicates_valid": int(len(bootstrap_excess)),
        "block_bootstrap_excess_ci95_low": quantile(bootstrap_excess, 0.025),
        "block_bootstrap_excess_ci95_high": quantile(bootstrap_excess, 0.975),
        "loco_min_excess": float(finite_loco.min()) if len(finite_loco) else float("nan"),
        "loco_max_excess": float(finite_loco.max()) if len(finite_loco) else float("nan"),
        "loco_positive_chromosomes": int((finite_loco > 0).sum()),
        "loco_chromosomes": int(len(finite_loco)),
    }
    return summary, loco_rows


def format_float(value: float | int) -> str:
    value = float(value)
    return "" if not math.isfinite(value) else f"{value:.12f}"


def summarize_groups(
    groups: Mapping[tuple[str, str, str], GroupStats],
    permutations: int,
    seed: int,
    bootstrap_replicates: int,
    bootstrap_seed: int,
    chromosomes: Sequence[str],
) -> tuple[list[dict], list[dict], list[dict]]:
    bootstrap_rng = np.random.default_rng(bootstrap_seed)
    bootstrap_indices = bootstrap_rng.integers(
        0,
        len(chromosomes),
        size=(bootstrap_replicates, len(chromosomes)),
        dtype=np.int16,
    )
    all_rows: list[dict] = []
    complexity_rows: list[dict] = []
    loco_rows: list[dict] = []
    for key in sorted(groups):
        group = groups[key]
        if group.denominator == 0:
            continue
        for metric in METRICS:
            null_rates = group.null_counts[metric].astype(float) / group.denominator
            observed_count = group.observed_success[metric]
            observed_rate = observed_count / group.denominator
            null_mean = float(null_rates.mean())
            bootstrap, metric_loco = bootstrap_and_loco(
                group,
                metric,
                chromosomes,
                bootstrap_indices,
                permutations,
            )
            loco_rows.extend(metric_loco)
            row = {
                "endpoint": group.endpoint,
                "comparison": group.comparison,
                "population": "fixed_observed_evaluable_regions",
                "k_bin": group.k_bin,
                "metric": metric,
                "denominator": group.denominator,
                "observed_n": observed_count,
                "observed_rate": observed_rate,
                "null_mean_n": float(group.null_counts[metric].mean()),
                "null_mean_rate": null_mean,
                "null_sd_rate": float(null_rates.std(ddof=1)) if permutations > 1 else 0.0,
                "null_q025_rate": quantile(null_rates, 0.025),
                "null_q975_rate": quantile(null_rates, 0.975),
                "observed_minus_null_mean": observed_rate - null_mean,
                "empirical_p_ge": (1 + int((null_rates >= observed_rate).sum()))
                / (permutations + 1),
                "permutations": permutations,
                "seed": seed,
                "null_method": "within_region_within_mapped_component_B_shared_site_label_permutation",
                "vaf_condition": (
                    "conditional_on_frozen_selected_top_set_no_VAF_rescoring"
                    if group.endpoint.startswith("vaf_")
                    else "not_applicable"
                ),
                "k_definition": (
                    "all_fixed_evaluable_regions"
                    if group.k_bin == "all"
                    else "max_shared_k_across_components_of_fixed_observed_mapping"
                ),
                **bootstrap,
            }
            if group.k_bin == "all":
                all_rows.append(row)
            else:
                complexity_rows.append(row)
    return all_rows, complexity_rows, loco_rows


def stringify_rows(rows: Sequence[dict], float_fields: Iterable[str]) -> list[dict]:
    fields = set(float_fields)
    return [
        {
            key: format_float(value) if key in fields and value != "" else value
            for key, value in row.items()
        }
        for row in rows
    ]


def main() -> None:
    args = parse_args()
    if args.permutations <= 0:
        raise SystemExit("--permutations must be positive")
    if args.bootstrap_replicates <= 0:
        raise SystemExit("--bootstrap-replicates must be positive")
    bootstrap_seed = args.bootstrap_seed if args.bootstrap_seed is not None else args.seed + 1
    checks: list[dict] = []

    records = read_jsonl_gz(args.input)
    add_check(checks, "fixed_region_population", "5720", str(len(records)), len(records) == 5720)
    match_ids = [str(record.get("match_id")) for record in records]
    add_check(
        checks,
        "unique_match_ids",
        str(len(records)),
        str(len(set(match_ids))),
        len(set(match_ids)) == len(records),
    )
    chromosomes = sorted(
        {str(record["chrom"]) for record in records},
        key=lambda value: int(value.removeprefix("chr")),
    )
    add_check(checks, "autosome_blocks", "22", str(len(chromosomes)), len(chromosomes) == 22)

    # Self-tests specifically target direction handling under vertex relabeling.
    add_check(
        checks,
        "relation_transform_identity",
        "FRP",
        transform_signature("FRP", 3, (0, 1, 2)),
        transform_signature("FRP", 3, (0, 1, 2)) == "FRP",
    )
    add_check(
        checks,
        "relation_transform_two_site_direction_flip",
        "R",
        transform_signature("F", 2, (1, 0)),
        transform_signature("F", 2, (1, 0)) == "R",
    )

    regions = build_region_plans(records, checks)
    groups = initialize_groups(regions, args.permutations, checks)
    run_null(
        regions,
        groups,
        args.permutations,
        args.seed,
        args.progress_every,
    )

    # Metric nesting is checked for every iteration and every non-empty group.
    nesting_pass = True
    for group in groups.values():
        if group.denominator == 0:
            continue
        exact = group.null_counts["exact"]
        nested = group.null_counts["robust_nested"]
        compatible = group.null_counts["any_compatible"]
        unique = group.null_counts["unique_exact"]
        nesting_pass &= bool(
            np.all(unique <= exact)
            and np.all(exact <= nested)
            and np.all(nested <= compatible)
            and np.all(compatible <= group.denominator)
        )
    add_check(checks, "null_metric_nesting_all_iterations", "PASS", "PASS" if nesting_pass else "FAIL", nesting_pass)
    add_check(
        checks,
        "permutation_count",
        str(args.permutations),
        str(args.permutations),
        args.permutations == len(next(iter(groups.values())).null_counts["exact"]),
    )
    add_check(
        checks,
        "frozen_VAF_top_set_condition",
        "no VAF rescoring",
        "conditional_on_frozen_selected_top_set_no_VAF_rescoring",
        True,
    )
    add_check(
        checks,
        "component_invariants",
        "k, Q size, candidate multiplicity, HP mapping preserved",
        "preserved by bijective label relabeling; source compact rows are read-only",
        True,
    )

    all_rows, complexity_rows, loco_rows = summarize_groups(
        groups,
        args.permutations,
        args.seed,
        args.bootstrap_replicates,
        bootstrap_seed,
        chromosomes,
    )
    add_check(
        checks,
        "main_metric_rows",
        str(len(ENDPOINTS) * len(COMPARISONS) * len(METRICS)),
        str(len(all_rows)),
        len(all_rows) == len(ENDPOINTS) * len(COMPARISONS) * len(METRICS),
    )
    add_check(
        checks,
        "complexity_strata_present",
        ",".join(K_BINS),
        ",".join(sorted({row["k_bin"] for row in complexity_rows})),
        set(K_BINS) == {row["k_bin"] for row in complexity_rows},
    )
    add_check(
        checks,
        "chromosome_block_bootstrap",
        str(args.bootstrap_replicates),
        str(min(row["block_bootstrap_replicates_valid"] for row in all_rows)),
        all(row["block_bootstrap_replicates_valid"] > 0 for row in all_rows),
    )
    add_check(
        checks,
        "leave_one_chromosome_out",
        f"{len(all_rows) * len(chromosomes)} main rows",
        f"{sum(row['k_bin'] == 'all' for row in loco_rows)} main rows",
        sum(row["k_bin"] == "all" for row in loco_rows)
        == len(all_rows) * len(chromosomes),
    )

    failed = [row for row in checks if row["status"] != "PASS"]
    if failed:
        raise RuntimeError(f"null validation checks failed: {canonical_json(failed)}")

    float_fields = {
        "observed_rate",
        "null_mean_n",
        "null_mean_rate",
        "null_sd_rate",
        "null_q025_rate",
        "null_q975_rate",
        "observed_minus_null_mean",
        "empirical_p_ge",
        "block_bootstrap_excess_ci95_low",
        "block_bootstrap_excess_ci95_high",
        "loco_min_excess",
        "loco_max_excess",
    }
    loco_float_fields = {"observed_rate", "null_mean_rate", "observed_minus_null_mean"}
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    null_tsv = output_dir / "hcc1395_topology_null.tsv"
    null_json = output_dir / "hcc1395_topology_null.json"
    complexity_tsv = output_dir / "hcc1395_topology_complexity_strata.tsv"
    loco_tsv = output_dir / "hcc1395_topology_null_loco.tsv"
    checks_tsv = output_dir / "hcc1395_topology_null_checks.tsv"

    result_fields = [
        "endpoint",
        "comparison",
        "population",
        "k_bin",
        "metric",
        "denominator",
        "observed_n",
        "observed_rate",
        "null_mean_n",
        "null_mean_rate",
        "null_sd_rate",
        "null_q025_rate",
        "null_q975_rate",
        "observed_minus_null_mean",
        "empirical_p_ge",
        "block_bootstrap_replicates_valid",
        "block_bootstrap_excess_ci95_low",
        "block_bootstrap_excess_ci95_high",
        "loco_min_excess",
        "loco_max_excess",
        "loco_positive_chromosomes",
        "loco_chromosomes",
        "permutations",
        "seed",
        "null_method",
        "vaf_condition",
        "k_definition",
    ]
    loco_fields = [
        "endpoint",
        "comparison",
        "k_bin",
        "metric",
        "held_out_chromosome",
        "denominator",
        "observed_rate",
        "null_mean_rate",
        "observed_minus_null_mean",
    ]
    write_tsv(null_tsv, stringify_rows(all_rows, float_fields), result_fields)
    write_tsv(complexity_tsv, stringify_rows(complexity_rows, float_fields), result_fields)
    write_tsv(loco_tsv, stringify_rows(loco_rows, loco_float_fields), loco_fields)
    write_tsv(checks_tsv, checks, ("check", "expected", "observed", "status"))

    document = {
        "schema_version": "1.0",
        "input": {
            "path": str(args.input.resolve()),
            "sha256": sha256(args.input),
            "regions": len(records),
            "chromosomes": chromosomes,
        },
        "method": {
            "null": "within-region, within-mapped-component B-side shared-site label permutation",
            "relation_vector": "F/R/P lexicographic site-pair vector with direction-aware vertex relabeling",
            "permutations": args.permutations,
            "seed": args.seed,
            "comparisons": list(COMPARISONS),
            "swap_tolerant_rule": "Jaccard, min directional coverage, intersection, category rank, identity tie-break",
            "endpoints": list(ENDPOINTS),
            "metrics": list(METRICS),
            "fixed_population": "observed-evaluable regions separately for each endpoint/comparison",
            "preserved": [
                "within-component shared k",
                "projected Q-set size",
                "pre-projection candidate multiplicity",
                "HP mapping",
            ],
            "vaf_condition": "conditional on frozen selected top set; VAF is not rescored",
            "complexity_k": "maximum shared k among mapped HP components of the fixed observed mapping",
            "block_bootstrap": "sample 22 chromosome blocks with replacement",
            "bootstrap_replicates": args.bootstrap_replicates,
            "bootstrap_seed": bootstrap_seed,
            "loco": "leave one autosome out from observed-minus-null-mean excess",
        },
        "claim_ceiling": (
            "Conditional technical reproducibility null only; not biological topology truth, "
            "not an independent validation, and not a calibrated VAF posterior."
        ),
        "results": all_rows,
        "complexity_strata": complexity_rows,
        "leave_one_chromosome_out": loco_rows,
        "checks": checks,
        "outputs": {
            # Keep the scientific receipt byte-identical across output roots.
            # The CLI receipt below still reports the concrete execution paths.
            "null_tsv": null_tsv.name,
            "complexity_tsv": complexity_tsv.name,
            "loco_tsv": loco_tsv.name,
            "checks_tsv": checks_tsv.name,
        },
    }
    write_json(null_json, document)
    print(
        canonical_json(
            {
                "status": "PASS",
                "regions": len(records),
                "permutations": args.permutations,
                "bootstrap_replicates": args.bootstrap_replicates,
                "checks": f"{len(checks)}/{len(checks)} PASS",
                "outputs": {
                    "null_tsv": str(null_tsv),
                    "null_json": str(null_json),
                    "complexity_tsv": str(complexity_tsv),
                    "loco_tsv": str(loco_tsv),
                    "checks_tsv": str(checks_tsv),
                },
            }
        )
    )


if __name__ == "__main__":
    main()
