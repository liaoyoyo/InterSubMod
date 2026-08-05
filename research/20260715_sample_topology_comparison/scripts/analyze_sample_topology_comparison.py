#!/usr/bin/env python3
"""Build a current-v5 seven-dataset topology comparison artifact.

The analysis keeps three mutually exclusive W_primary partitions separate:
structural determinacy, descriptive read-AF selection, and mutation-state
morphology.  HCC1395 and HCC1395_DORADO are treated as a technical/cross-
platform pair from one biological sample, never as biological replicates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SIDECAR_DIR = (
    REPO_ROOT
    / "research"
    / "20260715_layered_workstation_genome_topology_multiselect"
    / "data"
    / "current_v5_read_af_topology"
)
DEFAULT_UNKNOWN_K = (
    REPO_ROOT
    / "research"
    / "20260714_hcc1395_unknown_k_clone_state_consistency"
    / "results"
    / "unknown_k_consistency.json"
)
DEFAULT_EXACT_SIGNATURE = (
    REPO_ROOT
    / "research"
    / "20260715_sample_topology_comparison"
    / "artifacts"
    / "hcc1395_exact_signature_validation.json"
)
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parents[1] / "artifacts"

SAMPLE_ORDER = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
BIOLOGICAL_ID = {sample: sample for sample in SAMPLE_ORDER}
BIOLOGICAL_ID["HCC1395_DORADO"] = "HCC1395"
AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))

DIMENSIONS: dict[str, dict[str, Any]] = {
    "structural": {
        "summary_key": "structural_classes",
        "region_key": "structural_class",
        "label": "結構可辨識度",
        "categories": (
            "exact_and_topology_unique",
            "topology_unique_exact_multiple",
            "topology_multiple_exact_multiple",
            "incomplete",
        ),
    },
    "read_af_selection": {
        "summary_key": "selection_classes",
        "region_key": "selection_class",
        "label": "read-AF 第一順位",
        "categories": (
            "structural_exact_unique",
            "read_af_unique_first",
            "read_af_tied_same_topology",
            "read_af_tied_different_topology",
            "read_af_unavailable",
            "incomplete",
        ),
    },
    "morphology": {
        "summary_key": "morphology_classes",
        "region_key": "morphology_class",
        "label": "mutation-state morphology",
        "categories": (
            "single_no_within_hp_relation",
            "direct_chain",
            "sister_branch",
            "direct_and_sister",
            "unresolved",
        ),
    },
}
EVALUABLE_CATEGORIES = {
    "structural": (
        "exact_and_topology_unique",
        "topology_unique_exact_multiple",
        "topology_multiple_exact_multiple",
    ),
    "read_af_selection": (
        "structural_exact_unique",
        "read_af_unique_first",
        "read_af_tied_same_topology",
        "read_af_tied_different_topology",
    ),
    "morphology": (
        "single_no_within_hp_relation",
        "direct_chain",
        "sister_branch",
        "direct_and_sister",
    ),
}

OPERATIONAL_THRESHOLDS = {
    "tvd": {"close_max": 0.05, "moderate_max": 0.10},
    "region_jaccard": {"high_min": 0.80, "moderate_min": 0.60},
    "raw_agreement": {"high_min": 0.80, "moderate_min": 0.65},
    "kappa": {"substantial_min": 0.60, "moderate_min": 0.40},
    "hcc_relative_rank_max": 3,
    "note": "Pre-registered operational display bands; not universal field cutoffs.",
}


class ContractError(RuntimeError):
    """Raised when an input or output contract is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"cannot read JSON {path}: {exc}") from exc
    require(isinstance(payload, dict), f"JSON root is not an object: {path}")
    return payload


def normalized(counts: dict[str, int], categories: Iterable[str]) -> dict[str, float]:
    categories = tuple(categories)
    total = sum(int(counts.get(category, 0)) for category in categories)
    require(total > 0, "cannot normalize an empty composition")
    return {category: int(counts.get(category, 0)) / total for category in categories}


def total_variation(left: dict[str, float], right: dict[str, float]) -> float:
    keys = set(left) | set(right)
    return 0.5 * sum(abs(left.get(key, 0.0) - right.get(key, 0.0)) for key in keys)


def jensen_shannon_distance(left: dict[str, float], right: dict[str, float]) -> float:
    keys = set(left) | set(right)
    midpoint = {key: 0.5 * (left.get(key, 0.0) + right.get(key, 0.0)) for key in keys}

    def divergence(source: dict[str, float]) -> float:
        return sum(
            value * math.log2(value / midpoint[key])
            for key, value in source.items()
            if value > 0.0 and midpoint[key] > 0.0
        )

    return math.sqrt(0.5 * divergence(left) + 0.5 * divergence(right))


def hellinger_distance(left: dict[str, float], right: dict[str, float]) -> float:
    keys = set(left) | set(right)
    return math.sqrt(
        0.5
        * sum(
            (math.sqrt(left.get(key, 0.0)) - math.sqrt(right.get(key, 0.0))) ** 2
            for key in keys
        )
    )


def quantile(values: list[float], probability: float) -> float:
    require(values, "quantile requires at least one value")
    require(0.0 <= probability <= 1.0, "invalid quantile probability")
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def confusion_statistics(
    confusion: Counter[tuple[str, str]], categories: Iterable[str]
) -> dict[str, Any]:
    categories = tuple(categories)
    n = sum(confusion.values())
    require(n > 0, "confusion matrix is empty")
    row = {category: sum(confusion[(category, other)] for other in categories) for category in categories}
    column = {category: sum(confusion[(other, category)] for other in categories) for category in categories}
    diagonal = {category: confusion[(category, category)] for category in categories}
    agreement_count = sum(diagonal.values())
    agreement = agreement_count / n
    expected = sum(row[category] * column[category] for category in categories) / (n * n)
    kappa = (agreement - expected) / (1.0 - expected) if expected < 1.0 else 1.0
    left_profile = {category: row[category] / n for category in categories}
    right_profile = {category: column[category] / n for category in categories}
    category_jaccard = {
        category: (
            diagonal[category] / (row[category] + column[category] - diagonal[category])
            if row[category] + column[category] - diagonal[category] > 0
            else None
        )
        for category in categories
    }
    return {
        "n": n,
        "agreement_count": agreement_count,
        "raw_agreement": agreement,
        "expected_agreement": expected,
        "cohen_kappa": kappa,
        "matched_marginal_tvd": total_variation(left_profile, right_profile),
        "left_counts": row,
        "right_counts": column,
        "diagonal_counts": diagonal,
        "per_category_jaccard": category_jaccard,
        "matrix": {
            left: {right: confusion[(left, right)] for right in categories}
            for left in categories
        },
    }


def classify_tvd(value: float) -> str:
    if value <= OPERATIONAL_THRESHOLDS["tvd"]["close_max"]:
        return "close"
    if value <= OPERATIONAL_THRESHOLDS["tvd"]["moderate_max"]:
        return "moderate_difference"
    return "material_difference"


def classify_overlap(value: float) -> str:
    if value >= OPERATIONAL_THRESHOLDS["region_jaccard"]["high_min"]:
        return "high"
    if value >= OPERATIONAL_THRESHOLDS["region_jaccard"]["moderate_min"]:
        return "moderate"
    return "limited"


def classify_agreement(value: float) -> str:
    if value >= OPERATIONAL_THRESHOLDS["raw_agreement"]["high_min"]:
        return "high"
    if value >= OPERATIONAL_THRESHOLDS["raw_agreement"]["moderate_min"]:
        return "moderate"
    return "limited"


def classify_kappa(value: float) -> str:
    if value >= OPERATIONAL_THRESHOLDS["kappa"]["substantial_min"]:
        return "substantial"
    if value >= OPERATIONAL_THRESHOLDS["kappa"]["moderate_min"]:
        return "moderate"
    return "limited"


def bootstrap_confusion_by_chromosome(
    by_chromosome: dict[str, Counter[tuple[str, str]]],
    categories: Iterable[str],
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    categories = tuple(categories)
    require(iterations >= 100, "bootstrap iterations must be >=100")
    rng = random.Random(seed)
    agreement_values: list[float] = []
    kappa_values: list[float] = []
    tvd_values: list[float] = []
    for _ in range(iterations):
        sampled = rng.choices(AUTOSOMES, k=len(AUTOSOMES))
        aggregate: Counter[tuple[str, str]] = Counter()
        for chrom in sampled:
            aggregate.update(by_chromosome.get(chrom, Counter()))
        metrics = confusion_statistics(aggregate, categories)
        agreement_values.append(metrics["raw_agreement"])
        kappa_values.append(metrics["cohen_kappa"])
        tvd_values.append(metrics["matched_marginal_tvd"])
    return {
        "unit": "autosome",
        "n_blocks": len(AUTOSOMES),
        "iterations": iterations,
        "seed": seed,
        "raw_agreement_ci95": [quantile(agreement_values, 0.025), quantile(agreement_values, 0.975)],
        "cohen_kappa_ci95": [quantile(kappa_values, 0.025), quantile(kappa_values, 0.975)],
        "matched_marginal_tvd_ci95": [quantile(tvd_values, 0.025), quantile(tvd_values, 0.975)],
    }


def mean_profiles(profiles: list[dict[str, float]], categories: Iterable[str]) -> dict[str, float]:
    categories = tuple(categories)
    require(profiles, "cannot average zero profiles")
    return {
        category: sum(profile[category] for profile in profiles) / len(profiles)
        for category in categories
    }


def load_sidecars(sidecar_dir: Path) -> tuple[dict[str, Any], dict[str, dict[str, Any]], list[str]]:
    index_path = sidecar_dir / "current_v5_read_af_topology.index.json"
    index = load_json(index_path)
    require(index.get("schema_name") == "intersubmod.current_v5_read_af_topology_index", "sidecar index schema_name drift")
    require(index.get("schema_version") == "1.1.0", "sidecar index schema_version drift")
    require(index.get("scope") == "GRCh38 chr1-22 current canonical v5", "sidecar scope drift")
    require(index.get("dataset_count") == 7 and index.get("all_checks_pass") is True, "sidecar index is not 7/7 PASS")
    index_records = {record["sample"]: record for record in index.get("samples", [])}
    require(set(index_records) == set(SAMPLE_ORDER), "sidecar sample set drift")

    checks: list[str] = []
    datasets: dict[str, dict[str, Any]] = {}
    for sample in SAMPLE_ORDER:
        record = index_records[sample]
        source_path = Path(record["output"]).resolve()
        require(source_path.parent == sidecar_dir.resolve(), f"{sample}: sidecar escapes expected directory")
        require(source_path.is_file(), f"{sample}: missing sidecar")
        require(sha256_file(source_path) == record["output_sha256"], f"{sample}: sidecar hash drift")
        payload = load_json(source_path)
        require(payload.get("schema_name") == "intersubmod.current_v5_read_af_topology_sample", f"{sample}: schema_name drift")
        require(payload.get("schema_version") == "1.1.0", f"{sample}: schema_version drift")
        require(payload.get("sample") == sample, f"{sample}: payload sample mismatch")
        require(payload.get("summary", {}).get("all_checks_pass") is True, f"{sample}: summary checks failed")

        regions = payload.get("regions") or []
        region_map = {region["region"]: region for region in regions}
        require(len(region_map) == len(regions), f"{sample}: duplicate region keys")
        summary = payload["summary"]
        require(len(regions) == int(summary["W_tree"]), f"{sample}: region count/W_tree mismatch")
        primary = {
            key: region
            for key, region in region_map.items()
            if region.get("structural_class") != "no_primary_lineage"
        }
        require(len(primary) == int(summary["W_primary"]), f"{sample}: primary region/W_primary mismatch")
        require(
            all(region.get("primary_families") for region in primary.values()),
            f"{sample}: primary region lacks primary families",
        )
        require(
            {region.get("chrom") for region in regions} == set(AUTOSOMES),
            f"{sample}: chromosome scope is not chr1-22",
        )

        dimension_payload: dict[str, Any] = {}
        for dimension, contract in DIMENSIONS.items():
            summary_counts = {category: int(summary[contract["summary_key"]][category]) for category in contract["categories"]}
            observed = Counter(region[contract["region_key"]] for region in primary.values())
            require(
                {category: observed[category] for category in contract["categories"]} == summary_counts,
                f"{sample}: {dimension} region/summary partition mismatch",
            )
            require(sum(summary_counts.values()) == len(primary), f"{sample}: {dimension} does not conserve W_primary")
            dimension_payload[dimension] = {
                "counts": summary_counts,
                "proportions": normalized(summary_counts, contract["categories"]),
            }
            checks.append(f"{sample}:{dimension}:partition_conserves_W_primary")

        require(int(record["W_tree"]) == int(summary["W_tree"]), f"{sample}: index W_tree mismatch")
        require(int(record["W_primary"]) == int(summary["W_primary"]), f"{sample}: index W_primary mismatch")
        datasets[sample] = {
            "dataset": sample,
            "biological_id": BIOLOGICAL_ID[sample],
            "relationship": "technical_cross_platform_pair" if sample in {"HCC1395", "HCC1395_DORADO"} else "single_dataset",
            "W_tree": int(summary["W_tree"]),
            "W_primary": int(summary["W_primary"]),
            "no_primary": int(summary["no_primary"]),
            "dimensions": dimension_payload,
            "regions": region_map,
            "primary_regions": primary,
            "source_path": str(source_path),
            "source_sha256": record["output_sha256"],
        }
    return index, datasets, checks


def build_aggregates(datasets: dict[str, dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {
        "dataset_micro": {"W_primary": sum(datasets[sample]["W_primary"] for sample in SAMPLE_ORDER), "dimensions": {}},
        "dataset_macro": {"weighting": "equal dataset weight", "dimensions": {}},
        "biological_sample_macro": {
            "n_biological_samples": 6,
            "weighting": "average technical datasets within biological ID, then equal biological-sample weight",
            "members": {},
            "dimensions": {},
        },
    }
    biological_groups: dict[str, list[str]] = defaultdict(list)
    for sample in SAMPLE_ORDER:
        biological_groups[BIOLOGICAL_ID[sample]].append(sample)

    for dimension, contract in DIMENSIONS.items():
        categories = contract["categories"]
        micro_counts = {
            category: sum(datasets[sample]["dimensions"][dimension]["counts"][category] for sample in SAMPLE_ORDER)
            for category in categories
        }
        result["dataset_micro"]["dimensions"][dimension] = {
            "counts": micro_counts,
            "proportions": normalized(micro_counts, categories),
        }
        result["dataset_macro"]["dimensions"][dimension] = {
            "proportions": mean_profiles(
                [datasets[sample]["dimensions"][dimension]["proportions"] for sample in SAMPLE_ORDER],
                categories,
            )
        }
        biological_profiles = []
        for biological_id, members in sorted(biological_groups.items()):
            profile = mean_profiles(
                [datasets[sample]["dimensions"][dimension]["proportions"] for sample in members],
                categories,
            )
            biological_profiles.append(profile)
            member_record = result["biological_sample_macro"]["members"].setdefault(
                biological_id, {"datasets": members, "dimensions": {}}
            )
            member_record["dimensions"][dimension] = {"proportions": profile}
        result["biological_sample_macro"]["dimensions"][dimension] = {
            "proportions": mean_profiles(biological_profiles, categories)
        }
    return result


def build_pairwise(datasets: dict[str, dict[str, Any]]) -> dict[str, Any]:
    records = []
    for left, right in combinations(SAMPLE_ORDER, 2):
        dimensions = {}
        for dimension in DIMENSIONS:
            left_profile = datasets[left]["dimensions"][dimension]["proportions"]
            right_profile = datasets[right]["dimensions"][dimension]["proportions"]
            evaluable = EVALUABLE_CATEGORIES[dimension]
            left_evaluable = normalized(
                datasets[left]["dimensions"][dimension]["counts"], evaluable
            )
            right_evaluable = normalized(
                datasets[right]["dimensions"][dimension]["counts"], evaluable
            )
            dimensions[dimension] = {
                "tvd": total_variation(left_profile, right_profile),
                "jensen_shannon_distance": jensen_shannon_distance(left_profile, right_profile),
                "hellinger_distance": hellinger_distance(left_profile, right_profile),
                "conditional_evaluable_tvd": total_variation(left_evaluable, right_evaluable),
            }
        mean_tvd = sum(value["tvd"] for value in dimensions.values()) / len(dimensions)
        conditional_mean_tvd = sum(value["conditional_evaluable_tvd"] for value in dimensions.values()) / len(dimensions)
        records.append(
            {
                "left": left,
                "right": right,
                "same_biological_sample": BIOLOGICAL_ID[left] == BIOLOGICAL_ID[right],
                "dimensions": dimensions,
                "profile_mean_tvd": mean_tvd,
                "profile_mean_tvd_band": classify_tvd(mean_tvd),
                "conditional_evaluable_mean_tvd": conditional_mean_tvd,
                "conditional_evaluable_mean_tvd_band": classify_tvd(conditional_mean_tvd),
            }
        )
    records.sort(key=lambda record: (record["profile_mean_tvd"], record["left"], record["right"]))
    for rank, record in enumerate(records, start=1):
        record["rank_by_profile_mean_tvd"] = rank
    for dimension in DIMENSIONS:
        ordered = sorted(
            records,
            key=lambda record: (
                record["dimensions"][dimension]["tvd"],
                record["left"],
                record["right"],
            ),
        )
        for rank, record in enumerate(ordered, start=1):
            record["dimensions"][dimension]["tvd_rank_among_21_pairs"] = rank
        conditional_ordered = sorted(
            records,
            key=lambda record: (
                record["dimensions"][dimension]["conditional_evaluable_tvd"],
                record["left"],
                record["right"],
            ),
        )
        for rank, record in enumerate(conditional_ordered, start=1):
            record["dimensions"][dimension]["conditional_evaluable_tvd_rank_among_21_pairs"] = rank
    conditional_records = sorted(
        records,
        key=lambda record: (
            record["conditional_evaluable_mean_tvd"], record["left"], record["right"]
        ),
    )
    for rank, record in enumerate(conditional_records, start=1):
        record["rank_by_conditional_evaluable_mean_tvd"] = rank
    require(len(records) == 21, "pairwise comparison is not 21 unique pairs")

    matrices: dict[str, Any] = {}
    for matrix_dimension in (*DIMENSIONS.keys(), "profile_mean"):
        matrix = {sample: {other: 0.0 for other in SAMPLE_ORDER} for sample in SAMPLE_ORDER}
        for record in records:
            value = (
                record["profile_mean_tvd"]
                if matrix_dimension == "profile_mean"
                else record["dimensions"][matrix_dimension]["tvd"]
            )
            matrix[record["left"]][record["right"]] = value
            matrix[record["right"]][record["left"]] = value
        matrices[matrix_dimension] = matrix

    nearest = {}
    for sample in SAMPLE_ORDER:
        candidates = [
            record for record in records if sample in {record["left"], record["right"]}
        ]
        best = min(candidates, key=lambda record: record["profile_mean_tvd"])
        nearest[sample] = {
            "dataset": best["right"] if best["left"] == sample else best["left"],
            "profile_mean_tvd": best["profile_mean_tvd"],
        }
    return {
        "metric_contract": {
            "primary": "total variation distance on mutually exclusive W_primary category proportions",
            "range": [0.0, 1.0],
            "profile_mean_tvd": "unweighted arithmetic mean of structural, read-AF-selection, and morphology TVD; navigation only, not biological distance",
            "conditional_evaluable_mean_tvd": "same unweighted mean after structural incomplete, read-AF unavailable/incomplete, and morphology unresolved are excluded and each dimension renormalized",
        },
        "records": records,
        "matrices": {"order": list(SAMPLE_ORDER), "tvd": matrices},
        "nearest_neighbor_by_profile_mean_tvd": nearest,
    }


def build_hcc_profile_bootstrap(
    datasets: dict[str, dict[str, Any]], iterations: int, seed: int
) -> dict[str, Any]:
    require(iterations >= 100, "profile bootstrap iterations must be >=100")
    chromosome_counts: dict[str, dict[str, dict[str, Counter[str]]]] = {}
    for sample in SAMPLE_ORDER:
        chromosome_counts[sample] = {
            dimension: {chrom: Counter() for chrom in AUTOSOMES}
            for dimension in DIMENSIONS
        }
        for region in datasets[sample]["primary_regions"].values():
            chrom = region["chrom"]
            for dimension, contract in DIMENSIONS.items():
                chromosome_counts[sample][dimension][chrom][region[contract["region_key"]]] += 1

    def profiles_for_blocks(blocks: list[str], conditional: bool) -> dict[str, dict[str, dict[str, float]]]:
        profiles: dict[str, dict[str, dict[str, float]]] = {}
        for sample in SAMPLE_ORDER:
            profiles[sample] = {}
            for dimension, contract in DIMENSIONS.items():
                counts: Counter[str] = Counter()
                for chrom in blocks:
                    counts.update(chromosome_counts[sample][dimension][chrom])
                categories = EVALUABLE_CATEGORIES[dimension] if conditional else contract["categories"]
                profiles[sample][dimension] = normalized(counts, categories)
        return profiles

    def pair_mean_tvd(
        profiles: dict[str, dict[str, dict[str, float]]], left: str, right: str
    ) -> float:
        return sum(
            total_variation(profiles[left][dimension], profiles[right][dimension])
            for dimension in DIMENSIONS
        ) / len(DIMENSIONS)

    rng = random.Random(seed)
    full_values: list[float] = []
    conditional_values: list[float] = []
    full_ranks: list[int] = []
    conditional_ranks: list[int] = []
    for _ in range(iterations):
        blocks = rng.choices(AUTOSOMES, k=len(AUTOSOMES))
        full = profiles_for_blocks(blocks, conditional=False)
        conditional = profiles_for_blocks(blocks, conditional=True)
        full_pair_values = sorted(
            pair_mean_tvd(full, left, right)
            for left, right in combinations(SAMPLE_ORDER, 2)
        )
        conditional_pair_values = sorted(
            pair_mean_tvd(conditional, left, right)
            for left, right in combinations(SAMPLE_ORDER, 2)
        )
        full_hcc = pair_mean_tvd(full, "HCC1395", "HCC1395_DORADO")
        conditional_hcc = pair_mean_tvd(conditional, "HCC1395", "HCC1395_DORADO")
        full_values.append(full_hcc)
        conditional_values.append(conditional_hcc)
        full_ranks.append(1 + sum(value < full_hcc - 1e-15 for value in full_pair_values))
        conditional_ranks.append(1 + sum(value < conditional_hcc - 1e-15 for value in conditional_pair_values))

    loco_full = []
    loco_conditional = []
    for excluded in AUTOSOMES:
        blocks = [chrom for chrom in AUTOSOMES if chrom != excluded]
        loco_full.append(
            {
                "excluded": excluded,
                "profile_mean_tvd": pair_mean_tvd(
                    profiles_for_blocks(blocks, conditional=False),
                    "HCC1395",
                    "HCC1395_DORADO",
                ),
            }
        )
        loco_conditional.append(
            {
                "excluded": excluded,
                "profile_mean_tvd": pair_mean_tvd(
                    profiles_for_blocks(blocks, conditional=True),
                    "HCC1395",
                    "HCC1395_DORADO",
                ),
            }
        )
    return {
        "unit": "autosome",
        "n_blocks": len(AUTOSOMES),
        "iterations": iterations,
        "seed": seed,
        "full_profile_mean_tvd_ci95": [quantile(full_values, 0.025), quantile(full_values, 0.975)],
        "conditional_evaluable_mean_tvd_ci95": [
            quantile(conditional_values, 0.025),
            quantile(conditional_values, 0.975),
        ],
        "full_rank_le_2_fraction": sum(rank <= 2 for rank in full_ranks) / iterations,
        "full_rank_le_3_fraction": sum(rank <= 3 for rank in full_ranks) / iterations,
        "conditional_rank_le_2_fraction": sum(rank <= 2 for rank in conditional_ranks) / iterations,
        "conditional_rank_le_3_fraction": sum(rank <= 3 for rank in conditional_ranks) / iterations,
        "leave_one_chromosome_out": {
            "full": loco_full,
            "conditional_evaluable": loco_conditional,
            "full_range": [
                min(record["profile_mean_tvd"] for record in loco_full),
                max(record["profile_mean_tvd"] for record in loco_full),
            ],
            "conditional_evaluable_range": [
                min(record["profile_mean_tvd"] for record in loco_conditional),
                max(record["profile_mean_tvd"] for record in loco_conditional),
            ],
        },
    }


def build_hcc_pair(
    datasets: dict[str, dict[str, Any]],
    pairwise: dict[str, Any],
    bootstrap_iterations: int,
    seed: int,
) -> dict[str, Any]:
    left_name = "HCC1395"
    right_name = "HCC1395_DORADO"
    left = datasets[left_name]
    right = datasets[right_name]
    left_all = set(left["regions"])
    right_all = set(right["regions"])
    left_primary = set(left["primary_regions"])
    right_primary = set(right["primary_regions"])
    all_intersection = left_all & right_all
    primary_intersection = left_primary & right_primary
    require(primary_intersection, "HCC primary exact-coordinate intersection is empty")

    pair_record = next(
        record
        for record in pairwise["records"]
        if {record["left"], record["right"]} == {left_name, right_name}
    )
    dimension_results = {}
    matched_examples: list[dict[str, Any]] = []
    for dimension, contract in DIMENSIONS.items():
        categories = contract["categories"]
        by_chromosome: dict[str, Counter[tuple[str, str]]] = {chrom: Counter() for chrom in AUTOSOMES}
        confusion: Counter[tuple[str, str]] = Counter()
        for region_key in sorted(primary_intersection):
            left_region = left["primary_regions"][region_key]
            right_region = right["primary_regions"][region_key]
            labels = (left_region[contract["region_key"]], right_region[contract["region_key"]])
            confusion[labels] += 1
            by_chromosome[left_region["chrom"]][labels] += 1
            if labels[0] != labels[1] and len(matched_examples) < 12:
                matched_examples.append(
                    {
                        "region": region_key,
                        "chrom": left_region["chrom"],
                        "dimension": dimension,
                        "left_label": labels[0],
                        "right_label": labels[1],
                    }
                )
        metrics = confusion_statistics(confusion, categories)
        bootstrap = bootstrap_confusion_by_chromosome(
            by_chromosome,
            categories,
            bootstrap_iterations,
            seed + list(DIMENSIONS).index(dimension),
        )
        chromosome_breakdown = {}
        for chrom in AUTOSOMES:
            if sum(by_chromosome[chrom].values()) == 0:
                chromosome_breakdown[chrom] = {"n": 0, "raw_agreement": None, "cohen_kappa": None}
            else:
                chrom_metrics = confusion_statistics(by_chromosome[chrom], categories)
                chromosome_breakdown[chrom] = {
                    "n": chrom_metrics["n"],
                    "raw_agreement": chrom_metrics["raw_agreement"],
                    "cohen_kappa": chrom_metrics["cohen_kappa"],
                }
        metrics["raw_agreement_band"] = classify_agreement(metrics["raw_agreement"])
        metrics["cohen_kappa_band"] = classify_kappa(metrics["cohen_kappa"])
        metrics["bootstrap"] = bootstrap
        metrics["by_chromosome"] = chromosome_breakdown
        dimension_results[dimension] = metrics

    primary_jaccard = len(primary_intersection) / len(left_primary | right_primary)
    profile_rank = int(pair_record["rank_by_profile_mean_tvd"])
    rank_band = "top3_relative_close" if profile_rank <= OPERATIONAL_THRESHOLDS["hcc_relative_rank_max"] else "not_top3"
    profile_bootstrap = build_hcc_profile_bootstrap(datasets, bootstrap_iterations, seed + 100)
    return {
        "left": left_name,
        "right": right_name,
        "biological_id": "HCC1395",
        "relationship": "same biological sample; technical/cross-platform pair; not biological replication",
        "marginal_profiles": {
            "dimensions": pair_record["dimensions"],
            "profile_mean_tvd": pair_record["profile_mean_tvd"],
            "profile_mean_tvd_band": pair_record["profile_mean_tvd_band"],
            "rank_among_21_pairs": profile_rank,
            "rank_band": rank_band,
            "conditional_evaluable_mean_tvd": pair_record["conditional_evaluable_mean_tvd"],
            "conditional_evaluable_mean_tvd_band": pair_record["conditional_evaluable_mean_tvd_band"],
            "conditional_evaluable_rank_among_21_pairs": pair_record["rank_by_conditional_evaluable_mean_tvd"],
            "chromosome_block_bootstrap": profile_bootstrap,
        },
        "exact_coordinate_overlap": {
            "all_regions": {
                "left": len(left_all),
                "right": len(right_all),
                "intersection": len(all_intersection),
                "union": len(left_all | right_all),
                "jaccard": len(all_intersection) / len(left_all | right_all),
                "left_coverage": len(all_intersection) / len(left_all),
                "right_coverage": len(all_intersection) / len(right_all),
            },
            "primary_both": {
                "left": len(left_primary),
                "right": len(right_primary),
                "intersection": len(primary_intersection),
                "union": len(left_primary | right_primary),
                "jaccard": primary_jaccard,
                "jaccard_band": classify_overlap(primary_jaccard),
                "left_coverage": len(primary_intersection) / len(left_primary),
                "right_coverage": len(primary_intersection) / len(right_primary),
            },
        },
        "matched_primary_dimensions": dimension_results,
        "discordant_examples": matched_examples,
        "claim_ceiling": "partial cross-technical reproducibility of current-v5 regional mutation-state labels; not accuracy truth, biological replication, confirmed clones, or interchangeable datasets",
    }


def prior_current_v5_evidence(path: Path) -> dict[str, Any]:
    payload = load_json(path)
    run = payload.get("input_run") or {}
    require(run.get("run_id") == "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5", "unknown-K evidence is not current canonical v5")
    comparison = payload.get("cross_technical_comparison") or {}
    shared = comparison.get("shared_exact_allele_regions") or {}
    return {
        "source": str(path.resolve()),
        "source_sha256": sha256_file(path),
        "run_id": run["run_id"],
        "exact_retained_site_overlap": comparison["exact_retained_site_overlap"],
        "exact_allele_region_overlap": comparison["exact_allele_region_overlap"],
        "caller_af_concordance": comparison["caller_af_concordance"],
        "pooled_read_af_concordance": comparison["pooled_read_af_concordance"],
        "exact_allele_region_agreements": {
            "n": shared["n"],
            "region_determinacy": shared["region_determinacy_agreement"] / shared["n"],
            "hp_flip_invariant_mutation_k": shared["flip_invariant_mutation_k_signature_agreement"] / shared["n"],
            "pair_complete_shape": shared["flip_invariant_shape_signature_agreement_among_pair_complete"] / shared["both_pair_complete"],
            "pair_complete_shape_n": shared["both_pair_complete"],
        },
        "relationship": "independent current-v5 evidence layer; exact allele-set region key, not the sidecar exact-coordinate key",
    }


def exact_signature_evidence(path: Path, sidecar_index_path: Path) -> dict[str, Any]:
    payload = load_json(path)
    require(
        payload.get("schema_name") == "intersubmod.hcc1395_exact_signature_validation"
        and payload.get("schema_version") == "1.0.0",
        "HCC exact-signature schema drift",
    )
    require(payload.get("all_checks_pass") is True, "HCC exact-signature checks failed")
    provenance = payload.get("provenance") or {}
    require(
        Path(provenance.get("sidecar_index", "")).resolve() == sidecar_index_path.resolve()
        and provenance.get("sidecar_index_sha256") == sha256_file(sidecar_index_path),
        "HCC exact-signature sidecar-index binding drifted",
    )
    analysis_script = Path(provenance.get("analysis_script", "")).resolve()
    require(analysis_script.is_file(), "HCC exact-signature analysis script is missing")
    require(
        provenance.get("analysis_script_sha256") == sha256_file(analysis_script),
        "HCC exact-signature analysis script hash drifted",
    )
    return {
        "source": str(path.resolve()),
        "source_sha256": sha256_file(path),
        "analysis_script": str(analysis_script),
        "analysis_script_sha256": provenance["analysis_script_sha256"],
        "claim_ceiling": payload["claim_ceiling"],
        "region_universe": payload["region_universe"],
        "retained_ssnv_universe": payload["retained_ssnv_universe"],
        "internal_ssnv_pairing": payload["internal_ssnv_pairing"],
        "shape_agreement": payload["shape_agreement"],
        "exact_labeled_edge_agreement": payload["exact_labeled_edge_agreement"],
        "sample_contracts": payload["sample_contracts"],
        "all_checks_pass": True,
    }


def compact_dataset(dataset: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in dataset.items() if key not in {"regions", "primary_regions"}}


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_profiles_tsv(path: Path, datasets: dict[str, dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("dataset", "biological_id", "W_primary", "dimension", "category", "count", "proportion"))
        for sample in SAMPLE_ORDER:
            for dimension, contract in DIMENSIONS.items():
                values = datasets[sample]["dimensions"][dimension]
                for category in contract["categories"]:
                    writer.writerow(
                        (
                            sample,
                            BIOLOGICAL_ID[sample],
                            datasets[sample]["W_primary"],
                            dimension,
                            category,
                            values["counts"][category],
                            f'{values["proportions"][category]:.12f}',
                        )
                    )


def write_pairwise_tsv(path: Path, pairwise: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        columns = [
            "left",
            "right",
            "same_biological_sample",
            "rank",
            "profile_mean_tvd",
            "conditional_evaluable_rank",
            "conditional_evaluable_mean_tvd",
        ]
        for dimension in DIMENSIONS:
            columns.extend((f"{dimension}_tvd", f"{dimension}_js_distance", f"{dimension}_hellinger"))
        writer.writerow(columns)
        for record in sorted(pairwise["records"], key=lambda row: row["rank_by_profile_mean_tvd"]):
            row: list[Any] = [
                record["left"],
                record["right"],
                str(record["same_biological_sample"]).lower(),
                record["rank_by_profile_mean_tvd"],
                f'{record["profile_mean_tvd"]:.12f}',
                record["rank_by_conditional_evaluable_mean_tvd"],
                f'{record["conditional_evaluable_mean_tvd"]:.12f}',
            ]
            for dimension in DIMENSIONS:
                metrics = record["dimensions"][dimension]
                row.extend(
                    (
                        f'{metrics["tvd"]:.12f}',
                        f'{metrics["jensen_shannon_distance"]:.12f}',
                        f'{metrics["hellinger_distance"]:.12f}',
                    )
                )
            writer.writerow(row)


def write_hcc_chromosome_tsv(path: Path, hcc_pair: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("dimension", "chrom", "n", "raw_agreement", "cohen_kappa"))
        for dimension in DIMENSIONS:
            for chrom in AUTOSOMES:
                values = hcc_pair["matched_primary_dimensions"][dimension]["by_chromosome"][chrom]
                writer.writerow(
                    (
                        dimension,
                        chrom,
                        values["n"],
                        "" if values["raw_agreement"] is None else f'{values["raw_agreement"]:.12f}',
                        "" if values["cohen_kappa"] is None else f'{values["cohen_kappa"]:.12f}',
                    )
                )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sidecar-dir", type=Path, default=DEFAULT_SIDECAR_DIR)
    parser.add_argument("--unknown-k", type=Path, default=DEFAULT_UNKNOWN_K)
    parser.add_argument("--exact-signature", type=Path, default=DEFAULT_EXACT_SIGNATURE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--bootstrap-iterations", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260715)
    args = parser.parse_args()

    sidecar_dir = args.sidecar_dir.resolve()
    index, datasets, checks = load_sidecars(sidecar_dir)
    aggregates = build_aggregates(datasets)
    pairwise = build_pairwise(datasets)
    hcc_pair = build_hcc_pair(datasets, pairwise, args.bootstrap_iterations, args.seed)
    input_index_path = sidecar_dir / "current_v5_read_af_topology.index.json"
    prior_evidence = prior_current_v5_evidence(args.unknown_k.resolve())
    exact_evidence = exact_signature_evidence(args.exact_signature.resolve(), input_index_path)
    artifact = {
        "schema_name": "intersubmod.sample_topology_comparison",
        "schema_version": "1.0.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task_type": "B_comprehensive_validation",
        "scope": "7 datasets / 6 biological samples / GRCh38 chr1-22 / current canonical v5",
        "claim_ceiling": "regional mutation-state topology profiles and cross-technical reproducibility; not biological clone distance, accuracy truth, or independent biological replication",
        "definitions": {
            "dataset_micro": "raw W_primary counts summed across seven datasets; HCC1395 contributes two technical datasets",
            "biological_sample_macro": "technical datasets averaged within biological ID, then six biological IDs equally weighted",
            "TVD": "half the L1 distance between mutually exclusive category proportions; 0=same composition, 1=no shared category mass",
            "profile_mean_tvd": "unweighted mean of three dimension TVDs for navigation only",
            "exact_coordinate_primary_match": "identical chrom:start-end region key with a primary lineage in both datasets",
        },
        "operational_thresholds": OPERATIONAL_THRESHOLDS,
        "provenance": {
            "sidecar_index": str(input_index_path),
            "sidecar_index_sha256": sha256_file(input_index_path),
            "sidecar_generated_at": index["generated_at"],
            "current_summary": index["provenance"]["current_summary"],
            "current_summary_sha256": index["provenance"]["current_summary_sha256"],
            "run_root": index["provenance"]["run_root"],
            "analysis_script": str(Path(__file__).resolve()),
            "analysis_script_sha256": sha256_file(Path(__file__).resolve()),
        },
        "dimension_contracts": DIMENSIONS,
        "datasets": [compact_dataset(datasets[sample]) for sample in SAMPLE_ORDER],
        "aggregates": aggregates,
        "pairwise_composition": pairwise,
        "hcc1395_technical_pair": hcc_pair,
        "independent_current_v5_hcc_evidence": prior_evidence,
        "hcc_exact_signature_evidence": exact_evidence,
        "checks": {
            "dataset_count_is_7": len(datasets) == 7,
            "biological_sample_count_is_6": len(set(BIOLOGICAL_ID.values())) == 6,
            "all_partitions_conserve_W_primary": len(checks) == 21,
            "pair_count_is_21": len(pairwise["records"]) == 21,
            "pairwise_matrices_are_symmetric": all(
                abs(pairwise["matrices"]["tvd"][dimension][left][right] - pairwise["matrices"]["tvd"][dimension][right][left]) < 1e-15
                for dimension in (*DIMENSIONS.keys(), "profile_mean")
                for left in SAMPLE_ORDER
                for right in SAMPLE_ORDER
            ),
            "hcc_primary_confusions_conserve_matched_n": all(
                values["n"] == hcc_pair["exact_coordinate_overlap"]["primary_both"]["intersection"]
                for values in hcc_pair["matched_primary_dimensions"].values()
            ),
            "hcc_is_technical_pair_not_biological_replication": hcc_pair["relationship"].startswith("same biological sample"),
            "hcc_exact_signature_checks_pass": exact_evidence["all_checks_pass"] is True,
        },
    }
    artifact["all_checks_pass"] = all(artifact["checks"].values())
    require(artifact["all_checks_pass"], "one or more comparison checks failed")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_json = output_dir / "sample_topology_comparison.json"
    profiles_tsv = output_dir / "sample_topology_profiles.tsv"
    pairwise_tsv = output_dir / "pairwise_profile_distances.tsv"
    hcc_chromosome_tsv = output_dir / "hcc1395_matched_by_chromosome.tsv"
    write_json(output_json, artifact)
    write_profiles_tsv(profiles_tsv, datasets)
    write_pairwise_tsv(pairwise_tsv, pairwise)
    write_hcc_chromosome_tsv(hcc_chromosome_tsv, hcc_pair)

    receipt = {
        "schema_name": "intersubmod.sample_topology_comparison_validation_receipt",
        "schema_version": "1.0.0",
        "status": "PASS",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "checks": artifact["checks"],
        "outputs": {
            str(path): {"sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in (output_json, profiles_tsv, pairwise_tsv, hcc_chromosome_tsv)
        },
    }
    receipt_path = output_dir / "validation_receipt.json"
    write_json(receipt_path, receipt)

    structural = hcc_pair["matched_primary_dimensions"]["structural"]
    selection = hcc_pair["matched_primary_dimensions"]["read_af_selection"]
    morphology = hcc_pair["matched_primary_dimensions"]["morphology"]
    print(
        "PASS datasets=7 biological_samples=6 W_primary="
        f'{aggregates["dataset_micro"]["W_primary"]} pairs={len(pairwise["records"])}'
    )
    print(
        "HCC exact-primary="
        f'{hcc_pair["exact_coordinate_overlap"]["primary_both"]["intersection"]} '
        f'Jaccard={hcc_pair["exact_coordinate_overlap"]["primary_both"]["jaccard"]:.6f} '
        f'profile_mean_TVD={hcc_pair["marginal_profiles"]["profile_mean_tvd"]:.6f} '
        f'rank={hcc_pair["marginal_profiles"]["rank_among_21_pairs"]}/21'
    )
    print(
        "HCC agreement/kappa "
        f'structural={structural["raw_agreement"]:.6f}/{structural["cohen_kappa"]:.6f} '
        f'read_AF={selection["raw_agreement"]:.6f}/{selection["cohen_kappa"]:.6f} '
        f'morphology={morphology["raw_agreement"]:.6f}/{morphology["cohen_kappa"]:.6f}'
    )
    print(f"OUTPUT {output_json}")
    print(f"RECEIPT {receipt_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
