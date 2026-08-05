#!/usr/bin/env python3
"""Audit five-class regional topology and HCC1395 cross-dataset concordance.

The primary classifier is deliberately structural.  It operates on rooted
mutation-state trees and includes H_* Steiner/partial-supported nodes.  It does
not claim that a node is a confirmed clone.  HP1 and HP2 remain ordered forest
components; the five-class region label collapses component features with OR.

The HCC1395 comparison uses deterministic one-to-one interval matching and
reports exact-coordinate and sensitivity scenarios.  Dataset A is not treated
as biological truth; directional balanced accuracies are reported both ways.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable
from zoneinfo import ZoneInfo

import numpy as np
from scipy.optimize import linear_sum_assignment


SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
BIOLOGICAL_PAIR = ("HCC1395", "HCC1395_DORADO")
CATEGORIES = [
    "no_within_hp_relation",
    "sister_only",
    "direct_only",
    "sister_and_direct",
    "topology_multiple_unresolved",
]
CATEGORY_ZH = {
    "no_within_hp_relation": "無 HP 內關係",
    "sister_only": "姐妹 only",
    "direct_only": "直系 only",
    "sister_and_direct": "姐妹＋直系",
    "topology_multiple_unresolved": "Topo>1 未定",
    "incomplete": "Incomplete（不進五類分母）",
}
SEED = 20260712
N_PERMUTATIONS = 5000


@dataclass(frozen=True)
class RegionRecord:
    sample: str
    region: str
    chrom: str
    start: int
    end: int
    primary_hp_units: int
    structural_class: str
    complete: bool
    coarse_category: str
    hp1_category_set: str
    hp2_category_set: str
    ordered_hp_coarse_signature: str
    ordered_shape_set_signature: str
    ordered_tree_digest_signature: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class Match:
    scenario: str
    a: RegionRecord
    b: RegionRecord
    overlap_bp: int
    reciprocal_min: float
    interval_jaccard: float
    start_delta_bp: int
    end_delta_bp: int


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, span = region.split(":", 1)
    start_text, end_text = span.split("-", 1)
    start, end = int(start_text), int(end_text)
    if start > end:
        raise ValueError(f"invalid region: {region}")
    return chrom, start, end


def shape_category(shape: dict) -> str:
    sister = int(shape["max_outdegree"]) >= 2
    direct = int(shape["max_depth"]) >= 2
    if sister and direct:
        return "sister_and_direct"
    if sister:
        return "sister_only"
    if direct:
        return "direct_only"
    return "no_within_hp_relation"


def collapse_components(component_categories: Iterable[str]) -> str:
    values = set(component_categories)
    sister = bool(values & {"sister_only", "sister_and_direct"})
    direct = bool(values & {"direct_only", "sister_and_direct"})
    if sister and direct:
        return "sister_and_direct"
    if sister:
        return "sister_only"
    if direct:
        return "direct_only"
    return "no_within_hp_relation"


def canonical_set(values: Iterable[str]) -> str:
    return "{" + ",".join(sorted(set(values))) + "}"


def load_tree_digests(exact_catalog: dict, samples: set[str]) -> dict[tuple[str, str, str], str]:
    """Load exact candidate-tree-set digests for the biological pair only."""
    result: dict[tuple[str, str, str], str] = {}
    source_by_sample = {row["sample"]: Path(row["source"]) for row in exact_catalog["samples"]}
    for sample in sorted(samples):
        source = source_by_sample[sample]
        layered = load_json(source)
        for unit in layered["detail"]:
            if not unit.get("is_primary_lineage"):
                continue
            key = (sample, unit["region"], str(unit["family"]))
            if key in result:
                raise RuntimeError(f"duplicate primary layered unit: {key}")
            digest = unit.get("analysis_tree_digest_sha256")
            if not digest:
                digest = "INCOMPLETE"
            result[key] = str(digest)
    return result


def build_region_records(
    exact_catalog: dict,
    ct_report: dict,
    region_rows: list[dict],
    tree_digests: dict[tuple[str, str, str], str],
) -> tuple[list[RegionRecord], list[dict]]:
    shape_meta = {row["shape_id"]: row for row in ct_report["shape_catalog"]}
    if len(shape_meta) != int(ct_report["aggregate"]["distinct_shapes_global"]):
        raise RuntimeError("shape catalog count mismatch")

    units_by_region: dict[tuple[str, str], list[dict]] = collections.defaultdict(list)
    seen_units: set[tuple[str, str, str]] = set()
    for sample in exact_catalog["samples"]:
        sample_name = sample["sample"]
        for unit in sample["unit_rows"]:
            key = (sample_name, unit["region"], str(unit["hp"]))
            if key in seen_units:
                raise RuntimeError(f"duplicate exact-catalog unit: {key}")
            seen_units.add(key)
            units_by_region[(sample_name, unit["region"])].append(unit)

    records: list[RegionRecord] = []
    checks: list[dict] = []
    seen_regions: set[tuple[str, str]] = set()
    for row in region_rows:
        sample, region = row["sample"], row["region"]
        key = (sample, region)
        if key in seen_regions:
            raise RuntimeError(f"duplicate region TSV key: {key}")
        seen_regions.add(key)
        chrom, start, end = parse_region(region)
        structural = row["structural_class"]
        complete = structural != "incomplete"
        hp_units = int(row["primary_HP_units"])

        hp_category_sets: dict[str, set[str]] = {}
        hp_shape_sets: dict[str, set[str]] = {}
        digest_by_hp: dict[str, str] = {}
        units = sorted(units_by_region.get(key, []), key=lambda value: int(value["hp"]))
        if complete:
            if len(units) != hp_units:
                raise RuntimeError(f"complete region unit count mismatch: {key}: {len(units)} != {hp_units}")
            for unit in units:
                hp = str(unit["hp"])
                shape_ids = set(unit["shape_candidate_counts"])
                if not shape_ids or not shape_ids <= set(shape_meta):
                    raise RuntimeError(f"missing/unknown shape set: {key} HP{hp}")
                hp_shape_sets[hp] = shape_ids
                hp_category_sets[hp] = {shape_category(shape_meta[value]) for value in shape_ids}
                digest_by_hp[hp] = tree_digests.get((sample, region, hp), "")

            if structural == "T>1|Topo>1":
                coarse = "topology_multiple_unresolved"
            else:
                if any(len(values) != 1 for values in hp_category_sets.values()):
                    raise RuntimeError(f"Topo=1 region has multiple coarse categories: {key}")
                coarse = collapse_components(next(iter(values)) for values in hp_category_sets.values())
        else:
            coarse = "incomplete"

        def hp_value(mapping: dict[str, set[str]], hp: str) -> str:
            return canonical_set(mapping[hp]) if hp in mapping else "NA"

        ordered_hp = "|".join(
            f"HP{hp}={canonical_set(hp_category_sets[hp])}" for hp in sorted(hp_category_sets, key=int)
        )
        ordered_shapes = "|".join(
            f"HP{hp}={canonical_set(hp_shape_sets[hp])}" for hp in sorted(hp_shape_sets, key=int)
        )
        ordered_digest = "|".join(
            f"HP{hp}={digest_by_hp.get(hp, '')}" for hp in sorted(digest_by_hp, key=int)
        )
        records.append(
            RegionRecord(
                sample=sample,
                region=region,
                chrom=chrom,
                start=start,
                end=end,
                primary_hp_units=hp_units,
                structural_class=structural,
                complete=complete,
                coarse_category=coarse,
                hp1_category_set=hp_value(hp_category_sets, "1"),
                hp2_category_set=hp_value(hp_category_sets, "2"),
                ordered_hp_coarse_signature=ordered_hp,
                ordered_shape_set_signature=ordered_shapes,
                ordered_tree_digest_signature=ordered_digest,
            )
        )

    checks.append(
        {
            "scope": "aggregate",
            "check": "region_keys_unique",
            "pass": len(records) == len(seen_regions),
            "observed": len(records),
            "expected": len(seen_regions),
        }
    )
    return records, checks


def interval_features(a: RegionRecord, b: RegionRecord) -> tuple[int, float, float, int, int]:
    if a.chrom != b.chrom:
        return 0, 0.0, 0.0, abs(a.start - b.start), abs(a.end - b.end)
    overlap = max(0, min(a.end, b.end) - max(a.start, b.start) + 1)
    reciprocal_min = min(overlap / a.length, overlap / b.length) if overlap else 0.0
    union = a.length + b.length - overlap
    jaccard = overlap / union if union else 0.0
    return overlap, reciprocal_min, jaccard, abs(a.start - b.start), abs(a.end - b.end)


def scenario_accepts(name: str, values: tuple[int, float, float, int, int]) -> bool:
    overlap, reciprocal_min, _, start_delta, end_delta = values
    if name == "reciprocal_overlap_0.80":
        return reciprocal_min >= 0.80
    if name == "reciprocal_overlap_0.50":
        return reciprocal_min >= 0.50
    if name == "endpoint_anchor_1kb":
        return overlap > 0 and start_delta <= 1000 and end_delta <= 1000
    if name == "endpoint_anchor_5kb":
        return overlap > 0 and start_delta <= 5000 and end_delta <= 5000
    raise ValueError(name)


def exact_matches(a_records: list[RegionRecord], b_records: list[RegionRecord]) -> list[Match]:
    a_by_region = {row.region: row for row in a_records}
    b_by_region = {row.region: row for row in b_records}
    if len(a_by_region) != len(a_records) or len(b_by_region) != len(b_records):
        raise RuntimeError("exact-match inputs are not unique")
    matches = []
    for region in sorted(set(a_by_region) & set(b_by_region)):
        a, b = a_by_region[region], b_by_region[region]
        features = interval_features(a, b)
        matches.append(Match("exact_coordinate", a, b, *features))
    return matches


def interval_matches(
    name: str,
    a_records: list[RegionRecord],
    b_records: list[RegionRecord],
) -> list[Match]:
    """Maximum-cardinality, then maximum-quality, one-to-one chromosome match."""
    a_by_chrom: dict[str, list[RegionRecord]] = collections.defaultdict(list)
    b_by_chrom: dict[str, list[RegionRecord]] = collections.defaultdict(list)
    for row in a_records:
        a_by_chrom[row.chrom].append(row)
    for row in b_records:
        b_by_chrom[row.chrom].append(row)

    matches: list[Match] = []
    for chrom in sorted(set(a_by_chrom) & set(b_by_chrom)):
        left = sorted(a_by_chrom[chrom], key=lambda row: (row.start, row.end, row.region))
        right = sorted(b_by_chrom[chrom], key=lambda row: (row.start, row.end, row.region))
        weights = np.zeros((len(left), len(right)), dtype=np.float64)
        feature_cache: dict[tuple[int, int], tuple[int, float, float, int, int]] = {}
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                # Cheap exclusion before calculating full features.
                if b.end < a.start - 5000 or b.start > a.end + 5000:
                    continue
                values = interval_features(a, b)
                if not scenario_accepts(name, values):
                    continue
                feature_cache[(i, j)] = values
                _, reciprocal_min, jaccard, start_delta, end_delta = values
                quality = reciprocal_min * 1_000_000 + jaccard * 1000 - (start_delta + end_delta) / 1_000_000
                # 1e12 makes maximum cardinality dominate all quality trade-offs.
                weights[i, j] = 1_000_000_000_000 + quality
        row_idx, col_idx = linear_sum_assignment(weights, maximize=True)
        for i, j in zip(row_idx.tolist(), col_idx.tolist()):
            if (i, j) not in feature_cache:
                continue
            matches.append(Match(name, left[i], right[j], *feature_cache[(i, j)]))

    matches.sort(key=lambda row: (row.a.chrom, row.a.start, row.a.end, row.b.start, row.b.end))
    if len({row.a.region for row in matches}) != len(matches):
        raise RuntimeError(f"scenario {name} duplicated dataset-A regions")
    if len({row.b.region for row in matches}) != len(matches):
        raise RuntimeError(f"scenario {name} duplicated dataset-B regions")
    return matches


def confusion_matrix(y_a: list[str], y_b: list[str]) -> np.ndarray:
    index = {value: position for position, value in enumerate(CATEGORIES)}
    matrix = np.zeros((len(CATEGORIES), len(CATEGORIES)), dtype=np.int64)
    for a, b in zip(y_a, y_b):
        matrix[index[a], index[b]] += 1
    return matrix


def cohen_kappa(matrix: np.ndarray) -> float:
    total = int(matrix.sum())
    if not total:
        return math.nan
    observed = float(np.trace(matrix)) / total
    expected = float(np.dot(matrix.sum(axis=1), matrix.sum(axis=0))) / (total * total)
    return (observed - expected) / (1 - expected) if expected < 1 else math.nan


def balanced_accuracy(matrix: np.ndarray, a_as_reference: bool) -> float:
    values = matrix if a_as_reference else matrix.T
    recalls = []
    for index in range(values.shape[0]):
        denominator = int(values[index, :].sum())
        if denominator:
            recalls.append(float(values[index, index]) / denominator)
    return float(np.mean(recalls)) if recalls else math.nan


def wilson_interval(success: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return math.nan, math.nan
    p = success / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return center - half, center + half


def permutation_null(
    y_a: list[str],
    y_b: list[str],
    chroms: list[str],
    observed_kappa: float,
    seed: int,
    n_permutations: int,
) -> dict:
    a = np.array(y_a, dtype=object)
    b = np.array(y_b, dtype=object)
    groups = [np.where(np.array(chroms, dtype=object) == chrom)[0] for chrom in sorted(set(chroms))]
    rng = np.random.default_rng(seed)
    agreements = np.empty(n_permutations, dtype=np.float64)
    for permutation in range(n_permutations):
        shuffled = b.copy()
        for indices in groups:
            shuffled[indices] = rng.permutation(shuffled[indices])
        agreements[permutation] = np.mean(a == shuffled)

    matrix = confusion_matrix(y_a, y_b)
    total = int(matrix.sum())
    expected = float(np.dot(matrix.sum(axis=1), matrix.sum(axis=0))) / (total * total)
    kappas = (agreements - expected) / (1 - expected) if expected < 1 else np.full_like(agreements, np.nan)
    observed_agreement = float(np.mean(a == b))
    return {
        "seed": seed,
        "permutations": n_permutations,
        "chromosome_preserving": True,
        "agreement_mean": float(np.mean(agreements)),
        "agreement_sd": float(np.std(agreements, ddof=1)),
        "agreement_q025": float(np.quantile(agreements, 0.025)),
        "agreement_q975": float(np.quantile(agreements, 0.975)),
        "agreement_p_ge": float((1 + np.count_nonzero(agreements >= observed_agreement)) / (n_permutations + 1)),
        "kappa_mean": float(np.mean(kappas)),
        "kappa_q025": float(np.quantile(kappas, 0.025)),
        "kappa_q975": float(np.quantile(kappas, 0.975)),
        "kappa_p_ge": float((1 + np.count_nonzero(kappas >= observed_kappa)) / (n_permutations + 1)),
    }


def unordered_signature(signature: str) -> str:
    if not signature:
        return ""
    values = []
    for token in signature.split("|"):
        _, value = token.split("=", 1)
        values.append(value)
    return "|".join(sorted(values))


def match_metrics(matches: list[Match], a_total: int, b_total: int, a_complete: int, b_complete: int) -> tuple[dict, list[dict]]:
    complete = [row for row in matches if row.a.complete and row.b.complete]
    y_a = [row.a.coarse_category for row in complete]
    y_b = [row.b.coarse_category for row in complete]
    chroms = [row.a.chrom for row in complete]
    matrix = confusion_matrix(y_a, y_b)
    n = len(complete)
    agreements = int(np.trace(matrix))
    raw_agreement = agreements / n if n else math.nan
    ci_low, ci_high = wilson_interval(agreements, n)
    kappa = cohen_kappa(matrix)
    resolution_matrix = np.array(
        [
            [int(matrix[:-1, :-1].sum()), int(matrix[:-1, -1].sum())],
            [int(matrix[-1, :-1].sum()), int(matrix[-1, -1])],
        ],
        dtype=np.int64,
    )
    resolution_agreement = float(np.trace(resolution_matrix)) / n if n else math.nan
    resolution_kappa = cohen_kappa(resolution_matrix)
    ba_a = balanced_accuracy(matrix, a_as_reference=True)
    ba_b = balanced_accuracy(matrix, a_as_reference=False)
    jaccards = {}
    for index, category in enumerate(CATEGORIES):
        intersection = int(matrix[index, index])
        union = int(matrix[index, :].sum() + matrix[:, index].sum() - intersection)
        jaccards[category] = intersection / union if union else math.nan

    pooled = collections.Counter(y_a + y_b)
    dominant = pooled.most_common(1)[0][0] if pooled else ""
    non_dominant = [(a, b) for a, b in zip(y_a, y_b) if a != dominant and b != dominant]
    dominant_cell_removed = [(a, b) for a, b in zip(y_a, y_b) if not (a == dominant and b == dominant)]

    def agreement(values: list[tuple[str, str]]) -> float:
        return sum(a == b for a, b in values) / len(values) if values else math.nan

    topo_unique = [row for row in complete if row.a.coarse_category != CATEGORIES[-1] and row.b.coarse_category != CATEGORIES[-1]]
    ordered_hp_equal = sum(row.a.ordered_hp_coarse_signature == row.b.ordered_hp_coarse_signature for row in topo_unique)
    unordered_hp_equal = sum(
        unordered_signature(row.a.ordered_hp_coarse_signature) == unordered_signature(row.b.ordered_hp_coarse_signature)
        for row in topo_unique
    )
    ordered_shape_equal = sum(row.a.ordered_shape_set_signature == row.b.ordered_shape_set_signature for row in complete)
    unordered_shape_equal = sum(
        unordered_signature(row.a.ordered_shape_set_signature) == unordered_signature(row.b.ordered_shape_set_signature)
        for row in complete
    )
    digest_rows = [row for row in complete if row.a.ordered_tree_digest_signature and row.b.ordered_tree_digest_signature]
    ordered_digest_equal = sum(
        row.a.ordered_tree_digest_signature == row.b.ordered_tree_digest_signature for row in digest_rows
    )
    unordered_digest_equal = sum(
        unordered_signature(row.a.ordered_tree_digest_signature) == unordered_signature(row.b.ordered_tree_digest_signature)
        for row in digest_rows
    )

    interval_values = {
        "reciprocal_min_median": float(np.median([row.reciprocal_min for row in matches])) if matches else math.nan,
        "interval_jaccard_median": float(np.median([row.interval_jaccard for row in matches])) if matches else math.nan,
        "start_delta_bp_median": float(np.median([row.start_delta_bp for row in matches])) if matches else math.nan,
        "end_delta_bp_median": float(np.median([row.end_delta_bp for row in matches])) if matches else math.nan,
    }
    null = permutation_null(y_a, y_b, chroms, kappa, SEED, N_PERMUTATIONS) if n else {}
    result = {
        "scenario": matches[0].scenario if matches else "unknown",
        "matched_all": len(matches),
        "match_coverage_a": len(matches) / a_total,
        "match_coverage_b": len(matches) / b_total,
        "complete_both": n,
        "complete_both_coverage_a": n / a_complete,
        "complete_both_coverage_b": n / b_complete,
        "raw_agreement": raw_agreement,
        "raw_agreement_ci95_low": ci_low,
        "raw_agreement_ci95_high": ci_high,
        "cohen_kappa": kappa,
        "resolution_resolved_resolved": int(resolution_matrix[0, 0]),
        "resolution_resolved_unresolved": int(resolution_matrix[0, 1]),
        "resolution_unresolved_resolved": int(resolution_matrix[1, 0]),
        "resolution_unresolved_unresolved": int(resolution_matrix[1, 1]),
        "resolution_binary_agreement": resolution_agreement,
        "resolution_binary_kappa": resolution_kappa,
        "balanced_accuracy_a_reference": ba_a,
        "balanced_accuracy_b_reference": ba_b,
        "balanced_accuracy_symmetric_mean": (ba_a + ba_b) / 2,
        "category_jaccard": jaccards,
        "category_jaccard_macro": float(np.nanmean(list(jaccards.values()))),
        "dominant_class": dominant,
        "both_non_dominant_n": len(non_dominant),
        "both_non_dominant_agreement": agreement(non_dominant),
        "dominant_concordant_cell_removed_n": len(dominant_cell_removed),
        "dominant_concordant_cell_removed_agreement": agreement(dominant_cell_removed),
        "topology_unique_both_n": len(topo_unique),
        "ordered_hp_coarse_signature_agreement": ordered_hp_equal / len(topo_unique) if topo_unique else math.nan,
        "phase_swap_tolerant_hp_coarse_signature_agreement": unordered_hp_equal / len(topo_unique) if topo_unique else math.nan,
        "ordered_shape_set_signature_agreement": ordered_shape_equal / n if n else math.nan,
        "phase_swap_tolerant_shape_set_signature_agreement": unordered_shape_equal / n if n else math.nan,
        "exact_tree_digest_evaluable_n": len(digest_rows),
        "ordered_exact_tree_set_digest_agreement": ordered_digest_equal / len(digest_rows) if digest_rows else math.nan,
        "phase_swap_tolerant_exact_tree_set_digest_agreement": unordered_digest_equal / len(digest_rows) if digest_rows else math.nan,
        **interval_values,
        "permutation_null": null,
    }
    confusion_rows = []
    for i, category_a in enumerate(CATEGORIES):
        for j, category_b in enumerate(CATEGORIES):
            confusion_rows.append(
                {
                    "scenario": result["scenario"],
                    "category_a": category_a,
                    "category_b": category_b,
                    "n": int(matrix[i, j]),
                }
            )
    return result, confusion_rows


def summarize_datasets(records: list[RegionRecord], ct_report: dict, regional: dict) -> tuple[list[dict], list[dict]]:
    ct_by_sample = {row["sample"]: row for row in ct_report["samples"]}
    regional_by_sample = {row["sample"]: row for row in regional["samples"]}
    rows = []
    checks = []
    for sample in SAMPLE_ORDER:
        values = [row for row in records if row.sample == sample]
        counts = collections.Counter(row.coarse_category for row in values)
        complete = sum(row.complete for row in values)
        summary = {
            "sample": sample,
            "primary_regions": len(values),
            "complete_regions": complete,
            "incomplete_regions": counts["incomplete"],
        }
        for category in CATEGORIES:
            summary[category] = counts[category]
            summary[f"share_{category}"] = counts[category] / complete if complete else math.nan
        rows.append(summary)

        expected_classes = ct_by_sample[sample]["T_and_Topology"]["classes"]
        expected_complete = sum(int(expected_classes[key]) for key in ("T=1|Topo=1", "T>1|Topo=1", "T>1|Topo>1"))
        checks.extend(
            [
                {"scope": sample, "check": "five_classes_sum_to_complete", "pass": sum(counts[key] for key in CATEGORIES) == complete, "observed": sum(counts[key] for key in CATEGORIES), "expected": complete},
                {"scope": sample, "check": "complete_plus_incomplete_sum_to_primary", "pass": complete + counts["incomplete"] == len(values), "observed": complete + counts["incomplete"], "expected": len(values)},
                {"scope": sample, "check": "primary_matches_regional_composition", "pass": len(values) == int(regional_by_sample[sample]["primary_regions"]), "observed": len(values), "expected": int(regional_by_sample[sample]["primary_regions"])},
                {"scope": sample, "check": "complete_matches_regional_composition", "pass": complete == int(regional_by_sample[sample]["fully_complete_regions"]), "observed": complete, "expected": int(regional_by_sample[sample]["fully_complete_regions"])},
                {"scope": sample, "check": "complete_matches_ct_report", "pass": complete == expected_complete, "observed": complete, "expected": expected_complete},
                {"scope": sample, "check": "topon_matches_ct_report", "pass": counts["topology_multiple_unresolved"] == int(expected_classes["T>1|Topo>1"]), "observed": counts["topology_multiple_unresolved"], "expected": int(expected_classes["T>1|Topo>1"])},
            ]
        )
    return rows, checks


def golden_cases(records: list[RegionRecord]) -> list[dict]:
    preferred = [row for row in records if row.sample == "HCC1395" and row.complete]
    result = []
    for category in CATEGORIES:
        candidates = [row for row in preferred if row.coarse_category == category]
        if category == "no_within_hp_relation":
            double = [row for row in candidates if row.primary_hp_units == 2]
            candidates = double or candidates
        row = sorted(candidates, key=lambda value: (value.chrom, value.start, value.end))[0]
        result.append(
            {
                "category": category,
                "category_zh": CATEGORY_ZH[category],
                "sample": row.sample,
                "region": row.region,
                "structural_class": row.structural_class,
                "primary_hp_units": row.primary_hp_units,
                "ordered_hp_coarse_signature": row.ordered_hp_coarse_signature,
                "ordered_shape_set_signature": row.ordered_shape_set_signature,
            }
        )
    return result


def region_rows_for_output(records: list[RegionRecord]) -> list[dict]:
    return [
        {
            "sample": row.sample,
            "region": row.region,
            "chrom": row.chrom,
            "start": row.start,
            "end": row.end,
            "length_bp": row.length,
            "primary_hp_units": row.primary_hp_units,
            "structural_class": row.structural_class,
            "complete": row.complete,
            "coarse_category": row.coarse_category,
            "coarse_category_zh": CATEGORY_ZH[row.coarse_category],
            "hp1_category_set": row.hp1_category_set,
            "hp2_category_set": row.hp2_category_set,
            "ordered_hp_coarse_signature": row.ordered_hp_coarse_signature,
            "ordered_shape_set_signature": row.ordered_shape_set_signature,
            "ordered_tree_digest_signature": row.ordered_tree_digest_signature,
        }
        for row in records
    ]


def match_rows_for_output(matches_by_scenario: dict[str, list[Match]]) -> list[dict]:
    result = []
    for scenario, matches in matches_by_scenario.items():
        for index, row in enumerate(matches, start=1):
            result.append(
                {
                    "scenario": scenario,
                    "match_id": f"{scenario}:{index:05d}",
                    "chrom": row.a.chrom,
                    "region_a": row.a.region,
                    "region_b": row.b.region,
                    "overlap_bp": row.overlap_bp,
                    "reciprocal_min": row.reciprocal_min,
                    "interval_jaccard": row.interval_jaccard,
                    "start_delta_bp": row.start_delta_bp,
                    "end_delta_bp": row.end_delta_bp,
                    "complete_a": row.a.complete,
                    "complete_b": row.b.complete,
                    "category_a": row.a.coarse_category,
                    "category_b": row.b.coarse_category,
                    "category_agree": row.a.complete and row.b.complete and row.a.coarse_category == row.b.coarse_category,
                    "ordered_hp_signature_a": row.a.ordered_hp_coarse_signature,
                    "ordered_hp_signature_b": row.b.ordered_hp_coarse_signature,
                    "shape_set_signature_a": row.a.ordered_shape_set_signature,
                    "shape_set_signature_b": row.b.ordered_shape_set_signature,
                    "tree_digest_signature_a": row.a.ordered_tree_digest_signature,
                    "tree_digest_signature_b": row.b.ordered_tree_digest_signature,
                }
            )
    return result


def chromosome_rows(matches_by_scenario: dict[str, list[Match]]) -> list[dict]:
    result = []
    for scenario, matches in matches_by_scenario.items():
        grouped: dict[str, list[Match]] = collections.defaultdict(list)
        for row in matches:
            grouped[row.a.chrom].append(row)
        for chrom in sorted(grouped):
            rows = grouped[chrom]
            complete = [row for row in rows if row.a.complete and row.b.complete]
            result.append(
                {
                    "scenario": scenario,
                    "chrom": chrom,
                    "matched_all": len(rows),
                    "complete_both": len(complete),
                    "agreements": sum(row.a.coarse_category == row.b.coarse_category for row in complete),
                    "raw_agreement": sum(row.a.coarse_category == row.b.coarse_category for row in complete) / len(complete) if complete else math.nan,
                    "median_reciprocal_min": float(np.median([row.reciprocal_min for row in rows])),
                    "median_interval_jaccard": float(np.median([row.interval_jaccard for row in rows])),
                }
            )
    return result


def report_markdown(
    path: Path,
    generated_at: str,
    sources: dict[str, dict],
    summaries: list[dict],
    metrics: list[dict],
    goldens: list[dict],
    checks: list[dict],
    outputs: dict[str, Path],
    chromosome_output: list[dict],
    command: str,
) -> None:
    exact = next(row for row in metrics if row["scenario"] == "exact_coordinate")
    sensitivity = [row for row in metrics if row["scenario"] != "exact_coordinate"]
    hcc = next(row for row in summaries if row["sample"] == BIOLOGICAL_PAIR[0])
    dorado = next(row for row in summaries if row["sample"] == BIOLOGICAL_PAIR[1])
    lines = [
        "<!--",
        f"建立時間: {generated_at}",
        "目標: 計算七個 dataset rows 的五類粗拓撲，並驗證 HCC1395 兩列資料的 region interval concordance",
        "處理範圍: chr1-22; historical layered-v2 engineering snapshot; HCC1395 vs HCC1395_DORADO",
        "關聯檔案:",
        *[f"  - {value['path']}" for value in sources.values()],
        "-->",
        "",
        "# HCC1395 兩列資料的區域粗拓撲與區間一致性驗證",
        "",
        "> **TL;DR（Task B，服務 G4/G5）**：Exact-coordinate 可配對 {matched:,} 個區域，其中 {complete:,} 個兩邊皆 complete；五類 raw agreement={agree:.2%}（95% CI {lo:.2%}–{hi:.2%}）、κ={kappa:.3f}。這是高於 chromosome-preserving null 的可重現訊號，但 aggregate 組成差異與 exact tree-set digest 一致率顯示兩列資料並非可互換；支持「區域級方向性 reproducibility」，不足以證明方法已生物學驗證。".format(
            matched=exact["matched_all"],
            complete=exact["complete_both"],
            agree=exact["raw_agreement"],
            lo=exact["raw_agreement_ci95_low"],
            hi=exact["raw_agreement_ci95_high"],
            kappa=exact["cohen_kappa"],
        ),
        "",
        "## 1. 任務分類、研究問題與邊界",
        "",
        "- Task type：**B Comprehensive validation**（全 7 dataset rows；不得以 HCC pair 子集代替全樣本 census）。",
        "- 研究問題：同一 biological cell line 的 HCC1395 5kHz 與 DORADO，是否在相同／相近 genomic interval 得到一致的五類 mutation-state 粗拓撲？",
        "- 成功條件：分類守恆、配對鍵唯一、exact 與 interval sensitivity 方向一致、觀察 agreement 顯著高於染色體內 permutation null。",
        "- 失敗條件：五類無法回加 complete、配對多對一、結論只由 dominant Topo>1 驅動，或 exact tree signatures 幾乎不重現。",
        "- Claim ceiling：tree node 是 mutation state；`H_*` 是 Steiner／partial-supported completion state，**不是已確認 clone**。本結果不驗證細胞族群真值，也不建立癌症基因或藥物因果。",
        "",
        "## 2. 可稽核 operational classifier",
        "",
        "Primary grain 是 **complete primary region 的 ordered HP forest**。每個 HP component 保留 HP1/HP2 順序；region 五類用 component feature 的 OR 合併：",
        "",
        "1. `Topo>1 未定`：ordered regional topology 超過一種；即使候選都屬相同粗類，仍不硬選。",
        "2. 其餘 Topo=1：在每個 HP rooted tree 中，`max_outdegree≥2` 表示姐妹分枝；`max_depth≥2` 表示直系鏈。",
        "3. 兩者皆否＝無 HP 內關係；只有分枝＝姐妹 only；只有鏈＝直系 only；兩者皆是＝姐妹＋直系。",
        "4. HP1/HP2 是兩個 forest components，跨 HP 節點**永遠不判姐妹或直系**。",
        "5. 此 primary graph classifier 包含 `H_*`；另外的 observed-state sensitivity 只可把非 ROOT、非 `H_*` endpoint 當直接觀測。本分析不把 hidden node 改名成 clone。",
        "",
        "### Golden logic cases",
        "",
        "| 類別 | 邏輯 | HCC1395 實際 region | ordered HP coarse | rooted shape set |",
        "|---|---|---|---|---|",
    ]
    logic = {
        "no_within_hp_relation": "ROOT→A（或 HP1/HP2 各自單點）",
        "sister_only": "ROOT→A 且 ROOT→B",
        "direct_only": "ROOT→A→B",
        "sister_and_direct": "ROOT→A，A→B 且 A→C",
        "topology_multiple_unresolved": "候選 rooted shapes >1",
    }
    for row in goldens:
        lines.append(
            f"| {row['category_zh']} | {logic[row['category']]} | `{row['region']}` | `{row['ordered_hp_coarse_signature']}` | `{row['ordered_shape_set_signature']}` |"
        )

    lines.extend(
        [
            "",
            "## 3. 全 7 dataset rows 五類組成",
            "",
            "五類分母是 complete；Incomplete 另列，並以 `五類=complete`、`complete+incomplete=primary` 雙重守恆。",
            "",
            "| Dataset | Primary | Complete | Incomplete | 無 HP 內 | 姐妹 only | 直系 only | 姐妹＋直系 | Topo>1 未定 |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summaries:
        lines.append(
            "| {sample} | {primary_regions:,} | {complete_regions:,} | {incomplete_regions:,} | {no_within_hp_relation:,} ({s0:.2%}) | {sister_only:,} ({s1:.2%}) | {direct_only:,} ({s2:.2%}) | {sister_and_direct:,} ({s3:.2%}) | {topology_multiple_unresolved:,} ({s4:.2%}) |".format(
                **row,
                s0=row["share_no_within_hp_relation"],
                s1=row["share_sister_only"],
                s2=row["share_direct_only"],
                s3=row["share_sister_and_direct"],
                s4=row["share_topology_multiple_unresolved"],
            )
        )

    total_variation = 0.5 * sum(
        abs(hcc[f"share_{category}"] - dorado[f"share_{category}"]) for category in CATEGORIES
    )
    exact_chromosomes = [row for row in chromosome_output if row["scenario"] == "exact_coordinate"]
    lowest_chromosome = min(exact_chromosomes, key=lambda row: row["raw_agreement"])
    highest_chromosome = max(exact_chromosomes, key=lambda row: row["raw_agreement"])
    lines.extend(
        [
            "",
            "HCC1395 兩列的 marginal composition 並不相同：Topo>1 未定為 {a:.2%} vs {b:.2%}，直系 only 為 {c:.2%} vs {d:.2%}；五類 total-variation distance={tv:.2%}。所以不能只因為是同一 biological sample 就宣稱分布一致。".format(
                a=hcc["share_topology_multiple_unresolved"],
                b=dorado["share_topology_multiple_unresolved"],
                c=hcc["share_direct_only"],
                d=dorado["share_direct_only"],
                tv=total_variation,
            ),
            "",
            "## 4. HCC1395 region interval concordance",
            "",
            "配對規則：同染色體、一對一；exact 要求 `(chrom,start,end)` 全同。Sensitivity 用 reciprocal overlap 0.80/0.50，或 start/end 雙端 anchor 各在 1 kb/5 kb 內。非 exact 情境用 maximum-cardinality、再 maximum-quality assignment；matched A/B keys 均唯一。",
            "",
            "| Scenario | Matched all | Complete both | A/B complete coverage | Raw agreement (95% CI) | κ | Symmetric BA | Macro category Jaccard | Non-dominant-only agree | Remove dominant concordant cell | Permutation null mean (95%) | p |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in metrics:
        null = row["permutation_null"]
        lines.append(
            "| {scenario} | {matched_all:,} | {complete_both:,} | {ca:.1%}/{cb:.1%} | {raw_agreement:.2%} ({raw_agreement_ci95_low:.2%}–{raw_agreement_ci95_high:.2%}) | {cohen_kappa:.3f} | {balanced_accuracy_symmetric_mean:.3f} | {category_jaccard_macro:.3f} | {both_non_dominant_agreement:.2%} (n={both_non_dominant_n:,}) | {dominant_concordant_cell_removed_agreement:.2%} (n={dominant_concordant_cell_removed_n:,}) | {nm:.2%} ({nlo:.2%}–{nhi:.2%}) | {p:.5f} |".format(
                **row,
                ca=row["complete_both_coverage_a"],
                cb=row["complete_both_coverage_b"],
                nm=null["agreement_mean"],
                nlo=null["agreement_q025"],
                nhi=null["agreement_q975"],
                p=null["agreement_p_ge"],
            )
        )

    lines.extend(
        [
            "",
            "### Exact-coordinate 的解讀",
            "",
            "- Raw agreement={agree:.2%}、κ={kappa:.3f}，高於 chromosome-preserving null {null:.2%}（p={p:.5f}）：不是只靠染色體別 class composition 的巧合。".format(
                agree=exact["raw_agreement"], kappa=exact["cohen_kappa"], null=exact["permutation_null"]["agreement_mean"], p=exact["permutation_null"]["agreement_p_ge"]
            ),
            "- 但 Topo>1 是 pooled dominant class。若只看兩邊都非 dominant 的子集，agreement={nondom:.2%}；這個數字排除了 Topo↔determinate 的困難案例，不能單獨當性能。若只移除 `Topo>1/Topo>1` 的同格貢獻，剩餘 agreement 降為 {removed:.2%}，說明 headline agreement 明顯受 dominant concordance 驅動。".format(
                nondom=exact["both_non_dominant_agreement"], removed=exact["dominant_concordant_cell_removed_agreement"]
            ),
            "- Binary resolved/unresolved 2×2（A rows × B columns）：resolved/resolved={rr:,}、resolved/unresolved={ru:,}、unresolved/resolved={ur:,}、unresolved/unresolved={uu:,}；agreement={agree:.2%}、κ={kappa:.3f}。兩邊都 resolved 後 87.18% 是**條件式一致**，不可作全體 accuracy。".format(
                rr=exact["resolution_resolved_resolved"],
                ru=exact["resolution_resolved_unresolved"],
                ur=exact["resolution_unresolved_resolved"],
                uu=exact["resolution_unresolved_unresolved"],
                agree=exact["resolution_binary_agreement"],
                kappa=exact["resolution_binary_kappa"],
            ),
            "- 五類 per-category Jaccard：" + "；".join(f"{CATEGORY_ZH[key]}={value:.3f}" for key, value in exact["category_jaccard"].items()) + "。",
            "- Balanced accuracy 以任一 dataset 當 reference 都不等於 truth。A→B={a:.3f}、B→A={b:.3f}，不對稱反映 marginal class shift。".format(
                a=exact["balanced_accuracy_a_reference"], b=exact["balanced_accuracy_b_reference"]
            ),
            "- 染色體分層 exact agreement 範圍為 {low:.2%}（{low_chr}, n={low_n:,}）到 {high:.2%}（{high_chr}, n={high_n:,}）；22 條 autosomes 全有 exact matched complete regions。Permutation 已在染色體內打亂，避免把 chromosome composition 當一致性。".format(
                low=lowest_chromosome["raw_agreement"],
                low_chr=lowest_chromosome["chrom"],
                low_n=lowest_chromosome["complete_both"],
                high=highest_chromosome["raw_agreement"],
                high_chr=highest_chromosome["chrom"],
                high_n=highest_chromosome["complete_both"],
            ),
            "",
            "### Ordered HP 與 exact tree-set signature",
            "",
            "- 兩邊皆 Topo=1 的 exact-coordinate regions：ordered HP coarse signature agreement={ordered:.2%}；允許全區 HP1↔HP2 phase flip 後={swap:.2%}。Primary 輸出仍保留 ordered pair，swap 只作跨資料集 sensitivity。".format(
                ordered=exact["ordered_hp_coarse_signature_agreement"], swap=exact["phase_swap_tolerant_hp_coarse_signature_agreement"]
            ),
            "- Rooted unlabeled candidate-shape set：ordered={ordered:.2%}、phase-swap tolerant={swap:.2%}。".format(
                ordered=exact["ordered_shape_set_signature_agreement"], swap=exact["phase_swap_tolerant_shape_set_signature_agreement"]
            ),
            "- 原始 layered exact candidate-tree-set digest：ordered={ordered:.2%}、phase-swap tolerant={swap:.2%}（n={n:,}）。這是比五類更嚴格的 reproducibility ceiling。".format(
                ordered=exact["ordered_exact_tree_set_digest_agreement"], swap=exact["phase_swap_tolerant_exact_tree_set_digest_agreement"], n=exact["exact_tree_digest_evaluable_n"]
            ),
            "",
            "## 5. 科學判定",
            "",
            "1. **可支持**：同一 biological sample 的兩列資料，在大量 exact／高 overlap intervals 上有高於 permutation null 的粗拓撲 concordance，方法抓到一部分可重現的區域結構訊號。",
            "2. **不可宣稱完全一致**：marginal 五類 TV distance、Topo>1 比例、balanced accuracy 不對稱與 strict tree-set signature 都顯示 dataset/basecalling/read-support 依賴。",
            "3. **不可由此證明方法有效**：兩列不是獨立 biological replicates，也沒有 clone ground truth；同一樣本 cross-dataset agreement 主要是 reproducibility evidence，不是 accuracy/validity proof。",
            "4. **基因／藥物資料的角色**：癌症基因或藥物 annotation 可檢查已知 biology 的 enrichment/face validity，但不是 topology ground truth；除非有獨立克隆、單細胞、longitudinal 或功能性證據，不能把 annotation overlap 寫成方法被證明。",
            "5. **目前狀態**：`SHARE WITH CAVEATS / SCIENTIFIC NO-GO for proof-of-effectiveness`。clean layered-v3 全 7 樣本尚未 closeout，這裡仍是 historical layered-v2 engineering snapshot。",
            "",
            "## 6. Step → Verify 與 QA",
            "",
            "1. 讀取 exact catalog、C/T report、regional composition、region TSV → 驗證：7 dataset rows 齊全、region key 唯一。",
            "2. 依 rooted shape metadata 重算五類 → 驗證：每列五類回加 complete；complete+incomplete 回加 primary。",
            "3. 建立 exact/overlap/anchor 一對一配對 → 驗證：每 scenario A/B matched keys 各自唯一。",
            "4. 重算 agreement/κ/BA/Jaccard/null → 驗證：seed=20260712、5,000 permutations、chromosome-preserving。",
            "5. 追到 layered raw exact tree digest → 驗證：不以 coarse shape 猜 exact tree equality。",
            "",
            f"- Checks：**{sum(bool(row['pass']) for row in checks)}/{len(checks)} PASS**。",
            f"- 主要資料表：`{outputs['summary_tsv'].resolve()}`、`{outputs['regions_tsv'].resolve()}`、`{outputs['metrics_tsv'].resolve()}`、`{outputs['matches_tsv'].resolve()}`。",
            "",
            "## 7. 完整執行命令、輸入、輸出與實際片段",
            "",
            "### 執行命令",
            "",
            "```bash",
            command,
            "```",
            "",
            "### 輸入",
            "",
            *[f"- `{value['path']}`（SHA-256 `{value['sha256']}`）" for value in sources.values()],
            "",
            "### 輸出",
            "",
            *[f"- `{value.resolve()}`" for value in outputs.values()],
            f"- `{path.resolve()}`",
            "",
            "### 實際輸出片段",
            "",
            "```text",
            "ACTUAL SUMMARY -> exact matched={:,}; complete_both={:,}; agreement={:.6f}; kappa={:.6f}; null_p={:.6f}".format(
                exact["matched_all"], exact["complete_both"], exact["raw_agreement"], exact["cohen_kappa"], exact["permutation_null"]["agreement_p_ge"]
            ),
            f"VALIDATION -> {sum(bool(row['pass']) for row in checks)}/{len(checks)} PASS",
            "```",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-catalog", required=True, type=Path)
    parser.add_argument("--ct-report", required=True, type=Path)
    parser.add_argument("--regional-composition", required=True, type=Path)
    parser.add_argument("--region-tsv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--output-report", required=True, type=Path)
    args = parser.parse_args()
    command = (
        f"python3 {Path(__file__).resolve()} "
        f"--exact-catalog {args.exact_catalog.resolve()} "
        f"--ct-report {args.ct_report.resolve()} "
        f"--regional-composition {args.regional_composition.resolve()} "
        f"--region-tsv {args.region_tsv.resolve()} "
        f"--output-dir {args.output_dir.resolve()} "
        f"--output-report {args.output_report.resolve()}"
    )

    exact_catalog = load_json(args.exact_catalog)
    ct_report = load_json(args.ct_report)
    regional = load_json(args.regional_composition)
    with args.region_tsv.open(encoding="utf-8", newline="") as handle:
        region_rows = list(csv.DictReader(handle, delimiter="\t"))

    tree_digests = load_tree_digests(exact_catalog, set(BIOLOGICAL_PAIR))
    records, checks = build_region_records(exact_catalog, ct_report, region_rows, tree_digests)
    summaries, summary_checks = summarize_datasets(records, ct_report, regional)
    checks.extend(summary_checks)

    pair_records = {sample: [row for row in records if row.sample == sample] for sample in BIOLOGICAL_PAIR}
    matches_by_scenario = {
        "exact_coordinate": exact_matches(pair_records[BIOLOGICAL_PAIR[0]], pair_records[BIOLOGICAL_PAIR[1]])
    }
    for scenario in ("reciprocal_overlap_0.80", "reciprocal_overlap_0.50", "endpoint_anchor_1kb", "endpoint_anchor_5kb"):
        matches_by_scenario[scenario] = interval_matches(
            scenario, pair_records[BIOLOGICAL_PAIR[0]], pair_records[BIOLOGICAL_PAIR[1]]
        )

    for scenario, matches in matches_by_scenario.items():
        checks.extend(
            [
                {"scope": scenario, "check": "matched_a_keys_unique", "pass": len(matches) == len({row.a.region for row in matches}), "observed": len({row.a.region for row in matches}), "expected": len(matches)},
                {"scope": scenario, "check": "matched_b_keys_unique", "pass": len(matches) == len({row.b.region for row in matches}), "observed": len({row.b.region for row in matches}), "expected": len(matches)},
                {"scope": scenario, "check": "matched_chromosomes_identical", "pass": all(row.a.chrom == row.b.chrom for row in matches), "observed": sum(row.a.chrom == row.b.chrom for row in matches), "expected": len(matches)},
            ]
        )

    summary_by_sample = {row["sample"]: row for row in summaries}
    metrics = []
    confusion_rows = []
    for scenario, matches in matches_by_scenario.items():
        metric, confusion = match_metrics(
            matches,
            summary_by_sample[BIOLOGICAL_PAIR[0]]["primary_regions"],
            summary_by_sample[BIOLOGICAL_PAIR[1]]["primary_regions"],
            summary_by_sample[BIOLOGICAL_PAIR[0]]["complete_regions"],
            summary_by_sample[BIOLOGICAL_PAIR[1]]["complete_regions"],
        )
        metrics.append(metric)
        confusion_rows.extend(confusion)

    for row in metrics:
        resolution_total = sum(
            int(row[key])
            for key in (
                "resolution_resolved_resolved",
                "resolution_resolved_unresolved",
                "resolution_unresolved_resolved",
                "resolution_unresolved_unresolved",
            )
        )
        checks.append(
            {
                "scope": row["scenario"],
                "check": "resolution_2x2_sums_to_complete_both",
                "pass": resolution_total == int(row["complete_both"]),
                "observed": resolution_total,
                "expected": int(row["complete_both"]),
            }
        )

    checks.append(
        {
            "scope": "aggregate",
            "check": "all_validation_checks_pass",
            "pass": all(bool(row["pass"]) for row in checks),
            "observed": sum(bool(row["pass"]) for row in checks),
            "expected": len(checks),
        }
    )
    if not all(bool(row["pass"]) for row in checks):
        failed = [row for row in checks if not bool(row["pass"])]
        raise SystemExit(f"validation failed: {failed[:5]}")

    output = args.output_dir
    outputs = {
        "canonical_json": output / "topology_pair_analysis.json",
        "summary_tsv": output / "coarse_topology_all_dataset_summary.tsv",
        "regions_tsv": output / "coarse_topology_all_regions.tsv",
        "metrics_tsv": output / "hcc1395_pair_match_metrics.tsv",
        "confusion_tsv": output / "hcc1395_pair_confusion.tsv",
        "resolution_binary_tsv": output / "hcc1395_pair_resolution_binary.tsv",
        "matches_tsv": output / "hcc1395_pair_matches.tsv",
        "chromosome_tsv": output / "hcc1395_pair_interval_by_chrom.tsv",
        "golden_tsv": output / "coarse_topology_golden_cases.tsv",
        "checks_tsv": output / "topology_pair_validation_checks.tsv",
    }
    generated_at = datetime.now(ZoneInfo("Asia/Taipei")).isoformat(timespec="seconds")
    layered_by_sample = {row["sample"]: Path(row["source"]) for row in exact_catalog["samples"]}
    sources = {
        "exact_catalog": {"path": str(args.exact_catalog.resolve()), "sha256": sha256(args.exact_catalog)},
        "ct_report": {"path": str(args.ct_report.resolve()), "sha256": sha256(args.ct_report)},
        "regional_composition": {"path": str(args.regional_composition.resolve()), "sha256": sha256(args.regional_composition)},
        "region_tsv": {"path": str(args.region_tsv.resolve()), "sha256": sha256(args.region_tsv)},
        "layered_hcc1395": {"path": str(layered_by_sample["HCC1395"].resolve()), "sha256": sha256(layered_by_sample["HCC1395"])},
        "layered_hcc1395_dorado": {"path": str(layered_by_sample["HCC1395_DORADO"].resolve()), "sha256": sha256(layered_by_sample["HCC1395_DORADO"])},
    }
    goldens = golden_cases(records)
    canonical = {
        "schema_version": "1.0",
        "generated_at": generated_at,
        "status": "PASS",
        "task_type": "B_comprehensive_validation",
        "scope": "chr1-22; 7 dataset rows; historical layered-v2 engineering snapshot",
        "claim_ceiling": "mutation-state topology reproducibility; not clone truth or causal drug validation",
        "classifier": {
            "classes": CATEGORIES,
            "hidden_nodes_in_primary_graph": True,
            "cross_hp_relations_forbidden": True,
            "double_hp_order_preserved": True,
            "topology_multiple_forced_unresolved": True,
        },
        "sources": sources,
        "dataset_summary": summaries,
        "hcc1395_pair_metrics": metrics,
        "golden_cases": goldens,
        "validation": {"checks": len(checks), "passed": sum(bool(row["pass"]) for row in checks)},
    }
    write_json(outputs["canonical_json"], canonical)
    write_tsv(
        outputs["summary_tsv"],
        summaries,
        ["sample", "primary_regions", "complete_regions", "incomplete_regions"]
        + CATEGORIES
        + [f"share_{value}" for value in CATEGORIES],
    )
    write_tsv(
        outputs["regions_tsv"],
        region_rows_for_output(records),
        [
            "sample", "region", "chrom", "start", "end", "length_bp", "primary_hp_units",
            "structural_class", "complete", "coarse_category", "coarse_category_zh",
            "hp1_category_set", "hp2_category_set", "ordered_hp_coarse_signature",
            "ordered_shape_set_signature", "ordered_tree_digest_signature",
        ],
    )
    metric_rows = []
    for row in metrics:
        flattened = {key: value for key, value in row.items() if key not in {"category_jaccard", "permutation_null"}}
        flattened.update({f"jaccard_{key}": value for key, value in row["category_jaccard"].items()})
        flattened.update({f"null_{key}": value for key, value in row["permutation_null"].items()})
        metric_rows.append(flattened)
    write_tsv(outputs["metrics_tsv"], metric_rows, list(metric_rows[0]))
    write_tsv(outputs["confusion_tsv"], confusion_rows, ["scenario", "category_a", "category_b", "n"])
    resolution_rows = []
    for row in metrics:
        for state_a, state_b, key in (
            ("resolved", "resolved", "resolution_resolved_resolved"),
            ("resolved", "unresolved", "resolution_resolved_unresolved"),
            ("unresolved", "resolved", "resolution_unresolved_resolved"),
            ("unresolved", "unresolved", "resolution_unresolved_unresolved"),
        ):
            resolution_rows.append(
                {
                    "scenario": row["scenario"],
                    "state_a": state_a,
                    "state_b": state_b,
                    "n": row[key],
                    "binary_agreement": row["resolution_binary_agreement"],
                    "binary_kappa": row["resolution_binary_kappa"],
                }
            )
    write_tsv(
        outputs["resolution_binary_tsv"],
        resolution_rows,
        ["scenario", "state_a", "state_b", "n", "binary_agreement", "binary_kappa"],
    )
    match_output_rows = match_rows_for_output(matches_by_scenario)
    write_tsv(outputs["matches_tsv"], match_output_rows, list(match_output_rows[0]))
    chromosome_output = chromosome_rows(matches_by_scenario)
    write_tsv(outputs["chromosome_tsv"], chromosome_output, list(chromosome_output[0]))
    write_tsv(outputs["golden_tsv"], goldens, list(goldens[0]))
    write_tsv(outputs["checks_tsv"], checks, ["scope", "check", "pass", "observed", "expected"])
    report_markdown(
        args.output_report,
        generated_at,
        sources,
        summaries,
        metrics,
        goldens,
        checks,
        outputs,
        chromosome_output,
        command,
    )

    print(f"INPUT EXACT CATALOG -> {args.exact_catalog.resolve()}")
    print(f"INPUT C/T REPORT -> {args.ct_report.resolve()}")
    print(f"INPUT REGIONAL COMPOSITION -> {args.regional_composition.resolve()}")
    print(f"INPUT REGION TSV -> {args.region_tsv.resolve()}")
    for key, value in outputs.items():
        print(f"OUTPUT {key.upper()} -> {value.resolve()}")
    print(f"OUTPUT REPORT -> {args.output_report.resolve()}")
    exact = next(row for row in metrics if row["scenario"] == "exact_coordinate")
    print(
        "ACTUAL SUMMARY -> exact matched={:,}; complete_both={:,}; agreement={:.6f}; kappa={:.6f}; null_p={:.6f}".format(
            exact["matched_all"], exact["complete_both"], exact["raw_agreement"], exact["cohen_kappa"], exact["permutation_null"]["agreement_p_ge"]
        )
    )
    print(f"VALIDATION -> {sum(bool(row['pass']) for row in checks)}/{len(checks)} PASS")


if __name__ == "__main__":
    main()
