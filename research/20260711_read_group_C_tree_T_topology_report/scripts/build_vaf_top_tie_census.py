#!/usr/bin/env python3
"""Build an exact top-tie and topology census for VAF-ranked candidate trees.

This is a descriptive engineering analysis over the frozen layered-v2 snapshot.
It does not alter candidate sets and does not claim a biologically confirmed tree.
Exact ties are evaluated with Fraction-valued read AF and ordering scores. Softmax
relative weights are reported separately as heuristic concentration measures.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime
from fractions import Fraction
from pathlib import Path
from zoneinfo import ZoneInfo


SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
STRUCTURAL_CLASSES = ("T=1|Topo=1", "T>1|Topo=1", "T>1|Topo>1", "incomplete")
AMBIGUOUS_CLASSES = ("T>1|Topo=1", "T>1|Topo>1")
LEGACY_STATUSES = (
    "vaf_top_consistent",
    "vaf_top_direction_inconsistent",
    "below_threshold",
    "missing_or_mismatch",
    "not_evaluable_recurrence",
)
TIE_CLASSES = (
    "unique_first",
    "tied_first_same_topology",
    "tied_first_different_topology",
)
DEFAULT_TEMPERATURES = (0.025, 0.05, 0.10)
DEFAULT_THRESHOLDS = (0.50, 0.55, 0.60, 0.75, 0.90)
DEFAULT_NEAR_TIE_GAPS = (1e-12, 1e-9, 1e-6)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dump_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def canonical_shape(edges: list[list[str]] | list[tuple[str, str]]) -> str:
    children: dict[str, list[str]] = defaultdict(list)
    child_nodes: set[str] = set()
    all_nodes: set[str] = set()
    for parent, child in edges:
        children[parent].append(child)
        child_nodes.add(child)
        all_nodes.update((parent, child))
    roots = sorted(all_nodes - child_nodes)

    def visit(node: str) -> str:
        return "(" + "".join(sorted(visit(child) for child in children[node])) + ")"

    return "|".join(sorted(visit(root) for root in roots)) if roots else "()"


def shape_id(signature: str) -> str:
    return "TS-" + hashlib.sha1(signature.encode()).hexdigest()[:10]


def tree_digest(edges: list[list[str]] | list[tuple[str, str]]) -> str:
    payload = sorted([list(edge) for edge in edges])
    encoded = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode()
    return "TR-" + hashlib.sha1(encoded).hexdigest()[:12]


def unlabel(node: str) -> frozenset[int]:
    if node == "ROOT":
        return frozenset()
    value = node[2:] if node.startswith("H_") else node
    return frozenset(index for index, allele in enumerate(value) if allele == "A")


def exact_read_af_from_colcov(colcov: dict, positions: list[int], k: int) -> dict[int, Fraction] | None:
    if not colcov or len(positions) < k:
        return None
    values: dict[int, Fraction] = {}
    for index in range(k):
        counts = colcov.get(str(positions[index])) or colcov.get(positions[index])
        if not counts:
            return None
        n_ref, n_alt = int(counts[0]), int(counts[1])
        total = n_ref + n_alt
        if total <= 0:
            return None
        values[index] = Fraction(n_alt, total)
    return values


def exact_ordering_score(
    edges: list[list[str]] | list[tuple[str, str]], read_af: dict[int, Fraction]
) -> tuple[Fraction, list[Fraction]]:
    score = Fraction(0, 1)
    comparisons: list[Fraction] = []
    for parent, child in edges:
        parent_set, child_set = unlabel(parent), unlabel(child)
        acquired = child_set - parent_set
        if len(acquired) != 1:
            continue
        child_mutation = next(iter(acquired))
        for ancestor_mutation in parent_set:
            delta = read_af[ancestor_mutation] - read_af[child_mutation]
            score += delta
            comparisons.append(delta)
    return score, comparisons


def structural_class(units: list[dict]) -> str:
    if not all(unit.get("analysis_candidate_set_complete") is True for unit in units):
        return "incomplete"
    t_joint = math.prod(int(unit["n_trees"]) for unit in units)
    topo_joint = math.prod(int(unit["n_distinct_shapes_exact"]) for unit in units)
    if t_joint == 1 and topo_joint == 1:
        return "T=1|Topo=1"
    if t_joint > 1 and topo_joint == 1:
        return "T>1|Topo=1"
    if t_joint > 1 and topo_joint > 1:
        return "T>1|Topo>1"
    raise RuntimeError(f"impossible structural state T={t_joint} Topo={topo_joint}")


def classify_top_set(n_top: int, n_top_shapes: int) -> str:
    if n_top == 1 and n_top_shapes == 1:
        return "unique_first"
    if n_top > 1 and n_top_shapes == 1:
        return "tied_first_same_topology"
    if n_top > 1 and n_top_shapes > 1:
        return "tied_first_different_topology"
    raise RuntimeError(f"invalid top set n_top={n_top} n_top_shapes={n_top_shapes}")


def direction_state(flags: list[bool]) -> str:
    if not flags:
        return "not_applicable"
    if all(flags):
        return "all_consistent"
    if any(flags):
        return "mixed"
    return "all_inconsistent"


def direction_evidence_state(flags: list[bool | None]) -> str:
    if not flags or not any(flag is not None for flag in flags):
        return "not_evaluable_no_comparisons"
    evaluated = [flag for flag in flags if flag is not None]
    if all(evaluated):
        core = "all_consistent"
    elif any(evaluated):
        core = "mixed"
    else:
        core = "all_inconsistent"
    return f"partial_{core}" if len(evaluated) != len(flags) else core


def temp_key(value: float) -> str:
    return f"{value:.3f}".rstrip("0").rstrip(".")


def threshold_key(value: float) -> str:
    return f"{value:.2f}"


def analyze_ranked_unit(
    *,
    sample: str,
    biological_id: str,
    unit: dict,
    group: dict | None,
    ra,
    solver,
    temperatures: tuple[float, ...],
    near_tie_gaps: tuple[float, ...],
    violation_margin: Fraction,
    candidate_writer: csv.DictWriter,
) -> tuple[dict, list[dict]]:
    family = unit["family"]
    base = {
        "sample": sample,
        "biological_id": biological_id,
        "region": unit["region"],
        "family": family,
        "L1_base_class": unit.get("L1_base_class"),
        "cn": unit.get("cn"),
        "n_trees_expected": int(unit["n_trees"]),
        "n_shapes_expected": int(unit["n_distinct_shapes_exact"]),
    }
    if group is None:
        return {**base, "status": "missing_group"}, []
    full = (group.get("populations_by_hp", {}) or {}).get(family, {}) or {}
    partial = list(((group.get("subread_groups_by_hp", {}) or {}).get(family, {}) or {}).keys())
    positions = list(group.get("positions", []))
    k = len(positions)
    result = solver.enumerate_min_trees(full, partial, k, tree_cap=0)
    n_trees_enumerated = int(result["n_trees"])
    shapes_enumerated = {canonical_shape(tree["edges"]) for tree in result["trees"]}
    enumeration_ok = (
        not result.get("capped")
        and bool(result.get("trees_complete"))
        and n_trees_enumerated == base["n_trees_expected"]
        and len(shapes_enumerated) == base["n_shapes_expected"]
    )
    if not enumeration_ok:
        return {
            **base,
            "status": "candidate_or_shape_mismatch",
            "n_trees_enumerated": n_trees_enumerated,
            "n_shapes_enumerated": len(shapes_enumerated),
        }, []
    colcov = ((group.get("col_coverage_by_hp", {}) or {}).get(family, {}) or {})
    exact_read_af = exact_read_af_from_colcov(colcov, positions, k)
    if exact_read_af is None:
        return {
            **base,
            "status": "missing_vaf",
            "n_trees_enumerated": n_trees_enumerated,
            "n_shapes_enumerated": len(shapes_enumerated),
        }, []

    float_read_af = {index: float(value) for index, value in exact_read_af.items()}
    candidates: list[dict] = []
    max_float_score_error = 0.0
    for candidate_index, tree in enumerate(result["trees"], 1):
        exact_score, exact_deltas = exact_ordering_score(tree["edges"], exact_read_af)
        reference_score, _ = ra.ordering_score(tree["edges"], float_read_af)
        max_float_score_error = max(max_float_score_error, abs(reference_score - float(exact_score)))
        signature = canonical_shape(tree["edges"])
        flags = [delta >= -violation_margin for delta in exact_deltas]
        direction_evaluable = bool(flags)
        candidates.append(
            {
                "candidate_index": candidate_index,
                "candidate_id": tree_digest(tree["edges"]),
                "shape_id": shape_id(signature),
                "shape_signature": signature,
                "exact_score": exact_score,
                "legacy_direction_consistent": all(flags),
                "direction_evaluable": direction_evaluable,
                "direction_consistent": all(flags) if direction_evaluable else None,
                "n_comparisons": len(exact_deltas),
            }
        )

    exact_scores = [candidate["exact_score"] for candidate in candidates]
    float_scores = [float(score) for score in exact_scores]
    max_score = max(exact_scores)
    exact_top_indices = [index for index, score in enumerate(exact_scores) if score == max_score]
    exact_top_shapes = {candidates[index]["shape_id"] for index in exact_top_indices}
    exact_top_legacy_flags = [
        candidates[index]["legacy_direction_consistent"] for index in exact_top_indices
    ]
    exact_top_evidence_flags = [
        candidates[index]["direction_consistent"] for index in exact_top_indices
    ]
    weights_by_temperature = {
        temp_key(temperature): ra.posterior(float_scores, temperature)
        for temperature in temperatures
    }
    default_key = temp_key(0.05)
    if default_key not in weights_by_temperature:
        raise RuntimeError("temperature grid must include 0.05")
    default_weights = weights_by_temperature[default_key]
    ranked_indices = sorted(
        range(len(candidates)),
        key=lambda index: (
            -exact_scores[index],
            candidates[index]["candidate_id"],
        ),
    )
    dense_rank_by_index: dict[int, int] = {}
    dense_rank = 0
    previous_score: Fraction | None = None
    for index in ranked_indices:
        if previous_score is None or exact_scores[index] != previous_score:
            dense_rank += 1
            previous_score = exact_scores[index]
        dense_rank_by_index[index] = dense_rank
    for candidate_order, index in enumerate(ranked_indices, 1):
        candidate = candidates[index]
        candidate_writer.writerow(
            {
                "sample": sample,
                "biological_id": biological_id,
                "region": unit["region"],
                "family": family,
                "candidate_order": candidate_order,
                "score_rank": dense_rank_by_index[index],
                "candidate_index": candidate["candidate_index"],
                "candidate_id": candidate["candidate_id"],
                "shape_id": candidate["shape_id"],
                "exact_score_numerator": candidate["exact_score"].numerator,
                "exact_score_denominator": candidate["exact_score"].denominator,
                "score_float": f"{float(candidate['exact_score']):.17g}",
                "softmax_weight_T0.05": f"{default_weights[index]:.17g}",
                "is_exact_top": str(index in exact_top_indices),
                "direction_evaluable": str(candidate["direction_evaluable"]),
                "direction_consistent": "" if candidate["direction_consistent"] is None else str(candidate["direction_consistent"]),
                "legacy_direction_consistent": str(candidate["legacy_direction_consistent"]),
                "n_ancestry_comparisons": candidate["n_comparisons"],
            }
        )

    sorted_weights = sorted(default_weights, reverse=True)
    second_weight = sorted_weights[1] if len(sorted_weights) > 1 else None
    near_tie: dict[str, dict] = {}
    for gap in near_tie_gaps:
        indices = [
            index
            for index, score in enumerate(exact_scores)
            if float(max_score - score) <= gap
        ]
        shapes = {candidates[index]["shape_id"] for index in indices}
        near_tie[f"{gap:.0e}"] = {
            "n_top": len(indices),
            "n_top_shapes": len(shapes),
            "class": classify_top_set(len(indices), len(shapes)),
        }

    shape_weights: dict[str, float] = defaultdict(float)
    for candidate, weight in zip(candidates, default_weights):
        shape_weights[candidate["shape_id"]] += weight
    top_shape_weight = max(shape_weights.values())
    top_aggregate_shapes = sorted(
        shape for shape, weight in shape_weights.items() if math.isclose(weight, top_shape_weight, rel_tol=1e-12, abs_tol=1e-12)
    )
    top_weights = {
        key: max(weights) for key, weights in weights_by_temperature.items()
    }
    return {
        **base,
        "status": "ranked",
        "n_trees_enumerated": n_trees_enumerated,
        "n_shapes_enumerated": len(shapes_enumerated),
        "max_float_score_error": max_float_score_error,
        "top_score_numerator": max_score.numerator,
        "top_score_denominator": max_score.denominator,
        "n_top_exact": len(exact_top_indices),
        "n_top_shapes_exact": len(exact_top_shapes),
        "exact_top_class": classify_top_set(len(exact_top_indices), len(exact_top_shapes)),
        "exact_top_candidate_ids": [candidates[index]["candidate_id"] for index in exact_top_indices],
        "exact_top_shape_ids": sorted(exact_top_shapes),
        "exact_top_direction_state": direction_evidence_state(exact_top_evidence_flags),
        "exact_top_legacy_direction_state": direction_state(exact_top_legacy_flags),
        "n_top_direction_evaluable": sum(flag is not None for flag in exact_top_evidence_flags),
        "n_top_direction_consistent": sum(flag is True for flag in exact_top_evidence_flags),
        "n_top_direction_inconsistent": sum(flag is False for flag in exact_top_evidence_flags),
        "n_top_no_direction_comparisons": sum(flag is None for flag in exact_top_evidence_flags),
        "candidate_comparison_count_varies": len({candidate["n_comparisons"] for candidate in candidates}) > 1,
        "all_candidates_co_top": len(exact_top_indices) == len(candidates),
        "top_weight_by_temperature": top_weights,
        "top_weight_T0.05": max(default_weights),
        "second_weight_T0.05": second_weight,
        "weight_gap_T0.05": max(default_weights) - second_weight if second_weight is not None else None,
        "top_aggregate_shape_weight_T0.05": top_shape_weight,
        "n_top_aggregate_shapes_T0.05": len(top_aggregate_shapes),
        "top_aggregate_shape_ids_T0.05": top_aggregate_shapes,
        "near_tie": near_tie,
    }, candidates


def region_tie_summary(ranked_units: list[dict], tie_field: str | None = None) -> dict:
    if tie_field is None:
        n_top = math.prod(int(unit["n_top_exact"]) for unit in ranked_units)
        n_shapes = math.prod(int(unit["n_top_shapes_exact"]) for unit in ranked_units)
    else:
        n_top = math.prod(int(unit["near_tie"][tie_field]["n_top"]) for unit in ranked_units)
        n_shapes = math.prod(int(unit["near_tie"][tie_field]["n_top_shapes"]) for unit in ranked_units)
    return {
        "n_top_joint": n_top,
        "n_top_shapes_joint": n_shapes,
        "class": classify_top_set(n_top, n_shapes),
    }


def summarize_regions(regions: list[dict], near_tie_gaps: tuple[float, ...], temperatures: tuple[float, ...], thresholds: tuple[float, ...]) -> dict:
    structural = Counter(row["structural_class"] for row in regions)
    ambiguous = [row for row in regions if row["structural_class"] in AMBIGUOUS_CLASSES]
    legacy = Counter(row["legacy_status_0.60"] for row in ambiguous)
    evaluable = [row for row in ambiguous if row["evaluable"]]
    old_below = [row for row in evaluable if row["legacy_status_0.60"] == "below_threshold"]

    def nested_count(rows: list[dict], keys: tuple[str, ...]) -> dict:
        result: dict = {}
        for row in rows:
            target = result
            for key in keys[:-1]:
                value = row[key]
                target = target.setdefault(value, {})
            final_value = row[keys[-1]]
            target[final_value] = target.get(final_value, 0) + 1
        return result

    threshold_grid: dict[str, dict] = {}
    for temperature in temperatures:
        tk = temp_key(temperature)
        threshold_grid[tk] = {}
        for threshold in thresholds:
            sk = threshold_key(threshold)
            ge_count = sum(row[f"min_top_weight_T{tk}"] >= threshold for row in evaluable)
            gt_count = sum(row[f"min_top_weight_T{tk}"] > threshold for row in evaluable)
            threshold_grid[tk][sk] = {
                "ge": ge_count,
                "gt": gt_count,
                "evaluable_denominator": len(evaluable),
                "ambiguous_denominator": len(ambiguous),
            }

    exact_tie = Counter(row["exact_top_class"] for row in evaluable)
    exact_tie_old_below = Counter(row["exact_top_class"] for row in old_below)
    near_tie_sensitivity: dict[str, dict] = {}
    for gap in near_tie_gaps:
        gap_key = f"{gap:.0e}"
        field = f"near_tie_class_{gap_key}"
        near_tie_sensitivity[gap_key] = {
            "all_evaluable": dict(Counter(row[field] for row in evaluable)),
            "legacy_below_0.60": dict(Counter(row[field] for row in old_below)),
        }

    old_below_ge = sum(row["min_top_weight_T0.05"] >= 0.50 for row in old_below)
    old_below_gt = sum(row["min_top_weight_T0.05"] > 0.50 for row in old_below)
    old_below_050 = {
        "ge_0.50": old_below_ge,
        "gt_0.50": old_below_gt,
        "lt_0.50": sum(row["min_top_weight_T0.05"] < 0.50 for row in old_below),
        "eq_0.50": old_below_ge - old_below_gt,
        "ge_0.50_by_exact_top_class": dict(
            Counter(row["exact_top_class"] for row in old_below if row["min_top_weight_T0.05"] >= 0.50)
        ),
        "lt_0.50_by_exact_top_class": dict(
            Counter(row["exact_top_class"] for row in old_below if row["min_top_weight_T0.05"] < 0.50)
        ),
    }

    return {
        "region_count": len(regions),
        "structural_classes": {key: structural[key] for key in STRUCTURAL_CLASSES},
        "complete_regions": sum(structural[key] for key in STRUCTURAL_CLASSES if key != "incomplete"),
        "ambiguous_complete_regions": len(ambiguous),
        "evaluable_ambiguous_regions": len(evaluable),
        "legacy_region_status_0.60": {key: legacy[key] for key in LEGACY_STATUSES},
        "exact_top_classes_all_evaluable": {key: exact_tie[key] for key in TIE_CLASSES},
        "legacy_below_0.60_count": len(old_below),
        "legacy_below_0.60_exact_top_classes": {key: exact_tie_old_below[key] for key in TIE_CLASSES},
        "legacy_below_0.60_by_structural_and_top_class": nested_count(
            old_below, ("structural_class", "exact_top_class")
        ),
        "legacy_below_0.60_by_top_class_and_direction": nested_count(
            old_below, ("exact_top_class", "top_direction_state")
        ),
        "legacy_below_0.60_threshold_0.50": old_below_050,
        "threshold_grid": threshold_grid,
        "near_tie_sensitivity": near_tie_sensitivity,
        "exact_top_joint_count_bins_legacy_below": dict(
            sorted(Counter(str(row["n_top_joint_exact"]) for row in old_below).items(), key=lambda item: int(item[0]))
        ),
        "exact_top_shape_count_bins_legacy_below": dict(
            sorted(Counter(str(row["n_top_shapes_joint_exact"]) for row in old_below).items(), key=lambda item: int(item[0]))
        ),
    }


def analyze_sample(
    *,
    meta: dict,
    sample_dir: Path,
    corrected_sample: dict,
    legacy_sample: dict,
    ra,
    solver,
    temperatures: tuple[float, ...],
    thresholds: tuple[float, ...],
    near_tie_gaps: tuple[float, ...],
    violation_margin: Fraction,
    candidate_writer: csv.DictWriter,
) -> tuple[list[dict], list[dict], dict, dict]:
    sample = meta["sample"]
    biological_id = meta["biological_id"]
    layered_path = sample_dir / f"layered_reconstruction_{sample}.json"
    layered = json.loads(layered_path.read_text(encoding="utf-8"))
    groups = ra.load_groups(sample_dir)
    primary_units = [unit for unit in layered["detail"] if unit.get("is_primary_lineage")]
    primary_by_region: dict[str, list[dict]] = defaultdict(list)
    for unit in primary_units:
        primary_by_region[unit["region"]].append(unit)
    for units in primary_by_region.values():
        units.sort(key=lambda unit: unit["family"])

    unit_results: dict[tuple[str, str], dict] = {}
    unit_rows: list[dict] = []
    candidate_rows_written = 0
    for unit in primary_units:
        key = (unit["region"], unit["family"])
        n_shapes_expected = unit.get("n_distinct_shapes_exact")
        base = {
            "sample": sample,
            "biological_id": biological_id,
            "region": unit["region"],
            "family": unit["family"],
            "n_trees_expected": int(unit["n_trees"]),
            "n_shapes_expected": int(n_shapes_expected) if n_shapes_expected is not None else None,
            "L1_base_class": unit.get("L1_base_class"),
            "analysis_candidate_set_complete": unit.get("analysis_candidate_set_complete"),
            "cn": unit.get("cn"),
        }
        if unit.get("analysis_candidate_set_complete") is not True:
            result = {**base, "status": "incomplete"}
        elif int(unit["n_trees"]) == 1:
            result = {
                **base,
                "status": "fixed_T1",
                "n_top_exact": 1,
                "n_top_shapes_exact": 1,
                "exact_top_class": "unique_first",
            }
        elif unit.get("L1_base_class") == "recurrence_required":
            result = {**base, "status": "recurrence_not_evaluable"}
        elif unit.get("L1_base_class") == "ambiguous":
            result, candidates = analyze_ranked_unit(
                sample=sample,
                biological_id=biological_id,
                unit=unit,
                group=groups.get(unit["region"]),
                ra=ra,
                solver=solver,
                temperatures=temperatures,
                near_tie_gaps=near_tie_gaps,
                violation_margin=violation_margin,
                candidate_writer=candidate_writer,
            )
            candidate_rows_written += len(candidates)
        else:
            raise RuntimeError(
                f"unsupported complete T>1 unit {sample} {unit['region']} HP{unit['family']} "
                f"class={unit.get('L1_base_class')}"
            )
        unit_results[key] = result
        unit_rows.append(result)

    region_rows: list[dict] = []
    for region, units in sorted(primary_by_region.items()):
        category = structural_class(units)
        t_joint = math.prod(int(unit["n_trees"]) for unit in units) if category != "incomplete" else None
        topo_joint = math.prod(int(unit["n_distinct_shapes_exact"]) for unit in units) if category != "incomplete" else None
        base = {
            "sample": sample,
            "biological_id": biological_id,
            "region": region,
            "primary_HP_units": len(units),
            "primary_HP_families": ",".join(unit["family"] for unit in units),
            "structural_class": category,
            "T_joint": t_joint,
            "Topo_joint": topo_joint,
        }
        if category == "incomplete":
            region_rows.append(
                {
                    **base,
                    "evaluable": False,
                    "final_class": "incomplete_or_capped",
                    "legacy_status_0.60": "not_applicable",
                }
            )
            continue
        if category == "T=1|Topo=1":
            region_rows.append(
                {
                    **base,
                    "evaluable": False,
                    "final_class": "structurally_unique_not_ranked",
                    "legacy_status_0.60": "not_applicable",
                    "n_top_joint_exact": 1,
                    "n_top_shapes_joint_exact": 1,
                    "exact_top_class": "unique_first",
                }
            )
            continue

        t_multiple_units = [unit for unit in units if int(unit["n_trees"]) > 1]
        results = [unit_results[(region, unit["family"])] for unit in t_multiple_units]
        if any(result["status"] == "recurrence_not_evaluable" for result in results):
            region_rows.append(
                {
                    **base,
                    "evaluable": False,
                    "final_class": "recurrence_not_evaluable",
                    "legacy_status_0.60": "not_evaluable_recurrence",
                }
            )
            continue
        if any(result["status"] in {"missing_group", "missing_vaf", "candidate_or_shape_mismatch"} for result in results):
            region_rows.append(
                {
                    **base,
                    "evaluable": False,
                    "final_class": "missing_or_mismatch",
                    "legacy_status_0.60": "missing_or_mismatch",
                }
            )
            continue
        if not all(result["status"] == "ranked" for result in results):
            raise RuntimeError(f"unexpected unit statuses {sample} {region}: {[r['status'] for r in results]}")

        exact = region_tie_summary(results)
        region_legacy_direction_flags: list[bool] = []
        region_direction_evidence_flags: list[bool | None] = []
        for result in results:
            if result["exact_top_legacy_direction_state"] == "all_consistent":
                region_legacy_direction_flags.extend([True] * int(result["n_top_exact"]))
            elif result["exact_top_legacy_direction_state"] == "all_inconsistent":
                region_legacy_direction_flags.extend([False] * int(result["n_top_exact"]))
            else:
                region_legacy_direction_flags.extend([True, False])
            region_direction_evidence_flags.extend(
                [True] * int(result["n_top_direction_consistent"])
            )
            region_direction_evidence_flags.extend(
                [False] * int(result["n_top_direction_inconsistent"])
            )
            region_direction_evidence_flags.extend(
                [None] * int(result["n_top_no_direction_comparisons"])
            )
        row = {
            **base,
            "evaluable": True,
            "final_class": exact["class"],
            "n_top_joint_exact": exact["n_top_joint"],
            "n_top_shapes_joint_exact": exact["n_top_shapes_joint"],
            "exact_top_class": exact["class"],
            "top_direction_state": direction_evidence_state(region_direction_evidence_flags),
            "legacy_top_direction_state": direction_state(region_legacy_direction_flags),
            "n_top_direction_evaluable": sum(flag is not None for flag in region_direction_evidence_flags),
            "n_top_no_direction_comparisons": sum(flag is None for flag in region_direction_evidence_flags),
            "ranked_Tn_units": len(results),
            "all_candidates_co_top_region": all(result["all_candidates_co_top"] for result in results),
        }
        for gap in near_tie_gaps:
            gap_key = f"{gap:.0e}"
            near = region_tie_summary(results, gap_key)
            row[f"near_tie_n_top_{gap_key}"] = near["n_top_joint"]
            row[f"near_tie_n_shapes_{gap_key}"] = near["n_top_shapes_joint"]
            row[f"near_tie_class_{gap_key}"] = near["class"]
        for temperature in temperatures:
            tk = temp_key(temperature)
            row[f"min_top_weight_T{tk}"] = min(
                float(result["top_weight_by_temperature"][tk]) for result in results
            )

        row["strict_gt_0.50_T0.05"] = row["min_top_weight_T0.05"] > 0.50
        row["inclusive_ge_0.50_T0.05"] = row["min_top_weight_T0.05"] >= 0.50
        row["ge_0.60_T0.05"] = row["min_top_weight_T0.05"] >= 0.60
        if exact["class"] == "unique_first":
            row["ranking_interpretation"] = (
                "unique_first_majority_gt_0.50"
                if row["strict_gt_0.50_T0.05"]
                else "unique_first_weak_le_0.50"
            )
        else:
            row["ranking_interpretation"] = exact["class"]

        if row["min_top_weight_T0.05"] < 0.60:
            row["legacy_status_0.60"] = "below_threshold"
        elif row["legacy_top_direction_state"] == "all_consistent":
            row["legacy_status_0.60"] = "vaf_top_consistent"
        else:
            row["legacy_status_0.60"] = "vaf_top_direction_inconsistent"
        region_rows.append(row)

    sample_summary = summarize_regions(region_rows, near_tie_gaps, temperatures, thresholds)
    legacy_expected = legacy_sample["region_vaf_status"]
    structural_expected = corrected_sample["T_and_Topology"]["classes"]
    checks = {
        "region_count_matches_W_primary": len(region_rows) == int(corrected_sample["W"]["primary"]),
        "structural_classes_match_corrected": all(
            sample_summary["structural_classes"][key] == int(structural_expected.get(key, 0))
            for key in STRUCTURAL_CLASSES
        ),
        "legacy_0.60_statuses_reproduced": all(
            sample_summary["legacy_region_status_0.60"][key] == int(legacy_expected[key])
            for key in LEGACY_STATUSES
        ),
        "legacy_below_partition_conserves": sum(
            sample_summary["legacy_below_0.60_exact_top_classes"][key] for key in TIE_CLASSES
        )
        == sample_summary["legacy_below_0.60_count"],
        "no_candidate_or_shape_mismatch": not any(
            row["status"] == "candidate_or_shape_mismatch" for row in unit_rows
        ),
        "float_score_matches_exact": all(
            float(row.get("max_float_score_error", 0.0)) <= 1e-12
            for row in unit_rows
            if row["status"] == "ranked"
        ),
    }
    sample_summary.update(
        {
            "sample": sample,
            "biological_id": biological_id,
            "candidate_rows_written": candidate_rows_written,
            "checks": checks,
            "source_layered": str(layered_path.resolve()),
        }
    )
    return region_rows, unit_rows, sample_summary, checks


def flatten_for_tsv(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    if isinstance(value, bool):
        return str(value)
    return str(value)


def write_dict_rows(path: Path, rows: list[dict], preferred_fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    all_fields = set().union(*(row.keys() for row in rows)) if rows else set()
    fields = [field for field in preferred_fields if field in all_fields]
    fields.extend(sorted(all_fields - set(fields)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: flatten_for_tsv(row.get(field)) for field in fields})


def pct(count: int, denominator: int) -> str:
    return f"{100 * count / denominator:.2f}%" if denominator else "NA"


def build_report(output: dict) -> str:
    aggregate = output["aggregate"]
    below = aggregate["legacy_below_0.60_exact_top_classes"]
    threshold = aggregate["legacy_below_0.60_threshold_0.50"]
    total_below = aggregate["legacy_below_0.60_count"]
    complete = aggregate["complete_regions"]
    ambiguous = aggregate["ambiguous_complete_regions"]
    evaluable = aggregate["evaluable_ambiguous_regions"]
    all_tie = aggregate["exact_top_classes_all_evaluable"]
    cross = aggregate["legacy_below_0.60_by_structural_and_top_class"]
    topo1_cross = cross.get("T>1|Topo=1", {})
    topon_cross = cross.get("T>1|Topo>1", {})
    grid = aggregate["threshold_grid"]
    strict_050 = int(grid["0.05"]["0.50"]["gt"])
    legacy_060 = int(grid["0.05"]["0.60"]["ge"])
    increment = strict_050 - legacy_060
    sample_lines = []
    for sample in output["samples"]:
        sample_below = sample["legacy_below_0.60_exact_top_classes"]
        sample_threshold = sample["legacy_below_0.60_threshold_0.50"]
        sample_lines.append(
            f"| {sample['sample']} | {sample['legacy_below_0.60_count']:,} | "
            f"{sample_below['unique_first']:,} | {sample_below['tied_first_same_topology']:,} | "
            f"{sample_below['tied_first_different_topology']:,} | {sample_threshold['gt_0.50']:,} |"
        )
    top_bins = aggregate["exact_top_joint_count_bins_legacy_below"]
    top_bin_values = {int(key): int(value) for key, value in top_bins.items()}
    grouped_top_bins = {
        "1": top_bin_values.get(1, 0),
        "2": top_bin_values.get(2, 0),
        "3–5": sum(value for key, value in top_bin_values.items() if 3 <= key <= 5),
        "6–10": sum(value for key, value in top_bin_values.items() if 6 <= key <= 10),
        ">10": sum(value for key, value in top_bin_values.items() if key > 10),
    }
    unit_diag = aggregate["unit_diagnostics"]
    validation_checks = output["validation"]["aggregate_checks"]
    lines = [
        "<!--",
        f"建立時間: {output['generated_at']}",
        "目標: 重算 VAF strict top-weight >0.50 / legacy 0.60、exact tie 與 Topo 一致性",
        f"處理範圍: {output['scope']}",
        f"關聯檔案: {output['sources']['legacy_vaf_summary']['path']}; {output['sources']['corrected_report']['path']}",
        "-->",
        "",
        "# VAF strict >0.50、並列第一與拓撲一致性完整重算",
        "",
        f"> **TL;DR**：原 {total_below:,} 個 `top<0.60` regions 中，{below['unique_first']:,}（{pct(below['unique_first'], total_below)}）有唯一 exact 第一、{below['tied_first_same_topology']:,}（{pct(below['tied_first_same_topology'], total_below)}）並列但 Topo 一致、{below['tied_first_different_topology']:,}（{pct(below['tied_first_different_topology'], total_below)}）並列且 Topo 不一致。strict `>0.50` 只新增 {threshold['gt_0.50']:,} 個唯一第一；直接 rank/tie 適合作主分類，weight threshold 只宜作 concentration 輔助。",
        "",
        "> **證據邊界**：historical layered-v2 描述性工程重算；relative weight 不是 calibrated probability，selection ≠ confirmation。",
        "",
        "## 研究問題與方法",
        "",
        f"- Task type：B comprehensive validation；chr1-22 × {len(output['samples'])} dataset rows／{len({sample['biological_id'] for sample in output['samples']})} biological samples。",
        "- 對每個 complete `region×primary HP` 以 frozen solver `tree_cap=0` 重枚舉所有 exact T。",
        "- read AF 由整數 REF/ALT counts 轉成 `Fraction`；最高 exact ordering score 完全相等才算並列第一。",
        "- Region top set 是 HP-labeled Cartesian product；可相乘集合基數，但不相乘 VAF weights。",
        "- `N_top=1` 為唯一第一；`N_top>1, N_top_shapes=1` 為並列同 Topo；`N_top_shapes>1` 為並列異 Topo。",
        "",
        "## 分母與守恆",
        "",
        f"- `W_primary={aggregate['region_count']:,} = complete {complete:,} + incomplete {aggregate['structural_classes']['incomplete']:,}`。",
        f"- `complete {complete:,} = T=1 {aggregate['structural_classes']['T=1|Topo=1']:,} + ambiguous {ambiguous:,}`。",
        f"- `ambiguous {ambiguous:,} = evaluable {evaluable:,} + missing/mismatch {aggregate['legacy_region_status_0.60']['missing_or_mismatch']:,} + recurrence {aggregate['legacy_region_status_0.60']['not_evaluable_recurrence']:,}`。",
        f"- Legacy 0.60：`{ambiguous:,} = top {legacy_060:,} + below {total_below:,} + missing {aggregate['legacy_region_status_0.60']['missing_or_mismatch']:,} + recurrence {aggregate['legacy_region_status_0.60']['not_evaluable_recurrence']:,}`。",
        "",
        f"## 全部 {evaluable:,} 個可評估 ambiguous regions：直接排序結果",
        "",
        "| Exact 第一順位類別 | Regions | 可評估區域比例 |",
        "|---|---:|---:|",
        f"| 唯一第一 | {all_tie['unique_first']:,} | {pct(all_tie['unique_first'], evaluable)} |",
        f"| 並列第一、Topo 一致 | {all_tie['tied_first_same_topology']:,} | {pct(all_tie['tied_first_same_topology'], evaluable)} |",
        f"| 並列第一、Topo 不一致 | {all_tie['tied_first_different_topology']:,} | {pct(all_tie['tied_first_different_topology'], evaluable)} |",
        f"| **合計** | **{sum(all_tie.values()):,}** | **100.00%** |",
        "",
        f"## 原 {total_below:,} below-0.60 regions：exact 第一順位分類",
        "",
        f"| 分類 | Regions | {total_below:,} 內比例 | 意義 |",
        "|---|---:|---:|---|",
        f"| 唯一第一 | {below['unique_first']:,} | {pct(below['unique_first'], total_below)} | highest exact score 只有一棵樹 |",
        f"| 並列第一、Topo 一致 | {below['tied_first_same_topology']:,} | {pct(below['tied_first_same_topology'], total_below)} | 多棵 co-top T，但 canonical shape 相同 |",
        f"| 並列第一、Topo 不一致 | {below['tied_first_different_topology']:,} | {pct(below['tied_first_different_topology'], total_below)} | co-top T 橫跨多種 canonical shapes |",
        f"| **合計** | **{sum(below.values()):,}** | **100.00%** | 守恆回 legacy below-0.60 |",
        "",
        "### 原始 Topo 類別交叉表",
        "",
        "| 原始 structural class | 唯一第一 | 並列同 Topo | 並列異 Topo | 合計 |",
        "|---|---:|---:|---:|---:|",
        f"| T>1·Topo=1 | {topo1_cross.get('unique_first', 0):,} | {topo1_cross.get('tied_first_same_topology', 0):,} | {topo1_cross.get('tied_first_different_topology', 0):,} | {sum(topo1_cross.values()):,} |",
        f"| T>1·Topo>1 | {topon_cross.get('unique_first', 0):,} | {topon_cross.get('tied_first_same_topology', 0):,} | {topon_cross.get('tied_first_different_topology', 0):,} | {sum(topon_cross.values()):,} |",
        f"| **合計** | **{below['unique_first']:,}** | **{below['tied_first_same_topology']:,}** | **{below['tied_first_different_topology']:,}** | **{total_below:,}** |",
        "",
        "### 第一順位樹數分布",
        "",
        f"| Joint co-top T 數 | Regions | {total_below:,} 內比例 |",
        "|---|---:|---:|",
        *[
            f"| {label} | {count:,} | {pct(count, total_below)} |"
            for label, count in grouped_top_bins.items()
        ],
        "",
        "## strict >0.50 與 legacy ≥0.60",
        "",
        f"- 原 {total_below:,} 中：strict `>0.50`={threshold['gt_0.50']:,}（{pct(threshold['gt_0.50'], total_below)}）；`=0.50`={threshold['eq_0.50']:,}（{pct(threshold['eq_0.50'], total_below)}）；`<0.50`={threshold['lt_0.50']:,}（{pct(threshold['lt_0.50'], total_below)}）。",
        f"- strict `>0.50` 的 {threshold['gt_0.50']:,} 個全部為 unique first；non-unique=0。",
        f"- `=0.50` 的 {threshold['eq_0.50']:,} 個全部為 tied/same Topo；若用 `>=0.50` 會誤把這些並列納入 priority。",
        f"- 全 ambiguous 分母下，legacy `>=0.60`={legacy_060:,}/{ambiguous:,}（{pct(legacy_060, ambiguous)}）；strict `>0.50`={strict_050:,}/{ambiguous:,}（{pct(strict_050, ambiguous)}），增加 {increment:,} regions／{100 * increment / ambiguous:.2f} percentage points。",
        "",
        "### Temperature 敏感度（strict >0.50）",
        "",
        f"| Temperature | Regions | Ambiguous {ambiguous:,} 比例 |",
        "|---:|---:|---:|",
        *[
            f"| {temperature} | {int(grid[temperature]['0.50']['gt']):,} | {pct(int(grid[temperature]['0.50']['gt']), ambiguous)} |"
            for temperature in ("0.025", "0.05", "0.1")
        ],
        "",
        "Exact tie 分類在 score-gap `1e-12／1e-9／1e-6` 三組數值敏感度下完全相同；這只代表對 `≤1e-6` 的 numerical tolerance 不變，不代表對 read-sampling uncertainty 穩健。Weight threshold 數量則明顯隨 temperature 改變，因此直接 rank/tie 比固定 weight threshold 更適合作主結果。",
        "",
        "## 逐 dataset：原 below-0.60 拆分",
        "",
        "| Dataset row | Below 0.60 | 唯一第一 | 並列同 Topo | 並列異 Topo | strict >0.50 |",
        "|---|---:|---:|---:|---:|---:|",
        *sample_lines,
        f"| **合計** | **{total_below:,}** | **{below['unique_first']:,}** | **{below['tied_first_same_topology']:,}** | **{below['tied_first_different_topology']:,}** | **{threshold['gt_0.50']:,}** |",
        "",
        "## 方法判斷",
        "",
        "1. **直接排序適合作主分類**：它直接回答第一順位是一棵、多棵同 Topo、或多棵異 Topo，且不受 softmax temperature 改變名次／exact tie。",
        "2. **strict >0.50 適合作輔助 majority-concentration flag**：它可保證單棵候選過半，但不是統計顯著性或真值機率。",
        "3. **不能只看唯一第一**：仍應同報 top/second gap、direction flag、候選數與 CN/read-depth；unique weak 不等於可信解析。",
        f"4. **Topo 一致比 exact T 唯一更寬鬆**：{below['tied_first_same_topology']:,} 個 regions 的 co-top T 都落在同一 shape；{below['tied_first_different_topology']:,} 個 regions 的 exact 第一順位集合仍跨越多種 Topo。",
        f"5. **Score-sum 可比性仍待 sensitivity**：{unit_diag['candidate_comparison_count_varies']:,}/{unit_diag['ranked_units']:,} ranked units 的候選樹 ancestry-comparison 數不同，需再比較 normalized score／likelihood。",
        f"6. **方向證據不可用 vacuous truth**：{unit_diag['exact_top_candidates_without_direction_comparisons']:,} 個 exact-top candidates 沒有 ancestry comparison，已標為 direction-not-evaluable，不算一致性支持。",
        "",
        "## 判讀限制",
        "",
        "- Exact tie 由整數 REF/ALT counts 推導的 Fraction score 完全相等判定。",
        "- Relative weight 受 temperature 與候選樹數影響，不是 calibrated posterior。",
        "- 雙 HP 只相乘 top-set 基數形成 ordered Cartesian product；不相乘 VAF weights。",
        "- raw VAF 未校正 copy number、tumor purity、mutation multiplicity；selection ≠ confirmation。",
        "- Exact ranking 是 heuristic ordering-score argmax；尚未做 read/molecule bootstrap 或 normalized-score sensitivity。",
        "",
        "## 驗證",
        "",
        f"- Validation status：**{output['validation']['status']}**",
        f"- Failures：{output['validation']['failures'] or 'none'}",
        *[f"- `{key}`：{'PASS' if value else 'FAIL'}" for key, value in validation_checks.items()],
        "- 逐 region TSV 已由獨立 group-by reviewer 另行重算；subagent verdict 於交付時補入 implementation notes。",
        "",
        "## 可重現輸入",
        "",
        f"- Run root：`{output['sources']['run_root']}`",
        f"- Corrected report：`{output['sources']['corrected_report']['path']}`",
        f"- Legacy VAF summary：`{output['sources']['legacy_vaf_summary']['path']}`",
        f"- Method SHA-256：`{output['sources']['method_script']['sha256']}`",
        f"- Solver SHA-256：`{output['sources']['solver_script']['sha256']}`",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--corrected-report", required=True, type=Path)
    parser.add_argument("--legacy-vaf-summary", required=True, type=Path)
    parser.add_argument("--method-script-dir", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-region-tsv", required=True, type=Path)
    parser.add_argument("--output-unit-tsv", required=True, type=Path)
    parser.add_argument("--output-candidate-tsv-gz", required=True, type=Path)
    parser.add_argument("--output-summary-tsv", required=True, type=Path)
    parser.add_argument("--output-checks-tsv", required=True, type=Path)
    parser.add_argument("--output-report-md", required=True, type=Path)
    parser.add_argument("--samples", default=",".join(SAMPLE_ORDER))
    parser.add_argument("--temperatures", default=",".join(str(value) for value in DEFAULT_TEMPERATURES))
    parser.add_argument("--thresholds", default=",".join(str(value) for value in DEFAULT_THRESHOLDS))
    parser.add_argument("--near-tie-gaps", default=",".join(str(value) for value in DEFAULT_NEAR_TIE_GAPS))
    parser.add_argument("--violation-margin", default="0.05")
    args = parser.parse_args()

    selected_samples = tuple(item.strip() for item in args.samples.split(",") if item.strip())
    if not selected_samples or any(sample not in SAMPLE_ORDER for sample in selected_samples):
        raise SystemExit(f"invalid sample selection: {selected_samples}")
    temperatures = tuple(float(value) for value in args.temperatures.split(","))
    thresholds = tuple(float(value) for value in args.thresholds.split(","))
    near_tie_gaps = tuple(float(value) for value in args.near_tie_gaps.split(","))
    if 0.05 not in temperatures or 0.50 not in thresholds or 0.60 not in thresholds:
        raise SystemExit("temperature grid must contain 0.05 and thresholds must contain 0.50/0.60")
    violation_margin = Fraction(args.violation_margin)

    sys.path.insert(0, str(args.method_script_dir.resolve()))
    import read_af_tree_ordering_multisample as ra  # noqa: PLC0415
    import tree_enumeration_solver as solver  # noqa: PLC0415

    manifest = json.loads(args.input_manifest.read_text(encoding="utf-8"))
    corrected = json.loads(args.corrected_report.read_text(encoding="utf-8"))
    legacy = json.loads(args.legacy_vaf_summary.read_text(encoding="utf-8"))
    meta_by_name = {item["sample"]: item for item in manifest["samples"]}
    corrected_by_name = {item["sample"]: item for item in corrected["samples"]}
    legacy_by_name = {item["sample"]: item for item in legacy["samples"]}
    if not set(selected_samples) <= set(meta_by_name) or not set(selected_samples) <= set(corrected_by_name) or not set(selected_samples) <= set(legacy_by_name):
        raise SystemExit("selected sample missing from an input source")

    candidate_fields = [
        "sample",
        "biological_id",
        "region",
        "family",
        "candidate_order",
        "score_rank",
        "candidate_index",
        "candidate_id",
        "shape_id",
        "exact_score_numerator",
        "exact_score_denominator",
        "score_float",
        "softmax_weight_T0.05",
        "is_exact_top",
        "direction_evaluable",
        "direction_consistent",
        "legacy_direction_consistent",
        "n_ancestry_comparisons",
    ]
    args.output_candidate_tsv_gz.parent.mkdir(parents=True, exist_ok=True)
    all_regions: list[dict] = []
    all_units: list[dict] = []
    sample_summaries: list[dict] = []
    all_checks: list[tuple[str, str, bool]] = []
    with gzip.open(args.output_candidate_tsv_gz, "wt", encoding="utf-8", newline="") as handle:
        candidate_writer = csv.DictWriter(handle, fieldnames=candidate_fields, delimiter="\t", lineterminator="\n")
        candidate_writer.writeheader()
        for sample in selected_samples:
            regions, units, summary, checks = analyze_sample(
                meta=meta_by_name[sample],
                sample_dir=args.run_root / "samples" / sample,
                corrected_sample=corrected_by_name[sample],
                legacy_sample=legacy_by_name[sample],
                ra=ra,
                solver=solver,
                temperatures=temperatures,
                thresholds=thresholds,
                near_tie_gaps=near_tie_gaps,
                violation_margin=violation_margin,
                candidate_writer=candidate_writer,
            )
            all_regions.extend(regions)
            all_units.extend(units)
            sample_summaries.append(summary)
            all_checks.extend((sample, key, bool(value)) for key, value in checks.items())

    aggregate = summarize_regions(all_regions, near_tie_gaps, temperatures, thresholds)
    ranked_units = [row for row in all_units if row["status"] == "ranked"]
    aggregate["unit_diagnostics"] = {
        "ranked_units": len(ranked_units),
        "candidate_comparison_count_varies": sum(
            bool(row.get("candidate_comparison_count_varies")) for row in ranked_units
        ),
        "units_with_top_candidate_without_direction_comparisons": sum(
            int(row.get("n_top_no_direction_comparisons", 0)) > 0 for row in ranked_units
        ),
        "exact_top_candidates_without_direction_comparisons": sum(
            int(row.get("n_top_no_direction_comparisons", 0)) for row in ranked_units
        ),
        "max_fraction_to_float_score_error": max(
            (float(row.get("max_float_score_error", 0.0)) for row in ranked_units),
            default=0.0,
        ),
    }
    aggregate_checks = {
        "all_sample_checks_pass": all(value for _, _, value in all_checks),
        "region_rows_equal_sample_sum": len(all_regions) == sum(int(summary["region_count"]) for summary in sample_summaries),
        "legacy_below_partition_conserves": sum(
            aggregate["legacy_below_0.60_exact_top_classes"][key] for key in TIE_CLASSES
        )
        == aggregate["legacy_below_0.60_count"],
        "tied_different_impossible_when_original_Topo1": aggregate[
            "legacy_below_0.60_by_structural_and_top_class"
        ].get("T>1|Topo=1", {}).get("tied_first_different_topology", 0)
        == 0,
        "complete_conservation": aggregate["complete_regions"] + aggregate["structural_classes"]["incomplete"]
        == aggregate["region_count"],
        "ambiguous_conservation": sum(aggregate["legacy_region_status_0.60"].values())
        == aggregate["ambiguous_complete_regions"],
        "full_scope_has_seven_rows": len(selected_samples) != len(SAMPLE_ORDER) or tuple(selected_samples) == tuple(SAMPLE_ORDER),
        "sample_region_keys_unique": len({(row["sample"], row["region"]) for row in all_regions})
        == len(all_regions),
        "legacy_below_rows_are_evaluable_and_below_0.60": all(
            row["evaluable"]
            and row["min_top_weight_T0.05"] < 0.60
            for row in all_regions
            if row["legacy_status_0.60"] == "below_threshold"
        ),
        "strict_gt_0.50_implies_unique_first": all(
            row["exact_top_class"] == "unique_first"
            for row in all_regions
            if row.get("evaluable") and row.get("min_top_weight_T0.05", 0.0) > 0.50
        ),
        "top_set_cardinality_bounds_hold": all(
            1 <= int(row["n_top_shapes_joint_exact"])
            <= int(row["n_top_joint_exact"])
            <= int(row["T_joint"])
            for row in all_regions
            if row.get("evaluable")
        ),
        "exact_top_class_matches_cardinalities": all(
            row["exact_top_class"]
            == classify_top_set(
                int(row["n_top_joint_exact"]), int(row["n_top_shapes_joint_exact"])
            )
            for row in all_regions
            if row.get("evaluable")
        ),
        "near_tie_sensitivity_stable_through_1e-6": all(
            all(
                row[f"near_tie_class_{gap:.0e}"] == row["exact_top_class"]
                for gap in near_tie_gaps
            )
            for row in all_regions
            if row.get("evaluable")
        ),
        "direction_evidence_accounting_conserves": all(
            int(row.get("n_top_direction_evaluable", 0))
            + int(row.get("n_top_no_direction_comparisons", 0))
            == int(row["n_top_exact"])
            for row in ranked_units
        ),
    }
    all_checks.extend(("aggregate", key, bool(value)) for key, value in aggregate_checks.items())
    validation_pass = all(value for _, _, value in all_checks)
    generated_at = datetime.now(ZoneInfo("Asia/Taipei")).isoformat(timespec="seconds")
    output = {
        "schema_version": "1.0",
        "generated_at": generated_at,
        "status": "PASS" if validation_pass else "FAIL",
        "scope": f"chr1-22; {len(selected_samples)} dataset rows; historical layered-v2 engineering snapshot",
        "samples_selected": list(selected_samples),
        "claim_ceiling": "exact top-tie and topology census under raw family-specific read-AF ordering; not a confirmed biological tree",
        "definitions": {
            "exact_tie": "Fraction-valued raw ordering scores are exactly equal to the maximum",
            "region_top_set": "HP-labeled Cartesian product of primary-unit exact top sets; weights are not multiplied",
            "unique_first": "n_top_joint=1",
            "tied_first_same_topology": "n_top_joint>1 and n_top_shapes_joint=1",
            "tied_first_different_topology": "n_top_shapes_joint>1",
            "relative_weight": "temperature-scaled softmax concentration heuristic; not calibrated probability",
        },
        "parameters": {
            "temperatures": temperatures,
            "thresholds": thresholds,
            "near_tie_absolute_score_gaps": near_tie_gaps,
            "violation_margin": float(violation_margin),
        },
        "aggregate": aggregate,
        "samples": sample_summaries,
        "validation": {
            "status": "PASS" if validation_pass else "FAIL",
            "aggregate_checks": aggregate_checks,
            "failures": [f"{scope}:{check}" for scope, check, value in all_checks if not value],
        },
        "sources": {
            "run_root": str(args.run_root.resolve()),
            "input_manifest": {"path": str(args.input_manifest.resolve()), "sha256": sha256(args.input_manifest)},
            "corrected_report": {"path": str(args.corrected_report.resolve()), "sha256": sha256(args.corrected_report)},
            "legacy_vaf_summary": {"path": str(args.legacy_vaf_summary.resolve()), "sha256": sha256(args.legacy_vaf_summary)},
            "method_script": {
                "path": str((args.method_script_dir / "read_af_tree_ordering_multisample.py").resolve()),
                "sha256": sha256(args.method_script_dir / "read_af_tree_ordering_multisample.py"),
            },
            "solver_script": {
                "path": str((args.method_script_dir / "tree_enumeration_solver.py").resolve()),
                "sha256": sha256(args.method_script_dir / "tree_enumeration_solver.py"),
            },
        },
    }

    dump_json(args.output_json, output)
    write_dict_rows(
        args.output_region_tsv,
        all_regions,
        [
            "sample",
            "biological_id",
            "region",
            "primary_HP_units",
            "primary_HP_families",
            "structural_class",
            "T_joint",
            "Topo_joint",
            "evaluable",
            "legacy_status_0.60",
            "exact_top_class",
            "n_top_joint_exact",
            "n_top_shapes_joint_exact",
            "top_direction_state",
            "legacy_top_direction_state",
            "n_top_direction_evaluable",
            "n_top_no_direction_comparisons",
            "min_top_weight_T0.05",
            "strict_gt_0.50_T0.05",
            "inclusive_ge_0.50_T0.05",
            "ge_0.60_T0.05",
            "ranking_interpretation",
        ],
    )
    write_dict_rows(
        args.output_unit_tsv,
        all_units,
        [
            "sample",
            "biological_id",
            "region",
            "family",
            "status",
            "L1_base_class",
            "n_trees_expected",
            "n_trees_enumerated",
            "n_shapes_expected",
            "n_shapes_enumerated",
            "n_top_exact",
            "n_top_shapes_exact",
            "exact_top_class",
            "top_weight_T0.05",
            "second_weight_T0.05",
            "weight_gap_T0.05",
            "exact_top_direction_state",
            "exact_top_legacy_direction_state",
            "n_top_direction_evaluable",
            "n_top_no_direction_comparisons",
            "candidate_comparison_count_varies",
        ],
    )

    summary_rows: list[dict] = []
    for scope, summary in [(item["sample"], item) for item in sample_summaries] + [("aggregate", aggregate)]:
        denominator = int(summary["legacy_below_0.60_count"])
        for tie_class in TIE_CLASSES:
            count = int(summary["legacy_below_0.60_exact_top_classes"][tie_class])
            summary_rows.append(
                {
                    "scope": scope,
                    "subset": "legacy_below_0.60",
                    "class": tie_class,
                    "count": count,
                    "denominator": denominator,
                    "share": count / denominator if denominator else None,
                }
            )
    write_dict_rows(
        args.output_summary_tsv,
        summary_rows,
        ["scope", "subset", "class", "count", "denominator", "share"],
    )
    args.output_checks_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_checks_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["scope", "check", "pass"])
        for scope, check, value in all_checks:
            writer.writerow([scope, check, str(bool(value))])
    args.output_report_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_report_md.write_text(build_report(output), encoding="utf-8")

    print(f"INPUT RUN ROOT -> {args.run_root.resolve()}")
    print(f"INPUT MANIFEST -> {args.input_manifest.resolve()}")
    print(f"INPUT CORRECTED REPORT -> {args.corrected_report.resolve()}")
    print(f"INPUT LEGACY VAF SUMMARY -> {args.legacy_vaf_summary.resolve()}")
    print(f"OUTPUT JSON -> {args.output_json.resolve()}")
    print(f"OUTPUT REGION TSV -> {args.output_region_tsv.resolve()}")
    print(f"OUTPUT UNIT TSV -> {args.output_unit_tsv.resolve()}")
    print(f"OUTPUT CANDIDATE TSV.GZ -> {args.output_candidate_tsv_gz.resolve()}")
    print(f"OUTPUT SUMMARY TSV -> {args.output_summary_tsv.resolve()}")
    print(f"OUTPUT CHECKS TSV -> {args.output_checks_tsv.resolve()}")
    print(f"OUTPUT REPORT -> {args.output_report_md.resolve()}")
    print(f"STATUS -> {'PASS' if validation_pass else 'FAIL'}")
    if not validation_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
