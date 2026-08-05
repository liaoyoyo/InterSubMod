#!/usr/bin/env python3
"""Compare HCC1395 VAF-prioritized candidate forests across two datasets.

The analysis is intentionally descriptive.  Per-site read-derived AF is used
by the frozen upstream census as an ancestor >= descendant ordering heuristic.
The resulting exact-score maximizers and softmax weights are not calibrated
posterior probabilities and are not a biological tree truth set.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable


SAMPLE_A = "HCC1395"
SAMPLE_B = "HCC1395_DORADO"
SAMPLES = (SAMPLE_A, SAMPLE_B)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> Iterable[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def dump_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def truth(value: str | bool | None) -> bool:
    return value is True or value == "True"


def canonical_shape(edges: list[list[str]]) -> str:
    children: dict[str, list[str]] = collections.defaultdict(list)
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


def shape_id(edges: list[list[str]]) -> str:
    return "TS-" + hashlib.sha1(canonical_shape(edges).encode()).hexdigest()[:10]


def tree_id(edges: list[list[str]]) -> str:
    payload = sorted([list(edge) for edge in edges])
    encoded = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode()
    return "TR-" + hashlib.sha1(encoded).hexdigest()[:12]


def parse_labeled_set_signature(value: str) -> tuple[tuple[str, tuple[str, ...]], ...]:
    result = []
    for component in value.split("|") if value else []:
        label, members = component.split("=", 1)
        family = label.removeprefix("HP")
        values = tuple(sorted(item for item in members.strip("{}").split(",") if item))
        result.append((family, values))
    return tuple(sorted(result, key=lambda item: int(item[0])))


def ordered_signature(
    forest: tuple[tuple[str, tuple[int, ...], tuple[str, ...], tuple[str, ...], str], ...],
    field: str,
) -> str:
    index = {"candidate": 2, "shape": 3}[field]
    return "|".join(
        f"HP{family}@{','.join(map(str, positions))}=" + "{" + ",".join(component[index]) + "}"
        for family, positions, *_ in forest
        for component in [(family, positions, forest_item(forest, family, index), (), "")]
    )


def forest_item(
    forest: tuple[tuple[str, tuple[int, ...], tuple[str, ...], tuple[str, ...], str], ...],
    family: str,
    index: int,
) -> tuple[str, ...]:
    for component in forest:
        if component[0] == family:
            return component[index]
    raise KeyError(family)


def candidate_key(
    forest: tuple[tuple[str, tuple[int, ...], tuple[str, ...], tuple[str, ...], str], ...],
    swap_tolerant: bool,
) -> tuple:
    components = tuple((positions, candidates) for _, positions, candidates, _, _ in forest)
    if swap_tolerant:
        return tuple(sorted(components))
    return tuple((family, positions, candidates) for family, positions, candidates, _, _ in forest)


def shape_key(
    forest: tuple[tuple[str, tuple[int, ...], tuple[str, ...], tuple[str, ...], str], ...],
    swap_tolerant: bool,
) -> tuple:
    components = tuple(shapes for _, _, _, shapes, _ in forest)
    if swap_tolerant:
        return tuple(sorted(components))
    return tuple((family, shapes) for family, _, _, shapes, _ in forest)


def unlabeled_signature(value: str) -> tuple[str, ...]:
    return tuple(sorted(component.split("=", 1)[1] for component in value.split("|") if "=" in component))


def metric_row(
    family: str,
    population: str,
    comparison: str,
    numerator: int,
    denominator: int,
    interpretation: str,
) -> dict:
    return {
        "metric_family": family,
        "population": population,
        "comparison": comparison,
        "numerator": numerator,
        "denominator": denominator,
        "percent": f"{100 * numerator / denominator:.6f}" if denominator else "",
        "interpretation": interpretation,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vaf-units", type=Path, required=True)
    parser.add_argument("--vaf-candidates", type=Path, required=True)
    parser.add_argument("--vaf-regions", type=Path, required=True)
    parser.add_argument("--pair-matches", type=Path, required=True)
    parser.add_argument("--pair-analysis", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    pair_analysis = json.loads(args.pair_analysis.read_text(encoding="utf-8"))
    layered_paths = {
        SAMPLE_A: Path(pair_analysis["sources"]["layered_hcc1395"]["path"]),
        SAMPLE_B: Path(pair_analysis["sources"]["layered_hcc1395_dorado"]["path"]),
    }

    units: dict[tuple[str, str, str], dict[str, str]] = {}
    units_by_region: dict[tuple[str, str], list[dict[str, str]]] = collections.defaultdict(list)
    for row in read_tsv(args.vaf_units):
        if row["sample"] not in SAMPLES:
            continue
        key = (row["sample"], row["region"], row["family"])
        if key in units:
            raise RuntimeError(f"duplicate VAF unit: {key}")
        units[key] = row
        units_by_region[(row["sample"], row["region"])].append(row)
    for rows in units_by_region.values():
        rows.sort(key=lambda row: int(row["family"]))

    region_status: dict[tuple[str, str], dict[str, str]] = {}
    for row in read_tsv(args.vaf_regions):
        if row["sample"] not in SAMPLES:
            continue
        key = (row["sample"], row["region"])
        if key in region_status:
            raise RuntimeError(f"duplicate VAF region: {key}")
        region_status[key] = row

    top_candidates_from_detail: dict[tuple[str, str, str], set[str]] = collections.defaultdict(set)
    top_shapes_from_detail: dict[tuple[str, str, str], set[str]] = collections.defaultdict(set)
    candidate_rows = 0
    with gzip.open(args.vaf_candidates, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample"] not in SAMPLES:
                continue
            candidate_rows += 1
            is_top = truth(row["is_exact_top"])
            if is_top != (int(row["score_rank"]) == 1):
                raise RuntimeError(f"top/rank mismatch: {row['sample']} {row['region']} HP{row['family']}")
            if is_top:
                key = (row["sample"], row["region"], row["family"])
                top_candidates_from_detail[key].add(row["candidate_id"])
                top_shapes_from_detail[key].add(row["shape_id"])

    positions_by_region: dict[tuple[str, str], tuple[int, ...]] = {}
    mlhp_sources: list[dict] = []
    for sample in SAMPLES:
        sample_dir = args.run_root / "samples" / sample
        for path in sorted(sample_dir.glob("mlhp_part_*.json")):
            mlhp_sources.append({"sample": sample, "path": str(path.resolve()), "sha256": sha256(path)})
            document = json.loads(path.read_text(encoding="utf-8"))
            for group in document.get("groups", []):
                region = f"{group['chrom']}:{group['start']}-{group['end']}"
                key = (sample, region)
                value = tuple(int(position) for position in group.get("positions", []))
                if key in positions_by_region and positions_by_region[key] != value:
                    raise RuntimeError(f"duplicate group with discordant positions: {key}")
                positions_by_region[key] = value

    fixed: dict[tuple[str, str, str], tuple[str, str]] = {}
    layered_digests: dict[tuple[str, str, str], str] = {}
    for sample, path in layered_paths.items():
        document = json.loads(path.read_text(encoding="utf-8"))
        for unit in document["detail"]:
            if not unit.get("is_primary_lineage") or unit.get("analysis_candidate_set_complete") is not True:
                continue
            key = (sample, unit["region"], str(unit["family"]))
            layered_digests[key] = str(unit["analysis_tree_digest_sha256"])
            if int(unit["n_trees"]) != 1:
                continue
            if not unit.get("display_trees_complete") or len(unit.get("trees", [])) != 1:
                raise RuntimeError(f"fixed T1 tree was not materialized: {key}")
            edges = unit["trees"][0]["edges"]
            fixed[key] = (tree_id(edges), shape_id(edges))

    # Exact-coordinate rows are the population; complete-both is the primary denominator.
    exact_matches: list[dict[str, str]] = []
    for row in read_tsv(args.pair_matches):
        if row["scenario"] == "exact_coordinate":
            exact_matches.append(row)
    exact_complete = [row for row in exact_matches if truth(row["complete_a"]) and truth(row["complete_b"])]

    # Validate detailed exact-top rows against the frozen unit summary.
    ranked_units = 0
    for key, row in units.items():
        if row["status"] != "ranked":
            continue
        ranked_units += 1
        expected_candidates = set(json.loads(row["exact_top_candidate_ids"]))
        expected_shapes = set(json.loads(row["exact_top_shape_ids"]))
        if expected_candidates != top_candidates_from_detail.get(key, set()):
            raise RuntimeError(f"candidate exact-top mismatch: {key}")
        if expected_shapes != top_shapes_from_detail.get(key, set()):
            raise RuntimeError(f"shape exact-top mismatch: {key}")

    def build_exact_forest(sample: str, region: str):
        components = []
        rows = units_by_region.get((sample, region), [])
        for unit in rows:
            family = unit["family"]
            status = unit["status"]
            positions = positions_by_region.get((sample, region), ())
            if not positions:
                return None
            if status == "ranked":
                candidate_ids = tuple(sorted(json.loads(unit["exact_top_candidate_ids"])))
                shape_ids = tuple(sorted(json.loads(unit["exact_top_shape_ids"])))
                source = "vaf_exact_top"
            elif status == "fixed_T1":
                candidate_id, fixed_shape_id = fixed[(sample, region, family)]
                candidate_ids = (candidate_id,)
                shape_ids = (fixed_shape_id,)
                source = "structural_fixed_T1"
            else:
                return None
            components.append((family, positions, candidate_ids, shape_ids, source))
        return tuple(components) if components else None

    def build_selected_shape(sample: str, region: str, original_signature: str):
        status = region_status[(sample, region)]
        structural = status["structural_class"]
        if structural in {"T=1|Topo=1", "T>1|Topo=1"}:
            value = parse_labeled_set_signature(original_signature)
            if not value or any(len(shapes) != 1 for _, shapes in value):
                raise RuntimeError(f"structural Topo=1 did not have one shape per component: {(sample, region)}")
            return value, "structural_topo1"
        if (
            structural == "T>1|Topo>1"
            and truth(status["evaluable"])
            and int(status["n_top_shapes_joint_exact"]) == 1
        ):
            forest = build_exact_forest(sample, region)
            if forest is None or any(len(component[3]) != 1 for component in forest):
                raise RuntimeError(f"VAF-selected single shape could not be materialized: {(sample, region)}")
            return tuple((component[0], component[3]) for component in forest), "vaf_exact_top_shape"
        return None, "not_selectable"

    rows_out: list[dict] = []
    counters: collections.Counter[str] = collections.Counter()
    provenance_pairs: collections.Counter[tuple[str, str]] = collections.Counter()
    fixed_shape_components_checked = 0
    fixed_shape_component_mismatches = 0
    for match in exact_complete:
        region_a, region_b = match["region_a"], match["region_b"]
        if region_a != region_b:
            raise RuntimeError("exact-coordinate match had different region strings")
        status_a = region_status[(SAMPLE_A, region_a)]
        status_b = region_status[(SAMPLE_B, region_b)]
        forest_a = build_exact_forest(SAMPLE_A, region_a)
        forest_b = build_exact_forest(SAMPLE_B, region_b)
        exact_a = forest_a is not None
        exact_b = forest_b is not None
        both_exact = exact_a and exact_b
        uses_vaf_a = bool(forest_a and any(component[4] == "vaf_exact_top" for component in forest_a))
        uses_vaf_b = bool(forest_b and any(component[4] == "vaf_exact_top" for component in forest_b))
        both_vaf = both_exact and uses_vaf_a and uses_vaf_b
        unique_a = bool(forest_a and all(len(component[2]) == 1 for component in forest_a))
        unique_b = bool(forest_b and all(len(component[2]) == 1 for component in forest_b))
        both_unique = both_exact and unique_a and unique_b
        candidate_ordered_agree = (
            candidate_key(forest_a, False) == candidate_key(forest_b, False) if both_exact else None
        )
        candidate_swap_agree = (
            candidate_key(forest_a, True) == candidate_key(forest_b, True) if both_exact else None
        )

        selected_shape_a, shape_source_a = build_selected_shape(
            SAMPLE_A, region_a, match["shape_set_signature_a"]
        )
        selected_shape_b, shape_source_b = build_selected_shape(
            SAMPLE_B, region_b, match["shape_set_signature_b"]
        )

        # A fixed T=1 component is converted to the same canonical TS namespace
        # used by the structural catalog.  This prevents blank fixed-unit IDs
        # from being mistaken for missing shapes.
        for sample, region, original in (
            (SAMPLE_A, region_a, match["shape_set_signature_a"]),
            (SAMPLE_B, region_b, match["shape_set_signature_b"]),
        ):
            original_by_hp = dict(parse_labeled_set_signature(original))
            for unit in units_by_region[(sample, region)]:
                if unit["status"] != "fixed_T1":
                    continue
                fixed_shape_components_checked += 1
                generated_shape = fixed[(sample, region, unit["family"])][1]
                if generated_shape not in original_by_hp.get(unit["family"], ()):
                    fixed_shape_component_mismatches += 1
        provenance_pairs[(shape_source_a, shape_source_b)] += 1
        both_shape = selected_shape_a is not None and selected_shape_b is not None
        shape_ordered_agree = selected_shape_a == selected_shape_b if both_shape else None
        shape_swap_agree = (
            tuple(sorted(value for _, value in selected_shape_a))
            == tuple(sorted(value for _, value in selected_shape_b))
            if both_shape
            else None
        )

        both_actual_vaf_evaluable = truth(status_a["evaluable"]) and truth(status_b["evaluable"])
        both_vaf_single_shape = (
            both_actual_vaf_evaluable
            and forest_a is not None
            and forest_b is not None
            and all(len(component[3]) == 1 for component in forest_a)
            and all(len(component[3]) == 1 for component in forest_b)
        )
        strict_gt_050 = (
            both_actual_vaf_evaluable
            and truth(status_a["strict_gt_0.50_T0.05"])
            and truth(status_b["strict_gt_0.50_T0.05"])
        )
        ge_060 = (
            both_actual_vaf_evaluable
            and truth(status_a["ge_0.60_T0.05"])
            and truth(status_b["ge_0.60_T0.05"])
        )

        counters["exact_complete_pairs"] += 1
        counters["exact_selectable_a"] += exact_a
        counters["exact_selectable_b"] += exact_b
        counters["both_exact_selectable"] += both_exact
        counters["both_exact_ordered_agree"] += candidate_ordered_agree is True
        counters["both_exact_swap_agree"] += candidate_swap_agree is True
        counters["both_unique_exact"] += both_unique
        counters["both_unique_exact_ordered_agree"] += both_unique and candidate_ordered_agree is True
        counters["both_unique_exact_swap_agree"] += both_unique and candidate_swap_agree is True
        counters["both_actual_vaf_evaluable"] += both_actual_vaf_evaluable
        counters["both_actual_vaf_exact_ordered_agree"] += both_actual_vaf_evaluable and candidate_ordered_agree is True
        counters["both_actual_vaf_exact_swap_agree"] += both_actual_vaf_evaluable and candidate_swap_agree is True
        counters["both_actual_vaf_unique_exact"] += both_actual_vaf_evaluable and both_unique
        counters["both_actual_vaf_unique_ordered_agree"] += (
            both_actual_vaf_evaluable and both_unique and candidate_ordered_agree is True
        )
        counters["both_actual_vaf_unique_swap_agree"] += (
            both_actual_vaf_evaluable and both_unique and candidate_swap_agree is True
        )
        counters["shape_selectable_a"] += selected_shape_a is not None
        counters["shape_selectable_b"] += selected_shape_b is not None
        counters["both_shape_selectable"] += both_shape
        counters["both_shape_ordered_agree"] += shape_ordered_agree is True
        counters["both_shape_swap_agree"] += shape_swap_agree is True
        counters["both_vaf_shape_sources"] += (
            shape_source_a == "vaf_exact_top_shape" and shape_source_b == "vaf_exact_top_shape"
        )
        counters["both_vaf_shape_sources_ordered_agree"] += (
            shape_source_a == "vaf_exact_top_shape"
            and shape_source_b == "vaf_exact_top_shape"
            and shape_ordered_agree is True
        )
        counters["both_vaf_shape_sources_swap_agree"] += (
            shape_source_a == "vaf_exact_top_shape"
            and shape_source_b == "vaf_exact_top_shape"
            and shape_swap_agree is True
        )
        counters["both_actual_vaf_single_shape"] += both_vaf_single_shape
        counters["both_actual_vaf_single_shape_ordered_agree"] += (
            both_vaf_single_shape and shape_key(forest_a, False) == shape_key(forest_b, False)
        )
        counters["both_actual_vaf_single_shape_swap_agree"] += (
            both_vaf_single_shape and shape_key(forest_a, True) == shape_key(forest_b, True)
        )
        for gate, enabled in (("gt050", strict_gt_050), ("ge060", ge_060)):
            counters[f"both_{gate}"] += enabled
            counters[f"both_{gate}_candidate_ordered_agree"] += enabled and candidate_ordered_agree is True
            counters[f"both_{gate}_candidate_swap_agree"] += enabled and candidate_swap_agree is True
            counters[f"both_{gate}_shape_ordered_agree"] += (
                enabled and shape_key(forest_a, False) == shape_key(forest_b, False)
            )
            counters[f"both_{gate}_shape_swap_agree"] += (
                enabled and shape_key(forest_a, True) == shape_key(forest_b, True)
            )

        unranked_shape_ordered = match["shape_set_signature_a"] == match["shape_set_signature_b"]
        unranked_shape_swap = unlabeled_signature(match["shape_set_signature_a"]) == unlabeled_signature(
            match["shape_set_signature_b"]
        )
        unranked_tree_ordered = match["tree_digest_signature_a"] == match["tree_digest_signature_b"]
        unranked_tree_swap = unlabeled_signature(match["tree_digest_signature_a"]) == unlabeled_signature(
            match["tree_digest_signature_b"]
        )
        counters["unranked_shape_ordered_agree"] += unranked_shape_ordered
        counters["unranked_shape_swap_agree"] += unranked_shape_swap
        counters["unranked_tree_ordered_agree"] += unranked_tree_ordered
        counters["unranked_tree_swap_agree"] += unranked_tree_swap
        counters["vaf_subset_unranked_shape_ordered_agree"] += both_actual_vaf_evaluable and unranked_shape_ordered
        counters["vaf_subset_unranked_shape_swap_agree"] += both_actual_vaf_evaluable and unranked_shape_swap
        counters["vaf_subset_unranked_tree_ordered_agree"] += both_actual_vaf_evaluable and unranked_tree_ordered
        counters["vaf_subset_unranked_tree_swap_agree"] += both_actual_vaf_evaluable and unranked_tree_swap

        def forest_json(forest):
            if forest is None:
                return ""
            return json.dumps(
                [
                    {
                        "HP": family,
                        "positions": list(positions),
                        "candidate_ids": list(candidates),
                        "shape_ids": list(shapes),
                        "source": source,
                    }
                    for family, positions, candidates, shapes, source in forest
                ],
                ensure_ascii=False,
                separators=(",", ":"),
            )

        rows_out.append(
            {
                "match_id": match["match_id"],
                "region": region_a,
                "structural_class_a": status_a["structural_class"],
                "structural_class_b": status_b["structural_class"],
                "vaf_final_class_a": status_a["final_class"],
                "vaf_final_class_b": status_b["final_class"],
                "exact_forest_a": forest_json(forest_a),
                "exact_forest_b": forest_json(forest_b),
                "exact_selectable_a": exact_a,
                "exact_selectable_b": exact_b,
                "both_exact_selectable": both_exact,
                "uses_vaf_a": uses_vaf_a,
                "uses_vaf_b": uses_vaf_b,
                "candidate_top_combinations_a": math.prod(len(component[2]) for component in forest_a) if forest_a else "",
                "candidate_top_combinations_b": math.prod(len(component[2]) for component in forest_b) if forest_b else "",
                "unique_exact_a": unique_a,
                "unique_exact_b": unique_b,
                "candidate_ordered_agree": "" if candidate_ordered_agree is None else candidate_ordered_agree,
                "candidate_swap_tolerant_agree": "" if candidate_swap_agree is None else candidate_swap_agree,
                "shape_source_a": shape_source_a,
                "shape_source_b": shape_source_b,
                "selected_shape_a": "" if selected_shape_a is None else json.dumps(selected_shape_a, separators=(",", ":")),
                "selected_shape_b": "" if selected_shape_b is None else json.dumps(selected_shape_b, separators=(",", ":")),
                "both_shape_selectable": both_shape,
                "shape_ordered_agree": "" if shape_ordered_agree is None else shape_ordered_agree,
                "shape_swap_tolerant_agree": "" if shape_swap_agree is None else shape_swap_agree,
                "both_actual_vaf_evaluable": both_actual_vaf_evaluable,
                "both_actual_vaf_single_shape": both_vaf_single_shape,
                "both_strict_gt_0.50_T0.05": strict_gt_050,
                "both_ge_0.60_T0.05": ge_060,
                "min_top_weight_a_T0.05": status_a["min_top_weight_T0.05"],
                "min_top_weight_b_T0.05": status_b["min_top_weight_T0.05"],
            }
        )

    # Full-sample conservation for the official shape endpoint.
    sample_conservation = {}
    for sample in SAMPLES:
        complete = 0
        exact_selectable = 0
        unique_exact = 0
        structural_shape = 0
        vaf_shape = 0
        for (row_sample, region), status in region_status.items():
            if row_sample != sample or status["structural_class"] == "incomplete":
                continue
            complete += 1
            forest = build_exact_forest(sample, region)
            if forest is not None:
                exact_selectable += 1
                unique_exact += all(len(component[2]) == 1 for component in forest)
            if status["structural_class"] in {"T=1|Topo=1", "T>1|Topo=1"}:
                structural_shape += 1
            elif (
                status["structural_class"] == "T>1|Topo>1"
                and truth(status["evaluable"])
                and int(status["n_top_shapes_joint_exact"]) == 1
            ):
                vaf_shape += 1
        sample_conservation[sample] = {
            "complete_regions": complete,
            "exact_candidate_forest_selectable": exact_selectable,
            "unique_exact_candidate_forest": unique_exact,
            "structural_topo1_shape": structural_shape,
            "vaf_selected_topo_gt1_single_shape": vaf_shape,
            "structural_or_vaf_shape_selectable": structural_shape + vaf_shape,
        }

    metrics: list[dict] = []

    def add_pair(family, population, denominator_key, ordered_key, swap_key, interpretation):
        denominator = counters[denominator_key]
        metrics.append(metric_row(family, population, "ordered_HP", counters[ordered_key], denominator, interpretation))
        metrics.append(metric_row(family, population, "HP1_HP2_swap_tolerant", counters[swap_key], denominator, interpretation))

    add_pair(
        "unranked_full_candidate_shape_set",
        "all_exact_coordinate_complete_both",
        "exact_complete_pairs",
        "unranked_shape_ordered_agree",
        "unranked_shape_swap_agree",
        "Full feasible shape-set equality before VAF ranking; not a selected-tree endpoint.",
    )
    add_pair(
        "unranked_full_candidate_tree_set_digest",
        "all_exact_coordinate_complete_both",
        "exact_complete_pairs",
        "unranked_tree_ordered_agree",
        "unranked_tree_swap_agree",
        "Full feasible tree-set digest equality before VAF ranking.",
    )
    add_pair(
        "vaf_exact_first_candidate_set",
        "both_exact_candidate_forests_selectable",
        "both_exact_selectable",
        "both_exact_ordered_agree",
        "both_exact_swap_agree",
        "Exact VAF-score maximizer sets; ties retained; fixed T=1 components included.",
    )
    add_pair(
        "vaf_unique_exact_first_candidate",
        "both_unique_exact_candidate_forests",
        "both_unique_exact",
        "both_unique_exact_ordered_agree",
        "both_unique_exact_swap_agree",
        "Both sides have exactly one complete forest candidate after structural/VAF selection.",
    )
    add_pair(
        "vaf_exact_first_candidate_set",
        "both_sides_actually_use_vaf_ranking",
        "both_actual_vaf_evaluable",
        "both_actual_vaf_exact_ordered_agree",
        "both_actual_vaf_exact_swap_agree",
        "Both regions are VAF-evaluable; exact top candidate sets retain ties.",
    )
    add_pair(
        "vaf_unique_exact_first_candidate",
        "both_sides_actually_use_vaf_ranking_and_unique",
        "both_actual_vaf_unique_exact",
        "both_actual_vaf_unique_ordered_agree",
        "both_actual_vaf_unique_swap_agree",
        "Both VAF-evaluable sides select exactly one full forest candidate.",
    )
    add_pair(
        "selected_single_shape",
        "both_structural_or_vaf_shape_selectable",
        "both_shape_selectable",
        "both_shape_ordered_agree",
        "both_shape_swap_agree",
        "Original Topo=1 uses structural shape; original Topo>1 requires VAF-evaluable unique exact-top shape.",
    )
    add_pair(
        "selected_single_shape",
        "both_shapes_require_vaf_from_original_Topo_gt1",
        "both_vaf_shape_sources",
        "both_vaf_shape_sources_ordered_agree",
        "both_vaf_shape_sources_swap_agree",
        "Both original candidate shape sets had Topo>1 and VAF selected one shape.",
    )
    add_pair(
        "selected_single_shape",
        "both_sides_vaf_evaluable_and_single_shape",
        "both_actual_vaf_single_shape",
        "both_actual_vaf_single_shape_ordered_agree",
        "both_actual_vaf_single_shape_swap_agree",
        "Both regions use VAF exact ranking and each complete forest has one selected shape.",
    )
    for gate, label in (("gt050", "both_region_min_top_weight_strict_gt_0.50"), ("ge060", "both_region_min_top_weight_ge_0.60")):
        add_pair(
            "vaf_exact_first_candidate_set_sensitivity",
            label,
            f"both_{gate}",
            f"both_{gate}_candidate_ordered_agree",
            f"both_{gate}_candidate_swap_agree",
            "Relative softmax-weight gate at T=0.05; heuristic concentration, not calibrated probability.",
        )
        add_pair(
            "selected_single_shape_sensitivity",
            label,
            f"both_{gate}",
            f"both_{gate}_shape_ordered_agree",
            f"both_{gate}_shape_swap_agree",
            "The same gated subset compared at unlabeled rooted-shape level.",
        )

    checks = []

    def check(name: str, observed, expected) -> None:
        checks.append({"check": name, "pass": observed == expected, "observed": observed, "expected": expected})
        if observed != expected:
            raise RuntimeError(f"QA failed: {name}: {observed!r} != {expected!r}")

    check("exact_coordinate_match_count", len(exact_matches), 6252)
    check("exact_coordinate_complete_both_count", len(exact_complete), 5720)
    check("pair_ranked_unit_candidate_detail_join", ranked_units, len(top_candidates_from_detail))
    check("HCC1395_complete_region_conservation", sample_conservation[SAMPLE_A]["complete_regions"], 6940)
    check("HCC1395_DORADO_complete_region_conservation", sample_conservation[SAMPLE_B]["complete_regions"], 6750)
    check("HCC1395_shape_selectable_official", sample_conservation[SAMPLE_A]["structural_or_vaf_shape_selectable"], 6798)
    check("HCC1395_DORADO_shape_selectable_official", sample_conservation[SAMPLE_B]["structural_or_vaf_shape_selectable"], 6082)
    check("HCC1395_exact_candidate_forest_selectable", sample_conservation[SAMPLE_A]["exact_candidate_forest_selectable"], 6849)
    check("HCC1395_DORADO_exact_candidate_forest_selectable", sample_conservation[SAMPLE_B]["exact_candidate_forest_selectable"], 6250)
    check("HCC1395_unique_exact_candidate_forest", sample_conservation[SAMPLE_A]["unique_exact_candidate_forest"], 6395)
    check("HCC1395_DORADO_unique_exact_candidate_forest", sample_conservation[SAMPLE_B]["unique_exact_candidate_forest"], 5572)
    check("fixed_T1_shape_components_checked", fixed_shape_components_checked, 6554)
    check("fixed_T1_shape_component_mismatches", fixed_shape_component_mismatches, 0)
    check("pair_both_exact_candidate_forest_selectable", counters["both_exact_selectable"], 5347)
    check("pair_both_unique_exact_candidate_forest", counters["both_unique_exact"], 4602)
    check("pair_both_actual_vaf_evaluable", counters["both_actual_vaf_evaluable"], 3236)
    check("pair_both_actual_vaf_unique_exact", counters["both_actual_vaf_unique_exact"], 2543)
    check("pair_both_shape_selectable", counters["both_shape_selectable"], 5168)
    check("pair_both_vaf_shape_sources", counters["both_vaf_shape_sources"], 2089)
    check("region_output_rows", len(rows_out), len(exact_complete))

    fields = [
        "match_id",
        "region",
        "structural_class_a",
        "structural_class_b",
        "vaf_final_class_a",
        "vaf_final_class_b",
        "exact_forest_a",
        "exact_forest_b",
        "exact_selectable_a",
        "exact_selectable_b",
        "both_exact_selectable",
        "uses_vaf_a",
        "uses_vaf_b",
        "candidate_top_combinations_a",
        "candidate_top_combinations_b",
        "unique_exact_a",
        "unique_exact_b",
        "candidate_ordered_agree",
        "candidate_swap_tolerant_agree",
        "shape_source_a",
        "shape_source_b",
        "selected_shape_a",
        "selected_shape_b",
        "both_shape_selectable",
        "shape_ordered_agree",
        "shape_swap_tolerant_agree",
        "both_actual_vaf_evaluable",
        "both_actual_vaf_single_shape",
        "both_strict_gt_0.50_T0.05",
        "both_ge_0.60_T0.05",
        "min_top_weight_a_T0.05",
        "min_top_weight_b_T0.05",
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    regions_out = args.output_dir / "hcc1395_vaf_pair_regions.tsv"
    metrics_out = args.output_dir / "hcc1395_vaf_pair_metrics.tsv"
    checks_out = args.output_dir / "hcc1395_vaf_pair_checks.tsv"
    json_out = args.output_dir / "hcc1395_vaf_pair_summary.json"
    write_tsv(regions_out, rows_out, fields)
    write_tsv(
        metrics_out,
        metrics,
        ["metric_family", "population", "comparison", "numerator", "denominator", "percent", "interpretation"],
    )
    write_tsv(checks_out, checks, ["check", "pass", "observed", "expected"])

    source_paths = [args.vaf_units, args.vaf_candidates, args.vaf_regions, args.pair_matches, args.pair_analysis]
    summary = {
        "schema_version": "1.0",
        "analysis_date": "2026-07-12",
        "scope": "HCC1395 vs HCC1395_DORADO; chr1-22; exact-coordinate complete-both historical layered-v2 regions",
        "claim_ceiling": (
            "Same-read-AF ancestor-ordering heuristic inference only. Relative softmax weights are not calibrated "
            "probabilities; selected forests are not biological truth, clones, or independent validation."
        ),
        "definitions": {
            "candidate_exact_key": (
                "HP-labeled tuple of (ordered genomic positions, exact-score top candidate IDs); ties retained. "
                "The swap-tolerant key sorts the one/two HP components after dropping HP labels."
            ),
            "candidate_score": "sum over inferred ancestor-descendant pairs of read_AF(ancestor)-read_AF(descendant)",
            "relative_weight": "softmax(candidate ordering scores / temperature 0.05); heuristic concentration only",
            "shape_endpoint": (
                "Original Topo=1 is structurally selected without VAF. Original Topo>1 requires region "
                "evaluable=True and n_top_shapes_joint_exact=1 before the VAF-selected shape is materialized."
            ),
            "unranked_metric": "Full feasible candidate-set equality before VAF ranking; kept separate from top-set endpoints.",
        },
        "sources": [
            {"path": str(path.resolve()), "sha256": sha256(path)} for path in source_paths
        ]
        + [
            {"path": str(path.resolve()), "sha256": sha256(path), "role": "fixed_T1 tree materialization"}
            for path in layered_paths.values()
        ]
        + mlhp_sources,
        "source_row_counts": {
            "pair_vaf_units": len(units),
            "pair_ranked_units": ranked_units,
            "pair_candidate_rows": candidate_rows,
            "pair_group_regions": len(positions_by_region),
            "exact_coordinate_matches": len(exact_matches),
            "exact_coordinate_complete_both": len(exact_complete),
        },
        "sample_conservation": sample_conservation,
        "pair_counters": dict(sorted(counters.items())),
        "shape_provenance_pairs": [
            {"source_a": key[0], "source_b": key[1], "count": count}
            for key, count in sorted(provenance_pairs.items())
        ],
        "metrics": metrics,
        "validation": {"passed": sum(truth(row["pass"]) for row in checks), "total": len(checks), "checks": checks},
        "outputs": {
            "regions_tsv": str(regions_out.resolve()),
            "metrics_tsv": str(metrics_out.resolve()),
            "checks_tsv": str(checks_out.resolve()),
        },
    }
    dump_json(json_out, summary)
    print(f"INPUT VAF UNITS -> {args.vaf_units.resolve()}")
    print(f"INPUT VAF CANDIDATES -> {args.vaf_candidates.resolve()}")
    print(f"INPUT VAF REGIONS -> {args.vaf_regions.resolve()}")
    print(f"INPUT PAIR MATCHES -> {args.pair_matches.resolve()}")
    print(f"INPUT RUN ROOT -> {args.run_root.resolve()}")
    print(f"OUTPUT REGIONS -> {regions_out.resolve()}")
    print(f"OUTPUT METRICS -> {metrics_out.resolve()}")
    print(f"OUTPUT SUMMARY -> {json_out.resolve()}")
    print(f"QA -> {len(checks)}/{len(checks)} PASS")
    print(
        "CORE -> "
        f"exact selectable {counters['both_exact_selectable']}; unique {counters['both_unique_exact']}; "
        f"actual VAF {counters['both_actual_vaf_evaluable']}; shape selectable {counters['both_shape_selectable']}"
    )


if __name__ == "__main__":
    main()
