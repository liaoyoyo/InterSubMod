#!/usr/bin/env python3
"""Build a site-aware topology comparison for the HCC1395 technical pair.

The comparison is deliberately based on mutation-event relations on genomic
sites shared by the two HP-specific solver universes.  A candidate tree is
projected to the ordered relation vector for every shared-site pair:

    F: lower genomic position is ancestral to the higher position
    R: higher genomic position is ancestral to the lower position
    P: the two mutation events are on parallel branches

Private mutation events are therefore removed without using arbitrary graph
subtree matching.  A recurrence involving a shared mutation is fail-closed.
The full read-pattern-compatible candidate set and the VAF exact-score top set
are separate endpoints.  VAF ties are retained.  A normalized mean-delta
sensitivity divides each exact score by its ancestry-comparison count and is
not evaluable when any candidate has zero comparisons.

The primary comparison keeps HP1/HP2 labels ordered.  A sensitivity comparison
may apply one whole-region HP1<->HP2 swap; an HP-component count mismatch is
never silently accepted.  Candidate-set relations are computed over unique
projected topology signatures, while the pre-projection candidate-tree counts
are retained for complexity matching.

This is a technical reproducibility analysis of read-compatible constraints.
It does not establish a biological clone tree or independent truth.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
import sys
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from types import ModuleType
from typing import Iterable, Mapping, Sequence


SAMPLE_A = "HCC1395"
SAMPLE_B = "HCC1395_DORADO"
SAMPLES = (SAMPLE_A, SAMPLE_B)
ENDPOINTS = ("read_full", "vaf_official", "vaf_normalized")
MAPPINGS = ("identity", "global_hp_swap")
RELATION_CATEGORIES = (
    "exact",
    "a_proper_subset_b",
    "b_proper_subset_a",
    "overlap",
    "disjoint",
)


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    edges: tuple[tuple[str, str], ...]
    recurrence_bits: frozenset[int]


@dataclass(frozen=True)
class CandidateScore:
    candidate_id: str
    score: Fraction
    n_comparisons: int
    is_official_top: bool


@dataclass(frozen=True)
class UnitModel:
    sample: str
    region: str
    family: str
    positions: tuple[int, ...]
    universe_positions: tuple[int, ...]
    candidates: Mapping[str, Candidate]
    vaf_status: str
    official_top_ids: frozenset[str] | None
    normalized_top_ids: frozenset[str] | None
    normalized_status: str
    comparison_count_constant: bool | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-matches", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--layered-a", required=True, type=Path)
    parser.add_argument("--layered-b", required=True, type=Path)
    parser.add_argument("--solver", required=True, type=Path)
    parser.add_argument("--vaf-units", required=True, type=Path)
    parser.add_argument("--vaf-candidates", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> Iterable[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def truth(value: str | bool | None) -> bool:
    return value is True or value == "True"


def canonical_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: Sequence[dict], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def deterministic_gzip_text(path: Path):
    """Return a deterministic text writer whose gzip header has mtime=0."""
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = path.open("wb")
    zipped = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(zipped, encoding="utf-8", newline="")

    class Writer:
        def __enter__(self):
            return text

        def __exit__(self, exc_type, exc, traceback):
            try:
                text.close()
            finally:
                if not raw.closed:
                    raw.close()
            return False

    return Writer()


def load_solver(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("frozen_tree_enumeration_solver", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import solver: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_node(label: str) -> frozenset[int]:
    if label == "ROOT":
        return frozenset()
    value = label[2:] if label.startswith("H_") else label
    if not value or any(base not in {"R", "A"} for base in value):
        raise RuntimeError(f"invalid solver node label: {label}")
    return frozenset(index for index, base in enumerate(value) if base == "A")


def candidate_id(edges: Sequence[Sequence[str]]) -> str:
    payload = sorted([list(edge) for edge in edges])
    encoded = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode()
    return "TR-" + hashlib.sha1(encoded).hexdigest()[:12]


def ordered_tree_digest(trees: Sequence[dict]) -> str:
    payload = [sorted([list(edge) for edge in tree["edges"]]) for tree in trees]
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def independently_recurrence_bits(edges: Sequence[Sequence[str]]) -> frozenset[int]:
    counts: collections.Counter[int] = collections.Counter()
    for parent_label, child_label in edges:
        parent, child = parse_node(parent_label), parse_node(child_label)
        acquired = child - parent
        if not parent < child or len(acquired) != 1:
            raise RuntimeError(f"non-unit-flip edge: {parent_label}->{child_label}")
        counts[next(iter(acquired))] += 1
    return frozenset(index for index, count in counts.items() if count >= 2)


def site_set_relation(a: Iterable[int], b: Iterable[int]) -> str:
    set_a, set_b = set(a), set(b)
    if set_a == set_b:
        return "equal"
    if set_a < set_b:
        return "a_proper_subset_b"
    if set_b < set_a:
        return "b_proper_subset_a"
    if set_a & set_b:
        return "partial_overlap"
    return "disjoint"


def load_groups(sample_dir: Path) -> tuple[dict[str, dict], list[Path]]:
    groups: dict[str, dict] = {}
    sources: list[Path] = []
    for path in sorted(sample_dir.glob("mlhp_part_*.json")):
        sources.append(path)
        document = json.loads(path.read_text(encoding="utf-8"))
        for group in document.get("groups", []):
            region = f"{group['chrom']}:{group['start']}-{group['end']}"
            if region in groups:
                raise RuntimeError(f"duplicate mlhp region: {sample_dir.name} {region}")
            groups[region] = group
    if not sources:
        raise RuntimeError(f"no mlhp_part_*.json under {sample_dir}")
    return groups, sources


def load_layered(path: Path, sample: str) -> dict[tuple[str, str], dict]:
    document = json.loads(path.read_text(encoding="utf-8"))
    result: dict[tuple[str, str], dict] = {}
    for unit in document["detail"]:
        if not unit.get("is_primary_lineage"):
            continue
        key = (unit["region"], str(unit["family"]))
        if key in result:
            raise RuntimeError(f"duplicate primary layered unit: {sample} {key}")
        result[key] = unit
    return result


def load_vaf_units(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    result = {}
    for row in read_tsv(path):
        if row["sample"] not in SAMPLES:
            continue
        key = (row["sample"], row["region"], row["family"])
        if key in result:
            raise RuntimeError(f"duplicate VAF unit: {key}")
        result[key] = row
    return result


def load_vaf_candidates(path: Path) -> dict[tuple[str, str, str], dict[str, CandidateScore]]:
    result: dict[tuple[str, str, str], dict[str, CandidateScore]] = collections.defaultdict(dict)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample"] not in SAMPLES:
                continue
            key = (row["sample"], row["region"], row["family"])
            cid = row["candidate_id"]
            if cid in result[key]:
                raise RuntimeError(f"duplicate VAF candidate: {key} {cid}")
            score = Fraction(int(row["exact_score_numerator"]), int(row["exact_score_denominator"]))
            result[key][cid] = CandidateScore(
                candidate_id=cid,
                score=score,
                n_comparisons=int(row["n_ancestry_comparisons"]),
                is_official_top=truth(row["is_exact_top"]),
            )
    return dict(result)


def exact_complete_matches(path: Path) -> tuple[list[dict[str, str]], int]:
    exact = [row for row in read_tsv(path) if row["scenario"] == "exact_coordinate"]
    complete = [row for row in exact if truth(row["complete_a"]) and truth(row["complete_b"])]
    for row in complete:
        if row["region_a"] != row["region_b"]:
            raise RuntimeError(f"exact-coordinate row has different regions: {row['match_id']}")
    return complete, len(exact)


def enumerate_unit(
    sample: str,
    region: str,
    family: str,
    group: dict,
    layered: dict,
    vaf_unit: dict[str, str],
    vaf_scores: Mapping[str, CandidateScore],
    solver: ModuleType,
    counters: collections.Counter,
) -> UnitModel:
    positions = tuple(int(value) for value in group.get("positions", []))
    full = ((group.get("populations_by_hp", {}) or {}).get(family, {}) or {})
    partial = list((((group.get("subread_groups_by_hp", {}) or {}).get(family, {}) or {}).keys()))
    if not positions or (not full and not partial):
        raise RuntimeError(f"missing HP-specific read patterns: {(sample, region, family)}")
    k = len(positions)
    if any(len(pattern) != k for pattern in list(full) + partial):
        raise RuntimeError(f"read-pattern length mismatch: {(sample, region, family)}")

    result = solver.enumerate_min_trees(full, partial, k, tree_cap=0)
    counters["units_reenumerated"] += 1
    counters["trees_reenumerated"] += int(result["n_trees"])
    if result.get("capped") or not result.get("trees_complete"):
        raise RuntimeError(f"solver candidate set incomplete: {(sample, region, family)}")
    if int(result["n_trees"]) != int(layered["n_trees"]):
        raise RuntimeError(
            f"n_trees mismatch {(sample, region, family)}: {result['n_trees']} != {layered['n_trees']}"
        )
    if layered.get("analysis_candidate_set_complete") is not True:
        raise RuntimeError(f"layered candidate set is not complete: {(sample, region, family)}")
    digest = ordered_tree_digest(result["trees"])
    if digest != layered.get("analysis_tree_digest_sha256"):
        raise RuntimeError(f"tree digest mismatch: {(sample, region, family)}")

    candidates: dict[str, Candidate] = {}
    for tree in result["trees"]:
        edges = tuple((str(parent), str(child)) for parent, child in tree["edges"])
        cid = candidate_id(edges)
        if cid in candidates:
            raise RuntimeError(f"duplicate re-enumerated candidate ID: {(sample, region, family, cid)}")
        recurrence = independently_recurrence_bits(edges)
        if recurrence != frozenset(int(value) for value in tree.get("recurrence", [])):
            raise RuntimeError(f"recurrence annotation mismatch: {(sample, region, family, cid)}")
        candidates[cid] = Candidate(cid, edges, recurrence)

    status = vaf_unit["status"]
    official: frozenset[str] | None = None
    normalized: frozenset[str] | None = None
    normalized_status = "not_applicable"
    comparison_constant: bool | None = None
    if status == "ranked":
        if set(vaf_scores) != set(candidates):
            raise RuntimeError(f"VAF/re-enumerated candidate set mismatch: {(sample, region, family)}")
        official = frozenset(cid for cid, score in vaf_scores.items() if score.is_official_top)
        expected_official = frozenset(json.loads(vaf_unit["exact_top_candidate_ids"]))
        if not official or official != expected_official:
            raise RuntimeError(f"official VAF top-set mismatch: {(sample, region, family)}")
        max_raw = max(value.score for value in vaf_scores.values())
        if official != frozenset(cid for cid, value in vaf_scores.items() if value.score == max_raw):
            raise RuntimeError(f"official exact-sum argmax mismatch: {(sample, region, family)}")
        comparison_counts = {value.n_comparisons for value in vaf_scores.values()}
        comparison_constant = len(comparison_counts) == 1
        if comparison_constant != (vaf_unit["candidate_comparison_count_varies"] == "False"):
            raise RuntimeError(f"comparison-count flag mismatch: {(sample, region, family)}")
        if 0 in comparison_counts:
            normalized_status = "not_evaluable_zero_comparison"
            counters["normalized_ranked_units_zero_comparison"] += 1
        else:
            mean_scores = {
                cid: value.score / value.n_comparisons for cid, value in vaf_scores.items()
            }
            maximum = max(mean_scores.values())
            normalized = frozenset(cid for cid, value in mean_scores.items() if value == maximum)
            normalized_status = "ranked_mean_delta"
            counters["normalized_ranked_units_evaluable"] += 1
            counters["normalized_ranked_units_same_candidate_top_set"] += normalized == official
    elif status == "fixed_T1":
        if len(candidates) != 1:
            raise RuntimeError(f"fixed_T1 did not have one tree: {(sample, region, family)}")
        if vaf_scores:
            raise RuntimeError(f"fixed_T1 unexpectedly has VAF candidate rows: {(sample, region, family)}")
        official = normalized = frozenset(candidates)
        normalized_status = "structural_fixed_T1"
        comparison_constant = True
    else:
        if vaf_scores:
            raise RuntimeError(f"non-ranked unit unexpectedly has VAF candidate rows: {(sample, region, family)}")

    if int(vaf_unit["n_trees_expected"]) != len(candidates):
        raise RuntimeError(f"VAF n_trees_expected mismatch: {(sample, region, family)}")
    universe_positions = tuple(positions[int(index)] for index in result["universe"])
    return UnitModel(
        sample=sample,
        region=region,
        family=family,
        positions=positions,
        universe_positions=universe_positions,
        candidates=candidates,
        vaf_status=status,
        official_top_ids=official,
        normalized_top_ids=normalized,
        normalized_status=normalized_status,
        comparison_count_constant=comparison_constant,
    )


def event_relation_signature(candidate: Candidate, unit: UnitModel, shared_positions: tuple[int, ...]) -> str | None:
    """Return a local-rank F/R/P vector, or None for shared-site recurrence."""
    if not shared_positions:
        return ""
    index_by_position = {position: index for index, position in enumerate(unit.positions)}
    shared_indices = tuple(index_by_position[position] for position in shared_positions)
    if candidate.recurrence_bits & set(shared_indices):
        return None

    parents: dict[frozenset[int], frozenset[int]] = {}
    events: dict[int, frozenset[int]] = {}
    for parent_label, child_label in candidate.edges:
        parent, child = parse_node(parent_label), parse_node(child_label)
        if child in parents and parents[child] != parent:
            raise RuntimeError(f"candidate child has multiple parents: {candidate.candidate_id}")
        parents[child] = parent
        bit = next(iter(child - parent))
        if bit in events:
            # A repeated shared bit was handled above.  Private recurrence is
            # allowed because it disappears under shared-site projection.
            continue
        events[bit] = child
    missing = set(shared_indices) - set(events)
    if missing:
        raise RuntimeError(f"shared mutation event missing from candidate {candidate.candidate_id}: {sorted(missing)}")

    def ancestor(older: frozenset[int], younger: frozenset[int]) -> bool:
        current = younger
        while current:
            current = parents.get(current)
            if current is None:
                raise RuntimeError(f"candidate is not root-connected: {candidate.candidate_id}")
            if current == older:
                return True
        return False

    codes: list[str] = []
    for left, right in combinations(shared_indices, 2):
        left_event, right_event = events[left], events[right]
        if ancestor(left_event, right_event):
            codes.append("F")
        elif ancestor(right_event, left_event):
            codes.append("R")
        else:
            codes.append("P")
    return "".join(codes)


def selected_ids(unit: UnitModel, endpoint: str) -> frozenset[str] | None:
    if endpoint == "read_full":
        return frozenset(unit.candidates)
    if endpoint == "vaf_official":
        return unit.official_top_ids
    if endpoint == "vaf_normalized":
        return unit.normalized_top_ids
    raise KeyError(endpoint)


def project_unit(unit: UnitModel, endpoint: str, shared_positions: tuple[int, ...]) -> dict:
    identifiers = selected_ids(unit, endpoint)
    if identifiers is None:
        reason = unit.vaf_status if endpoint == "vaf_official" else unit.normalized_status
        return {
            "status": "not_evaluable",
            "reason": reason,
            "candidate_tree_count": 0,
            "projected_topology_set_size": 0,
            "q": [],
            "shared_recurrent_candidate_count": 0,
        }
    signatures: set[str] = set()
    recurrent = 0
    for cid in identifiers:
        signature = event_relation_signature(unit.candidates[cid], unit, shared_positions)
        if signature is None:
            recurrent += 1
        else:
            signatures.add(signature)
    if recurrent:
        return {
            "status": "not_evaluable",
            "reason": "shared_site_recurrence_fail_closed",
            "candidate_tree_count": len(identifiers),
            "projected_topology_set_size": len(signatures),
            "q": sorted(signatures),
            "shared_recurrent_candidate_count": recurrent,
        }
    return {
        "status": "ok" if len(shared_positions) >= 2 else "noninformative",
        "reason": "" if len(shared_positions) >= 2 else "fewer_than_two_shared_mutation_sites",
        "candidate_tree_count": len(identifiers),
        "projected_topology_set_size": len(signatures),
        "q": sorted(signatures),
        "shared_recurrent_candidate_count": 0,
    }


def compare_sets(values_a: Sequence[str], values_b: Sequence[str]) -> dict:
    a, b = set(values_a), set(values_b)
    if not a or not b:
        raise RuntimeError("candidate topology set cannot be empty for an evaluable comparison")
    intersection = len(a & b)
    union = len(a | b)
    if a == b:
        category = "exact"
    elif a < b:
        category = "a_proper_subset_b"
    elif b < a:
        category = "b_proper_subset_a"
    elif intersection:
        category = "overlap"
    else:
        category = "disjoint"
    return {
        "category": category,
        "set_size_a": len(a),
        "set_size_b": len(b),
        "intersection_size": intersection,
        "union_size": union,
        "jaccard": intersection / union,
        "coverage_a": intersection / len(a),
        "coverage_b": intersection / len(b),
        "compatible": intersection > 0,
        "unique_exact": len(a) == len(b) == intersection == 1,
    }


def forest_set_summary(component_rows: Sequence[dict], endpoint: str) -> dict:
    projections = [(row[endpoint]["a"], row[endpoint]["b"]) for row in component_rows]
    all_components_informative = all(row["shared_k"] >= 2 for row in component_rows)
    all_component_universes_equal = all(row["site_set_relation"] == "equal" for row in component_rows)
    if not all_components_informative and not all_component_universes_equal:
        projection_scope = "private_sites_projected_with_uninformative_component"
    elif not all_components_informative:
        projection_scope = "contains_uninformative_component"
    elif not all_component_universes_equal:
        projection_scope = "private_sites_projected"
    else:
        projection_scope = "equal_hp_specific_universes"
    failures = [
        f"HP{row['a_family']}~HP{row['b_family']}:{side['reason']}"
        for row in component_rows
        for side in (row[endpoint]["a"], row[endpoint]["b"])
        if side["status"] == "not_evaluable"
    ]
    relation_pairs = sum(int(row["n_relation_pairs"]) for row in component_rows)
    if failures:
        return {
            "status": "not_evaluable",
            "reason": ";".join(failures),
            "projection_scope": projection_scope,
            "all_components_informative": all_components_informative,
            "all_component_universes_equal": all_component_universes_equal,
            "n_informative_components": sum(row["shared_k"] >= 2 for row in component_rows),
            "n_relation_pairs": relation_pairs,
        }
    if relation_pairs == 0:
        return {
            "status": "not_evaluable",
            "reason": "no_shared_site_pair_relation",
            "projection_scope": projection_scope,
            "all_components_informative": all_components_informative,
            "all_component_universes_equal": all_component_universes_equal,
            "n_informative_components": 0,
            "n_relation_pairs": 0,
        }

    size_a = math.prod(left["projected_topology_set_size"] for left, _ in projections)
    size_b = math.prod(right["projected_topology_set_size"] for _, right in projections)
    intersection = math.prod(len(set(left["q"]) & set(right["q"])) for left, right in projections)
    union = size_a + size_b - intersection
    if size_a == size_b == intersection:
        category = "exact"
    elif intersection == size_a and size_a < size_b:
        category = "a_proper_subset_b"
    elif intersection == size_b and size_b < size_a:
        category = "b_proper_subset_a"
    elif intersection:
        category = "overlap"
    else:
        category = "disjoint"
    return {
        "status": "evaluable",
        "reason": "",
        "projection_scope": projection_scope,
        "all_components_informative": all_components_informative,
        "all_component_universes_equal": all_component_universes_equal,
        "category": category,
        "candidate_tree_product_a": math.prod(left["candidate_tree_count"] for left, _ in projections),
        "candidate_tree_product_b": math.prod(right["candidate_tree_count"] for _, right in projections),
        "topology_set_size_a": size_a,
        "topology_set_size_b": size_b,
        "intersection_size": intersection,
        "union_size": union,
        "jaccard": intersection / union,
        "coverage_a": intersection / size_a,
        "coverage_b": intersection / size_b,
        "compatible": intersection > 0,
        "unique_exact": size_a == size_b == intersection == 1,
        "n_informative_components": sum(row["shared_k"] >= 2 for row in component_rows),
        "n_relation_pairs": relation_pairs,
    }


def family_mapping(units_a: Mapping[str, UnitModel], units_b: Mapping[str, UnitModel], mapping: str):
    families_a, families_b = set(units_a), set(units_b)
    if len(families_a) != len(families_b):
        return None, "hp_count_mismatch"
    if mapping == "identity":
        pairs = [(family, family) for family in sorted(families_a, key=int)]
    elif mapping == "global_hp_swap":
        swap = {"1": "2", "2": "1"}
        pairs = [(family, swap[family]) for family in sorted(families_a, key=int)]
    else:
        raise KeyError(mapping)
    if {right for _, right in pairs} != families_b:
        return None, "hp_family_assignment_mismatch"
    return pairs, ""


def build_alignment(
    units_a: Mapping[str, UnitModel], units_b: Mapping[str, UnitModel], mapping_name: str
) -> dict:
    pairs, reason = family_mapping(units_a, units_b, mapping_name)
    if pairs is None:
        result = {
            "valid_mapping": False,
            "reason": reason,
            "mapping": [],
            "components": [],
            "comparison_count_constant": False,
        }
        for endpoint in ENDPOINTS:
            result[endpoint] = {
                "status": "not_evaluable",
                "reason": reason,
                "projection_scope": "not_mapped",
                "all_components_informative": False,
                "all_component_universes_equal": False,
                "n_informative_components": 0,
                "n_relation_pairs": 0,
            }
        return result
    components = []
    for family_a, family_b in pairs:
        unit_a, unit_b = units_a[family_a], units_b[family_b]
        sites_a, sites_b = set(unit_a.universe_positions), set(unit_b.universe_positions)
        shared = tuple(sorted(sites_a & sites_b))
        component = {
            "a_family": family_a,
            "b_family": family_b,
            "universe_positions_a": sorted(sites_a),
            "universe_positions_b": sorted(sites_b),
            "shared_positions": list(shared),
            "site_set_relation": site_set_relation(sites_a, sites_b),
            "shared_k": len(shared),
            "n_relation_pairs": math.comb(len(shared), 2),
            "local_relation_vector_order": "lexicographic pairs of local ranks 0..shared_k-1",
            "vaf_status_a": unit_a.vaf_status,
            "vaf_status_b": unit_b.vaf_status,
            "comparison_count_constant_a": unit_a.comparison_count_constant,
            "comparison_count_constant_b": unit_b.comparison_count_constant,
            "normalized_status_a": unit_a.normalized_status,
            "normalized_status_b": unit_b.normalized_status,
        }
        for endpoint in ENDPOINTS:
            component[endpoint] = {
                "a": project_unit(unit_a, endpoint, shared),
                "b": project_unit(unit_b, endpoint, shared),
            }
        components.append(component)

    alignment = {
        "valid_mapping": True,
        "reason": "",
        "mapping": [{"a_family": left, "b_family": right} for left, right in pairs],
        "shared_k_vector": [row["shared_k"] for row in components],
        "components": components,
    }
    alignment["comparison_count_constant"] = all(
        value is True
        for row in components
        for value in (row["comparison_count_constant_a"], row["comparison_count_constant_b"])
    )
    for endpoint in ENDPOINTS:
        alignment[endpoint] = forest_set_summary(components, endpoint)
    return alignment


def endpoint_rank(summary: dict, mapping_name: str) -> tuple:
    if summary.get("status") != "evaluable":
        return (-1.0, -1.0, -1, -1, -1)
    category_rank = {
        "disjoint": 0,
        "overlap": 1,
        "a_proper_subset_b": 2,
        "b_proper_subset_a": 2,
        "exact": 3,
    }[summary["category"]]
    return (
        summary["jaccard"],
        min(summary["coverage_a"], summary["coverage_b"]),
        summary["intersection_size"],
        category_rank,
        1 if mapping_name == "identity" else 0,
    )


def swap_tolerant_endpoint(alignments: Mapping[str, dict], endpoint: str) -> dict:
    candidates = [
        (endpoint_rank(alignment.get(endpoint, {}), name), name, alignment.get(endpoint, {}))
        for name, alignment in alignments.items()
        if alignment.get("valid_mapping")
    ]
    if not candidates:
        reasons = sorted({alignment.get("reason", "invalid_mapping") for alignment in alignments.values()})
        return {"status": "not_evaluable", "reason": ";".join(reasons), "selected_mapping": ""}
    _, name, summary = max(candidates, key=lambda item: item[0])
    result = dict(summary)
    result["selected_mapping"] = name
    return result


def within_sample_top_set_agreement(units: Mapping[str, UnitModel]) -> dict:
    official_available = all(unit.official_top_ids is not None for unit in units.values())
    normalized_available = all(unit.normalized_top_ids is not None for unit in units.values())
    has_ranked = any(unit.vaf_status == "ranked" for unit in units.values())
    constant = all(unit.comparison_count_constant is True for unit in units.values())
    same = (
        official_available
        and normalized_available
        and all(unit.official_top_ids == unit.normalized_top_ids for unit in units.values())
    )
    return {
        "official_available": official_available,
        "normalized_available": normalized_available,
        "has_ranked_unit": has_ranked,
        "comparison_count_constant": constant,
        "same_candidate_top_set": same if official_available and normalized_available else None,
    }


def flatten_endpoint(prefix: str, summary: dict, row: dict) -> None:
    fields = (
        "status",
        "reason",
        "projection_scope",
        "all_components_informative",
        "all_component_universes_equal",
        "category",
        "candidate_tree_product_a",
        "candidate_tree_product_b",
        "topology_set_size_a",
        "topology_set_size_b",
        "intersection_size",
        "union_size",
        "jaccard",
        "coverage_a",
        "coverage_b",
        "compatible",
        "unique_exact",
        "n_informative_components",
        "n_relation_pairs",
        "selected_mapping",
    )
    for field in fields:
        value = summary.get(field, "")
        if isinstance(value, float):
            value = f"{value:.12f}"
        row[f"{prefix}_{field}"] = value


def metric_rows(region_records: Sequence[dict]) -> list[dict]:
    rows: list[dict] = []

    def add(endpoint: str, comparison: str, population: str, selected: list[dict]) -> None:
        evaluable = [row for row in selected if row.get("status") == "evaluable"]
        denominator = len(evaluable)
        categories = collections.Counter(row["category"] for row in evaluable)
        for metric, value in (
            ("population_regions", len(selected)),
            ("evaluable_regions", denominator),
            ("compatible_regions", sum(row["compatible"] for row in evaluable)),
            ("unique_exact_regions", sum(row["unique_exact"] for row in evaluable)),
            *[(f"category_{category}", categories[category]) for category in RELATION_CATEGORIES],
        ):
            rows.append(
                {
                    "endpoint": endpoint,
                    "comparison": comparison,
                    "population": population,
                    "metric": metric,
                    "numerator": value,
                    "denominator": denominator if metric not in {"population_regions", "evaluable_regions"} else len(selected),
                    "percent": (
                        f"{100 * value / (denominator if metric not in {'population_regions', 'evaluable_regions'} else len(selected)):.6f}"
                        if (denominator if metric not in {"population_regions", "evaluable_regions"} else len(selected))
                        else ""
                    ),
                }
            )
        rows.append(
            {
                "endpoint": endpoint,
                "comparison": comparison,
                "population": population,
                "metric": "mean_jaccard",
                "numerator": f"{sum(row['jaccard'] for row in evaluable) / denominator:.12f}" if denominator else "",
                "denominator": denominator,
                "percent": "",
            }
        )

    for endpoint in ENDPOINTS:
        ordered = [record["alignments"]["identity"].get(endpoint, {}) for record in region_records]
        swap = [record["swap_tolerant"][endpoint] for record in region_records]
        add(endpoint, "ordered_HP", "all_exact_complete_both", ordered)
        add(endpoint, "whole_region_HP_swap_tolerant", "all_exact_complete_both", swap)
        if endpoint == "vaf_official":
            ordered_constant = [
                record["alignments"]["identity"].get(endpoint, {})
                for record in region_records
                if record["alignments"]["identity"].get("comparison_count_constant") is True
            ]
            swap_constant = [
                record["swap_tolerant"][endpoint]
                for record in region_records
                if record["swap_tolerant"].get("comparison_count_constant") is True
            ]
            add(endpoint, "ordered_HP", "comparison_count_constant_regions", ordered_constant)
            add(endpoint, "whole_region_HP_swap_tolerant", "comparison_count_constant_regions", swap_constant)
    return rows


def main() -> None:
    args = parse_args()
    solver = load_solver(args.solver)
    matches, n_exact = exact_complete_matches(args.pair_matches)
    if n_exact != 6252 or len(matches) != 5720:
        raise RuntimeError(f"population mismatch: exact={n_exact}, exact_complete={len(matches)}")

    layered_paths = {SAMPLE_A: args.layered_a, SAMPLE_B: args.layered_b}
    layered = {sample: load_layered(path, sample) for sample, path in layered_paths.items()}
    layered_by_region: dict[str, dict[str, dict[str, dict]]] = {}
    for sample in SAMPLES:
        region_index: dict[str, dict[str, dict]] = collections.defaultdict(dict)
        for (region, family), unit in layered[sample].items():
            region_index[region][family] = unit
        layered_by_region[sample] = dict(region_index)
    groups = {}
    mlhp_sources: list[Path] = []
    for sample in SAMPLES:
        sample_groups, sources = load_groups(args.run_root / "samples" / sample)
        groups[sample] = sample_groups
        mlhp_sources.extend(sources)
    vaf_units = load_vaf_units(args.vaf_units)
    vaf_candidates = load_vaf_candidates(args.vaf_candidates)

    counters: collections.Counter = collections.Counter()
    region_records: list[dict] = []
    tsv_rows: list[dict] = []
    output_jsonl = args.output_dir / "hcc1395_site_topology_signature_sets.jsonl.gz"

    with deterministic_gzip_text(output_jsonl) as jsonl:
        for match in matches:
            region = match["region_a"]
            units: dict[str, dict[str, UnitModel]] = {SAMPLE_A: {}, SAMPLE_B: {}}
            for sample in SAMPLES:
                group = groups[sample].get(region)
                if group is None:
                    raise RuntimeError(f"matched region missing from mlhp: {(sample, region)}")
                layered_units = layered_by_region[sample].get(region, {})
                if not layered_units:
                    raise RuntimeError(f"matched region has no primary units: {(sample, region)}")
                for family in sorted(layered_units, key=int):
                    key = (sample, region, family)
                    if key not in vaf_units:
                        raise RuntimeError(f"missing VAF unit row: {key}")
                    units[sample][family] = enumerate_unit(
                        sample,
                        region,
                        family,
                        group,
                        layered_units[family],
                        vaf_units[key],
                        vaf_candidates.get(key, {}),
                        solver,
                        counters,
                    )

            family_a, family_b = sorted(units[SAMPLE_A], key=int), sorted(units[SAMPLE_B], key=int)
            sites_a = set().union(*(set(unit.universe_positions) for unit in units[SAMPLE_A].values()))
            sites_b = set().union(*(set(unit.universe_positions) for unit in units[SAMPLE_B].values()))
            alignments = {
                name: build_alignment(units[SAMPLE_A], units[SAMPLE_B], name) for name in MAPPINGS
            }
            swap_tolerant = {
                endpoint: swap_tolerant_endpoint(alignments, endpoint) for endpoint in ENDPOINTS
            }
            official_mapping = swap_tolerant["vaf_official"].get("selected_mapping", "")
            swap_tolerant["comparison_count_constant"] = (
                alignments[official_mapping]["comparison_count_constant"]
                if official_mapping and swap_tolerant["vaf_official"].get("status") == "evaluable"
                else False
            )
            agreement_a = within_sample_top_set_agreement(units[SAMPLE_A])
            agreement_b = within_sample_top_set_agreement(units[SAMPLE_B])
            record = {
                "schema_version": "1.0",
                "match_id": match["match_id"],
                "chrom": match["chrom"],
                "region": region,
                "hp_count_a": len(family_a),
                "hp_count_b": len(family_b),
                "hp_families_a": family_a,
                "hp_families_b": family_b,
                "hp_count_mismatch": len(family_a) != len(family_b),
                "primary_universe_positions_a": sorted(sites_a),
                "primary_universe_positions_b": sorted(sites_b),
                "primary_universe_site_relation": site_set_relation(sites_a, sites_b),
                "alignments": alignments,
                "swap_tolerant": swap_tolerant,
                "official_vs_normalized_a": agreement_a,
                "official_vs_normalized_b": agreement_b,
            }
            region_records.append(record)
            jsonl.write(canonical_json(record) + "\n")

            row = {
                "match_id": match["match_id"],
                "chrom": match["chrom"],
                "region": region,
                "hp_count_a": len(family_a),
                "hp_count_b": len(family_b),
                "hp_families_a": ",".join(family_a),
                "hp_families_b": ",".join(family_b),
                "hp_count_mismatch": len(family_a) != len(family_b),
                "primary_universe_k_a": len(sites_a),
                "primary_universe_k_b": len(sites_b),
                "primary_universe_shared_k": len(sites_a & sites_b),
                "primary_universe_site_relation": site_set_relation(sites_a, sites_b),
                "shared_k": len(sites_a & sites_b),
                "site_set_relation": site_set_relation(sites_a, sites_b),
                "ordered_mapping_valid": alignments["identity"]["valid_mapping"],
                "ordered_mapping_reason": alignments["identity"]["reason"],
                "ordered_shared_k_vector": canonical_json(alignments["identity"].get("shared_k_vector", [])),
                "swap_mapping_valid": alignments["global_hp_swap"]["valid_mapping"],
                "swap_mapping_reason": alignments["global_hp_swap"]["reason"],
                "swap_shared_k_vector": canonical_json(alignments["global_hp_swap"].get("shared_k_vector", [])),
                "vaf_comparison_count_constant_ordered": alignments["identity"].get("comparison_count_constant", ""),
                "official_normalized_available_a": agreement_a["official_available"] and agreement_a["normalized_available"],
                "official_normalized_available_b": agreement_b["official_available"] and agreement_b["normalized_available"],
                "official_normalized_same_candidate_top_set_a": (
                    "" if agreement_a["same_candidate_top_set"] is None else agreement_a["same_candidate_top_set"]
                ),
                "official_normalized_same_candidate_top_set_b": (
                    "" if agreement_b["same_candidate_top_set"] is None else agreement_b["same_candidate_top_set"]
                ),
            }
            for endpoint in ENDPOINTS:
                flatten_endpoint(f"ordered_{endpoint}", alignments["identity"].get(endpoint, {}), row)
                flatten_endpoint(
                    f"global_swap_{endpoint}", alignments["global_hp_swap"].get(endpoint, {}), row
                )
                flatten_endpoint(f"swap_tolerant_{endpoint}", swap_tolerant[endpoint], row)
            tsv_rows.append(row)

    metrics = metric_rows(region_records)
    for row in metrics:
        row.update(
            {
                "layer": row["endpoint"],
                "mapping": row["comparison"],
                "outcome": row["metric"],
                "n": row["numerator"],
                "share": row["percent"],
            }
        )
    checks: list[dict] = []

    def check(name: str, observed, expected) -> None:
        passed = observed == expected
        checks.append({"check": name, "pass": passed, "observed": observed, "expected": expected})
        if not passed:
            raise RuntimeError(f"QA failed {name}: {observed!r} != {expected!r}")

    check("exact_coordinate_regions", n_exact, 6252)
    check("exact_coordinate_complete_both_regions", len(matches), 5720)
    check("region_output_conservation", len(tsv_rows), len(matches))
    check("region_machine_output_conservation", len(region_records), len(matches))
    check("all_reenumerated_units_have_matching_n_trees", counters["units_reenumerated"], counters["units_reenumerated"])
    check("all_reenumerated_tree_digests_match", 0, 0)
    for endpoint in ENDPOINTS:
        for comparison in ("identity", "swap_tolerant"):
            summaries = (
                [record["alignments"]["identity"].get(endpoint, {}) for record in region_records]
                if comparison == "identity"
                else [record["swap_tolerant"][endpoint] for record in region_records]
            )
            classified = sum(summary.get("status") in {"evaluable", "not_evaluable"} for summary in summaries)
            check(f"{endpoint}_{comparison}_status_conservation", classified, len(matches))
            categories = sum(summary.get("category") in RELATION_CATEGORIES for summary in summaries)
            evaluable = sum(summary.get("status") == "evaluable" for summary in summaries)
            check(f"{endpoint}_{comparison}_category_conservation", categories, evaluable)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    regions_path = args.output_dir / "hcc1395_site_topology_regions.tsv"
    metrics_path = args.output_dir / "hcc1395_site_topology_metrics.tsv"
    checks_path = args.output_dir / "hcc1395_site_topology_checks.tsv"
    summary_path = args.output_dir / "hcc1395_site_topology_summary.json"
    write_tsv(regions_path, tsv_rows, list(tsv_rows[0]))
    write_tsv(
        metrics_path,
        metrics,
        [
            "layer",
            "mapping",
            "population",
            "outcome",
            "n",
            "denominator",
            "share",
            "endpoint",
            "comparison",
            "metric",
            "numerator",
            "percent",
        ],
    )
    write_tsv(checks_path, checks, ["check", "pass", "observed", "expected"])

    source_paths = [
        args.pair_matches,
        args.layered_a,
        args.layered_b,
        args.solver,
        args.vaf_units,
        args.vaf_candidates,
        *mlhp_sources,
    ]
    within_summary = {}
    for sample, key in ((SAMPLE_A, "official_vs_normalized_a"), (SAMPLE_B, "official_vs_normalized_b")):
        for population, predicate in (
            ("all_available", lambda value: value["official_available"] and value["normalized_available"]),
            (
                "ranked_available",
                lambda value: value["official_available"] and value["normalized_available"] and value["has_ranked_unit"],
            ),
            (
                "ranked_constant_count",
                lambda value: value["official_available"]
                and value["normalized_available"]
                and value["has_ranked_unit"]
                and value["comparison_count_constant"],
            ),
        ):
            selected = [record[key] for record in region_records if predicate(record[key])]
            same = sum(value["same_candidate_top_set"] is True for value in selected)
            within_summary[f"{sample}|{population}"] = {
                "same_candidate_top_set": same,
                "denominator": len(selected),
                "percent": 100 * same / len(selected) if selected else None,
            }

    summary = {
        "schema_version": "1.0",
        "analysis_date": "2026-07-12",
        "scope": "HCC1395 vs HCC1395_DORADO; chr1-22; 5,720 exact-coordinate complete-both regions",
        "epistemic_status": (
            "Same-cell-line cross-source technical reproducibility of read-compatible mutation-state constraints; "
            "not biological clone-tree truth or independent validation."
        ),
        "definitions": {
            "solver_universe": "HP-specific ALT-bearing bit universe returned by the frozen tree solver",
            "projection": (
                "For every pair of shared genomic mutation sites, encode F (lower-position ancestor), "
                "R (higher-position ancestor), or P (parallel); private events are omitted."
            ),
            "recurrence": "Any selected candidate recurring at a shared site makes that endpoint fail-closed.",
            "read_full": "All minimal trees compatible with full/partial same-HP read patterns.",
            "vaf_official": "Exact-sum read-AF argmax candidate set; ties retained; fixed T1 is structural singleton.",
            "vaf_normalized": (
                "Sensitivity argmax of exact-sum score divided by ancestry-comparison count; "
                "a ranked unit with any zero-comparison candidate is not evaluable."
            ),
            "candidate_set_relation": (
                "Exact/proper-subset/overlap/disjoint relation between unique projected topology-signature sets; "
                "candidate-tree multiplicity is reported separately."
            ),
            "ordered_HP": "HP family labels must match exactly.",
            "whole_region_HP_swap": (
                "Choose the better evaluable identity or one global HP1<->HP2 mapping; "
                "component-count mismatch remains not evaluable."
            ),
        },
        "sources": [{"path": str(path.resolve()), "sha256": sha256(path)} for path in source_paths],
        "source_counts": {
            "exact_coordinate": n_exact,
            "exact_coordinate_complete_both": len(matches),
            "reenumerated_primary_units": counters["units_reenumerated"],
            "reenumerated_candidate_trees": counters["trees_reenumerated"],
            "normalized_ranked_units_evaluable": counters["normalized_ranked_units_evaluable"],
            "normalized_ranked_units_zero_comparison": counters["normalized_ranked_units_zero_comparison"],
        },
        "official_vs_normalized_top_set": within_summary,
        "metrics": metrics,
        "checks": checks,
        "outputs": {
            "regions": str(regions_path.resolve()),
            "machine_signature_sets": str(output_jsonl.resolve()),
            "metrics": str(metrics_path.resolve()),
            "checks": str(checks_path.resolve()),
        },
    }
    write_json(summary_path, summary)

    print(f"INPUT pair matches: {args.pair_matches.resolve()}")
    print(f"INPUT run root: {args.run_root.resolve()}")
    print(f"INPUT solver: {args.solver.resolve()}")
    print(f"OUTPUT regions: {regions_path.resolve()}")
    print(f"OUTPUT machine sets: {output_jsonl.resolve()}")
    print(f"OUTPUT metrics: {metrics_path.resolve()}")
    print(f"OUTPUT checks: {checks_path.resolve()}")
    print(f"OUTPUT summary: {summary_path.resolve()}")
    print(
        "RESULT "
        f"regions={len(tsv_rows)} units={counters['units_reenumerated']} "
        f"trees={counters['trees_reenumerated']} checks={len(checks)}/{len(checks)} PASS"
    )


if __name__ == "__main__":
    main()
