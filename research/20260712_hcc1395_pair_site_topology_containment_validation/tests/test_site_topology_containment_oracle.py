#!/usr/bin/env python3
"""Independent golden oracle for site-aware topology containment.

This test module intentionally does not import the production comparison
script.  It freezes the scientific contract in implementation-neutral terms:
mutation labels, pairwise ancestry/parallel relations, complete candidate
sets, and one whole-region HP mapping.  A later production adapter should be
tested against the same JSON vectors.
"""

from __future__ import annotations

import itertools
import json
import unittest
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable


GOLDEN_PATH = Path(__file__).with_name("golden_site_topology_cases.json")


class NotEvaluable(ValueError):
    """Raised when a fixture cannot support a fail-closed topology claim."""

    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason


def site_set_relation(a_sites: Iterable[str], b_sites: Iterable[str]) -> str:
    """Classify nesting of complete mutation-site sets before projection."""

    a = frozenset(a_sites)
    b = frozenset(b_sites)
    if a == b:
        return "exact"
    if a < b:
        return "a_subset_b"
    if b < a:
        return "b_subset_a"
    if a & b:
        return "partial_overlap"
    return "disjoint"


def set_relation(a_values: set[Any], b_values: set[Any]) -> str:
    """Classify two non-empty sets without converting overlap to equality."""

    if a_values == b_values:
        return "exact"
    if a_values < b_values:
        return "a_subset_b"
    if b_values < a_values:
        return "b_subset_a"
    if a_values & b_values:
        return "overlap"
    return "disjoint"


def _prepare_tree(tree: dict[str, Any]) -> dict[str, Any]:
    """Validate a rooted infinite-sites tree and return indexed structures."""

    if tree.get("recurrence"):
        raise NotEvaluable("recurrence_metadata_nonempty")

    declared_sites = frozenset(tree["sites"])
    if len(declared_sites) != len(tree["sites"]):
        raise NotEvaluable("duplicate_site_label")

    nodes: dict[str, dict[str, Any]] = {}
    for row in tree["nodes"]:
        node_id = row["id"]
        if node_id in nodes:
            raise NotEvaluable("duplicate_node_id")
        derived = frozenset(row["derived"])
        if not derived <= declared_sites:
            raise NotEvaluable("node_contains_undeclared_site")
        nodes[node_id] = {"derived": derived, "hidden": bool(row.get("hidden"))}

    parents: dict[str, str] = {}
    children: dict[str, list[str]] = defaultdict(list)
    for parent, child in tree["edges"]:
        if parent not in nodes or child not in nodes:
            raise NotEvaluable("edge_endpoint_missing")
        if child in parents:
            raise NotEvaluable("multiple_parents")
        parent_sites = nodes[parent]["derived"]
        child_sites = nodes[child]["derived"]
        if not parent_sites < child_sites:
            # Loss of a mutation is a recurrence/back-mutation violation; an
            # equal-state edge is also invalid in the unprojected input tree.
            raise NotEvaluable("nonmonotone_edge")
        parents[child] = parent
        children[parent].append(child)

    roots = [node_id for node_id in nodes if node_id not in parents]
    if len(roots) != 1:
        raise NotEvaluable("root_count_not_one")
    root = roots[0]
    if nodes[root]["derived"]:
        raise NotEvaluable("root_is_not_reference_state")

    seen: set[str] = set()
    active: set[str] = set()

    def visit(node_id: str) -> None:
        if node_id in active:
            raise NotEvaluable("cycle")
        if node_id in seen:
            return
        active.add(node_id)
        for child_id in children[node_id]:
            visit(child_id)
        active.remove(node_id)
        seen.add(node_id)

    visit(root)
    if seen != set(nodes):
        raise NotEvaluable("disconnected_tree")

    return {
        "sites": declared_sites,
        "nodes": nodes,
        "parents": parents,
        "children": children,
        "root": root,
        "edges": [tuple(edge) for edge in tree["edges"]],
    }


def topology_signature(
    tree: dict[str, Any], shared_sites: Iterable[str], *, ignore_hidden: bool
) -> tuple[Any, ...]:
    """Encode labeled pairwise mutation relations after private-site contraction.

    Each shared mutation must be introduced on exactly one edge.  Extra private
    mutations are absent from the signature, which contracts their unary path
    contribution while preserving ancestry among the retained mutations.
    """

    prepared = _prepare_tree(tree)
    shared = tuple(sorted(set(shared_sites)))
    if not set(shared) <= prepared["sites"]:
        raise NotEvaluable("shared_site_missing_from_tree")

    introduction_node: dict[str, str] = {}
    for site in shared:
        introductions = [
            child
            for parent, child in prepared["edges"]
            if site not in prepared["nodes"][parent]["derived"]
            and site in prepared["nodes"][child]["derived"]
        ]
        if not introductions:
            raise NotEvaluable("mutation_introduction_missing")
        if len(introductions) > 1:
            raise NotEvaluable("multiple_mutation_introductions")
        introduction_node[site] = introductions[0]

    def is_strict_ancestor(ancestor: str, descendant: str) -> bool:
        cursor = prepared["parents"].get(descendant)
        while cursor is not None:
            if cursor == ancestor:
                return True
            cursor = prepared["parents"].get(cursor)
        return False

    pair_relations: list[tuple[str, str, str]] = []
    for left, right in itertools.combinations(shared, 2):
        left_node = introduction_node[left]
        right_node = introduction_node[right]
        if left_node == right_node:
            relation = "same_event"
        elif is_strict_ancestor(left_node, right_node):
            relation = f"{left}>{right}"
        elif is_strict_ancestor(right_node, left_node):
            relation = f"{right}>{left}"
        else:
            relation = "parallel"
        pair_relations.append((left, right, relation))

    signature: tuple[Any, ...] = (shared, tuple(pair_relations))
    if not ignore_hidden:
        visibility = tuple(
            (site, prepared["nodes"][introduction_node[site]]["hidden"])
            for site in shared
        )
        signature += (visibility,)
    return signature


def compare_trees(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    """Compare two candidate trees with a directional site-containment gate."""

    a_sites = frozenset(a["sites"])
    b_sites = frozenset(b["sites"])
    shared = sorted(a_sites & b_sites)
    site_relation = site_set_relation(a_sites, b_sites)
    base = {"shared_sites": shared, "site_set_relation": site_relation}
    if len(shared) < 2:
        return {
            **base,
            "evaluable": False,
            "reason": "insufficient_shared_sites",
            "classification": "not_evaluable",
            "containment_eligible": False,
        }

    try:
        relaxed_a = topology_signature(a, shared, ignore_hidden=True)
        relaxed_b = topology_signature(b, shared, ignore_hidden=True)
        strict_a = topology_signature(a, shared, ignore_hidden=False)
        strict_b = topology_signature(b, shared, ignore_hidden=False)
    except NotEvaluable as error:
        return {
            **base,
            "evaluable": False,
            "reason": error.reason,
            "classification": "not_evaluable",
            "containment_eligible": False,
        }

    relaxed_equal = relaxed_a == relaxed_b
    strict_equal = strict_a == strict_b
    if not relaxed_equal:
        classification = "conflict_on_shared_sites"
        containment_eligible = False
    elif site_relation == "exact":
        classification = "exact_topology" if strict_equal else "relaxed_exact_topology"
        containment_eligible = True
    elif site_relation == "a_subset_b":
        classification = "a_substructure_of_b"
        containment_eligible = True
    elif site_relation == "b_subset_a":
        classification = "b_substructure_of_a"
        containment_eligible = True
    else:
        # Equal projection on a non-nested intersection is useful sensitivity
        # evidence, but it is not a directional subtree/containment result.
        classification = "shared_projection_only"
        containment_eligible = False

    return {
        **base,
        "evaluable": True,
        "strict_equal": strict_equal,
        "relaxed_equal": relaxed_equal,
        "classification": classification,
        "containment_eligible": containment_eligible,
    }


def compare_candidate_sets(
    a_refs: list[str], b_refs: list[str], trees: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    """Compare complete candidate spaces after one shared-site projection."""

    if not a_refs or not b_refs:
        return {"evaluable": False, "reason": "empty_candidate_set"}

    a_site_sets = {frozenset(trees[ref]["sites"]) for ref in a_refs}
    b_site_sets = {frozenset(trees[ref]["sites"]) for ref in b_refs}
    if len(a_site_sets) != 1 or len(b_site_sets) != 1:
        return {"evaluable": False, "reason": "inconsistent_candidate_site_set"}
    a_sites = next(iter(a_site_sets))
    b_sites = next(iter(b_site_sets))
    shared = sorted(a_sites & b_sites)
    if len(shared) < 2:
        return {"evaluable": False, "reason": "insufficient_shared_sites"}

    try:
        a_signatures = {
            topology_signature(trees[ref], shared, ignore_hidden=True) for ref in a_refs
        }
        b_signatures = {
            topology_signature(trees[ref], shared, ignore_hidden=True) for ref in b_refs
        }
    except NotEvaluable as error:
        return {"evaluable": False, "reason": error.reason}

    relation = set_relation(a_signatures, b_signatures)
    evidence = {
        "exact": "candidate_set_equal",
        "a_subset_b": "nested_candidate_space",
        "b_subset_a": "nested_candidate_space",
        "overlap": "possible_only",
        "disjoint": "conflict",
    }[relation]
    return {
        "evaluable": True,
        "relation": relation,
        "evidence": evidence,
        "site_set_relation": site_set_relation(a_sites, b_sites),
        "n_projected_signatures_a": len(a_signatures),
        "n_projected_signatures_b": len(b_signatures),
    }


def _mapping_result(
    a_hps: dict[str, list[str]],
    b_hps: dict[str, list[str]],
    mapping: list[tuple[str, str]],
    trees: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    component_results = [
        compare_candidate_sets(a_hps[a_hp], b_hps[b_hp], trees)
        for a_hp, b_hp in mapping
    ]
    if any(not result.get("evaluable") for result in component_results):
        return {"status": "not_evaluable", "components": component_results}
    relations = [result["relation"] for result in component_results]
    if "disjoint" in relations:
        status = "conflict"
    elif "overlap" in relations:
        status = "possible_only"
    else:
        # Exact and nested complete candidate spaces are robust under the
        # candidate-set contract; existential overlap is intentionally weaker.
        status = "robust_compatible"
    return {"status": status, "components": component_results}


def compare_regions(
    a_hps: dict[str, list[str]],
    b_hps: dict[str, list[str]],
    trees: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """Compare one ordered forest and one global HP1/HP2-swapped forest."""

    if len(a_hps) != len(b_hps):
        return {
            "ordered_status": "not_evaluable",
            "swap_status": "not_evaluable",
            "whole_region_status": "not_evaluable",
            "accepted_mapping": None,
            "reason": "hp_count_mismatch",
        }
    valid_hps = {"HP1", "HP2"}
    if not set(a_hps) <= valid_hps or not set(b_hps) <= valid_hps:
        return {
            "ordered_status": "not_evaluable",
            "swap_status": "not_evaluable",
            "whole_region_status": "not_evaluable",
            "accepted_mapping": None,
            "reason": "unsupported_hp_label",
        }

    if set(a_hps) == set(b_hps):
        ordered_mapping = [(hp, hp) for hp in sorted(a_hps)]
        ordered = _mapping_result(a_hps, b_hps, ordered_mapping, trees)
    else:
        ordered = {"status": "not_evaluable", "components": []}

    flipped = {"HP1": "HP2", "HP2": "HP1"}
    if {flipped[hp] for hp in a_hps} == set(b_hps):
        swap_mapping = [(hp, flipped[hp]) for hp in sorted(a_hps)]
        swapped = _mapping_result(a_hps, b_hps, swap_mapping, trees)
    else:
        swapped = {"status": "not_applicable", "components": []}

    accepted_mapping = None
    if ordered["status"] == "robust_compatible":
        whole_status = "robust_compatible"
        accepted_mapping = "ordered"
    elif swapped["status"] == "robust_compatible":
        whole_status = "robust_compatible"
        accepted_mapping = "global_swap"
    elif ordered["status"] == "possible_only":
        whole_status = "possible_only"
        accepted_mapping = "ordered"
    elif swapped["status"] == "possible_only":
        whole_status = "possible_only"
        accepted_mapping = "global_swap"
    elif ordered["status"] == "not_evaluable" and swapped["status"] in {
        "not_applicable",
        "not_evaluable",
    }:
        whole_status = "not_evaluable"
    else:
        whole_status = "conflict"

    return {
        "ordered_status": ordered["status"],
        "swap_status": swapped["status"],
        "whole_region_status": whole_status,
        "accepted_mapping": accepted_mapping,
    }


def exact_vaf_top(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Return every exact rational-score maximizer, preserving input order."""

    if not rows:
        raise NotEvaluable("empty_vaf_candidate_set")
    parsed = [(row, Fraction(row["score"])) for row in rows]
    maximum = max(score for _, score in parsed)
    return [row for row, score in parsed if score == maximum]


class GoldenOracleTests(unittest.TestCase):
    """Execute every frozen synthetic case without production imports."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.golden = json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))
        cls.trees = cls.golden["trees"]

    def assert_expected(self, observed: dict[str, Any], expected: dict[str, Any]) -> None:
        for key, value in expected.items():
            self.assertIn(key, observed, f"oracle omitted expected field {key}")
            self.assertEqual(value, observed[key], f"field mismatch: {key}")

    def test_schema_and_claim_contract(self) -> None:
        self.assertEqual("1.0", self.golden["schema_version"])
        self.assertEqual(2, self.golden["claim_contract"]["minimum_shared_sites"])
        self.assertIn("one ordered identity mapping", self.golden["claim_contract"]["hp_swap"])

    def test_tree_golden_cases(self) -> None:
        for case in self.golden["tree_cases"]:
            with self.subTest(case=case["id"]):
                observed = compare_trees(self.trees[case["a"]], self.trees[case["b"]])
                self.assert_expected(observed, case["expected"])

    def test_candidate_set_golden_cases(self) -> None:
        for case in self.golden["candidate_set_cases"]:
            with self.subTest(case=case["id"]):
                observed = compare_candidate_sets(case["a"], case["b"], self.trees)
                self.assert_expected(observed, case["expected"])

    def test_region_golden_cases(self) -> None:
        for case in self.golden["region_cases"]:
            with self.subTest(case=case["id"]):
                observed = compare_regions(case["a"], case["b"], self.trees)
                self.assert_expected(observed, case["expected"])

    def test_cherry_pick_case_would_fool_a_non_bijective_matcher(self) -> None:
        case = next(
            row
            for row in self.golden["region_cases"]
            if row["id"] == "per_component_hp_cherry_pick_is_forbidden"
        )
        # Deliberately invalid reference: each A component may reuse any B HP.
        naive_component_results = []
        for a_refs in case["a"].values():
            naive_component_results.append(
                any(
                    compare_candidate_sets(a_refs, b_refs, self.trees).get("relation")
                    != "disjoint"
                    for b_refs in case["b"].values()
                )
            )
        self.assertTrue(all(naive_component_results))
        self.assertEqual(
            "conflict",
            compare_regions(case["a"], case["b"], self.trees)["whole_region_status"],
        )

    def test_vaf_ties_are_retained_before_topology_comparison(self) -> None:
        for case in self.golden["vaf_cases"]:
            with self.subTest(case=case["id"]):
                top_a = exact_vaf_top(case["a"])
                top_b = exact_vaf_top(case["b"])
                observed = {
                    "top_ids_a": [row["id"] for row in top_a],
                    "top_ids_b": [row["id"] for row in top_b],
                }
                relation = compare_candidate_sets(
                    [row["tree"] for row in top_a],
                    [row["tree"] for row in top_b],
                    self.trees,
                )
                observed["top_signature_relation"] = relation["relation"]

                if "naive_first_signature_relation" in case["expected"]:
                    naive = compare_candidate_sets(
                        [top_a[0]["tree"]], [top_b[0]["tree"]], self.trees
                    )
                    observed["naive_first_signature_relation"] = naive["relation"]
                self.assert_expected(observed, case["expected"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
