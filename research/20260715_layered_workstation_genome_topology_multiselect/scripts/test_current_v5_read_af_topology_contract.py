#!/usr/bin/env python3
"""Contract tests for current-v5 read-AF topology primitives."""

from __future__ import annotations

import unittest
from fractions import Fraction

from build_current_v5_read_af_topology import (
    canonical_shape,
    exact_ordering_score,
    exact_read_af,
    morphology_class,
    shape_features,
)


class CurrentV5ReadAfTopologyContractTest(unittest.TestCase):
    def test_canonical_shape_is_invariant_to_sibling_and_edge_order(self) -> None:
        edges = [
            ["ROOT", "A"],
            ["ROOT", "B"],
            ["A", "AA"],
            ["B", "BA"],
            ["B", "BB"],
        ]
        reordered = [
            ["B", "BB"],
            ["ROOT", "B"],
            ["B", "BA"],
            ["ROOT", "A"],
            ["A", "AA"],
        ]

        self.assertEqual(canonical_shape(edges), canonical_shape(reordered))

    def test_synthetic_shape_features_and_morphology_classes(self) -> None:
        cases = {
            "single": {
                "edges": [["ROOT", "A"]],
                "max_depth": 1,
                "max_outdegree": 1,
                "morphology": "single_no_within_hp_relation",
            },
            "direct": {
                "edges": [["ROOT", "A"], ["A", "AB"]],
                "max_depth": 2,
                "max_outdegree": 1,
                "morphology": "direct_chain",
            },
            "sister": {
                "edges": [["ROOT", "A"], ["ROOT", "B"]],
                "max_depth": 1,
                "max_outdegree": 2,
                "morphology": "sister_branch",
            },
            "mixed": {
                "edges": [
                    ["ROOT", "A"],
                    ["A", "AB"],
                    ["A", "AC"],
                ],
                "max_depth": 2,
                "max_outdegree": 2,
                "morphology": "direct_and_sister",
            },
        }

        for name, case in cases.items():
            with self.subTest(name=name):
                features = shape_features(case["edges"])
                self.assertEqual(features["max_depth"], case["max_depth"])
                self.assertEqual(features["max_outdegree"], case["max_outdegree"])
                self.assertEqual(
                    morphology_class(
                        features["has_direct_chain"],
                        features["has_sister_branch"],
                    ),
                    case["morphology"],
                )

    def test_exact_read_af_returns_fraction_rows_and_zero_denominator_na(self) -> None:
        values, rows, reason = exact_read_af(
            {"101": [3, 1], 202: [2, 2]},
            [101, 202],
        )

        self.assertIsNone(reason)
        self.assertEqual(values, {0: Fraction(1, 4), 1: Fraction(1, 2)})
        self.assertEqual([row["read_af"]["fraction"] for row in rows], ["1/4", "1/2"])

        values, rows, reason = exact_read_af(
            {"101": [3, 1], "202": [0, 0]},
            [101, 202],
        )

        self.assertIsNone(values)
        self.assertEqual(len(rows), 1)
        self.assertEqual(reason, "zero denominator at 202")

    def test_exact_ordering_score_uses_fraction_and_acquired_label_state(self) -> None:
        # ARR -> ARA acquires site index 2.  The expected comparison therefore
        # uses AF[0] - AF[2], not AF[0] - AF[1].  H_ must not alter the state.
        read_af = {
            0: Fraction(5, 6),
            1: Fraction(1, 2),
            2: Fraction(1, 6),
        }
        score, comparisons = exact_ordering_score(
            [["ROOT", "H_ARR"], ["H_ARR", "H_ARA"]],
            read_af,
        )

        self.assertIsInstance(score, Fraction)
        self.assertEqual(score, Fraction(2, 3))
        self.assertEqual(comparisons, [Fraction(2, 3)])

    def test_independent_hp_components_collapse_by_or_without_cross_hp_edges(self) -> None:
        components = {
            "HP1": [
                ["ROOT", "HP1:A"],
                ["HP1:A", "HP1:AB"],
            ],
            "HP2": [
                ["ROOT", "HP2:C"],
                ["ROOT", "HP2:D"],
            ],
        }
        features = {family: shape_features(edges) for family, edges in components.items()}

        region_summary = {
            "has_direct_chain": any(
                item["has_direct_chain"] for item in features.values()
            ),
            "has_sister_branch": any(
                item["has_sister_branch"] for item in features.values()
            ),
        }
        region_summary["morphology_class"] = morphology_class(
            region_summary["has_direct_chain"],
            region_summary["has_sister_branch"],
        )

        self.assertTrue(features["HP1"]["has_direct_chain"])
        self.assertFalse(features["HP1"]["has_sister_branch"])
        self.assertFalse(features["HP2"]["has_direct_chain"])
        self.assertTrue(features["HP2"]["has_sister_branch"])
        self.assertEqual(region_summary["morphology_class"], "direct_and_sister")

        # Region aggregation contains only OR-reduced features.  It deliberately
        # has no combined edge list from which a cross-HP ancestry could arise.
        self.assertNotIn("edges", region_summary)
        for family, edges in components.items():
            for parent, child in edges:
                self.assertTrue(parent == "ROOT" or parent.startswith(family + ":"))
                self.assertTrue(child.startswith(family + ":"))


if __name__ == "__main__":
    unittest.main()
