#!/usr/bin/env python3
"""Unit and regression-contract tests for HCC exact-signature validation."""

from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("analyze_hcc_exact_signature.py")
SPEC = importlib.util.spec_from_file_location("hcc_exact_signature", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
ARTIFACT_PATH = MODULE_PATH.parents[1] / "artifacts" / "hcc1395_exact_signature_validation.json"


class SignatureUnitTest(unittest.TestCase):
    def test_shape_is_invariant_to_sibling_and_edge_order(self) -> None:
        first = [["ROOT", "H_AR"], ["ROOT", "H_RA"], ["H_AR", "AA"]]
        second = [["ROOT", "H_RA"], ["H_AR", "AA"], ["ROOT", "H_AR"]]
        self.assertEqual(MODULE.canonical_shape(first), MODULE.canonical_shape(second))

    def test_exact_edges_are_order_invariant_but_label_sensitive(self) -> None:
        first = [["ROOT", "H_AR"], ["H_AR", "AA"]]
        reordered = [["H_AR", "AA"], ["ROOT", "H_AR"]]
        relabeled = [["ROOT", "H_RA"], ["H_RA", "AA"]]
        self.assertEqual(MODULE.canonical_edges(first), MODULE.canonical_edges(reordered))
        self.assertNotEqual(MODULE.canonical_edges(first), MODULE.canonical_edges(relabeled))

    def test_primary_unit_multiset_ignores_hp_labels_and_preserves_count(self) -> None:
        hp1_hp2 = MODULE.region_multiset(["((()))", "(())"])
        hp2_hp1 = MODULE.region_multiset(["(())", "((()))"])
        one_unit = MODULE.region_multiset(["((()))"])
        self.assertEqual(hp1_hp2, hp2_hp1)
        self.assertNotEqual(hp1_hp2, one_unit)

    def test_fixed_unit_is_supplemented_from_source_tree(self) -> None:
        unit = {
            "analysis_candidate_set_complete": True,
            "capped": False,
            "verification_status": "full_pass",
            "verification_complete": True,
            "verify_pass": True,
            "analysis_trees_generated": 1,
            "n_trees": 1,
            "n_distinct_shapes_exact": 1,
            "trees": [{"edges": [["ROOT", "H_A"]]}],
        }
        result = MODULE.derive_unit_signatures(unit, None)
        self.assertEqual(result["shape_source"], "layered_reconstruction_fixed_unit")
        self.assertEqual(result["edge_source"], "layered_reconstruction_fixed_unit")
        self.assertIsNotNone(result["shape"])
        self.assertIsNotNone(result["exact_labeled_edges"])

    def test_tied_same_shape_is_not_an_exact_candidate(self) -> None:
        unit = {
            "analysis_candidate_set_complete": True,
            "capped": False,
            "verification_status": "full_pass",
            "verification_complete": True,
            "verify_pass": True,
            "analysis_trees_generated": 2,
            "n_trees": 2,
            "n_distinct_shapes_exact": 2,
            "trees": [{"edges": [["ROOT", "H_A"]]}, {"edges": [["ROOT", "H_R"]]}],
        }
        ranking = {
            "status": "ranked",
            "n_top_exact": 2,
            "n_top_shapes_exact": 1,
            "top_shape_representatives": [
                {"shape_signature": "(())", "edges": [["ROOT", "H_A"]]}
            ],
        }
        result = MODULE.derive_unit_signatures(unit, ranking)
        self.assertEqual(result["shape"], "(())")
        self.assertIsNone(result["exact_labeled_edges"])

    def test_cli_requires_all_provenance_inputs(self) -> None:
        parser = MODULE.build_parser()
        required = {
            action.dest
            for action in parser._actions
            if getattr(action, "required", False)
        }
        self.assertEqual(
            required,
            {
                "hcc_sidecar",
                "dorado_sidecar",
                "sidecar_index",
                "current_summary",
                "run_root",
                "input_manifest",
                "output",
            },
        )


class CurrentArtifactContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = json.loads(ARTIFACT_PATH.read_text(encoding="utf-8"))

    def test_artifact_schema_and_checks(self) -> None:
        self.assertEqual(
            self.artifact["schema_name"],
            "intersubmod.hcc1395_exact_signature_validation",
        )
        self.assertTrue(self.artifact["all_checks_pass"])
        self.assertTrue(all(self.artifact["checks"].values()))
        for sample in ("HCC1395", "HCC1395_DORADO"):
            self.assertTrue(self.artifact["sample_contracts"][sample]["all_checks_pass"])

    def test_current_v5_regression_counts(self) -> None:
        self.assertEqual(self.artifact["region_universe"]["intersection"], 6998)
        self.assertEqual(self.artifact["region_universe"]["union"], 9609)
        self.assertAlmostEqual(self.artifact["region_universe"]["jaccard"], 0.7282755749817879)
        same_sites = self.artifact["internal_ssnv_pairing"]["same_set_within_exact_common_region"]
        self.assertEqual((same_sites["numerator"], same_sites["denominator"]), (6682, 6998))

        shape = self.artifact["shape_agreement"]["same_internal_ssnv_both_shape_resolved"]
        self.assertEqual((shape["numerator"], shape["denominator"]), (3681, 5039))
        shape_units = self.artifact["shape_agreement"]["same_internal_ssnv_same_primary_unit_count"]
        self.assertEqual((shape_units["numerator"], shape_units["denominator"]), (3681, 4508))

        edges = self.artifact["exact_labeled_edge_agreement"]["same_internal_ssnv_both_candidate_unique"]
        self.assertEqual((edges["numerator"], edges["denominator"]), (1582, 4482))
        edge_units = self.artifact["exact_labeled_edge_agreement"]["same_internal_ssnv_same_primary_unit_count"]
        self.assertEqual((edge_units["numerator"], edge_units["denominator"]), (1582, 4005))

    def test_fixed_supplement_and_stored_signatures_are_validated(self) -> None:
        expected = {
            "HCC1395": (4217, 5897),
            "HCC1395_DORADO": (4033, 5595),
        }
        for sample, (fixed, checked) in expected.items():
            counts = self.artifact["sample_contracts"][sample]["counts"]
            self.assertEqual(counts["fixed_units_supplemented_from_reconstruction"], fixed)
            self.assertEqual(counts["stored_shape_signatures_checked"], checked)
            self.assertEqual(counts["stored_shape_signature_mismatches"], 0)
            self.assertEqual(counts["internal_ssnv_unit_mismatches"], 0)


if __name__ == "__main__":
    unittest.main()
