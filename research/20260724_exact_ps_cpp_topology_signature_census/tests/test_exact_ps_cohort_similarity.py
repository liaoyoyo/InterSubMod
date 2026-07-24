#!/usr/bin/env python3
"""Regression tests for the exact-PS all-seven cohort similarity sidecar."""

from __future__ import annotations

import importlib.util
import io
import json
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
PRODUCER_PATH = (
    REPO_ROOT
    / "research"
    / "20260724_exact_ps_cpp_topology_signature_census"
    / "scripts"
    / "build_exact_ps_cohort_similarity.py"
)
SPEC = importlib.util.spec_from_file_location(
    "exact_ps_cohort_similarity", PRODUCER_PATH
)
assert SPEC and SPEC.loader
PRODUCER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PRODUCER)


class ExactPsCohortSimilarityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = PRODUCER.build()
        cls.profiles = cls.result["sample_profiles"]
        cls.pairs = {
            row["pair_id"]: row for row in cls.result["pairwise_profiles"]
        }
        cls.hcc = cls.result["hcc1395_dorado_strict_comparison"]

    def test_schema_scope_formula_and_source_identities(self) -> None:
        self.assertEqual(
            self.result["schema_name"],
            "intersubmod.exact_ps_cohort_similarity",
        )
        self.assertEqual(self.result["schema_version"], "1.0.0")
        self.assertEqual(self.result["task_type"], "B_COMPREHENSIVE_VALIDATION")
        self.assertEqual(
            tuple(self.result["scope"]["datasets"]),
            PRODUCER.EXPECTED_SAMPLES,
        )
        self.assertEqual(
            tuple(self.result["scope"]["chromosomes"]),
            PRODUCER.EXPECTED_CHROMS,
        )
        self.assertFalse(self.result["scope"]["partial"])
        contract = self.result["metric_contract"]
        self.assertEqual(
            contract["formula"],
            "similarity = 1 - sqrt(JS_divergence_base2)",
        )
        self.assertEqual(
            tuple(contract["dimension_order"]), PRODUCER.DIMENSION_ORDER
        )
        self.assertAlmostEqual(
            sum(contract["dimension_weights"].values()), 1.0, places=15
        )
        self.assertEqual(
            contract["dimensions"]["status"]["denominator"], "all_groups"
        )
        self.assertEqual(
            contract["dimensions"]["resolution"]["denominator"],
            "ranked_units",
        )
        identities = self.result["source_identities"]
        self.assertEqual(set(identities["samples"]), set(PRODUCER.EXPECTED_SAMPLES))
        for sample in PRODUCER.EXPECTED_SAMPLES:
            self.assertEqual(
                set(identities["samples"][sample]),
                {
                    "pipeline_receipt",
                    "mlhp",
                    "mlhp_receipt",
                    "topology",
                    "topology_receipt",
                    "census",
                },
            )
            for value in identities["samples"][sample].values():
                self.assertEqual(len(value["sha256"]), 64)
                self.assertGreater(value["size_bytes"], 0)

    def test_all_seven_profile_denominators_and_partitions(self) -> None:
        self.assertEqual(
            sum(row["denominators"]["all_groups"] for row in self.profiles.values()),
            98955,
        )
        self.assertEqual(
            sum(
                row["denominators"]["ranked_units"]
                for row in self.profiles.values()
            ),
            71955,
        )
        self.assertEqual(
            self.profiles["HCC1395"]["denominators"],
            {"all_groups": 11590, "ranked_units": 9130},
        )
        self.assertEqual(
            self.profiles["HCC1395_DORADO"]["denominators"],
            {"all_groups": 6865, "ranked_units": 5308},
        )
        for sample, profile in self.profiles.items():
            for dimension in ("status", "resolution", "morphology", "active_k"):
                record = profile["dimensions"][dimension]
                self.assertEqual(
                    sum(record["counts"].values()), record["denominator"], sample
                )
                self.assertAlmostEqual(
                    sum(record["proportions"].values()), 1.0, places=12
                )
            hp = profile["dimensions"]["hp"]
            self.assertEqual(sum(hp["counts"].values()), hp["denominator"])
            self.assertAlmostEqual(sum(hp["comparison_vector"]), 1.0, places=12)
            self.assertLessEqual(
                hp["comparison_vector"][0], hp["comparison_vector"][1]
            )

    def test_exact_21_pair_matrix_and_key_similarity_values(self) -> None:
        self.assertEqual(len(self.pairs), 21)
        self.assertEqual(
            sum(
                row["relationship"] == "same_biological_sample_technical_pair"
                for row in self.pairs.values()
            ),
            1,
        )
        self.assertEqual(
            sum(
                row["relationship"]
                == "cross_biological_sample_distribution_only"
                for row in self.pairs.values()
            ),
            20,
        )
        expected_matrix = [
            [1.000000, 0.926381, 0.825172, 0.802188, 0.677016, 0.915468, 0.876908],
            [0.926381, 1.000000, 0.806075, 0.752867, 0.628762, 0.915269, 0.900082],
            [0.825172, 0.806075, 1.000000, 0.842168, 0.706436, 0.757463, 0.759409],
            [0.802188, 0.752867, 0.842168, 1.000000, 0.840625, 0.743807, 0.706301],
            [0.677016, 0.628762, 0.706436, 0.840625, 1.000000, 0.632636, 0.591388],
            [0.915468, 0.915269, 0.757463, 0.743807, 0.632636, 1.000000, 0.878931],
            [0.876908, 0.900082, 0.759409, 0.706301, 0.591388, 0.878931, 1.000000],
        ]
        observed = self.result["matrix"]["similarity"]
        self.assertEqual(len(observed), 7)
        for row_index, row in enumerate(observed):
            for column_index, value in enumerate(row):
                self.assertAlmostEqual(
                    value, expected_matrix[row_index][column_index], places=6
                )
                self.assertAlmostEqual(
                    value, observed[column_index][row_index], places=15
                )
        hcc = self.pairs["HCC1395__HCC1395_DORADO"]
        self.assertEqual(hcc["similarity_rank_of_21"], 1)
        self.assertAlmostEqual(hcc["profile_similarity"], 0.926381219637838, places=14)
        expected_components = {
            "status": 0.943738,
            "resolution": 0.944237,
            "morphology": 0.908943,
            "active_k": 0.888063,
            "hp_symmetric": 0.998437,
        }
        for dimension, expected in expected_components.items():
            self.assertAlmostEqual(
                hcc["components"][dimension]["similarity"], expected, places=6
            )

    def test_fixed_chromosome_block_bootstrap(self) -> None:
        bootstrap = self.result[
            "hcc1395_dorado_chromosome_block_bootstrap"
        ]
        self.assertEqual(bootstrap["seed"], 20260725)
        self.assertEqual(bootstrap["replicates"], 2000)
        self.assertAlmostEqual(
            bootstrap["similarity"]["median"], 0.925010, places=6
        )
        self.assertAlmostEqual(
            bootstrap["similarity"]["p2_5"], 0.912684, places=6
        )
        self.assertAlmostEqual(
            bootstrap["similarity"]["p97_5"], 0.934393, places=6
        )
        self.assertEqual(bootstrap["rank_of_21"]["point"], 1)
        self.assertEqual(bootstrap["rank_of_21"]["p97_5"], 3)
        self.assertEqual(
            bootstrap["rank_of_21"]["rank_counts"],
            {"1": 1334, "2": 397, "3": 241, "4": 28},
        )
        self.assertAlmostEqual(
            bootstrap["rank_of_21"]["rank_1_rate"], 0.667, places=12
        )
        self.assertAlmostEqual(
            bootstrap["rank_of_21"]["top_3_rate"], 0.986, places=12
        )
        self.assertLess(
            bootstrap["gap_vs_best_other_pair"]["p2_5"], 0.0
        )
        self.assertGreater(
            bootstrap["gap_vs_best_other_pair"]["p97_5"], 0.0
        )

    def test_hcc_unambiguous_strict_concordance_regression(self) -> None:
        overlap = self.hcc["ordered_locus_key_overlap"]
        self.assertEqual(overlap["sample_a_unique_keys"], 8854)
        self.assertEqual(overlap["sample_b_unique_keys"], 4922)
        self.assertEqual(overlap["intersection"], 2440)
        self.assertEqual(overlap["union"], 11336)
        self.assertAlmostEqual(overlap["jaccard"], 0.21524347212420608)
        self.assertEqual(overlap["unambiguous_1_to_1_keys"], 1002)
        self.assertEqual(overlap["ambiguous_common_keys_not_paired"], 1438)
        self.assertEqual(overlap["same_multiplicity_keys"], 1780)

        strict = self.hcc["unambiguous_comparison"]
        expected = {
            "same_status": (912, 1002),
            "both_ranked": (897, 1002),
            "same_active_loci_given_both_ranked": (805, 897),
            "same_resolution_class": (774, 805),
            "same_shape_signature_set": (717, 805),
            "same_weighted_shape_signature_set": (690, 805),
            "same_coarse_class_set": (727, 805),
            "both_unique_tree": (564, 805),
            "same_selected_labeled_tree": (443, 564),
            "same_selected_unlabeled_shape": (485, 564),
        }
        self.assertEqual(strict["pair_count"], 1002)
        for key, (numerator, denominator) in expected.items():
            self.assertEqual(strict[key]["numerator"], numerator, key)
            self.assertEqual(strict[key]["denominator"], denominator, key)
        self.assertEqual(len(self.hcc["unambiguous_pair_records"]), 1002)
        self.assertEqual(
            len(
                {
                    (row["chrom"], tuple(row["ordered_positions"]))
                    for row in self.hcc["unambiguous_pair_records"]
                }
            ),
            1002,
        )

    def test_hcc_label_invariant_multiset_regression(self) -> None:
        multiset = self.hcc["multiset_comparison"]
        self.assertEqual(multiset["same_multiplicity_keys"], 1780)
        self.assertEqual(multiset["all_ranked_same_multiplicity_keys"], 1243)
        self.assertEqual(
            (
                multiset["same_active_locus_multiset"]["numerator"],
                multiset["same_active_locus_multiset"]["denominator"],
            ),
            (1143, 1243),
        )
        self.assertEqual(
            (
                multiset["same_shape_signature_multiset"]["numerator"],
                multiset["same_shape_signature_multiset"]["denominator"],
            ),
            (1050, 1243),
        )
        self.assertEqual(
            (
                multiset["same_weighted_shape_signature_multiset"]["numerator"],
                multiset["same_weighted_shape_signature_multiset"]["denominator"],
            ),
            (1018, 1243),
        )
        self.assertEqual(multiset["all_unique_same_multiplicity_keys"], 941)
        self.assertEqual(
            (
                multiset["same_selected_tree_multiset"]["numerator"],
                multiset["same_selected_tree_multiset"]["denominator"],
            ),
            (728, 941),
        )

    def test_repeated_hp_key_is_never_arbitrarily_paired(self) -> None:
        repeated = ("chr1", (100, 200))
        unique = ("chr1", (300, 400))
        left = {
            repeated: ["left_hp1", "left_hp2"],
            unique: ["left_unique"],
        }
        right = {
            repeated: ["right_hp1", "right_hp2"],
            unique: ["right_unique"],
        }
        unambiguous, ambiguous = PRODUCER.partition_locus_keys(left, right)
        self.assertEqual(
            unambiguous,
            [(unique, "left_unique", "right_unique")],
        )
        self.assertEqual(len(ambiguous), 1)
        self.assertEqual(ambiguous[0]["positions"], [100, 200])
        self.assertEqual(
            ambiguous[0]["reason"],
            "repeated_ordered_locus_key_no_hp_pairing",
        )
        self.assertTrue(
            self.result["checks"][
                "hcc_repeated_ordered_locus_keys_not_arbitrarily_paired"
            ]
        )

    def test_write_once_and_existing_cli_target_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="exact-ps-similarity-test-") as raw:
            root = Path(raw)
            output = root / "cohort_similarity.json"
            PRODUCER.write_new_json(output, self.result)
            loaded = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(loaded["schema_name"], PRODUCER.SCHEMA_NAME)
            with self.assertRaises(FileExistsError):
                PRODUCER.write_new_json(output, self.result)
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                self.assertEqual(PRODUCER.main(["--output", str(output)]), 2)
            self.assertIn("output already exists", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
