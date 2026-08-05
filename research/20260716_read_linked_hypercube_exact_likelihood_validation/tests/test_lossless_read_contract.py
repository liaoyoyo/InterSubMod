#!/usr/bin/env python3
import pathlib
import random
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from lossless_read_contract import (  # noqa: E402
    EncodedPattern,
    SiteCall,
    aggregate_symbolic_patterns,
    component_digest,
    cut_bridge_support,
    linkage_components,
    molecule_id,
    molecule_patterns,
)


def call(molecule, pos, state, hp="1", ps="100"):
    return SiteCall("D", "chr1", pos, molecule, state, hp, ps)


class PatternContractTests(unittest.TestCase):
    def test_extended_reasons_preserved_but_solver_pattern_uses_x(self):
        encoded = EncodedPattern.from_codes(list("RAODSLX"))
        self.assertEqual(encoded.pattern, "RAXXXXX")
        self.assertEqual(encoded.ref_mask, 1)
        self.assertEqual(encoded.alt_mask, 2)
        self.assertEqual(encoded.n_free, 5)
        self.assertEqual(encoded.n_completions, 32)

    def test_all_x_is_one_symbolic_group_with_all_completions(self):
        encoded = EncodedPattern.from_codes(list("XXXX"))
        self.assertEqual(encoded.pattern, "XXXX")
        self.assertEqual(encoded.fixed_mask, 0)
        self.assertEqual(encoded.free_mask, 0b1111)
        self.assertEqual(encoded.n_completions, 16)
        self.assertTrue(all(encoded.compatible(vertex) for vertex in range(16)))

    def test_rax_has_two_completions_and_symbolic_compatibility(self):
        encoded = EncodedPattern.from_codes(list("RAX"))
        self.assertEqual(encoded.n_completions, 2)
        self.assertEqual({vertex for vertex in range(8) if encoded.compatible(vertex)}, {2, 6})

    def test_missing_site_rows_become_x_without_duplicate_molecules(self):
        rows = [call("m1", 10, "R"), call("m1", 30, "A"), call("m2", 20, "L")]
        patterns = molecule_patterns(rows, [10, 20, 30])
        self.assertEqual(patterns[("m1", "1", "100")].pattern, "RXA")
        self.assertEqual(patterns[("m2", "1", "100")].pattern, "XXX")
        aggregated = aggregate_symbolic_patterns(rows, [10, 20, 30])
        self.assertEqual(sum(row["n_molecules"] for row in aggregated), 2)

    def test_molecule_identity_is_dataset_scoped(self):
        self.assertEqual(molecule_id("D", None, "q"), molecule_id("D", ".", "q"))
        self.assertNotEqual(molecule_id("D1", None, "q"), molecule_id("D2", None, "q"))


class ReadLinkedBoundaryTests(unittest.TestCase):
    def test_abc_plus_cd_stitches_abcd_at_threshold_one(self):
        positions = [10, 20, 30, 40]
        rows = [
            call("read_abc", 10, "A"),
            call("read_abc", 20, "R"),
            call("read_abc", 30, "A"),
            call("read_cd", 30, "R"),
            call("read_cd", 40, "A"),
        ]
        support = cut_bridge_support(rows, positions)
        self.assertEqual(support, (1, 1, 1))
        self.assertEqual(linkage_components(positions, support, 1), ((10, 20, 30, 40),))
        self.assertEqual(linkage_components(positions, support, 2), ((10,), (20,), (30,), (40,)))

    def test_each_molecule_counts_once_per_cut(self):
        positions = [10, 20, 30]
        rows = [call("m1", 10, "R"), call("m1", 20, "A"), call("m1", 30, "A")]
        self.assertEqual(cut_bridge_support(rows, positions), (1, 1))

    def test_order_invariance(self):
        positions = [10, 20, 30]
        rows = [call("m1", 10, "R"), call("m1", 30, "A"), call("m2", 20, "R"), call("m2", 30, "A")]
        first = cut_bridge_support(rows, positions)
        random.Random(7).shuffle(rows)
        second = cut_bridge_support(rows, positions)
        self.assertEqual(first, second)
        self.assertEqual(component_digest(positions, first, 1), component_digest(positions, second, 1))

    def test_hp_specific_support(self):
        positions = [10, 20]
        rows = [call("m1", 10, "R", "1"), call("m1", 20, "A", "1"), call("m2", 10, "R", "2"), call("m2", 20, "A", "2")]
        self.assertEqual(cut_bridge_support(rows, positions), (2,))
        self.assertEqual(cut_bridge_support(rows, positions, hp_family="1"), (1,))

    def test_pooled_linkage_does_not_imply_within_hp_linkage(self):
        positions = [10, 20, 30]
        rows = [
            call("hp1_ab", 10, "R", "1"),
            call("hp1_ab", 20, "A", "1"),
            call("hp2_bc", 20, "R", "2"),
            call("hp2_bc", 30, "A", "2"),
        ]
        pooled = cut_bridge_support(rows, positions)
        hp1 = cut_bridge_support(rows, positions, hp_family="1")
        hp2 = cut_bridge_support(rows, positions, hp_family="2")
        self.assertEqual(pooled, (1, 1))
        self.assertEqual(hp1, (1, 0))
        self.assertEqual(hp2, (0, 1))
        self.assertEqual(linkage_components(positions, pooled, 1), ((10, 20, 30),))
        self.assertEqual(linkage_components(positions, hp1, 1), ((10, 20), (30,)))
        self.assertEqual(linkage_components(positions, hp2, 1), ((10,), (20, 30)))

    def test_components_reject_unsorted_or_duplicate_positions(self):
        with self.assertRaisesRegex(ValueError, "sorted and unique"):
            linkage_components([20, 10], [1], 1)
        with self.assertRaisesRegex(ValueError, "sorted and unique"):
            linkage_components([10, 10], [1], 1)


if __name__ == "__main__":
    unittest.main()
