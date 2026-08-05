#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import sys
import unittest


HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))

from methyl_edge_probe import (  # noqa: E402
    evaluate_edge_preference,
    make_reads,
    read_distance,
    synthetic_cases,
)


class MethylEdgeProbeTest(unittest.TestCase):
    def test_all_synthetic_controls(self):
        rows = synthetic_cases()
        self.assertEqual(7, len(rows))
        self.assertTrue(all(row["pass"] for row in rows), rows)

    def test_direction_reverses_when_10_and_01_profiles_swap(self):
        low = [0.1] * 12
        high = [0.9] * 12
        child = [0.15] * 12
        forward = evaluate_edge_preference(
            make_reads("10", low),
            make_reads("01", high),
            make_reads("11", child),
            min_reads_per_state=8,
            min_common_cpg=5,
            bootstrap_replicates=100,
        )
        reverse = evaluate_edge_preference(
            make_reads("10", high),
            make_reads("01", low),
            make_reads("11", child),
            min_reads_per_state=8,
            min_common_cpg=5,
            bootstrap_replicates=100,
        )
        self.assertEqual("MODEL_FAVORED_10_TO_11", forward.status)
        self.assertEqual("MODEL_FAVORED_01_TO_11", reverse.status)
        self.assertAlmostEqual(forward.delta_m, -reverse.delta_m, places=12)

    def test_uninformative_half_probability_is_invalid(self):
        left = make_reads("a", [0.5] * 10)[0]
        right = make_reads("b", [0.5] * 10)[0]
        self.assertIsNone(read_distance(left, right, min_common_cpg=5))


if __name__ == "__main__":
    unittest.main()
