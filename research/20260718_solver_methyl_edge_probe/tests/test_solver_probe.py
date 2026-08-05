#!/usr/bin/env python3
"""Regression tests for the isolated exact-solver probe."""

from __future__ import annotations

import itertools
import os
import pathlib
import random
import sys
import unittest


HERE = pathlib.Path(__file__).resolve().parent
PROBE_ROOT = HERE.parent
sys.path.insert(0, str(PROBE_ROOT / "scripts"))

from solver_probe import (  # noqa: E402
    ExactObligationBnb,
    brute_force_optimal,
    build_instance,
    current_scipy_enumerate,
    persistent_highs_enumerate,
    vertex_set_digest,
)


CURRENT_SOLVER = pathlib.Path(
    os.environ.get(
        "INTERSUBMOD_CURRENT_SOLVER",
        PROBE_ROOT.parent
        / "20260716_read_linked_hypercube_exact_likelihood_validation"
        / "scripts"
        / "hypercube_exact.py",
    )
)
HIGHSPY_PATH = pathlib.Path(
    os.environ.get("INTERSUBMOD_HIGHSPY_PATH", "/tmp/intersubmod_highspy_1_15_1")
)


def digests(k: int, vertex_sets) -> list[str]:
    return sorted(vertex_set_digest(k, values) for values in vertex_sets)


class SolverProbeTest(unittest.TestCase):
    def test_declared_small_oracles(self):
        cases = [
            ([], [], 2, 0, 1),
            (["AA"], [], 2, 1, 2),
            ([], ["AX", "XA"], 2, 2, 3),
            (["AAA"], [], 3, 2, 6),
            (["RRA", "AAA"], ["RAX"], 3, 1, 1),
            ([], ["AAA"], 3, 3, 6),
            (["AR", "RA", "AA"], [], 2, 0, 1),
            (["AAAA"], [], 4, 3, 24),
        ]
        for full, partial, k, expected_h, expected_sets in cases:
            with self.subTest(full=full, partial=partial, k=k):
                instance = build_instance(full, partial, k)
                brute = brute_force_optimal(instance)
                bnb = ExactObligationBnb(instance).run()
                self.assertTrue(brute["complete"])
                self.assertTrue(bnb.complete)
                self.assertEqual(expected_h, brute["objective"])
                self.assertEqual(expected_h, bnb.objective)
                self.assertEqual(expected_sets, len(brute["vertex_sets"]))
                self.assertEqual(
                    digests(k, brute["vertex_sets"]),
                    digests(k, bnb.vertex_sets),
                )

    def test_seeded_bnb_matches_independent_bruteforce(self):
        rng = random.Random(20260718)
        checked = 0
        for _ in range(80):
            k = rng.randint(1, 4)
            full_pool = [
                "".join(symbols)
                for symbols in itertools.product("RA", repeat=k)
            ]
            partial_pool = [
                "".join(symbols)
                for symbols in itertools.product("RAX", repeat=k)
                if "X" in symbols
            ]
            full = rng.sample(full_pool, rng.randint(0, min(2, len(full_pool))))
            partial = rng.sample(
                partial_pool,
                rng.randint(0, min(3, len(partial_pool))),
            )
            instance = build_instance(full, partial, k)
            brute = brute_force_optimal(instance, combination_cap=300_000)
            if not brute["complete"]:
                continue
            bnb = ExactObligationBnb(instance, time_limit_seconds=10).run()
            self.assertTrue(bnb.complete, (full, partial, k, bnb.status))
            self.assertEqual(brute["objective"], bnb.objective)
            self.assertEqual(
                digests(k, brute["vertex_sets"]),
                digests(k, bnb.vertex_sets),
            )
            checked += 1
        self.assertGreaterEqual(checked, 60)

    @unittest.skipUnless(HIGHSPY_PATH.exists(), "isolated highspy install not present")
    def test_persistent_highs_matches_bruteforce_and_builds_once(self):
        instance = build_instance(["AAAA"], [], 4)
        brute = brute_force_optimal(instance)
        result = persistent_highs_enumerate(
            instance,
            highspy_path=str(HIGHSPY_PATH),
            max_sets=100,
        )
        self.assertTrue(result["complete"])
        self.assertEqual(3, result["objective"])
        self.assertEqual(24, len(result["vertex_sets"]))
        self.assertEqual(1, result["model_builds"])
        self.assertEqual(25, result["solve_calls"])
        self.assertEqual(
            digests(4, brute["vertex_sets"]),
            digests(4, result["vertex_sets"]),
        )

    @unittest.skipUnless(
        HIGHSPY_PATH.exists() and CURRENT_SOLVER.exists(),
        "current solver or isolated highspy install not present",
    )
    def test_current_and_persistent_match_after_prepared_base_remediation(self):
        instance = build_instance(["AAAA"], [], 4)
        current = current_scipy_enumerate(
            str(CURRENT_SOLVER),
            ["AAAA"],
            [],
            4,
            max_sets=100,
            time_limit_seconds=10,
        )
        persistent = persistent_highs_enumerate(
            instance,
            highspy_path=str(HIGHSPY_PATH),
            max_sets=100,
        )
        self.assertTrue(current["complete"])
        self.assertTrue(persistent["complete"])
        # The historical probe observed 25 builds.  The current solver now
        # prepares the immutable base once; build count is diagnostic only.
        self.assertEqual(1, current["model_builds"])
        self.assertEqual(1, persistent["model_builds"])
        self.assertEqual(
            digests(4, current["vertex_sets"]),
            digests(4, persistent["vertex_sets"]),
        )


if __name__ == "__main__":
    unittest.main()
