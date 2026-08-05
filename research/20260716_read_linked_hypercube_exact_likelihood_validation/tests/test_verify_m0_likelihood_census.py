#!/usr/bin/env python3
import pathlib
import sys
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_m0_likelihood_census import (  # noqa: E402
    categorical_conservation,
    expected_selection_status,
    sha256_path,
    summarize_rows,
    verification_provenance,
)


def make_row(
    *,
    raw_t: int,
    n_vertex: int,
    best_vertex: int,
    top_edges: int,
    status: str,
    converged: bool = True,
):
    return {
        "dataset": "COLO829",
        "region": f"chr1:{raw_t}-{n_vertex}",
        "raw_tree_candidates_T": raw_t,
        "distinct_vertex_sets_V": n_vertex,
        "best_vertex_sets": best_vertex,
        "top_edge_variants": top_edges,
        "selection_status": status,
        "all_fits_converged": converged,
        "all_fits_monotone": True,
    }


class SelectionContractTests(unittest.TestCase):
    def test_all_status_branches(self):
        cases = [
            make_row(
                raw_t=1,
                n_vertex=1,
                best_vertex=1,
                top_edges=1,
                status="T1_CANDIDATE_UNIQUE",
            ),
            make_row(
                raw_t=3,
                n_vertex=1,
                best_vertex=1,
                top_edges=3,
                status="T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED",
            ),
            make_row(
                raw_t=4,
                n_vertex=4,
                best_vertex=2,
                top_edges=2,
                status="LIKELIHOOD_TIED_VERTEX_SETS",
            ),
            make_row(
                raw_t=4,
                n_vertex=2,
                best_vertex=1,
                top_edges=3,
                status="LIKELIHOOD_UNIQUE_VERTEX_EDGE_UNRESOLVED",
            ),
            make_row(
                raw_t=4,
                n_vertex=4,
                best_vertex=1,
                top_edges=1,
                status="LIKELIHOOD_UNIQUE_VERTEX_SINGLE_EDGE_CANDIDATE",
            ),
            make_row(
                raw_t=4,
                n_vertex=4,
                best_vertex=1,
                top_edges=1,
                status="RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE",
                converged=False,
            ),
        ]
        for row in cases:
            with self.subTest(status=row["selection_status"]):
                self.assertEqual(expected_selection_status(row), row["selection_status"])

    def test_nonconvergence_overrides_apparent_tie(self):
        row = make_row(
            raw_t=4,
            n_vertex=4,
            best_vertex=2,
            top_edges=2,
            status="RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE",
            converged=False,
        )
        self.assertEqual(expected_selection_status(row), row["selection_status"])


class AggregateTests(unittest.TestCase):
    def test_mutually_exclusive_partition(self):
        rows = [
            make_row(
                raw_t=1,
                n_vertex=1,
                best_vertex=1,
                top_edges=1,
                status="T1_CANDIDATE_UNIQUE",
            ),
            make_row(
                raw_t=4,
                n_vertex=4,
                best_vertex=2,
                top_edges=2,
                status="LIKELIHOOD_TIED_VERTEX_SETS",
            ),
            make_row(
                raw_t=4,
                n_vertex=4,
                best_vertex=1,
                top_edges=1,
                status="LIKELIHOOD_UNIQUE_VERTEX_SINGLE_EDGE_CANDIDATE",
            ),
            make_row(
                raw_t=4,
                n_vertex=4,
                best_vertex=1,
                top_edges=1,
                status="RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE",
                converged=False,
            ),
        ]
        conservation = categorical_conservation(rows)
        self.assertEqual(conservation["partition_sum"], 4)
        self.assertTrue(conservation["partition_conserves"])
        self.assertTrue(conservation["T_V_partition_conserves"])

        aggregate = summarize_rows(rows)
        self.assertEqual(aggregate["n_hp_lineage_units"], 4)
        self.assertEqual(sum(aggregate["selection_status_counts"].values()), 4)
        self.assertEqual(aggregate["n_optimizer_abstain"], 1)
        self.assertFalse(aggregate["all_fits_converged"])
        self.assertTrue(aggregate["all_fits_monotone"])


class ProvenanceBindingTests(unittest.TestCase):
    def test_required_byte_hashes_and_runtime_versions_are_bound(self):
        fixture = pathlib.Path(__file__).resolve()
        provenance = verification_provenance(fixture, fixture)
        byte_identity = provenance["byte_identity"]
        self.assertEqual(byte_identity["m0_receipt"]["sha256"], sha256_path(fixture))
        self.assertEqual(byte_identity["m0_rows_gzip"]["sha256"], sha256_path(fixture))
        self.assertIn("gzip container", byte_identity["m0_rows_gzip"]["hash_scope"])
        self.assertEqual(
            pathlib.Path(byte_identity["m0_census_script_on_disk_at_verification"]["path"]).name,
            "run_m0_likelihood_census.py",
        )
        self.assertEqual(
            pathlib.Path(byte_identity["independent_verifier"]["path"]).name,
            "verify_m0_likelihood_census.py",
        )
        runtime = provenance["verification_runtime"]
        self.assertTrue(runtime["python_version"])
        self.assertTrue(runtime["numpy_version"])
        self.assertTrue(runtime["scipy_version"])


if __name__ == "__main__":
    unittest.main()
