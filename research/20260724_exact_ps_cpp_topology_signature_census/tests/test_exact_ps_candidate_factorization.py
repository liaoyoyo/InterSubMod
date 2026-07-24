#!/usr/bin/env python3
"""Fast regression tests for the exact-PS candidate-factorization v2 bundle."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[3]
PRODUCER_PATH = (
    REPO_ROOT
    / "research"
    / "20260724_exact_ps_cpp_topology_signature_census"
    / "scripts"
    / "build_exact_ps_candidate_factorization.py"
)
SPEC = importlib.util.spec_from_file_location(
    "exact_ps_candidate_factorization", PRODUCER_PATH
)
assert SPEC and SPEC.loader
PRODUCER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PRODUCER)

V2_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/"
    "observation_workspaces/20260725_exact_ps_candidate_factorization/all7_v2"
)
RECEIPT_PATH = V2_ROOT / "receipt.json"
EXPECTED_SAMPLE_COUNTS = {
    "HCC1395": (9_130, 48_337, 48_337, 15_350),
    "HCC1395_DORADO": (5_308, 14_919, 14_919, 8_292),
    "COLO829": (10_757, 53_345, 53_345, 47_025),
    "H1437": (13_740, 190_549, 190_549, 123_734),
    "H2009": (23_128, 630_580, 630_580, 471_910),
    "HCC1937": (4_245, 15_527, 15_527, 6_674),
    "HCC1954": (5_647, 19_335, 19_335, 7_542),
}
EXPECTED_OUTPUT_IDENTITIES = {
    "HCC1395": (
        30_080_573,
        "508e4a93a6b5a69daee1312235a0f804e2e741b0c3850770f626b356b149faef",
    ),
    "HCC1395_DORADO": (
        13_267_689,
        "8738a7a77f3ae2dd6043d4b35e5f205f0e213aa6aacc85ed47a170f7ea9916ea",
    ),
    "COLO829": (
        34_218_826,
        "b33fec2c6693cc3ef446ef327a686f9fa579de7e445b41ffb1c1e9e89d4d83d2",
    ),
    "H1437": (
        83_921_963,
        "377f494d2284a18ccf0ccaecfe47e9f99c5948b79122552dea3308acb083371c",
    ),
    "H2009": (
        245_422_226,
        "002d50c41a307569b575d724f2e929e3db0aaafb179bd63298a4c912793aadc2",
    ),
    "HCC1937": (
        11_768_146,
        "a6f87f40c8a549200d0bc8e7c9df570f8b0a5c6c238907eea4cfd76ba1bc0bb2",
    ),
    "HCC1954": (
        15_121_829,
        "0fd3be9c6bac90c4a20e107e9cfab766504b6ea2a627909d65f971f9442ddef4",
    ),
}
EXPECTED_BINARY_IDENTITY = (
    377_432,
    "0fc76517c0a94987a6ec5ff1c6e4516bfbdf435f14d58429b3dfb189e5bb5819",
)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_hcc_target_and_hash(path: Path) -> tuple[dict[str, object], str]:
    """Locate the target while hashing HCC once, avoiding a second 30 MB scan."""
    digest = hashlib.sha256()
    target: dict[str, object] | None = None
    with path.open("rb") as handle:
        for raw_line in handle:
            digest.update(raw_line)
            if b'"group_index":7469' in raw_line:
                if target is not None:
                    raise AssertionError("duplicate HCC1395 group_index=7469")
                target = json.loads(raw_line)
    if target is None:
        raise AssertionError("missing HCC1395 group_index=7469")
    return target, digest.hexdigest()


class ExactPsCandidateFactorizationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))
        cls.hcc_path = Path(
            cls.receipt["samples"]["HCC1395"]["output"]["path"]
        )
        cls.target, cls.hcc_sha256 = load_hcc_target_and_hash(cls.hcc_path)

    def test_v2_receipt_scope_flags_denominators_and_totals(self) -> None:
        self.assertEqual(
            (
                self.receipt["schema_name"],
                self.receipt["schema_version"],
            ),
            PRODUCER.RECEIPT_SCHEMA,
        )
        self.assertEqual(
            self.receipt["task_type"], "B_COMPREHENSIVE_VALIDATION"
        )
        self.assertEqual(
            tuple(self.receipt["scope"]["datasets"]), PRODUCER.EXPECTED_SAMPLES
        )
        self.assertEqual(
            self.receipt["scope"]["chromosomes"],
            [f"chr{value}" for value in range(1, 23)],
        )
        self.assertTrue(self.receipt["scope"]["ranked_only"])
        self.assertTrue(self.receipt["scope"]["exact_ps_primary_hp"])
        self.assertEqual(
            self.receipt["scope"]["strict_endpoint_read_linkage_threshold"], 3
        )
        self.assertEqual(
            self.receipt["denominators"],
            {
                "all_exact_ps_hp_units": 98_955,
                "mutation_bearing_units": 85_941,
                "ranked_complete_units_included": 71_955,
                "resource_abstain_units_excluded": 10_717,
                "no_active_alt_units_excluded": 13_014,
            },
        )
        self.assertEqual(
            self.receipt["totals"],
            {
                "ranked_units": 71_955,
                "minimum_vertex_sets": 972_592,
                "minimum_trees": 972_592,
                "global_best_trees": 680_527,
            },
        )
        self.assertTrue(self.receipt["technical_all_pass"])
        self.assertTrue(self.receipt["all_pass"])
        self.assertFalse(self.receipt["validation_evidence_eligible"])
        self.assertFalse(
            self.receipt["all_mutation_bearing_families_complete"]
        )
        self.assertTrue(all(self.receipt["checks"].values()))

    def test_per_sample_exact_fixtures(self) -> None:
        self.assertEqual(
            set(self.receipt["samples"]), set(PRODUCER.EXPECTED_SAMPLES)
        )
        self.assertEqual(PRODUCER.EXPECTED_SAMPLE_COUNTS, EXPECTED_SAMPLE_COUNTS)
        fields = (
            "ranked_units",
            "minimum_vertex_sets",
            "minimum_trees",
            "global_best_trees",
        )
        for sample, expected in EXPECTED_SAMPLE_COUNTS.items():
            self.assertEqual(
                tuple(int(self.receipt["samples"][sample][key]) for key in fields),
                expected,
                sample,
            )

    def test_published_final_identity_specs_and_live_files(self) -> None:
        """Check all final specs without re-reading the whole 398 MB bundle.

        Every published file gets an exact frozen SHA fixture and a live size
        check. HCC1395 gets a live SHA check in the same scan used to locate the
        target row; the small executable also gets a live SHA check.
        """
        for sample, (expected_size, expected_sha256) in (
            EXPECTED_OUTPUT_IDENTITIES.items()
        ):
            spec = self.receipt["samples"][sample]["output"]
            path = Path(spec["path"])
            self.assertEqual(
                path,
                V2_ROOT / f"{sample}.candidate_factorization.jsonl",
                sample,
            )
            self.assertTrue(path.is_file(), path)
            self.assertEqual(path.stat().st_size, expected_size, sample)
            self.assertEqual(int(spec["size_bytes"]), expected_size, sample)
            self.assertEqual(spec["sha256"], expected_sha256, sample)

        binary = self.receipt["implementation"]["binary"]
        binary_path = Path(binary["path"])
        expected_size, expected_sha256 = EXPECTED_BINARY_IDENTITY
        self.assertEqual(
            binary_path,
            V2_ROOT
            / "implementation"
            / "exact_topology_candidate_factorization",
        )
        self.assertTrue(binary_path.is_file())
        self.assertEqual(binary_path.stat().st_size, expected_size)
        self.assertEqual(int(binary["size_bytes"]), expected_size)
        self.assertEqual(binary["sha256"], expected_sha256)
        self.assertEqual(sha256_path(binary_path), expected_sha256)
        self.assertEqual(
            self.hcc_sha256,
            EXPECTED_OUTPUT_IDENTITIES["HCC1395"][1],
        )

    def test_hcc1395_target_candidate_space_regression(self) -> None:
        self.assertEqual(self.target["group_index"], 7469)
        self.assertEqual(self.target["chrom"], "chr10")
        self.assertEqual(self.target["active_positions"], [87818272, 87840023])
        self.assertEqual(self.target["observed_vertices"], [0])
        self.assertEqual(
            [item["vertices"] for item in self.target["minimum_vertex_sets"]],
            [[0, 1, 2], [0, 1, 3], [0, 2, 3]],
        )
        self.assertEqual(self.target["minimum_vertex_set_count"], 3)
        self.assertEqual(
            self.target["all_minimum_tree_edge_incidence"],
            [[0, 1, "2"], [0, 2, "2"], [1, 3, "1"], [2, 3, "1"]],
        )
        self.assertEqual(
            self.target["global_best_edge_incidence"],
            [[0, 1, "1"], [1, 3, "1"]],
        )
        self.assertEqual(self.target["best_vertex_set_count"], 1)
        self.assertEqual(self.target["best_tree_tie_count"], "1")

    def validate_single_hcc_row(self, row: dict[str, object]) -> dict[str, object]:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "single.jsonl"
            path.write_text(
                json.dumps(row, separators=(",", ":")) + "\n",
                encoding="utf-8",
            )
            with patch.dict(
                PRODUCER.EXPECTED_SAMPLE_COUNTS,
                {"HCC1395": (1, 3, 3, 1)},
            ):
                return PRODUCER.validate_sample_output("HCC1395", path)

    def test_isolated_target_fixture_passes_generic_validator(self) -> None:
        self.assertEqual(
            self.validate_single_hcc_row(copy.deepcopy(self.target)),
            {
                "ranked_units": 1,
                "minimum_vertex_sets": 3,
                "minimum_trees": 3,
                "global_best_trees": 1,
            },
        )

    def test_validator_rejects_tampered_edge_incidence(self) -> None:
        row = copy.deepcopy(self.target)
        row["all_minimum_tree_edge_incidence"][0][2] = "3"
        with self.assertRaisesRegex(
            PRODUCER.FactorizationBuildError,
            "generic edge-incidence conservation failed",
        ):
            self.validate_single_hcc_row(row)

    def test_validator_rejects_best_parent_outside_legal_set(self) -> None:
        row = copy.deepcopy(self.target)
        row["minimum_vertex_sets"][0]["parent_factorization"][0][2] = [2]
        with self.assertRaisesRegex(
            PRODUCER.FactorizationBuildError,
            "invalid parent choices",
        ):
            self.validate_single_hcc_row(row)

    def test_validator_rejects_tampered_signature_tree_count(self) -> None:
        row = copy.deepcopy(self.target)
        row["global_best_signatures"][0]["tree_count"] = "2"
        with self.assertRaisesRegex(
            PRODUCER.FactorizationBuildError,
            "signature conservation failed",
        ):
            self.validate_single_hcc_row(row)

    def test_rename_no_replace_preserves_existing_target(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            target = root / "target"
            source.write_bytes(b"new candidate bundle")
            target.write_bytes(b"existing trusted bundle")
            before_sha256 = sha256_path(target)
            with self.assertRaisesRegex(
                PRODUCER.FactorizationBuildError,
                "atomic no-replace rename failed",
            ):
                PRODUCER.rename_no_replace(source, target)
            self.assertTrue(source.is_file())
            self.assertEqual(target.read_bytes(), b"existing trusted bundle")
            self.assertEqual(sha256_path(target), before_sha256)


if __name__ == "__main__":
    unittest.main()
