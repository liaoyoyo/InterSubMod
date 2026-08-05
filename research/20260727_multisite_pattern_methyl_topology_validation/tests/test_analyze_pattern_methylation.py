#!/usr/bin/env python3
"""Contract tests for exact-raw-HP pattern x methylation analysis."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "analyze_pattern_methylation.py"
)
SPEC = importlib.util.spec_from_file_location("analyze_pattern_methylation", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
analysis = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = analysis
SPEC.loader.exec_module(analysis)


class PatternMethylAnalysisTest(unittest.TestCase):
    def test_bernoulli_matches_expected_deterministic_disagreement(self) -> None:
        matrix = np.asarray([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]])
        distance = analysis.bernoulli_distance_matrix(matrix, min_common=3)
        self.assertAlmostEqual(distance[0, 1], 1.0)
        self.assertAlmostEqual(distance[0, 2], 0.0)
        self.assertAlmostEqual(distance[1, 2], 1.0)
        self.assertTrue(np.allclose(np.diag(distance), 0.0))

    def test_invalid_distance_stays_nan(self) -> None:
        matrix = np.asarray([[0.0, math.nan, math.nan], [1.0, math.nan, math.nan]])
        distance = analysis.bernoulli_distance_matrix(matrix, min_common=3)
        self.assertTrue(math.isnan(distance[0, 1]))
        self.assertNotEqual(distance[0, 1], 1.0)

    def test_exchangeable_indices_exclude_monomorphic_strata(self) -> None:
        labels = ["RR", "RR", "RA", "AA", "RA", "AA"]
        strata = ["mono", "mono", "mixed", "mixed", "mixed2", "mixed2"]
        indices = analysis.exchangeable_indices(labels, strata)
        self.assertEqual(indices.tolist(), [2, 3, 4, 5])

    def test_restricted_permutation_is_deterministic(self) -> None:
        labels = ["RR", "RA"] * 4
        strata = [f"s{index // 2}" for index in range(8)]
        matrix = np.asarray(
            [
                [0.0 if i == j else (0.1 if labels[i] == labels[j] else 0.9) for j in range(8)]
                for i in range(8)
            ]
        )
        first = analysis.restricted_permanova_permdisp(
            matrix, labels, strata, permutations=49, seed=7
        )
        second = analysis.restricted_permanova_permdisp(
            matrix, labels, strata, permutations=49, seed=7
        )
        self.assertEqual(first, second)
        self.assertGreater(first["r_squared"], 0.9)
        self.assertEqual(first["realized"], 49)

    def test_restricted_permutation_fails_without_exchangeability(self) -> None:
        matrix = np.ones((6, 6), dtype=float)
        np.fill_diagonal(matrix, 0.0)
        with self.assertRaisesRegex(ValueError, "no_exchangeable_labels"):
            analysis.restricted_permanova_permdisp(
                matrix,
                ["RR", "RR", "RR", "AA", "AA", "AA"],
                ["left"] * 3 + ["right"] * 3,
                permutations=9,
                seed=1,
            )

    def test_by_adjustment_uses_complete_family_not_adaptive_subset(self) -> None:
        complete_family = [
            {
                "pair_full4": "true",
                "k_ge_3": "false",
                "evaluation_status": "EVALUABLE",
                "permanova_p": 0.001 if index == 0 else 0.5,
            }
            for index in range(100)
        ]
        subset = [dict(complete_family[0])]
        analysis.adjust_p_values(complete_family)
        analysis.adjust_p_values(subset)
        self.assertAlmostEqual(subset[0]["q_by"], 0.001)
        self.assertAlmostEqual(
            complete_family[0]["q_by"], 0.518737751763962
        )
        self.assertGreater(complete_family[0]["q_by"], 0.05)

    def test_partial_pattern_cannot_project_to_topology(self) -> None:
        with self.assertRaisesRegex(ValueError, "partial pattern"):
            analysis.pattern_vertex("RX")

    def test_hamming_gt_one_is_pair_band_only(self) -> None:
        relation = analysis.topology_pair_relation(
            "RR",
            "AA",
            {
                "best_tree_unique": True,
                "representative_best_edges": [
                    {"parent_vertex": 0, "child_vertex": 1}
                ],
            },
        )
        self.assertEqual(relation, "PAIR_BAND_ONLY_HAMMING_GT1")

    def test_only_unique_best_hamming_one_edge_gets_projection(self) -> None:
        topology = {
            "best_tree_unique": True,
            "representative_best_edges": [{"parent_vertex": 0, "child_vertex": 1}],
        }
        self.assertEqual(
            analysis.topology_pair_relation("RR", "AR", topology),
            "HAMMING1_GLOBAL_BEST_UNANIMOUS",
        )
        topology["best_tree_unique"] = False
        self.assertEqual(
            analysis.topology_pair_relation("RR", "AR", topology),
            "HAMMING1_NOT_UNANIMOUS",
        )

    def test_geometry_smd_detects_large_read_length_imbalance(self) -> None:
        rows = []
        labels = []
        for index in range(10):
            rows.append({"mapq": "60", "start0": "0", "end0": str(100 + index)})
            labels.append("RR")
            rows.append({"mapq": "60", "start0": "0", "end0": str(1000 + index)})
            labels.append("AA")
        effect, feature = analysis.geometry_max_smd(rows, labels)
        self.assertGreater(effect, 0.5)
        self.assertTrue(feature.startswith("read_length:"))

    def test_unexpected_unit_failure_fails_closed(self) -> None:
        key = ("SYNTH", "chr1", "region-1", "1")
        count_row = {
            "dataset": "SYNTH",
            "chrom": "chr1",
            "region_id": "region-1",
            "unit_id": "unit-1",
            "phase_set": "100",
            "hp_family": "1",
            "hp_raw": "1",
            "active_positions": "100,200",
            "n_active_bits": "2",
            "n_total": "1",
            "n_complete": "1",
            "n_partial": "0",
            "state_count_json": '{"RA":1}',
            "partial_state_count_json": "{}",
        }
        rows = [
            {
                "dataset": "SYNTH",
                "chrom": "chr1",
                "region_id": "region-1",
                "unit_id": "unit-1",
                "phase_set": "100",
                "hp_family": "1",
                "hp_raw": "1",
                "qname_sha256": "a" * 64,
                "active_positions": "100,200",
                "n_active_bits": "2",
                "pattern": "RA",
                "complete_pattern": "true",
            }
        ]
        with mock.patch.object(
            analysis, "load_methyl_union", side_effect=RuntimeError("synthetic bug")
        ):
            with self.assertRaisesRegex(
                analysis.AnalysisContractError, "unexpected analysis failure"
            ):
                analysis.analyze_unit(
                    key,
                    rows,
                    count_row,
                    {},
                    {},
                    analysis.AnalysisConfig(permutations=9),
                )

    def test_formal_unit_must_match_topology_identity(self) -> None:
        key = ("SYNTH", "chr1", "region-1", "1-1")
        count = {
            "unit_id": "unit-1",
            "phase_set": "100",
            "hp_family": "1",
            "active_positions": "100,200",
        }
        topology = {
            "SYNTH\x1fregion-1": {
                "sample": "SYNTH",
                "chrom": "chr1",
                "unit_id": "unit-1",
                "phase_set": "100",
                "hp_family": "1",
                "active_positions": [100, 200],
            }
        }
        analysis.validate_formal_topology_bindings({key: count}, topology)
        topology["SYNTH\x1fregion-1"]["phase_set"] = "DRIFT"
        with self.assertRaisesRegex(
            analysis.AnalysisContractError, "identity mismatch"
        ):
            analysis.validate_formal_topology_bindings({key: count}, topology)

    def test_candidate_counts_and_identity_are_recomputed(self) -> None:
        key = ("SYNTH", "chr1", "region-1", "1")
        count = {
            "dataset": "SYNTH",
            "chrom": "chr1",
            "region_id": "region-1",
            "unit_id": "unit-1",
            "phase_set": "100",
            "hp_family": "1",
            "hp_raw": "1",
            "active_positions": "100,200",
            "n_active_bits": "2",
            "n_total": "2",
            "n_complete": "1",
            "n_partial": "1",
            "state_count_json": '{"RA":1}',
            "partial_state_count_json": '{"RX":1}',
        }
        rows = [
            {
                **{
                    field: count[field]
                    for field in (
                        "dataset",
                        "chrom",
                        "region_id",
                        "unit_id",
                        "phase_set",
                        "hp_family",
                        "hp_raw",
                        "active_positions",
                        "n_active_bits",
                    )
                },
                "qname_sha256": f"{index:064x}",
                "pattern": pattern,
                "complete_pattern": complete,
            }
            for index, pattern, complete in (
                (1, "RA", "true"),
                (2, "RX", "false"),
            )
        ]
        analysis.validate_candidate_unit(key, rows, count)
        rows[0]["phase_set"] = "DRIFT"
        with self.assertRaisesRegex(
            analysis.AnalysisContractError, "identity mismatch"
        ):
            analysis.validate_candidate_unit(key, rows, count)


class CommonBasisContractTest(unittest.TestCase):
    @staticmethod
    def _row(index: int, state: str) -> dict[str, str]:
        return {
            "qname_sha256": f"{index:064x}",
            "pattern": state,
            "complete_pattern": "true",
            "read_group": "RG1",
            "strand": "+",
        }

    def test_post_filter_loss_of_frozen_n5_state_fails_closed(self) -> None:
        states = ["AA"] * 5 + ["RA"] * 20 + ["RR"] * 25
        rows = [self._row(index + 1, state) for index, state in enumerate(states)]
        cpgs = tuple(range(1000, 1010))
        methyl = {
            row["qname_sha256"]: ({cpg: 0.2 for cpg in cpgs} if index else {})
            for index, row in enumerate(rows)
        }
        config = analysis.AnalysisConfig(
            min_common_cpg=10,
            marker_mask_bp=0,
        )
        with self.assertRaisesRegex(
            ValueError, "post_cpg_filter_formal_state_support_lost"
        ):
            analysis.build_common_basis(rows, methyl, (), config)

    def test_state_coverage_is_rechecked_after_read_filter(self) -> None:
        states = ["AA"] * 10 + ["RR"] * 40
        rows = [self._row(index + 1, state) for index, state in enumerate(states)]
        cpgs = tuple(range(1000, 1010))
        methyl: dict[str, dict[int, float]] = {}
        for index, row in enumerate(rows):
            values = {cpg: 0.2 for cpg in cpgs}
            if index == 0:
                values = {1000: 0.2}
            elif index in {8, 9}:
                values.pop(1000)
            methyl[row["qname_sha256"]] = values
        retained, _matrix, basis, counts = analysis.build_common_basis(
            rows,
            methyl,
            (),
            analysis.AnalysisConfig(min_common_cpg=1, marker_mask_bp=0),
        )
        self.assertNotIn(1000, basis)
        self.assertEqual(len(retained), 49)
        self.assertEqual(counts, {"AA": 9, "RR": 40})

    def test_exchangeability_state_loss_fails_closed(self) -> None:
        rows = [
            {
                **self._row(index + 1, "AA" if index < 5 else "RR"),
                "read_group": "AA_ONLY" if index < 4 else "MIXED",
            }
            for index in range(45)
        ]
        methyl = {
            row["qname_sha256"]: {
                cpg: 0.2 for cpg in range(1000, 1010)
            }
            for row in rows
        }
        with self.assertRaisesRegex(
            ValueError, "post_exchangeability_formal_state_support_lost"
        ):
            analysis.build_common_basis(
                rows,
                methyl,
                (),
                analysis.AnalysisConfig(marker_mask_bp=0),
            )


class CatalogBindingTest(unittest.TestCase):
    def test_catalog_artifact_hash_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reads = root / "reads.tsv"
            methyl = root / "methylation.csv"
            reads.write_text("read_id\tread_name\n0\tread-a\n", encoding="utf-8")
            methyl.write_text("read_id,1000\n0,0.2\n", encoding="utf-8")
            catalog = root / "catalog.tsv"

            def write_catalog() -> None:
                with catalog.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=[
                            "dataset",
                            "chrom",
                            "position1",
                            "status",
                            "reads_path",
                            "reads_size_bytes",
                            "reads_sha256",
                            "methylation_path",
                            "methylation_size_bytes",
                            "methylation_sha256",
                        ],
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerow(
                        {
                            "dataset": "SYNTH",
                            "chrom": "chr1",
                            "position1": 100,
                            "status": "PASS",
                            "reads_path": reads,
                            "reads_size_bytes": reads.stat().st_size,
                            "reads_sha256": analysis.sha256_file(reads),
                            "methylation_path": methyl,
                            "methylation_size_bytes": methyl.stat().st_size,
                            "methylation_sha256": analysis.sha256_file(methyl),
                        }
                    )

            write_catalog()
            analysis.load_catalog(catalog)
            reads.write_text(
                "read_id\tread_name\n0\tread-drift\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(
                analysis.AnalysisContractError, "mismatch"
            ):
                analysis.load_catalog(catalog)


class CandidateManifestTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.shard_dir = self.root / "candidate_read_join"
        self.shard_dir.mkdir()
        self.shards = []
        for index in range(2):
            path = self.shard_dir / f"SYNTH.chr{index + 1}.tsv.gz"
            path.write_bytes(f"shard-{index}\n".encode("utf-8"))
            self.shards.append(path)
        self.manifest = self.root / "candidate_read_join.manifest.tsv"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write_manifest(self, hashes: list[str] | None = None) -> None:
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                ["sort_ordinal", "relative_path", "sha256", "size_bytes"]
            )
            for index, path in enumerate(self.shards):
                writer.writerow(
                    [
                        index,
                        path.relative_to(self.root).as_posix(),
                        (
                            hashes[index]
                            if hashes is not None
                            else analysis.sha256_file(path)
                        ),
                        path.stat().st_size,
                    ]
                )

    def test_manifest_preserves_hash_bound_order(self) -> None:
        self._write_manifest()
        resolved = analysis.load_candidate_shards([], self.manifest)
        self.assertEqual(resolved, [path.resolve() for path in self.shards])

    def test_manifest_fails_closed_on_hash_drift(self) -> None:
        self._write_manifest(["0" * 64, analysis.sha256_file(self.shards[1])])
        with self.assertRaisesRegex(
            analysis.AnalysisContractError, "SHA-256 mismatch"
        ):
            analysis.load_candidate_shards([], self.manifest)

    def test_candidate_shard_same_size_drift_is_rejected_at_read(self) -> None:
        self._write_manifest()
        resolved = analysis.load_candidate_shards([], self.manifest)
        identities = {
            path.resolve(): analysis.file_identity(path) for path in resolved
        }
        original_size = self.shards[0].stat().st_size
        self.shards[0].write_bytes(b"Shard-0\n")
        self.assertEqual(self.shards[0].stat().st_size, original_size)
        with self.assertRaisesRegex(
            analysis.AnalysisContractError, "source identity drift"
        ):
            list(analysis.iter_candidate_groups(resolved, {}, identities))

    def test_unit_key_file_is_exact_and_rejects_duplicates(self) -> None:
        key_file = self.root / "unit_keys.tsv"
        key_file.write_text(
            "dataset\tchrom\tregion_id\thp_raw\n"
            "SYNTH\tchr1\tregion-1\t1-1\n",
            encoding="utf-8",
        )
        self.assertEqual(
            analysis.load_unit_keys(key_file),
            {("SYNTH", "chr1", "region-1", "1-1")},
        )
        key_file.write_text(
            "dataset\tchrom\tregion_id\thp_raw\n"
            "SYNTH\tchr1\tregion-1\t1-1\n"
            "SYNTH\tchr1\tregion-1\t1-1\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(
            analysis.AnalysisContractError, "duplicate key"
        ):
            analysis.load_unit_keys(key_file)

    def _write_adaptive_receipt(
        self,
        unit_keys: Path,
        screen: Path,
        *,
        screen_permutations: int = 999,
        refinement_permutations: int = 49999,
        all_pass: bool = True,
        selected_count: int = 1,
    ) -> Path:
        receipt = self.root / "adaptive_receipt.json"
        payload = {
            "schema_name": analysis.ADAPTIVE_SELECTION_SCHEMA_NAME,
            "schema_version": analysis.ADAPTIVE_SELECTION_SCHEMA_VERSION,
            "all_pass": all_pass,
            "contract": {
                "family": "CONFIRMATORY_FULL4_OR_LONG",
                "screen_permutations": screen_permutations,
                "screen_floor": 1.0 / (screen_permutations + 1),
                "refinement_permutations": refinement_permutations,
            },
            "counts": {"selected_for_refinement": selected_count},
            "inputs": {"screen_evidence": analysis.file_identity(screen)},
            "outputs": {"unit_keys": analysis.file_identity(unit_keys)},
        }
        receipt.write_text(
            json.dumps(payload, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return receipt

    def test_adaptive_receipt_binds_keys_screen_and_frozen_budget(self) -> None:
        unit_keys = self.root / "selected.tsv"
        unit_keys.write_text(
            "dataset\tchrom\tregion_id\thp_raw\n"
            "SYNTH\tchr1\tregion-1\t1-1\n",
            encoding="utf-8",
        )
        screen = self.root / "confirmatory.tsv.gz"
        screen.write_bytes(b"screen-evidence")
        receipt = self._write_adaptive_receipt(unit_keys, screen)
        sources = analysis.validate_adaptive_refinement_binding(
            unit_keys,
            receipt,
            analyzer_permutations=49999,
            selected_key_count=1,
        )
        self.assertEqual(sources["unit_key_file"], analysis.file_identity(unit_keys))
        self.assertEqual(sources["screen_evidence"], analysis.file_identity(screen))
        self.assertEqual(
            sources["unit_key_receipt"], analysis.file_identity(receipt)
        )

    def test_adaptive_receipt_rejects_budget_and_same_size_hash_drift(self) -> None:
        unit_keys = self.root / "selected.tsv"
        unit_keys.write_text(
            "dataset\tchrom\tregion_id\thp_raw\n"
            "SYNTH\tchr1\tregion-1\t1-1\n",
            encoding="utf-8",
        )
        screen = self.root / "confirmatory.tsv.gz"
        screen.write_bytes(b"screen-evidence")
        receipt = self._write_adaptive_receipt(
            unit_keys,
            screen,
            refinement_permutations=50000,
        )
        with self.assertRaisesRegex(
            analysis.AnalysisContractError, "refinement budget mismatch"
        ):
            analysis.validate_adaptive_refinement_binding(
                unit_keys,
                receipt,
                analyzer_permutations=49999,
                selected_key_count=1,
            )

        receipt = self._write_adaptive_receipt(unit_keys, screen)
        screen.write_bytes(b"Screen-evidence")
        self.assertEqual(screen.stat().st_size, len(b"screen-evidence"))
        with self.assertRaisesRegex(
            analysis.AnalysisContractError, "SHA-256 mismatch"
        ):
            analysis.validate_adaptive_refinement_binding(
                unit_keys,
                receipt,
                analyzer_permutations=49999,
                selected_key_count=1,
            )

    def test_unit_key_and_receipt_cli_arguments_are_atomic(self) -> None:
        with self.assertRaisesRegex(
            analysis.AnalysisContractError, "must be provided together"
        ):
            analysis.main(
                [
                    "--pattern-counts",
                    str(self.root / "counts.tsv"),
                    "--artifact-catalog",
                    str(self.root / "catalog.tsv"),
                    "--unit-key-file",
                    str(self.root / "selected.tsv"),
                    "--topology-root",
                    str(self.root / "topology"),
                    "--output-dir",
                    str(self.root / "output"),
                ]
            )


class MethylUnionTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.dataset = "SYNTH"
        self.chrom = "chr1"
        self.qnames = ["read-a", "read-b"]
        self.digests = [analysis.sha256_read_name(name) for name in self.qnames]

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _artifact(
        self, position: int, cpgs: list[int], values: list[list[float]]
    ) -> dict[str, str]:
        region = self.root / str(position)
        reads = region / "reads.tsv"
        methyl = region / "methylation.csv"
        region.mkdir(parents=True)
        with reads.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["read_id", "read_name"])
            for index, name in enumerate(self.qnames):
                writer.writerow([index, name])
        with methyl.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["read_id"] + cpgs)
            for index, row in enumerate(values):
                writer.writerow([index] + row)
        return {
            "reads_path": str(reads),
            "reads_size_bytes": str(reads.stat().st_size),
            "reads_sha256": analysis.sha256_file(reads),
            "methylation_path": str(methyl),
            "methylation_size_bytes": str(methyl.stat().st_size),
            "methylation_sha256": analysis.sha256_file(methyl),
        }

    def _rows(self) -> list[dict[str, str]]:
        return [
            {
                "dataset": self.dataset,
                "chrom": self.chrom,
                "qname_sha256": digest,
            }
            for digest in self.digests
        ]

    def test_union_deduplicates_identical_overlap(self) -> None:
        catalog = {
            (self.dataset, self.chrom, 1000): self._artifact(
                1000, [900, 1100], [[0.1, 0.2], [0.8, 0.7]]
            ),
            (self.dataset, self.chrom, 2000): self._artifact(
                2000, [1100, 2100], [[0.2, 0.3], [0.7, 0.6]]
            ),
        }
        union, audit = analysis.load_methyl_union(
            self._rows(), [1000, 2000], catalog
        )
        self.assertEqual(audit["overlap_conflicts"], 0)
        self.assertEqual(audit["join_fraction_min"], 1.0)
        self.assertEqual(sorted(union[self.digests[0]]), [900, 1100, 2100])

    def test_union_fails_on_overlap_conflict(self) -> None:
        catalog = {
            (self.dataset, self.chrom, 1000): self._artifact(
                1000, [1100], [[0.2], [0.7]]
            ),
            (self.dataset, self.chrom, 2000): self._artifact(
                2000, [1100], [[0.9], [0.7]]
            ),
        }
        with self.assertRaisesRegex(
            analysis.AnalysisContractError, "methylation value conflicts"
        ):
            analysis.load_methyl_union(self._rows(), [1000, 2000], catalog)

    def test_union_rejects_same_size_artifact_drift_at_use(self) -> None:
        record = self._artifact(
            1000, [1100], [[0.2], [0.7]]
        )
        methyl = Path(record["methylation_path"])
        original_size = methyl.stat().st_size
        methyl.write_text(
            methyl.read_text(encoding="utf-8").replace("0.2", "0.9"),
            encoding="utf-8",
        )
        self.assertEqual(methyl.stat().st_size, original_size)
        with self.assertRaisesRegex(
            analysis.AnalysisContractError, "source identity drift"
        ):
            analysis.load_methyl_union(
                self._rows(),
                [1000],
                {(self.dataset, self.chrom, 1000): record},
            )


if __name__ == "__main__":
    unittest.main()
