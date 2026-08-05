#!/usr/bin/env python3
"""Synthetic contract tests for the builder-ready pattern-methyl sidecar."""

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "build_pattern_methyl_sidecar.py"
)
SPEC = importlib.util.spec_from_file_location("build_pattern_methyl_sidecar", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
sidecar = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = sidecar
SPEC.loader.exec_module(sidecar)


class PatternMethylSidecarTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.evidence_path = self.root / "pattern_methyl_evidence.v1.tsv.gz"
        self.details_path = self.root / "pattern_methyl_details.v1.jsonl.gz"
        self.dataset = "D1"
        self.chrom = "chr1"
        self.region_id = "chr1|PS=10|HP=1|U_SYNTH:B0001"
        self.unit_id = "U_SYNTH"
        self.phase_set = "10"
        self.topology_rows = {
            f"D{index}": self._topology_row(f"D{index}") for index in range(1, 8)
        }
        self.topology_paths = [
            self.root / "topology" / f"D{index}.topology.jsonl"
            for index in range(1, 8)
        ]
        self.evidence_row = self._evidence_row()
        self.detail_payload = self._detail_payload()
        self._write_all_inputs()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _region_for_dataset(self, dataset: str) -> str:
        return self.region_id if dataset == self.dataset else f"chr1|PS=10|HP=1|U_{dataset}:B0001"

    def _topology_row(self, dataset: str) -> dict[str, object]:
        return {
            "schema_name": sidecar.TOPOLOGY_SCHEMA_NAME,
            "schema_version": sidecar.SCHEMA_VERSION,
            "sample": dataset,
            "chrom": self.chrom,
            "region_id": self._region_for_dataset(dataset),
            "unit_id": self.unit_id if dataset == self.dataset else f"U_{dataset}",
            "phase_set": self.phase_set,
            "hp_family": "1",
            "active_positions": [100, 200],
            "active_bit_count": 2,
            "original_bit_count": 2,
            "best_tree_unique": True,
            "best_tree_tie_count": "1",
            "minimum_vertex_set_count": 2,
            "best_vertex_set_count": 1,
            "total_tree_count": "2",
            "objective_h": 1,
            "search_nodes": 3,
            "af_basis": "synthetic frozen caller AF",
            "af_coverage": [
                {
                    "position": 100,
                    "ref_reads": 30,
                    "alt_reads": 10,
                    "fraction": "1/4",
                },
                {
                    "position": 200,
                    "ref_reads": 20,
                    "alt_reads": 20,
                    "fraction": "1/2",
                },
            ],
            "best_score_fraction": "1/4",
            "input_vaf_eligible": True,
            "representative_best_vertices": [
                {"vertex": 0, "label": "ROOT"},
                {"vertex": 1, "label": "H_AR"},
                {"vertex": 3, "label": "AA"},
            ],
            "representative_best_edges": [
                {
                    "parent_vertex": 0,
                    "child_vertex": 1,
                    "parent_label": "ROOT",
                    "child_label": "H_AR",
                    "acquired_active_bit": 0,
                    "acquired_position": 100,
                    "edge_score_fraction": "0/1",
                },
                {
                    "parent_vertex": 1,
                    "child_vertex": 3,
                    "parent_label": "H_AR",
                    "child_label": "AA",
                    "acquired_active_bit": 1,
                    "acquired_position": 200,
                    "edge_score_fraction": "1/4",
                },
            ],
            "representative_best_morphology": {
                "shape_signature": "((()))",
                "shape_sha256": "a" * 64,
            },
        }

    def _evidence_row(self) -> dict[str, str]:
        row = {field: "" for field in sidecar.EVIDENCE_FIELDS}
        row.update(
            {
                "schema_version": sidecar.SCHEMA_VERSION,
                "dataset": self.dataset,
                "chrom": self.chrom,
                "region_id": self.region_id,
                "unit_id": self.unit_id,
                "phase_set": self.phase_set,
                "hp_family": "1",
                "hp_raw": "1-1",
                "active_positions": "100,200",
                "n_active_bits": "2",
                "pair_full4": "true",
                "k_ge_3": "false",
                "input_n_complete": "40",
                "input_state_counts_json": '{"AA":10,"AR":10,"RA":10,"RR":10}',
                "analysis_n": "40",
                "analysis_state_counts_json": '{"AA":10,"AR":10,"RA":10,"RR":10}',
                "n_common_cpg": "3",
                "n_distal_cpg": "0",
                "qname_join_fraction_min": "1",
                "tile_overlap_conflicts": "0",
                "exchangeable_strata": "2",
                "exchangeable_n": "40",
                "permanova_pseudo_f": "12.5",
                "permanova_r2": "0.3",
                "permanova_p": "0.001",
                "permanova_permutations_requested": "999",
                "permanova_permutations_realized": "999",
                "permdisp_f": "0.4",
                "permdisp_p": "0.5",
                "best_pair": "RR|AR",
                "best_pair_hamming": "1",
                "best_pair_between_mean": "0.8",
                "best_pair_pooled_within_mean": "0.1",
                "best_pair_distance_contrast": "0.7",
                "best_pair_standardized_effect": "1.5",
                "best_pair_topology_relation": sidecar.HALO_RELATION,
                "max_geometry_smd": "0.2",
                "geometry_feature": "mapq:RR|AR",
                "all_states_n8": "true",
                "all_states_n10": "true",
                "equal_n_r2": "0.25",
                "equal_n_retention": "0.833333",
                "rarefaction_median_r2": "0.27",
                "rarefaction_retention": "0.9",
                "distal_r2": "",
                "distal_p": "",
                "distal_permutations_realized": "0",
                "distal_retention": "",
                "multiplicity_family": "CONFIRMATORY_FULL4_OR_LONG",
                "q_by": "0.01",
                "p_holm": "0.02",
                "assessment": "ROBUST_ASSOCIATION",
                "evaluation_status": "EVALUABLE",
                "invalid_reason": "",
            }
        )
        return row

    @staticmethod
    def _digest(index: int) -> str:
        return hashlib.sha256(f"read-{index}".encode("utf-8")).hexdigest()

    def _detail_payload(self) -> dict[str, object]:
        patterns = ["RR", "AR", "RA", "AA"]
        read_order = [
            {
                "qname_sha256": self._digest(index),
                "pattern": patterns[index % 4],
                "read_group": f"RG{index % 2}",
                "strand": "+" if index % 2 == 0 else "-",
            }
            for index in range(40)
        ]
        distance_matrix = [
            [
                0.0
                if first == second
                else 0.1
                if patterns[first % 4] == patterns[second % 4]
                else 0.8
                for second in range(40)
            ]
            for first in range(40)
        ]
        topology = self.topology_rows[self.dataset]

        def effect(
            first: str,
            second: str,
            between_mean: float,
            pooled_within_mean: float,
            standardized_effect: float,
        ) -> dict[str, object]:
            return {
                "first": first,
                "second": second,
                "hamming": sidecar.hamming_distance(first, second),
                "between_mean": between_mean,
                "pooled_within_mean": pooled_within_mean,
                "distance_contrast": between_mean - pooled_within_mean,
                "standardized_effect": standardized_effect,
                "topology_relation": sidecar.expected_topology_relation(
                    first, second, topology
                ),
            }

        return {
            "schema_name": sidecar.DETAIL_SCHEMA_NAME,
            "schema_version": sidecar.SCHEMA_VERSION,
            "dataset": self.dataset,
            "chrom": self.chrom,
            "region_id": self.region_id,
            "hp_raw": "1-1",
            "active_positions": [100, 200],
            "cpg_positions": [300, 400, 500],
            "state_counts": {"AA": 10, "AR": 10, "RA": 10, "RR": 10},
            "state_mean_profiles": {
                "AA": [0.9, None, 0.8],
                "AR": [0.8, 0.7, 0.9],
                "RA": [0.2, 0.3, 0.1],
                "RR": [0.1, 0.2, 0.1],
            },
            "pairwise_effects": [
                effect("RR", "AR", 0.8, 0.1, 1.5),
                effect("RR", "AA", 0.8, 0.15, 1.4),
                effect("AA", "AR", 0.7, 0.15, 1.2),
                effect("AA", "RA", 0.6, 0.15, 1.1),
                effect("AR", "RA", 0.5, 0.15, 0.9),
                effect("RA", "RR", 0.4, 0.15, 0.7),
            ],
            "topology": {
                field: topology[field] for field in sidecar.TOPOLOGY_SNAPSHOT_FIELDS
            },
            "balanced_n_per_state": 10,
            "rarefaction_r2": [0.25, 0.27],
            "read_order": read_order,
            "distance_matrix": distance_matrix,
            "assessment": "ROBUST_ASSOCIATION",
        }

    def _write_evidence(self) -> None:
        with gzip.open(
            self.evidence_path, "wt", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=sidecar.EVIDENCE_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow(self.evidence_row)

    def _write_details(self, payloads: list[dict[str, object]] | None = None) -> None:
        rows = [self.detail_payload] if payloads is None else payloads
        with gzip.open(self.details_path, "wt", encoding="utf-8", newline="\n") as handle:
            for payload in rows:
                handle.write(
                    json.dumps(
                        payload,
                        ensure_ascii=False,
                        sort_keys=True,
                        separators=(",", ":"),
                    )
                )
                handle.write("\n")

    def _write_topologies(self) -> None:
        for path in self.topology_paths:
            dataset = path.name.split(".")[0]
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(
                json.dumps(
                    self.topology_rows[dataset],
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                )
                + "\n",
                encoding="utf-8",
            )

    def _write_all_inputs(self) -> None:
        self._write_evidence()
        self._write_details()
        self._write_topologies()

    def _build(self, matrix_top_n: int = 1) -> dict[str, object]:
        return sidecar.build_sidecar(
            self.evidence_path,
            self.details_path,
            self.topology_paths,
            matrix_top_n=matrix_top_n,
        )

    @staticmethod
    def _only_bundle(payload: dict[str, object]) -> dict[str, object]:
        datasets = payload["datasets"]
        assert isinstance(datasets, list)
        regions = datasets[0]["regions"]
        return regions[0]["exact_raw_hp_bundles"][0]

    def test_builds_exact_raw_hp_overlay_with_strict_projection_contract(self) -> None:
        payload = self._build()
        self.assertEqual(payload["status"], "PASS")
        self.assertEqual(payload["summary"]["topology_authority_files"], 7)
        self.assertEqual(payload["summary"]["exact_raw_hp_bundles"], 1)
        bundle = self._only_bundle(payload)
        self.assertEqual(bundle["grain"]["hp_raw"], "1-1")
        self.assertEqual(bundle["grain"]["phase_set"], self.phase_set)
        self.assertEqual(
            {node["pattern"] for node in bundle["complete_node_projection"]},
            {"RR", "AR", "RA", "AA"},
        )
        self.assertEqual(len(bundle["edge_halos"]), 2)
        self.assertEqual(
            bundle["edge_halos"][0]["relation"], sidecar.HALO_RELATION
        )
        self.assertEqual(len(bundle["pair_bands"]), 2)
        self.assertEqual(
            bundle["pair_bands"][0]["topology_relation"],
            sidecar.HAMMING_GT1_RELATION,
        )
        self.assertTrue(bundle["detail"]["large_matrix_embedded"])
        self.assertIsNone(bundle["evidence"]["distal_r2"])
        self.assertIsNone(
            bundle["detail"]["top_case_payload"]["state_mean_profiles"]["AA"][1]
        )
        self.assertEqual(
            bundle["partial_subcube_evidence"]["source_status"],
            sidecar.PARTIAL_X_SOURCE_STATUS,
        )
        self.assertEqual(bundle["partial_subcube_evidence"]["states"], [])
        self.assertEqual(bundle["partial_subcube_evidence"]["pairs"], [])
        anchor = payload["datasets"][0]["regions"][0]["topology_anchor"]
        self.assertFalse(anchor["topology_counts_af_selected_tree_embedded"])
        self.assertNotIn("af_coverage", anchor)
        self.assertNotIn("representative_best_edges", anchor)
        self.assertEqual(
            payload["external_tables"]["full_evidence_tsv"]["sha256"],
            sidecar.sha256_file(self.evidence_path),
        )

    def test_unused_no_active_alt_topology_row_is_accepted(self) -> None:
        unused = dict(self.topology_rows["D7"])
        unused.update(
            {
                "region_id": "chr1|PS=20|HP=2|U_UNUSED:B0001",
                "unit_id": "U_UNUSED",
                "phase_set": "20",
                "hp_family": "2",
                "active_positions": [],
                "active_bit_count": 0,
                "best_tree_unique": None,
                "representative_best_edges": [],
                "unit_status": "no_active_alt",
            }
        )
        unused.pop("representative_best_vertices")
        unused.pop("representative_best_morphology")
        path = self.topology_paths[-1]
        path.write_text(
            json.dumps(
                unused,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n"
            + json.dumps(
                self.topology_rows["D7"],
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )
        payload = self._build()
        self.assertEqual(payload["status"], "PASS")
        binding = next(
            item
            for item in payload["sources"]["topology_jsonl"]
            if item["dataset"] == "D7"
        )
        self.assertEqual(binding["row_count"], 2)

    def test_evidence_linked_no_active_alt_topology_row_fails_closed(self) -> None:
        topology = self.topology_rows[self.dataset]
        topology.update(
            {
                "active_positions": [],
                "active_bit_count": 0,
                "best_tree_unique": None,
                "representative_best_edges": [],
                "unit_status": "no_active_alt",
            }
        )
        topology.pop("representative_best_vertices")
        topology.pop("representative_best_morphology")
        self._write_topologies()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError,
            "evidence-linked topology requires at least two active positions",
        ):
            self._build()

    def test_evidence_linked_unranked_topology_preserves_association_without_edges(
        self,
    ) -> None:
        topology = self.topology_rows[self.dataset]
        topology.update(
            {
                "best_tree_unique": None,
                "best_tree_tie_count": None,
                "representative_best_edges": [],
                "unit_status": "zero_denominator",
                "read_af_status": "zero_denominator",
            }
        )
        topology.pop("representative_best_vertices")
        topology.pop("representative_best_morphology")
        self.detail_payload["topology"] = sidecar.normalized_topology_snapshot(topology)
        for effect in self.detail_payload["pairwise_effects"]:
            effect["topology_relation"] = sidecar.expected_topology_relation(
                effect["first"], effect["second"], topology
            )
        self.evidence_row["best_pair_topology_relation"] = (
            sidecar.expected_topology_relation("RR", "AR", topology)
        )
        self._write_topologies()
        self._write_evidence()
        self._write_details()
        payload = self._build()
        bundle = self._only_bundle(payload)
        self.assertEqual(bundle["edge_halos"], [])
        self.assertTrue(bundle["pair_associations"])
        self.assertTrue(
            all(
                node["projection"]
                == "TOPOLOGY_VERTEX_REFERENCE_ONLY_SELECTED_EXEMPLAR_UNAVAILABLE"
                for node in bundle["complete_node_projection"]
            )
        )
        anchor = payload["datasets"][0]["regions"][0]["topology_anchor"]
        self.assertFalse(anchor["selected_exemplar_available"])
        self.assertEqual(anchor["authority_unit_status"], "zero_denominator")

    def test_input_x_is_rejected_without_hash_bound_partial_source(self) -> None:
        self.evidence_row["input_state_counts_json"] = (
            '{"AA":10,"AR":10,"RA":10,"RR":10,"RX":2}'
        )
        self._write_evidence()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError,
            "input_state_counts contains X but analyzer evidence is complete-R/A-only",
        ):
            self._build()

    def test_analysis_state_counts_with_x_fail_closed(self) -> None:
        self.evidence_row["analysis_state_counts_json"] = (
            '{"AA":10,"AR":10,"RR":10,"RX":10}'
        )
        self._write_evidence()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError, "analysis_state_counts contains X"
        ):
            self._build()

    def test_detail_read_order_with_x_fails_closed(self) -> None:
        self.detail_payload["read_order"][0]["pattern"] = "RX"
        self._write_details()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError, "detail read_order contains X"
        ):
            self._build()

    def test_detail_pair_effect_with_x_fails_closed(self) -> None:
        self.detail_payload["pairwise_effects"][0]["first"] = "RX"
        self._write_details()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError, "detail pair effect contains X"
        ):
            self._build()

    def test_empty_pairwise_effects_fail_closed(self) -> None:
        self.detail_payload["pairwise_effects"] = []
        self._write_details()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError,
            "not the exact unordered analysis-state pair set",
        ):
            self._build()

    def test_missing_pairwise_effect_fails_closed(self) -> None:
        self.detail_payload["pairwise_effects"].pop()
        self._write_details()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError,
            "not the exact unordered analysis-state pair set",
        ):
            self._build()

    def test_duplicate_unordered_pairwise_effect_fails_closed(self) -> None:
        duplicate = dict(self.detail_payload["pairwise_effects"][0])
        duplicate["first"], duplicate["second"] = (
            duplicate["second"],
            duplicate["first"],
        )
        self.detail_payload["pairwise_effects"].append(duplicate)
        self._write_details()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError, "duplicate unordered detail pair"
        ):
            self._build()

    def test_evidence_best_pair_mismatch_fails_closed(self) -> None:
        second = self.detail_payload["pairwise_effects"][1]
        self.evidence_row.update(
            {
                "best_pair": f"{second['first']}|{second['second']}",
                "best_pair_hamming": str(second["hamming"]),
                "best_pair_between_mean": str(second["between_mean"]),
                "best_pair_pooled_within_mean": str(second["pooled_within_mean"]),
                "best_pair_distance_contrast": str(second["distance_contrast"]),
                "best_pair_standardized_effect": str(second["standardized_effect"]),
                "best_pair_topology_relation": str(second["topology_relation"]),
            }
        )
        self._write_evidence()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError, "evidence/detail best_pair mismatch"
        ):
            self._build()

    def test_detail_metric_matching_tsv_ten_digit_rounding_is_accepted(self) -> None:
        self.detail_payload["pairwise_effects"][0]["between_mean"] = 0.80000000004
        self._write_details()
        payload = self._build()
        bundle = self._only_bundle(payload)
        best = next(
            item
            for item in bundle["pair_associations"]
            if {item["first"], item["second"]} == {"RR", "AR"}
        )
        self.assertEqual(best["between_mean"], 0.80000000004)

    def test_detail_metric_beyond_tsv_ten_digit_rounding_fails_closed(self) -> None:
        self.detail_payload["pairwise_effects"][0]["between_mean"] = 0.8000000006
        self._write_details()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError,
            "evidence/detail best_pair metric mismatch.*between_mean",
        ):
            self._build()

    def test_hamming_gt_one_cannot_claim_edge_halo(self) -> None:
        self.detail_payload["pairwise_effects"][1][
            "topology_relation"
        ] = sidecar.HALO_RELATION
        self._write_details()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError, "pair/topology relation mismatch"
        ):
            self._build()

    def test_topology_unanimity_is_recomputed_not_trusted(self) -> None:
        self.topology_rows[self.dataset]["best_tree_unique"] = False
        self.detail_payload["topology"]["best_tree_unique"] = False
        self._write_topologies()
        self._write_details()
        with self.assertRaisesRegex(
            sidecar.SidecarContractError, "pair/topology relation mismatch"
        ):
            self._build()

    def test_large_matrix_is_not_embedded_when_top_n_is_zero(self) -> None:
        payload = self._build(matrix_top_n=0)
        bundle = self._only_bundle(payload)
        self.assertFalse(bundle["detail"]["large_matrix_embedded"])
        self.assertEqual(
            bundle["detail"]["matrix_source_status"], "AVAILABLE_NOT_EMBEDDED"
        )
        self.assertNotIn("top_case_payload", bundle["detail"])
        self.assertEqual(payload["summary"]["matrix_candidates"], 1)
        self.assertEqual(payload["summary"]["matrices_embedded"], 0)

    def test_sidecar_and_receipt_are_byte_deterministic_and_hash_bound(self) -> None:
        payload = self._build()
        first_dir = self.root / "out-1"
        second_dir = self.root / "out-2"
        first_receipt = sidecar.write_outputs(payload, first_dir)
        second_receipt = sidecar.write_outputs(payload, second_dir)
        first_sidecar = first_dir / sidecar.OUTPUT_FILENAME
        second_sidecar = second_dir / sidecar.OUTPUT_FILENAME
        first_receipt_path = first_dir / sidecar.RECEIPT_FILENAME
        second_receipt_path = second_dir / sidecar.RECEIPT_FILENAME
        self.assertEqual(first_sidecar.read_bytes(), second_sidecar.read_bytes())
        self.assertEqual(
            first_receipt_path.read_bytes(), second_receipt_path.read_bytes()
        )
        self.assertEqual(first_sidecar.read_bytes()[4:8], b"\x00\x00\x00\x00")
        self.assertEqual(
            first_receipt["output"]["sha256"], sidecar.sha256_file(first_sidecar)
        )
        self.assertEqual(first_receipt, second_receipt)
        self.assertTrue(first_receipt["overlay_only"])
        self.assertEqual(
            first_receipt["topology_authority_unchanged"],
            {"counts": True, "af": True, "selected_tree": True},
        )
        self.assertTrue(all(first_receipt["checks"].values()))

    def test_failed_staging_files_are_archived_not_deleted(self) -> None:
        gzip_output = self.root / "failed-gzip" / "sidecar.json.gz"
        gzip_output.parent.mkdir()
        with mock.patch.object(
            sidecar.json, "dump", side_effect=RuntimeError("synthetic gzip failure")
        ):
            with self.assertRaisesRegex(RuntimeError, "synthetic gzip failure"):
                sidecar.deterministic_gzip_json(gzip_output, {"status": "PASS"})
        gzip_archive = gzip_output.parent / sidecar.FAILED_STAGING_DIRECTORY
        self.assertFalse(gzip_output.exists())
        self.assertEqual(len(list(gzip_archive.iterdir())), 1)

        text_output = self.root / "failed-text" / "receipt.json"
        text_output.parent.mkdir()
        with mock.patch.object(
            sidecar.os, "fsync", side_effect=OSError("synthetic fsync failure")
        ):
            with self.assertRaisesRegex(OSError, "synthetic fsync failure"):
                sidecar.atomic_write_text(text_output, '{"status":"PASS"}\n')
        text_archive = text_output.parent / sidecar.FAILED_STAGING_DIRECTORY
        self.assertFalse(text_output.exists())
        archived_text = list(text_archive.iterdir())
        self.assertEqual(len(archived_text), 1)
        self.assertEqual(archived_text[0].read_text(encoding="utf-8"), '{"status":"PASS"}\n')

    def test_cli_requires_exactly_seven_topology_authorities_and_writes_nothing(self) -> None:
        output_dir = self.root / "cli-fail"
        arguments = [
            "--evidence-tsv",
            str(self.evidence_path),
            "--details-jsonl",
            str(self.details_path),
        ]
        for path in self.topology_paths[:6]:
            arguments.extend(["--topology-jsonl", str(path)])
        arguments.extend(["--output-dir", str(output_dir)])
        self.assertEqual(sidecar.main(arguments), 2)
        self.assertFalse((output_dir / sidecar.OUTPUT_FILENAME).exists())
        self.assertFalse((output_dir / sidecar.RECEIPT_FILENAME).exists())

    def test_cli_success_writes_named_sidecar_and_receipt(self) -> None:
        output_dir = self.root / "cli-pass"
        arguments = [
            "--evidence-tsv",
            str(self.evidence_path),
            "--details-jsonl",
            str(self.details_path),
        ]
        for path in reversed(self.topology_paths):
            arguments.extend(["--topology-jsonl", str(path)])
        arguments.extend(
            [
                "--output-dir",
                str(output_dir),
                "--matrix-top-n",
                "1",
            ]
        )
        self.assertEqual(sidecar.main(arguments), 0)
        output_path = output_dir / sidecar.OUTPUT_FILENAME
        receipt_path = output_dir / sidecar.RECEIPT_FILENAME
        self.assertTrue(output_path.is_file())
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        self.assertEqual(receipt["status"], "PASS")
        self.assertTrue(receipt["overlay_only"])
        self.assertTrue(receipt["topology_authority_unchanged"]["counts"])
        self.assertTrue(receipt["topology_authority_unchanged"]["af"])
        self.assertTrue(receipt["topology_authority_unchanged"]["selected_tree"])
        self.assertEqual(receipt["output"]["sha256"], sidecar.sha256_file(output_path))

    def test_orphan_detail_record_fails_closed(self) -> None:
        orphan = dict(self.detail_payload)
        orphan["hp_raw"] = "2-1"
        self._write_details([self.detail_payload, orphan])
        with self.assertRaisesRegex(sidecar.SidecarContractError, "orphan detail row"):
            self._build()


if __name__ == "__main__":
    unittest.main()
