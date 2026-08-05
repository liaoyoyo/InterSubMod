#!/usr/bin/env python3
"""Synthetic tests for the multisite methyl artifact catalog."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts" / "build_artifact_catalog.py"
)
SPEC = importlib.util.spec_from_file_location("build_artifact_catalog", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
catalog = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = catalog
SPEC.loader.exec_module(catalog)


class ArtifactCatalogTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.root = self.base / "primary"
        self.dataset = "SYNTH"
        self.chrom = "chr1"
        self.position = 10001
        self.dataset_dir = self.root / self.dataset
        self.region = (
            self.dataset_dir
            / f"{self.dataset}{catalog.ANALYSIS_SUFFIX}"
            / self.chrom
            / f"{self.chrom}_{self.position}"
            / f"{self.chrom}_{self.position - 5000}_{self.position + 5000}"
        )
        self._write_artifacts(include_hp=True)
        self._write_receipt(pass_value=True)
        self.markers = self.base / "markers.tsv"
        self._write_markers([(self.dataset, self.chrom, self.position)])

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def _write_rows(path: Path, rows: list[list[object]], delimiter: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter=delimiter, lineterminator="\n")
            writer.writerows(rows)

    def _write_artifacts(self, *, include_hp: bool) -> None:
        reads_header = [
            "read_id",
            "read_name",
            "chr",
            "start",
            "end",
            "mapq",
        ]
        if include_hp:
            reads_header.append("hp")
        reads_header.extend(["alt_support", "is_tumor", "strand"])
        reads_rows: list[list[object]] = [reads_header]
        for index, hp_value in enumerate(["UNTRUSTED", "999", ""]):
            row: list[object] = [
                str(index),
                f"read-{index}",
                self.chrom,
                9000,
                11000,
                60,
            ]
            if include_hp:
                row.append(hp_value)
            row.extend(["ALT" if index else "REF", 1, "+" if index != 1 else "-"])
            reads_rows.append(row)
        self._write_rows(self.region / "reads" / "reads.tsv", reads_rows, "\t")

        self._write_rows(
            self.region / "methylation" / "methylation.csv",
            [
                ["read_id", "9900", "10100"],
                ["0", "0.1", "0.9"],
                ["1", "0.2", "NA"],
                ["2", "0.8", "0.7"],
            ],
            ",",
        )
        self._write_rows(
            self.region / "distance" / "BERNOULLI" / "matrix.csv",
            [
                ["read_id", "0", "1", "2"],
                ["0", "0", "0.2", "0.7"],
                ["1", "0.2", "0", "0.5"],
                ["2", "0.7", "0.5", "0"],
            ],
            ",",
        )

    def _write_receipt(self, *, pass_value: bool) -> None:
        self.dataset_dir.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema_name": "intersubmod.all_ssnv_site_run",
            "schema_version": "1.0.0",
            "sample": self.dataset,
            "command": [
                "/synthetic/inter_sub_mod",
                "-w",
                "5000",
                "--distance-metric",
                "BERNOULLI",
            ],
            "output_dir": str(self.dataset_dir.resolve()),
            "validation": {
                "expected_vcf_sites": 1,
                "reads_files": 1,
                "methylation_files": 1,
                "bernoulli_matrix_files": 1,
                "pass": pass_value,
            },
            "pass": pass_value,
        }
        (self.dataset_dir / "run_receipt.json").write_text(
            json.dumps(payload), encoding="utf-8"
        )

    def _write_markers(self, markers: list[tuple[str, str, int]]) -> None:
        with self.markers.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["dataset", "chrom", "position1"])
            writer.writerows(markers)

    def test_build_catalog_passes_and_hp_is_non_authoritative(self) -> None:
        document = catalog.build_catalog(
            self.markers,
            self.root,
            created_at_utc="2026-07-27T00:00:00+00:00",
        )
        self.assertTrue(document["pass"])
        self.assertEqual(document["status"], "PASS")
        self.assertEqual(document["summary"]["markers_pass"], 1)
        self.assertEqual(
            document["contract"]["reads_hp_authority"],
            "NOT_AUTHORITATIVE_AND_NOT_CONSUMED",
        )
        record = document["records"][0]
        self.assertEqual(record["status"], "PASS")
        self.assertEqual(record["reads_row_count"], 3)
        self.assertEqual(record["methylation_row_count"], 3)
        self.assertEqual(record["methylation_cpg_count"], 2)
        self.assertEqual(record["bernoulli_row_count"], 3)
        self.assertEqual(record["bernoulli_column_count"], 3)
        self.assertEqual(record["reads_hp_column_present"], 1)
        self.assertEqual(record["reads_hp_authoritative"], 0)
        self.assertEqual(len(record["reads_sha256"]), 64)
        self.assertEqual(len(record["artifact_bundle_sha256"]), 64)

    def test_catalog_does_not_require_hp_column(self) -> None:
        self._write_artifacts(include_hp=False)
        document = catalog.build_catalog(self.markers, self.root)
        self.assertTrue(document["pass"])
        self.assertEqual(document["records"][0]["reads_hp_column_present"], 0)

    def test_duplicate_input_markers_are_collapsed_with_count(self) -> None:
        self._write_markers(
            [
                (self.dataset, self.chrom, self.position),
                (self.dataset, self.chrom, self.position),
            ]
        )
        document = catalog.build_catalog(self.markers, self.root)
        self.assertTrue(document["pass"])
        self.assertEqual(document["source"]["input_rows"], 2)
        self.assertEqual(document["source"]["unique_markers"], 1)
        self.assertEqual(document["source"]["duplicate_rows_collapsed"], 1)
        self.assertEqual(document["records"][0]["input_occurrences"], 2)

    def test_missing_bernoulli_fails_closed_and_cli_emits_diagnostics(self) -> None:
        (self.region / "distance" / "BERNOULLI" / "matrix.csv").unlink()
        output_tsv = self.base / "catalog.tsv"
        output_json = self.base / "catalog.json"
        exit_code = catalog.main(
            [
                "--markers",
                str(self.markers),
                "--primary-root",
                str(self.root),
                "--output-tsv",
                str(output_tsv),
                "--output-json",
                str(output_json),
            ]
        )
        self.assertEqual(exit_code, 2)
        payload = json.loads(output_json.read_text(encoding="utf-8"))
        self.assertFalse(payload["pass"])
        self.assertEqual(payload["status"], "FAIL_CLOSED")
        self.assertEqual(payload["summary"]["markers_fail"], 1)
        self.assertEqual(payload["records"][0]["status"], "FAIL_ARTIFACT_MISSING")
        self.assertTrue(output_tsv.is_file())

    def test_methylation_read_id_mismatch_fails_closed(self) -> None:
        self._write_rows(
            self.region / "methylation" / "methylation.csv",
            [
                ["read_id", "9900", "10100"],
                ["0", "0.1", "0.9"],
                ["1", "0.2", "NA"],
                ["DIFFERENT", "0.8", "0.7"],
            ],
            ",",
        )
        document = catalog.build_catalog(self.markers, self.root)
        self.assertFalse(document["pass"])
        self.assertEqual(document["records"][0]["status"], "FAIL_READ_ID_MISMATCH")

    def test_invalid_receipt_fails_before_artifact_promotion(self) -> None:
        self._write_receipt(pass_value=False)
        document = catalog.build_catalog(self.markers, self.root)
        self.assertFalse(document["pass"])
        record = document["records"][0]
        self.assertEqual(record["status"], "FAIL_RECEIPT_INVALID")
        self.assertEqual(record["receipt_pass"], 0)
        self.assertEqual(record["reads_path"], "")

    def test_invalid_marker_header_raises_catalog_error(self) -> None:
        self.markers.write_text("dataset\tchrom\nSYNTH\tchr1\n", encoding="utf-8")
        with self.assertRaises(catalog.CatalogError) as context:
            catalog.build_catalog(self.markers, self.root)
        self.assertEqual(context.exception.code, "MARKER_TSV_INVALID")

    def test_write_catalog_round_trip(self) -> None:
        document = catalog.build_catalog(
            self.markers,
            self.root,
            created_at_utc="2026-07-27T00:00:00+00:00",
        )
        output_tsv = self.base / "out" / "catalog.v1.tsv"
        output_json = self.base / "out" / "catalog.v1.json"
        catalog.write_catalog(document, output_tsv, output_json)
        payload = json.loads(output_json.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "1.0.0")
        with output_tsv.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["status"], "PASS")
        self.assertEqual(rows[0]["methylation_cpg_count"], "2")


if __name__ == "__main__":
    unittest.main()
