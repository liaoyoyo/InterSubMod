#!/usr/bin/env python3
"""Synthetic tests for formal-marker BERNOULLI parity validation."""

from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "validate_bernoulli_parity.py"
)
SPEC = importlib.util.spec_from_file_location("validate_bernoulli_parity", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
parity = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = parity
SPEC.loader.exec_module(parity)


class BernoulliParityTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.methylation = self.root / "methylation.csv"
        self.bernoulli = self.root / "matrix.csv"
        self.catalog = self.root / "catalog.tsv"
        self._write_csv(
            self.methylation,
            [
                ["read_id", "101", "102", "103"],
                ["r0", "0", "0", "0"],
                ["r1", "1", "1", "1"],
                ["r2", "0", "0", "0"],
            ],
        )
        self._write_bernoulli("1")
        with self.catalog.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "dataset",
                    "chrom",
                    "position1",
                    "status",
                    "methylation_path",
                    "bernoulli_path",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow(
                {
                    "dataset": "SYNTH",
                    "chrom": "chr1",
                    "position1": "100",
                    "status": "PASS",
                    "methylation_path": self.methylation,
                    "bernoulli_path": self.bernoulli,
                }
            )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def _write_csv(path: Path, rows: list[list[str]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            csv.writer(handle, lineterminator="\n").writerows(rows)

    def _write_bernoulli(self, between: str) -> None:
        self._write_csv(
            self.bernoulli,
            [
                ["read_id", "r0", "r1", "r2"],
                ["r0", "0", between, "0"],
                ["r1", between, "0", between],
                ["r2", "0", between, "0"],
            ],
        )

    def test_every_cell_matches(self) -> None:
        summary = parity.execute(
            self.catalog,
            self.root / "parity.tsv",
            self.root / "parity.json",
            max_reads=3,
            tolerance=1e-6,
        )
        self.assertTrue(summary["all_pass"])
        self.assertEqual(summary["counts"]["markers_pass"], 1)
        self.assertEqual(summary["counts"]["pair_cells_checked"], 3)
        self.assertEqual(summary["max_absolute_error"], 0.0)

    def test_numeric_drift_is_reported(self) -> None:
        self._write_bernoulli("0.8")
        result = parity.validate_marker(
            {
                "dataset": "SYNTH",
                "chrom": "chr1",
                "position1": "100",
                "methylation_path": str(self.methylation),
                "bernoulli_path": str(self.bernoulli),
            },
            max_reads=3,
            tolerance=1e-6,
        )
        self.assertEqual(result["status"], "FAIL")
        self.assertAlmostEqual(result["max_absolute_error"], 0.2)


if __name__ == "__main__":
    unittest.main()
