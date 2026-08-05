#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts" / "select_formal_candidates.py"
)
SPEC = importlib.util.spec_from_file_location("select_formal_candidates", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
selection = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = selection
SPEC.loader.exec_module(selection)


class SelectFormalCandidatesTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.counts = self.root / "counts.tsv"
        self.markers = self.root / "markers.tsv"
        with self.counts.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "dataset",
                    "chrom",
                    "region_id",
                    "hp_raw",
                    "formal_n5",
                    "pair_full4",
                    "k_ge_3",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(
                [
                    {
                        "dataset": "S",
                        "chrom": "chr1",
                        "region_id": "r1",
                        "hp_raw": "1-1",
                        "formal_n5": "true",
                        "pair_full4": "true",
                        "k_ge_3": "false",
                    },
                    {
                        "dataset": "S",
                        "chrom": "chr1",
                        "region_id": "r2",
                        "hp_raw": "2",
                        "formal_n5": "false",
                        "pair_full4": "false",
                        "k_ge_3": "true",
                    },
                ]
            )
        with self.markers.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["dataset", "chrom", "region_id", "position1"])
            writer.writerow(["S", "chr1", "r1", 100])
            writer.writerow(["S", "chr1", "r1", 200])
            writer.writerow(["S", "chr1", "r2", 300])

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_selects_only_formal_region_markers(self) -> None:
        fields, formal, markers, summary = selection.select_candidates(
            self.counts, self.markers
        )
        self.assertIn("formal_n5", fields)
        self.assertEqual(len(formal), 1)
        self.assertEqual([row["position1"] for row in markers], ["100", "200"])
        self.assertEqual(summary["formal_units"], 1)
        self.assertEqual(summary["pair_full4"], 1)

    def test_missing_position1_fails_closed(self) -> None:
        self.markers.write_text(
            "dataset\tchrom\tregion_id\tposition\nS\tchr1\tr1\t100\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(selection.SelectionError, "position1"):
            selection.select_candidates(self.counts, self.markers)


if __name__ == "__main__":
    unittest.main()
