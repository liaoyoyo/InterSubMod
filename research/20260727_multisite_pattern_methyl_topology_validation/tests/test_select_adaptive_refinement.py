#!/usr/bin/env python3
"""Tests for pre-gated adaptive permutation selection."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "select_adaptive_refinement.py"
)
SPEC = importlib.util.spec_from_file_location("select_adaptive_refinement", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
selector = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = selector
SPEC.loader.exec_module(selector)


class AdaptiveRefinementSelectionTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.evidence = self.root / "screen.tsv"
        self.base = {
            "dataset": "SYNTH",
            "chrom": "chr1",
            "region_id": "region-1",
            "hp_raw": "1-1",
            "pair_full4": "true",
            "k_ge_3": "false",
            "evaluation_status": "EVALUABLE",
            "permanova_p": "0.001",
            "permanova_r2": "0.2",
            "permanova_permutations_requested": "999",
            "permanova_permutations_realized": "999",
            "permdisp_p": "0.2",
            "best_pair_distance_contrast": "0.15",
            "best_pair_standardized_effect": "0.8",
            "max_geometry_smd": "0.2",
            "all_states_n8": "true",
            "equal_n_retention": "0.8",
            "rarefaction_retention": "0.7",
            "multiplicity_family": "CONFIRMATORY_FULL4_OR_LONG",
        }

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write(self, rows: list[dict[str, str]]) -> None:
        with self.evidence.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(self.base),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

    def test_floor_candidate_passing_every_pre_gate_is_selected(self) -> None:
        self._write([self.base])
        rows, counts = selector.select_rows(
            self.evidence,
            screen_permutations=999,
            refinement_permutations=49999,
        )
        self.assertEqual(counts["selected_for_refinement"], 1)
        self.assertEqual(rows[0]["hp_raw"], "1-1")
        self.assertEqual(rows[0]["refinement_permutations"], 49999)

    def test_any_failed_pre_gate_excludes_candidate(self) -> None:
        failed = dict(self.base)
        failed["max_geometry_smd"] = "0.5"
        self._write([failed])
        rows, counts = selector.select_rows(
            self.evidence,
            screen_permutations=999,
            refinement_permutations=49999,
        )
        self.assertEqual(rows, [])
        self.assertEqual(counts["selected_for_refinement"], 0)

    def test_execute_receipt_binds_screen_and_unit_key_size_and_hash(self) -> None:
        self._write([self.base])
        output_tsv = self.root / "selected.tsv"
        output_json = self.root / "selected.receipt.json"
        selector.execute(
            self.evidence,
            output_tsv,
            output_json,
            screen_permutations=999,
            refinement_permutations=49999,
        )
        receipt = json.loads(output_json.read_text(encoding="utf-8"))
        self.assertEqual(
            receipt["inputs"]["screen_evidence"]["size_bytes"],
            self.evidence.stat().st_size,
        )
        self.assertEqual(
            receipt["inputs"]["screen_evidence"]["sha256"],
            selector.sha256_file(self.evidence),
        )
        self.assertEqual(
            receipt["outputs"]["unit_keys"]["size_bytes"],
            output_tsv.stat().st_size,
        )

    def test_non_frozen_permutation_budgets_fail_closed(self) -> None:
        self._write([self.base])
        with self.assertRaisesRegex(
            selector.RefinementSelectionError, "frozen budget 999"
        ):
            selector.select_rows(
                self.evidence,
                screen_permutations=998,
                refinement_permutations=49999,
            )
        with self.assertRaisesRegex(
            selector.RefinementSelectionError, "frozen budget 49999"
        ):
            selector.select_rows(
                self.evidence,
                screen_permutations=999,
                refinement_permutations=50000,
            )


if __name__ == "__main__":
    unittest.main()
