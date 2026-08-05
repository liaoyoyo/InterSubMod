#!/usr/bin/env python3
"""Tests for the source-backed current sSNV funnel receipt."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "scripts" / "build_current_funnel_receipt.py"
SPEC = importlib.util.spec_from_file_location("current_funnel", SCRIPT)
FUNNEL = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = FUNNEL
SPEC.loader.exec_module(FUNNEL)


def ledger(dataset: str, raw: int, tree: int, autosomal: int, retained: int) -> dict:
    excluded = raw - tree
    out_of_scope = tree - autosomal
    singleton = autosomal - retained
    return {
        "schema_version": "2.0",
        "sample": dataset,
        "pass": True,
        "longphase_input_contract": "clairs_raw_all",
        "tree_contract": "longphase_recalibrated_PASS",
        "raw_clairs_records": raw,
        "longphase_input_records": raw,
        "longphase_recalibrated_records": raw,
        "tree_input_records": tree,
        "autosomal_biallelic_snvs": autosomal,
        "branch_counts": {
            "excluded_by_longphase_filter": excluded,
            "max_snv_excluded": 0,
            "out_of_scope_non_autosomal": out_of_scope,
            "positional_singleton": singleton,
            "retained": retained,
        },
        "checks": {key: True for key in FUNNEL.REQUIRED_CHECKS},
        "duplicate_record_key_excess": {key: 0 for key in ("raw_clairs", "longphase_input", "longphase_recalibrated", "tree_input")},
    }


class CurrentFunnelReceiptTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.ledger_root = self.root / "canonical"
        samples = []
        for index, dataset in enumerate(FUNNEL.DATASETS, start=1):
            raw, tree, autosomal, retained = 100 + index, 90 + index, 80 + index, 50 + index
            directory = self.ledger_root / "samples" / dataset
            directory.mkdir(parents=True)
            (directory / f"ssnv_site_ledger_{dataset}.summary.json").write_text(
                json.dumps(ledger(dataset, raw, tree, autosomal, retained)) + "\n",
                encoding="utf-8",
            )
            samples.append(
                {
                    "sample": dataset,
                    "tree_input_records": tree,
                    "autosomal_biallelic_sSNV": autosomal,
                    "retained_sSNV": retained,
                }
            )
        aggregate = {
            "dataset_count": len(FUNNEL.DATASETS),
            "tree_input_records": sum(row["tree_input_records"] for row in samples),
            "autosomal_biallelic_sSNV": sum(row["autosomal_biallelic_sSNV"] for row in samples),
            "retained_sSNV": sum(row["retained_sSNV"] for row in samples),
        }
        self.canonical = self.root / "canonical.json"
        self.canonical.write_text(
            json.dumps({"all_pass": True, "canonical": {"aggregate": aggregate, "samples": samples}}) + "\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_builds_conserved_source_backed_aggregate(self) -> None:
        receipt = FUNNEL.build_receipt(self.canonical, self.ledger_root)
        self.assertTrue(receipt["all_pass"])
        self.assertEqual(receipt["scope"]["dataset_count"], 7)
        self.assertEqual(receipt["scope"]["biological_sample_count"], 6)
        self.assertEqual(receipt["aggregate"]["raw_clairs_records"], sum(100 + i for i in range(1, 8)))
        self.assertEqual(
            sum(receipt["aggregate"]["branch_counts"].values()),
            receipt["aggregate"]["raw_clairs_records"],
        )
        self.assertAlmostEqual(
            receipt["aggregate"]["relative_ratios"]["retained_over_autosomal"],
            sum(50 + i for i in range(1, 8)) / sum(80 + i for i in range(1, 8)),
        )
        self.assertEqual(len(receipt["inputs"]["site_ledger_summaries"]), 7)

    def test_rejects_failed_ledger_check(self) -> None:
        path = self.ledger_root / "samples" / FUNNEL.DATASETS[0] / f"ssnv_site_ledger_{FUNNEL.DATASETS[0]}.summary.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["checks"]["autosomal_snv_conservation"] = False
        path.write_text(json.dumps(payload) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(FUNNEL.FunnelReceiptError, "autosomal_snv_conservation"):
            FUNNEL.build_receipt(self.canonical, self.ledger_root)

    def test_rejects_canonical_mismatch(self) -> None:
        payload = json.loads(self.canonical.read_text(encoding="utf-8"))
        payload["canonical"]["samples"][0]["retained_sSNV"] += 1
        self.canonical.write_text(json.dumps(payload) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(FUNNEL.FunnelReceiptError, "retained 與 canonical 不符"):
            FUNNEL.build_receipt(self.canonical, self.ledger_root)


if __name__ == "__main__":
    unittest.main()
