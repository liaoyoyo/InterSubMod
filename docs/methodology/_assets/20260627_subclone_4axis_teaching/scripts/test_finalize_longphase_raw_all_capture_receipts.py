#!/usr/bin/env python3
"""Tests for restart-safe verification of an existing raw-all receipt closeout."""

from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

import build_longphase_raw_all_capture_receipt_v2 as builder
import finalize_longphase_raw_all_capture_receipts as finalizer
import test_build_longphase_raw_all_capture_receipt_v2 as receipt_fixture
import validate_layered_v3_inputs as contract


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


class ExistingCloseoutTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="raw-all-existing-closeout-")
        root = Path(self.temporary.name)
        fixture = receipt_fixture.make_fixture(root)
        template = builder.build_receipt(
            fixture["sample_dir"], fixture["manifest"], fixture["sample"], fixture["run_root"],
            fixture["binary"], ["normalizer", "--sample", fixture["sample"]],
        )
        self.run_root = fixture["run_root"]
        rows = []
        for sample in finalizer.EXPECTED_SAMPLES:
            receipt = copy.deepcopy(template)
            receipt["sample"] = sample
            path = self.run_root / "samples" / sample / "producer_capture_receipt_v2.json"
            write_json(path, receipt)
            rows.append({
                "sample": sample,
                "receipt": contract.artifact(path),
                "mapped_alignment_count": receipt["global_coordinate_counts"]["mapped_alignment_count"],
                "input_record_count": receipt["filter_transition_summary"]["input_record_count"],
                "rescued_nonpass_to_pass": receipt["filter_transition_summary"]["rescued_nonpass_to_pass"],
                "removed_pass_to_nonpass": receipt["filter_transition_summary"]["removed_pass_to_nonpass"],
                "persisted_bam": False,
            })
        self.closeout = self.run_root / "raw_all_receipt_closeout.json"
        write_json(self.closeout, {
            "schema_name": "intersubmod.longphase_raw_all_receipt_closeout",
            "schema_version": "1.0.0",
            "run_root": str(self.run_root),
            "dataset_count": 7,
            "all_pass": True,
            "builder": contract.artifact(Path(builder.__file__).resolve()),
            "receipt_schema": contract.artifact(builder.RECEIPT_SCHEMA.resolve()),
            "receipts": rows,
        })
        self.marker = self.run_root / "_RAW_ALL_RECEIPTS_SUCCESS"
        write_json(self.marker, {
            "status": "SUCCESS",
            "closeout": str(self.closeout),
            "closeout_sha256": contract.file_sha256(self.closeout),
            "dataset_count": 7,
        })

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_existing_closeout_verifies_without_sidecar_rescan(self) -> None:
        observed = finalizer.verify_existing_closeout(self.run_root, self.closeout, self.marker)
        self.assertTrue(observed["all_pass"])
        self.assertEqual(observed["dataset_count"], 7)

    def test_receipt_mutation_is_rejected(self) -> None:
        receipt_path = Path(json.loads(self.closeout.read_text())["receipts"][0]["receipt"]["path"])
        receipt_path.write_text(receipt_path.read_text() + " ", encoding="utf-8")
        with self.assertRaises(contract.ContractError) as caught:
            finalizer.verify_existing_closeout(self.run_root, self.closeout, self.marker)
        self.assertEqual(caught.exception.code, "E_HASH_MISMATCH")


if __name__ == "__main__":
    unittest.main(verbosity=2)
