#!/usr/bin/env python3
"""Fail-closed tests for the normalized-raw-all LongPhase-S producer."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[5]
METHOD_ROOT = Path(__file__).resolve().parent.parent
RUNNER = Path(__file__).resolve().parent / "run_longphase_raw_all_production_sidecars.sh"
MANIFEST = METHOD_ROOT / "data/longphase_raw_all_production_manifest_v2.json"
LEGACY_MANIFEST = METHOD_ROOT / "data/longphase_production_sidecar_manifest.json"
PATCH_RECEIPT = ROOT / "research/20260710_layered_reconstruction_v2/longphase_s_zero_read_patch_build_receipt.json"
EXPECTED_BINARY_SHA256 = "5ceba723d31c52f01202478de19952219371f0b7f9c136b8882cdd93026789b8"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class RawAllProductionContractTest(unittest.TestCase):
    def test_manifest_binds_exact_raw_all_universe_and_binary(self):
        value = json.loads(MANIFEST.read_text(encoding="utf-8"))
        self.assertEqual(value["dataset_count"], 7)
        self.assertEqual(value["biological_sample_count"], 6)
        self.assertEqual(value["longphase_binary"]["sha256"], EXPECTED_BINARY_SHA256)
        self.assertEqual(sum(item["caller_raw_scan"]["records"] for item in value["samples"]), 638259)
        self.assertEqual(
            sum(item["caller_raw_scan"]["filter_counts"]["LowQual"] for item in value["samples"]),
            70179,
        )
        self.assertTrue(all(item["caller_raw_scan"]["duplicate_record_key_excess"] == 0
                            for item in value["samples"]))

    def test_patch_receipt_binds_all_probe_evidence(self):
        value = json.loads(PATCH_RECEIPT.read_text(encoding="utf-8"))
        self.assertEqual(value["status"], "APPROVED_FOR_FAIL_CLOSED_7_DATASET_VALIDATION")
        self.assertEqual(value["patched_binary_sha256"], EXPECTED_BINARY_SHA256)
        probes = {
            "hcc1954_raw_all_success_sha256":
                ROOT / "research/20260710_layered_reconstruction_v2/probes/20260711_HCC1954_chr22_raw_all_longphase_contract_v1/raw_all_run/_SUCCESS",
            "hcc1954_original_vs_patch_regression_sha256":
                ROOT / "research/20260710_layered_reconstruction_v2/probes/20260711_HCC1954_chr22_raw_all_longphase_contract_v1/patch_regression_comparison.json",
            "hcc1395_known_crash_patch_success_sha256":
                ROOT / "research/20260710_layered_reconstruction_v2/probes/20260711_HCC1395_chrX_72880028_raw_all_crash_probe_v1/run_patched_v3/_SUCCESS",
            "hcc1395_whole_chrx_success_sha256":
                ROOT / "research/20260710_layered_reconstruction_v2/probes/20260711_HCC1395_chrX_raw_all_patched_contract_v1/run/_SUCCESS",
        }
        for field, path in probes.items():
            self.assertEqual(value["approval"][field], sha256(path))

    def test_legacy_pass_only_manifest_is_rejected_without_run_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_root = Path(tmp) / "must-not-exist"
            env = {**os.environ, "SM_LPS_MANIFEST": str(LEGACY_MANIFEST), "SM_RUN_ROOT": str(run_root)}
            result = subprocess.run(["bash", str(RUNNER)], env=env, text=True, capture_output=True, check=False)
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            self.assertFalse(run_root.exists())

    def test_missing_manifest_is_rejected_without_run_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_root = Path(tmp) / "must-not-exist"
            manifest = Path(tmp) / "missing.json"
            env = {**os.environ, "SM_LPS_MANIFEST": str(manifest), "SM_RUN_ROOT": str(run_root)}
            result = subprocess.run(["bash", str(RUNNER)], env=env, text=True, capture_output=True, check=False)
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            self.assertIn(f"ERROR manifest: {manifest}", result.stderr)
            self.assertFalse(run_root.exists())

    def test_wrong_binary_is_rejected_without_run_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_root = Path(tmp) / "must-not-exist"
            env = {
                **os.environ,
                "SM_LONGPHASE_S": "/bin/true",
                "SM_LPS_MANIFEST": str(MANIFEST),
                "SM_RUN_ROOT": str(run_root),
            }
            result = subprocess.run(["bash", str(RUNNER)], env=env, text=True, capture_output=True, check=False)
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            self.assertFalse(run_root.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
