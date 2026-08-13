#!/usr/bin/env python3
"""Regression tests for fail-closed public-claim governance."""

from __future__ import annotations

import copy
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[3]
PROJECT = ROOT / "research/20260813_public_docs_p0_correction"
SCRIPTS = PROJECT / "scripts"
VALIDATOR = SCRIPTS / "validate_claim_remediation_registry.py"
sys.path.insert(0, str(SCRIPTS))

from build_claim_remediation_registry import build  # noqa: E402
from claim_registry_contract import (  # noqa: E402
    C066_CANONICAL_WORDING,
    atomic_write_json,
    build_receipt_payload,
    canonical_json_bytes,
    canonical_object_sha256,
    sha256_bytes,
)
from github_about_evidence import validate_about_receipt  # noqa: E402


class ClaimRemediationRegistryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory()
        cls.temp = Path(cls.temporary.name)
        cls.registry = build()
        cls.registry_path = cls.temp / "registry.json"
        cls.receipt_path = cls.temp / "receipt.json"
        cls.write_pair(cls.registry, cls.registry_path, cls.receipt_path)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    @staticmethod
    def write_pair(registry: dict[str, object], registry_path: Path, receipt_path: Path) -> None:
        raw = canonical_json_bytes(registry)
        atomic_write_json(registry_path, registry)
        atomic_write_json(receipt_path, build_receipt_payload(registry, sha256_bytes(raw)))

    def run_validator(
        self, *args: str, registry: Path | None = None, receipt: Path | None = None
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                "python3",
                str(VALIDATOR),
                "--registry",
                str(registry or self.registry_path),
                "--receipt",
                str(receipt or self.receipt_path),
                *args,
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_00_registry_and_receipt_pass(self) -> None:
        result = self.run_validator()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('"public_source_files": 34', result.stdout)
        self.assertIn('"github_about_snapshot_files": 3', result.stdout)

    def test_missing_and_duplicate_claims_fail(self) -> None:
        self.assertNotEqual(self.run_validator("--simulate-drop-claim", "C001").returncode, 0)
        self.assertNotEqual(self.run_validator("--simulate-duplicate-claim", "C001").returncode, 0)

    def test_claim_semantic_drift_fails(self) -> None:
        for argument in (
            "--simulate-illegal-verdict",
            "--simulate-missing-evidence",
            "--simulate-unverified-upgrade",
            "--simulate-source-live-collapse",
            "--simulate-denominator-drift",
            "--simulate-capability-commit-drift",
        ):
            with self.subTest(argument=argument):
                self.assertNotEqual(self.run_validator(argument).returncode, 0)

    def test_fixed_top_level_contract_drift_fails(self) -> None:
        for argument in (
            "--simulate-schema-drift",
            "--simulate-generated-at-drift",
            "--simulate-task-type-drift",
            "--simulate-source-path-drift",
            "--simulate-allowed-verdicts-drift",
            "--simulate-counts-drift",
            "--simulate-gates-drift",
            "--simulate-p0-replay-drift",
            "--simulate-p0-guard-hash-drift",
            "--simulate-about-receipt-hash-drift",
        ):
            with self.subTest(argument=argument):
                self.assertNotEqual(self.run_validator(argument).returncode, 0)

    def test_recomputed_manifests_reject_synchronized_forgery(self) -> None:
        forged = copy.deepcopy(self.registry)
        forged["public_source_manifest"]["files"][0]["sha256"] = "0" * 64
        forged["public_source_manifest_sha256"] = canonical_object_sha256(forged["public_source_manifest"])
        registry_path = self.temp / "forged-registry.json"
        receipt_path = self.temp / "forged-receipt.json"
        self.write_pair(forged, registry_path, receipt_path)
        result = self.run_validator(registry=registry_path, receipt=receipt_path)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("public source manifest differs", result.stdout)

    def test_about_snapshot_manifest_recomputed(self) -> None:
        forged = copy.deepcopy(self.registry)
        forged["github_about_snapshot_manifest"]["files"][0]["sha256"] = "0" * 64
        forged["github_about_snapshot_manifest_sha256"] = canonical_object_sha256(
            forged["github_about_snapshot_manifest"]
        )
        registry_path = self.temp / "forged-about-registry.json"
        receipt_path = self.temp / "forged-about-receipt.json"
        self.write_pair(forged, registry_path, receipt_path)
        result = self.run_validator(registry=registry_path, receipt=receipt_path)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("snapshot manifest differs", result.stdout)

    def test_build_receipt_semantic_drift_fails(self) -> None:
        receipt = json.loads(self.receipt_path.read_text(encoding="utf-8"))
        receipt["gates"]["PUBLICATION_READY"]["status"] = "PASS"
        bad_receipt = self.temp / "bad-receipt.json"
        atomic_write_json(bad_receipt, receipt)
        result = self.run_validator(receipt=bad_receipt)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not exactly bind", result.stdout)

    def test_about_receipt_rejects_hash_valid_semantic_forgery(self) -> None:
        evidence_root = self.temp / "about-root"
        receipt = json.loads((PROJECT / "github_about_c108_receipt.json").read_text(encoding="utf-8"))
        for command in receipt["verification"]["commands"]:
            source = ROOT / command["response_snapshot"]
            destination = evidence_root / command["response_snapshot"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
            if command["command"] == "gh api -i repos/liaoyoyo/InterSubMod":
                payload = json.loads(destination.read_text(encoding="utf-8"))
                payload["description"] = "hash-valid but semantically false"
                destination.write_text(json.dumps(payload), encoding="utf-8")
                command["response_sha256"] = sha256_bytes(destination.read_bytes())
        receipt_path = evidence_root / "receipt.json"
        receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
        _receipt, errors = validate_about_receipt(evidence_root, receipt_path)
        self.assertTrue(any("snapshot description mismatch" in error for error in errors), errors)

    def test_c066_contract_separates_three_semantics(self) -> None:
        claim = next(item for item in self.registry["claims"] if item["claim_id"] == "C066")
        self.assertEqual(claim["remediated_claim"], C066_CANONICAL_WORDING)
        for anchor in ("raw Cramer's V", "reliability flag", "summary CramersV", "GlobalTest::passed_gate"):
            self.assertIn(anchor, claim["remediated_claim"])

    def test_sampled_storage_claims_are_bounded_validated_derived(self) -> None:
        by_id = {item["claim_id"]: item for item in self.registry["claims"]}
        for claim_id in ("C047", "C048", "C089"):
            with self.subTest(claim_id=claim_id):
                self.assertEqual(by_id[claim_id]["current_verdict"], "CONFIRMED_WITH_LIMITS")
                self.assertEqual(by_id[claim_id]["source_status"], "VALIDATED_DERIVED_WITH_LIMITS")
                self.assertIn("PARTIAL/NON_FINAL", by_id[claim_id]["remediated_claim"])
        self.assertEqual(
            self.registry["counts"]["by_current_verdict"],
            {"CONFIRMED": 69, "CONFIRMED_WITH_LIMITS": 69, "UNVERIFIED": 20},
        )

    def test_named_eligibility_flag_provenance_is_bounded_and_closed(self) -> None:
        claim = next(item for item in self.registry["claims"] if item["claim_id"] == "C114")
        self.assertEqual(claim["current_verdict"], "CONFIRMED_WITH_LIMITS")
        self.assertEqual(claim["source_status"], "FROZEN_CONFIRMED_WITH_LIMITS")
        for anchor in ("canonical/cohort_receipt.json", "summary/all7_summary.json", "not top-level"):
            self.assertIn(anchor, claim["remediated_claim"])

    def test_atomic_write_preserves_previous_file_on_replace_failure(self) -> None:
        target = self.temp / "atomic.json"
        target.write_bytes(b"old-bytes\n")
        with mock.patch("claim_registry_contract.os.replace", side_effect=OSError("simulated replace failure")):
            with self.assertRaises(OSError):
                atomic_write_json(target, {"new": True})
        self.assertEqual(target.read_bytes(), b"old-bytes\n")
        self.assertFalse(list(self.temp.glob(".atomic.json.*.tmp")))


if __name__ == "__main__":
    unittest.main()
