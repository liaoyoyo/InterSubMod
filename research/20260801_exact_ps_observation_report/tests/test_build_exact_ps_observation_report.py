#!/usr/bin/env python3
"""Contract tests for the authority-driven Exact-PS HTML builder."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[3]
BUILDER = REPO_ROOT / "scripts/analysis/build_exact_ps_observation_report.py"
HANDOFF = REPO_ROOT / "docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01"
MANIFEST = HANDOFF / "authority_manifest.json"
REGISTRY = HANDOFF / "denominator_registry.tsv"


class BuilderContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def run_builder(
        self,
        output: Path,
        *,
        manifest: Path = MANIFEST,
        registry: Path = REGISTRY,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                "--authority-manifest",
                str(manifest),
                "--denominator-registry",
                str(registry),
                "--output-dir",
                str(output),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_authority_build_emits_bound_report_and_receipt(self) -> None:
        output = self.root / "release"
        completed = self.run_builder(output)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        data_path = output / "report_data.json"
        html_path = output / "report.standalone.html"
        receipt_path = output / "report_build_receipt.json"
        self.assertTrue(data_path.is_file())
        self.assertTrue(html_path.is_file())
        self.assertTrue(receipt_path.is_file())
        data = json.loads(data_path.read_text())
        receipt = json.loads(receipt_path.read_text())
        self.assertEqual(
            data["schema_name"],
            "intersubmod.exact_ps_observation_report.data",
        )
        self.assertEqual(data["topology"]["one_exact_topology"], 63506)
        self.assertEqual(
            sum(data["topology"]["resolved_coarse_class"].values()),
            63511,
        )
        self.assertEqual(data["provenance"]["authority_artifact_count"], 13)
        self.assertEqual(data["provenance"]["strict_nested_bundle_count"], 9)
        self.assertTrue(receipt["all_pass"])
        for filename in (
            "report_data.json",
            "report.standalone.html",
            "report_build_receipt.json",
        ):
            sidecar = output / f"{filename}.sha256"
            target = output / filename
            self.assertTrue(sidecar.is_file())
            self.assertEqual(
                sidecar.read_text().split()[0],
                hashlib.sha256(target.read_bytes()).hexdigest(),
            )

    def test_rejects_authority_artifact_hash_drift(self) -> None:
        payload = json.loads(MANIFEST.read_text())
        payload["artifacts"][0]["sha256"] = "0" * 64
        manifest = self.root / "tampered_manifest.json"
        manifest.write_text(json.dumps(payload) + "\n")
        completed = self.run_builder(self.root / "release", manifest=manifest)
        self.assertEqual(completed.returncode, 2)
        self.assertIn("authority SHA mismatch", completed.stderr)

    def test_rejects_duplicate_authority_artifact_id(self) -> None:
        payload = json.loads(MANIFEST.read_text())
        payload["artifacts"][1]["artifact_id"] = payload["artifacts"][0][
            "artifact_id"
        ]
        manifest = self.root / "duplicate_manifest.json"
        manifest.write_text(json.dumps(payload) + "\n")
        completed = self.run_builder(self.root / "release", manifest=manifest)
        self.assertEqual(completed.returncode, 2)
        self.assertIn("duplicate authority artifact", completed.stderr)

    def test_rejects_denominator_registry_row_loss(self) -> None:
        rows = REGISTRY.read_text().splitlines()
        registry = self.root / "short_registry.tsv"
        registry.write_text("\n".join(rows[:-1]) + "\n")
        completed = self.run_builder(self.root / "release", registry=registry)
        self.assertEqual(completed.returncode, 2)
        self.assertIn("denominator row count drift", completed.stderr)

    def test_refuses_nonempty_output_directory(self) -> None:
        output = self.root / "release"
        output.mkdir()
        (output / "user_file.txt").write_text("preserve me\n")
        completed = self.run_builder(output)
        self.assertEqual(completed.returncode, 2)
        self.assertIn("output directory is not empty", completed.stderr)
        self.assertEqual((output / "user_file.txt").read_text(), "preserve me\n")


if __name__ == "__main__":
    unittest.main()
