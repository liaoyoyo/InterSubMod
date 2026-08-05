#!/usr/bin/env python3
"""Negative browser-contract tests for the Exact-PS observation report."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[3]
BUILDER = REPO_ROOT / "scripts/analysis/build_exact_ps_observation_report.py"
QA = (
    REPO_ROOT
    / "research/20260801_exact_ps_observation_report/scripts/"
    "qa_exact_ps_observation_report.py"
)
HANDOFF = REPO_ROOT / "docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01"
MANIFEST = HANDOFF / "authority_manifest.json"
REGISTRY = HANDOFF / "denominator_registry.tsv"


class QaNegativeContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tempdir = tempfile.TemporaryDirectory()
        cls.root = Path(cls.tempdir.name)
        cls.source = cls.root / "source"
        completed = subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                "--authority-manifest",
                str(MANIFEST),
                "--denominator-registry",
                str(REGISTRY),
                "--output-dir",
                str(cls.source),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tempdir.cleanup()

    def run_qa(
        self,
        html: Path,
        report_data: Path,
        output: Path,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(QA),
                "--html",
                str(html),
                "--report-data",
                str(report_data),
                "--output-dir",
                str(output),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_rejects_dom_metric_that_disagrees_with_report_data(self) -> None:
        source_html = (self.source / "report.standalone.html").read_text()
        old = 'data-metric="final_groups" data-value="98955"'
        self.assertIn(old, source_html)
        tampered_html = self.root / "tampered_metric.html"
        tampered_html.write_text(source_html.replace(old, old.replace("98955", "98954"), 1))
        output = self.root / "qa_metric"
        completed = self.run_qa(
            tampered_html,
            self.source / "report_data.json",
            output,
        )
        self.assertEqual(completed.returncode, 1)
        receipt = json.loads((output / "browser_qa_receipt.json").read_text())
        self.assertFalse(receipt["all_pass"])
        self.assertIn("metric DOM/data mismatch", receipt["failure"])

    def test_rejects_report_data_schema_version_drift(self) -> None:
        payload = json.loads((self.source / "report_data.json").read_text())
        payload["schema_version"] = "9.9.9"
        tampered_data = self.root / "tampered_report_data.json"
        tampered_data.write_text(json.dumps(payload) + "\n")
        output = self.root / "qa_schema"
        completed = self.run_qa(
            self.source / "report.standalone.html",
            tampered_data,
            output,
        )
        self.assertEqual(completed.returncode, 1)
        receipt = json.loads((output / "browser_qa_receipt.json").read_text())
        self.assertFalse(receipt["all_pass"])
        self.assertIn("report-data version mismatch", receipt["failure"])


if __name__ == "__main__":
    unittest.main()
