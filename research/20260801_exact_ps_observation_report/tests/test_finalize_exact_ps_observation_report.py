#!/usr/bin/env python3
"""Contract tests for the Exact-PS observation report finalizer."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[3]
FINALIZER = (
    Path(__file__).resolve().parents[1]
    / "scripts/finalize_exact_ps_observation_report.py"
)
BUILDER = REPO_ROOT / "scripts/analysis/build_exact_ps_observation_report.py"
HANDOFF = (
    REPO_ROOT
    / "docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01"
)
MANIFEST = HANDOFF / "authority_manifest.json"
REGISTRY = HANDOFF / "denominator_registry.tsv"


FAKE_QA_PASS = r"""
import argparse, hashlib, json
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument("--html", type=Path, required=True)
p.add_argument("--report-data", type=Path, required=True)
p.add_argument("--output-dir", type=Path, required=True)
p.add_argument("--chromium")
a = p.parse_args()
a.output_dir.mkdir(parents=True, exist_ok=True)
for name, payload in {
    "desktop.png": b"PNG-desktop",
    "laptop.png": b"PNG-laptop",
    "mobile.png": b"PNG-mobile",
    "narrow.png": b"PNG-narrow",
    "no_js.png": b"PNG-no-js",
    "report_A4.pdf": b"%PDF-1.7 fixture /Type /Page",
}.items():
    (a.output_dir / name).write_bytes(payload)
def ident(path):
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
checks = {
    "report_data_schema_pass": True,
    "report_data_contract_pass": True,
    "html_data_binding_pass": True,
    "sample_identity_and_binding_pass": True,
    "four_viewports_pass": True,
    "no_js_core_content_pass": True,
    "print_pdf_pass": True,
    "print_page_and_overflow_pass": True,
    "console_errors_zero": True,
    "page_errors_zero": True,
    "external_requests_zero": True,
    "sample_interaction_pass": True,
}
(a.output_dir / "browser_qa_receipt.json").write_text(json.dumps({
    "schema_name": "intersubmod.exact_ps_observation_report.browser_qa",
    "schema_version": "1.0.0",
    "inputs": {"html": ident(a.html), "report_data": ident(a.report_data)},
    "viewports": {
        key: {"observed": {"embedded_data_equal": True}}
        for key in ("desktop", "laptop", "mobile", "narrow")
    },
    "no_js": {"observed": {"sample_rows": 7}},
    "print": {"observed": {"pdf_page_count": 1}},
    "checks": checks,
    "all_pass": True,
}) + "\n")
"""

FAKE_QA_FORGED_SCHEMA = FAKE_QA_PASS.replace(
    "intersubmod.exact_ps_observation_report.browser_qa",
    "fixture.forged.browser_qa",
)

FAKE_QA_FAIL = r"""
raise SystemExit(2)
"""

FAKE_BUILDER_FORGED = r"""
import argparse, json
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument("--authority-manifest", type=Path, required=True)
p.add_argument("--denominator-registry", type=Path, required=True)
p.add_argument("--output-dir", type=Path, required=True)
a = p.parse_args()
a.output_dir.mkdir(parents=True, exist_ok=True)
(a.output_dir / "report_data.json").write_text("{}\n")
(a.output_dir / "report.standalone.html").write_text("<html></html>\n")
(a.output_dir / "report_build_receipt.json").write_text(json.dumps({
    "schema_name": "fixture.forged.build",
    "schema_version": "1.0.0",
    "all_pass": True,
}) + "\n")
"""


class FinalizerContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.qa = self.root / "qa.py"
        self.release_root = self.root / "releases"

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def run_finalizer(
        self,
        release_id: str,
        *,
        builder: Path = BUILDER,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(FINALIZER),
                "--authority-manifest",
                str(MANIFEST),
                "--denominator-registry",
                str(REGISTRY),
                "--release-root",
                str(self.release_root),
                "--release-id",
                release_id,
                "--builder",
                str(builder),
                "--qa-script",
                str(self.qa),
                "--allow-nonproduction-tools",
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_success_publishes_release_marker_last(self) -> None:
        self.qa.write_text(FAKE_QA_PASS)
        completed = self.run_finalizer("fixture_v1")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        release = self.release_root / "fixture_v1"
        receipt_path = release / "release_receipt.json"
        sidecar = release / "release_receipt.json.sha256"
        self.assertTrue(receipt_path.is_file())
        self.assertTrue(sidecar.is_file())
        receipt = json.loads(receipt_path.read_text())
        self.assertTrue(receipt["all_pass"])
        self.assertEqual(
            receipt["statuses"]["scientific_completeness_status"],
            "INCOMPLETE_WITH_ABSTAIN",
        )
        self.assertEqual(
            sidecar.read_text().split()[0],
            hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
        )
        self.assertTrue(
            (self.release_root / ".locks/fixture_v1.lock.json").is_file()
        )

    def test_qa_failure_preserves_failed_attempt_without_release(self) -> None:
        self.qa.write_text(FAKE_QA_FAIL)
        completed = self.run_finalizer("fixture_fail")
        self.assertEqual(completed.returncode, 2)
        self.assertFalse((self.release_root / "fixture_fail").exists())
        attempts = list((self.release_root / "failed_attempts").glob("fixture_fail.*"))
        self.assertEqual(len(attempts), 1)
        self.assertTrue((attempts[0] / "report_data.json").is_file())
        self.assertTrue((attempts[0] / "failure_receipt.json").is_file())
        self.assertTrue((attempts[0] / "release_reservation.json").is_file())

    def test_refuses_existing_release(self) -> None:
        self.qa.write_text(FAKE_QA_PASS)
        first = self.run_finalizer("fixture_once")
        self.assertEqual(first.returncode, 0, first.stderr)
        second = self.run_finalizer("fixture_once")
        self.assertEqual(second.returncode, 1)
        self.assertIn("release already exists", second.stderr)

    def test_rejects_forged_qa_receipt_schema(self) -> None:
        self.qa.write_text(FAKE_QA_FORGED_SCHEMA)
        completed = self.run_finalizer("fixture_forged_qa")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("QA receipt schema mismatch", completed.stderr)
        self.assertFalse((self.release_root / "fixture_forged_qa").exists())

    def test_rejects_forged_build_receipt(self) -> None:
        forged_builder = self.root / "forged_builder.py"
        forged_builder.write_text(FAKE_BUILDER_FORGED)
        self.qa.write_text(FAKE_QA_PASS)
        completed = self.run_finalizer(
            "fixture_forged_build",
            builder=forged_builder,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("builder omitted required output", completed.stderr)
        self.assertFalse((self.release_root / "fixture_forged_build").exists())


if __name__ == "__main__":
    unittest.main()
