#!/usr/bin/env python3
"""Adversarial tests for the immutable professor-report browser QA artifact."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "scripts" / "qa_validated_html_report.py"
SPEC = importlib.util.spec_from_file_location("validated_html_browser_qa", SCRIPT)
QA = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = QA
SPEC.loader.exec_module(QA)


def passing_checks() -> dict[str, bool]:
    return {name: True for name in QA.CHECK_NAMES}


def valid_metrics() -> dict[str, object]:
    return {
        "title": "教授版驗證報告",
        "status_text": "FINAL · VALIDATED",
        "h1_count": 1,
        "section_count": 12,
        "svg_count": 8,
        "details_count": 6,
        "desktop_widths": {"scroll": 1440, "client": 1440},
        "mobile_widths": {"scroll": 390, "client": 390},
        "missing_anchor_targets": [],
        "missing_local_files": [],
        "external_links": [],
        "remote_requests": [],
        "failed_requests": [],
        "console_errors": [],
        "page_errors": [],
    }


class BrowserQaArtifactTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.html = self.root / "report.html"
        self.html_payload = b'<!doctype html><html lang="zh-Hant"><p>exact bytes</p></html>\n'
        self.html.write_bytes(self.html_payload)
        self.producer = self.root / "qa_producer.py"
        self.producer.write_text("# frozen QA producer\n", encoding="utf-8")
        self.counter = 0

    def tearDown(self) -> None:
        self.temp.cleanup()

    def _make_artifact(self) -> Path:
        self.counter += 1
        output = self.root / f"qa_{self.counter}"
        candidate = self.root / f".qa_{self.counter}.staging.test"
        candidate.mkdir()
        (candidate / QA.OUTPUT_FILENAMES["desktop_screenshot"]).write_bytes(b"desktop-png")
        (candidate / QA.OUTPUT_FILENAMES["mobile_screenshot"]).write_bytes(b"mobile-png")
        (candidate / QA.OUTPUT_FILENAMES["print_pdf"]).write_bytes(b"%PDF-test")
        QA.seal_qa_artifact(
            candidate,
            output,
            html_path=self.html,
            producer_path=self.producer,
            expected_status="final",
            checks=passing_checks(),
            metrics=valid_metrics(),
            pre_html_identity=QA.file_identity(self.html),
            pre_producer_identity=QA.artifact_identity(self.producer),
            created_at_utc="2026-07-16T00:00:00Z",
        )
        return output

    def _receipt(self, output: Path) -> dict[str, object]:
        return json.loads((output / QA.RECEIPT_NAME).read_text(encoding="utf-8"))

    def _resign_tampered_receipt(self, output: Path, mutate) -> None:
        os.chmod(output, 0o755)
        receipt_path = output / QA.RECEIPT_NAME
        sidecar_path = output / QA.SIDECAR_NAME
        os.chmod(receipt_path, 0o644)
        document = json.loads(receipt_path.read_text(encoding="utf-8"))
        mutate(document)
        payload = json.dumps(document, ensure_ascii=False, allow_nan=False, indent=2).encode("utf-8") + b"\n"
        receipt_path.write_bytes(payload)
        os.chmod(receipt_path, 0o444)
        os.chmod(sidecar_path, 0o644)
        sidecar_path.write_text(
            f"{hashlib.sha256(payload).hexdigest()}  {QA.RECEIPT_NAME}\n",
            encoding="ascii",
        )
        os.chmod(sidecar_path, 0o444)
        os.chmod(output, 0o555)

    def test_success_binds_all_inputs_outputs_and_exact_modes(self) -> None:
        output = self._make_artifact()
        result = QA.verify_qa_artifact(
            output,
            html_path=self.html,
            producer_path=self.producer,
            expected_status="final",
        )
        receipt = self._receipt(output)

        self.assertTrue(result["all_pass"])
        self.assertEqual(receipt["schema_version"], "2.0.0")
        self.assertEqual(receipt["artifact_class"], QA.ARTIFACT_CLASS)
        self.assertEqual(set(receipt["inputs"]), {"html", "qa_producer"})
        self.assertEqual(set(receipt["outputs"]), set(QA.OUTPUT_FILENAMES))
        for row in (*receipt["inputs"].values(), *receipt["outputs"].values()):
            self.assertEqual(set(row), QA.IDENTITY_KEYS)
            self.assertEqual(len(row["sha256"]), 64)
            self.assertEqual(row["st_nlink"], 1)
        self.assertEqual(receipt["inputs"]["html"]["mode_octal"], "0444")
        self.assertEqual(stat.S_IMODE(self.html.stat().st_mode), 0o444)
        self.assertEqual(stat.S_IMODE(output.stat().st_mode), 0o555)
        self.assertEqual({path.name for path in output.iterdir()}, set(QA.EXACT_QA_FILES))
        for path in output.iterdir():
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o444)
            self.assertEqual(path.stat().st_nlink, 1)

    def test_html_or_rendered_bytes_tamper_is_rejected(self) -> None:
        html_output = self._make_artifact()
        os.chmod(self.html, 0o644)
        self.html.write_bytes(self.html_payload + b"<!-- tamper -->\n")
        os.chmod(self.html, 0o444)
        with self.assertRaisesRegex(QA.QaArtifactError, "HTML identity changed"):
            QA.verify_qa_artifact(html_output, html_path=self.html, producer_path=self.producer)

        # Restore the original HTML for a second independently sealed artifact.
        os.chmod(self.html, 0o644)
        self.html.write_bytes(self.html_payload)
        output = self._make_artifact()
        rendered = output / QA.OUTPUT_FILENAMES["desktop_screenshot"]
        os.chmod(rendered, 0o644)
        rendered.write_bytes(b"tampered screenshot")
        os.chmod(rendered, 0o444)
        with self.assertRaisesRegex(QA.QaArtifactError, "output identity changed"):
            QA.verify_qa_artifact(output, html_path=self.html, producer_path=self.producer)

    def test_receipt_or_sidecar_tamper_is_rejected(self) -> None:
        output = self._make_artifact()
        sidecar = output / QA.SIDECAR_NAME
        os.chmod(sidecar, 0o644)
        sidecar.write_text(f"{'0' * 64}  {QA.RECEIPT_NAME}\n", encoding="ascii")
        os.chmod(sidecar, 0o444)
        with self.assertRaisesRegex(QA.QaArtifactError, "sidecar mismatch"):
            QA.verify_qa_artifact(output, html_path=self.html, producer_path=self.producer)

        output2 = self._make_artifact()
        receipt = output2 / QA.RECEIPT_NAME
        os.chmod(receipt, 0o644)
        receipt.write_bytes(receipt.read_bytes() + b" ")
        os.chmod(receipt, 0o444)
        with self.assertRaisesRegex(QA.QaArtifactError, "sidecar mismatch"):
            QA.verify_qa_artifact(output2, html_path=self.html, producer_path=self.producer)

    def test_mode_and_hardlink_tamper_are_rejected(self) -> None:
        output = self._make_artifact()
        mobile = output / QA.OUTPUT_FILENAMES["mobile_screenshot"]
        os.chmod(mobile, 0o644)
        with self.assertRaisesRegex(QA.QaArtifactError, "mode is not 0444"):
            QA.verify_qa_artifact(output, html_path=self.html, producer_path=self.producer)

        output2 = self._make_artifact()
        pdf = output2 / QA.OUTPUT_FILENAMES["print_pdf"]
        os.link(pdf, self.root / "external_pdf_hardlink")
        with self.assertRaisesRegex(QA.QaArtifactError, "hardlink count is not one"):
            QA.verify_qa_artifact(output2, html_path=self.html, producer_path=self.producer)

    def test_extra_file_and_symlink_are_rejected(self) -> None:
        output = self._make_artifact()
        os.chmod(output, 0o755)
        (output / "unexpected.txt").write_text("extra\n", encoding="utf-8")
        os.chmod(output, 0o555)
        with self.assertRaisesRegex(QA.QaArtifactError, "exact five-file layout"):
            QA.verify_qa_artifact(output, html_path=self.html, producer_path=self.producer)

        output2 = self._make_artifact()
        os.chmod(output2, 0o755)
        mobile = output2 / QA.OUTPUT_FILENAMES["mobile_screenshot"]
        mobile.rename(self.root / "preserved-mobile-output")
        mobile.symlink_to(self.html)
        os.chmod(output2, 0o555)
        with self.assertRaisesRegex(QA.QaArtifactError, "not a regular file"):
            QA.verify_qa_artifact(output2, html_path=self.html, producer_path=self.producer)

    def test_resigned_path_swap_is_rejected_semantically(self) -> None:
        output = self._make_artifact()
        self._resign_tampered_receipt(
            output,
            lambda document: document["outputs"]["desktop_screenshot"].update(
                {"path": str(output / QA.OUTPUT_FILENAMES["mobile_screenshot"])}
            ),
        )
        with self.assertRaisesRegex(QA.QaArtifactError, "output path mismatch"):
            QA.verify_qa_artifact(output, html_path=self.html, producer_path=self.producer)

        output2 = self._make_artifact()
        alternate = self.root / "alternate.html"
        alternate.write_bytes(self.html_payload)
        os.chmod(alternate, 0o444)
        self._resign_tampered_receipt(
            output2,
            lambda document: (
                document["inputs"]["html"].update({"path": str(alternate)}),
                document.update({"html": str(alternate)}),
            ),
        )
        with self.assertRaisesRegex(QA.QaArtifactError, "HTML path differs"):
            QA.verify_qa_artifact(output2, html_path=self.html, producer_path=self.producer)

    def test_producer_drift_is_rejected(self) -> None:
        output = self._make_artifact()
        os.chmod(self.producer, 0o644)
        self.producer.write_text("# changed QA producer\n", encoding="utf-8")
        with self.assertRaisesRegex(QA.QaArtifactError, "producer identity changed"):
            QA.verify_qa_artifact(output, html_path=self.html, producer_path=self.producer)

    def _assert_sealing_identity_race_is_preserved(self, target: str) -> None:
        """Deterministically mutate an input after the last pre-seal check."""

        output = self.root / f"race_{target}"
        original_browser_qa = QA._run_browser_qa
        original_artifact_identity = QA.artifact_identity
        original_module_file = QA.__file__
        victim = self.html if target == "html" else self.producer
        calls = 0

        def fake_browser_qa(_html_path: Path, candidate: Path, _expected_status: str):
            for filename in QA.OUTPUT_FILENAMES.values():
                (candidate / filename).write_bytes(f"rendered {filename}\n".encode("utf-8"))
            return passing_checks(), valid_metrics()

        def racing_artifact_identity(path: Path, *, recorded_path: Path | None = None):
            nonlocal calls
            if Path(os.path.abspath(os.fspath(path))) == Path(os.path.abspath(os.fspath(victim))):
                calls += 1
                # run_browser_qa performs one initial capture, one post-browser
                # check, and one seal-entry check.  The fourth read is the
                # authoritative identity that used to be accepted without
                # comparing it back to the pre-QA identity.
                if calls == 4:
                    os.chmod(victim, 0o644)
                    victim.write_text(f"changed {target} after final check\n", encoding="utf-8")
                    if target == "html":
                        os.chmod(victim, 0o444)
            return original_artifact_identity(path, recorded_path=recorded_path)

        QA._run_browser_qa = fake_browser_qa
        QA.artifact_identity = racing_artifact_identity
        if target == "producer":
            QA.__file__ = str(self.producer)
        try:
            changed_label = "HTML" if target == "html" else "QA producer"
            with self.assertRaisesRegex(
                QA.QaArtifactError,
                f"{changed_label} changed between QA check",
            ):
                QA.run_browser_qa(self.html, output, "final")
        finally:
            QA._run_browser_qa = original_browser_qa
            QA.artifact_identity = original_artifact_identity
            QA.__file__ = original_module_file

        self.assertFalse(output.exists())
        preserved = list(self.root.glob(f"{output.name}.failed-staging.*"))
        self.assertEqual(len(preserved), 1)
        self.assertTrue((preserved[0] / "QA_NOT_PUBLISHED.txt").is_file())
        self.assertFalse((preserved[0] / QA.RECEIPT_NAME).exists())
        self.assertFalse((preserved[0] / QA.SIDECAR_NAME).exists())

    def test_html_sealing_identity_race_is_rejected_and_preserved(self) -> None:
        self._assert_sealing_identity_race_is_preserved("html")

    def test_producer_sealing_identity_race_is_rejected_and_preserved(self) -> None:
        self._assert_sealing_identity_race_is_preserved("producer")

    def test_preexisting_output_and_exclusive_receipt_fail_closed(self) -> None:
        output = self.root / "already_exists"
        output.mkdir()
        candidate = self.root / ".candidate"
        candidate.mkdir()
        for filename in QA.OUTPUT_FILENAMES.values():
            (candidate / filename).write_bytes(filename.encode("ascii"))
        with self.assertRaisesRegex(QA.QaArtifactError, "already exists"):
            QA.seal_qa_artifact(
                candidate,
                output,
                html_path=self.html,
                producer_path=self.producer,
                expected_status="final",
                checks=passing_checks(),
                metrics=valid_metrics(),
                pre_html_identity=QA.file_identity(self.html),
                pre_producer_identity=QA.artifact_identity(self.producer),
            )

        existing = self.root / "exclusive.json"
        QA._write_new_file(existing, b"first\n")
        with self.assertRaises(FileExistsError):
            QA._write_new_file(existing, b"second\n")
        self.assertEqual(existing.read_bytes(), b"first\n")

    def test_failed_staging_is_preserved_without_success_sidecar(self) -> None:
        output = self.root / "never_published"
        candidate = self.root / ".never_published.staging.test"
        candidate.mkdir()
        (candidate / "partial-output.png").write_bytes(b"diagnostic")
        (candidate / QA.RECEIPT_NAME).write_bytes(b'{"all_pass": true}\n')
        (candidate / QA.SIDECAR_NAME).write_bytes(b"staged-but-not-published\n")
        preserved = QA._preserve_failed_staging(candidate, output, "simulated browser failure")

        self.assertTrue(preserved.is_dir())
        self.assertIn("failed-staging", preserved.name)
        self.assertTrue((preserved / "QA_NOT_PUBLISHED.txt").is_file())
        self.assertFalse((preserved / QA.RECEIPT_NAME).exists())
        self.assertFalse((preserved / QA.SIDECAR_NAME).exists())
        self.assertTrue((preserved / f"UNAUTHENTICATED_{QA.RECEIPT_NAME}").is_file())
        self.assertTrue((preserved / f"UNAUTHENTICATED_{QA.SIDECAR_NAME}").is_file())
        self.assertFalse(output.exists())
        self.assertEqual(stat.S_IMODE(preserved.stat().st_mode), 0o555)

    def test_legacy_run_cli_and_verify_only_cli_both_parse(self) -> None:
        legacy = QA.parse_args([
            "--html", str(self.html), "--output-dir", str(self.root / "qa"), "--expect-status", "final",
        ])
        verify = QA.parse_args([
            "verify-only", "--html", str(self.html), "--output-dir", str(self.root / "qa"),
        ])
        explicit_run = QA.parse_args([
            "run", "--html", str(self.html), "--output-dir", str(self.root / "qa"), "--expect-status", "partial",
        ])

        self.assertEqual(legacy.command, "run")
        self.assertEqual(verify.command, "verify-only")
        self.assertEqual(explicit_run.command, "run")
        self.assertIsNone(verify.expect_status)

    def test_legacy_file_identity_remains_strict(self) -> None:
        identity = QA.file_identity(self.html)
        self.assertEqual(identity["path"], str(self.html.resolve()))
        self.assertEqual(identity["size_bytes"], len(self.html_payload))
        self.assertEqual(identity["sha256"], hashlib.sha256(self.html_payload).hexdigest())
        self.assertTrue(QA.file_identity_matches(self.html, identity))
        for malformed in (
            {"path": identity["path"], "sha256": identity["sha256"]},
            {**identity, "size_bytes": True},
            {**identity, "sha256": "0" * 63},
            {**identity, "unexpected": "field"},
        ):
            with self.subTest(malformed=malformed):
                self.assertFalse(QA.file_identity_matches(self.html, malformed))


if __name__ == "__main__":
    unittest.main()
