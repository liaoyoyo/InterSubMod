#!/usr/bin/env python3
"""Tests for canonical LongPhase-S tagged BAM immutability receipts."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parent / "audit_canonical_longphase_bam_immutability.py"
SAMPLES = (
    "HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954",
)


class CanonicalBamImmutabilityTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory(prefix="canonical-bam-audit-")
        self.root = Path(self.temp.name)
        entries = []
        for index, sample in enumerate(SAMPLES):
            bam = self.root / f"{sample}.bam"
            bam.write_bytes((f"{sample}:{index}\n" * 32).encode())
            entries.append({"sample": sample, "tumor_bam": str(bam)})
        self.manifest = self.root / "manifest.json"
        self.manifest.write_text(json.dumps({"samples": entries}) + "\n", encoding="utf-8")
        self.baseline = self.root / "baseline.json"

    def tearDown(self) -> None:
        self.temp.cleanup()

    def run_audit(self, output: Path, baseline: Path | None = None) -> subprocess.CompletedProcess[str]:
        argv = [sys.executable, str(SCRIPT), "--manifest", str(self.manifest), "--output", str(output)]
        if baseline is not None:
            argv.extend(["--baseline", str(baseline)])
        return subprocess.run(argv, text=True, capture_output=True, check=False)

    def test_unchanged_seven_bams_pass(self) -> None:
        capture = self.run_audit(self.baseline)
        self.assertEqual(capture.returncode, 0, capture.stderr)
        verification = self.root / "verification.json"
        checked = self.run_audit(verification, self.baseline)
        self.assertEqual(checked.returncode, 0, checked.stderr)
        result = json.loads(verification.read_text(encoding="utf-8"))
        self.assertTrue(result["all_match"])
        self.assertEqual(len(result["comparisons"]), 7)

    def test_mutated_bam_fails_and_preserves_evidence(self) -> None:
        self.assertEqual(self.run_audit(self.baseline).returncode, 0)
        target = self.root / "COLO829.bam"
        target.write_bytes(target.read_bytes() + b"mutation\n")
        verification = self.root / "mutation.json"
        checked = self.run_audit(verification, self.baseline)
        self.assertEqual(checked.returncode, 7)
        result = json.loads(verification.read_text(encoding="utf-8"))
        self.assertFalse(result["all_match"])
        changed = next(item for item in result["comparisons"] if item["sample"] == "COLO829")
        self.assertIn("size_bytes", changed["differing_fields"])

    def test_existing_output_is_not_overwritten(self) -> None:
        self.baseline.write_text("sentinel\n", encoding="utf-8")
        checked = self.run_audit(self.baseline)
        self.assertEqual(checked.returncode, 7)
        self.assertEqual(self.baseline.read_text(encoding="utf-8"), "sentinel\n")


if __name__ == "__main__":
    unittest.main()
