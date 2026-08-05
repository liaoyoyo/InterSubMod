#!/usr/bin/env python3
"""Tests for focal-ALT analysis orchestration helpers."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location(
    "analyze_focal_alt_multicluster", SCRIPTS / "analyze_focal_alt_multicluster.py"
)
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)


class AnalysisHelperTest(unittest.TestCase):
    def test_chromosome_sort_key_supports_sex_and_mitochondrial_contigs(self) -> None:
        chroms = ["chrM", "chr2", "chrX", "chr1", "chrY", "chr10"]
        self.assertEqual(
            sorted(chroms, key=M.chromosome_sort_key),
            ["chr1", "chr2", "chr10", "chrX", "chrY", "chrM"],
        )

    def test_unknown_contigs_sort_after_canonical_contigs(self) -> None:
        self.assertLess(M.chromosome_sort_key("chrM"), M.chromosome_sort_key("chrUn_KI270442v1"))


if __name__ == "__main__":
    unittest.main()
