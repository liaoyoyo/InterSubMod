#!/usr/bin/env python3
"""Regression tests for fail-closed LOH provenance in the Phase 1A Round 3 consumer."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from scripts.analysis.run_phase1a_round3_loh_feature_test import (  # noqa: E402
    prepare_loh_features,
)
from scripts.lib.verification_schema_contract import SchemaContractError  # noqa: E402


class Round3LohContractTest(unittest.TestCase):
    def test_canonical_field_is_the_model_feature_and_alias_is_validated(self):
        frame = pd.DataFrame(
            {
                "VerificationSchemaVersion": [2, 2],
                "Potential_LOH": [False, True],
                "LOH_Subtype_LegacyVC": ["None", "LOH_Strong"],
                "LOH_Subtype": ["None", "LOH_Strong"],
            }
        )
        prepared, metadata = prepare_loh_features(frame)
        self.assertEqual(prepared["Potential_LOH"].tolist(), ["False", "True"])
        self.assertEqual(
            prepared["LOH_Subtype_LegacyVC"].tolist(),
            ["None", "LOH_Strong"],
        )
        self.assertEqual(metadata["selection_field"], "LOH_Subtype_LegacyVC")
        self.assertEqual(metadata["schema_status"], "LOH_LEGACY_EXPLICIT")

        mismatched = frame.copy()
        mismatched.loc[1, "LOH_Subtype"] = "LOH_Noise"
        with self.assertRaisesRegex(SchemaContractError, "not an exact alias"):
            prepare_loh_features(mismatched)

    def test_historical_fallback_requires_explicit_authorization(self):
        historical = pd.DataFrame(
            {
                "Potential_LOH": [False, True],
                "LOH_Subtype": ["None", "LOH_Noise"],
            }
        )
        with self.assertRaisesRegex(SchemaContractError, "explicit unversioned fallback"):
            prepare_loh_features(historical)

        with self.assertWarnsRegex(UserWarning, "UNVERSIONED"):
            prepared, metadata = prepare_loh_features(
                historical,
                allow_unversioned_v1=True,
            )
        self.assertEqual(prepared["LOH_Subtype_LegacyVC"].tolist(), ["None", "LOH_Noise"])
        self.assertEqual(metadata["schema_status"], "UNVERSIONED_V1")

    def test_missing_or_inconsistent_values_never_become_none(self):
        missing = pd.DataFrame(
            {
                "Potential_LOH": [False],
                "LOH_Subtype_LegacyVC": [pd.NA],
                "LOH_Subtype": [pd.NA],
            }
        )
        with self.assertRaisesRegex(SchemaContractError, "<MISSING>"):
            prepare_loh_features(missing)

        inconsistent = pd.DataFrame(
            {
                "Potential_LOH": [False],
                "LOH_Subtype_LegacyVC": ["LOH_Weak"],
                "LOH_Subtype": ["LOH_Weak"],
            }
        )
        with self.assertRaisesRegex(SchemaContractError, "disagree"):
            prepare_loh_features(inconsistent)


if __name__ == "__main__":
    unittest.main()
