#!/usr/bin/env python3
"""Regression tests for schema guards in historical untracked consumers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.analysis import tp_fp_structure_label_association as association
from scripts.lib.verification_schema_contract import (
    UNKNOWN_LEGACY_CLASS,
    SchemaContractError,
)


def fixture_frame(classes, explicit_legacy=True):
    rows = []
    for index, legacy_class in enumerate(classes):
        row = {name: "false" for name in association.BOOL_COLS}
        row.update(
            {
                "RegionID": index,
                "VerificationClass": legacy_class,
                "GlobalP": 1.0,
                "HPFineNGroups": 1,
            }
        )
        if explicit_legacy:
            row["VerificationClass_Legacy"] = legacy_class
        rows.append(row)
    return pd.DataFrame(rows)


class HistoricalAssociationSchemaTest(unittest.TestCase):
    def write_pair(self, root, frame):
        tp = Path(root) / "tp.csv"
        fp = Path(root) / "fp.csv"
        frame.to_csv(tp, index=False)
        frame.to_csv(fp, index=False)
        return tp, fp

    def test_explicit_legacy_field_is_selected(self):
        with tempfile.TemporaryDirectory() as root:
            tp, fp = self.write_pair(root, fixture_frame(["Strong", "Noise"]))
            loaded, view = association.load(tp, fp)

        self.assertEqual(view.field, "VerificationClass_Legacy")
        self.assertEqual(
            loaded[association.VERIFICATION_SELECTED].tolist(),
            ["Strong", "Noise", "Strong", "Noise"],
        )

    def test_unversioned_fallback_requires_flag_and_preserves_unknown_bucket(self):
        frame = fixture_frame(["Strong", "FutureLegacy"], explicit_legacy=False)
        with tempfile.TemporaryDirectory() as root:
            tp, fp = self.write_pair(root, frame)
            with self.assertRaises(SchemaContractError):
                association.load(tp, fp)
            loaded, view = association.load(tp, fp, allow_unversioned_v1=True)

        self.assertEqual(view.schema_status, "UNVERSIONED_V1")
        self.assertEqual(view.unknown_counts, {"FutureLegacy": 2})
        self.assertEqual(
            loaded[association.VERIFICATION_SELECTED].tolist(),
            ["Strong", UNKNOWN_LEGACY_CLASS, "Strong", UNKNOWN_LEGACY_CLASS],
        )


if __name__ == "__main__":
    unittest.main()
