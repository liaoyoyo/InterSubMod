#!/usr/bin/env python3
"""Scoped regressions for Remaining-Pass-Historical consumers."""

from __future__ import annotations

import ast
import unittest
import warnings
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"

from scripts.lib.verification_schema_contract import (
    SchemaContractError,
    VERIFICATION_PROVENANCE_COLUMNS,
    extract_provenance_frame,
    select_loh_legacy,
)


PASS_FILES = (
    "export_phase1_manifest_shard.py",
    "export_phase1_read_training_table.py",
    "run_phase1a_read_classifier_benchmark.py",
    "run_to_support_feature_diagnostics.py",
)

HISTORICAL_FILES = (
    "20260423_B3_paired_obs18.py",
    "20260423_B5_colo829_s1_fold.py",
    "20260423_B7_loh_noise_signal.py",
    "20260423_s5_loh_af_cn_scatter.py",
)


def load_functions_from_ast(filename: str, names: tuple[str, ...], namespace: dict):
    source_path = ANALYSIS_DIR / filename
    tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
    selected = [
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in names
    ]
    if {node.name for node in selected} != set(names):
        raise AssertionError(f"Unable to find {names!r} in {filename}")
    future = ast.ImportFrom(
        module="__future__",
        names=[ast.alias(name="annotations")],
        level=0,
    )
    isolated = ast.Module(body=[future, *selected], type_ignores=[])
    ast.fix_missing_locations(isolated)
    exec(compile(isolated, str(source_path), "exec"), namespace)
    return tuple(namespace[name] for name in names)


def v2_provenance_fixture() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "VerificationSchemaVersion": [2, 2],
            "VerificationClass": ["Strong_Bidirectional", "ClusterFirstOnly"],
            "VerificationClass_V1_Deprecated": ["Strong", "Strong"],
            "VerificationClass_Legacy": ["Strong", "Subclone"],
            "LabelFirstSupport": ["true", "false"],
            "ClusterFirstSupport": ["true", "true"],
            "WithinHPSupport": ["true", "false"],
            "DispersionWarning": ["false", "true"],
            "EvidencePath": ["BIDIRECTIONAL", "CLUSTER_FIRST_ONLY"],
            "EvidenceDerivation": ["LIVE", "LIVE"],
            "LOH_Subtype_LegacyVC": ["LOH_Strong", "LOH_Subclone"],
            "LOH_Subtype": ["LOH_Strong", "LOH_Subclone"],
        }
    )


class ProvenancePassThroughTests(unittest.TestCase):
    def test_all_four_pass_consumers_require_consistent_v2(self) -> None:
        for filename in PASS_FILES:
            (attach,) = load_functions_from_ast(
                filename,
                ("attach_verification_provenance",),
                {
                    "pd": pd,
                    "VERIFICATION_PROVENANCE_COLUMNS": VERIFICATION_PROVENANCE_COLUMNS,
                    "extract_provenance_frame": extract_provenance_frame,
                },
            )
            with self.subTest(filename=filename, mode="valid"):
                selected = attach(v2_provenance_fixture())
                self.assertTrue(
                    set(VERIFICATION_PROVENANCE_COLUMNS).issubset(selected.columns)
                )
                self.assertEqual(
                    selected.attrs["verification_provenance"]["schema_status"],
                    "V2",
                )
            with self.subTest(filename=filename, mode="mixed_version"):
                mixed = v2_provenance_fixture()
                mixed.loc[1, "VerificationSchemaVersion"] = 1
                with self.assertRaises(SchemaContractError):
                    attach(mixed)
            with self.subTest(filename=filename, mode="missing_evidence"):
                missing = v2_provenance_fixture().drop(columns=["EvidencePath"])
                with self.assertRaises(SchemaContractError):
                    attach(missing)

    def test_prediction_rows_retain_every_provenance_field(self) -> None:
        (build_prediction_rows,) = load_functions_from_ast(
            "run_phase1a_read_classifier_benchmark.py",
            ("build_prediction_rows",),
            {"VERIFICATION_PROVENANCE_COLUMNS": VERIFICATION_PROVENANCE_COLUMNS},
        )
        frame = v2_provenance_fixture().iloc[[0]].assign(
            dataset_id="fixture",
            dataset_label="Fixture",
            dataset_role="discovery",
            harmonization_group="ONT|paired",
            region_key="chr1:10:A:T",
            read_id="read1",
            phase1a_read_label="ALT",
            truth_status="TP",
            PassedGating="true",
        )
        rows = build_prediction_rows("fixture", "holdout", frame, [1])
        self.assertEqual(len(rows), 1)
        for column in VERIFICATION_PROVENANCE_COLUMNS:
            self.assertEqual(rows[0][column], frame.iloc[0][column])


class HistoricalLohTests(unittest.TestCase):
    def test_all_four_historical_consumers_require_explicit_h1(self) -> None:
        h1 = pd.DataFrame({"LOH_Subtype": ["LOH_Strong", "None"]})
        for filename in HISTORICAL_FILES:
            (attach,) = load_functions_from_ast(
                filename,
                ("attach_loh_legacy_view",),
                {"pd": pd, "select_loh_legacy": select_loh_legacy},
            )
            with self.subTest(filename=filename, mode="canonical"):
                selected = attach(v2_provenance_fixture())
                self.assertEqual(
                    selected["_loh_subtype_legacy"].tolist(),
                    ["LOH_Strong", "LOH_Subclone"],
                )
            with self.subTest(filename=filename, mode="h1"):
                with self.assertRaises(SchemaContractError):
                    attach(h1)
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    selected = attach(h1, allow_unversioned_v1=True)
                self.assertTrue(caught)
                self.assertEqual(
                    selected.attrs["loh_schema_contract"]["schema_status"],
                    "UNVERSIONED_V1",
                )

    def test_historical_missing_unknown_and_alias_mismatch_fail(self) -> None:
        (attach,) = load_functions_from_ast(
            "20260423_B5_colo829_s1_fold.py",
            ("attach_loh_legacy_view",),
            {"pd": pd, "select_loh_legacy": select_loh_legacy},
        )
        missing = v2_provenance_fixture()
        missing.loc[1, "LOH_Subtype_LegacyVC"] = ""
        missing.loc[1, "LOH_Subtype"] = ""
        with self.assertRaisesRegex(SchemaContractError, "<MISSING>"):
            attach(missing)

        unknown = v2_provenance_fixture()
        unknown.loc[1, "LOH_Subtype_LegacyVC"] = "LOH_Unknown"
        unknown.loc[1, "LOH_Subtype"] = "LOH_Unknown"
        with self.assertRaisesRegex(SchemaContractError, "LOH_Unknown"):
            attach(unknown)

        mismatch = v2_provenance_fixture()
        mismatch.loc[1, "LOH_Subtype"] = "LOH_Weak"
        with self.assertRaisesRegex(SchemaContractError, "not an exact alias"):
            attach(mismatch)


class WiringTests(unittest.TestCase):
    def test_full_provenance_and_no_deprecated_loh_consumers(self) -> None:
        for filename in PASS_FILES:
            source = (ANALYSIS_DIR / filename).read_text(encoding="utf-8")
            with self.subTest(filename=filename):
                self.assertIn("VERIFICATION_PROVENANCE_COLUMNS", source)
                self.assertIn("extract_provenance_frame", source)
        for filename in HISTORICAL_FILES:
            source = (ANALYSIS_DIR / filename).read_text(encoding="utf-8")
            with self.subTest(filename=filename):
                self.assertNotRegex(source, r"\[[\"']LOH_Subtype[\"']\]")
                self.assertNotRegex(source, r"fillna\([\"']None[\"']\)")


if __name__ == "__main__":
    unittest.main(verbosity=2)
