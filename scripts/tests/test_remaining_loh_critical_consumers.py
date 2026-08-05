#!/usr/bin/env python3
"""Targeted schema regressions for Remaining-LOH-Critical consumers."""

from __future__ import annotations

import ast
import tempfile
import unittest
import warnings
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"

from scripts.lib.verification_schema_contract import (
    SchemaContractError,
    VERIFICATION_PROVENANCE_COLUMNS,
    select_current_view,
    select_legacy_view,
    select_loh_legacy,
)


TARGET_FILES = (
    "20260423_B4_S4_discrimination.py",
    "20260423_phase3_synthesis.py",
    "20260424_X5v2_corrected_S3S5.py",
    "build_allelesig_loh_study.py",
    "build_loh_feature_validation.py",
    "build_to_feature_study_Q6.py",
    "build_multilayer_hp_before_after_comparison.py",
    "20260424_X6_merge_caller_af_S3S5.py",
)

LOH_WRAPPER_FILES = TARGET_FILES[:6]


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


def v2_loh_fixture() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "VerificationSchemaVersion": [2, 2],
            "VerificationClass": ["Strong_Bidirectional", "ClusterFirstOnly"],
            "VerificationClass_Legacy": ["Strong", "Subclone"],
            "LOH_Subtype_LegacyVC": ["LOH_Strong", "LOH_Subclone"],
            "LOH_Subtype": ["LOH_Strong", "LOH_Subclone"],
        }
    )


class LohGuardTests(unittest.TestCase):
    def test_all_six_loh_wrappers_require_explicit_h1(self) -> None:
        h1 = pd.DataFrame({"LOH_Subtype": ["LOH_Strong", "None"]})
        for filename in LOH_WRAPPER_FILES:
            (attach,) = load_functions_from_ast(
                filename,
                ("attach_loh_legacy_view",),
                {"pd": pd, "select_loh_legacy": select_loh_legacy},
            )
            with self.subTest(filename=filename, mode="canonical"):
                selected = attach(v2_loh_fixture())
                self.assertEqual(
                    selected["_loh_subtype_legacy"].tolist(),
                    ["LOH_Strong", "LOH_Subclone"],
                )
                self.assertEqual(
                    selected.attrs["loh_schema_contract"]["selection_field"],
                    "LOH_Subtype_LegacyVC",
                )
            with self.subTest(filename=filename, mode="h1_gate"):
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

    def test_alias_missing_and_unknown_never_become_none(self) -> None:
        (attach,) = load_functions_from_ast(
            "20260423_phase3_synthesis.py",
            ("attach_loh_legacy_view",),
            {"pd": pd, "select_loh_legacy": select_loh_legacy},
        )
        mismatch = v2_loh_fixture()
        mismatch.loc[1, "LOH_Subtype"] = "LOH_Weak"
        with self.assertRaisesRegex(SchemaContractError, "not an exact alias"):
            attach(mismatch)

        missing = v2_loh_fixture()
        missing.loc[1, "LOH_Subtype_LegacyVC"] = ""
        missing.loc[1, "LOH_Subtype"] = ""
        with self.assertRaisesRegex(SchemaContractError, "<MISSING>"):
            attach(missing)

        unknown = v2_loh_fixture()
        unknown.loc[1, "LOH_Subtype_LegacyVC"] = "LOH_Unknown"
        unknown.loc[1, "LOH_Subtype"] = "LOH_Unknown"
        with self.assertRaisesRegex(SchemaContractError, "LOH_Unknown"):
            attach(unknown)


class PassThroughAndVersionTests(unittest.TestCase):
    def test_x6_retains_available_provenance_columns(self) -> None:
        (reader,) = load_functions_from_ast(
            "20260424_X6_merge_caller_af_S3S5.py",
            ("read_summary_with_loh_provenance",),
            {
                "Path": Path,
                "pd": pd,
                "VERIFICATION_PROVENANCE_COLUMNS": VERIFICATION_PROVENANCE_COLUMNS,
                "select_loh_legacy": select_loh_legacy,
            },
        )
        frame = v2_loh_fixture().assign(
            Chr=["chr1", "chr1"],
            Pos=[10, 20],
            Ref=["A", "C"],
            Alt=["T", "G"],
            HPFineNGroups=[2, 3],
            Coverage_Multiple=[1.0, 1.1],
            Diploid_Coverage_Used=[100, 100],
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "summary.csv"
            frame.to_csv(path, index=False)
            selected, metadata = reader(path)
        expected = {
            "VerificationSchemaVersion",
            "VerificationClass",
            "VerificationClass_Legacy",
            "LOH_Subtype_LegacyVC",
            "LOH_Subtype",
        }
        self.assertTrue(expected.issubset(selected.columns))
        self.assertTrue(expected.issubset(metadata["provenance_columns_passthrough"]))
        self.assertEqual(metadata["selection_field"], "LOH_Subtype_LegacyVC")

    def test_multilayer_requires_one_explicit_semantic_view(self) -> None:
        attach, load_and_merge = load_functions_from_ast(
            "build_multilayer_hp_before_after_comparison.py",
            ("attach_verification_view", "load_and_merge"),
            {
                "pd": pd,
                "SchemaContractError": SchemaContractError,
                "select_current_view": select_current_view,
                "select_legacy_view": select_legacy_view,
            },
        )
        h1 = pd.DataFrame({"VerificationClass": ["Strong", "Noise"]})
        with self.assertRaises(SchemaContractError):
            attach(h1, "current-v2")
        with self.assertRaises(SchemaContractError):
            attach(h1, "legacy")
        with warnings.catch_warnings(record=True):
            legacy = attach(h1, "legacy", allow_unversioned_v1=True)
        self.assertEqual(
            legacy["_verification_class_selected"].tolist(), ["Strong", "Noise"]
        )

        before = h1.assign(Chr=["chr1", "chr1"], Pos=[10, 20], Ref=["A", "C"], Alt=["T", "G"])
        after = v2_loh_fixture().assign(
            Chr=["chr1", "chr1"], Pos=[10, 20], Ref=["A", "C"], Alt=["T", "G"]
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            before_path = Path(temp_dir) / "before.csv"
            after_path = Path(temp_dir) / "after.csv"
            before.to_csv(before_path, index=False)
            after.to_csv(after_path, index=False)
            with self.assertRaises(SchemaContractError):
                load_and_merge(before_path, after_path, "fixture", "current-v2")
            with warnings.catch_warnings(record=True):
                merged, _, _, contract = load_and_merge(
                    before_path,
                    after_path,
                    "fixture",
                    "legacy",
                    allow_unversioned_v1=True,
                )
        self.assertIn("_verification_class_selected_before", merged.columns)
        self.assertIn("_verification_class_selected_after", merged.columns)
        self.assertEqual(contract["comparison_semantics"], "Legacy VerificationClass")
        self.assertFalse(contract["direct_v1_v2_current_mix"])


class StaticConsumerWiringTests(unittest.TestCase):
    def test_no_direct_deprecated_loh_or_raw_transition_consumers_remain(self) -> None:
        for filename in TARGET_FILES:
            source = (ANALYSIS_DIR / filename).read_text(encoding="utf-8")
            with self.subTest(filename=filename):
                self.assertNotRegex(source, r"\[[\"']LOH_Subtype[\"']\]")
                self.assertNotRegex(source, r"fillna\([\"']None[\"']\)")
        multilayer = (
            ANALYSIS_DIR / "build_multilayer_hp_before_after_comparison.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("VerificationClass_before", multilayer)
        self.assertNotIn("VerificationClass_after", multilayer)


if __name__ == "__main__":
    unittest.main(verbosity=2)
