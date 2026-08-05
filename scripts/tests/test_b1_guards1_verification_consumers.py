#!/usr/bin/env python3
"""Targeted fixture regressions for the B1-Guards-1 consumers."""

from __future__ import annotations

import ast
import importlib.util
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"
for path in (ANALYSIS_DIR, REPO_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from scripts.lib.verification_schema_contract import (
    SchemaContractError,
    select_current_view,
    select_legacy_view,
)


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ANALYSIS_DIR / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {filename}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_function_from_ast(filename: str, function_name: str, namespace: dict):
    source_path = ANALYSIS_DIR / filename
    tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
    function = next(
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == function_name
    )
    isolated = ast.Module(body=[function], type_ignores=[])
    ast.fix_missing_locations(isolated)
    exec(compile(isolated, str(source_path), "exec"), namespace)
    return namespace[function_name]


CANDIDATE = load_module("b1_candidate_rules", "analyze_candidate_rules.py")
ARI = load_module("b1_ari_typology", "ari_typology.py")
QUADRANT = load_module("b1_loh_quadrant", "build_loh_quadrant_explanation.py")
O03 = load_module("b1_o03", "build_observation_O03_loh_features_post_fix.py")
O15 = load_module("b1_o15", "build_observation_O15_loh_zone_metrics_hcc1395.py")
Q1Q5 = load_module("b1_q1q5", "build_to_feature_study_Q1Q5.py")

EXPANDED_ATTACH = load_function_from_ast(
    "build_loh_expanded_observation.py",
    "attach_historical_verification_view",
    {"pd": pd, "select_legacy_view": select_legacy_view},
)
POST_ATTACH = load_function_from_ast(
    "build_post_hp_fix_to_loh_investigation.py",
    "attach_current_verification_view",
    {
        "pd": pd,
        "SchemaContractError": SchemaContractError,
        "select_current_view": select_current_view,
    },
)


def v2_fixture(include_unknown: bool = False) -> pd.DataFrame:
    current = ["Strong_Bidirectional", "ClusterFirstOnly"]
    if include_unknown:
        current[1] = "UnexpectedClass"
    return pd.DataFrame(
        {
            "Chr": ["chr1", "chr1"],
            "Pos": ["10", "20"],
            "Ref": ["A", "C"],
            "Alt": ["T", "G"],
            "VerificationSchemaVersion": [2, 2],
            "VerificationClass": current,
            "VerificationClass_Legacy": ["Strong", "Subclone"],
            "LOH_Subtype_LegacyVC": ["LOH_Strong", "LOH_Subclone"],
            "LOH_Subtype": ["LOH_Strong", "LOH_Subclone"],
        }
    )


class CurrentViewConsumerTests(unittest.TestCase):
    def test_candidate_summary_keeps_unknown_bucket_and_metadata(self) -> None:
        frame = v2_fixture(include_unknown=True)
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "summary.csv"
            frame.to_csv(path, index=False)
            rows, metadata = CANDIDATE.parse_summary(path, "TP")
        selected = [row["_verification_class_current"] for row in rows.values()]
        self.assertEqual(selected, ["Strong_Bidirectional", "UnknownCurrentClass"])
        self.assertEqual(metadata["schema_status"], "V2")
        self.assertEqual(metadata["unknown_current_class_count"], 1)

    def test_dynamic_consumers_retain_unknown_bucket(self) -> None:
        frame = v2_fixture(include_unknown=True)
        for selector in (
            ARI.attach_current_verification_view,
            QUADRANT.attach_current_verification_view,
        ):
            with self.subTest(selector=selector.__module__):
                selected = selector(frame)
                self.assertEqual(
                    selected["_verification_class_current"].tolist(),
                    ["Strong_Bidirectional", "UnknownCurrentClass"],
                )

    def test_post_fix_alias_and_mixed_version_fail_closed(self) -> None:
        frame = v2_fixture()
        frame["verification_class"] = frame["VerificationClass"]
        selected = POST_ATTACH(frame)
        self.assertEqual(
            selected["_verification_class_current"].tolist(),
            ["Strong_Bidirectional", "ClusterFirstOnly"],
        )
        mismatched = frame.copy()
        mismatched.loc[1, "verification_class"] = "Strong_Bidirectional"
        with self.assertRaisesRegex(ValueError, "not an exact alias"):
            POST_ATTACH(mismatched)
        mixed = frame.copy()
        mixed.loc[1, "VerificationSchemaVersion"] = 1
        with self.assertRaises(ValueError):
            POST_ATTACH(mixed)


class LegacyAndLohConsumerTests(unittest.TestCase):
    def test_versioned_historical_consumers_read_canonical_legacy_fields(self) -> None:
        frame = v2_fixture()
        expanded = EXPANDED_ATTACH(frame)
        o15 = O15.attach_historical_schema_views(frame)
        self.assertEqual(
            expanded["_verification_class_legacy"].tolist(), ["Strong", "Subclone"]
        )
        self.assertEqual(
            o15.attrs["schema_metadata"]["verification_selection_field"],
            "VerificationClass_Legacy",
        )
        self.assertEqual(
            o15.attrs["schema_metadata"]["loh_selection_field"],
            "LOH_Subtype_LegacyVC",
        )

    def test_o03_uses_current_and_canonical_loh(self) -> None:
        selected = O03.attach_verification_loh_views(v2_fixture(include_unknown=True))
        self.assertEqual(
            selected["_verification_class_current"].tolist(),
            ["Strong_Bidirectional", "UnknownCurrentClass"],
        )
        self.assertEqual(
            selected["_loh_subtype_legacy"].tolist(),
            ["LOH_Strong", "LOH_Subclone"],
        )

    def test_historical_l4_requires_explicit_h1(self) -> None:
        historical = pd.DataFrame(
            {
                "VerificationClass": ["Strong", "Noise"],
                "LOH_Subtype": ["LOH_Strong", "LOH_Noise"],
            }
        )
        with self.assertRaises(ValueError):
            EXPANDED_ATTACH(historical)
        with self.assertRaises(ValueError):
            O15.attach_historical_schema_views(historical)
        with warnings.catch_warnings(record=True):
            expanded = EXPANDED_ATTACH(historical, allow_unversioned_v1=True)
            o15 = O15.attach_historical_schema_views(
                historical, allow_unversioned_v1=True
            )
        self.assertEqual(
            expanded["_verification_class_legacy"].tolist(), ["Strong", "Noise"]
        )
        self.assertEqual(
            o15["_loh_subtype_legacy"].tolist(), ["LOH_Strong", "LOH_Noise"]
        )

    def test_loh_alias_mismatch_and_missing_are_distinct(self) -> None:
        mismatch = v2_fixture()
        mismatch.loc[1, "LOH_Subtype"] = "LOH_Weak"
        with self.assertRaisesRegex(ValueError, "not an exact alias"):
            O03.attach_verification_loh_views(mismatch)

        missing = v2_fixture()
        missing.loc[1, "LOH_Subtype_LegacyVC"] = ""
        with self.assertRaisesRegex(ValueError, "<MISSING>"):
            O03.attach_verification_loh_views(missing)

        unknown = v2_fixture()
        unknown.loc[1, "LOH_Subtype_LegacyVC"] = "LOH_Unknown"
        unknown.loc[1, "LOH_Subtype"] = "LOH_Unknown"
        with self.assertRaisesRegex(ValueError, "LOH_Unknown"):
            O03.attach_verification_loh_views(unknown)

    def test_q1q5_explicit_raw_mode_and_canonical_mode(self) -> None:
        canonical = Q1Q5.attach_schema_views(v2_fixture(), "fixture")
        self.assertEqual(
            canonical.attrs["schema_metadata"]["loh_selection_field"],
            "LOH_Subtype_LegacyVC",
        )
        historical = pd.DataFrame(
            {
                "VerificationClass": ["Strong", "Noise"],
                "LOH_Subtype": ["LOH_Strong", "LOH_Noise"],
            }
        )
        with self.assertRaises(ValueError):
            Q1Q5.attach_schema_views(historical, "fixture")
        with warnings.catch_warnings(record=True):
            selected = Q1Q5.attach_schema_views(
                historical, "fixture", allow_unversioned_v1=True
            )
        self.assertEqual(
            selected.attrs["schema_metadata"]["verification_schema_status"],
            "UNVERSIONED",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
