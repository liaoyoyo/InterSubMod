#!/usr/bin/env python3
"""Focused regression tests for the B0-Clean-1 verification consumers."""

from __future__ import annotations

import importlib.util
import sys
import unittest
import warnings
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ANALYSIS_DIR / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {filename}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


CROSS_SAMPLE = load_module(
    "b0_cross_sample_gradient", "analyze_cross_sample_ism_af_gradient.py"
)
TP_FP = load_module(
    "b0_tp_fp_characterization", "analyze_to_tp_fp_characterization.py"
)
SCHEME = load_module(
    "b0_verification_scheme", "analyze_to_verification_scheme_adjustments.py"
)
ALLELE = load_module("b0_allele_deep_dive", "build_allele_deep_dive.py")
FEATURE = load_module("b0_feature_direction", "build_feature_direction_map.py")
NON_LOH = load_module("b0_non_loh", "build_non_loh_discrimination.py")
PROVENANCE = load_module(
    "b0_to_fp_provenance", "build_to_fp_provenance_analysis.py"
)
METHOD = load_module("b0_method_design", "validate_method_design.py")


def v2_fixture() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "VerificationSchemaVersion": [2, 2, 2, 2],
            "VerificationClass": [
                "Strong_Bidirectional",
                "ClusterFirstOnly",
                "LabelShift",
                "Noise_Uniform",
            ],
            "VerificationClass_V1_Deprecated": [
                "Strong",
                "Strong",
                "LabelShift",
                "Noise_Uniform",
            ],
            "VerificationClass_Legacy": ["Strong", "Subclone", "Weak", "Noise"],
            "LabelFirstSupport": ["true", "false", "true", "false"],
            "ClusterFirstSupport": ["true", "true", "false", "false"],
            "WithinHPSupport": ["false", "false", "false", "NA"],
            "DispersionWarning": ["false", "false", "false", "NA"],
            "EvidencePath": [
                "BIDIRECTIONAL",
                "CLUSTER_FIRST_ONLY",
                "LABEL_SHIFT",
                "NOISE_UNIFORM",
            ],
            "EvidenceDerivation": ["LIVE", "LIVE", "LIVE", "LIVE"],
        }
    )


class HistoricalLegacyConsumerTests(unittest.TestCase):
    def test_gradient_uses_explicit_historical_legacy_view(self) -> None:
        frame = pd.DataFrame({"VerificationClass": ["Strong", "Noise"]})
        with self.assertRaises(ValueError):
            CROSS_SAMPLE.attach_historical_legacy_view(frame)
        with warnings.catch_warnings(record=True):
            selected = CROSS_SAMPLE.attach_historical_legacy_view(
                frame, allow_unversioned_v1=True
            )
        self.assertEqual(
            selected["_verification_legacy_class"].tolist(), ["Strong", "Noise"]
        )
        self.assertEqual(
            selected["_verification_schema_status"].unique().tolist(),
            ["UNVERSIONED_V1"],
        )

    def test_characterization_rejects_nonlegacy_value(self) -> None:
        rows = [{"VerificationClass": "Strong_Bidirectional"}]
        with self.assertRaises(ValueError):
            TP_FP.attach_historical_legacy_records(rows, ["VerificationClass"])

    def test_versioned_input_reads_only_explicit_legacy_column(self) -> None:
        frame = v2_fixture()
        selected = CROSS_SAMPLE.attach_historical_legacy_view(frame)
        self.assertEqual(
            selected["_verification_legacy_class"].tolist(),
            ["Strong", "Subclone", "Weak", "Noise"],
        )
        self.assertEqual(
            selected["_verification_selection_field"].unique().tolist(),
            ["VerificationClass_Legacy"],
        )

    def test_unversioned_legacy_requires_explicit_authorization(self) -> None:
        rows = [{"VerificationClass": "Strong"}, {"VerificationClass": "Noise"}]
        with self.assertRaises(ValueError):
            TP_FP.attach_historical_legacy_records(rows, ["VerificationClass"])
        with warnings.catch_warnings(record=True):
            selected, metadata = TP_FP.attach_historical_legacy_records(
                rows,
                ["VerificationClass"],
                allow_unversioned_v1=True,
            )
        self.assertEqual(
            [row["VerificationClass_Legacy_Selected"] for row in selected],
            ["Strong", "Noise"],
        )
        self.assertEqual(metadata["verification_schema_status"], "UNVERSIONED_V1")


class EvidenceAndReferenceTests(unittest.TestCase):
    def test_scheme_keeps_v1_reference_and_typed_evidence(self) -> None:
        prepared = SCHEME.prepare_verification_frame(v2_fixture())
        self.assertEqual(
            prepared["_verification_class_v1_deprecated"].tolist(),
            ["Strong", "Strong", "LabelShift", "Noise_Uniform"],
        )
        self.assertEqual(
            prepared["_cluster_first_support"].tolist(), [True, True, False, False]
        )
        self.assertEqual(
            prepared["_verification_class_v2"].tolist(),
            [
                "Strong_Bidirectional",
                "ClusterFirstOnly",
                "LabelShift",
                "Noise_Uniform",
            ],
        )

    def test_scheme_missing_evidence_fails_closed(self) -> None:
        frame = v2_fixture().drop(columns=["ClusterFirstSupport"])
        with self.assertRaises(ValueError):
            SCHEME.prepare_verification_frame(frame)


class ExplicitTruthTests(unittest.TestCase):
    def test_all_truth_consumers_reject_verification_only_input(self) -> None:
        frame = pd.DataFrame({"VerificationClass": ["Strong", "Noise"]})
        for module in (ALLELE, FEATURE, NON_LOH):
            with self.subTest(module=module.__name__):
                with self.assertRaises(ValueError):
                    module.require_explicit_truth_labels(frame)

    def test_all_truth_consumers_preserve_explicit_binary_truth(self) -> None:
        frame = pd.DataFrame({"truth_label": ["TP", "FP"]})
        for module in (ALLELE, FEATURE, NON_LOH):
            with self.subTest(module=module.__name__):
                selected = module.require_explicit_truth_labels(frame)
                self.assertEqual(selected["truth_label"].tolist(), ["TP", "FP"])
                self.assertEqual(
                    selected["truth_schema_status"].unique().tolist(),
                    ["EXPLICIT_BINARY_TRUTH"],
                )


class ProvenanceRuleSelectionTests(unittest.TestCase):
    def test_rule_source_is_explicit_and_semantically_distinct(self) -> None:
        frame = v2_fixture()
        bidirectional = PROVENANCE.attach_verification_rule_view(
            frame, "current-bidirectional"
        )
        cluster = PROVENANCE.attach_verification_rule_view(
            frame, "cluster-first-evidence"
        )
        self.assertEqual(
            bidirectional["_verification_rule_match"].tolist(),
            [True, False, False, False],
        )
        self.assertEqual(
            cluster["_verification_rule_match"].tolist(),
            [True, True, False, False],
        )
        bidirectional_ids = {
            rule["rule_id"]
            for rule in PROVENANCE.build_candidate_rules(
                "fixture", "current-bidirectional"
            )
        }
        cluster_ids = {
            rule["rule_id"]
            for rule in PROVENANCE.build_candidate_rules(
                "fixture", "cluster-first-evidence"
            )
        }
        self.assertIn("strong_bidirectional_lowaf_highad", bidirectional_ids)
        self.assertIn("cluster_first_support_lowaf_highad", cluster_ids)
        self.assertNotIn("strong_lowaf_highad", bidirectional_ids | cluster_ids)

    def test_missing_summary_match_fails_closed(self) -> None:
        with self.assertRaises(ValueError):
            PROVENANCE.require_verification_rule_match({"variant_key": "chr1:1:A:T"})


class MethodDesignSelectionTests(unittest.TestCase):
    def test_legacy_and_evidence_modes_are_explicit(self) -> None:
        frame = v2_fixture()
        legacy = METHOD.select_cluster_class_frame(frame, "legacy")
        evidence = METHOD.select_cluster_class_frame(frame, "evidence")
        expected = ["Strong", "Subclone", "Weak", "Noise"]
        self.assertEqual(legacy["_cluster_class_selected"].tolist(), expected)
        self.assertEqual(evidence["_cluster_class_selected"].tolist(), expected)
        self.assertEqual(
            evidence["_cluster_class_selection_field"].unique().tolist(),
            ["LabelFirstSupport+ClusterFirstSupport"],
        )

    def test_unknown_current_class_fails_closed(self) -> None:
        frame = v2_fixture()
        frame.loc[0, "VerificationClass"] = "Strong"
        with self.assertRaises(ValueError):
            METHOD.select_cluster_class_frame(frame, "legacy")


if __name__ == "__main__":
    unittest.main(verbosity=2)
