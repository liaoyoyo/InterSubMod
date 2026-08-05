#!/usr/bin/env python3
"""Data-independent regressions for Batch B0 evidence consumers."""

from __future__ import annotations

import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from scripts.analysis import analyze_clairs_borderline_fn as clairs  # noqa: E402
from scripts.analysis import analyze_gq_methylation_rescue_matrix as gq_matrix  # noqa: E402
from scripts.analysis import analyze_longphase_rescue_with_methylation as longphase  # noqa: E402
from scripts.analysis import analyze_methylation_rescue_feature_space as feature_space  # noqa: E402
from scripts.analysis import build_phase2_annotation_layer as annotation_layer  # noqa: E402
from scripts.analysis import evaluate_rescue_with_methylation as evaluate_rescue  # noqa: E402
from scripts.analysis import fn_verdict_reclassify_t1 as fn_reclassify  # noqa: E402
from scripts.analysis import run_phase2_paired_model_feature_analysis as phase2  # noqa: E402
from scripts.lib.verification_schema_contract import SchemaContractError  # noqa: E402


def v2_fixture() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "VerificationSchemaVersion": [2, 2, 2, 2],
            "VerificationClass": [
                "Strong_Bidirectional",
                "ClusterFirstOnly",
                "LabelShift",
                "Noise_Uncorrelated",
            ],
            "VerificationClass_Legacy": ["Strong", "Subclone", "Weak", "Noise"],
            "LabelFirstSupport": ["true", "false", "true", "false"],
            "ClusterFirstSupport": ["true", "true", "false", "false"],
            "WithinHPSupport": ["false"] * 4,
            "DispersionWarning": ["false"] * 4,
            "EvidencePath": [
                "BIDIRECTIONAL",
                "CLUSTER_FIRST_ONLY",
                "LABEL_SHIFT",
                "NOISE_UNCORRELATED",
            ],
            "EvidenceDerivation": ["LIVE"] * 4,
            "source_scope": ["tp", "fp", "tp", "fp"],
        }
    )


class EvidenceAttachmentTest(unittest.TestCase):
    def test_phase2_v2_and_legacy_sources_are_equivalent(self):
        frame = v2_fixture()
        mask = frame["source_scope"].isin({"tp", "fp"})

        v2 = phase2.attach_verification_support(frame, "evidence-v2", mask)
        legacy = phase2.attach_verification_support(frame, "legacy", mask)

        self.assertEqual(v2["label_first_support"].tolist(), [True, False, True, False])
        self.assertEqual(v2["cluster_first_support"].tolist(), [True, True, False, False])
        self.assertEqual(v2["cluster_first_support"].tolist(), legacy["cluster_first_support"].tolist())
        self.assertEqual(v2["verification_support_source"].iloc[0], "LabelFirstSupport+ClusterFirstSupport")

    def test_phase2_missing_na_unknown_fail_closed(self):
        frame = v2_fixture()
        mask = pd.Series(True, index=frame.index)
        with self.assertRaisesRegex(SchemaContractError, "missing required columns"):
            phase2.attach_verification_support(
                frame.drop(columns="ClusterFirstSupport"),
                "evidence-v2",
                mask,
            )

        missing = frame.copy()
        missing.loc[0, "ClusterFirstSupport"] = "NA"
        with self.assertRaisesRegex(SchemaContractError, "<MISSING>"):
            phase2.attach_verification_support(missing, "evidence-v2", mask)

        unknown = frame.copy()
        unknown.loc[0, "VerificationClass"] = "FutureClass"
        with self.assertRaisesRegex(SchemaContractError, "unknown current classes"):
            phase2.attach_verification_support(unknown, "evidence-v2", mask)

        bad_legacy = frame.copy()
        bad_legacy.loc[0, "VerificationClass_Legacy"] = "Other"
        with self.assertRaisesRegex(SchemaContractError, "invalid values"):
            phase2.attach_verification_support(bad_legacy, "legacy", mask)

    def test_non_analyzed_rows_remain_na_not_false(self):
        frame = v2_fixture()
        mask = pd.Series([True, True, False, False], index=frame.index)

        gq = gq_matrix.attach_verification_support(frame, "evidence-v2", mask)
        states = feature_space.attach_verification_support(frame, "evidence-v2", mask)

        self.assertTrue(pd.isna(gq.loc[2, "cluster_first_support"]))
        self.assertTrue(pd.isna(states.loc[2, "VerificationSupportState"]))
        self.assertEqual(states.loc[:1, "VerificationSupportState"].tolist(), ["Strong", "Subclone"])


class DirectConsumerTest(unittest.TestCase):
    def test_longphase_rules_use_two_evidence_axes(self):
        prepared = longphase.attach_summary_support(v2_fixture(), "evidence-v2")
        rules = longphase.build_rules()

        self.assertTrue(rules["cluster_first_support"](prepared.iloc[0]))
        self.assertTrue(rules["cluster_first_support"](prepared.iloc[1]))
        self.assertTrue(rules["bidirectional_support"](prepared.iloc[0]))
        self.assertFalse(rules["bidirectional_support"](prepared.iloc[1]))
        self.assertTrue(rules["cluster_first_only_support"](prepared.iloc[1]))

    def test_evaluate_rescue_support_reads_cluster_evidence(self):
        prepared = evaluate_rescue.attach_summary_support(v2_fixture(), "evidence-v2")
        self.assertEqual(
            [evaluate_rescue.support_cluster_first(row) for _, row in prepared.iterrows()],
            [True, True, False, False],
        )

    def test_clairs_loader_and_rules_use_selected_evidence(self):
        frame = v2_fixture().assign(
            Chr=["chr1"] * 4,
            Pos=[1, 2, 3, 4],
            Ref=["A"] * 4,
            Alt=["C"] * 4,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "significance.csv"
            frame.to_csv(path, index=False)
            mapping, metadata = clairs.load_significance(path, "evidence-v2")

        self.assertEqual(metadata["support_source"], "LabelFirstSupport+ClusterFirstSupport")
        rows = list(mapping.values())
        rules = clairs._build_combined_rules()
        self.assertEqual(
            [rules["cluster_first_support"][1](row) for row in rows],
            [True, True, False, False],
        )
        self.assertEqual(
            [rules["gq15_and_any_verification_support"][1]({**row, "gq": 15}) for row in rows],
            [True, True, True, False],
        )

    def test_annotation_layer_uses_cluster_first_flag(self):
        frame = pd.DataFrame({"cluster_first_support_flag": [True, False]})
        self.assertEqual(
            annotation_layer.annotation_mask(frame, "cluster_first_support").tolist(),
            [True, False],
        )
        self.assertIn("cluster_first_support_flag", annotation_layer.ANNOTATED_FIELDS)
        for config in annotation_layer.ANNOTATION_CONFIG.values():
            self.assertNotIn("strong_subclone", config.weak_support_flag_features)


class FnReclassifyTest(unittest.TestCase):
    def test_v2_reads_evidence_and_unversioned_requires_flag(self):
        frame = v2_fixture()
        with tempfile.TemporaryDirectory() as tmpdir:
            v2_path = Path(tmpdir) / "v2.csv"
            old_path = Path(tmpdir) / "old.csv"
            old_explicit_path = Path(tmpdir) / "old_explicit.csv"
            frame.to_csv(v2_path, index=False)
            pd.DataFrame({"VerificationClass": ["Strong", "Subclone", "Weak", "Noise"]}).to_csv(
                old_path,
                index=False,
            )
            frame[["VerificationClass", "VerificationClass_Legacy"]].to_csv(
                old_explicit_path,
                index=False,
            )

            rows, metadata = fn_reclassify.load_original_support(v2_path)
            self.assertEqual(metadata["schema_status"], "V2_EVIDENCE")
            self.assertEqual(
                [row["_OriginalClusterFirstSupport"] for row in rows],
                [True, True, False, False],
            )
            with self.assertRaisesRegex(SchemaContractError, "--allow-unversioned-v1"):
                fn_reclassify.load_original_support(old_path)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                old_rows, old_metadata = fn_reclassify.load_original_support(
                    old_path,
                    allow_unversioned_v1=True,
                )
            self.assertEqual(old_metadata["schema_status"], "UNVERSIONED_V1")
            self.assertEqual(
                [row["_OriginalLabelFirstSupport"] for row in old_rows],
                [True, False, True, False],
            )
            self.assertTrue(any("UNVERSIONED" in str(item.message) for item in caught))

            with warnings.catch_warnings(record=True) as explicit_caught:
                warnings.simplefilter("always")
                explicit_rows, explicit_metadata = fn_reclassify.load_original_support(
                    old_explicit_path,
                    allow_unversioned_v1=True,
                )
            self.assertEqual(
                explicit_metadata["schema_status"],
                "UNVERSIONED_V1_EXPLICIT_LEGACY",
            )
            self.assertEqual(explicit_metadata["selection_field"], "VerificationClass_Legacy")
            self.assertEqual(
                [row["_OriginalClusterFirstSupport"] for row in explicit_rows],
                [True, True, False, False],
            )
            self.assertTrue(any("UNVERSIONED" in str(item.message) for item in explicit_caught))

    def test_audit_records_verification_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "v2.csv"
            v2_fixture().to_csv(path, index=False)
            result = fn_reclassify.audit(path)

        self.assertEqual(result["verification_provenance"]["schema_status"], "V2_EVIDENCE")
        self.assertEqual(result["T1i_label_rescue"]["noise_total"], 1)
        self.assertEqual(result["T1ii_full_upgrade_cluster_too"]["weak_total"], 1)


if __name__ == "__main__":
    unittest.main()
