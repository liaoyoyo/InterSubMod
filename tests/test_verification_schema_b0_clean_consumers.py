#!/usr/bin/env python3
"""Fixture regressions for the B0-Clean-2 verification schema consumers."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from scripts.analysis.build_cross_sample_methylation_observation_workspace import (  # noqa: E402
    apply_verification_view,
)
from scripts.analysis.build_observation_O06_verification_class import (  # noqa: E402
    apply_verification_contract as apply_o06_contract,
)
from scripts.analysis.build_observation_O10_read_level_methyl import (  # noqa: E402
    attach_legacy_verification_view,
)
from scripts.analysis.build_observation_O15b_loh_zone_metrics_cross_sample import (  # noqa: E402
    apply_historical_schema_contract,
)
from scripts.analysis.build_phase1a_head_to_head_baseline_table import (  # noqa: E402
    build_manifest_stats,
)
from scripts.analysis.fn_verdict_readback_audit import audit_csv  # noqa: E402
from scripts.analysis.verify_class_decision_tree_audit import (  # noqa: E402
    apply_verification_contract as apply_tree_contract,
    decision_tree_flow,
    enrich_row,
)
from scripts.lib.verification_schema_contract import (  # noqa: E402
    SchemaContractError,
    UNKNOWN_CURRENT_CLASS,
    UNKNOWN_LEGACY_CLASS,
)


def canonical_loh_columns(values):
    return {
        "LOH_Subtype_LegacyVC": values,
        "LOH_Subtype": values,
    }


def evidence_columns(cluster_values):
    n = len(cluster_values)
    return {
        "LabelFirstSupport": ["false"] * n,
        "ClusterFirstSupport": cluster_values,
        "WithinHPSupport": ["NA"] * n,
        "DispersionWarning": ["NA"] * n,
        "EvidencePath": ["CLUSTER_FIRST_ONLY"] * n,
        "EvidenceDerivation": ["LIVE"] * n,
    }


class ExplicitViewConsumerTest(unittest.TestCase):
    def test_cross_sample_current_unknown_is_retained_and_legacy_is_required(self):
        frame = pd.DataFrame(
            {
                "VerificationSchemaVersion": [2, 2],
                "VerificationClass": ["Strong_Bidirectional", "FutureClass"],
                "VerificationClass_Legacy": ["Strong", "Noise"],
            }
        )
        current, metadata = apply_verification_view(frame, "current")
        self.assertEqual(current["VerificationClass"].astype("object").tolist(), [
            "Strong_Bidirectional",
            UNKNOWN_CURRENT_CLASS,
        ])
        self.assertEqual(metadata["unknown_counts"], {"FutureClass": 1})
        self.assertEqual(current["VerificationClass_SourceValue"].tolist()[-1], "FutureClass")

        legacy, legacy_metadata = apply_verification_view(frame, "legacy")
        self.assertEqual(legacy_metadata["selection_field"], "VerificationClass_Legacy")
        self.assertEqual(legacy["VerificationClass"].astype("object").tolist(), ["Strong", "Noise"])
        with self.assertRaisesRegex(SchemaContractError, "VerificationClass_Legacy is missing"):
            apply_verification_view(frame.drop(columns="VerificationClass_Legacy"), "legacy")

    def test_o06_current_order_and_loh_alias_equality(self):
        frame = pd.DataFrame(
            {
                "VerificationSchemaVersion": [2, 2],
                "VerificationClass": ["ClusterFirstOnly", "FutureClass"],
                "VerificationClass_Legacy": ["Subclone", "Noise"],
                **canonical_loh_columns(["LOH_Subclone", "None"]),
            }
        )
        selected, metadata = apply_o06_contract(frame, "current")
        self.assertEqual(metadata["loh_selection_field"], "LOH_Subtype_LegacyVC")
        self.assertIn(UNKNOWN_CURRENT_CLASS, selected["VerificationClass"].cat.categories)
        self.assertEqual(selected["VerificationClass"].astype("object").iloc[1], UNKNOWN_CURRENT_CLASS)

        mismatch = frame.copy()
        mismatch.loc[1, "LOH_Subtype"] = "LOH_Noise"
        with self.assertRaisesRegex(SchemaContractError, "not an exact alias"):
            apply_o06_contract(mismatch, "current")


class HistoricalConsumerTest(unittest.TestCase):
    def test_o10_h1_requires_flag_and_counts_unknown(self):
        frame = pd.DataFrame({"VerificationClass": ["Strong", "Mystery"]})
        with self.assertRaisesRegex(SchemaContractError, "cannot infer legacy"):
            attach_legacy_verification_view(frame)

        selected, metadata = attach_legacy_verification_view(frame, allow_unversioned_v1=True)
        self.assertEqual(metadata["unknown_counts"], {"Mystery": 1})
        self.assertEqual(metadata["excluded_unknown_rows"], 1)
        self.assertTrue(pd.isna(selected["VerificationClass"].iloc[1]))
        self.assertEqual(selected["VerificationClass_SourceValue"].iloc[1], "Mystery")

    def test_o15_h1_unknown_bucket_and_canonical_loh_guard(self):
        historical = pd.DataFrame(
            {
                "VerificationClass": ["Noise", "Mystery"],
                "LOH_Subtype": ["LOH_Noise", "None"],
            }
        )
        selected, metadata = apply_historical_schema_contract(
            historical,
            allow_unversioned_v1=True,
        )
        self.assertEqual(metadata["unknown_bucket_count"], 1)
        self.assertEqual(selected["VerificationClass"].astype("object").iloc[1], UNKNOWN_LEGACY_CLASS)

        canonical = pd.DataFrame(
            {
                "VerificationSchemaVersion": [2],
                "VerificationClass_Legacy": ["Noise"],
                "VerificationClass": ["Noise_Uniform"],
                "LOH_Subtype_LegacyVC": ["LOH_Noise"],
                "LOH_Subtype": ["None"],
            }
        )
        with self.assertRaisesRegex(SchemaContractError, "not an exact alias"):
            apply_historical_schema_contract(canonical)

    def test_phase1a_explicit_current_does_not_emit_false_legacy_zeroes(self):
        manifest = pd.DataFrame(
            {
                "dataset_id": ["d1", "d1"],
                "platform": ["ONT", "ONT"],
                "mode": ["paired", "paired"],
                "truth_status": ["TP", "TP"],
                "PassedGating": [True, False],
                "Quality_Score": [10.0, 20.0],
                "PairwiseMedianDist": [0.1, 0.2],
                "VerificationSchemaVersion": [2, 2],
                "VerificationClass": ["Strong_Bidirectional", "FutureClass"],
            }
        )
        row = build_manifest_stats(manifest, "current").iloc[0]
        counts = json.loads(row["verification_class_counts_json"])
        self.assertEqual(row["verification_unknown_regions"], 1)
        self.assertEqual(counts[UNKNOWN_CURRENT_CLASS], 1)
        self.assertEqual(row["verification_noise_regions"], "")
        with self.assertRaisesRegex(SchemaContractError, "VerificationClass_Legacy is missing"):
            build_manifest_stats(manifest, "legacy")

    def test_fn_audit_h1_excludes_unknown_from_noise_weak_cohorts(self):
        frame = pd.DataFrame(
            {
                "VerificationClass": ["Noise", "Weak", "Mystery"],
                "PassedGating": [False, False, False],
                "ClusterPermanovaValid": [True, True, True],
                "ClusterPermanovaP": [0.01, 0.02, 0.03],
                "LabelHPPermanovaValid": [False, False, False],
                "LabelHPPermanovaP": [1.0, 1.0, 1.0],
                "LabelAllelePermanovaValid": [False, False, False],
                "LabelAllelePermanovaP": [1.0, 1.0, 1.0],
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "audit.csv"
            frame.to_csv(path, index=False)
            with self.assertRaises(SchemaContractError):
                audit_csv(path)
            result = audit_csv(path, allow_unversioned_v1=True)
        self.assertEqual(result["verification_unknown_excluded_rows"], 1)
        stage2 = result["stage2_verdict_ignored_permanova_FN1_FN10"]
        self.assertEqual(stage2["Noise"]["class_total"], 1)
        self.assertEqual(stage2["Weak"]["class_total"], 1)


class EvidenceAndPairSchemaTest(unittest.TestCase):
    def test_decision_tree_uses_typed_evidence_not_class_name(self):
        frame = pd.DataFrame(
            {
                "VerificationSchemaVersion": [2, 2],
                "VerificationClass": ["ClusterFirstOnly", "Noise_Uniform"],
                **evidence_columns(["false", "true"]),
            }
        )
        selected, _ = apply_tree_contract(frame)
        rows = [enrich_row(row) for row in selected.to_dict("records")]
        self.assertFalse(rows[0]["_cluster_sig"])
        self.assertTrue(rows[1]["_cluster_sig"])
        flow = decision_tree_flow(rows, [])
        self.assertEqual(flow["cluster_first_support_true"]["tp"], 1)
        self.assertEqual(flow["verification_classes"]["ClusterFirstOnly"]["tp"], 1)

        invalid = frame.copy()
        invalid.loc[0, "VerificationClass"] = "FutureClass"
        with self.assertRaisesRegex(SchemaContractError, "unknown VerificationClass"):
            apply_tree_contract(invalid)

    def test_s5_schema_check_rejects_mixed_versions(self):
        script = REPO_ROOT / "scripts" / "analysis" / "build_s5_downstream_impact_and_cv.py"
        valid = pd.DataFrame(
            {
                "mode": ["paired", "to"],
                "VerificationSchemaVersion": [2, 2],
                "VerificationClass": ["Strong_Bidirectional", "FutureClass"],
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            valid_path = Path(tmpdir) / "valid.tsv"
            invalid_path = Path(tmpdir) / "invalid.tsv"
            valid.to_csv(valid_path, sep="\t", index=False)
            valid.assign(VerificationSchemaVersion=[2, 1]).to_csv(
                invalid_path,
                sep="\t",
                index=False,
            )
            ok = subprocess.run(
                [sys.executable, str(script), "--schema-check-only", str(valid_path)],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            bad = subprocess.run(
                [sys.executable, str(script), "--schema-check-only", str(invalid_path)],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
        self.assertEqual(ok.returncode, 0, ok.stderr)
        self.assertIn(UNKNOWN_CURRENT_CLASS, ok.stdout)
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("expected VerificationSchemaVersion=2", bad.stderr)


if __name__ == "__main__":
    unittest.main()
