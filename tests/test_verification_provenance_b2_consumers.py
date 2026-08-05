#!/usr/bin/env python3
"""Focused provenance-contract regressions for the B2 consumer batch."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"
for path in (REPO_ROOT, ANALYSIS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from scripts.analysis import analyze_kism_vs_cn_perread as kism  # noqa: E402
from scripts.analysis import analyze_to_feature_recall as to_feature  # noqa: E402
from scripts.analysis import build_loh_round1_cross_sample_audit as loh_audit  # noqa: E402
from scripts.analysis import build_loh_round1_cross_sample_audit_v2_figures as loh_figures  # noqa: E402
from scripts.analysis import build_observation_O01_master_distribution as o01  # noqa: E402
from scripts.analysis import build_phase1_training_manifest as phase1  # noqa: E402
from scripts.analysis import build_phase1a_split_manifest as phase1a  # noqa: E402
from scripts.analysis import build_snv_methylation_association_analysis as snv  # noqa: E402
from scripts.lib.verification_schema_contract import (  # noqa: E402
    UNKNOWN_CURRENT_CLASS,
    VERIFICATION_PROVENANCE_COLUMNS,
    SchemaContractError,
)


def v2_frame(include_unknown: bool = True) -> pd.DataFrame:
    rows = [
        {
            "Chr": "chr1",
            "Pos": 101,
            "Ref": "A",
            "Alt": "C",
            "VerificationSchemaVersion": 2,
            "VerificationClass": "Strong_Bidirectional",
            "VerificationClass_V1_Deprecated": "Strong",
            "VerificationClass_Legacy": "Strong",
            "LabelFirstSupport": "true",
            "ClusterFirstSupport": "true",
            "WithinHPSupport": "NA",
            "DispersionWarning": "NA",
            "EvidencePath": "BIDIRECTIONAL",
            "EvidenceDerivation": "LEGACY_CLASS",
            "LOH_Subtype_LegacyVC": "LOH_Strong",
            "LOH_Subtype": "LOH_Strong",
        }
    ]
    if include_unknown:
        rows.append(
            {
                "Chr": "chr1",
                "Pos": 102,
                "Ref": "G",
                "Alt": "T",
                "VerificationSchemaVersion": 2,
                "VerificationClass": "FutureClass",
                "VerificationClass_V1_Deprecated": "Noise_Uniform",
                "VerificationClass_Legacy": "Noise",
                "LabelFirstSupport": "false",
                "ClusterFirstSupport": "false",
                "WithinHPSupport": "NA",
                "DispersionWarning": "NA",
                "EvidencePath": "NOISE_UNIFORM",
                "EvidenceDerivation": "LEGACY_CLASS",
                "LOH_Subtype_LegacyVC": "None",
                "LOH_Subtype": "None",
            }
        )
    return pd.DataFrame(rows)


def h1_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Chr": ["chr1"],
            "Pos": [101],
            "Ref": ["A"],
            "Alt": ["C"],
            "VerificationClass": ["Strong"],
            "LOH_Subtype": ["LOH_Strong"],
        }
    )


class PassThroughConsumerTests(unittest.TestCase):
    def test_all_dataframe_consumers_fail_closed_then_mark_h1(self):
        adapters = [
            kism.validate_summary_provenance,
            to_feature.validate_summary_provenance,
            loh_audit.attach_verification_contract,
            o01.attach_master_provenance,
            phase1.attach_input_provenance,
        ]
        historical = h1_frame()
        for adapter in adapters:
            with self.subTest(adapter=adapter.__module__):
                with self.assertRaises(SchemaContractError):
                    adapter(historical)
                selected = adapter(historical, allow_unversioned_v1=True)
                self.assertEqual(
                    selected["VerificationProvenanceStatus"].unique().tolist(),
                    ["UNVERSIONED_V1"],
                )
                self.assertEqual(selected["VerificationSchemaVersion"].iloc[0], "")
                self.assertIn("--allow-unversioned-v1", selected["VerificationProvenanceWarnings"].iloc[0])

    def test_pass_through_adapters_keep_all_fields_and_reject_mixed_version(self):
        canonical = v2_frame()
        adapters = [
            kism.validate_summary_provenance,
            to_feature.validate_summary_provenance,
            phase1.attach_input_provenance,
        ]
        for adapter in adapters:
            with self.subTest(adapter=adapter.__module__):
                selected = adapter(canonical)
                self.assertTrue(set(VERIFICATION_PROVENANCE_COLUMNS).issubset(selected.columns))
                self.assertEqual(selected["VerificationProvenanceStatus"].unique().tolist(), ["V2"])
                warning = json.loads(selected["VerificationProvenanceWarnings"].iloc[0])
                self.assertEqual(warning["unknown_current_counts"], {"FutureClass": 1})

                mixed = canonical.copy()
                mixed.loc[1, "VerificationSchemaVersion"] = 1
                with self.assertRaisesRegex(SchemaContractError, "expected VerificationSchemaVersion=2"):
                    adapter(mixed)

    def test_phase1a_passes_exact_fields_and_rejects_mixed_status(self):
        canonical = phase1.attach_input_provenance(v2_frame(include_unknown=False))
        self.assertEqual(phase1a.validate_manifest_provenance(canonical), "V2")
        rows = phase1a.build_split_rows(
            canonical.assign(
                dataset_id="d1",
                dataset_label="d1",
                sample="s1",
                platform="ONT",
                mode="paired",
                region_key="chr1:101:A:C",
                truth_status="TP",
                source_scope="tp",
            )
        )
        for field in phase1.PROVENANCE_EXPORT_FIELDS:
            self.assertEqual(rows[0][field], canonical.iloc[0][field])

        mixed = canonical.copy()
        mixed.loc[0, "VerificationProvenanceStatus"] = "UNVERSIONED_V1"
        with self.assertRaisesRegex(SchemaContractError, "cannot carry a schema version"):
            phase1a.validate_manifest_provenance(mixed, allow_unversioned_v1=True)

        historical = phase1.attach_input_provenance(h1_frame(), allow_unversioned_v1=True)
        with self.assertRaises(SchemaContractError):
            phase1a.validate_manifest_provenance(historical)
        self.assertEqual(
            phase1a.validate_manifest_provenance(historical, allow_unversioned_v1=True),
            "UNVERSIONED_V1",
        )


class CurrentAndLohViewTests(unittest.TestCase):
    def test_loh_audit_buckets_unknown_and_uses_canonical_alias(self):
        selected = loh_audit.attach_verification_contract(v2_frame())
        self.assertEqual(selected["verification_class"].iloc[1], UNKNOWN_CURRENT_CLASS)
        self.assertEqual(selected["loh_subtype"].tolist(), ["LOH_Strong", "None"])
        self.assertEqual(selected["LOHProvenanceSourceField"].unique().tolist(), ["LOH_Subtype_LegacyVC"])
        summary = loh_audit.summarize_verification_by_loh(
            selected.assign(
                sample="s1",
                sample_label="s1",
                mode="paired",
                mode_label="paired-pure",
                truth_label="TP",
            ),
            Path("synthetic.tsv.gz"),
            "2026-07-14",
        )
        self.assertTrue(set(VERIFICATION_PROVENANCE_COLUMNS).issubset(summary.columns))

        mismatch = v2_frame()
        mismatch.loc[0, "LOH_Subtype"] = "None"
        with self.assertRaisesRegex(SchemaContractError, "not an exact alias"):
            loh_audit.attach_verification_contract(mismatch)

    def test_o01_retains_raw_value_while_plot_view_uses_unknown_bucket(self):
        selected = o01.attach_master_provenance(v2_frame())
        self.assertEqual(selected["VerificationClass_SourceRaw"].iloc[1], "FutureClass")
        self.assertEqual(selected["VerificationClass"].iloc[1], UNKNOWN_CURRENT_CLASS)
        self.assertEqual(selected["LOH_Subtype"].tolist(), ["LOH_Strong", "None"])
        marked_h1 = phase1.attach_input_provenance(h1_frame(), allow_unversioned_v1=True)
        selected_h1 = o01.attach_master_provenance(marked_h1, allow_unversioned_v1=True)
        self.assertEqual(selected_h1["VerificationProvenanceStatus"].iloc[0], "UNVERSIONED_V1")

    def test_loh_figure_consumer_validates_derived_metadata(self):
        all_df = loh_audit.attach_verification_contract(v2_frame())
        verif_df = pd.DataFrame(
            {
                "verification_class": all_df["verification_class"],
                "loh_subtype": all_df["loh_subtype"],
                "VerificationProvenanceStatus": ["V2", "V2"],
                "VerificationSchemaVersion": [2, 2],
                "VerificationProvenanceSourceField": ["VerificationClass", "VerificationClass"],
                "LOHProvenanceSourceField": [
                    "LOH_Subtype_LegacyVC",
                    "LOH_Subtype_LegacyVC",
                ],
            }
        )
        receipt = loh_figures.validate_round1_provenance(verif_df, all_df)
        self.assertEqual(receipt["schema_status"], "V2")
        self.assertIn(UNKNOWN_CURRENT_CLASS, receipt["class_order"])

        verif_df.loc[1, "VerificationSchemaVersion"] = 1
        with self.assertRaisesRegex(SchemaContractError, "VerificationSchemaVersion=2"):
            loh_figures.validate_round1_provenance(verif_df, all_df)

        historical_all = pd.DataFrame(
            {"verification_class": ["Strong"], "loh_subtype": ["LOH_Strong"]}
        )
        historical_summary = historical_all.copy()
        with self.assertRaises(SchemaContractError):
            loh_figures.validate_round1_provenance(historical_summary, historical_all)
        receipt = loh_figures.validate_round1_provenance(
            historical_summary,
            historical_all,
            allow_unversioned_v1=True,
        )
        self.assertEqual(receipt["schema_status"], "UNVERSIONED_DERIVED")


class SnvSidecarLoaderTests(unittest.TestCase):
    @staticmethod
    def science_frame() -> pd.DataFrame:
        frame = v2_frame(include_unknown=False)
        defaults = {
            "RegionID": "r1",
            "NumReads": 20,
            "NumCpGs": 5,
            "HPMergedDelta": 0.2,
            "HPMergedP": 0.01,
            "HPMergedSig": True,
            "AlleleDelta": 0.2,
            "AlleleP": 0.01,
            "AlleleSig": True,
            "HPFineSig": True,
            "HPFineP": 0.01,
            "GlobalP": 0.01,
            "ClusterPermanovaP": 0.01,
            "LabelAllelePermanovaP": 0.01,
            "Significant": True,
            "PassedGating": True,
            "sample": "HCC1395",
            "mode": "paired",
            "truth_label": "TP",
            "Potential_LOH": True,
        }
        return frame.assign(**defaults)

    def test_snv_loader_carries_complete_provenance_and_rejects_h1_by_default(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            canonical_path = Path(temp_dir) / "canonical.tsv.gz"
            historical_path = Path(temp_dir) / "historical.tsv.gz"
            canonical = self.science_frame()
            canonical.to_csv(canonical_path, sep="\t", index=False, compression="gzip")
            loaded = snv.load_master_with_provenance(canonical_path)
            self.assertTrue(set(VERIFICATION_PROVENANCE_COLUMNS).issubset(loaded.columns))
            self.assertEqual(loaded["VerificationProvenanceStatus"].iloc[0], "V2")

            historical = canonical.drop(
                columns=[
                    field
                    for field in VERIFICATION_PROVENANCE_COLUMNS
                    if field not in {"VerificationClass", "LOH_Subtype"}
                ]
            )
            historical["VerificationClass"] = "Strong"
            historical.to_csv(historical_path, sep="\t", index=False, compression="gzip")
            with self.assertRaises(SchemaContractError):
                snv.load_master_with_provenance(historical_path)
            loaded_h1 = snv.load_master_with_provenance(
                historical_path,
                allow_unversioned_v1=True,
            )
            self.assertEqual(loaded_h1["VerificationProvenanceStatus"].iloc[0], "UNVERSIONED_V1")

            marked_path = Path(temp_dir) / "marked_h1.tsv.gz"
            marked_h1 = loh_audit.attach_verification_contract(
                historical,
                allow_unversioned_v1=True,
            )
            marked_h1.to_csv(marked_path, sep="\t", index=False, compression="gzip")
            loaded_marked_h1 = snv.load_master_with_provenance(
                marked_path,
                allow_unversioned_v1=True,
            )
            self.assertEqual(
                loaded_marked_h1["VerificationProvenanceStatus"].iloc[0],
                "UNVERSIONED_V1",
            )


if __name__ == "__main__":
    unittest.main()
