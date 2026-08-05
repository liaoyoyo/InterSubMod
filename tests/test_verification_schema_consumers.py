#!/usr/bin/env python3
"""Data-independent regression tests for verification schema consumers."""

from __future__ import annotations

import io
import tempfile
import unittest
import warnings
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from scripts.analysis.analyze_methyl_cluster_allele_cooccurrence import (
    select_replication_cohorts,
)
from scripts.analysis.compare_subclone_validation import derive_legacy_support_columns
from scripts.lib.verification_schema_contract import (
    CURRENT_CLASSES_V2,
    LEGACY_CLASSES,
    UNKNOWN_CURRENT_CLASS,
    VERIFICATION_PROVENANCE_COLUMNS,
    SchemaContractError,
    extract_provenance_frame,
    ordered_class_crosstab,
    read_evidence,
    select_current_view,
    select_legacy_view,
    select_loh_legacy,
    validate_region_strata,
)
from tools.compare_vcf_results import plot_verification_class, resolve_verification_view
from tools.find_verification_candidates import analyze_candidates, prepare_legacy_candidates


class CurrentAndLegacyViewTest(unittest.TestCase):
    def test_v2_all_classes_unknown_and_stable_order(self):
        raw = list(CURRENT_CLASSES_V2) + ["UnexpectedFutureClass"]
        df = pd.DataFrame(
            {
                "VerificationSchemaVersion": [2] * len(raw),
                "VerificationClass": raw,
            }
        )

        view = select_current_view(df)

        self.assertEqual(view.schema_status, "V2")
        self.assertEqual(view.categories, CURRENT_CLASSES_V2 + (UNKNOWN_CURRENT_CLASS,))
        self.assertEqual(view.values.iloc[-1], UNKNOWN_CURRENT_CLASS)
        self.assertEqual(view.unknown_counts, {"UnexpectedFutureClass": 1})

    def test_zero_frequency_classes_remain_in_compare_table(self):
        df = pd.DataFrame(
            {
                "VerificationSchemaVersion": [2, 2],
                "VerificationClass": ["Strong_Bidirectional", "FutureClass"],
                "Label": ["TP", "FP"],
            }
        )
        view = resolve_verification_view(df, requested_view="current")

        table = ordered_class_crosstab(df["Label"], view)

        self.assertEqual(list(table.columns), list(CURRENT_CLASSES_V2) + [UNKNOWN_CURRENT_CLASS])
        self.assertEqual(int(table["ClusterFirstOnly"].sum()), 0)
        self.assertEqual(int(table[UNKNOWN_CURRENT_CLASS].sum()), 1)

        with tempfile.TemporaryDirectory() as tmpdir:
            plotted = plot_verification_class(df, {"distribution": tmpdir}, view)
            self.assertTrue((Path(tmpdir) / "dist_verification_class.png").is_file())
        self.assertEqual(list(plotted.columns), list(CURRENT_CLASSES_V2) + [UNKNOWN_CURRENT_CLASS])

    def test_unversioned_current_is_raw_and_marked(self):
        df = pd.DataFrame({"VerificationClass": ["Strong", "UnmappedRaw"]})

        view = resolve_verification_view(df, requested_view="current")

        self.assertEqual(view.schema_status, "UNVERSIONED")
        self.assertEqual(view.field, "VerificationClass")
        self.assertEqual(view.values.tolist(), ["Strong", "UnmappedRaw"])
        self.assertIn("UNVERSIONED", view.warning_messages[0])

    def test_explicit_legacy_requires_four_valid_values(self):
        df = pd.DataFrame({"VerificationClass_Legacy": list(LEGACY_CLASSES)})
        view = select_legacy_view(df)
        self.assertEqual(view.values.tolist(), list(LEGACY_CLASSES))
        self.assertEqual(view.field, "VerificationClass_Legacy")

        with self.assertRaisesRegex(SchemaContractError, "invalid values"):
            select_legacy_view(pd.DataFrame({"VerificationClass_Legacy": ["Strong", "Other"]}))
        with self.assertRaisesRegex(SchemaContractError, "<MISSING>"):
            select_legacy_view(pd.DataFrame({"VerificationClass_Legacy": ["Strong", None]}))
        with self.assertRaisesRegex(SchemaContractError, "invalid values"):
            select_legacy_view(pd.DataFrame({"VerificationClass_Legacy": [" Strong"]}))
        with self.assertRaisesRegex(SchemaContractError, "is missing"):
            select_legacy_view(pd.DataFrame({"VerificationClass": ["Strong"]}))

    def test_h1_fallback_is_explicit_warns_and_counts_exclusions(self):
        df = pd.DataFrame({"VerificationClass": ["Strong", "Noise", "Mystery"]})
        with self.assertRaises(SchemaContractError):
            select_legacy_view(df)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            view = select_legacy_view(
                df,
                allow_unversioned_v1=True,
                unversioned_unknown_policy="exclude",
            )

        self.assertEqual(view.schema_status, "UNVERSIONED_V1")
        self.assertEqual(view.unknown_counts, {"Mystery": 1})
        self.assertTrue(pd.isna(view.values.iloc[2]))
        self.assertTrue(any("UNVERSIONED" in str(item.message) for item in caught))

        versioned = df.assign(VerificationSchemaVersion=2)
        with self.assertRaisesRegex(SchemaContractError, "fallback is forbidden"):
            select_legacy_view(versioned, allow_unversioned_v1=True)
        malformed_versioned = df.assign(VerificationSchemaVersion=pd.NA)
        with self.assertRaisesRegex(SchemaContractError, "fallback is forbidden"):
            select_legacy_view(malformed_versioned, allow_unversioned_v1=True)


class EvidenceLohRegionProvenanceTest(unittest.TestCase):
    def test_typed_evidence_booleans_and_na(self):
        df = pd.DataFrame(
            {
                "LabelFirstSupport": ["true", "false"],
                "ClusterFirstSupport": [True, False],
                "WithinHPSupport": ["NA", "true"],
                "DispersionWarning": [pd.NA, "false"],
                "EvidencePath": ["BIDIRECTIONAL", "NOISE_UNCORRELATED"],
                "EvidenceDerivation": ["LIVE", "LEGACY_CLASS"],
            }
        )

        evidence = read_evidence(df)

        self.assertEqual(str(evidence["LabelFirstSupport"].dtype), "boolean")
        self.assertEqual(evidence["LabelFirstSupport"].tolist(), [True, False])
        self.assertTrue(pd.isna(evidence["WithinHPSupport"].iloc[0]))
        self.assertTrue(pd.isna(evidence["DispersionWarning"].iloc[0]))

        invalid = df.copy()
        invalid.loc[0, "LabelFirstSupport"] = "TRUE"
        with self.assertRaisesRegex(SchemaContractError, "invalid boolean"):
            read_evidence(invalid)

    def test_loh_canonical_alias_and_h1(self):
        explicit = pd.DataFrame(
            {
                "LOH_Subtype_LegacyVC": ["None", "LOH_Subclone"],
                "LOH_Subtype": ["None", "LOH_Subclone"],
            }
        )
        view = select_loh_legacy(explicit)
        self.assertEqual(view.field, "LOH_Subtype_LegacyVC")

        mismatch = explicit.copy()
        mismatch.loc[1, "LOH_Subtype"] = "LOH_Strong"
        with self.assertRaisesRegex(SchemaContractError, "not an exact alias"):
            select_loh_legacy(mismatch)
        with self.assertRaisesRegex(SchemaContractError, "missing required columns"):
            select_loh_legacy(explicit.drop(columns="LOH_Subtype"))

        old = pd.DataFrame({"LOH_Subtype": ["None", "LOH_Weak"]})
        with self.assertRaises(SchemaContractError):
            select_loh_legacy(old)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fallback = select_loh_legacy(old, allow_unversioned_v1=True)
        self.assertEqual(fallback.schema_status, "UNVERSIONED_V1")
        self.assertTrue(any("UNVERSIONED" in str(item.message) for item in caught))

    @staticmethod
    def _status(status, assignment_count, occupied):
        return {
            "status": status,
            "reason": "OK" if status == "VALID" else "NO_ELIGIBLE_REGIONS",
            "schema_version": 1,
            "eligible_region_count": 1 if status == "VALID" else 0,
            "min_regions_required": 1,
            "assignment_count": assignment_count,
            "n_occupied_region_strata": occupied,
            "warning_count": 0,
            "generated_at": "2026-07-14T22:00:00+08:00",
        }

    def test_region_status_assignment_and_sentinels(self):
        valid = pd.DataFrame(
            {
                "RegionStratificationSchemaVersion": [1, 1],
                "RegionStratum_ID": [0, -1],
                "RegionStratum_Label": ["BaselineLowAsm", "Unassigned"],
                "RegionStratum_Reason": ["BASELINE_LOW_ASM", "INELIGIBLE_REGION"],
                "Subclone_ID": [0, -1],
            }
        )
        checked = validate_region_strata(valid, self._status("VALID", 1, 1))
        self.assertEqual(checked.assignment_count, 1)
        self.assertEqual(checked.n_occupied_region_strata, 1)

        insufficient = valid.assign(
            RegionStratum_ID=-1,
            RegionStratum_Label="Unassigned",
            RegionStratum_Reason="INSUFFICIENT_REGIONS",
            Subclone_ID=-1,
        )
        checked = validate_region_strata(
            insufficient,
            self._status("INSUFFICIENT_REGIONS", 0, 0),
        )
        self.assertEqual(checked.assignment_count, 0)

        invalid = insufficient.copy()
        invalid.loc[0, "RegionStratum_ID"] = 2
        invalid.loc[0, "Subclone_ID"] = 2
        invalid.loc[0, "RegionStratum_Label"] = "LohFlagged"
        invalid.loc[0, "RegionStratum_Reason"] = "LOH_FLAGGED"
        with self.assertRaisesRegex(SchemaContractError, "assigned stratum"):
            validate_region_strata(invalid, self._status("INSUFFICIENT_REGIONS", 1, 1))

        bad_status = self._status("VALID", 1, 1)
        bad_status["warning_count"] = 0.5
        with self.assertRaisesRegex(SchemaContractError, "non-negative integer"):
            validate_region_strata(valid, bad_status)

    def test_complete_v2_provenance_projection(self):
        df = pd.DataFrame(
            {
                "VerificationSchemaVersion": [2],
                "VerificationClass": ["Strong_Bidirectional"],
                "VerificationClass_V1_Deprecated": ["Strong"],
                "VerificationClass_Legacy": ["Strong"],
                "LabelFirstSupport": ["true"],
                "ClusterFirstSupport": ["true"],
                "WithinHPSupport": ["false"],
                "DispersionWarning": ["false"],
                "EvidencePath": ["BIDIRECTIONAL"],
                "EvidenceDerivation": ["LIVE"],
                "LOH_Subtype_LegacyVC": ["None"],
                "LOH_Subtype": ["None"],
            }
        )

        provenance = extract_provenance_frame(df)

        self.assertEqual(tuple(provenance.columns), VERIFICATION_PROVENANCE_COLUMNS)
        self.assertEqual(provenance.attrs["schema_status"], "V2")

        with self.assertRaisesRegex(SchemaContractError, "missing required columns"):
            extract_provenance_frame(df.drop(columns="EvidencePath"))

        inconsistent = df.copy()
        inconsistent.loc[0, "EvidencePath"] = "CLUSTER_FIRST_ONLY"
        with self.assertRaisesRegex(SchemaContractError, "EvidencePath"):
            extract_provenance_frame(inconsistent)

        offline = df.copy()
        offline.loc[0, "WithinHPSupport"] = "NA"
        offline.loc[0, "DispersionWarning"] = "NA"
        offline.loc[0, "EvidenceDerivation"] = "LEGACY_CLASS"
        self.assertEqual(extract_provenance_frame(offline).attrs["schema_status"], "V2")

        invalid_offline = offline.copy()
        invalid_offline.loc[0, "WithinHPSupport"] = "false"
        with self.assertRaisesRegex(SchemaContractError, "LEGACY_CLASS requires WithinHPSupport=NA"):
            extract_provenance_frame(invalid_offline)

        invalid_dispersion = offline.copy()
        invalid_dispersion.loc[0, "DispersionWarning"] = "false"
        with self.assertRaisesRegex(SchemaContractError, "LEGACY_CLASS requires DispersionWarning=NA"):
            extract_provenance_frame(invalid_dispersion)


class ConsumerContractTest(unittest.TestCase):
    def test_candidate_report_uses_legacy_field_and_hits_each_class(self):
        df = pd.DataFrame(
            {
                "VerificationClass": ["Noise_Uncorrelated"] * 4,
                "VerificationClass_Legacy": list(LEGACY_CLASSES),
                "Chr": ["chr1"] * 4,
                "Pos": [1, 2, 3, 4],
                "GlobalP": [0.1, 0.2, 0.3, 0.4],
                "LabelP": [0.1, 0.2, 0.3, 0.4],
                "LabelDelta": [0.1, 0.2, 0.3, 0.4],
            }
        )
        selected, view = prepare_legacy_candidates(df)
        self.assertEqual(view.field, "VerificationClass_Legacy")

        with tempfile.TemporaryDirectory() as tmpdir, redirect_stdout(io.StringIO()) as captured:
            analyze_candidates(selected, SimpleNamespace(output_dir=tmpdir, top_n=1))

        report = captured.getvalue()
        self.assertIn("# Legacy verification candidates", report)
        for legacy_class in LEGACY_CLASSES:
            self.assertIn(f"## Category: {legacy_class} (Total: 1)", report)

    def test_replication_selector_explicit_and_h1_unknown_policy(self):
        explicit = pd.DataFrame(
            {"VerificationClass_Legacy": ["Strong", "Noise", "Weak", "Subclone"]}
        )
        selected, metadata = select_replication_cohorts(explicit)
        self.assertEqual(selected["_VerificationClass_Legacy_Selected"].tolist(), ["Strong", "Noise"])
        self.assertEqual(metadata["selected_counts"], {"Strong": 1, "Noise": 1})
        self.assertEqual(metadata["selection_field"], "VerificationClass_Legacy")

        historical = pd.DataFrame({"VerificationClass": ["Strong", "Noise", "UnknownOld"]})
        with self.assertRaises(SchemaContractError):
            select_replication_cohorts(historical)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            selected, metadata = select_replication_cohorts(
                historical,
                allow_unversioned_v1=True,
            )
        self.assertEqual(len(selected), 2)
        self.assertEqual(metadata["excluded_unknown_count"], 1)
        self.assertEqual(metadata["unknown_counts"], {"UnknownOld": 1})
        self.assertTrue(any("UNVERSIONED" in str(item.message) for item in caught))

    def test_legacy_support_booleans_are_locked(self):
        df = pd.DataFrame({"VerificationClass_Legacy": list(LEGACY_CLASSES)})

        derived, view = derive_legacy_support_columns(df)

        self.assertEqual(view.field, "VerificationClass_Legacy")
        self.assertEqual(
            derived["LegacyVerificationSupport"].tolist(),
            [True, True, False, False],
        )
        self.assertEqual(
            derived["LegacyClusterFirstOnly"].tolist(),
            [False, True, False, False],
        )


if __name__ == "__main__":
    unittest.main()
