#!/usr/bin/env python3
"""Focused contract tests for the B1-Guards-2 downstream consumer migration."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"
for path in (REPO_ROOT, ANALYSIS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from scripts.analysis import compare_phase1a_model_errors as compare_phase1a  # noqa: E402
from scripts.analysis import extract_label_first_metrics as label_first  # noqa: E402
from scripts.analysis import pon_cross_sample_and_h2009_diagnosis as pon  # noqa: E402
from scripts.analysis import seqc2_cnv_cross_sample_and_rootcause as seqc2  # noqa: E402
from scripts.analysis import summarize_phase1a_prediction_failures as summarize_phase1a  # noqa: E402
from scripts.analysis import tp_fp_structure_label_figures as figures  # noqa: E402
from scripts.lib.verification_schema_contract import (  # noqa: E402
    LEGACY_CLASSES,
    UNKNOWN_CURRENT_CLASS,
    SchemaContractError,
)
from scripts.validation.metrics import compute_metrics  # noqa: E402


def v2_frame(classes=("Strong_Bidirectional",)) -> pd.DataFrame:
    rows = []
    for current in classes:
        if current == "Strong_Bidirectional":
            deprecated, legacy, path = "Strong", "Strong", "BIDIRECTIONAL"
            label_first_support, cluster_first_support = "true", "true"
        elif current == "ClusterFirstOnly":
            deprecated, legacy, path = "Strong", "Subclone", "CLUSTER_FIRST_ONLY"
            label_first_support, cluster_first_support = "false", "true"
        else:
            deprecated, legacy, path = current, "Weak", "NOISE_UNIFORM"
            label_first_support, cluster_first_support = "true", "false"
        rows.append(
            {
                "Chr": "chr1",
                "Pos": "100",
                "Ref": "A",
                "Alt": "C",
                "VerificationSchemaVersion": "2",
                "VerificationClass": current,
                "VerificationClass_V1_Deprecated": deprecated,
                "VerificationClass_Legacy": legacy,
                "LabelFirstSupport": label_first_support,
                "ClusterFirstSupport": cluster_first_support,
                "WithinHPSupport": "NA",
                "DispersionWarning": "NA",
                "EvidencePath": path,
                "EvidenceDerivation": "LEGACY_CLASS",
            }
        )
    return pd.DataFrame(rows)


class CurrentConsumerGuardTests(unittest.TestCase):
    def test_phase1_consumers_require_v2_provenance_and_bucket_unknown(self):
        frame = v2_frame(("Strong_Bidirectional", "FutureClass"))
        normalized, metadata = compare_phase1a.attach_current_taxonomy(frame)
        self.assertEqual(metadata["schema_status"], "V2")
        self.assertEqual(normalized.loc[1, "VerificationClass"], UNKNOWN_CURRENT_CLASS)
        self.assertEqual(metadata["unknown_counts"], {"FutureClass": 1})

        missing = frame.drop(columns=["EvidencePath"])
        with self.assertRaises(SchemaContractError):
            summarize_phase1a.attach_current_taxonomy(missing)

    def test_phase1_h1_requires_flag_and_warns_without_inventing_evidence(self):
        frame = pd.DataFrame({"VerificationClass": ["Strong", "FutureClass"]})
        with self.assertRaises(SchemaContractError):
            compare_phase1a.attach_current_taxonomy(frame)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            normalized, metadata = compare_phase1a.attach_current_taxonomy(
                frame, allow_unversioned_v1=True
            )
        self.assertTrue(caught)
        self.assertEqual(metadata["schema_status"], "UNVERSIONED_V1")
        self.assertEqual(normalized.loc[1, "VerificationClass"], UNKNOWN_CURRENT_CLASS)
        self.assertEqual(normalized.loc[0, "EvidencePath"], "")

    def test_label_first_uses_direct_evidence_and_preserves_serialization(self):
        rows = v2_frame().to_dict(orient="records")
        normalized, metadata = label_first.validate_label_first_schema(rows)
        self.assertEqual(metadata["schema_status"], "V2")
        self.assertEqual(normalized[0]["LabelFirstSupport"], "true")
        self.assertEqual(normalized[0]["EvidencePath"], "BIDIRECTIONAL")

        bad = v2_frame()
        bad.loc[0, "LabelFirstSupport"] = "Strong"
        with self.assertRaises(SchemaContractError):
            label_first.validate_label_first_schema(bad.to_dict(orient="records"))

    def test_pon_and_seqc2_reject_mixed_versions(self):
        frame = v2_frame(("Strong_Bidirectional", "ClusterFirstOnly"))
        frame.loc[1, "VerificationSchemaVersion"] = "1"
        rows = frame.to_dict(orient="records")
        with self.assertRaises(SchemaContractError):
            pon._validate_current_taxonomy(rows, list(frame.columns), False)
        with self.assertRaises(SchemaContractError):
            seqc2.normalize_current_taxonomy(frame)

    def test_label_first_cli_validates_before_writing_and_passes_evidence(self):
        script = REPO_ROOT / "scripts" / "analysis" / "extract_label_first_metrics.py"
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_path = root / "significance_summary.csv"
            output_path = root / "label_first.tsv"
            frame = v2_frame()
            frame["RegionID"] = ["chr1:100"]
            frame["Chr"] = ["chr1"]
            frame["Pos"] = ["100"]
            frame["Ref"] = ["A"]
            frame["Alt"] = ["C"]
            frame.to_csv(input_path, index=False)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--summary-csv",
                    str(input_path),
                    "--output-tsv",
                    str(output_path),
                ],
                cwd=REPO_ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            derived = pd.read_csv(output_path, sep="\t", dtype=str, keep_default_na=False)
            self.assertEqual(derived.loc[0, "LabelFirstSupport"], "true")
            self.assertEqual(derived.loc[0, "EvidencePath"], "BIDIRECTIONAL")

            invalid_output = root / "must_not_exist.tsv"
            frame.drop(columns=["EvidencePath"]).to_csv(input_path, index=False)
            rejected = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--summary-csv",
                    str(input_path),
                    "--output-tsv",
                    str(invalid_output),
                ],
                cwd=REPO_ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(rejected.returncode, 0)
            self.assertFalse(invalid_output.exists())


class MetricsAndFigureGuardTests(unittest.TestCase):
    def test_metrics_records_v2_evidence_and_unknown_bucket(self):
        frame = v2_frame(("Strong_Bidirectional", "FutureClass"))
        frame["Quality_Score"] = [0.8, 0.2]
        frame["HPFineNGroups"] = [1, 2]
        frame["DominantLabel"] = ["HP1", "HP2"]
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "significance_summary.csv"
            frame.to_csv(path, index=False)
            stats = compute_metrics.extract_significance_stats(path)
        self.assertEqual(stats["verification_taxonomy"]["schema_status"], "V2")
        self.assertEqual(stats["verification_class_dist"][UNKNOWN_CURRENT_CLASS], 1)
        self.assertEqual(stats["verification_evidence_dist"]["EvidencePath"]["BIDIRECTIONAL"], 1)
        self.assertEqual(stats["verification_evidence_dist"]["LabelFirstSupport"]["true"], 2)

    def test_metrics_h1_opt_in_and_cross_file_mixed_status_fail(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            v2_dir = root / "v2"
            h1_dir = root / "h1"
            v2_dir.mkdir()
            h1_dir.mkdir()
            v2_frame().to_csv(v2_dir / "significance_summary.csv", index=False)
            pd.DataFrame({"VerificationClass": ["Strong"]}).to_csv(
                h1_dir / "significance_summary.csv", index=False
            )
            with self.assertRaises(SchemaContractError):
                compute_metrics.build_metrics_bundle(root, allow_unversioned_v1=True)

    def test_metrics_cli_emits_bundle_only_after_v2_contract_passes(self):
        script = REPO_ROOT / "scripts" / "validation" / "metrics" / "compute_metrics.py"
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            run_dir = root / "run" / "intersubmod_tp"
            run_dir.mkdir(parents=True)
            v2_frame().to_csv(run_dir / "significance_summary.csv", index=False)
            output = root / "metrics_bundle.json"
            completed = subprocess.run(
                [str(script), "--experiment-dir", str(root / "run"), "--output", str(output)],
                cwd=REPO_ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            bundle = json.loads(output.read_text())
            self.assertEqual(bundle["verification_taxonomy"]["schema_status"], "V2")

    @staticmethod
    def _write_figure_assets(root: Path, status="LEGACY_EXPLICIT", old_key=False):
        key = "verdict_VerificationClass" if old_key else "legacy_VerificationClass"
        selection_field = "VerificationClass_Legacy" if status == "LEGACY_EXPLICIT" else "VerificationClass"
        warning_values = [] if status == "LEGACY_EXPLICIT" else ["UNVERSIONED upstream"]
        meta = {
            "verification_selection": {
                "selection_field": selection_field,
                "schema_status": status,
                "categories": list(LEGACY_CLASSES),
                "unknown_counts": {},
                "warnings": warning_values,
            },
            "verification_selection_values": list(LEGACY_CLASSES),
        }
        crosstabs = {key: {"cells": []}}
        confound = {"matched_subset": {"crosstabs": {key: {"cells": []}}}}
        (root / "run_meta.json").write_text(json.dumps(meta))
        (root / "crosstabs.json").write_text(json.dumps(crosstabs))
        (root / "confound_control.json").write_text(json.dumps(confound))

    def test_figures_require_explicit_upstream_taxonomy_and_new_dimension_key(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self._write_figure_assets(root)
            _, _, _, selection = figures.load_verified_assets(root)
            self.assertEqual(selection["selection_field"], "VerificationClass_Legacy")

            self._write_figure_assets(root, old_key=True)
            with self.assertRaises(SchemaContractError):
                figures.load_verified_assets(root)

    def test_figure_h1_requires_local_explicit_flag(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self._write_figure_assets(root, status="UNVERSIONED_V1")
            with self.assertRaises(SchemaContractError):
                figures.load_verified_assets(root)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                figures.load_verified_assets(root, allow_unversioned_v1=True)
            self.assertTrue(caught)


class FrozenGoldenGuardTests(unittest.TestCase):
    def test_regression_script_refuses_automatic_golden_rebuild(self):
        script = REPO_ROOT / "scripts" / "regression" / "regression_check.sh"
        completed = subprocess.run(
            [str(script), "--update-golden"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("golden artifacts are frozen", completed.stderr)
        self.assertNotIn("cp ", script.read_text())

    def test_stale_goldens_fail_before_cpp_execution(self):
        script = REPO_ROOT / "scripts" / "regression" / "regression_check.sh"
        completed = subprocess.run(
            [str(script)],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("stale/non-v2 header", completed.stderr)


if __name__ == "__main__":
    unittest.main()
