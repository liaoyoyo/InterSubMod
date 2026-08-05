#!/usr/bin/env python3
"""Regression tests for the seven-dataset backbone sensitivity contract."""

import importlib.util
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
SUMMARIZER = SCRIPT_DIR / "summarize_backbone_sensitivity.py"
REPORT = SCRIPT_DIR / "build_layered_v2_validation_report.py"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]


def load_module(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def view():
    return {
        "census": {
            "U1_sSNV_somatic_total": 10,
            "U1_sSNV_scope_chr1_22": 10,
            "n_regions": 2,
            "n_regions_with_primary_lineage": 2,
            "L1": {
                "n_primary_lineage_units": 2,
                "n_reference_only_controls": 1,
                "n_unresolved_H3_auxiliary": 0,
                "n_shared_H4_auxiliary": 0,
                "determinacy_primary_lineage": {"determined": 1},
            },
            "hp_multiplicity": {"2": 1},
            "region_determinacy": {"all_determined": 1},
            "funnel": {"L6_retained_sSNV": 2},
        }
    }


class BackboneSensitivityContractTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        producer = self.root / "producer"
        producer.mkdir()
        closeout = producer / "production_closeout.json"
        success = producer / "_SUCCESS"
        artifacts = producer / "artifacts.final.sha256"
        closeout.write_text("{}\n", encoding="utf-8")
        success.write_text("{}\n", encoding="utf-8")
        artifacts.write_text("synthetic\n", encoding="utf-8")
        digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
        self.manifest_path = self.root / "canonical.json"
        self.manifest = {
            "schema_version": "2.1",
            "task_type": "B_COMPREHENSIVE_VALIDATION",
            "dataset_count": 7,
            "tree_input_contract": "longphase_recalibrated_PASS",
            "tag_contract": {
                "truth_flags": False,
                "PS_preserved": True,
                "bam_identity_locked": True,
                "tree_backbone": "LongPhase-S _sc.vcf FILTER=PASS",
                "longphase_filtering_policy": "production_default_filter",
                "production_closeout": str(closeout),
                "production_closeout_sha256": digest(closeout),
                "production_success": str(success),
                "production_success_sha256": digest(success),
                "production_artifacts_manifest": str(artifacts),
                "production_artifacts_manifest_sha256": digest(artifacts),
            },
            "samples": [{"sample": sample} for sample in SAMPLES],
        }
        write_json(self.manifest_path, self.manifest)
        self.base = self.root / "base"
        self.clairs = self.root / "clairs"
        self._make_run(self.base, "longphase_recalibrated_PASS")
        self._make_run(self.clairs, "clairs_PASS_input")

    def tearDown(self):
        self.temp.cleanup()

    def _make_run(self, root, contract):
        root.mkdir()
        (root / "_SUCCESS").write_text("PASS\n", encoding="utf-8")
        write_json(root / "input_manifest.json", {
            "tree_input_contract": contract,
            "samples": [{"sample": sample} for sample in SAMPLES],
        })
        write_json(root / "verification_summary.json", {
            "all_pass": True,
            "n_pass": 7,
            "dataset_count": 7,
            "samples": [{"sample": sample} for sample in SAMPLES],
        })
        for index, sample in enumerate(SAMPLES, start=1):
            sample_dir = root / "samples" / sample
            write_json(sample_dir / f"layered_region_view_{sample}.json", view())
            write_json(sample_dir / f"layered_reconstruction_{sample}.json", {
                "analysis_contract": {
                    "PS": "preserved for phase-block QC; not used as a topology edge or lineage label"
                },
                "read_tag_census": {
                    "check_exact_sidecar_join": True,
                    "phase_set_region_counts": {"one": 1},
                    "n_regions_with_phase_set_census": 1,
                },
                "detail": [{
                    "chrom": "chr1", "start": 100 + index, "end": 200 + index,
                    "family": "1", "is_primary_lineage": True, "capped": False,
                    "analysis_tree_digest_sha256": "a" * 64,
                }]
            })
            write_json(sample_dir / f"ssnv_site_ledger_{sample}.summary.json", {
                "pass": True,
                "longphase_input_contract": "clairs_raw_all",
                "tree_contract": contract,
            })
            write_json(sample_dir / "mlhp_part_chr1.json", {
                "groups": [{"chrom": "chr1", "positions": [100 + index, 200 + index]}]
            })

    def _run(self, output):
        return subprocess.run([
            sys.executable, str(SUMMARIZER),
            "--base-run", str(self.base),
            "--clairs-run", str(self.clairs),
            "--input-manifest", str(self.manifest_path),
            "--output", str(output),
        ], text=True, capture_output=True, check=False)

    def test_full_seven_dataset_summary_passes_report_gate(self):
        output = self.root / "summary.json"
        proc = self._run(output)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        summary = json.loads(output.read_text(encoding="utf-8"))
        self.assertEqual(len(summary["comparisons"]), 7)
        self.assertEqual(summary["aggregate"]["verdict"], "robust_all_gates")
        self.assertEqual(summary["aggregate"]["min_primary_unit_key_jaccard"], 1.0)
        self.assertEqual(summary["aggregate"]["min_shared_topology_digest_concordance"], 1.0)

        verification = json.loads((self.base / "verification_summary.json").read_text(encoding="utf-8"))
        for item in verification["samples"]:
            item["metrics"] = {
                "read_tag_census": {"check_exact_sidecar_join": True},
                "site_ledger": {"pass": True},
            }
        report = load_module(REPORT)
        report.validate_authoritative_inputs(
            self.manifest,
            verification,
            {"all_forbidden_counts_zero": True, "samples": [{"sample": sample} for sample in SAMPLES]},
            {"all_candidate_sets_complete": True, "samples": [{"sample": sample} for sample in SAMPLES]},
            {"aggregate": {key: {"min_primary_unit_key_jaccard": 1.0,
                                  "min_shared_topology_digest_concordance": 1.0}
                           for key in ("mapq30", "baseq10", "minread4", "maxsnv6")},
             "rows": [{"variant": key, "sample": sample}
                      for key in ("mapq30", "baseq10", "minread4", "maxsnv6") for sample in SAMPLES]},
            summary,
        )

    def test_missing_verified_dataset_is_rejected(self):
        verification_path = self.clairs / "verification_summary.json"
        verification = json.loads(verification_path.read_text(encoding="utf-8"))
        verification["samples"].pop()
        write_json(verification_path, verification)
        output = self.root / "must_not_exist.json"
        proc = self._run(output)
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse(output.exists())
        self.assertIn("verification/sample set failed", proc.stderr)

    def test_report_rejects_noncanonical_backbone_summary(self):
        report = load_module(REPORT)
        verification = json.loads((self.base / "verification_summary.json").read_text(encoding="utf-8"))
        for item in verification["samples"]:
            item["metrics"] = {
                "read_tag_census": {"check_exact_sidecar_join": True},
                "site_ledger": {"pass": True},
            }
        malformed = {
            "schema_version": "2.1",
            "canonical": "ClairS PASS",
            "alternative": "LongPhase-S PASS",
            "comparisons": [],
        }
        with self.assertRaises(SystemExit):
            report.validate_authoritative_inputs(
                self.manifest,
                verification,
                {"all_forbidden_counts_zero": True, "samples": [{"sample": sample} for sample in SAMPLES]},
                {"all_candidate_sets_complete": True, "samples": [{"sample": sample} for sample in SAMPLES]},
                {"aggregate": {key: {"min_primary_unit_key_jaccard": 1.0,
                                      "min_shared_topology_digest_concordance": 1.0}
                               for key in ("mapq30", "baseq10", "minread4", "maxsnv6")},
                 "rows": [{"variant": key, "sample": sample}
                          for key in ("mapq30", "baseq10", "minread4", "maxsnv6") for sample in SAMPLES]},
                malformed,
            )

    def test_topology_digest_change_cannot_be_labeled_robust(self):
        path = self.clairs / "samples" / SAMPLES[0] / f"layered_reconstruction_{SAMPLES[0]}.json"
        document = json.loads(path.read_text(encoding="utf-8"))
        document["detail"][0]["analysis_tree_digest_sha256"] = "b" * 64
        write_json(path, document)
        output = self.root / "topology_changed.json"
        proc = self._run(output)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        summary = json.loads(output.read_text(encoding="utf-8"))
        self.assertEqual(summary["comparisons"][0]["shared_unit_topology_digest_concordance"], 0.0)
        self.assertEqual(summary["aggregate"]["verdict"], "backbone_sensitive")

    def test_canonical_v3_lock_and_sorted_verification_are_accepted(self):
        manifest_path = self.base / "input_manifest.json"
        manifest_path.rename(self.base / "input_manifest.v2-not-used.json")
        write_json(self.base / "frozen_input_lock.json", {
            "all_pass": True,
            "analysis_contract": {
                "task_type": "comprehensive_validation",
                "tree_input_contract": "longphase_s_recalibrated_FILTER_PASS",
                "production_filter_policy": {"name": "production_default_filter"},
            },
            "samples": [{
                "sample": sample,
                "somatic": {
                    "caller_raw_vcf": {}, "longphase_input_vcf": {},
                    "caller_pass_baseline_vcf": {"path": f"/{sample}/caller-pass.vcf.gz"},
                    "longphase_recalibrated_all_vcf": {"path": f"/{sample}/all.vcf.gz"},
                    "longphase_recalibrated_pass_vcf": {"path": f"/{sample}/pass.vcf.gz"},
                    "tree_vcf": {"path": f"/{sample}/pass.vcf.gz"},
                },
            } for sample in reversed(SAMPLES)],
        })
        verification_path = self.base / "verification_summary.json"
        verification = json.loads(verification_path.read_text(encoding="utf-8"))
        verification["samples"] = list(reversed(verification["samples"]))
        write_json(verification_path, verification)
        output = self.root / "v3_base.json"
        proc = self._run(output)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(json.loads(output.read_text())["aggregate"]["verdict"], "robust_all_gates")

    def test_sensitivity_v3_lock_is_accepted_only_with_clairs_selected_tree(self):
        manifest_path = self.clairs / "input_manifest.json"
        manifest_path.rename(self.clairs / "input_manifest.v2-not-used.json")
        samples = []
        for sample in reversed(SAMPLES):
            caller_pass = {"path": f"/{sample}/caller-pass.vcf.gz"}
            longphase_pass = {"path": f"/{sample}/longphase-pass.vcf.gz"}
            samples.append({
                "sample": sample,
                "somatic": {
                    "caller_raw_vcf": {},
                    "longphase_input_vcf": {},
                    "caller_pass_baseline_vcf": caller_pass,
                    "longphase_recalibrated_all_vcf": {"path": f"/{sample}/all.vcf.gz"},
                    "longphase_recalibrated_pass_vcf": longphase_pass,
                    "tree_vcf": caller_pass,
                },
            })
        write_json(self.clairs / "frozen_input_lock.json", {
            "all_pass": True,
            "analysis_contract": {
                "task_type": "backbone_sensitivity",
                "tree_input_contract": "clairs_FILTER_PASS_sensitivity",
                "production_filter_policy": {"name": "production_default_filter"},
            },
            "samples": samples,
        })
        verification_path = self.clairs / "verification_summary.json"
        verification = json.loads(verification_path.read_text(encoding="utf-8"))
        verification["samples"] = list(reversed(verification["samples"]))
        write_json(verification_path, verification)
        output = self.root / "v3_sensitivity.json"
        proc = self._run(output)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(json.loads(output.read_text())["aggregate"]["verdict"], "robust_all_gates")


if __name__ == "__main__":
    unittest.main()
