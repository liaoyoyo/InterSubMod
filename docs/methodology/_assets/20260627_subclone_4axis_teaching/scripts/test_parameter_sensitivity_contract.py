#!/usr/bin/env python3
"""Regression tests for full-scope parameter sensitivity and topology gates."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "summarize_parameter_sensitivity.py"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
VARIANTS = {
    "mapq30": {"params": {"MAPQ_MIN": 30}},
    "baseq10": {"params": {"BASEQ_MIN": 10}},
    "minread4": {"MINREAD": 4},
    "maxsnv6": {"params": {"MAX_SNV": 6}},
}


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def view():
    return {
        "census": {
            "n_regions": 2,
            "n_regions_with_primary_lineage": 2,
            "L1": {"n_primary_lineage_units": 2,
                   "determinacy_primary_lineage": {"determined": 1}},
            "hp_multiplicity": {"2": 1},
            "region_determinacy": {"all_determined": 1},
            "funnel": {"L6_retained_sSNV": 2},
        }
    }


class ParameterSensitivityContractTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.manifest = self.root / "canonical.json"
        self.manifest_doc = {
            "schema_version": "2.1",
            "tree_input_contract": "longphase_recalibrated_PASS",
            "samples": [{"sample": sample} for sample in SAMPLES],
        }
        write_json(self.manifest, self.manifest_doc)
        self.base = self.root / "base"
        self._make_run(self.base, None)
        self.parameter_root = self.root / "parameters"
        for variant, params in VARIANTS.items():
            run = self.parameter_root / variant / "run"
            self._make_run(run, params)

    def tearDown(self):
        self.temp.cleanup()

    def _make_run(self, root, params):
        root.mkdir(parents=True)
        (root / "_SUCCESS").write_text("PASS\n", encoding="utf-8")
        write_json(root / "input_manifest.json", self.manifest_doc)
        write_json(root / "verification_summary.json", {
            "all_pass": True, "n_pass": 7,
            "samples": [{"sample": sample} for sample in SAMPLES],
        })
        write_json(root / "params.json", params or {"params": {"MAPQ_MIN": 20}})
        for index, sample in enumerate(SAMPLES, start=1):
            sample_dir = root / "samples" / sample
            write_json(sample_dir / f"layered_region_view_{sample}.json", view())
            write_json(sample_dir / "mlhp_part_1.json", {
                "groups": [{"chrom": "chr1", "positions": [100 + index, 200 + index]}]
            })
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
            write_json(sample_dir / f"ssnv_site_ledger_{sample}.summary.json", {"pass": True})

    def run_summary(self, output):
        return subprocess.run([
            sys.executable, str(SCRIPT),
            "--base-run", str(self.base),
            "--parameter-root", str(self.parameter_root),
            "--input-manifest", str(self.manifest),
            "--output", str(output),
        ], text=True, capture_output=True, check=False)

    def test_all_four_full_scope_variants_are_robust_when_identical(self):
        output = self.root / "summary.json"
        proc = self.run_summary(output)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        summary = json.loads(output.read_text(encoding="utf-8"))
        self.assertEqual(set(summary["aggregate"]), set(VARIANTS))
        self.assertTrue(all(value["verdict"] == "robust_all_gates"
                            for value in summary["aggregate"].values()))

    def test_topology_change_forces_sensitive_verdict(self):
        path = self.parameter_root / "mapq30" / "run" / "samples" / SAMPLES[0] / f"layered_reconstruction_{SAMPLES[0]}.json"
        document = json.loads(path.read_text(encoding="utf-8"))
        document["detail"][0]["analysis_tree_digest_sha256"] = "b" * 64
        write_json(path, document)
        output = self.root / "topology.json"
        proc = self.run_summary(output)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        summary = json.loads(output.read_text(encoding="utf-8"))
        self.assertEqual(summary["aggregate"]["mapq30"]["verdict"], "sensitive")

    def test_missing_variant_success_marker_is_rejected(self):
        marker = self.parameter_root / "baseq10" / "run" / "_SUCCESS"
        marker.rename(marker.with_name("_SUCCESS.missing"))
        output = self.root / "must_not_exist.json"
        proc = self.run_summary(output)
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse(output.exists())

    def test_v3_canonical_base_lock_is_accepted(self):
        manifest_path = self.base / "input_manifest.json"
        manifest_path.rename(self.base / "input_manifest.v2-not-used.json")
        write_json(self.base / "frozen_input_lock.json", {
            "all_pass": True,
            "analysis_contract": {"production_filter_policy": {"name": "production_default_filter"}},
            "samples": [{
                "sample": sample,
                "somatic": {
                    "caller_raw_vcf": {}, "longphase_input_vcf": {},
                    "longphase_recalibrated_all_vcf": {"path": f"/{sample}/all.vcf.gz"},
                    "tree_vcf": {"path": f"/{sample}/pass.vcf.gz"},
                },
            } for sample in reversed(SAMPLES)],
        })
        verification_path = self.base / "verification_summary.json"
        verification = json.loads(verification_path.read_text(encoding="utf-8"))
        verification["samples"] = list(reversed(verification["samples"]))
        write_json(verification_path, verification)
        output = self.root / "v3_base.json"
        proc = self.run_summary(output)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(all(item["verdict"] == "robust_all_gates"
                            for item in json.loads(output.read_text())["aggregate"].values()))

    def test_v3_variant_reads_frozen_launch_analysis_params(self):
        run = self.parameter_root / "mapq30" / "run"
        params = run / "params.json"
        params.rename(run / "params.v2-not-used.json")
        write_json(run / "launch_receipt.json", {
            "extra": {"analysis_params": {"MAPQ_MIN": 30}}
        })
        output = self.root / "v3_variant.json"
        proc = self.run_summary(output)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(json.loads(output.read_text())["aggregate"]["mapq30"]["verdict"],
                         "robust_all_gates")


if __name__ == "__main__":
    unittest.main()
