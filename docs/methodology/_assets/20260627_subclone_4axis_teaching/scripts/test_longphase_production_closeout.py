#!/usr/bin/env python3
"""Synthetic lifecycle tests for LongPhase-S production closeout."""

import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pysam


SCRIPT_DIR = Path(__file__).resolve().parent
FINALIZER = SCRIPT_DIR / "finalize_longphase_production_sidecars.py"
PREPARER = SCRIPT_DIR / "prepare_clean_layered_manifest.py"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
BIOLOGICAL = {sample: ("HCC1395" if sample == "HCC1395_DORADO" else sample) for sample in SAMPLES}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def write_hash_manifest(path, artifacts):
    path.write_text("".join(f"{sha256(item)}  {item}\n" for item in artifacts), encoding="utf-8")


def write_vcf(path, position):
    plain = path.with_suffix("").with_suffix("")
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "##FILTER=<ID=PASS,Description=All filters passed>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t{position}\t.\tA\tC\t60\tPASS\t.\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True, csi=True)
    return plain


class ProductionCloseoutTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.tmp = Path(self.temp.name)
        self.root = self.tmp / "run"
        self.root.mkdir()
        self.binary = self.tmp / "longphase-s"
        self.binary.write_bytes(b"synthetic longphase binary\n")
        self.binary.chmod(0o755)
        self.binary_sha = sha256(self.binary)
        source = self.tmp / "producer.py"
        source.write_text("# frozen producer\n", encoding="utf-8")
        manifest = {"samples": []}
        verification_samples = []
        status_rows = []
        for index, name in enumerate(SAMPLES, start=1):
            input_dir = self.tmp / "inputs" / name
            input_dir.mkdir(parents=True)
            germline = input_dir / "germline.vcf.gz"
            germline.write_text("germline\n", encoding="utf-8")
            (input_dir / "germline.vcf.gz.csi").write_text("index\n", encoding="utf-8")
            normal = input_dir / "normal.bam"
            tumor = input_dir / "tumor.bam"
            if name == "HCC1395_DORADO":
                normal_target = input_dir / "normal.target.bam"
                tumor_target = input_dir / "tumor.target.bam"
                normal_target.write_text("normal target payload\n", encoding="utf-8")
                tumor_target.write_text("tumor target payload\n", encoding="utf-8")
                normal.symlink_to(normal_target.name)
                tumor.symlink_to(tumor_target.name)
            else:
                normal.write_text("normal\n", encoding="utf-8")
                tumor.write_text("tumor\n", encoding="utf-8")
            (input_dir / "normal.bam.bai").write_text("index\n", encoding="utf-8")
            (input_dir / "tumor.bam.bai").write_text("index\n", encoding="utf-8")
            reference = input_dir / "ref.fa"
            reference.write_text(">chr1\nA\n", encoding="utf-8")
            caller = input_dir / "caller.vcf.gz"
            write_vcf(caller, 100 + index)
            sample = {
                "sample": name,
                "biological_id": BIOLOGICAL[name],
                "germline_phased_vcf": str(germline),
                "normal_bam": str(normal),
                "tumor_bam": str(tumor),
                "reference": str(reference),
                "caller_raw_vcf": str(caller),
                "caller_pass_vcf": str(caller),
                "caller_pass_records": 1,
            }
            manifest["samples"].append(sample)
            wd = self.root / "samples" / name
            wd.mkdir(parents=True)
            inventory = wd / "input_files.tsv"
            with inventory.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["role", "path", "size_bytes", "mtime"])
                for role, artifact in (("germline", germline), ("normal_bam", normal),
                                       ("tumor_bam", tumor), ("reference", reference),
                                       ("caller_pass", caller)):
                    mtime = subprocess.run(["stat", "-c", "%y", str(artifact)], text=True,
                                           capture_output=True, check=True).stdout.strip()
                    writer.writerow([role, artifact, artifact.lstat().st_size, mtime])
            write_hash_manifest(wd / "input.sha256", [caller, germline, normal, tumor, reference])

            sidecar = wd / f"{name}.read_tags.tsv"
            sidecar.write_text(
                "#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n"
                f"chr1\t{99 + index}\t{100 + index}\tread{name}\t0\t60\t0123456789abcdef\t1-1\t10\n",
                encoding="utf-8",
            )
            sidecar_gz = Path(f"{sidecar}.gz")
            pysam.tabix_compress(str(sidecar), str(sidecar_gz), force=True)
            pysam.tabix_index(str(sidecar_gz), seq_col=0, start_col=1, end_col=2,
                              zerobased=True, force=True)
            vcf_all = wd / f"{name}.longphase_s.recalibrated.all.vcf.gz"
            vcf_pass = wd / f"{name}.longphase_s.recalibrated.pass.vcf.gz"
            write_vcf(vcf_all, 100 + index)
            write_vcf(vcf_pass, 100 + index)
            raw_vcf = wd / f"{name}_production_sc.vcf"
            raw_vcf.write_text("synthetic raw VCF receipt\n", encoding="utf-8")
            purity = wd / f"{name}_production_purity.out"
            purity.write_text("1.0\n", encoding="utf-8")
            fifo = wd / "consumed_tagged_bam.fifo"
            os.mkfifo(fifo)
            capture = {
                "pass": True,
                "rows_mapped": 1,
                "alignment_classes": {"total": 1, "primary": 1},
            }
            write_json(wd / "stream_capture_summary.json", capture)
            validation = {
                "pass": True,
                "region": "all",
                "sidecar": str(sidecar),
                "capture": capture,
                "HP_counts": {"1-1": 1},
                "HP_with_PS_counts": {"1-1": 1},
                "duplicate_exact_alignment_rows": 0,
                "duplicate_exact_alignment_conflicts": 0,
                "record_key_missing": 0,
                "record_key_extra": 0,
                "checks": {"truth_flags_absent": True, "recalibrated_preserves_all_input_keys": True},
            }
            write_json(wd / "sidecar_validation.json", validation)
            (wd / "command.sh.txt").write_text("longphase-s synthetic\n", encoding="utf-8")
            (wd / "longphase_s.production.log").write_text("complete\n", encoding="utf-8")
            (wd / "stream_capture.log").write_text("complete\n", encoding="utf-8")
            output_artifacts = [sidecar, sidecar_gz, vcf_all, vcf_pass, wd / "sidecar_validation.json"]
            write_hash_manifest(wd / "output.sha256", output_artifacts)
            verification_samples.append({"pass": True, "sidecar": str(sidecar)})
            status_rows.extend([
                ["2026-07-11T00:00:00+08:00", name, "production_tagging", "START", ""],
                ["2026-07-11T00:00:01+08:00", name, "production_tagging", "PASS", "rows=1"],
            ])

        self.manifest = self.tmp / "manifest.json"
        write_json(self.manifest, manifest)
        self.base_manifest = self.tmp / "base.json"
        write_json(self.base_manifest, {
            "samples": [{"sample": item["sample"], "biological_id": item["biological_id"],
                         "tumor_bam": item["tumor_bam"]} for item in manifest["samples"]]
        })
        write_json(self.root / "input_manifest.json", manifest)
        write_json(self.root / "verification_summary.json", {
            "dataset_count": 7, "n_pass": 7, "all_pass": True, "samples": verification_samples,
        })
        write_json(self.root / "params.json", {"truth_flags": False})
        write_json(self.root / "runtime_executable_receipt.json", {
            "binary": {"sha256": self.binary_sha},
            "all_active_proc_executables_match_path": True,
        })
        write_hash_manifest(self.root / "code.sha256", [self.manifest, source])
        with (self.root / "run_status.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["timestamp", "sample", "stage", "status", "detail"])
            writer.writerows(status_rows)
            writer.writerow(["2026-07-11T00:00:02+08:00", "ALL", "verify", "PASS", "7/7"])

    def tearDown(self):
        self.temp.cleanup()

    def run_closeout(self, expected_sha=None):
        return subprocess.run([
            sys.executable, str(FINALIZER),
            "--run-root", str(self.root),
            "--expected-manifest", str(self.manifest),
            "--longphase-binary", str(self.binary),
            "--expected-binary-sha256", expected_sha or self.binary_sha,
        ], text=True, capture_output=True, check=False)

    def test_success_marker_is_published_after_strict_closeout(self):
        proc = self.run_closeout()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((self.root / "closeout" / "production_closeout.json").is_file())
        self.assertTrue((self.root / "closeout" / "artifacts.final.sha256").is_file())
        self.assertTrue((self.root / "_SUCCESS").is_file())
        receipt = json.loads((self.root / "closeout" / "production_closeout.json").read_text())
        self.assertEqual(receipt["n_pass"], 7)
        self.assertTrue(all(item["record_key_missing"] == 0 for item in receipt["samples"]))
        for line in (self.root / "closeout" / "artifacts.final.sha256").read_text().splitlines():
            digest, raw_path = line.split(None, 1)
            artifact = Path(raw_path.strip())
            self.assertTrue(artifact.is_file(), artifact)
            self.assertEqual(sha256(artifact), digest)

    def test_input_mutation_fails_without_success_marker(self):
        caller = Path(json.loads(self.manifest.read_text())["samples"][0]["caller_pass_vcf"])
        caller.write_bytes(caller.read_bytes() + b"drift")
        proc = self.run_closeout()
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse((self.root / "_SUCCESS").exists())

    def test_binary_mismatch_fails_without_success_marker(self):
        proc = self.run_closeout("0" * 64)
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse((self.root / "_SUCCESS").exists())

    def test_clean_manifest_accepts_verified_closeout(self):
        self.assertEqual(self.run_closeout().returncode, 0)
        output = self.tmp / "clean.json"
        proc = subprocess.run([
            sys.executable, str(PREPARER),
            "--base-manifest", str(self.base_manifest),
            "--longphase-manifest", str(self.manifest),
            "--production-root", str(self.root),
            "--output", str(output),
        ], text=True, capture_output=True, check=False)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        clean = json.loads(output.read_text())
        self.assertEqual(clean["tree_input_contract"], "longphase_recalibrated_PASS")
        self.assertTrue(all(item["longphase_production_closeout"].endswith("production_closeout.json")
                            for item in clean["samples"]))

    def test_clean_manifest_rejects_unclosed_producer(self):
        output = self.tmp / "must_not_exist.json"
        proc = subprocess.run([
            sys.executable, str(PREPARER),
            "--base-manifest", str(self.base_manifest),
            "--longphase-manifest", str(self.manifest),
            "--production-root", str(self.root),
            "--output", str(output),
        ], text=True, capture_output=True, check=False)
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse(output.exists())
        self.assertIn("production root is incomplete", proc.stderr)

    def test_clean_manifest_rejects_post_closeout_artifact_drift(self):
        self.assertEqual(self.run_closeout().returncode, 0)
        sidecar = self.root / "samples" / SAMPLES[0] / f"{SAMPLES[0]}.read_tags.tsv.gz"
        sidecar.write_bytes(sidecar.read_bytes() + b"drift")
        output = self.tmp / "must_not_exist_after_drift.json"
        proc = subprocess.run([
            sys.executable, str(PREPARER),
            "--base-manifest", str(self.base_manifest),
            "--longphase-manifest", str(self.manifest),
            "--production-root", str(self.root),
            "--output", str(output),
        ], text=True, capture_output=True, check=False)
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse(output.exists())
        self.assertIn("changed after closeout", proc.stderr)


if __name__ == "__main__":
    unittest.main()
