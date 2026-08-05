#!/usr/bin/env python3
"""Fail-closed lifecycle tests for the canonical layered runner."""

import json
import hashlib
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import pysam


HERE = Path(__file__).resolve().parent
RUNNER = HERE / "run_layered_7samples_newbb.sh"
BASE_MANIFEST = HERE.parent / "data" / "layered_v2_input_manifest.json"


def clean_contract(manifest):
    manifest = dict(manifest)
    manifest["schema_version"] = "2.1"
    manifest["task_type"] = "B_COMPREHENSIVE_VALIDATION"
    manifest["tree_input_contract"] = "longphase_recalibrated_PASS"
    manifest["tag_contract"] = {
        "truth_flags": False,
        "PS_preserved": True,
        "bam_identity_locked": True,
        "tree_backbone": "LongPhase-S _sc.vcf FILTER=PASS",
        "longphase_filtering_policy": "production_default_filter",
    }
    return manifest


def run_fixture(root, manifest):
    manifest_path = root.parent / f"{root.name}.manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    env = os.environ.copy()
    env.update({"SM_INPUT_MANIFEST": str(manifest_path), "SM_RUN_ROOT": str(root)})
    return subprocess.run(["bash", str(RUNNER)], env=env, text=True, capture_output=True, check=False)


def write_vcf(path):
    plain = path.with_suffix("").with_suffix("")
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=ClairS\n"
        "##clairs_version=synthetic\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##FILTER=<ID=PASS,Description=All filters passed>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n"
        "chr1\t101\t.\tA\tC\t60\tPASS\t.\tGT\t0/1\t0/0\n"
        "chr1\t151\t.\tA\tG\t60\tPASS\t.\tGT\t0/1\t0/0\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def make_alignment(name, first_alt, second_alt):
    sequence = list("A" * 100)
    if first_alt:
        sequence[20] = "C"
    if second_alt:
        sequence[70] = "G"
    alignment = pysam.AlignedSegment()
    alignment.query_name = name
    alignment.query_sequence = "".join(sequence)
    alignment.flag = 0
    alignment.reference_id = 0
    alignment.reference_start = 80
    alignment.mapping_quality = 60
    alignment.cigarstring = "100M"
    alignment.query_qualities = pysam.qualitystring_to_array("I" * 100)
    return alignment


def build_success_fixture(parent):
    bam = parent / "tumor.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    specifications = []
    for label, first, second, hp in (
            ("rr1", False, False, "1-1"), ("ar1", True, False, "1-1"),
            ("aa1", True, True, "1-1"), ("rr2", False, False, "2-1")):
        specifications.extend((f"{label}_{index}", first, second, hp) for index in range(3))
    with pysam.AlignmentFile(str(bam), "wb", header=header) as output:
        for name, first, second, _hp in specifications:
            output.write(make_alignment(name, first, second))
    pysam.index(str(bam))

    sidecar_plain = parent / "read_tags.tsv"
    with pysam.AlignmentFile(str(bam), "rb") as source, sidecar_plain.open("w", encoding="utf-8") as output:
        output.write("#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n")
        for alignment, specification in zip(source.fetch(until_eof=True), specifications):
            hp = specification[3]
            digest = hashlib.blake2b(alignment.cigarstring.encode(), digest_size=8).hexdigest()
            output.write(
                f"chr1\t{alignment.reference_start}\t{alignment.reference_end}\t{alignment.query_name}\t"
                f"{alignment.flag}\t{alignment.mapping_quality}\t{digest}\t{hp}\t100\n"
            )
    sidecar = parent / "read_tags.tsv.gz"
    pysam.tabix_compress(str(sidecar_plain), str(sidecar), force=True)
    pysam.tabix_index(str(sidecar), seq_col=0, start_col=1, end_col=2, zerobased=True, force=True)

    raw = parent / "raw.vcf.gz"
    longphase_input = parent / "longphase_input.vcf.gz"
    recalibrated = parent / "longphase_sc_all.vcf.gz"
    tree = parent / "longphase_sc_pass.vcf.gz"
    for path in (raw, longphase_input, recalibrated, tree):
        write_vcf(path)

    validation = parent / "sidecar_validation.json"
    validation.write_text(json.dumps({
        "schema_version": "1.0",
        "region": "all",
        "sidecar": str(sidecar_plain),
        "HP_counts": {"1-1": 9, "2-1": 3},
        "duplicate_exact_alignment_rows": 0,
        "duplicate_exact_alignment_conflicts": 0,
        "checks": {
            "truth_flags_absent": True,
            "parser_count_matches_input": True,
            "sidecar_row_count_matches_capture": True,
            "tagged_count_matches_execution": True,
            "sidecar_no_unknown_HP": True,
            "sidecar_no_exact_identity_conflicts": True,
            "recalibrated_preserves_all_input_keys": True,
        },
        "pass": True,
    }), encoding="utf-8")
    inventory = parent / "input_files.tsv"
    stat_mtime = subprocess.check_output(["stat", "-c", "%y", str(bam)], text=True).strip()
    inventory.write_text(
        f"role\tpath\tsize_bytes\tmtime\n"
        f"tumor_bam\t{bam}\t{bam.lstat().st_size}\t{stat_mtime}\n",
        encoding="utf-8",
    )
    logical = bam.lstat()
    resolved = bam.resolve()
    target = resolved.stat()
    fingerprint = {
        "logical_path": str(bam),
        "logical_size_bytes": logical.st_size,
        "logical_mtime_ns": logical.st_mtime_ns,
        "resolved_path": str(resolved),
        "resolved_size_bytes": target.st_size,
        "resolved_mtime_ns": target.st_mtime_ns,
    }
    manifest = clean_contract(json.loads(BASE_MANIFEST.read_text(encoding="utf-8")))
    producer_dir = parent / "producer_closeout"
    producer_dir.mkdir()
    closeout_path = producer_dir / "production_closeout.json"
    artifacts_path = producer_dir / "artifacts.final.sha256"
    success_path = parent / "producer_SUCCESS.json"
    closeout_path.write_text(json.dumps({
        "status": "PASS", "dataset_count": 7, "n_pass": 7, "truth_flags": False,
        "tree_backbone": "LongPhase-S _sc.vcf FILTER=PASS",
    }), encoding="utf-8")
    artifacts_path.write_text("synthetic artifact manifest\n", encoding="utf-8")
    success_path.write_text(json.dumps({
        "status": "SUCCESS",
        "closeout_receipt": str(closeout_path),
        "closeout_receipt_sha256": hashlib.sha256(closeout_path.read_bytes()).hexdigest(),
        "artifacts_manifest": str(artifacts_path),
        "artifacts_manifest_sha256": hashlib.sha256(artifacts_path.read_bytes()).hexdigest(),
    }), encoding="utf-8")
    manifest["tag_contract"].update({
        "production_closeout": str(closeout_path),
        "production_closeout_sha256": hashlib.sha256(closeout_path.read_bytes()).hexdigest(),
        "production_success": str(success_path),
        "production_success_sha256": hashlib.sha256(success_path.read_bytes()).hexdigest(),
        "production_artifacts_manifest": str(artifacts_path),
        "production_artifacts_manifest_sha256": hashlib.sha256(artifacts_path.read_bytes()).hexdigest(),
    })
    for sample in manifest["samples"]:
        sample.update({
            "tumor_bam": str(bam),
            "tumor_bam_fingerprint": fingerprint,
            "normal_bam": str(bam),
            "caller_raw_vcf": str(raw),
            "longphase_input_vcf": str(longphase_input),
            "longphase_recalibrated_all_vcf": str(recalibrated),
            "longphase_recalibrated_pass_vcf": str(tree),
            "somatic_vcf": str(tree),
            "somatic_vcf_role": "LongPhase-S recalibrated FILTER=PASS; synthetic tree universe",
            "read_tag_sidecar": str(sidecar),
            "read_tag_sidecar_index": str(sidecar) + ".tbi",
            "read_tag_validation": str(validation),
            "longphase_input_inventory": str(inventory),
            "longphase_production_closeout": str(closeout_path),
            "longphase_tagging_scope": "production genome-wide; no truth VCF/BED flags",
            "cn_bed": None,
            "cn_source": "unavailable",
            "cn_int_gain": None,
            "cn_int_loss": None,
            "integration_json": None,
        })
    return manifest


class LayeredRunnerFailClosedTest(unittest.TestCase):
    def test_duplicate_dataset_binding_is_rejected_before_run_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            manifest = build_success_fixture(parent)
            manifest["samples"][-1]["sample"] = "HCC1395"
            root = parent / "duplicate_binding"
            result = run_fixture(root, manifest)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(root.exists())

    def test_tiny_clairs_backbone_sensitivity_is_explicit_and_runnable(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            manifest = build_success_fixture(parent)
            manifest["task_type"] = "B_BACKBONE_SENSITIVITY"
            manifest["tree_input_contract"] = "clairs_PASS_input"
            manifest["tag_contract"]["tree_backbone"] = "ClairS PASS sensitivity"
            for sample in manifest["samples"]:
                sample["somatic_vcf"] = sample["longphase_input_vcf"]
                sample["somatic_vcf_role"] = "ClairS PASS sensitivity; not canonical"
            root = parent / "sensitivity"
            result = run_fixture(root, manifest)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertTrue((root / "_SUCCESS").is_file())
            site_summary = json.loads(
                (root / "samples" / "HCC1395" / "ssnv_site_ledger_HCC1395.summary.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(site_summary["tree_contract"], "clairs_PASS_input")
            self.assertTrue(site_summary["checks"]["longphase_input_equals_tree_input"])

    def test_tiny_canonical_run_freezes_inputs_and_commits_success_last(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            manifest = build_success_fixture(parent)
            root = parent / "success"
            result = run_fixture(root, manifest)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertTrue((root / "_SUCCESS").is_file())
            self.assertFalse((root / "_FAILED").exists())
            verification = json.loads((root / "verification_summary.json").read_text(encoding="utf-8"))
            self.assertTrue(verification["all_pass"])
            self.assertEqual(verification["n_pass"], 7)
            receipt = json.loads((root / "launch_receipt.json").read_text(encoding="utf-8"))
            self.assertTrue(receipt["workers_read_frozen_manifest"])
            self.assertTrue(receipt["workers_read_source_bundle"])
            self.assertEqual(receipt["frozen_manifest"], str(root / "input_manifest.json"))
            hash_check = subprocess.run(
                ["sha256sum", "-c", str(root / "source_bundle.sha256")],
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(hash_check.returncode, 0, hash_check.stderr + hash_check.stdout)

    def test_empty_dataset_is_rejected_without_run_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            manifest = clean_contract(json.loads(BASE_MANIFEST.read_text(encoding="utf-8")))
            manifest.update({"dataset_count": 0, "biological_sample_count": 0, "samples": []})
            root = parent / "empty"
            result = run_fixture(root, manifest)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(root.exists())

    def test_missing_preflight_inputs_do_not_create_run_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            manifest = clean_contract(json.loads(BASE_MANIFEST.read_text(encoding="utf-8")))
            for sample in manifest["samples"]:
                sample.update({
                    "read_tag_sidecar": str(parent / "missing.tags.tsv.gz"),
                    "read_tag_sidecar_index": str(parent / "missing.tags.tsv.gz.tbi"),
                    "read_tag_validation": str(parent / "missing.validation.json"),
                    "caller_raw_vcf": str(parent / "missing.raw.vcf.gz"),
                    "longphase_input_vcf": str(parent / "missing.input.vcf.gz"),
                    "longphase_recalibrated_all_vcf": str(parent / "missing.recal.vcf.gz"),
                    "longphase_tagging_scope": "production genome-wide; no truth VCF/BED flags",
                })
            root = parent / "missing"
            result = run_fixture(root, manifest)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(root.exists())
            self.assertFalse((root / "_SUCCESS").exists())

    def test_existing_run_root_is_never_reused(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            manifest = clean_contract(json.loads(BASE_MANIFEST.read_text(encoding="utf-8")))
            root = parent / "existing"
            root.mkdir()
            sentinel = root / "sentinel.txt"
            sentinel.write_text("unchanged\n", encoding="utf-8")
            result = run_fixture(root, manifest)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "unchanged\n")
            self.assertFalse((root / "_SUCCESS").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
