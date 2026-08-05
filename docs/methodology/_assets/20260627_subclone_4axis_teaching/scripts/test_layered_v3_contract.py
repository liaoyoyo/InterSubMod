#!/usr/bin/env python3
"""Tiny adversarial tests for the layered v3 source-manifest and frozen-lock contract."""

from __future__ import annotations

import argparse
import contextlib
import copy
import gzip
import io
import json
import os
import platform
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path
from typing import Any

import prepare_clean_layered_manifest_v3 as prepare
import validate_layered_v3_inputs as validator


SCRIPT_DIR = Path(__file__).resolve().parent
SCHEMA_PATH = SCRIPT_DIR.parent / "schemas" / "layered_input_manifest_v3.schema.json"
LOCK_SCHEMA_PATH = SCRIPT_DIR.parent / "schemas" / "layered_input_lock_v1.schema.json"
PRODUCER_RECEIPT_SCHEMA_PATH = SCRIPT_DIR.parent / "schemas" / "longphase_raw_all_capture_receipt_v2.schema.json"


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_sidecar(path: Path, *, duplicate: bool = False) -> None:
    rows = [
        "chr1\t0\t10\treadA\t0\t60\t1111111111111111\t1\t100",
        "chr1\t20\t30\treadB\t16\t50\t2222222222222222\t2\t100",
    ]
    if duplicate:
        rows.insert(1, rows[0])
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n")
        handle.write("\n".join(rows) + "\n")


def make_method_receipts(root: Path) -> tuple[Path, Path, Path]:
    bounded_sidecar = root / "bounded.tsv.gz"
    write_sidecar(bounded_sidecar)
    bounded_index = Path(f"{bounded_sidecar}.tbi")
    bounded_index.write_bytes(b"tiny-tabix")
    payload_sha = validator.file_sha256(root / "seed.txt")
    real_receipt = root / "equivalence_probe.json"
    comparisons = {
        "payload_equal": True,
        "payload_sha256_equal": True,
        "calls_sha256_equal": True,
        "hp_sha256_equal": True,
        "ps_sha256_equal": True,
        "sidecar_exact_matches_all_exposures": True,
        "sidecar_missing_zero": True,
        "sidecar_conflicts_zero": True,
    }
    write_json(real_receipt, {
        "schema_version": "1.0",
        "scope": "PARTIAL_BOUNDED_REAL_DATA_PROBE",
        "claim_limit": "join equivalence only; historical tags are not clean production-tag evidence",
        "identity_schema": "coordinate_join_v1",
        "region": {"chrom": "chr1", "start": 1, "end": 30, "positions": [1, 30]},
        "sidecar": {
            "path": str(bounded_sidecar),
            "sha256": validator.file_sha256(bounded_sidecar),
            "index_path": str(bounded_index),
            "index_sha256": validator.file_sha256(bounded_index),
            "counts": {"rows": 2},
        },
        "direct": {"payload_sha256": payload_sha},
        "joined": {"payload_sha256": payload_sha},
        "comparisons": comparisons,
        "pass": True,
    })
    test_source = root / "join_tests.py"
    consumer_source = root / "consumer.py"
    test_source.write_text("# synthetic join tests\n", encoding="utf-8")
    consumer_source.write_text("# synthetic consumer\n", encoding="utf-8")
    synthetic_receipt = root / "synthetic_contract_receipt.json"
    write_json(synthetic_receipt, {
        "schema_version": "1.0",
        "scope": "SYNTHETIC_COORDINATE_JOIN_V1_CONTRACT",
        "claim_limit": "join edge-case behavior only; not production-data validation",
        "test_source": {"path": str(test_source), "sha256": validator.file_sha256(test_source)},
        "consumer_source": {"path": str(consumer_source), "sha256": validator.file_sha256(consumer_source)},
        "tests": [
            {"name": "test_exact_sidecar_matches_direct_tagged_bam", "pass": True},
            {"name": "test_missing_alignment_is_detected", "pass": True},
            {"name": "test_duplicate_identity_conflict_is_detected", "pass": True},
            {"name": "test_redundant_identical_sidecar_row_collapses_without_conflict", "pass": True},
            {"name": "test_duplicate_bam_identity_with_conflicting_allele_is_excluded", "pass": True},
        ],
        "tests_run": 5,
        "failures": 0,
        "errors": 0,
        "exit_code": 0,
        "pass": True,
    })
    return real_receipt, synthetic_receipt, test_source


def make_fixture(root: Path) -> tuple[dict[str, Any], Path]:
    root.mkdir(parents=True, exist_ok=True)
    (root / "seed.txt").write_text("stable-payload\n", encoding="utf-8")
    production_root = root / "production"
    executable = root / "longphase-s"
    runner_source = root / "run_longphase_raw_all.sh"
    capture_source = root / "capture.py"
    validator_source = root / "validate.py"
    filter_source = root / "filter_audit.py"
    patch_source = root / "zero_read.patch"
    patch_receipt_path = root / "patch_receipt.json"
    environment_lock = root / "environment.lock.json"
    code_hash_manifest = root / "code.sha256"
    normalizer_source = root / "normalizer.py"
    executable.write_bytes(b"tiny-longphase-binary")
    runner_source.write_text("# runner\n", encoding="utf-8")
    capture_source.write_text("# capture\n", encoding="utf-8")
    validator_source.write_text("# validator\n", encoding="utf-8")
    filter_source.write_text("# filter audit\n", encoding="utf-8")
    patch_source.write_text("diff --git a/x b/x\n", encoding="utf-8")
    write_json(patch_receipt_path, {
        "status": "APPROVED_FOR_FAIL_CLOSED_7_DATASET_VALIDATION",
        "patched_binary_sha256": validator.file_sha256(executable),
        "approval": {"scope": "FAIL_CLOSED_7_DATASET_VALIDATION_ONLY"},
        "patch": str(patch_source),
    })
    environment_lock.write_text('{"env":"tiny"}\n', encoding="utf-8")
    code_hash_manifest.write_text("tiny code manifest\n", encoding="utf-8")
    normalizer_source.write_text("# normalizer\n", encoding="utf-8")
    base_samples = []
    longphase_samples = []
    for sample, biological_id in validator.EXPECTED_DATASETS.items():
        source_dir = root / "source" / sample
        sample_dir = production_root / "samples" / sample
        source_dir.mkdir(parents=True)
        sample_dir.mkdir(parents=True)
        bam = source_dir / f"{sample}.bam"
        bam.write_bytes((sample.encode("ascii") + b"\0") * 32)
        bam_index = Path(f"{bam}.bai")
        bam_index.write_bytes(b"tiny-bai")
        storage = validator.storage_identity_v1(bam, bam_index)
        normal_bam = source_dir / f"{sample}.normal.bam"
        normal_bam.write_bytes((b"normal-" + sample.encode("ascii")) * 8)
        normal_index = Path(f"{normal_bam}.bai")
        normal_index.write_bytes(b"tiny-normal-bai")
        normal_storage = validator.storage_identity_v1(normal_bam, normal_index)
        germline = source_dir / f"{sample}.germline.phased.vcf.gz"
        germline.write_bytes(b"tiny-germline-vcf")
        Path(f"{germline}.tbi").write_bytes(b"tiny-germline-index")
        reference = source_dir / "reference.fa"
        reference.write_bytes(b">chr1\nACGT\n")
        reference_fai = Path(f"{reference}.fai")
        reference_fai.write_bytes(b"chr1\t4\t6\t4\t5\n")
        reference_storage = validator.storage_identity_v1(reference, reference_fai)
        vcfs: dict[str, Path] = {}
        for role in ("raw", "pass"):
            path = source_dir / f"{sample}.{role}.vcf.gz"
            path.write_bytes(f"{sample}-{role}\n".encode("ascii"))
            Path(f"{path}.tbi").write_bytes(b"tiny-vcf-index")
            vcfs[role] = path
        normalized = sample_dir / f"{sample}.clairs.raw_all.normalized.vcf.gz"
        normalized.write_bytes(f"{sample}-raw-normalized\n".encode("ascii"))
        Path(f"{normalized}.tbi").write_bytes(b"tiny-normalized-index")
        recal = sample_dir / f"{sample}.longphase_s.recalibrated.all.vcf.gz"
        recal.write_bytes(f"{sample}-recal\n".encode("ascii"))
        Path(f"{recal}.tbi").write_bytes(b"tiny-recal-index")
        recal_pass = sample_dir / f"{sample}.longphase_s.recalibrated.pass.vcf.gz"
        recal_pass.write_bytes(f"{sample}-recal-pass\n".encode("ascii"))
        Path(f"{recal_pass}.tbi").write_bytes(b"tiny-recal-pass-index")
        sidecar = sample_dir / f"{sample}.read_tags.tsv.gz"
        write_sidecar(sidecar)
        sidecar_index = Path(f"{sidecar}.tbi")
        sidecar_index.write_bytes(b"tiny-sidecar-index")
        argv = [
            str(executable), "somatic_haplotag", "-s", str(germline), "-b", str(normal_bam),
            "--tumor-snv-file", str(normalized), "--tumor-bam-file", str(bam), "-r", str(reference),
            "-t", "12", "-q", "20", "--tagSupplementary", "--output-somatic-vcf",
            "-o", str(sample_dir / f"{sample}_production"),
        ]
        validation = {
            "sample": sample,
            "region": "all",
            "pass": True,
            "duplicate_exact_alignment_rows": 0,
            "duplicate_exact_alignment_conflicts": 0,
            "capture": {"rows_mapped": 2},
            "checks": {
                "truth_flags_absent": True,
                "parser_count_matches_input": True,
                "sidecar_row_count_matches_capture": True,
                "tagged_count_matches_execution": True,
                "sidecar_no_unknown_HP": True,
                "sidecar_no_exact_identity_conflicts": True,
                "recalibrated_preserves_all_input_keys": True,
            },
        }
        native_validation_path = sample_dir / "sidecar_validation.json"
        write_json(native_validation_path, validation)
        stream_capture_path = sample_dir / "stream_capture_summary.json"
        write_json(stream_capture_path, {"rows_mapped": 2, "pass": True})
        normalization_document = {
            "input": {"record_count": 2},
            "output": {"record_count": 2},
            "filter_transition_counts": {"LowQual->LowQual": 1, "PASS->PASS": 1},
            "rescued_nonpass_to_pass": 0,
            "removed_pass_to_nonpass": 0,
            "pass": True,
        }
        transition_document = {
            "input": {"record_count": 2},
            "output": {"record_count": 2},
            "filter_transition_counts": {"LowQual->PASS": 1, "PASS->LowQual": 1},
            "rescued_nonpass_to_pass": 1,
            "removed_pass_to_nonpass": 1,
            "pass": True,
        }
        normalization_path = sample_dir / "normalization_audit.json"
        transition_path = sample_dir / "filter_transition_audit.json"
        sample_verification_path = sample_dir / "sample_verification.json"
        write_json(normalization_path, normalization_document)
        write_json(transition_path, transition_document)
        write_json(sample_verification_path, {
            "sample": sample,
            "expected_raw_records": 2,
            "normalization": normalization_document,
            "filter_transitions": transition_document,
            "pass": True,
        })
        consumed_fifo = sample_dir / "consumed_tagged_bam.fifo"
        os.mkfifo(consumed_fifo)
        producer_inputs = {
            "germline_phased_vcf": validator.indexed_artifact(germline),
            "normal_bam": {"path": str(normal_bam), "index_path": str(normal_index), "storage_identity_v1": normal_storage},
            "tumor_bam": {"path": str(bam), "index_path": str(bam_index), "storage_identity_v1": storage},
            "caller_raw_vcf": validator.indexed_artifact(vcfs["raw"]),
            "longphase_input_vcf": validator.indexed_artifact(normalized),
            "caller_pass_baseline_vcf": validator.indexed_artifact(vcfs["pass"]),
            "reference": {
                "path": str(reference),
                "fai_path": str(reference_fai),
                "storage_identity_v1": reference_storage,
            },
        }
        producer_receipt = {
            "schema_name": validator.RAW_ALL_PRODUCER_RECEIPT_SCHEMA,
            "schema_version": validator.RAW_ALL_PRODUCER_RECEIPT_VERSION,
            "sample": sample,
            "evidence_origin": "post_run_normalization_from_frozen_raw_all_execution_artifacts",
            "identity_schema": "coordinate_join_v1",
            "assurance": "bounded_coordinate_equivalence_not_global_payload_identity",
            "longphase_input_contract": validator.RAW_ALL_INPUT_CONTRACT,
            "tree_input_contract": validator.TREE_INPUT_CONTRACT,
            "duplicate_identity_policy": validator.REDUNDANT_IDENTITY_POLICY,
            "production_policy": validator.PRODUCTION_POLICY,
            "effective_options": validator.PRODUCTION_EFFECTIVE_OPTIONS,
            "command_argv": argv,
            "command_argv_sha256": validator.canonical_sha256(argv),
            "producer_inputs": producer_inputs,
            "producer_input_binding_sha256": validator.canonical_sha256(producer_inputs),
            "longphase_executable": {
                **validator.artifact(executable),
                "version": "longphase-s tiny-test",
            },
            "patch_evidence": {
                "approval_scope": "FAIL_CLOSED_7_DATASET_VALIDATION_ONLY",
                "patch_receipt": validator.artifact(patch_receipt_path),
                "source_patch": validator.artifact(patch_source),
            },
            "producer_code": {
                "runner": validator.artifact(runner_source),
                "capture": validator.artifact(capture_source),
                "validator": validator.artifact(validator_source),
                "filter_auditor": validator.artifact(filter_source),
            },
            "environment_lock": {
                "production_manifest": validator.artifact(environment_lock),
                "run_params": validator.artifact(environment_lock),
                "code_hash_manifest": validator.artifact(code_hash_manifest),
                "sample_input_inventory": validator.artifact(code_hash_manifest),
                "sample_input_hash_manifest": validator.artifact(code_hash_manifest),
                "sample_output_hash_manifest": validator.artifact(code_hash_manifest),
                "python_executable": validator.artifact(Path(sys.executable).resolve()),
                "python_version": sys.version,
                "platform": platform.platform(),
                "normalizer_source": validator.artifact(normalizer_source),
                "normalizer_argv": [str(normalizer_source), "--sample", sample],
            },
            "capture_outputs": {
                "sidecar": validator.artifact(sidecar),
                "sidecar_index": validator.artifact(sidecar_index),
                "native_validation": validator.artifact(native_validation_path),
                "stream_capture_summary": validator.artifact(stream_capture_path),
                "normalization_audit": validator.artifact(normalization_path),
                "filter_transition_audit": validator.artifact(transition_path),
                "sample_verification": validator.artifact(sample_verification_path),
                "longphase_recalibrated_all_vcf": validator.indexed_artifact(recal),
                "longphase_recalibrated_pass_vcf": validator.indexed_artifact(recal_pass),
            },
            "global_coordinate_counts": {
                "mapped_alignment_count": 2,
                "identity_unique_count": 2,
                "duplicate_count": 0,
                "conflict_count": 0,
            },
            "filter_transition_summary": {
                "input_record_count": 2,
                "output_record_count": 2,
                "rescued_nonpass_to_pass": 1,
                "removed_pass_to_nonpass": 1,
                "transition_counts": transition_document["filter_transition_counts"],
                "pass": True,
            },
            "bam_output_policy": {
                "transport": "named_fifo",
                "persisted_bam": False,
                "regular_bam_count": 0,
                "consumed_fifo_path": str(consumed_fifo),
                "is_fifo_at_closeout": True,
            },
        }
        write_json(sample_dir / "producer_capture_receipt_v2.json", producer_receipt)
        base_samples.append({
            "sample": sample,
            "biological_id": biological_id,
            "cn_bed": None,
            "cn_source": "unavailable",
            "cn_semantics": "missing; never interpreted neutral",
            "cn_int_gain": None,
            "cn_int_loss": None,
            "integration_json": None,
        })
        longphase_samples.append({
            "sample": sample,
            "tumor_bam": str(bam),
            "caller_raw_vcf": str(vcfs["raw"]),
            "caller_pass_vcf": str(vcfs["pass"]),
            "longphase_input_contract": validator.RAW_ALL_INPUT_CONTRACT,
            "platform": "ONT_TEST",
        })
    base_path = root / "base.json"
    longphase_path = root / "longphase.json"
    write_json(base_path, {
        "schema_version": "2.0",
        "dataset_count": 7,
        "biological_sample_count": 6,
        "samples": base_samples,
    })
    write_json(longphase_path, {"samples": longphase_samples})
    real_receipt, synthetic_receipt, test_source = make_method_receipts(root)
    manifest = prepare.build_manifest(
        base_path,
        longphase_path,
        production_root,
        real_receipt,
        synthetic_receipt,
        "tiny_layered_v3_fixture",
        created_at_utc="2026-07-11T00:00:00Z",
    )
    return manifest, test_source


class LayeredV3ContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="layered-v3-contract-")
        self.root = Path(self.temporary.name)
        self.manifest, self.synthetic_test_source = make_fixture(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_validator(self, manifest: dict[str, Any]) -> tuple[int, Path, Path, dict[str, Any] | None]:
        manifest_path = self.root / f"manifest-{len(list(self.root.glob('manifest-*.json')))}.json"
        lock_path = self.root / f"lock-{manifest_path.stem}.json"
        failure_path = self.root / f"failure-{manifest_path.stem}.json"
        write_json(manifest_path, manifest)
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            exit_code = validator.main([
                "--manifest", str(manifest_path),
                "--schema", str(SCHEMA_PATH),
                "--frozen-lock", str(lock_path),
                "--failure-report", str(failure_path),
            ])
        failure = json.loads(failure_path.read_text(encoding="utf-8")) if failure_path.exists() else None
        return exit_code, lock_path, failure_path, failure

    def assert_rejected(self, manifest: dict[str, Any], exit_code: int, stable_code: str) -> None:
        observed, lock, failure, document = self.run_validator(manifest)
        self.assertEqual(observed, exit_code)
        self.assertFalse(lock.exists())
        self.assertTrue(failure.is_file())
        self.assertEqual(document["error"]["code"], stable_code)
        self.assertFalse(document["valid_lock_written"])

    def rewrite_validation(self, manifest: dict[str, Any], index: int, validation: dict[str, Any]) -> None:
        item = manifest["samples"][index]
        validation_path = Path(item["read_tags"]["validation"]["path"])
        write_json(validation_path, validation)
        identity = validator.full_identity(validation_path)
        item["read_tags"]["validation"]["identity"] = identity
        item["read_tags"]["subject_binding"]["validation_sha256"] = identity["sha256"]
        manifest["production_summary"]["datasets"][index]["validation_sha256"] = identity["sha256"]

    def rewrite_producer_receipt(self, manifest: dict[str, Any], index: int, receipt: dict[str, Any]) -> None:
        item = manifest["samples"][index]
        path = Path(item["read_tags"]["producer_capture_receipt"]["path"])
        write_json(path, receipt)
        identity = validator.full_identity(path)
        item["read_tags"]["producer_capture_receipt"]["identity"] = identity
        item["read_tags"]["subject_binding"]["producer_capture_receipt_sha256"] = identity["sha256"]

    def test_valid_fixture_writes_only_frozen_lock(self) -> None:
        exit_code, lock_path, failure_path, _ = self.run_validator(self.manifest)
        self.assertEqual(exit_code, 0)
        self.assertTrue(lock_path.is_file())
        self.assertFalse(failure_path.exists())
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
        lock_sha = lock.pop("lock_sha256")
        self.assertEqual(validator.canonical_sha256(lock), lock_sha)
        self.assertEqual(len(lock["samples"]), 7)
        self.assertTrue(all(not item["alignment_payload"]["storage_identity_v1"]["is_full_content_hash"] for item in lock["samples"]))
        self.assertTrue(all(item["somatic"]["tree_vcf"] for item in lock["samples"]))
        self.assertTrue(all(
            set(item["read_tags"]["producer_evidence"]["producer_inputs"])
            == {
                "germline_phased_vcf", "normal_bam", "tumor_bam", "caller_raw_vcf",
                "longphase_input_vcf", "caller_pass_baseline_vcf", "reference",
            }
            for item in lock["samples"]
        ))
        unavailable = lock["samples"][0]["copy_number"]
        self.assertEqual(unavailable["unlisted_position_semantics"], "unavailable")
        self.assertTrue(unavailable["reason"])
        self.assertIsNone(unavailable["cn_bed"])

    def test_manifest_builder_scans_each_dataset_once(self) -> None:
        scan = validator.inspect_coordinate_sidecar
        with mock.patch.object(validator, "inspect_coordinate_sidecar", wraps=scan) as observed:
            make_fixture(self.root / "single-scan-fixture")
        self.assertEqual(observed.call_count, 7)

    def test_method_receipt_resolution_uses_explicit_repo_root_after_bundle_relocation(self) -> None:
        explicit_root = self.root / "InterSubMod"
        validator.set_method_receipt_repo_root(explicit_root)
        try:
            observed = validator._resolve_receipt_reference("InterSubMod/research/probe.json")
        finally:
            validator.set_method_receipt_repo_root(None)
        self.assertEqual(observed, explicit_root / "research/probe.json")

    def test_measured_cn_semantics_survive_frozen_lock(self) -> None:
        changed = copy.deepcopy(self.manifest)
        cn_bed = self.root / "measured-cn.bed"
        cn_bed.write_text("chr1\t0\t100\tgain\n", encoding="utf-8")
        changed["samples"][0]["copy_number"] = {
            "availability": "measured",
            "source": "TINY_REVIEWED_CN",
            "semantics": "non-neutral intervals; unlisted autosomal positions interpreted neutral",
            "coordinate_system": "0_based_half_open",
            "unlisted_position_semantics": "neutral",
            "allowed_states": ["gain", "loss", "loh", "neutral"],
            "overlap_policy": "forbid",
            "reason": None,
            "cn_bed": validator.artifact(cn_bed),
            "cn_int_gain": None,
            "cn_int_loss": None,
            "integration_json": None,
        }
        exit_code, lock_path, failure_path, _ = self.run_validator(changed)
        self.assertEqual(exit_code, 0)
        self.assertFalse(failure_path.exists())
        frozen = json.loads(lock_path.read_text(encoding="utf-8"))["samples"][0]["copy_number"]
        self.assertEqual(frozen["source"], "TINY_REVIEWED_CN")
        self.assertEqual(frozen["coordinate_system"], "0_based_half_open")
        self.assertEqual(frozen["unlisted_position_semantics"], "neutral")

    def test_clairs_pass_cannot_be_tree_backbone(self) -> None:
        broken = copy.deepcopy(self.manifest)
        broken["samples"][0]["somatic"]["backbone_role"] = "clairs_paired_filter_pass"
        self.assert_rejected(broken, 2, "E_SCHEMA")

    def test_explicit_clairs_pass_sensitivity_preserves_longphase_tag_binding(self) -> None:
        sensitivity = copy.deepcopy(self.manifest)
        sensitivity["analysis_contract"]["task_type"] = "backbone_sensitivity"
        sensitivity["analysis_contract"]["tree_input_contract"] = validator.SENSITIVITY_TREE_INPUT_CONTRACT
        for item in sensitivity["samples"]:
            somatic = item["somatic"]
            producer_bound_sha = item["read_tags"]["subject_binding"][
                "longphase_recalibrated_pass_vcf_sha256"
            ]
            somatic["backbone_role"] = "clairs_filter_pass_sensitivity"
            somatic["tree_vcf"] = copy.deepcopy(somatic["caller_pass_baseline_vcf"])
            self.assertEqual(
                producer_bound_sha,
                somatic["longphase_recalibrated_pass_vcf"]["identity"]["sha256"],
            )
        exit_code, lock_path, failure_path, _ = self.run_validator(sensitivity)
        self.assertEqual(exit_code, 0)
        self.assertFalse(failure_path.exists())
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
        self.assertEqual(lock["analysis_contract"]["task_type"], "backbone_sensitivity")
        self.assertTrue(all(
            item["somatic"]["tree_vcf"] == item["somatic"]["caller_pass_baseline_vcf"]
            and item["somatic"]["longphase_recalibrated_pass_vcf"]
            != item["somatic"]["longphase_recalibrated_all_vcf"]
            for item in lock["samples"]
        ))

    def test_zero_dataset_is_not_vacuous_pass(self) -> None:
        broken = copy.deepcopy(self.manifest)
        broken["dataset_count"] = 0
        broken["biological_sample_count"] = 0
        broken["samples"] = []
        self.assert_rejected(broken, 2, "E_SCHEMA_EMPTY")

    def test_duplicate_dataset_id_is_rejected(self) -> None:
        broken = copy.deepcopy(self.manifest)
        broken["samples"][1]["sample"] = broken["samples"][0]["sample"]
        self.assert_rejected(broken, 2, "E_SAMPLE_DUPLICATE")

    def test_unsafe_dataset_id_is_rejected(self) -> None:
        broken = copy.deepcopy(self.manifest)
        broken["samples"][0]["sample"] = "../escape"
        self.assert_rejected(broken, 2, "E_SAMPLE_ID_UNSAFE")

    def test_unknown_field_is_rejected(self) -> None:
        broken = copy.deepcopy(self.manifest)
        broken["read_tag_sidecr"] = "typo"
        self.assert_rejected(broken, 2, "E_UNKNOWN_FIELD")

    def test_noncomprehensive_scope_is_rejected(self) -> None:
        broken = copy.deepcopy(self.manifest)
        broken["analysis_contract"]["scope"]["contigs"] = validator.AUTOSOMES[:-1]
        self.assert_rejected(broken, 2, "E_SCOPE_NOT_COMPREHENSIVE")

    def test_six_of_seven_production_summary_is_rejected(self) -> None:
        broken = copy.deepcopy(self.manifest)
        broken["production_summary"]["passed_dataset_count"] = 6
        self.assert_rejected(broken, 2, "E_PRODUCTION_SUMMARY")

    def test_subject_swap_is_rejected(self) -> None:
        broken = copy.deepcopy(self.manifest)
        broken["samples"][0]["read_tags"]["subject_binding"]["sidecar_sha256"] = "0" * 64
        self.assert_rejected(broken, 5, "E_SIDECAR_SUBJECT_MISMATCH")

    def test_bam_mutation_breaks_storage_identity(self) -> None:
        bam = Path(self.manifest["samples"][0]["alignment_payload"]["path"])
        bam.write_bytes(bam.read_bytes() + b"mutation")
        self.assert_rejected(self.manifest, 3, "E_STORAGE_IDENTITY")

    def test_incomplete_native_validation_is_rejected(self) -> None:
        broken = copy.deepcopy(self.manifest)
        item = broken["samples"][0]
        validation_path = Path(item["read_tags"]["validation"]["path"])
        self.rewrite_validation(
            broken, 0, {"pass": True, "region": "all", "checks": {"truth_flags_absent": True}}
        )
        self.assert_rejected(broken, 5, "E_SIDECAR_VALIDATION")

    def test_truth_option_in_producer_argv_is_rejected(self) -> None:
        broken = copy.deepcopy(self.manifest)
        path = Path(broken["samples"][0]["read_tags"]["producer_capture_receipt"]["path"])
        receipt = json.loads(path.read_text(encoding="utf-8"))
        receipt["command_argv"].append("--truth-vcf=/tmp/truth.vcf.gz")
        receipt["command_argv_sha256"] = validator.canonical_sha256(
            receipt["command_argv"]
        )
        self.rewrite_producer_receipt(broken, 0, receipt)
        self.assert_rejected(broken, 4, "E_TRUTH_FLAG_PRESENT")

    def test_disable_filter_is_rejected(self) -> None:
        broken = copy.deepcopy(self.manifest)
        path = Path(broken["samples"][0]["read_tags"]["producer_capture_receipt"]["path"])
        receipt = json.loads(path.read_text(encoding="utf-8"))
        receipt["command_argv"].append("--disableFilter")
        receipt["command_argv_sha256"] = validator.canonical_sha256(
            receipt["command_argv"]
        )
        self.rewrite_producer_receipt(broken, 0, receipt)
        self.assert_rejected(broken, 4, "E_PRODUCTION_FILTER_POLICY")

    def test_effective_option_origin_is_required(self) -> None:
        broken = copy.deepcopy(self.manifest)
        path = Path(broken["samples"][0]["read_tags"]["producer_capture_receipt"]["path"])
        receipt = json.loads(path.read_text(encoding="utf-8"))
        receipt["effective_options"].pop("purity_read_assignment")
        self.rewrite_producer_receipt(broken, 0, receipt)
        self.assert_rejected(broken, 4, "E_PRODUCER_OPTION_MISMATCH")

    def test_global_duplicate_scan_cannot_be_overridden_by_validation(self) -> None:
        broken = copy.deepcopy(self.manifest)
        item = broken["samples"][0]
        sidecar_path = Path(item["read_tags"]["sidecar"]["path"])
        write_sidecar(sidecar_path, duplicate=True)
        sidecar_identity = validator.full_identity(sidecar_path)
        item["read_tags"]["sidecar"]["identity"] = sidecar_identity
        item["read_tags"]["subject_binding"]["sidecar_sha256"] = sidecar_identity["sha256"]
        receipt_path = Path(item["read_tags"]["producer_capture_receipt"]["path"])
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        receipt["capture_outputs"]["sidecar"]["identity"] = sidecar_identity
        receipt["global_coordinate_counts"]["mapped_alignment_count"] = 3
        receipt["global_coordinate_counts"]["identity_unique_count"] = 3
        item["read_tags"]["subject_binding"]["mapped_alignment_count"] = 3
        item["read_tags"]["subject_binding"]["identity_unique_count"] = 3
        self.rewrite_producer_receipt(broken, 0, receipt)
        self.assert_rejected(broken, 5, "E_SIDECAR_SUBJECT_MISMATCH")

    def test_coordinate_scan_rejects_decreasing_start(self) -> None:
        sidecar = self.root / "out_of_order.read_tags.tsv.gz"
        with gzip.open(sidecar, "wt", encoding="utf-8") as handle:
            handle.write("#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n")
            handle.write("chr1\t20\t30\treadB\t0\t60\t2222222222222222\t2\t100\n")
            handle.write("chr1\t0\t10\treadA\t0\t60\t1111111111111111\t1\t100\n")
        with self.assertRaises(validator.ContractError) as caught:
            validator.inspect_coordinate_sidecar(sidecar, "HCC1395")
        self.assertEqual(caught.exception.code, "E_SIDECAR_PARSE")

    def test_redundant_identical_tags_pass_with_explicit_collapse_evidence(self) -> None:
        observed = copy.deepcopy(self.manifest)
        item = observed["samples"][0]
        sidecar_path = Path(item["read_tags"]["sidecar"]["path"])
        write_sidecar(sidecar_path, duplicate=True)
        sidecar_identity = validator.full_identity(sidecar_path)
        item["read_tags"]["sidecar"]["identity"] = sidecar_identity
        item["read_tags"]["subject_binding"]["sidecar_sha256"] = sidecar_identity["sha256"]
        item["read_tags"]["subject_binding"].update({
            "mapped_alignment_count": 3,
            "identity_unique_count": 2,
            "duplicate_count": 1,
            "conflict_count": 0,
        })
        validation_path = Path(item["read_tags"]["validation"]["path"])
        native = json.loads(validation_path.read_text(encoding="utf-8"))
        native["capture"]["rows_mapped"] = 3
        native["duplicate_exact_alignment_rows"] = 1
        self.rewrite_validation(observed, 0, native)
        receipt_path = Path(item["read_tags"]["producer_capture_receipt"]["path"])
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        receipt["capture_outputs"]["sidecar"] = item["read_tags"]["sidecar"]
        receipt["capture_outputs"]["native_validation"] = item["read_tags"]["validation"]
        receipt["global_coordinate_counts"] = {
            "mapped_alignment_count": 3,
            "identity_unique_count": 2,
            "duplicate_count": 1,
            "conflict_count": 0,
        }
        self.rewrite_producer_receipt(observed, 0, receipt)

        exit_code, lock_path, failure_path, _ = self.run_validator(observed)
        self.assertEqual(exit_code, 0)
        self.assertTrue(lock_path.is_file())
        self.assertFalse(failure_path.exists())

    def test_synthetic_receipt_source_drift_is_rejected(self) -> None:
        self.synthetic_test_source.write_text("# drifted test source\n", encoding="utf-8")
        self.assert_rejected(self.manifest, 5, "E_JOIN_METHOD_RECEIPT_HASH")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(LayeredV3ContractTest)
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    print(stream.getvalue(), end="")
    if args.receipt:
        receipt = {
            "schema_name": "intersubmod.layered_v3_adversarial_test_receipt",
            "schema_version": "1.0.0",
            "tests_run": result.testsRun,
            "failures": len(result.failures),
            "errors": len(result.errors),
            "pass": result.wasSuccessful(),
            "test_source": {"path": str(Path(__file__).resolve()), "sha256": validator.file_sha256(Path(__file__).resolve())},
            "validator_source": {"path": str(Path(validator.__file__).resolve()), "sha256": validator.file_sha256(Path(validator.__file__).resolve())},
            "preparer_source": {"path": str(Path(prepare.__file__).resolve()), "sha256": validator.file_sha256(Path(prepare.__file__).resolve())},
            "manifest_schema": {"path": str(SCHEMA_PATH.resolve()), "sha256": validator.file_sha256(SCHEMA_PATH)},
            "lock_schema": {"path": str(LOCK_SCHEMA_PATH.resolve()), "sha256": validator.file_sha256(LOCK_SCHEMA_PATH)},
            "producer_receipt_schema": {
                "path": str(PRODUCER_RECEIPT_SCHEMA_PATH.resolve()),
                "sha256": validator.file_sha256(PRODUCER_RECEIPT_SCHEMA_PATH),
            },
        }
        validator.atomic_write_json(args.receipt, receipt)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
