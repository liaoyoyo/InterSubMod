#!/usr/bin/env python3
"""Tiny fail-closed tests for the raw-all LongPhase-S capture receipt v2."""

from __future__ import annotations

import datetime as dt
import gzip
import json
import os
import shlex
import tempfile
import unittest
from pathlib import Path
from typing import Any

import pysam

import audit_longphase_filter_transitions as filter_audit
import build_longphase_raw_all_capture_receipt_v2 as builder
import validate_layered_v3_inputs as contract


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_hashes(path: Path, artifacts: list[Path]) -> None:
    path.write_text(
        "".join(f"{contract.file_sha256(item)}  {item}\n" for item in artifacts), encoding="utf-8"
    )


def mtime_text(path: Path) -> str:
    observed = path.stat()
    whole = dt.datetime.fromtimestamp(observed.st_mtime, dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    return f"{whole}.{observed.st_mtime_ns % 1_000_000_000:09d} +0000"


def vcf_text(filters: tuple[str, ...]) -> str:
    records = (
        (10, "A", "C", filters[0], "0/1:3.5"),
        (20, "G", "T", filters[1], "1/1:20"),
    )
    body = "".join(
        f"chr1\t{pos}\t.\t{ref}\t{alt}\t10\t{label}\t.\tGT:GQ\t{sample}\n"
        for pos, ref, alt, label, sample in records
    )
    return (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        "##FILTER=<ID=LowQual,Description=low>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=Quality>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
        f"{body}"
    )


def write_bgzf_vcf(path: Path, text: str) -> None:
    with pysam.BGZFile(str(path), "wb") as handle:
        handle.write(text.encode("utf-8"))
    pysam.tabix_index(str(path), preset="vcf", force=True, csi=True)


def write_sidecars(plain: Path, compressed: Path, *, redundant: bool = False) -> None:
    first = "chr1\t0\t10\treadA\t0\t60\t1111111111111111\t1-1\t100\n"
    text = "#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n" + first
    if redundant:
        text += first
    text += "chr1\t20\t30\treadB\t16\t50\t2222222222222222\t2\t100\n"
    plain.write_text(text, encoding="utf-8")
    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(text)


def make_fixture(root: Path) -> dict[str, Any]:
    sample = "HCC1395"
    run_root = root / "run"
    sample_dir = run_root / "samples" / sample
    source_bundle = run_root / "source_bundle"
    inputs = root / "inputs"
    sample_dir.mkdir(parents=True)
    source_bundle.mkdir()
    inputs.mkdir()

    raw = inputs / "caller.raw.vcf.gz"
    baseline = inputs / "caller.pass.vcf.gz"
    germline = inputs / "germline.vcf.gz"
    normalized = sample_dir / f"{sample}.clairs.raw_all.normalized.vcf.gz"
    native_recal = sample_dir / f"{sample}_production_sc.vcf"
    recal_all = sample_dir / f"{sample}.longphase_s.recalibrated.all.vcf.gz"
    recal_pass = sample_dir / f"{sample}.longphase_s.recalibrated.pass.vcf.gz"
    write_bgzf_vcf(raw, vcf_text(("LowQual", "PASS")))
    write_bgzf_vcf(normalized, vcf_text(("LowQual", "PASS")))
    native_recal.write_text(vcf_text(("PASS", "LowQual")), encoding="utf-8")
    write_bgzf_vcf(recal_all, vcf_text(("PASS", "LowQual")))
    write_bgzf_vcf(recal_pass, vcf_text(("PASS", "LowQual")).replace(
        "chr1\t20\t.\tG\tT\t10\tLowQual\t.\tGT:GQ\t1/1:20\n", ""
    ))
    write_bgzf_vcf(baseline, vcf_text(("PASS", "PASS")))
    write_bgzf_vcf(germline, vcf_text(("PASS", "PASS")))

    normal = inputs / "normal.bam"
    tumor = inputs / "tumor.bam"
    for path, payload in ((normal, b"normal-bam"), (tumor, b"tumor-bam")):
        path.write_bytes(payload)
        Path(f"{path}.bai").write_bytes(b"tiny-bai")
    reference = inputs / "reference.fasta"
    reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text("chr1\t100\t6\t100\t101\n", encoding="utf-8")

    binary = source_bundle / "longphase-s"
    binary.write_text("#!/bin/sh\necho 'Version: 1.0.0'\n", encoding="utf-8")
    os.chmod(binary, 0o755)
    patch_file = source_bundle / "zero_read.patch"
    patch_file.write_text("diff --git a/x b/x\n", encoding="utf-8")
    patch_receipt = source_bundle / "patch_receipt.json"
    write_json(patch_receipt, {
        "status": "APPROVED_FOR_FAIL_CLOSED_7_DATASET_VALIDATION",
        "patched_binary_sha256": contract.file_sha256(binary),
        "patch": f"/original/{patch_file.name}",
        "approval": {"scope": "FAIL_CLOSED_7_DATASET_VALIDATION_ONLY"},
    })
    producer_sources = []
    for name in (
        "run_longphase_raw_all_production_sidecars.sh",
        "capture_longphase_tagged_bam_sidecar.py",
        "validate_streamed_longphase_sidecar.py",
        "audit_longphase_filter_transitions.py",
    ):
        source = source_bundle / name
        source.write_text(f"# {name}\n", encoding="utf-8")
        producer_sources.append(source)

    sample_meta = {
        "sample": sample,
        "biological_id": "HCC1395",
        "germline_phased_vcf": str(germline),
        "normal_bam": str(normal),
        "tumor_bam": str(tumor),
        "caller_raw_vcf": str(raw),
        "caller_raw_vcf_index": str(Path(f"{raw}.csi")),
        "caller_pass_vcf": str(baseline),
        "caller_pass_vcf_index": str(Path(f"{baseline}.csi")),
        "reference": str(reference),
        "longphase_input_contract": "normalized_ClairS_raw_all",
        "caller_raw_scan": {"records": 2},
    }
    manifest = run_root / "input_manifest.json"
    write_json(manifest, {
        "schema_version": "2.0",
        "dataset_count": 7,
        "biological_sample_count": 6,
        "longphase_binary": {
            "patch_receipt": f"/original/{patch_receipt.name}",
            "patch_receipt_sha256": contract.file_sha256(patch_receipt),
        },
        "samples": [sample_meta],
    })
    write_json(run_root / "params.json", {
        "run_root": str(run_root),
        "threads": 12,
        "parallel_samples": 1,
        "truth_flags": False,
        "mapq": 20,
        "tag_supplementary": True,
        "output_mode": "read_tag_sidecar",
        "longphase_input_contract": "normalized_ClairS_raw_all",
        "filter_contract": "bidirectional_recalibration",
    })

    command = [
        str(binary), "somatic_haplotag", "-s", str(germline), "-b", str(normal),
        "--tumor-snv-file", str(normalized), "--tumor-bam-file", str(tumor), "-r", str(reference),
        "-t", "12", "--tagSupplementary", "-q", "20", "--output-somatic-vcf",
        "-o", str(sample_dir / f"{sample}_production"),
    ]
    (sample_dir / "command.sh.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")
    inventory_paths = {
        "germline": germline,
        "normal_bam": normal,
        "tumor_bam": tumor,
        "reference": reference,
        "caller_raw": raw,
        "longphase_input": normalized,
        "caller_pass_baseline": baseline,
    }
    with (sample_dir / "input_files.tsv").open("w", encoding="utf-8") as handle:
        handle.write("role\tpath\tsize_bytes\tmtime\n")
        for role, path in inventory_paths.items():
            handle.write(f"{role}\t{path}\t{path.stat().st_size}\t{mtime_text(path)}\n")
    write_hashes(sample_dir / "input.sha256", [
        raw, Path(f"{raw}.csi"), normalized, Path(f"{normalized}.csi"),
        baseline, Path(f"{baseline}.csi"), germline, Path(f"{germline}.csi"),
        Path(f"{normal}.bai"), Path(f"{tumor}.bai"), binary,
    ])

    plain = sample_dir / f"{sample}.read_tags.tsv"
    sidecar = sample_dir / f"{sample}.read_tags.tsv.gz"
    write_sidecars(plain, sidecar)
    Path(f"{sidecar}.tbi").write_bytes(b"tiny-tabix")
    stream_summary = sample_dir / "stream_capture_summary.json"
    write_json(stream_summary, {
        "rows_mapped": 2,
        "alignment_classes": {"total": 2, "unmapped": 0},
        "pass": True,
    })
    normalization = filter_audit.audit(raw, normalized)
    transitions = filter_audit.audit(normalized, native_recal)
    normalization_path = sample_dir / "normalization_audit.json"
    transition_path = sample_dir / "filter_transition_audit.json"
    write_json(normalization_path, normalization)
    write_json(transition_path, transitions)
    checks = {name: True for name in builder.NATIVE_CHECKS}
    native = sample_dir / "sidecar_validation.json"
    write_json(native, {
        "schema_version": "1.0",
        "region": "all",
        "capture": {"rows_mapped": 2},
        "duplicate_exact_alignment_rows": 0,
        "duplicate_exact_alignment_conflicts": 0,
        "record_key_missing": 0,
        "record_key_extra": 0,
        "checks": checks,
        "pass": True,
    })
    sample_verification_path = sample_dir / "sample_verification.json"
    write_json(sample_verification_path, {
        "sample": sample,
        "expected_raw_records": 2,
        "zero_read_warnings": 0,
        "normalization": normalization,
        "filter_transitions": transitions,
        "sidecar": {"pass": True},
        "pass": True,
    })
    output_artifacts = [
        plain, sidecar, Path(f"{sidecar}.tbi"), recal_all, Path(f"{recal_all}.csi"),
        recal_pass, Path(f"{recal_pass}.csi"), normalization_path, transition_path,
        native, stream_summary, sample_verification_path,
    ]
    output_hash_manifest = sample_dir / "output.sha256"
    write_hashes(output_hash_manifest, output_artifacts)
    write_hashes(run_root / "code.sha256", [
        manifest, binary, patch_file, patch_receipt, *producer_sources,
    ])
    os.mkfifo(sample_dir / "consumed_tagged_bam.fifo")
    write_json(sample_dir / "_SUCCESS", {"status": "SUCCESS", "sample": sample})
    write_json(run_root / "_SUCCESS", {"status": "SUCCESS"})
    return {
        "sample": sample,
        "sample_dir": sample_dir,
        "run_root": run_root,
        "manifest": manifest,
        "binary": binary,
        "native_recal": native_recal,
        "recal_all": recal_all,
        "transition_path": transition_path,
        "sample_verification_path": sample_verification_path,
        "output_hash_manifest": output_hash_manifest,
        "output_artifacts": output_artifacts,
        "plain_sidecar": plain,
        "sidecar": sidecar,
        "native_validation": native,
        "stream_summary": stream_summary,
    }


class RawAllCaptureReceiptV2Test(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="raw-all-receipt-v2-")
        self.root = Path(self.temporary.name)
        self.fixture = make_fixture(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def build(self) -> dict[str, Any]:
        item = self.fixture
        return builder.build_receipt(
            item["sample_dir"], item["manifest"], item["sample"], item["run_root"], item["binary"],
            ["normalizer", "--sample", item["sample"]],
        )

    def assert_code(self, expected: str) -> None:
        with self.assertRaises(contract.ContractError) as caught:
            self.build()
        self.assertEqual(caught.exception.code, expected)

    def test_valid_receipt_binds_raw_all_transitions_and_fifo(self) -> None:
        receipt = self.build()
        self.assertEqual(receipt["schema_version"], "2.0.0")
        self.assertEqual(receipt["longphase_input_contract"], "normalized_ClairS_raw_all")
        self.assertEqual(receipt["filter_transition_summary"]["rescued_nonpass_to_pass"], 1)
        self.assertEqual(receipt["filter_transition_summary"]["removed_pass_to_nonpass"], 1)
        self.assertFalse(receipt["bam_output_policy"]["persisted_bam"])
        self.assertEqual(receipt["global_coordinate_counts"]["mapped_alignment_count"], 2)
        self.assertNotIn("dict", receipt["producer_inputs"]["reference"])
        self.assertEqual(
            receipt["duplicate_identity_policy"],
            "collapse_redundant_rows_with_identical_HP_PS",
        )
        contract.apply_json_schema(receipt, builder.RECEIPT_SCHEMA)

    def test_redundant_identical_tag_row_is_conserved_and_accepted(self) -> None:
        write_sidecars(self.fixture["plain_sidecar"], self.fixture["sidecar"], redundant=True)
        summary = json.loads(self.fixture["stream_summary"].read_text(encoding="utf-8"))
        summary["rows_mapped"] = 3
        summary["alignment_classes"]["total"] = 3
        write_json(self.fixture["stream_summary"], summary)
        native = json.loads(self.fixture["native_validation"].read_text(encoding="utf-8"))
        native["capture"]["rows_mapped"] = 3
        native["duplicate_exact_alignment_rows"] = 1
        write_json(self.fixture["native_validation"], native)
        write_hashes(self.fixture["output_hash_manifest"], self.fixture["output_artifacts"])
        receipt = self.build()
        self.assertEqual(
            receipt["global_coordinate_counts"],
            {
                "mapped_alignment_count": 3,
                "identity_unique_count": 2,
                "duplicate_count": 1,
                "conflict_count": 0,
            },
        )

    def test_persisted_bam_is_rejected(self) -> None:
        (self.fixture["sample_dir"] / "unexpected.bam").write_bytes(b"not allowed")
        self.assert_code("E_BAM_OUTPUT_POLICY")

    def test_command_using_clairs_pass_is_rejected(self) -> None:
        command_path = self.fixture["sample_dir"] / "command.sh.txt"
        argv = shlex.split(command_path.read_text(encoding="utf-8"))
        manifest = json.loads(self.fixture["manifest"].read_text(encoding="utf-8"))
        argv[argv.index("--tumor-snv-file") + 1] = manifest["samples"][0]["caller_pass_vcf"]
        command_path.write_text(shlex.join(argv) + "\n", encoding="utf-8")
        self.assert_code("E_INPUT_BINDING_MISMATCH")

    def test_stored_transition_path_must_bind_native_producer_vcf(self) -> None:
        transition = json.loads(self.fixture["transition_path"].read_text(encoding="utf-8"))
        transition["output"]["path"] = str(self.fixture["recal_all"])
        write_json(self.fixture["transition_path"], transition)
        verification = json.loads(
            self.fixture["sample_verification_path"].read_text(encoding="utf-8")
        )
        verification["filter_transitions"] = transition
        write_json(self.fixture["sample_verification_path"], verification)
        write_hashes(self.fixture["output_hash_manifest"], self.fixture["output_artifacts"])
        self.assert_code("E_FILTER_TRANSITION_AUDIT")

    def test_recalibrated_payload_mutation_is_rejected(self) -> None:
        write_bgzf_vcf(self.fixture["recal_all"], vcf_text(("PASS", "PASS")))
        self.assert_code("E_HASH_MISMATCH")


if __name__ == "__main__":
    unittest.main(verbosity=2)
