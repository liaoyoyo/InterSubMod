#!/usr/bin/env python3
"""Tiny tests for post-run LongPhase production receipt normalization."""

from __future__ import annotations

import contextlib
import datetime as dt
import gzip
import io
import json
import os
import shlex
import tempfile
import unittest
from pathlib import Path
from typing import Any

import build_longphase_production_capture_receipt_v3 as builder
import validate_layered_v3_inputs as contract


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_hashes(path: Path, artifacts: list[Path]) -> None:
    path.write_text(
        "".join(f"{contract.file_sha256(item)}  {item}\n" for item in artifacts), encoding="utf-8"
    )


def mtime_text(path: Path) -> str:
    # Match the production runner's GNU `stat -c` (logical path, no `-L`).
    stat = path.lstat()
    whole = dt.datetime.fromtimestamp(stat.st_mtime, dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    return f"{whole}.{stat.st_mtime_ns % 1_000_000_000:09d} +0000"


def write_sidecars(plain: Path, compressed: Path) -> None:
    text = (
        "#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n"
        "chr1\t0\t10\treadA\t0\t60\t1111111111111111\t1\t100\n"
        "chr1\t20\t30\treadB\t16\t50\t2222222222222222\t2\t100\n"
    )
    plain.write_text(text, encoding="utf-8")
    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(text)


def make_fixture(root: Path, *, complete: bool = True) -> dict[str, Any]:
    sample = "HCC1395"
    run_root = root / "run"
    sample_dir = run_root / "samples" / sample
    inputs = root / "inputs"
    code = root / "code"
    sample_dir.mkdir(parents=True)
    inputs.mkdir()
    code.mkdir()
    germline = inputs / "germline.vcf.gz"
    normal = inputs / "normal.bam"
    tumor = inputs / "tumor.bam"
    caller_pass = inputs / "caller.pass.vcf.gz"
    caller_raw = inputs / "caller.raw.vcf.gz"
    reference = inputs / "reference.fasta"
    for path, payload in (
        (germline, b"germline"), (normal, b"normal-bam"), (tumor, b"tumor-bam"),
        (caller_pass, b"caller-pass"), (caller_raw, b"caller-raw"), (reference, b">chr1\nACGT\n"),
    ):
        path.write_bytes(payload)
    for vcf in (germline, caller_pass, caller_raw):
        Path(f"{vcf}.csi").write_bytes(b"tiny-csi")
    for bam in (normal, tumor):
        Path(f"{bam}.bai").write_bytes(b"tiny-bai")
    Path(f"{reference}.fai").write_bytes(b"chr1\t4\t6\t4\t5\n")
    reference.with_suffix(".dict").write_bytes(b"@SQ\tSN:chr1\tLN:4\n")
    binary = root / "longphase-s"
    binary.write_text("#!/bin/sh\necho 'Version: 1.0.0'\n", encoding="utf-8")
    os.chmod(binary, 0o755)
    sample_meta = {
        "sample": sample,
        "biological_id": "HCC1395",
        "germline_phased_vcf": str(germline),
        "normal_bam": str(normal),
        "tumor_bam": str(tumor),
        "caller_pass_vcf": str(caller_pass),
        "caller_raw_vcf": str(caller_raw),
        "reference": str(reference),
    }
    manifest = run_root / "input_manifest.json"
    write_json(manifest, {
        "schema_version": "1.0", "dataset_count": 7, "biological_sample_count": 6, "samples": [sample_meta]
    })
    write_json(run_root / "params.json", {
        "run_root": str(run_root), "threads": 12, "parallel_samples": 1, "truth_flags": False,
        "mapq": 20, "tag_supplementary": True, "output_mode": "read_tag_sidecar",
    })
    producer_sources = []
    for name in (
        "run_longphase_production_sidecars.sh",
        "capture_longphase_tagged_bam_sidecar.py",
        "validate_streamed_longphase_sidecar.py",
    ):
        source = code / name
        source.write_text(f"# {name}\n", encoding="utf-8")
        producer_sources.append(source)
    write_hashes(run_root / "code.sha256", producer_sources)
    command = [
        str(binary), "somatic_haplotag", "-s", str(germline), "-b", str(normal),
        "--tumor-snv-file", str(caller_pass), "--tumor-bam-file", str(tumor), "-r", str(reference),
        "-t", "12", "--tagSupplementary", "-q", "20", "--output-somatic-vcf",
        "-o", str(sample_dir / f"{sample}_production"),
    ]
    (sample_dir / "command.sh.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")
    inventory_paths = {
        "germline": germline, "normal_bam": normal, "tumor_bam": tumor,
        "reference": reference, "caller_pass": caller_pass,
    }
    with (sample_dir / "input_files.tsv").open("w", encoding="utf-8") as handle:
        handle.write("role\tpath\tsize_bytes\tmtime\n")
        for role, path in inventory_paths.items():
            handle.write(f"{role}\t{path}\t{path.lstat().st_size}\t{mtime_text(path)}\n")
    write_hashes(sample_dir / "input.sha256", [
        caller_pass, Path(f"{caller_pass}.csi"), germline, Path(f"{germline}.csi"),
        Path(f"{normal}.bai"), Path(f"{tumor}.bai"),
    ])
    if complete:
        plain = sample_dir / f"{sample}.read_tags.tsv"
        sidecar = sample_dir / f"{sample}.read_tags.tsv.gz"
        write_sidecars(plain, sidecar)
        Path(f"{sidecar}.tbi").write_bytes(b"tiny-tabix")
        recal_all = sample_dir / f"{sample}.longphase_s.recalibrated.all.vcf.gz"
        recal_pass = sample_dir / f"{sample}.longphase_s.recalibrated.pass.vcf.gz"
        recal_all.write_bytes(b"recal-all")
        recal_pass.write_bytes(b"recal-pass")
        Path(f"{recal_all}.csi").write_bytes(b"tiny-csi")
        Path(f"{recal_pass}.csi").write_bytes(b"tiny-csi")
        checks = {name: True for name in builder.NATIVE_CHECKS}
        native = sample_dir / "sidecar_validation.json"
        write_json(native, {
            "schema_version": "1.0", "region": "all", "capture": {"rows_mapped": 2},
            "duplicate_exact_alignment_rows": 0, "duplicate_exact_alignment_conflicts": 0,
            "record_key_missing": 0, "record_key_extra": 0, "checks": checks, "pass": True,
        })
        write_hashes(sample_dir / "output.sha256", [plain, sidecar, recal_all, recal_pass, native])
    return {
        "sample": sample,
        "sample_dir": sample_dir,
        "run_root": run_root,
        "manifest": manifest,
        "binary": binary,
        "tumor": tumor,
        "recal_all": sample_dir / f"{sample}.longphase_s.recalibrated.all.vcf.gz",
    }


class ProductionCaptureReceiptBuilderTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="producer-receipt-v3-")
        self.root = Path(self.temporary.name)
        self.fixture = make_fixture(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def build(self, fixture: dict[str, Any] | None = None) -> dict[str, Any]:
        item = fixture or self.fixture
        return builder.build_receipt(
            item["sample_dir"], item["manifest"], item["sample"], item["run_root"], item["binary"],
            ["normalizer", "--sample", item["sample"]],
        )

    def assert_code(self, expected: str, fixture: dict[str, Any] | None = None) -> None:
        with self.assertRaises(contract.ContractError) as caught:
            self.build(fixture)
        self.assertEqual(caught.exception.code, expected)

    def test_valid_receipt_binds_full_producer_contract(self) -> None:
        receipt = self.build()
        self.assertEqual(set(receipt["producer_inputs"]), {
            "germline_phased_vcf", "normal_bam", "tumor_bam", "caller_pass_vcf", "reference"
        })
        self.assertEqual(receipt["global_coordinate_counts"]["mapped_alignment_count"], 2)
        self.assertFalse(receipt["producer_inputs"]["tumor_bam"]["storage_identity_v1"]["is_full_content_hash"])
        contract.apply_json_schema(receipt, builder.RECEIPT_SCHEMA)

    def test_cli_writes_atomic_receipt_and_never_overwrites(self) -> None:
        output = self.root / "receipt.json"
        failure = self.root / "failure.json"
        args = [
            "--sample-dir", str(self.fixture["sample_dir"]), "--production-manifest", str(self.fixture["manifest"]),
            "--sample", self.fixture["sample"], "--run-root", str(self.fixture["run_root"]),
            "--longphase-binary", str(self.fixture["binary"]), "--output", str(output),
            "--failure-report", str(failure),
        ]
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(builder.main(args), 0)
            self.assertEqual(builder.main(args), 7)
        self.assertTrue(output.is_file())
        self.assertFalse(failure.exists())

    def test_truth_token_is_rejected(self) -> None:
        path = self.fixture["sample_dir"] / "command.sh.txt"
        argv = shlex.split(path.read_text(encoding="utf-8"))
        argv.extend(["--truth-vcf", "/tmp/truth.vcf.gz"])
        path.write_text(shlex.join(argv) + "\n", encoding="utf-8")
        self.assert_code("E_TRUTH_FLAG_PRESENT")

    def test_manifest_input_swap_is_rejected(self) -> None:
        path = self.fixture["sample_dir"] / "command.sh.txt"
        argv = shlex.split(path.read_text(encoding="utf-8"))
        index = argv.index("--tumor-bam-file") + 1
        swapped = self.root / "swapped.bam"
        swapped.write_bytes(b"swapped")
        argv[index] = str(swapped)
        path.write_text(shlex.join(argv) + "\n", encoding="utf-8")
        self.assert_code("E_INPUT_BINDING_MISMATCH")

    def test_native_validation_missing_core_check_is_rejected(self) -> None:
        path = self.fixture["sample_dir"] / "sidecar_validation.json"
        native = json.loads(path.read_text(encoding="utf-8"))
        native["checks"].pop("truth_flags_absent")
        write_json(path, native)
        write_hashes(self.fixture["sample_dir"] / "output.sha256", [
            self.fixture["sample_dir"] / f"{self.fixture['sample']}.read_tags.tsv",
            self.fixture["sample_dir"] / f"{self.fixture['sample']}.read_tags.tsv.gz",
            self.fixture["sample_dir"] / f"{self.fixture['sample']}.longphase_s.recalibrated.all.vcf.gz",
            self.fixture["sample_dir"] / f"{self.fixture['sample']}.longphase_s.recalibrated.pass.vcf.gz",
            path,
        ])
        self.assert_code("E_SIDECAR_VALIDATION")

    def test_output_hash_mutation_is_rejected(self) -> None:
        self.fixture["recal_all"].write_bytes(b"mutated-after-output-sha")
        self.assert_code("E_HASH_MISMATCH")

    def test_symlink_bam_inventory_uses_logical_metadata_and_binds_target(self) -> None:
        manifest = json.loads(self.fixture["manifest"].read_text(encoding="utf-8"))
        meta = manifest["samples"][0]
        for key in ("normal_bam", "tumor_bam"):
            logical = Path(meta[key])
            payload = logical.read_bytes()
            target = logical.with_name(f"{logical.stem}.target{logical.suffix}")
            logical.unlink()
            target.write_bytes(payload + b"-target-payload")
            logical.symlink_to(target.name)

        inventory_paths = {
            "germline": Path(meta["germline_phased_vcf"]),
            "normal_bam": Path(meta["normal_bam"]),
            "tumor_bam": Path(meta["tumor_bam"]),
            "reference": Path(meta["reference"]),
            "caller_pass": Path(meta["caller_pass_vcf"]),
        }
        inventory = self.fixture["sample_dir"] / "input_files.tsv"
        with inventory.open("w", encoding="utf-8") as handle:
            handle.write("role\tpath\tsize_bytes\tmtime\n")
            for role, logical in inventory_paths.items():
                handle.write(
                    f"{role}\t{logical}\t{logical.lstat().st_size}\t{mtime_text(logical)}\n"
                )

        receipt = self.build()
        for role in ("normal_bam", "tumor_bam"):
            identity = receipt["producer_inputs"][role]["storage_identity_v1"]
            self.assertTrue(identity["logical_is_symlink"])
            self.assertNotEqual(identity["logical_size_bytes"], identity["size_bytes"])
            self.assertEqual(Path(identity["realpath"]), inventory_paths[role].resolve())

    def test_incomplete_active_sample_fails_before_large_identity_reads(self) -> None:
        incomplete = make_fixture(self.root / "incomplete", complete=False)
        self.assert_code("E_REQUIRED_ARTIFACT", incomplete)


if __name__ == "__main__":
    unittest.main(verbosity=2)
