from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_all_ssnv_intersubmod.py"
SPEC = importlib.util.spec_from_file_location("run_all_ssnv_intersubmod_reuse", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def make_inputs(tmp_path: Path) -> tuple[dict[str, object], Path, Path, Path, Path]:
    binary = tmp_path / "inter_sub_mod"
    reference = tmp_path / "reference.fa"
    bam = tmp_path / "tumor.bam"
    bai = tmp_path / "tumor.bam.bai"
    vcf = tmp_path / "input.vcf.gz"
    for path, content in (
        (binary, b"binary"),
        (reference, b">chr1\nA\n"),
        (Path(str(reference) + ".fai"), b"chr1\t1\t6\t1\t2\n"),
        (bam, b"bam"),
        (bai, b"bai"),
        (vcf, b"vcf"),
    ):
        path.write_bytes(content)
    entry: dict[str, object] = {
        "sample": "S1",
        "raw_alignment": {"path": str(bam)},
        "raw_alignment_index": {"path": str(bai)},
        "all_ssnv_vcf": {"path": str(vcf), "sha256": MODULE.sha256(vcf)},
        "latest_read_tag_sidecar": {"path": str(tmp_path / "tags.tsv.gz")},
        "counts": {"all_ssnv": 1},
    }
    output_root = tmp_path / "outputs"
    output_dir = output_root / "S1"
    output_dir.mkdir(parents=True)
    return entry, binary, reference, output_root, output_dir


def make_receipt(
    entry: dict[str, object],
    binary: Path,
    reference: Path,
    output_dir: Path,
    validation: dict[str, object],
) -> dict[str, object]:
    raw_alignment = entry["raw_alignment"]
    raw_index = entry["raw_alignment_index"]
    vcf = entry["all_ssnv_vcf"]
    assert isinstance(raw_alignment, dict)
    assert isinstance(raw_index, dict)
    assert isinstance(vcf, dict)
    return {
        "schema_name": "intersubmod.all_ssnv_site_run",
        "schema_version": "1.0.0",
        "sample": "S1",
        "pass": True,
        "exit_code": 0,
        "input_bam": raw_alignment["path"],
        "input_bai": raw_index["path"],
        "input_vcf": vcf["path"],
        "input_vcf_sha256": vcf["sha256"],
        "binary": str(binary),
        "binary_sha256": MODULE.sha256(binary),
        "reference": str(reference),
        "output_dir": str(output_dir),
        "reused_existing_validated_output": False,
        "validation": validation,
    }


def test_reused_receipt_provenance_checks_all_required_identities(tmp_path: Path) -> None:
    entry, binary, reference, _output_root, output_dir = make_inputs(tmp_path)
    receipt_path = output_dir / "run_receipt.json"
    validation = {"expected_vcf_sites": 1, "pass": True}
    receipt = make_receipt(entry, binary, reference, output_dir, validation)
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    assert MODULE.validate_reused_receipt_provenance(
        receipt_path, entry, binary, reference, output_dir
    ) == receipt

    mutations = {
        "sample": ("sample", "OTHER"),
        "input_vcf_sha256": ("input_vcf_sha256", "wrong"),
        "binary": ("binary", str(tmp_path / "other-binary")),
        "binary_sha256": ("binary_sha256", "wrong"),
        "reference": ("reference", str(tmp_path / "other-reference.fa")),
        "output_dir": ("output_dir", str(tmp_path / "other-output")),
        "original_receipt": ("reused_existing_validated_output", True),
    }
    for failed_gate, (field, value) in mutations.items():
        tampered = dict(receipt)
        tampered[field] = value
        receipt_path.write_text(json.dumps(tampered), encoding="utf-8")
        with pytest.raises(RuntimeError, match=failed_gate):
            MODULE.validate_reused_receipt_provenance(
                receipt_path, entry, binary, reference, output_dir
            )


def test_run_sample_reuses_original_receipt_without_reissuing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entry, binary, reference, output_root, output_dir = make_inputs(tmp_path)
    validation = {
        "expected_vcf_sites": 1,
        "summary_rows": 1,
        "summary_rows_not_emitted": 0,
        "reads_files": 1,
        "bernoulli_matrix_files": 1,
        "methylation_files": 1,
        "region_failures_in_log": 0,
        "region_stratification_status": "NOT_APPLICABLE_TUMOR_ONLY",
        "run_level_files": {
            "significance_summary": True,
            "significance_statistics": True,
            "region_stratification_status_tsv": True,
        },
        "pass": True,
    }
    receipt = make_receipt(entry, binary, reference, output_dir, validation)
    receipt_path = output_dir / "run_receipt.json"
    receipt_path.write_text(json.dumps(receipt, sort_keys=True), encoding="utf-8")
    original_bytes = receipt_path.read_bytes()
    original_mtime_ns = receipt_path.stat().st_mtime_ns
    monkeypatch.setattr(MODULE, "validate_output", lambda *_args: validation)
    monkeypatch.setattr(
        MODULE.subprocess,
        "run",
        lambda *_args, **_kwargs: pytest.fail("subprocess must not run for reused output"),
    )

    observed = MODULE.run_sample(
        entry,
        str(binary),
        str(reference),
        str(output_root),
        1,
    )

    assert observed == receipt
    assert receipt_path.read_bytes() == original_bytes
    assert receipt_path.stat().st_mtime_ns == original_mtime_ns
