from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_all_ssnv_output_completeness.py"
SPEC = importlib.util.spec_from_file_location("audit_all_ssnv_output_completeness", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_site_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "chrom", "pos", "ref", "alt"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_site_artifacts(sample_root: Path, chrom: str, pos: int) -> None:
    region = sample_root / "callset" / chrom / f"{chrom}_{pos}" / "window"
    for relative in (
        "reads/reads.tsv",
        "methylation/methylation.csv",
        "distance/BERNOULLI/matrix.csv",
    ):
        path = region / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("header\n", encoding="ascii")


def test_path_site_handles_each_artifact_depth() -> None:
    root = Path("/output/sample/sample-label/chr1/chr1_123/chr1_100_200")
    paths = [
        root / "reads/reads.tsv",
        root / "methylation/methylation.csv",
        root / "distance/BERNOULLI/matrix.csv",
    ]
    assert [MODULE.path_site(path) for path in paths] == [("chr1", 123)] * 3


def test_observed_keys_reports_duplicate_and_empty(tmp_path: Path) -> None:
    first = tmp_path / "a/chr2_50/window/reads/reads.tsv"
    second = tmp_path / "b/chr2_50/window/reads/reads.tsv"
    first.parent.mkdir(parents=True)
    second.parent.mkdir(parents=True)
    first.write_text("header\n", encoding="ascii")
    second.write_text("", encoding="ascii")
    expected = {("S1", "chr2", 50, "A", "G")}
    keys, duplicates, empty, unresolved = MODULE.observed_keys(
        [first, second], "S1", MODULE.index_expected_positions(expected)
    )
    assert keys == expected
    assert duplicates == [("S1", "chr2", 50, "A", "G")]
    assert empty == 1
    assert unresolved == []


def test_expected_keys_use_dataset_chrom_pos_ref_alt_and_reject_position_collision(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "sites.tsv.gz"
    write_site_manifest(
        manifest,
        [
            {"sample": "S1", "chrom": "chr1", "pos": 10, "ref": "a", "alt": "g"},
            {"sample": "S1", "chrom": "chr1", "pos": 20, "ref": "C", "alt": "T"},
        ],
    )
    keys = MODULE.expected_keys(manifest, "S1")
    assert keys == {
        ("S1", "chr1", 10, "A", "G"),
        ("S1", "chr1", 20, "C", "T"),
    }

    keys.add(("S1", "chr1", 10, "A", "T"))
    with pytest.raises(RuntimeError, match="do not encode REF/ALT"):
        MODULE.index_expected_positions(keys)


def test_receipt_audit_fails_closed_on_provenance_tampering(tmp_path: Path) -> None:
    sample_root = tmp_path / "outputs" / "S1"
    sample_root.mkdir(parents=True)
    binary = tmp_path / "inter_sub_mod"
    reference = tmp_path / "reference.fa"
    vcf = tmp_path / "input.vcf.gz"
    for path, content in ((binary, b"binary"), (reference, b">chr1\nA\n"), (vcf, b"vcf")):
        path.write_bytes(content)
    vcf_sha256 = MODULE.sha256(vcf)
    entry = {
        "sample": "S1",
        "all_ssnv_vcf": {"path": str(vcf), "sha256": vcf_sha256},
    }
    receipt = {
        "schema_name": "intersubmod.all_ssnv_site_run",
        "sample": "S1",
        "pass": True,
        "exit_code": 0,
        "input_vcf": str(vcf),
        "input_vcf_sha256": vcf_sha256,
        "binary": str(binary),
        "binary_sha256": MODULE.sha256(binary),
        "reference": str(reference),
        "output_dir": str(sample_root),
        "validation": {"expected_vcf_sites": 1, "pass": True},
    }
    assert MODULE.receipt_audit(receipt, entry, sample_root, 1)["pass"] is True

    mutations = {
        "sample": ("sample", "OTHER"),
        "input_vcf_sha256": ("input_vcf_sha256", "wrong"),
        "binary_sha256": ("binary_sha256", "wrong"),
        "reference_nonempty_absolute": ("reference", str(tmp_path / "missing.fa")),
        "output_dir": ("output_dir", str(tmp_path / "other-output")),
    }
    for failed_gate, (field, value) in mutations.items():
        tampered = json.loads(json.dumps(receipt))
        tampered[field] = value
        audit = MODULE.receipt_audit(tampered, entry, sample_root, 1)
        assert audit["pass"] is False
        assert audit["gates"][failed_gate] is False

    with_reference_hash = dict(receipt, reference_sha256="wrong")
    audit = MODULE.receipt_audit(with_reference_hash, entry, sample_root, 1)
    assert audit["pass"] is False
    assert audit["gates"]["reference_sha256"] is False

    vcf.write_bytes(b"drifted-vcf")
    audit = MODULE.receipt_audit(receipt, entry, sample_root, 1)
    assert audit["pass"] is False
    assert audit["gates"]["input_vcf_current_sha256"] is False


def test_audit_sample_reconciles_full_key_and_receipt_provenance(tmp_path: Path) -> None:
    output_root = tmp_path / "outputs"
    sample_root = output_root / "S1"
    sample_root.mkdir(parents=True)
    write_site_artifacts(sample_root, "chr1", 123)
    (sample_root / "inter_sub_mod.log").write_text("[INFO] complete\n", encoding="ascii")
    (sample_root / "region_stratification_status.tsv").write_text(
        "status\treason\nNOT_APPLICABLE_TUMOR_ONLY\tNORMAL_BAM_REQUIRED\n",
        encoding="ascii",
    )
    site_manifest = tmp_path / "sites.tsv.gz"
    write_site_manifest(
        site_manifest,
        [{"sample": "S1", "chrom": "chr1", "pos": 123, "ref": "A", "alt": "G"}],
    )
    binary = tmp_path / "inter_sub_mod"
    reference = tmp_path / "reference.fa"
    vcf = tmp_path / "input.vcf.gz"
    for path, content in ((binary, b"binary"), (reference, b">chr1\nA\n"), (vcf, b"vcf")):
        path.write_bytes(content)
    vcf_sha256 = MODULE.sha256(vcf)
    entry = {
        "sample": "S1",
        "site_manifest": {"path": str(site_manifest)},
        "all_ssnv_vcf": {"path": str(vcf), "sha256": vcf_sha256},
    }
    receipt = {
        "schema_name": "intersubmod.all_ssnv_site_run",
        "sample": "S1",
        "pass": True,
        "exit_code": 0,
        "input_vcf": str(vcf),
        "input_vcf_sha256": vcf_sha256,
        "binary": str(binary),
        "binary_sha256": MODULE.sha256(binary),
        "reference": str(reference),
        "output_dir": str(sample_root),
        "validation": {"expected_vcf_sites": 1, "pass": True},
    }
    (sample_root / "run_receipt.json").write_text(json.dumps(receipt), encoding="utf-8")

    result = MODULE.audit_sample(entry, output_root)
    assert result["pass"] is True
    assert result["dataset"] == "S1"
    assert result["receipt_audit"]["pass"] is True
    for role in ("reads", "methylation", "bernoulli"):
        assert result["artifacts"][role]["key_fields"] == [
            "dataset",
            "chrom",
            "pos",
            "ref",
            "alt",
        ]
