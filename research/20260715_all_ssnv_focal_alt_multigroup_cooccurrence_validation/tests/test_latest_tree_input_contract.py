from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_latest_tree_input_contract.py"
SPEC = importlib.util.spec_from_file_location("audit_latest_tree_input_contract", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_fixture(tmp_path: Path) -> Path:
    samples = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
    layered_rows = []
    focal_rows = []
    for index, sample in enumerate(samples):
        tree = tmp_path / f"{sample}.pass.vcf.gz"
        tree.write_bytes(f"{sample}-latest-pass".encode())
        focal_vcf = tmp_path / f"{sample}.autosomal.snv.vcf.gz"
        focal_vcf.write_bytes(f"{sample}-focal-subset".encode())
        raw_bam = tmp_path / f"{sample}.bam"
        sidecar = tmp_path / f"{sample}.read_tags.tsv.gz"
        raw_bam.write_bytes(b"bam")
        sidecar.write_bytes(b"tags")
        tree_record = {
            "path": str(tree),
            "identity": {"sha256": digest(tree)},
        }
        layered_rows.append(
            {
                "sample": sample,
                "alignment_payload": {"path": str(raw_bam), "embedded_hp_ps_policy": "ignore"},
                "somatic": {
                    "backbone_role": "longphase_s_recalibrated_filter_pass",
                    "tree_vcf": tree_record,
                    "longphase_recalibrated_pass_vcf": tree_record,
                },
                "read_tags": {
                    "authority": "external_coordinate_sidecar",
                    "fallback_policy": "forbidden",
                    "sidecar": {"path": str(sidecar)},
                },
            }
        )
        focal_rows.append(
            {
                "sample": sample,
                "biological_id": "HCC1395" if sample == "HCC1395_DORADO" else sample,
                "all_ssnv_vcf": {"path": str(focal_vcf), "sha256": digest(focal_vcf)},
                "raw_alignment": {"path": str(raw_bam)},
                "latest_read_tag_sidecar": {"path": str(sidecar)},
                "counts": {"all_ssnv": 469849 if index == 0 else 0},
                "checks": {"all_vcf_equals_latest_pass_ledger": True},
            }
        )
    layered = {
        "analysis_contract": {
            "tree_input_contract": "longphase_s_recalibrated_FILTER_PASS",
            "read_tag_mode": "external_sidecar",
            "require_exact_join": True,
            "embedded_tag_policy": "ignore",
        },
        "production_summary": {"all_pass": True, "passed_dataset_count": 7},
        "samples": layered_rows,
    }
    layered_path = tmp_path / "layered.json"
    layered_path.write_text(json.dumps(layered), encoding="utf-8")
    focal = {
        "pass": True,
        "scope": "synthetic",
        "totals": {"all_ssnv": 469849},
        "layered_root": {"path": str(layered_path), "sha256": digest(layered_path)},
        "samples": focal_rows,
    }
    focal_path = tmp_path / "focal.json"
    focal_path.write_text(json.dumps(focal), encoding="utf-8")
    return focal_path


def test_realpath_and_sha_contract_produces_seven_sample_receipt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    focal = build_fixture(tmp_path)
    output = tmp_path / "audit.json"
    monkeypatch.setattr(
        sys,
        "argv",
        ["audit", "--input-manifest", str(focal), "--output", str(output)],
    )
    MODULE.main()
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["pass"] is True
    assert len(payload["samples"]) == 7
    assert all(row["pass"] for row in payload["samples"])


def test_output_is_fail_closed_against_overwrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    focal = build_fixture(tmp_path)
    output = tmp_path / "audit.json"
    output.write_text("existing", encoding="ascii")
    monkeypatch.setattr(
        sys,
        "argv",
        ["audit", "--input-manifest", str(focal), "--output", str(output)],
    )
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        MODULE.main()
