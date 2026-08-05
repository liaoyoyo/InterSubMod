from __future__ import annotations

import copy
from dataclasses import replace
import hashlib
import importlib.util
import json
from pathlib import Path
import stat
import sys


TOPIC_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = TOPIC_ROOT / "scripts" / "verify_big7_input_binding.py"
PILOT_MANIFEST = TOPIC_ROOT / "input_contract" / "big7_hcc1395_pilot_manifest.json"
EXTRACTOR = (
    TOPIC_ROOT.parents[1]
    / "research"
    / "20260716_read_linked_hypercube_exact_likelihood_validation"
    / "scripts"
    / "extract_lossless_read_linkage.py"
)


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(name, None)
    return module


verifier = _load_module(SCRIPT, "verify_big7_input_binding_under_test")


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _full_artifact(path: Path) -> dict:
    payload = path.read_bytes()
    return {
        "path": str(path),
        "identity": {
            "policy": "full_sha256",
            "size_bytes": len(payload),
            "sha256": _sha256_bytes(payload),
        },
    }


def _sample_chunk(path: Path, label: str, offset: int, length: int):
    payload = path.read_bytes()[offset : offset + length]
    assert len(payload) == length
    return verifier.ChunkSpec(label, offset, length, _sha256_bytes(payload))


def _storage_digest(contract) -> str:
    value = {
        "policy": "storage_identity_v1",
        "assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
        "is_full_content_hash": False,
        "requested_path": str(contract.bam_path),
        "realpath": str(contract.bam_path),
        "logical_is_symlink": contract.logical_is_symlink,
        "logical_size_bytes": contract.logical_size_bytes,
        "logical_mtime_ns": contract.logical_mtime_ns,
        "st_dev": contract.st_dev,
        "st_ino": contract.st_ino,
        "size_bytes": contract.bam_size_bytes,
        "mtime_ns": contract.mtime_ns,
        "ctime_ns": contract.ctime_ns,
        "chunk_size_bytes": contract.chunk_size_bytes,
        "chunks": [
            {
                "label": chunk.label,
                "offset": chunk.offset,
                "length": chunk.length,
                "sha256": chunk.sha256,
            }
            for chunk in contract.chunks
        ],
        "index": {
            "path": str(contract.bai_path),
            "identity": {
                "policy": "full_sha256",
                "size_bytes": contract.bai_size_bytes,
                "sha256": contract.bai_sha256,
            },
        },
    }
    return verifier.canonical_sha256(value)


def _write_fake_samtools(path: Path, *, quickcheck_exit: int = 0) -> None:
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        "    print('samtools synthetic-1.0')\n"
        "    raise SystemExit(0)\n"
        f"raise SystemExit({quickcheck_exit})\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def _fixture(tmp_path: Path, *, quickcheck_exit: int = 0):
    bam = tmp_path / "input.bam"
    bai = tmp_path / "input.bam.bai"
    tree_vcf = tmp_path / "tree.vcf.gz"
    tree_vcf_index = tmp_path / "tree.vcf.gz.csi"
    sidecar = tmp_path / "read_tags.tsv.gz"
    sidecar_index = tmp_path / "read_tags.tsv.gz.tbi"
    bam.write_bytes(bytes(range(256)) * 8)
    bai.write_bytes(b"BAI\x01" + bytes(range(127)) * 3)
    tree_vcf.write_bytes(b"synthetic-tree-vcf")
    tree_vcf_index.write_bytes(b"synthetic-tree-vcf-index")
    sidecar.write_bytes(b"synthetic-sidecar")
    sidecar_index.write_bytes(b"synthetic-sidecar-index")

    source_alignment = {
        "path": str(tmp_path / "authority.bam"),
        "index_path": str(tmp_path / "authority.bam.bai"),
        "embedded_hp_ps_policy": "ignore",
        "identity_policy": "storage_identity_v1",
        "storage_identity_v1": {"authority_placeholder": True},
    }
    source_sample = {
        "sample": "HCC1395",
        "biological_id": "HCC1395",
        "platform": "ONT",
        "alignment_payload": source_alignment,
        "somatic": {
            "tree_vcf": {
                **_full_artifact(tree_vcf),
                "index": _full_artifact(tree_vcf_index),
            }
        },
        "read_tags": {
            "sidecar": _full_artifact(sidecar),
            "index": _full_artifact(sidecar_index),
            "subject_binding": {"alignment_payload_storage_identity_sha256": "authority-only"},
        },
        "marker": "all non-alignment fields must remain frozen",
    }
    source_manifest = tmp_path / "canonical_manifest.json"
    source_document = {"samples": [source_sample]}
    source_manifest.write_text(json.dumps(source_document, indent=2) + "\n", encoding="utf-8")

    target = bam.stat()
    logical = bam.lstat()
    chunk_size = 64
    middle = (target.st_size - chunk_size) // 2
    chunks = (
        _sample_chunk(bam, "first", 0, chunk_size),
        _sample_chunk(bam, "middle", middle, chunk_size),
        _sample_chunk(bam, "last", target.st_size - chunk_size, chunk_size),
    )
    contract = verifier.BindingContract(
        sample="HCC1395",
        biological_id="HCC1395",
        source_manifest_path=source_manifest,
        source_manifest_sha256=verifier.sha256_path(source_manifest),
        source_sample_semantic_sha256=verifier.canonical_sha256(source_sample),
        bam_path=bam,
        bai_path=bai,
        logical_is_symlink=stat.S_ISLNK(logical.st_mode),
        logical_size_bytes=logical.st_size,
        logical_mtime_ns=logical.st_mtime_ns,
        st_dev=target.st_dev,
        st_ino=target.st_ino,
        bam_size_bytes=target.st_size,
        mtime_ns=target.st_mtime_ns,
        ctime_ns=target.st_ctime_ns,
        chunk_size_bytes=chunk_size,
        chunks=chunks,
        bai_size_bytes=bai.stat().st_size,
        bai_sha256=verifier.sha256_path(bai),
        storage_identity_sha256="pending",
    )
    contract = replace(contract, storage_identity_sha256=_storage_digest(contract))

    pilot_sample = copy.deepcopy(source_sample)
    pilot_sample["alignment_payload"] = {
        "path": str(bam),
        "index_path": str(bai),
        "embedded_hp_ps_policy": "ignore",
        "identity_policy": "storage_identity_v1",
        "storage_identity_v1": verifier.expected_storage_identity(contract),
    }
    pilot_document = {
        "schema_name": verifier.PILOT_SCHEMA_NAME,
        "schema_version": verifier.PILOT_SCHEMA_VERSION,
        "manifest_id": "synthetic_hcc1395_partial",
        "created_at_utc": "2026-07-22T00:00:00Z",
        "task_type": "exploratory_pilot",
        "validation_evidence_eligible": False,
        "scope": {
            "claim_status": "PARTIAL",
            "dataset": "HCC1395",
            "biological_id": "HCC1395",
            "contigs": list(verifier.AUTOSOMES),
            "purpose": "exact_ps_x_hp_k_le_12_input_binding",
            "claim_limit": (
                "HCC1395 chr1-22 input binding only; not seven-sample, production, "
                "or paper-final validation"
            ),
        },
        "source_authority": {
            "manifest_path": str(source_manifest),
            "manifest_sha256": contract.source_manifest_sha256,
            "sample": "HCC1395",
            "sample_semantic_sha256": contract.source_sample_semantic_sha256,
            "alignment_transfer_assurance": (
                "same_size_plus_fixed_first_middle_last_chunk_sha256_plus_full_bai_sha256"
            ),
            "is_full_bam_content_hash": False,
        },
        "extractor_contract": copy.deepcopy(verifier.EXTRACTOR_CONTRACT),
        "dataset_count": 1,
        "biological_sample_count": 1,
        "samples": [pilot_sample],
    }
    pilot_manifest = tmp_path / "pilot_manifest.json"
    pilot_manifest.write_text(json.dumps(pilot_document, indent=2) + "\n", encoding="utf-8")
    fake_samtools = tmp_path / "samtools"
    _write_fake_samtools(fake_samtools, quickcheck_exit=quickcheck_exit)
    return {
        "bam": bam,
        "bai": bai,
        "tree_vcf": tree_vcf,
        "tree_vcf_index": tree_vcf_index,
        "sidecar": sidecar,
        "sidecar_index": sidecar_index,
        "contract": contract,
        "pilot_document": pilot_document,
        "pilot_manifest": pilot_manifest,
        "samtools": fake_samtools,
    }


def _rewrite_manifest(fixture: dict) -> None:
    fixture["pilot_manifest"].write_text(
        json.dumps(fixture["pilot_document"], indent=2) + "\n",
        encoding="utf-8",
    )


def test_exact_binding_passes_and_records_middle_offset(tmp_path: Path) -> None:
    fixture = _fixture(tmp_path)
    receipt = verifier.verify_binding(
        fixture["pilot_manifest"],
        contract=fixture["contract"],
        samtools=str(fixture["samtools"]),
    )
    assert receipt["all_pass"] is True
    assert receipt["failure"] is None
    assert receipt["binding_contract"]["bam"]["is_full_content_hash"] is False
    assert receipt["binding_contract"]["bam"]["chunks"][1]["offset"] == 992
    assert receipt["observed"]["bai"]["sha256"] == fixture["contract"].bai_sha256
    assert receipt["observed"]["samtools_quickcheck"]["returncode"] == 0
    authority = receipt["observed"]["frozen_authority_artifacts"]
    assert set(authority) == {
        "tree_vcf",
        "tree_vcf_index",
        "read_tag_sidecar",
        "read_tag_sidecar_index",
    }
    assert all(row["verification_level"] == "fresh_full_sha256" for row in authority.values())
    assert json.loads(json.dumps(receipt, allow_nan=False))["all_pass"] is True


def test_relative_declared_input_path_fails_closed(tmp_path: Path) -> None:
    fixture = _fixture(tmp_path)
    fixture["pilot_document"]["samples"][0]["somatic"]["tree_vcf"]["path"] = "relative.vcf.gz"
    _rewrite_manifest(fixture)
    receipt = verifier.verify_binding(
        fixture["pilot_manifest"],
        contract=fixture["contract"],
        samtools=str(fixture["samtools"]),
    )
    assert receipt["all_pass"] is False
    assert receipt["failure"]["code"] == "E_PATH_ABSOLUTE"


def test_frozen_vcf_or_sidecar_record_change_fails_closed(tmp_path: Path) -> None:
    fixture = _fixture(tmp_path)
    fixture["pilot_document"]["samples"][0]["somatic"]["tree_vcf"]["identity"]["sha256"] = "0" * 64
    _rewrite_manifest(fixture)
    receipt = verifier.verify_binding(
        fixture["pilot_manifest"],
        contract=fixture["contract"],
        samtools=str(fixture["samtools"]),
    )
    assert receipt["all_pass"] is False
    assert receipt["failure"]["code"] == "E_FROZEN_FIELDS"


def test_same_size_authority_artifact_content_change_fails_closed(tmp_path: Path) -> None:
    fixture = _fixture(tmp_path)
    with fixture["sidecar"].open("r+b") as handle:
        handle.seek(3)
        original = handle.read(1)
        handle.seek(3)
        handle.write(bytes([original[0] ^ 0xFF]))
    receipt = verifier.verify_binding(
        fixture["pilot_manifest"],
        contract=fixture["contract"],
        samtools=str(fixture["samtools"]),
    )
    assert fixture["sidecar"].stat().st_size == len(b"synthetic-sidecar")
    assert receipt["all_pass"] is False
    assert receipt["failure"]["code"] == "E_AUTHORITY_SHA256"


def test_middle_bam_chunk_change_fails_closed(tmp_path: Path) -> None:
    fixture = _fixture(tmp_path)
    with fixture["bam"].open("r+b") as handle:
        handle.seek(fixture["contract"].chunks[1].offset)
        original = handle.read(1)
        handle.seek(fixture["contract"].chunks[1].offset)
        handle.write(bytes([original[0] ^ 0xFF]))
    receipt = verifier.verify_binding(
        fixture["pilot_manifest"],
        contract=fixture["contract"],
        samtools=str(fixture["samtools"]),
    )
    assert receipt["all_pass"] is False
    assert receipt["failure"]["code"] == "E_BAM_IDENTITY"


def test_full_bai_hash_change_fails_closed(tmp_path: Path) -> None:
    fixture = _fixture(tmp_path)
    with fixture["bai"].open("r+b") as handle:
        handle.seek(8)
        original = handle.read(1)
        handle.seek(8)
        handle.write(bytes([original[0] ^ 0xFF]))
    receipt = verifier.verify_binding(
        fixture["pilot_manifest"],
        contract=fixture["contract"],
        samtools=str(fixture["samtools"]),
    )
    assert receipt["all_pass"] is False
    assert receipt["failure"]["code"] == "E_BAI_IDENTITY"


def test_samtools_quickcheck_failure_fails_closed(tmp_path: Path) -> None:
    fixture = _fixture(tmp_path, quickcheck_exit=1)
    receipt = verifier.verify_binding(
        fixture["pilot_manifest"],
        contract=fixture["contract"],
        samtools=str(fixture["samtools"]),
    )
    assert receipt["all_pass"] is False
    assert receipt["failure"]["code"] == "E_SAMTOOLS_QUICKCHECK"


def test_checked_in_pilot_manifest_is_accepted_by_current_extractor_loader() -> None:
    extractor = _load_module(EXTRACTOR, "m2_extractor_loader_compatibility_test")
    document = json.loads(PILOT_MANIFEST.read_text(encoding="utf-8"))
    entry = extractor.input_entry(document, "HCC1395")
    assert len(document["samples"]) == 1
    assert entry["sample"] == "HCC1395"
    assert entry["alignment_payload"]["path"] == str(verifier.PRODUCTION_CONTRACT.bam_path)
    assert entry["alignment_payload"]["index_path"] == str(verifier.PRODUCTION_CONTRACT.bai_path)
