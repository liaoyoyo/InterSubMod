from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import stat
from types import ModuleType

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    TOPIC_ROOT / "audit_notes" / "assemble_release_source_authority_v15.py"
)


def load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "assemble_release_source_authority_v15_test",
        MODULE_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def module() -> ModuleType:
    return load_module()


def artifact(token: str) -> dict[str, object]:
    return {
        "path": f"/tmp/{token}.json",
        "size_bytes": 1,
        "sha256": token * 64,
        "mode": "0o444",
    }


def test_reader_detects_in_place_and_path_replacement(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.txt"
    source.write_bytes(b"before\n")
    with module.BoundArtifactReader() as reader:
        reader.open(source, include_mode=True)
        source.write_bytes(b"after\n")
        with pytest.raises(module.AssemblyError, match="binding changed"):
            reader.require_paths_still_bound()

    source.write_bytes(b"original\n")
    displaced = tmp_path / "source.original.txt"
    with module.BoundArtifactReader() as reader:
        reader.open(source, include_mode=True)
        os.rename(source, displaced)
        source.write_bytes(b"replacement\n")
        with pytest.raises(module.AssemblyError, match="binding changed"):
            reader.require_paths_still_bound()


def test_reader_detects_directory_replacement(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    directory = tmp_path / "keys"
    displaced = tmp_path / "keys.original"
    directory.mkdir(mode=0o700)
    with module.BoundArtifactReader() as reader:
        reader.open_directory(directory)
        os.rename(directory, displaced)
        directory.mkdir(mode=0o700)
        with pytest.raises(module.AssemblyError, match="directory binding"):
            reader.require_paths_still_bound()


def test_exclusive_writer_is_readonly_and_no_clobber(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    output = tmp_path / "authority.json"
    data = module.encode_json({"pass": True})
    expected = module.predicted_readonly_record(output, data)
    assert module.write_exclusive_readonly(output, data) == expected
    assert stat.S_IMODE(output.stat().st_mode) == 0o444
    with pytest.raises(FileExistsError):
        module.write_exclusive_readonly(output, data)


def test_v13_has_no_preselected_head_or_source_digest(module: ModuleType) -> None:
    assert not hasattr(module, "EXPECTED_HEAD")
    assert not hasattr(module, "EXPECTED_SOURCE_SET_SHA256")
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert "release_source_authority_v12_bootstrap" not in source
    assert "with open(RELEASE_MODULE_PATH" not in source


def test_support_routes_are_exact_v13_v18(module: ModuleType) -> None:
    support = module.support_paths()
    assert support["authority_assembler"] == MODULE_PATH.resolve()
    assert support["review_normalizer"].name == (
        "normalize_external_source_review_v20.py"
    )
    assert support["external_review_runner"].name == (
        "run_external_claude_review_v22_attested.py"
    )
    assert support["attested_evidence_validator"].name == (
        "attested_release_evidence_v5.py"
    )


def test_approval_derives_signed_process_records(module: ModuleType) -> None:
    review_payload = {
        "reviewer_id": "45df7edf-681f-4a4c-b69f-d0ea2c3ef527",
        "reviewer_label": "External Reviewer A v19",
        "model": "claude-opus-4-8",
        "verdict": "APPROVE",
        "findings_closed": True,
        "review_scope": "Complete release source chain.",
        "reviewed_source_set_sha256": "a" * 64,
        "reviewed_git_head": "b" * 40,
        "reviewed_assembler": artifact("c"),
        "canonical_junit": {
            "artifact": artifact("d"),
            "counts": {
                "tests": 2,
                "failures": 0,
                "errors": 0,
                "skipped": 0,
            },
        },
    }
    verified = {
        "review_payload": review_payload,
        "attestation": artifact("e"),
        "signature": artifact("f"),
        "public_key": artifact("1"),
        "session_id": "7b543513-8f45-41af-a1df-02ec98f518d4",
        "process_records": {
            "command_record": artifact("2"),
            "request": artifact("3"),
            "envelope": artifact("4"),
            "stderr": artifact("5"),
        },
    }
    approval = module.approval_from_review(verified)
    assert approval["process_attestation"] == verified["attestation"]
    assert approval["process_attestation_signature"] == verified["signature"]
    assert approval["process_attestation_public_key"] == verified["public_key"]
    assert approval["process_records"] == verified["process_records"]


def test_execute_module_uses_bound_bytes(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    path = tmp_path / "module.py"
    path.write_text("VALUE = 2\n", encoding="utf-8")
    loaded = module.execute_module_from_bytes(
        "assembler_bound_module_test",
        path,
        b"VALUE = 1\n",
    )
    assert loaded.VALUE == 1
