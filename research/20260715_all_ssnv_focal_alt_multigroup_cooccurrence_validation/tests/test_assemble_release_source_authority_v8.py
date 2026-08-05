from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import stat
from types import ModuleType

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    TOPIC_ROOT / "audit_notes" / "assemble_release_source_authority_v8.py"
)


def load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "assemble_release_source_authority_v8", MODULE_PATH
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def module() -> ModuleType:
    return load_module()


def load_bound_normalizer(
    module: ModuleType,
    reader: object,
) -> ModuleType:
    _, data, record = reader.open(module.REVIEW_NORMALIZER)
    assert data is not None
    assert record["sha256"] == module.EXPECTED_REVIEW_NORMALIZER_SHA256
    return module.execute_module_from_bytes(
        "test_bound_normalizer",
        module.REVIEW_NORMALIZER,
        data,
    )


def load_bound_signer(
    module: ModuleType,
    reader: object,
) -> ModuleType:
    _, data, record = reader.open(module.SIGNER_SOURCE)
    assert data is not None
    assert record["sha256"] == module.EXPECTED_SIGNER_SOURCE_SHA256
    return module.execute_module_from_bytes(
        "test_bound_signer",
        module.SIGNER_SOURCE,
        data,
    )


def valid_review_payload(normalizer: ModuleType, token: str = "A") -> dict[str, object]:
    return {
        "schema_name": normalizer.EXPECTED_SCHEMA_NAME,
        "schema_version": normalizer.EXPECTED_SCHEMA_VERSION,
        "reviewer_label": f"External Reviewer {token} v15",
        "reviewer_id": (
            "45df7edf-681f-4a4c-b69f-d0ea2c3ef527"
            if token == "A"
            else "7b543513-8f45-41af-a1df-02ec98f518d4"
        ),
        "model": "claude-opus-4-8",
        "verdict": "APPROVE",
        "findings_closed": True,
        "f1_status": "RESOLVED_VERIFIED",
        "f2_status": "RESOLVED_VERIFIED",
        "reviewed_source_set_sha256": normalizer.EXPECTED_SOURCE_SET_SHA256,
        "reviewed_git_head": normalizer.EXPECTED_GIT_HEAD,
        "review_scope": "All 23 protected producers and release-chain contracts.",
        "blocking_findings": [],
        "nonblocking_findings": [],
        "evidence": {"source_roles": 23},
    }


def test_bound_reader_caches_one_descriptor_and_exact_bytes(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.txt"
    source.write_bytes(b"bound payload\n")
    with module.BoundArtifactReader() as reader:
        first_fd, first_data, first_record = reader.open(source)
        second_fd, second_data, second_record = reader.open(source)
        assert first_fd == second_fd
        assert first_data == second_data == b"bound payload\n"
        assert first_record == second_record
        reader.require_paths_still_bound()


def test_bound_reader_rejects_in_place_content_change(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.txt"
    source.write_bytes(b"before\n")
    with module.BoundArtifactReader() as reader:
        reader.open(source)
        source.write_bytes(b"after-content\n")
        with pytest.raises(module.AssemblyError, match="binding changed"):
            reader.require_paths_still_bound()


def test_bound_reader_rejects_path_replacement(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.txt"
    displaced = tmp_path / "source.original.txt"
    source.write_bytes(b"original\n")
    with module.BoundArtifactReader() as reader:
        reader.open(source)
        os.rename(source, displaced)
        source.write_bytes(b"replacement\n")
        with pytest.raises(module.AssemblyError, match="binding changed"):
            reader.require_paths_still_bound()


def test_directory_binding_detects_replacement(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    directory = tmp_path / "bound"
    displaced = tmp_path / "bound.original"
    directory.mkdir()
    with module.BoundArtifactReader() as reader:
        reader.open_directory(directory)
        os.rename(directory, displaced)
        directory.mkdir()
        with pytest.raises(module.AssemblyError, match="Directory binding changed"):
            reader.require_paths_still_bound()


def test_output_is_exclusive_readonly_and_matches_prediction(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    output = tmp_path / "authority.json"
    encoded = module.encode_json({"authorized": True})
    predicted = module.predicted_readonly_record(output, encoded)
    with module.BoundArtifactReader() as reader:
        _, observed_bytes, observed = reader.create_output(output, encoded)
        assert observed_bytes == encoded
        assert observed == predicted
        assert stat.S_IMODE(output.stat().st_mode) == 0o444
        reader.require_paths_still_bound()
        with pytest.raises(module.AssemblyError, match="already exists"):
            reader.create_output(output, encoded)


def test_junit_counts_require_attribute_element_agreement(
    module: ModuleType,
) -> None:
    valid = b"""\
<testsuites>
  <testsuite tests="2" failures="1" errors="0" skipped="0">
    <testcase name="pass" />
    <testcase name="fail"><failure /></testcase>
  </testsuite>
</testsuites>
"""
    invalid = valid.replace(b'tests="2"', b'tests="3"')
    assert module.parse_junit_counts(valid) == {
        "tests": 2,
        "failures": 1,
        "errors": 0,
        "skipped": 0,
    }
    with pytest.raises(module.AssemblyError, match="aggregate attributes"):
        module.parse_junit_counts(invalid)


def test_module_execution_uses_supplied_bound_bytes(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    path = tmp_path / "module.py"
    path.write_text("VALUE = 2\n", encoding="utf-8")
    loaded = module.execute_module_from_bytes(
        "bound_test_module",
        path,
        b"VALUE = 1\n",
    )
    assert loaded.VALUE == 1


def test_strict_review_uses_bound_bytes_and_accepts_clean_approval(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    review = tmp_path / "reviewer_A.json"
    with module.BoundArtifactReader() as reader:
        normalizer = load_bound_normalizer(module, reader)
        payload = valid_review_payload(normalizer)
        review.write_text(
            json.dumps(payload, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        review.chmod(0o444)
        observed_payload, record = module.strict_load_review(
            reader,
            review,
            normalizer,
            "A",
        )
        assert observed_payload == payload
        assert record["mode"] == "0o444"
        reader.require_paths_still_bound()


def test_strict_review_rejects_nonfinite_nested_value(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    review = tmp_path / "reviewer_A.json"
    with module.BoundArtifactReader() as reader:
        normalizer = load_bound_normalizer(module, reader)
        payload = valid_review_payload(normalizer)
        raw = json.dumps(payload, sort_keys=True).replace(
            '"source_roles": 23',
            '"source_roles": NaN',
        )
        review.write_text(raw + "\n", encoding="utf-8")
        review.chmod(0o444)
        with pytest.raises(module.AssemblyError, match="Strict review validation"):
            module.strict_load_review(reader, review, normalizer, "A")


def test_retired_v8_anchor_rejects_rotated_live_release_module(
    module: ModuleType,
) -> None:
    with module.BoundArtifactReader() as reader:
        _, release_bytes, release_record = reader.open(module.RELEASE_MODULE_PATH)
        assert release_bytes is not None
        assert release_record["mode"] == "0o444"
        release_module = module.execute_module_from_bytes(
            "test_bound_release_module",
            module.RELEASE_MODULE_PATH,
            release_bytes,
        )
        with pytest.raises(
            module.AssemblyError,
            match=(
                "Release-module (path binding drifted: PUBLIC_KEY_PATH|"
                "path binding drifted: AUTHORITY_PATH|"
                "constant drifted: PUBLIC_KEY_SHA256)"
            ),
        ):
            module.require_constant_bindings(release_module)


def test_source_contains_no_legacy_two_read_loader(module: ModuleType) -> None:
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert "importlib" not in source
    assert "read_text(" not in source
    assert "module.artifact" not in source


def test_bound_v4_signer_round_trips_exact_approval_identity(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    approval_path = tmp_path / "approval.json"
    encoded = module.encode_json({"approved": True})
    approval_record = module.predicted_readonly_record(approval_path, encoded)
    with module.BoundArtifactReader() as reader:
        signer = load_bound_signer(module, reader)
        instruction = signer.format_signing_instruction(approval_record)
        parsed = signer.parse_signing_instruction(
            instruction,
            expected_path=approval_path.resolve(strict=False),
        )

    assert parsed == approval_record
    assert instruction.startswith("SIGN {")
    assert json.loads(instruction.removeprefix("SIGN ")) == approval_record
