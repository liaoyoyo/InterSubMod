from __future__ import annotations

import importlib.util
import hashlib
import json
import os
from pathlib import Path
import stat
from types import ModuleType

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = TOPIC_ROOT / "audit_notes" / "normalize_external_source_review_v16.py"


def load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "normalize_external_source_review_v16", MODULE_PATH
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def module() -> ModuleType:
    return load_module()


@pytest.fixture()
def valid_payload(module: ModuleType) -> dict[str, object]:
    return {
        "schema_name": module.EXPECTED_SCHEMA_NAME,
        "schema_version": module.EXPECTED_SCHEMA_VERSION,
        "reviewer_label": "External Reviewer A v15",
        "reviewer_id": "45df7edf-681f-4a4c-b69f-d0ea2c3ef527",
        "model": "claude-opus-4-8",
        "verdict": "APPROVE",
        "findings_closed": True,
        "f1_status": "RESOLVED_VERIFIED",
        "f2_status": "RESOLVED_VERIFIED",
        "reviewed_source_set_sha256": module.EXPECTED_SOURCE_SET_SHA256,
        "reviewed_git_head": module.EXPECTED_GIT_HEAD,
        "review_scope": "All 23 protected producers and release-chain contracts.",
        "blocking_findings": [],
        "nonblocking_findings": [],
        "evidence": {"source_roles": 23},
    }


def test_duplicate_keys_are_rejected_at_nested_depth(module: ModuleType) -> None:
    with pytest.raises(module.ReviewValidationError, match="Duplicate JSON key"):
        module.parse_single_json_object('{"evidence":{"x":1,"x":2}}')


@pytest.mark.parametrize(
    "raw_text",
    (
        '{"evidence":NaN}',
        '{"evidence":{"score":Infinity}}',
        '{"evidence":[-Infinity]}',
    ),
)
def test_nonfinite_json_constants_are_rejected(
    module: ModuleType, raw_text: str
) -> None:
    with pytest.raises(module.ReviewValidationError, match="Non-finite"):
        module.parse_single_json_object(raw_text)


def test_trailing_prose_is_rejected(module: ModuleType) -> None:
    with pytest.raises(module.ReviewValidationError, match="trailing"):
        module.parse_single_json_object('{"x":1}\nreview complete')


def test_clean_approval_passes(
    module: ModuleType, valid_payload: dict[str, object]
) -> None:
    module.validate_review(valid_payload, reviewer_token="A", require_approve=True)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("reviewed_source_set_sha256", "0" * 64, "source-set digest"),
        ("reviewed_git_head", "0" * 40, "Git HEAD"),
        ("reviewer_id", "not-a-uuid", "valid UUID"),
        ("f1_status", "OPEN", "inconsistent"),
    ],
)
def test_wrong_release_binding_is_rejected(
    module: ModuleType,
    valid_payload: dict[str, object],
    field: str,
    value: object,
    message: str,
) -> None:
    valid_payload[field] = value
    with pytest.raises(module.ReviewValidationError, match=message):
        module.validate_review(valid_payload, reviewer_token="A", require_approve=True)


def test_unexpected_top_level_key_is_rejected(
    module: ModuleType, valid_payload: dict[str, object]
) -> None:
    valid_payload["unexpected"] = True
    with pytest.raises(module.ReviewValidationError, match="Top-level key mismatch"):
        module.validate_review(valid_payload, reviewer_token="A", require_approve=True)


def test_approval_with_blocker_is_rejected(
    module: ModuleType, valid_payload: dict[str, object]
) -> None:
    valid_payload["blocking_findings"] = [{"id": "F1"}]
    with pytest.raises(module.ReviewValidationError, match="inconsistent"):
        module.validate_review(valid_payload, reviewer_token="A", require_approve=True)


def test_request_changes_requires_blocker(
    module: ModuleType, valid_payload: dict[str, object]
) -> None:
    valid_payload["verdict"] = "REQUEST_CHANGES"
    valid_payload["findings_closed"] = False
    with pytest.raises(module.ReviewValidationError, match="blocking finding"):
        module.validate_review(valid_payload, reviewer_token="A", require_approve=False)


def test_wrong_reviewer_identity_is_rejected(
    module: ModuleType, valid_payload: dict[str, object]
) -> None:
    with pytest.raises(module.ReviewValidationError, match="Reviewer B"):
        module.validate_review(valid_payload, reviewer_token="B", require_approve=True)


def test_bound_output_success_is_canonical_readonly_and_content_verified(
    module: ModuleType,
    valid_payload: dict[str, object],
    tmp_path: Path,
) -> None:
    output = tmp_path / "review.normalized.json"
    expected = (
        json.dumps(valid_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")

    record = module.write_new_readonly_json(output, valid_payload)

    assert record == {
        "path": str(output.resolve()),
        "size_bytes": len(expected),
        "sha256": hashlib.sha256(expected).hexdigest(),
        "mode": "0o444",
    }
    assert output.read_bytes() == expected
    assert stat.S_IMODE(output.stat().st_mode) == 0o444


def test_existing_output_is_never_overwritten(
    module: ModuleType,
    valid_payload: dict[str, object],
    tmp_path: Path,
) -> None:
    output = tmp_path / "review.normalized.json"
    output.write_bytes(b"preexisting\n")
    original = output.read_bytes()

    with pytest.raises(FileExistsError):
        module.write_new_readonly_json(output, valid_payload)

    assert output.read_bytes() == original


def test_parent_path_replacement_during_publish_is_rejected(
    module: ModuleType,
    valid_payload: dict[str, object],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parent = tmp_path / "bound-parent"
    parent.mkdir()
    displaced = tmp_path / "bound-parent.displaced"
    output = parent / "review.normalized.json"
    real_fchmod = module.os.fchmod
    replaced = False

    def fchmod_then_replace_parent(fd: int, mode: int) -> None:
        nonlocal replaced
        real_fchmod(fd, mode)
        if not replaced and mode == 0o444:
            parent.rename(displaced)
            parent.mkdir()
            replaced = True

    monkeypatch.setattr(module.os, "fchmod", fchmod_then_replace_parent)
    with pytest.raises(module.ReviewValidationError, match="parent binding changed"):
        module.write_new_readonly_json(output, valid_payload)

    assert not output.exists()
    assert (displaced / output.name).is_file()


def test_output_entry_replacement_during_publish_is_rejected(
    module: ModuleType,
    valid_payload: dict[str, object],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "review.normalized.json"
    displaced = tmp_path / "review.normalized.displaced.json"
    real_fchmod = module.os.fchmod
    replaced = False

    def fchmod_then_replace_output(fd: int, mode: int) -> None:
        nonlocal replaced
        real_fchmod(fd, mode)
        if not replaced and mode == 0o444:
            output.rename(displaced)
            output.write_bytes(b"replacement\n")
            output.chmod(0o444)
            replaced = True

    monkeypatch.setattr(module.os, "fchmod", fchmod_then_replace_output)
    with pytest.raises(module.ReviewValidationError, match="output binding changed"):
        module.write_new_readonly_json(output, valid_payload)

    assert output.read_bytes() == b"replacement\n"
    assert displaced.is_file()


def test_same_inode_same_size_terminal_byte_mutation_is_rejected(
    module: ModuleType,
    valid_payload: dict[str, object],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "review.normalized.json"
    real_require_output_bound = module.require_output_bound
    calls = 0

    def require_then_mutate_same_inode(
        output_path: Path,
        **kwargs: object,
    ) -> os.stat_result:
        nonlocal calls
        observed = real_require_output_bound(output_path, **kwargs)
        calls += 1
        if calls == 2:
            descriptor = int(kwargs["output_descriptor"])
            os.pwrite(descriptor, b"x" * observed.st_size, 0)
            os.fsync(descriptor)
        return observed

    monkeypatch.setattr(
        module,
        "require_output_bound",
        require_then_mutate_same_inode,
    )
    with pytest.raises(
        module.ReviewValidationError,
        match="terminal bytes changed",
    ):
        module.write_new_readonly_json(output, valid_payload)

    assert calls == 2
    assert output.read_bytes() == b"x" * output.stat().st_size
