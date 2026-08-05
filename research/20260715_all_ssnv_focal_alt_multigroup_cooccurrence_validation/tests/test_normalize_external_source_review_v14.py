from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = TOPIC_ROOT / "audit_notes" / "normalize_external_source_review_v14.py"


def load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "normalize_external_source_review_v14", MODULE_PATH
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
        "reviewer_label": "External Reviewer A v14",
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
