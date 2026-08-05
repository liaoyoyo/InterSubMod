from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import ModuleType

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    TOPIC_ROOT
    / "audit_notes"
    / "extract_external_source_review_cli_envelope_v1.py"
)


def load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "extract_external_source_review_cli_envelope_v1", MODULE_PATH
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
def command_record() -> dict[str, object]:
    return {"reviewer": "A", "session_id": "session-a"}


@pytest.fixture()
def envelope() -> dict[str, object]:
    return {
        "type": "result",
        "subtype": "success",
        "is_error": False,
        "session_id": "session-a",
        "structured_output": {
            "reviewer_label": "External Reviewer A v14",
            "reviewer_id": "45df7edf-681f-4a4c-b69f-d0ea2c3ef527",
            "verdict": "APPROVE",
        },
        "permission_denials": [],
    }


def test_clean_bound_envelope_extracts_structured_object(
    module: ModuleType,
    envelope: dict[str, object],
    command_record: dict[str, object],
) -> None:
    observed = module.extract_structured_review(
        envelope,
        command_record,
        reviewer_token="A",
    )
    assert observed is envelope["structured_output"]


@pytest.mark.parametrize(
    "raw",
    (
        b'{"x":1,"x":2}',
        b'{"x":NaN}',
        b'{"x":1} trailing',
        b"[]",
    ),
)
def test_ambiguous_or_nonfinite_json_is_rejected(
    module: ModuleType,
    raw: bytes,
) -> None:
    with pytest.raises(module.EnvelopeExtractionError):
        module.parse_json_object(raw, label="test")


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("type", "assistant", "clean success"),
        ("subtype", "error", "clean success"),
        ("is_error", True, "clean success"),
        ("session_id", "other", "session"),
        ("structured_output", [], "structured_output"),
    ),
)
def test_wrong_envelope_state_is_rejected(
    module: ModuleType,
    envelope: dict[str, object],
    command_record: dict[str, object],
    field: str,
    value: object,
    message: str,
) -> None:
    envelope[field] = value
    with pytest.raises(module.EnvelopeExtractionError, match=message):
        module.extract_structured_review(
            envelope,
            command_record,
            reviewer_token="A",
        )


def test_wrong_reviewer_token_is_rejected(
    module: ModuleType,
    envelope: dict[str, object],
    command_record: dict[str, object],
) -> None:
    with pytest.raises(module.EnvelopeExtractionError, match="reviewer token"):
        module.extract_structured_review(
            envelope,
            command_record,
            reviewer_token="B",
        )


def test_output_is_exclusive_and_readonly(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    output = tmp_path / "payload.json"
    record = module.write_new_readonly_json(output, {"verdict": "APPROVE"})
    assert record["mode"] == "0o444"
    assert json.loads(output.read_text(encoding="utf-8")) == {"verdict": "APPROVE"}
    with pytest.raises(FileExistsError):
        module.write_new_readonly_json(output, {"verdict": "REQUEST_CHANGES"})
