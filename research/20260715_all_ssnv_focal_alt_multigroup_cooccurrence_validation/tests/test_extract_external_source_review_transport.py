from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import stat
from types import ModuleType

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    TOPIC_ROOT / "audit_notes" / "extract_external_source_review_transport.py"
)


def load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "extract_external_source_review_transport", MODULE_PATH
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def module() -> ModuleType:
    return load_module()


def test_clean_json_object_is_preserved(module: ModuleType) -> None:
    raw = '{"reviewer_id":"a","evidence":{"x":1}}\n'
    payload_text, payload, prefix = module.extract_eof_json_object(raw)
    assert payload_text == raw.rstrip()
    assert payload["evidence"] == {"x": 1}
    assert prefix == ""


def test_nonbraced_prose_prefix_is_auditable(module: ModuleType) -> None:
    raw = 'Review completed against digest `abc`.\n\n{"verdict":"APPROVE"}\n'
    payload_text, payload, prefix = module.extract_eof_json_object(raw)
    assert payload_text == '{"verdict":"APPROVE"}'
    assert payload == {"verdict": "APPROVE"}
    assert prefix == "Review completed against digest `abc`.\n\n"


@pytest.mark.parametrize(
    ("raw", "message"),
    [
        ("no object", "no JSON object"),
        ('prefix } text {"x":1}', "closing brace"),
        ('{"x":1}\\ntrailing', "trailing"),
        ('{"x":1}\\n{"y":2}', "trailing"),
        ('[{"x":1}]', "trailing"),
        ('{"x":NaN}', "Non-finite"),
        ('{"x":1,"x":2}', "Duplicate JSON key"),
    ],
)
def test_ambiguous_or_nonstandard_transport_is_rejected(
    module: ModuleType, raw: str, message: str
) -> None:
    with pytest.raises(module.TransportExtractionError, match=message):
        module.extract_eof_json_object(raw)


def test_output_is_readonly_and_cannot_be_overwritten(
    module: ModuleType, tmp_path: Path
) -> None:
    output = tmp_path / "payload.json"
    record = module.write_new_readonly_payload(output, '{"x":1}')
    assert json.loads(output.read_text(encoding="utf-8")) == {"x": 1}
    assert stat.S_IMODE(output.stat().st_mode) == 0o444
    assert record["mode"] == "0o444"
    with pytest.raises(FileExistsError):
        module.write_new_readonly_payload(output, '{"x":2}')
