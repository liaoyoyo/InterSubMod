from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = TOPIC_ROOT / "scripts" / "finalize_task_b_result_release.py"
SPEC = importlib.util.spec_from_file_location("finalize_task_b_result_release", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


@pytest.mark.parametrize(
    "payload",
    [
        b'{"pass":true,"pass":false}',
        b'{"value":NaN}',
    ],
)
def test_signed_receipt_json_is_strict(payload: bytes) -> None:
    with pytest.raises(MODULE.FinalReleaseError, match="invalid JSON"):
        MODULE.parse_json_bytes(payload, label="signed receipt")


def build_report_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> dict[str, Path]:
    paths = {
        "report_markdown": tmp_path / "report.md",
        "canonical_report_artifact": tmp_path / "artifact.json",
        "report_build_receipt": tmp_path / "report_build_receipt.json",
        "portable_html": tmp_path / "report.standalone.html",
        "portable_delivery_receipt": tmp_path / "portable_report_delivery_receipt.json",
        "portable_official_verification_screenshot": tmp_path / "official.png",
        "portable_desktop_screenshot": tmp_path / "desktop.png",
        "portable_mobile_screenshot": tmp_path / "mobile.png",
        "portable_desktop_qa": tmp_path / "desktop_qa.json",
        "portable_mobile_qa": tmp_path / "mobile_qa.json",
    }
    paths["report_markdown"].write_text("# report\n", encoding="utf-8")
    write_json(paths["canonical_report_artifact"], {"schema": "test"})
    paths["portable_html"].write_text("<!doctype html><title>report</title>\n", encoding="utf-8")
    for role in (
        "portable_official_verification_screenshot",
        "portable_desktop_screenshot",
        "portable_mobile_screenshot",
    ):
        paths[role].write_bytes(b"\x89PNG\r\n\x1a\nsynthetic")
    qa = {
        "pass": True,
        "overlapCount": 0,
        "consoleErrors": [],
        "documentScrollWidth": 100,
        "documentClientWidth": 100,
        "bodyScrollWidth": 100,
        "bodyClientWidth": 100,
    }
    write_json(paths["portable_desktop_qa"], qa)
    write_json(paths["portable_mobile_qa"], qa)
    write_json(
        paths["portable_delivery_receipt"],
        {
            "schema_name": "intersubmod.portable_report_delivery_receipt",
            "schema_version": "1.0.0",
            "status": "PASS",
            "pass": True,
            "artifact": MODULE.SOURCE_AUTHORITY.artifact(
                paths["canonical_report_artifact"]
            ),
            "output": MODULE.SOURCE_AUTHORITY.artifact(paths["portable_html"]),
        },
    )
    write_json(
        paths["report_build_receipt"],
        {
            "schema_name": "intersubmod.all_ssnv_report_build_receipt",
            "schema_version": "1.2.0",
            "task_type": "B_comprehensive_validation",
            "formal_task_b_release_eligible": True,
            "outputs": {
                "report_md": MODULE.SOURCE_AUTHORITY.artifact(
                    paths["report_markdown"]
                ),
                "artifact_json": MODULE.SOURCE_AUTHORITY.artifact(
                    paths["canonical_report_artifact"]
                ),
            },
            "code": {
                "report_builder": MODULE.SOURCE_AUTHORITY.artifact(
                    TOPIC_ROOT / "scripts" / "build_all_ssnv_report_artifact.py"
                ),
                "final_report_dataset_builder": MODULE.SOURCE_AUTHORITY.artifact(
                    Path(MODULE.BUILDER.__file__).resolve()
                ),
            },
            "pass": True,
        },
    )
    for path in paths.values():
        path.chmod(0o444)
    source_authority = {"authority_id": "TEST_ONLY_UNSIGNED_AUTHORITY", "pass": True}
    monkeypatch.setattr(MODULE, "REPORT_OUTPUTS", paths)
    monkeypatch.setattr(
        MODULE,
        "validate_release_signature_artifacts",
        lambda: {"source_authority": source_authority, "pass": True},
    )
    monkeypatch.setattr(
        MODULE.SOURCE_AUTHORITY,
        "validate_release_source_authority",
        lambda *_args, **_kwargs: source_authority,
    )
    return paths


def test_report_release_revalidates_all_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    paths = build_report_fixture(tmp_path, monkeypatch)
    source_authority, dataset_verification, outputs = MODULE.validate_report_release_outputs()
    assert source_authority["authority_id"] == "TEST_ONLY_UNSIGNED_AUTHORITY"
    assert dataset_verification["pass"] is True
    assert set(outputs) == set(paths)


def test_report_release_rejects_post_delivery_html_tamper(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = build_report_fixture(tmp_path, monkeypatch)
    html = paths["portable_html"]
    html.chmod(0o644)
    html.write_text(html.read_text(encoding="utf-8") + "<!-- tamper -->\n", encoding="utf-8")
    html.chmod(0o444)
    with pytest.raises(MODULE.FinalReleaseError, match="portable HTML artifact identity drift"):
        MODULE.validate_report_release_outputs()


def test_report_release_rejects_layout_qa_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = build_report_fixture(tmp_path, monkeypatch)
    qa_path = paths["portable_mobile_qa"]
    qa_path.chmod(0o644)
    qa = json.loads(qa_path.read_text())
    qa["overlapCount"] = 1
    qa["pass"] = False
    write_json(qa_path, qa)
    qa_path.chmod(0o444)
    with pytest.raises(MODULE.FinalReleaseError, match="portable_mobile_qa layout QA"):
        MODULE.validate_report_release_outputs()


def test_report_release_commands_are_isolated_and_report_specific() -> None:
    commands = (
        MODULE.create_command(),
        MODULE.verify_command(),
        MODULE.create_report_command(),
        MODULE.verify_report_command(),
    )
    assert all(command[:4] == MODULE.canonical_python_prefix() for command in commands)
    assert "--create-report" in commands[2]
    assert "--verify-report" in commands[3]
    assert MODULE.REPORT_RELEASE_RECEIPT != MODULE.RELEASE_RECEIPT
