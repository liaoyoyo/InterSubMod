from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "research/20260813_public_docs_p0_correction/scripts"
sys.path.insert(0, str(SCRIPTS))
RUNNER_PATH = SCRIPTS / "run_claim_validation.py"
SPEC = importlib.util.spec_from_file_location("run_claim_validation", RUNNER_PATH)
assert SPEC is not None and SPEC.loader is not None
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)

from claim_registry_contract import (  # noqa: E402
    ALLOWED_CURRENT_VERDICTS,
    BUILD_RECEIPT_SCHEMA_NAME,
    BUILD_RECEIPT_SCHEMA_VERSION,
    REGISTRY_ID,
    REGISTRY_SCHEMA_NAME,
    REGISTRY_SCHEMA_VERSION,
    atomic_write_json,
    build_receipt_payload,
    canonical_object_sha256,
    sha256_path,
)


def create_allowlisted_html(root: Path) -> list[Path]:
    paths = []
    for relative in RUNNER.EXPECTED_HTML_RELATIVE_PATHS:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("<!doctype html><title>contract</title>\n", encoding="utf-8")
        paths.append(path)
    return paths


def source_record(root: Path, path: Path) -> dict[str, str]:
    return {
        "path": str(path.relative_to(root)),
        "sha256": sha256_path(path),
        "git_state": "DIRTY",
    }


def make_browser_receipt(root: Path, receipt_path: Path, runner_path: Path) -> None:
    html_files = create_allowlisted_html(root)
    runner_path.parent.mkdir(parents=True, exist_ok=True)
    runner_path.write_text("# browser runner fixture\n", encoding="utf-8")
    html_records = [source_record(root, path) for path in html_files]
    runner_record = source_record(root, runner_path)
    cases = []
    for mode, contract in RUNNER.EXPECTED_BROWSER_MODES.items():
        for path in html_files:
            case: dict[str, Any] = {
                "path": str(path.relative_to(root)),
                "source_sha256": sha256_path(path),
                "mode": mode,
                **contract,
                "passed": True,
                "document_overflow_pass": True,
                "page_errors": [],
                "console_errors": [],
                "external_runtime_requests": [],
                "navigation_error": None,
                "print_pdf_size_bytes": None,
                "print_pdf_sha256": None,
            }
            if mode == "print":
                case["print_pdf_size_bytes"] = 2000
                case["print_pdf_sha256"] = "a" * 64
            cases.append(case)
    tested = [
        {**record, "errors": [], "passed": True}
        for record in html_records
    ]
    receipt = {
        "schema_version": 3,
        "verdict": "PASS_FOR_LOCAL_SOURCE_QA",
        "runner": {
            "path": str(runner_path.relative_to(root)),
            "sha256": sha256_path(runner_path),
            "python": "3.10.12",
        },
        "source_provenance": {
            "html": html_records,
            "standalone_svg": [],
            "runner": runner_record,
            "stable_during_run": True,
        },
        "tested_files": tested,
        "browser_cases": cases,
        "checks": {
            "browser_cases": {
                "expected": RUNNER.EXPECTED_BROWSER_CASES,
                "total": RUNNER.EXPECTED_BROWSER_CASES,
                "per_mode": {
                    mode: len(RUNNER.EXPECTED_HTML_RELATIVE_PATHS)
                    for mode in sorted(RUNNER.EXPECTED_BROWSER_MODES)
                },
                "passed": RUNNER.EXPECTED_BROWSER_CASES,
                "failed": 0,
            },
            "page_errors": 0,
            "console_errors": 0,
            "external_runtime_requests": {"attempted": 0},
            "source_provenance_stable_during_run": True,
        },
    }
    atomic_write_json(receipt_path, receipt)


def make_registry_bundle(root: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path, Path]:
    project = root / "research/20260813_public_docs_p0_correction"
    registry_path = project / "claim_remediation_registry.json"
    receipt_path = project / "claim_remediation_build_receipt.json"
    p0_path = project / "scripts/p0_claim_registry.json"
    inventory = root / "inventory.tsv"
    source_scope = root / "source_scope.tsv"
    about = root / "about.json"
    for path, value in (
        (p0_path, {"contract": "fixture"}),
        (about, {"contract": "fixture"}),
    ):
        atomic_write_json(path, value)
    inventory.write_text("claim_id\nC001\n", encoding="utf-8")
    source_scope.write_text("artifact_id\nA001\n", encoding="utf-8")
    fake_manifest = {
        "counts": {"files": 1},
        "files": [{"path": "README.md", "sha256": "b" * 64}],
    }
    monkeypatch.setattr(RUNNER, "build_public_source_manifest", lambda _root, _p0: fake_manifest)
    registry = {
        "schema_name": REGISTRY_SCHEMA_NAME,
        "schema_version": REGISTRY_SCHEMA_VERSION,
        "registry_id": REGISTRY_ID,
        "output_path": str(registry_path.relative_to(root)),
        "source_inventory": str(inventory.relative_to(root)),
        "source_inventory_sha256": sha256_path(inventory),
        "source_scope": str(source_scope.relative_to(root)),
        "source_scope_sha256": sha256_path(source_scope),
        "p0_guard_registry": str(p0_path.relative_to(root)),
        "p0_guard_registry_sha256": sha256_path(p0_path),
        "p0_guard_validation_summary": {"verdict": "PASS"},
        "github_about_receipt": str(about.relative_to(root)),
        "github_about_receipt_sha256": sha256_path(about),
        "github_about_snapshot_manifest": {"counts": {"files": 3}, "files": []},
        "github_about_snapshot_manifest_sha256": "c" * 64,
        "public_source_manifest": fake_manifest,
        "public_source_manifest_sha256": canonical_object_sha256(fake_manifest),
        "allowed_current_verdicts": ALLOWED_CURRENT_VERDICTS,
        "counts": {"claims": 158},
        "gates": {"RELEASE_READY": {"status": "BLOCKED"}},
    }
    atomic_write_json(registry_path, registry)
    registry_bytes = registry_path.read_bytes()
    receipt = build_receipt_payload(registry, RUNNER.sha256_path(registry_path))
    assert receipt["schema_name"] == BUILD_RECEIPT_SCHEMA_NAME
    assert receipt["schema_version"] == BUILD_RECEIPT_SCHEMA_VERSION
    atomic_write_json(receipt_path, receipt)
    assert RUNNER.validate_build_receipt(registry, registry_bytes, receipt) == []
    return registry_path, receipt_path, p0_path


def test_python_310_is_fail_closed() -> None:
    assert RUNNER.require_supported_python((3, 10, 0)) == "3.10.0"
    with pytest.raises(RuntimeError, match=r"Python >= 3\.10.*running 3\.9\.19"):
        RUNNER.require_supported_python((3, 9, 19))


def test_html_allowlist_rejects_missing_and_extra(tmp_path: Path) -> None:
    paths = create_allowlisted_html(tmp_path)
    assert RUNNER.resolve_html_allowlist(tmp_path) == paths
    extra = tmp_path / "docs/explain/spoof.standalone.html"
    extra.write_text("<!doctype html>\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="unexpected=.*spoof"):
        RUNNER.resolve_html_allowlist(tmp_path)
    extra.unlink()
    paths[0].unlink()
    with pytest.raises(RuntimeError, match="missing=.*01_background"):
        RUNNER.resolve_html_allowlist(tmp_path)


def test_claim_build_receipt_binds_current_source_manifest_and_hashes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    registry, receipt, p0 = make_registry_bundle(tmp_path, monkeypatch)
    result = RUNNER.validate_claim_registry_bundle(tmp_path, registry, receipt, p0)
    assert result["verdict"] == "PASS", result["errors"]

    (tmp_path / "source_scope.tsv").write_text("artifact_id\nDRIFT\n", encoding="utf-8")
    result = RUNNER.validate_claim_registry_bundle(tmp_path, registry, receipt, p0)
    assert result["verdict"] == "FAIL"
    assert any("source_scope_sha256" in error for error in result["errors"])


def test_browser_receipt_binds_exact_84_cases_and_current_sources(tmp_path: Path) -> None:
    receipt = tmp_path / "browser.json"
    runner = tmp_path / "research/20260813_public_docs_p0_correction/scripts/run_html_browser_qa.py"
    make_browser_receipt(tmp_path, receipt, runner)
    result = RUNNER.validate_browser_qa_receipt(tmp_path, receipt, runner)
    assert result["verdict"] == "PASS", result["errors"]
    assert result["browser_cases"] == RUNNER.EXPECTED_BROWSER_CASES == 84

    first_html = tmp_path / RUNNER.EXPECTED_HTML_RELATIVE_PATHS[0]
    first_html.write_text("<!doctype html><title>drift</title>\n", encoding="utf-8")
    result = RUNNER.validate_browser_qa_receipt(tmp_path, receipt, runner)
    assert result["verdict"] == "FAIL"
    assert any("hash is stale" in error for error in result["errors"])


def test_atomic_receipt_write_is_canonical_and_leaves_no_temp_file(tmp_path: Path) -> None:
    output = tmp_path / "validation_receipt.json"
    RUNNER.atomic_write_json(output, {"verdict": "PASS", "value": 1})
    assert json.loads(output.read_text(encoding="utf-8")) == {"verdict": "PASS", "value": 1}
    assert output.read_bytes().endswith(b"\n")
    assert list(tmp_path.glob(".validation_receipt.json.*.tmp")) == []
