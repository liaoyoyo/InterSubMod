from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
RUNNER_PATH = (
    ROOT
    / "research/20260813_public_docs_p0_correction/scripts/run_html_browser_qa.py"
)
SPEC = importlib.util.spec_from_file_location("run_html_browser_qa", RUNNER_PATH)
assert SPEC is not None and SPEC.loader is not None
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


def run_git(root: Path, *args: str) -> None:
    subprocess.run(["git", *args], cwd=root, check=True, capture_output=True)


def create_allowlisted_html(root: Path) -> None:
    for relative_path in RUNNER.EXPECTED_HTML_RELATIVE_PATHS:
        path = root / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("<!doctype html><title>test</title>", encoding="utf-8")


def test_resolve_html_allowlist_requires_exact_paths(tmp_path: Path) -> None:
    create_allowlisted_html(tmp_path)
    resolved = RUNNER.resolve_html_allowlist(tmp_path)
    assert len(resolved) == RUNNER.EXPECTED_HTML_COUNT == 21
    assert [path.relative_to(tmp_path) for path in resolved] == list(
        RUNNER.EXPECTED_HTML_RELATIVE_PATHS
    )

    extra = tmp_path / "docs/explain/17_spoof.standalone.html"
    extra.write_text("<!doctype html>", encoding="utf-8")
    with pytest.raises(RuntimeError, match=r"unexpected=.*17_spoof"):
        RUNNER.resolve_html_allowlist(tmp_path)


def test_resolve_html_allowlist_rejects_missing_path(tmp_path: Path) -> None:
    create_allowlisted_html(tmp_path)
    missing = tmp_path / RUNNER.EXPECTED_HTML_RELATIVE_PATHS[0]
    missing.unlink()
    with pytest.raises(RuntimeError, match=r"missing=.*01_background"):
        RUNNER.resolve_html_allowlist(tmp_path)


def test_source_state_is_derived_from_git_for_each_file(tmp_path: Path) -> None:
    run_git(tmp_path, "init", "-q")
    run_git(tmp_path, "config", "user.name", "QA Test")
    run_git(tmp_path, "config", "user.email", "qa@example.invalid")

    clean = tmp_path / "clean.html"
    dirty = tmp_path / "dirty.svg"
    untracked = tmp_path / "runner.py"
    clean.write_text("clean", encoding="utf-8")
    dirty.write_text("before", encoding="utf-8")
    run_git(tmp_path, "add", "clean.html", "dirty.svg")
    run_git(tmp_path, "commit", "-qm", "fixture")
    dirty.write_text("after", encoding="utf-8")
    untracked.write_text("runner", encoding="utf-8")

    records = [
        RUNNER.source_record(tmp_path, clean),
        RUNNER.source_record(tmp_path, dirty),
        RUNNER.source_record(tmp_path, untracked),
    ]
    assert [record["git_state"] for record in records] == [
        "CLEAN_CHECKOUT",
        "DIRTY",
        "UNTRACKED",
    ]
    assert all(set(record) == {"path", "sha256", "git_state"} for record in records)
    assert RUNNER.determine_source_state(records) == "UNTRACKED"
    assert RUNNER.determine_source_state(records[:2]) == "DIRTY"
    assert RUNNER.determine_source_state(records[:1]) == "CLEAN_CHECKOUT"


def test_source_state_cli_override_is_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", [str(RUNNER_PATH), "--source-state", "CLEAN_CHECKOUT"])
    with pytest.raises(SystemExit) as exc_info:
        RUNNER.parse_args()
    assert exc_info.value.code == 2
