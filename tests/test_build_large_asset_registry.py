import json
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "handoff" / "build_large_asset_registry.py"


def run(*args, cwd):
    return subprocess.run(args, cwd=cwd, check=True, capture_output=True, text=True)


def test_large_asset_registry_is_conservative_and_complete(tmp_path):
    run("git", "init", "-q", cwd=tmp_path)
    run("git", "config", "user.email", "test@example.invalid", cwd=tmp_path)
    run("git", "config", "user.name", "Test", cwd=tmp_path)
    (tmp_path / "small.txt").write_text("small", encoding="utf-8")
    (tmp_path / "large.html").write_text("x" * 41, encoding="utf-8")
    run("git", "add", ".", cwd=tmp_path)
    run("git", "commit", "-qm", "fixture", cwd=tmp_path)
    output = tmp_path / "registry.json"
    run(
        sys.executable,
        str(SCRIPT),
        "--repo-root",
        str(tmp_path),
        "--output",
        str(output),
        "--threshold-bytes",
        "40",
        "--verified-at",
        "2026-08-13T20:00:00+08:00",
        cwd=tmp_path,
    )
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["summary"]["tracked_assets_over_threshold"] == 1
    assert payload["source_state"] == "COMMIT_PINNED"
    assert payload["summary"]["verdict"] == "LARGE_ASSET_MIGRATION_BLOCKED"
    row = payload["assets"][0]
    assert row["public_uri"] == "large.html"
    assert row["evidence_status"] == "IN_PROGRESS"
    assert row["finality"] == "NON_FINAL"
    assert row["license"] == "NEEDS_REVIEW"
