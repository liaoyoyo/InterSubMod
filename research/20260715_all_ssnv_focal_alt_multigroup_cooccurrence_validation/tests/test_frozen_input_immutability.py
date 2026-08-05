from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_frozen_input_immutability.py"
SPEC = importlib.util.spec_from_file_location("audit_frozen_input_immutability", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def frozen(path: Path, include_hash: bool = True) -> dict[str, object]:
    record: dict[str, object] = {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "mtime_ns": path.stat().st_mtime_ns,
    }
    if include_hash:
        record["sha256"] = MODULE.sha256(path)
    return record


def test_hashed_artifact_passes_and_detects_content_drift(tmp_path: Path) -> None:
    path = tmp_path / "input.tsv"
    path.write_text("one\n", encoding="ascii")
    record = frozen(path)
    assert MODULE.audit_artifact("site_manifest", record)["pass"] is True
    path.write_text("two\n", encoding="ascii")
    assert MODULE.audit_artifact("site_manifest", record)["pass"] is False


def test_large_alignment_stat_policy_passes(tmp_path: Path) -> None:
    path = tmp_path / "input.bam"
    path.write_bytes(b"bam")
    result = MODULE.audit_artifact("raw_alignment", frozen(path, include_hash=False))
    assert result["hash_policy"] == "stat_only_for_large_alignment_artifact"
    assert result["pass"] is True


def test_unhashed_small_artifact_fails_closed(tmp_path: Path) -> None:
    path = tmp_path / "input.tsv"
    path.write_text("row\n", encoding="ascii")
    result = MODULE.audit_artifact("site_manifest", frozen(path, include_hash=False))
    assert result["checks"]["frozen_sha256_present"] is False
    assert result["pass"] is False
