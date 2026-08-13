#!/usr/bin/env python3
"""Build the portable workflow governance registry from Git index blobs.

The Git index is the authority on purpose: modified worktree bytes, untracked
files, and generated ``__pycache__`` files cannot change the registry.  A new
script must be added to the index before it can appear in this inventory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any, Mapping, Sequence


SCHEMA_VERSION = "1.0.0"
DEFAULT_GENERATED_AT = "2026-08-13T00:00:00+08:00"
CLASSIFICATIONS = ("SUPPORTED", "REPRODUCIBLE_LEGACY", "ARCHIVED")
REGULAR_GIT_MODES = frozenset({"100644", "100755"})
ARCHIVED_PREFIX = "scripts/analysis/legacy/"
UNSAFE_CLEANUP_PATH = "scripts/pipeline/steps/04_cleanup.sh"
DEFAULT_OUTPUT = Path(
    "docs/handoff/20260813_完整研究資料與軟體交接_01/registries/workflow_registry.json"
)

# This allowlist is intentionally narrow.  Adding a path is a public support
# decision and must be reviewed together with the tests and handoff docs.
SUPPORTED_ALLOWLIST = frozenset(
    {
        # Portable machine/site profile interface.
        "scripts/site/bootstrap",
        "scripts/site/doctor",
        "scripts/site/site_profile.py",
        # Portable pipeline interface and its minimum reviewed dependency set.
        "scripts/pipeline/config.sh",
        "scripts/pipeline/run_benchmark.sh",
        "scripts/pipeline/steps/00_prepare_germline.sh",
        "scripts/pipeline/steps/01_longphase_s.sh",
        "scripts/pipeline/steps/02_intersubmod.sh",
        "scripts/pipeline/steps/03_filter_analysis.py",
        "scripts/analysis/haplotag_qc.py",
        "scripts/analysis/extract_label_first_metrics.py",
        "scripts/lib/verification_schema_contract.py",
        # Public handoff registry/package and tiny DEMO validators/runners.
        "scripts/handoff/build_registries.py",
        "scripts/handoff/build_workflow_registry.py",
        "scripts/handoff/validate_handoff_package.py",
        "scripts/handoff/build_tiny_public_fixture.sh",
        "scripts/handoff/run_tiny_public_e2e.sh",
        "scripts/handoff/validate_tiny_public_e2e.py",
    }
)

PIPELINE_SUPPORTED_PATHS = frozenset(
    path
    for path in SUPPORTED_ALLOWLIST
    if path.startswith("scripts/pipeline/")
    or path in {
        "scripts/analysis/haplotag_qc.py",
        "scripts/analysis/extract_label_first_metrics.py",
        "scripts/lib/verification_schema_contract.py",
    }
)
SITE_SUPPORTED_PATHS = frozenset(path for path in SUPPORTED_ALLOWLIST if path.startswith("scripts/site/"))
TINY_SUPPORTED_PATHS = frozenset(
    {
        "scripts/handoff/build_tiny_public_fixture.sh",
        "scripts/handoff/run_tiny_public_e2e.sh",
        "scripts/handoff/validate_tiny_public_e2e.py",
    }
)

# This is a deliberately lexical count.  It finds literal POSIX-looking path
# tokens while avoiding URLs and variable-relative fragments such as ${ROOT}/x.
ABSOLUTE_PATH_TOKEN_RE = re.compile(
    rb"(?<![A-Za-z0-9_$}./:-])/(?!/)(?:[A-Za-z0-9._~+@%=-]+(?:/[A-Za-z0-9._~+@%={}$:-]+)*)"
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class RegistryError(RuntimeError):
    """Raised when the Git inventory or registry contract is invalid."""


@dataclass(frozen=True)
class GitTrackedFile:
    path: str
    mode: str
    object_id: str
    content: bytes


def run_git(repo_root: Path, arguments: Sequence[str], *, binary: bool = False) -> bytes | str:
    result = subprocess.run(
        ["git", "-C", str(repo_root), *arguments],
        check=False,
        capture_output=True,
        text=not binary,
    )
    if result.returncode != 0:
        stderr = result.stderr.decode("utf-8", "replace") if binary else result.stderr
        raise RegistryError(f"git {' '.join(arguments)} failed: {stderr.strip()}")
    return result.stdout


def validate_repo_root(repo_root: Path) -> Path:
    resolved = repo_root.resolve()
    top_level = str(run_git(resolved, ["rev-parse", "--show-toplevel"])).strip()
    if Path(top_level).resolve() != resolved:
        raise RegistryError(f"--repo-root must be the Git top level: {resolved}")
    return resolved


def is_pycache_path(path: str) -> bool:
    return "__pycache__" in PurePosixPath(path).parts or path.endswith((".pyc", ".pyo"))


def git_tracked_regular_files(repo_root: Path) -> list[GitTrackedFile]:
    """Return eligible stage-0 regular scripts using index blobs, not worktree bytes."""

    repo_root = validate_repo_root(repo_root)
    raw = run_git(repo_root, ["ls-files", "--stage", "-z", "--", "scripts"], binary=True)
    assert isinstance(raw, bytes)
    rows: list[tuple[str, str, str]] = []
    non_stage_zero: list[str] = []
    for item in raw.split(b"\0"):
        if not item:
            continue
        try:
            metadata, raw_path = item.split(b"\t", 1)
            mode, object_id, stage = metadata.decode("ascii").split()
            path = raw_path.decode("utf-8")
        except (UnicodeDecodeError, ValueError) as error:
            raise RegistryError("unable to parse git ls-files --stage output") from error
        if stage != "0":
            non_stage_zero.append(path)
            continue
        if mode not in REGULAR_GIT_MODES or is_pycache_path(path):
            continue
        rows.append((path, mode, object_id))
    if non_stage_zero:
        raise RegistryError(f"unmerged scripts entries are not registrable: {sorted(set(non_stage_zero))}")

    paths = [path for path, _, _ in rows]
    duplicates = sorted(path for path, count in Counter(paths).items() if count != 1)
    if duplicates:
        raise RegistryError(f"tracked script paths are not unique: {duplicates}")

    blob_cache: dict[str, bytes] = {}
    records: list[GitTrackedFile] = []
    for path, mode, object_id in sorted(rows):
        if object_id not in blob_cache:
            content = run_git(repo_root, ["cat-file", "blob", object_id], binary=True)
            assert isinstance(content, bytes)
            blob_cache[object_id] = content
        records.append(GitTrackedFile(path=path, mode=mode, object_id=object_id, content=blob_cache[object_id]))
    return records


def shebang_for(content: bytes) -> str | None:
    first_line = content.splitlines()[0] if content else b""
    if not first_line.startswith(b"#!"):
        return None
    return first_line.decode("utf-8", "replace")


def workflow_type(path: str, shebang: str | None) -> str:
    lower_shebang = (shebang or "").lower()
    suffix = PurePosixPath(path).suffix.lower()
    if "python" in lower_shebang or suffix == ".py":
        return "PYTHON"
    if any(shell in lower_shebang for shell in ("/bash", "/sh", "/zsh")) or suffix == ".sh":
        return "SHELL"
    if "node" in lower_shebang or suffix == ".js":
        return "JAVASCRIPT"
    if suffix == ".ts":
        return "TYPESCRIPT"
    if suffix == ".md":
        return "MARKDOWN"
    if suffix in {".tsv", ".csv"}:
        return "TABULAR_DATA"
    if suffix == ".json":
        return "JSON_DATA"
    if suffix in {".yaml", ".yml"}:
        return "YAML_DATA"
    if suffix in {".service", ".timer"}:
        return "SYSTEMD_UNIT"
    if suffix == ".css":
        return "STYLESHEET"
    return "OTHER_REGULAR_FILE"


def classification_for(path: str) -> str:
    if path.startswith(ARCHIVED_PREFIX):
        return "ARCHIVED"
    if path in SUPPORTED_ALLOWLIST:
        return "SUPPORTED"
    return "REPRODUCIBLE_LEGACY"


def reason_for(path: str, classification: str) -> str:
    if classification == "ARCHIVED":
        return "Located under scripts/analysis/legacy/** and intentionally retained as an archived research workflow."
    if path == UNSAFE_CLEANUP_PATH:
        return (
            "Retained for reproducibility but excluded from portable support because it directly deletes files; "
            "an archive-first replacement is required."
        )
    if classification == "REPRODUCIBLE_LEGACY":
        return (
            "Git-tracked workflow material retained for reproducibility, but outside the bounded bip7/bip8 "
            "portable public support allowlist."
        )
    if path in SITE_SUPPORTED_PATHS:
        return "Reviewed public interface for creating, loading, or validating a machine-local portable site profile."
    if path in TINY_SUPPORTED_PATHS:
        return "Reviewed public tiny-fixture builder, runner, or validator for portable DEMO execution."
    if path == "scripts/handoff/build_workflow_registry.py":
        return "Reviewed public builder for the deterministic Git-index workflow governance registry."
    if path == "scripts/handoff/build_registries.py":
        return "Reviewed public handoff registry builder with explicit authority and storage-root inputs."
    if path == "scripts/handoff/validate_handoff_package.py":
        return "Reviewed public fail-closed validator for the bounded handoff package contract."
    if path == "scripts/pipeline/run_benchmark.sh":
        return "Reviewed public portable benchmark orchestrator, bounded to validated site-profile execution."
    return "Minimum reviewed internal dependency of the supported portable benchmark interface."


def known_limits_for(path: str, classification: str, absolute_path_count: int) -> list[str]:
    if classification == "ARCHIVED":
        limits = [
            "ARCHIVED is a governance lifecycle label; the Git-tracked bytes remain available for historical audit.",
            "No bip7/bip8 portable execution or maintenance support is promised.",
        ]
    elif classification == "REPRODUCIBLE_LEGACY":
        limits = [
            "REPRODUCIBLE_LEGACY does not mean broken; it means bip7/bip8 portable support is not promised.",
            "Dependencies, side effects, and scientific validity require workflow-specific review before reuse.",
        ]
    else:
        limits = [
            "SUPPORTED is bounded to the documented public invocation and is not a scientific-result guarantee."
        ]
        if path in PIPELINE_SUPPORTED_PATHS:
            limits.append(
                "Portable execution requires a validated site profile; built-in /big* fallback paths are outside support."
            )
        if path in SITE_SUPPORTED_PATHS:
            limits.append("Operators must supply real machine paths and checksums; example placeholder hashes fail closed.")
        if path in TINY_SUPPORTED_PATHS:
            limits.append("The tiny fixture is DEMO/PARTIAL engineering evidence and cannot validate scientific claims.")
        if path == "scripts/handoff/build_registries.py":
            limits.append("Local authority, run-manifest, and storage roots must be supplied; large payloads are not bundled.")
        if path == "scripts/handoff/build_workflow_registry.py":
            limits.append("Only stage-0 Git index blobs are inventoried; dirty and untracked worktree bytes are ignored.")
        if path == "scripts/handoff/validate_handoff_package.py":
            limits.append("Package-contract validation does not rerun the underlying scientific analyses.")

    if path == UNSAFE_CLEANUP_PATH:
        limits.extend(
            [
                "Do not invoke this stage in the public portable workflow because it performs direct deletion.",
                "Use run_benchmark.sh with --skip-cleanup until an archive-first cleanup replacement is available.",
            ]
        )
    elif path == "scripts/pipeline/run_benchmark.sh":
        limits.append(
            "Portable runs must use --skip-cleanup until stage 04 is replaced by an archive-first implementation."
        )
    if absolute_path_count:
        limits.append(
            f"Lexical scan found {absolute_path_count} literal POSIX absolute-path token(s); review their runtime role."
        )
    return limits


def record_for(entry: GitTrackedFile) -> dict[str, Any]:
    shebang = shebang_for(entry.content)
    suffix = PurePosixPath(entry.path).suffix
    classification = classification_for(entry.path)
    absolute_path_count = len(ABSOLUTE_PATH_TOKEN_RE.findall(entry.content))
    return {
        "path": entry.path,
        "type": workflow_type(entry.path, shebang),
        "classification": classification,
        "reason": reason_for(entry.path, classification),
        "shebang": shebang,
        "suffix": suffix,
        "git_mode": entry.mode,
        "size_bytes": len(entry.content),
        "sha256": hashlib.sha256(entry.content).hexdigest(),
        "absolute_path_token_count": absolute_path_count,
        "known_limits": known_limits_for(entry.path, classification, absolute_path_count),
    }


def tree_content_sha256(records: Sequence[Mapping[str, Any]]) -> str:
    """Hash sorted path/mode/size/content tuples with unambiguous NUL framing."""

    digest = hashlib.sha256()
    for row in sorted(records, key=lambda item: str(item["path"])):
        fields = (row["path"], row["git_mode"], str(row["size_bytes"]), row["sha256"])
        for field in fields:
            digest.update(str(field).encode("utf-8"))
            digest.update(b"\0")
        digest.update(b"\n")
    return digest.hexdigest()


def records_sha256(records: Sequence[Mapping[str, Any]]) -> str:
    payload = json.dumps(records, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def validate_generated_at(value: str) -> str:
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise RegistryError(f"invalid --generated-at ISO-8601 value: {value}") from error
    if parsed.tzinfo is None:
        raise RegistryError("--generated-at must include an explicit UTC offset")
    return value


def build_registry(repo_root: Path, generated_at: str = DEFAULT_GENERATED_AT) -> dict[str, Any]:
    generated_at = validate_generated_at(generated_at)
    entries = git_tracked_regular_files(repo_root)
    records = [record_for(entry) for entry in entries]
    classification_counts = Counter(row["classification"] for row in records)
    type_counts = Counter(row["type"] for row in records)
    summary = {
        "total_records": len(records),
        "by_classification": {name: classification_counts[name] for name in CLASSIFICATIONS},
        "legacy_or_archived_records": classification_counts["REPRODUCIBLE_LEGACY"]
        + classification_counts["ARCHIVED"],
        "by_type": dict(sorted(type_counts.items())),
        "absolute_path_token_total": sum(row["absolute_path_token_count"] for row in records),
    }
    registry = {
        "schema_version": SCHEMA_VERSION,
        "registry_type": "git_tracked_workflow_governance",
        "generated_at": generated_at,
        "source": {
            "authority": "Git index stage-0 blobs",
            "scope": "Git-tracked regular files under scripts/**, excluding cache bytecode",
            "dirty_worktree_used": False,
            "untracked_files_used": False,
        },
        "classification_enum": list(CLASSIFICATIONS),
        "classification_meanings": {
            "SUPPORTED": "Bounded public bip7/bip8 portable interface or minimum reviewed dependency.",
            "REPRODUCIBLE_LEGACY": (
                "Retained for reproducibility; not known-broken, but no bip7/bip8 portable support promise."
            ),
            "ARCHIVED": "Historical scripts/analysis/legacy material retained for audit, not public execution.",
        },
        "portable_safety_boundary": {
            "cleanup_stage": UNSAFE_CLEANUP_PATH,
            "cleanup_classification": "REPRODUCIBLE_LEGACY",
            "required_pipeline_option": "--skip-cleanup",
            "replacement_gate": "archive-first cleanup implementation and validation",
        },
        "summary": summary,
        "tree_content_hash": {
            "algorithm": "sha256",
            "canonicalization": "sorted(path NUL git_mode NUL size_bytes NUL sha256 NUL LF)",
            "value": tree_content_sha256(records),
        },
        "records_hash": {
            "algorithm": "sha256",
            "canonicalization": "UTF-8 canonical JSON of records; generated_at excluded",
            "value": records_sha256(records),
        },
        "records": records,
    }
    validate_registry(registry)
    return registry


def validate_registry(registry: Mapping[str, Any]) -> None:
    records = registry.get("records")
    if not isinstance(records, list):
        raise RegistryError("records must be a list")
    paths = [row.get("path") for row in records]
    if len(paths) != len(set(paths)):
        raise RegistryError("records contain duplicate paths")
    if paths != sorted(paths):
        raise RegistryError("records must be sorted by path")
    if any(is_pycache_path(str(path)) for path in paths):
        raise RegistryError("cache bytecode must not appear in the registry")
    for row in records:
        if row.get("classification") not in CLASSIFICATIONS:
            raise RegistryError(f"invalid classification: {row.get('path')}")
        if row.get("classification") == "SUPPORTED" and row.get("path") not in SUPPORTED_ALLOWLIST:
            raise RegistryError(f"SUPPORTED path is outside allowlist: {row.get('path')}")
        if not isinstance(row.get("known_limits"), list) or not row["known_limits"]:
            raise RegistryError(f"known_limits must be a non-empty list: {row.get('path')}")
        if not SHA256_RE.fullmatch(str(row.get("sha256", ""))):
            raise RegistryError(f"invalid sha256: {row.get('path')}")
    if UNSAFE_CLEANUP_PATH in paths:
        cleanup = records[paths.index(UNSAFE_CLEANUP_PATH)]
        if cleanup["classification"] != "REPRODUCIBLE_LEGACY":
            raise RegistryError("unsafe cleanup stage must remain REPRODUCIBLE_LEGACY")

    summary = registry.get("summary", {})
    counts = Counter(row["classification"] for row in records)
    expected_counts = {name: counts[name] for name in CLASSIFICATIONS}
    if summary.get("total_records") != len(records) or summary.get("by_classification") != expected_counts:
        raise RegistryError("summary classification counts do not match records")
    if summary.get("legacy_or_archived_records") != counts["REPRODUCIBLE_LEGACY"] + counts["ARCHIVED"]:
        raise RegistryError("summary legacy/archive count does not match records")
    if summary.get("legacy_or_archived_records", 0) <= 200:
        raise RegistryError("legacy/archive inventory unexpectedly contains 200 or fewer records")
    if registry.get("tree_content_hash", {}).get("value") != tree_content_sha256(records):
        raise RegistryError("tree content hash mismatch")
    if registry.get("records_hash", {}).get("value") != records_sha256(records):
        raise RegistryError("records hash mismatch")


def render_registry(registry: Mapping[str, Any]) -> bytes:
    return (json.dumps(registry, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--generated-at", default=DEFAULT_GENERATED_AT)
    parser.add_argument("--check", action="store_true", help="Fail if --output is not the exact deterministic rebuild")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    repo_root = validate_repo_root(args.repo_root)
    output = args.output if args.output.is_absolute() else repo_root / args.output
    registry = build_registry(repo_root, args.generated_at)
    payload = render_registry(registry)

    print(f"[INPUT] repo_root={repo_root}")
    print("[INPUT] git_source=index_stage_0_blobs")
    print(f"[OUTPUT] workflow_registry={output.resolve()}")
    if args.check:
        if not output.is_file() or output.read_bytes() != payload:
            print("[RESULT] deterministic_rebuild=false", file=sys.stderr)
            return 1
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(payload)

    summary = registry["summary"]
    counts = summary["by_classification"]
    print(
        "[RESULT] "
        f"deterministic_rebuild=true records={summary['total_records']} supported={counts['SUPPORTED']} "
        f"reproducible_legacy={counts['REPRODUCIBLE_LEGACY']} archived={counts['ARCHIVED']} "
        f"tree_content_sha256={registry['tree_content_hash']['value']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
