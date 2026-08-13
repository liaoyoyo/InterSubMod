#!/usr/bin/env python3
"""Archive, repair, and verify release-tree symlinks and local settings.

This tool is intentionally narrow.  It knows the reviewed 2026-08-13
InterSubMod hygiene inventory and refuses to mutate a tree when an unexpected
tracked symlink is present.  Destructive worktree changes happen only after a
recoverable external archive has been written and verified.

The directory digests recorded for external derived assets hash metadata
(relative path, entry type, and byte size), not file contents.  They are useful
for inventory drift detection but are not scientific-content checksums.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Mapping, Sequence


SCHEMA_VERSION = "1.0.0"
DEFAULT_BASELINE = "ddd8909a838318d8a77969313e9561c8ff9d01c2"
DEFAULT_ARCHIVE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/release_archive/InterSubMod/"
    "research-handoff-2026.08.1/repo_hygiene"
)
DEFAULT_RECEIPT_REL = Path("research/20260813_complete_research_handoff/receipts")
SETTINGS_REL = ".claude/settings.local.json"
SETTINGS_TEMPLATE_REL = ".claude/settings.local.example.json"
PROC_SOURCE_PREFIX = "research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"


@dataclass(frozen=True)
class LinkRule:
    action: str
    expected_target: str
    replacement_target: str | None = None
    semantic_description: str | None = None
    availability: str | None = None
    claim_ceiling: str | None = None


REVIEWED_RULES: Mapping[str, LinkRule] = {
    "docs/methodology/_assets/20260618_subcluster_pilot/obs_ws/cpp_wg/figs": LinkRule(
        action="replace_with_pointer",
        expected_target=(
            "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/"
            "docs/methodology/_assets/20260618_subcluster_pilot/figs_cpp_wg_full"
        ),
        semantic_description="Whole-genome cpp_wg derived visualization directory.",
        availability="LOCAL_DERIVED_UNTRACKED",
        claim_ceiling="PARTIAL visualization evidence; not an authority artifact or biological validation.",
    ),
    "docs/methodology/_assets/20260618_subcluster_pilot/obs_ws/cpp_wg/figs_tn": LinkRule(
        action="replace_with_pointer",
        expected_target=(
            "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/"
            "docs/methodology/_assets/20260618_subcluster_pilot/figs_cpp_wg_full_tn"
        ),
        semantic_description="Whole-genome tumor-normal cpp_wg derived visualization directory.",
        availability="LOCAL_DERIVED_UNTRACKED",
        claim_ceiling="PARTIAL visualization evidence; not an authority artifact or biological validation.",
    ),
    "docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v3/figures": LinkRule(
        action="remove_broken",
        expected_target="../v2/figures",
        semantic_description="Broken historical relative link to a missing v2 figure directory.",
        availability="MISSING",
        claim_ceiling="HISTORICAL pointer only; no figure payload was present at the reviewed baseline.",
    ),
    "docs/reports/validated/2026/04/20260406_LOH雙定義交叉分析報告/figures/visual_inspection": LinkRule(
        action="replace_with_pointer",
        expected_target=(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
            "20260406_visual_inspection"
        ),
        semantic_description="LOH visual-inspection derived workspace.",
        availability="LOCAL_CANONICAL",
        claim_ceiling="Derived visual inspection only; consult the cited report and upstream evidence before reuse.",
    ),
    "docs/reports/validated/2026/04/purity_figures": LinkRule(
        action="replace_with_pointer",
        expected_target=(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
            "20260402_phasing_causal_chain/purity_figures"
        ),
        semantic_description="Purity-analysis derived figure bundle.",
        availability="LOCAL_CANONICAL",
        claim_ceiling="Derived visualization bundle; not an authority artifact or independent validation.",
    ),
    "research/tpfp_loh_af_kde_discrimination/scripts/step4_tpfp_discrimination.py": LinkRule(
        action="convert_to_relative",
        expected_target=(
            "/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/scripts/"
            "step4_tpfp_discrimination.py"
        ),
        replacement_target="../../ng_kde_rescaling/scripts/step4_tpfp_discrimination.py",
    ),
    "research/tpfp_loh_af_kde_discrimination/scripts/step4b_tpfp_tomode.py": LinkRule(
        action="convert_to_relative",
        expected_target=(
            "/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/scripts/"
            "step4b_tpfp_tomode.py"
        ),
        replacement_target="../../ng_kde_rescaling/scripts/step4b_tpfp_tomode.py",
    ),
    "research/tpfp_loh_af_kde_discrimination/scripts/step5_visualize_tpfp.py": LinkRule(
        action="convert_to_relative",
        expected_target=(
            "/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/scripts/"
            "step5_visualize_tpfp.py"
        ),
        replacement_target="../../ng_kde_rescaling/scripts/step5_visualize_tpfp.py",
    ),
    "research/tpfp_loh_af_kde_discrimination/scripts/step6_tpfp_detailed.py": LinkRule(
        action="convert_to_relative",
        expected_target=(
            "/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/scripts/"
            "step6_tpfp_detailed.py"
        ),
        replacement_target="../../ng_kde_rescaling/scripts/step6_tpfp_detailed.py",
    ),
    "research/tpfp_loh_af_kde_discrimination/scripts/step7_visualize_detailed.py": LinkRule(
        action="convert_to_relative",
        expected_target=(
            "/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/scripts/"
            "step7_visualize_detailed.py"
        ),
        replacement_target="../../ng_kde_rescaling/scripts/step7_visualize_detailed.py",
    ),
    "research/tpfp_loh_af_kde_discrimination/scripts/step8_multidim_panorama.py": LinkRule(
        action="convert_to_relative",
        expected_target=(
            "/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/scripts/"
            "step8_multidim_panorama.py"
        ),
        replacement_target="../../ng_kde_rescaling/scripts/step8_multidim_panorama.py",
    ),
    "research/tpfp_loh_af_kde_discrimination/scripts/utils_io.py": LinkRule(
        action="convert_to_relative",
        expected_target=(
            "/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/scripts/utils_io.py"
        ),
        replacement_target="../../ng_kde_rescaling/scripts/utils_io.py",
    ),
}


EXPECTED_COUNTS = {
    "tracked_symlinks": 512,
    "remove_proc_self_fd": 500,
    "remove_broken": 1,
    "convert_to_relative": 7,
    "replace_with_pointer": 4,
}


class HygieneError(RuntimeError):
    """Raised when a fail-closed hygiene invariant is not satisfied."""


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def run_git(repo: Path, args: Sequence[str], *, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=repo,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def git_text(repo: Path, args: Sequence[str]) -> str:
    return run_git(repo, args).stdout.decode("utf-8", "surrogateescape").strip()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def safe_repo_path(repo: Path, relative_path: str) -> Path:
    pure = PurePosixPath(relative_path)
    if pure.is_absolute() or not pure.parts or ".." in pure.parts:
        raise HygieneError(f"unsafe repository-relative path: {relative_path!r}")
    candidate = repo.joinpath(*pure.parts)
    parent = candidate.parent.resolve()
    root = repo.resolve()
    if parent != root and root not in parent.parents:
        raise HygieneError(f"path parent escapes repository: {relative_path!r}")
    return candidate


def write_json(path: Path, payload: Mapping[str, Any], *, mode: int = 0o644) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = (json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("wb") as handle:
        handle.write(encoded)
        handle.flush()
        os.fsync(handle.fileno())
    os.chmod(temporary, mode)
    os.replace(temporary, path)


def tracked_symlink_entries(repo: Path) -> list[dict[str, Any]]:
    raw = run_git(repo, ["ls-files", "-s", "-z"]).stdout
    records: list[dict[str, Any]] = []
    for item in raw.split(b"\0"):
        if not item:
            continue
        metadata, path_bytes = item.split(b"\t", 1)
        mode, blob_oid, stage = metadata.decode("ascii").split()
        if mode != "120000" or stage != "0":
            continue
        relative_path = path_bytes.decode("utf-8", "surrogateescape")
        link_path = safe_repo_path(repo, relative_path)
        if not link_path.is_symlink():
            raise HygieneError(
                f"tracked symlink is absent or already replaced in the worktree: {relative_path}"
            )
        target = os.readlink(link_path)
        records.append(
            {
                "path": relative_path,
                "git_blob_oid": blob_oid,
                "original_target": target,
                "original_target_is_absolute": os.path.isabs(target),
                "original_target_exists": link_path.exists(),
            }
        )
    return sorted(records, key=lambda record: record["path"])


def directory_metadata_inventory(root: Path) -> dict[str, Any]:
    """Return a content-light inventory without reading file payloads."""

    if not root.exists():
        return {
            "exists": False,
            "entry_type": "MISSING",
            "file_count": 0,
            "directory_count": 0,
            "symlink_count": 0,
            "total_file_bytes": 0,
            "directory_manifest_sha256": None,
            "hash_strategy": "not_applicable_missing_target",
        }
    if not root.is_dir():
        stat_result = root.stat()
        return {
            "exists": True,
            "entry_type": "FILE",
            "file_count": 1,
            "directory_count": 0,
            "symlink_count": 0,
            "total_file_bytes": stat_result.st_size,
            "directory_manifest_sha256": None,
            "content_sha256": sha256_file(root),
            "hash_strategy": "sha256_file_content",
        }

    digest = hashlib.sha256()
    file_count = 0
    directory_count = 0
    symlink_count = 0
    total_file_bytes = 0
    root_resolved = root.resolve()

    for current, directories, files in os.walk(root_resolved, followlinks=False):
        directories.sort()
        files.sort()
        current_path = Path(current)
        for name in list(directories):
            path = current_path / name
            relative = path.relative_to(root_resolved).as_posix()
            if path.is_symlink():
                symlink_count += 1
                directories.remove(name)
                line = f"L\t{relative}\t{os.readlink(path)}\n"
            else:
                directory_count += 1
                line = f"D\t{relative}\n"
            digest.update(line.encode("utf-8", "surrogateescape"))
        for name in files:
            path = current_path / name
            relative = path.relative_to(root_resolved).as_posix()
            if path.is_symlink():
                symlink_count += 1
                line = f"L\t{relative}\t{os.readlink(path)}\n"
            else:
                size = path.stat().st_size
                file_count += 1
                total_file_bytes += size
                line = f"F\t{relative}\t{size}\n"
            digest.update(line.encode("utf-8", "surrogateescape"))

    return {
        "exists": True,
        "entry_type": "DIRECTORY",
        "file_count": file_count,
        "directory_count": directory_count,
        "symlink_count": symlink_count,
        "total_file_bytes": total_file_bytes,
        "directory_manifest_sha256": digest.hexdigest(),
        "hash_strategy": (
            "metadata_sha256_v1 over sorted entry_type, relative_path, file_size_or_link_target; "
            "file contents are intentionally not read or hashed"
        ),
    }


def resolve_original_target(link_path: Path, target: str) -> Path:
    if os.path.isabs(target):
        return Path(target)
    return link_path.parent / target


def build_plan(
    repo: Path,
    *,
    rules: Mapping[str, LinkRule] = REVIEWED_RULES,
    expected_counts: Mapping[str, int] | None = EXPECTED_COUNTS,
    proc_source_prefix: str = PROC_SOURCE_PREFIX,
) -> list[dict[str, Any]]:
    entries = tracked_symlink_entries(repo)
    plan: list[dict[str, Any]] = []
    seen_rules: set[str] = set()

    for entry in entries:
        relative_path = entry["path"]
        original_target = entry["original_target"]
        link_path = safe_repo_path(repo, relative_path)
        rule = rules.get(relative_path)
        if rule is not None:
            seen_rules.add(relative_path)
            if original_target != rule.expected_target:
                raise HygieneError(
                    f"reviewed target changed for {relative_path}: "
                    f"expected {rule.expected_target!r}, got {original_target!r}"
                )
            operation = {**entry, **asdict(rule)}
        elif relative_path.startswith(proc_source_prefix) and original_target.startswith("/proc/self/fd/"):
            if entry["original_target_exists"]:
                raise HygieneError(f"/proc/self/fd link unexpectedly resolves: {relative_path}")
            operation = {
                **entry,
                "action": "remove_proc_self_fd",
                "replacement_target": None,
                "semantic_description": "Broken runtime file-descriptor link accidentally committed as output.",
                "availability": "INVALIDATED",
                "claim_ceiling": "No scientific evidence; runtime residue only.",
            }
        else:
            raise HygieneError(
                f"unreviewed tracked symlink refuses automatic cleanup: {relative_path} -> {original_target}"
            )

        target_path = resolve_original_target(link_path, original_target)
        if operation["action"] == "replace_with_pointer":
            operation["target_inventory"] = directory_metadata_inventory(target_path)
        elif operation["action"] == "convert_to_relative":
            replacement = operation["replacement_target"]
            if not replacement or os.path.isabs(replacement):
                raise HygieneError(f"invalid relative replacement for {relative_path}")
            replacement_path = link_path.parent / replacement
            if not replacement_path.is_file():
                raise HygieneError(f"reviewed repository target is missing: {replacement_path}")
            operation["replacement_target_sha256"] = sha256_file(replacement_path)
            operation["target_inventory"] = {
                "exists": True,
                "entry_type": "FILE",
                "file_count": 1,
                "directory_count": 0,
                "symlink_count": 0,
                "total_file_bytes": replacement_path.stat().st_size,
                "content_sha256": operation["replacement_target_sha256"],
                "hash_strategy": "sha256_file_content",
            }
        else:
            operation["target_inventory"] = directory_metadata_inventory(target_path)
        plan.append(operation)

    missing_rules = sorted(set(rules) - seen_rules)
    if missing_rules:
        raise HygieneError(f"reviewed symlink paths missing from tracked inventory: {missing_rules}")

    counts = Counter(operation["action"] for operation in plan)
    observed = {"tracked_symlinks": len(plan), **counts}
    if expected_counts is not None:
        for key, expected in expected_counts.items():
            actual = observed.get(key, 0)
            if actual != expected:
                raise HygieneError(f"inventory count mismatch for {key}: expected {expected}, got {actual}")
    return plan


def current_git_context(repo: Path, baseline: str) -> dict[str, Any]:
    head = git_text(repo, ["rev-parse", "HEAD"])
    baseline_full = git_text(repo, ["rev-parse", baseline])
    ancestor = run_git(repo, ["merge-base", "--is-ancestor", baseline_full, head], check=False)
    if ancestor.returncode != 0:
        raise HygieneError(f"baseline {baseline_full} is not an ancestor of HEAD {head}")
    status_raw = run_git(repo, ["status", "--porcelain=v1", "-z", "--untracked-files=all"]).stdout
    status_entries = [
        item.decode("utf-8", "surrogateescape") for item in status_raw.split(b"\0") if item
    ]
    return {
        "repo_root": str(repo.resolve()),
        "head": head,
        "baseline": baseline_full,
        "baseline_is_ancestor": True,
        "branch": git_text(repo, ["branch", "--show-current"]),
        "git_status_entries": status_entries,
    }


def settings_record(repo: Path) -> dict[str, Any]:
    path = safe_repo_path(repo, SETTINGS_REL)
    tracked = run_git(repo, ["ls-files", "--error-unmatch", "--", SETTINGS_REL], check=False)
    if tracked.returncode != 0 or not path.is_file() or path.is_symlink():
        raise HygieneError(f"expected tracked local settings file is missing: {SETTINGS_REL}")
    mode = stat.S_IMODE(path.stat().st_mode)
    return {
        "path": SETTINGS_REL,
        "tracked_in_index": True,
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "original_mode_octal": f"{mode:04o}",
        "archive_relative_path": "private/.claude/settings.local.json",
        "sensitive_payload_in_git_receipt": False,
    }


def make_before_receipt(repo: Path, baseline: str, plan: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    counts = Counter(operation["action"] for operation in plan)
    return {
        "schema_version": SCHEMA_VERSION,
        "receipt_type": "repo_hygiene_before",
        "generated_at": now_iso(),
        "hostname": platform.node(),
        "task_type": ["B_COMPREHENSIVE_VALIDATION", "D_EXTERNAL_HANDOFF"],
        "claim_ceiling": (
            "Repository hygiene only; this receipt does not validate or regenerate scientific results."
        ),
        "git": current_git_context(repo, baseline),
        "counts": {"tracked_symlinks": len(plan), **dict(sorted(counts.items()))},
        "settings": settings_record(repo),
        "symlink_operations": list(plan),
    }


def build_archive(
    repo: Path,
    archive_root: Path,
    before_receipt: Mapping[str, Any],
) -> tuple[Path, dict[str, Any]]:
    if archive_root.exists() or archive_root.is_symlink():
        raise HygieneError(f"archive target already exists; refusing overwrite: {archive_root}")
    archive_parent = archive_root.parent
    archive_parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".repo_hygiene-", dir=archive_parent))
    os.chmod(temporary, 0o700)
    try:
        archived_settings = temporary / before_receipt["settings"]["archive_relative_path"]
        archived_settings.parent.mkdir(parents=True, exist_ok=True)
        os.chmod(archived_settings.parent, 0o700)
        shutil.copyfile(safe_repo_path(repo, SETTINGS_REL), archived_settings)
        os.chmod(archived_settings, 0o600)
        archived_hash = sha256_file(archived_settings)
        if archived_hash != before_receipt["settings"]["sha256"]:
            raise HygieneError("archived settings checksum mismatch before mutation")

        manifest = {
            "schema_version": SCHEMA_VERSION,
            "manifest_type": "repo_hygiene_external_archive",
            "created_at": now_iso(),
            "hostname": platform.node(),
            "archive_root": str(archive_root.resolve(strict=False)),
            "source_repo": before_receipt["git"],
            "claim_ceiling": before_receipt["claim_ceiling"],
            "settings": before_receipt["settings"],
            "symlink_operations": before_receipt["symlink_operations"],
            "restore_contract": (
                "The archived settings payload plus per-link original target metadata can restore the "
                "reviewed worktree paths. External target payloads were not copied."
            ),
        }
        manifest_path = temporary / "archive_manifest.json"
        write_json(manifest_path, manifest, mode=0o600)
        os.replace(temporary, archive_root)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise

    final_manifest = archive_root / "archive_manifest.json"
    verified = verify_archive(archive_root)
    if not verified["pass"]:
        raise HygieneError(f"external archive verification failed: {verified['errors']}")
    return final_manifest, json.loads(final_manifest.read_text(encoding="utf-8"))


def verify_archive(archive_root: Path) -> dict[str, Any]:
    errors: list[str] = []
    manifest_path = archive_root / "archive_manifest.json"
    if not manifest_path.is_file():
        return {"pass": False, "errors": [f"missing archive manifest: {manifest_path}"]}
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return {"pass": False, "errors": [f"invalid archive manifest: {error}"]}

    operations = manifest.get("symlink_operations", [])
    paths = [operation.get("path") for operation in operations]
    if len(paths) != len(set(paths)):
        errors.append("duplicate symlink paths in archive manifest")
    settings = manifest.get("settings", {})
    archived_settings = archive_root / str(settings.get("archive_relative_path", ""))
    if not archived_settings.is_file():
        errors.append(f"missing archived settings payload: {archived_settings}")
    elif sha256_file(archived_settings) != settings.get("sha256"):
        errors.append("archived settings SHA-256 mismatch")
    if manifest.get("schema_version") != SCHEMA_VERSION:
        errors.append("unsupported archive schema version")
    return {
        "pass": not errors,
        "errors": errors,
        "archive_root": str(archive_root.resolve()),
        "manifest_path": str(manifest_path.resolve()),
        "manifest_sha256": sha256_file(manifest_path),
        "symlink_record_count": len(operations),
        "settings_payload_sha256": settings.get("sha256"),
    }


def pointer_payload(
    operation: Mapping[str, Any], source_head: str, *, verified_at: str | None = None
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_type": "EXTERNAL_DERIVED_DIRECTORY_POINTER",
        "semantic_description": operation.get("semantic_description"),
        "original_machine_location": operation["original_target"],
        "availability": operation.get("availability"),
        "source_head_when_reviewed": source_head,
        "verified_at": verified_at or now_iso(),
        "target_inventory": operation.get("target_inventory"),
        "claim_ceiling": operation.get("claim_ceiling"),
        "payload_policy": "The external payload is intentionally not copied into Git.",
    }


def pointer_readme(operation: Mapping[str, Any]) -> str:
    target = operation["original_target"]
    description = operation.get("semantic_description") or "External derived directory."
    claim_ceiling = operation.get("claim_ceiling") or "No claim beyond the linked provenance."
    return f"""<!--
建立時間: 2026-08-13
目標: 取代不可攜 absolute symlink，保留外部衍生資料的可追溯 pointer
處理範圍: {operation['path']}
關聯檔案: InterSubMod/research/20260813_complete_research_handoff/receipts/repo_hygiene_archive_manifest_20260813.json
驗證方式: 讀取同目錄 ARTIFACT_POINTER.json 的 path/count/size/metadata digest
證據等級: L3（provenance pointer；不含科研 payload）
-->

# 外部衍生資料指標

{description}

- 原始本機路徑：`{target}`
- 可用性：`{operation.get('availability')}`
- 科學 claim ceiling：{claim_ceiling}
- Git 政策：此目錄只保存 pointer 與 inventory；不複製大型衍生內容。

機器可讀的檔數、總大小及 metadata digest 請見 `ARTIFACT_POINTER.json`。該 digest
只涵蓋相對路徑、entry type 與檔案大小，不是逐檔內容 SHA-256。
"""


def apply_operations(repo: Path, manifest: Mapping[str, Any]) -> None:
    source_head = str(manifest["source_repo"]["head"])
    if git_text(repo, ["rev-parse", "HEAD"]) != source_head:
        raise HygieneError("HEAD changed after archive creation; refusing worktree mutation")

    for operation in manifest["symlink_operations"]:
        relative_path = operation["path"]
        path = safe_repo_path(repo, relative_path)
        if not path.is_symlink() or os.readlink(path) != operation["original_target"]:
            raise HygieneError(f"link changed after archive creation: {relative_path}")

    settings_path = safe_repo_path(repo, SETTINGS_REL)
    if sha256_file(settings_path) != manifest["settings"]["sha256"]:
        raise HygieneError("settings changed after archive creation")

    for operation in manifest["symlink_operations"]:
        path = safe_repo_path(repo, operation["path"])
        action = operation["action"]
        path.unlink()
        if action == "convert_to_relative":
            os.symlink(operation["replacement_target"], path)
        elif action == "replace_with_pointer":
            path.mkdir()
            write_json(
                path / "ARTIFACT_POINTER.json",
                pointer_payload(operation, source_head, verified_at=str(manifest["created_at"])),
            )
            (path / "README.md").write_text(pointer_readme(operation), encoding="utf-8")
        elif action in {"remove_proc_self_fd", "remove_broken"}:
            continue
        else:
            raise HygieneError(f"unsupported action after archive preflight: {action}")

    settings_path.unlink()


def versioned_worktree_symlinks(repo: Path) -> list[dict[str, Any]]:
    raw = run_git(repo, ["ls-files", "-z"]).stdout
    paths = [item.decode("utf-8", "surrogateescape") for item in raw.split(b"\0") if item]
    records: list[dict[str, Any]] = []
    for relative_path in paths:
        path = safe_repo_path(repo, relative_path)
        if path.is_symlink():
            target = os.readlink(path)
            records.append(
                {
                    "path": relative_path,
                    "target": target,
                    "target_is_absolute": os.path.isabs(target),
                    "target_exists": path.exists(),
                }
            )
    return sorted(records, key=lambda record: record["path"])


def verify_worktree(repo: Path, archive_root: Path, baseline: str) -> dict[str, Any]:
    archive_check = verify_archive(archive_root)
    manifest_path = archive_root / "archive_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    errors: list[str] = [] if archive_check["pass"] else list(archive_check["errors"])
    operations = manifest.get("symlink_operations", [])
    expected_relative_links = 0
    pointer_count = 0
    absent_count = 0

    for operation in operations:
        path = safe_repo_path(repo, operation["path"])
        action = operation["action"]
        if action == "convert_to_relative":
            expected_relative_links += 1
            if not path.is_symlink():
                errors.append(f"relative replacement is not a symlink: {operation['path']}")
            elif os.readlink(path) != operation["replacement_target"]:
                errors.append(f"relative replacement target mismatch: {operation['path']}")
            elif os.path.isabs(os.readlink(path)) or not path.exists():
                errors.append(f"relative replacement is absolute or broken: {operation['path']}")
        elif action == "replace_with_pointer":
            pointer_count += 1
            expected = {"ARTIFACT_POINTER.json", "README.md"}
            if not path.is_dir() or path.is_symlink():
                errors.append(f"pointer directory missing: {operation['path']}")
            elif {item.name for item in path.iterdir()} != expected:
                errors.append(f"pointer directory content mismatch: {operation['path']}")
        elif action in {"remove_proc_self_fd", "remove_broken"}:
            absent_count += 1
            if path.exists() or path.is_symlink():
                errors.append(f"removed link path still exists: {operation['path']}")

    settings_path = safe_repo_path(repo, SETTINGS_REL)
    template_path = safe_repo_path(repo, SETTINGS_TEMPLATE_REL)
    if settings_path.exists() or settings_path.is_symlink():
        errors.append(f"tracked local settings still present: {SETTINGS_REL}")
    if not template_path.is_file():
        errors.append(f"safe settings template missing: {SETTINGS_TEMPLATE_REL}")

    current_links = versioned_worktree_symlinks(repo)
    absolute_count = sum(record["target_is_absolute"] for record in current_links)
    broken_count = sum(not record["target_exists"] for record in current_links)
    if absolute_count:
        errors.append(f"versioned worktree still contains {absolute_count} absolute symlinks")
    if broken_count:
        errors.append(f"versioned worktree still contains {broken_count} broken symlinks")
    if len(current_links) != expected_relative_links:
        errors.append(
            f"expected {expected_relative_links} valid relative symlinks, found {len(current_links)}"
        )

    git_context = current_git_context(repo, baseline)
    # --no-index checks the rule even while the deletion is not staged yet.
    ignore_check = run_git(
        repo, ["check-ignore", "--no-index", "-q", "--", SETTINGS_REL], check=False
    )
    if ignore_check.returncode != 0:
        errors.append(
            "IGNORE_RULE_PENDING: .claude/settings.local.json is not ignored; parent integration owns .gitignore"
        )

    return {
        "schema_version": SCHEMA_VERSION,
        "receipt_type": "repo_hygiene_after",
        "verified_at": now_iso(),
        "hostname": platform.node(),
        "claim_ceiling": (
            "Repository hygiene only; this receipt does not validate or regenerate scientific results."
        ),
        "pass": not errors,
        "pass_excluding_parent_owned_ignore_rule": (
            not [error for error in errors if not error.startswith("IGNORE_RULE_PENDING:")]
        ),
        "errors": errors,
        "git": git_context,
        "archive_verification": archive_check,
        "operation_counts": {
            "valid_relative_symlinks": expected_relative_links,
            "pointer_directories": pointer_count,
            "removed_paths": absent_count,
        },
        "versioned_worktree_symlinks": current_links,
        "versioned_worktree_absolute_symlink_count": absolute_count,
        "versioned_worktree_broken_symlink_count": broken_count,
        "settings_local_present": settings_path.exists() or settings_path.is_symlink(),
        "settings_template_present": template_path.is_file(),
        "settings_local_ignore_rule_present": ignore_check.returncode == 0,
    }


def restore_from_archive(repo: Path, archive_root: Path, *, confirm: bool) -> dict[str, Any]:
    if not confirm:
        raise HygieneError("restore requires --confirm-restore")
    archive_check = verify_archive(archive_root)
    if not archive_check["pass"]:
        raise HygieneError(f"cannot restore an invalid archive: {archive_check['errors']}")
    manifest = json.loads((archive_root / "archive_manifest.json").read_text(encoding="utf-8"))

    restored_symlinks = 0
    already_original_symlinks = 0
    for operation in manifest["symlink_operations"]:
        path = safe_repo_path(repo, operation["path"])
        action = operation["action"]
        if path.is_symlink() and os.readlink(path) == operation["original_target"]:
            already_original_symlinks += 1
            continue
        if action == "convert_to_relative":
            if not path.is_symlink() or os.readlink(path) != operation["replacement_target"]:
                raise HygieneError(f"unsafe restore state for relative link: {operation['path']}")
            path.unlink()
        elif action == "replace_with_pointer":
            expected = {"ARTIFACT_POINTER.json", "README.md"}
            if not path.is_dir() or path.is_symlink():
                raise HygieneError(f"unsafe restore state for pointer directory: {operation['path']}")
            present = {item.name for item in path.iterdir()}
            if not present.issubset(expected):
                raise HygieneError(f"unexpected file in pointer directory: {operation['path']}")
            expected_pointer = (
                json.dumps(
                    pointer_payload(
                        operation,
                        str(manifest["source_repo"]["head"]),
                        verified_at=str(manifest["created_at"]),
                    ),
                    ensure_ascii=False,
                    indent=2,
                    sort_keys=True,
                )
                + "\n"
            ).encode("utf-8")
            expected_readme = pointer_readme(operation).encode("utf-8")
            if "ARTIFACT_POINTER.json" in present and sha256_file(
                path / "ARTIFACT_POINTER.json"
            ) != sha256_bytes(expected_pointer):
                raise HygieneError(f"pointer JSON changed; refusing restore: {operation['path']}")
            if "README.md" in present and sha256_file(path / "README.md") != sha256_bytes(
                expected_readme
            ):
                raise HygieneError(f"pointer README changed; refusing restore: {operation['path']}")
            for name in sorted(present):
                (path / name).unlink()
            path.rmdir()
        elif action in {"remove_proc_self_fd", "remove_broken"}:
            if path.exists() or path.is_symlink():
                raise HygieneError(f"unsafe restore state for removed link: {operation['path']}")
        os.symlink(operation["original_target"], path)
        restored_symlinks += 1

    settings = manifest["settings"]
    settings_path = safe_repo_path(repo, settings["path"])
    settings_restored = False
    if settings_path.exists() or settings_path.is_symlink():
        if not settings_path.is_file() or settings_path.is_symlink():
            raise HygieneError(f"unsafe existing settings restore target: {settings_path}")
        if sha256_file(settings_path) != settings["sha256"]:
            raise HygieneError(f"existing settings differ from archive: {settings_path}")
    else:
        archived_settings = archive_root / settings["archive_relative_path"]
        shutil.copyfile(archived_settings, settings_path)
        os.chmod(settings_path, int(settings["original_mode_octal"], 8))
        if sha256_file(settings_path) != settings["sha256"]:
            raise HygieneError("restored settings checksum mismatch")
        settings_restored = True
    return {
        "pass": True,
        "restored_at": now_iso(),
        "restored_symlinks": restored_symlinks,
        "already_original_symlinks": already_original_symlinks,
        "settings_restored": settings_restored,
        "restored_settings_sha256": settings["sha256"],
    }


def copy_archive_manifest_to_receipts(manifest_path: Path, receipt_path: Path) -> None:
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(manifest_path, receipt_path)
    if sha256_file(manifest_path) != sha256_file(receipt_path):
        raise HygieneError("Git receipt copy of archive manifest does not match external archive")


def command_scan(args: argparse.Namespace) -> int:
    repo = args.repo_root.resolve()
    print(f"[INPUT] repo_root={repo}")
    print(f"[INPUT] expected_baseline={args.expected_baseline}")
    plan = build_plan(repo)
    receipt = make_before_receipt(repo, args.expected_baseline, plan)
    if args.output_json:
        write_json(args.output_json, receipt)
        print(f"[OUTPUT] before_receipt={args.output_json.resolve()}")
    print(f"[RESULT] tracked_symlinks={len(plan)} actions={json.dumps(receipt['counts'], sort_keys=True)}")
    return 0


def command_apply(args: argparse.Namespace) -> int:
    repo = args.repo_root.resolve()
    archive_root = args.archive_root.resolve(strict=False)
    receipt_dir = (repo / args.receipt_dir).resolve()
    print(f"[INPUT] repo_root={repo}")
    print(f"[INPUT] expected_baseline={args.expected_baseline}")
    print(f"[INPUT] archive_root={archive_root}")
    print(f"[OUTPUT] receipt_dir={receipt_dir}")

    plan = build_plan(repo)
    before = make_before_receipt(repo, args.expected_baseline, plan)
    before_path = receipt_dir / "repo_hygiene_before_20260813.json"
    write_json(before_path, before)
    manifest_path, manifest = build_archive(repo, archive_root, before)
    archive_receipt = receipt_dir / "repo_hygiene_archive_manifest_20260813.json"
    copy_archive_manifest_to_receipts(manifest_path, archive_receipt)
    print(f"[ARCHIVE] verified_manifest={manifest_path}")
    print(f"[ARCHIVE] settings_sha256={manifest['settings']['sha256']}")

    apply_operations(repo, manifest)
    after = verify_worktree(repo, archive_root, args.expected_baseline)
    after_path = receipt_dir / "repo_hygiene_after_20260813.json"
    write_json(after_path, after)
    print(f"[OUTPUT] before_receipt={before_path}")
    print(f"[OUTPUT] archive_manifest_receipt={archive_receipt}")
    print(f"[OUTPUT] after_receipt={after_path}")
    print(
        "[RESULT] "
        f"pass={str(after['pass']).lower()} "
        f"pass_excluding_parent_owned_ignore_rule="
        f"{str(after['pass_excluding_parent_owned_ignore_rule']).lower()} "
        f"absolute={after['versioned_worktree_absolute_symlink_count']} "
        f"broken={after['versioned_worktree_broken_symlink_count']}"
    )
    for error in after["errors"]:
        print(f"[GATE] {error}")
    return 0 if after["pass_excluding_parent_owned_ignore_rule"] else 1


def command_verify(args: argparse.Namespace) -> int:
    repo = args.repo_root.resolve()
    print(f"[INPUT] repo_root={repo}")
    print(f"[INPUT] archive_root={args.archive_root.resolve(strict=False)}")
    receipt = verify_worktree(repo, args.archive_root, args.expected_baseline)
    if args.output_json:
        write_json(args.output_json, receipt)
        print(f"[OUTPUT] after_receipt={args.output_json.resolve()}")
    print(
        "[RESULT] "
        f"pass={str(receipt['pass']).lower()} "
        f"pass_excluding_parent_owned_ignore_rule="
        f"{str(receipt['pass_excluding_parent_owned_ignore_rule']).lower()}"
    )
    for error in receipt["errors"]:
        print(f"[GATE] {error}")
    return 0 if receipt["pass_excluding_parent_owned_ignore_rule"] else 1


def command_verify_archive(args: argparse.Namespace) -> int:
    print(f"[INPUT] archive_root={args.archive_root.resolve(strict=False)}")
    result = verify_archive(args.archive_root)
    print(f"[RESULT] {json.dumps(result, ensure_ascii=False, sort_keys=True)}")
    return 0 if result["pass"] else 1


def command_restore(args: argparse.Namespace) -> int:
    print(f"[INPUT] repo_root={args.repo_root.resolve()}")
    print(f"[INPUT] archive_root={args.archive_root.resolve(strict=False)}")
    result = restore_from_archive(args.repo_root.resolve(), args.archive_root, confirm=args.confirm_restore)
    print(f"[RESULT] {json.dumps(result, ensure_ascii=False, sort_keys=True)}")
    return 0


def add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--expected-baseline", default=DEFAULT_BASELINE)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    scan_parser = subparsers.add_parser("scan", help="Validate and report the reviewed pre-cleanup inventory.")
    add_common(scan_parser)
    scan_parser.add_argument("--output-json", type=Path)
    scan_parser.set_defaults(handler=command_scan)

    apply_parser = subparsers.add_parser("apply", help="Archive first, then apply the reviewed cleanup.")
    add_common(apply_parser)
    apply_parser.add_argument("--archive-root", type=Path, default=DEFAULT_ARCHIVE_ROOT)
    apply_parser.add_argument("--receipt-dir", type=Path, default=DEFAULT_RECEIPT_REL)
    apply_parser.set_defaults(handler=command_apply)

    verify_parser = subparsers.add_parser("verify", help="Verify the cleaned worktree and external archive.")
    add_common(verify_parser)
    verify_parser.add_argument("--archive-root", type=Path, default=DEFAULT_ARCHIVE_ROOT)
    verify_parser.add_argument("--output-json", type=Path)
    verify_parser.set_defaults(handler=command_verify)

    archive_parser = subparsers.add_parser("verify-archive", help="Verify archive recoverability metadata.")
    archive_parser.add_argument("--archive-root", type=Path, default=DEFAULT_ARCHIVE_ROOT)
    archive_parser.set_defaults(handler=command_verify_archive)

    restore_parser = subparsers.add_parser("restore", help="Restore reviewed paths from the external archive.")
    restore_parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    restore_parser.add_argument("--archive-root", type=Path, default=DEFAULT_ARCHIVE_ROOT)
    restore_parser.add_argument("--confirm-restore", action="store_true")
    restore_parser.set_defaults(handler=command_restore)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        return int(args.handler(args))
    except (HygieneError, OSError, subprocess.CalledProcessError, json.JSONDecodeError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
