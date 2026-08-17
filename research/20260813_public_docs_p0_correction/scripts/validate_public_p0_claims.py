#!/usr/bin/env python3
"""Fail-closed validator for the 2026-08-13 public-document P0 corrections."""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


ALLOWED_DISPOSITIONS = {
    "local_source_fixed_external_publish_required",
    "external_action_required",
}
EXPECTED_EXTERNAL_SURFACES = {"GitHub About", "default branch", "GitHub Wiki", "GitHub Pages"}


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--registry",
        type=Path,
        default=script_path.parent / "p0_claim_registry.json",
        help="Machine-readable P0 claim registry.",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=script_path.parents[3],
        help="InterSubMod repository root.",
    )
    parser.add_argument(
        "--simulate-drop-claim",
        metavar="CLAIM_ID",
        help="Diagnostic only: remove one registry claim in memory to prove fail-closed behavior.",
    )
    return parser.parse_args()


def normalize_document(path: Path) -> str:
    """Normalize Markdown/HTML enough for stable semantic anchor checks."""
    raw = path.read_text(encoding="utf-8")
    decoded = html.unescape(raw)
    without_tags = re.sub(r"<[^>]+>", " ", decoded)
    # Presentation markup must not let a false claim bypass semantic guards.
    # Remove Markdown emphasis/code delimiters while retaining underscores in
    # identifiers such as significance_summary.csv.
    without_markdown_markup = without_tags.replace("`", "").replace("*", "")
    without_markdown_markup = without_markdown_markup.replace("~~", "")
    return re.sub(r"\s+", " ", without_markdown_markup).strip()


def load_inventory_p0(repo_root: Path, source: str) -> set[str]:
    path = resolve_repo_relative_file(repo_root, source, "source inventory")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {row["claim_id"] for row in rows if row.get("priority") == "P0"}


def resolve_repo_relative_file(repo_root: Path, value: str, label: str) -> Path:
    """Resolve one declared input without allowing absolute or parent escapes."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty repo-relative path")
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"{label} must be repo-relative without '..': {value}")
    resolved_root = repo_root.resolve()
    resolved = (resolved_root / relative).resolve()
    try:
        resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(f"{label} escapes repository root: {value}") from exc
    if not resolved.is_file() or resolved.is_symlink():
        raise FileNotFoundError(f"{label} is not a non-symlink regular file: {value}")
    return resolved


def validate_regex(pattern: str, where: str, errors: list[str]) -> re.Pattern[str] | None:
    try:
        return re.compile(pattern, flags=re.IGNORECASE | re.MULTILINE)
    except re.error as exc:
        errors.append(f"{where}: invalid regex {pattern!r}: {exc}")
        return None


def expand_rule_targets(
    rule: dict[str, Any], repo_root: Path, location: str, errors: list[str]
) -> list[Path]:
    """Resolve an explicitly bounded target list without permitting paths outside the repo."""
    resolved: set[Path] = set()
    targets = rule.get("targets", [])
    target_globs = rule.get("target_globs", [])
    if not isinstance(targets, list) or not all(isinstance(item, str) and item.strip() for item in targets):
        errors.append(f"{location}: targets must be a string list")
        targets = []
    if not isinstance(target_globs, list) or not all(
        isinstance(item, str) and item.strip() for item in target_globs
    ):
        errors.append(f"{location}: target_globs must be a string list")
        target_globs = []
    if not targets and not target_globs:
        errors.append(f"{location}: at least one target or target_glob is required")
        return []

    def accept(path: Path, source: str) -> None:
        try:
            path.relative_to(repo_root)
        except ValueError:
            errors.append(f"{location}: target escapes repository root: {source}")
            return
        if not path.is_file():
            errors.append(f"{location}: target is not a regular file: {source}")
            return
        resolved.add(path)

    for target in targets:
        candidate = Path(target)
        if candidate.is_absolute() or ".." in candidate.parts:
            errors.append(f"{location}: target must be repo-relative without '..': {target}")
            continue
        accept((repo_root / candidate).resolve(), target)
    for pattern in target_globs:
        candidate = Path(pattern)
        if candidate.is_absolute() or ".." in candidate.parts:
            errors.append(f"{location}: glob must be repo-relative without '..': {pattern}")
            continue
        matches = [path.resolve() for path in repo_root.glob(pattern) if path.is_file()]
        if not matches:
            errors.append(f"{location}: target_glob matched no files: {pattern}")
            continue
        for path in matches:
            accept(path, pattern)
    return sorted(resolved)


def validate_pattern_list(value: Any, location: str, errors: list[str], *, required: bool) -> list[str]:
    if not isinstance(value, list) or (required and not value):
        errors.append(f"{location}: must be {'a non-empty' if required else 'a'} regex list")
        return []
    patterns: list[str] = []
    for index, pattern in enumerate(value):
        if not isinstance(pattern, str) or not pattern:
            errors.append(f"{location}[{index}]: must be a non-empty regex")
            continue
        if validate_regex(pattern, f"{location}[{index}]", errors) is not None:
            patterns.append(pattern)
    return patterns


def validate_cross_document_guards(
    registry: dict[str, Any], repo_root: Path
) -> tuple[list[str], dict[str, Any]]:
    """Validate the registry-owned cross-document and repository shipment contract."""
    errors: list[str] = []
    contract = registry.get("cross_document_guards")
    if not isinstance(contract, dict):
        return ["cross_document_guards must be an object"], {}

    public_globs = contract.get("public_source_globs")
    public_paths = expand_rule_targets(
        {"target_globs": public_globs}, repo_root, "cross_document_guards.public_source_globs", errors
    )
    raw_cache: dict[Path, str] = {}
    normalized_cache: dict[Path, str] = {}

    def text_for(path: Path, view: str) -> str:
        if view == "raw":
            if path not in raw_cache:
                raw_cache[path] = path.read_text(encoding="utf-8")
            return raw_cache[path]
        if view != "normalized":
            errors.append(f"unsupported content view {view!r}; expected raw or normalized")
        if path not in normalized_cache:
            normalized_cache[path] = normalize_document(path)
        return normalized_cache[path]

    checked_documents = 0
    required_anchors = 0
    forbidden_anchors = 0
    context_occurrences = 0
    guard_ids: list[str] = []

    for group_name in ("document_rules", "corpus_rules"):
        rules = contract.get(group_name)
        if not isinstance(rules, list) or not rules:
            errors.append(f"cross_document_guards.{group_name} must be a non-empty list")
            continue
        for index, rule in enumerate(rules):
            location = f"cross_document_guards.{group_name}[{index}]"
            if not isinstance(rule, dict):
                errors.append(f"{location}: rule must be an object")
                continue
            guard_id = rule.get("guard_id")
            if not isinstance(guard_id, str) or not guard_id.strip():
                errors.append(f"{location}: guard_id must be a non-empty string")
                guard_id = location
            guard_ids.append(guard_id)
            paths = expand_rule_targets(rule, repo_root, location, errors)
            required = validate_pattern_list(
                rule.get("required_all", []), f"{location}.required_all", errors, required=group_name == "document_rules"
            )
            forbidden = validate_pattern_list(
                rule.get("forbidden_any", []), f"{location}.forbidden_any", errors, required=True
            )
            view = rule.get("view", "normalized")
            for path in paths:
                checked_documents += 1
                text = text_for(path, view)
                relative = path.relative_to(repo_root)
                for pattern in required:
                    required_anchors += 1
                    if re.search(pattern, text, flags=re.IGNORECASE | re.MULTILINE) is None:
                        errors.append(f"{guard_id}: {relative} missing required cross-document anchor {pattern!r}")
                for pattern in forbidden:
                    forbidden_anchors += 1
                    if re.search(pattern, text, flags=re.IGNORECASE | re.MULTILINE) is not None:
                        errors.append(f"{guard_id}: {relative} matched forbidden cross-document claim {pattern!r}")

    context_rules = contract.get("context_rules")
    if not isinstance(context_rules, list) or not context_rules:
        errors.append("cross_document_guards.context_rules must be a non-empty list")
        context_rules = []
    for index, rule in enumerate(context_rules):
        location = f"cross_document_guards.context_rules[{index}]"
        if not isinstance(rule, dict):
            errors.append(f"{location}: rule must be an object")
            continue
        guard_id = rule.get("guard_id")
        if not isinstance(guard_id, str) or not guard_id.strip():
            errors.append(f"{location}: guard_id must be a non-empty string")
            guard_id = location
        guard_ids.append(guard_id)
        paths = expand_rule_targets(rule, repo_root, location, errors)
        anchor_value = rule.get("anchor_pattern")
        if not isinstance(anchor_value, str) or not anchor_value:
            errors.append(f"{location}.anchor_pattern must be a non-empty regex")
            continue
        anchor = validate_regex(anchor_value, f"{location}.anchor_pattern", errors)
        required = validate_pattern_list(
            rule.get("required_all"), f"{location}.required_all", errors, required=True
        )
        skip_patterns = validate_pattern_list(
            rule.get("skip_if_any", []), f"{location}.skip_if_any", errors, required=False
        )
        window_chars = rule.get("window_chars")
        if not isinstance(window_chars, int) or not 50 <= window_chars <= 5000:
            errors.append(f"{location}.window_chars must be an integer from 50 to 5000")
            continue
        minimum_occurrences = rule.get("minimum_occurrences", 0)
        if not isinstance(minimum_occurrences, int) or minimum_occurrences < 0:
            errors.append(f"{location}.minimum_occurrences must be a non-negative integer")
            continue
        if anchor is None:
            continue
        rule_occurrences = 0
        for path in paths:
            text = text_for(path, rule.get("view", "normalized"))
            relative = path.relative_to(repo_root)
            for match in anchor.finditer(text):
                rule_occurrences += 1
                context_occurrences += 1
                context = text[max(0, match.start() - window_chars) : min(len(text), match.end() + window_chars)]
                if any(re.search(pattern, context, flags=re.IGNORECASE | re.MULTILINE) for pattern in skip_patterns):
                    continue
                for pattern in required:
                    required_anchors += 1
                    if re.search(pattern, context, flags=re.IGNORECASE | re.MULTILINE) is None:
                        errors.append(
                            f"{guard_id}: {relative} occurrence lacks nearby semantic anchor {pattern!r}"
                        )
        if rule_occurrences < minimum_occurrences:
            errors.append(
                f"{guard_id}: anchor_pattern matched {rule_occurrences}; minimum is {minimum_occurrences}"
            )

    duplicate_guard_ids = sorted({item for item in guard_ids if guard_ids.count(item) > 1})
    if duplicate_guard_ids:
        errors.append("duplicate cross-document guard IDs: " + ", ".join(duplicate_guard_ids))

    required_paths = contract.get("required_tracked_paths")
    executable_paths = contract.get("required_executable_paths")
    if not isinstance(required_paths, list) or not required_paths or not all(
        isinstance(item, str) and item.strip() for item in required_paths
    ):
        errors.append("cross_document_guards.required_tracked_paths must be a non-empty string list")
        required_paths = []
    if not isinstance(executable_paths, list) or not all(
        isinstance(item, str) and item.strip() for item in executable_paths
    ):
        errors.append("cross_document_guards.required_executable_paths must be a string list")
        executable_paths = []
    if not set(executable_paths).issubset(set(required_paths)):
        errors.append("required_executable_paths must be a subset of required_tracked_paths")

    tracked_paths = 0
    for relative in required_paths:
        candidate = Path(relative)
        location = f"cross_document_guards.required_tracked_paths[{relative}]"
        if candidate.is_absolute() or ".." in candidate.parts:
            errors.append(f"{location}: path must be repo-relative without '..'")
            continue
        path = repo_root / candidate
        if not path.is_file() or path.is_symlink():
            errors.append(f"{location}: shipped fixture must be a non-symlink regular file")
            continue
        result = subprocess.run(
            ["git", "ls-files", "--stage", "--", relative],
            cwd=repo_root,
            text=True,
            capture_output=True,
            check=False,
        )
        rows = [line.split(maxsplit=3) for line in result.stdout.splitlines() if line.strip()]
        stage_zero = [row for row in rows if len(row) == 4 and row[2] == "0" and row[3] == relative]
        if result.returncode != 0 or len(stage_zero) != 1:
            errors.append(f"{location}: path is not exactly one stage-0 Git index entry")
            continue
        mode = stage_zero[0][0]
        index_oid = stage_zero[0][1]
        if mode not in {"100644", "100755"}:
            errors.append(f"{location}: unsupported Git mode {mode}; regular file required")
            continue
        working_oid = subprocess.run(
            ["git", "hash-object", "--", relative],
            cwd=repo_root,
            text=True,
            capture_output=True,
            check=False,
        )
        if working_oid.returncode != 0 or working_oid.stdout.strip() != index_oid:
            errors.append(f"{location}: working fixture bytes differ from the stage-0 Git blob")
            continue
        tracked_paths += 1
        if relative in executable_paths and (mode != "100755" or not os.access(path, os.X_OK)):
            errors.append(f"{location}: documented entrypoint must be tracked executable mode 100755")

    summary = {
        "guard_contract": "scripts/p0_claim_registry.json",
        "public_source_files": len(public_paths),
        "guard_ids": len(guard_ids),
        "checked_documents": checked_documents,
        "required_anchors": required_anchors,
        "forbidden_anchors": forbidden_anchors,
        "context_occurrences": context_occurrences,
        "required_tracked_paths": len(required_paths),
        "tracked_regular_paths": tracked_paths,
        "errors": len(errors),
        "verdict": "PASS" if not errors else "FAIL",
    }
    return errors, summary


def validate_registry(registry: dict[str, Any], repo_root: Path) -> tuple[list[str], dict[str, Any]]:
    errors: list[str] = []
    claims = registry.get("claims")
    if registry.get("schema_version") != 1:
        errors.append("registry schema_version must be exactly 1")
    if not isinstance(claims, list):
        return [*errors, "registry claims must be a list"], {}

    locked_public_state = registry.get("locked_public_state")
    if not isinstance(locked_public_state, list):
        errors.append("locked_public_state must be a list")
        locked_public_state = []
    locked_surfaces = {
        item.get("surface")
        for item in locked_public_state
        if isinstance(item, dict)
        and item.get("status") == "external_action_required"
        and isinstance(item.get("locked_at"), str)
        and item["locked_at"].strip()
    }
    if locked_surfaces != EXPECTED_EXTERNAL_SURFACES:
        errors.append(
            "locked_public_state must mark exactly these surfaces external_action_required: "
            + ", ".join(sorted(EXPECTED_EXTERNAL_SURFACES))
        )

    claim_ids = [item.get("claim_id") for item in claims if isinstance(item, dict)]
    duplicate_ids = sorted({item for item in claim_ids if claim_ids.count(item) > 1})
    if duplicate_ids:
        errors.append(f"duplicate claim IDs: {', '.join(duplicate_ids)}")

    try:
        inventory_ids = load_inventory_p0(repo_root, registry["source_inventory"])
    except (KeyError, OSError, csv.Error) as exc:
        return [*errors, f"cannot load source inventory: {exc}"], {}

    registry_ids = {item for item in claim_ids if isinstance(item, str)}
    missing = sorted(inventory_ids - registry_ids)
    extra = sorted(registry_ids - inventory_ids)
    if missing:
        errors.append(f"registry missing inventory P0 IDs: {', '.join(missing)}")
    if extra:
        errors.append(f"registry has non-inventory P0 IDs: {', '.join(extra)}")

    checked_targets = 0
    required_anchors = 0
    forbidden_anchors = 0
    external_actions = 0
    document_cache: dict[Path, str] = {}

    for index, claim in enumerate(claims):
        if not isinstance(claim, dict):
            errors.append(f"claims[{index}] must be an object")
            continue
        claim_id = claim.get("claim_id", f"claims[{index}]")
        disposition = claim.get("disposition")
        checks = claim.get("checks")
        actions = claim.get("external_actions")
        if disposition not in ALLOWED_DISPOSITIONS:
            errors.append(f"{claim_id}: invalid disposition {disposition!r}")
        if not isinstance(actions, list) or not actions or not all(isinstance(x, str) and x.strip() for x in actions):
            errors.append(f"{claim_id}: external_actions must be a non-empty string list")
        else:
            external_actions += len(actions)

        if disposition == "external_action_required":
            if checks not in ([], None):
                errors.append(f"{claim_id}: external-only claim must not pretend to have local checks")
            continue
        if not isinstance(checks, list) or not checks:
            errors.append(f"{claim_id}: local-source disposition requires non-empty checks")
            continue

        for check_index, check in enumerate(checks):
            location = f"{claim_id}.checks[{check_index}]"
            if not isinstance(check, dict):
                errors.append(f"{location}: check must be an object")
                continue
            target = check.get("target")
            required = check.get("required_all")
            forbidden = check.get("forbidden_any")
            if not isinstance(target, str) or not target.strip():
                errors.append(f"{location}: target must be a non-empty path")
                continue
            if not isinstance(required, list) or not required:
                errors.append(f"{location}: required_all must be non-empty")
                required = []
            if not isinstance(forbidden, list) or not forbidden:
                errors.append(f"{location}: forbidden_any must be non-empty")
                forbidden = []

            try:
                path = resolve_repo_relative_file(repo_root, target, f"{location}.target")
            except (ValueError, FileNotFoundError) as exc:
                errors.append(str(exc))
                continue
            checked_targets += 1
            if path not in document_cache:
                document_cache[path] = normalize_document(path)
            text = document_cache[path]

            for pattern_index, pattern in enumerate(required):
                if not isinstance(pattern, str) or not pattern:
                    errors.append(f"{location}.required_all[{pattern_index}]: must be a non-empty regex")
                    continue
                required_anchors += 1
                regex = validate_regex(pattern, f"{location}.required_all[{pattern_index}]", errors)
                if regex is not None and regex.search(text) is None:
                    errors.append(f"{location}: missing required anchor {pattern!r}")

            for pattern_index, pattern in enumerate(forbidden):
                if not isinstance(pattern, str) or not pattern:
                    errors.append(f"{location}.forbidden_any[{pattern_index}]: must be a non-empty regex")
                    continue
                forbidden_anchors += 1
                regex = validate_regex(pattern, f"{location}.forbidden_any[{pattern_index}]", errors)
                if regex is not None and regex.search(text) is not None:
                    errors.append(f"{location}: residual forbidden claim matched {pattern!r}")

    cross_errors, cross_summary = validate_cross_document_guards(registry, repo_root)
    errors.extend(f"cross-document guard: {item}" for item in cross_errors)

    summary = {
        "inventory_p0": len(inventory_ids),
        "registry_p0": len(registry_ids),
        "local_claims": sum(
            1 for item in claims if item.get("disposition") == "local_source_fixed_external_publish_required"
        ),
        "external_only_claims": sum(1 for item in claims if item.get("disposition") == "external_action_required"),
        "checked_target_rules": checked_targets,
        "required_anchors": required_anchors,
        "forbidden_anchors": forbidden_anchors,
        "external_actions": external_actions,
        "external_public_surfaces": len(locked_surfaces),
        "unique_documents": len(document_cache),
        "cross_document_guard": cross_summary,
        "errors": len(errors),
        "verdict": "PASS" if not errors else "FAIL",
    }
    return errors, summary


def main() -> int:
    args = parse_args()
    registry_path = args.registry.resolve()
    repo_root = args.repo_root.resolve()
    try:
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(json.dumps({"verdict": "FAIL", "errors": [f"cannot load registry: {exc}"]}, ensure_ascii=False, indent=2))
        return 1

    if args.simulate_drop_claim:
        claims = registry.get("claims")
        if isinstance(claims, list):
            registry["claims"] = [item for item in claims if item.get("claim_id") != args.simulate_drop_claim]

    errors, summary = validate_registry(registry, repo_root)
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    if errors:
        print("ERRORS:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
