#!/usr/bin/env python3
"""Fail-closed validator for the 2026-08-13 public-document P0 corrections."""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
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
    return re.sub(r"\s+", " ", without_tags).strip()


def load_inventory_p0(repo_root: Path, source: str) -> set[str]:
    path = repo_root / source
    if not path.is_file():
        raise FileNotFoundError(f"source inventory does not exist: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {row["claim_id"] for row in rows if row.get("priority") == "P0"}


def validate_regex(pattern: str, where: str, errors: list[str]) -> re.Pattern[str] | None:
    try:
        return re.compile(pattern, flags=re.IGNORECASE | re.MULTILINE)
    except re.error as exc:
        errors.append(f"{where}: invalid regex {pattern!r}: {exc}")
        return None


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

            path = repo_root / target
            if not path.is_file():
                errors.append(f"{location}: target does not exist: {path}")
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
