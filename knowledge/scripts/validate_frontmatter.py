#!/usr/bin/env python3
"""Validate YAML frontmatter of all .md files in the ISM KB.

Usage:
    python3 scripts/validate_frontmatter.py <kb_root>

Checks:
- Required fields present: id, name, description, status, last_verified,
  content_nature, doc_type, verified_scope, tags, canonical_paths
  (v0.5: content_nature renamed from source_type; backward compat warning emitted)
- Value constraints:
    status ∈ {active, archived, deprecated}
    content_nature ∈ {runtime-fact, frozen-decision, reference-summary,
                      reference, paper-derived, historical-note, postmortem}
    doc_type ∈ {reference, howto, explanation, tutorial}
    last_verified matches YYYY-MM-DD
- id is globally unique and starts with 'ism-kb-'
- canonical_paths list non-empty
"""

import argparse
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.stderr.write("ERROR: pyyaml not installed. Run: pip install pyyaml\n")
    sys.exit(2)


REQUIRED_FIELDS = {
    "id", "name", "description", "status", "last_verified",
    "content_nature", "doc_type", "verified_scope", "tags", "canonical_paths",
}
STATUS_VALUES = {"active", "archived", "deprecated"}
CONTENT_NATURE_VALUES = {
    "runtime-fact", "frozen-decision", "reference-summary",
    "reference", "paper-derived", "historical-note", "postmortem",
}
DOC_TYPE_VALUES = {"reference", "howto", "explanation", "tutorial"}
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
ID_RE = re.compile(r"^ism-kb-[0-9a-z\-]+$")

SKIP_FILES = {"README.md", "AGENT.md", "CHANGELOG.md"}


def extract_frontmatter(md_path):
    """Return frontmatter dict or None if not present."""
    text = md_path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        return None
    end = text.find("\n---\n", 4)
    if end == -1:
        return None
    block = text[4:end]
    return yaml.safe_load(block)


def validate_file(md_path, seen_ids, warnings):
    errors = []
    fm = extract_frontmatter(md_path)
    if fm is None:
        return [f"{md_path}: MISSING frontmatter block"]

    # v0.5 backward compat: accept source_type but warn
    if "source_type" in fm and "content_nature" not in fm:
        warnings.append(
            f"{md_path}: DEPRECATED field 'source_type' (rename to 'content_nature' per v0.5 schema)"
        )
        fm["content_nature"] = fm["source_type"]

    missing = REQUIRED_FIELDS - set(fm.keys())
    if missing:
        errors.append(f"{md_path}: MISSING fields {sorted(missing)}")

    if "id" in fm:
        if not isinstance(fm["id"], str) or not ID_RE.match(fm["id"]):
            errors.append(f"{md_path}: invalid id '{fm['id']}' (must match ism-kb-<group>-<topic>)")
        elif fm["id"] in seen_ids:
            errors.append(f"{md_path}: duplicate id '{fm['id']}' (also in {seen_ids[fm['id']]})")
        else:
            seen_ids[fm["id"]] = md_path

    if "status" in fm and fm["status"] not in STATUS_VALUES:
        errors.append(f"{md_path}: invalid status '{fm['status']}'")
    if "content_nature" in fm and fm["content_nature"] not in CONTENT_NATURE_VALUES:
        errors.append(f"{md_path}: invalid content_nature '{fm['content_nature']}'")
    if "doc_type" in fm and fm["doc_type"] not in DOC_TYPE_VALUES:
        errors.append(f"{md_path}: invalid doc_type '{fm['doc_type']}'")
    if "last_verified" in fm:
        lv = str(fm["last_verified"])
        if not DATE_RE.match(lv):
            errors.append(f"{md_path}: invalid last_verified '{lv}' (expected YYYY-MM-DD)")

    if "canonical_paths" in fm:
        cp = fm["canonical_paths"]
        if not isinstance(cp, list) or len(cp) == 0:
            errors.append(f"{md_path}: canonical_paths must be non-empty list")

    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path, help="KB root directory")
    args = parser.parse_args()

    root = args.root.resolve()
    if not root.is_dir():
        sys.stderr.write(f"ERROR: {root} is not a directory\n")
        sys.exit(2)

    seen_ids = {}
    all_errors = []
    warnings = []
    n_files = 0
    for md in sorted(root.rglob("*.md")):
        if md.name in SKIP_FILES and md.parent == root:
            continue
        n_files += 1
        errs = validate_file(md, seen_ids, warnings)
        all_errors.extend(errs)

    if warnings:
        print(f"WARN: {len(warnings)} deprecation warning(s):")
        for w in warnings:
            print(f"  - {w}")
        print()

    if all_errors:
        print(f"FAIL: {len(all_errors)} error(s) in {n_files} file(s):")
        for e in all_errors:
            print(f"  - {e}")
        sys.exit(1)
    print(f"PASS: {n_files} file(s) validated, {len(seen_ids)} unique ids")


if __name__ == "__main__":
    main()
