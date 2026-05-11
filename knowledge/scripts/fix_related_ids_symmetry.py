#!/usr/bin/env python3
"""Auto-fix bidirectional symmetry by adding reverse related_ids.

Rule: If A lists B, ensure B lists A (unless A is an index file).

Usage:
    python3 scripts/fix_related_ids_symmetry.py <kb_root>
"""

import argparse
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.stderr.write("ERROR: pyyaml not installed\n")
    sys.exit(2)


SKIP_FILES = {"README.md", "AGENT.md", "CHANGELOG.md"}
FM_RE = re.compile(r"^---\n(.*?)\n---\n", re.DOTALL)
RELATED_BLOCK_RE = re.compile(
    r"(related_ids:\s*(?:\n  -\s*.+)+)",
    re.MULTILINE,
)
RELATED_EMPTY_RE = re.compile(r"related_ids:\s*\[\s*\]", re.MULTILINE)


def extract_frontmatter(text):
    m = FM_RE.match(text)
    if not m:
        return None, None
    try:
        fm = yaml.safe_load(m.group(1))
    except yaml.YAMLError:
        return None, None
    return fm, m.end()


def insert_related_id(text, new_id):
    """Insert new_id into the related_ids block, keep YAML sorted as-written."""
    m = RELATED_BLOCK_RE.search(text)
    if m:
        block = m.group(1)
        if new_id in block:
            return text
        indent_line = f"\n  - {new_id}"
        new_block = block + indent_line
        return text[: m.start(1)] + new_block + text[m.end(1) :]
    m2 = RELATED_EMPTY_RE.search(text)
    if m2:
        return text[: m2.start()] + f"related_ids:\n  - {new_id}" + text[m2.end() :]
    fm_end = FM_RE.match(text)
    if fm_end:
        insert_at = fm_end.end(1)
        return (
            text[:insert_at]
            + f"\nrelated_ids:\n  - {new_id}"
            + text[insert_at:]
        )
    return text


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    id_to_path = {}
    all_files = []
    for md in sorted(root.rglob("*.md")):
        if md.name in SKIP_FILES and md.parent == root:
            continue
        text = md.read_text(encoding="utf-8")
        fm, _ = extract_frontmatter(text)
        if not fm or "id" not in fm:
            continue
        all_files.append((md, fm, text))
        id_to_path[fm["id"]] = md

    fixes = []
    for md, fm, _ in all_files:
        doc_id = fm["id"]
        if md.name == "00_index.md":
            continue
        for rel_id in fm.get("related_ids", []) or []:
            if rel_id not in id_to_path:
                continue
            target = id_to_path[rel_id]
            target_text = target.read_text(encoding="utf-8")
            target_fm, _ = extract_frontmatter(target_text)
            if not target_fm:
                continue
            target_related = target_fm.get("related_ids", []) or []
            if doc_id not in target_related:
                fixes.append((target, doc_id))

    print(f"Found {len(fixes)} missing reverse reference(s)")
    if args.dry_run:
        for target, missing_id in fixes:
            print(f"  DRY: {target.relative_to(root)} needs += {missing_id}")
        return

    touched = {}
    for target, missing_id in fixes:
        text = target.read_text(encoding="utf-8")
        new_text = insert_related_id(text, missing_id)
        if new_text != text:
            target.write_text(new_text, encoding="utf-8")
            touched.setdefault(target, []).append(missing_id)

    for path, added in touched.items():
        print(f"  FIXED: {path.relative_to(root)} += {added}")
    print(f"Total files modified: {len(touched)}")


if __name__ == "__main__":
    main()
