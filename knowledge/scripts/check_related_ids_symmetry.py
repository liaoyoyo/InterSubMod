#!/usr/bin/env python3
"""Check bidirectional symmetry of related_ids in the ISM KB.

Rule: If A's frontmatter lists B in related_ids, then B's frontmatter must list A.
Exception: index files (00_index.md) are allowed to list children unidirectionally.

Usage:
    python3 scripts/check_related_ids_symmetry.py <kb_root>
"""

import argparse
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.stderr.write("ERROR: pyyaml not installed. Run: pip install pyyaml\n")
    sys.exit(2)


SKIP_FILES = {"README.md", "AGENT.md", "CHANGELOG.md"}


def extract_frontmatter(md_path):
    text = md_path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        return None
    end = text.find("\n---\n", 4)
    if end == -1:
        return None
    return yaml.safe_load(text[4:end])


def is_index(md_path):
    return md_path.name == "00_index.md"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()

    root = args.root.resolve()
    if not root.is_dir():
        sys.stderr.write(f"ERROR: {root} is not a directory\n")
        sys.exit(2)

    id_to_path = {}
    relations = {}

    for md in sorted(root.rglob("*.md")):
        if md.name in SKIP_FILES and md.parent == root:
            continue
        fm = extract_frontmatter(md)
        if not fm or "id" not in fm:
            continue
        doc_id = fm["id"]
        id_to_path[doc_id] = md
        relations[doc_id] = {
            "related": list(fm.get("related_ids", []) or []),
            "is_index": is_index(md),
        }

    errors = []
    for doc_id, info in relations.items():
        for rel_id in info["related"]:
            if rel_id not in relations:
                errors.append(
                    f"DANGLING: {doc_id} → {rel_id} (target doesn't exist)"
                )
                continue
            target_info = relations[rel_id]
            if doc_id not in target_info["related"]:
                if info["is_index"]:
                    continue
                errors.append(
                    f"ASYMMETRIC: {doc_id} → {rel_id}, but {rel_id} does not list {doc_id}"
                )

    if errors:
        print(f"FAIL: {len(errors)} issue(s):")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    print(f"PASS: {len(id_to_path)} document(s) checked, all related_ids symmetric")


if __name__ == "__main__":
    main()
