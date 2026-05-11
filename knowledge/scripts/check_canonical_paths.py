#!/usr/bin/env python3
"""Verify canonical_paths in frontmatter point to existing files.

Usage:
    python3 scripts/check_canonical_paths.py <kb_root>
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()

    root = args.root.resolve()
    if not root.is_dir():
        sys.stderr.write(f"ERROR: {root} is not a directory\n")
        sys.exit(2)

    errors = []
    n_files = 0
    for md in sorted(root.rglob("*.md")):
        if md.name in SKIP_FILES and md.parent == root:
            continue
        n_files += 1
        fm = extract_frontmatter(md)
        if not fm:
            continue
        cps = fm.get("canonical_paths", [])
        rel_md = md.relative_to(root)
        for cp in cps:
            target = root / cp
            if not target.exists():
                errors.append(f"BROKEN: {rel_md} → canonical_paths: {cp} (file not found)")
            elif target.resolve() != md.resolve():
                errors.append(
                    f"MISMATCH: {rel_md} has canonical_paths={cp} "
                    f"but file is actually at {rel_md}"
                )

    if errors:
        print(f"FAIL: {len(errors)} issue(s):")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    print(f"PASS: {n_files} file(s) checked, all canonical_paths resolve correctly")


if __name__ == "__main__":
    main()
