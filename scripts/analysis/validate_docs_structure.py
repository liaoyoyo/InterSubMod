#!/usr/bin/env python3
"""Validate docs naming, metadata, and obvious broken relative links."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"

DATE_NAME = re.compile(r"^\d{8}_.+_\d{2}\.md$")
FIXED = {"README.md", "CURRENT_FOCUS.md", "INDEX.md"}
LINK_RE = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
FENCE_RE = re.compile(r"```.*?```", re.DOTALL)


def is_archive_deep(p: Path) -> bool:
    rel = p.relative_to(DOCS)
    return str(rel).startswith("archive/deep/")


def check_names(md_files: list[Path]) -> list[str]:
    issues: list[str] = []
    for f in md_files:
        if is_archive_deep(f):
            continue
        name = f.name
        if name in FIXED:
            continue
        if f.parent.name in {"architecture", "modules"} and name.endswith(".md") and "_" in name:
            # architecture 長期文件允許 snake_case
            continue
        if not DATE_NAME.match(name):
            issues.append(f"[NAME] {f.relative_to(ROOT)}")
    return issues


def check_metadata(md_files: list[Path]) -> list[str]:
    issues: list[str] = []
    for f in md_files:
        if is_archive_deep(f):
            continue
        txt = f.read_text(errors="ignore").lstrip()
        if not txt.startswith("<!--"):
            issues.append(f"[META] {f.relative_to(ROOT)}")
    return issues


def check_links(md_files: list[Path]) -> list[str]:
    issues: list[str] = []
    for md in md_files:
        txt = md.read_text(errors="ignore")
        # Skip fenced code blocks to avoid false positives (e.g. C++ lambda `[&](...)`).
        txt = FENCE_RE.sub("", txt)
        for m in LINK_RE.finditer(txt):
            link = m.group(1).strip()
            if not link or link.startswith(("http://", "https://", "#", "mailto:")):
                continue
            if link.startswith("file://"):
                issues.append(f"[LINK:file://] {md.relative_to(ROOT)} -> {link}")
                continue
            target_link = link.split("#", 1)[0]
            if target_link.startswith("/"):
                target = Path(target_link)
            else:
                target = (md.parent / target_link).resolve()
            if not target.exists():
                issues.append(f"[LINK:missing] {md.relative_to(ROOT)} -> {link}")
    return issues


def main() -> int:
    md_files = sorted(DOCS.rglob("*.md"))
    name_issues = check_names(md_files)
    meta_issues = check_metadata(md_files)
    link_issues = check_links(md_files)

    print(f"Total markdown files: {len(md_files)}")
    print(f"Name issues: {len(name_issues)}")
    print(f"Metadata issues: {len(meta_issues)}")
    print(f"Link issues: {len(link_issues)}")

    for item in (name_issues + meta_issues + link_issues)[:200]:
        print(item)

    return 1 if (name_issues or meta_issues or link_issues) else 0


if __name__ == "__main__":
    raise SystemExit(main())
