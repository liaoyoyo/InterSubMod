#!/usr/bin/env python3
"""Check Markdown files for broken relative links and image references.

Scans each given .md file for:
  - Image refs:  ![alt](path)
  - Plain links: [text](path)

Skips URLs (http/https/mailto/ftp/#anchor) and absolute paths. Resolves each
relative path against the containing .md's parent directory and reports any
that do not exist on disk.

Advisory-only: always exits 0 so PostToolUse hooks never block Edit/Write.
Findings are printed to stdout with [LinkCheck] prefix for hook visibility.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

# Matches both ![alt](path) and [text](path). Group 1 = is-image flag, 2 = label, 3 = target.
# Target captured non-greedily up to first unescaped ')'; strips optional "title".
LINK_PATTERN = re.compile(
    r'(!?)\[([^\]]*)\]\(\s*<?([^)\s<>]+)>?(?:\s+"[^"]*")?\s*\)'
)

SKIP_SCHEMES = ("http://", "https://", "mailto:", "ftp://", "ftps://", "tel:", "data:")


def iter_links(md_path: Path):
    """Yield (lineno, is_image, label, target) for each link in the markdown file."""
    try:
        text = md_path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError) as exc:
        print(f"[LinkCheck] WARN cannot read {md_path}: {exc}", file=sys.stderr)
        return
    for lineno, line in enumerate(text.splitlines(), start=1):
        # Cheap skip: lines in fenced code blocks still get scanned — acceptable false-positive
        # rate is lower than missing a real broken link inside a prose paragraph.
        for m in LINK_PATTERN.finditer(line):
            yield lineno, bool(m.group(1)), m.group(2), m.group(3)


def is_external(target: str) -> bool:
    target_lower = target.lower()
    if target.startswith("#"):
        return True
    return any(target_lower.startswith(s) for s in SKIP_SCHEMES)


def check_file(md_path: Path) -> list[str]:
    """Return list of human-readable broken-link messages."""
    findings: list[str] = []
    md_dir = md_path.parent.resolve()
    for lineno, is_image, label, target in iter_links(md_path):
        if is_external(target):
            continue
        # Drop fragment (#section) and query (?q=...) for filesystem resolution
        fs_target = target.split("#", 1)[0].split("?", 1)[0]
        if not fs_target:
            continue
        # Absolute paths: check as-is
        target_path = Path(fs_target)
        if not target_path.is_absolute():
            target_path = (md_dir / fs_target).resolve()
        if not target_path.exists():
            kind = "image" if is_image else "link"
            label_snippet = (label[:40] + "…") if len(label) > 40 else label
            findings.append(
                f"[LinkCheck] {md_path}:{lineno} broken {kind} → "
                f"'{target}' (resolved: {target_path}) [{label_snippet!r}]"
            )
    return findings


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: check_md_links.py <file.md> [<file.md> ...]", file=sys.stderr)
        return 0
    all_findings: list[str] = []
    for arg in argv[1:]:
        p = Path(arg)
        if not p.exists() or not p.is_file() or p.suffix.lower() != ".md":
            continue
        all_findings.extend(check_file(p))
    if all_findings:
        for msg in all_findings:
            print(msg)
        print(f"[LinkCheck] Found {len(all_findings)} broken reference(s) — advisory only.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
