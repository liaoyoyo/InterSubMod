#!/usr/bin/env python3
"""List stale KB files (last_verified > 90d) or update a specific file's date.

Usage:
    # List stale files
    python3 scripts/refresh_last_verified.py <kb_root> --list-stale

    # Update single file's last_verified to today
    python3 scripts/refresh_last_verified.py <kb_root> --file <path.md>
"""

import argparse
import re
import sys
from datetime import date, datetime
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.stderr.write("ERROR: pyyaml not installed. Run: pip install pyyaml\n")
    sys.exit(2)


SKIP_FILES = {"README.md", "AGENT.md", "CHANGELOG.md"}
STALE_DAYS = 90
CRITICAL_DAYS = 180

DATE_LINE_RE = re.compile(r"^(last_verified:\s*)(.+)$", re.MULTILINE)


def extract_frontmatter(md_path):
    text = md_path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        return None, None, None
    end = text.find("\n---\n", 4)
    if end == -1:
        return None, None, None
    block = text[4:end]
    return yaml.safe_load(block), text, end


def classify(days):
    if days < STALE_DAYS:
        return "verified"
    if days < CRITICAL_DAYS:
        return "needs-recheck"
    return "stale"


def list_stale(root: Path):
    today = date.today()
    stats = {"verified": 0, "needs-recheck": 0, "stale": 0}
    entries = []
    for md in sorted(root.rglob("*.md")):
        if md.name in SKIP_FILES and md.parent == root:
            continue
        fm, _, _ = extract_frontmatter(md)
        if not fm or "last_verified" not in fm:
            continue
        lv = str(fm["last_verified"])
        try:
            dt = datetime.strptime(lv, "%Y-%m-%d").date()
        except ValueError:
            continue
        days = (today - dt).days
        cls = classify(days)
        stats[cls] += 1
        rel = md.relative_to(root)
        entries.append((days, cls, rel, lv))

    entries.sort(reverse=True)
    print(f"=== Freshness Report ({today}) ===")
    for cnt_key in ("verified", "needs-recheck", "stale"):
        print(f"  {cnt_key}: {stats[cnt_key]}")
    print()
    print("--- Files needing attention ---")
    for days, cls, rel, lv in entries:
        if cls == "verified":
            continue
        flag = "⚠" if cls == "needs-recheck" else "🔴"
        print(f"  {flag} {cls} ({days}d): {rel} (last_verified: {lv})")


def refresh_file(root: Path, rel_path: str):
    md = root / rel_path
    if not md.exists():
        sys.stderr.write(f"ERROR: {md} not found\n")
        sys.exit(2)
    text = md.read_text(encoding="utf-8")
    today_str = date.today().isoformat()
    new_text, count = DATE_LINE_RE.subn(
        lambda m: f"{m.group(1)}{today_str}",
        text,
        count=1,
    )
    if count == 0:
        sys.stderr.write(f"ERROR: no last_verified line found in {md}\n")
        sys.exit(2)
    md.write_text(new_text, encoding="utf-8")
    print(f"Updated {rel_path}: last_verified → {today_str}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--list-stale", action="store_true")
    parser.add_argument("--file", type=str, help="Relative path under root to refresh")
    args = parser.parse_args()

    root = args.root.resolve()
    if not root.is_dir():
        sys.stderr.write(f"ERROR: {root} is not a directory\n")
        sys.exit(2)

    if args.list_stale:
        list_stale(root)
    elif args.file:
        refresh_file(root, args.file)
    else:
        list_stale(root)


if __name__ == "__main__":
    main()
