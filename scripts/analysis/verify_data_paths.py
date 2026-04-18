#!/usr/bin/env python3
"""Verify data path integrity across InterSubMod project.

Checks:
  1. master_run_manifest.tsv: canonical_run_path and archive_run_path existence
  2. docs/reports/**/*.md: absolute paths (/big7_disk/, /big8_disk/, /bip8_disk/) existence
  3. observation_common.py MASTER_DATASET_PATH existence
  4. Output workspace directories existence

Outputs a TSV report to stdout (redirect to file as needed).

Usage:
    python3 scripts/analysis/verify_data_paths.py [--output path.tsv]
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
MANIFEST_PATH = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/master_run_manifest.tsv"
)
MASTER_DATASET_PATH = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/"
    "observation_workspaces/20260327_loh_round1_cross_sample_audit/"
    "all_region_rows.tsv.gz"
)
OBSERVATION_WS_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces"
)
DOCS_REPORTS_DIR = PROJECT_ROOT / "docs" / "reports"

# Regex to match absolute paths in markdown files
ABS_PATH_RE = re.compile(
    r"(/big[78]_disk/[^\s\)\"'>`\u3000\uff08\uff09\uff1a]+|"
    r"/bip8_disk/[^\s\)\"'>`\u3000\uff08\uff09\uff1a]+)"
)

# Characters to strip from end of matched paths
TRAIL_CHARS = "`）：。，、；*]|"


def classify_path(p: Path) -> str:
    """Classify a path as file, dir, symlink, or missing."""
    if not p.exists():
        if p.is_symlink():
            return "broken_symlink"
        return "missing"
    if p.is_symlink():
        return "symlink"
    if p.is_dir():
        return "dir"
    return "file"


def check_manifest(rows: list[dict]) -> None:
    """Check paths from master_run_manifest.tsv."""
    for col in ("canonical_run_path", "archive_run_path"):
        for row in rows:
            path_str = row.get(col, "").strip()
            if not path_str:
                continue
            p = Path(path_str)
            ptype = classify_path(p)
            yield {
                "path": path_str,
                "source_file": str(MANIFEST_PATH),
                "source_column": col,
                "exists": ptype != "missing" and ptype != "broken_symlink",
                "type": ptype,
                "note": f"sample={row.get('sample', '?')}, run_id={row.get('run_id', '?')}",
            }


def check_report_paths() -> None:
    """Scan docs/reports/**/*.md for absolute paths and verify existence."""
    seen = set()
    for md_file in sorted(DOCS_REPORTS_DIR.rglob("*.md")):
        try:
            content = md_file.read_text(encoding="utf-8", errors="replace")
        except (OSError, UnicodeDecodeError):
            continue
        for match in ABS_PATH_RE.finditer(content):
            raw = match.group(1).rstrip(".,;:)" + TRAIL_CHARS)
            # Strip markdown link syntax leakage: path](other_path
            if "](" in raw:
                raw = raw.split("](")[0]
            # Skip paths that are clearly template placeholders
            if "..." in raw or "{" in raw or len(raw) < 15:
                continue
            key = (str(md_file), raw)
            if key in seen:
                continue
            seen.add(key)
            p = Path(raw)
            ptype = classify_path(p)
            yield {
                "path": raw,
                "source_file": str(md_file.relative_to(PROJECT_ROOT)),
                "source_column": "markdown_content",
                "exists": ptype != "missing" and ptype != "broken_symlink",
                "type": ptype,
                "note": "",
            }


def check_master_dataset() -> None:
    """Verify the master dataset file exists."""
    ptype = classify_path(MASTER_DATASET_PATH)
    yield {
        "path": str(MASTER_DATASET_PATH),
        "source_file": "scripts/analysis/observation_common.py:41",
        "source_column": "MASTER_DATASET_PATH",
        "exists": ptype != "missing" and ptype != "broken_symlink",
        "type": ptype,
        "note": "master dataset for all O-series observations",
    }


def check_workspace_readmes() -> None:
    """Check which observation workspaces have README.md."""
    if not OBSERVATION_WS_ROOT.is_dir():
        return
    for d in sorted(OBSERVATION_WS_ROOT.iterdir()):
        if not d.is_dir() or d.name.startswith("."):
            continue
        readme = d / "README.md"
        has_readme = readme.exists()
        yield {
            "path": str(d),
            "source_file": "observation_workspaces/",
            "source_column": "workspace_dir",
            "exists": True,
            "type": "dir",
            "note": f"README={'yes' if has_readme else 'MISSING'}",
        }


def main():
    parser = argparse.ArgumentParser(description="Verify InterSubMod data paths")
    parser.add_argument("--output", "-o", help="Output TSV file path")
    args = parser.parse_args()

    out_f = open(args.output, "w", newline="") if args.output else sys.stdout
    writer = csv.DictWriter(
        out_f,
        fieldnames=["path", "source_file", "source_column", "exists", "type", "note"],
        delimiter="\t",
    )
    writer.writeheader()

    counts = {"total": 0, "ok": 0, "broken": 0}

    # 1. Manifest paths
    if MANIFEST_PATH.exists():
        with open(MANIFEST_PATH) as f:
            reader = csv.DictReader(f, delimiter="\t")
            manifest_rows = list(reader)
        for row in check_manifest(manifest_rows):
            writer.writerow(row)
            counts["total"] += 1
            counts["ok" if row["exists"] else "broken"] += 1

    # 2. Master dataset
    for row in check_master_dataset():
        writer.writerow(row)
        counts["total"] += 1
        counts["ok" if row["exists"] else "broken"] += 1

    # 3. Report paths
    for row in check_report_paths():
        writer.writerow(row)
        counts["total"] += 1
        counts["ok" if row["exists"] else "broken"] += 1

    # 4. Workspace READMEs
    for row in check_workspace_readmes():
        writer.writerow(row)
        counts["total"] += 1
        counts["ok" if row["exists"] else "broken"] += 1

    if out_f is not sys.stdout:
        out_f.close()

    # Summary to stderr
    print(f"\n=== Path Verification Summary ===", file=sys.stderr)
    print(f"Total checked: {counts['total']}", file=sys.stderr)
    print(f"OK: {counts['ok']}", file=sys.stderr)
    print(f"Broken: {counts['broken']}", file=sys.stderr)
    if counts["broken"] > 0:
        print(f"⚠ {counts['broken']} broken paths found!", file=sys.stderr)
    else:
        print("✓ All paths valid", file=sys.stderr)

    return 1 if counts["broken"] > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
