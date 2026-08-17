#!/usr/bin/env python3
"""Structural, link, and inline-SVG accessibility QA for InterSubMod Pages."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from urllib.parse import unquote, urlsplit

from bs4 import BeautifulSoup


EXPECTED_FILES = [
    "index.standalone.html",
    "01_background-glossary.standalone.html",
    "02_ism-core.standalone.html",
    "03_methylation-read-filter.standalone.html",
    "04_subclone-reconstruction-chr2-18M.standalone.html",
    "05_subclone-correction-audit-chr2-18M.standalone.html",
    "06_ism-subclone-pipeline-concept.standalone.html",
    "07_subclone-judgment-workstation-chr2-18M.standalone.html",
    "08_subclone-logic-chain-chr2-18M.standalone.html",
    "09_three-stats-division-of-labor.standalone.html",
    "10_ism-cpp-vs-chr2-subclone-capability.standalone.html",
    "11_system-map-overview.standalone.html",
    "12_intersubmod-io.standalone.html",
    "13_longlineage-io.standalone.html",
    "14_upstream-data.standalone.html",
    "15_python-html-layer.standalone.html",
    "16_how-to-run.standalone.html",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pages-dir", type=Path, default=Path("docs/explain"))
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--simulate-drop-svg-desc",
        action="store_true",
        help="Negative-control probe: treat the first SVG description as absent.",
    )
    return parser.parse_args()


def load_soup(path: Path) -> BeautifulSoup:
    return BeautifulSoup(path.read_text(encoding="utf-8"), "html.parser")


def local_link_error(source: Path, href: str, pages_dir: Path) -> str | None:
    if not href or href.startswith(("http://", "https://", "mailto:", "tel:", "data:")):
        return None
    if href.startswith("javascript:"):
        return f"unsafe javascript href: {href}"
    parts = urlsplit(href)
    path_part = unquote(parts.path)
    fragment = unquote(parts.fragment)
    if not path_part:
        target = source
    elif path_part.startswith("/"):
        # Repository-root links cannot be resolved reliably from the Pages subdirectory.
        return None
    else:
        target = (source.parent / path_part).resolve()
        if pages_dir.resolve() not in target.parents and target != pages_dir.resolve():
            return None
    if not target.exists():
        return f"missing local target: {href}"
    if fragment and target.is_file() and target.suffix.lower() in {".html", ".htm"}:
        target_soup = load_soup(target)
        if target_soup.find(id=fragment) is None:
            return f"missing fragment #{fragment} in {target.name}"
    return None


def main() -> int:
    args = parse_args()
    pages_dir_display = args.pages_dir.as_posix()
    pages_dir = args.pages_dir.resolve()
    errors: list[str] = []
    warnings: list[str] = []
    page_rows: list[dict[str, object]] = []
    total_svg = 0
    simulated = False

    actual_names = sorted(path.name for path in pages_dir.glob("*.standalone.html"))
    if set(actual_names) != set(EXPECTED_FILES):
        errors.append(
            f"page inventory mismatch: missing={sorted(set(EXPECTED_FILES)-set(actual_names))}, "
            f"extra={sorted(set(actual_names)-set(EXPECTED_FILES))}"
        )

    for name in EXPECTED_FILES:
        path = pages_dir / name
        row_errors: list[str] = []
        row_warnings: list[str] = []
        if not path.is_file():
            row_errors.append("file missing")
            page_rows.append(
                {"page": name, "svg": 0, "errors": row_errors, "warnings": row_warnings}
            )
            continue
        text = path.read_text(encoding="utf-8")
        soup = BeautifulSoup(text, "html.parser")
        if "<!doctype html" not in text.lower():
            row_warnings.append("legacy HTML fragment: missing doctype")
        if soup.html is None or not soup.html.get("lang"):
            row_warnings.append("legacy HTML fragment: missing explicit html lang")
        if soup.title is None or not soup.title.get_text(strip=True):
            row_errors.append("missing non-empty title")
        if soup.find("h1") is None:
            row_errors.append("missing h1")

        ids = [tag.get("id") for tag in soup.find_all(attrs={"id": True})]
        duplicates = sorted({item for item in ids if ids.count(item) > 1})
        if duplicates:
            row_errors.append(f"duplicate ids: {duplicates}")

        svgs = soup.find_all("svg")
        total_svg += len(svgs)
        for index, svg in enumerate(svgs, start=1):
            label = f"svg[{index}]"
            title = svg.find("title")
            desc = svg.find("desc")
            if args.simulate_drop_svg_desc and not simulated:
                desc = None
                simulated = True
            if svg.get("role") != "img":
                row_errors.append(f"{label}: role must be img")
            if title is None or not title.get_text(" ", strip=True):
                row_errors.append(f"{label}: missing non-empty title")
            if desc is None or not desc.get_text(" ", strip=True):
                row_errors.append(f"{label}: missing non-empty desc")
            title_id = title.get("id") if title else None
            desc_id = desc.get("id") if desc else None
            aria_tokens = (svg.get("aria-labelledby") or "").split()
            if not title_id or not desc_id:
                row_errors.append(f"{label}: title/desc require ids")
            elif title_id not in aria_tokens or desc_id not in aria_tokens:
                row_errors.append(f"{label}: aria-labelledby must reference title and desc ids")
            # BeautifulSoup's html parser normalizes viewBox to lowercase.
            if not (svg.get("viewBox") or svg.get("viewbox")):
                row_errors.append(f"{label}: missing viewBox")

        for anchor in soup.find_all("a", href=True):
            error = local_link_error(path, anchor.get("href", ""), pages_dir)
            if error:
                row_errors.append(error)

        errors.extend(f"{name}: {error}" for error in row_errors)
        warnings.extend(f"{name}: {warning}" for warning in row_warnings)
        page_rows.append(
            {"page": name, "svg": len(svgs), "errors": row_errors, "warnings": row_warnings}
        )

    if total_svg != 37:
        errors.append(f"inline SVG total must be 37, observed {total_svg}")

    receipt = {
        "pages_dir": pages_dir_display,
        "expected_pages": 17,
        "observed_pages": len(actual_names),
        "expected_inline_svg": 37,
        "observed_inline_svg": total_svg,
        "pages": page_rows,
        "errors": errors,
        "warnings": warnings,
        "verdict": "PASS" if not errors else "FAIL",
    }
    rendered = json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
