#!/usr/bin/env python3
"""Run reproducible Chromium QA for the standalone public documentation.

This runner deliberately produces no screenshots, PDFs, or browser cache in the
repository.  Print PDFs are rendered in memory only, and the JSON receipt is
written atomically through a temporary file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional
from urllib.parse import urlparse

from lxml import etree
from playwright.sync_api import Browser, Page, Route, sync_playwright


RUNNER_RELATIVE_PATH = Path(
    "research/20260813_public_docs_p0_correction/scripts/run_html_browser_qa.py"
)
DEFAULT_RECEIPT_RELATIVE_PATH = Path(
    "research/20260813_public_docs_p0_correction/html_browser_qa_receipt.json"
)
EXPECTED_HTML_COUNT = 21
EXPECTED_HTML_RELATIVE_PATHS = (
    Path("docs/explain/01_background-glossary.standalone.html"),
    Path("docs/explain/02_ism-core.standalone.html"),
    Path("docs/explain/03_methylation-read-filter.standalone.html"),
    Path("docs/explain/04_subclone-reconstruction-chr2-18M.standalone.html"),
    Path("docs/explain/05_subclone-correction-audit-chr2-18M.standalone.html"),
    Path("docs/explain/06_ism-subclone-pipeline-concept.standalone.html"),
    Path("docs/explain/07_subclone-judgment-workstation-chr2-18M.standalone.html"),
    Path("docs/explain/08_subclone-logic-chain-chr2-18M.standalone.html"),
    Path("docs/explain/09_three-stats-division-of-labor.standalone.html"),
    Path("docs/explain/10_ism-cpp-vs-chr2-subclone-capability.standalone.html"),
    Path("docs/explain/11_system-map-overview.standalone.html"),
    Path("docs/explain/12_intersubmod-io.standalone.html"),
    Path("docs/explain/13_longlineage-io.standalone.html"),
    Path("docs/explain/14_upstream-data.standalone.html"),
    Path("docs/explain/15_python-html-layer.standalone.html"),
    Path("docs/explain/16_how-to-run.standalone.html"),
    Path("docs/explain/index.standalone.html"),
    Path(
        "docs/handoff/20260813_完整研究資料與軟體交接_01/"
        "20260813_完整研究交接總覽_01.standalone.html"
    ),
    Path(
        "docs/handoff/20260805_系統交接與驗收_01/"
        "00_交接總覽與驗收.standalone.html"
    ),
    Path(
        "docs/handoff/20260805_系統交接與驗收_01/"
        "05_repo整理與可觀察性盤點.standalone.html"
    ),
    Path("docs/handoff/20260806_兩repo端到端串接可行性盤點_01.standalone.html"),
)
SOURCE_STATES = ("CLEAN_CHECKOUT", "DIRTY", "UNTRACKED")
CLAIM_CEILING = (
    "This receipt validates local source parsing, Chromium runtime safety, "
    "document-level responsive overflow, absence of attempted HTTP(S) runtime "
    "requests, and in-memory Chromium print rendering. It does not validate "
    "scientific claim correctness, accessibility conformance, non-Chromium "
    "browsers, physical printers, publication state, or release readiness."
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate the exact 21-page public/handoff standalone HTML allowlist in desktop, "
            "mobile, no-JavaScript, and print Chromium modes."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
        help="InterSubMod repository root (default: inferred from this script).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Receipt path (default: "
            "research/20260813_public_docs_p0_correction/"
            "html_browser_qa_receipt.json under --root)."
        ),
    )
    return parser.parse_args()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_output(root: Path, *args: str) -> Optional[str]:
    completed = subprocess.run(
        ["git", *args],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def git_path_state(root: Path, path: Path) -> str:
    """Return the repository-derived state for one existing source file."""
    try:
        relative_path = path.relative_to(root)
    except ValueError as exc:
        raise RuntimeError(f"Source is outside repository root: {path}") from exc

    completed = subprocess.run(
        [
            "git",
            "status",
            "--porcelain=v1",
            "--untracked-files=all",
            "--",
            str(relative_path),
        ],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or "git status failed"
        raise RuntimeError(
            f"Cannot derive Git state for {relative_path}: {detail}"
        )
    status_lines = [line for line in completed.stdout.splitlines() if line]
    if any(line.startswith("?? ") for line in status_lines):
        return "UNTRACKED"
    if status_lines:
        return "DIRTY"

    tracked = subprocess.run(
        ["git", "ls-files", "--error-unmatch", "--", str(relative_path)],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    return "CLEAN_CHECKOUT" if tracked.returncode == 0 else "UNTRACKED"


def source_record(root: Path, path: Path) -> Dict[str, str]:
    return {
        "path": str(path.relative_to(root)),
        "sha256": sha256_path(path),
        "git_state": git_path_state(root, path),
    }


def determine_source_state(source_records: Iterable[Dict[str, str]]) -> str:
    states = {record["git_state"] for record in source_records}
    unknown = states.difference(SOURCE_STATES)
    if unknown:
        raise RuntimeError(f"Unexpected source Git state(s): {sorted(unknown)}")
    if "UNTRACKED" in states:
        return "UNTRACKED"
    if "DIRTY" in states:
        return "DIRTY"
    return "CLEAN_CHECKOUT"


def resolve_html_allowlist(root: Path) -> List[Path]:
    expected = {root / relative for relative in EXPECTED_HTML_RELATIVE_PATHS}
    discovered = set((root / "docs" / "explain").glob("*.standalone.html"))
    discovered.update((root / "docs" / "handoff").rglob("*.standalone.html"))
    missing = sorted(expected.difference(discovered))
    unexpected = sorted(discovered.difference(expected))
    if missing or unexpected:
        details: List[str] = []
        if missing:
            details.append(
                "missing=" + ",".join(str(path.relative_to(root)) for path in missing)
            )
        if unexpected:
            details.append(
                "unexpected="
                + ",".join(str(path.relative_to(root)) for path in unexpected)
            )
        raise RuntimeError("Standalone HTML allowlist mismatch: " + "; ".join(details))
    return [root / relative for relative in EXPECTED_HTML_RELATIVE_PATHS]


def strict_xml_parser() -> etree.XMLParser:
    return etree.XMLParser(
        recover=False,
        no_network=True,
        resolve_entities=False,
        load_dtd=False,
    )


def parse_html_files(root: Path, html_files: Iterable[Path]) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for path in html_files:
        parser = etree.HTMLParser(recover=False, no_network=True)
        errors: List[str] = []
        try:
            etree.parse(str(path), parser)
        except (OSError, etree.XMLSyntaxError) as exc:
            errors.append(str(exc))
        errors.extend(
            str(entry)
            for entry in parser.error_log
            if entry.level_name in {"ERROR", "FATAL"}
        )
        results.append(
            {
                "path": str(path.relative_to(root)),
                "sha256": sha256_path(path),
                "errors": errors,
                "passed": not errors,
            }
        )
    return results


def parse_standalone_svg_files(
    root: Path, svg_files: Iterable[Path]
) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for path in svg_files:
        errors: List[str] = []
        try:
            etree.parse(str(path), strict_xml_parser())
        except (OSError, etree.XMLSyntaxError) as exc:
            errors.append(str(exc))
        results.append(
            {
                "path": str(path.relative_to(root)),
                "sha256": sha256_path(path),
                "errors": errors,
                "passed": not errors,
            }
        )
    return results


def parse_inline_svg_blocks(root: Path, html_files: Iterable[Path]) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for path in html_files:
        source = path.read_text(encoding="utf-8")
        blocks = re.findall(r"<svg\b.*?</svg>", source, flags=re.IGNORECASE | re.DOTALL)
        for index, block in enumerate(blocks, start=1):
            errors: List[str] = []
            try:
                etree.fromstring(block.encode("utf-8"), strict_xml_parser())
            except etree.XMLSyntaxError as exc:
                errors.append(str(exc))
            results.append(
                {
                    "path": f"{path.relative_to(root)}#inline-svg-{index}",
                    "sha256": hashlib.sha256(block.encode("utf-8")).hexdigest(),
                    "errors": errors,
                    "passed": not errors,
                }
            )
    return results


MODES: List[Dict[str, Any]] = [
    {
        "name": "desktop",
        "viewport": {"width": 1440, "height": 900},
        "javascript_enabled": True,
        "print_media": False,
    },
    {
        "name": "mobile",
        "viewport": {"width": 390, "height": 844},
        "javascript_enabled": True,
        "print_media": False,
    },
    {
        "name": "no_js",
        "viewport": {"width": 390, "height": 844},
        "javascript_enabled": False,
        "print_media": False,
    },
    {
        "name": "print",
        "viewport": {"width": 1440, "height": 900},
        "javascript_enabled": True,
        "print_media": True,
    },
]


def install_external_request_block(page: Page, external_requests: List[str]) -> None:
    def observe_request(request: Any) -> None:
        if urlparse(request.url).scheme in {"http", "https"}:
            external_requests.append(request.url)

    def block_external_route(route: Route) -> None:
        if urlparse(route.request.url).scheme in {"http", "https"}:
            route.abort("blockedbyclient")
        else:
            route.continue_()

    page.on("request", observe_request)
    page.route("**/*", block_external_route)


def run_browser_case(
    browser: Browser,
    root: Path,
    path: Path,
    mode: Dict[str, Any],
) -> Dict[str, Any]:
    context = browser.new_context(
        viewport=mode["viewport"],
        java_script_enabled=mode["javascript_enabled"],
    )
    page = context.new_page()
    page_errors: List[str] = []
    console_errors: List[str] = []
    external_requests: List[str] = []
    page.on("pageerror", lambda exc: page_errors.append(str(exc)))
    page.on(
        "console",
        lambda message: (
            console_errors.append(message.text) if message.type == "error" else None
        ),
    )
    install_external_request_block(page, external_requests)

    dimensions: Optional[Dict[str, int]] = None
    offenders: List[Dict[str, Any]] = []
    navigation_error: Optional[str] = None
    print_pdf_size_bytes: Optional[int] = None
    print_pdf_sha256: Optional[str] = None

    try:
        page.goto(path.as_uri(), wait_until="load", timeout=30_000)
        if mode["javascript_enabled"]:
            page.wait_for_load_state("networkidle", timeout=30_000)
        if mode["print_media"]:
            page.emulate_media(media="print")

        dimensions = page.evaluate(
            """() => ({
                clientWidth: document.documentElement.clientWidth,
                scrollWidth: document.documentElement.scrollWidth,
                bodyClientWidth: document.body.clientWidth,
                bodyScrollWidth: document.body.scrollWidth
            })"""
        )
        if (
            dimensions["scrollWidth"] > dimensions["clientWidth"]
            or dimensions["bodyScrollWidth"] > dimensions["clientWidth"]
        ):
            offenders = page.evaluate(
                """() => Array.from(document.querySelectorAll('*'))
                  .filter((element) => {
                    const rect = element.getBoundingClientRect();
                    return rect.right > document.documentElement.clientWidth + 0.5
                      || rect.left < -0.5;
                  })
                  .slice(0, 20)
                  .map((element) => {
                    const rect = element.getBoundingClientRect();
                    return {
                      tag: element.tagName,
                      id: element.id,
                      cls: typeof element.className === 'string' ? element.className : '',
                      left: Math.round(rect.left * 10) / 10,
                      right: Math.round(rect.right * 10) / 10,
                      scrollWidth: element.scrollWidth,
                      clientWidth: element.clientWidth
                    };
                  })"""
            )
        if mode["print_media"]:
            pdf = page.pdf(format="A4", print_background=True)
            print_pdf_size_bytes = len(pdf)
            print_pdf_sha256 = hashlib.sha256(pdf).hexdigest()
    except Exception as exc:  # Playwright exposes several runtime exception types.
        navigation_error = str(exc)
    finally:
        context.close()

    overflow_pass = bool(
        dimensions
        and dimensions["scrollWidth"] <= dimensions["clientWidth"]
        and dimensions["bodyScrollWidth"] <= dimensions["clientWidth"]
    )
    print_pass = bool(
        not mode["print_media"]
        or (
            print_pdf_size_bytes is not None
            and print_pdf_size_bytes > 1_000
            and print_pdf_sha256
        )
    )
    passed = bool(
        navigation_error is None
        and overflow_pass
        and not page_errors
        and not console_errors
        and not external_requests
        and print_pass
    )
    return {
        "path": str(path.relative_to(root)),
        "source_sha256": sha256_path(path),
        "mode": mode["name"],
        "viewport": mode["viewport"],
        "javascript_enabled": mode["javascript_enabled"],
        "print_media": mode["print_media"],
        "dimensions": dimensions,
        "document_overflow_pass": overflow_pass,
        "overflow_offenders": offenders,
        "page_errors": page_errors,
        "console_errors": console_errors,
        "external_runtime_requests": external_requests,
        "print_pdf_size_bytes": print_pdf_size_bytes,
        "print_pdf_sha256": print_pdf_sha256,
        "navigation_error": navigation_error,
        "passed": passed,
    }


def summarize_mode(mode: Dict[str, Any], cases: List[Dict[str, Any]]) -> Dict[str, Any]:
    selected = [case for case in cases if case["mode"] == mode["name"]]
    print_sizes = [
        case["print_pdf_size_bytes"]
        for case in selected
        if case["print_pdf_size_bytes"] is not None
    ]
    return {
        "mode": mode["name"],
        "viewport": mode["viewport"],
        "javascript_enabled": mode["javascript_enabled"],
        "print_media": mode["print_media"],
        "cases": len(selected),
        "passed": sum(1 for case in selected if case["passed"]),
        "failed": sum(1 for case in selected if not case["passed"]),
        "page_errors": sum(len(case["page_errors"]) for case in selected),
        "console_errors": sum(len(case["console_errors"]) for case in selected),
        "overflow_failures": sum(
            1 for case in selected if not case["document_overflow_pass"]
        ),
        "external_runtime_requests": sum(
            len(case["external_runtime_requests"]) for case in selected
        ),
        "print_pdf_rendered_in_memory": len(print_sizes),
        "minimum_print_pdf_size_bytes": min(print_sizes) if print_sizes else None,
        "maximum_print_pdf_size_bytes": max(print_sizes) if print_sizes else None,
    }


def atomic_write_json(output: Path, payload: Dict[str, Any]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, output)
        temporary_path = None
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    output = (
        args.output.resolve()
        if args.output is not None
        else root / DEFAULT_RECEIPT_RELATIVE_PATH
    )
    runner_path = Path(__file__).resolve()
    expected_runner = root / RUNNER_RELATIVE_PATH
    if runner_path != expected_runner:
        raise SystemExit(
            f"Runner must execute from {expected_runner}; actual path: {runner_path}"
        )

    try:
        html_files = resolve_html_allowlist(root)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    if len(html_files) != EXPECTED_HTML_COUNT:
        raise SystemExit("Internal error: HTML allowlist count is not exactly 21")
    standalone_svg_files = sorted((root / "docs").rglob("*.svg"))
    source_paths = [*html_files, *standalone_svg_files, runner_path]
    try:
        source_provenance = [source_record(root, path) for path in source_paths]
        source_state = determine_source_state(source_provenance)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    provenance_by_path = {
        record["path"]: record for record in source_provenance
    }

    html_parse = parse_html_files(root, html_files)
    standalone_svg_parse = parse_standalone_svg_files(root, standalone_svg_files)
    inline_svg_parse = parse_inline_svg_blocks(root, html_files)
    for record in html_parse + standalone_svg_parse:
        record["git_state"] = provenance_by_path[record["path"]]["git_state"]

    cases: List[Dict[str, Any]] = []
    browser_version: Optional[str] = None
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        browser_version = browser.version
        try:
            for mode in MODES:
                for path in html_files:
                    cases.append(run_browser_case(browser, root, path, mode))
        finally:
            browser.close()

    mode_results = [summarize_mode(mode, cases) for mode in MODES]
    case_counts = Counter(case["mode"] for case in cases)
    parse_pass = all(item["passed"] for item in html_parse)
    svg_pass = all(
        item["passed"] for item in standalone_svg_parse + inline_svg_parse
    )
    browser_pass = len(cases) == EXPECTED_HTML_COUNT * len(MODES) and all(
        case["passed"] for case in cases
    )
    try:
        final_source_provenance = [source_record(root, path) for path in source_paths]
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    sources_stable = source_provenance == final_source_provenance
    all_pass = parse_pass and svg_pass and browser_pass and sources_stable
    verified_at = datetime.now().astimezone().isoformat(timespec="seconds")
    source_state_counts = Counter(
        record["git_state"] for record in source_provenance
    )

    receipt: Dict[str, Any] = {
        "schema_version": 3,
        "verified_at": verified_at,
        "task_type": ["B_COMPREHENSIVE_VALIDATION", "D_EXTERNAL_HANDOFF"],
        "scope": "LOCAL_SOURCE_HTML_BROWSER_QA",
        "worktree": str(root),
        "tested_git_head": git_output(root, "rev-parse", "HEAD"),
        "source_state": source_state,
        "source_state_counts": dict(sorted(source_state_counts.items())),
        "source_provenance": {
            "html": source_provenance[: len(html_files)],
            "standalone_svg": source_provenance[
                len(html_files) : len(html_files) + len(standalone_svg_files)
            ],
            "runner": source_provenance[-1],
            "stable_during_run": sources_stable,
        },
        "runner": {
            "path": str(runner_path.relative_to(root)),
            "sha256": sha256_path(runner_path),
            "git_state": provenance_by_path[str(runner_path.relative_to(root))][
                "git_state"
            ],
            "python": sys.version.split()[0],
        },
        "command": [
            "python3",
            str(runner_path.relative_to(root)),
            "--root",
            str(root),
            "--output",
            str(output),
        ],
        "browser": {
            "engine": "Chromium",
            "version": browser_version,
            "headless": True,
        },
        "tested_files": html_parse,
        "parsing": {
            "html": {
                "files": len(html_parse),
                "errors": sum(len(item["errors"]) for item in html_parse),
                "verdict": "PASS" if parse_pass else "FAIL",
            },
            "standalone_svg_xml": {
                "files": len(standalone_svg_parse),
                "errors": sum(
                    len(item["errors"]) for item in standalone_svg_parse
                ),
                "verdict": (
                    "PASS"
                    if all(item["passed"] for item in standalone_svg_parse)
                    else "FAIL"
                ),
                "tested_files": standalone_svg_parse,
            },
            "inline_svg_xml": {
                "blocks": len(inline_svg_parse),
                "errors": sum(len(item["errors"]) for item in inline_svg_parse),
                "verdict": (
                    "PASS"
                    if all(item["passed"] for item in inline_svg_parse)
                    else "FAIL"
                ),
                "tested_blocks": inline_svg_parse,
            },
        },
        "test_matrix": mode_results,
        "browser_cases": cases,
        "checks": {
            "browser_cases": {
                "expected": EXPECTED_HTML_COUNT * len(MODES),
                "total": len(cases),
                "per_mode": dict(sorted(case_counts.items())),
                "passed": sum(1 for case in cases if case["passed"]),
                "failed": sum(1 for case in cases if not case["passed"]),
            },
            "document_overflow": {
                "rule": (
                    "document.documentElement.scrollWidth <= clientWidth and "
                    "document.body.scrollWidth <= document.documentElement.clientWidth"
                ),
                "failed": sum(
                    1 for case in cases if not case["document_overflow_pass"]
                ),
                "component_horizontal_scroll_allowed": True,
            },
            "page_errors": sum(len(case["page_errors"]) for case in cases),
            "console_errors": sum(len(case["console_errors"]) for case in cases),
            "external_runtime_requests": {
                "attempted": sum(
                    len(case["external_runtime_requests"]) for case in cases
                ),
                "policy": (
                    "HTTP(S) requests are observed, blocked before network access, "
                    "and fail the case."
                ),
            },
            "print_pdf": {
                "rendered_in_memory": sum(
                    1 for case in cases if case["print_pdf_size_bytes"] is not None
                ),
                "persisted_files": 0,
                "verdict": (
                    "PASS"
                    if next(
                        result for result in mode_results if result["mode"] == "print"
                    )["failed"]
                    == 0
                    else "FAIL"
                ),
            },
            "screenshots_written": 0,
            "repository_browser_cache_written": 0,
            "source_provenance_stable_during_run": sources_stable,
        },
        "claim_ceiling": CLAIM_CEILING,
        "verdict": "PASS_FOR_LOCAL_SOURCE_QA" if all_pass else "FAIL_CLOSED",
        "publication_status": (
            "BLOCKED_DEFAULT_BRANCH_WIKI_PAGES_NOT_PUBLISHED_OR_REFETCHED__"
            "ABOUT_C108_HAS_SEPARATE_LIVE_RECEIPT"
        ),
        "release_status": "BLOCKED",
        "known_limits": [
            "Only Chromium is exercised; Firefox and WebKit are outside this receipt.",
            "Print validation covers in-memory Chromium PDF bytes, not a physical printer or every platform print driver.",
            "No screenshots or manual visual-review artifacts are produced by this runner.",
            "Scientific claim correctness remains governed by the claim registry and frozen authority bundle.",
            "Accessibility conformance and assistive-technology testing are outside this receipt.",
            "Default branch, GitHub Wiki, Pages, and public URLs are not published or re-fetched by this local runner.",
        ],
    }
    atomic_write_json(output, receipt)

    summary = {
        "verdict": receipt["verdict"],
        "html_files": len(html_parse),
        "standalone_svg_files": len(standalone_svg_parse),
        "inline_svg_blocks": len(inline_svg_parse),
        "browser_cases": len(cases),
        "browser_passed": receipt["checks"]["browser_cases"]["passed"],
        "browser_failed": receipt["checks"]["browser_cases"]["failed"],
        "overflow_failed": receipt["checks"]["document_overflow"]["failed"],
        "page_errors": receipt["checks"]["page_errors"],
        "console_errors": receipt["checks"]["console_errors"],
        "external_runtime_requests": receipt["checks"]["external_runtime_requests"][
            "attempted"
        ],
        "receipt": str(output),
    }
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
