#!/usr/bin/env python3
"""Runtime UI smoke tests for the layered and archived topology workstations.

The validator opens the generated, self-contained HTML files directly with
``file://`` URLs.  It writes its machine-readable report to stdout and only
writes PNG files when ``--screenshots-dir`` is explicitly provided.

Exit codes:
    0: every validation check passed
    1: one or more UI validation checks failed
    2: validator configuration, dependency, or browser startup failed
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple
from urllib.parse import unquote, urljoin, urlparse

try:
    from playwright.sync_api import Browser, Locator, Page, Playwright, sync_playwright
except ImportError as exc:  # Reported as structured JSON by main().
    Browser = Locator = Page = Playwright = Any  # type: ignore[misc,assignment]
    sync_playwright = None  # type: ignore[assignment]
    PLAYWRIGHT_IMPORT_ERROR: Optional[str] = str(exc)
else:
    PLAYWRIGHT_IMPORT_ERROR = None


SCRIPT_PATH = Path(__file__).resolve()
ASSETS_ROOT = SCRIPT_PATH.parents[2]
WORKSTATION_DIRS = {
    "layered": ASSETS_ROOT / "layered_workstation",
    "topology": ASSETS_ROOT / "topology_workstation",
}
SAMPLES: Tuple[str, ...] = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
QUICK_SAMPLE = "HCC1395"
VIEWPORTS = {
    "desktop": {"width": 1440, "height": 1000},
    "mobile": {"width": 390, "height": 844},
}


def utc_now() -> str:
    """Return a compact, timezone-aware timestamp."""

    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def json_safe(value: Any) -> Any:
    """Convert common Playwright return values to JSON-safe values."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return str(value)


def add_check(
    run: Dict[str, Any],
    name: str,
    passed: bool,
    *,
    expected: Any,
    actual: Any,
    details: Optional[Dict[str, Any]] = None,
) -> None:
    """Append one deterministic assertion to a page run."""

    check: Dict[str, Any] = {
        "name": name,
        "status": "pass" if passed else "fail",
        "expected": json_safe(expected),
        "actual": json_safe(actual),
    }
    if details:
        check["details"] = json_safe(details)
    run["checks"].append(check)


def document_specs(mode: str) -> List[Dict[str, Any]]:
    """Build the fixed audit scope; do not infer coverage from index links."""

    samples: Sequence[str] = SAMPLES if mode == "all" else (QUICK_SAMPLE,)
    specs: List[Dict[str, Any]] = []
    for workstation, directory in WORKSTATION_DIRS.items():
        specs.append(
            {
                "workstation": workstation,
                "kind": "index",
                "sample": None,
                "path": directory / "index.html",
            }
        )
        specs.extend(
            {
                "workstation": workstation,
                "kind": "sample",
                "sample": sample,
                "path": directory / f"{sample}.html",
            }
            for sample in samples
        )
    return specs


def check_index_links(page: Page, run: Dict[str, Any], index_path: Path) -> None:
    """Verify every local index link and the seven expected sample links."""

    anchors = page.eval_on_selector_all(
        "a[href]",
        """elements => elements.map(a => ({
            href: a.getAttribute('href') || '',
            text: (a.innerText || a.getAttribute('aria-label') || '').trim()
        }))""",
    )
    local_links: List[Dict[str, Any]] = []
    missing_links: List[Dict[str, Any]] = []
    found_samples = set()
    base_url = index_path.resolve().as_uri()

    for anchor in anchors:
        href = str(anchor.get("href", "")).strip()
        if not href or href.startswith("#"):
            continue
        parsed = urlparse(urljoin(base_url, href))
        if parsed.scheme in {"http", "https", "mailto", "javascript", "data"}:
            continue
        if parsed.scheme not in {"", "file"}:
            continue
        target = Path(unquote(parsed.path)).resolve()
        record = {
            "href": href,
            "text": anchor.get("text", ""),
            "target": str(target),
            "exists": target.exists(),
        }
        local_links.append(record)
        if not target.exists():
            missing_links.append(record)
        if target.parent == index_path.parent and target.suffix.lower() == ".html":
            found_samples.add(target.stem)

    add_check(
        run,
        "index_local_links_exist",
        not missing_links,
        expected={"missing": 0},
        actual={"checked": len(local_links), "missing": len(missing_links)},
        details={"missing_links": missing_links[:20]},
    )
    missing_samples = sorted(set(SAMPLES) - found_samples)
    add_check(
        run,
        "index_expected_sample_links",
        not missing_samples,
        expected={"samples": list(SAMPLES)},
        actual={
            "found": sorted(found_samples.intersection(SAMPLES)),
            "missing": missing_samples,
        },
    )


def check_mobile_body_overflow(page: Page, run: Dict[str, Any]) -> None:
    """Reject page-level horizontal overflow while allowing local scrollers."""

    metrics = page.evaluate(
        """() => {
            const root = document.documentElement;
            const body = document.body;
            const width = window.innerWidth;
            const overflow = Math.max(root.scrollWidth, body ? body.scrollWidth : 0) - width;
            const offenders = Array.from(document.querySelectorAll('body *'))
                .filter(el => {
                    const style = getComputedStyle(el);
                    const rect = el.getBoundingClientRect();
                    return style.display !== 'none' && style.visibility !== 'hidden' &&
                        rect.width > 0 && rect.right > width + 1;
                })
                .slice(0, 8)
                .map(el => {
                    const rect = el.getBoundingClientRect();
                    return {
                        tag: el.tagName.toLowerCase(),
                        id: el.id || null,
                        class: typeof el.className === 'string' ? el.className : null,
                        right: Math.round(rect.right * 10) / 10,
                        width: Math.round(rect.width * 10) / 10
                    };
                });
            return {
                viewport_width: width,
                root_client_width: root.clientWidth,
                root_scroll_width: root.scrollWidth,
                body_scroll_width: body ? body.scrollWidth : null,
                overflow_px: Math.max(0, Math.round(overflow * 10) / 10),
                offenders
            };
        }"""
    )
    passed = float(metrics.get("overflow_px", 0)) <= 1
    add_check(
        run,
        "layered_index_mobile_body_overflow",
        passed,
        expected={"overflow_px_max": 1},
        actual=metrics,
    )


def check_layered_home_link(page: Page, run: Dict[str, Any]) -> None:
    """Require a visible, top-of-page route back to the sample index."""

    match = page.evaluate(
        """() => {
            const targetPath = new URL('index.html', location.href).pathname;
            const visible = el => {
                const rect = el.getBoundingClientRect();
                const style = getComputedStyle(el);
                return style.display !== 'none' && style.visibility !== 'hidden' &&
                    rect.width > 0 && rect.height > 0 && rect.top < Math.max(600, innerHeight);
            };
            const links = Array.from(document.querySelectorAll('a[href]'));
            const found = links.find(a => {
                let sameIndex = false;
                try { sameIndex = new URL(a.getAttribute('href'), location.href).pathname === targetPath; }
                catch (_) { return false; }
                const name = (a.innerText || a.getAttribute('aria-label') || a.title || '').trim();
                return sameIndex && visible(a) && /首頁|總覽|home|工作站/i.test(name);
            });
            if (!found) return null;
            const rect = found.getBoundingClientRect();
            return {
                text: (found.innerText || found.getAttribute('aria-label') || found.title || '').trim(),
                href: found.getAttribute('href'),
                top: Math.round(rect.top * 10) / 10,
                container: found.closest('nav,[aria-label*="breadcrumb" i],.breadcrumb,[class*="crumb" i],header')?.tagName.toLowerCase() || null
            };
        }"""
    )
    add_check(
        run,
        "layered_sample_breadcrumb_or_home",
        match is not None,
        expected="visible top-of-page link to layered index with a clear home/overview name",
        actual=match,
    )


def check_topology_archive_banner(page: Page, run: Dict[str, Any]) -> None:
    """Require direct legacy sample URLs to announce archive status and canonical UI."""

    banner = page.evaluate(
        """() => {
            const statusRe = /封存|已停用|deprecated|archived|archive|唯讀|legacy/i;
            const selectors = [
                '[class*="archive" i]', '[id*="archive" i]',
                '[class*="legacy" i]', '[id*="legacy" i]',
                '[class*="deprecat" i]', '[id*="deprecat" i]',
                '[class*="banner" i]', '[role="alert"]', '[data-status]'
            ].join(',');
            const visible = el => {
                const rect = el.getBoundingClientRect();
                const style = getComputedStyle(el);
                return style.display !== 'none' && style.visibility !== 'hidden' &&
                    rect.width > 0 && rect.height > 0 && rect.top < Math.max(700, innerHeight);
            };
            const candidates = Array.from(document.querySelectorAll(selectors));
            const found = candidates.find(el => visible(el) && statusRe.test(el.innerText || ''));
            if (!found) return null;
            const canonical = Array.from(found.querySelectorAll('a[href]')).find(a => {
                try {
                    return new URL(a.getAttribute('href'), location.href).pathname
                        .includes('/layered_workstation/index.html');
                } catch (_) { return false; }
            });
            const rect = found.getBoundingClientRect();
            return {
                tag: found.tagName.toLowerCase(),
                id: found.id || null,
                class: typeof found.className === 'string' ? found.className : null,
                text: (found.innerText || '').trim().slice(0, 500),
                top: Math.round(rect.top * 10) / 10,
                canonical_href: canonical ? canonical.getAttribute('href') : null
            };
        }"""
    )
    add_check(
        run,
        "topology_direct_sample_archive_banner",
        banner is not None,
        expected="visible archive/deprecated banner near the top of the direct sample page",
        actual=banner,
    )
    add_check(
        run,
        "topology_archive_banner_canonical_link",
        bool(banner and banner.get("canonical_href")),
        expected="archive banner links to layered_workstation/index.html",
        actual=banner.get("canonical_href") if banner else None,
    )


def choose_region_row(page: Page, offset: int = 0) -> Tuple[Locator, Dict[str, Any]]:
    """Choose a low-complexity region row to keep interaction smoke fast."""

    metadata = page.evaluate(
        """offset => {
            const rows = Array.from(document.querySelectorAll('#list .row'));
            let preferred = rows
                .map((el, index) => ({el, index, text: (el.innerText || '').trim()}))
                .filter(item => /(^|\D)2\s*sSNV/i.test(item.text));
            if (!preferred.length) {
                preferred = rows.map((el, index) => ({el, index, text: (el.innerText || '').trim()}));
            }
            const item = preferred[Math.min(offset, Math.max(0, preferred.length - 1))];
            if (!item) return null;
            const labelEl = item.el.querySelector('.mono,b');
            return {
                dom_index: item.index,
                data_i: item.el.getAttribute('data-i'),
                region: (labelEl ? labelEl.textContent : item.text.split(/\s+/)[0]).trim(),
                tag: item.el.tagName.toLowerCase(),
                role: item.el.getAttribute('role'),
                tabindex: item.el.getAttribute('tabindex'),
                tab_index_value: item.el.tabIndex
            };
        }""",
        offset,
    )
    if metadata is None:
        raise LookupError("No rendered #list .row element was found")
    return page.locator("#list .row").nth(int(metadata["dom_index"])), metadata


def region_selection_state(page: Page, metadata: Dict[str, Any]) -> Dict[str, Any]:
    """Read selected-row and detail evidence without assuming a tag name."""

    return page.evaluate(
        """meta => {
            const rows = Array.from(document.querySelectorAll('#list .row'));
            const row = rows.find(el => el.getAttribute('data-i') === String(meta.data_i));
            const detail = document.getElementById('detail');
            const detailHeading = detail ? detail.querySelector('h3') : null;
            const active = document.activeElement;
            return {
                row_exists: !!row,
                selected: !!row && (row.classList.contains('sel') || row.getAttribute('aria-selected') === 'true'),
                active: active === row,
                active_data_i: active ? active.getAttribute('data-i') : null,
                detail_heading_contains_region: !!detailHeading &&
                    (detailHeading.innerText || '').includes(meta.region),
                detail_text_head: detail ? (detail.innerText || '').trim().slice(0, 240) : null
            };
        }""",
        metadata,
    )


def check_detail_in_mobile_view(page: Page, run: Dict[str, Any]) -> None:
    """Require the selected detail heading to enter the mobile viewport."""

    visibility = page.evaluate(
        """() => {
            const detail = document.getElementById('detail');
            const target = detail ? (detail.querySelector('h3') || detail.firstElementChild || detail) : null;
            if (!target) return null;
            const rect = target.getBoundingClientRect();
            const visibleHeight = Math.max(0, Math.min(rect.bottom, innerHeight) - Math.max(rect.top, 0));
            return {
                target: target.tagName.toLowerCase(),
                viewport_height: innerHeight,
                top: Math.round(rect.top * 10) / 10,
                bottom: Math.round(rect.bottom * 10) / 10,
                visible_height: Math.round(visibleHeight * 10) / 10,
                in_view: rect.top < innerHeight - 32 && rect.bottom > 32 && visibleHeight >= 20
            };
        }"""
    )
    add_check(
        run,
        "mobile_region_selection_reveals_detail",
        bool(visibility and visibility.get("in_view")),
        expected="selected region detail heading has at least 20 px visible in the mobile viewport",
        actual=visibility,
    )


def check_region_interaction(page: Page, run: Dict[str, Any], viewport_name: str, timeout_ms: int) -> None:
    """Exercise click and Enter-key selection for one rendered region row."""

    try:
        page.wait_for_selector("#list .row", state="visible", timeout=timeout_ms)
        row_count = page.locator("#list .row").count()
    except Exception as exc:
        add_check(
            run,
            "region_rows_rendered",
            False,
            expected="at least one visible #list .row",
            actual={"error": str(exc)},
        )
        return

    add_check(
        run,
        "region_rows_rendered",
        row_count > 0,
        expected={"minimum": 1},
        actual={"count": row_count},
    )
    if row_count < 1:
        return

    click_row, click_meta = choose_region_row(page, 0)
    click_error: Optional[str] = None
    try:
        click_row.click(timeout=timeout_ms)
        page.wait_for_timeout(800 if viewport_name == "mobile" else 120)
        click_state = region_selection_state(page, click_meta)
    except Exception as exc:
        click_error = str(exc)
        click_state = {"error": click_error}
    click_passed = bool(
        not click_error
        and click_state.get("selected")
        and click_state.get("detail_heading_contains_region")
    )
    add_check(
        run,
        "region_row_click_selection",
        click_passed,
        expected="click marks the same row selected and renders its region in #detail",
        actual={"row": click_meta, "state": click_state},
    )

    if viewport_name == "mobile":
        check_detail_in_mobile_view(page, run)

    key_row, key_meta = choose_region_row(page, 1 if row_count > 1 else 0)
    key_error: Optional[str] = None
    try:
        key_row.focus(timeout=timeout_ms)
        focused_before_key = page.evaluate(
            "meta => document.activeElement?.getAttribute('data-i') === String(meta.data_i)",
            key_meta,
        )
        key_row.press("Enter", timeout=timeout_ms)
        page.wait_for_timeout(180)
        key_state = region_selection_state(page, key_meta)
        key_state["focused_before_key"] = focused_before_key
    except Exception as exc:
        key_error = str(exc)
        key_state = {"error": key_error, "focused_before_key": False}
    key_passed = bool(
        not key_error
        and key_state.get("focused_before_key")
        and key_state.get("selected")
        and key_state.get("detail_heading_contains_region")
    )
    add_check(
        run,
        "region_row_keyboard_selection",
        key_passed,
        expected="focusable row activates with Enter, becomes selected, and renders its region in #detail",
        actual={"row": key_meta, "state": key_state},
    )


def screenshot_path_for(
    screenshots_dir: Path,
    workstation: str,
    kind: str,
    sample: Optional[str],
    viewport_name: str,
) -> Path:
    """Return a stable, filesystem-safe screenshot path."""

    identity = sample if sample else kind
    safe_identity = re.sub(r"[^A-Za-z0-9_.-]+", "_", identity)
    return screenshots_dir / f"{workstation}__{safe_identity}__{viewport_name}.png"


def should_capture_screenshot(mode: str, run: Dict[str, Any]) -> bool:
    """Decide screenshot capture after all validation checks are known."""

    if mode == "all":
        return True
    if mode == "failures":
        return any(check["status"] != "pass" for check in run["checks"])
    return False


def audit_page(
    browser: Browser,
    spec: Dict[str, Any],
    viewport_name: str,
    timeout_ms: int,
    screenshot_mode: str,
    screenshots_dir: Optional[Path],
) -> Dict[str, Any]:
    """Run one page at one viewport and return a self-contained result."""

    started = time.monotonic()
    path = Path(spec["path"]).resolve()
    run: Dict[str, Any] = {
        "workstation": spec["workstation"],
        "kind": spec["kind"],
        "sample": spec["sample"],
        "viewport": viewport_name,
        "viewport_size": VIEWPORTS[viewport_name],
        "input_path": str(path),
        "url": path.as_uri(),
        "checks": [],
        "screenshot": None,
    }
    console_errors: List[Dict[str, Any]] = []
    page_errors: List[str] = []
    request_failures: List[Dict[str, Any]] = []
    context = browser.new_context(
        viewport=VIEWPORTS[viewport_name],
        color_scheme="light",
        locale="zh-TW",
        device_scale_factor=1,
    )
    page = context.new_page()
    page.set_default_timeout(timeout_ms)

    def on_console(message: Any) -> None:
        if message.type == "error":
            console_errors.append(
                {
                    "text": message.text,
                    "location": json_safe(message.location),
                }
            )

    def on_page_error(error: Any) -> None:
        page_errors.append(str(error))

    def on_request_failed(request: Any) -> None:
        request_failures.append(
            {
                "url": request.url,
                "failure": json_safe(request.failure),
            }
        )

    page.on("console", on_console)
    page.on("pageerror", on_page_error)
    page.on("requestfailed", on_request_failed)

    navigation_ok = False
    try:
        if not path.is_file():
            raise FileNotFoundError(f"Input HTML does not exist: {path}")
        page.goto(path.as_uri(), wait_until="load", timeout=timeout_ms)
        page.wait_for_load_state("networkidle", timeout=timeout_ms)
        navigation_ok = True
        add_check(
            run,
            "document_load",
            True,
            expected="file:// document reaches load and networkidle",
            actual={"ready_state": page.evaluate("document.readyState")},
        )
    except Exception as exc:
        add_check(
            run,
            "document_load",
            False,
            expected="file:// document reaches load and networkidle",
            actual={"error": str(exc)},
        )

    if navigation_ok:
        try:
            title = page.title().strip()
            add_check(
                run,
                "document_title_nonempty",
                bool(title),
                expected="non-empty document title",
                actual=title,
            )
            if spec["kind"] == "sample":
                sample = str(spec["sample"])
                heading = page.locator("h1").first.inner_text(timeout=timeout_ms).strip()
                add_check(
                    run,
                    "title_sample_match",
                    sample.lower() in title.lower(),
                    expected={"sample_in_title": sample},
                    actual=title,
                )
                add_check(
                    run,
                    "heading_sample_match",
                    sample.lower() in heading.lower(),
                    expected={"sample_in_h1": sample},
                    actual=heading,
                )
                if spec["workstation"] == "layered":
                    check_layered_home_link(page, run)
                else:
                    check_topology_archive_banner(page, run)
                check_region_interaction(page, run, viewport_name, timeout_ms)
            else:
                check_index_links(page, run, path)
                if spec["workstation"] == "layered" and viewport_name == "mobile":
                    check_mobile_body_overflow(page, run)
        except Exception as exc:
            add_check(
                run,
                "validator_page_assertions",
                False,
                expected="page-specific assertions complete without validator exception",
                actual={"error": str(exc)},
            )

    page.wait_for_timeout(100)
    add_check(
        run,
        "console_errors",
        not console_errors,
        expected={"count": 0},
        actual={"count": len(console_errors), "errors": console_errors[:20]},
    )
    add_check(
        run,
        "page_errors",
        not page_errors,
        expected={"count": 0},
        actual={"count": len(page_errors), "errors": page_errors[:20]},
    )
    add_check(
        run,
        "request_failures",
        not request_failures,
        expected={"count": 0},
        actual={"count": len(request_failures), "failures": request_failures[:20]},
    )

    if screenshots_dir and should_capture_screenshot(screenshot_mode, run):
        target = screenshot_path_for(
            screenshots_dir,
            spec["workstation"],
            spec["kind"],
            spec["sample"],
            viewport_name,
        )
        try:
            page.screenshot(path=str(target), full_page=False)
            run["screenshot"] = str(target)
        except Exception as exc:
            add_check(
                run,
                "screenshot_capture",
                False,
                expected={"path": str(target)},
                actual={"error": str(exc)},
            )

    run["duration_seconds"] = round(time.monotonic() - started, 3)
    context.close()
    return run


def summarize(report: Dict[str, Any]) -> None:
    """Attach stable aggregate counts and failing-check coordinates."""

    runs = report["runs"]
    checks = [check for run in runs for check in run["checks"]]
    failed = [check for check in checks if check["status"] != "pass"]
    failing_coordinates: List[Dict[str, Any]] = []
    for run in runs:
        names = [check["name"] for check in run["checks"] if check["status"] != "pass"]
        if names:
            failing_coordinates.append(
                {
                    "workstation": run["workstation"],
                    "kind": run["kind"],
                    "sample": run["sample"],
                    "viewport": run["viewport"],
                    "failed_checks": names,
                }
            )
    documents = {
        (run["workstation"], run["kind"], run["sample"])
        for run in runs
    }
    report["summary"] = {
        "status": "pass" if not failed else "fail",
        "documents_tested": len(documents),
        "page_runs": len(runs),
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failed),
        "checks_failed": len(failed),
        "page_runs_failed": len(failing_coordinates),
        "failing_coordinates": failing_coordinates,
    }


def fatal_report(mode: str, message: str, started_at: str) -> Dict[str, Any]:
    """Create structured output for a validator-level failure."""

    return {
        "schema_version": "1.0",
        "tool": SCRIPT_PATH.name,
        "mode": mode,
        "started_at": started_at,
        "finished_at": utc_now(),
        "scope": {
            "partial": mode != "all",
            "documents_expected_full": 16,
            "input_roots": {key: str(value) for key, value in WORKSTATION_DIRS.items()},
        },
        "browser": None,
        "runs": [],
        "summary": {
            "status": "error",
            "documents_tested": 0,
            "page_runs": 0,
            "checks_total": 0,
            "checks_passed": 0,
            "checks_failed": 0,
            "page_runs_failed": 0,
            "failing_coordinates": [],
        },
        "fatal_errors": [message],
        "exit_code": 2,
    }


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI without touching browser state."""

    parser = argparse.ArgumentParser(
        description=(
            "Use headless Chromium to smoke-test generated layered/topology workstation "
            "HTML over file:// and emit one JSON report to stdout."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Modes: --quick audits both indices plus HCC1395 in desktop/mobile; "
            "--all audits both indices and all 7+7 sample pages in desktop/mobile. "
            "Exit codes: 0 pass, 1 validation failure, 2 runner/configuration failure."
        ),
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--quick",
        action="store_true",
        help="Run the representative 4-document / 8-page-run partial smoke.",
    )
    mode.add_argument(
        "--all",
        action="store_true",
        help="Run the complete 16-document / 32-page-run smoke.",
    )
    parser.add_argument(
        "--screenshots-dir",
        type=Path,
        help=(
            "Write viewport PNG screenshots here. Supplying this option enables "
            "--screenshot-mode=all unless that mode is set explicitly."
        ),
    )
    parser.add_argument(
        "--screenshot-mode",
        choices=("off", "failures", "all"),
        default=None,
        help="Capture no screenshots, only failing page runs, or every page run.",
    )
    parser.add_argument(
        "--timeout-ms",
        type=int,
        default=60_000,
        help="Per-navigation, selector, and interaction timeout in milliseconds.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the selected smoke suite and emit exactly one JSON document."""

    parser = build_parser()
    args = parser.parse_args(argv)
    mode = "all" if args.all else "quick"
    started_at = utc_now()
    screenshot_mode = args.screenshot_mode or ("all" if args.screenshots_dir else "off")

    if args.timeout_ms <= 0:
        report = fatal_report(mode, "--timeout-ms must be greater than zero", started_at)
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 2
    if screenshot_mode != "off" and not args.screenshots_dir:
        report = fatal_report(
            mode,
            "--screenshot-mode requires --screenshots-dir unless the mode is off",
            started_at,
        )
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 2
    if PLAYWRIGHT_IMPORT_ERROR:
        report = fatal_report(
            mode,
            f"Python Playwright import failed: {PLAYWRIGHT_IMPORT_ERROR}",
            started_at,
        )
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 2

    screenshots_dir: Optional[Path] = None
    if args.screenshots_dir:
        try:
            screenshots_dir = args.screenshots_dir.expanduser().resolve()
            screenshots_dir.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            report = fatal_report(mode, f"Cannot create screenshots directory: {exc}", started_at)
            print(json.dumps(report, ensure_ascii=False, indent=2))
            return 2

    specs = document_specs(mode)
    report: Dict[str, Any] = {
        "schema_version": "1.0",
        "tool": SCRIPT_PATH.name,
        "mode": mode,
        "started_at": started_at,
        "finished_at": None,
        "scope": {
            "partial": mode != "all",
            "documents_expected_full": 16,
            "documents_selected": len(specs),
            "page_runs_selected": len(specs) * len(VIEWPORTS),
            "samples_expected": list(SAMPLES),
            "input_roots": {key: str(value) for key, value in WORKSTATION_DIRS.items()},
            "viewports": VIEWPORTS,
        },
        "screenshot": {
            "mode": screenshot_mode,
            "output_dir": str(screenshots_dir) if screenshots_dir else None,
            "full_page": False,
        },
        "browser": None,
        "runs": [],
        "fatal_errors": [],
    }

    playwright: Optional[Playwright] = None
    browser: Optional[Browser] = None
    try:
        playwright = sync_playwright().start()  # type: ignore[union-attr]
        browser = playwright.chromium.launch(
            headless=True,
            args=["--allow-file-access-from-files"],
        )
        report["browser"] = {
            "engine": "chromium",
            "version": browser.version,
            "headless": True,
        }
        for spec in specs:
            for viewport_name in VIEWPORTS:
                report["runs"].append(
                    audit_page(
                        browser,
                        spec,
                        viewport_name,
                        args.timeout_ms,
                        screenshot_mode,
                        screenshots_dir,
                    )
                )
    except Exception as exc:
        report["fatal_errors"].append(str(exc))
    finally:
        if browser is not None:
            try:
                browser.close()
            except Exception as exc:
                report["fatal_errors"].append(f"Browser close failed: {exc}")
        if playwright is not None:
            try:
                playwright.stop()
            except Exception as exc:
                report["fatal_errors"].append(f"Playwright stop failed: {exc}")

    summarize(report)
    if report["fatal_errors"]:
        report["summary"]["status"] = "error"
        exit_code = 2
    elif report["summary"]["checks_failed"]:
        exit_code = 1
    else:
        exit_code = 0
    report["finished_at"] = utc_now()
    report["exit_code"] = exit_code
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
